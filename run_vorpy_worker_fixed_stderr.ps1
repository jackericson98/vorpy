param(
    [string]$DataRoot = "X:\Jack\ProtSA\p53tet_1pet\solvated_new",
    [string]$VorpyRoot = (Get-Location).Path,
    [int]$StartIndex = 0,
    [int]$EndIndex = 1999,
    [int]$PollSeconds = 60,
    [int]$MinFileAgeSeconds = 10
)

$ErrorActionPreference = "Stop"

# Expected snapshots: 100, 200, ..., 1000 ps
$ExpectedFrames = 100..1000 | Where-Object { $_ % 100 -eq 0 } | ForEach-Object {
    "frame_{0:D4}ps.pdb" -f $_
}

$TotalSystems = $EndIndex - $StartIndex + 1
$TotalFrames = $TotalSystems * $ExpectedFrames.Count

$WorkerRoot = Join-Path $DataRoot "_vorpy_worker"
$LogDir = Join-Path $WorkerRoot "logs"
$StatusDir = Join-Path $WorkerRoot "status"

New-Item -ItemType Directory -Force -Path $WorkerRoot, $LogDir, $StatusDir | Out-Null

$MasterLog = Join-Path $LogDir "worker.log"
$SummaryCsv = Join-Path $WorkerRoot "status.csv"

function Write-WorkerLog {
    param([string]$Message)

    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $line = "[$stamp] $Message"
    Write-Host $line
    Add-Content -LiteralPath $MasterLog -Value $line
}

function Get-SystemPaths {
    param([int]$Index)

    $name = "n$Index"
    $systemDir = Join-Path $DataRoot $name
    $frameDir = Join-Path $systemDir "analysis_frames"
    $systemStatusDir = Join-Path $StatusDir $name

    [PSCustomObject]@{
        Name = $name
        SystemDir = $systemDir
        FrameDir = $frameDir
        StatusDir = $systemStatusDir
    }
}

function Test-SystemReady {
    param($Paths)

    if (-not (Test-Path -LiteralPath $Paths.FrameDir -PathType Container)) {
        return $false
    }

    foreach ($frameName in $ExpectedFrames) {
        $framePath = Join-Path $Paths.FrameDir $frameName
        if (-not (Test-Path -LiteralPath $framePath -PathType Leaf)) {
            return $false
        }

        # Avoid grabbing a file while another process is still writing it.
        try {
            $item = Get-Item -LiteralPath $framePath
            if ($item.Length -le 0) {
                return $false
            }

            # Treat very recently modified files as still being produced/copied.
            $ageSeconds = ((Get-Date) - $item.LastWriteTime).TotalSeconds
            if ($ageSeconds -lt $MinFileAgeSeconds) {
                return $false
            }
        }
        catch {
            return $false
        }
    }

    return $true
}

function Get-DoneMarker {
    param($Paths, [string]$FrameName)

    $base = [System.IO.Path]::GetFileNameWithoutExtension($FrameName)
    return Join-Path $Paths.StatusDir "$base.done"
}

function Test-FrameDone {
    param($Paths, [string]$FrameName)

    $marker = Get-DoneMarker -Paths $Paths -FrameName $FrameName
    return Test-Path -LiteralPath $marker -PathType Leaf
}

function Get-CompletedFrameCount {
    $count = 0

    for ($i = $StartIndex; $i -le $EndIndex; $i++) {
        $paths = Get-SystemPaths -Index $i

        foreach ($frameName in $ExpectedFrames) {
            if (Test-FrameDone -Paths $paths -FrameName $frameName) {
                $count++
            }
        }
    }

    return $count
}

function Write-StatusCsv {
    $rows = foreach ($i in $StartIndex..$EndIndex) {
        $paths = Get-SystemPaths -Index $i
        $ready = Test-SystemReady -Paths $paths

        $done = 0
        foreach ($frameName in $ExpectedFrames) {
            if (Test-FrameDone -Paths $paths -FrameName $frameName) {
                $done++
            }
        }

        [PSCustomObject]@{
            system = $paths.Name
            ready = $ready
            frames_done = $done
            frames_total = $ExpectedFrames.Count
            complete = ($done -eq $ExpectedFrames.Count)
            frame_directory = $paths.FrameDir
        }
    }

    $rows | Export-Csv -LiteralPath $SummaryCsv -NoTypeInformation
}

function Invoke-VorpyFrame {
    param(
        $Paths,
        [string]$FrameName,
        [int]$SystemIndex,
        [int]$CompletedBefore
    )

    New-Item -ItemType Directory -Force -Path $Paths.StatusDir | Out-Null

    $framePath = Join-Path $Paths.FrameDir $FrameName
    $frameBase = [System.IO.Path]::GetFileNameWithoutExtension($FrameName)
    $doneMarker = Get-DoneMarker -Paths $Paths -FrameName $FrameName
    $failMarker = Join-Path $Paths.StatusDir "$frameBase.failed"
    $frameLog = Join-Path $LogDir ("{0}_{1}.log" -f $Paths.Name, $frameBase)

    if (Test-Path -LiteralPath $doneMarker) {
        return $true
    }

    if (Test-Path -LiteralPath $failMarker) {
        Remove-Item -LiteralPath $failMarker -Force
    }

    $overallNumber = $CompletedBefore + 1
    $overallPercent = if ($TotalFrames -gt 0) {
        [math]::Round(($CompletedBefore / $TotalFrames) * 100, 2)
    } else {
        100
    }

    Write-WorkerLog ""
    Write-WorkerLog ("=" * 78)
    Write-WorkerLog "Starting $($Paths.Name) / $FrameName"
    Write-WorkerLog "Overall completed: $CompletedBefore / $TotalFrames frames ($overallPercent%)"
    Write-WorkerLog "Input:  $framePath"
    Write-WorkerLog "Output: $($Paths.FrameDir)"
    Write-WorkerLog ("=" * 78)

    Push-Location $VorpyRoot
    try {
        # Windows PowerShell 5.1 can promote a native program's stderr output
        # (including an ordinary Python traceback) into NativeCommandError when
        # $ErrorActionPreference is "Stop". Temporarily allow native stderr so
        # the complete VorPy output is displayed/logged, and use LASTEXITCODE
        # as the authoritative success/failure signal.
        $savedErrorActionPreference = $ErrorActionPreference
        $ErrorActionPreference = "Continue"

        try {
            # User-requested VorPy invocation:
            # python vorpy [FILE] -s mv 5 -e small and info -e dir [FILE DIR]
            #
            # & preserves argument boundaries correctly even when paths contain spaces.
            & python vorpy $framePath -s mv 5 -e small and info -e dir $Paths.FrameDir 2>&1 |
                Tee-Object -FilePath $frameLog

            $exitCode = $LASTEXITCODE
        }
        finally {
            $ErrorActionPreference = $savedErrorActionPreference
        }
    }
    catch {
        $_ | Out-String | Tee-Object -FilePath $frameLog -Append | Out-Host
        $exitCode = 1
    }
    finally {
        Pop-Location
    }

    if ($exitCode -eq 0) {
        $finished = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
        @(
            "system=$($Paths.Name)"
            "frame=$FrameName"
            "input=$framePath"
            "output_dir=$($Paths.FrameDir)"
            "finished=$finished"
            "exit_code=0"
        ) | Set-Content -LiteralPath $doneMarker

        if (Test-Path -LiteralPath $failMarker) {
            Remove-Item -LiteralPath $failMarker -Force
        }

        $completedNow = $CompletedBefore + 1
        $pctNow = [math]::Round(($completedNow / $TotalFrames) * 100, 2)
        Write-WorkerLog "SUCCESS: $($Paths.Name) / $FrameName"
        Write-WorkerLog "Overall completed: $completedNow / $TotalFrames frames ($pctNow%)"
        return $true
    }
    else {
        @(
            "system=$($Paths.Name)"
            "frame=$FrameName"
            "input=$framePath"
            "failed=$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
            "exit_code=$exitCode"
            "log=$frameLog"
        ) | Set-Content -LiteralPath $failMarker

        Write-WorkerLog "FAILED: $($Paths.Name) / $FrameName (exit code $exitCode)"
        Write-WorkerLog "Failure log: $frameLog"
        return $false
    }
}


# ---------------------------------------------------------------------------
# Startup checks
# ---------------------------------------------------------------------------

if (-not (Test-Path -LiteralPath $DataRoot -PathType Container)) {
    throw "Data root does not exist or is not accessible: $DataRoot"
}

if (-not (Test-Path -LiteralPath $VorpyRoot -PathType Container)) {
    throw "VorPy root does not exist: $VorpyRoot"
}

$vorpyEntry = Join-Path $VorpyRoot "vorpy"
if (-not (Test-Path -LiteralPath $vorpyEntry)) {
    Write-WorkerLog "WARNING: '$vorpyEntry' was not found."
    Write-WorkerLog "Make sure -VorpyRoot points to the repository directory where 'python vorpy ...' works."
}

try {
    $pythonVersion = & python --version 2>&1
    Write-WorkerLog "Python: $pythonVersion"
}
catch {
    throw "Could not execute 'python'. Activate the correct environment before starting this worker."
}

Write-WorkerLog "VorPy worker starting."
Write-WorkerLog "Data root:   $DataRoot"
Write-WorkerLog "VorPy root:  $VorpyRoot"
Write-WorkerLog "Range:       n$StartIndex through n$EndIndex"
Write-WorkerLog "Frames/run:  $($ExpectedFrames.Count)"
Write-WorkerLog "Total work:  $TotalFrames frames"
Write-WorkerLog "Poll period: $PollSeconds seconds"
Write-WorkerLog "Minimum frame age: $MinFileAgeSeconds seconds"
Write-WorkerLog "Status CSV:  $SummaryCsv"
Write-WorkerLog "Master log:  $MasterLog"

# ---------------------------------------------------------------------------
# Worker loop
# ---------------------------------------------------------------------------

while ($true) {
    $madeProgress = $false
    $readySystems = 0
    $waitingSystems = 0
    $completeSystems = 0
    $scanCompletedFrames = Get-CompletedFrameCount

    for ($i = $StartIndex; $i -le $EndIndex; $i++) {
        $paths = Get-SystemPaths -Index $i

        $doneCount = 0
        foreach ($frameName in $ExpectedFrames) {
            if (Test-FrameDone -Paths $paths -FrameName $frameName) {
                $doneCount++
            }
        }

        if ($doneCount -eq $ExpectedFrames.Count) {
            $completeSystems++
            continue
        }

        if (-not (Test-SystemReady -Paths $paths)) {
            $waitingSystems++
            continue
        }

        $readySystems++
        Write-WorkerLog "READY: $($paths.Name) has all $($ExpectedFrames.Count) analysis frames."

        foreach ($frameName in $ExpectedFrames) {
            if (Test-FrameDone -Paths $paths -FrameName $frameName) {
                continue
            }

            $success = Invoke-VorpyFrame `
                -Paths $paths `
                -FrameName $frameName `
                -SystemIndex $i `
                -CompletedBefore $scanCompletedFrames

            if ($success) {
                $madeProgress = $true
                $scanCompletedFrames++
            }
            else {
                # Do not hammer a broken frame repeatedly in the same scan.
                # It will be retried on the next scan.
                Write-WorkerLog "Moving on after failure; this frame will be retried during the next scan."
            }
        }
    }

    Write-StatusCsv

    $completedFrames = $scanCompletedFrames
    $completedPct = [math]::Round(($completedFrames / $TotalFrames) * 100, 2)

    Write-WorkerLog ""
    Write-WorkerLog "SCAN COMPLETE"
    Write-WorkerLog "Frames complete:  $completedFrames / $TotalFrames ($completedPct%)"
    Write-WorkerLog "Systems complete: $completeSystems / $TotalSystems"
    Write-WorkerLog "Systems ready:    $readySystems"
    Write-WorkerLog "Systems waiting:  $waitingSystems"

    if ($completedFrames -ge $TotalFrames) {
        Write-WorkerLog "ALL VORPY RUNS COMPLETE. Worker exiting."
        break
    }

    Write-WorkerLog "Rescanning in $PollSeconds seconds. Press Ctrl+C to stop safely."
    Start-Sleep -Seconds $PollSeconds
}
