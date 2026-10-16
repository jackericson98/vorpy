[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string] $FilePath,
    [string] $Arguments = '',
    [Parameter(Mandatory = $true)]
    [string] $OutputRoot,
    [Parameter(Mandatory = $true)]
    [string] $ArchivePath,
    [string] $LogDirectory = '',
    [int] $MaxWallSeconds = 1800,
    [int] $MaxCpuSeconds = 1800,
    [long] $MaxPeakRssBytes = 4294967296,
    [long] $MaxTotalOutputBytes = 1073741824,
    [int] $MaxArchiveMembers = 100000,
    [int] $PollMilliseconds = 1000,
    [switch] $RequireArchive
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Resolve-AbsolutePath([string] $PathValue) {
    if ([IO.Path]::IsPathRooted($PathValue)) {
        return [IO.Path]::GetFullPath($PathValue)
    }
    return [IO.Path]::GetFullPath((Join-Path (Get-Location) $PathValue))
}

function Get-OutputBytes([string] $Root) {
    try {
        try {
            $rootItem = Get-Item -LiteralPath $Root -ErrorAction Stop
        } catch {
            $category = $_.CategoryInfo.Category
            if ($category -eq 'ObjectNotFound' -or
                $_.Exception -is [IO.FileNotFoundException] -or
                $_.Exception -is [IO.DirectoryNotFoundException]) {
                return [long] 0
            }
            throw
        }
        if (-not $rootItem.PSIsContainer) {
            throw "Output root is not a directory: $Root"
        }

        [long] $sum = 0
        $files = @(Get-ChildItem -LiteralPath $Root -File -Recurse -ErrorAction Stop)
        foreach ($file in $files) {
            [long] $length = $file.Length
            if ($length -lt 0 -or $sum -gt ([long]::MaxValue - $length)) {
                throw "Output size exceeds Int64 capacity at: $($file.FullName)"
            }
            $sum += $length
        }
        return $sum
    } catch {
        throw "Output-size scan failed for '$Root': $($_.Exception.Message)"
    }
}

function Get-ArchiveMemberCount([string] $PathValue) {
    if (-not (Test-Path -LiteralPath $PathValue -PathType Leaf)) {
        return $null
    }
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $archive = $null
    try {
        $archive = [IO.Compression.ZipFile]::OpenRead($PathValue)
        return [int] $archive.Entries.Count
    } catch [IO.InvalidDataException] {
        return $null
    } catch [IO.IOException] {
        return $null
    } finally {
        if ($null -ne $archive) {
            $archive.Dispose()
        }
    }
}

function Get-ProcessTreeIds([int] $RootId) {
    $snapshot = @(Get-CimInstance -ClassName Win32_Process -ErrorAction Stop)
    $seen = New-Object 'System.Collections.Generic.HashSet[int]'
    $queue = New-Object 'System.Collections.Generic.Queue[int]'
    $ordered = New-Object 'System.Collections.Generic.List[int]'
    $null = $seen.Add($RootId)
    $queue.Enqueue($RootId)
    $ordered.Add($RootId)
    while ($queue.Count -gt 0) {
        $parentId = $queue.Dequeue()
        foreach ($candidate in $snapshot) {
            if ([int]$candidate.ParentProcessId -ne $parentId) {
                continue
            }
            $childId = [int]$candidate.ProcessId
            if ($seen.Add($childId)) {
                $queue.Enqueue($childId)
                $ordered.Add($childId)
            }
        }
    }
    return $ordered.ToArray()
}

function Stop-ProcessTree([int] $RootId) {
    $treeIds = @($RootId)
    $treeError = $null
    try {
        $treeIds = @(Get-ProcessTreeIds $RootId)
    } catch {
        $treeError = $_.Exception.Message
    }

    $taskkillExit = $null
    $taskkillOutput = @()
    try {
        $taskkillOutput = @(& taskkill.exe /PID $RootId /T /F 2>&1)
        $taskkillExit = $LASTEXITCODE
    } catch {
        $taskkillOutput = @($_.Exception.Message)
    }

    if ($null -eq $taskkillExit -or $taskkillExit -ne 0) {
        foreach ($processId in ($treeIds | Sort-Object -Descending)) {
            try {
                Stop-Process -Id $processId -Force -ErrorAction SilentlyContinue
            } catch {
                # The process may already have exited after taskkill's tree walk.
            }
        }
    }

    return [PSCustomObject]@{
        tree_ids = $treeIds
        tree_error = $treeError
        method = 'taskkill_tree_with_stop_process_fallback'
        taskkill_exit_code = $taskkillExit
        taskkill_output = (($taskkillOutput | ForEach-Object { $_.ToString() }) -join ' ')
    }
}

function Write-Summary([hashtable] $Summary, [string] $SummaryPath) {
    try {
        $json = [PSCustomObject]$Summary | ConvertTo-Json -Compress -Depth 8
        [IO.File]::WriteAllText($SummaryPath, $json)
    } catch {
        $Summary['summary_write_error'] = $_.Exception.Message
        $json = [PSCustomObject]$Summary | ConvertTo-Json -Compress -Depth 8
    }
    Write-Output $json
}

$resolvedFile = Resolve-AbsolutePath $FilePath
$resolvedOutput = Resolve-AbsolutePath $OutputRoot
$resolvedArchive = Resolve-AbsolutePath $ArchivePath
$resolvedTimingSidecar = $resolvedArchive + '.timing.json'
$resolvedLogDirectory = if ([string]::IsNullOrWhiteSpace($LogDirectory)) {
    $resolvedOutput + '.watchdog'
} else {
    Resolve-AbsolutePath $LogDirectory
}
$resolvedStdout = Join-Path $resolvedLogDirectory 'stdout.log'
$resolvedStderr = Join-Path $resolvedLogDirectory 'stderr.log'
$resolvedSummary = Join-Path $resolvedLogDirectory 'summary.json'

if ([string]::IsNullOrWhiteSpace($Arguments) -or
    $Arguments -notmatch '(^|\s)--fail-if-exists(?=\s|$)') {
    throw 'Watchdog requires the approved --fail-if-exists reservation flag.'
}

if (-not (Test-Path -LiteralPath $resolvedFile -PathType Leaf)) {
    throw "Watchdog executable does not exist: $resolvedFile"
}
if (Test-Path -LiteralPath $resolvedOutput) {
    throw "Refusing to overwrite an existing output root: $resolvedOutput"
}
if (Test-Path -LiteralPath $resolvedArchive) {
    throw "Refusing to overwrite an existing archive: $resolvedArchive"
}
if (Test-Path -LiteralPath $resolvedTimingSidecar) {
    throw "Refusing to overwrite an existing archive timing sidecar: $resolvedTimingSidecar"
}
if (Test-Path -LiteralPath $resolvedLogDirectory) {
    throw "Refusing to overwrite existing watchdog log directory: $resolvedLogDirectory"
}
if ($MaxWallSeconds -le 0 -or $MaxCpuSeconds -le 0 -or $MaxPeakRssBytes -le 0 -or
    $MaxTotalOutputBytes -le 0 -or $MaxArchiveMembers -le 0 -or $PollMilliseconds -le 0) {
    throw 'All watchdog limits must be positive.'
}

$startedAt = [DateTime]::UtcNow
$started = [Diagnostics.Stopwatch]::StartNew()
[IO.Directory]::CreateDirectory($resolvedLogDirectory) | Out-Null
$processStartError = $null
$process = $null
$launcherPath = if ([string]::IsNullOrWhiteSpace($env:ComSpec)) { 'cmd.exe' } else { $env:ComSpec }
$launcherCommand = '""' + $resolvedFile + '" ' + $Arguments +
    ' > "' + $resolvedStdout + '" 2> "' + $resolvedStderr +
    '" & exit /b !ERRORLEVEL!"'
$launcherArguments = @('/d', '/v:on', '/s', '/c', $launcherCommand)
try {
    $process = Start-Process -FilePath $launcherPath -ArgumentList $launcherArguments `
        -WorkingDirectory (Get-Location).Path -PassThru
} catch {
    $processStartError = $_.Exception.Message
}
if ($null -eq $process) {
    $started.Stop()
    $summary = @{
        status = 'FAILED_LAUNCH'
        failure = 'process_launch_failed'
        trigger = $null
        exit_code = $null
        launch_error = $processStartError
        started_at_utc = $startedAt.ToString('o')
        ended_at_utc = [DateTime]::UtcNow.ToString('o')
        wall_seconds = $started.Elapsed.TotalSeconds
        stdout_path = $resolvedStdout
        stderr_path = $resolvedStderr
        summary_path = $resolvedSummary
        output_root = $resolvedOutput
        archive_path = $resolvedArchive
        launcher_path = $launcherPath
        launcher_arguments = $launcherArguments
        strict_kernel_limits = $false
    }
    Write-Summary $summary $resolvedSummary
    exit 11
}
$trigger = $null
$termination = $null
$treeIds = @($process.Id)
$treeError = $null
$lastRss = [long] 0
$lastOutput = [long] 0
$lastCpu = [double] 0
$lastMembers = $null
$outputScanError = $null
$processMetricsError = $null

try {
    while (-not $process.HasExited) {
        $process.Refresh()
        try {
            $treeIds = @(Get-ProcessTreeIds $process.Id)
            $treeError = $null
        } catch {
            $treeError = $_.Exception.Message
        }
        try {
            $cpuTime = $process.TotalProcessorTime
            if ($null -eq $cpuTime) {
                throw 'Process CPU time is unavailable.'
            }
            $lastCpu = [double] $cpuTime.TotalSeconds
            $lastRss = [long] $process.WorkingSet64
        } catch {
            $processMetricsError = $_.Exception.Message
            $lastCpu = [double] 0
            $lastRss = [long] 0
        }
        try {
            $lastOutput = Get-OutputBytes $resolvedOutput
        } catch {
            $outputScanError = $_.Exception.Message
            $termination = Stop-ProcessTree $process.Id
            break
        }
        $lastMembers = Get-ArchiveMemberCount $resolvedArchive

        if ($started.Elapsed.TotalSeconds -ge $MaxWallSeconds) {
            $trigger = 'wall_seconds'
        } elseif ($lastCpu -ge $MaxCpuSeconds) {
            $trigger = 'cpu_seconds'
        } elseif ($lastRss -ge $MaxPeakRssBytes) {
            $trigger = 'peak_rss_bytes'
        } elseif ($lastOutput -ge $MaxTotalOutputBytes) {
            $trigger = 'total_output_bytes'
        } elseif ($null -ne $lastMembers -and $lastMembers -gt $MaxArchiveMembers) {
            $trigger = 'archive_member_count'
        }

        if ($null -ne $trigger) {
            $termination = Stop-ProcessTree $process.Id
            try {
                $process.WaitForExit(5000)
            } catch {
                # The process object may already be invalid after taskkill.
            }
            break
        }
        Start-Sleep -Milliseconds $PollMilliseconds
    }
    if ($process.HasExited) {
        try {
            $process.WaitForExit(5000) | Out-Null
        } catch {
            $processMetricsError = $_.Exception.Message
        }
        $process.Refresh()
        try {
            $cpuTime = $process.TotalProcessorTime
            if ($null -eq $cpuTime) {
                throw 'Process CPU time is unavailable.'
            }
            $lastCpu = [double] $cpuTime.TotalSeconds
            $lastRss = [long] $process.WorkingSet64
        } catch {
            $processMetricsError = $_.Exception.Message
        }
    }
} finally {
    $started.Stop()
}

try {
    $lastOutput = Get-OutputBytes $resolvedOutput
} catch {
    if ($null -eq $outputScanError) {
        $outputScanError = $_.Exception.Message
    }
}
$lastMembers = Get-ArchiveMemberCount $resolvedArchive
$streamCaptureError = $null
$stdoutBytes = $null
$stderrBytes = $null
try {
    $stdoutBytes = [long](Get-Item -LiteralPath $resolvedStdout -ErrorAction Stop).Length
    $stderrBytes = [long](Get-Item -LiteralPath $resolvedStderr -ErrorAction Stop).Length
} catch {
    $streamCaptureError = $_.Exception.Message
}
$exitCode = $null
if ($process.HasExited) {
    try {
        $process.WaitForExit(5000) | Out-Null
        $process.Refresh()
        $exitCode = [int] $process.ExitCode
    } catch {
        $processMetricsError = $_.Exception.Message
    }
}
$status = 'PROCESS_EXITED'
$failure = $null
if ($null -ne $outputScanError) {
    $status = 'FAILED_OUTPUT_SCAN'
    $failure = $outputScanError
} elseif ($null -ne $trigger) {
    $status = 'ABORTED_LIMIT'
    $failure = "limit:$trigger"
} elseif ($null -ne $exitCode -and $exitCode -ne 0) {
    $status = 'PROCESS_FAILED'
    $failure = "process_exit:$exitCode"
} elseif ($RequireArchive -and $null -eq $lastMembers) {
    $status = 'FAILED_ARCHIVE_MISSING_OR_UNREADABLE'
    $failure = 'required_archive_missing_or_unreadable'
} elseif ($null -ne $lastMembers -and $lastMembers -gt $MaxArchiveMembers) {
    $status = 'FAILED_ARCHIVE_MEMBER_LIMIT'
    $failure = "archive_members:$lastMembers>$MaxArchiveMembers"
}
$remainingTreeIds = @()
try {
    if (-not $process.HasExited) {
        $remainingTreeIds = @(Get-ProcessTreeIds $process.Id)
    }
} catch {
    if ($null -eq $treeError) {
        $treeError = $_.Exception.Message
    }
}
$summary = @{
    status = $status
    failure = $failure
    trigger = $trigger
    exit_code = $exitCode
    started_at_utc = $startedAt.ToString('o')
    ended_at_utc = [DateTime]::UtcNow.ToString('o')
    wall_seconds = $started.Elapsed.TotalSeconds
    cpu_seconds = $lastCpu
    peak_rss_bytes = $lastRss
    total_output_bytes = $lastOutput
    archive_member_count = $lastMembers
    output_scan_error = $outputScanError
    launch_error = $processStartError
    stdout_path = $resolvedStdout
    stderr_path = $resolvedStderr
    summary_path = $resolvedSummary
    stdout_bytes = $stdoutBytes
    stderr_bytes = $stderrBytes
    stream_capture_error = $streamCaptureError
    process_metrics_error = $processMetricsError
    output_root = $resolvedOutput
    archive_path = $resolvedArchive
    launcher_path = $launcherPath
    launcher_arguments = $launcherArguments
    child_process_ids = @($treeIds | Where-Object { $_ -ne $process.Id })
    child_process_count = @($treeIds | Where-Object { $_ -ne $process.Id }).Count
    remaining_process_ids = $remainingTreeIds
    child_process_handling = if ($null -eq $treeError) { 'observed_and_tree_terminated' } else { 'UNRESOLVED' }
    child_process_error = $treeError
    termination_method = if ($null -ne $termination) { $termination.method } else { $null }
    termination_taskkill_exit_code = if ($null -ne $termination) { $termination.taskkill_exit_code } else { $null }
    archive_member_enforcement = 'post_finalize_validation'
    enforcement_scope = 'primary_process_polling'
    strict_kernel_limits = $false
}
Write-Summary $summary $resolvedSummary

if ($null -ne $trigger) {
    exit 20
}
if ($null -ne $outputScanError) {
    exit 22
}
if ($null -ne $exitCode -and $exitCode -ne 0) {
    exit 10
}
if ($RequireArchive -and $null -eq $lastMembers) {
    [Console]::Error.WriteLine('Required archive is missing or is not a finalized readable ZIP archive.')
    exit 23
}
if ($null -ne $lastMembers -and $lastMembers -gt $MaxArchiveMembers) {
    [Console]::Error.WriteLine("Archive member limit exceeded: $lastMembers > $MaxArchiveMembers.")
    exit 21
}
exit 0
