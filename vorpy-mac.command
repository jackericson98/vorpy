#!/bin/sh
# Double-clickable macOS launcher. Always run relative to this checkout.
set -u

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
venv_dir="$repo_dir/.venv"
python_cmd="$venv_dir/bin/python"
ready_file="$venv_dir/.vorpy-gui-ready"
launch_log="$venv_dir/vorpy-launch.log"
cd "$repo_dir"

fail_setup() {
    printf '\n%s\n' "$1" >&2
    printf '%s\n' 'Review the installation error above, then try again.' >&2
    printf '%s' 'Press Return to close...'
    read -r _answer
    exit 1
}

if [ ! -x "$python_cmd" ]; then
    if ! command -v python3 >/dev/null 2>&1; then
        printf '%s\n' 'Python 3 was not found. Install Python 3 from https://python.org and try again.' >&2
        printf '%s' 'Press Return to close...'
        read -r _answer
        exit 1
    fi
    printf '%s\n' '============================================================'
    printf '%s\n' 'VorPy first-time setup'
    printf '%s\n' '============================================================'
    printf '%s\n\n' 'Creating a local Python environment. This normally happens only once.'
    python3 -m venv "$venv_dir" || fail_setup "Could not create VorPy's Python environment."
fi

if [ ! -f "$ready_file" ]; then
    "$python_cmd" -c "import importlib.util, sys; names = ('PySide6', 'pyvista', 'pyvistaqt', 'vorpy.workbench'); sys.exit(0 if all(importlib.util.find_spec(name) for name in names) else 1)" >/dev/null 2>&1
    dependencies_ready=$?
    if [ "$dependencies_ready" -ne 0 ]; then
        "$python_cmd" "$repo_dir/vorpy/workbench/bootstrap.py" || fail_setup "VorPy's dependencies could not be installed."
    fi
    : > "$ready_file"
    printf '\n%s\n' 'Installation complete. Starting VorPy...'
fi

# Detach the GUI so Terminal can close after setup. Startup errors are retained
# in the launch log instead of leaving a terminal window attached to VorPy.
nohup "$python_cmd" -m vorpy.workbench "$@" >"$launch_log" 2>&1 &
exit 0
