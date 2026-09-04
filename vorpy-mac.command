#!/bin/sh
# Double-clickable macOS launcher. Always run relative to this checkout.
set -u

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$repo_dir"

if [ -x "$repo_dir/.venv/bin/python" ]; then
    python_cmd="$repo_dir/.venv/bin/python"
elif command -v python3 >/dev/null 2>&1; then
    python_cmd=python3
else
    printf '%s\n' 'VorPy needs Python 3. Install Python, then run:' >&2
    printf '  %s\n' 'python3 -m pip install -e ".[gui]"' >&2
    printf '%s' 'Press Return to close...'
    read -r _answer
    exit 1
fi

if "$python_cmd" -m vorpy.workbench "$@"; then
    status=0
else
    status=$?
fi
if [ "$status" -ne 0 ]; then
    printf '\nVorPy exited with an error. Install this checkout with:\n' >&2
    printf '  %s\n' 'python3 -m pip install -e ".[gui]"' >&2
    printf '%s' 'Press Return to close...'
    read -r _answer
fi
exit "$status"
