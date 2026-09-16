"""Main entry point for VorPy."""

from __future__ import annotations

import sys
from pathlib import Path


# Directory execution (python3 vorpy ...) puts the package directory on
# sys.path; absolute imports need its parent to select the local checkout.
if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def main() -> None:
    if len(sys.argv) == 1:
        from vorpy.workbench.__main__ import main as gui_main

        raise SystemExit(gui_main())

    from vorpy.src.command.vpy_cmnd import Command

    Command().run()


if __name__ == "__main__":
    main()
