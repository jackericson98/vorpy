"""Main entry point for VorPy."""

from __future__ import annotations

import sys


def main() -> None:
    if len(sys.argv) == 1:
        from vorpy.workbench.__main__ import main as gui_main

        raise SystemExit(gui_main())

    from vorpy.src.command.vpy_cmnd import Command

    Command().run()


if __name__ == "__main__":
    main()
