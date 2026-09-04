from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication

from vorpy.workbench.ui.main_window import MainWindow
from vorpy.workbench.ui.theme import STUDIO_STYLESHEET


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv if argv is None else argv
    app = QApplication(arguments)
    app.setApplicationName("VorPy")
    app.setOrganizationName("Ericson Labs")
    app.setDesktopFileName("vorpy")
    icon_path = Path(__file__).resolve().parent / "assets" / "VorpyIcon_transparent.png"
    app.setWindowIcon(QIcon(str(icon_path)))
    app.setStyle("Fusion")
    app.setStyleSheet(STUDIO_STYLESHEET)
    window = MainWindow()
    window.setWindowIcon(app.windowIcon())
    window.show()
    if len(arguments) > 1:
        selected_path = Path(arguments[1]).expanduser().resolve()
        QTimer.singleShot(0, lambda: window.load_path(selected_path))
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
