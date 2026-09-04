from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen


def _configure_windows_taskbar() -> None:
    """Give an unpackaged Windows launch its own taskbar identity."""
    if sys.platform != "win32":
        return
    try:
        import ctypes

        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
            "EricsonLabs.VorPy.Workbench"
        )
    except (AttributeError, OSError):
        # The GUI can still run if an unusual Windows environment does not
        # expose the taskbar API.
        pass


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv if argv is None else argv
    _configure_windows_taskbar()
    app = QApplication(arguments)
    app.setApplicationName("VorPy")
    app.setOrganizationName("Ericson Labs")
    app.setDesktopFileName("vorpy")
    icon_path = Path(__file__).resolve().parent / "assets" / "VorpyIcon_transparent.png"
    app.setWindowIcon(QIcon(str(icon_path)))
    app.setStyle("Fusion")

    splash = QSplashScreen(QPixmap(str(icon_path)))
    splash.showMessage(
        "Starting VorPy Analysis Studio...",
        Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter,
        Qt.GlobalColor.white,
    )
    splash.show()
    app.processEvents()

    # These modules load the 3D visualization stack and are intentionally
    # imported after the splash screen is visible.
    from vorpy.workbench.ui.main_window import MainWindow
    from vorpy.workbench.ui.theme import STUDIO_STYLESHEET

    app.setStyleSheet(STUDIO_STYLESHEET)
    window = MainWindow()
    window.setWindowIcon(app.windowIcon())
    window.show()
    splash.finish(window)
    if len(arguments) > 1:
        selected_path = Path(arguments[1]).expanduser().resolve()
        QTimer.singleShot(0, lambda: window.load_path(selected_path))
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
