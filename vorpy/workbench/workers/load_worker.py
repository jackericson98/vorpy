"""Read structures off the GUI thread; VTK and widgets stay on the GUI thread."""
from threading import Event

from PySide6.QtCore import QThread, Signal


class LoadWorker(QThread):
    progress = Signal(str, int)
    preview = Signal(object)

    def __init__(self, reader, parent=None):
        super().__init__(parent)
        self.reader = reader
        self.result = None
        self.error = None
        self.preview_displayed = Event()

    def run(self):
        try:
            self.result = self.reader(self.progress.emit, self._publish_preview)
        except Exception as error:
            self.error = error

    def _publish_preview(self, result):
        self.preview.emit(result)
        # Give the GUI first use of the CPU; then continue solvent work while
        # the user navigates the molecule. The GUI acknowledges even on error.
        self.preview_displayed.wait()
