"""Background preloading for PDB trajectory frames."""
from PySide6.QtCore import QThread, Signal

from vorpy.workbench.services.structure_loader import load_pdb


class TrajectoryPreloadWorker(QThread):
    frame_loaded = Signal(object)
    failed = Signal(str)

    def __init__(self, source, frame_ranges, current_frame, parent=None):
        super().__init__(parent)
        self.source = source
        self.frame_ranges = frame_ranges
        self.current_frame = current_frame

    def run(self):
        for number in range(1, len(self.frame_ranges) + 1):
            if number == self.current_frame or self.isInterruptionRequested():
                continue
            try:
                result = load_pdb(self.source, frame_index=number,
                                  frame_ranges=self.frame_ranges)
                self.frame_loaded.emit(result)
            except Exception as error:  # noqa: BLE001
                self.failed.emit(f"Frame {number}: {error}")
