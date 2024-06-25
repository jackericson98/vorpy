import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
                             QLineEdit, QLabel, QComboBox, QCheckBox, QGridLayout, QFrame,
                             QFileDialog, QListWidget)
from PyQt5.QtCore import Qt
from Visualize.GUIs.sphere_viewer import GLWidget


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Vorpy')
        self.init_ui()

    def init_ui(self):
        # Main layout
        main_layout = QVBoxLayout()

        # Top area: Load Files and Settings
        top_layout = QHBoxLayout()
        top_layout.addWidget(self.load_files_section())
        top_layout.addWidget(self.settings_section())

        # Bottom area: Viewer and Exports
        bottom_layout = QHBoxLayout()
        bottom_layout.addWidget(self.viewer_section())
        bottom_layout.addWidget(self.exports_section())

        # Add top and bottom layouts to main layout
        main_layout.addLayout(top_layout)
        main_layout.addLayout(bottom_layout)

        # Set main widget
        main_widget = QWidget()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

    def load_files_section(self):
        section = FileDropWidget()
        layout = QVBoxLayout()
        section_title = QLabel("Load Files")
        section_title.setAlignment(Qt.AlignCenter)
        layout.addWidget(section_title)

        self.file_list = QListWidget()
        layout.addWidget(self.file_list)

        browse_button = QPushButton('Browse Files')
        browse_button.clicked.connect(self.browse_files)
        layout.addWidget(browse_button)

        section.setLayout(layout)
        return section

    def settings_section(self):
        section = QWidget()
        layout = QGridLayout()

        layout.addWidget(QLabel('Surface Resolution:'), 0, 0)
        surface_resolution = QLineEdit('0.2')
        layout.addWidget(surface_resolution, 0, 1)

        layout.addWidget(QLabel('Max Vertex:'), 1, 0)
        max_vertex = QLineEdit('5')
        layout.addWidget(max_vertex, 1, 1)

        layout.addWidget(QLabel('Box Multiplier:'), 2, 0)
        box_multiplier = QLineEdit()
        layout.addWidget(box_multiplier, 2, 1)

        layout.addWidget(QLabel('Network Type:'), 3, 0)
        network_type = QComboBox()
        network_type.addItems(['Additively Weighted', 'Power', 'Primitive'])
        layout.addWidget(network_type, 3, 1)

        change_radii_button = QPushButton('Change Atomic Radii')
        layout.addWidget(change_radii_button, 4, 0, 1, 2)

        save_surfaces = QCheckBox('Save Surfaces')
        layout.addWidget(save_surfaces, 5, 0, 1, 2)

        section.setLayout(layout)
        return section

    def viewer_section(self):
        section = QWidget()
        layout = QVBoxLayout()
        layout.addWidget(GLWidget())
        section.setLayout(layout)
        return section

    def exports_section(self):
        section = QWidget()
        layout = QVBoxLayout()

        export_pdb = QCheckBox('PDB file')
        export_logs = QCheckBox('Logs File')
        export_shell = QCheckBox('Shell')
        export_set_atoms = QCheckBox('Set Atoms File')
        export_vertices = QCheckBox('Vertices File')
        export_edges = QCheckBox('Edges File')
        export_folder_button = QPushButton('Select Export Folder')
        export_folder_button.clicked.connect(self.select_export_folder)

        layout.addWidget(export_pdb)
        layout.addWidget(export_logs)
        layout.addWidget(export_shell)
        layout.addWidget(export_set_atoms)
        layout.addWidget(export_vertices)
        layout.addWidget(export_edges)
        layout.addWidget(export_folder_button)

        section.setLayout(layout)
        return section

    def browse_files(self):
        file_dialog = QFileDialog(self)
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        if file_dialog.exec_():
            file_paths = file_dialog.selectedFiles()
            for path in file_paths:
                self.file_list.addItem(path)

    def select_export_folder(self):
        folder_path = QFileDialog.getExistingDirectory(self, 'Select Folder')
        print("Export folder selected:", folder_path)


class FileDropWidget(QFrame):
    def __init__(self):
        super().__init__()
        self.setFrameShape(QFrame.StyledPanel)
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            self.parent().file_list.addItem(file_path)
        event.acceptProposedAction()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
