"""Reusable presentation containers; scientific state remains in the workbench."""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QScrollArea, QToolButton, QVBoxLayout, QWidget, QSizePolicy,
)


def action_button(action, text=None):
    button = QToolButton()
    button.setDefaultAction(action)
    if text:
        button.setText(text)
        action.changed.connect(lambda: button.setText(text))
    button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
    return button


def scroll_panel(widget):
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.NoFrame)
    scroll.setWidget(widget)
    return scroll


class WorkflowSidebar(QWidget):
    def __init__(self, tabs):
        super().__init__()
        self.setObjectName("workflowPanel")
        self.setMinimumWidth(240)
        self.setMaximumWidth(340)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(tabs)


class ViewerPanel(QWidget):
    def __init__(self, viewer, actions):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 0)
        layout.setSpacing(4)
        bar = QFrame()
        bar.setObjectName("viewerBar")
        tools = QHBoxLayout(bar)
        tools.setContentsMargins(6, 3, 6, 3)
        tools.setSpacing(3)
        for action, label in actions:
            tools.addWidget(action_button(action, label))
        tools.addStretch()
        bar.setToolTip("Drag to rotate · middle-drag to pan · wheel for depth clipping · Shift-click to toggle selection")
        layout.addWidget(bar)
        layout.addWidget(viewer, 1)
        self.metadata = QLabel("No structure loaded")
        self.metadata.setObjectName("sectionLabel")
        self.metadata.setWordWrap(True)
        layout.addWidget(self.metadata)
        self.selection_label = QLabel("Residue selection active")
        self.selection_label.setObjectName("sectionLabel")
        layout.addWidget(self.selection_label)


class ViewInspector(QWidget):
    def __init__(self, tabs):
        super().__init__()
        self.setMinimumWidth(360)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 6)
        title = QLabel("View Settings")
        title.setObjectName("viewSettingsTitle")
        layout.addWidget(title)
        layout.addWidget(tabs, 1)


class ResultsInspector(QWidget):
    def __init__(self, content):
        super().__init__()
        self.setObjectName("resultsInspector")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 0, 8, 4)
        layout.setSpacing(2)
        self.toggle = QToolButton()
        self.toggle.setText("Results / Analysis")
        self.toggle.setCheckable(True)
        self.toggle.setChecked(True)
        self.toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle.setArrowType(Qt.DownArrow)
        layout.addWidget(self.toggle, alignment=Qt.AlignLeft)
        content.setMinimumHeight(0)
        content.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Ignored)
        self.content = content
        layout.addWidget(content, 1)
        self.toggle.toggled.connect(self.set_expanded)

    def set_expanded(self, expanded):
        self.content.setVisible(expanded)
        self.toggle.setArrowType(Qt.DownArrow if expanded else Qt.RightArrow)
        self.setMaximumHeight(16777215 if expanded else 34)
        if self.toggle.isChecked() != expanded:
            self.toggle.setChecked(expanded)
