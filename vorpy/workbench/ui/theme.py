"""Application-wide styling for the dark analysis-studio interface."""

STUDIO_STYLESHEET = """
QWidget { color: #d9e2ec; background: #10161d; font-family: "Segoe UI", "Inter", sans-serif; font-size: 10pt; }
QMainWindow, QMenuBar, QMenu, QStatusBar { background: #0b1117; }
QMenuBar { border-bottom: 1px solid #27313d; padding: 2px; }
QMenuBar::item { padding: 5px 10px; border-radius: 4px; }
QMenuBar::item:selected, QMenu::item:selected { background: #25304a; color: #ffffff; }
QMenu { border: 1px solid #303b48; padding: 5px; }
QMenu::item { padding: 6px 28px 6px 10px; }
QToolBar { background: #101820; border: none; border-bottom: 1px solid #27313d; spacing: 4px; padding: 5px 8px; }
QToolBar QToolButton { padding: 6px 9px; border-radius: 5px; }
QToolButton, QPushButton { background: #18222d; border: 1px solid #33404e; border-radius: 5px; padding: 6px 10px; }
QToolButton:hover, QPushButton:hover { background: #233141; border-color: #596b7e; }
QToolButton:pressed, QPushButton:pressed { background: #2c385f; }
QToolButton:checked { background: #5546b8; border-color: #7869e3; color: white; }
QToolButton#railButton { border: none; border-radius: 7px; padding: 9px 3px; min-height: 49px; color: #b8c5d2; }
QToolButton#railButton:hover { background: #1a2530; color: white; }
QToolButton#railButton:checked { background: #5546b8; color: white; }
QFrame#toolRail { background: #0b1117; border-right: 1px solid #27313d; }
QWidget#workflowPanel { background: #111820; border-right: 1px solid #27313d; }
QLabel#brandMark { background: transparent; border: none; padding: 0; }
QFrame#viewerBar, QFrame#metricCard { background: #141d26; border: 1px solid #2b3743; border-radius: 6px; }
QFrame#metricCard { border-top: 2px solid #6857d9; }
QLabel#metricTitle { color: #8fa0b2; font-size: 9pt; }
QLabel#metricValue { color: #f3f7fb; font-size: 18pt; font-weight: 600; }
QLabel#sectionLabel { color: #8fa0b2; font-size: 9pt; font-weight: 600; }
QTabWidget::pane { border: 1px solid #2b3743; background: #111820; }
QTabBar::tab { background: #111820; color: #9cabba; padding: 9px 14px; border-bottom: 2px solid transparent; }
QTabBar::tab:selected { color: #ffffff; border-bottom-color: #806df0; }
QTabBar::tab:hover { color: #ffffff; }
QGroupBox { border: 1px solid #2d3945; border-radius: 6px; margin-top: 12px; padding: 12px 8px 8px 8px; font-weight: 600; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; color: #aebdca; }
QTreeWidget, QTableWidget, QComboBox, QSpinBox { background: #111820; alternate-background-color: #141e28; border: 1px solid #2d3945; border-radius: 4px; selection-background-color: #4f449e; selection-color: white; }
QTreeWidget::item { padding: 5px; }
QHeaderView::section { background: #18222d; color: #9eafbf; border: none; border-right: 1px solid #2d3945; border-bottom: 1px solid #2d3945; padding: 6px; }
QComboBox, QSpinBox { padding: 4px 7px; min-height: 20px; }
QComboBox::drop-down { border: none; width: 22px; }
QCheckBox { spacing: 6px; }
QCheckBox::indicator { width: 15px; height: 15px; }
QCheckBox::indicator:unchecked { background: #0d141b; border: 1px solid #516170; border-radius: 3px; }
QCheckBox::indicator:checked { background: #6b59df; border: 1px solid #8b7af2; border-radius: 3px; }
QSlider::groove:horizontal { height: 4px; background: #2c3742; border-radius: 2px; }
QSlider::handle:horizontal { width: 13px; margin: -5px 0; background: #806df0; border-radius: 6px; }
QProgressBar { background: #17212b; border: 1px solid #2c3945; border-radius: 4px; text-align: center; min-height: 14px; }
QProgressBar::chunk { background: #6654d7; border-radius: 3px; }
QSplitter::handle { background: #27313d; }
QSplitter::handle:horizontal { width: 1px; }
QSplitter::handle:vertical { height: 1px; }
QStatusBar { border-top: 1px solid #27313d; color: #9aabba; }
QScrollBar:vertical { width: 10px; background: #10161d; }
QScrollBar::handle:vertical { background: #364453; border-radius: 5px; min-height: 24px; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
"""


