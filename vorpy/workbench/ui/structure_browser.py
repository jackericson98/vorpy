"""Virtual selection list: allocate data per group, widgets only for visible rows."""
from PySide6.QtCore import QAbstractListModel, QModelIndex, QSortFilterProxyModel, Qt, Signal
from PySide6.QtWidgets import QListView


class SelectionModel(QAbstractListModel):
    checked = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.entries = []
        self.selection = set()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.entries)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or not 0 <= index.row() < len(self.entries):
            return None
        label, indices = self.entries[index.row()]
        if role == Qt.DisplayRole:
            return label
        if role == Qt.UserRole:
            return indices
        if role == Qt.CheckStateRole:
            count = len(self.selection.intersection(indices))
            return Qt.Checked if count == len(indices) else Qt.PartiallyChecked if count else Qt.Unchecked
        return None

    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable

    def setData(self, index, value, role=Qt.EditRole):
        if role != Qt.CheckStateRole or not index.isValid():
            return False
        indices = self.entries[index.row()][1]
        if value in (Qt.Checked, Qt.Checked.value):
            self.selection.update(indices)
        else:
            self.selection.difference_update(indices)
        self.dataChanged.emit(index, index, [Qt.CheckStateRole])
        self.checked.emit(index.row())
        return True

    def reset_entries(self, entries, selection):
        self.beginResetModel()
        self.entries = entries
        self.selection = set(selection)
        self.endResetModel()


class BrowserItem:
    """Stable Python handle for a source row, including filtered-out rows."""
    def __init__(self, browser, row):
        self.browser, self.row = browser, row

    def data(self, role):
        model = self.browser.source_model
        return model.data(model.index(self.row, 0), role)

    def text(self):
        return self.data(Qt.DisplayRole)

    def checkState(self):
        return self.data(Qt.CheckStateRole)

    def setCheckState(self, state):
        model = self.browser.source_model
        model.setData(model.index(self.row, 0), state, Qt.CheckStateRole)

    def isHidden(self):
        return self.browser.query not in self.text().casefold()


class StructureBrowser(QListView):
    itemChanged = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.query = ''
        self.source_model = SelectionModel(self)
        self.proxy = QSortFilterProxyModel(self)
        self.proxy.setSourceModel(self.source_model)
        self.proxy.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self.setModel(self.proxy)
        self.setUniformItemSizes(True)
        self.source_model.checked.connect(lambda row: self.itemChanged.emit(self.item(row)))

    def set_entries(self, entries, selection):
        self.source_model.reset_entries(entries, selection)

    def set_selection(self, selection):
        model = self.source_model
        if model.selection != selection:
            model.selection = set(selection)
            if model.entries:
                model.dataChanged.emit(model.index(0, 0), model.index(len(model.entries) - 1, 0),
                                       [Qt.CheckStateRole])

    def set_query(self, query):
        self.query = query.strip().casefold()
        self.proxy.setFilterFixedString(self.query)

    def clear(self):
        self.set_entries([], set())

    def count(self):
        return len(self.source_model.entries)

    def item(self, row):
        return BrowserItem(self, row)
