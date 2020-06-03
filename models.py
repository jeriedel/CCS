import pandas as pd
import numpy as np
from PyQt5 import QtCore

class PandasModel(QtCore.QAbstractTableModel):
    def __init__(self, data, parent=None):
         super(PandasModel, self).__init__()
         self._data = data.round(4)
         self.header_labels = [str(i+1) for i in range(data.shape[1])]

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self.header_labels[section]
        return QtCore.QAbstractTableModel.headerData(self, section, orientation, role)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(str(
                self._data.values[index.row()][index.column()]))

class ArrayModel(QtCore.QAbstractTableModel):
    def __init__(self, data, parent=None):
        super(ArrayModel, self).__init__()

        if not isinstance(data, np.ndarray):
            raise Exception
        else:
            self._data = data.reshape((1, len(data)))
            self._data = np.around(self._data, 4)
            self.header_labels = [str(i+1) for i in range(self._data.shape[1])]

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self.header_labels[section]
        return QtCore.QAbstractTableModel.headerData(self, section, orientation, role)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(str(
                self._data[index.row()][index.column()]))

class CcsModel(ArrayModel):
    def __init__(self, data, parent=None):
        super().__init__(data)

class FitModel(ArrayModel):
    def __init__(self, data, parent=None):
        super().__init__(data)
