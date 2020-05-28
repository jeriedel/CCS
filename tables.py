from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import QPalette, QRegExpValidator


class TableResults(QFrame):
    def __init__(self, *args, **kwargs):
        super(TableResults, self).__init__(*args, **kwargs)
        self.coeffs_table   = QTableView()
        self.ccs_table      = QTableView()
        self.fit_table      = QTableView()

        coeffs_header   = self.coeffs_table.horizontalHeader() 
        ccs_header      = self.ccs_table.horizontalHeader() 
        fit_header      = self.fit_table.horizontalHeader() 

        coeffs_header.setSectionResizeMode(QHeaderView.Stretch)
        ccs_header.setSectionResizeMode(QHeaderView.Stretch)
        fit_header.setSectionResizeMode(QHeaderView.Stretch)

        self.__init_tables()
        self.setStyleSheet(open("QSS/tables.qss").read())

    def __init_tables(self):
        table_area_layout   = QVBoxLayout(self)
        table_area_layout.setContentsMargins(0,0,0,0)

        label_coeffs = QLabel("Coefficients") 
        label_ccs    = QLabel("CCS results") 
        label_fit    = QLabel("Fit results") 

        table_area_layout.addWidget(label_coeffs)
        table_area_layout.addWidget(self.coeffs_table)
        table_area_layout.addWidget(label_ccs)
        table_area_layout.addWidget(self.ccs_table)
        table_area_layout.addWidget(label_fit)
        table_area_layout.addWidget(self.fit_table)