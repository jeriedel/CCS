from PyQt5.QtWidgets import *
from PyQt5.QtCore import  *
import numpy as np

ERROR_MESSAGE = """The provided filename is not valid. Please select a new filename"""

class ExportDialog(QDialog):
    def __init__(self, parent, plot=True, figure=None, data=None, *args, **kwargs):
        super(ExportDialog, self).__init__(parent, *args, **kwargs)
        if plot and figure is not None:
            self.plot = True
            self.figure = figure
            self.output = None
            self.__render_ui_plot()
        else:
            self.plot       = False
            self.data       = data
            self.coeffFile  = None
            self.ccsFile    = None
            self.statsFile  = None
            self.__render_ui_tables()

    def __render_ui_plot(self):
        self.setWindowTitle("Export Plot")
        mainLayout = QVBoxLayout(self)
        exportArea = QFrame()
        buttonArea = QFrame()

        exportLayout = QFormLayout(exportArea)
        buttonLayout = QHBoxLayout(buttonArea)

        exportLayout.setFormAlignment(Qt.AlignLeft)
        exportLayout.setLabelAlignment(Qt.AlignLeft)

        group = QFrame()
        groupLayout = QHBoxLayout(group)
        group.setContentsMargins(0,0,0,0)
        groupLayout.setContentsMargins(0,0,0,0)

        browse = QPushButton()
        browse.setText("Browse")
        self.file   = QLineEdit()
        groupLayout.addWidget(self.file)
        groupLayout.addWidget(browse)

        browse.clicked.connect(self.openFileDialog)
        self.file.textChanged.connect(self.outputFileChanged)

        quit = QPushButton("Quit")
        export = QPushButton("Export")

        quit.clicked.connect(self.close)
        export.clicked.connect(self.startExport)

        exportLayout.addRow(QLabel("Export File:"), group)

        buttonLayout.addWidget(quit)
        buttonLayout.addWidget(export)

        mainLayout.addWidget(exportArea)
        mainLayout.addWidget(buttonArea)

    def __render_ui_tables(self):
        self.setWindowTitle("Export Tables")
        mainLayout = QVBoxLayout(self)
        exportArea = QFrame()
        buttonArea = QFrame()

        exportLayout = QFormLayout(exportArea)
        buttonLayout = QHBoxLayout(buttonArea)

        exportLayout.setContentsMargins(0,0,0,0)

        exportLayout.setFormAlignment(Qt.AlignLeft)
        exportLayout.setLabelAlignment(Qt.AlignLeft)

        frameCoeff  = QFrame()
        frameCCS    = QFrame()
        frameStats  = QFrame()

        frameCoeff.setContentsMargins(0,0,0,0)
        frameCCS.setContentsMargins(0,0,0,0)
        frameStats.setContentsMargins(0,0,0,0)

        layoutCoeff = QHBoxLayout(frameCoeff)
        layoutCCS   = QHBoxLayout(frameCCS)
        layoutStats = QHBoxLayout(frameStats)

        layoutCoeff.setContentsMargins(0,0,0,0)
        layoutCCS.setContentsMargins(0,0,0,0)
        layoutStats.setContentsMargins(0,0,0,0)

        self.coeffFileEdit   = QLineEdit()
        self.ccsFileEdit     = QLineEdit()
        self.statsFileEdit   = QLineEdit()

        browseCoeff          = QPushButton("Browse")
        browseCCS            = QPushButton("Browse")
        browseStats          = QPushButton("Browse")

        browseCoeff.clicked.connect(lambda: self.openFileDialog(key='table', to='coeff'))
        browseCCS.clicked.connect(lambda: self.openFileDialog(key='table', to='ccs'))
        browseStats.clicked.connect(lambda: self.openFileDialog(key='table', to='stats'))

        layoutCoeff.addWidget(self.coeffFileEdit)
        layoutCoeff.addWidget(browseCoeff)
        layoutCCS.addWidget(self.ccsFileEdit)
        layoutCCS.addWidget(browseCCS)
        layoutStats.addWidget(self.statsFileEdit)
        layoutStats.addWidget(browseStats)

        exportLayout.addRow(QLabel("Export Coeff. Table:"), frameCoeff)
        exportLayout.addRow(QLabel("Export CCS Table:"), frameCCS)
        exportLayout.addRow(QLabel("Export Stats Table:"), frameStats)

        quit = QPushButton("Quit")
        export = QPushButton("Export")

        quit.clicked.connect(self.close)
        export.clicked.connect(self.startExport)

        buttonLayout.addWidget(quit)
        buttonLayout.addWidget(export)

        mainLayout.addWidget(exportArea)
        mainLayout.addWidget(buttonArea)

    def openFileDialog(self, key='plot', to='coff'):
        if key == 'plot':
            self.output, filter = QFileDialog.getSaveFileName(self, "Output File", QDir.homePath())
            if not self.output == '':
                self.file.setText(self.output)
            else:
                self.output = None
        else:
            file, filter = QFileDialog.getSaveFileName(self, "Output File", QDir.homePath())
            if not file == '':
                if to == 'coeff':
                    self.coeffFile = file
                    self.coeffFileEdit.setText(file)
                elif to == 'ccs':
                    self.ccsFile = file
                    self.ccsFileEdit.setText(file)
                else:
                    self.statsFile = file
                    self.statsFileEdit.setText(file)

    def outputFileChanged(self, text):
        if not text == "":
            self.output = text
        else:
            self.output = None

    def startExport(self):
        if self.plot:
            if self.output is not None:
                self.figure.savefig(self.output)
                self.close()
            else:
                QMessageBox.warning(self, "Form Error", ERROR_MESSAGE, QMessageBox.Ok)
        else:
            coeff, ccs, stats = self.data
            if self.coeffFile is not None:
                coeff.to_csv(self.coeffFile, sep='\t')
            if self.ccsFile is not None:
                ccs = np.asarray(ccs)
                if len(ccs) == 2:
                    ccs = ccs.reshape((-1, 2))
                    np.savetxt(self.ccsFile, ccs, header="CCS\terror", delimiter='\t')
                else:
                    ccs = ccs.reshape((-1, 3))
                    np.savetxt(self.ccsFile, ccs, header="CCS\tCCS corr.\terror", delimiter='\t')
            if self.statsFile is not None:
                stats = np.asarray(stats)
                stats = np.reshape(stats, (-1, 3))
                np.savetxt(self.statsFile, stats, header="R\tP\tstd. error", delimiter="\t")

            self.close()
