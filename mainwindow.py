from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from overview import Overview
from canvas import Canvas
from tables import TableResults
from export import ExportDialog
import ccs.synapt as synapt
import ccs.imob as imob
from models import PandasModel, CcsModel, FitModel
import numpy as np
import pandas as pd

class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        # Set window size
        self.setGeometry(1200, 800, 1200, 800)
        # Set window title
        self.setWindowTitle("CCS-Calculator")
        # Create menu
        self.__set_menu()
        # initialize UI
        self.__init_ui()
        # default mode
        self.project = 'synapt'
        self.setStyleSheet("QSplitter::handle {background: #6d6d6d}")

        # Private variables
        self.header_file = None

    def __init_ui(self):
        self.__render_synapt_dialog()

    def __set_menu(self):
        # Define all actions in menu
        quit_app   = QAction('&Quit', self)
        synapt     = QAction('&Synapt', self)
        imob       = QAction('&Imob', self)
        drift      = QAction('&Drift Time Scan', self)
        ms         = QAction('&Visualize MS', self)
        export_plt = QAction('&Export Plot', self)
        export_tab = QAction('&Export Table', self)

        # Add signals
        quit_app.triggered.connect(qApp.quit)
        synapt.triggered.connect(lambda: self.__handle_project_change('synapt'))
        imob.triggered.connect(lambda: self.__handle_project_change('imob'))
        drift.triggered.connect(lambda: self.__handle_project_change('drift'))
        ms.triggered.connect(lambda : self.__handle_project_change('ms'))
        export_plt.triggered.connect(lambda: self.__export_plot())
        export_tab.triggered.connect(lambda: self.__export_table())

        # Add menu tabs
        main_menu = self.menuBar()
        ccs_menu  = main_menu.addMenu('&CSS')
        project   = main_menu.addMenu('&Project')
        export    = main_menu.addMenu('&Export')
        # Add actions
        ccs_menu.addAction(quit_app)
        project.addAction(synapt)
        project.addAction(imob)
        project.addAction(drift)
        project.addAction(ms)
        export.addAction(export_plt)
        export.addAction(export_tab)

    def __plot_results(self, imms_data, coeff, results):
        if self.project == 'synapt':
            self.graph.clear_plot()
            self.__plot_fit_synapt(imms_data, coeff)
            self.__plot_ccs_synapt(results)
            self.__populate_tables_synapt(coeff, results)
        elif self.project == 'imob':
            self.graph.clear_plot()
            self.__plot_fit_imob(imms_data, coeff)
            self.__plot_ccs_imob(results)
            self.__populate_tables_imob(coeff, results)
    
    def __plot_results_drift(self, res):
            self.graph.clear_plot(method='drift')
            self.graph.plot_drift_time_scope(res)
        
    def __plot_results_ms(self, data):
        self.graph.clear_plot(method='ms')
        self.graph.plot_ms(data)

    def __plot_fit_synapt(self, imms_data, coeffs):
        for key,  coeff in coeffs.items():
            temp = imms_data[key]
            x = temp['Time'].to_numpy()
            y = temp['Intensity'].to_numpy()
            fitted = synapt.normal_distribution(x, coeff)
            self.graph.plot_fit(x, y)
            self.graph.plot_fit(x, fitted, method='dashed')

    def __plot_fit_imob(self, imms_data, coeffs):
        time = imms_data.index
        for col in imms_data:
            intensity = np.asarray(imms_data[col])
            self.graph.plot_fit(time, intensity)

        coeffs = coeffs.to_numpy()

        for coeff in coeffs: 
            fitted = imob.normal_distribution(time, coeff)
            self.graph.plot_fit(time, fitted, method='dashed')

    def __plot_ccs_synapt(self, results):
        for exp, regress, ccs in results:
            exp_x, exp_y = exp
            fitted_x    = regress.voltage
            fitted_y    = regress.fitted
            self.graph.plot_ccs(exp_x, exp_y, method="scatter")
            self.graph.plot_ccs(fitted_x, fitted_y, method="plot")

    def __plot_ccs_imob(self, results):
        for exp, regress, ccs in results:
            exp_x, exp_y = exp
            fitted_x    = regress.voltage
            fitted_y    = regress.fitted
            self.graph.plot_ccs(exp_x, exp_y, method="scatter")
            self.graph.plot_ccs(fitted_x, fitted_y, method="plot")

    def __populate_tables_synapt(self, coeffs, results):
        _, regress, ccs_results = results[0]
        _, _, R, P, std_err = regress

        ccs_results = np.asarray(ccs_results)
        fit_params  = np.asarray([R, P, std_err])

        ccs_model    = CcsModel(ccs_results)
        ccs_model.header_labels = ['ccs', 'error']

        fit_model    = FitModel(fit_params) 
        fit_model.header_labels = ['R', 'P', 'std err.']

        coeffs = pd.DataFrame.from_dict(coeffs, orient='index', columns=['A', 'mu', 'sigma'])
        coeff_model = PandasModel(coeffs)
        coeff_model.header_labels = ['A', 'mu', 'sigma']

        self.table_area.coeffs_table.setModel(coeff_model)
        self.table_area.ccs_table.setModel(ccs_model)
        self.table_area.fit_table.setModel(fit_model)
        self.table_area.show()

    def __populate_tables_imob(self, coeffs, results):
        _, regress, ccs_results = results[0]
        _, _, slope, intercept, R, P, std_err = regress

        ccs_results = np.asarray(ccs_results)
        fit_params  = np.asarray([slope, intercept, R, P, std_err])

        pandas_model = PandasModel(coeffs)
        pandas_model.header_labels = ['A', 'mu', 'sigma', 'off']

        ccs_model    = CcsModel(ccs_results)
        ccs_model.header_labels = ['ccs', 'ccs corr.', 'error']

        fit_model    = FitModel(fit_params) 
        fit_model.header_labels = ['slope', 'intercept', 'R', 'P', 'std err.']

        self.table_area.coeffs_table.setModel(pandas_model)
        self.table_area.ccs_table.setModel(ccs_model)
        self.table_area.fit_table.setModel(fit_model)
        self.table_area.show()

    def __handle_project_change(self, project):
        self.setCentralWidget(None)

        if project == 'synapt':
            self.project = 'synapt'
            self.__render_synapt_dialog()
        elif project == 'imob':
            self.project = 'imob'
            self.__render_imob_dialog()
        elif project == 'drift':
            self.project == 'drift'
            self.__render_drifttime_dialog()
        elif project == 'ms':
            self.project = 'ms'
            self.__render_ms_dialog()

    def __render_synapt_dialog(self):
        central     = QFrame()
        hbox        = QHBoxLayout(central)
        hbox.setContentsMargins(0,0,0,0)

        # Define components visible in main layout
        self.params         = Overview('synapt')
        self.graph          = Canvas()

        # sublayout for results
        results_area        = QFrame()
        results_area_layout = QHBoxLayout(results_area)
        results_area_layout.setContentsMargins(0,0,0,0)

        self.table_area     = TableResults()
        self.table_area.hide()

        results_splitter = QSplitter(Qt.Horizontal)
        results_splitter.addWidget(self.graph)
        results_splitter.addWidget(self.table_area)

        results_area_layout.addWidget(results_splitter)

        # Add splitter for main layout
        splitter         = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.params)
        splitter.addWidget(results_area)

        splitter.setStretchFactor(1,10)
        hbox.addWidget(splitter)

        self.setCentralWidget(central)

        self.params.ccs_synapt_success.connect(self.__plot_results) 

        self.show()

    def __render_imob_dialog(self):
        central     = QFrame()
        hbox        = QHBoxLayout(central)
        hbox.setContentsMargins(0,0,0,0)

        # Define components visible in main layout
        self.params         = Overview('imob')
        self.graph          = Canvas()

        # sublayout for results
        results_area     = QFrame()
        results_area_layout = QHBoxLayout(results_area)
        results_area_layout.setContentsMargins(0,0,0,0)

        self.table_area     = TableResults()
        self.table_area.hide()

        results_splitter = QSplitter(Qt.Horizontal)
        results_splitter.addWidget(self.graph)
        results_splitter.addWidget(self.table_area)

        results_area_layout.addWidget(results_splitter)

        # Add splitter for main layout
        splitter         = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.params)
        splitter.addWidget(results_area)

        splitter.setStretchFactor(1,10)
        hbox.addWidget(splitter)

        self.setCentralWidget(central)

        self.params.ccs_imob_success.connect(self.__plot_results) 

        self.show()

    def __render_drifttime_dialog(self):
        central     = QFrame()
        hbox        = QHBoxLayout(central)
        hbox.setContentsMargins(0,0,0,0)

        # Define components visible in main layout
        self.params         = Overview('drift')
        self.graph          = Canvas(project='drift')

        # Add splitter for main layout
        splitter         = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.params)
        splitter.addWidget(self.graph)

        splitter.setStretchFactor(1,10)
        hbox.addWidget(splitter)

        self.setCentralWidget(central)
        
        self.params.drift_success.connect(self.__plot_results_drift)
        self.show()
        
    def __render_ms_dialog(self):
        central = QFrame()
        hbox = QHBoxLayout(central)
        hbox.setContentsMargins(0, 0, 0, 0)

        # Define components visible in main layout
        self.params = Overview('ms')
        self.graph = Canvas(project='ms')

        # Add splitter for main layout
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.params)
        splitter.addWidget(self.graph)

        splitter.setStretchFactor(1, 10)
        hbox.addWidget(splitter)

        self.setCentralWidget(central)
        self.params.ms_success.connect(self.__plot_results_ms)
        self.show()

    def __export_plot(self):
        dialog = ExportDialog(self, figure=self.graph.figure, plot=True)
        dialog.exec()

    def __export_table(self):
        try:
            coeffs = self.table_area.coeffs_table.model().raw
            ccs    = self.table_area.ccs_table.model().raw
            stats  = self.table_area.fit_table.model().raw
            tdata = (coeffs, ccs, stats)
            dialog = ExportDialog(self, plot=False, data=tdata)
            dialog.exec()
        except:
            print("Currently, there are no tables for plotting available.")

