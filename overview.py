from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import QRegExpValidator

import os
import re
from functools import partial
import pandas as pd
import numpy as np
import ccs.synapt as synapt
import ccs.imob as imob
import ccs.drifttime as drifttime
from helper.fileloader import FileLoader
from helper.correction import  correction
from error import InvalidForm, NotDirectory

class Overview(QFrame):
    ccs_synapt_success = pyqtSignal(dict, dict, list)
    ccs_imob_success   = pyqtSignal(pd.DataFrame, pd.DataFrame, list)
    drift_success      = pyqtSignal(tuple)
    ms_success         = pyqtSignal(pd.DataFrame)

    def __init__(self, project, *args, **kwargs):
        super(Overview, self).__init__(*args, **kwargs)
        
        if project == 'synapt':
            self.__set_default_values_synapt()
            self.__init_ui_synapt()
            self.project = 'synapt'
        elif project == 'imob':
            self.__set_default_values_imob()
            self.__init_ui_imob()
            self.project = 'imob'
        elif project == 'drift':
            self.__set_default_values_drift()
            self.__init_ui_drift()
            self.project = 'drift'
        
        elif project == 'ms':
            self.__init_ui_ms()
            self.__set_default_values_ms()
            self.project = 'ms'

        self.setStyleSheet(open("qss/overview.qss").read())

    def __set_default_values_synapt(self):
        self.base_dir   = None
        self.pat_imms   = r".*imms.txt"
        self.ext_file   = r"_extern.inf"
        self.log_file   = r"Logfile.txt"
        self.head_file  = r"_HEADER.txt"
        self.mass       = ""
        self.charge     = ""
        self.gaussian   = "1"
    
    def __set_default_values_imob(self):
        self.h5file     = "/Users/Jerome/Downloads/Rohdaten/test.h5"
        self.dev        = "190"
        self.mz         = None
        self.charge     = "1"
        self.ld         = "161.2"
        self.gas        = 'He'
        self.gaussian   = "1"
        self.selecting  = False
        self.selections = ""
    
    def __set_default_values_drift(self):
        self.h5file         = None
        self.minmz          = "50.0"
        self.maxmz          = "3000.0"
        self.scaling        = "0.1"
        self.drift_gates    = None
        self.mass_gates     = None
        self.project        = "0"
        
    def __set_default_values_ms(self):
        self.ftype  = "XY File"
        self.file   = None
        self.xcol   = "1"
        self.ycol   = "2"

    def __init_ui_synapt(self):
        vbox        = QVBoxLayout(self)
        header      = QWidget() 
        hbox        = QHBoxLayout(header)
        title       = QLabel("Project Selection Synapt")

        header.setStyleSheet(open("qss/header.qss").read())
        hbox.setAlignment(title, Qt.AlignCenter)
        hbox.addWidget(title)

        form        = QWidget()
        form_layout = QFormLayout()
        form_layout.setFormAlignment(Qt.AlignLeft)
        form_layout.setLabelAlignment(Qt.AlignLeft)
        base_dir    = QLabel("Base directory")
        pat_imms    = QLabel("Pattern imms files")
        ext_file    = QLabel("External file")
        log_file    = QLabel("Log file")
        head_file   = QLabel("Header file")
        mass        = QLabel("Mass")
        charge      = QLabel("Charge")
        gaussian    = QLabel("No. Gaussians")

        self.base_edit      = QLineEdit()
        self.browse         = QPushButton("Browse")
        browse_group        = QFrame()
        browse_group_lay    = QHBoxLayout(browse_group)

        self.pat_edit       = QLineEdit()
        self.ext_edit       = QLineEdit()
        self.log_edit       = QLineEdit()
        self.header_edit    = QLineEdit()
        self.mass_edit      = QLineEdit()
        self.charge_edit    = QLineEdit()
        self.gaussian_edit  = QLineEdit()

        self.browse.setStyleSheet("""
            font-size: 12pt;
            min-height: 20px;
            min-width: 100px;
        """) 
        browse_group_lay.setContentsMargins(0,0,0,0)
        browse_group_lay.addWidget(self.base_edit)
        browse_group_lay.addWidget(self.browse)

        self.browse.clicked.connect(partial(self.__create_file_dialog, True))
        self.base_edit.textChanged.connect(self.__base_dir_changed)
        self.pat_edit.textChanged.connect(self.__pat_imms_changed)
        self.ext_edit.textChanged.connect(self.__ext_file_changed)
        self.log_edit.textChanged.connect(self.__log_file_changed)
        self.header_edit.textChanged.connect(self.__head_file_changed)
        self.mass_edit.textChanged.connect(self.__mass_changed)
        self.charge_edit.textChanged.connect(self.__charge_changed)
        self.gaussian_edit.textChanged.connect(self.__gaussian_changed)

        self.pat_edit.setText(".*imms.txt")
        self.ext_edit.setText("_extern.inf")
        self.log_edit.setText("Logfile.txt")
        self.header_edit.setText("_HEADER.TXT")
        self.gaussian_edit.setText("1")

        float_rgx   = QRegExp("\d+\.*\d*")
        integer_rgx = QRegExp("\d+")

        validator_mass      = QRegExpValidator(float_rgx, self)
        validator_charge    = QRegExpValidator(integer_rgx, self)
        validator_gaussian  = QRegExpValidator(integer_rgx, self)

        self.mass_edit.setValidator(validator_mass)
        self.charge_edit.setValidator(validator_charge)
        self.gaussian_edit.setValidator(validator_gaussian)

        form.setLayout(form_layout)
        form_layout.addRow(base_dir,    browse_group)
        form_layout.addRow(pat_imms,    self.pat_edit)
        form_layout.addRow(ext_file,    self.ext_edit)
        form_layout.addRow(log_file,    self.log_edit)
        form_layout.addRow(head_file,   self.header_edit)
        form_layout.addRow(mass,        self.mass_edit)
        form_layout.addRow(charge,      self.charge_edit)
        form_layout.addRow(gaussian,    self.gaussian_edit)

        button_run = QPushButton("Run")
        button_run.clicked.connect(self.__run_calculation)

        vbox.addWidget(header)
        vbox.addWidget(form)
        vbox.addWidget(button_run)
        self.show()

    def __init_ui_imob(self):
        vbox        = QVBoxLayout(self)
        header      = QWidget() 
        hbox        = QHBoxLayout(header)
        title       = QLabel("Project Selection Imob")

        header.setStyleSheet(open("qss/header.qss").read())
        hbox.setAlignment(title, Qt.AlignCenter)
        hbox.addWidget(title)

        form        = QWidget()
        form_layout = QFormLayout()
        form_layout.setFormAlignment(Qt.AlignLeft)
        form_layout.setLabelAlignment(Qt.AlignLeft)

        h5file      = QLabel("HDF File")
        dev         = QLabel("Drift End Voltage")
        mz          = QLabel("m/z")
        charge      = QLabel("Charge")
        ld          = QLabel("Length drift tube")
        gas         = QLabel("Drift gas")
        gaussian    = QLabel("No. Gaussians")
        select      = QLabel("Select Peaks")
        selection   = QLabel("Provide selection")

        self.file_edit      = QLineEdit()
        self.browse         = QPushButton("Browse")
        browse_group        = QFrame()
        browse_group_lay    = QHBoxLayout(browse_group)

        self.dev_edit       = QLineEdit()
        self.mz_edit        = QLineEdit()
        self.charge_edit    = QLineEdit()
        self.ld_edit        = QLineEdit()
        self.gas_edit       = QLineEdit()
        self.gaussian_edit  = QLineEdit()
        self.select_cb      = QCheckBox()
        self.select_edit    = QLineEdit()

        self.browse.setStyleSheet("""
            font-size: 12pt;
            min-height: 20px;
            min-width: 100px;
        """)

        browse_group_lay.setContentsMargins(0,0,0,0)
        browse_group_lay.addWidget(self.file_edit)
        browse_group_lay.addWidget(self.browse)

        self.browse.clicked.connect(partial(self.__create_file_dialog, False, True))
        self.file_edit.textChanged.connect(self.__file_changed)
        self.dev_edit.textChanged.connect(self.__dev_changed)
        self.mz_edit.textChanged.connect(self.__mz_value_changed)
        self.charge_edit.textChanged.connect(self.__charge_changed)
        self.ld_edit.textChanged.connect(self.__ld_value_changed)
        self.gas_edit.textChanged.connect(self.__gas_changed)
        self.gaussian_edit.textChanged.connect(self.__gaussian_changed)
        self.select_cb.stateChanged.connect(self.__handle_select)
        self.select_edit.textChanged.connect(self.__selection_changed)

        self.dev_edit.setText("190")
        self.charge_edit.setText("1")
        self.ld_edit.setText("161.2")
        self.gas_edit.setText("He")
        self.gaussian_edit.setText("1")

        float_rgx   = QRegExp("\d+\.*\d*")
        integer_rgx = QRegExp("\d+")
        gas_rgx     = QRegExp("N2|n2|He|he")

        validator_dev       = QRegExpValidator(float_rgx, self)
        validator_mz        = QRegExpValidator(float_rgx, self)
        validator_charge    = QRegExpValidator(integer_rgx, self)
        validator_ld        = QRegExpValidator(float_rgx, self)
        validator_gaussian  = QRegExpValidator(integer_rgx, self)
        validator_gas       = QRegExpValidator(gas_rgx, self)

        self.dev_edit.setValidator(validator_dev)
        self.mz_edit.setValidator(validator_mz)
        self.charge_edit.setValidator(validator_charge)
        self.ld_edit.setValidator(validator_ld)
        self.gas_edit.setValidator(validator_gas)
        self.gaussian_edit.setValidator(validator_gaussian)

        form.setLayout(form_layout)
        form_layout.addRow(h5file,   browse_group)
        form_layout.addRow(dev,      self.dev_edit)
        form_layout.addRow(mz,       self.mz_edit)
        form_layout.addRow(charge,   self.charge_edit)
        form_layout.addRow(ld,       self.ld_edit)
        form_layout.addRow(gas,      self.gas_edit)
        form_layout.addRow(gaussian, self.gaussian_edit)
        form_layout.addRow(select,   self.select_cb)
        form_layout.addRow(selection, self.select_edit)

        button_run = QPushButton("Run")
        button_run.clicked.connect(self.__run_calculation)

        vbox.addWidget(header)
        vbox.addWidget(form)
        vbox.addWidget(button_run)
        self.show()

    def __init_ui_drift(self):
        vbox        = QVBoxLayout(self)
        header      = QWidget() 
        hbox        = QHBoxLayout(header)
        title       = QLabel("Drift Time Scan")

        header.setStyleSheet(open("qss/header.qss").read())
        hbox.setAlignment(title, Qt.AlignCenter)
        hbox.addWidget(title)

        form        = QWidget()
        form_layout = QFormLayout()
        form_layout.setFormAlignment(Qt.AlignLeft)
        form_layout.setLabelAlignment(Qt.AlignLeft)

        h5file      = QLabel("HDF File")
        minmz       = QLabel("Min m/z")
        maxmz       = QLabel("Max m/z")
        scale       = QLabel("Scaling")
        drift_gates = QLabel("Drift time gates")
        mass_gates  = QLabel("Mass gates")

        self.file_edit      = QLineEdit()
        self.browse         = QPushButton("Browse")
        browse_group        = QFrame()
        browse_group_lay    = QHBoxLayout(browse_group)

        self.browse.setStyleSheet("""
            font-size: 12pt;
            min-height: 20px;
            min-width: 100px;
        """) 
        browse_group_lay.setContentsMargins(0,0,0,0)
        browse_group_lay.addWidget(self.file_edit)
        browse_group_lay.addWidget(self.browse)

        self.minmz_edit     = QLineEdit()
        self.maxmz_edit     = QLineEdit()
        self.scale_edit     = QLineEdit()
        self.drift_gates_e  = QLineEdit()
        self.mass_gates_e   = QLineEdit()

        self.browse.clicked.connect(partial(self.__create_file_dialog, False, True))
        self.file_edit.textChanged.connect(self.__file_changed)
        self.minmz_edit.textChanged.connect(self.__minmz_changed)
        self.maxmz_edit.textChanged.connect(self.__maxmz_changed)
        self.scale_edit.textChanged.connect(self.__scale_changed)
        self.drift_gates_e.textChanged.connect(self.__drift_gates_changed)
        self.mass_gates_e.textChanged.connect(self.__mass_gates_changed)

        self.minmz_edit.setText("50.0")
        self.maxmz_edit.setText("3000.0")
        self.scale_edit.setText("0.1")

        float_rgx   = QRegExp("\d+\.*\d*")
        integer_rgx = QRegExp("\d+")

        float_validator = QRegExpValidator(float_rgx, self)
        int_validator   = QRegExpValidator(integer_rgx, self)

        self.minmz_edit.setValidator(float_validator)
        self.maxmz_edit.setValidator(float_validator)
        self.scale_edit.setValidator(float_validator)
        
        form.setLayout(form_layout)
        form_layout.addRow(h5file,      browse_group)
        form_layout.addRow(minmz,       self.minmz_edit) 
        form_layout.addRow(maxmz,       self.maxmz_edit)
        form_layout.addRow(scale,       self.scale_edit)
        form_layout.addRow(drift_gates, self.drift_gates_e)
        form_layout.addRow(mass_gates,  self.mass_gates_e)

        button_run = QPushButton("Run")
        button_run.clicked.connect(self.__run_calculation)

        vbox.addWidget(header)
        vbox.addWidget(form)
        vbox.addWidget(button_run)
        
    def __init_ui_ms(self):
        vbox        = QVBoxLayout(self)
        header      = QWidget()
        hbox        = QHBoxLayout(header)
        title       = QLabel("Visualize MS spectrum")
    
        header.setStyleSheet(open("qss/header.qss").read())
        hbox.setAlignment(title, Qt.AlignCenter)
        hbox.addWidget(title)
    
        form        = QWidget()
        form_layout = QFormLayout(form)
        form_layout.setFormAlignment(Qt.AlignLeft)
        form_layout.setLabelAlignment(Qt.AlignLeft)
        
        file_type   = QLabel("File Type")
        file        = QLabel("File")
        col_x       = QLabel("Col. No. Wavenumber")
        col_y       = QLabel("Col. No. Intensity")
        
        self.cb_file        = QComboBox()
        self.file_edit      = QLineEdit()
        self.browse         = QPushButton("Browse")
        self.col_no_x       = QLineEdit()
        self.col_no_y       = QLineEdit()

        browse_group        = QFrame()
        browse_group_lay    = QHBoxLayout(browse_group)
        browse_group_lay.setContentsMargins(0,0,0,0)
        self.browse.setStyleSheet("""
            font-size: 12pt;
            min-height: 20px;
            min-width: 100px;
        """)
        
        browse_group_lay.addWidget(self.file_edit)
        browse_group_lay.addWidget(self.browse)
        
        self.cb_file.addItem("XY File")
        self.cb_file.addItem("XY File with correction")
        self.cb_file.addItem("Imob File")
        self.col_no_x.setText("1")
        self.col_no_y.setText("2")
        
        IntValidator = QRegExpValidator(QRegExp(r"[0-9]+"), self)
        self.col_no_x.setValidator(IntValidator)
        self.col_no_y.setValidator(IntValidator)

        form_layout.addRow(file_type, self.cb_file)
        form_layout.addRow(file, browse_group)
        form_layout.addRow(col_x, self.col_no_x)
        form_layout.addRow(col_y, self.col_no_y)
        
        self.cb_file.currentTextChanged.connect(self.__file_type_changed)
        self.browse.clicked.connect(partial(self.__create_file_dialog, False, True))
        self.file_edit.textChanged.connect(self.__file_changed)
        self.col_no_x.textChanged.connect(self.__xcol_changed)
        self.col_no_y.textChanged.connect(self.__ycol_changed)

        button_run = QPushButton("Run")
        button_run.clicked.connect(self.__run_calculation)
        
        vbox.addWidget(header)
        vbox.addWidget(form)
        vbox.addWidget(button_run)
        
    def __create_file_dialog(self, dir_only=False, file_only=False):
        home = os.path.expanduser('~')
        if dir_only:
            directory = QFileDialog.getExistingDirectory(self, 'Select directory', home)
            self.base_edit.setText(directory)
        elif file_only:
            h5file = QFileDialog.getOpenFileName(self, 'Select file', home, "Files (*.h5 *.hdf *.txt *.dat)")[0]
            self.file_edit.setText(h5file)
        else:
            return 

    def __select_base_dir(self, old, new):
        if self.base_dir is None or "":
            directory = self.__create_file_dialog(dir_only=True)
            self.base_dir = directory
            self.base_edit.setText(directory)
        else:
            self.base_dir = self.base_edit.text() 
    
    def __base_dir_changed(self, text):
        if self.base_dir == "":
            self.base_dir = None
        else:
            self.base_dir = text
    
    def __pat_imms_changed(self, text):
        self.pat_imms = text

    def __ext_file_changed(self, text):
        self.ext_file = text

    def __log_file_changed(self, text):
        self.log_file = text

    def __head_file_changed(self, text):
        self.head_file = text

    def __mass_changed(self, text):
        self.mass = text

    def __charge_changed(self, text):
        self.charge = text
    
    def __gaussian_changed(self, text):
        self.gaussian = text

    def __populate_results_table(self, ccs):
        print(ccs)
        
    def __file_type_changed(self, text):
        self.ftype = text

    def __file_changed(self, text):
        if not self.project == 'ms':
            self.h5file = text
        else:
            self.file = text

    def __dev_changed(self, text):
        self.dev = text

    def __mz_value_changed(self, text):
        self.mz = text

    def __ld_value_changed(self, text):
        self.ld = text

    def __gas_changed(self, text):
        self.gas = text

    def __minmz_changed(self, text):
        self.minmz = self.minmz_edit.text()

    def __maxmz_changed(self, text):
        self.maxmz = self.maxmz_edit.text()

    def __scale_changed(self, text):
        self.scaling = self.scale_edit.text()

    def __drift_gates_changed(self, text):
        self.drift_gates = self.drift_gates_e.text()

    def __mass_gates_changed(self, text):
        self.mass_gates = self.mass_gates_e.text()
        
    def __xcol_changed(self, text):
        self.xcol = text
        
    def __ycol_changed(self, text):
        self.ycol = text

    def __handle_select(self, state):
        if state > 1:
            self.selecting = True
            return
        self.selecting = False

    def __selection_changed(self, text):
        self.selections = text

    def __check_form_input(self):
        if self.project == 'synapt':
            self.__check_form_input_synapt()
        elif self.project == 'imob':
            self.__check_form_input_imob()
        elif self.project == 'drift':
            self.__check_form_input_drift()
        elif self.project == 'ms':
            self.__check_form_input_ms()


    def __check_form_input_synapt(self):
        try:
            path = os.path.expanduser(self.base_dir)
        except:
            raise InvalidForm("Please provide a valid base directory")

        if self.base_dir is None or self.base_dir == "":
            raise InvalidForm("Base directory input is not valid")

        if not os.path.isdir(path):
            raise NotDirectory("Input for base directory does not correspond to an existing directory") 

        if self.pat_imms == "":
            raise InvalidForm("IMMS file pattern is not valid")

        if self.head_file == "":
            raise InvalidForm("Header file pattern is not valid")

        if self.log_file == "":
            raise InvalidForm("Logfile pattern is not valid")

        if self.ext_file == "":
            raise InvalidForm("Extern file pattern is not valid")

        if self.mass == "":
            raise InvalidForm("Please provide a valid mass")

        if self.charge == "":
            raise InvalidForm("Please provide a valid charge")

    def __check_form_input_imob(self):
        try:
            path = os.path.expanduser(self.h5file)
            if not os.path.isfile(path):
                raise Exception

            self.dev        = float(self.dev)
            self.mz         = float(self.mz)
            self.charge     = float(self.charge)
            self.ld         = float(self.ld)
            self.gaussian   = float(self.gaussian)

            if self.selecting:
                if self.selections != "":
                    self.selections = self.__convert2array(self.selections)
                    if self.selections == 0:
                        raise ValueError
        except ValueError:
            raise InvalidForm("Input for range is not valid")

        except:
            raise InvalidForm("Please provide a valid base directory")

        if self.mz == 0:
            raise InvalidForm

        if self.charge == 0:
            raise InvalidForm

        if self.ld      == 0:
            raise InvalidForm

        if self.gaussian == 0 or self.gaussian > 1:
            raise InvalidForm
        
    def __check_form_input_drift(self):
        try:
            if not os.path.isfile(self.h5file):
                raise Exception
            
            self.minmz          = float(self.minmz)
            self.maxmz          = float(self.maxmz)
            self.scaling        = float(self.scaling)
            
            if self.mass_gates is not None:
                mass_gates = self.mass_gates.split(',')
                mass_gates = [float(x) for x in mass_gates]
                mass_gates = np.asarray(mass_gates)
                mass_gates = np.reshape(mass_gates, (-1, 2))
                
                self.mass_gates = mass_gates
            
            if self.drift_gates is not None:
                drift_gates = self.drift_gates.split(',')
                drift_gates = [float(x) for x in mass_gates]
                drift_gates = np.asarray(drift_gates)
                drift_gates = np.reshape(drift_gates, (-1, 2))
                
                self.drift_gates = drift_gates
                
        except Exception as e:
            raise InvalidForm
    
    def __check_form_input_ms(self):
        try:
            if not os.path.isfile(self.file):
                raise Exception
            
            if re.match("XY", self.ftype):
                if self.xcol is None:
                    raise Exception
                if self.ycol is None:
                    raise Exception
                
                self.xcol = int(self.xcol)
                self.ycol = int(self.ycol)
                
        except Exception:
            raise InvalidForm("")

    def __convert2array(self, input):
        numbers = input.split(',')
        rgx     = re.compile(r"\d+-\d+")
        rgxInt  = re.compile(r"\d+")

        selected = []
        for element in numbers:
            if re.match(rgx, element):
                peaks = element.split('-')
                sub   = list(range(int(peaks[0]), int(peaks[1]) + 1, 1))
                selected += sub
            elif re.match(rgxInt, element):
                peak = int(element)
                selected.append(peak)
        return selected

    def __run_calculation(self):
        try:
            self.__check_form_input()
            if self.project == 'synapt':
                imms_data, coeff, results = synapt.analysis_interface(
                    float(self.mass),
                    int(self.charge),
                    int(self.gaussian),
                    base=self.base_dir
                    )
                self.ccs_synapt_success.emit(imms_data, coeff, results)

            elif self.project == 'imob':
                if not self.selecting:
                    imms_data, coeff, results = imob.analysis_interface(
                        self.h5file,
                        float(self.dev),
                        float(self.mz),
                        float(self.charge),
                        float(self.ld),
                        self.gas,
                        )
                else:
                    imms_data, coeff, results = imob.analysis_interface(
                        self.h5file,
                        float(self.dev),
                        float(self.mz),
                        float(self.charge),
                        float(self.ld),
                        self.gas,
                        traces=self.selections
                    )
                self.ccs_imob_success.emit(imms_data, coeff, results)

            elif self.project == 'drift':
                data = drifttime.analysis_interface(
                    self.h5file,
                    self.minmz,
                    self.maxmz,
                    self.scaling,
                    self.drift_gates,
                    self.mass_gates)
                self.drift_success.emit(data)
            
            elif self.project == 'ms':
                if self.ftype == "XY File":
                    data = FileLoader.getFileAsDataFrame(self.file, cols=[self.xcol - 1, self.ycol - 1])
                    self.ms_success.emit(data)
                elif self.ftype == "XY File with correction":
                    data = FileLoader.getFileAsDataFrame(self.file, cols=[self.xcol - 1, self.ycol - 1])
                    data.iloc[:,0] = data.iloc[:,0].apply(correction)
                    self.ms_success.emit(data)
                else:
                    pass
                
        except InvalidForm as form_error:
            print(form_error)

        except NotDirectory as dir_error:
            print(dir_error)

        except Exception as e:
            print(e)