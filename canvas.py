from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
from ccs.drifttime import get_sub_atd, get_sub_ms

gray                                = '#757575'
mpl.rcParams['text.color']          = gray 
mpl.rcParams['text.color']          = gray
mpl.rcParams['axes.labelcolor']     = gray 
mpl.rcParams['xtick.color']         = gray 
mpl.rcParams['ytick.color']         = gray 
mpl.rcParams['axes.edgecolor']      = gray

class Canvas(QFrame):
    def __init__(self, project='ccs'):
        super().__init__()

        if project == 'ccs':
            self.__init_ui_ccs()
        elif project == 'drift':
            self.__init_ui_drift()

        self.setStyleSheet(open("QSS/canvas.qss").read())
    
    def __init_ui_ccs(self):
        vbox = QVBoxLayout(self)
        self.figure = Figure(constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.axes   = self.canvas.figure.subplots(nrows=2, ncols=1)
        self.__set_default_plot_settings_ccs()

        vbox.addWidget(self.canvas)

    def __init_ui_drift(self):
        self.vbox = QVBoxLayout(self)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)

        grid        = gridspec.GridSpec(ncols=3, nrows=3, figure=self.figure) 
        top_ax      = self.figure.add_subplot(grid[0, 0:2])
        middle_ax   = self.figure.add_subplot(grid[1:, 0:2])
        right_ax    = self.figure.add_subplot(grid[1:, 2])

        self.axes   = [top_ax, middle_ax, right_ax]

        self.__set_default_plot_settings_drift()
        self.vbox.addWidget(self.canvas)

    def __set_default_plot_settings_drift(self):
        for ax in self.axes:
            ax.patch.set_facecolor('white')
            ax.patch.set_alpha(0)

        self.figure.patch.set_facecolor("None")
        self.canvas.setStyleSheet("background-color:transparent")

    def __set_default_plot_settings_ccs(self):
        for ax in self.axes:
            ax.patch.set_facecolor('white')
            ax.patch.set_alpha(0)

        self.figure.patch.set_facecolor("None")
        self.canvas.setStyleSheet("background-color:transparent")
        self.axes[0].set_xlabel("Drift time [ms]")
        self.axes[0].set_ylabel("Intensity [a.u]")
        self.axes[0].set_title("Extracted Peaks")

        self.axes[1].set_xlabel("$ Drift Voltage^{-1} [V]$")
        self.axes[1].set_ylabel("Drift time [ms]")

    def clear_plot(self, method='ccs'):
        self.axes[0].cla()
        self.axes[1].cla()
        
        if method == 'ccs':
            self.__set_default_plot_settings_ccs()
        else:
            try:
                self.colorbar.remove()
            except Exception:
                pass
            
            self.axes[2].cla()
            self.__set_default_plot_settings_drift()

    def plot_fit(self, x, y, method='full'):
        if method == "full":
            self.axes[0].plot(x,y, gid='exp')

        elif method == 'dashed':
            self.axes[0].plot(x,y, '--', gid='fit')
        
        self.canvas.draw()

    def plot_ccs(self, x, y, method="plot"):
        if method == "plot":
            self.axes[1].plot(x,y, gid='reg')

        elif method == "scatter":
            self.axes[1].scatter(x,y, gid='exp')

        self.canvas.draw()
        
    def plot_drift_time_scope(self, res):
        mass, sliced, traces, acc, atd, scale, mgates, dgates = res
        nrows, ncols = traces.shape
        traces.insert(ncols, 'Mass', mass)
        minmass = sliced[0]
        maxmass = sliced[-1]
        indices = np.where(np.logical_and(mass > minmass, mass < maxmass))
        traces = traces.iloc[indices[0], :]

        sub_ms = []
        sub_atd = []
        if mgates is not None:
            sub_ms = get_sub_ms(traces, mgates)

        if dgates is not None:
            sub_atd = get_sub_atd(atd, dgates)

        traces = traces.set_index('Mass')
        masses = traces.index
        
        time = atd['Time']
        time = np.asarray([time for i in range(len(traces))])
        mass = np.asarray([np.repeat(i, ncols) for i in masses])
        count = traces.to_numpy()
        absmax = np.amax(count)
        
        mass_data = traces.sum(axis=1)

        mass_spectrum, heat_plot, atd_plot = tuple(self.axes)
        
        cmap = 'viridis'
        im = heat_plot.pcolormesh(mass, time, count, cmap=cmap, vmin=0, vmax=scale * absmax, label='Intensities')
        heat_plot.tick_params(labelright=True, labeltop=True, top=True, right=True, which='both')
        heat_plot.xaxis.set_minor_locator(MultipleLocator(25))
        heat_plot.yaxis.set_minor_locator(MultipleLocator(1))
        heat_plot.set_xlabel('m/z')
        heat_plot.set_ylabel('Drift Time [ms]')
        
        mass_spectrum.plot(mass_data.index, mass_data, color=(0, 0, 0))
        mass_spectrum.xaxis.set_ticks([])
        mass_spectrum.yaxis.set_ticks([])
        mass_spectrum.axis('off')
        
        atd_plot.plot(acc, time[0], color=(0, 0, 0))
        atd_plot.xaxis.set_ticks([])
        atd_plot.yaxis.set_ticks([])
        atd_plot.axis('off')

        for df in sub_ms:
            mass_spectrum.plot(df.index, df)

        for df in sub_atd:
            rgb = [np.random.uniform(0, 1) for i in range(4)]
            atd_plot.plot(df, df.index, color=tuple(rgb))

        self.colorbar = self.figure.colorbar(im, ax=atd_plot, use_gridspec=True)
        self.canvas.draw()
