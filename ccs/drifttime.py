import os
import numpy as np
import pandas as pd
import re
import h5py
from collections import OrderedDict, namedtuple
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

def get_sub_atd(atd, dgates):
    sub_atd = []
    time = atd['Time']
    for time_range in dgates:
        mintime = time_range[0]
        maxtime = time_range[1]
        indices = np.where(np.logical_and(time > mintime, time < maxtime))
        df      = atd.iloc[indices[0], :] 
        df      = df.set_index('Time')
        sub_atd.append(df)

    return sub_atd

def get_sub_ms(traces, mgates):
    sub_ms = []
    mass= traces['Mass']
    for mass_range in mgates:
        minmz   = mass_range[0] 
        maxmz   = mass_range[1] 
        indices = np.where(np.logical_and(mass > minmz, mass < maxmz))
        df      = traces.iloc[indices[0], :]
        df      = df.set_index('Mass')
        df      = df.sum(axis=1)
        sub_ms.append(df)

    return sub_ms    

def create_heat_map(mass, sliced, traces, acc, atd, scale=0.1, mgates=None, dgates=None):
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

    time    = atd['Time']
    time    = np.asarray([time for i in range(len(traces))])
    mass    = np.asarray([np.repeat(i, ncols) for i in masses])
    count   = traces.to_numpy()
    absmax = np.amax(count)

    mass_data = traces.sum(axis=1)

    fig             = plt.figure(figsize=(12,7))    
    grid            = GridSpec(ncols=3, nrows=3, figure=fig) 
    mass_spectrum   = fig.add_subplot(grid[0, 0:2])
    heat_plot       = fig.add_subplot(grid[1:, 0:2])
    atd_plot        = fig.add_subplot(grid[1:, 2])

    cmap = 'hot'
    im = heat_plot.pcolormesh(mass, time, count, cmap=cmap, vmin=0, vmax=scale * absmax, label='Intensities')
    heat_plot.tick_params(labelright=True, labeltop=True, top=True, right=True, which='both')
    heat_plot.xaxis.set_minor_locator(MultipleLocator(25))
    heat_plot.yaxis.set_minor_locator(MultipleLocator(1))
    heat_plot.set_xlabel('m/z')
    heat_plot.set_ylabel('Drift Time [ms]')

    mass_spectrum.plot(mass_data.index, mass_data, color=(0,0,0))
    mass_spectrum.xaxis.set_ticks([])
    mass_spectrum.yaxis.set_ticks([])
    mass_spectrum.axis('off')
    atd_plot.plot(acc, time[0], color=(0,0,0))
    atd_plot.xaxis.set_ticks([])
    atd_plot.yaxis.set_ticks([])
    atd_plot.axis('off')

    for df in sub_ms: 
        mass_spectrum.plot(df.index, df)

    for df in sub_atd:
        rgb = [np.random.uniform(0,1) for i in range(4)]
        atd_plot.plot(df, df.index, color=tuple(rgb))

    fig.colorbar(im, ax=atd_plot)
    plt.tight_layout()
    plt.show()

def calculate_atd(h5file, acc):
    scanparams = h5file.get('/Scanparameter')[()]
    start, stop, step, _ = scanparams
    time = [ 1000 * (float(start) + i * float(step)) for i in range(len(acc))] 

    atd = np.column_stack((time, acc))
    atd = pd.DataFrame(data=atd, columns=['Time', 'Accumulated'])
    return atd

def accumulate_traces(traces):
    traces = traces.transpose()
    acc    = traces.sum(axis=1)
    return acc

def load_traces(h5file):
    traces = h5file.get('Traces/PXIe5160/Chan_2')
    d_data = {}
    int_rgx = re.compile('\d+')
    for i, trace in enumerate(traces):
        dataset = h5file.get('Traces/PXIe5160/Chan_2/%s' % trace)[()]
        number  = int(re.findall(int_rgx, trace)[0])
        d_data[number] = dataset

    d_data = OrderedDict(sorted(d_data.items()))
    traces = pd.DataFrame.from_dict(d_data)
    column_names = ['Trace %i' %i for i in range(traces.shape[1])]
    traces.columns = column_names
    traces = traces.apply(np.negative)
    return traces

def create_sliced_array(regress, minmz, maxmz):
    indices = np.where(np.logical_and(regress > minmz, regress < maxmz))
    return regress[indices]

def create_mass_data(h5file, slope, intercept):
    record_length = h5file.get('/PXIe-5160_Sets/Horizontal/Min Record Length')[()]
    sample_rate   = h5file.get('/PXIe-5160_Sets/Horizontal/Min Sample Rate')[()]
    time          = np.asarray([i / float(sample_rate) for i in range(int(record_length))])
    return (slope * time + intercept)**2

def get_calibration_params(h5file):
    slope     = h5file.get('/PXIe-5160_Sets/Mass/Slope2')[()]
    intercept = h5file.get('/PXIe-5160_Sets/Mass/Intercept2')[()]
    return slope, intercept

def load_hdf_file(path2file):
    path2file = os.path.expanduser(path2file)
    return h5py.File(path2file, 'r')

def analysis_interface(path2file, minmz=50, maxmz=3000, scaling=0.1, drift_gates=None, mass_gates=None):
    """
        Path to HDF-file
    minmz: float
        Smallest m/z value displayed
    maxmz: float
        Largest m/z value displayed
    scale: float
        Scaling factor for visualization of plot
    drift_gates: str
        Path to file 
    mass_gates: str
        Path to file
    prc: int
        ?
    """
    h5file              = load_hdf_file(path2file)
    slope, intercept    = get_calibration_params(h5file)
    mass                = create_mass_data(h5file, slope, intercept)
    sliced              = create_sliced_array(mass, minmz, maxmz)
    traces              = load_traces(h5file)
    acc                 = accumulate_traces(traces)
    atd                 = calculate_atd(h5file, acc)
    
    Results = namedtuple('Results', [ 'Masses', 'Sliced', 'Traces', 'Acc', 'ATD', 'Scaling', 'dgates', 'mgates'])
    res = Results(mass, sliced, traces, acc, atd, scaling, drift_gates, mass_gates)
    return res
    
def main():
    analysis_interface(
        '~/Downloads/Drift Time Scan Skript/New Data/newtest.h5',
        minmz=200,
        maxmz=800,
        drift_gates=[[8, 10]],
        mass_gates=[[500, 550]])
    
if __name__ == '__main__':
    main()
