import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import pandas as pd
import re
from collections import namedtuple, OrderedDict
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.signal import lfilter

def optimize_coefficients_gaussian(data, N):
    coeff = {}
    x = data.index
    for data_set in data.iteritems():
        set_name, y = data_set
        y = np.asarray(y)
        A = np.max(y)
        median  = x[np.argmax(y)]
        idx     = np.where(y >= 0.5*A)[0]
        fwhm    = x[idx[-1]] - x[idx[0]]
        sigma   = fwhm / 2*np.sqrt(2*np.log(2))
        off     = np.mean(y)
        params  = [[A, median, sigma, off] for i in range(N)]
        coeff[set_name], cov = curve_fit(normal_distribution, x, y, p0=params)

    coeff = pd.DataFrame.from_dict(coeff).transpose()
    return coeff

def normal_distribution(x, *args):
    trial_func = lambda x, A, mu, sigma, off: A * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-1 * (x - mu)**2 / ( 2 * sigma ** 2)) + off
    accumulated = np.zeros(len(x))
    params = np.asarray(args)
    params = np.reshape(params, (-1, 4))
    
    for set in params:
        y = trial_func(x, set[0], set[1], set[2], set[3])
        accumulated += y

    return accumulated 

def cancel_noise(series):
    n = 15
    b = [1 / n] * n
    a = 1
    series = lfilter(b, a, series)
    return series

def load_hdf_file(path):
    path   = os.path.expanduser(path)
    h5file = h5.File(path, mode='r')
    return h5file

def get_xy_data(h5file):
    df        = pd.DataFrame()
    xdata     = pd.Series(data=calculate_x_points(h5file))

    data = {0: xdata}

    data_sets = h5file["/Traces/PXI6361/Chan_0"]

    rgx = re.compile(r'\d+')
    for i, data_set in enumerate(data_sets):
        selection   = "/Traces/PXI6361/Chan_0/%s" % data_set
        number      = int(re.findall(rgx, data_set)[0])
        ydata       = np.asarray(h5file[selection])
        data[1+ number] = ydata

    ordered = OrderedDict(sorted(data.items()))
    df = pd.DataFrame.from_dict(ordered)
    df = df.rename(columns={0: 'Time'})
    df = df.set_index('Time')
    df = df.apply(np.negative)
    return df

def calculate_x_points(h5file):
    record_length = np.asarray(h5file["/PXI6361_Sets/Horizontal/Record Length"])
    sample_rate   = np.asarray(h5file["/PXI6361_Sets/Horizontal/Sample Rate"])

    return [1000 * time / sample_rate for time in range(record_length)]

def extract_voltage_temp_pressure(h5file, dev=190):
    environment = namedtuple('Environment', ['Voltage', 'Temperature', 'Pressure'])
    h5file['Scanparameter']

    start   = h5file['Scanparameter']['Start']
    stop    = h5file['Scanparameter']['Stop'] 
    step    = h5file['Scanparameter']['Step'] 

    voltage = np.arange(start, stop + step, step) - dev 

    temperature = h5file['/Traces/EMB-PT100/Chan_0/DataSet_0']
    temperature = np.asarray(temperature)[0]
    pressure = h5file['Traces/PR4000/Chan_1/DataSet_0']
    pressure = np.asarray(pressure)[0]

    return environment(voltage, temperature, pressure)

def fit_gaussian(data):
    coeff   = optimize_coefficients_gaussian(data, 1)
    return coeff

def calculate_ccs(data, coeff, environment_conditions, measurement_params):
    Experimental = namedtuple('Experimental', ['voltage', 'drift_time'])
    Regress = namedtuple('Regress', ['voltage', 'fitted', 'slope', 'intercept','R', 'P', 'std_err'])
    CCS = namedtuple('CCS', ['ccs', 'ccs_corr', 'error'])

    voltage, temperature, pressure = environment_conditions
    mz, charge, ld, gas = measurement_params

    try:
        if ( gas == 'He'):
            drift_gas = 4
        elif ( gas == 'N2'):
            drift_gas = 28
        else:
            raise Exception
    except:
        sys.exit("Unknown gas")

    temperature += 273.15
    pressure     = 0.7501 * pressure 

    drift_times = np.asarray(coeff.iloc[:,1])

    if len(drift_times) < len(voltage):
        voltage = voltage[:len(drift_times)]

    inv_voltage =  1 / voltage

    slope, intercept, r_value, p_value, std_err = linregress(inv_voltage, drift_times)
    neutral_mass = mz * charge
    K            = ld**2 / (slope/1000) 
    K0           = K*(273.15/temperature)*(pressure/760)

    reduced_mass = (neutral_mass * drift_gas ) / (neutral_mass + drift_gas)
    ccs          = (charge/K0)*((1/(reduced_mass * temperature ))**0.5)*18495.88486
    ccs_corr     = 0.965 * ccs
    error        = (std_err / slope ) * ccs

    exp = Experimental(inv_voltage, drift_times)
    reg = Regress(inv_voltage, slope*inv_voltage + intercept, slope, intercept, r_value, p_value, std_err)
    ccs = CCS(ccs, ccs_corr, error)

    return [(exp, reg, ccs)]

def visualize_raw_data(data):
    x = data.index

    fig = plt.figure()
    columns = data.columns

    for i in range(data.shape[1]):
        y = data.iloc[:,i]
        plt.plot(x, y, label=columns[i])

    plt.legend() 
    plt.show()

def visualize_gauss_fit(data, coeff):
    time = np.asarray(data.index)
    columns = data.columns

    for i in range(data.shape[1]):
        y = data.iloc[:,i]
        plt.plot(time, y, label=columns[i])
    
    for data_set, row in coeff.iterrows():
        coefficients = np.asarray(row)
        y = normal_distribution(time, coefficients)
        plt.plot(time,y)

    plt.show()

def visualize_reg_fit(x, my, ry):
    plt.plot(x, my)
    plt.scatter(x, ry)
    plt.show()

def analysis_interface(path2file, dev, mz, z, ld, gas, traces=None):
    h5file      = load_hdf_file(path2file)
    data        = get_xy_data(h5file)
    data        = data.apply(cancel_noise)

    if traces is not None:
        try:
            data = data.iloc[:, np.array(traces) - 1]
        except:
            print("Input range is not valid")

    measurement_params = namedtuple('Experimental', ['mz', 'z', 'DriftTube', 'Gas'])
    measurement_params = measurement_params(mz, z, ld, gas)
    environment_conditions = extract_voltage_temp_pressure(h5file, dev)
    coeff   = fit_gaussian(data)
    results = calculate_ccs(data, coeff, environment_conditions, measurement_params)
    return data, coeff, results
