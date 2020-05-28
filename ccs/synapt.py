#! python3
import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import numpy as np
from collections import namedtuple

def normal_distribution(x, *args):
    """
    Accumulated sum of N gaussian functions depending on the number of
    input elements in the list

    params: list
        List of parameters for normal distribution. 
        The normal distribution has three variables, A, mu and sigma.
        For a n-dimensional generalization, provide a list with 3*N parameters.
        N is the number of gaussian functions used for approximation
    """
    # Textbook normal distribution
    trial_func = lambda x, A, mu, sigma: A * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-1 * (x - mu)**2 / ( 2 * sigma ** 2))
    
    #Accumulated sum of values from the different gaussian functions
    accumulated = np.zeros(len(x))
    # converts args to numpy array 
    params = np.asarray(args)
    params = np.reshape(params, (-1, 3))
    
    # generate y-data for multiple gaussians depending on input size
    for set in params:
        y = trial_func(x, set[0], set[1], set[2])
        accumulated += y

    return accumulated 

def optimize_coefficients_gaussian(data, N, maxfev=100000):
    """
    Calculate optimized coefficients based on the number of gaussian functions
    N: int
        Number of gaussian functions used for deconvolution
    data: dict
        Dictionary of dataframes, where each key corresponds to a new voltage

    Returns:
    coeff: dict
        Dictionary of coefficients saved for each CCS key
    """
    coeff = {}

    for key, dataset in data.items():
        x = dataset['Time'].to_numpy()
        y = dataset['Intensity'].to_numpy()
        max = np.max(y)
        median = np.median(x)
        params = [[max, median, 8] for i in range(N)]
    
        coeff[key], cov = curve_fit(normal_distribution, x, y, p0=params, maxfev=maxfev)

    return coeff

def get_data_of_imms_files(files, mass):
    """
    Get data of imms files based on mass provided.
    Attention: Mass has to be preselected currently, so that it corresponds to the maximum
    """
    data = {}
    for key, file in files.items():
        df = pd.read_csv(file, header=None, sep=" ").round(2)
        df.columns = ['Mass', 'Time', 'Intensity', 'None']
        df = df.drop(['None'], axis=1).set_index('Mass')
        select_mass = df.loc[(df.index == mass)]
        data[key] =  select_mass

    return data

def get_voltage_from_ext_files(files):
    get_float = lambda line: float(re.findall('[+, -]*[\d\.\d]+', line)[0])

    voltages = {}
    for key, file in files.items():
        with open(file, encoding="latin-1") as f:
            for line in f:
                if re.match('ADC Pusher Frequency', line):
                    value = get_float(line)
                    pusher_interval = 0.00025 + value/1000

                if re.match('Helium Cell DC', line):
                    helium_cell_dc = get_float(line)

                if re.match('Helium Exit', line):
                    helium_exit = get_float(line)

                if re.match('IMSBias', line):
                    ims_bias = get_float(line)

                if re.match('Transfer DC Entrance', line):
                    transfer_dc_entrance = get_float(line)
            
            voltages[key] = (helium_cell_dc + helium_exit + ims_bias - transfer_dc_entrance) * (1 - 2 / 170.0)

    return voltages, pusher_interval

def get_all_imms_files(base, file):
    """
    Get the directory of all imms files
    dir: string
        Directory in which searching is started
    file:
        Pattern which is used for searching files
    """
    directory = os.path.expanduser(base)
    rgx = re.compile(r'' + file)
    paths = {} 
    for root, dirs, files in os.walk(directory):
        for file in files:
            if re.search(rgx, file):
                # create path to file
                path = os.path.join(root, file)
                # Extract the CCS key
                ccs = re.findall('CCS\d+', file)
                # dictionary with CCS key and path to the corresponding file
                paths[ccs[0]] = path

    return paths

def get_all_ext_and_head_files(dir, ext_file, head_file):
    """
    Get the directory of all ext files
    dir: string
        Directory in which searching is started
    ext_file:
        Pattern which is used for searching external files
    head_file:
        Pattern which is used for searching header files
    """
    dir = os.path.expanduser(dir)
    ext_rgx = re.compile(r'' + ext_file)
    head_rgx = re.compile(r"" + head_file)
    ext_paths = {} 
    head_paths = {}

    for root, dirs, files in os.walk(dir):
        for file in files:
            if re.search(ext_rgx, file):
                # create path to file
                path = os.path.join(root, file)
                # Extract the CCS key
                ccs = re.findall('CCS\d+', root)
                # dictionary with CCS key and path to the corresponding file
                ext_paths[ccs[0]] = os.path.join(root, file)

            if re.search(head_rgx, file):
                # create path to file
                path = os.path.join(root, file)
                # Extract the CCS key
                ccs = re.findall('CCS\d+', root)
                # dictionary with CCS key and path to the corresponding file
                head_paths[ccs[0]] = os.path.join(root, file)

    return ext_paths, head_paths

def get_temp_drift_pressure(dir, log, head):
    """
    Function for extracting data from the log file based on the 
    acquired time saved in the header file

    dir: string
        Base directory for starting search for dir file
    log: string
        Pattern used for matching the log file
    head: dict
        Dictionary of header files from which the acquired time 
        gets extracted

    Returns:
    Temp, Drift gas, Pressure
    """
    dir = os.path.expanduser(dir)
    log_rgx = re.compile(r'' + log)
    log_file = ''
    
    # Find the log file
    for root, dirs, files in os.walk(dir):
        for file in files:
            if re.search(log_rgx, file):
                log_file = os.path.join(root, file)

    # Extract the acquired data from header file
    environment = {}
    for key, file in head.items():
        date_time = extract_acquired_time(file)
        environment[key] = get_data_from_log_file(date_time, log_file)

    return environment

def get_data_from_log_file(date, log):
    """
    Extracts temperature, drift gas and pressure from log file
    """
    date_rgx  = re.compile(r'' + date)
    with open(log) as f:
        for line in f:
            if re.match(date_rgx, line):
                match = re.search('(?P<p>[\d,\.]+) (?P<gas>N2|He) (?P<T>[\d, \.]+)', line) 
                break
    
    return (float(match['p']), match['gas'], 273.15 + float(match['T']))

def get_ccs(mass, charge, voltage, pusher, coeffs, N, env):
    """
    Function for calculating the collision cross section
    mass: float
        Mass of the selected peak
    charge: int
        Charge of the selected peak
    voltage: dict
        Dictionary of calculated voltages for the different CCSs
    pusher: float
        For conversion of drift times
    coeff: dict
        Dictionary of coefficients used for fitting the Peaks
    N: int
        Number of gaussian functions used for fit
    env: dict
        Dictionary of environment variables such as pressure, drift gas and temperature
    """
    # Number of voltage points
    points = len(voltage)
    # drift times 
    """
    --> number of gaussian functions 
    |
    |
    |
    number of voltage points
    """
    drift_times = np.zeros((points, N))
    # voltage of each measurement
    voltages = np.zeros(points)
    # k
    K = np.zeros(N)
    K0 = np.zeros(N)
    # number of ccs values is equal to number of fitted gaussians
    ccs = np.zeros(N)
    # error resulting from ccs fit in percent
    error = np.zeros(N)
    # construct data structure for voltages and coefficients
    for i, (key, V) in enumerate(voltage.items()):
        drift_times[i] = pusher * coeffs[key][1::3]
        voltages[i] = 1 / V 

    
    # each col of the coefficient array will represent a stand-alone ccs calculation
    # iterate over cols
    Experimental = namedtuple('Experimental', ['voltage', 'drift_time'])
    Regress = namedtuple('Regress', ['voltage', 'fitted', 'R', 'P', 'std_err'])
    CCS = namedtuple('CCS', ['ccs', 'error'])
    results = []
    for i in range(drift_times.shape[1]):
        key, (p, gas, T) = env.popitem()

        if gas == 'He':
            drift_mass = 4
        elif gas == 'N2':
            drift_mass = 28
        else:
            sys.exit("Unknown gas")

        slope, intercept, r_value, p_value, std_err = linregress(voltages, drift_times[:,i])
        neutral_mass = mass * charge
        K[i] = (25.05**2)/(slope/1000) 
        K0[i] =  K[i]*(273.15/T)*(p/760)
        reduced_mass = (neutral_mass * drift_mass) / (neutral_mass + drift_mass)
        ccs[i] = (charge/K0[i])*((1/(reduced_mass *T ))**0.5)*18495.88486
        error[i] = (std_err / slope ) * ccs[i]

        exp = Experimental(voltages, drift_times[:,i])
        reg = Regress(voltages, slope*voltages + intercept, r_value, p_value, std_err)
        ccs = CCS(ccs[i], error[i])

        results.append((exp, reg, ccs))

    return results

def extract_acquired_time(file):
    """
    Extract acquired data and time from history file
    """
    date_rgx = re.compile(r'\$\$ Acquired Date')
    time_rgx = re.compile(r'\$\$ Acquired Time')
    date = ''
    time = ''
    with open(file) as f:
        for line in f:
            if re.match(date_rgx, line):
                match = re.search(r'(?P<date>[0-9]+[-, A-Z, a-z]+[0-9]+)', line)
                date = match['date']

            if re.match(time_rgx, line):
                match = re.search(r'(?P<time>[0-9]+\:[0-9]+)', line)
                time = match['time']

    return date + ' ' + time

def visualize_fitting_results(imms_data, coeffs):
    """
    Visualize fitting results based on extracted imms_data
    and coefficients obtained from fitting

    imms_data: dict
        Dictionary of CCS keys and the corresponding data frame

    coeffs: dict
        Dictionary of CCS keys and the corresponding coefficients
    """

    fig = plt.figure()
    for key,  coeff in coeffs.items():
        temp = imms_data[key]
        x = temp['Time'].to_numpy()
        y = temp['Intensity'].to_numpy()
        fitted = normal_distribution(x, coeff)
        plt.plot(x, y)
        plt.plot(x, fitted, label='Fitted')

    plt.legend()
    plt.show()

def analysis_interface(mass, charge, gauss=1, base=".", pat_imms='.*imms.txt',
                ext_file="_extern.inf", head_file="_HEADER.TXT", log_file="Logfile.txt"):

    imms_files = get_all_imms_files(base, pat_imms)
    ext_files, head_files = get_all_ext_and_head_files(base, ext_file, head_file)
    imms_data = get_data_of_imms_files(imms_files, mass)
    # get coefficients
    coeff = optimize_coefficients_gaussian(imms_data, gauss)
    # voltage calculated from data in ext files
    voltage, pusher_interval = get_voltage_from_ext_files(ext_files)
    # get temp, drift_gas and pressure from log file
    env_var = get_temp_drift_pressure(base, log_file, head_files) 
    results = get_ccs(mass, charge, voltage, pusher_interval, coeff, gauss, env_var)
    return imms_data, coeff, results
