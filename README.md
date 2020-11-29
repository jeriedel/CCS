# Collision Cross Section Calculator

GUI application for calculating collision cross section parameters based on the Mason-Shamp equation. Currently, the application supports HDF files and
files outputted from the CDCReader.exe provided by Waters. Additionally, the application provides utility functions for visualizing drift time scans as
contour plots and quick visualization of mass spectra provided as text files in a XY format.

# Running the application 

At the moment, the application is not bundled into a package or provided as a standalone application. For running the application, please make sure that 
all package requirements are met. For installing all requirements, please run the following command 

$ pip3 install -r requirements.txt

To start the application, please execute main.py

$ python3 main.py
