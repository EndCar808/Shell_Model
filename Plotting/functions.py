#!/usr/bin/env python    

## Author: Enda Carroll
## Date: Sept 20222
## Info: General python functions for analysing
#        Solver data


#######################
##  Library Imports  ##
#######################
import numpy as np
import h5py
import sys
import os
import re
from numba import njit
import pyfftw
from collections.abc import Iterable
from itertools import zip_longest
from subprocess import Popen, PIPE

#################################
## Colour Printing to Terminal ##
#################################
class tc:
    H         = '\033[95m'
    B         = '\033[94m'
    C         = '\033[96m'
    G         = '\033[92m'
    Y         = '\033[93m'
    R         = '\033[91m'
    Rst       = '\033[0m'
    Bold      = '\033[1m'
    Underline = '\033[4m'

#####################################
##       DATA FILE FUNCTIONS       ##
#####################################
def sim_data(input_dir, method = "default"):

    """
    Reads in the system parameters from the simulation provided in SimulationDetails.txt.

    Input Parameters:
        input_dir : string
                    - If method == "defualt" is True then this will be the path to
                    the input folder. if not then this will be the input folder
        method    : string
                    - Determines whether the data is to be read in from a file or
                    from an input folder
    """

    ## Define a data class 
    class SimulationData:

        """
        Class for the system parameters.
        """

        ## Initialize class
        def __init__(self, N = 0, nu = 0.0, eta = 0.0, a = 0.0, b = 0.0, t0 = 0.0, T = 0.0, ndata = 0, u0 = "N_SCALING", dt = 0.):
            self.N         = int(N)
            self.nu        = float(nu)
            self.eta       = float(eta)
            self.a         = float(a)
            self.b         = float(b)
            self.t0        = float(t0)
            self.dt        = float(dt)
            self.T         = float(T)
            self.ndata     = int(ndata)
            self.u0        = str(u0)
            self.dt        = float(dt)


    ## Create instance of class
    data = SimulationData()

    if method == "default":
        ## Read in simulation data from file
        with open(input_dir + "SimulationDetails.txt") as f:

            ## Loop through lines and parse data
            for line in f.readlines():

                ## Parse viscosity
                if line.startswith('Viscosity'):
                    data.nu = float(line.split()[-1])

                ## Parse diffusivity
                if line.startswith('Magnetic Diffusivity'):
                    data.eta = float(line.split()[-1])

                ## Parse velocity spectrum slope
                if line.startswith('Velocity Spect Slope'):
                    data.a = float(line.split()[-1])

                ## Parse magnetic spectrum slope
                if line.startswith('Magnetic Spect Slope'):
                    data.b = float(line.split()[-1])

                ## Parse number of modes/shells
                if line.startswith('Fourier Modes'):
                    data.N = int(line.split()[-1])

                ## Parse the start and end time
                if line.startswith('Time Range'):
                    data.t0 = float(line.split()[-3].lstrip('['))
                    data.T  = float(line.split()[-1].rstrip(']'))

                ## Parse number of saving steps
                if line.startswith('Total Saving Steps'):
                    data.ndata = int(line.split()[-1]) + 1 # plus 1 to include initial condition

                ## Parse the initial condition
                if line.startswith('Initial Conditions'):
                    data.u0 = str(line.split()[-1])

                ## Parse the timestep
                if line.startswith('Finishing Timestep'):
                    data.dt = float(line.split()[-1])
    else:
        for term in input_dir.split('_'):

            ## Parse Viscosity
            if term.startswith("NU"):
                data.nu = float(term.split('[')[-1].rstrip(']'))

            ## Parse Diffusivity
            if term.startswith("ETA"):
                data.eta = float(term.split('[')[-1].rstrip(']'))

            ## Parse velocity spectrum slope
            if term.startswith("ALPHA"):
                data.a = float(term.split('[')[-1].rstrip(']'))

            ## Parse magnetic spectrum slope
            if term.startswith("BETA"):
                data.b = float(term.split('[')[-1].rstrip(']'))

            ## Parse Number of collocation points & Fourier modes
            if term.startswith("N["):
                data.N = int(term.split('[')[-1].rstrip(']'))

            ## Parse Time range
            if term.startswith('T['):
                data.t0 = float(term.split('-')[0].lstrip('T['))
                data.dt = float(term.split('-')[1])
                data.T  = float(term.split('-')[-1].rstrip(']'))

            ## Parse initial condition
            if term.startswith('u0'):
                data.u0 = str(term.split('[')[-1])
            if not term.startswith('u0') and term.endswith('].h5'):
                data.u0 = data.u0 + '_' + str(term.split(']')[0])


    return data


def import_data(input_file, sim_data, method = "default"):

    """
    Reads in run data from main HDF5 file.

    input_dir : string
                - If method == "defualt" is True then this will be the path to
               the input folder. if not then this will be the input folder
    method    : string
                - Determines whether the data is to be read in from a file or
               from an input folder
    sim_data  : class
                - object containing the simulation parameters
    """
    
    ## Depending on the output mmode of the solver the input files will be named differently
    if method == "default":
        in_file = input_file + "Main_HDF_Data.h5"
    else:
        in_file = input_file

    ## Define a data class for the solver data
    class SolverData:

        """
        Class for the run data.
        """
        def __init__(self):
            ## Open file and read in the data
            with h5py.File(in_file, 'r') as f:
                ## Allocate global arrays
                if 'VelModes' in list(f.keys()):
                    self.u     = f["VelModes"][:, :]
                if 'VelAmps' in list(f.keys()):
                    self.a_n   = f["VelAmps"][:, :]
                if 'VelPhases' in list(f.keys()):
                    self.phi_n = f["VelPhases"][:, :]
                if 'MagModes' in list(f.keys()):
                    self.b     = f["MagModes"][:, :]
                if 'MagAmps' in list(f.keys()):
                    self.b_n   = f["MagAmps"][:, :]
                if 'MagPhases' in list(f.keys()):
                    self.psi_n = f["MagPhases"][:, :]
                if 'Time' in list(f.keys()):
                    self.time  = f["Time"][:]
                if 'k' in list(f.keys()):
                    self.k     = f["k"][:]
                ## Allocate system measure arrays
                if 'TotalEnergy' in list(f.keys()):
                    self.tot_enrg      = f["TotalEnergy"][:]
                if 'TotalHelicity' in list(f.keys()):
                    self.tot_hel       = f["TotalHelicity"][:]
                if 'TotalCrossHelicity' in list(f.keys()):
                    self.tot_cross_hel = f["TotalCrossHelicity"][:]

    ## Create instance of data class
    data = SolverData()

    return data