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



#################################
##          MISC               ##
#################################
def compute_pdf_from_hist(counts, ranges, normed = False):

    """
    Computes the PDF of input data from a histogram of that data

    Input Parameters:
        counts : array
                - Input array of the histogram counts
        ranges : array
                - Input array of the bin ranges
        normed : bool
                - Indicates if the PDF is to be normalized by the variances or not
    """

    ## Get only non_zero counts
    non_zero_indx = np.argwhere(counts != 0)

    ## Compute the bin info for plotting the pdf
    bin_centres = (ranges[1:] + ranges[:-1]) * 0.5
    bin_width   = ranges[1] - ranges[0]
    bin_centres = bin_centres[non_zero_indx]

    ## Compute the PDF
    bin_counts = counts[non_zero_indx]
    pdf = bin_counts / (np.sum(bin_counts) * bin_width)

    if normed is True:
        ## Compute the variance from the PDF data
        std         = np.sqrt(np.sum(pdf * bin_centres**2 * bin_width))
        pdf         *= std
        bin_centres /= std
        bin_width   /= std

    return pdf, bin_centres


def compute_pdf(data, bin_lims = None, nbins = 1000, normed = False):

    """
    Computes the PDF of input data from a histogram of that data

    Input Parameters:
        data : array
                - Input array of the data used to compute the PDF
        nbins : int
                - Determines the number of bins to use in the histogram
        normed : bool
                - Indicates if the PDF is to be normalized by the variances or not
    """

    ## Compute the histogram of the data
    hist, edges = np.histogram(data, bins = nbins, range = bin_lims)

    ## Compute the bin info for plotting the pdf
    bin_centres = (edges[1:] + edges[:-1]) * 0.5
    bin_width   = edges[1] - edges[0]

    ## Compute the PDF
    pdf = hist / (np.sum(hist) * bin_width)

    if normed is True:
        ## Compute the variance from the PDF data
        std         = np.sqrt(np.sum(pdf * bin_centres**2 * bin_width))
        pdf         *= std
        bin_centres /= std
        bin_width   /= std

    return pdf, bin_centres

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
        def __init__(self, N = 0, nu = 0.0, eta = 0.0, a = 0.0, b = 0.0, k0 = 1.0, eps = 0.5, eps_m = 1.0/3.0, Lambda = 2.0, t0 = 0.0, T = 0.0, ndata = 0, forc = "NONE", forc_k = 0, forc_scale = 1.0, u0 = "N_SCALING", dt = 0.):
            self.N          = int(N)
            self.nu         = float(nu)
            self.eta        = float(eta)
            self.a          = float(a)
            self.b          = float(b)
            self.t0         = float(t0)
            self.dt         = float(dt)
            self.T          = float(T)
            self.ndata      = int(ndata)
            self.u0         = str(u0)
            self.dt         = float(dt)
            self.k0         = float(k0)
            self.Lambda     = float(Lambda)
            self.eps        = float(eps)
            self.eps_m      = float(eps_m)
            self.forc       = str(forc)
            self.forc_k     = int(forc_k)
            self.forc_scale = str(forc_scale)

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

                ## Parse the shell wavenumber prefactor
                if line.startswith('Shell Wavenumber Prefactor'):
                    data.k0 = float(line.split()[-1])

                ## Parse the intershell ratio
                if line.startswith('Intershell Ratio'):
                    data.Lambda = float(line.split()[-1])

                ## Parse the forcing type
                if line.startswith('Forcing Type'):
                    data.forc = str(line.split()[-1])

                ## Parse the forcing wavenumber
                if line.startswith('Forcing Wavenumber'):
                    data.forc_k = int(line.split()[-1])

                ## Parse the forcing scale val
                if line.startswith('Forcing Scale Val'):
                    data.forc_scale = float(line.split()[-1])
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
                data.t0 = float(term.split(',')[0].lstrip('T['))
                data.dt = float(term.split(',')[1])
                data.T  = float(term.split(',')[-1].rstrip(']'))

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
                if 'VelModes' not in list(f.keys()) and ('VelAmps' in list(f.keys()) and 'VelPhases' in list(f.keys())):
                    self.u     = f["VelAmps"][:, :] * np.exp(1j * f["VelPhases"][:, :])
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
                if 'EnergyFlux' in list(f.keys()):
                    self.enrg_flux = f["EnergyFlux"][:, :]
                if 'VelEnergyDiss' in list(f.keys()):
                    self.enrg_diss_u = f["VelEnergyDiss"][:, :]
                if 'VelEnergyInput' in list(f.keys()):
                    self.enrg_input_u = f["VelEnergyInput"][:, :]
                if 'MagEnergyDiss' in list(f.keys()):
                    self.enrg_diss_b = f["MagEnergyDiss"][:, :]
                if 'MagEnergyInput' in list(f.keys()):
                    self.enrg_input_b = f["MagEnergyInput"][:, :]
                if 'Time' in list(f.keys()):
                    self.time  = f["Time"][:]
                if 'k' in list(f.keys()):
                    self.k     = f["k"][:]
                ## Read in spectra
                if 'EnergySpectrum' in list(f.keys()):
                    self.enrg_spect = f["EnergySpectrum"][:]
                if 'DissipationSpectrum' in list(f.keys()):
                    self.diss_spect = f["DissipationSpectrum"][:]
                ## Allocate system measure arrays
                if 'TotalEnergy' in list(f.keys()):
                    self.tot_enrg      = f["TotalEnergy"][:]
                if 'TotalVelocityHelicity' in list(f.keys()):
                    self.tot_hel_u       = f["TotalVelocityHelicity"][:]
                if 'TotalMagneticHelicity' in list(f.keys()):
                    self.tot_hel_b       = f["TotalMagneticHelicity"][:]
                if 'TotalCrossHelicity' in list(f.keys()):
                    self.tot_cross_hel = f["TotalCrossHelicity"][:]
                ## Allocate stats arrays
                if 'RealVelHist_Counts' in list(f.keys()):
                    self.vel_hist_counts = f["RealVelHist_Counts"][:, :]
                if 'RealVelHist_Ranges' in list(f.keys()):
                    self.vel_hist_ranges = f["RealVelHist_Ranges"][:, :]
                if 'VelStats' in list(f.keys()):
                    self.vel_stats = f["VelStats"][:, :]
                if 'MagStats' in list(f.keys()):
                    self.mag_stats = f["MagStats"][:, :]
                if 'StructureFunctionVel' in list(f.keys()):
                    self.vel_str_func = f["StructureFunctionVel"][:, :]
                if 'StructureFunctionVelFlux' in list(f.keys()):
                    self.vel_flux_str_func = f["StructureFunctionVelFlux"][:, :]
                if 'StructureFunctionMag' in list(f.keys()):
                    self.mag_str_func = f["StructureFunctionMag"][:, :]
                if 'StructureFunctionMagFlux' in list(f.keys()):
                    self.mag_flux_str_func = f["StructureFunctionMagFlux"][:, :]


    ## Create instance of data class
    data = SolverData()

    return data



def import_stats_data(input_file, sim_data, method = "default"):

    """
    Reads in stats data from stats HDF5 file.

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
        in_file = input_file + "Stats_HDF_Data.h5"
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
                ## Allocate stats arrays
                if 'RealVelHist_Counts' in list(f.keys()):
                    self.vel_hist_counts = f["RealVelHist_Counts"][:, :]
                if 'RealVelHist_Ranges' in list(f.keys()):
                    self.vel_hist_ranges = f["RealVelHist_Ranges"][:, :]
                if 'VelStats' in list(f.keys()):
                    self.vel_stats = f["VelStats"][:, :]
                if 'MagStats' in list(f.keys()):
                    self.mag_stats = f["MagStats"][:, :]
                if 'StructureFunctionVel' in list(f.keys()):
                    self.vel_str_func = f["StructureFunctionVel"][:, :]
                if 'StructureFunctionVelFlux' in list(f.keys()):
                    self.vel_flux_str_func = f["StructureFunctionVelFlux"][:, :]
                if 'StructureFunctionVelFluxAbs' in list(f.keys()):
                    self.vel_flux_str_func_abs = f["StructureFunctionVelFluxAbs"][:, :]
                if 'StructureFunctionMag' in list(f.keys()):
                    self.mag_str_func = f["StructureFunctionMag"][:, :]
                if 'StructureFunctionMagFlux' in list(f.keys()):
                    self.mag_flux_str_func = f["StructureFunctionMagFlux"][:, :]
                if 'StructureFunctionMagFluxAbs' in list(f.keys()):
                    self.mag_flux_str_func_abs = f["StructureFunctionMagFluxAbs"][:, :]
                ## Read in final state of system
                if 'VelModes' in list(f.keys()):
                    self.u = f["VelModes"][:]
                if 'VelAmps' in list(f.keys()):
                    self.a_n = f["VelAmps"][:]
                if 'VelPhases' in list(f.keys()):
                    self.phi_n = f["VelPhases"][:]
                if 'MagModes' in list(f.keys()):
                    self.b = f["MagModes"][:]
                if 'MagAmps' in list(f.keys()):
                    self.b_n = f["MagAmps"][:]
                if 'MagPhases' in list(f.keys()):
                    self.psi_n = f["MagPhases"][:]



    ## Create instance of data class
    data = SolverData()

    return data


def import_sys_msr_data(input_file, sim_data, method = "default"):

    """
    Reads in system measures from system measures HDF5 file.

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
        in_file = input_file + "System_Measure_HDF_Data.h5"
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
                if 'Time' in list(f.keys()):
                    self.time  = f["Time"][:]
                if 'k' in list(f.keys()):
                    self.k     = f["k"][:]
                ## Allocate system measure arrays
                if 'TotalEnergy' in list(f.keys()):
                    self.tot_enrg      = f["TotalEnergy"][:]
                if 'TotalDissipation' in list(f.keys()):
                    self.tot_diss      = f["TotalDissipation"][:]
                if 'TotalVelocityHelicity' in list(f.keys()):
                    self.tot_hel_u       = f["TotalVelocityHelicity"][:]
                if 'TotalMagneticHelicity' in list(f.keys()):
                    self.tot_hel_b       = f["TotalMagneticHelicity"][:]
                if 'TotalCrossHelicity' in list(f.keys()):
                    self.tot_cross_hel = f["TotalCrossHelicity"][:]
                if 'CharacteristicVel' in list(f.keys()):
                    self.u_charact = f["CharacteristicVel"][:]
                if 'IntegralLengthScale' in list(f.keys()):
                    self.int_scale = f["IntegralLengthScale"][:]
                if 'KolmogorovLengthScale' in list(f.keys()):
                    self.kolm_scale = f["KolmogorovLengthScale"][:]
                if 'ReynoldsNo' in list(f.keys()):
                    self.reynolds_no = f["ReynoldsNo"][:]
                if 'TaylorMicroScale' in list(f.keys()):
                    self.taylor_micro_scale = f["TaylorMicroScale"][:]
                if 'VelocityForcing' in list(f.keys()):
                    self.forcing_u = f["VelocityForcing"][:]

    ## Create instance of data class
    data = SolverData()

    return data