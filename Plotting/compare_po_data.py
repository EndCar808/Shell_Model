#!/usr/bin/env python    

## Author: Enda Carroll
## Date: Sept 2021
## Info: Script to compare solver results with decaying turbulence papers
#        Solver data

#######################
##  Library Imports  ##
#######################
import numpy as np
import h5py
import sys
import os
from numba import njit
import matplotlib as mpl
# mpl.use('TkAgg') # Use this backend for displaying plots in window
mpl.use('Agg') # Use this backend for writing plots to file
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif']  = 'Computer Modern Roman'
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import getopt
from itertools import zip_longest
import multiprocessing as mprocs
import time as TIME
from subprocess import Popen, PIPE, run
from matplotlib.pyplot import cm
from functions import tc, sim_data, import_data, compute_pdf


###############################
##       FUNCTION DEFS       ##
###############################
def parse_cml(argv):

    """
    Parses command line arguments
    """

    ## Create arguments class
    class cmd_args:

        """
        Class for command line arguments
        """

        def __init__(self, in_dir = None, out_dir = None, in_file = None, num_procs = 20, parallel = False, plotting = False, PO = False):
            self.in_dir         = in_dir
            self.out_dir_phases = out_dir
            self.out_dir_triads = out_dir
            self.in_file        = out_dir
            self.plotting       = plotting
            self.tag = "None"
            self.PO = PO


    ## Initialize class
    cargs = cmd_args()

    try:
        ## Gather command line arguments
        opts, args = getopt.getopt(argv, "i:o:f:t:", ["plot", "po"])
    except:
        print("[" + tc.R + "ERROR" + tc.Rst + "] ---> Incorrect Command Line Arguements.")
        sys.exit()

    ## Parse command line args
    for opt, arg in opts:

        if opt in ['-i']:
            ## Read input directory
            cargs.in_dir = str(arg)
            print("\nInput Folder: " + tc.C + "{}".format(cargs.in_dir) + tc.Rst)

            cargs.out_dir = str(arg)
            print("Output Folder: " + tc.C + "{}".format(cargs.out_dir) + tc.Rst)

        if opt in ['-o']:
            ## Read output directory
            cargs.out_dir = str(arg)
            print("Output Folder: " + tc.C + "{}".format(cargs.out_dir) + tc.Rst)

        elif opt in ['-f']:
            ## Read input directory
            cargs.in_file = str(arg)
            print("Input Post Processing File: " + tc.C + "{}".format(cargs.in_file) + tc.Rst)

        elif opt in ['--plot']:
            ## Read in plotting indicator
            cargs.plotting = True

        elif opt in ['--po']:
            ## Read in plotting indicator
            cargs.PO = True

        elif opt in ['-t']:
            cargs.tag = str(arg)

    return cargs

######################
##       MAIN       ##
######################
if __name__ == '__main__':
    # -------------------------------------
    # # --------- Parse Commnad Line
    # -------------------------------------
    cmdargs = parse_cml(sys.argv[1:])
    if cmdargs.in_file is None:
        method = "default"
        data_file_path = cmdargs.in_dir
    else: 
        method = "file"
        data_file_path = cmdargs.in_dir + cmdargs.in_file
    # -----------------------------------------
    # # --------  Read In data
    # -----------------------------------------
    ## Read in simulation parameters
    sys_vars = sim_data(data_file_path, method)

    ## Read in solver data
    run_data = import_data(data_file_path, sys_vars, method)

    # -----------------------------------------
    # # --------  Load in Sasha's PO Data
    # -----------------------------------------
    Y = np.loadtxt("./Sashas_Work/codes/Updated_Codes/Phase_only_"+ str(10) + "_" + str(1e-3)+ "_" + str(0.4)+".txt")
    Phi = Y[0:sys_vars.N] %(2*np.pi)
    Psi = Y[sys_vars.N:2*sys_vars.N] %(2*np.pi)

    Y = np.loadtxt("./Sashas_Work/codes/Data_Temp/Complex_"+ str(12) + "_" + str(1e-3) + ".txt", dtype=complex)
    U = Y[0:sys_vars.N]
    B = Y[sys_vars.N:2*sys_vars.N]
   
    # -----------------------------------------
    # # --------  Plot Data
    # -----------------------------------------
    # -----------------------------------------
    # # --------  Plot Data
    # -----------------------------------------
    if cmdargs.plotting is True:
        if cmdargs.PO is True:
            ##-------------- Plot Phases Tseries
            i = 5
            fig = plt.figure(figsize = (16, 8))
            gs  = GridSpec(1, 2)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.plot(run_data.phi_n[:, i] % 2.0 * np.pi, '--', label = "Mine")
            ax1.plot(np.transpose(Phi)[:, i], '.-', label = "Sashas")
            ax1.set_ylabel(r"$\phi_n$")
            ax1.set_xlabel(r'$t$')
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.legend()

            ax2 = fig.add_subplot(gs[0, 1])
            ax2.plot(run_data.psi_n[:, i] % 2.0 * np.pi, '--', label = "Mine")
            ax2.plot(np.transpose(Psi)[:, i], '.-', label = "Sashas")
            ax2.set_ylabel(r"$\psi_n$")
            ax2.set_xlabel(r'$t$')
            ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax2.legend()

            plt.savefig(cmdargs.out_dir + "/PO_Compare.png")
            plt.close()

        if cmdargs.PO is False:
            ##-------------- Plot Velocity and Magnetic Tseries
            i = -1
            fig = plt.figure(figsize = (16, 8))
            gs  = GridSpec(2, 2)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.plot(np.real(run_data.u[:, i]), '--', label = "Mine Real")
            ax1.plot(np.real(np.transpose(U)[:, i]), '.-', label = "Sashas Real")
            ax1.set_ylabel(r"$u$")
            ax1.set_xlabel(r'$t$')
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.legend()

            ax2 = fig.add_subplot(gs[0, 1])
            ax2.plot(np.real(run_data.b[:, i]), '--', label = "Mine Real")
            ax2.plot(np.real(np.transpose(B)[:, i]), '.-', label = "Sashas Real")
            ax2.set_ylabel(r"$b$")
            ax2.set_xlabel(r'$t$')
            ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax2.legend()

            ax1 = fig.add_subplot(gs[1, 0])
            ax1.plot(np.imag(run_data.u[:, i]), '--', label = "Mine Imag")
            ax1.plot(np.imag(np.transpose(U)[:, i]), '.-', label = "Sashas Imag")
            ax1.set_ylabel(r"$u$")
            ax1.set_xlabel(r'$t$')
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.legend()

            ax2 = fig.add_subplot(gs[1, 1])
            ax2.plot(np.imag(run_data.b[:, i]), '--', label = "Mine Imag")
            ax2.plot(np.imag(np.transpose(B)[:, i]), '.-', label = "Sashas Imag")
            ax2.set_ylabel(r"$b$")
            ax2.set_xlabel(r'$t$')
            ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax2.legend()

            plt.savefig(cmdargs.out_dir + "/UB_Compare.png")
            plt.close()


            ##-------------- Plot Error Tseries
            i = -1
            fig = plt.figure(figsize = (16, 8))
            gs  = GridSpec(1, 2)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.plot(np.absolute(np.real(run_data.u[:-1, i]) - np.real(np.transpose(U)[:, i])) / np.absolute(np.real(run_data.u[:-1, i])), '--', label = "Mine Real")
            ax1.set_ylabel(r"Error")
            ax1.set_xlabel(r'$t$')
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.legend()

            ax2 = fig.add_subplot(gs[0, 1])
            ax2.plot(np.absolute(np.real(run_data.b[:-1, i]) - np.real(np.transpose(B)[:, i])) / np.absolute(np.real(run_data.b[:-1, i])), '--', label = "Mine Real")
            ax2.set_ylabel(r"Error")
            ax2.set_xlabel(r'$t$')
            ax2.set_yscale('log')
            ax2.set_xscale('log')
            ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax2.legend()

            plt.savefig(cmdargs.out_dir + "/UB_Compare.png")
            plt.close()
