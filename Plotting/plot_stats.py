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

        def __init__(self, in_dir = None, out_dir = None, in_file = None, num_procs = 20, parallel = False, plotting = False, video = False, triads = False, phases = False, triad_type = None, full = False):
            self.in_dir         = in_dir
            self.out_dir_phases = out_dir
            self.out_dir_triads = out_dir
            self.in_file        = out_dir
            self.plotting       = plotting
            self.tag = "None"


    ## Initialize class
    cargs = cmd_args()

    try:
        ## Gather command line arguments
        opts, args = getopt.getopt(argv, "i:o:f:t:", ["plot"])
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
    # # --------  Plot Data
    # -----------------------------------------
    if cmdargs.plotting is True:
        ##-------------- Plot PDFs
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 4)
        for j, i in enumerate([1, 5, 10, 15]):
            ax1 = fig.add_subplot(gs[0, j])
            var = np.sqrt(np.mean(np.real(run_data.u[:, i])**2))
            pdf, centres = compute_pdf(np.real(run_data.u[:, i]) / var, nbins = 100)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i))
            ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
            ax1.set_ylabel(r"PDF")
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.set_yscale('log')
            ax1.legend()

        plt.savefig(cmdargs.out_dir + "/RealVel_PDF.png")
        plt.close()

        if hasattr(run_data, "b"):
            fig = plt.figure(figsize = (16, 8))
            gs  = GridSpec(1, 4)
            for j, i in enumerate([1, 5, 10, 15]):
                ax2 = fig.add_subplot(gs[0, j])
                var = np.sqrt(np.mean(np.real(run_data.b[:, i])**2))
                pdf, centres = compute_pdf(np.real(run_data.b[:, i]) / var, nbins = 100)
                ax2.plot(centres, pdf, label = "$n = {}$".format(i + 1))
                ax2.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
                ax2.set_ylabel(r"PDF")
                ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                ax2.set_yscale('log')
                ax2.legend()

            plt.savefig(cmdargs.out_dir + "/RealMag_PDF.png")
            plt.close()


        ##---------------- Plot Structure funcitons
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, "StructureFunctionMag"):
            gs  = GridSpec(1, 2)
        else:   
            gs  = GridSpec(1, 1)

        indx_shift = 2

        k_n = sys_vars.k0 * (sys_vars.Lambda**(np.arange(1, sys_vars.N + 1)))
        # print(k_n)

        ax1 = fig.add_subplot(gs[0, 0])
        for i in [4 - indx_shift, 6 - indx_shift]:
            ax1.plot(k_n, run_data.vel_str_func[:, i], label = "$p = {}$".format(i + indx_shift))
        ax1.set_xlabel(r"$k_n$")
        ax1.set_ylabel(r"$S_p(k_n)$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_xlim(1, k_n[-1])
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.legend()

        plt.savefig(cmdargs.out_dir + "/VelStrFunc.png")
        plt.close()


        ##---------------- Plot Flux Structure funcitons
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, "StructureFunctionMagFlux"):
            gs  = GridSpec(2, 2)
        else:
            gs  = GridSpec(1, 2)

        indx_shift = 2
        ax1 = fig.add_subplot(gs[0, 0])
        for i in [4 - indx_shift, 6 - indx_shift]:
            ax1.plot(k_n, run_data.vel_flux_str_func[:, i, 0], label = "$p = {}$".format(i + indx_shift))
        ax1.set_xlabel(r"$k_n$")
        ax1.set_ylabel(r"$S_p^{E}(k_n)$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_xlim(1, k_n[-1])
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.legend()


        print(run_data.vel_flux_str_func[:, :, 1])
        ax2 = fig.add_subplot(gs[0, 1])
        indx_shift = 2
        for i in [4 - indx_shift, 6 - indx_shift]:
            ax2.plot(k_n, run_data.vel_flux_str_func[:, i, 1], label = "$p = {}$".format(i + indx_shift))
        ax2.set_xlabel(r"$k_n$")
        ax2.set_ylabel(r"$S_p^{H}(k_n)$")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2.set_xlim(1, k_n[-1])
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.legend()

        plt.savefig(cmdargs.out_dir + "/VelFluxStrFunc.png")
        plt.close()