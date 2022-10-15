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
from functions import tc, sim_data, import_data, compute_pdf, compute_pdf_from_hist, import_stats_data, import_sys_msr_data
np.seterr(divide = 'ignore')


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

        def __init__(self, in_dir = None, out_dir = None, stats_file = None, plotting = False):
            self.in_dir         = in_dir
            self.out_dir_stats  = stats_file
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


    ## Make output folder for plots
    cmdargs.out_dir_stats = cmdargs.out_dir + "STATS_PLOTS/"
    if os.path.isdir(cmdargs.out_dir_stats) != True:
        print("Making folder:" + tc.C + " STATS_PLOTS/" + tc.Rst)
        os.mkdir(cmdargs.out_dir_stats)
    # -----------------------------------------
    # # --------  Read In data
    # -----------------------------------------
    ## Read in simulation parameters
    sys_vars = sim_data(data_file_path, method)

    ## Read in solver data
    run_data = import_data(data_file_path, sys_vars, method)

    ## Read in stats data
    stats_data = import_stats_data(data_file_path, sys_vars, method)

    ## Read in sys_msr data
    sys_msr_data = import_sys_msr_data(data_file_path, sys_vars, method)

    # -----------------------------------------
    # # --------  Plot Data
    # -----------------------------------------
    if cmdargs.plotting is True:
        ##-------------- Plot PDFs
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 4)
        for j, i in enumerate([1, 5, 10, 15]):
            ax1 = fig.add_subplot(gs[0, j])
            pdf, centres = compute_pdf(np.real(run_data.u[:, i]), nbins = 50, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i))
            ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
            ax1.set_ylabel(r"PDF")
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.set_yscale('log')
            ax1.legend()

        plt.savefig(cmdargs.out_dir_stats + "RealVel_PDF.png")
        plt.close()

        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([5, 10, 15]):
            pdf, centres = compute_pdf(np.real(run_data.u[:, i]), nbins = 50, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i))
            ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
            ax1.set_ylabel(r"PDF")
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.set_yscale('log')
            ax1.legend()

        plt.savefig(cmdargs.out_dir_stats + "RealVel_PDF_InOne.png", bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([5, 10, 15]):
            pdf, centres = compute_pdf_from_hist(stats_data.vel_hist_counts[i, :], stats_data.vel_hist_ranges[i, :], normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i))
            ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
            ax1.set_ylabel(r"PDF")
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.set_yscale('log')
            ax1.legend()

        plt.savefig(cmdargs.out_dir_stats + "RealVel_PDF_InOne_C_data.png")
        plt.close()

        if hasattr(run_data, "b"):
            fig = plt.figure(figsize = (16, 8))
            gs  = GridSpec(1, 4)
            for j, i in enumerate([1, 5, 10, 15]):
                ax2 = fig.add_subplot(gs[0, j])
                pdf, centres = compute_pdf(np.real(run_data.b[:, i]), nbins = 50, normed = True)
                ax2.plot(centres, pdf, label = "$n = {}$".format(i + 1))
                ax2.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
                ax2.set_ylabel(r"PDF")
                ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                ax2.set_yscale('log')
                ax2.legend()

            plt.savefig(cmdargs.out_dir_stats + "RealMag_PDF.png", bbox_inches='tight')
            plt.close()


        ##---------------- Plot Structure funcitons
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, "StructureFunctionMag"):
            gs  = GridSpec(1, 2)
        else:   
            gs  = GridSpec(1, 1)

        indx_shift = 1

        ax1 = fig.add_subplot(gs[0, 0])
        for i in [4 - indx_shift, 6 - indx_shift]:
            ax1.plot(np.log2(sys_msr_data.k), np.log2(stats_data.vel_str_func[:, i]), label = "$p = {}$".format(i + indx_shift))
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (S_p(k_n))$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.legend()
        plt.savefig(cmdargs.out_dir_stats + "VelStrFunc.png")
        plt.close()

        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, "StructureFunctionMag"):
            gs  = GridSpec(1, 2)
        else:   
            gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(stats_data.vel_str_func.shape[-1]):
            ax1.plot(np.log2(sys_msr_data.k), np.log2(stats_data.vel_str_func[:, i]), label = "$p = {}$".format(i + indx_shift))
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (S_p(k_n))$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.legend()
        plt.savefig(cmdargs.out_dir_stats + "VelStrFunc_All.png", bbox_inches='tight')
        plt.close()


        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, "StructureFunctionMag"):
            gs  = GridSpec(1, 2)
        else:   
            gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(stats_data.vel_str_func.shape[-1]):
            p, = ax1.plot(np.log2(sys_msr_data.k), np.log2((sys_msr_data.k**((i + indx_shift) / 3.0)) * stats_data.vel_str_func[:, i]) + i * 10, label = "$p = {}$".format(i + indx_shift))
            ax1.plot(np.log2(sys_msr_data.k), np.log2((sys_msr_data.k**(0))) + i * 10, '--', color = p.get_color())
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (k_n^{p/3} S_p(k_n))$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.legend()
        plt.savefig(cmdargs.out_dir_stats + "PreMult_VelStrFunc_All.png", bbox_inches='tight')
        plt.close()


        ##---------------- Plot Flux Structure funcitons
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, "StructureFunctionMagFlux"):
            gs  = GridSpec(2, 2)
        else:
            gs  = GridSpec(1, 2)

        indx_shift = 1

        ax1 = fig.add_subplot(gs[0, 0])
        for i in [4 - indx_shift, 6 - indx_shift]:
            ax1.plot(np.log2(sys_msr_data.k), np.log2(stats_data.vel_flux_str_func[:, i, 0]), label = "$p = {}$".format(i + indx_shift))
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (S_p^{E}(k_n))$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_title("Energy Flux Structure Function")
        ax1.legend()
        ax2 = fig.add_subplot(gs[0, 1])
        for i in [4 - indx_shift, 6 - indx_shift]:
            ax2.plot(np.log2(sys_msr_data.k), np.log2(np.absolute(stats_data.vel_flux_str_func[:, i, 1])), label = "$p = {}$".format(i + indx_shift))
        ax2.set_xlabel(r"$log_2 (k_n)$")
        ax2.set_ylabel(r"$log_2 (|S_p^{H}(k_n)|)$")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2.set_title("Helicity Flux Structure Function")
        ax2.legend()
        plt.savefig(cmdargs.out_dir_stats + "VelFluxStrFunc.png", bbox_inches='tight')
        plt.close()



        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(2, 2, hspace = 0.4)
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            ax1.plot(np.log2(sys_msr_data.k), np.log2(np.absolute(stats_data.vel_flux_str_func[:, i, 0])), label = "$p = {}$".format(i + indx_shift))
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (|S_p^{E}(k_n)|)$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_title("Energy Flux Structure Function")
        ax1.legend()
        print()
        print()
        ax2 = fig.add_subplot(gs[0, 1])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            ax2.plot(np.log2(sys_msr_data.k), np.log2(np.absolute(stats_data.vel_flux_str_func[:, i, 1])), label = "$p = {}$".format(i + indx_shift))
        ax2.set_xlabel(r"$log_2 (k_n)$")
        ax2.set_ylabel(r"$log_2 (|S_p^{H}(k_n)|)$")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2.set_title("Helicity Flux Structure Function")
        ax2.legend()

        ax1 = fig.add_subplot(gs[1, 0])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            ax1.plot(np.log2(sys_msr_data.k), np.log2(stats_data.vel_flux_str_func_abs[:, i, 0]), label = "$p = {}$".format(i + indx_shift))
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (S_p^{E, abs}(k_n))$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_title("Absolute Energy Flux Structure Function")
        ax1.legend()
        ax2 = fig.add_subplot(gs[1, 1])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            ax2.plot(np.log2(sys_msr_data.k), np.log2(stats_data.vel_flux_str_func_abs[:, i, 1]), label = "$p = {}$".format(i + indx_shift))
        ax2.set_xlabel(r"$log_2 (k_n)$")
        ax2.set_ylabel(r"$log_2 (S_p^{H, abs}(k_n))$")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2.set_title("Absolute Helicity Flux Structure Function")
        ax2.legend()
        plt.savefig(cmdargs.out_dir_stats + "VelFluxStrFunc_All.png", bbox_inches='tight')
        plt.close()


        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(2, 2, hspace = 0.4)
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            p, = ax1.plot(np.log2(sys_msr_data.k), np.log2((sys_msr_data.k**((i + indx_shift) / 3.0)) * np.absolute(stats_data.vel_flux_str_func[:, i, 0])) + i * 10, label = "$p = {}$".format(i + indx_shift))
            ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.k**0) + i * 10, '--', color = p.get_color())
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (k_n^{p / 3}|S_p^{E}(k_n)|)$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_title("Energy Flux Structure Function")
        ax1.legend()
        ax2 = fig.add_subplot(gs[0, 1])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            p, = ax2.plot(np.log2(sys_msr_data.k), np.log2((sys_msr_data.k**((i + indx_shift) / 3.0)) * np.absolute(stats_data.vel_flux_str_func[:, i, 1])) + i * 10, label = "$p = {}$".format(i + indx_shift))
            ax2.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.k**0) + i * 10, '--', color = p.get_color())
        ax2.set_xlabel(r"$log_2 (k_n)$")
        ax2.set_ylabel(r"$log_2 (k_n^{p / 3}|S_p^{H}(k_n))|$")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2.set_title("Helicity Flux Structure Function")
        ax2.legend()
        ax1 = fig.add_subplot(gs[1, 0])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            p, = ax1.plot(np.log2(sys_msr_data.k), np.log2((sys_msr_data.k**((i + indx_shift) / 3.0)) * stats_data.vel_flux_str_func_abs[:, i, 0]) + i * 10, label = "$p = {}$".format(i + indx_shift))
            ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.k**0) + i * 10, '--', color = p.get_color())
        ax1.set_xlabel(r"$log_2 (k_n)$")
        ax1.set_ylabel(r"$log_2 (k_n^{p / 3}S_p^{E, abs}(k_n))$")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_title("Absolute Energy Flux Structure Function")
        ax1.legend()
        ax2 = fig.add_subplot(gs[1, 1])
        for i in range(stats_data.vel_flux_str_func.shape[1]):
            p, = ax2.plot(np.log2(sys_msr_data.k), np.log2((sys_msr_data.k**((i + indx_shift) / 3.0)) * stats_data.vel_flux_str_func_abs[:, i, 1]) + i * 10, label = "$p = {}$".format(i + indx_shift))
            ax2.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.k**0) + i * 10, '--', color = p.get_color())
        ax2.set_xlabel(r"$log_2 (k_n)$")
        ax2.set_ylabel(r"$log_2 (k_n^{p / 3}S_p^{H, abs}(k_n))$")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2.set_title("Absolute Helicity Flux Structure Function")
        ax2.legend()
        plt.savefig(cmdargs.out_dir_stats + "PreMult_VelFluxStrFunc_All.png", bbox_inches='tight')
        plt.close()


        # with h5py.File(cmdargs.out_dir + "/Stats_HDF_Data.hdf5", "w") as f:
        #     f.create_dataset("RealVelHist_Counts", data = stats_data.vel_hist_counts)
        #     f.create_dataset("RealVelHist_Ranges", data = stats_data.vel_hist_ranges)
        #     f.create_dataset("VelStats", data = stats_data.vel_stats)
        #     f.create_dataset("StructureFunctionVel", data = stats_data.vel_str_func)
        #     f.create_dataset("StructureFunctionVelFlux", data = stats_data.vel_flux_str_func)
        #     