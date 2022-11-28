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
from functions import tc, sim_data, import_data, compute_pdf, import_stats_data, import_sys_msr_data, compute_pdf_from_hist
from plot_functions import plot_anomalous_exponent, plot_str_funcs_with_slope
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

        def __init__(self, in_dir = None, out_dir = None, stats_dir = None, plotting = False):
            self.in_dir         = in_dir
            self.out_dir_stats  = out_dir
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
    # run_data = import_data(data_file_path, sys_vars, method)

    ## Read in stats data
    stats_data = import_stats_data(data_file_path, sys_vars, method)

    ## Read in sys_msr data
    sys_msr_data = import_sys_msr_data(data_file_path, sys_vars, method)

    ## Make output folder for Run Info plots
    cmdargs.out_dir_stats = cmdargs.out_dir + "STATS/"
    if os.path.isdir(cmdargs.out_dir_stats) != True:
        print("Making folder:" + tc.C + " STATS/" + tc.Rst)
        os.mkdir(cmdargs.out_dir_stats)


    # -----------------------------------------
    # # --------  Plot Strucutre Function
    # -----------------------------------------
    ## --------  Structure function with fit
    zeta_p, ns_zeta_p = plot_str_funcs_with_slope(cmdargs.out_dir_stats + "VelStrFunc_Fit.png", sys_msr_data.k, stats_data.vel_str_func, np.arange(4, 16))

    ## --------  Plot Anomalous Exponent
    p = np.arange(2, stats_data.vel_str_func.shape[-1] + 1)
    plot_anomalous_exponent(cmdargs.out_dir_stats + "Vel_Anonalous_Exponent_Zeta_p.png", p, zeta_p[1:], label_str = r"Velocity; Shell Modell")

    if hasattr(stats_data, "mag_str_func"):
        ## --------  Structure function with fit
        mag_zeta_p, ns_zeta_p = plot_str_funcs_with_slope(cmdargs.out_dir_stats + "MagStrFunc_Fit.png", sys_msr_data.k, stats_data.mag_str_func, np.arange(4, 16))

        ## --------  Plot Anomalous Exponent
        p = np.arange(2, stats_data.mag_str_func.shape[-1] + 1)
        plot_anomalous_exponent(cmdargs.out_dir_stats + "Mag_Anonalous_Exponent_Zeta_p.png", p, mag_zeta_p[1:], label_str = r"Magnetic; Shell Modell")


    # -----------------------------------------
    # # --------  Velocity Structure Functions w/ Fit & Anomolous Exponent
    # -----------------------------------------
    ## --------  Structure function with fit
    enrg_flux_zeta_p, ns_zeta_p = plot_str_funcs_with_slope(cmdargs.out_dir_stats + "VelEnergyFluxAbsStrFunc_Fit.png", sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, :, 0], np.arange(4, 16))

    ## --------  Plot Anomalous Exponent
    plot_anomalous_exponent(cmdargs.out_dir_stats + "VelEnergyFluxAbs_Anonalous_Exponent_Zeta_p.png", p, enrg_flux_zeta_p[1:], label_str = r"Velocity Energy Flux; Shell Modell")

    ## --------  Structure function with fit
    hel_flux_zeta_p, ns_zeta_p = plot_str_funcs_with_slope(cmdargs.out_dir_stats + "VelHelicityFluxAbsStrFunc_Fit.png", sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, :, 1], np.arange(4, 16))

    ## --------  Plot Anomalous Exponent
    plot_anomalous_exponent(cmdargs.out_dir_stats + "VelHelicityFluxAbs_Anonalous_Exponent_Zeta_p.png", p, hel_flux_zeta_p[1:], label_str = r"Velocity HelicityFlux; Shell Modell")
    
    # -----------------------------------------
    # # --------  Velocity Anomolous Expoenent Combined
    # -----------------------------------------
    fig   = plt.figure(figsize = (16, 8))
    gs    = GridSpec(1, 1)
    ax1   = fig.add_subplot(gs[0, 0])
    mark_style = ['o','s','^','x','D','p']
    ax1.plot(p, zeta_p[1:], marker = mark_style[0], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Shell Model")
    ax1.plot(p, enrg_flux_zeta_p[1:], marker = mark_style[0], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Energy Flux Shell Model")
    ax1.plot(p, hel_flux_zeta_p[1:], marker = mark_style[0], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Helicity lux Shell Model")
    ax1.plot(p, ns_zeta_p, marker = mark_style[1], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Navier Stokes")
    ax1.plot(p, p / 3, 'k--', label = "K41")
    # ax1.set_xlim(2.0, 6.0)
    ax1.set_xlabel(r"$p$")
    ax1.set_ylabel(r"$\zeta_p$")
    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
    ax1.legend()
    plt.savefig(cmdargs.out_dir_stats + "VelCombined_Anonalous_Exponent_Zeta_p.png", bbox_inches='tight')
    plt.close()



    # -----------------------------------------
    # # --------  Plot PDFs
    # -----------------------------------------
    if hasattr(stats_data, "vel_hist_counts"):
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

        plt.savefig(cmdargs.out_dir_stats + "VelReal_PDF_InOne_C_data.png")
        plt.close()

    if hasattr(stats_data, "mag_hist_counts"):
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([5, 10, 15]):
            pdf, centres = compute_pdf_from_hist(stats_data.mag_hist_counts[i, :], stats_data.mag_hist_ranges[i, :], normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i))
            ax1.set_xlabel(r"$\Re b_n / \langle (\Re b_n)^2 \rangle^{1/2}$")
            ax1.set_ylabel(r"PDF")
            ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ax1.set_yscale('log')
            ax1.legend()

        plt.savefig(cmdargs.out_dir_stats + "VelReal_PDF_InOne_C_data.png")
        plt.close()
