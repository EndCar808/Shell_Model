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
from functions import tc, sim_data, import_data, compute_pdf, import_stats_data, import_sys_msr_data


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
    # sys_msr_data = import_sys_msr_data(data_file_path, sys_vars, method)

    ## Make output folder for Run Info plots
    cmdargs.out_dir_stats = cmdargs.out_dir + "STATS/"
    if os.path.isdir(cmdargs.out_dir_stats) != True:
        print("Making folder:" + tc.C + " STATS/" + tc.Rst)
        os.mkdir(cmdargs.out_dir_stats)

    # k = sys_msr_data.k
    k = 0.05 * (2**np.arange(1, 26))

    indx_shift = 1

    inert_lim_low  = 4
    inert_lim_high = 18

    fig   = plt.figure(figsize = (16, 8))
    gs    = GridSpec(1, 1)
    ax1   = fig.add_subplot(gs[0, 0])
    x0     = 0.15 
    y0     = 0.15
    width  = 0.3
    height = 0.2
    ax1in = fig.add_axes([x0, y0, width, height])
    for i in range(stats_data.vel_str_func.shape[-1]):
        p, = ax1.plot(np.log2(k), np.log2(stats_data.vel_str_func[:, i]) + i * 10, '.-', label = "$p = {}$".format(i + indx_shift))
        pfit_info  = np.polyfit(np.log2(k[inert_lim_low:inert_lim_high]), np.log2(stats_data.vel_str_func[inert_lim_low:inert_lim_high, i]) + i * 10, 1)
        pfit_slope = pfit_info[0]
        pfit_c     = pfit_info[1]
        print(i +indx_shift, -(i +indx_shift) / 3, pfit_slope, pfit_c)
        ax1.plot(np.log2(k[inert_lim_low:inert_lim_high]), np.log2(k[inert_lim_low:inert_lim_high])*pfit_slope + pfit_c + 1.5, '--', color = p.get_color())
        local_deriv = np.gradient(stats_data.vel_str_func[:, i])
        ax1in.plot(np.log2(k), local_deriv, color = p.get_color())
        ax1in.set_ylabel(r"$\zeta_p$")
    ax1.set_xlabel(r"$log_2 (k_n)$")
    ax1.set_ylabel(r"$log_2 (S_p(k_n))$")
    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
    ax1in.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
    ax1.legend()
    ## Add insert
    plt.savefig(cmdargs.out_dir_stats + "VelStrFunc_Fit.png", bbox_inches='tight')
    plt.close()