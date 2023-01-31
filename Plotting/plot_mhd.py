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
from functions import tc, sim_data, import_data, compute_pdf, import_stats_data, import_sys_msr_data, import_phase_sync_data, compute_pdf_from_hist, compute_pdf


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

        def __init__(self, in_dir = None, out_dir = None, info_dir = None, plotting = False):
            self.in_dir         = in_dir
            self.out_dir_info   = out_dir
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


    ## Make output folder for plots
    cmdargs.out_dir_MHD = cmdargs.out_dir + "MHD_PLOTS/"
    if os.path.isdir(cmdargs.out_dir_MHD) != True:
        print("Making folder:" + tc.C + " MHD_PLOTS/" + tc.Rst)
        os.mkdir(cmdargs.out_dir_MHD)
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

        ## Plot time averaged amplitude of Z_plus and Z_minus
        if hasattr(run_data, 'z_plus'):
            amp_z_plus_avg      = run_data.z_plus_avg
            amp_z_plus_avg_long = run_data.z_plus_avg
        else:
            ## Get Z_plus
            z_plus = run_data.u[:, :] + run_data.b[:, :]

            ## Get averages
            averaging_indx = sys_vars.ndata // 3
            amp_z_plus_avg      = np.mean(np.absolute(z_plus[:averaging_indx, :]), axis = 0)
            amp_z_plus_avg_long = np.mean(np.absolute(z_plus[:, :]), axis = 0)

        if hasattr(run_data, 'z_minus'):
            amp_z_minus_avg      = run_data.z_minus_avg
            amp_z_minus_avg_long = run_data.z_minus_avg
        else:
            ## Get Z_minus
            z_minus = run_data.u[:, :] - run_data.b[:, :]

            ## Get averages
            averaging_indx = sys_vars.ndata // 3
            amp_z_minus_avg      = np.mean(np.absolute(z_minus[:averaging_indx, :]), axis = 0)
            amp_z_minus_avg_long = np.mean(np.absolute(z_minus[:, :]), axis = 0)
        
        slope = -1/3
        intercept = 0

        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(np.log(sys_msr_data.k), np.log(amp_z_plus_avg), label = "Short Average Window")
        ax1.plot(np.log(sys_msr_data.k), np.log(amp_z_plus_avg_long), 'r', label = "Long Average Window")
        ax1.plot(np.log(sys_msr_data.k), np.log(sys_msr_data.k) * slope + intercept, 'b--', label = "$k^{-a}, a = 1/3$")
        ax1.set_xlabel(r"$\log k_n$")
        ax1.set_ylabel(r"$\log \langle |z^{-}|\rangle$")
        ax1.legend()
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(np.log(sys_msr_data.k), np.log(amp_z_minus_avg), label = "Short Average Window")
        ax2.plot(np.log(sys_msr_data.k), np.log(amp_z_minus_avg), 'r', label = "Long Average Window")
        ax2.plot(np.log(sys_msr_data.k), np.log(sys_msr_data.k) * slope + intercept, 'b--', label = "$k^{-a}, a = 1/3$")
        ax2.set_xlabel(r"$\log k_n$")
        ax2.set_ylabel(r"$\log \langle |z^{+}|\rangle$")
        ax2.legend()
        plt.savefig(cmdargs.out_dir_MHD + "TimeAveraged_Z_Both.png", bbox_inches='tight')
        plt.close()


        ## Plot Reduced Cross Helicity
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.t, run_data.tot_cross_hel / run_data.tot_enrg, label = "Reduced Cross Helicity")
        ax1.set_xlabel(r"Time")
        ax1.set_ylabel(r"$H_{c} / E$")
        ax1.legend()

        plt.savefig(cmdargs.out_dir_MHD + "ReducedCrossHelicity.png", bbox_inches='tight')
        plt.close()




        ## Plot Time averaged nonlinear term -> should scale like k^{-1}
        if not hasattr(run_data, "z_plus"):
            z_plus = run_data.u[:, :] + run_data.b[:, :]
        else:
            z_plus = run_data.z_plus[:, :]
        if not hasattr(run_data, "z_minus"):
            z_minus = run_data.u[:, :] - run_data.b[:, :]
        else:
            z_minus = run_data.z_minus[:, :]

        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.k, np.multiply(np.multiply(z_plus, z_plus), z_minus), label = "Time Averaged Nonlinear Term")
        ax1.plot(sys_msr_data.k, sys_msr_data.k * -1, label = "k^{-a}, a = 1")
        ax1.set_xlabel(r"$k$")
        ax1.set_ylabel(r"$\langle z^{+}z^{+}z^{-} \rangle$")
        ax1.legend()

        plt.savefig(cmdargs.out_dir_MHD + "TimeAveraged_NonlinearTerm.png", bbox_inches='tight')
        plt.close()

        ## Plot The magnetic field structure function
        inert_range = np.arange(2, 13 + 1)
        dns_zeta_p, ns_zeta_p = plot_str_funcs_with_slope(cmdargs.out_dir_MHD + "Magnetic_StrFunc.png", sys_msr_data.k, stats_data.mag_str_func[:, :3], inert_range, insert_fig = False, scaling = 'loge')