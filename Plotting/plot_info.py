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

    ## Make output folder for Run Info plots
    cmdargs.out_dir_info = cmdargs.out_dir + "RUN_INFO_PLOTS/"
    if os.path.isdir(cmdargs.out_dir_info) != True:
        print("Making folder:" + tc.C + " RUN_INFO_PLOTS/" + tc.Rst)
        os.mkdir(cmdargs.out_dir_info)
    # -----------------------------------------
    # # --------  Compute Post data
    # -----------------------------------------
    u_sqr_av = np.mean(np.absolute(run_data.u)**2, axis = 0)
    for i in range(1, 100, 10):
        # print(run_data.tot_enrg[i], 0.5 * np.sum(np.absolute(run_data.u[i, :])**2), run_data.tot_enrg[i]/np.sum(np.absolute(run_data.u[i, :])**2))
        print(2 * sys_msr_data.tot_enrg[i], np.sum(np.absolute(run_data.u[i, :] * np.conjugate(run_data.u[i, :]))))
    print()
    
    print()
    # -----------------------------------------
    # # --------  Plot Data
    # -----------------------------------------
    if cmdargs.plotting is True:
        ##-------------- Plot System Measures
        fig = plt.figure(figsize = (32, 8))
        gs  = GridSpec(2, 4, hspace = 0.35)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.time, sys_msr_data.tot_diss / sys_msr_data.tot_diss[0] - 1)
        ax1.set_xlabel(r"$t$")
        ax1.set_yscale('log')
        ax1.set_title(r"Relative Dissipation")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2 = fig.add_subplot(gs[1, 0])
        ax2.plot(sys_msr_data.time, sys_msr_data.tot_diss)
        ax2.set_xlabel(r"$t$")
        ax1.set_yscale('log')
        ax2.set_title(r"Total Dissipation")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax3 = fig.add_subplot(gs[0, 1])
        ax3.plot(sys_msr_data.time, sys_msr_data.u_charact)
        ax3.set_xlabel(r"$t$")
        ax3.set_title(r"Characterisitic Velocity")
        ax3.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax4 = fig.add_subplot(gs[0, 2])
        ax4.plot(sys_msr_data.time, sys_msr_data.int_scale)
        ax4.set_xlabel(r"$t$")
        ax4.set_title(r"Integral Length Scale")
        ax4.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax4 = fig.add_subplot(gs[0, 3])
        ax4.plot(sys_msr_data.time, sys_msr_data.taylor_micro_scale)
        ax4.set_xlabel(r"$t$")
        ax4.set_title(r"Taylor Micro Scale")
        ax4.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.plot(sys_msr_data.time, sys_msr_data.reynolds_no)
        ax4.set_xlabel(r"$t$")
        ax4.set_title(r"Reynolds No.")
        ax4.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax4 = fig.add_subplot(gs[1, 2])
        ax4.plot(sys_msr_data.time, sys_msr_data.int_scale / sys_msr_data.u_charact, label = "$\ell / U_{rms}$")
        ax4.plot(sys_msr_data.time, 1 / sys_msr_data.tot_diss, label = "$k / \epsilon$")
        ax4.plot(sys_msr_data.time, sys_msr_data.u_charact / sys_msr_data.tot_diss, label = "$U_{rms} / \epsilon$")
        ax4.plot(sys_msr_data.time, 1 / sys_msr_data.u_charact, label = "$1 / U_{rms}$") ## shouled be 1 / (U k)
        ax4.set_xlabel(r"$t$")
        ax4.set_title(r"Eddy Turn Over Time")
        ax4.set_yscale('log')
        ax4.legend()
        ax4.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax4 = fig.add_subplot(gs[1, 3])
        ax4.plot(sys_msr_data.time, sys_msr_data.int_scale / sys_msr_data.u_charact, label = "$\ell / U_{rms}$")
        ax4.plot(sys_msr_data.time, 1 / sys_msr_data.tot_diss, label = "$k / \epsilon$")
        ax4.plot(sys_msr_data.time, sys_msr_data.u_charact / sys_msr_data.tot_diss, label = "$U_{rms} / \epsilon$")
        ax4.plot(sys_msr_data.time, 1 / sys_msr_data.u_charact, label = "$1 / U_{rms}$")  ## shouled be 1 / (U k)
        ax4.set_xlabel(r"$t$")
        ax4.set_title(r"Eddy Turn Over Time")
        ax4.legend()
        ax4.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)


        plt.savefig(cmdargs.out_dir_info + "System_Measures.png", bbox_inches='tight')
        plt.close()

        ##-------------- Plot Conserved quntities
        fig = plt.figure(figsize = (32, 8))
        if hasattr(run_data, 'b'):
            gs  = GridSpec(2, 4)
        else:
            gs  = GridSpec(2, 2)
        ## Plot the relative energy
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.time, sys_msr_data.tot_enrg / sys_msr_data.tot_enrg[0] - 1)
        ax1.set_xlabel(r"$t$")
        ax1.set_title(r"Relative Energy")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ## Plot the relative helicity
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(sys_msr_data.time, 1 - sys_msr_data.tot_hel_u / sys_msr_data.tot_hel_u[0])
        ax2.set_xlabel(r"$t$")
        ax2.set_title(r"Relative Velocity Helicity")
        ax2.set_yscale('symlog')
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ## Plot the relative energy
        ax1 = fig.add_subplot(gs[1, 0])
        ax1.plot(sys_msr_data.time, sys_msr_data.tot_enrg)
        ax1.set_xlabel(r"$t$")
        ax1.set_title(r"Total Energy")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ## Plot the relative helicity
        ax2 = fig.add_subplot(gs[1, 1])
        ax2.plot(sys_msr_data.time, sys_msr_data.tot_hel_u)
        ax2.set_xlabel(r"$t$")
        ax2.set_title(r"Total Velocity Helicity")
        ax2.set_yscale('symlog')
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        if hasattr(run_data, 'b'):
            ## Plot the relative helicity
            ax2 = fig.add_subplot(gs[0, 2])
            ax2.plot(sys_msr_data.time, 1 - sys_msr_data.tot_hel_b / sys_msr_data.tot_hel_b[0])
            ax2.set_xlabel(r"$t$")
            ax2.set_title(r"Relative Magnetic Helicity")
            ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ## Plot the relative cross helicity
            ax3 = fig.add_subplot(gs[0, 3])
            ax3.plot(sys_msr_data.time, 1 - sys_msr_data.tot_cross_hel / sys_msr_data.tot_cross_hel[0])
            ax3.set_xlabel(r"$t$")
            ax3.set_title(r"Relative Cross Helicity")
            ax3.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ## Plot the relative helicity
            ax2 = fig.add_subplot(gs[1, 2])
            ax2.plot(sys_msr_data.time, sys_msr_data.tot_hel_b)
            ax2.set_xlabel(r"$t$")
            ax2.set_title(r"Total Magnetic Helicity")
            ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ## Plot the relative cross helicity
            ax3 = fig.add_subplot(gs[1, 3])
            ax3.plot(sys_msr_data.time, sys_msr_data.tot_cross_hel)
            ax3.set_xlabel(r"$t$")
            ax3.set_title(r"Total Cross Helicity")
            ax3.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)

        plt.savefig(cmdargs.out_dir_info + "Quadratic_Invariants.png", bbox_inches='tight')
        plt.close()

        ##-------------- Plot the Spectra
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        for i in [0, 100, 500, -2]:
            ax1.plot(sys_msr_data.k, run_data.enrg_spect[i, :], label = "Iter = {}".format(i))
        ax1.set_xlabel("$k_n$")
        ax1.set_ylabel("$E$")
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.legend()
        ax2 = fig.add_subplot(gs[0, 1])
        for i in [0, 100, 500, -2]:
            ax2.plot(sys_msr_data.k, run_data.diss_spect[i, :], label = "Iter = {}".format(i))
        ax2.set_xlabel("$k_n$")
        ax2.set_ylabel("$\epsilon$")
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2.legend()
        plt.savefig(cmdargs.out_dir_info + "Spectra.png", bbox_inches='tight')
        plt.close()



        ##-------------- Plot Amplitudes over time
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, 'b'):
            gs  = GridSpec(2, 2)
        else: 
            gs  = GridSpec(1, 2)
        ## Plot the velocity amplitudes
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(sys_vars.N):
            ax1.plot(sys_msr_data.time, np.absolute(run_data.u[:, i]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$A_n$")
        ## Plot the velocity amplitudes on log scale
        ax1 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N):
            ax1.plot(sys_msr_data.time, np.absolute(run_data.u[:, i]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$A_n$")
        ax1.set_yscale('log')
        if hasattr(run_data, 'b'):
            ## Plot the magnetic amplitudes
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N):
                ax2.plot(sys_msr_data.time, np.absolute(run_data.b[:, i]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$B_n$")
            ## Plot the magnetic amplitudes on log scale
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N):
                ax2.plot(sys_msr_data.time, np.absolute(run_data.b[:, i]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$B_n$")
            ax2.set_yscale('log')

        plt.savefig(cmdargs.out_dir_info + "Amplitudes_Tseries.png", bbox_inches='tight')
        plt.close()


        ##-------------- Plot Phases over time
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, 'b'):
            gs  = GridSpec(2, 2)
        else: 
            gs  = GridSpec(1, 2)        
        ## Plot the velocity amplitudes
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(sys_vars.N//2):
            ax1.plot(sys_msr_data.time, np.angle(run_data.u[:, i]) % 2.0 * np.pi)
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("Low N")
        ## Plot the velocity amplitudes
        ax2.set_title("Low N")
        ax1 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N//2, sys_vars.N):
            ax1.plot(sys_msr_data.time, np.angle(run_data.u[:, i]) % 2.0 * np.pi)
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("High N")

        if hasattr(run_data, 'b'):
            ## Plot the magnetic phases
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N//2):
                ax2.plot(sys_msr_data.time, np.angle(run_data.b[:, i]) % 2.0 * np.pi)
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ## Plot the magnetic phases
            ax2 = fig.add_subplot(gs[1, 1])
            for i in range(sys_vars.N//2, sys_vars.N):
                ax2.plot(sys_msr_data.time, np.angle(run_data.b[:, i]) % 2.0 * np.pi)
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ax2.set_title("High N")

        plt.savefig(cmdargs.out_dir_info + "Phases_Tseries.png", bbox_inches='tight')
        plt.close()


        ##-------------- Plot Phases difference phi_n - phi_n+3 over time
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, 'b'):
            gs  = GridSpec(2, 2)
        else: 
            gs  = GridSpec(1, 2)     
        ## Plot the velocity amplitudes
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(sys_vars.N//2):
            ax1.plot(sys_msr_data.time, np.angle(run_data.u[:, i]) - np.angle(run_data.u[:, i + 3]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("Low N")
        ## Plot the velocity amplitudes
        ax2.set_title("Low N")
        ax1 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N//2, sys_vars.N - 3):
            ax1.plot(sys_msr_data.time, np.angle(run_data.u[:, i]) - np.angle(run_data.u[:, i + 3]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("High N")
        if hasattr(run_data, 'b'):
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N//2):
                ax2.plot(sys_msr_data.time, np.angle(run_data.b[:, i]) - np.angle(run_data.b[:, i + 3]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[1, 1])
            for i in range(sys_vars.N//2, sys_vars.N - 3):
                ax2.plot(sys_msr_data.time, np.angle(run_data.b[:, i]) - np.angle(run_data.b[:, i + 3]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ax2.set_title("High N")

        plt.savefig(cmdargs.out_dir_info + "PhaseDifference_Tseries.png", bbox_inches='tight')
        plt.close()



        ##-------------- Compute triads
        T_ppp = np.zeros((run_data.u.shape[0], run_data.u.shape[1] - 2))
        if hasattr(run_data, 'b'):
            T_pss = np.zeros((run_data.u.shape[0], run_data.u.shape[1] - 2))
            T_sps = np.zeros((run_data.u.shape[0], run_data.u.shape[1] - 2))
            T_ssp  = np.zeros((run_data.u.shape[0], run_data.u.shape[1] - 2))
        for i in range(1, sys_vars.N - 2):
            T_ppp[:, i - 1] = np.mod(np.angle(run_data.u[:, i]) + np.angle(run_data.u[:, i + 1]) + np.angle(run_data.u[:, i + 2]), 2.0 * np.pi)
            if hasattr(run_data, 'b'):
                T_pss[:, i - 1] = np.mod(np.angle(run_data.u[:, i]) + np.angle(run_data.b[:, i + 1]) + np.angle(run_data.b[:, i + 2]), 2.0 * np.pi)
                T_sps[:, i - 1] = np.mod(np.angle(run_data.b[:, i]) + np.angle(run_data.u[:, i + 1]) + np.angle(run_data.b[:, i + 2]), 2.0 * np.pi)
                T_ssp[:, i - 1] = np.mod(np.angle(run_data.b[:, i]) + np.angle(run_data.b[:, i + 1]) + np.angle(run_data.u[:, i + 2]), 2.0 * np.pi)


        ##-------------- Plot Phases difference phi_n - phi_n+3 over time
        fig = plt.figure(figsize = (16, 8))
        if hasattr(run_data, 'b'):
            gs  = GridSpec(2, 2)
        else:
            gs  = GridSpec(1, 1)
        t_i = 4.
        t_f = 4.5
        ## Plot the velocity amplitudes
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(sys_vars.N - 2):
            ax1.plot(sys_msr_data.time, T_ppp[:, i])
        ax1.set_xlabel("$t$")
        ax1.set_title("$\phi_n + \phi_{n + 1} + \phi_{n + 2}$")
        # ax1.set_xlim(t_i, t_f)

        if hasattr(run_data, 'b'):
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[0, 1])
            for i in range(sys_vars.N - 2):
                ax2.plot(sys_msr_data.time, T_pss[:, i])
            ax2.set_xlabel("$t$")
            ax2.set_title("$\phi_n + \psi_{n + 1} + \psi_{n + 2}$")
            # ax2.set_xlim(t_i, t_f)
            ## Plot the velocity amplitudes
            ax1 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N - 2):
                ax1.plot(sys_msr_data.time, T_sps[:, i])
            ax1.set_xlabel("$t$")
            ax1.set_title("$\psi_n + \phi_{n + 1} + \psi_{n + 2}$")
            # ax1.set_xlim(t_i, t_f)
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[1, 1])
            for i in range(sys_vars.N - 2):
                ax2.plot(sys_msr_data.time, T_ssp[:, i])
            ax2.set_xlabel("$t$")
            ax2.set_title("$\psi_n + \psi_{n + 1} + \phi_{n + 2}$")
            # ax2.set_xlim(t_i, t_f)

        plt.savefig(cmdargs.out_dir_info + "Triads_Tseries.png", bbox_inches='tight')
        plt.close()

        ##-------------- Plot Triad T_ppp over time
        fig = plt.figure(figsize = (16, 16))
        gs  = GridSpec(5, 4)
        ## Plot the velocity amplitudes
        for i in range(5):
            for j in range(4):
                if i * 4 + j < 17:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf(T_ppp[:, i * 4 + j], bin_lims = (0.0, 2.0 * np.pi), normed = False)
                    ax1.plot(centres, pdf, label = "$Tppp({})$".format(i * 4 + j))
                    ax1.set_xlabel("$\phi_n + \phi_{n + 1} + \phi_{n + 2}$")
                    ax1.set_ylabel("PDF")
                    # ax1.set_yscale("log")
                    ax1.legend()
                    ax1.set_xlim(0, 2.0*np.pi)

        plt.savefig(cmdargs.out_dir_info + "TriadPDF_T_ppp.png", bbox_inches='tight')
        plt.close()

        if hasattr(run_data, 'b'):
            ##-------------- Plot Triad T_pss over time
            fig = plt.figure(figsize = (16, 16))
            gs  = GridSpec(5, 4)
            ## Plot the velocity amplitudes
            for i in range(5):
                for j in range(4):
                    if i * 4 + j < 17:
                        ax1 = fig.add_subplot(gs[i, j])
                        pdf, centres = compute_pdf(T_pss[:, i * 4 + j], bin_lims = (0.0, 2.0 * np.pi), normed = False)
                        ax1.plot(centres, pdf, label = "$Tpss({})$".format(i * 4 + j))
                        ax1.set_xlabel("$\phi_n + \psi_{n + 1} + \psi_{n + 2}$")
                        ax1.set_ylabel("PDF")
                        # ax1.set_yscale("log")
                        ax1.legend()
                        ax1.set_xlim(0, 2.0*np.pi)

            plt.savefig(cmdargs.out_dir_info + "TriadPDF_T_pss.png", bbox_inches='tight')
            plt.close()

            ##-------------- Plot Triad T_sps over time
            fig = plt.figure(figsize = (16, 16))
            gs  = GridSpec(5, 4)
            ## Plot the velocity amplitudes
            for i in range(5):
                for j in range(4):
                    if i * 4 + j < 17:
                        ax1 = fig.add_subplot(gs[i, j])
                        pdf, centres = compute_pdf(T_sps[:, i * 4 + j], bin_lims = (0.0, 2.0 * np.pi), normed = False)
                        ax1.plot(centres, pdf, label = "$Tsps({})$".format(i * 4 + j))
                        ax1.set_xlabel("$\psi_n + \phi_{n + 1} + \psi_{n + 2}$")
                        ax1.set_ylabel("PDF")
                        # ax1.set_yscale("log")
                        ax1.legend()
                        ax1.set_xlim(0, 2.0*np.pi)

            plt.savefig(cmdargs.out_dir_info + "TriadPDF_T_sps.png", bbox_inches='tight')
            plt.close()

            ##-------------- Plot Triad T_ssp over time
            fig = plt.figure(figsize = (16, 16))
            gs  = GridSpec(5, 4)
            ## Plot the velocity amplitudes
            for i in range(5):
                for j in range(4):
                    if i * 4 + j < 17:
                        ax1 = fig.add_subplot(gs[i, j])
                        pdf, centres = compute_pdf(T_ssp[:, i * 4 + j], bin_lims = (0.0, 2.0 * np.pi), normed = False)
                        ax1.plot(centres, pdf, label = "$Tssp({})$".format(i * 4 + j))
                        ax1.set_xlabel("$\psi_n + \psi_{n + 1} + \phi_{n + 2}$")
                        ax1.set_ylabel("PDF")
                        # ax1.set_yscale("log")
                        ax1.legend()
                        ax1.set_xlim(0, 2.0*np.pi)

            plt.savefig(cmdargs.out_dir_info + "TriadPDF_T_ssp.png", bbox_inches='tight')
            plt.close()



        ##-------------- Plot The Total Flux & Diss over time
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 3)
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(sys_vars.N):
            ax1.plot(sys_msr_data.time, run_data.enrg_flux[:, i])
        ax1.set_xlabel("$t$")
        ax1.set_title("Energ Flux")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N):
            ax2.plot(sys_msr_data.time, run_data.enrg_diss_u[:, i])
        ax2.set_xlabel("$t$")
        ax2.set_title("Velocity Energy Dissipation")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1 = fig.add_subplot(gs[0, 2])
        for i in range(sys_vars.N):
            ax1.plot(sys_msr_data.time, run_data.enrg_input_u[:, i])
        ax1.set_xlabel("$t$")
        ax1.set_title("Velocity Energ Input")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)

        plt.savefig(cmdargs.out_dir_info + "EnergyFlux_EnergyDiss_EnergyInput_Tseries.png", bbox_inches='tight')
        plt.close()


        ##-------------- Plot The Total Flux & Diss over time
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(1, 3)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.time, np.sum(run_data.enrg_flux[:, :], axis = -1))
        ax1.set_xlabel("$t$")
        ax1.set_title("Energ Flux")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(sys_msr_data.time, np.sum(run_data.enrg_diss_u[:, :], axis = -1))
        ax2.set_xlabel("$t$")
        ax2.set_title("Velocity Energy Dissipation")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1 = fig.add_subplot(gs[0, 2])
        ax1.plot(sys_msr_data.time, np.sum(run_data.enrg_input_u[:, :], axis = -1))
        ax1.set_xlabel("$t$")
        ax1.set_title("Velocity Energ Input")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)

        print(run_data.enrg_input_u[:, 0])
        print(run_data.enrg_diss_u[:, 0])
        print(run_data.enrg_flux[:, 0])

        plt.savefig(cmdargs.out_dir_info + "Totals_EnergyFlux_EnergyDiss_EnergyInput_Tseries.png", bbox_inches='tight')
        plt.close()
