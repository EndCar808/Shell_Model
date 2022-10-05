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
        ##-------------- Plot Conserved quntities
        fig = plt.figure(figsize = (32, 8))
        if hasattr(run_data, 'b'):
            gs  = GridSpec(1, 4)
        else:
            gs  = GridSpec(1, 2)
        ## Plot the relative energy
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(run_data.time, run_data.tot_enrg / run_data.tot_enrg[0] - 1)
        ax1.set_xlabel(r"$t$")
        ax1.set_title(r"Total Energy")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ## Plot the relative helicity
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(run_data.time, 1 - run_data.tot_hel_u / run_data.tot_hel_u[0])
        ax2.set_xlabel(r"$t$")
        ax2.set_title(r"Total Velocity Helicity")
        ax2.set_yscale('symlog')
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        if hasattr(run_data, 'b'):
            ## Plot the relative helicity
            ax2 = fig.add_subplot(gs[0, 2])
            ax2.plot(run_data.time, 1 - run_data.tot_hel_b / run_data.tot_hel_b[0])
            ax2.set_xlabel(r"$t$")
            ax2.set_title(r"Total Magnetic Helicity")
            ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
            ## Plot the relative cross helicity
            ax3 = fig.add_subplot(gs[0, 3])
            ax3.plot(run_data.time, 1 - run_data.tot_cross_hel / run_data.tot_cross_hel[0])
            ax3.set_xlabel(r"$t$")
            ax3.set_title(r"Total Cross Helicity")
            ax3.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)

        plt.savefig(cmdargs.out_dir + "Quadratic_Invariants.png")
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
            ax1.plot(run_data.time, np.absolute(run_data.u[:, i]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$A_n$")
        ## Plot the velocity amplitudes on log scale
        ax1 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N):
            ax1.plot(run_data.time, np.absolute(run_data.u[:, i]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$A_n$")
        ax1.set_yscale('log')
        if hasattr(run_data, 'b'):
            ## Plot the magnetic amplitudes
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N):
                ax2.plot(run_data.time, np.absolute(run_data.b[:, i]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$B_n$")
            ## Plot the magnetic amplitudes on log scale
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N):
                ax2.plot(run_data.time, np.absolute(run_data.b[:, i]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$B_n$")
            ax2.set_yscale('log')

        plt.savefig(cmdargs.out_dir + "Amplitudes_Tseries.png")
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
            ax1.plot(run_data.time, np.angle(run_data.u[:, i]) % 2.0 * np.pi)
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("Low N")
        ## Plot the velocity amplitudes
        ax2.set_title("Low N")
        ax1 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N//2, sys_vars.N):
            ax1.plot(run_data.time, np.angle(run_data.u[:, i]) % 2.0 * np.pi)
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("High N")

        if hasattr(run_data, 'b'):
            ## Plot the magnetic phases
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N//2):
                ax2.plot(run_data.time, np.angle(run_data.b[:, i]) % 2.0 * np.pi)
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ## Plot the magnetic phases
            ax2 = fig.add_subplot(gs[1, 1])
            for i in range(sys_vars.N//2, sys_vars.N):
                ax2.plot(run_data.time, np.angle(run_data.b[:, i]) % 2.0 * np.pi)
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ax2.set_title("High N")

        plt.savefig(cmdargs.out_dir + "Phases_Tseries.png")
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
            ax1.plot(run_data.time, np.angle(run_data.u[:, i]) - np.angle(run_data.u[:, i + 3]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("Low N")
        ## Plot the velocity amplitudes
        ax2.set_title("Low N")
        ax1 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N//2, sys_vars.N - 3):
            ax1.plot(run_data.time, np.angle(run_data.u[:, i]) - np.angle(run_data.u[:, i + 3]))
        ax1.set_xlabel("$t$")
        ax1.set_ylabel("$\phi_n$")
        ax1.set_title("High N")
        if hasattr(run_data, 'b'):
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N//2):
                ax2.plot(run_data.time, np.angle(run_data.b[:, i]) - np.angle(run_data.b[:, i + 3]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[1, 1])
            for i in range(sys_vars.N//2, sys_vars.N - 3):
                ax2.plot(run_data.time, np.angle(run_data.b[:, i]) - np.angle(run_data.b[:, i + 3]))
            ax2.set_xlabel("$t$")
            ax2.set_ylabel("$\psi_n$")
            ax2.set_title("High N")

        plt.savefig(cmdargs.out_dir + "PhaseDifference_Tseries.png")
        plt.close()



        ##-------------- Compute triads
        T_ppp = np.zeros((run_data.u.shape[0], run_data.u.shape[1] - 3))
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
            ax1.plot(run_data.time, T_ppp[:, i])
        ax1.set_xlabel("$t$")
        ax1.set_title("$\phi_n + \phi_{n + 1} + \phi_{n + 2}$")
        ax1.set_xlim(t_i, t_f)

        if hasattr(run_data, 'b'):
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[0, 1])
            for i in range(sys_vars.N - 2):
                ax2.plot(run_data.time, T_pss[:, i])
            ax2.set_xlabel("$t$")
            ax2.set_title("$\phi_n + \psi_{n + 1} + \psi_{n + 2}$")
            ax2.set_xlim(t_i, t_f)
            ## Plot the velocity amplitudes
            ax1 = fig.add_subplot(gs[1, 0])
            for i in range(sys_vars.N - 2):
                ax1.plot(run_data.time, T_sps[:, i])
            ax1.set_xlabel("$t$")
            ax1.set_title("$\psi_n + \phi_{n + 1} + \psi_{n + 2}$")
            ax1.set_xlim(t_i, t_f)
            ## Plot the velocity amplitudes
            ax2 = fig.add_subplot(gs[1, 1])
            for i in range(sys_vars.N - 2):
                ax2.plot(run_data.time, T_ssp[:, i])
            ax2.set_xlabel("$t$")
            ax2.set_title("$\psi_n + \psi_{n + 1} + \phi_{n + 2}$")
            ax2.set_xlim(t_i, t_f)

        plt.savefig(cmdargs.out_dir + "Triads_Tseries.png")
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

        plt.savefig(cmdargs.out_dir + "TriadPDF_T_ppp.png")
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

            plt.savefig(cmdargs.out_dir + "TriadPDF_T_pss.png")
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

            plt.savefig(cmdargs.out_dir + "TriadPDF_T_sps.png")
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

            plt.savefig(cmdargs.out_dir + "TriadPDF_T_ssp.png")
            plt.close()



        ##-------------- Plot The Total Flux & Diss over time
        fig = plt.figure(figsize = (16, 8))
        gs  = GridSpec(2, 2)
        ## Plot the velocity amplitudes
        ax1 = fig.add_subplot(gs[0, 0])
        for i in range(sys_vars.N):
            ax1.plot(run_data.time, run_data.enrg_flux[:, i])
        ax1.set_xlabel("$t$")
        ax1.set_title("Energ Flux")
        ax1.set_yscale('symlog')
        ax2 = fig.add_subplot(gs[0, 1])
        for i in range(sys_vars.N):
            ax2.plot(run_data.time, run_data.enrg_diss[:, i])
        ax2.set_xlabel("$t$")
        ax2.set_title("Energy Dissipation")
        ax2.set_yscale('symlog')
        ax1 = fig.add_subplot(gs[1, 0])

        for i in range(sys_vars.N):
            ax1.plot(run_data.time, run_data.enrg_flux[:, i])
        ax1.set_xlabel("$t$")
        ax1.set_title("Energ Flux")
        ax2 = fig.add_subplot(gs[1, 1])
        for i in range(sys_vars.N):
            ax2.plot(run_data.time, run_data.enrg_diss[:, i])
        ax2.set_xlabel("$t$")
        ax2.set_title("Energy Dissipation")


        plt.savefig(cmdargs.out_dir + "EnergyFlux_EnergyDiss_Tseries.png")
        plt.close()
