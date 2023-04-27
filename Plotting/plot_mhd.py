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
from functions import tc, sim_data, import_data, compute_pdf, import_stats_data, import_sys_msr_data, import_phase_sync_data, compute_pdf_from_hist, compute_pdf, get_els_flux_field, get_elsassar_sf
from plot_functions import plot_anomalous_exponent, plot_str_funcs_with_slope, phase_only_space_time, plot_str_func_with_anom_scaling

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

    ## Read in sys_msr data
    phase_sync = import_phase_sync_data(data_file_path, sys_vars, method)

    # -----------------------------------------
    # # --------  Plot Data
    # -----------------------------------------
    if cmdargs.plotting is True:
        ## Plotting variables
        fig_size = (12.5, 0.5 * 12.5)
        fig_format = 'png'

        ##------------------------- Get Elsassar variables
        ## Plot time averaged amplitude of Z_plus and Z_minus
        if hasattr(run_data, 'z_plus'):
            amp_z_plus_avg      = run_data.z_plus_avg
        else:
            ## Get Z_plus
            z_plus = run_data.u[:, :] + run_data.b[:, :]
        if hasattr(run_data, 'z_minus'):
            amp_z_minus_avg      = run_data.z_minus_avg
        else:
            ## Get Z_minus
            z_minus = run_data.u[:, :] - run_data.b[:, :]

        # NOTE: Plots that start with "VALID_" are to be compared with the this paper: https://arxiv.org/pdf/nlin/0107032.pdf


        # --------------------------------------------------
        # # --------  Plot Magnetic Structure Function
        # ---------------------------------------------------
        ## Plot The magnetic field structure function
        inert_lim_low  = 3
        inert_lim_high = 12
        inert_range = np.arange(inert_lim_low - 1, (inert_lim_high - 1) + 1) ## -1 to get the correct shell, +1 to include that shell
        dns_zeta_p, ns_zeta_p, dns_zetap_p_res = plot_str_funcs_with_slope(cmdargs.out_dir_MHD + "VALID_Magnetic_StrFunc.png", sys_msr_data.k, stats_data.mag_str_func[:, :3], inert_range, insert_fig = False, scaling = 'loge', fig_size = fig_size)

        # --------------------------------------------------
        # # --------  Plot Time Averaged Elsassar Variables -- Kolmogorov scaling k_n^{-1/3}
        # ---------------------------------------------------     
        # This plot should show the time averaged variables over different windows are the same unlike the
        # non exponentially correlated forcing case where the Kolmogorov is a transient that dies out over long times in Z_n^-
        # but becomes more steep over time in Z_n^+
        slope = -1/3
        intercept = -0.5
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.k, np.mean(np.absolute(z_minus[:sys_vars.ndata//3, :]) , axis = 0), label = r"$\tau = T / 3$")
        ax1.plot(sys_msr_data.k, np.mean(np.absolute(z_minus[:2*sys_vars.ndata//3, :]) , axis = 0), label = r"$\tau = 2 T / 3$")
        ax1.plot(sys_msr_data.k, np.mean(np.absolute(z_minus[:, :]) , axis = 0), label = r"$\tau = T$")
        ax1.plot(sys_msr_data.k[inert_lim_low:inert_lim_high + 1], sys_msr_data.k[inert_lim_low:inert_lim_high + 1] ** slope + intercept, 'k--', label = "$k^{-1/3}$")
        ax1.set_xlabel(r"$k_n$")
        ax1.set_ylabel(r"$\langle |z^{-}|\rangle_\tau$")
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.legend()
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(sys_msr_data.k, np.mean(np.absolute(z_plus[:sys_vars.ndata//3, :]) , axis = 0), label = r"$\tau = T / 3$")
        ax2.plot(sys_msr_data.k, np.mean(np.absolute(z_plus[:2*sys_vars.ndata//3, :]) , axis = 0), label = r"$\tau = 2 T / 3$")
        ax2.plot(sys_msr_data.k, np.mean(np.absolute(z_plus[:, :]) , axis = 0), label = r"$\tau = T$")
        ax2.plot(sys_msr_data.k[inert_lim_low:inert_lim_high + 1], sys_msr_data.k[inert_lim_low:inert_lim_high + 1] ** slope + intercept, 'k--', label = "$k^{-1/3}$")
        ax2.set_xlabel(r"$k_n$")
        ax2.set_ylabel(r"$\langle |z^{+}|\rangle_\tau$")
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.legend()
        plt.savefig(cmdargs.out_dir_MHD + "VALID_TimeAveraged_Z_Both" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        # --------------------------------------------------
        # # --------  Plot Time Averaged Field Variables -- Kolmogorov scaling k_n^{-1/3}
        # ---------------------------------------------------
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.k, np.mean(np.absolute(run_data.u[:, :]), axis=0))
        ax1.plot(sys_msr_data.k[inert_lim_low:inert_lim_high + 1], sys_msr_data.k[inert_lim_low:inert_lim_high + 1] ** slope + intercept, 'k--', label = "$k^{-1/3}$")
        ax1.set_xlabel(r"$ k_n$")
        ax1.set_ylabel(r"$ \langle |u_n| \rangle$")
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.legend()
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(sys_msr_data.k, np.mean(np.absolute(run_data.b[:, :]), axis=0))
        ax2.plot(sys_msr_data.k[inert_lim_low:inert_lim_high + 1], sys_msr_data.k[inert_lim_low:inert_lim_high + 1] ** slope + intercept, 'k--', label = "$k^{-1/3}$")
        ax2.set_xlabel(r"$ k_n$")
        ax2.set_ylabel(r"$ \langle |b_n| \rangle$")
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.legend()
        plt.savefig(cmdargs.out_dir_MHD + "VALID_TimeAveraged_U_and_B" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        # --------------------------------------------------
        # # --------  Plot Tseries of Reduced Cross Helicity
        # ---------------------------------------------------
        ## Plot Reduced Cross Helicity
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.time, sys_msr_data.tot_cross_hel / sys_msr_data.tot_enrg, label = "Reduced Cross Helicity")
        ax1.set_xlabel(r"Time")
        ax1.set_ylabel(r"$\mathcal{H}_{c} / \mathcal{K}$")
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax1.set_ylim(-1.01, 1.01)
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "VALID_ReducedCrossHelicity" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        # --------------------------------------------------
        # # --------  Plot Scaling in Elsassar Variables Flux
        # ---------------------------------------------------
        ## Plot Time averaged nonlinear term -> should scale like k^{-1}
        ## Get mixed transfer term
        t_plus, t_minus = get_els_flux_field(z_plus, z_minus, sys_vars.EPS, sys_vars.EPS_M, sys_vars.Lambda)

        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.k, np.mean(np.absolute(np.imag(t_plus)), axis = 0))
        ax1.plot(sys_msr_data.k[inert_lim_low:inert_lim_high + 1], sys_msr_data.k[inert_lim_low:inert_lim_high + 1] ** (-1) + intercept, 'k--', label = "$k^{-1}$")
        ax1.set_xlabel(r"$k$")
        ax1.set_ylabel(r"$\langle z^{+}z^{+}z^{-} \rangle$")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax1.legend()
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(sys_msr_data.k, np.mean(np.absolute(np.imag(t_minus)), axis = 0))
        ax2.plot(sys_msr_data.k[inert_lim_low:inert_lim_high + 1], sys_msr_data.k[inert_lim_low:inert_lim_high + 1] ** (-1) + intercept, 'k--', label = "$k^{-1}$")
        ax2.set_xlabel(r"$k$")
        ax2.set_ylabel(r"$\langle z^{-}z^{-}z^{+} \rangle$")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax2.legend()
        plt.savefig(cmdargs.out_dir_MHD + "VALID_TimeAveraged_ElsNonlinearTerm" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()


        # print(np.mod(np.angle(t_plus) + 2.0 * np.pi, 2.0 * np.pi), np.mod(np.angle(t_minus) + 2.0 * np.pi, 2.0 * np.pi))


        
        # --------------------------------------------------
        # # --------  Plot System Measures
        # ---------------------------------------------------
        fig = plt.figure(figsize = (32, 8))
        gs  = GridSpec(2, 2, hspace = 0.35)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.time, sys_msr_data.tot_enrg)
        ax1.set_xlabel(r"$t$")
        ax1.set_title(r"Total Energy")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(sys_msr_data.time, sys_msr_data.tot_cross_hel)
        ax2.set_xlabel(r"$t$")
        ax2.set_title(r"Total Cross Helicity")
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax3 = fig.add_subplot(gs[1, 0])
        ax3.plot(sys_msr_data.time, sys_msr_data.tot_hel_b)
        ax3.set_xlabel(r"$t$")
        ax3.set_title(r"Total MagneticHelicity")
        ax3.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.plot(sys_msr_data.time, sys_msr_data.tot_diss_u, label = "$\epsilon_u$")
        if hasattr(run_data, 'b'):
            ax4.plot(sys_msr_data.time, sys_msr_data.tot_diss_b, label = "$\epsilon_b$")
        ax4.set_xlabel(r"$t$")
        ax4.set_yscale('log')
        ax4.set_title(r"Total Dissipation")
        ax4.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax4.legend()
        fig.savefig(cmdargs.out_dir_MHD  + "System_Measures" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()


        # --------------------------------------------------
        # # --------  Plot Field PDFs
        # ---------------------------------------------------
        ###---------- Individual Phases
        print("Plotting --- Phases")
        num_bins = 1000
        norm_hist = False

        input_data = [np.mod(np.angle(run_data.u) + 2.0 * np.pi, 2.0 * np.pi), np.mod(np.angle(run_data.b) + 2.0 * np.pi, 2.0 * np.pi), np.mod(np.angle(z_plus) + 2.0 * np.pi, 2.0 * np.pi), np.mod(np.angle(z_minus) + 2.0 * np.pi, 2.0 * np.pi), np.mod(np.angle(t_plus) + 2.0 * np.pi, 2.0 * np.pi), np.mod(np.angle(t_minus) + 2.0 * np.pi, 2.0 * np.pi)]
        figure_names = [r"U_Phases", r"B_Phases", r"Zplus_Phases", r"Zminus_Phases", r"Tplus_Phases", r"Tminus_Phases"]
        data_labels  = [r"$\arg\{ u_n \}$", r"$\arg\{ b_n \}$", r"$\arg\{ Z_n^+ \}$", r"$\arg\{ Z_n^- \}$", r"$\arg\{ T_n^+ \}$", r"$\arg\{ T_n^- \}$"]

        for in_data, fig_name, data_labs in zip(input_data, figure_names, data_labels):
            fig = plt.figure(figsize=(24, 24))
            gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
            pdf_angles = []
            for i in range(5):
              for j in range(5):
                  indx = i * 5 + j
                  if indx < in_data.shape[-1]:
                      ax1 = fig.add_subplot(gs[i, j])
                      pdf, centres = compute_pdf(in_data, nbins=num_bins, normed=norm_hist)
                      pdf_angles.append(pdf)
                      p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
                      ax1.set_xlabel(data_labs)
                      ax1.set_xlim(0, 2.0*np.pi)
                      ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
                      ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
                      ax1.set_ylabel(r"PDF")
                      ax1.set_yscale('log')
                      ax1.set_title("n = {}".format(indx + 1))
            fig.savefig(cmdargs.out_dir_MHD  + fig_name + "_PDFs" + "." + fig_format, format=fig_format, bbox_inches='tight')
            plt.close()


        ###---------- Individual Amps
        print("Plotting --- Amps")
        num_bins = 1000
        norm_hist = False

        input_data = [np.absolute(run_data.u), np.absolute(run_data.b), np.absolute(z_plus), np.absolute(z_minus), np.absolute(t_plus), np.absolute(t_minus)]
        figure_names = [r"U_Amps", r"B_Amps", r"Zplus_Amps", r"Zminus_Amps", r"Tplus_Amps", r"Tminus_Amps"]
        data_labels  = [r"$| u_n |$", r"$| b_n |$", r"$| Z_n^+ |$", r"$| Z_n^- |$", r"$| T_n^+ |$", r"$| T_n^- |$"]

        for in_data, fig_name, data_labs in zip(input_data, figure_names, data_labels):
            fig = plt.figure(figsize=(24, 24))
            gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
            pdf_angles = []
            for i in range(5):
              for j in range(5):
                  indx = i * 5 + j
                  if indx < in_data.shape[-1]:
                      ax1 = fig.add_subplot(gs[i, j])
                      pdf, centres = compute_pdf(in_data, nbins=num_bins, normed=norm_hist)
                      pdf_angles.append(pdf)
                      p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
                      ax1.set_xlabel(data_labs)
                      ax1.set_ylabel(r"PDF")
                      ax1.set_yscale('log')
                      ax1.set_title("n = {}".format(indx + 1))
            fig.savefig(cmdargs.out_dir_MHD  + fig_name + "_PDFs" + "." + fig_format, format=fig_format, bbox_inches='tight')
            plt.close()


        # --------------------------------------------------
        # # --------  Plot Triad PDFs
        # ---------------------------------------------------
        print("Plotting --- Triads")
        ## Plot the velocity triad pdfs
        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                if i * 5 + j < phase_sync.num_triads:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf_from_hist(phase_sync.vel_triad_hist_counts[i * 5 + j, :], phase_sync.vel_triad_hist_ranges[:], remove_zeros = False)
                    ax1.plot(centres, pdf, label = "$Tppp({})$".format(i * 5 + j + 1))
                    ax1.set_xlabel("$\phi_n + \phi_{n + 1} + \phi_{n + 2}$")
                    ax1.set_ylabel("PDF")
                    ax1.set_yscale("log")
                    ax1.legend()
                    ax1.set_xlim(0, 2.0*np.pi)
                    ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
                    ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
        plt.savefig(cmdargs.out_dir_MHD + "U_Triad_PDF" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        ## Plot the Magnetic triad pdfs
        for tri_type, type_str in enumerate(["Tpss", "Tsps", "Tssp"]):
            fig = plt.figure(figsize=(24, 24))
            gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
            for i in range(5):
                for j in range(5):
                    if i * 5 + j < phase_sync.num_triads:
                        ax1 = fig.add_subplot(gs[i, j])
                        pdf, centres = compute_pdf_from_hist(phase_sync.mag_triad_hist_counts[i * 5 + j, :, tri_type], phase_sync.mag_triad_hist_ranges[:], remove_zeros = False)
                        ax1.plot(centres, pdf, label = r"${}({})$".format(type_str, i * 5 + j + 1))
                        ax1.set_xlabel(type_str)
                        ax1.set_ylabel("PDF")
                        ax1.set_yscale("log")
                        ax1.legend()
                        ax1.set_xlim(0, 2.0*np.pi)
                        ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
                        ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
            plt.savefig(cmdargs.out_dir_MHD + "B_Triad_PDF_" + type_str + "." + fig_format, format=fig_format, bbox_inches='tight')
            plt.close()

        # --------------------------------------------------
        # # --------  Plot Phase Diff
        # ---------------------------------------------------
        print("Plotting --- Phase Diffs")
        ## Plot the velocity phase differences
        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                if i * 5 + j < phase_sync.num_phase_diffs:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf_from_hist(phase_sync.vel_phase_diff_hist_counts[i * 5 + j, :], phase_sync.vel_phase_diff_hist_ranges[:], remove_zeros = False)
                    ax1.plot(centres, pdf, label = "$Diffs({})$".format(i * 5 + j + 1))
                    ax1.set_xlabel("$\phi_n - \phi_{n + 3}$")
                    ax1.set_ylabel("PDF")
                    ax1.set_yscale("log")
                    ax1.legend()
                    ax1.set_xlim(0, 2.0*np.pi)
                    ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
                    ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
        plt.savefig(cmdargs.out_dir_MHD + "U_PhaseDiffs_PDF" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        ## Plot the Magnetic phase differences
        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                if i * 5 + j < phase_sync.num_phase_diffs:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf_from_hist(phase_sync.mag_phase_diff_hist_counts[i * 5 + j, :], phase_sync.mag_phase_diff_hist_ranges[:], remove_zeros = False)
                    ax1.plot(centres, pdf, label = "$Diffs({})$".format(i * 5 + j + 1))
                    ax1.set_xlabel("$\psi_n - \psi_{n + 3}$")
                    ax1.set_ylabel("PDF")
                    ax1.set_yscale("log")
                    ax1.legend()
                    ax1.set_xlim(0, 2.0*np.pi)
                    ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
                    ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
        plt.savefig(cmdargs.out_dir_MHD + "B_PhaseDiffs_PDF" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()


        # --------------------------------------------------
        # # --------  Plot PDFs of Real Part
        # ---------------------------------------------------
        print("Plotting --- Real Part PDFs")
        ##----------------- Plot Velocity
        nbins = 100
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            if hasattr(stats_data, "vel_hist_counts"):
                pdf, centres = compute_pdf_from_hist(stats_data.vel_hist_counts[i, :], stats_data.vel_hist_ranges[i, :], normed = True)
                ax1.set_title("C Data")
            else:
                if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
                    pdf, centres = compute_pdf(np.real(run_data.a_n[:, i] * np.exp(1j * run_data.phi_n[:, i])), nbins = nbins, normed = True)                   
                else:
                    pdf, centres = compute_pdf(np.real(run_data.u[:, i]), nbins = nbins, normed = True)
            if i == -1:
                ax1.plot(centres, pdf, label = "$n = {}$".format(sys_vars.N))    
            else:
                ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1))    
        ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "U_Real_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                indx = i * 5 + j
                if indx < sys_vars.N:
                    ax1 = fig.add_subplot(gs[i, j])
                    if hasattr(stats_data, "vel_hist_counts"):
                        pdf, centres = compute_pdf_from_hist(stats_data.vel_hist_counts[indx, :], stats_data.vel_hist_ranges[indx, :], normed = False)
                        ax1.set_title("C Data")
                    else:
                        if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
                            pdf, centres = compute_pdf(np.real(run_data.a_n[:, indx] * np.exp(1j * run_data.phi_n[:, indx])), nbins = nbins, normed = False)                 
                        else:
                            pdf, centres = compute_pdf(np.real(run_data.u[:, indx]), nbins = nbins, normed = False)
                    ax1.plot(centres, pdf, label = "$n = {}$".format(indx + 1)) 
                    ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
                    ax1.set_ylabel(r"PDF")
                    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                    ax1.set_yscale('log')
                    ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "U_Real_PDF_All" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        ##----------------- Plot Magnetic
        nbins = 100
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            if hasattr(stats_data, "mag_hist_counts"):
                pdf, centres = compute_pdf_from_hist(stats_data.mag_hist_counts[i, :], stats_data.mag_hist_ranges[i, :], normed = True)
                ax1.set_title("C Data")
            else:
                if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
                    pdf, centres = compute_pdf(np.real(run_data.b_n[:, i] * np.exp(1j * run_data.psi_n[:, i])), nbins = nbins, normed = True)                   
                else:
                    pdf, centres = compute_pdf(np.real(run_data.b[:, i]), nbins = nbins, normed = True)
            if i == -1:
                ax1.plot(centres, pdf, label = "$n = {}$".format(sys_vars.N))    
            else:
                ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1))    
        ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "U_Real_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                indx = i * 5 + j
                if indx < sys_vars.N:
                    ax1 = fig.add_subplot(gs[i, j])
                    if hasattr(stats_data, "mag_hist_counts"):
                        pdf, centres = compute_pdf_from_hist(stats_data.mag_hist_counts[indx, :], stats_data.mag_hist_ranges[indx, :], normed = False)
                        ax1.set_title("C Data")
                    else:
                        if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
                            pdf, centres = compute_pdf(np.real(run_data.b_n[:, indx] * np.exp(1j * run_data.psi_n[:, indx])), nbins = nbins, normed = False)                 
                        else:
                            pdf, centres = compute_pdf(np.real(run_data.b[:, indx]), nbins = nbins, normed = False)
                    ax1.plot(centres, pdf, label = "$n = {}$".format(indx + 1)) 
                    ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
                    ax1.set_ylabel(r"PDF")
                    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                    ax1.set_yscale('log')
                    ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "U_Real_PDF_All" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        ##-------------------- Elsassar Variables
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            pdf, centres = compute_pdf(np.real(z_plus[:, i]), nbins = nbins, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1)) 
        ax1.set_xlabel(r"$\Re Z_n^+ / \langle (\Re Z_n^+)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Zplus_Real_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            pdf, centres = compute_pdf(np.real(z_minus[:, i]), nbins = nbins, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1)) 
        ax1.set_xlabel(r"$\Re Z_n^- / \langle (\Re Z_n^-)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Zminus_Real_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            pdf, centres = compute_pdf(np.real(t_plus[:, i]), nbins = nbins, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1)) 
        ax1.set_xlabel(r"$\Re T_n^+ / \langle (\Re T_n^+)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tplus_Real_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            pdf, centres = compute_pdf(np.real(t_minus[:, i]), nbins = nbins, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1)) 
        ax1.set_xlabel(r"$\Re T_n^- / \langle (\Re T_n^-)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tminus_Real_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                indx = i * 5 + j
                if indx < sys_vars.N:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf(np.real(t_plus[:, indx]), nbins = nbins, normed = False)
                    ax1.plot(centres, pdf, label = "$n = {}$".format(indx + 1)) 
                    ax1.set_xlabel(r"$\Re T_n^+ / \langle (\Re T_n^+)^2 \rangle^{1/2}$")
                    ax1.set_ylabel(r"PDF")
                    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                    ax1.set_yscale('log')
                    ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tplus_Real_PDF_All" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                indx = i * 5 + j
                if indx < sys_vars.N:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf(np.real(t_minus[:, indx]), nbins = nbins, normed = False)
                    ax1.plot(centres, pdf, label = "$n = {}$".format(indx + 1)) 
                    ax1.set_xlabel(r"$\Re T_n^- / \langle (\Re T_n^-)^2 \rangle^{1/2}$")
                    ax1.set_ylabel(r"PDF")
                    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                    ax1.set_yscale('log')
                    ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tminus_Real_PDF_All" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            pdf, centres = compute_pdf(np.imag(t_plus[:, i]), nbins = nbins, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1)) 
        ax1.set_xlabel(r"$\Im T_n^+ / \langle (\Im T_n^+)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tplus_Imag_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([6, 11, 16, 20]):
            pdf, centres = compute_pdf(np.imag(t_minus[:, i]), nbins = nbins, normed = True)
            ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1)) 
        ax1.set_xlabel(r"$\Im T_n^- / \langle (\Im T_n^-)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tminus_Imag_PDF_InOne" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                indx = i * 5 + j
                if indx < sys_vars.N:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf(np.imag(t_plus[:, indx]), nbins = nbins, normed = False)
                    ax1.plot(centres, pdf, label = "$n = {}$".format(indx + 1)) 
                    ax1.set_xlabel(r"$\Im T_n^+ / \langle (\Im T_n^+)^2 \rangle^{1/2}$")
                    ax1.set_ylabel(r"PDF")
                    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                    ax1.set_yscale('log')
                    ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tplus_Imag_PDF_All" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(24, 24))
        gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
        for i in range(5):
            for j in range(5):
                indx = i * 5 + j
                if indx < sys_vars.N:
                    ax1 = fig.add_subplot(gs[i, j])
                    pdf, centres = compute_pdf(np.imag(t_minus[:, indx]), nbins = nbins, normed = False)
                    ax1.plot(centres, pdf, label = "$n = {}$".format(indx + 1)) 
                    ax1.set_xlabel(r"$\Im T_n^- / \langle (\Im T_n^-)^2 \rangle^{1/2}$")
                    ax1.set_ylabel(r"PDF")
                    ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
                    ax1.set_yscale('log')
                    ax1.legend()
        plt.savefig(cmdargs.out_dir_MHD + "Tminus_Imag_PDF_All" + "." + fig_format, format=fig_format, bbox_inches='tight')
        plt.close()



        # --------------------------------------------------
        # # --------  Plot Spacetime
        # ---------------------------------------------------
        print("Plotting --- Spacetime")
        ##-------- Velocity
        phase_only_space_time(cmdargs.out_dir_MHD + "U_SpaceTime_Phases" + ".png", np.mod(np.angle(run_data.u) + 2.0*np.pi, 2.0*np.pi), sys_msr_data.time, sys_vars.N, r"$\phi_n$")
        phase_only_space_time(cmdargs.out_dir_MHD + "U_SpaceTime_Triads" + ".png", phase_sync.vel_triads, sys_msr_data.time, phase_sync.num_triads, r"$\varphi_{n, n + 1, n + 2}$")
        phase_only_space_time(cmdargs.out_dir_MHD + "U_SpaceTime_PhaseDiffs" + ".png", phase_sync.vel_phase_diffs, sys_msr_data.time, phase_sync.num_phase_diffs, r"$\phi_n - \phi_{n + 3}$")

        ##-------- Magnetic
        phase_only_space_time(cmdargs.out_dir_MHD + "B_SpaceTime_Phases" + ".png", np.mod(np.angle(run_data.b) + 2.0*np.pi, 2.0*np.pi), sys_msr_data.time, sys_vars.N, r"$\psi_n$")
        phase_only_space_time(cmdargs.out_dir_MHD + "B_SpaceTime_PhaseDiffs" + ".png", phase_sync.mag_phase_diffs, sys_msr_data.time, phase_sync.num_phase_diffs, r"$\psi_n - \psi_{n + 3}$")
        for tri_type, type_str in enumerate(["Tpss", "Tsps", "Tssp"]):
            phase_only_space_time(cmdargs.out_dir_MHD + "B_SpaceTime_Triads_" + type_str + ".png", phase_sync.mag_triads[:, tri_type, :], sys_msr_data.time, phase_sync.num_triads, type_str)

        ## Elsassar Variables
        phase_only_space_time(cmdargs.out_dir_MHD + "Zplus_SpaceTime_Phases" + ".png", np.mod(np.angle(z_plus) + 2.0*np.pi, 2.0*np.pi), sys_msr_data.time, sys_vars.N, r"$\arg Z_n^+$")
        phase_only_space_time(cmdargs.out_dir_MHD + "Zminus_SpaceTime_Phases" + ".png", np.mod(np.angle(z_minus) + 2.0*np.pi, 2.0*np.pi), sys_msr_data.time, sys_vars.N, r"$\arg Z_n^-$")
        phase_only_space_time(cmdargs.out_dir_MHD + "Tplus_SpaceTime_Phases" + ".png", np.mod(np.angle(t_plus) + 2.0*np.pi, 2.0*np.pi), sys_msr_data.time, sys_vars.N, r"$\arg T_n^+$")
        phase_only_space_time(cmdargs.out_dir_MHD + "Tminus_SpaceTime_Phases" + ".png", np.mod(np.angle(t_minus) + 2.0*np.pi, 2.0*np.pi), sys_msr_data.time, sys_vars.N, r"$\arg T_n^-$")




        # --------------------------------------------------
        # # --------  Plot Structure Funcs + Anomalous Scaling
        # ---------------------------------------------------
        print("Plotting --- Structure Funcs")

        inert_lim_low  = 3
        inert_lim_high = 12
        range_lims     = [inert_lim_low, inert_lim_high]
        log_func = 'loge'
        k           = sys_msr_data.k

        ########################### Amp
        ##----------- Velocity Amp 
        str_funcs   = stats_data.vel_str_func[:, :] / stats_data.num_stats_steps
        inert_range = range_lims
        print("\nVel")
        vel_zeta_p, ns_zeta_p, vel_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_U" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

        ##----------- Magnetic Amp 
        str_funcs   = stats_data.mag_str_func[:, :] / stats_data.num_stats_steps
        inert_range = range_lims
        print("\nMag")
        mag_zeta_p, ns_zeta_p, mag_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_B" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

        ############################ Tripple Prod
        ##----------- Velocity Tripple Prod
        str_funcs   = stats_data.vel_trip_prod_str_func_abs[:, :] / stats_data.num_stats_steps
        inert_range = range_lims
        print("\nVel Trip")
        vel_trip_zeta_p, ns_zeta_p, vel_trip_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_U_TripProd" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

        ##----------- Magnetic Tripple Prod
        str_funcs   = stats_data.mag_trip_prod_str_func_abs[:, :] / stats_data.num_stats_steps
        inert_range = range_lims
        print("\nMag Trip")
        mag_trip_zeta_p, ns_zeta_p, mag_trip_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_B_TripProd" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))


        ############################ Flux Terms
        ##----------- Velocity Energy Flux
        str_funcs   = stats_data.vel_flux_str_func_abs[:, :, 0] / stats_data.num_stats_steps
        inert_range = range_lims
        print("\nVel Flux")
        vel_flux_zeta_p, ns_zeta_p, vel_flux_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_U_Flux" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

        ##----------- Magnetic Energy Flux
        str_funcs   = stats_data.mag_flux_str_func_abs[:, :, 0] / stats_data.num_stats_steps
        inert_range = range_lims
        print("\nMag Flux")
        mag_flux_zeta_p, ns_zeta_p, mag_flux_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_B_Flux" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))


        ############################ Elsassar Variables
        ## Get the Elsassar Str Funcs
        zp_sf, zm_sf, tp_sf, tm_sf = get_elsassar_sf(z_plus, z_minus, np.imag(t_plus), np.imag(t_minus))

        ##----------- Zplus
        str_funcs   = zp_sf
        inert_range = range_lims
        print("\nZplus")
        zp_zeta_p, ns_zeta_p, zp_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_Zplus" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

        ##----------- Zminus
        str_funcs   = zm_sf
        inert_range = range_lims
        print("\nZminus")
        zm_zeta_p, ns_zeta_p, zm_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_Zminus" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

        ##----------- Tplus
        str_funcs   = tp_sf
        inert_range = range_lims
        print("\nTplus")
        tp_zeta_p, ns_zeta_p, tp_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_Tplus" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

        ##----------- Tminus
        str_funcs   = tm_sf
        inert_range = range_lims
        print("\nTminus")
        tm_zeta_p, ns_zeta_p, tm_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_MHD + "SF_Tminus" + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))



        ###### ---------------- Plot results of fit to screen and text file
        with open(cmdargs.out_dir_MHD + "ComputedStrFuncSlopes.txt", 'w') as f:
            print("\n\nPower\tp/3\t\tNS\t\tVel\t\t\tMag\t\t     Vel Trip\t\t   Mag Trip\t\t   Vel Flux\t\t   Mag Flux\t\t   Zplus\t\t\tZminus\t\t\tTplus\t\t\tTminus")
            f.write("\n\nPower\tp/3\t\tNS\t\tVel\t\t\tMag\t\t     Vel Trip\t\t   Mag Trip\t\t   Vel Flux\t\t   Mag Flux\t\t   Zplus\t\t\tZminus\t\t\tTplus\t\t\tTminus\n")
            for i in range(len(vel_zeta_p)):
                if i >= 1 and i <= len(ns_zeta_p):
                    print(" {}\t {:1.3f} \t{:1.3f} \t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}".format(i + 1, -(i + 1) / 3, -ns_zeta_p[i - 1], 
                                                                                                        vel_zeta_p[i], vel_zeta_p_resid[i],
                                                                                                        mag_zeta_p[i], mag_zeta_p_resid[i],
                                                                                                        vel_trip_zeta_p[i], vel_trip_zeta_p_resid[i],
                                                                                                        mag_trip_zeta_p[i], mag_trip_zeta_p_resid[i],
                                                                                                        vel_flux_zeta_p[i], vel_flux_zeta_p_resid[i],
                                                                                                        mag_flux_zeta_p[i], mag_flux_zeta_p_resid[i],
                                                                                                        zp_zeta_p[i], zp_zeta_p_resid[i],
                                                                                                        zm_zeta_p[i], zm_zeta_p_resid[i],
                                                                                                        tp_zeta_p[i], tp_zeta_p_resid[i],
                                                                                                        tm_zeta_p[i], tm_zeta_p_resid[i]
                                                                                                        ))
                    f.write(" {}\t {:1.3f} \t{:1.3f} \t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, -ns_zeta_p[i - 1], 
                                                                                                        vel_zeta_p[i], vel_zeta_p_resid[i],
                                                                                                        mag_zeta_p[i], mag_zeta_p_resid[i],
                                                                                                        vel_trip_zeta_p[i], vel_trip_zeta_p_resid[i],
                                                                                                        mag_trip_zeta_p[i], mag_trip_zeta_p_resid[i],
                                                                                                        vel_flux_zeta_p[i], vel_flux_zeta_p_resid[i],
                                                                                                        mag_flux_zeta_p[i], mag_flux_zeta_p_resid[i],
                                                                                                        zp_zeta_p[i], zp_zeta_p_resid[i],
                                                                                                        zm_zeta_p[i], zm_zeta_p_resid[i],
                                                                                                        tp_zeta_p[i], tp_zeta_p_resid[i],
                                                                                                        tm_zeta_p[i], tm_zeta_p_resid[i]
                                                                                                        ))
                else:
                    print(" {}\t {:1.3f} \t{:1.3f} \t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}".format(i + 1, -(i + 1) / 3, 0.000, 
                                                                                                        vel_zeta_p[i], vel_zeta_p_resid[i],
                                                                                                        mag_zeta_p[i], mag_zeta_p_resid[i],
                                                                                                        vel_trip_zeta_p[i], vel_trip_zeta_p_resid[i],
                                                                                                        mag_trip_zeta_p[i], mag_trip_zeta_p_resid[i],
                                                                                                        vel_flux_zeta_p[i], vel_flux_zeta_p_resid[i],
                                                                                                        mag_flux_zeta_p[i], mag_flux_zeta_p_resid[i],
                                                                                                        zp_zeta_p[i], zp_zeta_p_resid[i],
                                                                                                        zm_zeta_p[i], zm_zeta_p_resid[i],
                                                                                                        tp_zeta_p[i], tp_zeta_p_resid[i],
                                                                                                        tm_zeta_p[i], tm_zeta_p_resid[i]
                                                                                                        ))
                    f.write(" {}\t {:1.3f} \t{:1.3f} \t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\t {:1.4f} +/-{:0.3f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, 0.000, 
                                                                                                        vel_zeta_p[i], vel_zeta_p_resid[i],
                                                                                                        mag_zeta_p[i], mag_zeta_p_resid[i],
                                                                                                        vel_trip_zeta_p[i], vel_trip_zeta_p_resid[i],
                                                                                                        mag_trip_zeta_p[i], mag_trip_zeta_p_resid[i],
                                                                                                        vel_flux_zeta_p[i], vel_flux_zeta_p_resid[i],
                                                                                                        mag_flux_zeta_p[i], mag_flux_zeta_p_resid[i],
                                                                                                        zp_zeta_p[i], zp_zeta_p_resid[i],
                                                                                                        zm_zeta_p[i], zm_zeta_p_resid[i],
                                                                                                        tp_zeta_p[i], tp_zeta_p_resid[i],
                                                                                                        tm_zeta_p[i], tm_zeta_p_resid[i]
                                                                                                        ))