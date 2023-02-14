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
from functions import tc, sim_data, import_data, compute_pdf, import_stats_data, import_sys_msr_data, import_phase_sync_data, compute_pdf_from_hist, compute_pdf, parse_cml, slope_fit
from plot_functions import plot_anomalous_exponent, plot_str_funcs_with_slope, phase_only_space_time

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
	cmdargs.out_dir_HD = cmdargs.out_dir + "HD_PLOTS/"
	if os.path.isdir(cmdargs.out_dir_HD) != True:
		print("Making folder:" + tc.C + " HD_PLOTS/" + tc.Rst)
		os.mkdir(cmdargs.out_dir_HD)
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

		fig_size = (10, 8)
		fig_file_type = ".png"
		num_pow = 6
		# num_pow = stats_data.vel_str_func.shape[-1]

		###################################
		##  PO Model Dynamics Diagnostic
		###################################
		if sys_vars.model_type == "PHASEONLY":
			phases = run_data.phi_n
		else:
			phases = np.mod(np.angle(run_data.u) + 2.0*np.pi, 2.0*np.pi)
		# Individual phases
		phase_only_space_time(cmdargs.out_dir_HD + "Phases_SpaceTime_Dynamics_Diagnostic_Phases" + fig_file_type, phases, sys_msr_data.time, sys_vars.N, r"$\phi_n$")

		# ## Read in sys_msr data
		phase_sync = import_phase_sync_data(data_file_path, sys_vars, method)
		
		## Triads
		phase_only_space_time(cmdargs.out_dir_HD + "Phases_SpaceTime_Dynamics_Diagnostic_Vel_Triads" + fig_file_type, phase_sync.vel_triads, sys_msr_data.time, phase_sync.num_triads, r"$\varphi_{n ,n + 1}^{n + 2}$")

		## Phase Differences
		phase_only_space_time(cmdargs.out_dir_HD + "Phases_SpaceTime_Dynamics_Diagnostic_Vel_PhaseDiffs" + fig_file_type, phase_sync.vel_phase_diffs, sys_msr_data.time, phase_sync.num_phase_diffs, r"$\phi_n - \phi_{n + 3}$")

		## Plot the velocity triad pdfs
		fig = plt.figure(figsize = (16, 16))
		gs  = GridSpec(5, 5, wspace = 0.25, hspace = 0.25)
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

		plt.savefig(cmdargs.out_dir_HD + "Phases_Vel_Triad_PDF.png", bbox_inches='tight')
		plt.close()

		## Plot the velocity triad phase differences
		fig = plt.figure(figsize = (16, 16))
		gs  = GridSpec(5, 5, wspace = 0.25, hspace = 0.25)
		for i in range(5):
			for j in range(5):
				if i * 5 + j < phase_sync.num_phase_diffs:
					ax1 = fig.add_subplot(gs[i, j])
					pdf, centres = compute_pdf_from_hist(phase_sync.vel_phase_diff_hist_counts[i * 5 + j, :], phase_sync.vel_phase_diff_hist_ranges[:], remove_zeros = False)
					ax1.plot(centres, pdf, label = "$Tppp({})$".format(i * 5 + j + 1))
					ax1.set_xlabel("$\phi_n - \phi_{n + 3}$")
					ax1.set_ylabel("PDF")
					ax1.set_yscale("log")
					ax1.legend()
					ax1.set_xlim(0, 2.0*np.pi)
					ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
					ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])

		plt.savefig(cmdargs.out_dir_HD + "Phases_Vel_PhaseDifferences_PDF.png", bbox_inches='tight')
		plt.close()

		###################################
		##  PDFs
		###################################
		fig = plt.figure(figsize = fig_size)
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		for j, i in enumerate([6, 11, 16, 20]):
			if hasattr(stats_data, "vel_hist_counts"):
				pdf, centres = compute_pdf_from_hist(stats_data.vel_hist_counts[i, :], stats_data.vel_hist_ranges[i, :], normed = True)
				ax1.set_title("C Data")
			else:
				if sys_vars.model_type == "PHASEONLY":
					pdf, centres = compute_pdf(np.real(run_data.a_n[:, i] * np.exp(1j * run_data.phi_n[:, i])), nbins = 500, normed = True)					
				else:
					pdf, centres = compute_pdf(np.real(run_data.u[:, i]), nbins = 500, normed = True)
			if i == -1:
				ax1.plot(centres, pdf, label = "$n = {}$".format(sys_vars.N))    
			else:
				ax1.plot(centres, pdf, label = "$n = {}$".format(i - 1))    
		ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
		ax1.set_ylabel(r"PDF")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.set_yscale('log')
		# ax1.set_title(sys_vars.model_type)
		ax1.legend()

		plt.savefig(cmdargs.out_dir_HD + "VelReal_PDF_InOne" + fig_file_type)
		plt.close()

		###################################
		##  Structure Functions
		###################################
		inertial_range = np.arange(7, 13 + 1)
		p = np.arange(2, num_pow + 1)

		ax_scale = "log2"

		## --------  Structure function with fit
		print("Velocity Structure Func")
		zeta_p, ns_zeta_p, zeta_p_res = plot_str_funcs_with_slope(cmdargs.out_dir_HD + "VelStrFunc_Fit" + fig_file_type, sys_msr_data.k, stats_data.vel_str_func[:, :num_pow], inertial_range, ax_scale, fig_size = fig_size)
		with open(cmdargs.out_dir_HD + "ComputedStrFuncSlopes.txt", 'w') as f:
			f.write("Velocity Structure Func\n")
			f.write("Power\tp/3\t\tDNS Slope\tNS\n")
			for i in range(num_pow):
				if i >= 1:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f} \t{:1.3f}\n".format(i + 1, -(i + 1) / 3, zeta_p[i], zeta_p_res[i], -ns_zeta_p[i - 1]))
				else:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, zeta_p[i], zeta_p_res[i]))

		print(len(p), len(zeta_p), len(ns_zeta_p))
		## --------  Plot Anomalous Exponent
		plot_anomalous_exponent(cmdargs.out_dir_HD + "Vel_Anonalous_Exponent_Zeta_p" + fig_file_type, p, zeta_p[1:], label_str = r"Velocity; Shell Modell", fig_size = fig_size)


		## --------  Velocity Tripple Product
		if hasattr(stats_data, "vel_trip_prod_str_func_abs"):
			print("Velocity Tripple Product Structure Func")
			trip_prod_zeta_p, ns_zeta_p, trip_prod_zeta_p_res = plot_str_funcs_with_slope(cmdargs.out_dir_HD + "TrippleProd_VelStrFunc_Fit" + fig_file_type, sys_msr_data.k, stats_data.vel_trip_prod_str_func_abs[:, :num_pow], inertial_range, ax_scale, fig_size = fig_size)
			with open(cmdargs.out_dir_HD + "ComputedStrFuncSlopes.txt", 'a') as f:
				f.write("Velocity Tripple Product Structure Func\n")
				f.write("Power\tp/3\t\tDNS Slope\tNS\n")
				for i in range(num_pow):
					if i >= 1:
						f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f} \t{:1.3f}\n".format(i + 1, -(i + 1) / 3, trip_prod_zeta_p[i], trip_prod_zeta_p_res[i], -ns_zeta_p[i - 1]))
					else:
						f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, trip_prod_zeta_p[i], trip_prod_zeta_p_res[i]))

		## --------  Plot Anomalous Exponent
		plot_anomalous_exponent(cmdargs.out_dir_HD + "TrippleProd_Vel_Anonalous_Exponent_Zeta_p" + fig_file_type, p, trip_prod_zeta_p[1:], label_str = r"Tripple Prod Velocity; Shell Modell", fig_size = fig_size)

		## --------  Structure function with fit
		print("Velocity Energy Structure Func")
		enrg_flux_zeta_p, ns_zeta_p, enrg_flux_zeta_p_res = plot_str_funcs_with_slope(cmdargs.out_dir_HD + "VelEnergyFluxAbsStrFunc_Fit" + fig_file_type, sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, :num_pow, 0], inertial_range, ax_scale, fig_size = fig_size)
		with open(cmdargs.out_dir_HD + "ComputedStrFuncSlopes.txt", 'a') as f:
			f.write("Velocity Energy Structure Func\n")
			f.write("Power\tp/3\t\tDNS Slope\tNS\n")
			for i in range(num_pow):
				if i >= 1:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f} \t{:1.3f}\n".format(i + 1, -(i + 1) / 3, enrg_flux_zeta_p[i], enrg_flux_zeta_p_res[i], -ns_zeta_p[i - 1]))
				else:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, enrg_flux_zeta_p[i], enrg_flux_zeta_p_res[i]))

		## --------  Plot Anomalous Exponent
		plot_anomalous_exponent(cmdargs.out_dir_HD + "VelEnergyFluxAbs_Anonalous_Exponent_Zeta_p" + fig_file_type, p, enrg_flux_zeta_p[1:], label_str = r"Velocity Energy Flux; Shell Modell", fig_size = fig_size)

		## --------  Structure function with fit
		print("Velocity Helicity Structure Func")
		hel_flux_zeta_p, ns_zeta_p, hel_flux_zeta_p_res = plot_str_funcs_with_slope(cmdargs.out_dir_HD + "VelHelicityFluxAbsStrFunc_Fit" + fig_file_type, sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, :num_pow, 1], inertial_range, ax_scale, fig_size = fig_size)
		with open(cmdargs.out_dir_HD + "ComputedStrFuncSlopes.txt", 'a') as f:
			f.write("Velocity Helicity Structure Func\n")
			f.write("Power\tp/3\t\tDNS Slope\tNS\n")
			for i in range(num_pow):
				if i >= 1:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f} \t{:1.3f}\n".format(i + 1, -(i + 1) / 3, hel_flux_zeta_p[i], hel_flux_zeta_p_res[i], -ns_zeta_p[i - 1]))
				else:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, hel_flux_zeta_p[i], hel_flux_zeta_p_res[i]))


		## --------  Plot Anomalous Exponent
		plot_anomalous_exponent(cmdargs.out_dir_HD + "VelHelicityFluxAbs_Anonalous_Exponent_Zeta_p" + fig_file_type, p, hel_flux_zeta_p[1:], label_str = r"Velocity HelicityFlux; Shell Modell", fig_size = fig_size)


		###################################
		##  Anomalous Exponent
		###################################
		fig   = plt.figure(figsize = fig_size)
		gs    = GridSpec(1, 1)
		ax1   = fig.add_subplot(gs[0, 0])
		mark_style = ['o','s','^','x','D','p']
		ax1.plot(p, zeta_p[1:], marker = mark_style[0], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Shell Model")
		if hasattr(stats_data, "vel_trip_prod_str_func_abs"):
			ax1.plot(p, trip_prod_zeta_p[1:], marker = mark_style[1], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Tripple Prod Shell Model")
		ax1.plot(p, enrg_flux_zeta_p[1:], marker = mark_style[2], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Energy Flux Shell Model")
		ax1.plot(p, hel_flux_zeta_p[1:], marker = mark_style[3], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Helicity lux Shell Model")
		ax1.plot(p, ns_zeta_p[:len(zeta_p[1:])], marker = mark_style[4], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Navier Stokes")
		ax1.plot(p, p/9.0 + 2.0 * (1.0 - np.power(2.0/3.0, p/3.0)), marker = mark_style[5], label = "She-Leveque Model; $p/9 + 2(1 - (2/3)^{p/3})$")
		ax1.plot(p, p / 3, 'k--', label = "K41")
		ax1.set_xlim(2.0, num_pow)
		ax1.set_xlabel(r"$p$")
		ax1.set_ylabel(r"$\zeta_p$")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.legend()
		plt.savefig(cmdargs.out_dir_HD + "VelCombined_Anonalous_Exponent_Zeta_p" + fig_file_type, bbox_inches='tight')
		plt.close()

		###################################
		##  Time Averaged Flux Spectra
		###################################
		fig = plt.figure(figsize = fig_size)
		gs  = GridSpec(1, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		tmp1 = np.copy(sys_msr_data.enrg_flux_t_avg)
		tmp2 = np.copy(sys_msr_data.enrg_flux_t_avg)
		tmp1[tmp1 < 0] = 0.0
		tmp2[tmp2 > 0] = 0.0
		ax1.plot(sys_msr_data.k, tmp1 * sys_msr_data.k**1, '.')
		ax1.plot(sys_msr_data.k, np.absolute(tmp2) * sys_msr_data.k**1, '.')
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$k_n \mathcal{E}_n$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax2 = fig.add_subplot(gs[0, 1])
		p, = ax2.plot(sys_msr_data.k, np.absolute(sys_msr_data.enrg_flux_t_avg))
		slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(sys_msr_data.enrg_flux_t_avg)), inertial_range[0], inertial_range[-1])
		ax2.plot(sys_msr_data.k, sys_msr_data.k ** (slope) * np.exp(c), '--', label = "$slope = {}$".format(np.around(slope, 3)), color = p.get_color())
		ax2.plot(sys_msr_data.k, sys_msr_data.k ** (-1), 'k--', label = "$k^{-1}$")
		ax2.set_xlabel(r"$k_n$")
		ax2.set_ylabel(r"$\mathcal{E}_n$")
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_HD + "Time_Averaged_Energy_Flux_Spectra" + fig_file_type, bbox_inches='tight')
		plt.close()


		if hasattr(sys_msr_data, "kin_hel_flux_t_avg"):
			tmp1 = np.copy(sys_msr_data.kin_hel_flux_t_avg)
			tmp2 = np.copy(sys_msr_data.kin_hel_flux_t_avg)
			tmp1[tmp1 < 0] = 0.0
			tmp2[tmp2 > 0] = 0.0
			fig = plt.figure(figsize = fig_size)
			gs  = GridSpec(1, 2)
			ax1 = fig.add_subplot(gs[0, 0])
			ax1.plot(sys_msr_data.k, tmp1 * sys_msr_data.k**2, '.')
			ax1.plot(sys_msr_data.k, np.absolute(tmp2) * sys_msr_data.k**2, '.')
			ax1.set_xlabel(r"$k_n$")
			ax1.set_ylabel(r"$k_n \mathcal{E}_n$")
			ax1.set_xscale("log")
			ax1.set_yscale("log")
			ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
			ax2 = fig.add_subplot(gs[0, 1])
			p, = ax2.plot(sys_msr_data.k, np.absolute(sys_msr_data.kin_hel_flux_t_avg))
			slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(sys_msr_data.kin_hel_flux_t_avg)), inertial_range[0], inertial_range[-1])
			ax2.plot(sys_msr_data.k, sys_msr_data.k ** (slope) * np.exp(c), '--', label = "$slope = {}$".format(np.around(slope, 3)), color = p.get_color())
			ax2.plot(sys_msr_data.k, sys_msr_data.k ** (-2), 'k--', label = "$k^{-2}$")
			ax2.set_xlabel(r"$k_n$")
			ax2.set_ylabel(r"$\mathcal{E}_n$")
			ax2.set_xscale("log")
			ax2.set_yscale("log")
			ax2.legend()
			ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
			plt.savefig(cmdargs.out_dir_HD + "Time_Averaged_Kinetic_Helicity_Flux_Spectra" + fig_file_type, bbox_inches='tight')
			plt.close()


		###################################
		##  Time Averaged Spectra
		###################################
		fig = plt.figure(figsize = fig_size)
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(sys_msr_data.k, sys_msr_data.enrg_spec_t_avg, label = "Time Averaged Energy Spectrum")
		slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(sys_msr_data.enrg_spec_t_avg)), inertial_range[0], inertial_range[-1])
		ax1.plot(sys_msr_data.k, sys_msr_data.k ** (slope) * np.exp(c), '--', label = "$slope = {}$".format(np.around(slope, 3)), color = p.get_color())
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$\mathcal{E}_n$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.legend()
		plt.savefig(cmdargs.out_dir_HD + "Time_Averaged_Energy_Spectra" + fig_file_type, bbox_inches='tight')
		plt.close()

		###################################
		##  Time Averaged Amplitudes
		###################################
		fig = plt.figure(figsize = fig_size)
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(sys_msr_data.k, sys_msr_data.a_n_t_avg, label = "Time Averaged Amplitudes")
		slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(sys_msr_data.a_n_t_avg)), inertial_range[0], inertial_range[-1])
		ax1.plot(sys_msr_data.k, sys_msr_data.k ** (slope) * np.exp(c), '--', label = "$slope = {}$".format(np.around(slope, 3)), color = p.get_color())
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$\langle a_n \rangle$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.legend()
		plt.savefig(cmdargs.out_dir_HD + "Time_Averaged_Amplitudes" + fig_file_type, bbox_inches='tight')
		plt.close()

		###################################
		##  Time Series of Energy Dissipation
		###################################
		fig = plt.figure(figsize = fig_size)
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(sys_msr_data.time, sys_msr_data.tot_diss_u, label = "Total Energy Dissipation")
		ax1.set_xlabel(r"Time")
		# ax1.set_xlim(sys_msr_data.time[0], sys_msr_data.time[-1])
		ax1.set_ylim(bottom = 0)
		ax1.set_ylabel(r"$\epsilon(t)$")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.legend()
		plt.savefig(cmdargs.out_dir_HD + "Energy_Dissipation" + fig_file_type, bbox_inches='tight')
		plt.close()


		###################################
		##  Total Fluxes Time Series
		###################################
		fig = plt.figure(figsize = fig_size)
		gs  = GridSpec(1, 2, wspace = 0.35)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(sys_msr_data.time, sys_msr_data.tot_vel_enrg_flux, label = "Total Energy Flux")
		ax1.set_xlabel(r"Time")
		# ax1.set_xlim(sys_msr_data.time[0], sys_msr_data.time[-1])
		ax1.set_ylabel(r"$\Pi^{\mathcal{E}}(t)$")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.legend()
		ax2 = fig.add_subplot(gs[0, 1])
		ax2.plot(sys_msr_data.time, sys_msr_data.tot_kin_hel_flux, label = "Total Kinetic Helicity Flux")
		ax2.set_xlabel(r"Time")
		# ax2.set_xlim(sys_msr_data.time[0], sys_msr_data.time[-1])
		ax2.set_ylabel(r"$\Pi^{\mathcal{H}}(t)$")
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax2.legend()
		plt.savefig(cmdargs.out_dir_HD + "Total_Fluxes" + fig_file_type, bbox_inches='tight')
		plt.close()
