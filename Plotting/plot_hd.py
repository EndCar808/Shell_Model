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
from functions import tc, sim_data, import_data, compute_pdf, import_stats_data, import_sys_msr_data, import_phase_sync_data, compute_pdf_from_hist, compute_pdf, parse_cml, slope_fit, get_vel_field_flux
from plot_functions import plot_anomalous_exponent, plot_str_funcs_with_slope, phase_only_space_time, plot_spectrum

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
		fig_format = "png"
		# num_pow = 6
		num_pow = stats_data.vel_str_func.shape[-1]

		## Plot The structure function
        inert_lim_low  = 3
        inert_lim_high = 12
        range_lims     = [inert_lim_low, inert_lim_high]
        log_func = 'loge'
        k           = sys_msr_data.k
		inert_range = range_lims

		print("Getting Flux Values")
		if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
			u_n = run_data.a_n * np.exp(1j * run_data.phi_n)
		else:
			u_n = run_data.u
		trip_prod, dub_prod, hel_flux, enrg_flux = get_vel_field_flux(u_n, N, delta, l)

		####---------------------------------------- Plot Flux Field Values
		input_data = [trip_prod, dub_prod, hel_flux, enrg_flux]
		figure_names = [r"TripProd", r"DubProd", r"HelFlux", r"EnergyFlux"]
		data_labels  = [r"u_{n + 2}u_{n + 1}u_{n}", r"u_{n}u_{n + 3}^*", r"\Pi_n^{\mathcal{H}^u}", r"\Pi_n^{\mathcal{K}^u}"]
		for in_data, fig_name, data_labs in zip(input_data, figure_names, data_labels):
			#-------------------- 1D - PDF - Angle
			num_bins = 100
			norm_hist = False
			fig = plt.figure(figsize=(24, 24))
			gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
			pdf_angles = []
			for i in range(5):
				for j in range(5):
					indx = i * 5 + j
					if indx < in_data.shape[-1]:
						ax1 = fig.add_subplot(gs[i, j])
						pdf, centres = compute_pdf(np.mod(np.angle(in_data[:, indx]) + 2.0 * np.pi, 2.0 * np.pi), nbins=num_bins, normed=norm_hist)
						pdf_angles.append(pdf)
						p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
						ax1.set_xlabel(r"$ \arg \left\{" +  data_labs + r" \right\}$")
						ax1.set_xlim(0, 2.0 * np.pi)
						ax1.set_ylabel(r"PDF")
						ax1.set_yscale('log')
						ax1.set_title("n = {}".format(indx + 1))
			fig.savefig(cmdargs.out_dir_HD + fig_name + "_1D_PDF_Angles" + fig_format, bbox_inches='tight')
			plt.close()

			#-------------------- 1D - PDF - Amp
			fig = plt.figure(figsize=(24, 24))
			gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
			for i in range(5):
				for j in range(5):
					indx = i * 5 + j
					if indx < in_data.shape[-1]:
						ax1 = fig.add_subplot(gs[i, j])
						pdf, centres = compute_pdf(np.absolute(in_data[:, indx]), nbins=num_bins, normed=norm_hist)
						p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
						ax1.set_xlabel(r"$ \left|" +  data_labs + r" \right|$")
						ax1.set_ylabel(r"PDF")
						ax1.set_yscale('log')
						ax1.set_title("n = {}".format(indx + 1))
			fig.savefig(cmdargs.out_dir_HD + fig_name + "_1D_PDF_Amp" + fig_format, bbox_inches='tight')
			plt.close()

			## 1D - PDF - Imag
			num_bins = 100
			norm_hist = False
			fig = plt.figure(figsize=(24, 24))
			gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
			for i in range(5):
				for j in range(5):
					indx = i * 5 + j
					if indx < in_data.shape[-1]:
						ax1 = fig.add_subplot(gs[i, j])
						pdf, centres = compute_pdf(np.imag(in_data[:, indx]), nbins=num_bins, normed=norm_hist)
						p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
						ax1.set_xlabel(r"$ \Im \left\{" +  data_labs + r" \right\}$")
						ax1.set_ylabel(r"PDF")
						ax1.set_yscale('log')
						ax1.set_title("n = {}".format(indx + 1))
			fig.savefig(cmdargs.out_dir_HD + fig_name + "_1D_PDF_Imag" + fig_format, bbox_inches='tight')
			plt.close()

			## 1D - PDF - Real
			fig = plt.figure(figsize=(24, 24))
			gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
			for i in range(5):
				for j in range(5):
					indx = i * 5 + j
					if indx < in_data.shape[-1]:
						ax1 = fig.add_subplot(gs[i, j])
						pdf, centres = compute_pdf(np.real(in_data[:, indx]), nbins=num_bins, normed=norm_hist)
						p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
						ax1.set_xlabel(r"$ \Re \left\{" +  data_labs + r" \right\}$")
						ax1.set_ylabel(r"PDF")
						ax1.set_yscale('log')
						ax1.set_title("n = {}".format(indx + 1))
			fig.savefig(cmdargs.out_dir_HD + fig_name + "_1D_PDF_Real" + fig_format, bbox_inches='tight')
			plt.close()

			#---------------------- Test SF Independence
			fig = plt.figure(figsize=(15, 6))
			gs  = GridSpec(1, 3, hspace=0.4, wspace=0.3)
			ax1 = fig.add_subplot(gs[0, 0])
			for p in range(1, 6 + 1):
				ax1.plot(k, np.mean(np.power(np.absolute(in_data), p/3.), axis=0), label=r"$Abs$; $p = {}$".format(p))
			ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
			ax1.set_yscale('log')
			ax1.set_xscale('log')
			ax1.legend()
			ax1 = fig.add_subplot(gs[0, 1])
			for p in range(1, 6 + 1):
				ax1.plot(k, np.mean(np.power(np.sign(np.sin(np.angle(in_data))), p) * np.power(np.absolute(np.sin(np.angle(in_data))), p/3.), axis=0), label=r"$\sin \arg$; $p = {}$".format(p))
			ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
			ax1.set_xscale('log')
			ax1.legend()
			ax1 = fig.add_subplot(gs[0, 2])
			slope = []
			for p in range(1, 6 + 1):
				p_info, = ax1.plot(k[:-6], np.mean(np.power(np.absolute(np.sin(np.angle(in_data))), p/3.), axis=0)[:-6], label=r"$\arg $; $p = {}$".format(p), marker = 'o')
				# ax1.plot(k, np.mean(np.power(np.absolute(np.cos(np.angle(in_data))), p/3.), axis=0), label=r"$\arg $; $p = {}$".format(p))
				# ax1.plot(k, np.mean( np.power(np.sign(np.angle(in_data)), p) * np.power(np.absolute(np.angle(in_data)), p/3.), axis=0), label=r"$\arg $; $p = {}$".format(p))
				x = k
				y = np.mean(np.power(np.absolute(np.sin(np.angle(in_data))), p/3.), axis=0)
				s, c, res = slope_fit(np.log10(x), np.log10(y), 2, 11)
				x = k
				y = np.mean(np.power(np.absolute(in_data), p/3.), axis=0)
				s_amp, c, res_amp = slope_fit(np.log10(x), np.log10(y), 2, 11)
				print(fig_name, p, s, res, s_amp, res_amp)
				# ax1.plot(x[2:11], x[2:11] ** slope + c, ':', color = p_info.get_color())

			ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
			ax1.set_xscale('log')
			ax1.set_yscale('log')
			ax1.legend()
			plt.savefig(cmdargs.out_dir_HD + fig_name + "_Test_Independence_SF" + fig_format, bbox_inches='tight')
			plt.close()


			#----------------------
			# Time Averaged Spectra
			#----------------------
			fig = plt.figure(figsize=fig_size)
			gs = GridSpec(2, 2, hspace=0.4, wspace=0.3)
			ax1 = fig.add_subplot(gs[0, 0])
			ax2 = fig.add_subplot(gs[0, 1])
			ax3 = fig.add_subplot(gs[1, 0])
			ax4 = fig.add_subplot(gs[1, 1])

			# Plot the time averaged tripple product
			data = np.absolute(np.mean(np.imag(in_data), axis=0))
			xlab = r"$\log k_n$"
			ylab = r"$\log \left|\Im \left\{" + data_labs + r"\right\}\right|$"
			plot_spectrum(fig, ax1, data, k, xlab, ylab)
			ax1.set_title("Time Averaged Imaginary part")

			# Plot the time averaged tripple product time k_n
			data = np.absolute(np.mean(np.imag(in_data * k), axis=0))
			ylab = r"$\log \left|k_n \Im \left\{" + data_labs + r"\right\}\right|$"
			plot_spectrum(fig, ax2, data, k, xlab, ylab)
			ax2.set_title("Time Averaged Imaginary part time $k_n$")

			# Plot the time averaged tripple product
			data = np.absolute(np.mean(np.real(in_data), axis=0))
			xlab = r"$\log k_n$"
			ylab = r"$\log \left|\Re \left\{" + data_labs + r"\right\}\right|$"
			plot_spectrum(fig, ax3, data, k, xlab, ylab)
			ax3.set_title("Time Averaged real part")

			# Plot the time averaged tripple product time k_n
			data = np.absolute(np.mean(np.real(in_data *  k), axis=0))
			ylab = r"$\log \left|k_n \Re \left\{" + data_labs + r"\right\}\right|$"
			plot_spectrum(fig, ax4, data, k, xlab, ylab)
			ax4.set_title("Time Averaged Real part time $k_n$")

			# Save figure
			plt.suptitle(fig_name)
			fig.savefig(cmdargs.out_dir_HD + fig_name + "_Time_Averaged_Spectrum" + fig_format, bbox_inches='tight')
			plt.close()

			fig = plt.figure(figsize=fig_size)
			gs = GridSpec(1, 1, hspace=0.4, wspace=0.5)
			ax1 = fig.add_subplot(gs[0, 0])
			ax1.scatter(np.absolute(in_data).flatten(), np.sin(np.angle(in_data).flatten()))
			ax2.set_xlabel(r"$\left|" + data_labs + r"\right|$")
			ax2.set_ylabel(r"$\sin \arg \left\{ " + data_labs + r" \right\}$")
			plt.suptitle(fig_name)
			fig.savefig(cmdargs.out_dir_HD + fig_name + "_Scatter" + fig_format, bbox_inches='tight')
			plt.close()

		###----------------------------------------- Plot SF
		input_data = [stats_data.vel_str_func[:, :], stats_data.vel_trip_prod_str_func_abs[:, :], stats_data.vel_flux_str_func_abs[:, :, 1], stats_data.vel_flux_str_func_abs[:, :, 0]]
		figure_names = [r"SF_U", r"SF_Trip", r"SF_HFlux", r"SF_EFlux"]
		for in_data, fig_name in zip(input_data, figure_names):
			str_funcs   = in_data / stats_data.num_stats_steps
			vel_zeta_p, ns_zeta_p, vel_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_HD + fig_name + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))



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
		ax2.plot(sys_msr_data.time, sys_msr_data.tot_hel_u)
		ax2.set_xlabel(r"$t$")
		ax2.set_title(r"Total Kinetic Helicity")
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax4 = fig.add_subplot(gs[1, 0])
		ax4.plot(sys_msr_data.time, sys_msr_data.tot_diss_u, label = "$\epsilon_u$")
		ax4.set_xlabel(r"$t$")
		ax4.set_yscale('log')
		ax4.set_title(r"Total Dissipation")
		ax4.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		fig.savefig(cmdargs.out_dir_HD  + "System_Measures" + "." + fig_format, format=fig_format, bbox_inches='tight')
		plt.close()


		# --------------------------------------------------
		# # --------  Plot Field PDFs
		# ---------------------------------------------------
		###---------- Individual Phases
		print("Plotting --- Phases")
		num_bins = 1000
		norm_hist = False

		in_data = np.mod(np.angle(run_data.u) + 2.0 * np.pi, 2.0 * np.pi)
		fig_name = r"U_Phases"
		data_labs = r"$\arg\{ u_n \}$"

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
		fig.savefig(cmdargs.out_dir_HD  + fig_name + "_PDFs" + "." + fig_format, format=fig_format, bbox_inches='tight')
		plt.close()

		###---------- Individual Amps
		in_data = np.absolute(run_data.u)
		fig_name = r"U_Amps"
		data_labs = r"$| u_n |$"
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
		fig.savefig(cmdargs.out_dir_HD  + fig_name + "_PDFs" + "." + fig_format, format=fig_format, bbox_inches='tight')
		plt.close()

		###################################
		##  PO Model Dynamics Diagnostic
		###################################
		if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
			phases = run_data.phi_n
		else:
			phases = np.mod(np.angle(run_data.u) + 2.0*np.pi, 2.0*np.pi)

		if sys_vars.model_type == "PO" or sys_vars.model_type == "AO" or sys_vars.model_type == "FULL":
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
		# ax1.set_title(sys_vars.model_type)
		ax1.legend()

		plt.savefig(cmdargs.out_dir_HD + "VelReal_PDF_InOne" + fig_file_type)
		plt.close()

		fig = plt.figure(figsize = (16, 16))
		gs  = GridSpec(5, 5, wspace = 0.25, hspace = 0.25)
		for i in range(5):
			for j in range(5):
				indx = i * 5 + j
				if indx < sys_vars.N:
					ax1 = fig.add_subplot(gs[i, j])
					if hasattr(stats_data, "vel_hist_counts"):
						pdf, centres = compute_pdf_from_hist(stats_data.vel_hist_counts[indx, :], stats_data.vel_hist_ranges[indx, :], normed = True)
						ax1.set_title("C Data")
					else:
						if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
							pdf, centres = compute_pdf(np.real(run_data.a_n[:, indx] * np.exp(1j * run_data.phi_n[:, indx])), nbins = nbins, normed = True)					
						else:
							pdf, centres = compute_pdf(np.real(run_data.u[:, indx]), nbins = nbins, normed = True)
					ax1.plot(centres, pdf, label = "$n = {}$".format(indx + 1)) 
					ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
					ax1.set_ylabel(r"PDF")
					ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
					ax1.set_yscale('log')
					ax1.legend()

		plt.savefig(cmdargs.out_dir_HD + "VelReal_PDF_All" + fig_file_type, bbox_inches = 'tight')
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
		print(len(p), len(zeta_p), len(ns_zeta_p))

		with open(cmdargs.out_dir_HD + "ComputedStrFuncSlopes.txt", 'w') as f:
			f.write("Velocity Structure Func\n")
			f.write("Power\tp/3\t\tDNS Slope\tNS\n")
			for i in range(num_pow):
				if i >= 1 and i < len(ns_zeta_p):
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f} \t{:1.3f}\n".format(i + 1, -(i + 1) / 3, zeta_p[i], zeta_p_res[i], -ns_zeta_p[i - 1]))
				else:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, zeta_p[i], zeta_p_res[i]))

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
					if i >= 1 and i < len(ns_zeta_p):
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
				if i >= 1 and i < len(ns_zeta_p):
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
				if i >= 1 and i < len(ns_zeta_p):
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f} \t{:1.3f}\n".format(i + 1, -(i + 1) / 3, hel_flux_zeta_p[i], hel_flux_zeta_p_res[i], -ns_zeta_p[i - 1]))
				else:
					f.write(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f}\n".format(i + 1, -(i + 1) / 3, hel_flux_zeta_p[i], hel_flux_zeta_p_res[i]))


		## --------  Plot Anomalous Exponent
		plot_anomalous_exponent(cmdargs.out_dir_HD + "VelHelicityFluxAbs_Anonalous_Exponent_Zeta_p" + fig_file_type, p, hel_flux_zeta_p[1:], label_str = r"Velocity HelicityFlux; Shell Modell", fig_size = fig_size)


		# Plotting variables
		fig_size = (10, 6)
		fig_format = ".png"

		## Normal Structure Functions
		str_funcs = [stats_data.vel_str_func[:, :num_pow], stats_data.vel_trip_prod_str_func_abs[:, :num_pow], stats_data.vel_flux_str_func_abs[:, :num_pow, 0], stats_data.vel_flux_str_func_abs[:, :num_pow, 1]]
		msg = ["Velocity Modes", "Tripple Product", "Energy Flux", "Helicity Flux"]
		filename = ["Vel", "Trip_Prod", "Enrg_Flux", "Hel_Flux"]
		for sf, m, name in zip(str_funcs, msg, filename):

			## Pre Multiplied Str Funcs 
			log_func    = np.log2
			k_n         = sys_msr_data.k[:sf.shape[0]]
			outdir_path = cmdargs.out_dir_HD + name + "_PreMult_SF" + fig_format

			fig = plt.figure(figsize=fig_size)
			gs  = GridSpec(1, 1)
			ax1 = fig.add_subplot(gs[0, 0])
			for i in range(sf.shape[1]):
			    p, = ax1.plot(log_func(k_n), log_func((k_n**((i + 1) / 3.0)) * sf[:, i]) + i * 10, label = "$p = {}$".format(i + 1))
			    ax1.plot(log_func(k_n), log_func(k_n**0) + i * 10, '--', color = p.get_color())
			ax1.set_xlabel(r"$\ln (k_n)$")
			ax1.set_ylabel(r"$\ln (k_n^{p / 3}|S_p(k_n)|)$")
			ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
			ax1.set_title("Pre Multiplied " + m)
			ax1.legend()
			plt.savefig(outdir_path, bbox_inches='tight')
			plt.close()

			## ESS Plots Combined
			log_func    = np.log2
			outdir_path = cmdargs.out_dir_HD + name + "_ESS_SF" + fig_format
			fig = plt.figure(figsize=fig_size)
			gs  = GridSpec(1, 1)
			ax1 = fig.add_subplot(gs[0, 0])
			for i in range(sf.shape[1]):
				if i == 2: 
					continue
				p, = ax1.plot(log_func(sf[:, 2]), log_func(sf[:, i]), label = "$p = {}$".format(i + 1))
				# ax1.plot(log_func(sf[:, 2]), log_func(she_leveque(sf[:, i])), ':', color=p.get_color())
			ax1.set_xlabel(r"$|S_3(k_n)|$")
			ax1.set_ylabel(r"$|S_p(k_n)|)$")
			ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
			ax1.set_title("ESS " + m)
			ax1.legend()
			plt.savefig(outdir_path, bbox_inches='tight')
			plt.close()

			## Local Slopes
			log_func    = np.log2
			outdir_path = cmdargs.out_dir_HD + name + "_LocalSlope_SF" + fig_format
			fig      = plt.figure(figsize=fig_size)
			gs       = GridSpec(1, 2)
			ax1    = fig.add_subplot(gs[0, 0])
			for i in range(sf.shape[1]):
				x      = log_func(sys_msr_data.k[:len(sf[:, i])])
				y      = log_func(sf[:, i])
				d_str_func = np.diff(y, n = 1)
				d_k        = np.diff(x, n = 1)
				derivs     = d_str_func / d_k
				# z1     = np.hstack((y[0], y[:-1]))
				# z2     = np.hstack((y[1:], y[-1]))
				# dx1    = np.hstack((0, np.diff(x)))
				# dx2    = np.hstack((np.diff(x), 0))
				# derivs = (z2-z1) / (dx2+dx1)
				ax1.plot(log_func(sys_msr_data.k[:len(derivs)]), derivs, label=r"$p = {}$".format(i + 1))
			ax1.set_xlabel(r"$k_n$")
			ax1.set_ylabel(r"$\zeta_p$")
			ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
			ax1.set_title("Local Slope " + m)
			ax1.legend()
			ax1    = fig.add_subplot(gs[0, 1])
			for i in range(sf.shape[1]):
				if i == 2: 
					continue
				x      = log_func(sf[:, 2])
				y      = log_func(sf[:, i])
				d_str_func = np.diff(y, n = 1)
				d_k        = np.diff(x, n = 1)
				derivs     = d_str_func / d_k
				# z1     = np.hstack((y[0], y[:-1]))
				# z2     = np.hstack((y[1:], y[-1]))
				# dx1    = np.hstack((0, np.diff(x)))
				# dx2    = np.hstack((np.diff(x), 0))
				# derivs = (z2-z1) / (dx2+dx1)
				ax1.plot(log_func(sf[:len(derivs), 2]), derivs, label=r"$p = {}$".format(i + 1))
			ax1.set_xlabel(r"$k_n$")
			ax1.set_ylabel(r"$\zeta_p$")
			ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
			ax1.set_title("ESS Local Slope " + m)
			ax1.legend()
			plt.savefig(outdir_path, bbox_inches='tight')
			plt.close()

		###################################
		##  Structure Function Contribution Breakdown
		###################################
		fig = plt.figure(figsize=(15, 6))
		gs  = GridSpec(1, 3, hspace=0.4, wspace=0.3)
		ax1 = fig.add_subplot(gs[0, 0])
		for p in range(1, 6 + 1):
			ax1.plot(sys_msr_data.k[:-2], np.mean(np.power(np.absolute(run_data.enrg_flux[:, :-2]), p/3.), axis=0), label=r"$Abs$; $p = {}$".format(p))
		ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
		ax1.set_yscale('log')
		ax1.set_xscale('log')
		ax1.legend()
		ax1 = fig.add_subplot(gs[0, 1])
		for p in range(1, 6 + 1):
			ax1.plot(sys_msr_data.k[:-2], np.mean(np.power(np.sign(np.sin(np.angle(run_data.enrg_flux[:, :-2]))), p) * np.power(np.absolute(np.sin(np.angle(run_data.enrg_flux[:, :-2]))), p/3.), axis=0), label=r"$\sin \arg$; $p = {}$".format(p))
		ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
		ax1.set_xscale('log')
		ax1.set_yscale('log')
		ax1.legend()
		ax1 = fig.add_subplot(gs[0, 2])
		print(run_data.enrg_flux[:, :-2])
		for p in range(1, 6 + 1):
			func = np.mean(np.power(np.absolute(np.sin(np.angle(run_data.enrg_flux[:, :-2]))), p/3.), axis=0)
			ax1.plot(sys_msr_data.k[:-2], func / func[0], label=r"$\arg $; $p = {}$".format(p))
			# ax1.plot(k, np.mean( np.power(np.sign(np.angle(run_data.enrg_flux[:, :])), p) * np.power(np.absolute(np.angle(run_data.enrg_flux[:, :])), p/3.), axis=0), label=r"$\arg $; $p = {}$".format(p))
		ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
		ax1.set_xscale('log')
		ax1.set_yscale('log')
		ax1.legend()
		plt.savefig(cmdargs.out_dir_HD + "EnrgyFlux" + "_Test_Independence_SF" + ".png", bbox_inches='tight')
		plt.close()

		###################################
		##  Anomalous Exponent
		###################################
		inertial_range = np.arange(7, 13 + 1)
		p = np.arange(2, num_pow + 1)

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
