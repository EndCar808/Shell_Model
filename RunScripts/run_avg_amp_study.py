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
# mpl.use('Agg') # Use this backend for writing plots to file
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
# from Plotting.functions import parse_cml, tc, sim_data, import_data, import_stats_data, import_sys_msr_data
sys.path.append('/home/enda/PhD/Shell_Model/Plotting')
from functions import parse_cml, tc, sim_data, import_data, import_stats_data, import_sys_msr_data, compute_pdf, compute_u_flux, compute_str_func

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
	cmdargs.out_dir_AVGAMP = cmdargs.out_dir + "AVGAMP_PLOTS/"
	if os.path.isdir(cmdargs.out_dir_AVGAMP) != True:
		print("Making folder:" + tc.C + " AVGAMP_PLOTS/" + tc.Rst)
		os.mkdir(cmdargs.out_dir_AVGAMP)

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
	# # --------  Generate Averaged Amp Field
	# -----------------------------------------
	if hasattr(sys_msr_data, 'a_n_t_avg'):
		## Test the average of amplitudes is correct
		avg_amp = np.mean(np.absolute(run_data.u), axis = 0)
		
		print(np.allclose(avg_amp, sys_msr_data.a_n_t_avg))
		
		plt.figure()
		plt.plot(avg_amp, label = "python")
		plt.plot(sys_msr_data.a_n_t_avg, label = "C")
		plt.xscale('log')
		plt.yscale('log')
		plt.legend()
		plt.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "Test_Avg_Amp.png", bbox_inches='tight')
		plt.close()

		## Create average amp field
		u_avg_amp_all = sys_msr_data.a_n_t_avg * np.exp(1j * np.angle(run_data.u))
		u_avg_amp     = avg_amp * np.exp(1j * np.angle(run_data.u))
		print(u_avg_amp.shape)

		## Compute PDFs
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		for j, i in enumerate([5 - 1, 10 - 1, 15 - 1, 20 - 1]):
			pdf, centres 		 = compute_pdf(np.real(u_avg_amp[:, i]), nbins = 500, normed = False)
			pdf_all, centres_all = compute_pdf(np.real(u_avg_amp_all[:, i]), nbins = 500, normed = False)
			if i == -1:
				ax1.plot(centres, pdf, label = "$n = {}$".format(sys_vars.N))    
				ax1.plot(centres_all, pdf_all, label = "All; $n = {}$".format(sys_vars.N))    
			else:
				ax1.plot(centres, pdf, label = "$n = {}$".format(i + 1))    
				ax1.plot(centres_all, pdf_all, label = "All; $n = {}$".format(i + 1))    
		ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
		ax1.set_ylabel(r"PDF")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.set_yscale('log')
		ax1.legend()
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_VelReal_PDF_InOne.png")
		plt.close()

		## Loop through time and compute quantities for testing
		avg_u_enrg_flux          = np.zeros((sys_vars.N,))
		avg_u_hel_flux           = np.zeros((sys_vars.N,))
		str_func_avg_u           = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
		str_func_avg_u_enrg_flux = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
		str_func_avg_u_hel_flux  = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
		avg_u_enrg_flux_all          = np.zeros((sys_vars.N,))
		avg_u_hel_flux_all           = np.zeros((sys_vars.N,))
		str_func_avg_u_all           = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
		str_func_avg_u_enrg_flux_all = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
		str_func_avg_u_hel_flux_all  = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
		num_t_data = sys_vars.ndata
		for t in range(num_t_data):
			print(t)

			## Get fluxes
			tmp_enrg, tmp_hel = compute_u_flux(u_avg_amp[t, :], sys_vars.EPS, sys_vars.Lambda)
			avg_u_enrg_flux   += tmp_enrg
			avg_u_hel_flux    += tmp_hel
			tmp_enrg, tmp_hel   = compute_u_flux(u_avg_amp_all[t, :], sys_vars.EPS, sys_vars.Lambda)
			avg_u_enrg_flux_all += tmp_enrg
			avg_u_hel_flux_all  += tmp_hel

			## Get Structure Functions
			tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_avg_amp[t, :], stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
			str_func_avg_u                     += tmp_u
			str_func_avg_u_enrg_flux           += tmp_enrg_flux
			str_func_avg_u_hel_flux            += tmp_hel_flux
			tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_avg_amp_all[t, :], stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
			str_func_avg_u_all                 += tmp_u
			str_func_avg_u_enrg_flux_all       += tmp_enrg_flux
			str_func_avg_u_hel_flux_all        += tmp_hel_flux

		## Average
		avg_u_enrg_flux          /= num_t_data
		avg_u_hel_flux           /= num_t_data
		str_func_avg_u           /= num_t_data
		str_func_avg_u_enrg_flux /= num_t_data
		str_func_avg_u_hel_flux  /= num_t_data
		avg_u_enrg_flux_all          /= num_t_data
		avg_u_hel_flux_all           /= num_t_data
		str_func_avg_u_all           /= num_t_data
		str_func_avg_u_enrg_flux_all /= num_t_data
		str_func_avg_u_hel_flux_all  /= num_t_data
		mark_style = ['o','s','^','x','D','p', '>', '<']

		## Plot Time Averaged Flux Scaling
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(sys_msr_data.k, np.absolute(sys_msr_data.enrg_flux_t_avg), label = "Pre-Multiplied Time Averaged Energy Flux Spectrum")
		ax1.plot(sys_msr_data.k, sys_msr_data.k * np.absolute(avg_u_enrg_flux), marker = mark_style[1], label = "Averaged Amps")
		ax1.plot(sys_msr_data.k, sys_msr_data.k * np.absolute(avg_u_enrg_flux_all), marker = mark_style[1], label = "All Averaged Amps")
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$k_n \mathcal{E}_n$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax2 = fig.add_subplot(gs[0, 1])
		ax2.plot(sys_msr_data.k, np.absolute(sys_msr_data.enrg_flux_t_avg) / sys_msr_data.k, label = "Time Averaged Energy Flux Spectrum")
		ax2.plot(sys_msr_data.k, np.absolute(avg_u_enrg_flux), marker = mark_style[1], label = "Averaged Amps")
		ax2.plot(sys_msr_data.k, np.absolute(avg_u_enrg_flux_all), marker = mark_style[1], label = "All Averaged Amps")
		ax2.plot(sys_msr_data.k, sys_msr_data.k ** (-1), 'k--', label = "$k^{-1}$")
		ax2.set_xlabel(r"$k_n$")
		ax2.set_ylabel(r"$\mathcal{E}_n$")
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Time_Averaged_Energy_Flux.png", bbox_inches='tight')
		plt.close()


		## Plot Str Funcs - Velocity
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		p = 2
		ax1.plot(sys_msr_data.k, stats_data.vel_str_func[:, p - 1], label = "Full Field - p = {}".format(p))
		ax1.plot(sys_msr_data.k, str_func_avg_u[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
		ax1.plot(sys_msr_data.k, str_func_avg_u_all[:, p - 1], marker = mark_style[1], label = "All Averaged Amps - p = {}".format(p))
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$\ln \mathcal{S}_p (|u_n|)$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.set_title("2nd Order Velocity Str Func")
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		p = 6
		ax2 = fig.add_subplot(gs[0, 1])
		ax2.plot(sys_msr_data.k, stats_data.vel_str_func[:, p - 1], label = "Full Field - p = {}".format(p))
		ax2.plot(sys_msr_data.k, str_func_avg_u[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
		ax2.plot(sys_msr_data.k, str_func_avg_u_all[:, p - 1], marker = mark_style[1], label = "All Averaged Amps - p = {}".format(p))
		ax2.set_xlabel(r"$k_n$")
		ax2.set_ylabel(r"$\ln \mathcal{S}_p (|u_n|)$")
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		ax2.set_title("6nd Order Velocity Str Func")
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Velocity.png", bbox_inches='tight')
		plt.close()


		## Plot Str Funcs - Energy Flux
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		p = 2
		ax1.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 0], label = "Full Field - p = {}".format(p))
		ax1.plot(sys_msr_data.k, str_func_avg_u_enrg_flux[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
		ax1.plot(sys_msr_data.k, str_func_avg_u_enrg_flux_all[:, p - 1], marker = mark_style[1], label = "All Averaged Amps - p = {}".format(p))
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.set_title("2nd Order Energy Flux Str Func")
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		p = 6
		ax2 = fig.add_subplot(gs[0, 1])
		ax2.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 0], label = "Full Field - p = {}".format(p))
		ax2.plot(sys_msr_data.k, str_func_avg_u_enrg_flux[:, p - 1], marker = mark_style[1], label = "All Averaged Amps - p = {}".format(p))
		ax2.plot(sys_msr_data.k, str_func_avg_u_enrg_flux_all[:, p - 1], marker = mark_style[1], label = "All Averaged Amps - p = {}".format(p))
		ax2.set_xlabel(r"$k_n$")
		ax2.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		ax2.set_title("6nd Order Energy Flux Str Func")
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Energy_Flux.png", bbox_inches='tight')
		plt.close()


		## Plot Str Funcs - Helicity Flux
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		p = 2
		ax1.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 1], label = "Full Field - p = {}".format(p))
		ax1.plot(sys_msr_data.k, str_func_avg_u_hel_flux[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
		ax1.plot(sys_msr_data.k, str_func_avg_u_hel_flux_all[:, p - 1], marker = mark_style[1], label = "All Averaged Amps - p = {}".format(p))
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.set_title("2nd Order Helicity Flux Str Func")
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		p = 6
		ax2 = fig.add_subplot(gs[0, 1])
		ax2.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 1], label = "Full Field - p = {}".format(p))
		ax2.plot(sys_msr_data.k, str_func_avg_u_hel_flux[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
		ax2.plot(sys_msr_data.k, str_func_avg_u_hel_flux_all[:, p - 1], marker = mark_style[1], label = "All Averaged Amps - p = {}".format(p))
		ax2.set_xlabel(r"$k_n$")
		ax2.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		ax2.set_title("6nd Order Helicity Flux Str Func")
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Helicity_Flux.png", bbox_inches='tight')
		plt.close()