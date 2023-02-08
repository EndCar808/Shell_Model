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


def slope_fit(x, y, low, high):

	poly_output = np.polyfit(x[low:high], y[low:high], 1, full = True)
	pfit_info   = poly_output[0]
	poly_resid  = poly_output[1][0]
	pfit_slope  = pfit_info[0]
	pfit_c      = pfit_info[1]

	return pfit_slope, pfit_c, poly_resid
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
		mark_style = ['o','s','^','x','D','p', '>', '<']

		# shell 6 - 19
		inert_lim_low  = 11
		inert_lim_high = 16

		num_pow = stats_data.vel_str_func.shape[-1]


		## Test the average of amplitudes is correct
		avg_amp = np.mean(np.absolute(run_data.u), axis = 0)
		
		print(np.allclose(avg_amp, sys_msr_data.a_n_t_avg), avg_amp - sys_msr_data.a_n_t_avg)
		
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(sys_msr_data.k, avg_amp, label = "python")
		ax1.plot(sys_msr_data.k, sys_msr_data.a_n_t_avg, label = "C")
		ax1.set_xscale('log')
		ax1.set_yscale('log')
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "Test_Avg_Amp.png", bbox_inches='tight')
		plt.close()

		## Spectra
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.enrg_spec_t_avg), label = "$E_n$", marker = mark_style[0], markevery = 1)
		poly_output = np.polyfit(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.enrg_spec_t_avg[inert_lim_low:inert_lim_high]), 1, full = True)
		pfit_info   = poly_output[0]
		poly_resid  = poly_output[1][0]
		pfit_slope  = pfit_info[0]
		pfit_c      = pfit_info[1]
		ax1.plot(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high])*pfit_slope + pfit_c, 'k--', label = "$k^-{}$".format(np.absolute(pfit_slope)))
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Spectra.png", bbox_inches='tight')
		plt.close()


		## Amp
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.a_n_t_avg), label = "$a_n$", marker = mark_style[0], markevery = 1)
		poly_output = np.polyfit(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.a_n_t_avg[inert_lim_low:inert_lim_high]), 1, full = True)
		pfit_info   = poly_output[0]
		poly_resid  = poly_output[1][0]
		pfit_slope  = pfit_info[0]
		pfit_c      = pfit_info[1]
		ax1.plot(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high])*pfit_slope + pfit_c, 'k--', label = "$k^-{}$".format(np.absolute(pfit_slope)))
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Amplitudes.png", bbox_inches='tight')
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
				p, = ax1.plot(centres, pdf, label = "$n = {}$".format(sys_vars.N), marker = mark_style[0])    
				ax1.plot(centres_all, pdf_all, label = "All; $n = {}$".format(sys_vars.N), color = p.get_color(), marker = mark_style[1])    
			else:
				p, = ax1.plot(centres, pdf, label = "$n = {}$".format(i + 1), marker = mark_style[0])    
				ax1.plot(centres_all, pdf_all, label = "All; $n = {}$".format(i + 1), color = p.get_color(), marker = mark_style[1])    
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
		u_pad           = np.zeros((sys_vars.N + 4, ), dtype = "complex128")
		u_pad_all       = np.zeros((sys_vars.N + 4, ), dtype = "complex128")

		num_t_data = sys_vars.ndata
		for t in range(num_t_data):
			print(t)

			## Get padded data
			u_pad[2    :sys_vars.N + 2] = u_avg_amp[t, :]
			u_pad_all[2:sys_vars.N + 2] = u_avg_amp_all[t, :]

			## Get fluxes
			tmp_enrg, tmp_hel = compute_u_flux(u_pad, sys_vars.N, sys_vars.EPS, sys_vars.Lambda)
			avg_u_enrg_flux   += tmp_enrg
			avg_u_hel_flux    += tmp_hel
			tmp_enrg, tmp_hel   = compute_u_flux(u_pad_all, sys_vars.N, sys_vars.EPS, sys_vars.Lambda)
			avg_u_enrg_flux_all += tmp_enrg
			avg_u_hel_flux_all  += tmp_hel

			## Get Structure Functions
			tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_pad, sys_vars.N, stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
			str_func_avg_u                     += tmp_u
			str_func_avg_u_enrg_flux           += tmp_enrg_flux
			str_func_avg_u_hel_flux            += tmp_hel_flux
			tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_pad_all, sys_vars.N, stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
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
		

		## Plot Time Averaged Energy Flux Scaling
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(sys_msr_data.k, sys_msr_data.k * np.absolute(sys_msr_data.enrg_flux_t_avg), label = "Pre-Multiplied Time Averaged Energy Flux Spectrum")
		ax1.plot(sys_msr_data.k, sys_msr_data.k * np.absolute(avg_u_enrg_flux), marker = mark_style[1], label = "Averaged Amps")
		ax1.plot(sys_msr_data.k, sys_msr_data.k * np.absolute(avg_u_enrg_flux_all), marker = mark_style[2], label = "All Averaged Amps")
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$k_n \mathcal{E}_n$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax2 = fig.add_subplot(gs[0, 1])
		ax2.plot(sys_msr_data.k, np.absolute(sys_msr_data.enrg_flux_t_avg), label = "Time Averaged Energy Flux Spectrum")
		ax2.plot(sys_msr_data.k, np.absolute(avg_u_enrg_flux), marker = mark_style[1], label = "Averaged Amps")
		slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(avg_u_enrg_flux)), inert_lim_low, inert_lim_high)
		ax2.plot(sys_msr_data.k, np.absolute(avg_u_enrg_flux_all), marker = mark_style[2], label = "All Averaged Amps")
		ax2.plot(sys_msr_data.k, sys_msr_data.k ** (-1), 'k--', label = "$k^{-1}$")
		ax2.plot(sys_msr_data.k, sys_msr_data.k ** (slope), '--', label = "{}".format(slope))
		ax2.set_xlabel(r"$k_n$")
		ax2.set_ylabel(r"$\mathcal{E}_n$")
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Time_Averaged_Energy_Flux.png", bbox_inches='tight')
		plt.close()


		## Plot Time Averaged Helicity Flux Scaling
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 2)
		ax1 = fig.add_subplot(gs[0, 0])
		tmp = sys_msr_data.hel_flux_t_avg
		print(tmp)
		tmp1 = np.copy(sys_msr_data.hel_flux_t_avg)
		tmp2 = np.copy(sys_msr_data.hel_flux_t_avg)
		tmp1[tmp1 > 0] = 0.0
		tmp2[tmp2 < 0] = 0.0
		ax1.plot(sys_msr_data.k, np.absolute(tmp1), '.', label = "Neg")
		ax1.plot(sys_msr_data.k, tmp2, 'o', label = "Pos")
		# ax1.plot(sys_msr_data.k, np.absolute(avg_u_hel_flux), marker = mark_style[1], label = "Averaged Amps")
		# ax1.plot(sys_msr_data.k, np.absolute(avg_u_hel_flux_all), marker = mark_style[2], label = "All Averaged Amps")
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$k_n \mathcal{E}_n$")
		ax1.set_xscale("log")
		ax1.set_yscale("log")
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax2 = fig.add_subplot(gs[0, 1])
		tmp1 = np.copy(avg_u_hel_flux)
		tmp2 = np.copy(avg_u_hel_flux)
		tmp1[tmp1 > 0] = 0.0
		tmp2[tmp2 < 0] = 0.0
		ax2.plot(sys_msr_data.k, np.absolute(tmp1), '.', label = "Neg")
		ax2.plot(sys_msr_data.k, tmp2, '.', label = "Pos")
		ax2.plot(sys_msr_data.k, np.absolute(sys_msr_data.hel_flux_t_avg), label = "Time Averaged Helicity Flux Spectrum")
		# ax2.plot(sys_msr_data.k, np.absolute(avg_u_hel_flux_all), marker = mark_style[2], label = "All Averaged Amps")
		ax2.plot(sys_msr_data.k, sys_msr_data.k ** (-2), 'k--', label = "$k^{-2}$")
		ax2.set_xlabel(r"$k_n$")
		ax2.set_ylabel(r"$\mathcal{E}_n$")
		ax2.set_xscale("log")
		ax2.set_yscale("log")
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Time_Averaged_Helicity_Flux.png", bbox_inches='tight')
		plt.close()


		## Plot Str Funcs - Velocity
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(2, 3)
		for i in range(2):
			for j in range(3):
				ax1 = fig.add_subplot(gs[i, j])
				p = i * 3 + j + 1
				ax1.plot(sys_msr_data.k, stats_data.vel_str_func[:, p - 1] / sys_vars.ndata, label = "Full Field - p = {}".format(p))
				ax1.plot(sys_msr_data.k, str_func_avg_u[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
				ax1.plot(sys_msr_data.k, str_func_avg_u_all[:, p - 1], marker = mark_style[2], label = "All Averaged Amps - p = {}".format(p))
				ax1.set_xlabel(r"$k_n$")
				ax1.set_ylabel(r"$\ln \mathcal{S}_p (|u_n|)$")
				ax1.set_xscale("log")
				ax1.set_yscale("log")
				ax1.set_title("{}nd Order Velocity Str Func".format(p))
				ax1.legend()
				ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Velocity.png", bbox_inches='tight')
		plt.close()

		## Plot Str Funcs - Energy Flux
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(2, 3)
		for i in range(2):
			for j in range(3):
				ax1 = fig.add_subplot(gs[i, j])
				p = i * 3 + j + 1
				ax1.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 0] / sys_vars.ndata, label = "Full Field - p = {}".format(p))
				ax1.plot(sys_msr_data.k, str_func_avg_u_enrg_flux[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
				ax1.plot(sys_msr_data.k, str_func_avg_u_enrg_flux_all[:, p - 1], marker = mark_style[2], label = "All Averaged Amps - p = {}".format(p))
				ax1.set_xlabel(r"$k_n$")
				ax1.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
				ax1.set_xscale("log")
				ax1.set_yscale("log")
				ax1.set_title("{}nd Order Energy Flux Str Func".format(p))
				ax1.legend()
				ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Energy_Flux.png", bbox_inches='tight')
		plt.close()
		
		## Plot Str Funcs - Helicity Flux
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(2, 3)
		for i in range(2):
			for j in range(3):
				ax1 = fig.add_subplot(gs[i, j])
				p = i * 3 + j + 1
				ax1.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 1] / sys_vars.ndata, label = "Full Field - p = {}".format(p))
				ax1.plot(sys_msr_data.k, str_func_avg_u_hel_flux[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
				ax1.plot(sys_msr_data.k, str_func_avg_u_hel_flux_all[:, p - 1], marker = mark_style[2], label = "All Averaged Amps - p = {}".format(p))
				ax1.set_xlabel(r"$k_n$")
				ax1.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
				ax1.set_xscale("log")
				ax1.set_yscale("log")
				ax1.set_title("{}nd Order Helicity Flux Str Func".format(p))
				ax1.legend()
				ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Helicity_Flux.png", bbox_inches='tight')
		plt.close()
		





		## Get adjusted time averaged amplitudes
		a_n = sys_msr_data.a_n_t_avg

		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		a_n_slope, a_n_c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(a_n), inert_lim_low, inert_lim_high)
		ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.a_n_t_avg), label = "$a_n \sim k^{}$".format(np.around(a_n_slope, 3)), marker = None, markevery = 1)
		for beta in [0., 0.3, 0.6, 0.9, 1.2, 1.5]:
			a_n_adjust = a_n * (sys_msr_data.k ** (np.absolute(a_n_slope) - beta))
			p, = ax1.plot(np.log2(sys_msr_data.k), np.log2(a_n_adjust), marker = mark_style[1], markevery = 1)
			slope, c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(a_n_adjust), inert_lim_low, inert_lim_high)
			ax1.plot(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high])*slope + c, '--', label = r"$k^a \quad a = $ {:0.3f}".format(slope), color = p.get_color())
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_All_Slopes.png", bbox_inches='tight')
		plt.close()






		## Get adjusted time averaged amplitudes
		a_n          = sys_msr_data.a_n_t_avg
		num_t_data   = sys_vars.ndata
		phases       = np.angle(run_data.u)
		alpha_slopes = np.arange(0.0, 1.51, 0.01)

		## alocate memory
		u_pad = np.zeros((sys_vars.N + 4, ), dtype = "complex128")
		enrg_flux_slope          = np.zeros((len(alpha_slopes),))
		hel_flux_slope           = np.zeros((len(alpha_slopes),))
		str_func_u_slope         = np.zeros((len(alpha_slopes), num_pow))
		str_func_enrg_flux_slope = np.zeros((len(alpha_slopes), num_pow))
		str_func_hel_flux_slope  = np.zeros((len(alpha_slopes), num_pow))

		for i, alpha in enumerate(alpha_slopes):
			print(i)

			## Alocate memory
			avg_u_enrg_flux          = np.zeros((sys_vars.N,))
			avg_u_hel_flux           = np.zeros((sys_vars.N,))
			str_func_avg_u           = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
			str_func_avg_u_enrg_flux = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
			str_func_avg_u_hel_flux  = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))

			## Create average amp field
			a_n_adjust = a_n * (sys_msr_data.k ** (np.absolute(a_n_slope) - alpha))
			u_avg_amp  = a_n_adjust * np.exp(1j * phases)
			
			## Loop through time and compute quantities for testing
			for t in range(num_t_data):

				## Get padded data
				u_pad[2    :sys_vars.N + 2] = u_avg_amp[t, :]

				## Get fluxes
				tmp_enrg, tmp_hel = compute_u_flux(u_pad, sys_vars.N, sys_vars.EPS, sys_vars.Lambda)
				avg_u_enrg_flux   += tmp_enrg
				avg_u_hel_flux    += tmp_hel

				## Get Structure Functions
				tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_pad, sys_vars.N, stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
				str_func_avg_u                     += tmp_u
				str_func_avg_u_enrg_flux           += tmp_enrg_flux
				str_func_avg_u_hel_flux            += tmp_hel_flux

			## Average
			avg_u_enrg_flux          /= num_t_data
			avg_u_hel_flux           /= num_t_data
			str_func_avg_u           /= num_t_data
			str_func_avg_u_enrg_flux /= num_t_data
			str_func_avg_u_hel_flux  /= num_t_data

			## Get scaling slopes
			enrg_flux_slope[i], c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(np.absolute(avg_u_enrg_flux)), inert_lim_low, inert_lim_high)
			hel_flux_slope[i], c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(np.absolute(avg_u_hel_flux)), inert_lim_low, inert_lim_high)
			for p in range(num_pow):
				str_func_enrg_flux_slope[i, p], c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(str_func_avg_u_hel_flux[:, p]), inert_lim_low, inert_lim_high)
				str_func_hel_flux_slope[i, p], c, _  = slope_fit(np.log2(sys_msr_data.k), np.log2(str_func_avg_u_enrg_flux[:, p]), inert_lim_low, inert_lim_high)
				str_func_u_slope[i, p], c, _         = slope_fit(np.log2(sys_msr_data.k), np.log2(str_func_avg_u[:, p]), inert_lim_low, inert_lim_high)


		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 1)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(alpha_slopes, enrg_flux_slope, label = "Energy Flux")
		ax1.plot(alpha_slopes, hel_flux_slope, label = "Helicity Flux")
		ax1.set_xlabel(r"$\alpha$")
		ax1.set_ylabel(r"Slopes")
		ax1.legend()
		# ax1.set_xlim(alpha_slopes[0], alpha_slopes[-1])
		plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_EnergyFlux_Slopes_vs_Alpha.png", bbox_inches='tight')
		plt.close()