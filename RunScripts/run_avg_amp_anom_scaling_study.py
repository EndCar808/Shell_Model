#!/usr/bin/env python

# Author: Enda Carroll
# Date: Sept 2021
# Info: Script to compare solver results with decaying turbulence papers
#       Solver data

#######################
#  Library Imports  ##
#######################
import numpy as np
import h5py
import sys
import os
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
sys.path.append('/home/enda/PhD/Shell_Model/Plotting')
from functions import parse_cml, tc, sim_data, compute_field_values, compute_str_func_field_values, slope_fit, she_leveque
from plot_functions import plot_str_funcs_with_slope, plot_anomalous_exponent
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import cm

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Computer Modern Roman'

np.seterr(divide = 'ignore') 

######################
#       MAIN       ##
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

	# Make output folder for plots
	cmdargs.out_dir_ANOM_SCALING = cmdargs.out_dir + "ANOM_SCALING_PLOTS/"
	if os.path.isdir(cmdargs.out_dir_ANOM_SCALING) != True:
	    print("Making folder:" + tc.C + " ANOM_SCALING_PLOTS/" + tc.Rst)
	    os.mkdir(cmdargs.out_dir_ANOM_SCALING)



	# -----------------------------------------
	# # --------  Read In data
	# -----------------------------------------
	# Read in simulation parameters
	sys_vars = sim_data(data_file_path, method)

	with h5py.File(cmdargs.in_dir + "Main_HDF_Data.h5", 'r') as in_file:

		if "VelModes" in list(in_file.keys()):
			u = in_file["VelModes"][:, :]

	with h5py.File(cmdargs.in_dir + "System_Measure_HDF_Data.h5", 'r') as in_file:

		if "k" in list(in_file.keys()):
			k = in_file["k"][:]
		if "TimeAveragedVelocityAmplitudes" in list(in_file.keys()):
			a_n_t_avg = in_file["TimeAveragedVelocityAmplitudes"][:]


	# -----------------------------------------
	# # --------  Compute the Needed Data
	# -----------------------------------------
	# Get the number of time data points
	num_t_data = sys_vars.ndata
	num_pow    = 6

	u_sf                 = np.zeros((sys_vars.N, num_pow))
	trip_prod_sf         = np.zeros((sys_vars.N, num_pow))
	doub_prod_sf         = np.zeros((sys_vars.N - 1, num_pow))
	hel_flux_sf          = np.zeros((sys_vars.N, num_pow))
	enrg_flux_sf         = np.zeros((sys_vars.N, num_pow))
	avg_amp_u_sf         = np.zeros((sys_vars.N, num_pow))
	avg_amp_trip_prod_sf = np.zeros((sys_vars.N, num_pow))
	avg_amp_doub_prod_sf = np.zeros((sys_vars.N - 1, num_pow))
	avg_amp_hel_flux_sf  = np.zeros((sys_vars.N, num_pow))
	avg_amp_enrg_flux_sf = np.zeros((sys_vars.N, num_pow))

	avg_amp   = np.mean(np.absolute(u[:, :]), axis=0)

	
	# Loop through simulation and compute field values
	for t in range(num_t_data):
		# print(t)

		# Get padded field 
		u_pad = np.pad(u[t, :], 2, "constant")

		# Get field values
		trip_prod, _, doub_prod, hel_flux, enrg_flux = compute_field_values(u_pad, sys_vars.N, sys_vars.eps, sys_vars.Lambda)
		u_tmp, trip_tmp, dub_tmp, enrg_tmp, hel_tmp  = compute_str_func_field_values(u_pad, trip_prod, doub_prod, enrg_flux, hel_flux, sys_vars.N, num_pow)
		u_sf[:, :]         += u_tmp
		trip_prod_sf[:, :] += trip_tmp
		doub_prod_sf[:, :] += dub_tmp
		enrg_flux_sf[:, :] += enrg_tmp
		hel_flux_sf[:, :]  += hel_tmp

		# Get padded field 
		# u_avg_amp = avg_amp * np.exp(1j * np.angle(u[t, :]))
		u_avg_amp = a_n_t_avg * np.exp(1j * np.angle(u[t, :]))
		u_pad_avg = np.pad(u_avg_amp, 2, "constant")


		# Get field values
		trip_prod, _, doub_prod, hel_flux, enrg_flux = compute_field_values(u_pad_avg, sys_vars.N, sys_vars.eps, sys_vars.Lambda)
		u_tmp, trip_tmp, dub_tmp, enrg_tmp, hel_tmp  = compute_str_func_field_values(u_pad_avg, trip_prod, doub_prod, enrg_flux, hel_flux, sys_vars.N, num_pow)
		avg_amp_u_sf[:, :]         += u_tmp
		avg_amp_trip_prod_sf[:, :] += trip_tmp
		avg_amp_doub_prod_sf[:, :] += dub_tmp
		avg_amp_enrg_flux_sf[:, :] += enrg_tmp
		avg_amp_hel_flux_sf[:, :]  += hel_tmp

	## Normalize all
	u_sf[                :, :] /= num_t_data
	trip_prod_sf[        :, :] /= num_t_data
	doub_prod_sf[        :, :] /= num_t_data
	enrg_flux_sf[        :, :] /= num_t_data
	hel_flux_sf[         :, :] /= num_t_data
	avg_amp_u_sf[        :, :] /= num_t_data
	avg_amp_trip_prod_sf[:, :] /= num_t_data
	avg_amp_doub_prod_sf[:, :] /= num_t_data
	avg_amp_enrg_flux_sf[:, :] /= num_t_data
	avg_amp_hel_flux_sf[ :, :] /= num_t_data

	# -----------------------------------------
	# # --------  Plot Data
	# -----------------------------------------
	# Plotting variables
	fig_size = (10, 6)
	fig_format = ".png"

	## Normal Structure Functions
	str_funcs = [u_sf, trip_prod_sf, doub_prod_sf[:-2, :], enrg_flux_sf, hel_flux_sf]
	msg = ["Velocity Modes", "Tripple Product", "Double Product", "Energy Flux", "Helicity Flux"]
	filename = ["Vel", "Trip_Prod", "Dub_Prod", "Enrg_Flux", "Hel_Flux"]
	# ranges = [np.arange(6, 18)]
	for sf, m, name in zip(str_funcs, msg, filename):

		## Plot SFs with fit
		print(m)
		inert_range = np.arange(3, 16)
		str_funcs   = sf
		scale       = 'log2'
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_SF_withFit" + fig_format
		zeta_p, ns_zeta_p, zeta_p_resid = plot_str_funcs_with_slope(outdir_path, k, str_funcs, inert_range, insert_fig=True, scaling=scale, fig_size=fig_size)

		## Anomalous Exponent
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_Anomalous_Exp" + fig_format
		plot_anomalous_exponent(outdir_path, np.arange(1, num_pow + 1), zeta_p, ns_zeta_p, label_str=m, fig_size=fig_size)


		## Pre Multiplied Str Funcs 
		log_func = np.log2
		k_n = k[:sf.shape[0]]
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_PreMult_SF" + fig_format

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
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_ESS_SF" + fig_format
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
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_LocalSlope_SF" + fig_format
		fig      = plt.figure(figsize=fig_size)
		gs       = GridSpec(1, 2)
		# local_deriv = np.concatenate((local_deriv, [(log_func(str_funcs[-1, i]) - 2.0 * log_func(str_funcs[-2, i]) + log_func(str_funcs[-3, i])) / (log_func(k[-1]) - log_func(k[-2]))**2, (log_func(str_funcs[-1, i]) - 2.0 * log_func(str_funcs[-2, i]) + log_func(str_funcs[-3, i])) / (log_func(k[-1]) - log_func(k[-2]))**2]))
		ax1    = fig.add_subplot(gs[0, 0])
		for i in range(sf.shape[1]):
			x      = log_func(k[:len(sf[:, i])])
			y      = log_func(sf[:, i])
			d_str_func = np.diff(y, n = 1)
			d_k        = np.diff(x, n = 1)
			derivs     = d_str_func / d_k
			# z1     = np.hstack((y[0], y[:-1]))
			# z2     = np.hstack((y[1:], y[-1]))
			# dx1    = np.hstack((0, np.diff(x)))
			# dx2    = np.hstack((np.diff(x), 0))
			# derivs = (z2-z1) / (dx2+dx1)
			ax1.plot(x, derivs, label=r"$p = {}$".format(i + 1))
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
			ax1.plot(x, derivs, label=r"$p = {}$".format(i + 1))
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$\zeta_p$")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.set_title("ESS Local Slope " + m)
		ax1.legend()
		plt.savefig(outdir_path, bbox_inches='tight')
		plt.close()

	## Normal Structure Functions
	str_funcs = [avg_amp_u_sf, avg_amp_trip_prod_sf, avg_amp_doub_prod_sf[:-2, :], avg_amp_enrg_flux_sf, avg_amp_hel_flux_sf]
	msg = ["Avg Amp Velocity Modes", "Avg Amp Tripple Product", "Avg Amp Double Product", "Avg Amp Energy Flux", "Avg Amp Helicity Flux"]
	filename = ["AvgAmp_Vel", "AvgAmp_Trip_Prod", "AvgAmp_Dub_Prod", "AvgAmp_Enrg_Flux", "AvgAmp_Hel_Flux"]
	# ranges = [np.arange(6, 18)]
	for sf, m, name in zip(str_funcs, msg, filename):

		## Plot SFs with fit
		print(m)
		inert_range = np.arange(3, 16)
		str_funcs   = sf
		scale       = 'loge'
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_SF_withFit" + fig_format
		zeta_p, ns_zeta_p, zeta_p_resid = plot_str_funcs_with_slope(outdir_path, k, str_funcs, inert_range, insert_fig=True, scaling=scale, fig_size=fig_size)

		## Anomalous Exponent
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_Anomalous_Exp" + fig_format
		plot_anomalous_exponent(outdir_path, np.arange(1, num_pow + 1), zeta_p, ns_zeta_p, label_str=m, fig_size=fig_size)

		# Get the slope of the str function slopes as a function of p
		slope, c, res = slope_fit(np.arange(1, num_pow + 1), zeta_p, 0, len(zeta_p))
		print(slope, res)

		## Pre Multiplied Str Funcs 
		log_func = np.log
		k_n = k[:sf.shape[0]]
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_PreMult_SF" + fig_format

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
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_ESS_SF" + fig_format
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
		outdir_path = cmdargs.out_dir_ANOM_SCALING + name + "_LocalSlope_SF" + fig_format
		fig      = plt.figure(figsize=fig_size)
		gs       = GridSpec(1, 2)
		# local_deriv = np.concatenate((local_deriv, [(log_func(str_funcs[-1, i]) - 2.0 * log_func(str_funcs[-2, i]) + log_func(str_funcs[-3, i])) / (log_func(k[-1]) - log_func(k[-2]))**2, (log_func(str_funcs[-1, i]) - 2.0 * log_func(str_funcs[-2, i]) + log_func(str_funcs[-3, i])) / (log_func(k[-1]) - log_func(k[-2]))**2]))
		ax1    = fig.add_subplot(gs[0, 0])
		for i in range(sf.shape[1]):
			x      = log_func(k[:len(sf[:, i])])
			y      = log_func(sf[:, i])
			d_str_func = np.diff(y, n = 1)
			d_k        = np.diff(x, n = 1)
			derivs     = d_str_func / d_k
			# z1     = np.hstack((y[0], y[:-1]))
			# z2     = np.hstack((y[1:], y[-1]))
			# dx1    = np.hstack((0, np.diff(x)))
			# dx2    = np.hstack((np.diff(x), 0))
			# derivs = (z2-z1) / (dx2+dx1)
			ax1.plot(log_func(k[:len(derivs)]), derivs, label=r"$p = {}$".format(i + 1))
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
			ax1.plot(log_func(k[:len(derivs)]), derivs, label=r"$p = {}$".format(i + 1))
		ax1.set_xlabel(r"$k_n$")
		ax1.set_ylabel(r"$\zeta_p$")
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax1.set_title("ESS Local Slope " + m)
		ax1.legend()
		plt.savefig(outdir_path, bbox_inches='tight')
		plt.close()

