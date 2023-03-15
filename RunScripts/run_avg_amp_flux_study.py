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
from functions import parse_cml, tc, sim_data, import_data, import_stats_data, import_sys_msr_data, compute_pdf, compute_u_flux, compute_str_func, compute_field_values, slope_fit
from plot_functions import plot_spectrum
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
	cmdargs.out_dir_AVGFLUX = cmdargs.out_dir + "AVGFLUX_PLOTS/"
	if os.path.isdir(cmdargs.out_dir_AVGFLUX) != True:
	    print("Making folder:" + tc.C + " AVGFLUX_PLOTS/" + tc.Rst)
	    os.mkdir(cmdargs.out_dir_AVGFLUX)

	# -----------------------------------------
	# # --------  Read In data
	# -----------------------------------------
	# Read in simulation parameters
	sys_vars = sim_data(data_file_path, method)

	# # Read in solver data
	# run_data = import_data(data_file_path, sys_vars, method)

	# # Read in stats data
	# stats_data = import_stats_data(data_file_path, sys_vars, method)

	# # Read in sys_msr data
	# sys_msr_data = import_sys_msr_data(data_file_path, sys_vars, method)

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

	# Allocate memory
	trip_prod_alt = np.ones((num_t_data, sys_vars.N)) * 1j
	trip_prod     = np.ones((num_t_data, sys_vars.N)) * 1j
	doub_prod     = np.ones((num_t_data, sys_vars.N)) * 1j
	hel_flux      = np.ones((num_t_data, sys_vars.N)) * 1j
	enrg_flux     = np.ones((num_t_data, sys_vars.N)) * 1j

	avg_trip_prod     = np.ones((num_t_data, sys_vars.N)) * 1j
	avg_doub_prod     = np.ones((num_t_data, sys_vars.N)) * 1j
	avg_hel_flux      = np.ones((num_t_data, sys_vars.N)) * 1j
	avg_enrg_flux     = np.ones((num_t_data, sys_vars.N)) * 1j

	# Loop through simulation and compute field values
	for t in range(num_t_data):
		# print(t)

		# Get padded field 
		u_pad = np.pad(u[t, :], 2, "constant")

		# Get field values
		trip_prod[t, :], trip_prod_alt[t, :], tmp_dbl_prod, hel_flux[t, :], enrg_flux[t, :] = compute_field_values(u_pad, sys_vars.N, sys_vars.eps, sys_vars.Lambda)
		doub_prod[t, :sys_vars.N - 1] = tmp_dbl_prod

		# Get padded field with time averaged amplitudes 
		u_pad = np.pad(u[t, :] * (a_n_t_avg / np.absolute(u[t, :])), 2, "constant")

		# Get field values
		avg_trip_prod[t, :], _, tmp_dbl_prod, avg_hel_flux[t, :], avg_enrg_flux[t, :] = compute_field_values(u_pad, sys_vars.N, sys_vars.eps, sys_vars.Lambda)
		avg_doub_prod[t, :sys_vars.N - 1] = tmp_dbl_prod

	# -----------------------------------------
	# # --------  Plot Data
	# -----------------------------------------
	# Plotting variables
	fig_size = (10, 6)
	fig_format = ".png"

	# Data
	input_data = [doub_prod, trip_prod, u, avg_trip_prod] #, k * trip_prod, trip_prod_alt, doub_prod, enrg_flux, hel_flux] #, avg_trip_prod, k * avg_trip_prod, avg_doub_prod, avg_enrg_flux, avg_hel_flux]
	figure_names = ["DoubleProd", "TripProd", "u", "TimeAvgTripProd"] #, "kTripProd", "TripProdAlt", "DoubleProd", "EnrgFlux", "HelFlux"] #, "TimeAvgTripProd", "TimeAvgkTripProd", "TimeAvgDoubleProd", "TimeAvgEnrgFlux", "TimeAvgHelFlux"]
	data_labels = [r"u_{n}u_{n + 3}^{*}", r"u_{n + 2}u_{n + 1}u_{n}", r"u_{n}", r"u_{n + 2}u_{n + 1}u_{n}"] #, r"k_n u_{n + 2}u_{n + 1}u_{n}", r"(1 - \delta) / \lambda u_{n + 2}u_{n + 1}u_{n}", r"u_{n}u_{n + 3}^{*}", r"\Pi_n^{\mathcal{E}}", r"\Pi_n^{\mathcal{H}}"] #, r"u_{n + 2}u_{n + 1}u_{n}", r"k_n u_{n + 2}u_{n + 1}u_{n}", r"u_{n}u_{n + 3}^{*}", r"\Pi_n^{\mathcal{E}}", r"\Pi_n^{\mathcal{H}}"]

	# Loop through data
	for in_data, fig_name, data_labs in zip(input_data, figure_names, data_labels):

		cmdargs.out_dir_AVGFLUX = cmdargs.out_dir + "AVGFLUX_PLOTS/" + fig_name + "/"
		if os.path.isdir(cmdargs.out_dir_AVGFLUX) != True:
		    print("Making folder:" + tc.C + " AVGFLUX_PLOTS/" + tc.Rst)
		    os.mkdir(cmdargs.out_dir_AVGFLUX)

		#----------------------
		# 1D PDF
		#----------------------
		## 1D - PDF - Angle
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
		fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1D_PDF_Angles" + fig_format, bbox_inches='tight')
		plt.close()

		num_bins = 100
		norm_hist = False
		fig = plt.figure(figsize=(24, 24))
		gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		for i in range(5):
			for j in range(5):
				indx = i * 5 + j
				if indx < in_data.shape[-1] - 1:
					ax1 = fig.add_subplot(gs[i, j])
					p, = ax1.plot(centres, pdf_angles[indx] - 2.0 * pdf_angles[indx + 1], label="$n = {}$".format(indx + 1))    
					ax1.set_xlabel(r"$ \arg \left\{" +  data_labs + r" \right\}$")
					ax1.set_xlim(0, 2.0 * np.pi)
					ax1.set_ylabel(r"PDF")
					# ax1.set_yscale('log')
					ax1.set_title("n = {}".format(indx + 1))
		fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1D_PDF_Angles_Diff" + fig_format, bbox_inches='tight')
		plt.close()

		# ## 1D - PDF - Amp
		# fig = plt.figure(figsize=(24, 24))
		# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		# for i in range(5):
		# 	for j in range(5):
		# 		indx = i * 5 + j
		# 		if indx < in_data.shape[-1]:
		# 			ax1 = fig.add_subplot(gs[i, j])
		# 			pdf, centres = compute_pdf(np.absolute(in_data[:, indx]), nbins=num_bins, normed=norm_hist)
		# 			p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
		# 			ax1.set_xlabel(r"$ \left|" +  data_labs + r" \right|$")
		# 			ax1.set_ylabel(r"PDF")
		# 			ax1.set_yscale('log')
		# 			ax1.set_title("n = {}".format(indx + 1))
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1D_PDF_Amp" + fig_format, bbox_inches='tight')
		# plt.close()

		# ## 1D - PDF - Imag
		# num_bins = 100
		# norm_hist = False
		# fig = plt.figure(figsize=(24, 24))
		# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		# for i in range(5):
		# 	for j in range(5):
		# 		indx = i * 5 + j
		# 		if indx < in_data.shape[-1]:
		# 			ax1 = fig.add_subplot(gs[i, j])
		# 			pdf, centres = compute_pdf(np.imag(in_data[:, indx]), nbins=num_bins, normed=norm_hist)
		# 			p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
		# 			ax1.set_xlabel(r"$ \Im \left\{" +  data_labs + r" \right\}$")
		# 			ax1.set_ylabel(r"PDF")
		# 			ax1.set_yscale('log')
		# 			ax1.set_title("n = {}".format(indx + 1))
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1D_PDF_Imag" + fig_format, bbox_inches='tight')
		# plt.close()

		# ## 1D - PDF - Real
		# fig = plt.figure(figsize=(24, 24))
		# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		# for i in range(5):
		# 	for j in range(5):
		# 		indx = i * 5 + j
		# 		if indx < in_data.shape[-1]:
		# 			ax1 = fig.add_subplot(gs[i, j])
		# 			pdf, centres = compute_pdf(np.real(in_data[:, indx]), nbins=num_bins, normed=norm_hist)
		# 			p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
		# 			ax1.set_xlabel(r"$ \Re \left\{" +  data_labs + r" \right\}$")
		# 			ax1.set_ylabel(r"PDF")
		# 			ax1.set_yscale('log')
		# 			ax1.set_title("n = {}".format(indx + 1))
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1D_PDF_Real" + fig_format, bbox_inches='tight')
		# plt.close()

		# #----------------------
		# # 1D Tseries
		# #----------------------
		# ## 1D Tseries - Amps
		# fig = plt.figure(figsize=(24, 24))
		# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		# for i in range(5):
		# 	for j in range(5):
		# 		indx = i * 5 + j
		# 		if indx < in_data.shape[-1]:
		# 			ax1 = fig.add_subplot(gs[i, j])
		# 			ax1.plot(np.absolute(in_data[:, indx]))
		# 			ax1.set_xlabel(r"Time")
		# 			ax1.set_ylabel(r"Amp")
		# 			ax1.set_title("n = {}".format(indx + 1))
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DTSeries_Abs" + fig_format, bbox_inches='tight')
		# plt.close()

		# ## 1D Tseries - Phases
		# fig = plt.figure(figsize=(24, 24))
		# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		# for i in range(5):
		# 	for j in range(5):
		# 		indx = i * 5 + j
		# 		if indx < in_data.shape[-1]:
		# 			ax1 = fig.add_subplot(gs[i, j])
		# 			ax1.plot(np.unwrap(np.angle(in_data[:, indx])))
		# 			ax1.set_xlabel(r"Time")
		# 			ax1.set_ylabel(r"Arg")
		# 			ax1.set_title("n = {}".format(indx + 1))
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DTSeries_Arg" + fig_format, bbox_inches='tight')
		# plt.close()

		# ## 1D Tseries - Re
		# fig = plt.figure(figsize=(24, 24))
		# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		# for i in range(5):
		# 	for j in range(5):
		# 		indx = i * 5 + j
		# 		if indx < in_data.shape[-1]:
		# 			ax1 = fig.add_subplot(gs[i, j])
		# 			ax1.plot(np.real(in_data[:, indx]))
		# 			ax1.set_xlabel(r"Time")
		# 			ax1.set_ylabel(r"$\Re$")
		# 			ax1.set_title("n = {}".format(indx + 1))
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DTSeries_Real" + fig_format, bbox_inches='tight')
		# plt.close()

		# ## 1D Tseries - Im
		# fig = plt.figure(figsize=(24, 24))
		# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
		# for i in range(5):
		# 	for j in range(5):
		# 		indx = i * 5 + j
		# 		if indx < in_data.shape[-1]:
		# 			ax1 = fig.add_subplot(gs[i, j])
		# 			ax1.plot(np.imag(in_data[:, indx]))
		# 			ax1.set_xlabel(r"Time")
		# 			ax1.set_ylabel(r"$\Im$")
		# 			ax1.set_title("n = {}".format(indx + 1))
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DTSeries_Imag" + fig_format, bbox_inches='tight')
		# plt.close()

		# #----------------------
		# # Time Averaged Spectra
		# #----------------------
		# fig = plt.figure(figsize=fig_size)
		# gs = GridSpec(2, 2, hspace=0.4, wspace=0.3)
		# ax1 = fig.add_subplot(gs[0, 0])
		# ax2 = fig.add_subplot(gs[0, 1])
		# ax3 = fig.add_subplot(gs[1, 0])
		# ax4 = fig.add_subplot(gs[1, 1])

		# # Plot the time averaged tripple product
		# data = np.absolute(np.mean(np.imag(in_data), axis=0))
		# xlab = r"$\log k_n$"
		# ylab = r"$\log \left|\Im \left\{" + data_labs + r"\right\}\right|$"
		# plot_spectrum(fig, ax1, data, k, xlab, ylab)
		# ax1.set_title("Time Averaged Imaginary part")

		# # Plot the time averaged tripple product time k_n
		# data = np.absolute(np.mean(np.imag(in_data * k), axis=0))
		# ylab = r"$\log \left|k_n \Im \left\{" + data_labs + r"\right\}\right|$"
		# plot_spectrum(fig, ax2, data, k, xlab, ylab)
		# ax2.set_title("Time Averaged Imaginary part time $k_n$")

		# # Plot the time averaged tripple product
		# data = np.absolute(np.mean(np.real(in_data), axis=0))
		# xlab = r"$\log k_n$"
		# ylab = r"$\log \left|\Re \left\{" + data_labs + r"\right\}\right|$"
		# plot_spectrum(fig, ax3, data, k, xlab, ylab)
		# ax3.set_title("Time Averaged real part")

		# # Plot the time averaged tripple product time k_n
		# data = np.absolute(np.mean(np.real(in_data *  k), axis=0))
		# ylab = r"$\log \left|k_n \Re \left\{" + data_labs + r"\right\}\right|$"
		# plot_spectrum(fig, ax4, data, k, xlab, ylab)
		# ax4.set_title("Time Averaged Real part time $k_n$")

		# # Save figure
		# plt.suptitle(fig_name)
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_Time_Averaged_Spectrum" + fig_format, bbox_inches='tight')
		# plt.close()

		# fig = plt.figure(figsize=fig_size)
		# gs = GridSpec(1, 1, hspace=0.4, wspace=0.5)
		# ax1 = fig.add_subplot(gs[0, 0])
		# ax1.scatter(np.absolute(in_data).flatten(), np.sin(np.angle(in_data).flatten()))
		# ax2.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# ax2.set_ylabel(r"$\sin \arg \left\{ " + data_labs + r" \right\}$")
		# plt.suptitle(fig_name)
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_Scatter" + fig_format, bbox_inches='tight')
		# plt.close()

		# #----------------------
		# # 2D Histogram
		# #----------------------
		# num_bins = 500
		# fig = plt.figure(figsize=fig_size)
		# gs = GridSpec(1, 2, hspace=0.4, wspace=0.5)

		# # Real vs Imag 
		# ax1 = fig.add_subplot(gs[0, 0])
		# hist, xedges, yedges = np.histogram2d(np.real(in_data).flatten(), np.imag(in_data).flatten(), bins=num_bins, density=True)
		# ax1.set_xlabel(r"$\Re \left\{" + data_labs + r"\right\}$")
		# ax1.set_ylabel(r"$\Im \left\{" + data_labs + r"\right\}$")
		# ax1.set_title(r"Real vs Imaginary")
		# im1 = ax1.imshow(np.rot90(hist), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		# div1  = make_axes_locatable(ax1)
		# ax1.set_xlim(-10, 10)
		# ax1.set_ylim(-10, 10)
		# cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
		# cb1   = plt.colorbar(im1, cax = cbax1)
		# cb1.set_label("PDF")

		# # Abs vs arg 
		# ax2 = fig.add_subplot(gs[0, 1])
		# hist, xedges, yedges = np.histogram2d(np.absolute(in_data).flatten(), np.mod(np.angle(in_data) + 2.0 * np.pi, 2.0 * np.pi).flatten(), bins=num_bins, density=True)
		# ax2.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# ax2.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# ax2.set_title(r"Abs vs Arg")
		# im2 = ax2.imshow(np.rot90(hist, k=1), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		# ax2.set_ylim(0, 2.0 * np.pi)
		# ax2.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# ax2.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# div2  = make_axes_locatable(ax2)
		# cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
		# cb2   = plt.colorbar(im2, cax = cbax2)
		# cb2.set_label("PDF")

		# # Save figure
		# plt.suptitle(fig_name)
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist" + fig_format, bbox_inches='tight')
		# plt.close()

		# #----------------------------------------
		# # 2D Histogram - k_n Pre multiplied Data
		# #----------------------------------------
		# fig = plt.figure(figsize=fig_size)
		# gs = GridSpec(1, 2, hspace=0.4, wspace=0.5)

		# # Real vs Imag 
		# ax1 = fig.add_subplot(gs[0, 0])
		# hist, xedges, yedges = np.histogram2d(np.real(k * in_data).flatten(), np.imag(k * in_data).flatten(), bins=num_bins, density=True)
		# ax1.set_xlabel(r"$\Re \left\{ k_n " + data_labs + r"\right\}$")
		# ax1.set_ylabel(r"$\Im \left\{ k_n " + data_labs + r"\right\}$")
		# ax1.set_title(r"Real vs Imaginary")
		# im1 = ax1.imshow(hist.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		# div1  = make_axes_locatable(ax1)
		# cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
		# cb1   = plt.colorbar(im1, cax = cbax1)
		# cb1.set_label("PDF")

		# # Abs vs arg 
		# ax2 = fig.add_subplot(gs[0, 1])
		# hist, xedges, yedges = np.histogram2d(np.absolute(k * in_data).flatten(), np.mod(np.angle(k * in_data) + 2.0 * np.pi, 2.0 * np.pi).flatten(), bins=num_bins, density=True)
		# ax2.set_ylim(0, 2.0 * np.pi)
		# ax2.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# ax2.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# ax2.set_xlabel(r"$\left| k_n " + data_labs + r"\right|$")
		# ax2.set_ylabel(r"$\arg \left\{  k_n " + data_labs + r" \right\}$")
		# ax2.set_title(r"Abs vs Arg")
		# im2 = ax2.imshow(np.rot90(hist.T), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		# ax2.set_ylim(0, 2.0 * np.pi)
		# ax2.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# ax2.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# div2  = make_axes_locatable(ax2)
		# cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
		# cb2   = plt.colorbar(im2, cax = cbax2)
		# cb2.set_label("PDF")

		# # Save figure
		# plt.suptitle(fig_name)
		# fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist_knPreMult" + fig_format, bbox_inches='tight')
		# plt.close()



		# #----------------------
		# # 1D Histograms Im and Re
		# #----------------------
		# num_bins = 500
		# norm_hist = True
		# shells = [1 - 1, 5 - 1, 15 - 1]

		# fig = plt.figure(figsize=(10, 6))
		# gs = GridSpec(1, 2, hspace=0.4, wspace=0.3)

		# # Real Part
		# ax1 = fig.add_subplot(gs[0, 0])
		# for j, i in enumerate(shells):
		#     pdf, centres = compute_pdf(np.real(in_data[:, i]), nbins=num_bins, normed=norm_hist)
		#     if i == -1:
		#         p, = ax1.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		#     else:
		#         p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		# ax1.set_xlabel(r"$\Re \left\{" +  data_labs + r" \right\} / \langle (\Re \left\{" +  data_labs + r" \right\})^2 \rangle^{1/2}$")
		# ax1.set_ylabel(r"PDF")
		# ax1.set_yscale('log')
		# ax1.legend()
		# ax1.set_title(r"Real Part")
		# ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# # Imag Part
		# ax2 = fig.add_subplot(gs[0, 1])
		# for j, i in enumerate(shells):
		#     pdf, centres = compute_pdf(np.imag(in_data[:, i]), nbins=num_bins, normed=norm_hist)
		#     if i == -1:
		#         p, = ax2.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		#     else:
		#         p, = ax2.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		# ax2.set_xlabel(r"$\Im \left\{" +  data_labs + r" \right\} / \langle (\Im \left\{" +  data_labs + r" \right\})^2 \rangle^{1/2}$")
		# ax2.set_ylabel(r"PDF")
		# ax2.set_yscale('log')
		# ax2.legend()
		# ax2.set_title(r"Imag Part")
		# ax2.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# # Save fig
		# plt.suptitle(fig_name)
		# plt.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DHist_RealImag" + fig_format, bbox_inches='tight')
		# plt.close()


		# #----------------------
		# # 1D Histograms Abs, Arg and sin(Arg)
		# #----------------------
		# num_bins = 500
		# norm_hist = True
		# shells = [1 - 1, 5 - 1, 15 - 1]
		# fig = plt.figure(figsize=(10, 6))
		# gs = GridSpec(1, 3, hspace=0.4, wspace=0.3)

		# # Abs Part
		# ax3 = fig.add_subplot(gs[0, 0])
		# for j, i in enumerate(shells):
		#     pdf, centres = compute_pdf(np.absolute(in_data[:, i]), nbins=num_bins, normed=norm_hist)
		#     if i == -1:
		#         p, = ax3.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		#     else:
		#         p, = ax3.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		# ax3.set_xlabel(r"$ \left|" +  data_labs + r" \right| / \langle ( \left|" +  data_labs + r" \right|)^2 \rangle^{1/2}$")
		# ax3.set_ylabel(r"PDF")
		# ax3.set_yscale('log')
		# ax3.legend()
		# ax3.set_title(r"Absolute Part")
		# ax3.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# # Angle Part
		# ax4 = fig.add_subplot(gs[0, 1])
		# for j, i in enumerate(shells):
		#     pdf, centres = compute_pdf(np.mod(np.angle(in_data[:, i]) + 2.0 * np.pi, 2.0*np.pi), nbins=num_bins, normed=False, bin_lims=[0.0, 2.0 * np.pi])
		#     if i == -1:
		#         p, = ax4.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		#     else:
		#         p, = ax4.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		# ax4.set_xlim(0, 2.0*np.pi)
		# ax4.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# ax4.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# ax4.set_xlabel(r"$\arg \left\{" +  data_labs + r" \right\} $")
		# ax4.set_ylabel(r"PDF")
		# ax4.set_yscale('log')
		# ax4.legend()
		# ax4.set_title(r"Phase Part")
		# ax4.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# # Sin (Angle) Part
		# ax4 = fig.add_subplot(gs[0, 2])
		# for j, i in enumerate(shells):
		#     pdf, centres = compute_pdf(np.sin(np.angle(in_data[:, i])), nbins=num_bins, normed=False)
		#     if i == -1:
		#         p, = ax4.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		#     else:
		#         p, = ax4.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		# ax4.set_xlabel(r"$\sin \arg \left\{" +  data_labs + r" \right\} $")
		# ax4.set_ylabel(r"PDF")
		# ax4.set_yscale('log')
		# ax4.legend()
		# ax4.set_title(r"sin(Phase) Part")
		# ax4.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# # Save fig
		# plt.suptitle(fig_name)
		# plt.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DHist_AbsAngle" + fig_format, bbox_inches='tight')
		# plt.close()

		# #----------------------
		# # PDF of Arg
		# #----------------------
		# bin_limits = [0.0 - 2.0 * np.pi / num_bins, 2.0 * np.pi + 2.0 * np.pi/num_bins]
		# fig = plt.figure(figsize=(10, 6))
		# gs  = GridSpec(1, 1, hspace=0.4, wspace=0.3)
		# ax1 = fig.add_subplot(gs[0, 0])
		# pdf, centres = compute_pdf(np.mod(np.angle(in_data[:, :].flatten()) + 2.0 * np.pi, 2.0*np.pi), nbins=num_bins, normed=False, bin_lims=bin_limits)
		# p,           = ax1.plot(centres, pdf)
		# # ax1.set_xlim(0, 2.0*np.pi)
		# ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# ax1.set_xlabel(r"$\arg \left\{" + data_labs + r" \right\} $")
		# ax1.set_ylabel(r"PDF")
		# ax1.set_yscale('log')
		# ax1.set_title(r"Phase Part")
		# ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
		# plt.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_Test_PDF_PhasesAll" + fig_format, bbox_inches='tight')
		# plt.close()

		# #----------------------
		# # Time Average of Abs and sin(Arg)
		# #----------------------
		# fig = plt.figure(figsize=(15, 6))
		# gs  = GridSpec(1, 3, hspace=0.4, wspace=0.3)
		# ax1 = fig.add_subplot(gs[0, 0])
		# ax1.plot(k * np.mean(np.absolute(in_data), axis=0), label=r"$Abs$")
		# ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
		# ax1.set_yscale('log')
		# ax1.legend()
		# ax1 = fig.add_subplot(gs[0, 1])
		# ax1.plot(np.mean(np.sin(np.angle(in_data)), axis=0), label=r"$\sin \arg $")
		# ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
		# # ax1.set_yscale('log')
		# ax1.legend()
		# ax1 = fig.add_subplot(gs[0, 2])
		# ax1.plot(np.mean(np.angle(in_data), axis=0), label=r"$\arg $")
		# ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
		# # ax1.set_yscale('log')
		# ax1.legend()
		# plt.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_Test_Average" + fig_format, bbox_inches='tight')
		# plt.close()

		#----------------------
		# Test SF Independence
		#----------------------
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
		plt.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_Test_Independence_SF" + fig_format, bbox_inches='tight')
		plt.close()


		# #----------------------
		# # Independence of Abs and Angle
		# #----------------------
		# num_bins = 250

		# shells = np.arange(in_data.shape[-1])
		# # shells = [1 - 1, 5 - 1, 10 - 1, 15 - 1, 20 - 1]

		# density_flag = True
		# aspect_flag = "auto"
		# c_map = mpl.colors.ListedColormap(cm.magma.colors)
		# c_map_norm = mpl.colors.LogNorm()

		# # Norms
		# diff_sqr_error_msr = np.zeros((len(shells), ))
		# hellinger_msr      = np.zeros((len(shells), ))
		# q_2norm_msr        = np.zeros((len(shells), ))
		# q_1norm_msr        = np.zeros((len(shells), ))
		# kl_div_msr         = np.zeros((len(shells), ))
		# for s, n in enumerate(shells):
		# 	# if not np.all(np.imag(in_data[:, n])):
		# 	# 	print(np.all(np.imag(in_data[:, n])), np.imag(in_data[:, n]))
		# 	if np.all(np.imag(in_data[:, n])):
		# 		## Comparison of 2D Distributions
		# 		fig = plt.figure(figsize=(16, 22))
		# 		gs = GridSpec(3, 2, hspace=0.4, wspace=0.5)

		# 		print(fig_name, n)

		# 		# Abs vs arg 2D (joint) distribution
		# 		ax1 = fig.add_subplot(gs[0, 0])
		# 		x = np.absolute(in_data[:, n])
		# 		y = np.mod(np.angle(in_data[:, n]) + 2.0 * np.pi, 2.0 * np.pi)
		# 		hist, xedges, yedges = np.histogram2d(x, y, bins=(np.linspace(x.min(), x.max(), num_bins + 1), np.linspace(0.0, 2.0 * np.pi, num_bins + 1)), density=density_flag)
		# 		ax1.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax1.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		ax1.set_title(r"Abs and Arg")
		# 		im1 = ax1.imshow(np.rot90(hist, k=1), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax1.set_ylim(0, 2.0 * np.pi)
		# 		ax1.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# 		ax1.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# 		div1  = make_axes_locatable(ax1)
		# 		margaxr = div1.append_axes("right", size = "15%", pad = 0.05)
		# 		pdf, centres = compute_pdf(np.mod(np.angle(in_data[:, n]) + 2.0 * np.pi, 2.0*np.pi), nbins=num_bins, normed=False, bin_lims=bin_limits)
		# 		p,           = margaxr.plot(pdf, centres, label="$n = {}$".format(n + 1))
		# 		margaxr.set_xscale('log')
		# 		# margaxr.set_ylim(centres[0], centres[-1])
		# 		margaxr.set_xticks([])
		# 		margaxr.set_xticklabels([])
		# 		margaxr.set_yticks([])
		# 		margaxr.set_yticklabels([])
		# 		margaxr.spines['bottom'].set_visible(False)
		# 		margaxr.spines['top'].set_visible(False)
		# 		margaxr.spines['right'].set_visible(False)
		# 		cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
		# 		cb2   = plt.colorbar(im1, cax = cbax1)
		# 		cb2.set_label("PDF")
		# 		margaxt = div1.append_axes("top", size = "15%", pad = 0.05)
		# 		pdf, centres = compute_pdf(np.absolute(in_data[:, n]), nbins=num_bins, normed=False)
		# 		p,           = margaxt.plot(centres, pdf, label="$n = {}$".format(n + 1))
		# 		margaxt.set_yscale('log')
		# 		# margaxt.set_xlim(centres[0], centres[-1])
		# 		margaxt.set_xticks([])
		# 		margaxt.set_xticklabels([])
		# 		margaxt.set_yticks([])
		# 		margaxt.set_yticklabels([])
		# 		margaxt.spines['left'].set_visible(False)
		# 		margaxt.spines['top'].set_visible(False)
		# 		margaxt.spines['right'].set_visible(False)
				
		# 		# Product of the marginal 1d distributions
		# 		ax2 = fig.add_subplot(gs[0, 1])
		# 		ax2.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax2.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		ax2.set_title(r"Marginal Abs * Marginal Arg")
		# 		marg_abs_pdf = np.sum(hist, axis=0)	* (yedges[1] - yedges[0])
		# 		marg_angle_pdf = np.sum(hist, axis=1) * (xedges[1] - xedges[0])
		# 		pdf_prod_data = np.outer(marg_angle_pdf, marg_abs_pdf)
		# 		im2 = ax2.imshow(np.rot90(pdf_prod_data), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax2.set_ylim(0, 2.0 * np.pi)
		# 		ax2.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# 		ax2.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# 		div2  = make_axes_locatable(ax2)
		# 		cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
		# 		cb2   = plt.colorbar(im2, cax = cbax2)
		# 		cb2.set_label("PDF")

		# 		# Error between the two
		# 		dx = xedges[1] - xedges[0]
		# 		dy = yedges[1] - yedges[0]
		# 		hist_2d = np.rot90(hist, k=1) 
		# 		marg_data = np.rot90(pdf_prod_data)

		# 		# Compute 2d distances
		# 		nx, ny = hist_2d.shape
		# 		diff_sqr_error_dist = np.zeros((nx, ny))
		# 		hellinger_dist      = np.zeros((nx, ny))
		# 		q_dist              = np.zeros((nx, ny))
		# 		kl_div_dist         = np.zeros((nx, ny))
		# 		for i in range(nx):
		# 			for j in range(ny):
		# 				diff_sqr_error_dist[i, j] = (hist_2d[i, j] - marg_data[i, j])**2

		# 				hellinger_dist[i, j]      = (np.sqrt(hist_2d[i, j]) - np.sqrt(marg_data[i, j]))**2

		# 				if (hist_2d[i, j] + marg_data[i, j]) == 0.0:
		# 					q_dist[i, j] = 0.0
		# 				else:
		# 					q_dist[i, j] = (hist_2d[i, j] - marg_data[i, j]) / (hist_2d[i, j] + marg_data[i, j])

		# 				if marg_data[i, j] == 0.0 or hist_2d[i, j] / marg_data[i, j] == 0.0:
		# 					kl_div_dist[i, j] = 0.0
		# 				else:
		# 					kl_div_dist[i, j] = hist_2d[i, j] * np.log(hist_2d[i, j] / marg_data[i, j])

		# 		# Compute measures
		# 		diff_sqr_error_msr[s] = np.sum(diff_sqr_error_dist) * dy * dx
		# 		hellinger_msr[s]      = np.sqrt(0.5 * np.sum(hellinger_dist)* dx * dy)
		# 		q_2norm_msr[s]        = np.linalg.norm(hist_2d - marg_data) / (np.linalg.norm(hist_2d) + np.linalg.norm(marg_data))
		# 		q_1norm_msr[s]        = np.linalg.norm(hist_2d - marg_data, ord=1) / (np.linalg.norm(hist_2d, ord=1) + np.linalg.norm(marg_data, ord=1))
		# 		kl_div_msr[s]         = np.sum(kl_div_dist) * dx * dy


		# 		ax3 = fig.add_subplot(gs[1, 0])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(hellinger_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax3.set_ylim(0, 2.0 * np.pi)
		# 		ax3.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# 		ax3.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"Helling Distance: {:1.6f}".format(hellinger_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"Hellinger Distance")

		# 		ax3 = fig.add_subplot(gs[1, 1])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(q_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax3.set_ylim(0, 2.0 * np.pi)
		# 		ax3.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# 		ax3.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"Q: 1norm {:1.6f} - 2norm {:1.6f}".format(q_1norm_msr[s], q_2norm_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"Q")

		# 		ax3 = fig.add_subplot(gs[2, 0])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(kl_div_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax3.set_ylim(0, 2.0 * np.pi)
		# 		ax3.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# 		ax3.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"KL Divergence: {:1.6f}".format(kl_div_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"KL Divergence")

		# 		ax3 = fig.add_subplot(gs[2, 1])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(diff_sqr_error_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax3.set_ylim(0, 2.0 * np.pi)
		# 		ax3.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		# 		ax3.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"Diff Sqrd Error: {:1.6f}".format(diff_sqr_error_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"Diff Sqrd Error")

		# 		print("Amp - Phase -- Diff Sqr Error:\t{:1.6f}".format(diff_sqr_error_msr[s]))
		# 		print("Amp - Phase -- Hellinger Measure:\t{:1.6f}".format(hellinger_msr[s]))	
		# 		print("Amp - Phase -- Qoutient Error 1norm:\t{:1.6f}".format(q_1norm_msr[s]))	
		# 		print("Amp - Phase -- Qoutient Error 2norm:\t{:1.6f}".format(q_2norm_msr[s]))	
		# 		print("Amp - Phase -- KL Divergence:\t{:1.6f}".format(kl_div_msr[s]))

		# 		# Save figure
		# 		plt.suptitle(fig_name + " $n = {} $".format(n + 1))
		# 		fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist_Comparison_n{}".format(n) + fig_format, bbox_inches='tight')
		# 		plt.close()


		# 		## Comparison of 2D Distributions
		# 		fig = plt.figure(figsize=(16, 22))
		# 		gs = GridSpec(3, 2, hspace=0.4, wspace=0.5)

		# 		print(fig_name, n)

		# 		# Abs vs arg 2D (joint) distribution
		# 		ax1 = fig.add_subplot(gs[0, 0])
		# 		x = np.absolute(in_data[:, n])
		# 		y = np.sin(np.mod(np.angle(in_data[:, n]) + 2.0 * np.pi, 2.0 * np.pi))
		# 		hist, xedges, yedges = np.histogram2d(x, y, bins=(np.linspace(x.min(), x.max(), num_bins + 1), np.linspace(-1.0, 1.0, num_bins + 1)), density=density_flag)
		# 		ax1.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax1.set_ylabel(r"$\sin \arg \left\{ " + data_labs + r" \right\}$")
		# 		ax1.set_title(r"Abs and sinArg")
		# 		im1 = ax1.imshow(np.rot90(hist, k=1), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax1.set_ylim(-1.0, 1.0)
		# 		div1  = make_axes_locatable(ax1)
		# 		margaxr = div1.append_axes("right", size = "15%", pad = 0.05)
		# 		pdf, centres = compute_pdf(np.sin(np.angle(in_data[:, n])), nbins=num_bins, normed=False, bin_lims=[-1.0, 1.0])
		# 		p,           = margaxr.plot(pdf, centres, label="$n = {}$".format(n + 1))
		# 		margaxr.set_xscale('log')
		# 		margaxr.set_ylim(centres[0], centres[-1])
		# 		margaxr.set_xticks([])
		# 		margaxr.set_xticklabels([])
		# 		margaxr.set_yticks([])
		# 		margaxr.set_yticklabels([])
		# 		margaxr.spines['bottom'].set_visible(False)
		# 		margaxr.spines['top'].set_visible(False)
		# 		margaxr.spines['right'].set_visible(False)
		# 		cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
		# 		cb2   = plt.colorbar(im1, cax = cbax1)
		# 		cb2.set_label("PDF")
		# 		margaxt = div1.append_axes("top", size = "15%", pad = 0.05)
		# 		pdf, centres = compute_pdf(np.absolute(in_data[:, n]), nbins=num_bins, normed=False)
		# 		p,           = margaxt.plot(centres, pdf, label="$n = {}$".format(n + 1))
		# 		margaxt.set_yscale('log')
		# 		margaxt.set_xlim(centres[0], centres[-1])
		# 		margaxt.set_xticks([])
		# 		margaxt.set_xticklabels([])
		# 		margaxt.set_yticks([])
		# 		margaxt.set_yticklabels([])
		# 		margaxt.spines['left'].set_visible(False)
		# 		margaxt.spines['top'].set_visible(False)
		# 		margaxt.spines['right'].set_visible(False)
				
		# 		# Product of the marginal 1d distributions
		# 		ax2 = fig.add_subplot(gs[0, 1])
		# 		ax2.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax2.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		ax2.set_title(r"Abs * Arg")
		# 		marg_abs_pdf = np.sum(hist, axis=0)	* (yedges[1] - yedges[0])
		# 		marg_angle_pdf = np.sum(hist, axis=1) * (xedges[1] - xedges[0])
		# 		pdf_prod_data = np.outer(marg_angle_pdf, marg_abs_pdf)
		# 		im2 = ax2.imshow(np.flipud(pdf_prod_data), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		ax2.set_ylim(-1.0, 1.0)
		# 		div2  = make_axes_locatable(ax2)
		# 		cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
		# 		cb2   = plt.colorbar(im2, cax = cbax2)
		# 		cb2.set_label("PDF")

		# 		# Error between the two
		# 		dx = xedges[1] - xedges[0]
		# 		dy = yedges[1] - yedges[0]
		# 		hist_2d = np.rot90(hist, k=1) 
		# 		marg_data = np.rot90(pdf_prod_data)

		# 		# Compute 2d distances
		# 		nx, ny = hist_2d.shape
		# 		diff_sqr_error_dist = np.zeros((nx, ny))
		# 		hellinger_dist      = np.zeros((nx, ny))
		# 		q_dist              = np.zeros((nx, ny))
		# 		kl_div_dist         = np.zeros((nx, ny))
		# 		for i in range(nx):
		# 			for j in range(ny):
		# 				diff_sqr_error_dist[i, j] = (hist_2d[i, j] - marg_data[i, j])**2

		# 				hellinger_dist[i, j]      = (np.sqrt(hist_2d[i, j]) - np.sqrt(marg_data[i, j]))**2

		# 				if (hist_2d[i, j] + marg_data[i, j]) == 0.0:
		# 					q_dist[i, j] = 0.0
		# 				else:
		# 					q_dist[i, j] = (hist_2d[i, j] - marg_data[i, j]) / (hist_2d[i, j] + marg_data[i, j])

		# 				if marg_data[i, j] == 0.0 or hist_2d[i, j] / marg_data[i, j] == 0.0:
		# 					kl_div_dist[i, j] = 0.0
		# 				else:
		# 					kl_div_dist[i, j] = hist_2d[i, j] * np.log(hist_2d[i, j] / marg_data[i, j])

		# 		# Compute measures
		# 		diff_sqr_error_msr[s] = np.sum(diff_sqr_error_dist) * dy * dx
		# 		hellinger_msr[s]      = np.sqrt(0.5 * np.sum(hellinger_dist)* dx * dy)
		# 		q_2norm_msr[s]        = np.linalg.norm(hist_2d - marg_data) / (np.linalg.norm(hist_2d) + np.linalg.norm(marg_data))
		# 		q_1norm_msr[s]        = np.linalg.norm(hist_2d - marg_data, ord=1) / (np.linalg.norm(hist_2d, ord=1) + np.linalg.norm(marg_data, ord=1))
		# 		kl_div_msr[s]         = np.sum(kl_div_dist) * dx * dy

		# 		ax3 = fig.add_subplot(gs[1, 0])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(hellinger_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"Helling Distance: {:1.6f}".format(hellinger_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"Hellinger Distance")

		# 		ax3 = fig.add_subplot(gs[1, 1])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(q_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"Q: 1norm {:1.6f} - 2norm {:1.6f}".format(q_1norm_msr[s], q_2norm_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"Q")

		# 		ax3 = fig.add_subplot(gs[2, 0])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(kl_div_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"KL Divergence: {:1.6f}".format(kl_div_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"KL Divergence")

		# 		ax3 = fig.add_subplot(gs[2, 1])
		# 		ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 		ax3.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		# 		im3 = ax3.imshow(diff_sqr_error_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 		div3  = make_axes_locatable(ax3)
		# 		ax3.set_title(r"Diff Sqrd Error: {:1.6f}".format(diff_sqr_error_msr[s]))
		# 		cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 		cb3   = plt.colorbar(im3, cax = cbax3)
		# 		cb3.set_label(r"Diff Sqrd Error")

		# 		print("Amp - Sin(Phase) --- Diff Sqr Error:\t{:1.6f}".format(diff_sqr_error_msr[s]))
		# 		print("Amp - Sin(Phase) --- Hellinger Measure:\t{:1.6f}".format(hellinger_msr[s]))	
		# 		print("Amp - Sin(Phase) --- Qoutient Error 1norm:\t{:1.6f}".format(q_1norm_msr[s]))	
		# 		print("Amp - Sin(Phase) --- Qoutient Error 2norm:\t{:1.6f}".format(q_2norm_msr[s]))	
		# 		print("Amp - Sin(Phase) --- KL Divergence:\t{:1.6f}".format(kl_div_msr[s]))

		# 		# Save figure
		# 		plt.suptitle(fig_name + " $n = {} $".format(n + 1))
		# 		fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist_ComparisonSinArg_n{}".format(n) + fig_format, bbox_inches='tight')
		# 		plt.close()


		# if fig_name == "u":
		# 	num_bins   = 250
		# 	separation = 4
		# 	shells     = np.arange(in_data.shape[-1] - separation)
		# 	# shells = [1 - 1, 5 - 1, 10 - 1, 15 - 1, 20 - 1]

		# 	density_flag = True
		# 	aspect_flag = "auto"
		# 	c_map = mpl.colors.ListedColormap(cm.magma.colors)
		# 	c_map_norm = mpl.colors.LogNorm()

		# 	# Norms
		# 	diff_sqr_error_msr = np.zeros((len(shells), ))
		# 	hellinger_msr      = np.zeros((len(shells), ))
		# 	q_2norm_msr        = np.zeros((len(shells), ))
		# 	q_1norm_msr        = np.zeros((len(shells), ))
		# 	kl_div_msr         = np.zeros((len(shells), ))
		# 	for s, n in enumerate(shells):
		# 		# if not np.all(np.imag(in_data[:, n])):
		# 		# 	print(np.all(np.imag(in_data[:, n])), np.imag(in_data[:, n]))
		# 		if np.all(np.imag(in_data[:, n])):
		# 			## Comparison of 2D Distributions
		# 			fig = plt.figure(figsize=(16, 22))
		# 			gs = GridSpec(3, 2, hspace=0.4, wspace=0.5)

		# 			print(fig_name, n)

		# 			# Abs vs arg 2D (joint) distribution
		# 			ax1 = fig.add_subplot(gs[0, 0])
		# 			x = np.absolute(in_data[:, n])
		# 			y = np.absolute(in_data[:, n + 2])
		# 			hist, xedges, yedges = np.histogram2d(x, y, bins=(np.linspace(x.min(), x.max(), num_bins + 1), np.linspace(y.min(), y.max(), num_bins + 1)), density=density_flag)
		# 			ax1.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax1.set_ylabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax1.set_title(r"Abs $u_{n}$ and Abs $u_{n + 1}$")
		# 			im1 = ax1.imshow(np.rot90(hist, k=1), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 			div1  = make_axes_locatable(ax1)
		# 			margaxr = div1.append_axes("right", size = "15%", pad = 0.05)
		# 			pdf, centres = compute_pdf(y, nbins=num_bins, normed=False)
		# 			p,           = margaxr.plot(pdf, centres, label="$n = {}$".format(n + 1))
		# 			margaxr.set_xscale('log')
		# 			# margaxr.set_ylim(centres[0], centres[-1])
		# 			margaxr.set_xticks([])
		# 			margaxr.set_xticklabels([])
		# 			margaxr.set_yticks([])
		# 			margaxr.set_yticklabels([])
		# 			margaxr.spines['bottom'].set_visible(False)
		# 			margaxr.spines['top'].set_visible(False)
		# 			margaxr.spines['right'].set_visible(False)
		# 			cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
		# 			cb2   = plt.colorbar(im1, cax = cbax1)
		# 			cb2.set_label("PDF")
		# 			margaxt = div1.append_axes("top", size = "15%", pad = 0.05)
		# 			pdf, centres = compute_pdf(x, nbins=num_bins, normed=False)
		# 			p,           = margaxt.plot(centres, pdf, label="$n = {}$".format(n + 1))
		# 			margaxt.set_yscale('log')
		# 			margaxt.set_xticks([])
		# 			margaxt.set_xticklabels([])
		# 			margaxt.set_yticks([])
		# 			margaxt.set_yticklabels([])
		# 			margaxt.spines['left'].set_visible(False)
		# 			margaxt.spines['top'].set_visible(False)
		# 			margaxt.spines['right'].set_visible(False)
					
		# 			# Product of the marginal 1d distributions
		# 			ax2 = fig.add_subplot(gs[0, 1])
		# 			ax1.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax1.set_ylabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax2.set_title(r"Marginal Abs $u_n$ * Marginal Abs $u_{n + 1}$")
		# 			marg_abs_pdf = np.sum(hist, axis=0)	* (yedges[1] - yedges[0])
		# 			marg_angle_pdf = np.sum(hist, axis=1) * (xedges[1] - xedges[0])
		# 			pdf_prod_data = np.outer(marg_angle_pdf, marg_abs_pdf)
		# 			im2 = ax2.imshow(np.rot90(pdf_prod_data), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 			div2  = make_axes_locatable(ax2)
		# 			cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
		# 			cb2   = plt.colorbar(im2, cax = cbax2)
		# 			cb2.set_label("PDF")

		# 			# Error between the two
		# 			dx = xedges[1] - xedges[0]
		# 			dy = yedges[1] - yedges[0]
		# 			hist_2d = np.rot90(hist, k=1) 
		# 			marg_data = np.rot90(pdf_prod_data)

		# 			# Compute 2d distances
		# 			nx, ny = hist_2d.shape
		# 			diff_sqr_error_dist = np.zeros((nx, ny))
		# 			hellinger_dist      = np.zeros((nx, ny))
		# 			q_dist              = np.zeros((nx, ny))
		# 			kl_div_dist         = np.zeros((nx, ny))
		# 			for i in range(nx):
		# 				for j in range(ny):
		# 					diff_sqr_error_dist[i, j] = (hist_2d[i, j] - marg_data[i, j])**2

		# 					hellinger_dist[i, j]      = (np.sqrt(hist_2d[i, j]) - np.sqrt(marg_data[i, j]))**2

		# 					if (hist_2d[i, j] + marg_data[i, j]) == 0.0:
		# 						q_dist[i, j] = 0.0
		# 					else:
		# 						q_dist[i, j] = (hist_2d[i, j] - marg_data[i, j]) / (hist_2d[i, j] + marg_data[i, j])

		# 					if marg_data[i, j] == 0.0 or hist_2d[i, j] / marg_data[i, j] == 0.0:
		# 						kl_div_dist[i, j] = 0.0
		# 					else:
		# 						kl_div_dist[i, j] = hist_2d[i, j] * np.log(hist_2d[i, j] / marg_data[i, j])

		# 			# Compute measures
		# 			diff_sqr_error_msr[s] = np.sum(diff_sqr_error_dist) * dy * dx
		# 			hellinger_msr[s]      = np.sqrt(0.5 * np.sum(hellinger_dist)* dx * dy)
		# 			q_2norm_msr[s]        = np.linalg.norm(hist_2d - marg_data) / (np.linalg.norm(hist_2d) + np.linalg.norm(marg_data))
		# 			q_1norm_msr[s]        = np.linalg.norm(hist_2d - marg_data, ord=1) / (np.linalg.norm(hist_2d, ord=1) + np.linalg.norm(marg_data, ord=1))
		# 			kl_div_msr[s]         = np.sum(kl_div_dist) * dx * dy


		# 			ax3 = fig.add_subplot(gs[1, 0])
		# 			ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
		# 			im3 = ax3.imshow(hellinger_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 			div3  = make_axes_locatable(ax3)
		# 			ax3.set_title(r"Helling Distance: {:1.6f}".format(hellinger_msr[s]))
		# 			cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 			cb3   = plt.colorbar(im3, cax = cbax3)
		# 			cb3.set_label(r"Hellinger Distance")

		# 			ax3 = fig.add_subplot(gs[1, 1])
		# 			ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
		# 			im3 = ax3.imshow(q_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 			div3  = make_axes_locatable(ax3)
		# 			ax3.set_title(r"Q: 1norm {:1.6f} - 2norm {:1.6f}".format(q_1norm_msr[s], q_2norm_msr[s]))
		# 			cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 			cb3   = plt.colorbar(im3, cax = cbax3)
		# 			cb3.set_label(r"Q")

		# 			ax3 = fig.add_subplot(gs[2, 0])
		# 			ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
		# 			im3 = ax3.imshow(kl_div_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 			div3  = make_axes_locatable(ax3)
		# 			ax3.set_title(r"KL Divergence: {:1.6f}".format(kl_div_msr[s]))
		# 			cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 			cb3   = plt.colorbar(im3, cax = cbax3)
		# 			cb3.set_label(r"KL Divergence")

		# 			ax3 = fig.add_subplot(gs[2, 1])
		# 			ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		# 			ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
		# 			im3 = ax3.imshow(diff_sqr_error_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
		# 			div3  = make_axes_locatable(ax3)
		# 			ax3.set_title(r"Diff Sqrd Error: {:1.6f}".format(diff_sqr_error_msr[s]))
		# 			cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
		# 			cb3   = plt.colorbar(im3, cax = cbax3)
		# 			cb3.set_label(r"Diff Sqrd Error")

		# 			print("Amp - Phase -- Diff Sqr Error:\t{:1.6f}".format(diff_sqr_error_msr[s]))
		# 			print("Amp - Phase -- Hellinger Measure:\t{:1.6f}".format(hellinger_msr[s]))	
		# 			print("Amp - Phase -- Qoutient Error 1norm:\t{:1.6f}".format(q_1norm_msr[s]))	
		# 			print("Amp - Phase -- Qoutient Error 2norm:\t{:1.6f}".format(q_2norm_msr[s]))	
		# 			print("Amp - Phase -- KL Divergence:\t{:1.6f}".format(kl_div_msr[s]))

		# 			# Save figure
		# 			plt.suptitle(fig_name + " $n = {} $".format(n + 1))
		# 			fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist_Comparison_Un{}_n{}".format(separation, n) + fig_format, bbox_inches='tight')
		# 			plt.close()


		if fig_name == "TripProd":
			for sep in range(2, 4):
				separation = sep
				shells     = np.arange(in_data.shape[-1] - separation)
				num_bins   = 250
				# shells = [1 - 1, 5 - 1, 10 - 1, 15 - 1, 20 - 1]

				density_flag = True
				aspect_flag = "auto"
				c_map = mpl.colors.ListedColormap(cm.magma.colors)
				c_map_norm = mpl.colors.LogNorm()

				# Norms
				diff_sqr_error_msr = np.zeros((len(shells), ))
				hellinger_msr      = np.zeros((len(shells), ))
				q_2norm_msr        = np.zeros((len(shells), ))
				q_1norm_msr        = np.zeros((len(shells), ))
				kl_div_msr         = np.zeros((len(shells), ))
				for s, n in enumerate(shells):
					# if not np.all(np.imag(in_data[:, n])):
					# 	print(np.all(np.imag(in_data[:, n])), np.imag(in_data[:, n]))
					if np.all(np.imag(in_data[:, n])):
						## Comparison of 2D Distributions
						fig = plt.figure(figsize=(16, 22))
						gs = GridSpec(3, 2, hspace=0.4, wspace=0.5)

						print(fig_name, n)

						# Abs vs arg 2D (joint) distribution
						ax1 = fig.add_subplot(gs[0, 0])
						x = np.angle(in_data[:, n])
						y = np.angle(in_data[:, n + separation])
						hist, xedges, yedges = np.histogram2d(x, y, bins=(np.linspace(x.min(), x.max(), num_bins + 1), np.linspace(y.min(), y.max(), num_bins + 1)), density=density_flag)
						ax1.set_xlabel(r"$\left|" + data_labs + r"\right|$")
						ax1.set_ylabel(r"$\left|" + data_labs + r"\right|$")
						ax1.set_title(r"Abs $u_{n}$ and Abs $u_{n + 1}$")
						im1 = ax1.imshow(np.rot90(hist, k=1), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
						div1  = make_axes_locatable(ax1)
						margaxr = div1.append_axes("right", size = "15%", pad = 0.05)
						pdf, centres = compute_pdf(y, nbins=num_bins, normed=False)
						p,           = margaxr.plot(pdf, centres, label="$n = {}$".format(n + 1))
						margaxr.set_xscale('log')
						# margaxr.set_ylim(centres[0], centres[-1])
						margaxr.set_xticks([])
						margaxr.set_xticklabels([])
						margaxr.set_yticks([])
						margaxr.set_yticklabels([])
						margaxr.spines['bottom'].set_visible(False)
						margaxr.spines['top'].set_visible(False)
						margaxr.spines['right'].set_visible(False)
						cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
						cb2   = plt.colorbar(im1, cax = cbax1)
						cb2.set_label("PDF")
						margaxt = div1.append_axes("top", size = "15%", pad = 0.05)
						pdf, centres = compute_pdf(x, nbins=num_bins, normed=False)
						p,           = margaxt.plot(centres, pdf, label="$n = {}$".format(n + 1))
						margaxt.set_yscale('log')
						margaxt.set_xticks([])
						margaxt.set_xticklabels([])
						margaxt.set_yticks([])
						margaxt.set_yticklabels([])
						margaxt.spines['left'].set_visible(False)
						margaxt.spines['top'].set_visible(False)
						margaxt.spines['right'].set_visible(False)
						
						# Product of the marginal 1d distributions
						ax2 = fig.add_subplot(gs[0, 1])
						ax1.set_xlabel(r"$\left|" + data_labs + r"\right|$")
						ax1.set_ylabel(r"$\left|" + data_labs + r"\right|$")
						ax2.set_title(r"Marginal Abs $u_n$ * Marginal Abs $u_{n + 1}$")
						marg_abs_pdf = np.sum(hist, axis=0)	* (yedges[1] - yedges[0])
						marg_angle_pdf = np.sum(hist, axis=1) * (xedges[1] - xedges[0])
						pdf_prod_data = np.outer(marg_angle_pdf, marg_abs_pdf)
						im2 = ax2.imshow(np.rot90(pdf_prod_data), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
						div2  = make_axes_locatable(ax2)
						cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
						cb2   = plt.colorbar(im2, cax = cbax2)
						cb2.set_label("PDF")

						# Error between the two
						dx = xedges[1] - xedges[0]
						dy = yedges[1] - yedges[0]
						hist_2d = np.rot90(hist, k=1) 
						marg_data = np.rot90(pdf_prod_data)

						# Compute 2d distances
						nx, ny = hist_2d.shape
						diff_sqr_error_dist = np.zeros((nx, ny))
						hellinger_dist      = np.zeros((nx, ny))
						q_dist              = np.zeros((nx, ny))
						kl_div_dist         = np.zeros((nx, ny))
						for i in range(nx):
							for j in range(ny):
								diff_sqr_error_dist[i, j] = (hist_2d[i, j] - marg_data[i, j])**2

								hellinger_dist[i, j]      = (np.sqrt(hist_2d[i, j]) - np.sqrt(marg_data[i, j]))**2

								if (hist_2d[i, j] + marg_data[i, j]) == 0.0:
									q_dist[i, j] = 0.0
								else:
									q_dist[i, j] = (hist_2d[i, j] - marg_data[i, j]) / (hist_2d[i, j] + marg_data[i, j])

								if marg_data[i, j] == 0.0 or hist_2d[i, j] / marg_data[i, j] == 0.0:
									kl_div_dist[i, j] = 0.0
								else:
									kl_div_dist[i, j] = hist_2d[i, j] * np.log(hist_2d[i, j] / marg_data[i, j])

						# Compute measures
						diff_sqr_error_msr[s] = np.sum(diff_sqr_error_dist) * dy * dx
						hellinger_msr[s]      = np.sqrt(0.5 * np.sum(hellinger_dist)* dx * dy)
						q_2norm_msr[s]        = np.linalg.norm(hist_2d - marg_data) / (np.linalg.norm(hist_2d) + np.linalg.norm(marg_data))
						q_1norm_msr[s]        = np.linalg.norm(hist_2d - marg_data, ord=1) / (np.linalg.norm(hist_2d, ord=1) + np.linalg.norm(marg_data, ord=1))
						kl_div_msr[s]         = np.sum(kl_div_dist) * dx * dy


						ax3 = fig.add_subplot(gs[1, 0])
						ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
						ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
						im3 = ax3.imshow(hellinger_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
						div3  = make_axes_locatable(ax3)
						ax3.set_title(r"Helling Distance: {:1.6f}".format(hellinger_msr[s]))
						cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
						cb3   = plt.colorbar(im3, cax = cbax3)
						cb3.set_label(r"Hellinger Distance")

						ax3 = fig.add_subplot(gs[1, 1])
						ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
						ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
						im3 = ax3.imshow(q_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
						div3  = make_axes_locatable(ax3)
						ax3.set_title(r"Q: 1norm {:1.6f} - 2norm {:1.6f}".format(q_1norm_msr[s], q_2norm_msr[s]))
						cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
						cb3   = plt.colorbar(im3, cax = cbax3)
						cb3.set_label(r"Q")

						ax3 = fig.add_subplot(gs[2, 0])
						ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
						ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
						im3 = ax3.imshow(kl_div_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
						div3  = make_axes_locatable(ax3)
						ax3.set_title(r"KL Divergence: {:1.6f}".format(kl_div_msr[s]))
						cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
						cb3   = plt.colorbar(im3, cax = cbax3)
						cb3.set_label(r"KL Divergence")

						ax3 = fig.add_subplot(gs[2, 1])
						ax3.set_xlabel(r"$\left|" + data_labs + r"\right|$")
						ax3.set_ylabel(r"$\left|" + data_labs + r"\right|$")
						im3 = ax3.imshow(diff_sqr_error_dist, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect=aspect_flag, cmap=c_map, norm=c_map_norm)
						div3  = make_axes_locatable(ax3)
						ax3.set_title(r"Diff Sqrd Error: {:1.6f}".format(diff_sqr_error_msr[s]))
						cbax3 = div3.append_axes("right", size = "5%", pad = 0.05)
						cb3   = plt.colorbar(im3, cax = cbax3)
						cb3.set_label(r"Diff Sqrd Error")

						print("Amp - Phase -- Diff Sqr Error:\t{:1.6f}".format(diff_sqr_error_msr[s]))
						print("Amp - Phase -- Hellinger Measure:\t{:1.6f}".format(hellinger_msr[s]))	
						print("Amp - Phase -- Qoutient Error 1norm:\t{:1.6f}".format(q_1norm_msr[s]))	
						print("Amp - Phase -- Qoutient Error 2norm:\t{:1.6f}".format(q_2norm_msr[s]))	
						print("Amp - Phase -- KL Divergence:\t{:1.6f}".format(kl_div_msr[s]))

						# Save figure
						plt.suptitle(fig_name + " $n = {} $".format(n + 1))
						fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist_Comparison_Angle(Triad)n{}_n{}".format(separation, n) + fig_format, bbox_inches='tight')
						plt.close()



