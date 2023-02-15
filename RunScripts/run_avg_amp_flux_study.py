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
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
sys.path.append('/home/enda/PhD/Shell_Model/Plotting')
from functions import parse_cml, tc, sim_data, import_data, import_stats_data, import_sys_msr_data, compute_pdf, compute_u_flux, compute_str_func, compute_field_values
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
		
	# -----------------------------------------
	# # --------  Compute the Needed Data
	# -----------------------------------------
	# Get the number of time data points
	num_t_data = sys_vars.ndata

	# Allocate memory
	trip_prod     = np.ones((num_t_data, sys_vars.N)) * 1j
	trip_prod_alt = np.ones((num_t_data, sys_vars.N)) * 1j
	doub_prod     = np.ones((num_t_data, sys_vars.N)) * 1j
	hel_flux      = np.ones((num_t_data, sys_vars.N)) * 1j
	enrg_flux     = np.ones((num_t_data, sys_vars.N)) * 1j

	# Loop through simulation and compute field values
	for t in range(num_t_data):
		# print(t)

		# Get padded field 
		u_pad = np.pad(u[t, :], 2, "constant")

		# Get field values
		trip_prod[t, :], trip_prod_alt[t, :], tmp_dbl_prod, hel_flux[t, :], enrg_flux[t, :] = compute_field_values(u_pad, sys_vars.N, sys_vars.eps, sys_vars.Lambda)
		doub_prod[t, :sys_vars.N - 1] = tmp_dbl_prod

	# -----------------------------------------
	# # --------  Plot Data
	# -----------------------------------------
	# Plotting variables
	fig_size = (10, 6)
	fig_format = ".png"

	# Data
	input_data = [trip_prod, trip_prod_alt, doub_prod, enrg_flux, hel_flux]
	figure_names = ["TripProd", "TripProdAlt", "DoubleProd", "EnrgFlux", "HelFlux"]
	data_labels = [r"u_{n + 2}u_{n + 1}u_{n}", r"(1 - \delta) / \lambda u_{n + 2}u_{n + 1}u_{n}", r"u_{n}u_{n + 3}^{*}", r"\Pi_n^{\mathcal{E}}", r"\Pi_n^{\mathcal{H}}"]

	# Loop through data
	for in_data, fig_name, data_labs in zip(input_data, figure_names, data_labels):
		
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
		fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_Time_Averaged_Spectrum" + fig_format, bbox_inches='tight')
		plt.close()


		#----------------------
		# 2D Histogram
		#----------------------
		num_bins = 500
		fig = plt.figure(figsize=fig_size)
		gs = GridSpec(1, 2, hspace=0.4, wspace=0.5)

		# Real vs Imag 
		ax1 = fig.add_subplot(gs[0, 0])
		hist, xedges, yedges = np.histogram2d(np.real(in_data).flatten(), np.imag(in_data).flatten(), bins=num_bins, density=True)
		ax1.set_xlabel(r"$\Re \left\{" + data_labs + r"\right\}$")
		ax1.set_ylabel(r"$\Im \left\{" + data_labs + r"\right\}$")
		ax1.set_title(r"Real vs Imaginary")
		im1 = ax1.imshow(hist.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		div1  = make_axes_locatable(ax1)
		cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
		cb1   = plt.colorbar(im1, cax = cbax1)
		cb1.set_label("PDF")

		# Abs vs arg 
		ax2 = fig.add_subplot(gs[0, 1])
		hist, xedges, yedges = np.histogram2d(np.absolute(in_data).flatten(), np.mod(np.angle(in_data) + 2.0 * np.pi, 2.0 * np.pi).flatten(), bins=num_bins, density=True)
		ax2.set_xlabel(r"$\left|" + data_labs + r"\right|$")
		ax2.set_ylabel(r"$\arg \left\{ " + data_labs + r" \right\}$")
		ax2.set_title(r"Abs vs Arg")
		im2 = ax2.imshow(np.rot90(hist, k=1), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		ax2.set_ylim(0, 2.0 * np.pi)
		ax2.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		ax2.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		div2  = make_axes_locatable(ax2)
		cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
		cb2   = plt.colorbar(im2, cax = cbax2)
		cb2.set_label("PDF")

		# Save figure
		plt.suptitle(fig_name)
		fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist" + fig_format, bbox_inches='tight')
		plt.close()


		#----------------------------------------
		# 2D Histogram - k_n Pre multiplied Data
		#----------------------------------------
		fig = plt.figure(figsize=fig_size)
		gs = GridSpec(1, 2, hspace=0.4, wspace=0.5)

		# Real vs Imag 
		ax1 = fig.add_subplot(gs[0, 0])
		hist, xedges, yedges = np.histogram2d(np.real(k * in_data).flatten(), np.imag(k * in_data).flatten(), bins=num_bins, density=True)
		ax1.set_xlabel(r"$\Re \left\{ k_n " + data_labs + r"\right\}$")
		ax1.set_ylabel(r"$\Im \left\{ k_n " + data_labs + r"\right\}$")
		ax1.set_title(r"Real vs Imaginary")
		im1 = ax1.imshow(hist.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		div1  = make_axes_locatable(ax1)
		cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
		cb1   = plt.colorbar(im1, cax = cbax1)
		cb1.set_label("PDF")

		# Abs vs arg 
		ax2 = fig.add_subplot(gs[0, 1])
		hist, xedges, yedges = np.histogram2d(np.absolute(k * in_data).flatten(), np.mod(np.angle(k * in_data) + 2.0 * np.pi, 2.0 * np.pi).flatten(), bins=num_bins, density=True)
		ax2.set_ylim(0, 2.0 * np.pi)
		ax2.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		ax2.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		ax2.set_xlabel(r"$\left| k_n " + data_labs + r"\right|$")
		ax2.set_ylabel(r"$\arg \left\{  k_n " + data_labs + r" \right\}$")
		ax2.set_title(r"Abs vs Arg")
		im2 = ax2.imshow(np.rot90(hist.T), extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto", cmap=mpl.colors.ListedColormap(cm.magma.colors), norm=mpl.colors.LogNorm())
		ax2.set_ylim(0, 2.0 * np.pi)
		ax2.set_yticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		ax2.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		div2  = make_axes_locatable(ax2)
		cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
		cb2   = plt.colorbar(im2, cax = cbax2)
		cb2.set_label("PDF")

		# Save figure
		plt.suptitle(fig_name)
		fig.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_2DHist_knPreMult" + fig_format, bbox_inches='tight')
		plt.close()


		#----------------------
		# 1D Histograms
		#----------------------
		num_bins = 250
		norm_hist = False
		shells = [1 - 1, 5 - 1, 20 - 1]
		fig = plt.figure(figsize=(10, 6))
		gs = GridSpec(1, 2, hspace=0.4, wspace=0.3)

		# Real Part
		ax1 = fig.add_subplot(gs[0, 0])
		for j, i in enumerate(shells):
		    pdf, centres = compute_pdf(np.real(in_data[:, i]), nbins=num_bins, normed=norm_hist)
		    if i == -1:
		        p, = ax1.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		    else:
		        p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		ax1.set_xlabel(r"$\Re \left\{" +  data_labs + r" \right\} / \langle (\Re \left\{" +  data_labs + r" \right\})^2 \rangle^{1/2}$")
		ax1.set_ylabel(r"PDF")
		ax1.set_yscale('log')
		ax1.legend()
		ax1.set_title(r"Real Part")
		ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# Imag Part
		ax2 = fig.add_subplot(gs[0, 1])
		for j, i in enumerate(shells):
		    pdf, centres = compute_pdf(np.imag(in_data[:, i]), nbins=num_bins, normed=norm_hist)
		    if i == -1:
		        p, = ax2.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		    else:
		        p, = ax2.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		ax2.set_xlabel(r"$\Im \left\{" +  data_labs + r" \right\} / \langle (\Im \left\{" +  data_labs + r" \right\})^2 \rangle^{1/2}$")
		ax2.set_ylabel(r"PDF")
		ax2.set_yscale('log')
		ax2.legend()
		ax2.set_title(r"Imag Part")
		ax2.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# Save fig
		plt.suptitle(fig_name)
		plt.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DHist_RealImag" + fig_format, bbox_inches='tight')
		plt.close()
		
		num_bins = 250
		norm_hist = True
		shells = [1 - 1, 5 - 1, 20 - 1]
		fig = plt.figure(figsize=(10, 6))
		gs = GridSpec(1, 3, hspace=0.4, wspace=0.3)

		# Abs Part
		ax3 = fig.add_subplot(gs[0, 0])
		for j, i in enumerate(shells):
		    pdf, centres = compute_pdf(np.absolute(in_data[:, i]), nbins=num_bins, normed=norm_hist)
		    if i == -1:
		        p, = ax3.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		    else:
		        p, = ax3.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		ax3.set_xlabel(r"$ \left|" +  data_labs + r" \right| / \langle ( \left|" +  data_labs + r" \right|)^2 \rangle^{1/2}$")
		ax3.set_ylabel(r"PDF")
		ax3.set_yscale('log')
		ax3.legend()
		ax3.set_title(r"Absolute Part")
		ax3.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# Angle Part
		ax4 = fig.add_subplot(gs[0, 1])
		for j, i in enumerate(shells):
		    pdf, centres = compute_pdf(np.mod(np.angle(in_data[:, i]) + 2.0 * np.pi, 2.0*np.pi), nbins=num_bins, normed=False, bin_lims=[0.0, 2.0 * np.pi])
		    if i == -1:
		        p, = ax4.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		    else:
		        p, = ax4.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		ax4.set_xlim(0, 2.0*np.pi)
		ax4.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
		ax4.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
		ax4.set_xlabel(r"$\arg \left\{" +  data_labs + r" \right\} $")
		ax4.set_ylabel(r"PDF")
		ax4.set_yscale('log')
		ax4.legend()
		ax4.set_title(r"Phase Part")
		ax4.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# Sin (Angle) Part
		ax4 = fig.add_subplot(gs[0, 2])
		for j, i in enumerate(shells):
		    pdf, centres = compute_pdf(np.sin(np.mod(np.angle(in_data[:, i]) + 2.0 * np.pi, 2.0*np.pi)), nbins=num_bins, normed=False, bin_lims=[-1, 1])
		    if i == -1:
		        p, = ax4.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
		    else:
		        p, = ax4.plot(centres, pdf, label="$n = {}$".format(i + 1))    
		ax4.set_xlabel(r"$\sin \arg \left\{" +  data_labs + r" \right\} $")
		ax4.set_ylabel(r"PDF")
		ax4.set_yscale('log')
		ax4.legend()
		ax4.set_title(r"sin(Phase) Part")
		ax4.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)

		# Save fig
		plt.suptitle(fig_name)
		plt.savefig(cmdargs.out_dir_AVGFLUX + fig_name + "_1DHist_AbsAngle" + fig_format, bbox_inches='tight')
		plt.close()



fig = plt.figure(figsize=(10, 6))
gs = GridSpec(1, 1, hspace=0.4, wspace=0.3)
ax1 = fig.add_subplot(gs[0, 0])

# Plot the time averaged tripple product
data = np.absolute(np.mean(np.imag(trip_prod), axis=0))
xlab = r"$\log k_n$"
ylab = r"$\log \left|\Im \left\{" + data_labels[0] + r"\right\}\right|$"
plot_spectrum(fig, ax1, data, k, xlab, ylab)
data = np.mean(np.absolute(trip_prod), axis=0) * np.absolute(np.mean(np.sin(np.angle(trip_prod)), axis=0))
xlab = r"$\log k_n$"
ylab = r"$\log \langle \left|" + data_labels[0] + r"\right| \rangle \langle \sin \arg\left\{" + data_labels[0] + r"\right\}\rangle$"
plot_spectrum(fig, ax1, data, k, xlab, ylab)
plt.suptitle(figure_names[0])
plt.savefig(cmdargs.out_dir_AVGFLUX + figure_names[0] + "_Test_TimeAverage_Indepence" + fig_format, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=(10, 6))
gs = GridSpec(1, 1, hspace=0.4, wspace=0.3)
ax1 = fig.add_subplot(gs[0, 0])

# Plot the time averaged tripple product
data = np.mean(np.imag(enrg_flux), axis=0)
ax1.plot(data, label=r"Full")
data = np.mean(np.absolute(enrg_flux), axis=0) * np.mean(np.sin(np.angle(enrg_flux)), axis=0)
ax1.plot(data, label=r"Prod")
plt.suptitle(figure_names[-2])
plt.savefig(cmdargs.out_dir_AVGFLUX + figure_names[-2] + "_Test_TimeAverage_Indepence" + fig_format, bbox_inches='tight')
plt.close()




fig = plt.figure(figsize=(10, 6))
gs = GridSpec(1, 1, hspace=0.4, wspace=0.3)
ax1 = fig.add_subplot(gs[0, 0])
pdf, centres = compute_pdf(np.mod(np.angle(enrg_flux[:, :].flatten()) + 2.0 * np.pi, 2.0*np.pi), nbins=num_bins, normed=False, bin_lims=[0.0 - 2.0 * np.pi / num_bins, 2.0 * np.pi + 2.0 * np.pi/num_bins])
p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))
# ax1.set_xlim(0, 2.0*np.pi)
ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
ax1.set_xlabel(r"$\arg \left\{" + data_labels[0] + r" \right\} $")
ax1.set_ylabel(r"PDF")
ax1.set_yscale('log')
ax1.legend()
ax1.set_title(r"Phase Part")
ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
plt.savefig(cmdargs.out_dir_AVGFLUX + figure_names[0] + "_Test_PDF_TriadsAll" + fig_format, bbox_inches='tight')
plt.close()