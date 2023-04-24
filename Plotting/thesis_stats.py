import matplotlib as mpl
# mpl.use('TkAgg') # Use this backend for displaying plots in window
mpl.use('Agg') # Use this backend for writing plots to file
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import h5py
import sys
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import cm 
import numpy as np
import itertools
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc, rcParams
from functions import tc, sim_data, import_data, import_stats_data, import_sys_msr_data, import_phase_sync_data, compute_pdf_from_hist, compute_pdf, slope_fit
from plot_functions import plot_anomalous_exponent, plot_str_funcs_with_slope, phase_only_space_time
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['font.serif']  = 'Computer Modern Roman'
rcParams['font.size']   = 12.5 # def: 10.0
## Lines
rcParams['lines.linewidth']  = 1.5 # def: 1.5
rcParams['lines.markersize'] = 5.0 # def: 6
rcParams['lines.markersize'] = 5.0 # def: 6
## Grid lines
rcParams['grid.color']     = 'k'
rcParams['grid.linestyle'] = ':'
rcParams['grid.linewidth'] = 0.5
rcParams['grid.alpha']     = 0.8
##Ticks
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top']       = True
rcParams['ytick.right']     = True
## Figsize
textwidth = 12.25
rcParams['figure.figsize'] = [textwidth, 0.4 * textwidth]

inset_lab_size = 10
label_size = 16

plot_colours = plt.cm.YlGnBu(np.linspace(0.3, 1.0, 10))

fig_size_nx1 = (0.7 * textwidth, 1.2 * textwidth)
fig_size_2x2 = (textwidth, 0.8 * textwidth)

my_magma = matplotlib.colors.ListedColormap(cm.magma.colors[::-1])
my_magma.set_under(color = "white")
my_cmap = matplotlib.colors.ListedColormap(cm.YlGnBu(np.arange(0,cm.YlGnBu.N)))
my_cmap.set_under(color = "white")


fig_format = 'pdf'

data_dir = str(sys.argv[1])

## ------ Input and Output Dirs
input_dir  = data_dir 
output_dir = "/home/ecarroll/PhD/Shell_Model/Data/Thesis/" + "Plots/"
data_file_path = input_dir

# -----------------------------------------
# # --------  Read In data
# -----------------------------------------
## Read in simulation parameters
sys_vars = sim_data(data_file_path, "default")

## Read in solver data
run_data = import_data(data_file_path, sys_vars, "default")

## Read in stats data
stats_data = import_stats_data(data_file_path, sys_vars, "default")

## Read in sys_msr data
sys_msr_data = import_sys_msr_data(data_file_path, sys_vars, "default")


# -----------------------------------------
# # --------  Plot Flux Scaling
# -----------------------------------------
fig = plt.figure()
gs  = GridSpec(1, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(sys_msr_data.k / sys_vars.k0, np.absolute(sys_msr_data.enrg_flux_t_avg), color = plot_colours[-3])
ax1.plot(sys_msr_data.k[2:18] / sys_vars.k0, sys_msr_data.k[2:18] ** (-1), '--', color=plot_colours[-1], label=r"$k_n^{-1}$")
ax1.set_xlabel(r"$k_n / k_0$", fontsize=label_size)
ax1.set_ylabel(r"$\left\langle\Pi^{\mathcal{E}}_n\right\rangle$", fontsize=label_size)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()
ax1.grid()
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(sys_msr_data.k / sys_vars.k0, np.absolute(sys_msr_data.kin_hel_flux_t_avg), color = plot_colours[-3])
ax2.plot(sys_msr_data.k[2:11] / sys_vars.k0,  sys_msr_data.k[2:11] ** (-2), '--', color=plot_colours[-1], label=r"$k_n^{-2}$")
ax2.set_xlabel(r"$k_n / k_0$", fontsize=label_size)
ax2.set_ylabel(r"$\left\langle\Pi^{\mathcal{H}}_n\right\rangle$", fontsize=label_size)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()
ax2.grid()
plt.savefig(output_dir + "Shell_FluxScaling" + "." + fig_format, format = fig_format, bbox_inches='tight', dpi=1200)
plt.close()


# -----------------------------------------
# # --------  Plot Shell Tseries
# -----------------------------------------
short_file = "/home/ecarroll/PhD/Shell_Model/Data/Thesis/HD-INTFACRK4-FULL_N[25]_T[0.0,0.0001,100.000]_NU[5e-07]_ALPHA[1.500]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[N_SCALING]_TAG[Short]"
with h5py.File(short_file + "/Main_HDF_Data.h5", 'r') as infile:
	short_u = infile["VelModes"][:, :]
with h5py.File(short_file + "/System_Measure_HDF_Data.h5", 'r') as infile:
	short_time = infile["Time"][:]

fig = plt.figure()
gs  = GridSpec(3, 1)
ax1 = fig.add_subplot(gs[0, 0])
t1 = 5
data = np.real(short_u[:, t1 - 1])
ax1.plot(short_time[:], data / np.sqrt(np.mean(data**2)), label=r"$n = {}$".format(t1), color=plot_colours[-3])
# ax1.set_ylabel(r"$\Re \left\{u_n \right\}$", fontsize = label_size)
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_xlim(short_time[0], short_time[-1])
ax1.grid()
ax1.legend(loc='upper right', prop={'size':10})
ax2 = fig.add_subplot(gs[1, 0])
t1 = 10
data = np.real(short_u[:, t1 - 1])
ax2.plot(short_time[:], data / np.sqrt(np.mean(data**2)), label=r"$n = {}$".format(t1), color=plot_colours[-3])
ax2.set_ylabel(r"$\Re \left\{u_n \right\}$", fontsize = label_size)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_xlim(short_time[0], short_time[-1])
ax2.grid()
ax2.legend(loc='upper right', prop={'size':10})
ax3 = fig.add_subplot(gs[2, 0])
t1 = 15
data = np.real(short_u[:, t1 - 1])
ax3.plot(short_time[:], data / np.sqrt(np.mean(data**2)), label=r"$n = {}$".format(t1), color=plot_colours[-3])
ax3.set_xlabel(r"$t$", fontsize = label_size)
# ax3.set_ylabel(r"$\Re \left\{u_n \right\}$", fontsize = label_size)
ax3.set_yticklabels([])
ax3.set_xlim(short_time[0], short_time[-1])
ax3.grid()
ax3.legend(loc='upper right', prop={'size':10})
plt.savefig(output_dir + "Shell_Tsereis" + "." + fig_format, format = fig_format, bbox_inches='tight', dpi=1200)
plt.close()


# -----------------------------------------
# # --------  Plot PDFs
# -----------------------------------------
nbins = 100
fig = plt.figure()
gs  = GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0, 0])
indxs = [5, 10, 15, 20]
for j, i in enumerate(indxs):
	if hasattr(stats_data, "vel_hist_counts"):
		pdf, centres = compute_pdf_from_hist(stats_data.vel_hist_counts[i, :], stats_data.vel_hist_ranges[i, :], normed = True)
	else:
		if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
			pdf, centres = compute_pdf(np.real(run_data.a_n[:, i] * np.exp(1j * run_data.phi_n[:, i])), nbins = nbins, normed = True)					
		else:
			pdf, centres = compute_pdf(np.real(run_data.u[:, i]), nbins = nbins, normed = True)    
	ax1.plot(centres, pdf, label = "$n = {}$".format(i), color=plot_colours[(j + 5 - len(indxs)) * 2])    
ax1.set_xlabel(r"$\Re \left\{u_n\right\} / \langle \Re\left\{ u_n \right\}^2 \rangle^{1/2}$")
ax1.set_ylabel(r"PDF")
ax1.grid()
ax1.set_yscale('log')
ax1.legend()
plt.savefig(output_dir + "Shell_PDFs" + "." + fig_format, format = fig_format, bbox_inches='tight', dpi=1200)
plt.close()
