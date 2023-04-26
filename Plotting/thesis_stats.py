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
from numba import njit
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc, rcParams
from functions import tc, sim_data, import_data, import_stats_data, import_sys_msr_data, import_phase_sync_data, compute_pdf_from_hist, compute_pdf, slope_fit, compute_u_flux
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


@njit
def get_fluxes(u, delta, l):

	num_t, num_shell = u.shape

	enrg_flux = np.zeros(u.shape)
	hel_flux  = np.zeros(u.shape)
	for t in range(num_t):
		enrg_flux[t, :], hel_flux[t, :] = compute_u_flux(u[t, :], num_shell, delta, l)

	return enrg_flux, hel_flux


def plot_sf(fig, ax, p_range, k, str_funcs, inert_range, ylab, insert_fig = True, scaling = 'loge', ax_pos=[0.15, 0.15, 0.15, 0.15]):

	if insert_fig:
		axin = fig.add_axes(ax_pos)

	ns_zeta_p    = [0.7, 1, 1.27, 1.53, 1.78]
	zeta_p       = []
	zeta_p_resid = []

	inert_lim_low  = inert_range[0]
	inert_lim_high = inert_range[-1]
	print("Inertial Range Shell No.s: {} = {}".format(inert_lim_low + 1, inert_lim_high + 1))

	if scaling == 'log2':
		log_func = np.log2
		xlabel_Str = r"$\log_2 k_n/k_0 $"
		ylabel_Str = r"$\log_2 " + ylab + " $"
	elif scaling == 'log10':
		log_func = np.log10
		xlabel_Str = r"$\log_10 k_n/k_0 $"
		ylabel_Str = r"$\log_10 " + ylab + " $"
	else:
		log_func = np.log
		xlabel_Str = r"$\ln k_n/k_0$"
		ylabel_Str = r"$\ln " + ylab + " $"

	## Loop over structure functions and plot
	print("Power\tp/3\t\tDNS Slope\tNS")
	for i in p_range:
	    
		## Plot strucure function
		p, = ax.plot(log_func(k[:str_funcs.shape[0]]), log_func(str_funcs[:, i]), label = "$p = {}$".format(int(i)), color=plot_colours[(i - 2) * 2])

		## Find polynomial fit and plot
		poly_output = np.polyfit(log_func(k[inert_lim_low:inert_lim_high]), log_func(str_funcs[inert_lim_low:inert_lim_high, i]), 1, full = True)
		pfit_info   = poly_output[0]
		poly_resid  = poly_output[1][0]
		pfit_slope  = pfit_info[0]
		pfit_c      = pfit_info[1]
		zeta_p.append(np.absolute(pfit_slope))
		zeta_p_resid.append(poly_resid)
		ax.plot(log_func(k[inert_lim_low:inert_lim_high]), log_func(k[inert_lim_low:inert_lim_high])*pfit_slope + pfit_c, '--', color = p.get_color())

		if i >= 1 and i <= len(ns_zeta_p):
			print(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f} \t{:1.3f}".format(i, -(i) / 3, pfit_slope, poly_resid, -ns_zeta_p[i - 2]))
		else:
			print(" {}\t {:1.4f} \t {:1.4f} +/-{:0.3f}".format(i, -(i) / 3, pfit_slope, poly_resid))

		if insert_fig:
			# ## Compute the local derivative and plot in insert
			# d_str_func  = np.diff(log_func(str_funcs[:, i]), n = 2)
			# d_k         = np.diff(log_func(k), n = 2)
			# print(d_str_func)
			# print(d_k)
			# local_deriv = d_str_func / d_k
			# print(local_deriv)
			# local_deriv = np.concatenate((local_deriv, [(log_func(str_funcs[-1, i]) - 2.0 * log_func(str_funcs[-2, i]) + log_func(str_funcs[-3, i])) / (log_func(k[-1]) - log_func(k[-2]))**2, (log_func(str_funcs[-1, i]) - 2.0 * log_func(str_funcs[-2, i]) + log_func(str_funcs[-3, i])) / (log_func(k[-1]) - log_func(k[-2]))**2]))
			# print()
			local_deriv = np.gradient(log_func(str_funcs[:, i]))
			axin.plot(log_func(k[:len(local_deriv)]), local_deriv, color = p.get_color())
			axin.set_ylabel(r"$\zeta_p$", labelpad = -40)
			axin.set_xlabel(xlabel_Str, labelpad = -30)
			axin.set_ylim(bottom=-5)
			axin.tick_params(axis='both', which='major', labelsize=inset_lab_size)

	## Set axis labels 
	ax.set_xlabel(xlabel_Str)
	ax.set_ylabel(ylabel_Str)
	ax.grid()
	axin.grid()
	ax.legend(loc='upper right', prop={'size':10})

	return zeta_p, zeta_p_resid


def plot_anom_scaling(fig, ax, p, zeta_p, ns_zeta_p, label_str):
	ax.plot(p, zeta_p, marker = 'o', markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = label_str, color=plot_colours[-1])
	ax.plot(np.arange(2, 6 + 1), ns_zeta_p, marker = '.', markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Navier Stokes", color=plot_colours[-3])
	ax.plot(p, p / 3, 'k--', label = r"$p/3$")
	ax.set_xlabel(r"$p$", fontsize=label_size)
	ax.set_ylabel(label_str, fontsize=label_size)
	ax.grid()
	ax.legend(loc = 'upper left', prop={'size':10})
	ax.set_xlim(2.0 - 0.05, 6.0 + 0.05)


data_dir = str(sys.argv[1])

## ------ Input and Output Dirs
input_dir  = data_dir 
output_dir = "/home/enda/PhD/Shell_Model/Data/Thesis/" + "Plots/"
data_file_path = input_dir

# -----------------------------------------
# # --------  Read In data
# -----------------------------------------
print("Reading In Data")
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
print("Ploting Flux Scaling")
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
print("Plotting Tseries")
short_file = "/home/enda/PhD/Shell_Model/Data/Thesis/HD-INTFACRK4-FULL_N[25]_T[0.0,0.0001,100.000]_NU[5e-07]_ALPHA[1.500]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[N_SCALING]_TAG[Short]"
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

# # -----------------------------------------
# # # --------  Plot Spacetime Energy Flux
# # -----------------------------------------
# enrg_flux, hel_flux = get_fluxes(short_u, sys_vars.EPS, sys_vars.Lambda)
# print("Plotting Flux SpaceTime")
# fig = plt.figure(figsize = (21, 12))
# gs  = GridSpec(2, 1)
# ax1 = fig.add_subplot(gs[0, 0])
# im1 = ax1.imshow(np.rot90(sys_msr_data.k[np.newaxis, :] * enrg_flux[:, :], k =-1), extent=(short_time[0], short_time[-1], sys_vars.N, 1), aspect = 'auto', cmap = my_magma)
# ax1.set_xlabel(r"$t$")
# ax1.set_ylabel(r"$n$")
# # ax1.set_ylim(1.0, N)
# div1  = make_axes_locatable(ax1)
# cbax1 = div1.append_axes("right", size = "5%", pad = 0.05)
# cb1   = plt.colorbar(im1, cax = cbax1)
# cb1.set_label(r"$\Pi_n^{\mathcal{K}}$")

# ax2 = fig.add_subplot(gs[1, 0])
# im2 = ax2.imshow(np.rot90(sys_msr_data.k[np.newaxis, :] * hel_flux[:, :], k=-1), aspect = 'auto', cmap = my_magma)
# ax2.set_xlabel(r"$t$")
# ax2.set_ylabel(r"$n$")
# # ax2.set_ylim(1.0, N)
# div2  = make_axes_locatable(ax2)
# cbax2 = div2.append_axes("right", size = "5%", pad = 0.05)
# cb2   = plt.colorbar(im2, cax = cbax2)
# cb2.set_label(r"$\Pi_n^{\mathcal{H}}$")
# plt.savefig(output_dir + "Shell_SpacetimeFlux" + "." + fig_format, format = fig_format, bbox_inches='tight', dpi=1200)
# plt.close()



# -----------------------------------------
# # --------  Plot PDFs
# -----------------------------------------
print("Plotting PDFs")
nbins = 100
fig = plt.figure()
gs  = GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0, 0])
indxs = [5, 10, 15, 20]
for j, i in enumerate(indxs):
	if hasattr(stats_data, "vel_hist_counts"):
		pdf, centres = compute_pdf_from_hist(stats_data.vel_hist_counts[i, :], stats_data.vel_hist_ranges[i, :], normed = True, remove_zeros=True)
	else:
		if sys_vars.model_type == "PO" or sys_vars.model_type == "AO":
			pdf, centres = compute_pdf(np.real(run_data.a_n[:, i] * np.exp(1j * run_data.phi_n[:, i])), nbins = nbins, normed = True, remove_zeros=True)					
		else:
			pdf, centres = compute_pdf(np.real(run_data.u[:, i]), nbins = nbins, normed = True, remove_zeros=True)    
	ax1.plot(centres, pdf, label = "$n = {}$".format(i), color=plot_colours[(j + 5 - len(indxs)) * 2])    
ax1.set_xlabel(r"$\Re \left\{u_n\right\} / \langle \Re\left\{ u_n \right\}^2 \rangle^{1/2}$", fontsize=label_size)
ax1.set_yscale('log')
ax1.set_ylabel(r"PDF", fontsize=label_size)
ax1.grid()
ax1.legend()
plt.savefig(output_dir + "Shell_PDFs" + "." + fig_format, format = fig_format, bbox_inches='tight', dpi=1200)
plt.close()



# -----------------------------------------
# # --------  Plot Structure Functions
# -----------------------------------------
p_range = np.arange(2.0, 6.0 + 1, dtype=np.int64)
zeta_p = [0.3729, 0.7055, 0.9963, 1.245, 1.4475, 1.6218]
ns_zeta_p    = [0.7, 1, 1.27, 1.53, 1.78]

k = sys_msr_data.k/sys_vars.k0
inert_range = [3 - 1, 12 - 1]

fig = plt.figure(figsize=fig_size_2x2)
gs  = GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0, 0])
zeta_p, zeta_p_resi = plot_sf(fig, ax1, p_range, k, stats_data.vel_str_func[:, :] / stats_data.num_stats_steps, inert_range, r"\mathcal{S}_p^{u}", insert_fig = True, scaling = 'loge', ax_pos=[0.15, 0.55, 0.15, 0.15])
ax1 = fig.add_subplot(gs[0, 1])
plot_anom_scaling(fig, ax1, p_range, zeta_p[:], ns_zeta_p, r"$\zeta_p^{u}$")
ax1 = fig.add_subplot(gs[1, 0])
zeta_p, zeta_p_resi = plot_sf(fig, ax1, p_range, k, stats_data.vel_flux_str_func_abs[:, :, 0] / stats_data.num_stats_steps, inert_range, r"\mathcal{S}_p^{\Pi^{\mathcal{K}^{u}}}", insert_fig = True, scaling = 'loge')
ax1 = fig.add_subplot(gs[1, 1])
plot_anom_scaling(fig, ax1, p_range, zeta_p[:], ns_zeta_p, r"$\zeta_p^{\Pi^{\mathcal{K}^{u}}}$")
plt.savefig(output_dir + "Shell_SFs" + "." + fig_format, format = fig_format, bbox_inches='tight', dpi=1200)
plt.close()







print("\n\nFinished\n\n")