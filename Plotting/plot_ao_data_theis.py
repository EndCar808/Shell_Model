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
from matplotlib.pyplot import cm 
from functions import tc, sim_data, import_data, import_stats_data, import_sys_msr_data, import_phase_sync_data, compute_pdf_from_hist, compute_pdf, slope_fit, compute_u_flux
from plot_functions import plot_anomalous_exponent, plot_str_funcs_with_slope, phase_only_space_time
from plot_functions import plot_str_func_with_anom_scaling
from thesis_stats import plot_sf, plot_anom_scaling


from functions import compute_field_values
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

my_cmap_rev = matplotlib.colors.ListedColormap(cm.YlGnBu(np.arange(0,cm.YlGnBu.N))[::-1])
my_cmap_rev.set_under(color = "white")


fig_format = 'pdf'



jitter_dir = "/home/enda/PhD/Shell_Model/Data/Thesis/AO_FxdPhase/HD-INTFACRK4-AOFXDPHASE_N[25]_T[0.0,5e-05,10000.000]_NU[5e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[AO_ALGND_PHASE]_TAG[Jitter-Test]/"

plot_dir = "/home/enda/PhD/Shell_Model/Data/Thesis/AO_Thesis_Plots/"
sys_vars = sim_data(jitter_dir, "Default")

with h5py.File(jitter_dir + "/Stats_HDF_Data.h5") as infile:
	sf = infile["StructureFunctionVelFluxAbs"][:, :, 0]

with h5py.File(jitter_dir + "/System_Measure_HDF_Data.h5") as infile:
	k = infile["k"][:]

log_func = np.log
p_range = np.arange(2, 6.0 + 1, dtype=np.int64)
zeta_p = [0.3729, 0.7055, 0.9963, 1.245, 1.4475, 1.6218]
ns_zeta_p    = [0.7, 1, 1.27, 1.53, 1.78]

kk = k/sys_vars.k0
inert_range = [3, 12]

fig = plt.figure()
gs  = GridSpec(1, 2)
ax1 = fig.add_subplot(gs[0, 0])
jit_zeta_p, zeta_p_resi = plot_sf(fig, ax1, p_range, kk, sf[:, 1:], inert_range, r"\mathcal{S}_p^{\Pi^{\mathcal{K}}}", insert_fig = True, scaling = 'loge')
# for idx, i in enumerate(p_range):
# 	p, = ax1.plot(log_func(k[:sf.shape[0]]), log_func(sf[:, idx]), label = "$p = {}$".format(int(i)), color=plot_colours[(i - 2) * 2])

ax1 = fig.add_subplot(gs[0, 1])
plot_anom_scaling(fig, ax1, p_range, jit_zeta_p[:], ns_zeta_p, r"$\zeta_p^{\Pi^{\mathcal{K}}}$")
plt.savefig(plot_dir + "AO_Shell_SFs" + "." + 'pdf', format = 'pdf', bbox_inches='tight', dpi=1200)
plt.close()






# repl_indx = [1, 10, 100, 1000, 10000]
# new_sf = []
# for repl in repl_indx:
# 	print(repl)

# 	indir = "/home/enda/PhD/Shell_Model/Data/Thesis/AO_FxdPhase/HD-INTFACRK4-AOFXDPHASE_N[25]_T[200.0,0.0001,1000.000]_NU[5e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[AO_INPUT_PHASE_REPLACE]_TAG[NEW-REPL-Test-{}]/".format(repl)

# 	with h5py.File(indir + "System_Measure_HDF_Data.h5", 'r') as infile:
# 		k = infile["k"][:]

# 	with h5py.File(indir + "Stats_HDF_Data.h5", 'r') as infile:
# 		num_stats_steps = infile["NumStatsSteps"][0]

# 		enrg_flux_str = infile["StructureFunctionVelFluxAbs"][:, :, 0]

# 		str_funcs   = enrg_flux_str / num_stats_steps
# 		vel_zeta_p, ns_zeta_p, vel_zeta_p_resid = plot_str_func_with_anom_scaling(plot_dir + "EFlux_Repl_{}".format(repl) + "." + "png", k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))
# 		new_sf.append(vel_zeta_p)


# p=np.arange(2, 6 + 1)
# plt.figure()
# # plt.plot(p, ns_zeta_p, label=r"NS")
# for i in range(len(new_sf)):
# 	plt.plot(p, new_sf[i][1:] / new_sf[i][2], label=r"$m = {}$".format(repl_indx[i]), color=plot_colours[i*2])
# plt.plot(p, jit_zeta_p[1:] / jit_zeta_p[2], label=r"$m = \infty$", color='grey')
# plt.plot(p, p/3, 'k--', label=r"K41")
# plt.xlabel("$p$", fontsize=label_size)
# plt.xlim(2-0.1, 6+0.1)
# plt.ylabel("$\zeta_p / \zeta_3$", fontsize=label_size)
# plt.legend()
# plt.grid()
# plt.savefig(plot_dir + "Anomalous_REPL.pdf", bbox_inches='tight', dpi=1200)
