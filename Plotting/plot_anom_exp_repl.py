import matplotlib as mpl
# mpl.use('TkAgg') # Use this backend for displaying plots in window
mpl.use('Agg') # Use this backend for writing plots to file
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import h5py
import sys
import os
import numpy as np
from plot_functions import plot_str_func_with_anom_scaling
import matplotlib.pyplot as plt

repl_indx = [1, 10, 25, 100, 500, 1000, 10000]

## Plot The structure function
inert_lim_low  = 3
inert_lim_high = 18
range_lims     = [inert_lim_low, inert_lim_high]
log_func = 'loge'
inert_range = range_lims


outdir = "/home/enda/PhD/Shell_Model/Data/Thesis/Plots/"

sf = []

for repl in repl_indx:

	indir = "/home/enda/PhD/Shell_Model/Data/Thesis/AO_FxdPhase/HD-INTFACRK4-AOFXDPHASE_N[25]_T[200.0,0.0001,1000.000]_NU[5e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[AO_INPUT_PHASE_REPLACE]_TAG[NEW-REPL-Test-{}]/".format(repl)

	with h5py.File(indir + "System_Measure_HDF_Data.h5", 'r') as infile:
		k = infile["k"][:]

	with h5py.File(indir + "Stats_HDF_Data.h5", 'r') as infile:
		num_stats_steps = infile["NumStatsSteps"][0]

		enrg_flux_str = infile["StructureFunctionVelFluxAbs"][:, :, 0]

		str_funcs   = enrg_flux_str / num_stats_steps
		vel_zeta_p, ns_zeta_p, vel_zeta_p_resid = plot_str_func_with_anom_scaling(outdir + "EFlux_Repl_{}".format(repl) + "." + "png", k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))
		sf.append(vel_zeta_p)


p=np.arange(2, 6 + 1)
plt.figure()
plt.plot(p, p/3, 'k--', label=r"K41")
plt.plot(p, ns_zeta_p, label=r"NS")
for i in range(len(sf)):
	plt.plot(p, sf[i][1:] / sf[i][2], label=r"$m = {}$".format(repl_indx[i]))
plt.xlabel("$p$")
plt.ylabel("$\zeta_p / \zeta_3$")
plt.legend()
plt.grid()
plt.savefig(outdir + "Anomalous_REPL")
