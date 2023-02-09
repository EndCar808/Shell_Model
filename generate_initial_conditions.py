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
from Plotting.functions import tc, parse_cml, slope_fit
######################
##       MAIN       ##
######################
if __name__ == '__main__':

	# -------------------------------------
	# # --------- Parse Commnad Line
	# -------------------------------------
	cmdargs = parse_cml(sys.argv[1:])


	# -------------------------------------
	# # --------- Read In Time Averaged Amps
	# -------------------------------------
	with h5py.File(cmdargs.in_dir + "System_Measure_HDF_Data.h5", 'r') as in_file:

		if "TimeAveragedVelocityAmplitudes" in list(in_file.keys()):
			a_n_t_avg = in_file["TimeAveragedVelocityAmplitudes"][:]
		else:
			print("No Time Averaged Amplitudes in Data Directory!")
			sys.exit()
		if "k" in list(in_file.keys()):
			k = in_file["k"][:]
		
	# -------------------------------------
	# # --------- Compute the various ICs
	# -------------------------------------
	alpha_slopes = np.arange(0.0, 1.51, 0.1)

	inert_lim_low  = 5
	inert_lim_high = 16

	fig = plt.figure(figsize = (16, 8))
	gs  = GridSpec(1, 1)
	ax1 = fig.add_subplot(gs[0, 0])
	a_n_slope, a_n_c, _ = slope_fit(np.log2(k), np.log2(a_n_t_avg), inert_lim_low, inert_lim_high)
	ax1.plot(np.log2(k), np.log2(a_n_t_avg), label = "$a_n \sim k^{}$".format(np.around(a_n_slope, 3)), marker = None, markevery = 1)
	for beta in alpha_slopes:
		a_n_adjust = a_n_t_avg * (k ** (np.absolute(a_n_slope) - beta))
		p, = ax1.plot(np.log2(k), np.log2(a_n_adjust), marker = '.', markevery = 1)
		slope, c, _ = slope_fit(np.log2(k), np.log2(a_n_adjust), inert_lim_low, inert_lim_high)
		ax1.plot(np.log2(k[inert_lim_low:inert_lim_high]), np.log2(k[inert_lim_low:inert_lim_high])*slope + c, '--', label = r"$k^a \quad a = $ ${:0.3f}$".format(slope), color = p.get_color())
	ax1.legend()
	ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
	plt.savefig(cmdargs.out_dir + "AvgAmps_All_Slopes.png", bbox_inches='tight')
	plt.close()

	# -------------------------------------
	# # --------- Write Data to Output Directory
	# -------------------------------------
	## Write the original
	with h5py.File(cmdargs.out_dir + "InitialData_ALPHA[{:1.3f}].h5".format(np.absolute(a_n_slope)), 'w') as out_file:
		out_file.create_dataset("VelAmps", data = a_n_t_avg)


	## Write the slope adjusted initial amplitudes
	for a in alpha_slopes:
		a_n_adjust = a_n_t_avg * (k ** (np.absolute(a_n_slope) - a))
		with h5py.File(cmdargs.out_dir + "InitialData_ALPHA[{:1.3f}].h5".format(a), 'w') as out_file:
			out_file.create_dataset("VelAmps", data = a_n_adjust)