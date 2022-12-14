#!/usr/bin/env python    

## Author: Enda Carroll
## Date: Sept 20222
## Info: General python functions for analysing
#        Solver data


#######################
##  Library Imports  ##
#######################
import numpy as np
import h5py
import sys
import os
import re
from numba import njit
import pyfftw
from collections.abc import Iterable
from itertools import zip_longest
from subprocess import Popen, PIPE
import matplotlib as mpl
# mpl.use('TkAgg') # Use this backend for displaying plots in window
mpl.use('Agg') # Use this backend for writing plots to file
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif']  = 'Computer Modern Roman'
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import cm

#################################
## Colour Printing to Terminal ##
#################################
class tc:
    H         = '\033[95m'
    B         = '\033[94m'
    C         = '\033[96m'
    G         = '\033[92m'
    Y         = '\033[93m'
    R         = '\033[91m'
    Rst       = '\033[0m'
    Bold      = '\033[1m'
    Underline = '\033[4m'



###################################
##          STATS PLOT           ##
###################################
def plot_str_funcs_with_slope(outdir_path, k, str_funcs, inert_range):

	## Set up figure
	fig   = plt.figure(figsize = (16, 8))
	gs    = GridSpec(1, 1)
	ax1   = fig.add_subplot(gs[0, 0])

	## Set up for insert figure
	x0     = 0.15 
	y0     = 0.15
	width  = 0.3
	height = 0.2
	## Add insert
	ax1in = fig.add_axes([x0, y0, width, height])

	## Slopes of structure functions
	ns_zeta_p = [0.72, 1, 1.273, 1.534, 1.786]
	zeta_p    = []

	## Get inertial range limits
	inert_lim_low  = inert_range[0]
	inert_lim_high = inert_range[-1]

	## Marker style list
	mark_style = ['o','s','^','x','D','p']

	## Loop over structure functions and plot
	print("Power\tp/3\t\tDNS Slope\tNS")
	for i in range(str_funcs.shape[-1]):
	    
	    ## Plot strucure function
	    p, = ax1.plot(np.log2(k), np.log2(str_funcs[:, i]) + i * 10, label = "$p = {}$".format(i + 1), marker = mark_style[i], markerfacecolor = 'None', markersize = 5.0, markevery = 1)
	    
	    ## Find polynomial fit and plot
	    pfit_info  = np.polyfit(np.log2(k[inert_lim_low:inert_lim_high]), np.log2(str_funcs[inert_lim_low:inert_lim_high, i]) + i * 10, 1)
	    pfit_slope = pfit_info[0]
	    pfit_c     = pfit_info[1]
	    zeta_p.append(np.absolute(pfit_slope))
	    
	    ## Plot slope ddata to screen
	    if i >= 1:
	        print(" {}\t {:1.4f} \t {:1.4f} \t{:1.3f}".format(i + 1, -(i + 1) / 3, pfit_slope, -ns_zeta_p[i - 1]))
	    else:
	        print(" {}\t {:1.4f} \t {:1.4f}".format(i + 1, -(i + 1) / 3, pfit_slope))

	    ## Plot slopes computed from llinear fit    
	    ax1.plot(np.log2(k[inert_lim_low:inert_lim_high]), np.log2(k[inert_lim_low:inert_lim_high])*pfit_slope + pfit_c + 1.5, '--', color = p.get_color())
	    
	    ## Compute the local derivative and plot in insert
	    d_str_func  = np.diff(np.log2(str_funcs[:, i]))
	    d_k         = np.diff(np.log2(k))
	    local_deriv = d_str_func / d_k
	    local_deriv = np.concatenate((local_deriv, [(np.log2(str_funcs[-1, i]) - np.log2(str_funcs[-2, i])) / (np.log2(k[-1] - np.log2(k[-2])))]))
	    
	    ## Plot local slopes
	    ax1in.plot(np.log2(k), local_deriv, color = p.get_color(), marker = mark_style[i], markerfacecolor = 'None', markersize = 5.0, markevery = 2)
	    ax1in.set_ylabel(r"$\zeta_p$", labelpad = -40)
	    ax1in.set_xlabel(r"$log2(k_n)$", labelpad = -30)

	## Set axis labels 
	ax1.set_xlabel(r"$log_2 (k_n)$")
	ax1.set_ylabel(r"$log_2 (S_p(k_n))$")

	## Set grids
	ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
	ax1in.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)

	## Turn on axes
	ax1.legend()

	## Save figure
	plt.savefig(outdir_path, bbox_inches='tight')
	plt.close()


	return zeta_p, ns_zeta_p


def plot_anomalous_exponent(outdir_path, p, zeta_p, ns_zeta_p = None, label_str = "GOY Model"):

	## Setup figure
	fig   = plt.figure(figsize = (16, 8))
	gs    = GridSpec(1, 1)
	ax1   = fig.add_subplot(gs[0, 0])

	## Marker style list
	mark_style = ['o','s','^','x','D','p']

	## Get p range if not provided
	if p is None:
		p = np.arange(2, len(zeta_p) + 3)

	if ns_zeta_p is None:
		ns_zeta_p = [0.72, 1, 1.273, 1.534, 1.786]


	## Plot the structure function slopes
	ax1.plot(p, zeta_p, marker = mark_style[0], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = label_str)
	
	## Plot the NS structure function slopes
	ax1.plot(p, ns_zeta_p, marker = mark_style[1], markerfacecolor = 'None', markersize = 5.0, markevery = 1, label = "Navier Stokes")

	## Plot the She-Leveque prediction
	ax1.plot(p, p/9.0 + 2.0 * (1.0 - np.power(2.0/3.0, p/3.0)), label = "She-Leveque Model; $p/9 + 2(1 - (2/3)^{p/3})$")

	## Plot the slopes according to K41 theory
	ax1.plot(p, p / 3, 'k--', label = "K41")

	## Plot formatting
	# ax1.set_xlim(2.0, 6.0)
	ax1.set_xlabel(r"$p$")
	ax1.set_ylabel(r"$\zeta_p$")
	ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
	ax1.legend()

	## Save figure
	plt.savefig(outdir_path, bbox_inches='tight')
	plt.close()