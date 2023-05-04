import numpy as np
import h5py
import sys
import os
from numba import njit
import matplotlib as mpl
# mpl.use('TkAgg') # Use this backend for displaying plots in window
mpl.use('Agg') # Use this backend for writing plots to file
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif']  = 'Computer Modern Roman'
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import getopt
from functions import tc, sim_data, import_stats_data
from plot_functions import plot_str_func_with_anom_scaling

np.seterr(divide = 'ignore') 
###############################
##       FUNCTION DEFS       ##
###############################
def parse_cml(argv):

    """
    Parses command line arguments
    """

    ## Create arguments class
    class cmd_args:

        """
        Class for command line arguments
        """

        def __init__(self, in_dir = None, out_dir = None, info_dir = None, plotting = False):
            self.in_dir         = in_dir
            self.out_dir_info   = out_dir
            self.in_file        = out_dir
            self.plotting       = plotting
            self.tag = "None"


    ## Initialize class
    cargs = cmd_args()

    try:
        ## Gather command line arguments
        opts, args = getopt.getopt(argv, "i:o:f:t:", ["plot"])
    except:
        print("[" + tc.R + "ERROR" + tc.Rst + "] ---> Incorrect Command Line Arguements.")
        sys.exit()

    ## Parse command line args
    for opt, arg in opts:

        if opt in ['-i']:
            ## Read input directory
            cargs.in_dir = str(arg)
            print("\nInput Folder: " + tc.C + "{}".format(cargs.in_dir) + tc.Rst)

            cargs.out_dir = str(arg)
            print("Output Folder: " + tc.C + "{}".format(cargs.out_dir) + tc.Rst)

        if opt in ['-o']:
            ## Read output directory
            cargs.out_dir = str(arg)
            print("Output Folder: " + tc.C + "{}".format(cargs.out_dir) + tc.Rst)

        elif opt in ['-f']:
            ## Read input directory
            cargs.in_file = str(arg)
            print("Input Post Processing File: " + tc.C + "{}".format(cargs.in_file) + tc.Rst)

        elif opt in ['--plot']:
            ## Read in plotting indicator
            cargs.plotting = True

        elif opt in ['-t']:
            cargs.tag = str(arg)

    return cargs

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


	cmdargs.out_dir_HD = cmdargs.out_dir + "HD_PLOTS/"
	if os.path.isdir(cmdargs.out_dir_HD) != True:
		print("Making folder:" + tc.C + " HD_PLOTS/" + tc.Rst)
		os.mkdir(cmdargs.out_dir_HD)
	# -----------------------------------------
	# # --------  Read In data
	# -----------------------------------------
	## Read in simulation parameters
	sys_vars = sim_data(data_file_path, method)

	k = np.zeros((sys_vars.N, ))
	for i in range(sys_vars.N):
		k[i] = sys_vars.k0 * sys_vars.Lambda**(i)

	## Read in stats data
	stats_data = import_stats_data(data_file_path, sys_vars, method)

	fig_size = (10, 8)
	fig_file_type = ".png"
	fig_format = "png"
	# num_pow = 6
	num_pow = stats_data.vel_str_func.shape[-1]

	## Plot The structure function
	inert_lim_low  = 3
	inert_lim_high = 18
	range_lims     = [inert_lim_low, inert_lim_high]
	log_func = 'loge'
	inert_range = range_lims

	###----------------------------------------- Plot SF
	input_data = [stats_data.vel_str_func[:, :], stats_data.vel_trip_prod_str_func_abs[:, :], stats_data.vel_flux_str_func_abs[:, :, 1], stats_data.vel_flux_str_func_abs[:, :, 0]]
	figure_names = [r"SF_U", r"SF_Trip", r"SF_HFlux", r"SF_EFlux"]
	for in_data, fig_name in zip(input_data, figure_names):
		print("SF: {}".format(fig_name))
		str_funcs   = in_data / stats_data.num_stats_steps
		vel_zeta_p, ns_zeta_p, vel_zeta_p_resid = plot_str_func_with_anom_scaling(cmdargs.out_dir_HD + fig_name + "." + fig_format, k, str_funcs, inert_range, insert_fig = True, scaling = log_func, fig_size = (16, 8))

	## Read in Velocities
	with h5py.File(cmdargs.in_dir + "Main_HDF_Data.h5", 'r') as in_file:
		num_t, num_shell = in_file["VelModes"].shape
		
		start=0
		stride=1000
		count=(num_t) // stride
		block=1
		
		u = in_file["VelModes"][h5py.MultiBlockSlice(start, stride, count, block), :]
		num_t_slice, num_shell_slice = u.shape


	tot_energy = np.sum(np.absolute(u)**2, axis=-1)
	tot_diss   = np.sum((k[np.newaxis, :]**2 * np.absolute(u)**2), axis=-1)
	eddy_turn  = 1.0 / (k[0] * np.absolute(u[:, 0]))
	
	time = np.linspace(sys_vars.T * 0.2, sys_vars.T, num_t_slice)

	fig = plt.figure(figsize = (32, 8))
	gs  = GridSpec(1, 3, hspace = 0.35)
	ax1 = fig.add_subplot(gs[0, 0])
	ax1.plot(time, tot_energy)
	ax1.set_xlabel(r"$t$")
	ax1.set_title(r"Total Energy")
	ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
	ax1 = fig.add_subplot(gs[0, 1])
	ax1.plot(time, tot_diss)
	ax1.set_xlabel(r"$t$")
	ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
	ax1.set_title(r"Total Dissipation")
	ax1 = fig.add_subplot(gs[0, 2])
	ax1.plot(time, eddy_turn, label=r"Mean = {}".format(np.mean(eddy_turn)))
	ax1.set_xlabel(r"$t$")
	ax1.set_title(r"Eddy Turnover")
	ax1.legend()
	ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
	fig.savefig(cmdargs.out_dir_HD  + "System_Measures" + "." + fig_format, format=fig_format, bbox_inches='tight')
	plt.close()