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
# from Plotting.functions import parse_cml, tc, sim_data, import_data, import_stats_data, import_sys_msr_data
sys.path.append('/home/enda/PhD/Shell_Model/Plotting')
from functions import tc
#################################
##          MISC               ##
#################################
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
            self.solver = False
            self.tag = "None"
            self.cmd_only = False
            self.print = True


    ## Initialize class
    cargs = cmd_args()

    try:
        ## Gather command line arguments
        opts, args = getopt.getopt(argv, "i:o:f:t:", ["plot", "solver"])
    except Exception as e:
        print("[" + tc.R + "ERROR" + tc.Rst + "] ---> Incorrect Command Line Arguements.")
        print(e)
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

        elif opt in ['--solver']:
            ## Read in solver indicator
            cargs.solver = True

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

	## Make output folder for plots
	cmdargs.out_dir_AVGAMP_SLOPE = cmdargs.out_dir + "AVGAMP_SLOPE_PLOTS/"
	if os.path.isdir(cmdargs.out_dir_AVGAMP_SLOPE) != True:
		print("Making folder:" + tc.C + " AVGAMP_SLOPE_PLOTS/" + tc.Rst)
		os.mkdir(cmdargs.out_dir_AVGAMP_SLOPE)


	## Solver variables
	executable = "Solver/bin/solver_phase_only"
	output_dir = cmdargs.out_dir
	collect_data = False
	## Main variables
	n                = 25
	t0               = 0
	t                = 1e2
	trans_iters      = 1
	trans_iters_frac = 0.01
	h                = 1e-4
	k0               = 0.05
	Lambda           = 2.0
	ep               = 0.5
	ep_m             = 0.0
	v                = 5e-7
	eta              = 0.
	u0               = "RANDOM"
	s_tag            = "Avg-Amp-Slope-Study"
	forcing          = "DELTA"
	force_k          = 1
	force_scale      = 0.1
	save             = 500
	## Slopes
	alpha = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
	beta  = [0.]
	# Extras
	cfl_cond       = 0
	c              = 0.9
	step_type      = 0
	hypervisc      = 0
	hypervisc_pow  = 2
	ekmn_hypo_diff = 0
	ekmn_hypo_pow  = -2





	#########################
	##      RUN SOLVER     ##
	#########################
	if cmdargs.solver:

	    ## Get the number of processes to launch
	    proc_limit = len(alpha)
	    print("\nNumber of Solver Processes Created = [" + tc.C + "{}".format(proc_limit) + tc.Rst + "]\n")

	    # Create output objects to store process error and output
	    if collect_data:
	        solver_output = []
	        solver_error  = []

	    ## Generate command list 
	    cmd_list = [["{} -o {} -n {} -s {:3.5f} -e {:3.5f} -T {} -T {} -c {} -c {:1.6f} -h {:1.16f} -h {} -a {:1.10f} -b {:1.10f} -w {:1.3f} -w {:1.3f} -y {:1.16f} -y {:1.16f} -v {:g} -v {} -v {:1.1f} -d {:g} -d {} -d {:1.1f} -i {} -t {} -f {} -f {} -f {:1.3f} -p {}".format(
	                                                                                                                                                executable, 
	                                                                                                                                                output_dir,
	                                                                                                                                                n,
	                                                                                                                                                t0, t, trans_iters, trans_iters_frac,
	                                                                                                                                                cfl_cond, c, 
	                                                                                                                                                h, step_type,
	                                                                                                                                                a,
	                                                                                                                                                b,
	                                                                                                                                                k0, Lambda,
	                                                                                                                                                ep, ep_m, 
	                                                                                                                                                v, hypervisc, hypervisc_pow, 
	                                                                                                                                                eta, ekmn_hypo_diff, ekmn_hypo_pow,
	                                                                                                                                                u0, 
	                                                                                                                                                s_tag, 
	                                                                                                                                                forcing, force_k, force_scale, 
	                                                                                                                                                save)] for a in alpha for b in beta]

	    if cmdargs.cmd_only:
	        print(tc.C + "\nSolver Commands:\n" + tc.Rst)
	        for c in cmd_list:
	            print(c)
	            print()
	    else:
	        ## Create grouped iterable of subprocess calls to Popen() - see grouper recipe in itertools
	        groups = [(Popen(cmd, shell = True, stdout = PIPE, stdin = PIPE, stderr = PIPE, universal_newlines = True) for cmd in cmd_list)] * proc_limit 

	        ## Loop through grouped iterable
	        for processes in zip_longest(*groups): 
	            for proc in filter(None, processes): # filters out 'None' fill values if proc_limit does not divide evenly into cmd_list
	                ## Print command to screen
	                print("\nExecuting the following command:\n\n\t" + tc.C + "{}\n".format(proc.args[0]) + tc.Rst)
	                
	                if cmdargs.print:
	                    ## Print output to terminal as it comes
	                    for line in proc.stdout:
	                        sys.stdout.write(line)

	                # Communicate with process to retrive output and error
	                [run_CodeOutput, run_CodeErr] = proc.communicate()

	                # Append to output and error objects
	                if collect_data:
	                    solver_output.append(run_CodeOutput)
	                    solver_error.append(run_CodeErr)
	                
	                ## Print both to screen
	                print(run_CodeOutput)
	                print(run_CodeErr)

	                ## Wait until all finished
	                proc.wait()

	        if collect_data:
	            # Get data and time
	            now = datetime.now()
	            d_t = now.strftime("%d%b%Y_%H:%M:%S")

	            # Write output to file
	            with open(par_runs_output_dir + "par_run_solver_output_{}_{}.txt".format(cmdargs.init_file.lstrip('InitFiles/').rstrip(".ini"), d_t), "w") as file:
	                for item in solver_output:
	                    file.write("%s\n" % item)

	            # Write error to file
	            with open(par_runs_output_dir + "par_run_solver_error_{}_{}.txt".format(cmdargs.init_file.lstrip('InitFiles/').rstrip(".ini"), d_t), "w") as file:
	                for i, item in enumerate(solver_error):
	                    file.write("%s\n" % cmd_list[i])
	                    file.write("%s\n" % item)


