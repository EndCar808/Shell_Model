#!/usr/bin/env python    
#######################
##  LIBRARY IMPORTS  ##
#######################
from configparser import ConfigParser
import numpy as np
import os
import sys
import getopt
import distutils.util as utils
from matplotlib.gridspec import GridSpec
from datetime import datetime
from collections.abc import Iterable
from itertools import zip_longest
from subprocess import Popen, PIPE
from Plotting.functions import tc, sim_data, import_data
import matplotlib.pyplot as plt
#########################
##  READ COMMAND LINE  ##
#########################
def parse_cml(argv):

	"""
	Parses command line arguments
	"""

	## Create arguments class
	class cmd_args:
		"""
		Class for command line arguments
		"""

		def __init__(self, solver_tag = None, in_dir = "Data/Test/", err = False, compare = False):
			self.in_dir     = in_dir
			self.out_dir     = in_dir
			self.solver_tag = solver_tag
			self.err 		= err
			self.compare 	= compare

	## Initialize class
	cargs = cmd_args()

	try:
		## Gather command line arguments
		opts, args = getopt.getopt(argv, "i:t:", ["tag", "err", "compare"])
	except Exception as e:
		print("[" + tc.R + "ERROR" + tc.Rst + "] ---> Incorrect Command Line Arguements.")
		print(e)
		sys.exit()

	## Parse command line args
	for opt, arg in opts:

		if opt in ['-t']:
			## Read in solver tag to search for test data
			cargs.solver_tag = str(arg)

		if opt in ['--err']:
			## Read in solver tag to search for test data
			cargs.err = True

		if opt in ['--compare']:
			## Read in solver tag to search for test data
			cargs.compare = True

	return cargs

######################
##       MAIN       ##
######################
if __name__ == '__main__':
	##########################
	##  PARSE COMMAND LINE  ##
	##########################
	cmdargs = parse_cml(sys.argv[1:])

	#####################
	##  GET FILE LIST  ##
	#####################
	## Read in list of folders
	test_dirs = []
	for folder in os.listdir(cmdargs.in_dir):
		if folder.split('_')[-1].lstrip('TAG[').rstrip(']') == str(cmdargs.solver_tag):
			for k in folder.split('_'):
				if k.startswith('T['):
					test_dirs.append(str(k.split(',')[1]) + "_" + folder)

	## Sort folders
	test_dirs.sort(reverse = True, key = lambda x: float(x.split("_")[0]))
	# for dirs in test_dirs:
	# 	print(dirs)

	print("\nNo. of Data Directories: {}\n".format(len(test_dirs)))

	print("Measure Accuracy: {}\nCompare Solver: {}\n".format(cmdargs.err, cmdargs.compare))
	########################
	##  MEASURE ACCURACY  ##
	########################
	if cmdargs.err:
		## Read in Data
		steps = []
		time  = []
		data  = []
		for f in test_dirs:
			input_dir = cmdargs.in_dir + f.split("_", 1)[-1] + "/"

			## Read in simulation parameters
			sys_vars = sim_data(input_dir, "default")

			## Read in solver data
			run_data = import_data(input_dir, sys_vars, "default")

			time.append(run_data.time)
			data.append(run_data.u)
			steps.append(sys_vars.dt)

		## Check timestamps match to ensure we are comparing similar data
		for i in range(1, len(time)):
			if not np.allclose(time[i], time[i - 1]):
				print("Time for step {} and {} do not match".format(steps[i], steps[i - 1]))
				break

		## Compute error
		err_at_T = []
		for i in range(1, len(time)):
			err_at_T.append(np.sum(np.abs(data[i][-1,:]-data[i-1][-1,:])**2, axis=0))
			print(np.sum(np.abs(data[i][-1,:]-data[i-1][-1,:])**2, axis=0))

		plot_steps = steps[1:]
		## Get fit of the data
		p = np.polyfit(np.log(plot_steps[:]), np.log(err_at_T[:]), 1)

		## Plot data
		plt.figure()
		plt.plot(plot_steps, err_at_T, 'o')
		plt.plot(plot_steps[:], np.exp(p[1]) * plot_steps[:]**p[0], '--', color='orangered',label="Error $\propto$ dt^{:.2f} (Best fit)".format(p[0]))
		plt.plot(plot_steps[:], np.array(plot_steps[:])**round(p[0]), "--", color = 'black', label = "$h^{}$".format(int(round(p[0]))))
		plt.xscale('log')
		plt.yscale('log')
		plt.ylabel(r"Error at T")
		plt.xlabel(r"dt")
		plt.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		plt.legend()

		plt.savefig(cmdargs.in_dir + "/TestResults/Error_vs_dt_" + cmdargs.solver_tag + ".png")
		plt.close()

	########################
	##  COMPARE SOLVERS  ##
	########################
	## Read in data
	data        = []
	time        = []
	solver_type = []
	if cmdargs.compare:
		for f in test_dirs:
			input_dir = cmdargs.in_dir + f.split("_", 1)[-1] + "/"

			solv_type = str(input_dir.split("_")[2].split("-")[1])

			## Read in simulation parameters
			sys_vars = sim_data(input_dir, "default")

			## Read in solver data
			run_data = import_data(input_dir, sys_vars, "default")

			time.append(run_data.time)
			data.append(run_data.u)
			solver_type.append(solv_type)


		## Plot data
		fig = plt.figure(figsize = (16, 8))
		gs  = GridSpec(1, 3)
		ax1 = fig.add_subplot(gs[0, 0])
		ax1.plot(np.linalg.norm(np.real(data[0][:, :]) - np.real(data[1][:, :]), axis = -1), '--', label = "Error Between {} & {}".format(solver_type[0], solver_type[1]))
		ax1.set_ylabel(r"Error")
		ax1.set_xlabel(r'Iterations')
		ax1.set_yscale('log')
		ax1.set_xscale('log')
		ax1.legend()
		ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax2 = fig.add_subplot(gs[0, 1])
		ax2.plot(np.linalg.norm(np.real(data[0][:, :]) - np.real(data[2][:, :]), axis = -1), '--', label = "Error Between {} & {}".format(solver_type[0], solver_type[2]))
		ax2.set_ylabel(r"Error")
		ax2.set_xlabel(r'Iterations')
		ax2.set_yscale('log')
		ax2.set_xscale('log')
		ax2.legend()
		ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
		ax3 = fig.add_subplot(gs[0, 2])
		ax3.plot(np.linalg.norm(np.real(data[1][:, :]) - np.real(data[2][:, :]), axis = -1), '--', label = "Error Between {} & {}".format(solver_type[1], solver_type[2]))
		ax3.set_ylabel(r"Error")
		ax3.set_xlabel(r'Iterations')
		ax3.set_yscale('log')
		ax3.set_xscale('log')
		ax3.legend()
		ax3.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)

		plt.savefig(cmdargs.out_dir + "/TestResults/Compare_Solvers_" + cmdargs.solver_tag + ".png")
		plt.close()