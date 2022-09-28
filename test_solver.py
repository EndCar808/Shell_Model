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

		def __init__(self, solver_tag = None, in_dir = "Data/Test/"):
			self.in_dir     = in_dir
			self.solver_tag = solver_tag

	## Initialize class
	cargs = cmd_args()

	try:
		## Gather command line arguments
		opts, args = getopt.getopt(argv, "i:t:", ["tag"])
	except Exception as e:
		print("[" + tc.R + "ERROR" + tc.Rst + "] ---> Incorrect Command Line Arguements.")
		print(e)
		sys.exit()

	## Parse command line args
	for opt, arg in opts:

		if opt in ['-i']:
			## Read in input directory
			cargs.in_dir = str(arg)
			print("Input Directory: " + tc.C + cargs.in_dir + tc.Rst)

		if opt in ['-t']:
			## Read in solver tag to search for test data
			cargs.solver_tag = str(arg)

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

	#####################
	##  READ IN DATA  ##
	#####################
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

	#################
	##  TEST DATA  ##
	#################
	## Check timestamps match to ensure we are comparing similar data
	for i in range(1, len(time)):
		if not np.allclose(time[i], time[i - 1]):
			print("Time for step {} and {} do not match".format(steps[i], steps[i - 1]))
			break

	## Compute error
	err = []
	for i in range(1, len(time)):
		err.append(np.sum(np.abs(data[i][-1,:]-data[i-1][-1,:])**2, axis=0))
		print(np.sum(np.abs(data[i][-1,:]-data[i-1][-1,:])**2, axis=0))

	plot_steps = steps[1:]
	## Get fit of the data
	p = np.polyfit(np.log(plot_steps[:-2]), np.log(err[:-2]), 1)

	#################
	##  PLOT DATA  ##
	#################
	plt.figure()
	plt.plot(plot_steps, err, 'o')
	plt.plot(plot_steps[:-2], np.exp(p[1]) * plot_steps[:-2]**p[0], '--', color='orangered',label="Erreur $\propto$ dt^{:.2f}".format(p[0]))
	plt.xscale('log')
	plt.yscale('log')
	plt.ylabel(r"Error at T")
	plt.xlabel(r"dt")
	plt.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
	plt.legend()

	plt.savefig(cmdargs.in_dir + "/Test_Solver_Error_vs_dt.png")
	plt.close()
