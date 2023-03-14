#!/usr/bin/env python

# Author: Enda Carroll
# Date: Sept 2021
# Info: Script to compare solver results with decaying turbulence papers
#       Solver data

#######################
#  Library Imports  ##
#######################
import numpy as np
import h5py as h5
import sys
import os
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
sys.path.append('/home/enda/PhD/Shell_Model/Plotting')
from functions import parse_cml, tc, sim_data, compute_pdf
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import cm
from subprocess import Popen, PIPE, run
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Computer Modern Roman'
import concurrent.futures as cf
np.seterr(divide = 'ignore') 


def plot_amp_snap(outdir, a_n, ymin, ymax, title, t):

	print(t)

	fig = plt.figure(figsize=(12, 8))
	gs  = GridSpec(1, 1, hspace=0.4, wspace=0.5)
	ax1 = fig.add_subplot(gs[0, 0])
	ax1.plot(k, a_n)
	ax1.set_ylim(ymin, ymax)
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xlabel(r"$k_n$")
	ax1.set_ylabel(r"$a_n$")
	ax1.set_title(title)
	fig.savefig(outdir + "AMP_SNAP_{:05d}.png".format(t), bbox_inches='tight')
	plt.close()

######################
#       MAIN       ##
######################
if __name__ == '__main__':
	# -------------------------------------
	# # --------- Parse Commnad Line
	# -------------------------------------
	cmdargs = parse_cml(sys.argv[1:])
	method = "default"
	data_file_path = cmdargs.in_dir

	# -----------------------------------------
	# # --------  Read In data
	# -----------------------------------------
	# Read in simulation parameters
	sys_vars = sim_data(data_file_path, method)



	with h5.File(cmdargs.in_dir + "Main_HDF_Data.h5", 'r') as in_file:

		if "VelModes" in list(in_file.keys()):
			u = in_file["VelModes"][:, :]
	
	# with h5.File(cmdargs.in_dir + "System_Measure_HDF_Data.h5", 'r') as in_file:

	# 	if "k" in list(in_file.keys()):
	# 		k = in_file["k"][:]
	# 	if "TimeAveragedVelocityAmplitudes" in list(in_file.keys()):
	# 		a_n_t_avg = in_file["TimeAveragedVelocityAmplitudes"][:]

	# cmdargs.out_dir_AMP = cmdargs.out_dir + "AMP_PLOTS/"
	# if os.path.isdir(cmdargs.out_dir_AMP) != True:
	#     print("Making folder:" + tc.C + " AMP_PLOTS/" + tc.Rst)
	#     os.mkdir(cmdargs.out_dir_AMP)
	# cmdargs.out_dir_AMP_AMP = cmdargs.out_dir_AMP + "AMP_SNAPS/"
	# if os.path.isdir(cmdargs.out_dir_AMP_AMP) != True:
	#     print("Making folder:" + tc.C + " AMP_SNAPS/" + tc.Rst)
	#     os.mkdir(cmdargs.out_dir_AMP_AMP)


	# # -----------------------------------------
	# # # --------  Plot data
	# # -----------------------------------------
	# num_bins  = 100
	# norm_hist = False

	# ## 1D - PDF - Amp
	# fig = plt.figure(figsize=(24, 24))
	# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
	# for i in range(5):
	# 	for j in range(5):
	# 		indx = i * 5 + j
	# 		if indx < u.shape[-1]:
	# 			ax1 = fig.add_subplot(gs[i, j])
	# 			pdf, centres = compute_pdf(np.absolute(u[:, indx]), nbins=num_bins, normed=norm_hist, remove_zeros=False)
	# 			p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
	# 			ax1.set_xlabel(r"$ \left|" +  "u_n" + r" \right|$")
	# 			ax1.set_ylabel(r"PDF")
	# 			ax1.set_yscale('log')
	# 			ax1.set_title("n = {}".format(indx + 1))
	# fig.savefig(cmdargs.out_dir_AMP + "VelModes" + "_1D_PDF_Amp" + ".png", bbox_inches='tight')
	# plt.close()

	# ## 1D - CDF - Amp
	# fig = plt.figure(figsize=(24, 24))
	# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
	# cdf = []
	# centres_p = []
	# for i in range(5):
	# 	for j in range(5):
	# 		indx = i * 5 + j
	# 		if indx < u.shape[-1]:
	# 			ax1 = fig.add_subplot(gs[i, j])
	# 			pdf, centres = compute_pdf(np.absolute(u[:, indx]), nbins=num_bins, normed=norm_hist, remove_zeros=False)
	# 			centres_p.append(centres)
	# 			cdf.append(np.cumsum(pdf) * (centres[1] - centres[0]))
	# 			p, = ax1.plot(centres, np.cumsum(pdf) * (centres[1] - centres[0]), label="$n = {}$".format(indx + 1))  
	# 			ax1.plot(np.sort(np.absolute(u[:, indx])), np.arange(len(np.absolute(u[:, indx]))) / len(np.absolute(u[:, indx])), ':', label=r"Direct")  
	# 			ax1.set_xlabel(r"$ \left|" +  "u_n" + r" \right|$")
	# 			ax1.set_ylabel(r"CDF")
	# 			ax1.set_title("n = {}".format(indx + 1))
	# 			ax1.set_ylim(0-.1, 1.0+.1)
	# 			ax1.legend()
	# fig.savefig(cmdargs.out_dir_AMP + "VelModes" + "_1D_CDF_Amp" + ".png", bbox_inches='tight')
	# plt.close()


	with h5.File("/home/enda/PhD/Shell_Model/Data/InputAmplitudes/" + "Amp_Data_N[{}].h5".format(sys_vars.N), "w") as out_f:
		for t in range(sys_vars.ndata):
			out_f.create_dataset("Amp_t{}".format(t), data = np.absolute(u[t, :]))


	# with h5.File("/home/enda/PhD/Shell_Model/Data/RandomAmplitudes/" + "Amp_CDF_Data_N[{}].h5".format(sys_vars.N), "w") as out_f:
	# 	for n in range(sys_vars.N):
	# 		out_f.create_dataset("X_values_{}".format(n), data = centres_p[n])
	# 		out_f.create_dataset("CDF_values_{}".format(n), data = cdf[n])


	# ## Focus on a_1
	# cdf_a_1   = cdf[0]
	# n_samples = int(1e5)
	# tmp_pdf   = np.empty_like(pdf)
	# for i in range(n_samples):
	# 	uni = np.random.rand()
	# 	bin_indx = np.argmin(cdf_a_1 <= uni)
	# 	tmp_pdf[bin_indx] += 1
	# tmp_pdf = tmp_pdf / np.sum(tmp_pdf)

	# fig = plt.figure(figsize=(12, 8))
	# gs  = GridSpec(1, 2, hspace=0.4, wspace=0.5)
	# ax1 = fig.add_subplot(gs[0, 0])
	# pdf, centres = compute_pdf(np.absolute(u[:, 0]), nbins=num_bins, normed=norm_hist, remove_zeros=False)
	# p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
	# p, = ax1.plot(centres, tmp_pdf / (centres[1] - centres[0]), label="Empirical")    
	# ax1.set_xlabel(r"$ \left|" +  "u_n" + r" \right|$")
	# ax1.set_ylabel(r"PDF")
	# ax1.set_yscale('log')
	# ax1.set_title("Emperical Vs PDF")
	# ax1.legend()
	# ax1 = fig.add_subplot(gs[0, 1])
	# p, = ax1.plot(centres, cdf_a_1, label="$n = {}$".format(indx + 1))  
	# uni  = np.random.rand()
	# indx = np.argmin(cdf_a_1 <= uni)
	# print(uni, centres[indx])
	# ax1.axhline(y = uni, xmin = centres.min(), xmax = centres.max())
	# ax1.axvline(x = centres[indx], ymin = 0.0, ymax = 1.0)
	# ax1.plot(np.sort(np.absolute(u[:, 0])), np.arange(len(np.absolute(u[:, 0]))) / len(np.absolute(u[:, 0])), label=r"Direct")  

	# ax1.set_xlabel(r"$ \left|" +  "u_n" + r" \right|$")
	# ax1.set_ylabel(r"CDF")
	# ax1.set_title("n = {}".format(indx + 1))
	# ax1.set_ylim(0-.1, 1.0+.1)
	# ax1.legend()
	# fig.savefig(cmdargs.out_dir_AMP + "VelModes" + "_1D_PDF_a_1" + ".png", bbox_inches='tight')
	# plt.close()


	# ## Simulate The amplitudes
	# rnd_a_n_t_avg = np.zeros((sys_vars.N, ))
	# a_n           = np.zeros((sys_vars.ndata, sys_vars.N))	
	# print(sys_vars.ndata)

	# # fig = plt.figure(figsize=(12, 8))
	# # gs  = GridSpec(1, 1, hspace=0.4, wspace=0.5)
	# # ax1 = fig.add_subplot(gs[0, 0])
	# for t in range(sys_vars.ndata):

	# 	for n in range(sys_vars.N):
	# 		# Get CDF and centres
	# 		x_n, cdf_n = centres_p[n], cdf[n]
	# 		# x_n, cdf_n = np.sort(np.absolute(u[:, n])), np.arange(len(np.absolute(u[:, n]))) / len(np.absolute(u[:, n]))

	# 		# Get random sample
	# 		uni      = np.random.rand()
	# 		bin_indx = np.argmin(cdf_n <= uni)

	# 		# Set a_n to randomly sampled variate
	# 		a_n[t, n] = x_n[bin_indx]

	# 		# Time average it
	# 		rnd_a_n_t_avg[n] += a_n[t, n]

	# # 	ax1.plot(k, a_n[t, :], 'r', alpha = 0.15)
	# # rnd_a_n_t_avg /= sys_vars.ndata
	# # ax1.plot(k, rnd_a_n_t_avg, 'k')
	# # ax1.set_xscale('log')
	# # ax1.set_yscale('log')
	# # fig.savefig(cmdargs.out_dir_AMP + "VelModes" + "_Sim_an_TimeAvg_All" + ".png", bbox_inches='tight')
	# # plt.close()

	# fig = plt.figure(figsize=(12, 8))
	# gs  = GridSpec(1, 1, hspace=0.4, wspace=0.5)
	# ax1 = fig.add_subplot(gs[0, 0])
	# p, = ax1.plot(k, rnd_a_n_t_avg, label="Randomly Generated")
	# p, = ax1.plot(k, a_n_t_avg, ':', label="Full Model")
	# ax1.set_xscale('log')
	# ax1.set_yscale('log')
	# ax1.legend()
	# fig.savefig(cmdargs.out_dir_AMP + "VelModes" + "_Sim_an_TimeAvg" + ".png", bbox_inches='tight')
	# plt.close()

	# num_threads = 25
	
	# ## Create process Pool executer
	# executor = cf.ProcessPoolExecutor(num_threads)

	# ## Submit jobs to the executor
	# # data = np.absolute(u)
	# data = a_n
	# step = 100
	# ymin = data[:step:, :].min()
	# ymax = data[:step:, :].max()
	# print(ymin, ymax)
	# futures = [executor.submit(plot_amp_snap, cmdargs.out_dir_AMP_AMP, data[t, :], ymin, ymax, "Random Amplitudes", t//step) for t in range(0, sys_vars.ndata, step)]

	# ## Wait until these jobs are finished
	# cf.wait(futures)
	# # for f in futures:
	# #     print(f.result())


	# framesPerSec = 15
	# inputFile    = cmdargs.out_dir_AMP_AMP + "AMP_SNAP_%05d.png"
	# videoName    = cmdargs.out_dir_AMP_AMP + "Rand_Amp_Evolution.mp4"
	# cmd = "ffmpeg -y -r {} -f image2 -s 1920x1080 -i {} -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec libx264 -crf 25 -pix_fmt yuv420p {}".format(framesPerSec, inputFile, videoName)

	# process = Popen(cmd, shell = True, stdout = PIPE, stdin = PIPE, universal_newlines = True)
	# [runCodeOutput, runCodeErr] = process.communicate()
	# print(runCodeOutput)
	# print(runCodeErr)
	# process.wait()

	# ## Prin summary of timmings to screen
	# print("\n" + tc.Y + "Finished making video..." + tc.Rst)
	# print("Video Location: " + tc.C + videoName + tc.Rst + "\n")


	# ## Remove the generated snaps after video is created
	# run("cd {};".format(cmdargs.out_dir_AMP_AMP) + "rm {};".format("./*.png") + "cd -;", shell = True)
