#######################
##  Library Imports  ##
#######################
import numpy as np
import h5py
import sys
from numba import njit
import matplotlib as mpl
# mpl.use('TkAgg') # Use this backend for displaying plots in window
# mpl.use('Agg') # Use this backend for writing plots to file
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif']  = 'Computer Modern Roman'
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Plotting.functions import parse_cml, slope_fit
from scipy.optimize import curve_fit, least_squares
from Plotting.functions import compute_pdf


@njit
def set_phases(N, jitt):	
	phi = np.zeros((N + 4, ))
	for i in range(0, N):
		n = i + 2
		if n == 2:
			phi[n] = np.pi/4.0
		elif n == 3:
			phi[n] = 2.0 * np.pi * np.random.rand()
		elif n < N + 2:
			jitter = jitt[i]
			phi[n] = (-np.pi/2 + jitter) - phi[n - 1] - phi[n - 2]
	return phi[2:N+2]

@njit
def get_triads(phi):
	num_shell = phi.shape[-1]

	tri = np.zeros((num_shell - 2))
	diff = np.zeros((num_shell - 3))
	for i in range(num_shell - 2):
		tri[i] = phi[i] + phi[i + 1] + phi[i + 2]
		if i < num_shell - 3:
			diff[i] = phi[i] - phi[i + 3]
	return tri, diff

@njit
def get_triads_phase_diffs(num_t, N, triads_sd):

	phi        = np.zeros((num_t, N))
	triads     = np.zeros((num_t, N - 2))
	phase_diff = np.zeros((num_t, N - 3))
	tmp_jit = np.zeros((N - 2, ))
	for t in range(num_t):

		# phase_jitt = np.random.uniform(low=-1.0, high=1.0, size=N)
		for i in range(N - 2):
			tmp_jit[i] = ((1.0 - (-1.0)) * np.random.rand() + (-1.0)) * triads_sd[i]
		# phase_jitt = np.concatenate([tmp_jit, [0.0, 0.0]])
		phase_jitt = tmp_jit

		phi[t, :] = set_phases(N, phase_jitt)
		triads[t, :], phase_diff[t, :] = get_triads(phi[t, :])
	return phi, triads, phase_diff

######################
##       MAIN       ##
######################
if __name__ == '__main__':

	# -------------------------------------
	# # --------- Parse Commnad Line
	# -------------------------------------
	cmdargs = parse_cml(sys.argv[1:])

	# -------------------------------------
	# # --------- Read In Data
	# -------------------------------------
	with h5py.File(cmdargs.in_dir + "Main_HDF_Data.h5", 'r') as in_file:
		u = in_file["VelModes"][:, :]
		num_t, num_shell = u.shape
	with h5py.File(cmdargs.in_dir + "System_Measure_HDF_Data.h5", 'r') as in_file:
		time = in_file["Time"][:]

	# tag = "Long"

	for t in range(11):
		print(time[t])
		for n in range(num_shell):
			print("u[{},{}]: {:1.16f} {:1.16f}\ta: {:1.16f}\tp: {:1.16f}".format(t, n, np.real(u[t, n]), np.imag(u[t, n]), np.absolute(u[t, n]), np.angle(u[t, n])))
		print()

	# with h5py.File("/home/enda/PhD/Shell_Model/Data/InputPhases/" + "Phase_Data_N[{}]_T[{}]_TAG[{}].h5".format(num_shell, 5e3, tag), "w") as out_f:
	# 	for t in range(num_t):
	# 		if t == 0:
	# 			out_f.create_dataset("Amp_t{}".format(t), data = np.absolute(u[t, :]))
	# 		out_f.create_dataset("Phase_t{}".format(t), data = np.angle(u[t, :]))

	# with h5py.File(cmdargs.in_dir + "Phase_Sync_HDF_Data.h5", 'r') as in_file:
	# 	vel_triads     = in_file["VelTriads"][:, :]
	# 	vel_phase_diff = in_file["VelPhaseDifferences"][:, :]

	# phases     = np.angle(u)
	# phases_avg = np.mean(phases, axis=0)

	# triads_avg = np.mean(vel_triads, axis=0)
	# triads_sd  = np.std(vel_triads, axis=0)


	# out_dir = "/home/enda/PhD/Shell_Model/Data/InitialConditions/TimeAveragedPhases/"
	# # with h5py.File(out_dir + "InitPhases.h5", 'w') as out_file:
	# # 	out_file.create_dataset("Vel_Phases_Avg", data=phases_avg)
	# # 	out_file.create_dataset("Triads_Avg", data=np.concatenate([triads_sd, [0.0, 0.0]]))
	# # 	out_file.create_dataset("Triads_Std", data=np.concatenate([triads_sd, [0.0, 0.0]]))

	# fig   = plt.figure(figsize = (21, 9))
	# gs    = GridSpec(1, 3)
	# ax1   = fig.add_subplot(gs[0, 0])
	# ax1.plot(phases_avg)
	# ax2   = fig.add_subplot(gs[0, 1])
	# ax2.plot(triads_avg)
	# ax2   = fig.add_subplot(gs[0, 2])
	# ax2.plot(triads_sd)
	# plt.savefig(out_dir + "AvgPhases.png")
	# plt.close()

	# num_bins  = 100
	# norm_hist = False
	# input_data   = [phases, vel_triads, vel_phase_diff]
	# figure_names = [r"Vel_Phases", r"Vel_Triads", r"Vel_PhaseDiff"]
	# data_labels  = [r"u_n", r"u_{n + 2}u_{n + 1}u_{n}", r"u_{n}u_{n + 3}^*"]
	# for in_data, fig_name, data_labs in zip(input_data, figure_names, data_labels):
	# 	fig        = plt.figure(figsize=(24, 24))
	# 	gs         = GridSpec(5, 5, hspace=0.4, wspace=0.5)
	# 	for i in range(5):
	# 		for j in range(5):
	# 			indx = i * 5 + j
	# 			if indx < in_data.shape[-1]:
	# 				ax1 = fig.add_subplot(gs[i, j])
	# 				pdf, centres = compute_pdf(np.mod(in_data[:, indx] + 2.0 * np.pi, 2.0 * np.pi), bin_lims=(0.0, 2.0 * np.pi), nbins=num_bins, normed=norm_hist)
	# 				p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
	# 				ax1.set_xlabel(r"$ \arg \left\{" +  data_labs + r" \right\}$")
	# 				ax1.set_xlim(0 - 0.1, 2.0*np.pi + 0.1)
	# 				ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
	# 				ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
	# 				ax1.set_ylabel(r"PDF")
	# 				ax1.set_yscale('log')
	# 				ax1.set_title("n = {}".format(indx + 1))
	# 	fig.savefig(out_dir + fig_name + "_1D_PDF_Angles" + "." + 'png', format = 'png', bbox_inches='tight')
	# 	plt.close()


	# N          = 25
	# num_triads = N - 2
	# num_diff   = N - 3

	# phi, triads, phase_diffs = get_triads_phase_diffs(num_t, N, triads_sd)

	# num_bins  = 100
	# norm_hist = False

	# input_data   = [phi, triads, phase_diffs]
	# figure_names = [r"Phases", r"Triads", r"PhaseDiff"]
	# data_labels  = [r"u_n", r"u_{n + 2}u_{n + 1}u_{n}", r"u_{n}u_{n + 3}^*"]
	# for in_data, fig_name, data_labs in zip(input_data, figure_names, data_labels):
	# 	fig        = plt.figure(figsize=(24, 24))
	# 	gs         = GridSpec(5, 5, hspace=0.4, wspace=0.5)
	# 	for i in range(5):
	# 		for j in range(5):
	# 			indx = i * 5 + j
	# 			if indx < in_data.shape[-1]:
	# 				ax1 = fig.add_subplot(gs[i, j])
	# 				pdf, centres = compute_pdf(np.mod(in_data[:, indx] + 2.0 * np.pi, 2.0 * np.pi), bin_lims=(0.0, 2.0 * np.pi), nbins=num_bins, normed=norm_hist)
	# 				p, = ax1.plot(centres, pdf, label="$n = {}$".format(indx + 1))    
	# 				ax1.set_xlabel(r"$ \arg \left\{" +  data_labs + r" \right\}$")
	# 				ax1.set_xlim(0 - 0.1, 2.0*np.pi + 0.1)
	# 				ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
	# 				ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
	# 				ax1.set_ylabel(r"PDF")
	# 				ax1.set_yscale('log')
	# 				ax1.set_title("n = {}".format(indx + 1))
	# 	fig.savefig(out_dir + fig_name + "_1D_PDF_Angles" + "." + 'png', format = 'png', bbox_inches='tight')
	# 	plt.close()


	# #################################################
	# #
	# #
	# #	Estimate PDF via CDF
	# #
	# #
	# #################################################
	# num_bins  = 100
	# norm_hist = False

	# ## 1D - CDF - Amp
	# fig = plt.figure(figsize=(24, 24))
	# gs = GridSpec(5, 5, hspace=0.4, wspace=0.5)
	# cdf = []
	# centres_p = []
	# for i in range(5):
	# 	for j in range(5):
	# 		indx = i * 5 + j
	# 		print(indx)
	# 		if indx < vel_triads.shape[-1]:
	# 			ax1 = fig.add_subplot(gs[i, j])
	# 			pdf, centres = compute_pdf(vel_triads[:, indx], nbins=num_bins, normed=norm_hist, remove_zeros=False)
	# 			centres_p.append(centres)
	# 			cdf.append(np.cumsum(pdf) * (centres[1] - centres[0]))
	# 			p, = ax1.plot(centres, np.cumsum(pdf) * (centres[1] - centres[0]), label="$n = {}$".format(indx + 1))  
	# 			ax1.plot(np.sort(vel_triads[:, indx]), np.arange(len(vel_triads[:, indx])) / len(vel_triads[:, indx]), ':', label=r"Direct")  
	# 			ax1.set_xlabel(r"$ \varphi_{n, n + 1, n + 2}$")
	# 			ax1.set_ylabel(r"CDF")
	# 			ax1.set_title("n = {}".format(indx + 1))
	# 			ax1.set_ylim(0-.1, 1.0+.1)
	# 			ax1.legend()
	# fig.savefig(out_dir + "Triads" + "_1D_CDF_Amp" + ".png", bbox_inches='tight')
	# plt.close()


	# # with h5.File("/home/enda/PhD/Shell_Model/Data/RandomPhases/" + "TPhase_CDF_Data_N[{}].h5".format(sys_vars.N), "w") as out_f:
	# # 	for n in range(sys_vars.N):
	# # 		out_f.create_dataset("X_values_{}".format(n), data = centres_p[n])
	# # 		out_f.create_dataset("CDF_values_{}".format(n), data = cdf[n])


	# ## Focus on a_1
	# cdf_indx = -1
	# cdf_a_1   = cdf[cdf_indx]
	# print(cdf_a_1.shape, len(cdf))
	# n_samples = int(1e5)
	# tmp_pdf   = np.zeros((num_bins, ))
	# for i in range(n_samples):
	# 	uni = np.random.rand()
	# 	bin_indx = np.argmin(cdf_a_1 <= uni)
	# 	tmp_pdf[bin_indx] += 1
	# tmp_pdf = tmp_pdf / np.sum(tmp_pdf)

	# fig = plt.figure(figsize=(12, 8))
	# gs  = GridSpec(1, 2, hspace=0.4, wspace=0.5)
	# ax1 = fig.add_subplot(gs[0, 0])
	# pdf, centres = compute_pdf(vel_triads[:, cdf_indx], nbins=num_bins, normed=norm_hist, remove_zeros=False)
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
	# ax1.axhline(y = uni, xmin = centres.min(), xmax = centres.max(), linestyle=':', color='r')
	# ax1.axvline(x = centres[indx], ymin = 0.0, ymax = 1.0, linestyle=':', color='r')
	# # ax1.plot(np.sort(vel_triads[:, cdf_indx]), np.arange(len(vel_triads[:, cdf_indx])) / len(vel_triads[:, cdf_indx]), label=r"Direct")  

	# ax1.set_xlabel(r"$ \left|" +  "u_n" + r" \right|$")
	# ax1.set_ylabel(r"CDF")
	# ax1.set_title("n = {}".format(indx + 1))
	# ax1.set_ylim(0-.1, 1.0+.1)
	# ax1.legend()
	# fig.savefig(out_dir + "Triads" + "_1D_PDF_a_1" + ".png", bbox_inches='tight')
	# plt.close()







