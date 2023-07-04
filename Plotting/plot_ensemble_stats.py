import numpy as np
import h5py
import sys
import os
from numba import njit
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import cm 
import time as TIME
from functions import compute_pdf, sim_data
from plot_functions import plot_str_func_with_anom_scaling, plot_PDF_InOne
from thesis_stats import plot_anom_scaling, plot_sf
from matplotlib import rc, rcParams
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


# --------------------------------------------------
# # --------  CMD Line & Data Directories
# --------------------------------------------------
## Get script name and check for input directories
script_name  = sys.argv[0]
out_dir_name = sys.argv[1]
tag_name     = sys.argv[2]
print("\nFile name: %s " % script_name)
print("\nTag name: %s " % tag_name)
print("\nOutdir name: %s " % out_dir_name)
if (len(sys.argv) == 1):
	print("No Input file specified, Error.\n")
	sys.exit()
data_no = 0

if not(os.path.isdir(out_dir_name)):
	print("Not a valid output dir: {}\n".format(out_dir_name))
	sys.exit()

## Initialize array to hold data directories
out_dir_data = np.array(['' for _ in range(len(sys.argv)-3)], dtype='object')
fig_format  = 'pdf'

## Get input data directories
for i in range(0,len(sys.argv)-3):
	if (len(sys.argv) > data_no+1):
		out_dir_data[data_no] = str(sys.argv[data_no+3])
		data_no += 1
print()

## Print directories to screen
dud_files = []
dud_file_indx = []
for i in range(0, data_no):
	print("%dth Input file: %s" % (i,out_dir_data[i]))
	if not (os.path.isfile(out_dir_data[i]+"/Stats_HDF_Data.h5")) or not os.path.isfile(out_dir_data[i]+"/SimulationDetails.txt"): 
		print("Cannot open %d th input file!" % (i))
		dud_files.append(out_dir_data[i])
		dud_file_indx.append(i)
print()

## Print and remove dud files
for i in range(len(dud_files)):
	print("Dud file no {}: {}".format(dud_file_indx[i], dud_files[i]))
print()

## Print true files list
true_data_no      = data_no - len(dud_files)
true_out_dir_data = []
dud_file          = False
for i in range(data_no):
	for j in range(len(dud_files)):
		if i == dud_file_indx[j]:
			dud_file = True

	if not dud_file:
		true_out_dir_data.append(out_dir_data[i])
	else:
		dud_file = False
	
true_out_dir_data = np.array(true_out_dir_data)	
for i in range(true_data_no):
	print("True File no. {}: {}".format(i, true_out_dir_data[i]))

print()


# plot_dir="/home/enda/PhD/Shell_Model/Data/Thesis/Plots/ReplaceData/"
plot_dir=out_dir_name

# --------------------------------------------------
# # --------  Initialize Data Arrays
# --------------------------------------------------
## Open first data directory to sizes of stats arrays for initialization
with h5py.File(true_out_dir_data[-1]+"/Stats_HDF_Data.h5",'r') as HDFfileRuntime:
	## Get the number of stats steps
	num_stats_steps = HDFfileRuntime["NumStatsSteps"][:]

	## Get the HD Stats Data
	vel_sf               = HDFfileRuntime["StructureFunctionVel"][:, :]
	num_shell, num_pow = vel_sf.shape

## Open first data directory to sizes of field arrays for initialization
with h5py.File(true_out_dir_data[-1]+"/Main_HDF_Data.h5",'r') as HDFfileRuntime:
	## Get the Velocity Modes
	u = HDFfileRuntime["VelModes"][:, :]
	num_t_steps = u.shape[0]

## Open first data directory to sizes of field arrays for initialization
with h5py.File(true_out_dir_data[-1]+"/System_Measure_HDF_Data.h5",'r') as HDFfileRuntime:
	k = HDFfileRuntime["k"][:]
	time = HDFfileRuntime["Time"][:]

## Get sim data
sys_vars = sim_data(true_out_dir_data[-1] + "/", "default")

## Initialize Arrays
Tot_Enrg_all 	 = np.empty([true_data_no, num_t_steps], dtype=float, order='C')
Vel_Modes_all    = np.empty([true_data_no, num_t_steps, num_shell], dtype=np.complex128, order='C')
Vel_SF_all       = np.empty([true_data_no, num_shell, num_pow], dtype=float, order='C')
Vel_Trip_SF_all  = np.empty([true_data_no, num_shell, num_pow], dtype=float, order='C')
Vel_EFlux_SF_all = np.empty([true_data_no, num_shell, num_pow], dtype=float, order='C')
Vel_HFlux_SF_all = np.empty([true_data_no, num_shell, num_pow], dtype=float, order='C')

Triads_all = np.empty([true_data_no, num_t_steps, 23], dtype=float, order='C')

# --------------------------------------------------
# # --------  Read In Data Arrays
# --------------------------------------------------
print("\nReading In Data")
skipped=0
for i in range(0, true_data_no):
	with h5py.File(true_out_dir_data[i]+"/Stats_HDF_Data.h5",'r') as HDFfileRuntime:
		## Get the number of stats steps
		num_stats_steps = HDFfileRuntime["NumStatsSteps"][:]

	with h5py.File(true_out_dir_data[i]+"/System_Measure_HDF_Data.h5",'r') as HDFfileRuntime:
		tmp = HDFfileRuntime["TotalEnergy"][:]

	if np.any(tmp > 75):
		skipped+=1
		continue

	with h5py.File(true_out_dir_data[i]+"/System_Measure_HDF_Data.h5",'r') as HDFfileRuntime:
		Tot_Enrg_all[i, :]     = HDFfileRuntime["TotalEnergy"][:]


		## Get the HD Stats Data
		Vel_SF_all[i, :, :]       = HDFfileRuntime["StructureFunctionVel"][:, :] / num_stats_steps
		Vel_Trip_SF_all[i, :, :]  = HDFfileRuntime["StructureFunctionTrippleProdVelAbs"][:, :] / num_stats_steps
		Vel_EFlux_SF_all[i, :, :] = HDFfileRuntime["StructureFunctionVelFluxAbs"][:, :, 0] / num_stats_steps
		Vel_HFlux_SF_all[i, :, :] = HDFfileRuntime["StructureFunctionVelFluxAbs"][:, :, 1] / num_stats_steps

	with h5py.File(true_out_dir_data[i]+"/Main_HDF_Data.h5",'r') as HDFfileRuntime:
		## Get the Velocity Modes
		Vel_Modes_all[i, :, :] = HDFfileRuntime["VelModes"][:, :]



	with h5py.File(true_out_dir_data[i]+"/Phase_Sync_HDF_Data.h5",'r') as HDFfileRuntime:
		Triads_all[i, :, :] = HDFfileRuntime["VelTriads"][:, :]

# --------------------------------------------------
# # --------  Plotting Data
# --------------------------------------------------
print("Plotting Data")

###------------------ Plot Total Energy Time Series 
fig = plt.figure(figsize = (12, 8))
gs  = GridSpec(1, 1, hspace = 0.35)
ax1 = fig.add_subplot(gs[0, 0])
for i in range(0, true_data_no, 10):
	ax1.plot(time, Tot_Enrg_all[i, :], label=r"i = {}".format(i))
ax1.set_xlabel(r"$t$")
ax1.set_title(r"Total Energy")
ax1.legend()
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
plt.suptitle("True No. of Data Dirs: {}/".format(true_data_no, data_no))
fig.savefig(plot_dir + "Ens_TotalEnergy_Tag[{}]".format(tag_name) + "." + fig_format, format=fig_format, bbox_inches='tight')
plt.close()

###------------------ Plot Structure Functions
inert_lim_low  = 3
inert_lim_high = 12
range_lims     = [inert_lim_low, inert_lim_high]
log_func = 'loge'

sf = np.mean(Vel_SF_all, axis=0)
print("SF --- Vel")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel_Tag[{}]".format(tag_name) + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))
sf = np.mean(Vel_Trip_SF_all, axis=0)
print("SF --- Vel Trip Product")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel_TripProd_Tag[{}]".format(tag_name) + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))
sf = np.mean(Vel_EFlux_SF_all, axis=0)
print("SF --- Vel Energy Flux")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel_EnergyFlux_Tag[{}]".format(tag_name) + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))
sf = np.mean(Vel_HFlux_SF_all, axis=0)
print("SF --- Vel Hel Flux")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel_HelFlux_Tag[{}]".format(tag_name) + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))


###------------------ Plot PDFs
plot_PDF_InOne(plot_dir + "EnsPDFs_Vel_Real_Tag[{}]".format(tag_name) + "." + fig_format, np.real(Vel_Modes_all), lab="\Re u_n", remove_zeros=True)





# --------------------------------------------------
# # --------  Plotting Data for Thesis
# --------------------------------------------------
###------------ Plot SF
p_range = np.arange(2.0, 6.0 + 1, dtype=np.int64)
zeta_p = [0.3729, 0.7055, 0.9963, 1.245, 1.4475, 1.6218]
ns_zeta_p    = [0.7, 1, 1.27, 1.53, 1.78]

kk = k/sys_vars.k0
inert_range = [3, 12]

fig = plt.figure(figsize=fig_size_2x2)
gs  = GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0, 0])
zeta_p, zeta_p_resi = plot_sf(fig, ax1, p_range, kk, np.mean(Vel_SF_all, axis=0)[:, 1:], inert_range, r"\mathcal{S}_p^{u}", insert_fig = True, scaling = 'loge', ax_pos=[0.15, 0.55, 0.15, 0.15])
ax1 = fig.add_subplot(gs[0, 1])
plot_anom_scaling(fig, ax1, p_range, zeta_p[:], ns_zeta_p, r"$\zeta_p^{u}$")
ax1 = fig.add_subplot(gs[1, 0])
zeta_p, zeta_p_resi = plot_sf(fig, ax1, p_range, kk, np.mean(Vel_EFlux_SF_all, axis=0)[:, 1:], inert_range, r"\mathcal{S}_p^{\Pi^{\mathcal{K}^{u}}}", insert_fig = True, scaling = 'loge')
ax1 = fig.add_subplot(gs[1, 1])
plot_anom_scaling(fig, ax1, p_range, zeta_p[:], ns_zeta_p, r"$\zeta_p^{\Pi^{\mathcal{K}^{u}}}$")
plt.savefig(plot_dir + "Shell_SFs" + "." + 'pdf', format = 'pdf', bbox_inches='tight', dpi=1200)
plt.close()


####------------- Plot PDFs
nbins = 1000
fig = plt.figure()
gs  = GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0, 0])
indxs = [5, 10, 15, 20]
for j, i in enumerate(indxs):
	pdf, centres = compute_pdf(np.real(Vel_Modes_all[:, :, i]).flatten(), nbins = nbins, normed = True, remove_zeros=True)    
	ax1.plot(centres, pdf, label = "$n = {}$".format(i), color=plot_colours[(j + 5 - len(indxs)) * 2])    
ax1.set_xlabel(r"$\Re \left\{u_n\right\} / \langle \Re\left\{ u_n \right\}^2 \rangle^{1/2}$", fontsize=label_size)
ax1.set_yscale('log')
ax1.set_ylabel(r"PDF", fontsize=label_size)
ax1.set_ylim(bottom=1e-5)
ax1.grid()
ax1.legend()
plt.savefig(plot_dir + "Shell_PDFs" + "." + 'pdf', format = 'pdf', bbox_inches='tight', dpi=1200)
plt.close()


###------------ Plot SF
p_range = np.arange(2.0, 6.0 + 1, dtype=np.int64)
zeta_p = [0.3729, 0.7055, 0.9963, 1.245, 1.4475, 1.6218]
ns_zeta_p    = [0.7, 1, 1.27, 1.53, 1.78]

kk = k/sys_vars.k0
inert_range = [3, 12]

fig = plt.figure()
gs  = GridSpec(1, 2)
ax1 = fig.add_subplot(gs[0, 0])
zeta_p, zeta_p_resi = plot_sf(fig, ax1, p_range, kk, np.mean(Vel_EFlux_SF_all, axis=0)[:, 1:], inert_range, r"\mathcal{S}_p^{\Pi^{\mathcal{K}}}", insert_fig = True, scaling = 'loge')
ax1 = fig.add_subplot(gs[0, 1])
plot_anom_scaling(fig, ax1, p_range, zeta_p[:], ns_zeta_p, r"$\zeta_p^{\Pi^{\mathcal{K}}}$")
plt.savefig(plot_dir + "PO_Shell_SFs" + "." + 'pdf', format = 'pdf', bbox_inches='tight', dpi=1200)
plt.close()

fig = plt.figure()
gs  = GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(kk, np.mean(np.absolute(Vel_Modes_all), axis=(0, 1)), color = plot_colours[-3])
ax1.plot(kk[3:15], kk[3:15]**(-0.3333), 'k--', label=r"$k_n^{-1/3}$")
ax1.set_xlabel(r"$k_n / k_0$")
ax1.set_ylabel(r"$a_n$")
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend()
ax1.grid()
plt.savefig(plot_dir + "PO_Amps" + "." + 'pdf', format = 'pdf', bbox_inches='tight', dpi=1200)
plt.close()

fig = plt.figure(figsize = (16, 16))
gs  = GridSpec(5, 5, wspace = 0.35, hspace = 0.25)
for i in range(5):
	for j in range(5):
		if i * 5 + j < 23:
			ax1 = fig.add_subplot(gs[i, j])
			pdf, ranges = np.histogram(Triads_all[:, :, i * 5 + j].flatten(), bins=1000, range=(0.0, 2.0*np.pi), density=True)
			centres = (ranges[1:] + ranges[:-1]) * 0.5
			ax1.plot(centres, pdf, label = "$n = {}$".format(i * 5 + j + 1), color=plot_colours[-3])
			ax1.set_xlabel("$\phi_n + \phi_{n + 1} + \phi_{n + 2}$")
			ax1.set_ylabel("PDF")
			ax1.set_yscale("log")
			ax1.legend()
			ax1.set_xlim(0, 2.0*np.pi)
			# ax1.set_xlim(-np.pi, np.pi)
			ax1.set_xticks([0.0, np.pi/2.0, np.pi, 1.5*np.pi, 2.0 * np.pi])
			ax1.set_xticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2 \pi$"])
			# ax1.set_xticklabels([-np.pi, -np.pi/2.0, 0.0, np.pi/2, np.pi])
			# ax1.set_xticklabels([r"$-\pi$", r"$\frac{-\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"])
			ax1.grid()
plt.savefig(plot_dir + "PO_Phases_Vel_Triad_PDF" + "." + fig_format, format = fig_format, bbox_inches='tight', dpi=1200)
plt.close()


print("True: {} - Skipped: {}".format(true_data_no, skipped))

print("\n\nFinshied\n\n")