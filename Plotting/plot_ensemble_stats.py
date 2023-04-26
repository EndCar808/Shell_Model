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
import time as TIME
from functions import compute_pdf, sim_data
from plot_functions import plot_str_func_with_anom_scaling, plot_PDF_InOne
# --------------------------------------------------
# # --------  CMD Line & Data Directories
# --------------------------------------------------
## Get script name and check for input directories
script_name       = sys.argv[0]
print("\nFile name: %s " % script_name)
if (len(sys.argv) == 1):
	print("No Input file specified, Error.\n")
	sys.exit()
data_no = 0

## Initialize array to hold data directories
out_dir_data = np.array(['' for _ in range(len(sys.argv)-1)], dtype='object')
fig_format  = 'png'

## Get input data directories
for i in range(0,len(sys.argv)-1):
	if (len(sys.argv) > data_no+1):
		out_dir_data[data_no] = str(sys.argv[data_no+1])
		data_no += 1
print()
## Print directories to screen
for i in range(0, data_no):
	print("%dth Input file: %s" % (i,out_dir_data[i]))
	if ( not (os.path.isfile(out_dir_data[i]+"/Stats_HDF_Data.h5")) ) : 
		print("Cannot open %d th input file, Error, will now exit!\n" % (i))
		sys.exit()	
print()

plot_dir="/home/enda/PhD/Shell_Model/Data/Thesis/Plots/"

# --------------------------------------------------
# # --------  Initialize Data Arrays
# --------------------------------------------------
## Open first data directory to sizes of stats arrays for initialization
with h5py.File(out_dir_data[0]+"/Stats_HDF_Data.h5",'r') as HDFfileRuntime:
	## Get the number of stats steps
	num_stats_steps = HDFfileRuntime["NumStatsSteps"][:]

	## Get the HD Stats Data
	vel_sf               = HDFfileRuntime["StructureFunctionVel"][:, :]
	num_shell, num_pow = vel_sf.shape

## Open first data directory to sizes of field arrays for initialization
with h5py.File(out_dir_data[0]+"/Main_HDF_Data.h5",'r') as HDFfileRuntime:
	## Get the Velocity Modes
	u = HDFfileRuntime["VelModes"][:, :]
	num_t_steps = u.shape[0]

## Open first data directory to sizes of field arrays for initialization
with h5py.File(out_dir_data[0]+"/System_Measure_HDF_Data.h5",'r') as HDFfileRuntime:
	k = HDFfileRuntime["k"][:]
	time = HDFfileRuntime["Time"][:]

## Get sim data
sys_vars = sim_data(out_dir_data[0] + "/", "default")

## Initialize Arrays
Tot_Enrg_all 	 = np.empty([data_no, num_t_steps], dtype=float, order='C')
Vel_Modse_all    = np.empty([data_no, num_t_steps, num_shell], dtype=np.complex128, order='C')
Vel_SF_all       = np.empty([data_no, num_shell, num_pow], dtype=float, order='C')
Vel_Trip_SF_all  = np.empty([data_no, num_shell, num_pow], dtype=float, order='C')
Vel_EFlux_SF_all = np.empty([data_no, num_shell, num_pow], dtype=float, order='C')
Vel_HFlux_SF_all = np.empty([data_no, num_shell, num_pow], dtype=float, order='C')



# --------------------------------------------------
# # --------  Read In Data Arrays
# --------------------------------------------------
print("\nReading In Data")
for i in range(0, data_no):
	with h5py.File(out_dir_data[0]+"/Stats_HDF_Data.h5",'r') as HDFfileRuntime:
		## Get the number of stats steps
		num_stats_steps = HDFfileRuntime["NumStatsSteps"][:]

		## Get the HD Stats Data
		Vel_SF_all[i, :, :]       = HDFfileRuntime["StructureFunctionVel"][:, :] / num_stats_steps
		Vel_Trip_SF_all[i, :, :]  = HDFfileRuntime["StructureFunctionTrippleProdVelAbs"][:, :] / num_stats_steps
		Vel_EFlux_SF_all[i, :, :] = HDFfileRuntime["StructureFunctionVelFluxAbs"][:, :, 0] / num_stats_steps
		Vel_HFlux_SF_all[i, :, :] = HDFfileRuntime["StructureFunctionVelFluxAbs"][:, :, 1] / num_stats_steps

	with h5py.File(out_dir_data[0]+"/Main_HDF_Data.h5",'r') as HDFfileRuntime:
		## Get the Velocity Modes
		Vel_Modse_all[i, :, :] = HDFfileRuntime["VelModes"][:, :]

	with h5py.File(out_dir_data[0]+"/System_Measure_HDF_Data.h5",'r') as HDFfileRuntime:
		Tot_Enrg_all[i, :]     = HDFfileRuntime["TotalEnergy"][:]



# --------------------------------------------------
# # --------  Plotting Data
# --------------------------------------------------
print("Plotting Data")

###------------------ Plot Total Energy Time Series 
fig = plt.figure(figsize = (12, 8))
gs  = GridSpec(1, 1, hspace = 0.35)
ax1 = fig.add_subplot(gs[0, 0])
for i in range(data_no):
	ax1.plot(time, Tot_Enrg_all[i, :], label=r"i = {}".format(i))
ax1.set_xlabel(r"$t$")
ax1.set_title(r"Total Energy")
ax1.legend()
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
fig.savefig(plot_dir + "Ens_TotalEnergy" + "." + fig_format, format=fig_format, bbox_inches='tight')
plt.close()

###------------------ Plot Structure Functions
inert_lim_low  = 3
inert_lim_high = 12
range_lims     = [inert_lim_low, inert_lim_high]
log_func = 'loge'

sf = np.mean(Vel_SF_all, axis=0)
print("SF --- Vel")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel" + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))
sf = np.mean(Vel_Trip_SF_all, axis=0)
print("SF --- Vel Trip Product")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel_TripProd" + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))
sf = np.mean(Vel_EFlux_SF_all, axis=0)
print("SF --- Vel Energy Flux")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel_EnergyFlux" + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))
sf = np.mean(Vel_HFlux_SF_all, axis=0)
print("SF --- Vel Hel Flux")
plot_str_func_with_anom_scaling(plot_dir + "EnsSF_Vel_HelFlux" + "." + fig_format, k/sys_vars.k0, sf, range_lims, insert_fig = True, scaling = 'loge', fig_size = (21, 9))


###------------------ Plot PDFs
plot_PDF_InOne(plot_dir + "EnsPDFs_Vel_Real" + "." + fig_format, np.real(Vel_Modse_all), lab="\Re u_n")






print("\n\nFinshied\n\n")