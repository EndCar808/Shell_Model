import numpy as np
import os 
import h5py
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib import rc, rcParams
from functions import slope_fit
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
label_size     = 16

fig_size_1row = (textwidth, 0.4 * textwidth)
fig_size_2row = (textwidth, 0.8 * textwidth)

plot_colours = plt.cm.YlGnBu(np.linspace(0.3, 1.0, 10))

# my_cmap = mpl.colors.ListedColormap(cm.magma.colors[::-1])
my_cmap = matplotlib.colors.ListedColormap(cm.YlGnBu(np.arange(0,cm.YlGnBu.N)))
my_cmap.set_under(color = "white")


out_dir="/home/enda/PhD/Shell_Model/Data/Convergence_Study/"

t_steps = np.array([1e-4*2**n for n in range(-3,11)])[::-1]
print(t_steps)

## Load error data
data = np.load("/home/enda/PhD/Shell_Model/Data/Convergence_Study/Err.npy")[::-1]
err_inf = np.zeros((len(t_steps)))
for i in range(1, len(t_steps)):
	# diff = np.absolute(data[i - 1, 10] - data[i, 10])
	# err_inf[i] = diff
	# print(diff)
	tmp_err = np.sum(np.absolute(data[i - 1, :])**2) - np.sum(np.absolute(data[i, :])**2)
	# tmp_err = np.sum(np.absolute(data[i - 1, :] - data[i, :])**2)
	# tmp_err = np.amax(np.absolute(data[i - 1, :] - data[i, :]))
	err = np.absolute(tmp_err)
	err_inf[i] = err
	print(err)


# plt.plot(np.log2(t_steps[1:]), np.log2(err_inf[1:]))
# plt.plot(np.log2(t_steps[1:]), np.log2(t_steps[1:]**4))
fig = plt.figure(figsize=(0.8 * textwidth, 0.4 * textwidth))
gs  = GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0, 0])
plt.plot(t_steps[1:8], err_inf[1:8], '.-')
plt.plot(t_steps[2:7], t_steps[2:7]**4 * 10 **3, 'k--', label = r"$h^4$" )
plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{\infty}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
plt.savefig(out_dir + "Shell_Solver_Error" + "." + "pdf", format="pdf", bbox_inches = 'tight', dpi=1200)
plt.close()












## Get the timesteps
dt0 = 0.01
t_steps_list = []
for i in range(10):
    # print(dt0 / 2**i)
    t_steps_list.append(dt0 / 2**i)

t_steps = np.array(t_steps_list)
print(t_steps)
dt_strings = ['0.01', '0.005', '0.0025', '0.00125', '0.000625', '0.0003125', '0.00015625', '7.8125e-05', '3.90625e-05', '1.95313e-05']

N = 19
data = np.zeros((len(dt_strings), N))

out_dir = "/home/enda/PhD/Shell_Model/Data/Convergence_Study/"



## Read in data
u_n_all = []
time_all = []
infile = "/Main_HDF_Data.h5"
for i, dt in enumerate(dt_strings):

    # in_dir = "/home/enda/PhD/Shell_Model/Data/Convergence_Study/Full/HD-INTFACRK4-FULL_N[19]_T[0.0,{},1.000]_SMP[1,1,0]_NU[0]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[NONE,1,0.100]_u0[N_SCALING]_TAG[Err-{}]".format(dt, i + 1)
    in_dir = "/home/enda/PhD/Shell_Model/Data/Convergence_Study/Full/HD-INTFACRK4-FULL_N[19]_T[0.0,{},1.000]_SMP[1,1,0]_NU[1e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[NONE,1,0.100]_u0[N_SCALING]_TAG[Err-{}]".format(dt, i + 1)
    print(in_dir + infile)

    with h5py.File(in_dir + infile, "r") as in_file:
        u_n_last = in_file["VelModes"][-1, :]
        u_n_all.append(u_n_last)




full_err_inf = np.zeros((len(dt_strings)))
full_err_l2 = np.zeros((len(dt_strings)))
full_err_l2_alt = np.zeros((len(dt_strings)))
for i in range(1, len(dt_strings)):
    tmp_err_alt = np.sum(np.absolute(u_n_all[i - 1][:])**2) - np.sum(np.absolute(u_n_all[i][:])**2)
    err_alt = np.absolute(tmp_err_alt)
    full_err_l2_alt[i] = err_alt

    tmp_err_l2 = np.sum(np.absolute(u_n_all[i - 1][:] - u_n_all[i][:])**2)
    err_l2 = np.absolute(tmp_err_l2)
    full_err_l2[i]  = err_l2

    tmp_err_inf = np.amax(np.absolute(u_n_all[i - 1][:] - u_n_all[i][:]))
    # tmp_err_inf = np.amax(np.absolute(u_n_all[i - 1][:])) - np.amax(np.absolute(u_n_all[i][:]))
    err_inf = np.absolute(tmp_err_inf)
    full_err_inf[i] = err_inf

    # print(err_inf, err_l2, err_alt)

fig = plt.figure(figsize=(textwidth, 0.4 * textwidth))
gs  = GridSpec(1, 3)
ax1 = fig.add_subplot(gs[0, 0])
plt.plot(t_steps[1:], full_err_inf[1:], '.-')
plt.plot(t_steps[1:], t_steps[1:]**4, 'k--', label = r"$h^4$" )
plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{\infty}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
ax1 = fig.add_subplot(gs[0, 1])
plt.plot(t_steps[1:], full_err_l2[1:], '.-')
plt.plot(t_steps[1:], t_steps[1:]**4, 'k--', label = r"$h^4$" )
plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{2}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
ax1 = fig.add_subplot(gs[0, 2])
plt.plot(t_steps[1:], t_steps[1:]**4 *10*8, 'k--', label = r"$h^4$" )
plt.plot(t_steps[1:], full_err_l2_alt[1:], '.-')
plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{2}^{alt}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
plt.savefig(out_dir + "NS_Solver_Error_Full_All" + "." + "pdf", format="pdf", bbox_inches = 'tight', dpi=1200)
plt.close()

fig = plt.figure(figsize=(0.8 * textwidth, 0.4 * textwidth))
gs  = GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0, 0])
plt.plot(t_steps[1:-3], full_err_inf[1:-3], '.-')
plt.plot(t_steps[2:-4], t_steps[2:-4]**4 * 10**6 , 'k--', label = r"$h^4$" )
# plt.plot(t_steps[1:-3], full_err_l2[1:-3], '.-')
# plt.plot(t_steps[2:-4], t_steps[2:-4]**4 * 10**-3 , 'k--', label = r"$h^4$" )

print("Visc Linf")
pfit_slope, pfit_c, poly_resid = slope_fit(np.log(t_steps[1:-3]), np.log(full_err_inf[1:-3]), 1, -1)
print(pfit_slope, pfit_c, poly_resid)
# plt.plot(t_steps[1:-3], np.exp(pfit_c)*t_steps[1:-3]**(pfit_slope) , 'r--', label = r"$Linf slope={}$".format(pfit_slope) )
print("Visc L2")
pfit_slope, pfit_c, poly_resid = slope_fit(np.log(t_steps[1:-3]), np.log(full_err_l2[1:-3]), 1, -1)
print(pfit_slope, pfit_c, poly_resid)
# plt.plot(t_steps[2:-4], np.exp(pfit_c)*t_steps[2:-4]**(pfit_slope) , 'b--', label = r"$L2 slope={}$".format(pfit_slope) )

plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{\infty}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
plt.savefig(out_dir + "NS_Solver_Error_Full_Visc" + "." + "pdf", format="pdf", bbox_inches = 'tight', dpi=1200)
plt.close()




#---------------------------------------------------------------------------------



u_n_all = []
time_all = []
infile = "/Main_HDF_Data.h5"
for i, dt in enumerate(dt_strings):

    in_dir = "/home/enda/PhD/Shell_Model/Data/Convergence_Study/Full/HD-INTFACRK4-FULL_N[19]_T[0.0,{},1.000]_SMP[1,1,0]_NU[0]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[NONE,1,0.100]_u0[N_SCALING]_TAG[Err-{}]".format(dt, i + 1)
    # in_dir = "/home/enda/PhD/Shell_Model/Data/Convergence_Study/Full/HD-INTFACRK4-FULL_N[19]_T[0.0,{},1.000]_SMP[1,1,0]_NU[1e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[NONE,1,0.100]_u0[N_SCALING]_TAG[Err-{}]".format(dt, i + 1)
    print(in_dir + infile)

    with h5py.File(in_dir + infile, "r") as in_file:
        u_n_last = in_file["VelModes"][-1, :]
        u_n_all.append(u_n_last)




no_vis_full_err_inf = np.zeros((len(dt_strings)))
no_vis_full_err_l2 = np.zeros((len(dt_strings)))
no_vis_full_err_l2_alt = np.zeros((len(dt_strings)))
for i in range(1, len(dt_strings)):
    tmp_err_alt = np.sum(np.absolute(u_n_all[i - 1][:])**2) - np.sum(np.absolute(u_n_all[i][:])**2)
    err_alt = np.absolute(tmp_err_alt)
    no_vis_full_err_l2_alt[i] = err_alt

    tmp_err_l2 = np.sum(np.absolute(u_n_all[i - 1][:] - u_n_all[i][:])**2)
    err_l2 = np.absolute(tmp_err_l2)
    no_vis_full_err_l2[i]  = err_l2

    tmp_err_inf = np.amax(np.absolute(u_n_all[i - 1][:] - u_n_all[i][:]))
    # tmp_err_inf = np.amax(np.absolute(u_n_all[i - 1][:])) - np.amax(np.absolute(u_n_all[i][:]))
    err_inf = np.absolute(tmp_err_inf)
    no_vis_full_err_inf[i] = err_inf

    # print(err_inf, err_l2, err_alt)


fig = plt.figure(figsize=(0.8 * textwidth, 0.4 * textwidth))
gs  = GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0, 0])
plt.plot(t_steps[1:-3], no_vis_full_err_inf[1:-3], '.-')
plt.plot(t_steps[2:-4], t_steps[2:-4]**4 * 10**6 , 'k--', label = r"$h^4$" )
# plt.plot(t_steps[1:-3], no_vis_full_err_l2[1:-3], '.-')
# plt.plot(t_steps[2:-4], t_steps[2:-4]**4 * 10**-3 , 'k--', label = r"$h^4$" )

print("No Visc Linf")
pfit_slope, pfit_c, poly_resid = slope_fit(np.log(t_steps[1:-3]), np.log(no_vis_full_err_inf[1:-3]), 1, -1)
print(pfit_slope, pfit_c, poly_resid)
# plt.plot(t_steps[1:-3], np.exp(pfit_c)*t_steps[1:-3]**(pfit_slope) , 'r--', label = r"$Linf slope={}$".format(pfit_slope) )
print("No Visc L2")
pfit_slope, pfit_c, poly_resid = slope_fit(np.log(t_steps[1:-3]), np.log(no_vis_full_err_l2[1:-3]), 1, -1)
print(pfit_slope, pfit_c, poly_resid)
# plt.plot(t_steps[2:-4], np.exp(pfit_c)*t_steps[2:-4]**(pfit_slope) , 'b--', label = r"$L2 slope={}$".format(pfit_slope) )

plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{\infty}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
plt.savefig(out_dir + "NS_Solver_Error_Full_NoVisc" + "." + "pdf", format="pdf", bbox_inches = 'tight', dpi=1200)
plt.close()




fig = plt.figure(figsize=(0.8 * textwidth, 0.4 * textwidth))
gs  = GridSpec(1, 2, wspace=0.25)
ax1 = fig.add_subplot(gs[0, 0])
plt.plot(t_steps[1:-3], no_vis_full_err_inf[1:-3], '.-', label = r"$\nu = 0$" )
plt.plot(t_steps[2:-4], t_steps[2:-4]**4 * 10**5 , 'k--', label = r"$h^4$" )
plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{\infty}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
ax1 = fig.add_subplot(gs[0, 1])
plt.plot(t_steps[1:-3], full_err_inf[1:-3], '.-', label = r"$\nu \neq 0$")
plt.plot(t_steps[2:-4], t_steps[2:-4]**4 * 10**4 , 'k--', label = r"$h^4$" )
plt.xscale('log')
plt.yscale('log')
ax1.set_xlabel(r"$h$", fontsize=label_size)
ax1.set_ylabel(r"$\left\|  \epsilon_h \right\|_{\infty}$", fontsize=label_size)
ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
ax1.legend()
plt.savefig(out_dir + "NS_Solver_Error_Both" + "." + "pdf", format="pdf", bbox_inches = 'tight', dpi=1200)
plt.close()