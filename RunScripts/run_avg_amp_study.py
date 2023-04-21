#!/usr/bin/env python

# Author: Enda Carroll
# Date: Sept 2021
# Info: Script to compare solver results with decaying turbulence papers
#       Solver data

#######################
#  Library Imports  ##
#######################
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
sys.path.append('/home/enda/PhD/Shell_Model/Plotting')
from functions import parse_cml, tc, sim_data, import_data, import_stats_data, import_sys_msr_data, compute_pdf, compute_u_flux, compute_str_func
# from Plotting.functions import parse_cml, tc, sim_data, import_data, import_stats_data, import_sys_msr_data
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Computer Modern Roman'

from scipy.optimize import curve_fit

np.seterr(divide = 'ignore') 


def slope_fit(x, y, low, high):
    """Fits Slope over the inertial range."""
    poly_output = np.polyfit(x[low:high], y[low:high], 1, full=True)
    pfit_info = poly_output[0]
    poly_resid = poly_output[1][0]
    pfit_slope = pfit_info[0]
    pfit_c = pfit_info[1]

    return pfit_slope, pfit_c, poly_resid

def spec_func(k, a, b, c, d):
    return a * k**(-b) * np.exp(-c * k) + d


######################
#       MAIN       ##
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

    # Make output folder for plots
    cmdargs.out_dir_AVGAMP = cmdargs.out_dir + "AVGAMP_PLOTS/"
    if os.path.isdir(cmdargs.out_dir_AVGAMP) != True:
        print("Making folder:" + tc.C + " AVGAMP_PLOTS/" + tc.Rst)
        os.mkdir(cmdargs.out_dir_AVGAMP)

    # -----------------------------------------
    # # --------  Read In data
    # -----------------------------------------
    # Read in simulation parameters
    sys_vars = sim_data(data_file_path, method)

    # Read in solver data
    run_data = import_data(data_file_path, sys_vars, method)

    # Read in stats data
    stats_data = import_stats_data(data_file_path, sys_vars, method)

    # Read in sys_msr data
    sys_msr_data = import_sys_msr_data(data_file_path, sys_vars, method)

    # -----------------------------------------
    # # --------  Generate Averaged Amp Field
    # -----------------------------------------
    if hasattr(sys_msr_data, 'a_n_t_avg'):
        fig_size = (10, 6)

        mark_style = ['o', 's', '^', 'x', 'D', 'p', '>', '<']

        # shell 6 - 19
        inert_lim_low = 11
        inert_lim_high = 16

        num_pow = stats_data.vel_str_func.shape[-1]

        # Test the average of amplitudes is correct
        avg_amp = np.mean(np.absolute(run_data.u), axis=0)

        print(np.allclose(avg_amp, sys_msr_data.a_n_t_avg), avg_amp - sys_msr_data.a_n_t_avg)

        # Testing the computation of the time averaged amps
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.k, avg_amp, label="python")
        ax1.plot(sys_msr_data.k, sys_msr_data.a_n_t_avg, label="C")
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.legend()
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "Test_Avg_Amp.png", bbox_inches='tight')
        plt.close()

        # Testing fitting 
        popt, pcov = curve_fit(spec_func, sys_msr_data.k, sys_msr_data.a_n_t_avg, maxfev=50000)
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(sys_msr_data.k, sys_msr_data.a_n_t_avg, '.', label="C")
        ax1.plot(sys_msr_data.k, spec_func(sys_msr_data.k, *popt), label=r"$a k^{-b}e^{-c k} + d$;" + r"$ a = {}, b = {}, c = {}, d = {}$".format(np.around(popt[0], 6), np.around(popt[1], 6), np.around(popt[2], 6), np.around(popt[3], 6)))
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.legend()
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "Test_Avg_Amp.png", bbox_inches='tight')
        plt.close()

        # Energy Spectra
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.enrg_spec_t_avg), label="$E_n$", marker = mark_style[0], markevery = 1)
        poly_output = np.polyfit(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.enrg_spec_t_avg[inert_lim_low:inert_lim_high]), 1, full = True)
        pfit_info = poly_output[0]
        poly_resid = poly_output[1][0]
        pfit_slope = pfit_info[0]
        pfit_c = pfit_info[1]
        ax1.plot(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high])*pfit_slope + pfit_c, 'k--', label = "$k^-{}$".format(np.absolute(pfit_slope)))
        ax1.legend()
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Spectra.png", bbox_inches='tight')
        plt.close()

        # Amp
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.a_n_t_avg), label="$a_n$", marker=mark_style[0], markevery=1)
        poly_output = np.polyfit(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.a_n_t_avg[inert_lim_low:inert_lim_high]), 1, full = True)
        pfit_info = poly_output[0]
        poly_resid = poly_output[1][0]
        pfit_slope = pfit_info[0]
        pfit_c = pfit_info[1]
        ax1.plot(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high])*pfit_slope + pfit_c, 'k--', label = "$k^-{}$".format(np.absolute(pfit_slope)))
        ax1.legend()
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Amplitudes.png", bbox_inches='tight')
        plt.close()

        # Create average amp field
        u_avg_amp_all = sys_msr_data.a_n_t_avg * np.exp(1j * np.angle(run_data.u))
        u_avg_amp = avg_amp * np.exp(1j * np.angle(run_data.u))

        # Compute PDFs from the full solver data
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        for j, i in enumerate([1 - 1]):
            pdf, centres = compute_pdf(np.real(run_data.u[:, i]), nbins=500, normed=True)
            if i == -1:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
            else:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))    
        ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_AVGAMP + "VelReal_PDF_InOne.png")
        plt.close()

        # Compute the PDFs with the time averaged amplitudes
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        shells = [1 - 1, 5 - 1, 10 - 1, 15 - 1, 20 - 1]
        to_normalize = True
        for j, i in enumerate(shells):
            pdf, centres = compute_pdf(np.real(u_avg_amp[:, i]), nbins=500, normed=to_normalize)
            if i == -1:
            	p, = ax1.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
            else:
            	p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))
            if to_normalize:
                ax1.plot(centres, 1/ (np.pi * np.sqrt( sys_msr_data.a_n_t_avg[i]**2 / np.mean(np.real(u_avg_amp[:, i])**2) - centres**2)), marker='.', label="Analytical PDF", color = p.get_color())
                ax1.set_xlabel(r"$\Re u_n / \langle (\Re u_n)^2 \rangle^{1/2}$")
            else:
                ax1.plot(centres, 1/ (np.pi * np.sqrt( sys_msr_data.a_n_t_avg[i]**2 - centres**2)), marker='.', label="Analytical PDF", color = p.get_color())
                ax1.set_xlabel(r"$\Re u_n $")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_VelReal_PDF_InOne.png")
        plt.close()

        # Compute the PDF of the time varying amplitudes
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        shells = [1 - 1, 5 - 1, 10 - 1, 15 - 1, 20 - 1]
        to_normalize = True
        for j, i in enumerate(shells):
            pdf, centres = compute_pdf(np.absolute(run_data.u[:, i]), nbins=500, normed=to_normalize)
            if i == -1:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
            else:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))
        if to_normalize:
            ax1.set_xlabel(r"$ |u_n| / \langle ( |u_n|)^2 \rangle^{1/2}$")
        else:
            ax1.set_xlabel(r"$ |u_n|$")            
        ax1.set_ylabel(r"PDF")
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_VelAbsolute_PDF_InOne.png")
        plt.close()

        # Plot the PDFs of the cosine of the Phases
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        shells = [1 - 1, 5 - 1, 10 - 1, 15 - 1, 20 - 1]
        to_normalize = True
        for j, i in enumerate(shells):
            pdf, centres = compute_pdf(np.cos(np.angle((run_data.u[:, i]))), nbins=500, normed=to_normalize)
            if i == -1:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
            else:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))
        if to_normalize:   
            ax1.set_xlabel(r"$ \cos(\phi_n) / \langle ( \cos(\phi_n))^2 \rangle^{1/2}$") 
        else:
            ax1.set_xlabel(r"$ \cos(\phi_n) $")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_CosArgUn_PDF_InOne.png")
        plt.close()
        
        # Plot the PDFs of the Phases
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        shells = [1 - 1, 5 -1, 10 - 1, 15 - 1, 20 - 1]
        to_normalize = True
        for j, i in enumerate(shells):
            pdf, centres = compute_pdf(np.angle((run_data.u[:, i])), nbins=500, normed=to_normalize)
            if i == -1:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(sys_vars.N))    
            else:
                p, = ax1.plot(centres, pdf, label="$n = {}$".format(i + 1))
        if to_normalize:   
            ax1.set_xlabel(r"$ \phi_n  / \langle \phi_n^2 \rangle^{1/2}$") 
        else:
            ax1.set_xlabel(r"$ \phi_n $")
        ax1.set_ylabel(r"PDF")
        ax1.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        ax1.set_yscale('log')
        ax1.legend()
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Phases_PDF_InOne.png")
        plt.close()
        
        # Loop through time and compute quantities for testing
        avg_u_enrg_flux = np.zeros((sys_vars.N,))
        avg_u_hel_flux = np.zeros((sys_vars.N,))
        str_func_avg_u = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        str_func_avg_u_enrg_flux = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        str_func_avg_u_hel_flux = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        u_pad = np.zeros((sys_vars.N + 4, ), dtype="complex128")
        # avg_u_enrg_flux_all = np.zeros((sys_vars.N,))
        # avg_u_hel_flux_all = np.zeros((sys_vars.N,))
        # str_func_avg_u_all = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        # str_func_avg_u_enrg_flux_all = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        # str_func_avg_u_hel_flux_all = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        # u_pad_all = np.zeros((sys_vars.N + 4, ), dtype="complex128")

        num_t_data = sys_vars.ndata
        for t in range(num_t_data):
            print(t)

            # Get padded data
            u_pad[2:sys_vars.N + 2] = u_avg_amp[t, :]
            # u_pad_all[2:sys_vars.N + 2] = u_avg_amp_all[t, :]

            # Get fluxes
            tmp_enrg, tmp_hel = compute_u_flux(u_pad, sys_vars.N, sys_vars.EPS, sys_vars.Lambda)
            avg_u_enrg_flux += tmp_enrg
            avg_u_hel_flux += tmp_hel
            # tmp_enrg, tmp_hel = compute_u_flux(u_pad_all, sys_vars.N, sys_vars.EPS, sys_vars.Lambda)
            # avg_u_enrg_flux_all += tmp_enrg
            # avg_u_hel_flux_all += tmp_hel

            # Get Structure Functions
            tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_pad, sys_vars.N, stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
            str_func_avg_u += tmp_u
            str_func_avg_u_enrg_flux += tmp_enrg_flux
            str_func_avg_u_hel_flux += tmp_hel_flux
            # tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_pad_all, sys_vars.N, stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
            # str_func_avg_u_all += tmp_u
            # str_func_avg_u_enrg_flux_all += tmp_enrg_flux
            # str_func_avg_u_hel_flux_all += tmp_hel_flux

        # Average
        avg_u_enrg_flux          /= num_t_data
        avg_u_hel_flux           /= num_t_data
        str_func_avg_u           /= num_t_data
        str_func_avg_u_enrg_flux /= num_t_data
        str_func_avg_u_hel_flux  /= num_t_data
        # avg_u_enrg_flux_all          /= num_t_data
        # avg_u_hel_flux_all           /= num_t_data
        # str_func_avg_u_all           /= num_t_data
        # str_func_avg_u_enrg_flux_all /= num_t_data
        # str_func_avg_u_hel_flux_all  /= num_t_data

        # Plot Time Averaged Energy Flux Scaling
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        tmp1 = np.copy(sys_msr_data.enrg_flux_t_avg)
        tmp2 = np.copy(sys_msr_data.enrg_flux_t_avg)
        tmp1[tmp1 > 0] = 0.0
        tmp2[tmp2 < 0] = 0.0
        ax1.plot(sys_msr_data.k, np.absolute(tmp1), '.', label = "Neg - Solver Data")
        ax1.plot(sys_msr_data.k, tmp2, '.', label = "Pos - Solver Data")
        tmp1 = np.copy(avg_u_enrg_flux)
        tmp2 = np.copy(avg_u_enrg_flux)
        tmp1[tmp1 > 0] = 0.0
        tmp2[tmp2 < 0] = 0.0
        ax1.plot(sys_msr_data.k, np.absolute(tmp1), '*', label = "Neg - Averaged Amps")
        ax1.plot(sys_msr_data.k, tmp2, '*', label = "Pos - Averaged Amps")
        ax1.set_xlabel(r"$k_n$")
        ax1.set_ylabel(r"$\mathcal{E}_n$")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.legend()
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2 = fig.add_subplot(gs[0, 1])
        p, = ax2.plot(sys_msr_data.k, np.absolute(sys_msr_data.enrg_flux_t_avg), label = "Time Averaged Energy Flux")
        p1, = ax2.plot(sys_msr_data.k, np.absolute(avg_u_enrg_flux), marker = mark_style[1], label = "Averaged Amps")
        slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(avg_u_enrg_flux)), inert_lim_low, inert_lim_high)
        ax2.plot(sys_msr_data.k, sys_msr_data.k ** (-1), '--', label = "$k^{-1}$", color = p.get_color())
        ax2.plot(sys_msr_data.k, sys_msr_data.k ** (slope), '--', label = "$k^{}$".format(np.around(slope, 3)), color = p1.get_color())
        ax2.set_xlabel(r"$k_n$")
        ax2.set_ylabel(r"$\mathcal{E}_n$")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.legend()
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Time_Averaged_Energy_Flux.png", bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 1)
        ax2 = fig.add_subplot(gs[0, 0])
        y = sys_msr_data.k * np.absolute(sys_msr_data.enrg_flux_t_avg)
        yy = y[inert_lim_low:inert_lim_high] - np.mean(y[inert_lim_low:inert_lim_high])
        delta_n = yy[1:]
        delta_n1 = yy[:-1]
        p, = ax2.plot(sys_msr_data.k[inert_lim_low:inert_lim_high-1], delta_n - delta_n1, label = "Time Averaged Energy Flux")
        slope, c, resid = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(sys_msr_data.enrg_flux_t_avg)), inert_lim_low, inert_lim_high)
        # ax2.plot(sys_msr_data.k, sys_msr_data.k ** (slope), '--', label = "$k^{}$".format(np.around(slope, 3)), color = p1.get_color())
        print(slope, c, resid)
        ax2.set_xlabel(r"$k_n$")
        ax2.set_ylabel(r"$\mathcal{E}_n$")
        ax2.set_xscale("log")
        # ax2.set_yscale("log")
        ax2.legend()
        ax2.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Time_Averaged_kn_Energy_Flux.png", bbox_inches='tight')
        plt.close()

        # Plot Time Averaged Helicity Flux Scaling
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        tmp1 = np.copy(sys_msr_data.kin_hel_flux_t_avg)
        tmp2 = np.copy(sys_msr_data.kin_hel_flux_t_avg)
        tmp1[tmp1 > 0] = 0.0
        tmp2[tmp2 < 0] = 0.0
        ax1.plot(sys_msr_data.k, np.absolute(tmp1), '.', label = "Neg - Solver Data")
        ax1.plot(sys_msr_data.k, tmp2, '.', label = "Pos - Solver Data")
        tmp1 = np.copy(avg_u_hel_flux)
        tmp2 = np.copy(avg_u_hel_flux)
        tmp1[tmp1 > 0] = 0.0
        tmp2[tmp2 < 0] = 0.0
        ax1.plot(sys_msr_data.k, np.absolute(tmp1), '*', label = "Neg - Averaged Amps")
        ax1.plot(sys_msr_data.k, tmp2, '*', label = "Pos - Averaged Amps")
        ax1.set_xlabel(r"$k_n$")
        ax1.set_ylabel(r"$\mathcal{E}_n$")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.legend()
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        ax2 = fig.add_subplot(gs[0, 1])
        tmp1 = np.copy(avg_u_hel_flux)
        tmp2 = np.copy(avg_u_hel_flux)
        tmp1[tmp1 > 0] = 0.0
        tmp2[tmp2 < 0] = 0.0
        p1, = ax2.plot(sys_msr_data.k, np.absolute(sys_msr_data.kin_hel_flux_t_avg), label = "Time Averaged Helicity Flux")
        p, = ax2.plot(sys_msr_data.k, np.absolute(avg_u_hel_flux), label = "Avg Amps")
        slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(np.absolute(avg_u_hel_flux)), 1, 8)
        ax2.plot(sys_msr_data.k, sys_msr_data.k ** (-2), '--', label = "$k^{-2}$", color = p.get_color())
        ax2.plot(sys_msr_data.k, sys_msr_data.k ** (slope), '--', label = "$k^{}$".format(np.around(slope, 3)), color = p1.get_color())
        ax2.set_xlabel(r"$k_n$")
        ax2.set_ylabel(r"$\mathcal{E}_n$")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.legend()
        ax2.grid(which="both", axis="both", color='k', linestyle=":", linewidth=0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Time_Averaged_Helicity_Flux.png", bbox_inches='tight')
        plt.close()

        # Plot Str Funcs - Velocity
        fig = plt.figure(figsize=fig_size)
        gs = GridSpec(2, 3)
        for i in range(2):
            for j in range(3):
                ax1 = fig.add_subplot(gs[i, j])
                p = i * 3 + j + 1
                ax1.plot(sys_msr_data.k, stats_data.vel_str_func[:, p - 1] / sys_vars.ndata, label = "Full Field - p = {}".format(p))
                ax1.plot(sys_msr_data.k, str_func_avg_u[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
                # ax1.plot(sys_msr_data.k, str_func_avg_u_all[:, p - 1], marker = mark_style[2], label = "All Averaged Amps - p = {}".format(p))
                ax1.set_xlabel(r"$k_n$")
                ax1.set_ylabel(r"$\ln \mathcal{S}_p (|u_n|)$")
                ax1.set_xscale("log")
                ax1.set_yscale("log")
                ax1.set_title("{}nd Order Velocity Str Func".format(p))
                ax1.legend()
                ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Velocity.png", bbox_inches='tight')
        plt.close()
    	

        # # delta u_n comparison
        # delta_u_n_abs = np.absolute(run_data.u) - sys_msr_data.a_n_t_avg
        # delta_u_n_abs_sqrd = delta_u_n_abs**2
        # delta_u_n_abs_sqrd_avg = np.mean(delta_u_n_abs_sqrd, axis=0)
        # delta_u_n_abs_sqrd_from_str_func = np.absolute(stats_data.vel_str_func[:, 1] / sys_vars.ndata - str_func_avg_u[:, 1])
        # fig = plt.figure(figsize=fig_size)
        # gs = GridSpec(1, 1)
        # ax1 = fig.add_subplot(gs[0, 0])
        # ax1.plot(sys_msr_data.k, delta_u_n_abs_sqrd_avg, marker='*', label=r"From Data")
        # ax1.plot(sys_msr_data.k, delta_u_n_abs_sqrd_from_str_func, label=r"From Struc Func; p = 2")
        # slope, c, _ = slope_fit(np.log(sys_msr_data.k), np.log(delta_u_n_abs_sqrd_avg), inert_lim_low, inert_lim_high)
        # ax2.plot(sys_msr_data.k, np.absolute(avg_u_enrg_flux_all), marker = mark_style[2], label = "Slope; $p = {}$".format(np.around(slope)))
        # ax1.set_xlabel(r"$k_n$")
        # ax1.set_ylabel(r"$\ln \mathcal{S}_p (|u_n|)$")
        # ax1.set_xscale("log")
        # ax1.set_yscale("log")
        # ax1.set_title(r"$\langle (\delta u_n)^2\rangle$ vs $\langle |u_n|^2\rangle - \langle |u_n|\rangle^2$")
        # ax1.legend()
        # ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        # plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_DeltaUn_StrFunc_Velocity.png", bbox_inches='tight')
        # plt.close()

        # # Plot the difference between the str functions and the time aveeraged amplitude ones
        # fig = plt.figure(figsize=fig_size)
        # gs = GridSpec(2, 3)
        # for i in range(2):
        #     for j in range(3):
        #         ax1 = fig.add_subplot(gs[i, j])
        #         p = i * 3 + j + 1
        #         ax1.plot(sys_msr_data.k, np.absolute(stats_data.vel_str_func[:, p - 1] / sys_vars.ndata - str_func_avg_u[:, p - 1]), label = "Differnce - p = {}".format(p))
        #         if p >= 2:
        #             approx_2nd_order = (p * (p - 1) / 2) * sys_msr_data.a_n_t_avg**(p - 2) * delta_u_n_abs_sqrd_from_str_func
        #             approx_3nd_order = approx_2nd_order + (p * (p -1) * (p -2)/ 3) * sys_msr_data.a_n_t_avg**(p - 3) * np.mean(delta_u_n_abs**3, axis=0)
        #             approx_4nd_order = approx_3nd_order + (p * (p -1) * (p -2) * (p - 3)/ 24) * sys_msr_data.a_n_t_avg**(p - 4) * np.mean(delta_u_n_abs**4, axis=0)
        #             ax1.plot(sys_msr_data.k, approx_2nd_order, label=r"Approx for Diff 2nd order".format(p))
        #             ax1.plot(sys_msr_data.k, approx_3nd_order, label=r"Approx for Diff 3rd order".format(p))
        #             ax1.plot(sys_msr_data.k, approx_4nd_order, label=r"Approx for Diff 4rd order".format(p))
        #         ax1.set_xlabel(r"$k_n$")
        #         ax1.set_ylabel(r"$\ln \mathcal{S}_p (|u_n|)$")
        #         ax1.set_xscale("log")
        #         ax1.set_yscale("log")
        #         ax1.set_title("{}nd Order Velocity Str Func".format(p))
        #         ax1.legend()
        #         ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        # plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_Difference_StrFunc_Velocity.png", bbox_inches='tight')
        # plt.close()

        # Plot Str Funcs - Energy Flux
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(2, 3)
        for i in range(2):
            for j in range(3):
                ax1 = fig.add_subplot(gs[i, j])
                p = i * 3 + j + 1
                ax1.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 0] / sys_vars.ndata, label = "Full Field - p = {}".format(p))
                ax1.plot(sys_msr_data.k, str_func_avg_u_enrg_flux[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
                # ax1.plot(sys_msr_data.k, str_func_avg_u_enrg_flux_all[:, p - 1], marker = mark_style[2], label = "All Averaged Amps - p = {}".format(p))
                ax1.set_xlabel(r"$k_n$")
                ax1.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
                ax1.set_xscale("log")
                ax1.set_yscale("log")
                ax1.set_title("{}nd Order Energy Flux Str Func".format(p))
                ax1.legend()
                ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Energy_Flux.png", bbox_inches='tight')
        plt.close()

        # Plot the difference in Str Funcs - Energy Flux
        fig = plt.figure(figsize = fig_size)
        gs  = GridSpec(2, 3)
        for i in range(2):
            for j in range(3):
                ax1 = fig.add_subplot(gs[i, j])
                p = i * 3 + j + 1
                ax1.plot(sys_msr_data.k, np.absolute(stats_data.vel_flux_str_func_abs[:, p - 1, 0] / sys_vars.ndata - str_func_avg_u_enrg_flux[:, p - 1]), label = "Difference - p = {}".format(p))
                # flux_approx_1nd_order = 
                # ax1.plot(sys_msr_data.k, flux_approx_1st_order, marker = mark_style[1], label = "Approx - p = {}".format(p))
                ax1.set_xlabel(r"$k_n$")
                ax1.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
                ax1.set_xscale("log")
                ax1.set_yscale("log")
                ax1.set_title("{}nd Order Energy Flux Str Func".format(p))
                ax1.legend()
                ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Energy_Flux.png", bbox_inches='tight')
        plt.close()

        # Plot Str Funcs - Helicity Flux
        fig = plt.figure(figsize = fig_size)
        gs = GridSpec(2, 3)
        for i in range(2):
            for j in range(3):
                ax1 = fig.add_subplot(gs[i, j])
                p = i * 3 + j + 1
                ax1.plot(sys_msr_data.k, stats_data.vel_flux_str_func_abs[:, p - 1, 1] / sys_vars.ndata, label = "Full Field - p = {}".format(p))
                ax1.plot(sys_msr_data.k, str_func_avg_u_hel_flux[:, p - 1], marker = mark_style[1], label = "Averaged Amps - p = {}".format(p))
                # ax1.plot(sys_msr_data.k, str_func_avg_u_hel_flux_all[:, p - 1], marker = mark_style[2], label = "All Averaged Amps - p = {}".format(p))
                ax1.set_xlabel(r"$k_n$")
                ax1.set_ylabel(r"$\ln \mathcal{S}_p (\Pi)$")
                ax1.set_xscale("log")
                ax1.set_yscale("log")
                ax1.set_title("{}nd Order Helicity Flux Str Func".format(p))
                ax1.legend()
                ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_StrFunc_Helicity_Flux.png", bbox_inches='tight')
        plt.close()

        # Get adjusted time averaged amplitudes
        a_n = sys_msr_data.a_n_t_avg

        fig = plt.figure(figsize = fig_size)
        gs = GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        a_n_slope, a_n_c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(a_n), inert_lim_low, inert_lim_high)
        ax1.plot(np.log2(sys_msr_data.k), np.log2(sys_msr_data.a_n_t_avg), label = "$a_n \sim k^{}$".format(np.around(a_n_slope, 3)), marker = None, markevery = 1)
        for beta in [0., 0.3, 0.6, 0.9, 1.2, 1.5]:
            a_n_adjust = a_n * (sys_msr_data.k ** (np.absolute(a_n_slope) - beta))
            p, = ax1.plot(np.log2(sys_msr_data.k), np.log2(a_n_adjust), marker = mark_style[1], markevery = 1)
            slope, c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(a_n_adjust), inert_lim_low, inert_lim_high)
            ax1.plot(np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high]), np.log2(sys_msr_data.k[inert_lim_low:inert_lim_high])*slope + c, '--', label = r"$k^a \quad a = $ {:0.3f}".format(slope), color = p.get_color())
        ax1.legend()
        ax1.grid(which = "both", axis = "both", color = 'k', linestyle = ":", linewidth = 0.5)
        plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_All_Slopes.png", bbox_inches='tight')
        plt.close()

        # # Get adjusted time averaged amplitudes
        # a_n          = sys_msr_data.a_n_t_avg
        # num_t_data   = sys_vars.ndata
        # phases       = np.angle(run_data.u)
        # alpha_slopes = np.arange(0.0, 1.51, 0.01)

        # # alocate memory
        # u_pad = np.zeros((sys_vars.N + 4, ), dtype = "complex128")
        # enrg_flux_slope          = np.zeros((len(alpha_slopes),))
        # hel_flux_slope           = np.zeros((len(alpha_slopes),))
        # str_func_u_slope         = np.zeros((len(alpha_slopes), num_pow))
        # str_func_enrg_flux_slope = np.zeros((len(alpha_slopes), num_pow))
        # str_func_hel_flux_slope  = np.zeros((len(alpha_slopes), num_pow))

        # for i, alpha in enumerate(alpha_slopes):
        #     print(i)

        #     # Alocate memory
        #     avg_u_enrg_flux          = np.zeros((sys_vars.N,))
        #     avg_u_hel_flux           = np.zeros((sys_vars.N,))
        #     str_func_avg_u           = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        #     str_func_avg_u_enrg_flux = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))
        #     str_func_avg_u_hel_flux  = np.zeros((sys_vars.N, stats_data.vel_str_func.shape[-1]))

        #     # Create average amp field
        #     a_n_adjust = a_n * (sys_msr_data.k ** (np.absolute(a_n_slope) - alpha))
        #     u_avg_amp  = a_n_adjust * np.exp(1j * phases)

        #     # Loop through time and compute quantities for testing
        #     for t in range(num_t_data):

        #         # Get padded data
        #         u_pad[2    :sys_vars.N + 2] = u_avg_amp[t, :]

        #         # Get fluxes
        #         tmp_enrg, tmp_hel = compute_u_flux(u_pad, sys_vars.N, sys_vars.EPS, sys_vars.Lambda)
        #         avg_u_enrg_flux   += tmp_enrg
        #         avg_u_hel_flux    += tmp_hel

        #         # Get Structure Functions
        #         tmp_u, tmp_enrg_flux, tmp_hel_flux = compute_str_func(u_pad, sys_vars.N, stats_data.vel_str_func.shape[-1], sys_vars.EPS, sys_vars.Lambda)
        #         str_func_avg_u                     += tmp_u
        #         str_func_avg_u_enrg_flux           += tmp_enrg_flux
        #         str_func_avg_u_hel_flux            += tmp_hel_flux

        #     # Average
        #     avg_u_enrg_flux          /= num_t_data
        #     avg_u_hel_flux           /= num_t_data
        #     str_func_avg_u           /= num_t_data
        #     str_func_avg_u_enrg_flux /= num_t_data
        #     str_func_avg_u_hel_flux  /= num_t_data

        #     # Get scaling slopes
        #     enrg_flux_slope[i], c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(np.absolute(avg_u_enrg_flux)), inert_lim_low, inert_lim_high)
        #     hel_flux_slope[i], c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(np.absolute(avg_u_hel_flux)), inert_lim_low, inert_lim_high)
        #     for p in range(num_pow):
        #         str_func_enrg_flux_slope[i, p], c, _ = slope_fit(np.log2(sys_msr_data.k), np.log2(str_func_avg_u_hel_flux[:, p]), inert_lim_low, inert_lim_high)
        #         str_func_hel_flux_slope[i, p], c, _  = slope_fit(np.log2(sys_msr_data.k), np.log2(str_func_avg_u_enrg_flux[:, p]), inert_lim_low, inert_lim_high)
        #         str_func_u_slope[i, p], c, _         = slope_fit(np.log2(sys_msr_data.k), np.log2(str_func_avg_u[:, p]), inert_lim_low, inert_lim_high)


        # fig = plt.figure(figsize = fig_size)
        # gs  = GridSpec(1, 1)
        # ax1 = fig.add_subplot(gs[0, 0])
        # ax1.plot(alpha_slopes, enrg_flux_slope, label = "Energy Flux")
        # ax1.plot(alpha_slopes, hel_flux_slope, label = "Helicity Flux")
        # ax1.set_xlabel(r"$\alpha$")
        # ax1.set_ylabel(r"Slopes")
        # ax1.legend()
        # # ax1.set_xlim(alpha_slopes[0], alpha_slopes[-1])
        # plt.savefig(cmdargs.out_dir_AVGAMP + "AvgAmps_EnergyFlux_Slopes_vs_Alpha.png", bbox_inches='tight')
        # plt.close()