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
from Plotting.functions import tc
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

        def __init__(self, init_file = None, cmd_only = False, print_updates = False):
            self.init_file = init_file
            self.cmd_only  = cmd_only
            self.print     = print_updates
            
    ## Initialize class
    cargs = cmd_args()

    # print(getopt.getopt(argv, "i:c:", ["cmdonly"]))
    try:
        ## Gather command line arguments
        opts, args = getopt.getopt(argv, "i:c:", ["cmdonly", "print"])
    except Exception as e:
        print("[" + tc.R + "ERROR" + tc.Rst + "] ---> Incorrect Command Line Arguements.")
        print(e)
        sys.exit()

    ## Parse command line args
    for opt, arg in opts:

        if opt in ['-i']:
            ## Read in config file
            cargs.init_file = str(arg)
            print("\nInput configuration file: " + tc.C + cargs.init_file + tc.Rst)

            if not os.path.isfile(cargs.init_file):
                print("[" + tc.R + "ERROR" + tc.Rst + "] ---> File Does not exist, double check input file path.")
                sys.exit()

        if opt in ['--cmdonly']:
            ## Read in indicator to print out commands to terminal only
            cargs.cmd_only = True

        if opt in ['--print']:
            ## Read in printing indicator
            cargs.print = True

    return cargs

######################
##       MAIN       ##
######################
if __name__ == '__main__':
    ##########################
    ##  PARSE COMMAND LINE  ##
    ##########################
    cmdargs = parse_cml(sys.argv[1:])


    ##########################
    ##  DEFAULT PARAMETERS  ##
    ##########################
    ## Space variables
    N = 19

    ## System parameters
    nu             = 1e-7
    eta            = 1e-7
    hypervisc      = False
    hypervisc_pow  = 2.0 
    ekmn_hypo_diff = False
    ekmn_hypo_pow  = -2.0
    eps            = 0.5
    esp_m          = 1.0/3.0
    k0             = 1.0
    Lambda         = 2.0
    alpha          = 1.4
    beta           = 1.4

    ## Time parameters
    t0          = 0.0
    T           = 12.0
    dt          = 1e-3
    step_type   = False
    cfl_cond    = False
    trans_iters = False
    cfl         = 0.9
    
    ## Solver parameters
    ic          = "N_SCALING"
    forcing     = "NONE"
    force_k     = 0
    force_scale = 1.0
    save_every  = 1

    ## Directory/File parameters
    input_dir       = "NONE"
    output_dir      = "./Data/Tmp/"
    file_only_mode  = False
    solver_tag      = "Test"
    system_tag      = "[MAG_HYDRO_INTFACRK4_FULL]"
    post_input_dir  = output_dir
    post_output_dir = output_dir

    ## Job parameters
    executable                  = "Solver/bin/solver"
    plot_options                = ""
    post_options                = ""
    plotting                    = False
    solver                      = False
    postprocessing              = False
    collect_data                = False
    solver_procs                = 1
    num_solver_job_threads      = 1
    num_postprocess_job_threads = 1
    num_plotting_job_threads    = 1

    #########################
    ##  PARSE CONFIG FILE  ##
    #########################
    ## Create parser instance
    parser = ConfigParser()

    ## Read in config file
    parser.read(cmdargs.init_file)

    ## Create list objects
    N          = []
    nu         = []
    eta        = []
    ic         = []
    T          = []
    dt         = []
    solver_tag = []
    alpha      = []
    beta       = []
    eps        = []
    eps_m      = []
    cfl        = []
    save_every = []
    input_file = []

    ## Parse input parameters
    for section in parser.sections():
        if section in ['SYSTEM']:
            if 'n' in parser[section]:
                for n in parser[section]['n'].lstrip('[').rstrip(']').split(', '):
                    N.append(int(n))
            if 'viscosity' in parser[section]:
                for n in parser[section]['viscosity'].lstrip('[').rstrip(']').split(', '):
                    nu.append(float(n))
            if 'diffusivity' in parser[section]:
                for n in parser[section]['diffusivity'].lstrip('[').rstrip(']').split(', '):
                    eta.append(float(n))
            if 'hyperviscosity' in parser[section]:
                hypervisc = int(parser[section]['hyperviscosity'] == 'True')
            if 'hypo_diffusion' in parser[section]:
                ekmn_hypo_diff = int(parser[section]['hypo_diffusion'] == 'True')
            if 'hyperviscosity_pow' in parser[section]:
                hypervisc_pow = float(parser[section]['hyperviscosity_pow'])
            if 'hypo_diffusion_pow' in parser[section]:
                ekmn_hypo_pow = float(parser[section]['hypo_diffusion_pow'])
            if 'vel_spect_slope' in parser[section]:
                for n in parser[section]['vel_spect_slope'].lstrip('[').rstrip(']').split(', '):
                    alpha.append(float(n))
            if 'mag_spect_slope' in parser[section]:
                for n in parser[section]['mag_spect_slope'].lstrip('[').rstrip(']').split(', '):
                    beta.append(float(n))
            if 'vel_interact_coeff' in parser[section]:
                for n in parser[section]['vel_interact_coeff'].lstrip('[').rstrip(']').split(', '):
                    eps.append(float(n))
            if 'mag_interact_coeff' in parser[section]:
                for n in parser[section]['mag_interact_coeff'].lstrip('[').rstrip(']').split(', '):
                    eps_m.append(float(n))
            if 'shell_k_prefac' in parser[section]:
                k0 = float(parser[section]['shell_k_prefac'])
            if 'shell_k_lambda' in parser[section]:
                lam = float(parser[section]['shell_k_lambda'])
        if section in ['SOLVER']:
            if 'initial_condition' in parser[section]:
                for n in parser[section]['initial_condition'].lstrip('[').rstrip(']').split(', '):
                    ic.append(str(n.lstrip('"').rstrip('"')))
            if 'forcing' in parser[section]:
                forcing = str(parser[section]['forcing'])
            if 'forcing_wavenumber' in parser[section]:
                force_k = int(float(parser[section]['forcing_wavenumber']))
            if 'forcing_scale' in parser[section]:
                force_scale = float(parser[section]['forcing_scale'])
            if 'save_data_every' in parser[section]:
                for n in parser[section]['save_data_every'].lstrip('[').rstrip(']').split(', '):
                    save_every.append(int(n))
            if 'stats_data_every' in parser[section]:
                stats_data_every = int(parser[section]['stats_data_every'])
            if 'replace_data_every' in parser[section]:
                replace_data_every = int(parser[section]['replace_data_every'])
        if section in ['TIME']:
            if 'end_time' in parser[section]:
                for n in parser[section]['end_time'].lstrip('[').rstrip(']').split(', '):
                    T.append(float(n))
            if 'timestep' in parser[section]:
                for n in parser[section]['timestep'].lstrip('[').rstrip(']').split(', '):
                    dt.append(float(n))
            if 'cfl' in parser[section]:
                for n in parser[section]['cfl'].lstrip('[').rstrip(']').split(', '):
                    cfl.append(float(parser[section]['cfl']))
            if 'start_time' in parser[section]:
                t0 = float(parser[section]['start_time'])
            if 'cfl_cond' in parser[section]:
                cfl_cond = int(parser[section]['cfl_cond'] == 'True')
            if 'trans_iters' in parser[section]:
                trans_iters = int(parser[section]['trans_iters'] == 'True')
            if 'trans_iters_frac' in parser[section]:
                trans_iters_frac = float(parser[section]['trans_iters_frac'])
            if 'adaptive_step_type' in parser[section]:
                step_type = int(parser[section]['adaptive_step_type'] == 'True')
        if section in ['DIRECTORIES']:
            if 'solver_input_dir' in parser[section]:
                for n in parser[section]['solver_input_dir'].lstrip('[').rstrip(']').split(', '):
                    input_file.append(str(n))
            if 'solver_input_str' in parser[section]:
                solver_input_str = str(parser[section]['solver_input_str'])
            if 'solver_input_param' in parser[section]:
                solver_input_param = str(parser[section]['solver_input_param'])
            if 'solver_output_dir' in parser[section]:
                output_dir = str(parser[section]['solver_output_dir'])
            if 'solver_tag' in parser[section]:
                for n in parser[section]['solver_tag'].lstrip('[').rstrip(']').split(', '):
                    solver_tag.append(n)
            if 'post_input_dir' in parser[section]:
                post_input_dir = str(parser[section]['post_input_dir'])
            if 'post_output_dir' in parser[section]:
                post_output_dir = str(parser[section]['post_output_dir'])
            if 'solver_file_only_mode' in parser[section]:
                file_only_mode = bool(utils.strtobool(parser[section]['solver_file_only_mode']))
            if 'system_tag' in parser[section]:
                system_tag = str(parser[section]['system_tag'])
        if section in ['JOB']:
            if 'executable' in parser[section]:
                executable = str(parser[section]['executable'])
            if 'plotting' in parser[section]:
                plotting = str(parser[section]['plotting'])
            if 'plot_script' in parser[section]:
                plot_script = str(parser[section]['plot_script'])
            if 'plot_options' in parser[section]:
                plot_options = str(parser[section]['plot_options'])
            if 'post_options' in parser[section]:
                post_options = str(parser[section]['post_options'])
            if 'call_solver' in parser[section]:
                solver = bool(utils.strtobool(parser[section]['call_solver']))
            if 'call_plotting' in parser[section]:
                plotting = bool(utils.strtobool(parser[section]['call_plotting']))
            if 'call_postprocessing' in parser[section]:
                postprocessing = bool(utils.strtobool(parser[section]['call_postprocessing']))
            if 'solver_procs' in parser[section]:
                solver_procs = int(parser[section]['solver_procs'])
            if 'collect_data' in parser[section]:
                collect_data = bool(utils.strtobool(parser[section]['collect_data']))
            if 'num_solver_job_threads' in parser[section]:
                num_solver_job_threads = int(parser[section]['num_solver_job_threads'])
            if 'num_postprocess_job_threads' in parser[section]:
                num_postprocess_job_threads = int(parser[section]['num_postprocess_job_threads'])
            if 'num_plotting_job_threads' in parser[section]:
                num_plotting_job_threads = int(parser[section]['num_plotting_job_threads'])

    ## Get the path to the runs output directory

    par_runs_output_dir = os.path.split(output_dir)[0]
    print(par_runs_output_dir)
    par_runs_output_dir += '/ParallelRunsDump/'
    if os.path.isdir(par_runs_output_dir) != True:
        print("Making folder:" + tc.C + " ParallelRunsDump/" + tc.Rst)
        os.mkdir(par_runs_output_dir)

    #########################
    ##      RUN SOLVER     ##
    #########################
    if solver:
        

        ## Get the number of processes to launch
        proc_limit = num_solver_job_threads
        print("\nNumber of Solver Processes Created = [" + tc.C + "{}".format(proc_limit) + tc.Rst + "]\n")

        # Create output objects to store process error and output
        if collect_data:
            solver_output = []
            solver_error  = []

        ## Generate command list 
        cmd_list = [["{} -o {} -n {} -s {:3.5f} -e {:3.5f} -T {} -T {} -c {} -c {:1.6f} -h {:1.16f} -h {} -a {:1.10f} -b {:1.10f} -w {:1.3f} -w {:1.3f} -y {:1.16f} -y {:1.16f} -v {:g} -v {} -v {:1.1f} -d {:g} -d {} -d {:1.1f} -i {} -t {} -f {} -f {} -f {:1.3f} -p {} -p {} -p {} -z {} -z {} -z {}".format(
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
                                                                                                                                                    et, ekmn_hypo_diff, ekmn_hypo_pow,
                                                                                                                                                    u0, 
                                                                                                                                                    s_tag, 
                                                                                                                                                    forcing, force_k, force_scale, 
                                                                                                                                                    save, stats_data_every, replace_data_every,
                                                                                                                                                    in_file, solver_input_str, solver_input_param)] for n in N for t in T for h, save in zip(dt, save_every) for u0 in ic for v in nu for et in eta for ep in eps for a in alpha for b in beta for ep_m in eps_m for c in cfl for s_tag in solver_tag for in_file in input_file]

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
                    print("Output:\n\n\t")
                    print(run_CodeOutput)
                    print("Error:\n\n\t")
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

    ###########################
    ##      RUN PLOTTING     ##
    ###########################
    if plotting:
    
        ## Get the number of processes to launch
        proc_limit = num_plotting_job_threads
        print("\nNumber of Post Processing Processes Created = [" + tc.C + "{}".format(proc_limit) + tc.Rst + "]\n")

        # Create output objects to store process error and output
        if collect_data:
            plot_output = []
            plot_error  = []
            
        
        ## Generate command list
        if "_mag_hydro" in executable or "_elsassar_mhd" in executable:
            cmd_list = [["python3 {} -i {} {}".format(
                                                plot_script, 
                                                post_input_dir + "_N[{}]_T[{:1.1f},{:g},{:1.3f}]_NU[{:g}]_ETA[{:g}]_ALPHA[{:1.3f}]_BETA[{:1.3f}]_K[{:1.3f},{:1.3f}]_EPS[{:1.2f},{:1.2f}]_FORC[{},{},{:1.3f}]_u0[{}]_TAG[{}]/".format(n, t0, h, t, v, e, a, b, k0, lam, ep, ep_m, forcing, force_k, force_scale, u0, s_tag), 
                                                plot_options)] for n in N for h in dt for t in T for e in eta for v in nu for a in alpha for b in beta for ep in eps for ep_m in eps_m for u0 in ic for s_tag in solver_tag]
        else:
            cmd_list = [["python3 {} -i {} {}".format(
                                                plot_script, 
                                                post_input_dir + "_N[{}]_T[{:1.1f},{:g},{:1.3f}]_NU[{:g}]_ALPHA[{:1.3f}]_K[{:1.3f},{:1.3f}]_EPS[{:1.2f}]_FORC[{},{},{:1.3f}]_u0[{}]_TAG[{}]/".format(n, t0, h, t, v, a, k0, lam, ep, forcing, force_k, force_scale, u0, s_tag), 
                                                plot_options)] for n in N for h in dt for t in T for v in nu for a in alpha for ep in eps for u0 in ic for s_tag in solver_tag]

        if cmdargs.cmd_only:
            print(tc.C + "\nPlotting Commands:\n" + tc.Rst)
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

                    ## Print output to terminal as it comes
                    for line in proc.stdout:
                        sys.stdout.write(line)
                    
                    # Communicate with process to retrive output and error
                    [run_CodeOutput, run_CodeErr] = proc.communicate()

                    # Append to output and error objects
                    if collect_data:
                        plot_output.append(run_CodeOutput)
                        plot_error.append(run_CodeErr)
                    
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
                with open(par_runs_output_dir + "par_run_plot_output_{}_{}.txt".format(cmdargs.init_file.lstrip('InitFiles/').rstrip(".ini"), d_t), "w") as file:
                    for item in plot_output:
                        file.write("%s\n" % item)

                # Write error to file
                with open(par_runs_output_dir + "par_run_plot_error_{}_{}.txt".format(cmdargs.init_file.lstrip('InitFiles/').rstrip(".ini"), d_t), "w") as file:
                    for i, item in enumerate(plot_error):
                        file.write("%s\n" % cmd_list[i])
                        file.write("%s\n" % item)
                    