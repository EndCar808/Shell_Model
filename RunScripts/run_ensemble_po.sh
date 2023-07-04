#!/bin/bash

N=25

dt=0.0001
T=1000000.000
trans_frac=0.1

nu=5e-07

forcing="NONE"
force_k=1
force_scale=0.100

u0="PO_PLAWEXP_RND"
data_dir="./Data/Thesis/EnsembleRuns/POStatsEnsemble/"

save_every=10000
stats_every=25000

# ############################## Run Solver
# for tag in {0..20}
# do
# 	solv_cmd="Solver/bin/solver_phase_only -o $data_dir -n $N -s 0.00000 -e $T -T 1 -T $trans_frac -c 0 -c 0.900000 -h $dt -h 0 -a 0.3333333333333333 -b 0.3333333333333333 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v $nu -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i $u0 -t $tag -f $forcing -f $force_k -f $force_scale -p $save_every -p $stats_every -p 10000 -z NONE"
# 	# echo $solv_cmd
# 	$solv_cmd &
# done

# wait


# ############################## Run Plotting
# for tag in {0..15}
# do
# 	plot_cmd="python3 Plotting/plot_hd.py -i "$data_dir"HD-INTFACRK4-FULL_N[$N]_T[0.0,"$dt","$T"]_SMP["$save_every","$stats_every",10000]_NU["$nu"]_ALPHA[1.500]_K[0.050,2.000]_EPS[0.50]_FORC["$forcing","$force_k","$force_scale"]_u0[$u0]_TAG[$tag]/ --plot"
# 	echo $plot_cmd
# 	$plot_cmd &
# done

# wait

############################## Run Ensemble Plotting Script
ens_cmd="python3 Plotting/plot_ensemble_stats.py /home/enda/PhD/Shell_Model/Data/Thesis/Plots/POStats/ POModel "
for tag in {0..20}
do
	input_data_dir=""$data_dir"HD-RK4-PO_N[$N]_T[0.0,"$dt","$T"]_SMP["$save_every","$stats_every",10000]_NU["$nu"]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC["$forcing","$force_k","$force_scale"]_u0[$u0]_TAG[$tag] "
	ens_cmd+=$input_data_dir
done
echo $ens_cmd
$ens_cmd

