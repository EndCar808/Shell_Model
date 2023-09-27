#!/bin/bash

dt=5e-05
t0=200.0
T=2000.000

save_every=10000
stats_every=100000

# for batch in {0..12} 
# do
# 	for id in {1..50} ## batchsize
# 	do
# 		dense_data_cmd="Solver/bin/solver_hydro_dense -o ./Data/DenseData/ -n 25 -s 0.00000 -e $T -T 1 -T 0.1 -c 0 -c 0.900000 -h $dt -h 0 -a 1.5000000000 -b 1.5000000000 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 5e-07 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i N_SCALING_RAND -t "$id" -f DELTA -f 1 -f 0.100 -p 1 -p 500000 -p 500000 -z NONE -z NONE -z 0"
# 		# echo $dense_data_cmd
# 		$dense_data_cmd &
# 	done

# 	wait	

# 	for r in 1 10 100 1000 10000 100000 1000000 10000000
# 	do
# 		for id in {1..50} ## batchsize
# 		do
# 		batch_size=50
# 		batch_id=`expr $batch \* $batch_size + $id`

# 		solver_repl_cmd="Solver/bin/solver_amp_only_fxd_phase -o ./Data/AOReplaceData/ -n 25 -s $t0 -e $T -T 0 -T 0.1 -c 0 -c 0.900000 -h $dt -h 0 -a 0.3333000000 -b 0.3333000000 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 5e-07 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i AO_INPUT_PHASE_REPLACE -t $batch_id-$r -f DELTA -f 1 -f 0.100 -p "$save_every" -p "$stats_every" -p "$r" -z ./Data/DenseData/HD-INTFACRK4-FULL_N[25]_T[0.0,"$dt","$T"]_SMP[1,500000,500000]_NU[5e-07]_ALPHA[1.500]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[N_SCALING_RAND]_TAG["$id"]/Main_HDF_Data.h5 -z NONE -z 0"
# 		# echo $solver_repl_cmd
# 		$solver_repl_cmd &
# 		# echo $batch_id $r
# 		done

# 		wait
# 	done

# 	wait
# done

for r in 1 10 100 1000 10000 100000 1000000 10000000
	do
	
	plot_ens_cmd="python3 Plotting/plot_ensemble_stats.py /home/enda/PhD/Shell_Model/Data/Thesis/Plots/ReplaceData/FullData/Repl-"$r"/ EnsREPL-"$r" "
	for batch in {0..12} 
	do
		for id in {1..50} ## batchsize
		do
			batch_size=50
			batch_id=`expr $batch \* $batch_size + $id`

			plot_ens_cmd+="/home/enda/PhD/Shell_Model/Data/AOReplaceData/HD-INTFACRK4-AOFXDPHASE_N[25]_T["$t0","$dt","$T"]_SMP["$save_every","$stats_every","$r"]_NU[5e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[AO_INPUT_PHASE_REPLACE]_TAG["$batch_id"-"$r"] "
		done
	done
	echo $plot_ens_cmd
	$plot_ens_cmd &
done

wait