#!/bin/bash

dt=5e-05
t0=500.0
T=5000.000

save_every=5000
stats_every=50000

for batch in {0..1} 
do
	for id in {1..50} ## batchsize
	do
		dense_data_cmd="Solver/bin/solver_phase_only_dense -o ./Data/DenseData_PO/ -n 25 -s 0.00000 -e $T -T 1 -T 0.1 -c 0 -c 0.900000 -h $dt -h 0 -a 1.5000000000 -b 1.5000000000 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 5e-07 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i PO_PLAWEXP_RND -t "$id" -f NONE -f 0 -f 0.100 -p 1 -p 500000 -p 500000 -z NONE -z NONE -z 0"
		# echo $dense_data_cmd
		$dense_data_cmd &
	done

	wait	

	for r in 1 10 100 1000 10000 100000 1000000 10000000
	do
		for id in {1..50} ## batchsize
		do
		batch_size=50
		batch_id=`expr $batch \* $batch_size + $id`

		solver_repl_cmd="Solver/bin/solver_amp_only_fxd_phase -o ./Data/AOReplaceData_PO/ -n 25 -s $t0 -e $T -T 0 -T 0.2 -c 0 -c 0.900000 -h $dt -h 0 -a 0.3333000000 -b 0.3333000000 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 5e-07 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i AO_INPUT_PHASE_REPLACE -t $batch_id-$r -f NONE -f 0 -f 0.100 -p "$save_every" -p "$stats_every" -p "$r" -z ./Data/DenseData_PO/HD-RK4-PO_N[25]_T[0.0,"$dt","$T"]_SMP[1,500000,500000]_NU[5e-07]_ALPHA[1.500]_K[0.050,2.000]_EPS[0.50]_FORC[NONE,0,0.100]_u0[PO_PLAWEXP_RND]_TAG["$id"]/Main_HDF_Data.h5 -z NONE -z 0"
		# echo $solver_repl_cmd
		$solver_repl_cmd &
		# echo $batch_id $r
		done

		wait
	done

	wait
done

for r in 1 10 100 1000 10000 100000 1000000 10000000
	do
	
	plot_ens_cmd="python3 Plotting/plot_ensemble_stats.py /home/enda/PhD/Shell_Model/Data/Thesis/Plots/ReplaceData/POData/Repl-"$r"/ POEnsREPL-"$r" "
	for batch in {0..1} 
	do
		for id in {1..50} ## batchsize
		do
			batch_size=50
			batch_id=`expr $batch \* $batch_size + $id`

			plot_ens_cmd+="/home/enda/PhD/Shell_Model/Data/AOReplaceData_PO/HD-INTFACRK4-AOFXDPHASE_N[25]_T["$t0","$dt","$T"]_SMP["$save_every","$stats_every","$r"]_NU[5e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[NONE,0,0.100]_u0[AO_INPUT_PHASE_REPLACE]_TAG["$batch_id"-"$r"] "
		done
	done
	echo $plot_ens_cmd
	$plot_ens_cmd &
done

wait