#!/bin/bash


data_every=1000
stats_every=50000

repl=1

for repl in 1 2 5 10 25 50 100 1000 10000
do
	python_ens_cmd="python3 Plotting/plot_ensemble_stats.py Repl-"$repl" "
	for r in {0..4}
	do
		############################## Run Solver
		for i in {0..50}
		do
			num1=50
			num2=$r
			num3=$i
			tag=`expr $num1 \* $num2 + $num3`

			# solv_cmd="Solver/bin/solver_amp_only_fxd_phase -o ./Data/AOReplaceData/ -n 25 -s 200.00000 -e 1000.00000 -T 0 -T 0.2 -c 0 -c 0.900000 -h 0.0001000000000000 -h 0 -a 0.3333000000 -b 0.3333000000 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 5e-07 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i AO_INPUT_PHASE_REPLACE -t EnsREPL-"$tag"-"$repl" -f DELTA -f 1 -f 0.100 -p "$data_every" -p "$stats_every" -p "$repl" -z ./Data/DenseData/HD-INTFACRK4-FULL_N[25]_T[0.0,0.0001,1000.000]_SMP[1,50000,1000]_NU[5e-07]_ALPHA[1.500]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[N_SCALING_RAND]_TAG["$tag"]/Main_HDF_Data.h5 -z NONE -z 0"
			# echo $solv_cmd
			# $solv_cmd &

			python_ens_cmd+="./Data/AOReplaceData/HD-INTFACRK4-AOFXDPHASE_N[25]_T[200.0,0.0001,1000.000]_SMP[10000,10000,"$repl"]_NU[5e-07]_ALPHA[0.333]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[AO_INPUT_PHASE_REPLACE]_TAG[EnsREPL-"$tag"-"$repl"] "
		done

		wait

	done
	echo $python_ens_cmd
	$python_ens_cmd

	wait
done
