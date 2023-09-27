#!/bin/bash


plot_dense_data_cmd="python3 Plotting/plot_dense_ensemble_data.py "
for r in {0..0}
do
	for i in {0..5}
	do
		num1=50
		num2=$r
		num3=$i
		tag=`expr $num1 \* $num2 + $num3`

		cmd="Solver/bin/solver_hydro_IF -o ./Data/DenseData/ -n 22 -s 0.00000 -e 5000.00000 -T 1 -T 0.2 -c 0 -c 0.900000 -h 0.0001000000000000 -h 0 -a 0.33333333 -b 0.33333333 -w 0.062 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 5e-07 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i N_SCALING_RAND -t $tag -f DELTA -f 1 -f 0.005 -p 1 -p 50000 -p 1000 -z NONE -z NONE -z 0"
		# echo $cmd
		$cmd &
	done

	wait 

	# for i in {0..50}
	# do
	# 	num1=50
	# 	num2=$r
	# 	num3=$i
	# 	tag=`expr $num1 \* $num2 + $num3`

	# 	plot_cmd="python3 Plotting/plot_dense_data.py -i ./Data/DenseData/HD-INTFACRK4-FULL_N[22]_T[0.0,0.0001,1000.000]_SMP[1,50000,1000]_NU[5e-07]_ALPHA[0.333]_K[0.062,2.000]_EPS[0.50]_FORC[DELTA,1,0.005]_u0[N_SCALING_RAND]_TAG[$tag]/ --plot"
	# 	$plot_cmd &
	# done

	# wait 

	for i in {0..5}
	do
		num1=50
		num2=$r
		num3=$i
		tag=`expr $num1 \* $num2 + $num3`

		plot_dense_data_cmd+="./Data/DenseData/HD-INTFACRK4-FULL_N[22]_T[0.0,0.0001,5000.000]_SMP[1,50000,1000]_NU[5e-07]_ALPHA[0.333]_K[0.062,2.000]_EPS[0.50]_FORC[DELTA,1,0.005]_u0[N_SCALING_RAND]_TAG[$tag] "
	done	
done

echo $plot_dense_data_cmd
$plot_dense_data_cmd

