#!/bin/bash

t_steps=(0.01 0.005 0.0025 0.00125 0.000625 0.0003125 0.00015625 0.000078125 0.0000390625 0.00001953125)
n=0

for dt in "${t_steps[@]}"
do 
	let n++
	echo "$n" "$dt"
	cmd="Solver/bin/solver_hydro_IF -o /home/enda/PhD/Shell_Model/Data/Convergence_Study/Full/ -n 19 -s 0.00000 -e 1.0 -T 0 -T 0.0 -c 0 -c 0.900000 -h "$dt" -h 0 -a 0.3333333333 -b 0.3333333333 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 0.0000001 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i N_SCALING -t Err-"$n" -f NONE -f 1 -f 0.100 -p 1 -p 1 -p 0 -z NONE -z NONE -z 0"
	# cmd="Solver/bin/solver_hydro_IF -o /home/enda/PhD/Shell_Model/Data/Convergence_Study/Full/ -n 19 -s 0.00000 -e 1.0 -T 0 -T 0.0 -c 0 -c 0.900000 -h "$dt" -h 0 -a 0.3333333333 -b 0.3333333333 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 0.000 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i N_SCALING -t Err-"$n" -f NONE -f 1 -f 0.100 -p 1 -p 1 -p 0 -z NONE -z NONE -z 0"
	echo "$cmd"
	$cmd &
done

wait

# n=0
# for dt in "${t_steps[@]}"
# do 
# 	let n++
# 	echo "$n" "$dt"
# 	cmd="Solver/bin/solver_amp_only_fxd_phase -o /home/enda/PhD/Shell_Model/Data/Convergence_Study/AmpOnly/ -n 19 -s 0.00000 -e 1.0000 -T 1 -T 0.1 -c 0 -c 0.900000 -h "$dt" -h 0 -a 0.3333333333 -b 0.3333333333 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 0.0005 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i N_SCALING -t Err-"$n" -f DELTA -f 1 -f 0.100 -p 1 -p 1 -p 0 -z NONE -z NONE -z 0"
# 	echo "$cmd"
# 	$cmd &
# done

# wait


# n=0
# for dt in "${t_steps[@]}"
# do 
# 	let n++
# 	echo "$n" "$dt"
# 	cmd="Solver/bin/solver_phase_only -o /home/enda/PhD/Shell_Model/Data/Convergence_Study/PhaseOnly/ -n 19 -s 0.00000 -e 1.0000 -T 1 -T 0.1 -c 0 -c 0.900000 -h "$dt" -h 0 -a 0.3333333333 -b 0.3333333333 -w 0.050 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 0.0005 -v 0 -v 2.0 -d 0 -d 0 -d -2.0 -i N_SCALING -t Err-"$n" -f DELTA -f 1 -f 0.100 -p 1 -p 1 -p 0 -z NONE -z NONE -z 0"
# 	echo "$cmd"
# 	$cmd &
# done

# wait