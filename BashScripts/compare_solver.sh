#!/bin/bash

rk_command="Solver/bin/solver_hydro_RK -o ./Data/Test/ -n 10 -s 0.00000 -e 10.0 -T 0 -c 0 -c 0.900000 -h 0.00005 -h 0 -a 1.5000000000 -b 1.5000000000 -w 1.000 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 0.0000001 -v 0 -v 2.0 -d 0.0000000000 -d 0 -d -2.0 -i N_SCALING -t Compare-Solver -f NONE -f 0 -f 0.100 -p 50"
if_command=" Solver/bin/solver_hydro_IF -o ./Data/Test/ -n 10 -s 0.00000 -e 10.0 -T 0 -c 0 -c 0.900000 -h 0.00005 -h 0 -a 1.5000000000 -b 1.5000000000 -w 1.000 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 0.0000001 -v 0 -v 2.0 -d 0.0000000000 -d 0 -d -2.0 -i N_SCALING -t Compare-Solver -f NONE -f 0 -f 0.100 -p 50"
ab_command="Solver/bin/solver_hydro_AB -o ./Data/Test/ -n 10 -s 0.00000 -e 10.0 -T 0 -c 0 -c 0.900000 -h 0.00005 -h 0 -a 1.5000000000 -b 1.5000000000 -w 1.000 -w 2.000 -y 0.5000000000000000 -y 0.3333333333333333 -v 0.0000001 -v 0 -v 2.0 -d 0.0000000000 -d 0 -d -2.0 -i N_SCALING -t Compare-Solver -f NONE -f 0 -f 0.100 -p 50"


echo -e "\nCommand run:\n\t \033[1;36m $rk_command \033[0m"
$rk_command

echo -e "\nCommand run:\n\t \033[1;36m $if_command \033[0m"
$if_command

echo -e "\nCommand run:\n\t \033[1;36m $ab_command \033[0m"
$ab_command

test_command="python3 test_solver.py -t Compare-Solver --compare"

echo -e "\nCommand run:\n\t \033[1;36m $test_command \033[0m"
$test_command