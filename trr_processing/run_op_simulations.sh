#!/bin/bash

#-
#- 	    Osmotic Pressure Simulations
#- 		run_op_simulations 
#- Last updated: Thu Jun 26 19:09:59 EDT 2014
#-
##
## Usage: ./run_op_simulations [options] <residue> <water> ...
##
##	-h	Show help options.
##	-v	Print version info
##
## DEPENDENCIES: In oder to be able to run this script you need to have "run_mpi_op.sh".
##
## This script runs multiple simulations of different concentrations.
## Each simulation first sets up, minimizes, and equilibrates the system.
## After the system is equilibrated the production run script will be created and sent to the cluster.
##
## In order to run this script you need to specify the name of the pdb file to be used. Omit the file handle.
## The second argument needed to run is the water model to be used for simulation.
## During the original runs either TIP3P or TIP4P would be used.
## 

help=$(grep "^##" "${BASH_SOURCE[0]}" | cut -c 4-)
version=$(grep "^#-" "${BASH_SOURCE[0]}" | cut -c 4-)

opt_h() 
{
	echo "$help"
}
opt_v() 
{
	echo "$version"
}

while getopts "hv" opt; do
	eval "opt_$opt"
	exit
done

# Input
residue=$1
water=$2

# Checks for the right amount of input parameters
if [ $# -ne 2 ]; then
	if [ "$residue" -ne "-h"] || [ "$residue" -ne "-v" ]; then
		echo "ERROR: not enough arguments given."
		echo "Try ./all_run_simul.sh -h for help."
		echo "Or ./all_run_simul.sh -v for information about the script."
		exit
	fi
fi

# Creates directories, retrives the simulation script and sends simulation to cluster, then repeats

for i in 5 15 25 35 44; do
	for j in 1 2 3 4 5; do
		cp -r /data/disk02/alarcj/data_base_stampede/osmotic_pressure_run_july12 ${i}nmg_${j}run
		cd ${i}nmg_${j}run
		sbatch run_mpi_op.sh $residue $water $i
		cd ..
	done
done

