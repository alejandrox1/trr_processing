#!/bin/bash

#	Updated: Mon Jun 22 10:40:58 EDT 2015
#	Last updated: Fri Jun 19 15:49:05 EDT 2015
# Specify number of nodes, name of run, box size for simulation (preferably >= 6.5).
#
# In terra the best set up is:
# ./system_setup.sh 2 32 <name> 6.5
#
 
curdir=$(pwd)

struct=arg
nnodes=1
ppn=2
name=GROMACS
size=6.5


cd ../
dir=simul_pro
mkdir ${dir}
cd ${dir}

for j in 1; do

	for i in 4; do
		
		cp -r ${curdir} pi_${i}_${j};
		cd pi_${i}_${j};	

		sbatch --nodes=${nnodes} run_trpcage.sh ${size} ${i} ${nnodes} ${ppn} ${name} ${struct} 

		cd ../
	done
done
