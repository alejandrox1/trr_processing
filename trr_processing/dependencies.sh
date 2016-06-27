#!/bin/bash

# <residue> should be the name of the cosolvent (excluding the file handle).
# <water> should correspond to the desired water type (e.g. tip3p, tip4p).

# Input
residue=$1
water=$2

echo
echo "Did you sibmit the right structure and wate type?"
echo "...Are you sure?"
echo
echo "Here is what I have..."
echo "Residue: " ${residue}
echo "Water: " ${water}
echo

for i in 5 15 25 35 44; do
        for j in 1 2 3 4 5; do
		cp -r /scratch/03561/alarcj/data_base_stampede/osmotic_pressure_april15 ${i}nmg_${j}run
		cd ${i}nmg_${j}run

		sbatch run_stampede_op.sh ${residue} ${water} ${i} > last_submitted.txt
		id=`awk 'END {print $NF}' last_submitted.txt`
		echo
		echo "This should match up:"
		echo ${id}
		echo	
	
		mv run_simul_stampede.sh ${i}nmg_op_${residue}.sh
		sbatch --dependency=afterok:${id} ${i}nmg_op_${residue}.sh 
		
		cd ../
	done
done
