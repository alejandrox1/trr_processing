#!/bin/bash

file1=$1
file2=$2
nnodes=$3
ppn=$4
name=$5

echo " Script is checking whether the MD simulation ran as expected..." > README.txt
echo " Looking for files $file1 and $file2 ..." >> README.txt
if [ -a $file1 ] && [ -a $file2 ]; then
        echo
	echo " So far, we are good!  " >> README.txt

        if [ "$file1" == "min.gro" ]; then
                sbatch --nodes=${nnodes} -J ${name} terra_trpcage.sh ${ppn}

        elif [ "$file1" == "simul.gro" ]; then
                sbatch --nodes=${nnodes} -J ${name} denature.sh ${ppn}

	elif [ "$file1" == "denature.gro" ]; then
		#sbatch production.sh
		continue
	fi
else
        echo " oops, something went wrong..." >> README.txt
fi

