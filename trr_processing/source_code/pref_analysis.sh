#!/bin/bash

# Compile Analysis Code
g++ functions_pref.cpp info.cpp pref_analysis.cpp -o pre -I/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/ -L/home/alarcj/exe/gromacs-4.0.7_flatbottom/exec/lib -lgmx -lm

./pre -s simul.tpr -f simul.xtc 
./pre -s production_run.tpr -f production_run.xtc 


for conc in 1 4 8 12 16 20; do
	sbatch -J ${i}_${struct} fast_autocorrelation.sh ${struct} ${i} folded
        sbatch -J ${i}_${struct} fast_autocorrelation.sh ${struct} ${i} unfolded
done

#-----------------------------------------------------------------------------------------#
conc=$2		# Concentration
flag=$3		# folded || unfolded

touch PROG_${conc}_${flag}.out

for i in {1..5}; do
	cd pi_${conc}_${i}

	topology=tpr.top
        N_co=`tail -4 ${topology} | awk 'FNR == 1 {print $2}'`
	N_wat=`tail -3 ${topology} | awk 'FNR == 1 {print $2}'`    
	
	if [ "$flag" == "folded" ]; then
		if  [ ! -d "f_analysis" ]; then
                	mkdir f_analysis
		fi

		cd f_analysis
		for bin in 1 2 4 8 16 32 64 128 256 512 1000; do
                	time ./pre -s ../simul.tpr -f ../simul.xtc -skip 10 -block $bin -Nco $N_co -Sol $N_wat  > PROG_${conc}_${flag}.out 2>&1      
		done
		cd ../

	elif [ "$flag" == "unfolded" ]; then
		if [ ! -d "u_analysis" ]; then
			mkdir u_analysis
                fi

		cd u_analysis
		for bin in 1 2 4 8 16 32 64 128 256 512 1000; do
                	./pre -s ../production_run.tpr -f ../production_run.xtc -skip 10 -block $bin -Nco $N_co -Sol $N_wat  > PROG_${conc}_${flag}.out 2>&1
		done
		cd ../
	fi
	cd ../
done

