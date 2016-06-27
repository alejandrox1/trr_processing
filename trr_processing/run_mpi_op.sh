#!/bin/bash

#SBATCH --exclusive           # Individual nodes
#SBATCH --nodes=1             # Total number of mpi tasks requested
#SBATCH -t 2-0                # Run time (hh:mm:ss)

#-
#-			run_mpi_op.sh
#-		Osmotic Pressure RUNing script
#-		Made: June 16, 2014
#-		Last Updated: Wed Apr 22 14:46:14 EDT 2015 
#- Following the method used by Luo and Roux DOI: 10.1021/jz900079w. 
#- Contributions (many) from Dr. Akansha Saxena and Dr. Christopher Neale.
#-
##
## Usage: sbatch run_mpi_op.sh <structure> <water> <number of cosolvents>
##
##	-h	Show help options.
##	-v	Print version info.
##
## This script set up, minimizes, and equilibrates the system required to calculate the osmotic pressure
## of a given input structure (originally amino acids).
## The script is runned by "all_simul_run.sh" 
## NOTE: Input structure does not need the file handle.
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

module load mpi/openmpi
gmxdir=/home/alarcj/exe/gromacs-4.0.7_flatbottom/exec/bin

input=$1
output=$input
water=$2
L="3.0"
Nmg=1


# Generate central chamber with correct topology, using OPLSAA/L
if [ "$input" == "pro" ]; then
	if [ "$water" == "tip3p" ]; then
		cp /data/disk02/alarcj/data_base_stampede/topologies/proline/proline.top topol.top
		cp /data/disk02/alarcj/data_base_stampede/topologies/proline/proline.gro pro.gro
		cp /data/disk02/alarcj/data_base_stampede/topologies/proline/posre.itp posre.itp
	else
		if [ $water == "tip4p" ];then
			cp /data/disk02/alarcj/data_base_stampede/topologies/pro_tip4p/proline.top topol.top
			cp /data/disk02/alarcj/data_base_stampede/topologies/pro_tip4p/proline.gro pro.gro
			cp /data/disk02/alarcj/data_base_stampede/topologies/pro_tip4p/posre.itp posre.itp
		fi
	fi
else
	echo -e "5\n1\n1\n" | $gmxdir/pdb2gmx -f ${input}.pdb -water $water -o ${output}.gro -ignh -p topol.top -ter > parameters.txt #Parameters for GLY,TYR  
fi

# Add N numer of compounds to the central chmaber, N = Nmg. Molecules should be limited to 0.5 nm from the edges of the box to be set up later
$gmxdir/genbox -ci ${output}.gro -nmol $Nmg -box `echo "$L - 0.5" | bc -l` `echo "$L - 0.5" | bc -l` `echo "$L - 0.5" | bc -l` -p topol.top -o ${Nmg}_${output}.gro

# Modify topology in order to have the correct the number of molecules
if [ "$input" == "pro" ]; then
	sed -i "$ c Protein_chain_X \t ${Nmg}" topol.top
elif [ "$input" == "arg" ] || [ "$input" == "glu" ]; then
	sed -i "$ c Protein_X \t ${Nmg}" topol.top
else
	sed -i "$ c Protein \t ${Nmg}" topol.top
fi

# Defining central chamber's dimensions
$gmxdir/editconf -f ${Nmg}_${output}.gro -o ${output}_box.gro -box $L $L $L 
# Solvating the central chamber
if [ "$water" == "tip3p" ]; then
	$gmxdir/genbox -cp ${output}_box -cs spc216.gro -p topol.top -o ${output}_solvated.gro
else
	if [ $water == "tip4p" ];then
		$gmxdir/genbox -cp ${output}_box -cs $water -p topol.top -o ${output}_solvated.gro
	fi
fi



# Minimization of central chamber
$gmxdir/grompp -f chamber_min.mdp -c ${output}_solvated.gro -p topol.top -o chamber_min.tpr 
$gmxdir/mdrun -v -deffnm chamber_min 


echo -e "0\n" | $gmxdir/trjconv -f chamber_min.trr -s chamber_min.tpr -pbc mol -o chamber_min.xtc
