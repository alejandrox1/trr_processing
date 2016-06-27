#!/bin/bash

#SBATCH --exclusive           # Individual nodes
#SBATCH -t 2-0                # Run time (hh:mm:ss)
#SBATCH -J osmotic
#SBATCH -o slurm.%j.out
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -A TG-MCB130178 

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

module load intel
module load impi
gmxdir=/home1/03561/alarcj/gromacs/exe/gromacs-4.0.7_flatbottom/exec/bin

input=$1
output=$input
water=$2
L="3.9"  #3.9
Nmg=$3


# Check wether the correct input is being given
if [ $# -ne 3 ]; then
	if [ "$input" -ne "-h" ] || [ "$residue" -ne "-v" ]; then
		echo "ERROR: must run this script with three arguments."
		echo "\$1 = input structure."
		echo "\$2 = output name. No \".gro\""
		echo "\$3 = Number of residues."
		exit
	fi
fi


# Generate central chamber with correct topology, using OPLSAA/L
if [ "$input" == "pro" ]; then
	cp /scratch/03561/alarcj/topologies/proline/proline.top topol.top
	cp /scratch/03561/alarcj/topologies/proline/proline.gro pro.gro
	cp /scratch/03561/alarcj/topologies/proline/posre.itp posre.itp
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
ibrun $gmxdir/mdrun -v -deffnm chamber_min 



# Create second chamber, fill with water and aling it with the central chamber
if [ "$water" == "tip3p" ];then
	$gmxdir/genbox -cs spc216.gro -box $L $L $L -o water_${L}nm.gro
else
	if [ "$water" == "tip4p" ];then
		$gmxdir/genbox -cs $water -box $L $L $L -o water_${L}nm.gro
	fi
fi
$gmxdir/editconf -f water_${L}nm.gro -o water_${L}nm_top.gro -box $L $L $L -center 0 0 `echo "scale=2; $L + $L/2" |bc -l`
$gmxdir/editconf -f chamber_min.gro -o central_shifted.gro -box $L $L $L -center 0 0 `echo "scale=2; $L/2" | bc -l`



# Merging the two chambers and fixing topology
./combine2gro.sh central_shifted.gro water_${L}nm_top.gro system.gro

if [ "$input" == "arg" ]; then
	$gmxdir/grompp -f min.mdp -c system.gro -p topol.top -o sys_ions.tpr >> parameters.txt
	echo -e "12\n" | $gmxdir/genion -s sys_ions.tpr -o sys_ions.gro -p topol.top -pname NA+ -nname CL- -nn $Nmg >> parameters.txt
	# Minimization
	$gmxdir/grompp -f min.mdp -c sys_ions.gro -p topol.top -o min.tpr
	ibrun $gmxdir/mdrun_mpi -np 16 -v -deffnm min 
	echo -e "0\n" | $gmxdir/trjconv -f min.trr -s min.tpr -pbc mol -o min.xtc
	
elif [ "$input" == "glu" ]; then
	$gmxdir/grompp -f min.mdp -c system.gro -p topol.top -o sys_ions.tpr >> parameters.txt
        echo -e "12\n" | $gmxdir/genion -s sys_ions.tpr -o sys_ions.gro -p topol.top -pname NA+ -nname CL- -np $Nmg >> parameters.txt
        # Minimization
        $gmxdir/grompp -f min.mdp -c sys_ions.gro -p topol.top -o min.tpr
        ibrun $gmxdir/mdrun_mpi -np 16 -v -deffnm min 
	echo -e "0\n" | $gmxdir/trjconv -f min.trr -s min.tpr -pbc mol -o min.xtc
        
else
	$gmxdir/grompp -f min.mdp -c system.gro -p topol.top -o min.tpr >> parameters.txt
	ibrun $gmxdir/mdrun_mpi -np 16 -v -deffnm min 
	echo -e "0\n" | $gmxdir/trjconv -f min.trr -s min.tpr -pbc mol -o min.xtc
        
fi

# Canonical ensemble equilibration
echo -e "0\nq\n" | $gmxdir/make_ndx -f min.gro -o index.ndx
$gmxdir/grompp -f nvt.mdp -c min.gro -n index.ndx -p topol.top -o nvt.tpr
ibrun $gmxdir/mdrun_mpi -np 16 -v -deffnm nvt 
echo -e "0\n" | $gmxdir/trjconv -f nvt.trr -s nvt.tpr -pbc mol -o nvt.xtc


# Generate mdp, index file for production run
./make_pull_mdp.sh simul.mdp nvt.gro $L $Nmg

rm \#*
rm *trr

# Production run
#mv run_simul_stampede.sh ${Nmg}nmg_op_${input}.sh
#sbatch ${Nmg}nmg_op_${input}.sh $L $Nmg 


