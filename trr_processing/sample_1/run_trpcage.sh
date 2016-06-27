#!/bin/bash

#SBATCH --exclusive           # Individual nodes

#####SBATCH --nodes=1             # Total number of mpi tasks requested
######SBATCH -J TRP_CAGE

#-
#-                      run_trpcage.sh
#-              TC10B on solution with ARG
#-              Made: Setptember 2014
#-              Last Updated: Thu Jun 11 09:17:45 EDT 2015
#- 
##
## Usage: sbatch run_trpcage.sh <Number of ARGinines>
##
##      -h      Show help options.
##      -v      Print version info.
##
## This script set up, minimizes, and equilibrates the system required later 
## for preferential coefficent interaction calculations.
## This script will run a series of simulations that will proceed equilibrate 
## the trp-cage mini-protein in both the folded and unfolded states.
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
module load gromacs/407
L=$1   #"4.7"
Nmg=$2
nnodes=$3
ppn=$4
name=$5


#echo "\n"
#echo "Please input the number of cosolutes to be added: "
#read Nmg


# Check wether the correct input is being given
if [ $# -ne 1 ]; then
        if [ "$input" -ne "-h" ] || [ "$residue" -ne "-v" ]; then
                echo "ERROR: must run this script with an argument."
                echo "\$1 = Number of cosolvents."
                exit
        fi
fi



#-----------------------------------------------------------------------------------------------#
# Set protein in box
echo -e "5\n1\n1\n" | pdb2gmx_s -f 1L2Y.pdb -water tip3p -o trpcage.gro -p topol.top -ignh -ter
mv posre.itp posre_trp.itp
editconf_s -f trpcage.gro -o trpcage_box.gro -box $L $L $L -c -d 2.0 -bt cubic 
echo -e "5\n1\n1\n" | pdb2gmx_s -f arg.pdb -water tip3p -o arg.gro -ignh -ter -p topol_arg.top
mv posre.itp posre_arg.itp
genbox_s -cp trpcage_box.gro -ci arg.gro -nmol $Nmg -box `echo "$L - 0.25" | bc -l` `echo "$L - 0.25" | bc -l` `echo "$L - 0.25" | bc -l` -p topol.top -o trp_arg.gro

# Modify topology file
sed -i 's:Protein_A:tc5b:g' topol.top
sed -i "$ c tc5b \t 1" topol.top
sed -i 's:Protein_X:arg:g' topol_arg.top
sed -i "\$aarg \t ${Nmg}" topol.top
#sed -i "$ c arg \t 32" topol_arg.top

rm tpr.top	# It doesnt work if the output already exists
./py_itp_maker.py topol.top tpr.top tpr.itp
./py_itp_maker.py topol_arg.top arg.top arg.itp 

sed -i '/#include "ions.itp"/a #include "tpr.itp"' tpr.top
sed -i '/#include "tpr.itp"/a #include "arg.itp"' tpr.top
sed -i '1s/^/[ moleculetype ] /' tpr.itp
sed -i '1s/^/[ moleculetype ] /' arg.itp
sed -i 's:posre.itp:posre_trp.itp:g' tpr.itp
sed -i 's:posre.itp:posre_arg.itp:g' arg.itp

# Solvate simulation box
genbox_s -cp trp_arg.gro -cs spc216.gro -p tpr.top -o cage_solvated.gro
# Neutrilize the system
grompp_s -f min.mdp -c cage_solvated.gro -p tpr.top -o cage_ions.tpr
echo -e "12\n" | genion_s -s cage_ions.tpr -o cage_ions.gro -p tpr.top -pname NA+ -nname CL- -nn `echo "$Nmg +1" | bc -l`
# Minimization
grompp_s -f min.mdp -c cage_ions.gro -p tpr.top -o min.tpr
mpirun -np ${ppn} mdrun_sm -np ${ppn} -v -deffnm min -append -cpi min.cpt -cpo min.cpt


echo -e "0\n" | trjconv_s -f min.trr -s min.tpr -pbc mol -o min.xtc

