#!/bin/bash

#SBATCH --exclusive           # Individual nodes
#SBATCH -t 2-0                # Run time (hh:mm:ss)
#SBATCH -J Preferential
#SBATCH -o slurm.%j.out
#SBATCH -p normal
#SBATCH -N 2
#SBATCH -n 32
#SBATCH -A TG-MCB130178 

module load intel
module load impi
gmxdir=/home1/03561/alarcj/gromacs/exe/gromacs-4.0.7_flatbottom/exec/bin

nnodes=$1
ppn=$2
name=$3

echo -e "0\nq\n" | $gmxdir/make_ndx -f min.gro -o index.ndx
$gmxdir/grompp -f simul_prep.mdp -c min.gro -n index.ndx -p tpr.top -o simul_prep.tpr 
ibrun $gmxdir/mdrun_mpi -np ${ppn} -v -deffnm simul_prep -append -cpi simul_prep.cpt -cpo simul_prep.cpt -maxh 47.95

$gmxdir/grompp -f simul.mdp -c simul_prep.gro -n index.ndx -p tpr.top -o simul.tpr
ibrun $gmxdir/mdrun_mpi -np ${ppn} -v -deffnm simul -append -cpi simul.cpt -cpo simul.cpt -maxh 47.95


echo -e "0\n" | $gmxdir/trjconv -f simul_prep.trr -s simul_prep.tpr -pbc mol -o simul_prep.xtc
echo -e "0\n" | $gmxdir/trjconv -f simul.trr -s simul.tpr -pbc mol -o simul.xtc

rm *trr

#mkdir pre_folded
#cd pre_folded
#mv ../simul.xtc .
#mv ../simul.tpr .
#cp ../tpr.top .
#echo -e "0\n" | trjconv_s -f simul.xtc -s simul.tpr -skip 100 -sep -o step_folded.gro
#mv simul.xtc ../
#mv simul.tpr ../
#cd ../


#./run.sh simul.gro simul.tpr ${nnodes} ${ppn} ${name}
