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

ppn=$1

$gmxdir/grompp -f production_prep.mdp -c denature.gro -n index.ndx -p tpr.top -o production_prep.tpr 
ibrun $gmxdir/mdrun_mpi -np ${ppn} -v -deffnm production_prep -append -cpi production_prep.cpt -cpo production_prep.cpt -maxh 47.95

$gmxdir/grompp -f production.mdp -c production_prep.gro -p tpr.top -n index.ndx -o production_run.tpr
ibrun $gmxdi/mdrun_mpi -np ${ppn} -v -deffnm production_run -append -cpi production_run.cpt -cpo production_run.cpt -maxh 47.95

