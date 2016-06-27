#!/bin/bash

#SBATCH --exclusive           # Individual nodes
#SBATCH -t 2-0                # Run time (hh:mm:ss)
#SBATCH -J pressure
#SBATCH -o slurm.%j.out
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -A TG-MCB130178 

module load intel
module load impi
gmxdir=/home1/03561/alarcj/gromacs/exe/gromacs-4.0.7_flatbottom/exec/bin


# Production run
$gmxdir/grompp -f this.mdp -c nvt.gro -p topol.top -o simul_equil.tpr -n this.ndx -maxwarn 1
ibrun $gmxdir/mdrun_mpi -np 16 -v -deffnm simul_equil -px pullx.xvg -pf pullf.xvg -append -cpi simul_equil.cpt -cpo simul_equil.cpt -maxh 47.95

echo -e "0\n" | $gmxdir/trjconv -f simul_equil.trr -s simul_equil.tpr -pbc mol -o simul_equil.xtc

rm -rf \#* step*.pdb
rm *trr
