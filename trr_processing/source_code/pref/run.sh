#!/bin/bash

#-
#-                      run.sh
#-      Preferential Interaction Coefficient
#-      Made: Wed Jun 1 2015i
#-	Updated: Wed Jun 17 15:20:09 EDT 2015
#-      Previous Update: Wed Jun  3 10:36:32 EDT 2015
#- 	Contributions (many) from Ryan Krafnick.
#-
##
## Usage: ./run.sh <threads> <topology> <protein structure> <cosolvent structure>
##	i.e. ./run.sh 3 tpr.top trpcage.pdb arg.pdb 
##
##      -h      Show help options.
##      -v      Print version info.
##
## This script will compile all necessary files for calculating the preferential interaction
## coefficient.
##
## DEPEDENCIES: coef.cpp  concentration.cpp  file_pro.cpp  INFO.cpp  pref_co.cpp  spatial.cpp
##
## OUTPUT: pic_calculation.dat
##
## NOTES: For input file structures use gro files instead of pdbs. (Beware of -ter)
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


cp ../trpcage.gro .
cp ../arg.gro .

topology=tpr.top						# tpr.top
protein=trpcage.gro
cosolvent=arg.gro

N_co=`tail -4 ${topology} | awk 'FNR == 1 {print $2}'`
protein_size=`awk 'NR==2' ${protein}`
co_size=`awk 'NR==2' ${cosolvent}`

mv trpcage.gro ../
mv arg.gro ../

echo ${N_co}
tail -4 ${topology}
time ./calc -n ${N_co} -p ${protein_size} -c ${co_size} 
