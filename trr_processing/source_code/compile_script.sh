#!/bin/bash

filef=$1
flag=$2
gmxdir=/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/

g++ functions_pref.cpp ${filef} -o pre -I/home/alarcj/exe/gromacs-4.0.7_flatbottom/include/ -L/home/alarcj/exe/gromacs-4.0.7_flatbottom/exec/lib -lgmx -lm 

sleep 1

mv preferential_val1.dat val1.dat

if [ "$flag" == "less" ]; then
	./a.out -s ../sample_1/min.tpr -f ../sample_1/min.xtc | less
else
	./pre -s ../sample_1/min.tpr -f ../sample_1/min.xtc -p topol.top 
fi

#~ Wed Aug 19 15:34:29 EDT 2015
#
# Including the #ifdef __cplusplus directives? in cpp_readf.cpp solved the previous problem g++ had at the 
# time of linking the libraries.
