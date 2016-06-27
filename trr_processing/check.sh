#!/bin/bash

touch simulation_results.txt
for i in 1 5 10 20 25 30 40 45 50 55 60 65; do
        for j in 1 2 3 4 5; do
                find ${i}nmg_${j}run/pullf.xvg >> simulation_results.txt 2>&1
        done
done

