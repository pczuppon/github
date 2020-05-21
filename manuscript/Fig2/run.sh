#!/bin/bash

g++ -std=c++11 sim_virus.cpp `pkg-config --libs gsl`

for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
    for j in 0 1
    do
	for k in 0 1
	do
	    qsub script.sh "$i" "$j" "$k"
        done
    done
done
