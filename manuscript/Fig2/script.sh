#!/bin/bash

#$ -t 2
#$ -N origin
#$ -S /bin/bash
#$ -cwd

./a.out $2 $3 $1
