#!/bin/bash
#SBATCH --partition=pdlabs
#SBATCH --time=1:00

declare -r exec=$1
declare -i n=$2
declare -i k=$3

./sequential $n $k $exec.txt
