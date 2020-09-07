#!/bin/bash
#SBATCH --job-name=temp
#SBATCH --nodes=1
#SBATCH --partition=pdlabs
#SBATCH --gres=gpu:1
#SBATCH --time=1:00

declare -r exec=$1
declare -i n=$2
declare -i k=$3

if [ "$exec" == "v1" ]
then
	./v1 $n $k $exec.txt
elif [ "$exec" == "v2" ]
then
    ./v2 $n $k $exec.txt
elif [ "$exec" == "v3" ]
then
	./v3 $n $k $exec.txt
else
	echo First argument must be "seq" or "v1" or "v2" or "v3"
fi
