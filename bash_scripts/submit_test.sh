#!/bin/bash

module load gcc/7.3.0 cuda/10.0.130

declare -r exec=$1

if [ "$exec" = "seq" ]
then
	make sequential
elif [ "$exec" == "v1" ]
then
	make v1
elif [ "$exec" == "v2" ]
then
    make v2
elif [ "$exec" == "v3" ]
then
	make v3
else
	echo First argument must be "seq" or "v1" or "v2" or "v3"
fi

declare -i min_n=1000
declare -i max_n=10000
declare -i inc_n=$min_n

declare -i min_k=50
declare -i max_k=100
declare -i inc_k=$min_k

declare -i n
declare -i k

for ((nc=$min_n; nc<=$max_n; nc += $inc_n))
do
	n=nc
	for ((kc=$min_k; kc<=$max_k; kc += $inc_k))
	do
		k=kc
		for ((ic=1; ic<=1; ic += 1))
        do
			if [ "$exec" = "seq" ]
			then
				./seq.sh $exec $n $k
			else
				./cuda.sh $exec $n $k
			fi
		done
	done
done
