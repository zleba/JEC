#!/bin/bash
N=10

add=$(pwd)

for per in {B..H}
do
    res=$(./GetEntries.py /nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat/merged/jets${per}.root |& tail -1 |grep -o "[0-9]*$")
    ((N =  res / 3000000 + 1))
    for ((i=0; i < N; ++i))
    do
        #./run.sh $N $i $per
        #break
        qsub -e $add/logs/${per}_${N}_${i}.err  -o $add/logs/${per}_${N}_${i}.out -N  ${per}_${N}_${i} run.sh $N $i $per
    done
    #break
done
