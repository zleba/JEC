export LD_LIBRARY_PATH_STORED=$LD_LIBRARY_PATH

export MYDIR=$PWD



#!/bin/bash
N=10

add=$(pwd)

for per in {B..H}
do
    res=$(./GetEntries.py /nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat/merged/jets${per}.root |& tail -1 |grep -o "[0-9]*$")
    ((N =  res / 1500000 + 1))

    echo $N $per 
    export NJOBS=$N
    export PERIOD=$per
    condor_submit jobs.submit
done
