#!/bin/zsh

for i in `seq 0 6`
do
    echo $i
    sed "s/MYID/$i/" crabConfig.py > CrConf${i}.py
    #crab submit -c CrConf${i}.py
    python CrConf${i}.py

done
