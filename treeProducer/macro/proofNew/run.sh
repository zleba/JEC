#!/bin/bash
#$ -V
#$ -cwd
##$ -e /afs/desy.de/user/z/zlebcr/cms/CMSSW_9_3_0/src/KKousour/TopAnalysis/test/farm/err/jetsB_001.err
##$ -o /afs/desy.de/user/z/zlebcr/cms/CMSSW_9_3_0/src/KKousour/TopAnalysis/test/farm/out/jetsB_001.out
#$ -l h_vmem=2G

root -l -b -q "runProof.C($1,$2,'$3')"
