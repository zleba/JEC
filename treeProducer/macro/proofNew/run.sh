#!/bin/zsh
#$ -V
#$ -cwd
#$ -l h_vmem=2G

#arguments
#nDivisions nNow period
echo $PWD
echo $SHELL
#module use -a /afs/desy.de/group/cms/modulefiles/ && module load cmssw
#cmsenv

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_STORED:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
root-config --version

root -l -b -q "runProof.C($1,$2,'$3')"
