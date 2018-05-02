#!/bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=2G

#arguments
#nDivisions nNow period

root -l -b -q "runProof.C($1,$2,'$3')"
