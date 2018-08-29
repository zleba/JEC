#!/bin/bash
inFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/ak4study/runsAll.root
outFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/histos/histoPuppi.root
parallel filler $inFile $outFile
#filler $inFile $outFile 1000

