#!/bin/bash
inFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/ak4study/runs18All.root
outFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/histos/histo18Puppi.root
#submit filler $inFile $outFile
parallel filler $inFile $outFile
#filler $inFile $outFile 1000

