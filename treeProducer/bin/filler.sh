#!/bin/bash
inFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/ntuples/runs4.root
outFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/histos/histo4AK8CHS.root
parallel filler $inFile $outFile
#filler $inFile $outFile 1000

