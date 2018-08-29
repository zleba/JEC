#!/bin/bash
inDir=/pnfs/desy.de/cms/tier2/store/user/rzlebcik/JetHT/puppi3_JetHT_Run2016C-07Aug17-v1
outFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/runTest.root
#processor $inDir $outFile 100
#exit
#parallel processor $inDir $outFile -j5
#submit processor $inDir $outFile

for per in {B..H}
do
    inDir=`ls -d /pnfs/desy.de/cms/tier2/store/user/rzlebcik/JetHT/puppi3_JetHT_Run2016$per*`
    outFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/ak4study/run3${per}.root
    echo $inDir
    submit processor $inDir $outFile &
done
