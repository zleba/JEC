#!/bin/bash
inDir=/pnfs/desy.de/cms/tier2/store/user/rzlebcik/JetHT/puppi3_JetHT_Run2016C-07Aug17-v1
outFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/runTest.root
#processor $inDir $outFile 100
#parallel processor $inDir $outFile -j5
#exit
#submit processor $inDir $outFile

for per in {B..H}
do
    inDir=`ls -d /pnfs/desy.de/cms/tier2/store/user/rzlebcik/JetHT/puppi3_JetHT_Run2016$per*`
    outFile=/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/run2${per}.root
    echo $inDir
    submit processor $inDir $outFile &
done
