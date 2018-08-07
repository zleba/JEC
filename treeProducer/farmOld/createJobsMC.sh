ADDRESS=`pwd`
listADDRESS=$ADDRESS/fileLists/QCD_TuneCUETP8M1_13TeV_pythia8
DUST=/nfs/dust/cms/user/zlebcr/JEC/ntuplesMC

step=3

#rm sub/*.sub
mcTags=`for i in  $listADDRESS/*.txt; do basename $i | sed 's/\.txt//'; done`

for run in $mcTags
do
    file=${run}.txt
    nFiles=`cat $listADDRESS/$file | wc -l`

    #echo $file $nFiles
    #continue
    for i in `seq  1 $step $nFiles`
    do
        name=jets${run}_$i
        subFile=sub/${name}.sub
        echo "#!/bin/bash" > $subFile
        echo "#$ -V" >> $subFile
        echo "#$ -e $ADDRESS/err/${name}.err" >> $subFile
        echo "#$ -o $ADDRESS/out/${name}.out" >> $subFile
        echo "#$ -l h_vmem=2G" >> $subFile

        #echo 'cd $TMPDIR' >> $subFile
        echo 'pwd' >> $subFile
        #echo "cp $ADDRESS/flatData-TTJets-cfg.py ." >> $subFile

        #echo "cmsRun $ADDRESS/flatData-TTJets-cfg.py  listFile=$listADDRESS/$file  startFile=$i  nFiles=$step outputFile=$DUST/${name}.root" >> $subFile 
        echo "cmsRun $ADDRESS/flatData-new.py  listFile=$listADDRESS/$file  startFile=$i  nFiles=$step outputFile=$DUST/${name}.root" >> $subFile 
        #echo "cp *.root $DUST/${name}.root" >> $subFile 

        chmod u+x $subFile

    done

done
