#!/usr/local/bin/python2.7

import os
curDir = os.getcwd()
listADDRESS = curDir + '/fileLists/Aug17'
#listADDRESS = curDir + '/fileLists/QCD_TuneCUETP8M1_13TeV_pythia8'

#outDir = '/nfs/dust/cms/user/zlebcr/JEC/ntuplesJCunc'
outDir = '/nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat'
pemFile = '/afs/desy.de/user/z/zlebcr/proxy.pem'

step=6

Tags = filter(lambda x: x.endswith('.txt'), os.listdir(listADDRESS))
Tags = [s.replace('.txt','') for s in Tags]
print Tags


for run in Tags:
    File= run + '.txt'
    nFiles = sum(1 for line in open(listADDRESS +'/' + File ))

    for i in range( (nFiles-1) / step + 1):
        #print nFiles, i+1, step*i+1
        name='jets{}_{:0>3}'.format(run, i+1)
        subFile='sub/' + name + '.sub'

        steps = ''
        fileTags = ''
        for k in range(1, step+1):
            if (i*step+k) > nFiles: break
            steps += str(i*step+k) + ' '
            fileTags += str(k) + ' '
        #steps    = ' '.join(str(i*step+k) for k in range(1,step+1))
        #fileTags = ' '.join(str(k) for k in range(1,step+1))
        #print steps
        #continue


        fSub = open(subFile, "w")
        fSub.writelines('\n'.join([       
        "#!/bin/bash",
        "#$ -V",
        "#$ -e {}/logs/{}.err".format(curDir, name),
        "#$ -o {}/logs/{}.out".format(curDir, name),
        "#$ -l h_vmem=3G",
        "#$ -l s_vmem=3G",
        "pwd",
        "export X509_USER_PROXY="+pemFile,
#        "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/desy.de/user/z/zlebcr/cms/CMSSW_9_3_0/src/KKousour/TopAnalysis/plugins/",
         "cmsRun "+curDir+"/flatData-new.py  listFile="+listADDRESS+"/"+File+
         " startFile="+str(i*step+1)+"  nFiles="+str(step)+
         " outputFile="+ outDir+"/"+name+".root",
#        "time parallel cmsRun "+curDir+"/flatData-new.py  listFile="+listADDRESS+"/"+File+
#        " startFile={1}  nFiles="+str(1)+
#        " outputFile="+ outDir+"/"+name+"_{2}.root  \"2>\" "+ curDir +"/out/"+name+"_{2}.out  ::: " + steps + " :::+ " + fileTags,
        ]))
        fSub.close()

        import stat
        os.chmod(subFile, stat.S_IRWXU + stat.S_IRGRP)


