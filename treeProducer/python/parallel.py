#!/usr/bin/env python


#nProc = 19
#tmpDir = '/tmp'
import os
from subprocess import check_output
#tmpDir = '/nfs/dust/cms/user/'+os.environ["USER"]+ '/tmp/'
tmpDir = check_output(["getTempDir"]).rstrip()

def setEnv(cmd, exDir):
    myEx   = check_output(["which", cmd]).rstrip()
    myLibs = check_output("ldd  " +  myEx + " | grep $CMSSW_BASE | awk '{print $3}'" , shell=True).splitlines()

    from shutil import copy
    #print 'RADEK', myEx, exDir
    copy(myEx, exDir)
    for f in myLibs:
        copy(f, exDir)
    #print 'After', myEx, exDir
    
    import os
    my_env = os.environ.copy()
    my_env["LD_LIBRARY_PATH"] = exDir+":" + my_env["LD_LIBRARY_PATH"]
    my_env["PATH"] = exDir+":" + my_env["PATH"]

    return my_env



def parallel(cmd, In, Out, nProc=0, *argv):
    cmd = cmd.strip()

    import multiprocessing
    if nProc==0:
        nProc = multiprocessing.cpu_count()


    from subprocess import Popen, call 
    import os.path
    from random import randint
    tag = 'Parallel_'+ cmd + str(randint(1, 10000000))
    outDir = tmpDir + '/' + tag
    os.mkdir(outDir)
    myEnv = setEnv(cmd,outDir)
    OutTmp = [outDir + '/'+str(i)+'.root' for i in range(nProc)]

    argvStr = [str(s) for s in argv]

    print "Parallel run start : " + tag

    procs = []
    for i in range(nProc):
        shellCMD = [cmd, In, OutTmp[i]] + argvStr +  [str(nProc), str(i)]
        print shellCMD
        procs.append(Popen(shellCMD, env=myEnv))

    for p in procs: p.wait()

    print "Parallel run done"
    #Check that files are there
    for f in OutTmp:
        if not os.path.isfile(f):
            print 'File '+f+' does not exist'
            os.exit(1);

    print "Merging rootFiles : start"
    call(['hadd','-f', Out] + OutTmp)
    print "Merging rootFiles : end"

