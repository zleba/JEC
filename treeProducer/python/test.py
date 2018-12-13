#!/usr/bin/env python
from dataSets import *

def runTest(testSet):
    inputSets = GetFiles(testSet, 1)
    print inputSets, len(inputSets)
    inputSet = inputSets[0]
    #return
    output="output_"+Short(testSet).replace('/','_')+".root"
    from subprocess import call
    call(['cmsRun', 'treeProducerCrab.py', 'inputFiles='+inputSet, 'maxEvents=3000', 'outputFile='+output])
    #root://cms-xrd-global.cern.ch/

runTest(data16[0]);
