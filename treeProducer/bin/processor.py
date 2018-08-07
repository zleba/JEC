#!/usr/bin/env python
from JEC.treeProducer.parallel import *

inDir='/pnfs/desy.de/cms/tier2/store/user/rzlebcik/JetHT/puppi3_JetHT_Run2016C-07Aug17-v1'
outFile='/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/runTest.root'
parallel('processor', inDir, outFile)
