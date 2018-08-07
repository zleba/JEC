#!/usr/bin/env python
from JEC.treeProducer.parallel import *

inDir='/pnfs/desy.de/cms/tier2/store/user/rzlebcik/JetHT/puppi2_JetHT_Run2016C-07Aug17-v1'
outFile='/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/runB.root'
parallel('processor', inDir, outFile)
