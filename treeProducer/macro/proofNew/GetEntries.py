#!/usr/bin/env python
from ROOT import TFile, TTree

import sys

f = TFile.Open(sys.argv[1])
tr = f.Get('ak4/events')
print tr.GetEntries()
