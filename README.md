# JEC

## Installation
```
scram p CMSSW_9_3_0
cd CMSSW_9_3_0/src
git clone git@github.com:zleba/JEC.git
cd JEC/treeProducer
scram b
```
Note that file `treeProducer/plugins/JEC.h` contains hardcoded location to the directory with JEC tables, you may replace it with the current directory.

## Production of the nTuples
The `treeProducer.cc` is in `plugins`. To run this code and produce ntuples, we use scripts in `farm` folder.
The python config file for `cmsRun` is called `farm/treeProducer.py`.
1) Fill the current directory into variable `currDir` of `treeProducer.py`
2) Try test run by `cmsRun treeProducer.py`
3) Run `./createJobs.py` to produce submit files for 2016 data, before running please fill the output dirrectory, where the nTuples will be stored, and path to the `proxy.pem` file.
4) Submit all files in `sub` directory
5) After all runs are finished, merge each period into single `root`-file using `hadd`

## Analyse the events
Currently all files fulfiling single- or double-jet trigger condition are sotred.
For each such event both AK4CHS and AK4PUPPI jets are kept.
The simple macro to run over the events and fill histogram with matching property is located in `macro/proofNew/matching.C`.
