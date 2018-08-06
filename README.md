# JEC

## Installation
```
scram p CMSSW_9_3_0
cd CMSSW_9_3_0/src
git clone git@github.com:zleba/JEC.git
cd JEC/treeProducer
cmsenv
scram b
```

## Production of the nTuples
The `treeProducer.cc` is in `plugins` and contains code to produce user trees (ntuples).
To see the jet properties currently stored in the ntuple, see the file `interface/QCDjet.h`

The scripts related to the running of the code producing nTuples are stored in `python` directory.
Especially consider these files:
`treeProducerCrab.py` -- the configuration for the `cmsRun` command
`dataSets.py` -- file containing list of data sets which can be currently used
`crabConfig.py` -- script to run nTuple production on grid
`test.py` -- script to run short test (runs cmsRun with proper setup)

Before running, please init your grid-certificate and crab:
```
voms-proxy-init --voms cms  -valid 192:00
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

As a next step, I recommend to run test locally with few events.
```
./test.py
```

If the test run was successful, one can proceed to run on grid over the whole data sample
```
./crabConfig.py
```
Here notice the first lines of the crabConfig script
```
fileTag = 'puppi2_'
from dataSets  import *
sets =  data16
```
the `fileTag` is the nameTag to distinguish between different runs of the ntupliser.
The list of data sets are imported from `dataSets.py` file and consequently the script submit all 2016 data sets to grid.
To run over 2016+2017, one needs to specify `sets = data16 + data17`.


When the run is finished, the root files with ntuples should be stored in the following directory
```
/pnfs/desy.de/cms/tier2/store/user/$USER/JetHT
```

Since cram produce grid pack from the whole content of the JEC directory, please don't put there any large files (e.g. heavy root files), otherwise the grid pack will be heavy and the transfer to the node would be slow.


## Analysing the nTuples - General remarks
The nTuples from the previous step, where the production time is several hours are supposed to contain stable content. 
Whereas analysis specific stuff, changing from time to time is calculated fast at the nTuple level.
That includes for example JECs.

In general the nTuplisers are processed using short compiled programs which are supposed to have the following arguments.
```
myAnalyser input output extraPars nSlices iSlice
```
`myAnalyser` - some name of the program to analyse ntuples
`input`  - is the input root file with nTuple or input directory containing several root files
`output` - is single output root file (containing histograms and/or nTuple)
`extraPars` - these are possible users specific extra parameters (not mandatory)
`nSlices`  - number of divisions of the whole input sample (these divisions have equal statistics), default value is 1 = no division
`iSlice`   - current slice index, can have value between 0..nSlices-1, default value is 0 = run on the first slice

These kind of programs can be used either to fill histograms from nTuple or to modify/filter events in the given nTuple.
This is inspired by [fold](https://en.wikipedia.org/wiki/Fold_(higher-order_function))/[map](https://en.wikipedia.org/wiki/Map_(higher-order_function))/[filter](https://en.wikipedia.org/wiki/Filter_(higher-order_function)) paradigm of the functional programing.

For short test run, one can use something like:
```
myAnalyser input output extraPars 1000
```
which will run only on the first slice from 1000 (=only on few events).
The motivation for such argument structure is to allow simple parallelization, as each slice can run independently and, in the next step, the outputs have to be merged using `hadd`.

The parallelization is done in the automatic way, by putting an extra command in front of given program, it can look like:
```
parallel myAnalyser input output extraPars
submit   myAnalyser input output extraPars
```
These commands always analyse the whole input and save the result to output.
The `parallel` command run the program parallel on cores of the current working machine, whereas `submit` command submit the task to HTcondor batch system.
The intermediate results from runs over the slices are stored in the temp directory and can be investigated for debugging purposes. 
The overall speed depends on the given working machine and complexity of the task.
In general for shorter tasks `paralell` is faster as there is no waiting time.
When `submit` command is used, it takes some time for jobs to start running, but there can be more of them than in `parallel`.


##Adding JEC and calculating trigger weight (nTuple -> nTuple)
The evaluation of the JEC correction by `JECcorrector` is quite slow, therefore, before the histograms are filled, the intermediate step is introduced to correct the jet pTs for the JEC.





The python config file for `cmsRun` is called `farm/treeProducer.py`.
1) Fill the current directory into variable `currDir` of `treeProducer.py`
2) Try test run by `cmsRun treeProducer.py`
3) Run `./createJobs.py` to produce submit files for 2016 data, before running please fill the output directory, where the nTuples will be stored, and path to the valid `proxy.pem` file.
4) Submit all files in `sub` directory by `for i in *.sub; do qsub $i; done`
5) After all runs are finished, merge each period into single `root`-file using `hadd`

## Filling the histograms with matching properties
Currently all files fulfilling single- or double-jet trigger condition are sorted.
For each such event both AK4CHS and AK4PUPPI jets are kept.
The simple macro to run over the events and fill histogram with matching property is located in `macro/proofNew/matching.C`.
The JEC are calculated on-the-fly.
The path of the input data (from previous step) is filled in the `runProof.C`.

To run test, try `./run.sh 30 0 B`, where `B` is the data period, `30` number of divisions and `0` is the acuall piece, here it can be number from `0` to `29`.

If test is successful, run`./runAll.sh` to submit jobs to farm for whole data sample (should be finished in ~10 minutes).
The resulting histogram files in `histos` folder needs to be merged into single `root`-file by `hadd`.

## Plotting histograms from the root-files with histograms
Currently there is simple plotting macro in the `macro/proofNew` directory called `plotter.C`, to simplify the code, some functions from `plottingHelper` library are employed, more information on:
```https://github.com/zleba/PlottingHelper```
for plotting itself do `root -l -b -q plotter.C`, check the input and output file which are hard-coded in the beginning of `void plotter()` function.


## TODO
Move to `ht_condor` batch system.
