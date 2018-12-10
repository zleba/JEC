# JEC

## Installation
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
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
myAnalyser input output nSlices iSlice
```
`myAnalyser` - some name of the program to analyse ntuples
`input`  - is the input root file with nTuple or input directory containing several root files
`output` - is single output root file (containing histograms and/or nTuple)
`nSlices`  - number of divisions of the whole input sample (these divisions have equal statistics), default value is 1 = no division
`iSlice`   - current slice index, can have value between 0..nSlices-1, default value is 0 = run on the first slice

These kind of programs can be used either to fill histograms from nTuple or to modify/filter events in the given nTuple.
This is inspired by [fold](https://en.wikipedia.org/wiki/Fold_(higher-order_function))/[map](https://en.wikipedia.org/wiki/Map_(higher-order_function))/[filter](https://en.wikipedia.org/wiki/Filter_(higher-order_function)) paradigm of the functional programing.

For short test run, one can use something like:
```
myAnalyser input output 1000
```
which will run only on the first slice from 1000 (=only on few events).
The motivation for such argument structure is to allow simple parallelization, as each slice can run independently and, in the next step, the outputs have to be merged using `hadd`.

The parallelization is done in the automatic way, by putting an extra command in front of given program, it can look like:
```
parallel myAnalyser input output 
submit   myAnalyser input output 
```
These commands always analyse the whole input and save the result to output.
The `parallel` command run the program parallel on cores of the current working machine, whereas `submit` command submit the task to HTcondor batch system.
The intermediate results from runs over the slices are stored in the temp directory and can be investigated for debugging purposes. 
To change address of the temp directory, see the script
```
scripts/getTempDir
```
The overall speed depends on the given working machine and complexity of the task.
In general for shorter tasks `paralell` is faster as there is no waiting time.
When `submit` command is used, it takes some time for jobs to start running, but there can be more of them than in `parallel`.


## Adding JEC and calculating trigger weight (nTuple &rightarrow; nTuple)
The evaluation of the JEC correction by `FactorizedJetCorrector` (interface coded in file `bin/JEC.h`) is quite slow, therefore, before the histograms are filled, the intermediate step is introduced to correct the jet pTs for the JEC effects.

By default all levels of JEC are applied, one can change it in the code.

In addition to the application of the JEC, at this step the event weight is calculated in accordance with the trigger efficiency intervals based on leading (or tagged) jet pt.
The events were the corresponding trigger element was not fired are removed.
Run the pre-analysis using:
```
submit processor inputDir output.root
```
Here inputDir is the directory containing ntuples from crab run (many root files), and the output.root is single root file containing result.
The example of the run can be found in short shell script.
```
processor.sh
```
Processing of 2016 data nTuples typically takes ~10 minutes on DESY htcondor based batch system(submit command).

## Filling histograms from pre-processed ntuple (nTuple &rightarrow; histos)
In the next step the histograms can be filled from nTupes obtained in previous step.
To do so run
```
parallel filler input.root histos.root
```
Where input.root is the address of the input nTuple file and histos.root is address of the output file containing filled histograms.
Example, how to call such analysis code is given in.
```
filler.sh
```

## Plotting histograms from the root-files with histograms (histos &rightarrow; pdf)
Currently there is simple plotting macro in the `plotting` directory called `plotter.C`, to simplify the code, some functions from `plottingHelper` library are employed, more information on:
```https://github.com/zleba/PlottingHelper```
for plotting itself do `root -l -b -q plotter.C`, check the input and output file which are hard-coded in the beginning of `void plotter()` function.

