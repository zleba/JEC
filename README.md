# JEC

## Installation
```
scram p CMSSW_9_3_0
cd CMSSW_9_3_0/src
git clone git@github.com:zleba/JEC.git
cd JEC/treeProducer
scram b
```
Note that file `treeProducer/plugins/JEC.h` contains hard-coded location to the directory with JEC tables, you may replace it with the current directory. Only V7 corrections are included in this repository.

## Production of the nTuples
The `treeProducer.cc` is in `plugins`. To run this code and produce ntuples, we use scripts in `farm` folder.
The python config file for `cmsRun` is called `farm/treeProducer.py`.
1) Fill the current directory into variable `currDir` of `treeProducer.py`
2) Try test run by `cmsRun treeProducer.py`
3) Run `./createJobs.py` to produce submit files for 2016 data, before running please fill the output directory, where the nTuples will be stored, and path to the `proxy.pem` file.
4) Submit all files in `sub` directory
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
