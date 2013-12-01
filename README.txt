HWW(2l2v) CMSSW Release Validation
contact: Jason Gran (jgran@physics.ucsb.edu)

https://github.com/jgran/HWWValidation



--------------------------------------------
To install and run
--------------------------------------------

cd CMSSW_5_3_11/src/ (or whatever release you are using)

git clone https://github.com/jgran/HWWValidation

cd HWWValidation/HWWBase

cmsenv

scram b -j10

cmsRun python/HWWValidation_cfg.py



This produces HWWBase/cutflow.txt and files in HWWBase/cut/ containing the run, lumi, and event# for each event passing each cut.



--------------------------------------------
To do
--------------------------------------------

-Test in other CMSSW releases. Currently only tested in CMSSW_5_3_11.
-Change default input file location, currently only runs on the uaf at UCSD.
-Add options to config file for input file, maxEvents, etc.
-Add utility to compare events from different releases. 
