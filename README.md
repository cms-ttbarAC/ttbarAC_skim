# ttbarAC_skim
Skim code setup for producing Ntuples for ttbar charge asymmetry analysis


Instructions for producing EDM Ntuples with BEST tagger inputs and then flat trees for analysis:

Check out CMSSW release:
```
cmsrel CMSSW_9_0_1
cd CMSSW_9_0_1/src
cmsenv
```

Check out BESTProducer code:
```
git clone git@github.com:justinrpilot/BESTAnalysis -b 90x_prod
git clone git@github.com:cms-ttbarAC/ttbarAC_skim
scram b
```

To run EDM production:

```
cd ttbarAC_skim/ttbarAC_skim/test
voms-proxy-init
cmsRun runSkim.py
```

To submit crab jobs:

Edit ```crab_*.py``` with a new dataset name.
Make sure to change your username in the stageout area (/store/group/lpctop/YOUR_USERNAME/)

2016 Dataset names:
```
'/JetHT/Run2016B-03Feb2017_ver2-v2/MINIAOD' 
'/JetHT/Run2016C-03Feb2017-v1/MINIAOD' 
'/JetHT/Run2016D-03Feb2017-v1/MINIAOD' 
'/JetHT/Run2016E-03Feb2017-v1/MINIAOD' 
'/JetHT/Run2016F-03Feb2017-v1/MINIAOD' 
'/JetHT/Run2016G-03Feb2017-v1/MINIAOD' 
'/JetHT/Run2016H-03Feb2017_ver2-v1/MINIAOD' 
```

```crab submit -c crab.py --dryrun``` (remove dryrun for actual submission)

After crab jobs are finished, make a text file with the filenames (starting with /store/group/...):

```
eosls /store/group/lpctop/.../ | grep "ana" > files_out.txt
```


