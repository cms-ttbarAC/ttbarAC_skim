# ttbarAC_skim

Skim code setup for producing Ntuples for ttbar charge asymmetry analysis.  
The framework is setup for running at the LPC, and saving ntuples to the LPC Tier3.


## Getting Started

Check out CMSSW release:
```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git cms-init
```

Next, checkout [BESTProducer](https://github.com/cms-ttbarAC/BESTAnalysis.git) and other necessary packages.  
There are two distinct sets of instructions depending on which dataset you will be processing.  These only differ for the EGamma dependencies:

```
### 2016 ###
# VID (Electron ID)
git cms-addpkg RecoEgamma/ElectronIdentification
git cms-addpkg PhysicsTools/SelectorUtils


### 2017 ###
# Updated Egamma IDs (creates RecoEgamma and PhysicsTools directories)
# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP


### COMMON ###

# BEST (top vs qcd tagging)
git clone https://github.com/cms-ttbarAC/BESTAnalysis.git -b 94X
# LWTNN (running BEST in c++)
mkdir lwtnn
git clone https://github.com/demarley/lwtnn.git -b CMSSW_8_0_X-compatible lwtnn/lwtnn

# Analysis code -- now linked with this repository
mkdir Analysis
git clone https://github.com/cms-ttbarAC/CyMiniAna.git Analysis/CyMiniAna
git clone https://github.com/cms-ttbarAC/ttbarAC_skim

scram b -j8
```

## Ntuple Production

Instructions for producing flat Ntuples with BEST tagger inputs.

```
cd ttbarAC_skim/ttbarAC_skim/test
voms-proxy-init -voms cms
source /cvmfs/cms.cern.ch/crab3/crab.csh
```

To submit crab jobs, the script `runSkim.py` is processed using the `cmsRun` executable.

There are two options for submitting crab jobs, (1) submit many jobs at once or (2) submit jobs manually one at a time.

### Multiple Submission

Write a text file that contains the datasets you would like to process, see `test/crab_datasets-data.txt` and `test/crab_datasets-mc.txt` as examples (or you can use these files).  
Then, execute the script with your datsets text file as a command-line argument: 
```
python crab-submit-multiple.py <datasets.txt> <year>
```
If no command-line argument is given, a default file & year will be chosen.
This script will loop through the different datasets and submit crab jobs for each one.
_Note: the text file needs to be written such that each sample has a 'nickname' that can be used to create the crab directory.  The example text files include this structure -- replicate it if you write your own file! Example:_
`ttbar /TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM`
_`ttbar` is the nickname, and `/TT.../MINIAODSIM` is the sample_

### Single Submission
_The following method is not recommended.  Instead, if you want to only submit one sample, use the method described in "Multiple Submission" and fill your text file with only 1 dataset._

Edit `crab_*.py` with a new dataset name to reflect the sample you are processing, e.g., `crab_SMttbar.py`.  
Make sure to the output directory is pointed to `/store/group/lpctop/ttbarAC/ttbarAC_skim/`.  
If you are processing data, modify the argument `config.JobType.pyCfgParams = ['isMC=True']` to `config.JobType.pyCfgParams = ['isMC=False']`.  
To submit the CRAB jobs, enter the command:
```
crab submit -c crab_*.py --dryrun
crab proceed
```
You can also submit the crab jobs directly by removing `--dryrun`.

For either option, monitor the CRAB jobs and resubmit any that fail.


## Previous setup (not recommended)

### EDM Production:

Instructions for producing EDM Ntuples with BEST tagger inputs and then flat trees for analysis.

```
cd ttbarAC_skim/ttbarAC_skim/test
voms-proxy-init -voms cms
source /cvmfs/cms.cern.ch/crab3/crab.csh
```

To submit crab jobs, the script `runSkim_data.py` (`runSkim.py` has been modified too much) is processed for MC (data) using the `cmsRun` executable.

First, edit `crab_*.py` with a new dataset name to reflect the sample you are processing, e.g., `crab_SMttbar.py`.  
Make sure to the output directory is pointed to `/store/group/lpctop/edmNtuples/`.  
To submit the CRAB jobs, enter the command:
```
crab submit -c crab.py --dryrun
```
either remove `--dryrun` to submit the jobs, or enter `crab proceed` after running the above command.  
Monitor the CRAB jobs and resubmit any failed jobs.


### Flat Ntuple Production

After crab jobs are finished, make a text file with the filenames (starting with `/store/group/...`):

```
xrdfsls -u /store/group/lpctop/ttbarAC/edmNtuples/.../ | grep "ana" > textfile.txt
```
If there are multiple EDM files in this directory, change `ana` to a specific keyword in the desired files.  
Then, setup the condor scripts to submit jobs to the LPC condor system that process the EDM ntuples into flat ntuples.  
The condor scripts to check are `submit.sh`, `analyze.sh`, and `treeMaker_fwlite.py`.  
There shouldn't be any changes necessary in `treeMaker_fwlite.py` or `submit.sh`, but `analyze.sh` may need changes for the following two lines:
```
output_path="/store/group/lpctop/ttbarAC/flatNtuples/"
...
python treeMaker_fwlite.py --files ${file} --isMC 1
```
where `output_path` may need to be changed, and `isMC` needs to be 0 for data samples.  

To submit the jobs, enter:

```
./submit.sh files.txt <process>
```

Where `<process>` is the name of the process, e.g., `SMttbar`.  
Note: 
> this string gets added to the ROOT output files, 
> so you use that same string later with `do_hadd.sh` 
> to grab only the files from a certain process and `hadd` them

Monitor these jobs using `condor_q ${USER}` or using [landscape](https://landscape.fnal.gov/lpc/dashboard/db/lpc-summary?orgId=1).  
Once the jobs are finished, flat ntuples are available for analysis using [CyMiniAna](https://github.com/cms-ttbarAC/CyMiniAna).

## Questions
Contact the authors.
