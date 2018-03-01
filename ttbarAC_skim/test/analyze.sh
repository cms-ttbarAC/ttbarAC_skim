#!/bin/bash


file=$1
sample=$2
process=$3
output_path="/store/group/lpctop/ttbarAC/flatNtuples/"

cd ${_CONDOR_SCRATCH_DIR}
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a tcsh script, use .csh instead of .sh
export SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 project CMSSW CMSSW_9_0_1`
cp treeMaker_fwlite.py BEST_mlp.pkl BEST_scaler.pkl CMSSW_9_0_1/src/.
cd CMSSW_9_0_1/src/
eval `scramv1 runtime -sh`
ls
echo ${file}


python treeMaker_fwlite.py --files ${file}


mv ttbarAC_outtree.root ttbarAC_tree_${sample}_${process}.root
xrdcp -f ttbarAC_tree_${sample}_${process}.root root://cmseos.fnal.gov/${output_path}
