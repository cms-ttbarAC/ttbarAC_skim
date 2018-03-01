#!/bin/bash

#do hadd of the vlq files
eval `scramv1 runtime -sh`

eospath="/store/user/pilot/ttbarACFiles/"

for argument in "$@"
do
    echo ${argument};
    hadd -f root://cmseos.fnal.gov/${eospath}/ttbarAC_${argument}_merged.root `xrdfs root://cmseos.fnal.gov ls -u ${eospath} | grep ${argument} | grep "ttbarAC"`; 
done
