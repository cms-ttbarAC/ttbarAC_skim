"""
Created 3 April

Script to submit multiple CRAB jobs
 https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial#1_CRAB_configuration_file_to_run

To run:
  cd test/
  python crab-submit-multiple.py <datasets>

where <datasets> is a text file that contains the datasets you want to process.
See 'test/crab-datasets-mc.txt' for an example.
If no argument is provided, a default option is selected
"""
import os
import sys
from multiprocessing import Process
import Analysis.CyMiniAna.util as util



def main(input_datasets="crab-datasets.txt"):

    from CRABClient.UserUtilities import config
    config = config()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    input_datasets = util.file2list(input_datasets)


    def submit(config):
        try:
            crabCommand('submit', config=config, dryrun=True)   # add 'dryrun=False' as argument for testing
            print ' Executed crabCommand() '
        except HTTPException, hte:
            print ' ERROR :: Cannot execute command! '
            print hte.headers


    for id in input_datasets:
        if (not id or id.startswith('#')): continue

        name,dataset = id.split(" ")
        isMC = True
        if name.startswith("SingleElectron") or name.startswith("SingleMuon") or name.startswith("JetHT"):
            isMC = False

        primary_dataset = dataset.split('/')[1]

        print '  --> Added {0}'.format(id)
        print '      - isMC = {0}'.format(isMC)
        print '      - primary dataset = {0}'.format(primary_dataset)

        # General
        config.General.requestName = 'ttbarAC_'+name
        config.General.workArea    = 'crab_'+name
        config.General.transferOutputs = True
        config.General.transferLogs    = True

        # JobType
        config.JobType.pluginName  = 'Analysis'
        config.JobType.psetName    = 'runSkim.py'
        config.JobType.allowUndistributedCMSSW = True
        config.JobType.inputFiles = ['BEST_mlp.json','metadataFile.txt','JERDatabase','JECDatabase']
        if isMC:
            config.JobType.pyCfgParams = ['isMC=True','sampleName={0}'.format(primary_dataset)]
        else:
            config.JobType.pyCfgParams = ['isMC=False','sampleName={0}'.format(primary_dataset)]

        # Data
        config.Data.splitting     = 'Automatic'
        #config.Data.unitsPerJob   = 10 if isMC else 3
        config.Data.outLFNDirBase = '/store/group/lpctop/ttbarAC/ttbarAC_skim_v0.3/'
        config.Data.publication   = False
        config.Data.inputDataset  = dataset
        if not isMC: config.Data.lumiMask = 'goldenJSON_2016.txt'

        # Site
        config.Site.storageSite   = "T3_US_FNALLPC"

        print '\n Configuration :'
        print config
        try :
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
        except :
            print ' ERROR :: Not submitted!'


if __name__=='__main__':
    if len(sys.argv)>1:
        main(sys.argv[1])    # pass datasets file as command line argument
    else:
        main("crab-datasets-test.txt")

## THE END ##
