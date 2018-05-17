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



def main(input_datasets="crab-datasets.txt",year='2016'):

    from CRABClient.UserUtilities import config
    config = config()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    input_datasets = util.file2list(input_datasets)


    def submit(config):
        try:
            crabCommand('submit', config=config, dryrun=False)   # add 'dryrun=False' as argument for testing
            print ' Executed crabCommand() '
        except HTTPException, hte:
            print ' ERROR :: Cannot execute command! '
            print hte.headers


    for id in input_datasets:
        if (not id or id.startswith('#')): continue

        name,dataset = id.split(" ")
        isMC = not (name.startswith("SingleElectron") or name.startswith("SingleMuon") or name.startswith("JetHT"))  # check data names

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
        config.JobType.inputFiles  = ['BEST_mlp.json','metadataFile.txt','JERDatabase','JECDatabase']
        config.JobType.pyCfgParams = ['isMC={0}'.format(isMC),'sampleName={0}'.format(primary_dataset),'year={0}'.format(year)]
        config.JobType.allowUndistributedCMSSW = True

        # Data
        #config.Data.splitting     = 'Automatic'
        config.Data.splitting     = 'FileBased'
        config.Data.unitsPerJob   = 5 if isMC else 3
        config.Data.outLFNDirBase = '/store/group/lpctop/ttbarAC/ttbarAC_skim_v0.4/'
        config.Data.publication   = False
        config.Data.inputDataset  = dataset
        if not isMC: 
            config.Data.lumiMask = 'goldenJSON_{0}.txt'.format(year)

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
    try:
        main(sys.argv[1],sys.argv[2])    # pass datasets file and year as command line arguments
    except IndexError:
        print " Not enough arguments! "
        try:
            main(sys.argv[1])
        except IndexError:
            main("crab-datasets-test.txt")

## THE END ##
