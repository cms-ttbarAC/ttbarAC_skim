from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'ttbarAC_dataH'
config.General.workArea = 'crab_dataH'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runSkim.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['isMC=False']

config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2016H-03Feb2017_ver2-v1/MINIAOD'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.publication = False
config.Data.lumiMask = 'goldenJSON_2016.txt'
config.Data.outLFNDirBase = '/store/group/lpctlbsm/'


config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
