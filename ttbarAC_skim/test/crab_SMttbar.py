from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'ttbarAC_ttbar'
config.General.workArea = 'crab_ttbar'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runSkim.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/MINIAODSIM' 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False

#config.Data.outLFNDirBase = '/store/user/lpctlbsm/pilot/'
config.Data.outLFNDirBase = '/store/group/lpctop/ttbarAC/edmNtuples/'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
