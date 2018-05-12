"""
Created:      --
Last Updated:  2 April 2018

Justin Pilot
UC Davis

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
---

Primary script for running ttbarAC_skim
over MiniAOD samples

To run:
  cmsRun runSkim.py
"""
import json
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

options = VarParsing('analysis')
options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Simulated data sample" )
options.register('sampleName',"",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of simulated/real data sample" )
options.parseArguments()


process = cms.Process("ttbarACskim")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        '/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/000C6E52-8BEC-E611-B3FF-0025905C42FE.root'
#        '/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/00000/E602EB22-B93A-E611-9FAC-0025904C6564.root'
#        '/store/mc/RunIISpring16MiniAODv2/ZprimeToTT_M-3000_W-30_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/1A405524-3E3D-E611-B103-047D7BD6DDB2.root'
#         '/store/user/demarley/ttbar/example.root'
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'
	)
)
process.TFileService = cms.Service("TFileService",fileName = cms.string("ttbarAC_outtree.root"))

## VID Electrons
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


## BEST
process.BESTProducer = cms.EDProducer('BESTProducer',
	pdgIDforMatch = cms.int32(6),
	NNtargetX = cms.int32(1),
	NNtargetY = cms.int32(1),
	isMC = cms.int32(options.isMC),
	doMatch = cms.int32(0)
)

## PHYSICS OBJECTS
process.selectedMuons = cms.EDFilter('PATMuonSelector',
    src = cms.InputTag('slimmedMuons'),
    cut = cms.string('pt > 40.0 && abs(eta) < 2.4')
)

process.selectedElectrons = cms.EDFilter('PATElectronSelector',
    src = cms.InputTag('slimmedElectrons'),
    cut = cms.string('pt > 40.0 && abs(eta) < 2.4')
)

process.selectedAK4Jets = cms.EDFilter('PATJetSelector',
    src = cms.InputTag('slimmedJets'),
    cut = cms.string('pt > 30.0 && abs(eta) < 2.4')
)
"""
updateJetCollection(
    process,
    labelName = "DeepCSV",
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = [
      'pfDeepCSVJetTags:probb',
      'pfDeepCSVJetTags:probbb',
      'pfDeepCSVJetTags:probc',
      'pfDeepCSVJetTags:probcc',
      'pfDeepCSVJetTags:probudsg',
      ] ## to add discriminators
)
"""
process.selectedMET = cms.EDFilter('PATMETSelector',
    src = cms.InputTag('slimmedMETs'),
    cut = cms.string('pt > -999.9'),
)

## GENERATOR PARTICLES
process.selectedGenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
        'drop *',
        'keep status == 3',
        'keep status >= 20 && status <= 100',
        'keep abs(pdgId) == 6 && status >= 20 && status <= 40',
        'keep abs(pdgId) >= 1 && abs(pdgId) <= 5 && status <= 100',
        'keep abs(pdgId) >= 11 && abs(pdgId) <= 18 && status <= 100',
        'keep abs(pdgId) == 23 && status >= 20 && status <= 40',
        'keep abs(pdgId) == 24 && status >= 20 && status <= 100',
        'keep abs(pdgId) == 25 && status >= 20 && status <= 40',
        'keep numberOfMothers() == 1 && abs(mother().pdgId()) == 6 && status >= 20 && status <= 40',
        'keep numberOfMothers() >= 1 && abs(mother().pdgId()) == 24 && status >= 20 && status <= 100',
    )
)


## EVENT SAVER FLAT NTUPLE
## Get the sampleName (primary dataset)
if not options.sampleName:
    try:
        sampleName = process.source.fileNames[0]
    except IndexError:
        sampleName = options.sampleName
else:
    sampleName = options.sampleName

if sampleName.startswith("/store/"):
    sampleName = sampleName.split("/")[4]  # try to get the primary dataset

process.tree = cms.EDAnalyzer("EventSaverFlatNtuple",
    isMC = cms.bool(options.isMC),
    sampleName = cms.string(sampleName),
    metadataFile = cms.string("metadataFile.txt"),
    elIdFullInfoMap_Loose  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    elIdFullInfoMap_Medium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    elIdFullInfoMap_Tight  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    #elIdFullInfoMap_HEEP   = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"), 
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

## PROCESS
#  -- 'process.selectedElectrons*' removed -- need original collection to use VID
if options.isMC:
    process.p = cms.Path(
        process.BESTProducer*
        process.selectedMuons*
        process.selectedAK4Jets*
        process.selectedMET*
        process.selectedGenParticles*
        process.egmGsfElectronIDSequence*
        process.tree
    )
else:
    process.p = cms.Path(
        process.BESTProducer*
        process.selectedMuons*
        process.selectedAK4Jets*
        process.selectedMET*
        process.egmGsfElectronIDSequence*
        process.tree
    )


## THE END 

"""
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("ana_out.root"),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*BEST*_*_*',
								      'keep *_*_*_*ttbarAC*',
								      'drop *_*_*pfCand*_*',
								      'drop *_*triggerResult*_*_*',
                        					      'drop *_*_*genJets*_*',
			                                              'drop *_*_*tagInfos*_*',
			                                              'drop *_*_*caloTowers*_*'
                                                                      ) 
                               )
process.outpath = cms.EndPath(process.out)
"""
