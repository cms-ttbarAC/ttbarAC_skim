import FWCore.ParameterSet.Config as cms

process = cms.Process("ttbarACskim")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/00000/E602EB22-B93A-E611-9FAC-0025904C6564.root'
#	'/store/mc/RunIISpring16MiniAODv2/ZprimeToTT_M-3000_W-30_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/70000/1A405524-3E3D-E611-B103-047D7BD6DDB2.root'
	)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.BESTProducer = cms.EDProducer('BESTProducer',
	pdgIDforMatch = cms.int32(6),
	NNtargetX = cms.int32(1),
	NNtargetY = cms.int32(1),
	isMC = cms.int32(1),
	doMatch = cms.int32(0)

)

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
    cut = cms.string('pt > 50.0 && abs(eta) < 2.4')
)

process.selectedMET = cms.EDFilter('PATMETSelector',
    src = cms.InputTag('slimmedMETs'),
    cut = cms.string('pt > -999.9'),
)

#process.selectedGenParticles = cms.EDFilter("CandPtrSelector",
#    src = cms.InputTag('packedGenParticles'),
#    cut = cms.string( 'abs(pdgId) == 6 || abs(pdgId) ==  15 || abs(pdgId) == 23 || abs(pdgId) ==  24 || abs(pdgId) == 25'),
#    filter = cms.bool(False)
#)

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
								      #, 'keep *_goodPatJetsCATopTagPF_*_*'
                     
			                                                 #, 'keep recoPFJets_*_*_*'
                                                                      ) 
                               )


process.outpath = cms.EndPath(process.out)

process.p = cms.Path(process.BESTProducer*
		     process.selectedMuons*
		     process.selectedElectrons*
		     process.selectedAK4Jets*
		     process.selectedMET*
		     process.selectedGenParticles
)


