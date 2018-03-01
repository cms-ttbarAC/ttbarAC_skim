from optparse import OptionParser
import numpy as np
from sklearn import svm, metrics, preprocessing
from sklearn.externals import joblib
import warnings
import os
warnings.simplefilter("ignore")


parser = OptionParser()

parser.add_option('--files', type='string', action='store', dest='files', help='Input Files')
parser.add_option('--outname', type='string', action='store',default='ttbarAC_outtree.root', dest='outname',help='Name of output file')
parser.add_option('--maxevents', type='int', action='store',default=-1,dest='maxevents',help='Number of events to run. -1 is all events')


(options, args) = parser.parse_args()
argv = []

import ROOT
import sys, copy
from array import array
from DataFormats.FWLite import Events, Handle


f = ROOT.TFile(options.outname, "RECREATE")
f.cd()

eventTree = ROOT.TTree('eventVars', 'eventVars')
maxObjects = 5
vBESTprob_t = ROOT.vector('float')()
vBESTprob_W = ROOT.vector('float')()
vBESTprob_Z = ROOT.vector('float')()
vBESTprob_H = ROOT.vector('float')()
vBESTprob_j = ROOT.vector('float')()
vBESTclass = ROOT.vector('float')()
eventTree.Branch( 'BESTProb_t', vBESTprob_t)
eventTree.Branch( 'BESTProb_W', vBESTprob_W)
eventTree.Branch( 'BESTProb_Z', vBESTprob_Z)
eventTree.Branch( 'BESTProb_H', vBESTprob_H)
eventTree.Branch( 'BESTProb_j', vBESTprob_j)
eventTree.Branch( 'BESTclass', vBESTclass)
vAK8pt = ROOT.vector('float')()
vAK8eta = ROOT.vector('float')()
vAK8phi = ROOT.vector('float')()
vAK8mass = ROOT.vector('float')()
vAK8SDmass = ROOT.vector('float')()
vAK8tau1 = ROOT.vector('float')()
vAK8tau2 = ROOT.vector('float')()
vAK8tau3 = ROOT.vector('float')()
vAK8charge = ROOT.vector('float')()
vAK8bDiscSubjet1 = ROOT.vector('float')()
vAK8bDiscSubjet2 = ROOT.vector('float')()
vAK8ChargeSubjet1 = ROOT.vector('float')()
vAK8ChargeSubjet2 = ROOT.vector('float')()
eventTree.Branch( 'AK8pt', vAK8pt)
eventTree.Branch( 'AK8eta', vAK8eta)
eventTree.Branch( 'AK8phi', vAK8phi)
eventTree.Branch( 'AK8mass', vAK8mass)
eventTree.Branch( 'AK8SDmass', vAK8SDmass)
eventTree.Branch( 'AK8tau1', vAK8tau1)
eventTree.Branch( 'AK8tau2', vAK8tau2)
eventTree.Branch( 'AK8tau3', vAK8tau3)
eventTree.Branch( 'AK8charge', vAK8charge)
eventTree.Branch( 'AK8bDiscSubjet1', vAK8bDiscSubjet1)
eventTree.Branch( 'AK8bDiscSubjet2', vAK8bDiscSubjet2)
eventTree.Branch( 'AK8ChargeSubjet1', vAK8ChargeSubjet1)
eventTree.Branch( 'AK8ChargeSubjet2', vAK8ChargeSubjet2)

HTak8 = array('f', [-999.])
eventTree.Branch( 'HTak8', HTak8, 'HTak8/F')

vAK4pt = ROOT.vector('float')()
vAK4eta = ROOT.vector('float')()
vAK4phi = ROOT.vector('float')()
vAK4mass = ROOT.vector('float')()
vAK4bDisc = ROOT.vector('float')()
eventTree.Branch( 'AK4pt', vAK4pt)
eventTree.Branch( 'AK4eta', vAK4eta)
eventTree.Branch( 'AK4phi', vAK4phi)
eventTree.Branch( 'AK4mass', vAK4mass)
eventTree.Branch( 'AK4bDisc', vAK4bDisc)

vELpt = ROOT.vector('float')()
vELeta = ROOT.vector('float')()
vELphi = ROOT.vector('float')()
vELenergy = ROOT.vector('float')()
vELcharge = ROOT.vector('float')()
vELiso = ROOT.vector('float')()
vELid = ROOT.vector('float')()
eventTree.Branch( 'ELpt', vELpt)
eventTree.Branch( 'ELeta', vELeta)
eventTree.Branch( 'ELphi', vELphi)
eventTree.Branch( 'ELenergy', vELenergy)
eventTree.Branch( 'ELcharge', vELcharge)
eventTree.Branch( 'ELiso', vELiso)
eventTree.Branch( 'ELid', vELid)
vMUpt = ROOT.vector('float')()
vMUeta = ROOT.vector('float')()
vMUphi = ROOT.vector('float')()
vMUenergy = ROOT.vector('float')()
vMUcharge = ROOT.vector('float')()
vMUlooseID = ROOT.vector('float')()
vMUcorrIso = ROOT.vector('float')()
eventTree.Branch( 'MUpt', vMUpt)
eventTree.Branch( 'MUeta', vMUeta)
eventTree.Branch( 'MUphi', vMUphi)
eventTree.Branch( 'MUenergy', vMUenergy)
eventTree.Branch( 'MUcharge', vMUcharge)
eventTree.Branch( 'MUlooseID', vMUlooseID)
eventTree.Branch( 'MUcorrIso', vMUcorrIso)

METpt = array('f', [-999.])
METphi = array('f', [-999.])
eventTree.Branch( 'METpt', METpt, 'METpt/F')
eventTree.Branch( 'METphi', METphi, 'METphi/F')


eventNum = array('i', [-999])
runNum = array('i', [-999])
lumiNum = array('i', [-999])
eventTree.Branch( 'eventNum', eventNum, 'eventNum/i')
eventTree.Branch( 'lumiNum', lumiNum, 'lumiNum/i')
eventTree.Branch( 'runNum', runNum, 'runNum/i')



jetsHandle = Handle("std::vector<pat::Jet>")
jetsLabel = ("BESTProducer", "savedJets", "ttbarACskim")
AK4jetsHandle = Handle("std::vector<pat::Jet>")
AK4jetsLabel = ("selectedAK4Jets", "", "ttbarACskim")
muonsHandle = Handle("std::vector<pat::Muon>")
muonsLabel = ("selectedMuons", "", "ttbarACskim")
electronsHandle = Handle("std::vector<pat::Electron>")
electronsLabel = ("selectedElectrons", "", "ttbarACskim")
metLabel = ("selectedMET", "", "ttbarACskim")
metHandle = Handle("std::vector<pat::MET>")


nnLabels = []
nnLabels.append( ("BESTProducer", "FWmoment1H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment1W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment1Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment1top", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment2H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment2W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment2Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment2top", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment3H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment3W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment3Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment3top", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment4H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment4W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment4Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "FWmoment4top", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "SDmass", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "aplanarityH", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "aplanarityW", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "aplanarityZ", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "aplanaritytop", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "bDisc", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "bDisc1", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "bDisc2", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "et", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "eta", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "isotropyH", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "isotropyW", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "isotropyZ", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "isotropytop", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "q", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "qsubjet0", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "qsubjet1", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sphericityH", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sphericityW", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sphericityZ", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sphericitytop", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPH", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPW", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPZ", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPtop", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPzH", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPzW", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPzZ", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "sumPztop", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "tau21", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "tau32", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "thrustH", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "thrustW", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "thrustZ", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "thrusttop", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m12H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m23H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m13H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m1234H", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m12W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m23W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m13W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m1234W", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m12Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m23Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m13Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m1234Z", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m12top", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m23top", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m13top", "ttbarACskim") )
nnLabels.append( ("BESTProducer", "m1234top", "ttbarACskim") )
nnHandles = [ Handle("vector<float>") ] * len(nnLabels)




filelist = options.files
files= []
nevents = 0
#s = options.files
s = 'root://cmsxrootd.fnal.gov/' + options.files
files.append(s)
print 'Added ' + s

mlp = joblib.load('BEST_mlp.pkl')
scaler = joblib.load('BEST_scaler.pkl')



nevents_file = 0
for ifile in files:
	print 'Processing', ifile

	if options.maxevents > 0 and nevents > options.maxevents :
        	break

	events = Events( ifile )


	products = {}


	for event in events:

		vBESTprob_t.clear()
		vBESTprob_W.clear()
		vBESTprob_Z.clear()
		vBESTprob_H.clear()
		vBESTprob_j.clear()
		vBESTclass.clear()	
		vAK8pt.clear()
		vAK8eta.clear()
		vAK8phi.clear()
		vAK8mass.clear()
		vAK8SDmass.clear()
		vAK8tau1.clear()
		vAK8tau2.clear()
		vAK8tau3.clear()
		vAK8charge.clear()
		vAK8bDiscSubjet1.clear()
		vAK8bDiscSubjet2.clear()
		vAK8ChargeSubjet1.clear()
		vAK8ChargeSubjet2.clear()
		vAK4pt.clear()
		vAK4eta.clear()
		vAK4phi.clear()
		vAK4mass.clear()
		vAK4bDisc.clear()
		vELpt.clear()
		vELeta.clear()
		vELphi.clear()
		vELenergy.clear()
		vELiso.clear()
		vELid.clear()
		vELcharge.clear()
		vMUpt.clear()
		vMUeta.clear()
		vMUphi.clear()
		vMUenergy.clear()
		vMUcharge.clear()
		vMUlooseID.clear()
	
		eventNum[0] = event.eventAuxiliary().event()
		lumiNum[0] = event.eventAuxiliary().luminosityBlock()
		runNum[0] = event.eventAuxiliary().run()
	
		if options.maxevents > 0 and nevents > options.maxevents :
        		break
		nevents += 1
		nevents_file += 1
		if (nevents_file % 1000 == 0):
			print nevents, ' Completed'


		event.getByLabel(jetsLabel, jetsHandle)
		if not jetsHandle.isValid():
			continue
		jets = jetsHandle.product()
		event.getByLabel(AK4jetsLabel, AK4jetsHandle)
		AK4jets = AK4jetsHandle.product()
		
		event.getByLabel(muonsLabel, muonsHandle)
		event.getByLabel(electronsLabel, electronsHandle)
		event.getByLabel(metLabel, metHandle)
		electrons = electronsHandle.product()
		muons = muonsHandle.product()
		mets = metHandle.product()

		METpt[0] = mets[0].pt()
		METphi[0] = mets[0].phi()

		for i in xrange(len(nnHandles)):
			event.getByLabel((nnLabels[i]), nnHandles[i])
			products[nnLabels[i][1]] = nnHandles[i].product() 

		ht_ak8 = 0.0

		if len(jets) < 1:
			continue
	
		for j in xrange(len(jets)):

			nhf = jets[j].neutralHadronEnergy() / jets[j].energy()
                	nef = jets[j].neutralEmEnergy() / jets[j].energy()
                	chf = jets[j].chargedHadronEnergy() / jets[j].energy()
                	cef = jets[j].chargedEmEnergy() / jets[j].energy()
                	nconstituents = jets[j].numberOfDaughters()
                	nch = jets[j].chargedMultiplicity()
                	goodJet = \
                  		nhf < 0.99 and \
                  		nef < 0.99 and \
                  		chf > 0.00 and \
                  		cef < 0.99 and \
                  		nconstituents > 1 and \
                  		nch > 0

                	if not goodJet :
				continue




			pzOverp_top = products['sumPztop'][j] / (products['sumPtop'][j] + 0.01) 
			pzOverp_W = products['sumPzW'][j] / (products['sumPW'][j] + 0.01) 
			pzOverp_Z = products['sumPzZ'][j] / (products['sumPZ'][j] + 0.01)
			pzOverp_H = products['sumPzH'][j] / (products['sumPH'][j] + 0.01)


			nnArray = np.array( [  products['SDmass'][j], products['tau32'][j], products['tau21'][j], products['FWmoment1top'][j], products['FWmoment2top'][j], products['FWmoment3top'][j], products['FWmoment4top'][j], products['isotropytop'][j], products['aplanaritytop'][j], products['sphericitytop'][j], products['thrusttop'][j], products['FWmoment1W'][j], products['FWmoment2W'][j], products['FWmoment3W'][j], products['FWmoment4W'][j], products['isotropyW'][j], products['aplanarityW'][j], products['sphericityW'][j], products['thrustW'][j], products['FWmoment1Z'][j], products['FWmoment2Z'][j], products['FWmoment3Z'][j], products['FWmoment4Z'][j], products['isotropyZ'][j], products['aplanarityZ'][j], products['sphericityZ'][j], products['thrustZ'][j], products['FWmoment1H'][j], products['FWmoment2H'][j], products['FWmoment3H'][j], products['FWmoment4H'][j], products['isotropyH'][j], products['aplanarityH'][j], products['sphericityH'][j], products['thrustH'][j], products['bDisc'][j], products['bDisc1'][j], products['bDisc2'][j], products['q'][j], products['m12W'][j], products['m13W'][j], products['m23W'][j], products['m1234W'][j], products['m12Z'][j], products['m13Z'][j], products['m23Z'][j], products['m1234Z'][j], products['m12top'][j], products['m13top'][j], products['m23top'][j], products['m1234top'][j], products['m12H'][j], products['m13H'][j], products['m23H'][j], products['m1234H'][j], pzOverp_top, pzOverp_W, pzOverp_Z, pzOverp_H  ])

			if np.isnan(np.min(nnArray)):
				continue

		 	vAK8pt.push_back(jets[j].pt())
			vAK8eta.push_back(jets[j].eta())
			vAK8phi.push_back(jets[j].phi())
			vAK8mass.push_back(jets[j].mass())
			vAK8SDmass.push_back( products['SDmass'][j] )

			vAK8tau1.push_back( jets[j].userFloat('NjettinessAK8:tau1') )
			vAK8tau2.push_back( jets[j].userFloat('NjettinessAK8:tau2') )
			vAK8tau3.push_back( jets[j].userFloat('NjettinessAK8:tau3') )
			vAK8charge.push_back(  products['q'][j] )
		
			vAK8bDiscSubjet1.push_back( products['bDisc1'][j] )
			vAK8bDiscSubjet2.push_back( products['bDisc2'][j] )
			vAK8ChargeSubjet1.push_back( products['qsubjet0'][j] )
			vAK8ChargeSubjet2.push_back( products['qsubjet1'][j] )

			ht_ak8 += jets[j].pt()

			nnArray_transformed = scaler.transform( nnArray )



			best_bin = mlp.predict( nnArray_transformed )

			best_probs = mlp.predict_proba( nnArray_transformed )


			vBESTclass.push_back(best_bin[0])	
			vBESTprob_t.push_back(best_probs[0][0])
			vBESTprob_W.push_back(best_probs[0][1])
			vBESTprob_Z.push_back(best_probs[0][2])
			vBESTprob_H.push_back(best_probs[0][3])
			vBESTprob_j.push_back(best_probs[0][4])
			
		HTak8[0] = ht_ak8
		j = 0
		for j in xrange(len(AK4jets)):


			vAK4pt.push_back(AK4jets[j].pt())
			vAK4eta.push_back(AK4jets[j].eta())
			vAK4phi.push_back(AK4jets[j].phi())
			vAK4mass.push_back(AK4jets[j].mass())
			vAK4bDisc.push_back(AK4jets[j].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"))


		for e in xrange(len(electrons)):

			
			vELpt.push_back(electrons[e].pt())
			vELeta.push_back(electrons[e].eta())
			vELphi.push_back(electrons[e].phi())
			vELenergy.push_back(electrons[e].energy())
			vELcharge.push_back(electrons[e].charge())
			vELiso.push_back( (electrons[e].trackIso() + electrons[e].caloIso()) / electrons[e].pt() )
			vELid.push_back( electrons[e].userFloat('ElectronMVAEstimatorRun2Spring15Trig25nsV1Values') )

		for m in xrange(len(muons)):
	
			vMUpt.push_back(muons[m].pt())
			vMUeta.push_back(muons[m].eta())
			vMUphi.push_back(muons[m].phi())
			vMUenergy.push_back(muons[m].energy())
			vMUcharge.push_back(muons[m].charge())
			vMUlooseID.push_back(muons[m].isLooseMuon())
			chPt = muons[m].pfIsolationR04().sumChargedHadronPt
			nhEt = muons[m].pfIsolationR04().sumNeutralHadronEt
			phEt = muons[m].pfIsolationR04().sumPhotonEt
			puPt = muons[m].pfIsolationR04().sumPUPt

			corrCombRelIso = (chPt + max(0.0, nhEt + phEt - 0.5*puPt) ) / muons[m].pt()
			vMUcorrIso.push_back(corrCombRelIso)


		eventTree.Fill()



f.cd()
f.Write()
f.Close()
