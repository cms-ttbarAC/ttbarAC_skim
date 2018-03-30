"""
Created:        --
Last Updated:   22 March 2018

Justin Pilot
UC Davis
-----

Convert EDM Ntuples (MiniAOD) into flat ntuples

Use with batch scripts (analyze.sh)
To test locally:
$ python treeMaker_fwlite.py --files <name_of_file> --maxevents N --outname <output_name> --isMC <is_mc>
"""
from optparse import OptionParser
import numpy as np
from sklearn import svm, metrics, preprocessing
from sklearn.externals import joblib
import warnings
import os
warnings.simplefilter("ignore")

# parse command-line arguments
parser = OptionParser()

parser.add_option('--files', 
          type='string', action='store', 
          dest='files', 
          help='Input File')
parser.add_option('--outname', 
          type='string', action='store',
          default='ttbarAC_outtree.root',
          dest='outname',
          help='Name of output file')
parser.add_option('--maxevents', 
          type='int',    action='store',
          default=-1,
          dest='maxevents',
          help='Number of events to run. -1 is all events')
parser.add_option('--isMC', 
          type='int',    action='store', 
          default=1, 
          dest='isMC',
          help='is it MC?')

(options, args) = parser.parse_args()
argv = []


import ROOT
import sys, copy
from array import array
from DataFormats.FWLite import Events, Handle
from RecoEgamma.EgammaTools import *


def ak8Quality(jet):
    """Quality cut on AK8 jet"""
    energy = jet.energy()
    nhf = jet.neutralHadronEnergy() / jet.energy()
    nef = jet.neutralEmEnergy() / jet.energy()
    chf = jet.chargedHadronEnergy() / jet.energy()
    cef = jet.chargedEmEnergy() / jet.energy()
    nconstituents = jet.numberOfDaughters()
    nch = jet.chargedMultiplicity()
    goodJet = \
        nhf < 0.99 and \
        nef < 0.99 and \
        chf > 0.00 and \
        cef < 0.99 and \
        nconstituents > 1 and \
        nch > 0

    return goodJet


def checkTopDecay(daughter):
    """Check the decay type of top quark by checking W decay products"""
    w_child = daughter.daughter(0)
    if abs(w_child.pdgId())==24 and w_child.numberOfDaughters()>0:
        # w decays to itself, check its decays
        mode = w_child.daughter(0).pdgId()
    else:
        # w decays normally
        mode = w_child.pdgId()
    return (abs(mode)<10)


def effectiveArea(eta):
    """Effective area for electron ID
       - The following values refer to EA for cone 0.3 and fixedGridRhoFastjetAll. 
       - They are valid for electrons only, different EA are available for muons.
    """
    effArea = 0.;
    if(np.fabs(eta)>0.0 and np.fabs(eta)<=1.0): effArea = 0.1752;
    elif(np.fabs(eta)>1.0 and np.fabs(eta)<=1.479): effArea = 0.1862;
    elif(np.fabs(eta)>1.479 and np.fabs(eta)<=2.0): effArea = 0.1411;
    elif(np.fabs(eta)>2.0 and np.fabs(eta)<=2.2): effArea = 0.1534;
    elif(np.fabs(eta)>2.2 and np.fabs(eta)<=2.3): effArea = 0.1903;
    elif(np.fabs(eta)>2.3 and np.fabs(eta)<=2.4): effArea = 0.2243;
    elif(np.fabs(eta)>2.4 and np.fabs(eta)<=2.5): effArea = 0.2687;

    return effArea;


def passIDWP(isEB, dEtaIn, dPhiIn, full5x5, hoe, relIso, ooemoop, missHits):
    """Check the electron values against working point thresholds
       Isolation cut is not relative isolation -- no pT used.

       Values accessed on 29 March 2018 from 
       https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
    """
    passVETO  = False
    passLOOSE = False
    passMED   = False
    passTIGHT = False

    dEtaIn  = np.fabs(dEtaIn)
    dPhiIn  = np.fabs(dPhiIn)
    ooemoop = np.fabs(ooemoop)

    if(isEB):
        passVETO  = (dEtaIn < 0.00749) and (dPhiIn < 0.228)  and (full5x5 < 0.0115)  and (hoe < 0.356)  and (relIso < 0.175)  and (ooemoop < 0.299)  and (missHits <= 2)
        passLOOSE = (dEtaIn < 0.00477) and (dPhiIn < 0.222)  and (full5x5 < 0.011)   and (hoe < 0.298)  and (relIso < 0.0994) and (ooemoop < 0.241)  and (missHits <= 1)
        passMED   = (dEtaIn < 0.00311) and (dPhiIn < 0.103)  and (full5x5 < 0.00998) and (hoe < 0.253)  and (relIso < 0.0695) and (ooemoop < 0.134)  and (missHits <= 1)
        passTIGHT = (dEtaIn < 0.00308) and (dPhiIn < 0.0816) and (full5x5 < 0.00998) and (hoe < 0.0414) and (relIso < 0.0588) and (ooemoop < 0.0129) and (missHits <= 1)
    else:
        passVETO  = (dEtaIn < 0.00895) and (dPhiIn < 0.213)  and (full5x5 < 0.037)  and (hoe < 0.211)  and (relIso < 0.159)  and (ooemoop < 0.15)   and (missHits <= 3)
        passLOOSE = (dEtaIn < 0.00868) and (dPhiIn < 0.213)  and (full5x5 < 0.0314) and (hoe < 0.101)  and (relIso < 0.107)  and (ooemoop < 0.14)   and (missHits <= 1)
        passMED   = (dEtaIn < 0.00609) and (dPhiIn < 0.045)  and (full5x5 < 0.0298) and (hoe < 0.0878) and (relIso < 0.0821) and (ooemoop < 0.13)   and (missHits <= 1)
        passTIGHT = (dEtaIn < 0.00605) and (dPhiIn < 0.0394) and (full5x5 < 0.0292) and (hoe < 0.0641) and (relIso < 0.0571) and (ooemoop < 0.0129) and (missHits <= 1)

    return {'isVETO':passVETO,'isLoose':passLOOSE,'isMedium':passMED,'isTight':passTIGHT}


def electronID(el,rho=None,conversions=None,beamspot=None):
    """Get the ID working points for electrons
       https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
       future iterations will use VID: https://twiki.cern.ch/twiki/bin/view/Main/VIDTutorial2017
         https://github.com/monoTopHelicityWG/B2GTTbar/blob/master/test/run_B2GTTbarTreeMaker_RSG_Toolbox.py#L52-L64

       using: https://github.com/cmsb2g/B2GAnaFW/blob/66ef0aa0dc5de724ececfb43f27346430c01350b/src/ElectronUserData.cc#L216
       loose,medium,tight + barrel/endcap
    """
    pfIso   = el.pfIsolationVariables();
    dEtaIn  = el.deltaEtaSuperClusterTrackAtVtx();
    dPhiIn  = el.deltaPhiSuperClusterTrackAtVtx();
    full5x5 = el.full5x5_sigmaIetaIeta();
    hoe     = el.hadronicOverEm();
    absiso  = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    relIsoWithDBeta_ = absiso/el.pt();

    if rho is not None:
        effArea   = effectiveArea(el.eta());
        absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * effArea );
    else:
        absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    relIso = absiso/el.pt();

    ooEmooP = 0.0
    if el.ecalEnergy() < 1e-8:
      ooEmooP = 999
    elif not np.isfinite(el.ecalEnergy()):
      ooEmooP = 998
    else:
      ooEmooP = np.fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() )

    missHits     = el.gsfTrack().hitPattern().numberOfLostHits(1)
    if beamspot is not None and conversions is not None: 
        hasMatchConv = ConversionTools.hasMatchedConversion(el, conversions, beamspot.position())
    else:
        hasMatchConv = False

    if not hasMatchConv:
        electronIDs = passIDWP(el.isEB(), dEtaIn, dPhiIn, full5x5, hoe, relIso, ooEmooP, missHits);
    else:
        electronIDs = {'isVETO':False,'isLoose':False,'isMedium':False,'isTight':False}

    return electronIDs


print "Creating output file "+options.outname
f = ROOT.TFile(options.outname, "RECREATE")
f.cd()


# Setup cutflow
cutflow = ROOT.TH1D('cutflow', 'cutflow',4,0,4)
cutflow.GetXaxis().SetBinLabel(1,"INITIAL");
cutflow.GetXaxis().SetBinLabel(2,"AK8_VALID");
cutflow.GetXaxis().SetBinLabel(3,"NUM_AK8");
cutflow.GetXaxis().SetBinLabel(4,"AK8_LEAD_PT");


# Setup TTree
eventTree   = ROOT.TTree('eventVars', 'eventVars')

vBESTprob_t = ROOT.vector('float')()
vBESTprob_W = ROOT.vector('float')()
vBESTprob_Z = ROOT.vector('float')()
vBESTprob_H = ROOT.vector('float')()
vBESTprob_j = ROOT.vector('float')()
vBESTclass  = ROOT.vector('float')()
eventTree.Branch( 'BESTProb_t', vBESTprob_t)
eventTree.Branch( 'BESTProb_W', vBESTprob_W)
eventTree.Branch( 'BESTProb_Z', vBESTprob_Z)
eventTree.Branch( 'BESTProb_H', vBESTprob_H)
eventTree.Branch( 'BESTProb_j', vBESTprob_j)
eventTree.Branch( 'BESTclass', vBESTclass)
vAK8pt   = ROOT.vector('float')()
vAK8eta  = ROOT.vector('float')()
vAK8phi  = ROOT.vector('float')()
vAK8mass = ROOT.vector('float')()
vAK8SDmass = ROOT.vector('float')()
vAK8tau1 = ROOT.vector('float')()
vAK8tau2 = ROOT.vector('float')()
vAK8tau3 = ROOT.vector('float')()
vAK8charge = ROOT.vector('float')()
vAK8bDiscSubjet1  = ROOT.vector('float')()
vAK8bDiscSubjet2  = ROOT.vector('float')()
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

vAK4pt    = ROOT.vector('float')()
vAK4eta   = ROOT.vector('float')()
vAK4phi   = ROOT.vector('float')()
vAK4mass  = ROOT.vector('float')()
vAK4bDisc = ROOT.vector('float')()
eventTree.Branch( 'AK4pt', vAK4pt)
eventTree.Branch( 'AK4eta', vAK4eta)
eventTree.Branch( 'AK4phi', vAK4phi)
eventTree.Branch( 'AK4mass', vAK4mass)
eventTree.Branch( 'AK4bDisc', vAK4bDisc)

vELpt  = ROOT.vector('float')()
vELeta = ROOT.vector('float')()
vELphi = ROOT.vector('float')()
vELenergy = ROOT.vector('float')()
vELcharge = ROOT.vector('float')()
vELiso = ROOT.vector('float')()
vELidLoose  = ROOT.vector('float')()
vELidMedium = ROOT.vector('float')()
vELidTight  = ROOT.vector('float')()
eventTree.Branch( 'ELpt', vELpt)
eventTree.Branch( 'ELeta', vELeta)
eventTree.Branch( 'ELphi', vELphi)
eventTree.Branch( 'ELenergy', vELenergy)
eventTree.Branch( 'ELcharge', vELcharge)
eventTree.Branch( 'ELiso', vELiso)
eventTree.Branch( 'ELid', vELidLoose)
eventTree.Branch( 'ELid', vELidMedium)
eventTree.Branch( 'ELid', vELidTight)
vMUpt  = ROOT.vector('float')()
vMUeta = ROOT.vector('float')()
vMUphi = ROOT.vector('float')()
vMUenergy  = ROOT.vector('float')()
vMUcharge  = ROOT.vector('float')()
vMUlooseID  = ROOT.vector('float')()
vMUmediumID = ROOT.vector('float')()
vMUcorrIso  = ROOT.vector('float')()
eventTree.Branch( 'MUpt', vMUpt)
eventTree.Branch( 'MUeta', vMUeta)
eventTree.Branch( 'MUphi', vMUphi)
eventTree.Branch( 'MUenergy', vMUenergy)
eventTree.Branch( 'MUcharge', vMUcharge)
eventTree.Branch( 'MUlooseID', vMUlooseID)
eventTree.Branch( 'MUmediumID',vMUmediumID)
eventTree.Branch( 'MUcorrIso', vMUcorrIso)

vGENpt  = ROOT.vector('float')()
vGENeta = ROOT.vector('float')()
vGENphi = ROOT.vector('float')()
vGENenergy   = ROOT.vector('float')()
vGENid       = ROOT.vector('int')()
vGENstatus   = ROOT.vector('int')()
vGENisHadTop = ROOT.vector('int')()
eventTree.Branch( 'GENpt', vGENpt )
eventTree.Branch( 'GENeta', vGENeta )
eventTree.Branch( 'GENphi', vGENphi )
eventTree.Branch( 'GENenergy', vGENenergy )
eventTree.Branch( 'GENid', vGENid )
eventTree.Branch( 'GENstatus', vGENstatus )
eventTree.Branch( 'GENisHadTop', vGENisHadTop )

METpt  = array('f', [-999.])
METphi = array('f', [-999.])
eventTree.Branch( 'METpt', METpt, 'METpt/F')
eventTree.Branch( 'METphi', METphi, 'METphi/F')


eventNum = array('L', [0])
runNum   = array('L', [0])
lumiNum  = array('L', [0])
eventTree.Branch( 'eventNum', eventNum, 'eventNum/L')
eventTree.Branch( 'lumiNum', lumiNum, 'lumiNum/L')
eventTree.Branch( 'runNum', runNum, 'runNum/L')


# Setup access to data in EDM file
jetsHandle = Handle("std::vector<pat::Jet>")
jetsLabel  = ("BESTProducer", "savedJets", "ttbarACskim")
AK4jetsHandle = Handle("std::vector<pat::Jet>")
AK4jetsLabel  = ("selectedAK4Jets", "", "ttbarACskim")
muonsHandle = Handle("std::vector<pat::Muon>")
muonsLabel  = ("selectedMuons", "", "ttbarACskim")
electronsHandle = Handle("std::vector<pat::Electron>")
electronsLabel  = ("selectedElectrons", "", "ttbarACskim")
metHandle  = Handle("std::vector<pat::MET>")
metLabel   = ("selectedMET", "", "ttbarACskim")
genPHandle = Handle("std::vector<reco::GenParticle>")
genPLabel  = ("selectedGenParticles", "", "ttbarACskim")

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
nevents  = 0
#s = options.files
s = 'root://cmsxrootd.fnal.gov/' + options.files

files = [s]
print 'Added ' + s

# BEST
mlp    = joblib.load('BEST_mlp.pkl')
scaler = joblib.load('BEST_scaler.pkl')

nevents_file = 0
for ifile in files:
    print 'Processing', ifile

    if options.maxevents > 0 and nevents > options.maxevents :
        break

    events   = Events( ifile )
    products = {}

    # Loop over events
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
        vELidLoose.clear()
        vELidMedium.clear()
        vELidTight.clear()
        vELcharge.clear()
        vMUpt.clear()
        vMUeta.clear()
        vMUphi.clear()
        vMUenergy.clear()
        vMUcharge.clear()
        vMUlooseID.clear()
        vMUmediumID.clear()
        vGENpt.clear()
        vGENeta.clear()
        vGENphi.clear()
        vGENenergy.clear()
        vGENid.clear()
        vGENstatus.clear()
        vGENisHadTop.clear()
    
        if options.maxevents > 0 and nevents > options.maxevents :
            break
        nevents += 1
        nevents_file += 1
        if (nevents_file % 1000 == 0):
            print nevents, ' Completed'

        cutflow.Fill(0.5)   # "INITIAL"

        try:
            eventNum[0] = event.eventAuxiliary().event()
        except OverflowError:
            print " OVERFLOW ERROR AT ",nevents
            print event.eventAuxiliary().event()
            break
        lumiNum[0] = event.eventAuxiliary().luminosityBlock()
        runNum[0]  = event.eventAuxiliary().run()

        event.getByLabel(jetsLabel, jetsHandle)
        if not jetsHandle.isValid():
            continue
        cutflow.Fill(1.5)   # "AK8_VALID"

        jets = jetsHandle.product()

        event.getByLabel(AK4jetsLabel, AK4jetsHandle)
        AK4jets = AK4jetsHandle.product()

        if options.isMC:
            event.getByLabel(genPLabel, genPHandle)
            genParticles = genPHandle.product()

            for particle in genParticles:
                if abs(particle.pdgId()) == 6 and particle.numberOfDaughters() == 2:
                    vGENpt.push_back(particle.pt())
                    vGENeta.push_back(particle.eta())
                    vGENphi.push_back(particle.phi())
                    vGENenergy.push_back(particle.energy())

                    vGENid.push_back(particle.pdgId())
                    vGENstatus.push_back(particle.status())                

                    daughter1 = particle.daughter(0)
                    daughter2 = particle.daughter(1)
                    if abs(daughter1.pdgId()) == 24:
                        vGENisHadTop.push_back( checkTopDecay(daughter1) )

                    if abs(daughter2.pdgId()) == 24:
                        vGENisHadTop.push_back( checkTopDecay(daughter2) )
                # end check top decay products
            # end loop over generator particles
        # end if isMC

        event.getByLabel(muonsLabel, muonsHandle)
        event.getByLabel(electronsLabel, electronsHandle)
        electrons = electronsHandle.product()
        muons     = muonsHandle.product()

        event.getByLabel(metLabel, metHandle)
        mets      = metHandle.product()
        METpt[0]  = mets[0].pt()
        METphi[0] = mets[0].phi()

        for i,nnHandle in enumerate(nnHandles):
            event.getByLabel(nnLabels[i], nnHandle)
            products[nnLabels[i][1]] = nnHandle.product() 

        # Minimum requirements on AK8 jets
        if len(jets) < 1:
            continue
        cutflow.Fill(2.5)   # "NUM_AK8"

        if jets[0].pt() < 350.0 or products['SDmass'][0]<20. or not ak8Quality(jets[0]):
            continue
        cutflow.Fill(3.5)   # "AK8_LEAD_PT"

        ht_ak8 = 0.0
        for j,jet in enumerate(jets):

            # Minimum pT
            if jet.pt() < 300.0 or products['SDmass'][j]<20.:
                continue

            # Quality cut
            goodJet = ak8Quality(jet)
            if not goodJet :
                continue


            # Predict BEST
            pzOverp_top = products['sumPztop'][j] / (products['sumPtop'][j] + 0.01) 
            pzOverp_W = products['sumPzW'][j] / (products['sumPW'][j] + 0.01) 
            pzOverp_Z = products['sumPzZ'][j] / (products['sumPZ'][j] + 0.01)
            pzOverp_H = products['sumPzH'][j] / (products['sumPH'][j] + 0.01)

            nnArray = np.array( [   products['SDmass'][j], products['tau32'][j], products['tau21'][j], 
                                    products['FWmoment1top'][j], products['FWmoment2top'][j], products['FWmoment3top'][j], products['FWmoment4top'][j], 
                                    products['isotropytop'][j], products['aplanaritytop'][j], products['sphericitytop'][j], products['thrusttop'][j], 
                                    products['FWmoment1W'][j], products['FWmoment2W'][j], products['FWmoment3W'][j], products['FWmoment4W'][j], 
                                    products['isotropyW'][j], products['aplanarityW'][j], products['sphericityW'][j], products['thrustW'][j], 
                                    products['FWmoment1Z'][j], products['FWmoment2Z'][j], products['FWmoment3Z'][j], products['FWmoment4Z'][j], 
                                    products['isotropyZ'][j], products['aplanarityZ'][j], products['sphericityZ'][j], products['thrustZ'][j], 
                                    products['FWmoment1H'][j], products['FWmoment2H'][j], products['FWmoment3H'][j], products['FWmoment4H'][j], 
                                    products['isotropyH'][j], products['aplanarityH'][j], products['sphericityH'][j], products['thrustH'][j], 
                                    products['bDisc'][j], products['bDisc1'][j], products['bDisc2'][j], products['q'][j], products['m12W'][j], 
                                    products['m13W'][j], products['m23W'][j], products['m1234W'][j], products['m12Z'][j], products['m13Z'][j], 
                                    products['m23Z'][j], products['m1234Z'][j], products['m12top'][j], products['m13top'][j], products['m23top'][j], 
                                    products['m1234top'][j], products['m12H'][j], products['m13H'][j], products['m23H'][j], products['m1234H'][j], 
                                    pzOverp_top, pzOverp_W, pzOverp_Z, pzOverp_H  ])

            if np.isnan(np.min(nnArray)):
                # this cut isn't included above, so it's possible to have events
                # that don't have any good AK8 jets...not sure np.isNaN ever occurs
                continue

            nnArray_transformed = scaler.transform( nnArray )
            best_bin   = mlp.predict( nnArray_transformed )
            best_probs = mlp.predict_proba( nnArray_transformed )

            # Set values
            vAK8pt.push_back(jet.pt())
            vAK8eta.push_back(jet.eta())
            vAK8phi.push_back(jet.phi())
            vAK8mass.push_back(jet.mass())
            vAK8SDmass.push_back( products['SDmass'][j] )

            vAK8tau1.push_back( jet.userFloat('NjettinessAK8:tau1') )
            vAK8tau2.push_back( jet.userFloat('NjettinessAK8:tau2') )
            vAK8tau3.push_back( jet.userFloat('NjettinessAK8:tau3') )
            vAK8charge.push_back(  products['q'][j] )
        
            vAK8bDiscSubjet1.push_back( products['bDisc1'][j] )
            vAK8bDiscSubjet2.push_back( products['bDisc2'][j] )
            vAK8ChargeSubjet1.push_back( products['qsubjet0'][j] )
            vAK8ChargeSubjet2.push_back( products['qsubjet1'][j] )

            vBESTclass.push_back(best_bin[0])    
            vBESTprob_t.push_back(best_probs[0][0])
            vBESTprob_W.push_back(best_probs[0][1])
            vBESTprob_Z.push_back(best_probs[0][2])
            vBESTprob_H.push_back(best_probs[0][3])
            vBESTprob_j.push_back(best_probs[0][4])

            ht_ak8 += jet.pt()
            
        HTak8[0] = ht_ak8


        # AK4 Jets
        for j,jet in enumerate(AK4jets):
            vAK4pt.push_back(jet.pt())
            vAK4eta.push_back(jet.eta())
            vAK4phi.push_back(jet.phi())
            vAK4mass.push_back(jet.mass())
            vAK4bDisc.push_back(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"))

        # Electrons
        # cut-based ID (want tight)
        for e,el in enumerate(electrons):
            vELpt.push_back(el.pt())
            vELeta.push_back(el.eta())
            vELphi.push_back(el.phi())
            vELenergy.push_back(el.energy())
            vELcharge.push_back(el.charge())
            # ID
            ids = electronID(el)
            vELidLoose.push_back(ids['isLoose'])
            vELidMedium.push_back(ids['isMedium'])
            vELidTight.push_back(ids['isTight'])
            # ISO
            vELiso.push_back( (el.trackIso() + el.caloIso()) / el.pt() )

        # Muons
        # cut-based ID (want medium)
        for m,mu in enumerate(muons):
            vMUpt.push_back(mu.pt())
            vMUeta.push_back(mu.eta())
            vMUphi.push_back(mu.phi())
            vMUenergy.push_back(mu.energy())
            vMUcharge.push_back(mu.charge())
            # ID
            vMUlooseID.push_back(mu.isLooseMuon())
            vMUmediumID.push_back(mu.isMediumMuon())
            #vMUtightID.push_back(mu.isTightMuon())   # requires vertex
            # ISO
            chPt = mu.pfIsolationR04().sumChargedHadronPt
            nhEt = mu.pfIsolationR04().sumNeutralHadronEt
            phEt = mu.pfIsolationR04().sumPhotonEt
            puPt = mu.pfIsolationR04().sumPUPt

            corrCombRelIso = (chPt + max(0.0, nhEt + phEt - 0.5*puPt) ) / mu.pt()
            vMUcorrIso.push_back(corrCombRelIso)

        # Fill Tree
        eventTree.Fill()

    # end loop over events
# end loop over files

# Write and close file
f.cd()
f.Write()
f.Close()


## THE END ##
