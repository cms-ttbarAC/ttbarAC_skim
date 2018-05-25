#ifndef EVENTSAVERFLATNTUPLE_H
#define EVENTSAVERFLATNTUPLE_H

#include <string>
#include <map>
#include <vector>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ActiveAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/contrib/Njettiness.hh>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH2D.h"

#include "lwtnn/lwtnn/interface/LightweightNeuralNetwork.hh"
#include "lwtnn/lwtnn/interface/parse_json.hh"

#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/physicsObjects.h"


struct LeptonKey {
    TLorentzVector p4;
    std::vector<int> keys;
};

struct CleanJet : Jet {
    TLorentzVector uncorrP4;   // jet without JECs
    TLorentzVector originalP4; // nominal jet with JECs before cleaning (saved to update MET)
    bool isLoose;              // ID of jet (same as nominal if no cleaning)

    // PF information (for jet ID)
    float neutralHadronEnergy;
    float neutralEmEnergy;
    float chargedHadronEnergy;
    float chargedEmEnergy;
    float chargedMultiplicity;
    int numberOfDaughters;
};



class EventSaverFlatNtuple : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:

    EventSaverFlatNtuple( const edm::ParameterSet & );
    virtual ~EventSaverFlatNtuple();

    void initialize_branches();
    bool checkTopDecay(const reco::Candidate& daughter) const;
    bool passAK8( const pat::Jet& j, const float& SDmass) const;
    bool passAK4( const CleanJet& j) const;

    bool jetID( const pat::Jet& j ) const;
    bool jetID( const CleanJet& j ) const;
    CleanJet leptonJetCleaning( const pat::Jet& j );
    void updateMET( const CleanJet& j, const std::vector<TLorentzVector>& leptons4cleaning );
    void applyJEC( CleanJet& j, const float& area ) const;

    std::vector<float> charge( const reco::Jet& jet, const std::vector<float> kappas, const unsigned int first_idx=0 ) const;
    std::vector<float> getTau( unsigned int N, const reco::Jet& ij ) const;
    int findPartonIndex( const std::vector<reco::GenParticle>& items, const reco::Candidate& item ) const;

  private:

    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    TTree* m_ttree;                // physics information
    TTree* m_metadata_ttree;       // metadata
    TH1D* m_hist_cutflow;          // cutflow (for bookkeeping)

    // histograms of truth-level information
    TH1D* m_hist_truth_dy;
    TH2D* m_hist_truth_mtt_dy;
    TH2D* m_hist_truth_pttt_dy;
    TH2D* m_hist_truth_beta_dy;
    TH2D* m_hist_truth_ytt_dy;


    std::auto_ptr<fastjet::contrib::Njettiness> m_nsub;

    lwt::LightweightNeuralNetwork* m_lwtnn;
    std::map<std::string,std::vector<double>> m_BEST_products;

    std::string m_name;
    bool m_isMC;
    int m_year;
    float m_analysisAK4Pt;

    std::string t_sampleName;
    std::string t_metadataFile;

    std::map<std::string,Sample> m_mapOfSamples;     // map of Sample structs
    std::map<std::string, float> m_XSections;        // map sample name to XSection
    std::map<std::string, float> m_KFactors;         // map sample name to KFactor
    std::map<std::string, float> m_sumOfMCWeights;   // map sample name to sum of weights

    std::map<std::string,boost::shared_ptr<FactorizedJetCorrector>> m_ak4_jec;     // JECs
    JME::JetResolutionScaleFactor m_resolution_ak4sf;
    JME::JetResolutionScaleFactor m_resolution_ak8sf;

    std::vector<unsigned int> m_goodIDs = {6,15,23,24,25,7000001,8000001,9900113,9900213};   // top, bosons, plus some BSM (VLQ, Zprime, & Wprime I think)

    // Tokens
    edm::EDGetTokenT<std::vector<pat::Muon>> t_muons;
    edm::EDGetTokenT<edm::View<pat::Electron>> t_electrons;
    edm::EDGetTokenT<std::vector<pat::Jet>> t_jets;
    edm::EDGetTokenT<std::vector<pat::Jet>> t_ljets;
    edm::EDGetTokenT<reco::GenJetCollection> t_truth_ljets;
    edm::EDGetTokenT<pat::METCollection> t_met;

    edm::EDGetTokenT<int> t_runno;
    edm::EDGetTokenT<int> t_npv;
    edm::EDGetTokenT<int> t_npuTrue;
    edm::EDGetTokenT<int> t_evtno;
    edm::EDGetTokenT<int> t_lumisec;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> t_genEvtInfoProd;
    edm::EDGetTokenT<double> t_rho;
    edm::EDGetTokenT<std::vector<reco::Vertex>> t_vertices;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> t_pileup;
    edm::EDGetTokenT<reco::BeamSpot> t_beamspot;
    edm::EDGetTokenT<reco::ConversionCollection> t_conversions;
    edm::EDGetTokenT<edm::TriggerResults> t_triggerBits;
    edm::EDGetTokenT<edm::TriggerResults> t_METFilter;

    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> t_elIdFullInfoMap_Loose;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> t_elIdFullInfoMap_Medium;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> t_elIdFullInfoMap_Tight;
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> t_elIdFullInfoMap_HEEP;

    std::vector< edm::EDGetTokenT<std::vector<float>> > t_BEST_products;

    // Handles
    edm::Handle<pat::MuonCollection> m_muons;
//    edm::Handle<std::vector<pat::Electron>> m_electrons;
    edm::Handle<edm::View<pat::Electron>> m_electrons;
    edm::Handle<pat::JetCollection> m_jets;    // AK4
    edm::Handle<pat::JetCollection> m_ljets;   // AK8
    edm::Handle<pat::METCollection> m_met;
    edm::Handle<reco::GenJetCollection> m_truth_ljets; // AK8

    edm::Handle<std::vector<reco::Vertex> > h_vertices;
    edm::Handle<std::vector<PileupSummaryInfo> > h_pileup;
    edm::Handle<double> h_rho;
    edm::Handle<int> h_runno;
    edm::Handle<int> h_npv;
    edm::Handle<int> h_npuTrue;
    edm::Handle<int> h_evtno;
    edm::Handle<int> h_lumisec;
    edm::Handle<std::vector<reco::GenParticle>> h_genEvtInfoProd;
    edm::Handle<edm::TriggerResults> h_triggerBits;
    edm::Handle<edm::TriggerResults> h_METFilter;

    edm::Handle<edm::ValueMap<vid::CutFlowResult> > h_cutflow_elId_Loose;
    edm::Handle<edm::ValueMap<vid::CutFlowResult> > h_cutflow_elId_Medium;
    edm::Handle<edm::ValueMap<vid::CutFlowResult> > h_cutflow_elId_Tight;
    edm::Handle<edm::ValueMap<vid::CutFlowResult> > h_cutflow_elId_HEEP;

    // Branches
    std::string m_sampleName;
    float m_xsection;
    float m_kfactor;
    float m_sumOfWeights;
    unsigned int m_NEvents;

    std::vector<float> m_jet_pt;
    std::vector<float> m_jet_eta;
    std::vector<float> m_jet_phi;
    std::vector<float> m_jet_mass;
    std::vector<float> m_jet_area;
    std::vector<float> m_jet_deepCSV;
    std::vector<float> m_jet_bdisc;
    std::vector<float> m_jet_charge;
    std::vector<int> m_jet_ID_loose;
    std::vector<int> m_jet_ID_medium;
    std::vector<int> m_jet_ID_tight;
    std::vector<int> m_jet_ID_tightlepveto;
    std::vector<int> m_jet_true_flavor;
    std::vector<float> m_jet_uncorrPt;
    std::vector<float> m_jet_uncorrE;
    std::vector<float> m_jet_jerSF;
    std::vector<float> m_jet_jerSF_UP;
    std::vector<float> m_jet_jerSF_DOWN;

    std::vector<float> m_ljet_pt;
    std::vector<float> m_ljet_eta;
    std::vector<float> m_ljet_phi;
    std::vector<float> m_ljet_mass;
    std::vector<float> m_ljet_area;
    std::vector<float> m_ljet_tau1;
    std::vector<float> m_ljet_tau2;
    std::vector<float> m_ljet_tau3;
    std::vector<float> m_ljet_BEST_t;
    std::vector<float> m_ljet_BEST_w;
    std::vector<float> m_ljet_BEST_z;
    std::vector<float> m_ljet_BEST_h;
    std::vector<float> m_ljet_BEST_j;
    std::vector<int> m_ljet_BEST_class;
    std::vector<float> m_ljet_charge;
    std::vector<float> m_ljet_chargeSD;
    std::vector<float> m_ljet_charge3;
    std::vector<float> m_ljet_charge10;
    std::vector<float> m_ljet_SDmass;
    std::vector<int> m_ljet_ID_loose;
    std::vector<int> m_ljet_ID_medium;
    std::vector<int> m_ljet_ID_tight;
    std::vector<int> m_ljet_ID_tightlepveto;
    std::vector<float> m_ljet_subjet0_pt;
    std::vector<float> m_ljet_subjet0_mass;
    std::vector<float> m_ljet_subjet0_tau1;
    std::vector<float> m_ljet_subjet0_tau2;
    std::vector<float> m_ljet_subjet0_tau3;
    std::vector<float> m_ljet_subjet0_bdisc;
    std::vector<float> m_ljet_subjet0_deepCSV;
    std::vector<float> m_ljet_subjet0_charge;
    std::vector<float> m_ljet_subjet0_charge3;
    std::vector<float> m_ljet_subjet0_charge10;
    std::vector<float> m_ljet_subjet1_pt;
    std::vector<float> m_ljet_subjet1_mass;
    std::vector<float> m_ljet_subjet1_tau1;
    std::vector<float> m_ljet_subjet1_tau2;
    std::vector<float> m_ljet_subjet1_tau3;
    std::vector<float> m_ljet_subjet1_bdisc;
    std::vector<float> m_ljet_subjet1_deepCSV;
    std::vector<float> m_ljet_subjet1_charge;
    std::vector<float> m_ljet_subjet1_charge3;
    std::vector<float> m_ljet_subjet1_charge10;
    std::vector<float> m_ljet_uncorrPt;
    std::vector<float> m_ljet_uncorrE;
    std::vector<float> m_ljet_jerSF;
    std::vector<float> m_ljet_jerSF_UP;
    std::vector<float> m_ljet_jerSF_DOWN;

    std::vector<float> m_el_pt;
    std::vector<float> m_el_eta;
    std::vector<float> m_el_phi;
    std::vector<float> m_el_e;
    std::vector<int> m_el_iso;
    std::vector<float> m_el_charge;
    std::vector<unsigned int> m_el_ID_loose;
    std::vector<unsigned int> m_el_ID_medium;
    std::vector<unsigned int> m_el_ID_tight;
    std::vector<unsigned int> m_el_ID_HEEP;
    std::vector<unsigned int> m_el_ID_looseNoIso;
    std::vector<unsigned int> m_el_ID_mediumNoIso;
    std::vector<unsigned int> m_el_ID_tightNoIso;
    std::vector<unsigned int> m_el_ID_HEEPNoIso;
    std::vector<int> m_el_reco;
    std::vector<float> m_el_SF_ID;
    std::vector<float> m_el_SF_reco;
    std::vector<float> m_el_SF_ID_UP;
    std::vector<float> m_el_SF_reco_UP;
    std::vector<float> m_el_SF_ID_DN;
    std::vector<float> m_el_SF_reco_DN;
    std::vector<LeptonKey> m_electronKeys;

    std::vector<float> m_mu_pt;
    std::vector<float> m_mu_eta;
    std::vector<float> m_mu_phi;
    std::vector<float> m_mu_e;
    std::vector<float> m_mu_iso;
    std::vector<float> m_mu_charge;
    std::vector<unsigned int> m_mu_ID_loose;
    std::vector<unsigned int> m_mu_ID_medium;
    std::vector<unsigned int> m_mu_ID_tight;
    std::vector<float> m_mu_SF_ID;
    std::vector<float> m_mu_SF_ISO;
    std::vector<float> m_mu_SF_trigger;
    std::vector<float> m_mu_SF_track;
    std::vector<float> m_mu_SF_ID_UP;
    std::vector<float> m_mu_SF_ISO_UP;
    std::vector<float> m_mu_SF_trigger_UP;
    std::vector<float> m_mu_SF_track_UP;
    std::vector<float> m_mu_SF_ID_DN;
    std::vector<float> m_mu_SF_ISO_DN;
    std::vector<float> m_mu_SF_trigger_DN;
    std::vector<float> m_mu_SF_track_DN;
    std::vector<LeptonKey> m_muonKeys;

    float m_met_met;
    float m_met_phi;
    float m_HTAK8;
    float m_HTAK4;
    float m_HT_branch;
    float m_ST_branch;

    unsigned int m_runNumber;
    unsigned long long m_eventNumber;
    unsigned int m_lumiblock;
    float m_rho;
    int m_nGoodVtx;
    int m_LHA_PDF_ID;
    int m_nIsoTrk;
    int m_true_pileup;
    unsigned int m_npv;

    std::map<std::string,unsigned int> m_triggerBits;
    std::map<std::string,unsigned int> m_filterBits;

    float m_weight_mc;
    float m_weight_btag;
    float m_weight_pileup;
    float m_weight_jet_jer;
    float m_weight_ljet_jer;

    float m_met_met_sf;
    float m_met_phi_sf;
    std::vector<float> m_jet_jec;
    std::vector<float> m_jet_jer_up;
    std::vector<float> m_jet_jer_down;
    std::vector<float> m_ljet_jec;
    std::vector<float> m_ljet_jer_up;
    std::vector<float> m_ljet_jer_down;

    std::vector<float> m_mc_pt;
    std::vector<float> m_mc_eta;
    std::vector<float> m_mc_phi;
    std::vector<float> m_mc_e;
    std::vector<int> m_mc_pdgId;
    std::vector<int> m_mc_status;
    std::vector<int> m_mc_parent_idx;
    std::vector<int> m_mc_child0_idx;
    std::vector<int> m_mc_child1_idx;
    std::vector<int> m_mc_isHadTop;

    std::vector<float> m_truth_jet_pt;
    std::vector<float> m_truth_jet_eta;
    std::vector<float> m_truth_jet_phi;
    std::vector<float> m_truth_jet_e;
    std::vector<float> m_truth_jet_charge;

    std::vector<float> m_truth_ljet_pt;
    std::vector<float> m_truth_ljet_eta;
    std::vector<float> m_truth_ljet_phi;
    std::vector<float> m_truth_ljet_mass;
    std::vector<float> m_truth_ljet_charge;
    std::vector<float> m_truth_ljet_tau1;
    std::vector<float> m_truth_ljet_tau2;
    std::vector<float> m_truth_ljet_tau3;
    std::vector<float> m_truth_ljet_SDmass;

    // Loose jet ID
    float m_nhfLoose = 0.99;
    float m_nefLoose = 0.99;
    float m_chfLoose = 0.00;
    float m_cefLoose = 0.99;
    int m_nconstitLoose = 1;
    float m_nchLoose = 0;

    // Triggers
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
    // Boosted diff x-sec (l+jets); Zprime->ttbar (0-lepton)
    std::vector<std::string> m_triggers = {
        "HLT_Mu40_Eta2P1_PFJet200_PFJet50",                 // 2016
        "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50",    // 2016
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",            // 2017
        "HLT_Mu50",
        "HLT_TkMu50",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT",
	"HLT_PFHT800",
        "HLT_PFHT900",
        "HLT_PFHT700TrimMass50",
        "HLT_AK8PFJet450",
        "HLT_PFJet360TrimMass30"
    };

    // MET Filters
    
    std::vector<std::string> m_filters = {
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_goodVertices",
        "Flag_eeBadScFilter",
        "Flag_BadPFMuonFilter",
        "Flag_BadChargedCandidateFilter",
        "Flag_globalTightHalo2016Filter" // (data?)
    };


    // BEST variables
    std::vector<std::string> m_BEST_variables = {
        "FWmoment1H",
        "FWmoment1W",
        "FWmoment1Z",
        "FWmoment1top",
        "FWmoment2H",
        "FWmoment2W",
        "FWmoment2Z",
        "FWmoment2top",
        "FWmoment3H",
        "FWmoment3W",
        "FWmoment3Z",
        "FWmoment3top",
        "FWmoment4H",
        "FWmoment4W",
        "FWmoment4Z",
        "FWmoment4top",
        "SDmass",
        "aplanarityH",
        "aplanarityW",
        "aplanarityZ",
        "aplanaritytop",
        "bDisc",
        "bDisc1",
        "bDisc2",
        "et",
        "eta",
        "isotropyH",
        "isotropyW",
        "isotropyZ",
        "isotropytop",
        "q",
        "qsubjet0",
        "qsubjet1",
        "sphericityH",
        "sphericityW",
        "sphericityZ",
        "sphericitytop",
        "sumPH",
        "sumPW",
        "sumPZ",
        "sumPtop",
        "sumPzH",
        "sumPzW",
        "sumPzZ",
        "sumPztop",
        "tau21",
        "tau32",
        "thrustH",
        "thrustW",
        "thrustZ",
        "thrusttop",
        "m12H",
        "m23H",
        "m13H",
        "m1234H",
        "m12W",
        "m23W",
        "m13W",
        "m1234W",
        "m12Z",
        "m23Z",
        "m13Z",
        "m1234Z",
        "m12top",
        "m23top",
        "m13top",
        "m1234top"
    };

};

#endif
