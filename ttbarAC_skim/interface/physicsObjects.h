/* 
   Physics objects to be used in analyses
   This structure allows the Event class
   and other classes to access these objects
   without circular inclusion (which breaks!)
*/
#ifndef PHYSICSOBJECTS_H_
#define PHYSICSOBJECTS_H_

#include "TLorentzVector.h"
#include <map>
#include <string>

// base object (consistent reference to TLorentzVector)
struct CmaBase {
    TLorentzVector p4;
    int isGood;
};


// Truth information
struct Parton : CmaBase {
    int pdgId;
    int index;       // index in vector of truth partons
    int decayIdx;    // index in truth record
    int parent_ref;  // index in truth vector of parent
    int parent_idx;  // index in truth record of parent
    int containment;
    int top_index;   // index of parton in the truth_top vector

    // Heavy Object Booleans
    bool isTop;
    bool isW;
    // Lepton Booleans
    bool isLepton;
    bool isTau;
    bool isElectron;
    bool isMuon;
    bool isNeutrino;
    // Quark Booleans
    bool isQuark;
    bool isBottom;
    bool isLight;
};


struct TruthTop {
    // collect indices in truth_partons vector of top parton info
    bool isTop;
    bool isAntiTop;
    int Top;                    // access the truth top in vector of partons
    int W;                      // access the truth W in vector of partons
    int bottom;                 // access the truth b in vector of partons
    std::vector<int> Wdecays;   // for storing W daughters
    std::vector<int> daughters; // for storing non-W/bottom daughters

    bool isHadronic;  // W decays to quarks
    bool isLeptonic;  // W decays to leptons
};



// Struct for jets
// -- common to all types of jets
struct Jet : CmaBase{
    float bdisc;
    std::map<std::string, bool> isbtagged;
    float charge;

    int index;       // index in vector of jets
    float radius;    // radius of jet (for truth-matching in DeltaR)

    bool isHadTop;   // matched to hadronically-decaying truth top parton
    int truth_jet;   // index in vector of truth jets that is matched to this jet
    int containment; // level of containment for partons matched to jet
    int matchId;     // keep track if jet is matched to top or anti-top
    std::vector<int> truth_partons;  // vector containing partons that are truth-matched to jet
};

struct Ljet : Jet{
    // extra ljet attributes
    float tau1;
    float tau2;
    float tau3;
    float tau21;
    float tau32;
    float softDropMass;

    float BEST_t;
    float BEST_w;
    float BEST_z;
    float BEST_h;
    float BEST_j;
    float BEST_class;

    float charge;
    std::vector<Jet> subjets;
    float subjet0_charge;
    float subjet0_bdisc;
    float subjet1_charge;
    float subjet1_bdisc;

    int target;
    std::map<std::string, double> features;  // store features in map to easily access later
    std::map<std::string, double> dnn;       // store full dnn results
};



// Extra lepton attributes
struct Lepton : CmaBase{
    // common to electrons and muons
    int charge;
    bool isElectron;
    bool isMuon;
    int index;       // index in vector of leptons

    float iso;
    float id;
    float loose;
    float medium;
    float tight;
};

struct Electron : Lepton{
    // extra electron attributes
    Electron() {
        isElectron = true;
        isMuon     = false;
    }
};
struct Muon : Lepton{
    // extra muon attributes
    Muon() {
        isElectron = false;
        isMuon     = true;
    }
};

struct Neutrino : CmaBase{
    // extra neutrino attributes
};

struct MET : CmaBase{
    // extra MET attributs
};



// Top Quark Definitions
struct Ttbar0L : CmaBase{
    // 2 FatJets
    std::vector<Ljet> ljets;
};
struct Ttbar1L : CmaBase{
    // 1 FatJet
    Ljet ljet;
    // 1 LepTop
    Lepton lepton;
    Neutrino neutrino;
    Jet jet;
};
struct Ttbar2L : CmaBase{
    // 2 LepTops -- split by top/antitop (lepton charge)
    Lepton lepton_t;
    Neutrino neutrino_t;
    Jet jet_t;

    Lepton lepton_tbar;
    Neutrino neutrino_tbar;
    Jet jet_tbar;
};

// ************************************ //
// dilepton ttbar System -- needs to be cleaned-up (21 March 2018)
struct Top : CmaBase{
    // Define a top quark
    bool isTop;
    bool isAntiTop;
};

struct LepTop : Top{
    // Leptonically-decaying top quark aspects
    Jet jet;
    Lepton lepton;
    Neutrino neutrino;

    float weight;        // reconstruction weight
    float weight_ES;     // reconstruction weight (EventShape)
    float weight_tt;     // reconstruction weight (ttbar XS method)
};

struct DileptonReco {
    // struct of information needed for neutrino reconstruction (AMWT)
    Lepton lepton_pos;      // positively charged lepton
    Lepton lepton_neg;      // negatively charged lepton
    TVector2 met;           // MET (stored as TVector2 for convenience)
    std::vector<Jet> jets;  // jets in event
    std::vector<Jet> bjets; // two 'b'-jets in ttbar decay
    Jet bJet;
    Jet bbarJet;
    size_t bJet_index, bbarJet_index;
};

struct TtbarDilepton : DileptonReco {
    /* dilepton ttbar system for dileptonTtbarReco setup
       formerly: Struct_KinematicReconstruction
       incorporated some attributes of the "TopSolution" struct
    */
    Neutrino neutrino;
    Neutrino neutrinoBar;

    TLorentzVector Wplus;
    TLorentzVector Wminus;

    LepTop top;
    LepTop topBar;
    TLorentzVector ttbar;

    double recMtop;     // top mass used in the reconstruction (needed if doing massLoop)
    double weight;      // weight of the solution
    int ntags;          // number of b-tags

    double dR;
    double dN;
    double x1;
    double x2;
    double mtt;         // invariant mass of the ttbar system
    bool isNoSmearSol;  // if smearing, this bool identifies solutions found w/o smearing

    // Enumeration for all defined weights of a kinematic reconstruction solution
    enum WeightType{defaultForMethod, neutrinoEnergy, averagedSumSmearings_mlb, undefinedWeight};
    std::map<WeightType,double> mapOfWeights;

    void set_ttbar(){
        ttbar = top.p4 + topBar.p4;
    }
    void set_W(){
        Wplus  = lepton_pos.p4 + neutrino.p4;
        Wminus = lepton_neg.p4 + neutrinoBar.p4;
    }
};




// ------------------------ // 
struct Sample {
    // Struct to contain sample information (processing the input file)
    std::string primaryDataset;
    float XSection;
    float KFactor;
    float sumOfWeights;
    unsigned int NEvents;
};


#endif
