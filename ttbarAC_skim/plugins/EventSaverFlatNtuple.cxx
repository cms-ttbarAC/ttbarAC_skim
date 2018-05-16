/*
Created:        29 March    2018
Last Updated:   10 May      2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Write data to flat ntuple
*/
#include "ttbarAC_skim/ttbarAC_skim/interface/EventSaverFlatNtuple.h"

using namespace edm;

EventSaverFlatNtuple::EventSaverFlatNtuple( const ParameterSet & cfg ) :
  m_lwtnn(nullptr),
  t_sampleName(cfg.getParameter<std::string>("sampleName")),
  t_metadataFile(cfg.getParameter<std::string>("metadataFile")),
  t_muons(consumes<pat::MuonCollection>(edm::InputTag("selectedMuons", "", "ttbarACskim"))),
  t_electrons(consumes<edm::View<pat::Electron>>(edm::InputTag("slimmedElectrons","",""))),
  t_jets(consumes<pat::JetCollection>(edm::InputTag("selectedAK4Jets", "", "ttbarACskim"))),
  t_ljets(consumes<pat::JetCollection>(edm::InputTag("BESTProducer", "savedJets", "ttbarACskim"))),
  t_met(consumes<pat::METCollection>(edm::InputTag("selectedMET", "", "ttbarACskim"))),
  t_rho(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  t_vertices(consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"))),
  t_pileup(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
  t_beamspot(consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"))),
  t_conversions(consumes<reco::ConversionCollection>(edm::InputTag("reducedEgamma:reducedConversions"))),
  t_METFilter(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO"))),
  t_elIdFullInfoMap_Loose(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_Loose"))),
  t_elIdFullInfoMap_Medium(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_Medium"))),
  t_elIdFullInfoMap_Tight(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_Tight"))){
  //t_elIdFullInfoMap_HEEP(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_HEEP"))){
    m_isMC = cfg.getParameter<bool>("isMC");           // filling truth branches

    t_BEST_products.clear();
    for (const auto& name : m_BEST_variables){
        edm::EDGetTokenT<std::vector<float>> tmp_token = consumes<std::vector<float>>(edm::InputTag("BESTProducer",name,"ttbarACskim"));
        t_BEST_products.push_back( tmp_token );
    }

    if (m_isMC){
        t_truth_ljets = consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"));
        t_genEvtInfoProd = consumes<std::vector<reco::GenParticle>>(edm::InputTag("selectedGenParticles", "", "ttbarACskim"));
    }

    bool reHLT = (t_sampleName.find("reHLT")!=std::string::npos) || (t_sampleName.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8")!=std::string::npos);
    std::string hlt = reHLT ? "HLT2" : "HLT";
    t_triggerBits = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", hlt));


    // Make output TTrees
    edm::Service<TFileService> fs;
    m_ttree = fs->make<TTree>("eventVars","eventVars");
    m_metadata_ttree = fs->make<TTree>("metadata","metadata");

    m_hist_cutflow = fs->make<TH1D>( "cutflow","cutflow",3,0,3);
    m_hist_cutflow->GetXaxis()->SetBinLabel(1,"INITIAL");
    m_hist_cutflow->GetXaxis()->SetBinLabel(2,"PRIMARYVTX");
    m_hist_cutflow->GetXaxis()->SetBinLabel(3,"AK8JETS");

    // TH1D
    m_hist_truth_dy = fs->make<TH1D>( "truth_dy","truth_dy",2000,-10,10);  // (unbounded, likely between -5,5)

    // TH2D (x,y)
    m_hist_truth_mtt_dy  = fs->make<TH2D>( "truth_mtt_dy", "truth_mtt_dy", 6000,0,6000, 2000,-10,10);   // (unbounded, likely < 5000)
    m_hist_truth_pttt_dy = fs->make<TH2D>( "truth_pttt_dy","truth_pttt_dy",1000,0,1000, 2000,-10,10);   // (unbounded, likely < 500)
    m_hist_truth_beta_dy = fs->make<TH2D>( "truth_beta_dy","truth_beta_dy",1000,0,1,    2000,-10,10);   // (0,1)
    m_hist_truth_ytt_dy  = fs->make<TH2D>( "truth_ytt_dy", "truth_ytt_dy", 2000,0,10,   2000,-10,10);   // (unbounded, likely between -5,5 -> absolute value)

    initialize_branches();


    // options set by config
    cma::setVerboseLevel("WARNING");
    m_sampleName = t_sampleName;
    if (m_isMC){
        m_mapOfSamples.clear();
        cma::getSampleWeights( t_metadataFile,m_mapOfSamples );
    }
    bool sampInFile = (m_mapOfSamples.find(m_sampleName)!=m_mapOfSamples.end());
    std::cout << " SAMPLE NAME " << m_sampleName << ": Found = " << sampInFile << std::endl;

    // N-subjettiness for softdrop subjets
    fastjet::contrib::NormalizedMeasure nsub_normalizedMeasure(1.0,0.8);
    fastjet::contrib::OnePass_KT_Axes nsub_onepass_kt_axes;
    fastjet::contrib::MeasureDefinition const* measureDef = nullptr;
    fastjet::contrib::AxesDefinition const* axesDef = nullptr;

    measureDef = &nsub_normalizedMeasure;
    axesDef    = &nsub_onepass_kt_axes;

    m_nsub = std::auto_ptr<fastjet::contrib::Njettiness> ( new fastjet::contrib::Njettiness( *axesDef, *measureDef ) );

    // Lightweight NN interface with BEST
    std::string dnnFile("");

    try {
        edm::FileInPath lwtnn_cfg("ttbarAC_skim/ttbarAC_skim/test/BEST_mlp.json");
        dnnFile = lwtnn_cfg.fullPath();
    }
    catch (...) {
        dnnFile = "BEST_mlp.json";
    }

    std::ifstream input_cfg( dnnFile );                     // original: "data/BEST_mlp.json"
    lwt::JSONConfig lwcfg = lwt::parse_json( input_cfg );
    m_lwtnn = new lwt::LightweightNeuralNetwork(lwcfg.inputs, lwcfg.layers, lwcfg.outputs);

    // New JEC for cleaned-jets
    //   Guide: https://github.com/dmajumder/VLQAna/blob/CMSSW_8_0_X_NewB2GAnaFW/src/JetMaker.cc
    //   -> https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Jet_Energy_Corrections_in_Run2
    //      Not updating to newest recommendation (as of 11 May) because there is no GT listed
    //   JECDatabase/textFiles/
    //    Summer16_23Sep2016 + BCDV4_DATA || EFV4_DATA || GV4_DATA || HV4_DATA || V4_MC
    m_ak4_jec.clear();
    std::string JECpathbase("JECDatabase/textFiles/");
    std::string JECdate("Summer16_23Sep2016");
    std::vector<std::string> eras = {"V4_MC","BCDV4_DATA","EFV4_DATA","GV4_DATA","HV4_DATA"};
    for (const auto& era : eras){
        std::vector<JetCorrectorParameters> vPar;  

        std::string JECpath = JECpathbase+JECdate+era;
        std::vector<std::string> JECpayloads = {
            JECdate+era+"_L1FastJet_AK4PFchs.txt",
            JECdate+era+"_L2Relative_AK4PFchs.txt",
            JECdate+era+"_L3Absolute_AK4PFchs.txt"    // JECpath+"/"+
        };

        // extra correction for data
        if (era.find("_MC")==std::string::npos)
            JECpayloads.push_back(JECdate+era+"_L2L3Residual_AK4PFchs.txt");

        for (const auto& payload : JECpayloads){
            std::string jec_file(JECpath+"/"+payload);
            JetCorrectorParameters pars(jec_file);
            vPar.push_back(pars) ; 
        }

        m_ak4_jec[era] = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
    } // end loop over data eras and MC

    std::string ak4jersfFile("JERDatabase/textFiles/Summer16_25nsV1_MC_SF_AK4PFchs.txt");
    std::string ak8jersfFile("JERDatabase/textFiles/Summer16_25nsV1_MC_SF_AK8PFchs.txt");

    m_resolution_ak4sf = JME::JetResolutionScaleFactor(ak4jersfFile);
    m_resolution_ak8sf = JME::JetResolutionScaleFactor(ak8jersfFile);
}


EventSaverFlatNtuple::~EventSaverFlatNtuple() {
    delete m_lwtnn;
}



void EventSaverFlatNtuple::beginJob(){
    /* Begin Job */
    return;
}


void EventSaverFlatNtuple::analyze( const edm::Event& event, const edm::EventSetup& setup) {
    /* Fill TTree 
       This is the function to modify / inherit for analysis-specific purposes
    */
    // Load data 
    event.getByToken( t_electrons, m_electrons );
    event.getByToken( t_muons, m_muons );
    event.getByToken( t_jets, m_jets );
    event.getByToken( t_met, m_met );
    event.getByToken( t_rho, h_rho );
    event.getByToken( t_vertices,h_vertices );
    event.getByToken( t_pileup,h_pileup );

    event.getByToken( t_triggerBits,h_triggerBits );
    event.getByToken( t_METFilter,  h_METFilter );

    event.getByToken( t_elIdFullInfoMap_Loose,  h_cutflow_elId_Loose );
    event.getByToken( t_elIdFullInfoMap_Medium, h_cutflow_elId_Medium );
    event.getByToken( t_elIdFullInfoMap_Tight,  h_cutflow_elId_Tight );
    //event.getByToken( t_elIdFullInfoMap_HEEP,   h_cutflow_elId_HEEP );

    // Start with filling cutflow
    m_hist_cutflow->Fill(0.5);  // INITIAL

    // Check generator-level information first
    // Fill ttbar truth distributions (for unfolding) 
    m_mc_pt.clear();
    m_mc_eta.clear();
    m_mc_phi.clear();
    m_mc_e.clear();
    m_mc_pdgId.clear();
    m_mc_status.clear();
    m_mc_parent_idx.clear();
    m_mc_child0_idx.clear();
    m_mc_child1_idx.clear();
    m_mc_isHadTop.clear();

    if (m_isMC){
        event.getByToken( t_truth_ljets, m_truth_ljets );
        event.getByToken( t_genEvtInfoProd,h_genEvtInfoProd );

        // create temporary vector of genParticles to store only a few particles
        std::vector<reco::GenParticle> genCollection_tmp;
        std::unique_ptr<std::vector<reco::GenParticle> > genCollection( new std::vector<reco::GenParticle> (*h_genEvtInfoProd) );
        for (unsigned int j=0, size=genCollection->size(); j<size; j++){
            reco::GenParticle particle = genCollection->at(j);

            unsigned int absPdgId = std::abs( particle.pdgId() );
            int parent_pdgId(0);
            if (particle.numberOfMothers()>0 && particle.mother(0)!=nullptr)
                parent_pdgId = particle.mother(0)->pdgId();

            // Check that this particle has a PDGID of interest, or that its parent does
            if ( std::find(m_goodIDs.begin(), m_goodIDs.end(), absPdgId) == m_goodIDs.end() &&
                 std::find(m_goodIDs.begin(), m_goodIDs.end(), std::abs(parent_pdgId)) == m_goodIDs.end() )
                continue;

            genCollection_tmp.push_back(particle);
        }

        TLorentzVector top;
        TLorentzVector antitop;
        // Now loop over the slimmed container of gen particles to save them and references to parent/children
        for (const auto& particle : genCollection_tmp){
            m_mc_pt.push_back(particle.pt());
            m_mc_eta.push_back(particle.eta());
            m_mc_phi.push_back(particle.phi());
            m_mc_e.push_back(particle.energy());
            m_mc_pdgId.push_back(particle.pdgId());
            m_mc_status.push_back(particle.status());

            // save the index (in the 'goodMCs' vector) of the parent/child (if they exist)
            int parent_idx(-1);
            int child0_idx(-1);
            int child1_idx(-1);

            if (particle.numberOfMothers()>0){
                auto parent = particle.mother(0);
                parent_idx  = findPartonIndex(genCollection_tmp,*parent);
            }

            if (particle.numberOfDaughters()>0){
                auto child0 = particle.daughter(0);
                child0_idx  = findPartonIndex(genCollection_tmp,*child0);
                if (particle.numberOfDaughters()>1){
                    auto child1 = particle.daughter(1);
                    child1_idx  = findPartonIndex(genCollection_tmp,*child1);
                }
            }

            m_mc_parent_idx.push_back( parent_idx );
            m_mc_child0_idx.push_back( child0_idx );
            m_mc_child1_idx.push_back( child1_idx );

            unsigned int isHadTop(0);
            if (std::abs(particle.pdgId())==6 && particle.numberOfDaughters()==2){
                auto* daughter1 = particle.daughter(0);
                auto* daughter2 = particle.daughter(1);

                if (std::abs(daughter1->pdgId()) == 24)
                    isHadTop = checkTopDecay(*daughter1);
                if (std::abs(daughter2->pdgId()) == 24)
                    isHadTop = checkTopDecay(*daughter2);

                if (particle.pdgId()>0)
                    top.SetPtEtaPhiE( particle.pt(), particle.eta(), particle.phi(), particle.energy() );
                else
                    antitop.SetPtEtaPhiE( particle.pt(), particle.eta(), particle.phi(), particle.energy() );

            }

            m_mc_isHadTop.push_back( isHadTop );
        } //  end loop over slimmed collection of gen particles

        if (top.Pt()>1 && antitop.Pt()>1){
            TLorentzVector ttbar = top+antitop;
            double dy   = std::abs(top.Rapidity()) - std::abs(antitop.Rapidity());
            double mtt  = ttbar.M();
            double pttt = ttbar.Pt();
            double ytt  = std::abs(ttbar.Rapidity());
            double beta = std::abs(top.Pz() + antitop.Pz()) / (top.E() + antitop.E());

            m_hist_truth_dy->Fill(dy);
            m_hist_truth_mtt_dy->Fill(mtt,dy);
            m_hist_truth_pttt_dy->Fill(pttt,dy);
            m_hist_truth_beta_dy->Fill(beta,dy);
            m_hist_truth_ytt_dy->Fill(ytt,dy);
         }
    } // end if isMC


    if (h_vertices->empty()) return;               // skip the event if no PV found
    m_hist_cutflow->Fill(1.5);   // Primary Vertex


    // May not have the BEST products
    // -- https://twiki.cern.ch/twiki/bin/view/Main/CMSSWCheatSheet#Does_a_certain_product_exist_in
    try {event.getByToken( t_ljets, m_ljets );}
    catch( cms::Exception& ex ) {std::cout << " > Large-R Jets not found " << std::endl;}
    if (!m_ljets.isValid()) {
        std::cout << " > Product not valid: Large-R Jets " << std::endl;
        return;
    }


    // Set branch values
    m_runNumber   = event.id().run();
    m_eventNumber = event.id().event();
    m_lumiblock   = event.id().luminosityBlock();
    m_rho         = *h_rho;

    // -- primary vertex information
    m_npv = h_vertices->size();                    // primary vertices
    const reco::Vertex &PV = h_vertices->front();  // save PV for tight muon ID

    // -- pileup
    m_true_pileup = 0;
    if (h_pileup.isValid()) { // protection for data
        for(std::vector<PileupSummaryInfo>::const_iterator iPV=h_pileup->begin(); iPV!=h_pileup->end(); ++iPV) {
            if(iPV->getBunchCrossing()==0) {
                m_true_pileup = iPV->getTrueNumInteractions();  
                break;
            }
        }
    }

    // Trigger Bits (save the decision and apply offline)
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
    const edm::TriggerNames& names = event.triggerNames(*h_triggerBits);
    for (unsigned int i=0, size=h_triggerBits->size(); i<size; ++i) {
        const std::string name(names.triggerName(i));
        for (const auto& trig : m_triggers){
            std::size_t found = name.find( trig );
            if (found==std::string::npos) continue;
            m_triggerBits.at(trig) = h_triggerBits->accept(i);
        }
    }


    // MET Filters (save the decision and apply offline)
    //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    edm::TriggerNames const& filterNames = event.triggerNames(*h_METFilter);
    unsigned int nMETfilters = h_METFilter->size();

    // Check the filters (that we care about) passed
    for (unsigned int i=0; i<nMETfilters; i++){
        const std::string name(filterNames.triggerName(i));
        if ( std::find( m_filters.begin(), m_filters.end(), name ) == m_filters.end() ) continue;

        m_filterBits.at(name) = h_METFilter->accept(i);
    }

    // Missing Transverse Momentum
    m_met_met = (*m_met.product())[0].pt();
    m_met_phi = (*m_met.product())[0].phi();


    // AK8 jets
    m_ljet_pt.clear();
    m_ljet_eta.clear();
    m_ljet_phi.clear();
    m_ljet_mass.clear();
    m_ljet_area.clear();
    m_ljet_tau1.clear();
    m_ljet_tau2.clear();
    m_ljet_tau3.clear();
    m_ljet_charge.clear();
    m_ljet_chargeSD.clear();
    m_ljet_charge3.clear();
    m_ljet_charge10.clear();
    m_ljet_SDmass.clear();
    m_ljet_BEST_t.clear();
    m_ljet_BEST_w.clear();
    m_ljet_BEST_z.clear();
    m_ljet_BEST_h.clear();
    m_ljet_BEST_j.clear();
    m_ljet_BEST_class.clear();
    m_ljet_subjet0_pt.clear();
    m_ljet_subjet0_mass.clear();
    m_ljet_subjet0_tau1.clear();
    m_ljet_subjet0_tau2.clear();
    m_ljet_subjet0_tau3.clear();
    m_ljet_subjet0_bdisc.clear();
    m_ljet_subjet0_deepCSV.clear();
    m_ljet_subjet0_charge.clear();
    m_ljet_subjet0_charge3.clear();
    m_ljet_subjet0_charge10.clear();
    m_ljet_subjet1_pt.clear();
    m_ljet_subjet1_mass.clear();
    m_ljet_subjet1_tau1.clear();
    m_ljet_subjet1_tau2.clear();
    m_ljet_subjet1_tau3.clear();
    m_ljet_subjet1_bdisc.clear();
    m_ljet_subjet1_deepCSV.clear();
    m_ljet_subjet1_charge.clear();
    m_ljet_subjet1_charge3.clear();
    m_ljet_subjet1_charge10.clear();
    m_ljet_uncorrPt.clear();
    m_ljet_uncorrE.clear();
    m_ljet_jerSF.clear();
    m_ljet_jerSF_UP.clear();
    m_ljet_jerSF_DOWN.clear();
    m_HTAK8 = 0;

    m_BEST_products.clear();
    for (unsigned int i=0,size=m_BEST_variables.size(); i<size; i++){
        edm::Handle<std::vector<float>> h_tmp;
        event.getByToken( t_BEST_products.at(i), h_tmp );

        std::vector<float> ftmp = *h_tmp.product();
        std::vector<double> dtmp(ftmp.begin(),ftmp.end());
        m_BEST_products[m_BEST_variables[i]] = dtmp;
    }

    unsigned int nj(0);
    for (const auto& ljet : *m_ljets.product()){
        // Check jet
        float SDmass = m_BEST_products.at("SDmass")[nj];
        bool pass = passAK8( ljet,SDmass );
        if (!pass){
            nj++;
            continue;
        }

        m_ljet_pt.push_back(  ljet.pt());
        m_ljet_eta.push_back( ljet.eta());
        m_ljet_phi.push_back( ljet.phi());
        m_ljet_mass.push_back(ljet.mass());
        m_ljet_area.push_back(ljet.jetArea());

        m_ljet_tau1.push_back(ljet.userFloat("NjettinessAK8:tau1"));
        m_ljet_tau2.push_back(ljet.userFloat("NjettinessAK8:tau2"));
        m_ljet_tau3.push_back(ljet.userFloat("NjettinessAK8:tau3"));

        // BESTProducer
        m_ljet_charge.push_back( m_BEST_products.at("q")[nj] );
        m_ljet_SDmass.push_back( SDmass );

        // BEST score
        std::map<std::string,double> BEST_products;
        for (auto& x : m_BEST_products)
            BEST_products[x.first] = x.second[nj];

        std::vector<std::string> jetNames = {"top","W","Z","H","jet"};
        for (unsigned int pp=0,size=jetNames.size(); pp<size; pp++){
            std::string jetName = jetNames[pp];

            float sumPz = BEST_products["sumPz"+jetName];
            float sumP  = BEST_products["sumP"+jetName];
            BEST_products["pzOverp_"+jetName] =  ( sumPz / (sumP + 0.0001) ); // not used for 'jet'
        }

        std::map<std::string,double> NNresults = m_lwtnn->compute(BEST_products);
        std::vector<double> values{ NNresults["dnn_qcd"],   NNresults["dnn_top"],
                                    NNresults["dnn_higgs"], NNresults["dnn_z"], NNresults["dnn_w"] };

        unsigned int particleID(0);
        float max_value(-1.0);
        for (unsigned int pid=0,size=values.size();pid<size;pid++){
            if (values.at(pid) > max_value){
                max_value  = values.at(pid);
                particleID = pid;
            }
        }

        m_ljet_BEST_t.push_back( NNresults.at("dnn_top") );
        m_ljet_BEST_w.push_back( NNresults.at("dnn_w") );
        m_ljet_BEST_z.push_back( NNresults.at("dnn_z") );
        m_ljet_BEST_h.push_back( NNresults.at("dnn_higgs") );
        m_ljet_BEST_j.push_back( NNresults.at("dnn_qcd") );
        m_ljet_BEST_class.push_back( particleID );

        // Subjet information (soft drop produces 2 subjets)
        auto subjets = ljet.subjets();
        auto subjet0 = subjets.at(0);
        auto subjet1 = subjets.at(1);


        // Calculate some extra jet charges
        std::vector<float> kappas = {0.3,1.0,0.6};                    // nominal = 0.6, re-calculate to compare

        std::vector<float> subjet0_q = charge( *subjet0, kappas );    // first subjet
        std::vector<float> subjet1_q = charge( *subjet1, kappas );    // second subjet
        std::vector<float> ljet_q    = charge( ljet, kappas, 2 );     // non-softdrop constituents

        // -- jet charge due to only the constituents in the soft drop subjet
        float softDropOnly_q(0.0);
        softDropOnly_q = subjet0_q.at(2) + subjet1_q.at(2);
        TLorentzVector softdrop0;
        softdrop0.SetPtEtaPhiE( subjet0->pt(), subjet0->eta(), subjet0->phi(), subjet0->energy() );
        TLorentzVector softdrop1;
        softdrop1.SetPtEtaPhiE( subjet1->pt(), subjet1->eta(), subjet1->phi(), subjet1->energy() );
        float softdrop_pt_w = pow( (softdrop0+softdrop1).Pt(), 0.6);
        softDropOnly_q /= softdrop_pt_w;
        m_ljet_chargeSD.push_back( softDropOnly_q );

        // Normalize pt-weighted sum by the weighted pt of the jet
        for (unsigned int k=0,size=kappas.size(); k<size; k++){
            ljet_q.at(k)    += subjet0_q.at(k) + subjet1_q.at(k);     // include soft-drop components

            ljet_q.at(k)    /= pow( ljet.pt(),kappas.at(k) );
            subjet0_q.at(k) /= pow( subjet0->pt(),kappas.at(k) );
            subjet1_q.at(k) /= pow( subjet1->pt(),kappas.at(k) );
        }
        m_ljet_charge3.push_back( subjet0_q.at(0) );
        m_ljet_charge10.push_back( subjet1_q.at(1) );


//        double deepCSV0_b  = subjet0->bDiscriminator("pfDeepCSVJetTags:probb");
//        double deepCSV0_bb = subjet0->bDiscriminator("pfDeepCSVJetTags:probbb");

        m_ljet_subjet0_bdisc.push_back(  m_BEST_products["bDisc1"][nj] );
        m_ljet_subjet0_charge.push_back( m_BEST_products["qsubjet0"][nj] );
        //m_ljet_subjet0_deepCSV.push_back(deepCSV0_b+deepCSV0_bb);
        m_ljet_subjet0_deepCSV.push_back(subjet0->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") );
        m_ljet_subjet0_pt.push_back(   subjet0->pt() );
        m_ljet_subjet0_mass.push_back( subjet0->mass() );
        m_ljet_subjet0_charge3.push_back(subjet0_q.at(0) );   // jet charge with kappa=0.3
        m_ljet_subjet0_charge10.push_back(subjet0_q.at(1) );  // jet charge with kappa=1.0
        // add some substructure to distinguish W/b subjets
        std::vector<float> taus = getTau(3,*subjet0);
        m_ljet_subjet0_tau1.push_back( taus.at(0) );    // subjet0->userFloat("NjettinessAK8PFCHSSoftDropSubjets:tau1")
        m_ljet_subjet0_tau2.push_back( taus.at(1) );    // subjet0->userFloat("NjettinessAK8PFCHSSoftDropSubjets:tau2")
        m_ljet_subjet0_tau3.push_back( taus.at(2) );    // subjet0->userFloat("NjettinessAK8PFCHSSoftDropSubjets:tau3")


//        double deepCSV1_b  = subjet1->bDiscriminator("pfDeepCSVJetTags:probb");
//        double deepCSV1_bb = subjet1->bDiscriminator("pfDeepCSVJetTags:probbb");

        m_ljet_subjet1_bdisc.push_back(  m_BEST_products["bDisc2"][nj] );
        m_ljet_subjet1_charge.push_back( m_BEST_products["qsubjet1"][nj] );
        //m_ljet_subjet1_deepCSV.push_back(deepCSV1_b+deepCSV1_bb);
        m_ljet_subjet1_deepCSV.push_back(subjet1->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") );
        m_ljet_subjet1_pt.push_back(   subjet1->pt() );
        m_ljet_subjet1_mass.push_back( subjet1->mass() );
        m_ljet_subjet1_charge3.push_back( subjet1_q.at(0) );  // jet charge with kappa=0.3
        m_ljet_subjet1_charge10.push_back( subjet1_q.at(1) ); // jet charge with kappa=1.0
        // add some substructure to distinguish W/b subjets
        taus = getTau(3,*subjet1);
        m_ljet_subjet1_tau1.push_back( taus.at(0) );    // subjet1->userFloat("NjettinessAK8PFCHSSoftDropSubjets:tau1")
        m_ljet_subjet1_tau2.push_back( taus.at(1) );    // subjet2->userFloat("NjettinessAK8PFCHSSoftDropSubjets:tau2")
        m_ljet_subjet1_tau3.push_back( taus.at(2) );    // subjet3->userFloat("NjettinessAK8PFCHSSoftDropSubjets:tau3")

        // save 'uncorrected' information in case we need to redo the JECs offline
        reco::Candidate::LorentzVector uncorrJet = ljet.correctedP4(0);
        m_ljet_uncorrPt.push_back(uncorrJet.pt());
        m_ljet_uncorrE.push_back(uncorrJet.energy());

        JME::JetParameters parameters = {{JME::Binning::JetEta, ljet.eta()}, {JME::Binning::Rho, m_rho}};
        m_ljet_jerSF.push_back( m_resolution_ak8sf.getScaleFactor(parameters) );
        m_ljet_jerSF_UP.push_back( m_resolution_ak8sf.getScaleFactor(parameters, Variation::UP) );
        m_ljet_jerSF_DOWN.push_back( m_resolution_ak8sf.getScaleFactor(parameters, Variation::DOWN) );

        m_HTAK8+=ljet.pt();
        nj++;
    } // end loop over AK8

    if (m_ljet_pt.size()<1)
        return;
    m_hist_cutflow->Fill(2.5);    // AK8Jets


    // Truth AK8
    m_truth_ljet_pt.clear();
    m_truth_ljet_eta.clear();
    m_truth_ljet_phi.clear();
    m_truth_ljet_mass.clear();
    m_truth_ljet_tau1.clear();
    m_truth_ljet_tau2.clear();
    m_truth_ljet_tau3.clear();
    m_truth_ljet_SDmass.clear();

    if (m_isMC){
        for (const auto& ljet : *m_truth_ljets.product()){
            m_truth_ljet_pt.push_back(  ljet.pt());
            m_truth_ljet_eta.push_back( ljet.eta());
            m_truth_ljet_phi.push_back( ljet.phi());
            m_truth_ljet_mass.push_back(ljet.mass());

//            m_truth_ljet_SDmass.push_back(ljet.userFloat("ak8PFJetsCHSSoftDropMass"));
            std::vector<float> taus = getTau(3,ljet);
            m_truth_ljet_tau1.push_back( taus.at(0) );  // "NjettinessAK8CHS:tau1"
            m_truth_ljet_tau2.push_back( taus.at(1) );  // "NjettinessAK8CHS:tau2"
            m_truth_ljet_tau3.push_back( taus.at(2) );  // "NjettinessAK8CHS:tau3"
        } // end loop over truth AK8
    } 

    // Leptons
    m_el_pt.clear();
    m_el_eta.clear();
    m_el_phi.clear();
    m_el_e.clear();
    //m_el_iso.clear();
    m_el_charge.clear();
    m_el_ID_loose.clear();
    m_el_ID_medium.clear();
    m_el_ID_tight.clear();
    //m_el_ID_HEEP.clear();
    m_el_ID_looseNoIso.clear();
    m_el_ID_mediumNoIso.clear();
    m_el_ID_tightNoIso.clear();
    m_electronKeys.clear();

    for (size_t i=0, size=m_electrons->size(); i<size; ++i){   
        const auto el = m_electrons->ptrAt(i);          // easier if we use ptrs for the id

        // ID
        vid::CutFlowResult idLoose  = (*h_cutflow_elId_Loose)[el];    // includes isolation requirment
        vid::CutFlowResult idMedium = (*h_cutflow_elId_Medium)[el];
        vid::CutFlowResult idTight  = (*h_cutflow_elId_Tight)[el];

        // ID -- no isolation
        // https://twiki.cern.ch/twiki/bin/view/Main/VIDTutorial2017#Performing_N_1_Studies
        vid::CutFlowResult idLooseNoIso  = idLoose.getCutFlowResultMasking(7);   // 7 = GsfEleEffAreaPFIsoCut_0
        vid::CutFlowResult idMediumNoIso = idMedium.getCutFlowResultMasking(7);
        vid::CutFlowResult idTightNoIso  = idTight.getCutFlowResultMasking(7);

        // selector: pt > 40.0 && abs(eta) < 2.4; require at least one of the IDs (no iso.) to fire
        if (el->pt() < 40. || std::abs(el->eta())>2.4) continue;
        if (!idLooseNoIso.cutFlowPassed() && !idMediumNoIso.cutFlowPassed() && !idTightNoIso.cutFlowPassed()) continue;

        // save information
        m_el_pt.push_back( el->pt() );
        m_el_eta.push_back(el->eta() );
        m_el_phi.push_back(el->phi() );
        m_el_e.push_back(  el->energy() );

        m_el_charge.push_back( el->charge() );
        //m_el_iso.push_back( (el->trackIso() + el->caloIso()) / el->pt() );

        m_el_ID_loose.push_back( idLoose.cutFlowPassed() );
        m_el_ID_medium.push_back(idMedium.cutFlowPassed() );
        m_el_ID_tight.push_back( idTight.cutFlowPassed() );
        //m_el_ID_HEEP.push_back(  idHEEP.cutFlowPassed() );

        m_el_ID_looseNoIso.push_back( idLooseNoIso.cutFlowPassed() );
        m_el_ID_mediumNoIso.push_back(idMediumNoIso.cutFlowPassed() );
        m_el_ID_tightNoIso.push_back( idTightNoIso.cutFlowPassed() );

        // lepton cleaning -- selector: TightID (no iso); pt > 50.0; abs(eta) < 2.4
        if (el->pt()>50. && std::abs(el->eta())<2.4 && idTightNoIso.cutFlowPassed() ){
            LeptonKey el_key;
            el_key.p4.SetPtEtaPhiE( el->pt(), el->eta(), el->phi(), el->energy() );
            el_key.keys.clear();
            for (unsigned int kk=0, size=el->numberOfSourceCandidatePtrs(); kk<size; kk++)
                el_key.keys.push_back( el->sourceCandidatePtr(kk).key() ); // el->originalObjectRef().key()

            m_electronKeys.push_back( el_key );
        }
    } // end loop over electrons

    m_mu_pt.clear();
    m_mu_eta.clear();
    m_mu_phi.clear();
    m_mu_e.clear();
    m_mu_iso.clear();
    m_mu_charge.clear();
    m_mu_ID_loose.clear();
    m_mu_ID_medium.clear();
    m_mu_ID_tight.clear();
    m_muonKeys.clear();

    for (const auto& mu : *m_muons.product()){
        // ID
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Muon_Identification
        bool isLoose( muon::isLooseMuon(mu) );
        bool isTight( muon::isTightMuon(mu,PV) );

        // Medium ID  -- 2016 data B-F vs G-H
        // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=27#MediumID2016_to_be_used_with_Run
        bool goodGlob   = mu.isGlobalMuon() && 
                          mu.globalTrack()->normalizedChi2() < 3 && 
                          mu.combinedQuality().chi2LocalPosition < 12 && 
                          mu.combinedQuality().trkKink < 20;
        bool isMedium   = muon::isMediumMuon(mu);
        bool isMediumBF = muon::isLooseMuon(mu) && 
                          mu.innerTrack()->validFraction() > 0.49 && 
                          muon::segmentCompatibility(mu) > (goodGlob ? 0.303 : 0.451); 
        bool isMediumMuon = (event.isRealData() && m_runNumber <= 278808) ? isMediumBF : isMedium;

        if (!isLoose && !isMediumMuon && !isTight) continue;  // require at least 1 ID to fire

        // save information
        m_mu_pt.push_back(  mu.pt() );
        m_mu_eta.push_back( mu.eta() );
        m_mu_phi.push_back( mu.phi() );
        m_mu_e.push_back(   mu.energy() );
        m_mu_charge.push_back( mu.charge() );

        // isolation "ISO04"
        float chPt = mu.pfIsolationR04().sumChargedHadronPt;
        float nhPt = mu.pfIsolationR04().sumNeutralHadronEt;
        float phPt = mu.pfIsolationR04().sumPhotonEt;
        float puPt = mu.pfIsolationR04().sumPUPt;
        float iso04 = (chPt+TMath::Max(0.,nhPt+phPt-0.5*puPt))/mu.pt();
        m_mu_iso.push_back( iso04 );

        m_mu_ID_loose.push_back(  isLoose );
        m_mu_ID_medium.push_back( isMediumMuon );
        m_mu_ID_tight.push_back(  isTight );

        // lepton cleaning -- selector: MediumID (no iso); pt > 50.0; abs(eta) < 2.4
        if (mu.pt()>50. && std::abs(mu.eta())<2.4 && isMediumMuon ) {
            LeptonKey mu_key;
            mu_key.p4.SetPtEtaPhiE( mu.pt(), mu.eta(), mu.phi(), mu.energy() );
            mu_key.keys.clear();
            for (unsigned int kk=0, size=mu.numberOfSourceCandidatePtrs(); kk<size; kk++)
                mu_key.keys.push_back( mu.sourceCandidatePtr(kk).key() );

            m_muonKeys.push_back( mu_key );
        }
    } // end loop over muons


    // AK4 jets
    m_jet_pt.clear();
    m_jet_eta.clear();
    m_jet_phi.clear();
    m_jet_mass.clear();
    m_jet_area.clear();
    m_jet_bdisc.clear();
    m_jet_deepCSV.clear();
    m_jet_uncorrPt.clear();
    m_jet_uncorrE.clear();
    m_jet_jerSF.clear();
    m_jet_jerSF_UP.clear();
    m_jet_jerSF_DOWN.clear();

    m_HTAK4 = 0;

    for (const auto& jet : *m_jets.product()){
        // lepton-cleaning
        // look at mediumID muons & tightID electrons (both w/o isolation)
        // re-check the 'loose' requirement
        CleanJet cjet = leptonJetCleaning( jet );

        if (!passAK4(cjet)) continue;                       // check kinematics (cjet should be the same as jet if not cleaned)

        // save variables -- if the jet isn't cleaned, the cleanJet.p4 matches the jet 4-vector
        m_jet_pt.push_back(   cjet.p4.Pt() );      // jet.pt()
        m_jet_eta.push_back(  cjet.p4.Eta() );     // jet.eta()
        m_jet_phi.push_back(  cjet.p4.Phi() );     // jet.phi()
        m_jet_mass.push_back( cjet.p4.M() );       // jet.mass()

        m_jet_bdisc.push_back(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
        //double deepCSV_b  = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        //double deepCSV_bb = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        //m_jet_deepCSV.push_back( deepCSV_b+deepCSV_bb );   // b+bb
        m_jet_deepCSV.push_back( jet.bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") );

        m_jet_area.push_back( jet.jetArea() );
        m_jet_uncorrPt.push_back( cjet.uncorrP4.Pt());   // same as jet.correctedP4(0)
        m_jet_uncorrE.push_back(  cjet.uncorrP4.E());    // same as jet.correctedP4(0)

        JME::JetParameters parameters = {{JME::Binning::JetEta, cjet.p4.Eta()}, {JME::Binning::Rho, m_rho}};
        m_jet_jerSF.push_back( m_resolution_ak4sf.getScaleFactor(parameters) );
        m_jet_jerSF_UP.push_back( m_resolution_ak4sf.getScaleFactor(parameters, Variation::UP) );
        m_jet_jerSF_DOWN.push_back( m_resolution_ak4sf.getScaleFactor(parameters, Variation::DOWN) );

        m_HTAK4 += cjet.p4.Pt();
    }

    // Fill output tree
    m_ttree->Fill();

    return;
}


void EventSaverFlatNtuple::endJob(){
    /* End of job 
       - Fill the metadata tree (only 1 "event")
    */
    m_xsection = 1;
    m_kfactor  = 1;
    m_NEvents  = 1;
    m_sumOfWeights = 1;

    if (m_isMC && m_mapOfSamples.find(m_sampleName)!=m_mapOfSamples.end()){
        Sample ss = m_mapOfSamples.at(m_sampleName);
        m_xsection = ss.XSection;
        m_kfactor  = ss.KFactor;
        m_sumOfWeights = ss.sumOfWeights;
        m_NEvents = ss.NEvents;
    }

    m_metadata_ttree->Fill();

    return;
}


void EventSaverFlatNtuple::initialize_branches(){
    /* Setup the output trees */
    // Metadata
    m_metadata_ttree->Branch("primaryDataset", &m_sampleName);                      // string
    m_metadata_ttree->Branch("xsection",       &m_xsection,     "xsection/F");      // float
    m_metadata_ttree->Branch("kfactor",        &m_kfactor,      "kfactor/F");       // float
    m_metadata_ttree->Branch("sumOfWeights",   &m_sumOfWeights, "sumOfWeights/F");  // float
    m_metadata_ttree->Branch("NEvents",        &m_NEvents,      "NEvents/i");       // uint

    // Event information
    //    Don't save flags/triggers (already used in MiniAOD!)
    m_ttree->Branch("runNumber",   &m_runNumber,   "runNumber/i");      // uint
    m_ttree->Branch("eventNumber", &m_eventNumber, "eventNumber/l");    // ulong
    m_ttree->Branch("lumiblock",   &m_lumiblock,   "lumiblock/i");      // uint
    m_ttree->Branch("rho",         &m_rho,         "rho/F");            // float
    m_ttree->Branch("npv",         &m_npv,         "npv/i");            // uint
    m_ttree->Branch("true_pileup", &m_true_pileup, "true_pileup/i");    // uint

    // Triggers
    for (const auto& name : m_triggers){
        m_triggerBits[name] = 0;
        m_ttree->Branch(name.c_str(),  &m_triggerBits[name], (name+"/i").c_str());  // uint
    }

    // Filters
    for (const auto& name : m_filters){
        m_filterBits[name] = 0;
        m_ttree->Branch(name.c_str(),  &m_filterBits[name],  (name+"/i").c_str());  // uint
    }

    // -- AK8 Jets
    m_ttree->Branch("AK8pt",     &m_ljet_pt);     // vector of floats
    m_ttree->Branch("AK8eta",    &m_ljet_eta);    // vector of floats
    m_ttree->Branch("AK8phi",    &m_ljet_phi);    // vector of floats
    m_ttree->Branch("AK8mass",   &m_ljet_mass);   // vector of floats
    m_ttree->Branch("AK8area",   &m_ljet_area);   // vector of floats
    m_ttree->Branch("AK8charge", &m_ljet_charge); // vector of floats
    m_ttree->Branch("AK8chargeSD", &m_ljet_chargeSD); // vector of floats
    m_ttree->Branch("AK8charge3",  &m_ljet_charge3);  // vector of floats
    m_ttree->Branch("AK8charge10", &m_ljet_charge10); // vector of floats
    m_ttree->Branch("AK8SDmass", &m_ljet_SDmass); // vector of floats
    m_ttree->Branch("AK8tau1",   &m_ljet_tau1);   // vector of floats
    m_ttree->Branch("AK8tau2",   &m_ljet_tau2);   // vector of floats
    m_ttree->Branch("AK8tau3",   &m_ljet_tau3);   // vector of floats
    m_ttree->Branch("AK8subjet0pt",     &m_ljet_subjet0_pt);     // vector of floats
    m_ttree->Branch("AK8subjet0mass",   &m_ljet_subjet0_mass);   // vector of floats
    m_ttree->Branch("AK8subjet0tau1",   &m_ljet_subjet0_tau1);   // vector of floats
    m_ttree->Branch("AK8subjet0tau2",   &m_ljet_subjet0_tau2);   // vector of floats
    m_ttree->Branch("AK8subjet0tau3",   &m_ljet_subjet0_tau3);   // vector of floats
    m_ttree->Branch("AK8subjet0bDisc",  &m_ljet_subjet0_bdisc);  // vector of floats
    m_ttree->Branch("AK8subjet0deepCSV",&m_ljet_subjet0_deepCSV);// vector of floats
    m_ttree->Branch("AK8subjet0charge", &m_ljet_subjet0_charge); // vector of floats
    m_ttree->Branch("AK8subjet0charge3",  &m_ljet_subjet0_charge3);  // vector of floats
    m_ttree->Branch("AK8subjet0charge10", &m_ljet_subjet0_charge10); // vector of floats
    m_ttree->Branch("AK8subjet1pt",     &m_ljet_subjet1_pt);     // vector of floats
    m_ttree->Branch("AK8subjet1mass",   &m_ljet_subjet1_mass);   // vector of floats
    m_ttree->Branch("AK8subjet1tau1",   &m_ljet_subjet1_tau1);   // vector of floats
    m_ttree->Branch("AK8subjet1tau2",   &m_ljet_subjet1_tau2);   // vector of floats
    m_ttree->Branch("AK8subjet1tau3",   &m_ljet_subjet1_tau3);   // vector of floats
    m_ttree->Branch("AK8subjet1bDisc",  &m_ljet_subjet1_bdisc);  // vector of floats
    m_ttree->Branch("AK8subjet1deepCSV",&m_ljet_subjet1_deepCSV);// vector of floats
    m_ttree->Branch("AK8subjet1charge", &m_ljet_subjet1_charge); // vector of floats
    m_ttree->Branch("AK8subjet1charge3",  &m_ljet_subjet1_charge3);  // vector of floats
    m_ttree->Branch("AK8subjet1charge10", &m_ljet_subjet1_charge10); // vector of floats
    m_ttree->Branch("AK8BEST_t", &m_ljet_BEST_t); // vector of floats
    m_ttree->Branch("AK8BEST_w", &m_ljet_BEST_w); // vector of floats
    m_ttree->Branch("AK8BEST_z", &m_ljet_BEST_z); // vector of floats
    m_ttree->Branch("AK8BEST_h", &m_ljet_BEST_h); // vector of floats
    m_ttree->Branch("AK8BEST_j", &m_ljet_BEST_j); // vector of floats
    m_ttree->Branch("AK8BEST_class", &m_ljet_BEST_class); // vector of floats
    m_ttree->Branch("AK8uncorrPt", &m_ljet_uncorrPt);     // vector of floats
    m_ttree->Branch("AK8uncorrE",  &m_ljet_uncorrE);      // vector of floats
    m_ttree->Branch("AK8jerSF",    &m_jet_jerSF);         // vector of floats
    m_ttree->Branch("AK8jerSF_UP",   &m_jet_jerSF_UP);    // vector of floats
    m_ttree->Branch("AK8jerSF_DOWN", &m_jet_jerSF_DOWN);  // vector of floats

    m_ttree->Branch("AK8truth_pt",     &m_truth_ljet_pt);     // vector of floats
    m_ttree->Branch("AK8truth_eta",    &m_truth_ljet_eta);    // vector of floats
    m_ttree->Branch("AK8truth_phi",    &m_truth_ljet_phi);    // vector of floats
    m_ttree->Branch("AK8truth_mass",   &m_truth_ljet_mass);   // vector of floats
    m_ttree->Branch("AK8truth_charge", &m_truth_ljet_charge); // vector of floats
    //m_ttree->Branch("AK8truth_SDmass", &m_truth_ljet_SDmass); // vector of floats  // not available right now (19 April 2018)
    m_ttree->Branch("AK8truth_tau1",   &m_truth_ljet_tau1);   // vector of floats
    m_ttree->Branch("AK8truth_tau2",   &m_truth_ljet_tau2);   // vector of floats
    m_ttree->Branch("AK8truth_tau3",   &m_truth_ljet_tau3);   // vector of floats

    // Physics Objects
    // -- AK4 Jets
    m_ttree->Branch("AK4pt",   &m_jet_pt);     // vector of floats
    m_ttree->Branch("AK4eta",  &m_jet_eta);    // vector of floats
    m_ttree->Branch("AK4phi",  &m_jet_phi);    // vector of floats
    m_ttree->Branch("AK4mass", &m_jet_mass);   // vector of floats
    m_ttree->Branch("AK4area", &m_jet_area);   // vector of floats
    m_ttree->Branch("AK4bDisc",&m_jet_bdisc);  // vector of floats
    m_ttree->Branch("AK4deepCSV",  &m_jet_deepCSV);   // vector of floats
    m_ttree->Branch("AK4uncorrPt", &m_jet_uncorrPt);  // vector of floats
    m_ttree->Branch("AK4uncorrE",  &m_jet_uncorrE);   // vector of floats
    m_ttree->Branch("AK4jerSF",    &m_jet_jerSF);     // vector of floats
    m_ttree->Branch("AK4jerSF_UP",   &m_jet_jerSF_UP);   // vector of floats
    m_ttree->Branch("AK4jerSF_DOWN", &m_jet_jerSF_DOWN); // vector of floats

    // -- Leptons (electrons & muons)
    m_ttree->Branch("ELpt",    &m_el_pt);     // vector of floats
    m_ttree->Branch("ELeta",   &m_el_eta);    // vector of floats
    m_ttree->Branch("ELphi",   &m_el_phi);    // vector of floats
    m_ttree->Branch("ELenergy",&m_el_e);      // vector of floats
    //m_ttree->Branch("ELiso",   &m_el_iso);    // vector of ints
    m_ttree->Branch("ELcharge",   &m_el_charge);        // vector of floats
    m_ttree->Branch("ELlooseID",  &m_el_ID_loose);      // vector of ints
    m_ttree->Branch("ELmediumID", &m_el_ID_medium);     // vector of ints
    m_ttree->Branch("ELtightID",  &m_el_ID_tight);      // vector of ints
    m_ttree->Branch("ELlooseIDnoIso",  &m_el_ID_looseNoIso);      // vector of ints
    m_ttree->Branch("ELmediumIDnoIso", &m_el_ID_mediumNoIso);     // vector of ints
    m_ttree->Branch("ELtightIDnoIso",  &m_el_ID_tightNoIso);      // vector of ints

    m_ttree->Branch("MUpt",    &m_mu_pt);     // vector of floats
    m_ttree->Branch("MUeta",   &m_mu_eta);    // vector of floats
    m_ttree->Branch("MUphi",   &m_mu_phi);    // vector of floats
    m_ttree->Branch("MUenergy",&m_mu_e);      // vector of floats
    m_ttree->Branch("MUcharge",   &m_mu_charge);        // vector of floats
    m_ttree->Branch("MUcorrIso",  &m_mu_iso);           // vector of floats
    m_ttree->Branch("MUlooseID",  &m_mu_ID_loose);      // vector of ints
    m_ttree->Branch("MUmediumID", &m_mu_ID_medium);     // vector of ints
    m_ttree->Branch("MUtightID",  &m_mu_ID_tight);      // vector of ints

    // MET, HT
    m_ttree->Branch("METpt",  &m_met_met, "METpt/F");  // float
    m_ttree->Branch("METphi", &m_met_phi, "METphi/F"); // float
    m_ttree->Branch("HTak8",  &m_HTAK8,   "HTak8/F");  // float
    m_ttree->Branch("HTak4",  &m_HTAK4,   "HTak4/F");  // float

    // Misc.
    // -- Generator-level information -- saving TOP only (30 March 2018)
    m_ttree->Branch("GENpt",       &m_mc_pt);       // vector of floats
    m_ttree->Branch("GENeta",      &m_mc_eta);      // vector of floats
    m_ttree->Branch("GENphi",      &m_mc_phi);      // vector of floats
    m_ttree->Branch("GENenergy",   &m_mc_e);        // vector of floats
    m_ttree->Branch("GENid",       &m_mc_pdgId);    // vector of ints
    m_ttree->Branch("GENstatus",   &m_mc_status);   // vector of ints
    m_ttree->Branch("GENparent_idx", &m_mc_parent_idx);   // vector of ints
    m_ttree->Branch("GENchild0_idx", &m_mc_child0_idx);   // vector of ints
    m_ttree->Branch("GENchild1_idx", &m_mc_child1_idx);   // vector of ints
    m_ttree->Branch("GENisHadTop", &m_mc_isHadTop); // vector of ints

    return;
}


bool EventSaverFlatNtuple::checkTopDecay(const reco::Candidate& daughter) const{
    /* Check the decay type of top quark by checking W decay products */
    auto w_child = daughter.daughter(0);
    int mode(0);

    if (std::abs(w_child->pdgId())==24 and w_child->numberOfDaughters()>0)
        mode = w_child->daughter(0)->pdgId();     // w decays to itself, check its decays
    else
        mode = w_child->pdgId();                  // w decays normally

    return (std::abs(mode)<10);
}


bool EventSaverFlatNtuple::passAK4( const CleanJet& cjet ) const {
    /* Check if small-R jet passes basic cuts (after cleaning!) */
    bool pass = (cjet.isLoose && cjet.p4.Pt()>15 && std::abs( cjet.p4.Eta() )<2.4);
    return pass;
}


bool EventSaverFlatNtuple::passAK8( const pat::Jet& j, const float& SDmass) const{
    /* Check if large-R jet passes basic cuts */
    bool goodJet = jetID(j);
    bool pass    = (j.pt() > 350. && goodJet && SDmass>20);
 
    return pass;
}


bool EventSaverFlatNtuple::jetID( const CleanJet& j ) const {
    /* Check Jet ID (loose) for ~cleaned~ AK4
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016
    */
    float inv_energy = 1. / j.uncorrP4.E();            // different jet energy after cleaning
    float nhf = j.neutralHadronEnergy * inv_energy;    // same as original jet
    float nef = j.neutralEmEnergy     * inv_energy;    // same as original jet
    float chf = j.chargedHadronEnergy * inv_energy;    // same as original jet
    float cef = j.chargedEmEnergy     * inv_energy;    // different if electron cleaned
    float nch = j.chargedMultiplicity;                 // different for muon & electron cleaning
    int nconstituents = j.numberOfDaughters;           // different for muon & electron cleaning

    bool goodJet = 
        nhf < m_nhfLoose &&
        nef < m_nefLoose &&
        chf > m_chfLoose &&
        cef < m_cefLoose &&
        nconstituents > m_nconstitLoose &&
        nch > m_nchLoose;

    return goodJet;
}


bool EventSaverFlatNtuple::jetID( const pat::Jet& j ) const {
    /* Check Jet ID (loose) 
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016
    */
    float nhf = j.neutralHadronEnergyFraction();
    float nef = j.neutralEmEnergyFraction();
    float chf = j.chargedHadronEnergyFraction();
    float cef = j.chargedEmEnergyFraction();
    float nch = j.chargedMultiplicity();
    int nconstituents = nch + j.neutralMultiplicity(); 

    bool goodJet = 
        nhf < m_nhfLoose &&
        nef < m_nefLoose &&
        chf > m_chfLoose &&
        cef < m_cefLoose &&
        nconstituents > m_nconstitLoose &&
        nch > m_nchLoose;

    return goodJet;
}


std::vector<float> EventSaverFlatNtuple::charge( const reco::Jet& jet, const std::vector<float> kappas, const unsigned int first_idx ) const {
    /* Calculate the jet charge
       Nest the kappa loop inside the jet consistituents loop
       to calculate jet charge for multiple kappa values -- might be less resource
       intensive than looping over the jet constituents multiple times...?
       !! DO NOT CALL THIS DIRECTLY FOR AN AK8 JET, IT WILL NOT ACCESS ALL DAUGHTERS !!
    */
    unsigned int size = jet.numberOfDaughters();
    unsigned int n_kappas = kappas.size();

    std::vector<float> charges;
    charges.resize( n_kappas,0.0 );

    for (unsigned int jd=first_idx; jd<size; jd++){
        const reco::CandidatePtr& dp = jet.daughterPtr(jd);
        if ( !(dp.isNonnull() && dp.isAvailable()) ) continue;    // quality check on constituent

        float child_q  = dp->charge();
        float child_pt = dp->pt();

        // sum the pt-weighted charges for different kappas
        for (unsigned int kap=0; kap<n_kappas; kap++)
            charges.at(kap) += child_q*pow(child_pt,kappas.at(kap));
    }
    // don't normalize by the total jet pt here, do that after calling this function!

    return charges;
}


std::vector<float> EventSaverFlatNtuple::getTau( unsigned int N, const reco::Jet& ij ) const {
    /* Calculate n-subjettiness for soft drop subjets.
       Loop over constituents of the soft drop subjet (just access the daughters) & convert to pseudojets.
       Calculate FJparticles only once, then calculate different taus (don't call this function repeatedly for the same jet!)
       Following existing setup in CMSSW for 94X (see links below).
       NB: It appears this can be done with the toolbox, but I'm not sure what else comes with that

        https://github.com/cms-jet/JetToolbox/blob/jetToolbox_94X/python/jetToolbox_cff.py#L798
        https://github.com/cms-sw/cmssw/blob/master/RecoJets/JetProducers/python/nJettinessAdder_cfi.py
        https://github.com/cms-sw/cmssw/blob/master/RecoJets/JetProducers/plugins/NjettinessAdder.cc
    */
    std::vector<float> taus;
    taus.resize(N,0.0);

    std::vector<fastjet::PseudoJet> FJparticles;
    for (unsigned k=0,size=ij.numberOfDaughters(); k<size; ++k){
        const reco::CandidatePtr & dp = ij.daughterPtr(k);
        if ( dp.isNonnull() && dp.isAvailable() )
            FJparticles.push_back( fastjet::PseudoJet( dp->px(), dp->py(), dp->pz(), dp->energy() ) );
    }

    // Calculate n-subjettiness (N=1,2,3,...)
    for (unsigned int nt=1; nt<N+1; nt++)
        taus.at(nt-1) = m_nsub->getTau(nt,FJparticles);

    return taus;
}


int EventSaverFlatNtuple::findPartonIndex( const std::vector<reco::GenParticle>& items, const reco::Candidate& item ) const{
    /* loop over particles (in the decay chain before this particle) to get parent/children
       -- Easier than comparing attributes?
    */
    int p0_idx(-1);
    for (unsigned int p0=0,size=items.size(); p0<size; p0++){
        // compare pdgId, status, charge
        if ( items.at(p0).charge() == item.charge() &&
             items.at(p0).status() == item.status() &&
             items.at(p0).pdgId()  == item.pdgId() ){
            p0_idx = p0;
            break;
        }
    } // end loop over truth particles

    return p0_idx;
}


CleanJet EventSaverFlatNtuple::leptonJetCleaning( const pat::Jet& j ) {
    /* Clean jet of candidate leptons 
       > If any of the lepton keys are in the AK4, 
         subtract the lepton four vector from the AK4
    */
    reco::Candidate::LorentzVector uncorrJet = j.correctedP4(0);

    CleanJet cleanJet;
    cleanJet.p4.SetPtEtaPhiE( j.pt(), j.eta(), j.phi(), j.energy() );  // initialize to nominal jet
    cleanJet.uncorrP4.SetPtEtaPhiE( uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), uncorrJet.energy() );     // set to 4-vector WITHOUT JECs
    cleanJet.isLoose  = jetID(j);                 // nominal jet ID

    // save the lepton p4 if it is used in cleaning for MET re-calculation
    std::vector<TLorentzVector> cleaningLeptons;  // may be more than one, but not likely (?)

    // check constituents of this jet
    bool clean(false);
    float chargedEMEnergy_clean(0.0);             // count the charged EM energy for cleaned constituents
    float chargedMultiplicity_clean(0.0);         // count the charged multiplicity for cleaned constituents

    auto constituents = j.daughterPtrVector();
    for ( auto & constituent : constituents ) {
        int key = constituent.key();

        // check muons
        for (const auto mukey : m_muonKeys){
            for (const auto& k : mukey.keys){
                if (key == k){
                    clean = true;
                    cleanJet.uncorrP4 -= mukey.p4;
                    chargedMultiplicity_clean++;
                    cleaningLeptons.push_back( mukey.p4 );
                    break; // only subtract 4-vector if at least one of the keys match
                }
            }
        } // end loop over muon candidates

        // check electrons
        for (const auto& elkey : m_electronKeys){
            for (const auto& k : elkey.keys){
                if (key == k){
                    clean = true;
                    cleanJet.uncorrP4 -= elkey.p4;
                    chargedEMEnergy_clean += elkey.p4.E();
                    chargedMultiplicity_clean++;
                    cleaningLeptons.push_back( elkey.p4 );
                    break; // only subtract 4-vector if at least one of the keys match
                }
            }
        } // end loop over electron candidates
    } // end loop over jet constituents

    if (clean) {
        // recalculate 'isLoose' -- update PF information of jet
        cleanJet.neutralHadronEnergy = j.neutralHadronEnergy();
        cleanJet.neutralEmEnergy     = j.neutralEmEnergy();
        cleanJet.chargedHadronEnergy = j.chargedHadronEnergy();
        cleanJet.chargedEmEnergy     = j.chargedEmEnergy()     - chargedEMEnergy_clean;         // electrons
        cleanJet.chargedMultiplicity = j.chargedMultiplicity() - chargedMultiplicity_clean;     // electrons+muons
        cleanJet.numberOfDaughters   = cleanJet.chargedMultiplicity + j.neutralMultiplicity();  // electrons+muons

        cleanJet.isLoose = jetID(cleanJet);   // uses uncorrected jet energy

        // re-apply JECs
        applyJEC( cleanJet, j.jetArea() );  // update p4 attribute

        // update MET
        updateMET( cleanJet, cleaningLeptons );
    }

    return cleanJet;
}


void EventSaverFlatNtuple::updateMET( const CleanJet& j, const std::vector<TLorentzVector>& leptons4cleaning ){
    /* Update the MET based on lepton-jet cleaning */
    TLorentzVector met;
    met.SetPtEtaPhiE( m_met_met, 0., m_met_phi, 0.);

    met += j.originalP4;
    met -= j.p4;

    for (const auto& l : leptons4cleaning)
        met -= l;

    m_met_met = met.Pt();    // new corrected result
    m_met_phi = met.Phi();   // new corrected result

    return;
}


void EventSaverFlatNtuple::applyJEC( CleanJet& j, const float& area ) const {
    /* Apply JECs to jet (nominally for jet cleaning) 
        > IOV BCD: [1,276811]        corresponds to Summer16_23Sep2016BCDV4_DATA (For Runs B/C/D)
        > IOV EF:  [276831,278801]   corresponds to Summer16_23Sep2016EFV4_DATA (For Runs E/early F)
        > IOV G:   [278802,280385]   corresponds to Summer16_23Sep2016GV4_DATA (For Runs lateF/G)
        > IOV H:   [280919,Infinity] corresponds to Summer16_23Sep2016HV4_DATA (For Run H)
    */
    std::string key("");
    if (m_isMC) key = "V4_MC";
    else if (m_runNumber>1 && m_runNumber<=276811) key = "BCDV4_DATA";
    else if (m_runNumber>=276831 && m_runNumber<=278801) key = "EFV4_DATA";
    else if (m_runNumber>=278802 && m_runNumber<=280385) key = "GV4_DATA";
    else if (m_runNumber>=280919) key = "HV4_DATA";

    m_ak4_jec.at(key)->setJetPt(  j.uncorrP4.Pt() );
    m_ak4_jec.at(key)->setJetEta( j.uncorrP4.Eta() );
    m_ak4_jec.at(key)->setJetE(   j.uncorrP4.Energy() );
    m_ak4_jec.at(key)->setJetA(area );
    m_ak4_jec.at(key)->setRho( m_rho );
    m_ak4_jec.at(key)->setNPV( m_npv );

    float JECfactor = m_ak4_jec.at(key)->getCorrection();
    j.originalP4.SetPtEtaPhiE( j.p4.Pt(), j.p4.Eta(), j.p4.Phi(), j.p4.E() );
    j.p4 = j.uncorrP4 * JECfactor;                  // scale 4-vector by correction

    return;
}


DEFINE_FWK_MODULE(EventSaverFlatNtuple);
// THE END
