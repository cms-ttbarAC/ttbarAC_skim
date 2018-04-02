/*
Created:        29 March    2018
Last Updated:   30 March    2018

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
  t_electrons(consumes<std::vector<pat::Electron>>(edm::InputTag("selectedElectrons", "", "ttbarACskim"))),
  t_electrons_orig(consumes<edm::View<pat::Electron>>(edm::InputTag("slimmedElectrons"))),
  t_jets(consumes<pat::JetCollection>(edm::InputTag("selectedAK4Jets", "", "ttbarACskim"))),
  t_ljets(consumes<pat::JetCollection>(edm::InputTag("BESTProducer", "savedJets", "ttbarACskim"))),
  t_truth_ljets(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"))),
  t_met(consumes<pat::METCollection>(edm::InputTag("selectedMET", "", "ttbarACskim"))),
  t_genEvtInfoProd(consumes<std::vector<reco::GenParticle>>(edm::InputTag("selectedGenParticles", "", "ttbarACskim"))),
  t_rho(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  t_vertices(consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"))),
  t_pileup(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
  t_beamspot(consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"))),
  t_conversions(consumes<reco::ConversionCollection>(edm::InputTag("reducedEgamma:reducedConversions"))),
  t_elIdFullInfoMap_Loose(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_Loose"))),
  t_elIdFullInfoMap_Medium(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_Medium"))),
  t_elIdFullInfoMap_Tight(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_Tight"))){
  //t_elIdFullInfoMap_HEEP(consumes<edm::ValueMap<vid::CutFlowResult>>(cfg.getParameter<edm::InputTag>("elIdFullInfoMap_HEEP"))){
    t_BEST_products.clear();
    for (const auto& name : m_BEST_variables){
        edm::EDGetTokenT<std::vector<float>> tmp_token = consumes<std::vector<float>>(edm::InputTag("BESTProducer",name,"ttbarACskim"));
        t_BEST_products.push_back( tmp_token );
    }

    // Make output TTrees
    edm::Service<TFileService> fs;
    m_ttree = fs->make<TTree>("eventVars","eventVars");
    m_metadata_ttree = fs->make<TTree>("metadata","metadata");

    m_hist_cutflow = fs->make<TH1D>( "cutflow","cutflow",3,0,3);
    m_hist_cutflow->GetXaxis()->SetBinLabel(1,"INITIAL");
    m_hist_cutflow->GetXaxis()->SetBinLabel(2,"PRIMARYVTX");
    m_hist_cutflow->GetXaxis()->SetBinLabel(3,"AK8JETS");

    initialize_branches();

    // options set by config
    m_isMC = cfg.getParameter<bool>("isMC");           // filling truth branches

    cma::setVerboseLevel("WARNING");
    m_sampleName = t_sampleName;
    m_XSections.clear();
    m_KFactors.clear();
    m_sumOfMCWeights.clear();
    m_NEvents.clear();
    if (m_isMC){
        m_mapOfSamples.clear();
        std::cout << " Get sample weights " << std::endl;
        cma::getSampleWeights( t_metadataFile,m_mapOfSamples );
    }

    // Lightweight NN interface with BEST
    std::string dnnFile("../test/BEST_mlp.json");
    std::ifstream input_cfg( dnnFile );                     // original: "data/BEST_mlp.json"
    lwt::JSONConfig lwcfg = lwt::parse_json( input_cfg );
    m_lwtnn = new lwt::LightweightNeuralNetwork(lwcfg.inputs, lwcfg.layers, lwcfg.outputs);
}


EventSaverFlatNtuple::~EventSaverFlatNtuple() {
    delete m_lwtnn;
}



void EventSaverFlatNtuple::beginJob(){
    /* Begin Job */
    return;
}


void EventSaverFlatNtuple::analyze( const edm::Event& event, const edm::EventSetup& ) {
    /* Fill TTree 
       This is the function to modify / inherit for analysis-specific purposes
    */
    // Load data 
    event.getByToken( t_electrons, m_electrons );
    event.getByToken( t_muons, m_muons );
    event.getByToken( t_jets, m_jets );
    //event.getByToken( t_truth_jets, m_truth_jets );
    event.getByToken( t_truth_ljets, m_truth_ljets );
    event.getByToken( t_met, m_met );
    event.getByToken( t_rho, h_rho );
    event.getByToken( t_genEvtInfoProd,h_genEvtInfoProd );
    event.getByToken( t_vertices,h_vertices );
    event.getByToken( t_pileup,h_pileup );

    //event.getByToken( t_beamspot,h_beamspot );
    //event.getByToken( t_conversions,h_conversions );
    event.getByToken( t_elIdFullInfoMap_Loose,  h_cutflow_elId_Loose );
    event.getByToken( t_elIdFullInfoMap_Medium, h_cutflow_elId_Medium );
    event.getByToken( t_elIdFullInfoMap_Tight,  h_cutflow_elId_Tight );
//    event.getByToken( t_elIdFullInfoMap_HEEP,   h_cutflow_elId_HEEP );

    m_hist_cutflow->Fill(0.5);  // INITIAL

    // Set branch values
    m_runNumber   = event.id().run();
    m_eventNumber = event.id().event();
    m_lumiblock   = event.id().luminosityBlock();
    m_rho         = *h_rho;

    // -- primary vertex information
    m_npv = h_vertices->size();                    // primary vertices
    if (h_vertices->empty()) return;               // skip the event if no PV found
    const reco::Vertex &PV = h_vertices->front();  // save PV for tight muon ID

    m_hist_cutflow->Fill(1.5);   // Primary Vertex

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

    // AK4 jets
    m_jet_pt.clear();
    m_jet_eta.clear();
    m_jet_phi.clear();
    m_jet_mass.clear();
    m_jet_bdisc.clear();

    for (const auto& jet : *m_jets.product()){
        m_jet_pt.push_back(   jet.pt() );
        m_jet_eta.push_back(  jet.eta() );
        m_jet_phi.push_back(  jet.phi() );
        m_jet_mass.push_back( jet.mass() );
        m_jet_bdisc.push_back(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
    }

    // AK8 jets
    m_ljet_pt.clear();
    m_ljet_eta.clear();
    m_ljet_phi.clear();
    m_ljet_mass.clear();
    m_ljet_tau1.clear();
    m_ljet_tau2.clear();
    m_ljet_tau3.clear();
    m_ljet_charge.clear();
    m_ljet_SDmass.clear();
    m_ljet_BEST_t.clear();
    m_ljet_BEST_w.clear();
    m_ljet_BEST_z.clear();
    m_ljet_BEST_h.clear();
    m_ljet_BEST_j.clear();
    m_ljet_BEST_class.clear();
    m_ljet_subjet0_pt.clear();
    m_ljet_subjet0_mass.clear();
    m_ljet_subjet0_bdisc.clear();
    m_ljet_subjet0_charge.clear();
    m_ljet_subjet1_pt.clear();
    m_ljet_subjet1_mass.clear();
    m_ljet_subjet1_bdisc.clear();
    m_ljet_subjet1_charge.clear();
    m_HTAK8 = 0;

    // May not have the BEST products
    // -- https://twiki.cern.ch/twiki/bin/view/Main/CMSSWCheatSheet#Does_a_certain_product_exist_in
    try {event.getByToken( t_ljets, m_ljets );}
    catch( cms::Exception& ex ) {std::cout << "Large-R Jets not found" << std::endl;}
    if (!m_ljets.isValid()) {
        std::cout << " Product not valid: Large-R Jets (none remaining after BEST) " << std::endl;
        return;
    }
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
        bool pass = passAK8( ljet,nj,m_BEST_products.at("AK8SDmass")[nj] );
        if (!pass){
            std::cout << " AK8 failed; pT = " << ljet.pt() << "; sdmass = " << m_BEST_products.at("AK8SDmass")[nj] << std::endl;
            nj++;
            continue;
        }

        m_ljet_pt.push_back(  ljet.pt());
        m_ljet_eta.push_back( ljet.eta());
        m_ljet_phi.push_back( ljet.phi());
        m_ljet_mass.push_back(ljet.mass());

        m_ljet_tau1.push_back(ljet.userFloat("NjettinessAK8:tau1"));
        m_ljet_tau2.push_back(ljet.userFloat("NjettinessAK8:tau2"));
        m_ljet_tau3.push_back(ljet.userFloat("NjettinessAK8:tau3"));

        // BESTProducer
        m_ljet_charge.push_back( m_BEST_products.at("AK8charge")[nj] );
        m_ljet_SDmass.push_back( m_BEST_products.at("AK8SDmass")[nj] );

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

        std::cout << " LWTNN " << std::endl;
        std::map<std::string,double> NNresults = m_lwtnn->compute(BEST_products);
        std::vector<double> values{ NNresults["dnn_qcd"],   NNresults["dnn_top"],
                                    NNresults["dnn_higgs"], NNresults["dnn_z"], NNresults["dnn_w"] };
        std::cout << " end LWTNN " << std::endl;

        unsigned int particleID(0);
        float max_value(-1.0);
        for (unsigned int pid=0,size=values.size();pid<size;pid++){
            if (values.at(pid) > max_value){
                max_value  = values.at(pid);
                particleID = pid;
            }
        }

        std::cout << " set LWTNN " << std::endl;
        m_ljet_BEST_t.push_back( NNresults.at("dnn_top") );
        m_ljet_BEST_w.push_back( NNresults.at("dnn_w") );
        m_ljet_BEST_z.push_back( NNresults.at("dnn_z") );
        m_ljet_BEST_h.push_back( NNresults.at("dnn_higgs") );
        m_ljet_BEST_j.push_back( NNresults.at("dnn_qcd") );
        m_ljet_BEST_class.push_back( particleID );

        std::cout << " set subjets " << std::endl;
        for (const auto& x : m_BEST_products)
            std::cout << x.first << std::endl;
        m_ljet_subjet0_pt.push_back(     m_BEST_products["AK8subjet0pT"][nj]);
        m_ljet_subjet0_mass.push_back(   m_BEST_products["AK8subjet0mass"][nj]);
        m_ljet_subjet0_bdisc.push_back(  m_BEST_products["AK8subjet0bDisc"][nj] );
        m_ljet_subjet0_charge.push_back( m_BEST_products["AK8subjet0charge"][nj] );
        m_ljet_subjet1_pt.push_back(     m_BEST_products["AK8subjet1pT"][nj]);
        m_ljet_subjet1_mass.push_back(   m_BEST_products["AK8subjet1mass"][nj]);
        m_ljet_subjet1_bdisc.push_back(  m_BEST_products["AK8subjet1bDisc"][nj] );
        m_ljet_subjet1_charge.push_back( m_BEST_products["AK8subjet1charge"][nj] );

        m_HTAK8+=ljet.pt();
        nj++;
        std::cout << " next " << std::endl;
    } // end loop over AK8
    std::cout << " END AK8 " << std::endl;

    if (m_ljet_pt.size()<1) return;
    m_hist_cutflow->Fill(2.5);    // AK8Jets



    // Leptons
    edm::Handle<edm::View<pat::Electron>> electrons; //Collection
    event.getByToken( t_electrons_orig, electrons );
    for (size_t i = 0; i < electrons->size(); ++i){   
        const auto el = electrons->ptrAt(i);          // easier if we use ptrs for the id
        if (el->pt()<40 || fabs(el->eta())>2.4 ) continue;

        std::cout << i << ": " << el->pt() << std::endl;

        vid::CutFlowResult idLoose  = (*h_cutflow_elId_Loose)[el];
        vid::CutFlowResult idMedium = (*h_cutflow_elId_Medium)[el];
        vid::CutFlowResult idTight  = (*h_cutflow_elId_Tight)[el];
//        vid::CutFlowResult idHEEP   = (*h_cutflow_elId_HEEP)[el];
    }

    m_el_pt.clear();
    m_el_eta.clear();
    m_el_phi.clear();
    m_el_e.clear();
    m_el_iso.clear();
    m_el_charge.clear();
    m_el_ID_loose.clear();
    m_el_ID_medium.clear();
    m_el_ID_tight.clear();
    //m_el_ID_HEEP.clear();

    unsigned int i(0);
    for (const auto& el : *m_electrons.product()){
        std::cout << i << ": " << el.pt() << std::endl;

        m_el_pt.push_back( el.pt() );
        m_el_eta.push_back(el.eta() );
        m_el_phi.push_back(el.phi() );
        m_el_e.push_back(  el.energy() );

        m_el_charge.push_back( el.charge() );
        m_el_iso.push_back( (el.trackIso() + el.caloIso()) / el.pt() );

        // ID
//        m_el_ID_loose.push_back( idLoose.cutFlowPassed() );
//        m_el_ID_medium.push_back(idMedium.cutFlowPassed() );
//        m_el_ID_tight.push_back( idTight.cutFlowPassed() );
//        m_el_ID_HEEP.push_back(  idHEEP.cutFlowPassed() );
        i++;
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

    for (const auto& mu : *m_muons.product()){
        m_mu_pt.push_back(  mu.pt() );
        m_mu_eta.push_back( mu.eta() );
        m_mu_phi.push_back( mu.phi() );
        m_mu_e.push_back(   mu.energy() );
        m_mu_charge.push_back( mu.charge() );

        // ISO04
        float chPt = mu.pfIsolationR04().sumChargedHadronPt;
        float nhPt = mu.pfIsolationR04().sumNeutralHadronEt;
        float phPt = mu.pfIsolationR04().sumPhotonEt;
        float puPt = mu.pfIsolationR04().sumPUPt;
        float iso04 = (chPt+TMath::Max(0.,nhPt+phPt-0.5*puPt))/mu.pt();
        m_mu_iso.push_back( iso04 );

        // Loose ID
        m_mu_ID_loose.push_back( muon::isLooseMuon(mu) );

        // Medium ID  -- 2016 data B-F vs G-H
        bool goodGlob   = mu.isGlobalMuon() && 
                          mu.globalTrack()->normalizedChi2() < 3 && 
                          mu.combinedQuality().chi2LocalPosition < 12 && 
                          mu.combinedQuality().trkKink < 20;
        bool isMedium   = muon::isMediumMuon(mu);
        bool isMediumBF = muon::isLooseMuon(mu) && 
                          mu.innerTrack()->validFraction() > 0.49 && 
                          muon::segmentCompatibility(mu) > (goodGlob ? 0.303 : 0.451); 
        bool isMediumMuon = (event.isRealData() && m_runNumber <= 278808) ? isMediumBF : isMedium;
        m_mu_ID_medium.push_back( isMediumMuon );

        m_mu_ID_tight.push_back( mu.isTightMuon(PV) );
    } // end loop over muons

    m_met_met = (*m_met.product())[0].pt();
    m_met_phi = (*m_met.product())[0].phi();

    if (m_isMC){
        m_mc_pt.clear();
        m_mc_eta.clear();
        m_mc_phi.clear();
        m_mc_e.clear();
        m_mc_pdgId.clear();
        m_mc_status.clear();
        m_mc_isHadTop.clear();

        for (const auto& particle: *h_genEvtInfoProd.product()){
            if (std::abs(particle.pdgId())==6 && particle.numberOfDaughters()==2){
                m_mc_pt.push_back(particle.pt());
                m_mc_eta.push_back(particle.eta());
                m_mc_phi.push_back(particle.phi());
                m_mc_e.push_back(particle.energy());
                m_mc_pdgId.push_back(particle.pdgId());
                m_mc_status.push_back(particle.status());

                auto* daughter1 = particle.daughter(0);
                auto* daughter2 = particle.daughter(1);
                if (std::abs(daughter1->pdgId()) == 24)
                    m_mc_isHadTop.push_back( checkTopDecay(*daughter1) );
                if (std::abs(daughter2->pdgId()) == 24)
                    m_mc_isHadTop.push_back( checkTopDecay(*daughter2) );
            }
        }
    } // end if isMC

    // Fill output tree
    std::cout << " Fill tree " << std::endl;
    m_ttree->Fill();

    return;
}


void EventSaverFlatNtuple::endJob(){
    /* End of job 
       - Fill the metadata tree (only 1 "event")
    */
    m_xsection = 1;
    m_kfactor  = 1;
    m_sumOfWeights = 1;
    if (m_isMC && m_mapOfSamples.find(m_sampleName)!=m_mapOfSamples.end()){
        Sample ss = m_mapOfSamples.at(m_sampleName);
        m_xsection = ss.XSection;
        m_kfactor  = ss.KFactor;
        m_sumOfWeights = ss.sumOfWeights;
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

    // Physics Objects
    // -- AK4 Jets
    m_ttree->Branch("AK4pt",   &m_jet_pt);     // vector of floats
    m_ttree->Branch("AK4eta",  &m_jet_eta);    // vector of floats
    m_ttree->Branch("AK4phi",  &m_jet_phi);    // vector of floats
    m_ttree->Branch("AK4mass", &m_jet_mass);   // vector of floats
    m_ttree->Branch("AK4bDisc",&m_jet_bdisc);  // vector of floats

    // -- AK8 Jets
    m_ttree->Branch("AK8pt",     &m_ljet_pt);     // vector of floats
    m_ttree->Branch("AK8eta",    &m_ljet_eta);    // vector of floats
    m_ttree->Branch("AK8phi",    &m_ljet_phi);    // vector of floats
    m_ttree->Branch("AK8mass",   &m_ljet_mass);   // vector of floats
    m_ttree->Branch("AK8charge", &m_ljet_charge); // vector of floats
    m_ttree->Branch("AK8SDmass", &m_ljet_SDmass); // vector of floats
    m_ttree->Branch("AK8tau1",   &m_ljet_tau1);   // vector of floats
    m_ttree->Branch("AK8tau2",   &m_ljet_tau2);   // vector of floats
    m_ttree->Branch("AK8tau3",   &m_ljet_tau3);   // vector of floats
    m_ttree->Branch("AK8PtSubjet1",     &m_ljet_subjet0_pt);     // vector of floats
    m_ttree->Branch("AK8MassSubjet1",   &m_ljet_subjet0_mass);   // vector of floats
    m_ttree->Branch("AK8bDiscSubjet1",  &m_ljet_subjet0_bdisc);  // vector of floats
    m_ttree->Branch("AK8ChargeSubjet1", &m_ljet_subjet0_charge); // vector of floats
    m_ttree->Branch("AK8PtSubjet2",     &m_ljet_subjet1_pt);     // vector of floats
    m_ttree->Branch("AK8MassSubjet2",   &m_ljet_subjet1_mass);   // vector of floats
    m_ttree->Branch("AK8bDiscSubjet2",  &m_ljet_subjet1_bdisc);  // vector of floats
    m_ttree->Branch("AK8ChargeSubjet2", &m_ljet_subjet1_charge); // vector of floats

    m_ttree->Branch("BESTProb_t", &m_ljet_BEST_t); // vector of floats
    m_ttree->Branch("BESTProb_w", &m_ljet_BEST_w); // vector of floats
    m_ttree->Branch("BESTProb_z", &m_ljet_BEST_z); // vector of floats
    m_ttree->Branch("BESTProb_h", &m_ljet_BEST_h); // vector of floats
    m_ttree->Branch("BESTProb_j", &m_ljet_BEST_j); // vector of floats
    m_ttree->Branch("BESTProb_class", &m_ljet_BEST_class); // vector of floats

    // -- Leptons (electrons & muons)
    m_ttree->Branch("ELpt",    &m_el_pt);     // vector of floats
    m_ttree->Branch("ELeta",   &m_el_eta);    // vector of floats
    m_ttree->Branch("ELphi",   &m_el_phi);    // vector of floats
    m_ttree->Branch("ELenergy",&m_el_e);      // vector of floats
    m_ttree->Branch("ELiso",   &m_el_iso);    // vector of ints
    m_ttree->Branch("ELcharge",   &m_el_charge);        // vector of floats
    m_ttree->Branch("ELlooseID",  &m_el_ID_loose);      // vector of ints
    m_ttree->Branch("ELmediumID", &m_el_ID_medium);     // vector of ints
    m_ttree->Branch("ELtightID",  &m_el_ID_tight);      // vector of ints

    m_ttree->Branch("MUpt",    &m_mu_pt);     // vector of floats
    m_ttree->Branch("MUeta",   &m_mu_eta);    // vector of floats
    m_ttree->Branch("MUphi",   &m_mu_phi);    // vector of floats
    m_ttree->Branch("MUenergy",&m_mu_e);      // vector of floats
    m_ttree->Branch("MUcharge",   &m_mu_charge);        // vector of floats
    m_ttree->Branch("MUcorrIso",  &m_mu_iso);           // vector of ints
    m_ttree->Branch("MUlooseID",  &m_mu_ID_loose);      // vector of ints
    m_ttree->Branch("MUmediumID", &m_mu_ID_medium);     // vector of ints
    m_ttree->Branch("MUtightID",  &m_mu_ID_tight);      // vector of ints

    // MET, HT
    m_ttree->Branch("METpt",  &m_met_met, "METpt/F");  // float
    m_ttree->Branch("METphi", &m_met_phi, "METphi/F"); // float
    m_ttree->Branch("HTak8",  &m_HTAK8,   "HTak8/F");  // float

    // Event information
    //    Don't save flags/triggers (already used in MiniAOD!)
    m_ttree->Branch("runNumber",   &m_runNumber,   "runNumber/i");      // uint
    m_ttree->Branch("eventNumber", &m_eventNumber, "eventNumber/l");    // ulong
    m_ttree->Branch("lumiblock",   &m_lumiblock,   "lumiblock/i");      // uint
    m_ttree->Branch("rho",         &m_rho,         "rho/F");            // float
    m_ttree->Branch("npv",         &m_npv,         "npv/I");            // int
    m_ttree->Branch("true_pileup", &m_true_pileup, "true_pileup/I");    // int

    // Misc.
    // -- Generator-level information -- saving TOP only (30 March 2018)
    m_ttree->Branch("GENpt",       &m_mc_pt);       // vector of floats
    m_ttree->Branch("GENeta",      &m_mc_eta);      // vector of floats
    m_ttree->Branch("GENphi",      &m_mc_phi);      // vector of floats
    m_ttree->Branch("GENenergy",   &m_mc_e);        // vector of floats
    m_ttree->Branch("GENid",       &m_mc_pdgId);    // vector of ints
    m_ttree->Branch("GENstatus",   &m_mc_status);   // vector of ints
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


bool EventSaverFlatNtuple::passAK8( const pat::Jet& j, const int index, const float& SDmass) const{
    /* Check if large-R jet passes basic cuts */
    bool pass(false);
    double pt_cut = (index==0) ? 350. : 300;

    float energy = j.energy();
    float nhf = j.neutralHadronEnergy() / energy;
    float nef = j.neutralEmEnergy() / energy;
    float chf = j.chargedHadronEnergy() / energy;
    float cef = j.chargedEmEnergy() / energy;
    float nch = j.chargedMultiplicity();
    int nconstituents = j.numberOfDaughters();
    bool goodJet = 
        nhf < 0.99 && 
        nef < 0.99 && 
        chf > 0.00 && 
        cef < 0.99 && 
        nconstituents > 1 &&
        nch > 0;

    if (j.pt() > pt_cut && goodJet && SDmass>20) pass = true;

    return pass;
}


DEFINE_FWK_MODULE(EventSaverFlatNtuple);
// THE END
