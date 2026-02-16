#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"

//---------------------------------------------------------------------
//  constructors and destructor
//---------------------------------------------------------------------

// constructor
V0Analyzer::V0Analyzer(const edm::ParameterSet& iConfig):
    beamSpotToken(              consumes<reco::BeamSpot>(                        iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
    vtxToken(                   consumes<std::vector<reco::Vertex>>(             iConfig.getParameter<edm::InputTag>("vertices"))),
    genEventInfoToken(          consumes<GenEventInfoProduct>(                   iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    lheEventInfoToken(          consumes<LHEEventProduct>(                       iConfig.getParameter<edm::InputTag>("lheEventInfo"))),
    pileUpToken(                consumes<std::vector<PileupSummaryInfo>>(        iConfig.getParameter<edm::InputTag>("pileUpInfo"))),
    genParticleToken(           consumes<reco::GenParticleCollection>(           iConfig.getParameter<edm::InputTag>("genParticles"))),
    particleLevelPhotonsToken(  consumes<reco::GenParticleCollection>(           iConfig.getParameter<edm::InputTag>("particleLevelPhotons"))),
    particleLevelLeptonsToken(  consumes<reco::GenJetCollection>(                iConfig.getParameter<edm::InputTag>("particleLevelLeptons"))),
    particleLevelJetsToken(     consumes<reco::GenJetCollection>(                iConfig.getParameter<edm::InputTag>("particleLevelJets"))),
    genJetsToken(               consumes<reco::GenJetCollection>(                iConfig.getParameter<edm::InputTag>("particleLevelJets"))),
   particleLevelMetsToken(     consumes<reco::METCollection>(                   iConfig.getParameter<edm::InputTag>("particleLevelMets"))),
    muonToken(                  consumes<std::vector<pat::Muon>>(                iConfig.getParameter<edm::InputTag>("muons"))),
    eleToken(                   consumes<std::vector<pat::Electron>>(            iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken(                   consumes<std::vector<pat::Tau>>(                 iConfig.getParameter<edm::InputTag>("taus"))),
    tauGenJetsToken(            consumes<std::vector<reco::GenJet>>(             iConfig.getParameter<edm::InputTag>("tauGenJets"))),
    photonToken(                consumes<std::vector<pat::Photon>>(              iConfig.getParameter<edm::InputTag>("photons"))),
    packedCandidatesToken(      consumes<std::vector<pat::PackedCandidate>>(     iConfig.getParameter<edm::InputTag>("packedCandidates"))),
    lostTracksToken(            consumes<std::vector<pat::PackedCandidate>>(     iConfig.getParameter<edm::InputTag>("lostTracks"))),
    //rhoToken(                   consumes<double>(                                iConfig.getParameter<edm::InputTag>("rho"))),
    metToken(                   consumes<std::vector<pat::MET>>(                 iConfig.getParameter<edm::InputTag>("met"))),
    metPuppiToken(              consumes<std::vector<pat::MET>>(                 iConfig.getParameter<edm::InputTag>("metPuppi"))),
    jetToken(                   consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jets"))),
    jetPuppiToken(              consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
    jetSmearedToken(            consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jetsSmeared"))),
    jetSmearedUpToken(          consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jetsSmearedUp"))),
    jetSmearedDownToken(        consumes<std::vector<pat::Jet>>(                 iConfig.getParameter<edm::InputTag>("jetsSmearedDown"))),
    recoResultsPrimaryToken(    consumes<edm::TriggerResults>(                   iConfig.getParameter<edm::InputTag>("recoResultsPrimary"))),
    recoResultsSecondaryToken(  consumes<edm::TriggerResults>(                   iConfig.getParameter<edm::InputTag>("recoResultsSecondary"))),
    triggerToken(               consumes<edm::TriggerResults>(                   iConfig.getParameter<edm::InputTag>("triggers"))),
    prescalesToken(             consumes<pat::PackedTriggerPrescales>(           iConfig.getParameter<edm::InputTag>("prescales"))),
    trigObjToken(               consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    //prefireWeightToken(         consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProb"))),
    //prefireWeightUpToken(       consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
    //prefireWeightDownToken(     consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
    //prefireWeightMuonToken(     consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbMuon"))),
    //prefireWeightMuonUpToken(   consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbMuonUp"))),
    //prefireWeightMuonDownToken( consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbMuonDown"))),
    //prefireWeightECALToken(     consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbECAL"))),
    //prefireWeightECALUpToken(   consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbECALUp"))),
    //prefireWeightECALDownToken( consumes<double>(                                edm::InputTag("prefiringweight:nonPrefiringProbECALDon"))),
    bFieldToken_(               esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    skim(                                                                        iConfig.getUntrackedParameter<std::string>("skim")),
    sampleIsData(                                                                iConfig.getUntrackedParameter<bool>("isData")),
    sampleIs2016preVFP(                                                          iConfig.getUntrackedParameter<bool>("is2016preVFP")),
    sampleIs2017(                                                                iConfig.getUntrackedParameter<bool>("is2017")),
    sampleIs2018(                                                                iConfig.getUntrackedParameter<bool>("is2018")),
    sampleIs2022(                                                                iConfig.getUntrackedParameter<bool>("is2022")),
    sampleIs2022EE(                                                              iConfig.getUntrackedParameter<bool>("is2022EE")),
    sampleIs2023(                                                                iConfig.getUntrackedParameter<bool>("is2023")),
    sampleIs2023BPix(                                                            iConfig.getUntrackedParameter<bool>("is2023BPix")),
    sampleIs2024(                                                                iConfig.getUntrackedParameter<bool>("is2024")),
    sampleIs2025(                                                                iConfig.getUntrackedParameter<bool>("is2025")),
    storeLheParticles(                                                           iConfig.getUntrackedParameter<bool>("storeLheParticles")),
    storeGenParticles(                                                           iConfig.getUntrackedParameter<bool>("storeGenParticles")),
    storeParticleLevel(                                                          iConfig.getUntrackedParameter<bool>("storeParticleLevel")),
    storeJecSourcesAll(                                                          iConfig.getUntrackedParameter<bool>("storeJecSourcesAll")),
    storeJecSourcesGrouped(                                                      iConfig.getUntrackedParameter<bool>("storeJecSourcesGrouped")),
    storeAllTauID(                                                               iConfig.getUntrackedParameter<bool>("storeAllTauID")),
    //storePrefireComponents(                                                      iConfig.getUntrackedParameter<bool>("storePrefireComponents")),
    storeJetSubstructure(                                                        iConfig.getUntrackedParameter<bool>("storeJetSubstructure"))
{
    // In case of run-2 data apply ecal filter
    if( is2017() || is2018() ) ecalBadCalibFilterToken = consumes<bool>(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));
    if( is2022() || is2022EE() || is2023() || is2023BPix()() ) ecalBadCalibFilterToken = consumes<bool>(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));
    
    // Create new sub-analyzer objects
    triggerAnalyzer             = new TriggerAnalyzer(      iConfig, this);
    lheAnalyzer                 = new LheAnalyzer(          iConfig, this);
    genAnalyzer                 = new GenAnalyzer(          iConfig, this);
    particleLevelAnalyzer       = new ParticleLevelAnalyzer(iConfig, this);
    jetAnalyzer                 = new JetAnalyzer(          iConfig, this);
    leptonAnalyzer              = new LeptonAnalyzer(       iConfig, this);
    k0Analyzer                  = new K0Analyzer(           iConfig, this);
}

// destructor
V0Analyzer::~V0Analyzer() 
{
    // Delete sub-analyzer objects
    delete triggerAnalyzer;
    delete lheAnalyzer;
    delete genAnalyzer;
    delete particleLevelAnalyzer;
    delete jetAnalyzer;
    delete leptonAnalyzer;
    delete k0Analyzer;
}

//---------------------------------------------------------------------
// General Functions 
//---------------------------------------------------------------------

//  method called for each lumi block 
void V0Analyzer::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup)
{
    _lumiBlock = (unsigned long) iLumi.id().luminosityBlock();
}

//  method fills 'descriptions' with the allowed parameters for the module
void V0Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

}

//---------------------------------------------------------------------
//  method called for each run
//---------------------------------------------------------------------
void V0Analyzer::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
    _runNb                      = (unsigned long) iRun.id().run();
    triggerAnalyzer->reIndex    = true;
}

//---------------------------------------------------------------------
//  method called for each Job:
//          -method called once each job just before starting event loop 
//---------------------------------------------------------------------

void V0Analyzer::beginJob() 
{
    //Initialize tree with event info
    outputTree          = fs->make<TTree>(   "blackJackAndHookersTree", "blackJackAndHookersTree");
    nVertices           = fs->make<TH1D>(    "nVertices",               "Number of vertices", 120, 0, 120);

    //Set general branches of the outputTree
    outputTree->Branch( "_runNb",            &_runNb,                   "_runNb/l");
    outputTree->Branch( "_lumiBlock",        &_lumiBlock,               "_lumiBlock/l");
    outputTree->Branch( "_eventNb",          &_eventNb,                 "_eventNb/l");
    outputTree->Branch( "_nVertex",          &_nVertex,                 "_nVertex/i");

    // Run-2 branches
    if( isRun2())
    {
        outputTree->Branch("_is2017",        &sampleIs2017,             "_is2017/O");
        outputTree->Branch("_is2018",        &sampleIs2018,             "_is2018/O");
        outputTree->Branch("_is2016preVFP",  &sampleIs2016preVFP,       "_is2016preVFP/O");
    }
    else {  // Run-3 branches
        outputTree->Branch("_is2022",        &sampleIs2022,             "_is2022/O");
        outputTree->Branch("_is2022EE",      &sampleIs2022EE,           "_is2022EE/O");
        outputTree->Branch("_is2023",        &sampleIs2023,             "_is2023/O");
        outputTree->Branch("_is2023BPix",    &sampleIs2023BPix,         "_is2023BPix/O");
        outputTree->Branch("_is2024",        &sampleIs2024,             "_is2024/O");
        outputTree->Branch("_is2025",        &sampleIs2025,             "_is2025/O");
    }

    // Set Prefiring branches for MC samples
    /*
    if( isMC())
    {
        outputTree->Branch("_prefireWeight",                &_prefireWeight,        "_prefireWeight/F");
        outputTree->Branch("_prefireWeightUp",              &_prefireWeightUp,      "_prefireWeightUp/F");
        outputTree->Branch("_prefireWeightDown",            &_prefireWeightDown,    "_prefireWeightDown/F");

        if (storePrefireComponents)
        {
            outputTree->Branch("_prefireWeightMuon",        &_prefireWeightMuon,    "_prefireWeightMuon/F");
            outputTree->Branch("_prefireWeightMuonUp",      &_prefireWeightMuonUp,  "_prefireWeightMuonUp/F");
            outputTree->Branch("_prefireWeightMuonDown",    &_prefireWeightMuonDown,"_prefireWeightMuonDown/F");
            outputTree->Branch("_prefireWeightECAL",        &_prefireWeightECAL,    "_prefireWeightECAL/F");
            outputTree->Branch("_prefireWeightECALUp",      &_prefireWeightECALUp,  "_prefireWeightECALUp/F");
            outputTree->Branch("_prefireWeightECALDown",    &_prefireWeightECALDown,"_prefireWeightECALDown/F");
        }
    }*/

    //Initialize MC analyzers 
    if( isMC() )                        lheAnalyzer->beginJob(outputTree, fs);
    if( isMC() )                        genAnalyzer->beginJob(outputTree, fs);
    if( isMC() && storeParticleLevel)   particleLevelAnalyzer->beginJob(outputTree);

    //Initialize other analyzers 
    triggerAnalyzer->beginJob(          outputTree);
    leptonAnalyzer->beginJob(           outputTree);
    jetAnalyzer->beginJob(              outputTree);
    k0Analyzer->beginJob(               outputTree);

    _runNb = 0;

}

//---------------------------------------------------------------------
//  Analyzer:
//          -method called for each event  
//---------------------------------------------------------------------
void V0Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
    using namespace edm;
    auto                                                        vertices = getHandle(iEvent, vtxToken);

    if(isMC())                                                  lheAnalyzer->analyze(iEvent);

    _nVertex                                                    = vertices->size();
    nVertices->Fill(                                            _nVertex, lheAnalyzer->getWeight());

    bool                                                        applySkim;
    if(isMC() && storeParticleLevel)                            applySkim = !particleLevelAnalyzer->analyze(iEvent);
    else                                                        applySkim = true;

    if(_nVertex == 0)                                           return;
    if(!leptonAnalyzer->analyze(iEvent, *(vertices->begin())) 
            and applySkim)                                      return;
    if(!jetAnalyzer->analyze(iEvent) 
            and applySkim)                                      return;
    if( isMC() )                                                genAnalyzer->analyze(iEvent);
    triggerAnalyzer->analyze(                                   iEvent);
    k0Analyzer->analyze(                                        iEvent, iSetup, *(vertices->begin()));

    _eventNb                                                    = (unsigned long) iEvent.id().event();

    // Calculate Prefiring for MC samples
    /*
    if(isMC())
    {
        _prefireWeight                                          = *(getHandle(iEvent, prefireWeightToken));
        _prefireWeightUp                                        = *(getHandle(iEvent, prefireWeightUpToken));
        _prefireWeightDown                                      = *(getHandle(iEvent, prefireWeightDownToken));

        if (storePrefireComponents)
        {
            _prefireWeightMuon                                  = *(getHandle(iEvent, prefireWeightMuonToken));
            _prefireWeightMuonUp                                = *(getHandle(iEvent, prefireWeightMuonUpToken));
            _prefireWeightMuonDown                              = *(getHandle(iEvent, prefireWeightMuonDownToken));
            _prefireWeightECAL                                  = *(getHandle(iEvent, prefireWeightECALToken));
            _prefireWeightECALUp                                = *(getHandle(iEvent, prefireWeightECALUpToken));
            _prefireWeightECALDown                              = *(getHandle(iEvent, prefireWeightECALDownToken));

        }
    }
    */

    // Fill output TTree
    outputTree->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0Analyzer);
