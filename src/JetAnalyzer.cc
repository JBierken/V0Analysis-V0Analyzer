//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

//include c++ library classes
#include <cmath>
#include <set>
#include <sstream>

//include other parts of framework
#include "V0Analysis/V0Analyzer/src/JetAnalyzer.h"

namespace {
    // Tagger labels (AK4 PUPPI, MiniAOD)
    const std::string kDeepJet      = "pfDeepFlavourJetTags:";
    const std::string kPNet         = "pfParticleNetFromMiniAODAK4PuppiCentralJetTags:";
    const std::string kPNetDisc     = "pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:";
    const std::string kPNetFwdDisc  = "pfParticleNetFromMiniAODAK4PuppiForwardDiscriminatorsJetTags:";
    const std::string kRobustParT   = "pfParticleTransformerAK4DiscriminatorsJetTags:";
    const std::string kUParT        = "pfUnifiedParticleTransformerAK4DiscriminatorsJetTags:";

    // Returns the tagger score, or -1 if the label is absent (bDiscriminator gives -1000)
    // or the tagger did not evaluate the jet (negative value).
    double tagScore(const pat::Jet& jet, const std::string& label){
        const double value = jet.bDiscriminator(label);
        return (std::isfinite(value) && value >= 0.) ? value : -1.;
    }
}

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* myAnalyzer):
    myAnalyzer(myAnalyzer)
{
    // Jet ID for PUPPI jets. For Run 3 JME recommends the "2022FG/2023CD" definition for
    // the 22Sep2023 ReReco of 2022 and for 2023; it is also used here for 2024, 2025 and 2026.
    // (The "runsBCDEprompt" version only applies to the 2022 PromptReco, which is not used.)
    const PFJetIDSelectionFunctor::Version_t version = myAnalyzer->isRun3()
        ? PFJetIDSelectionFunctor::RUN3PUPPIruns2022FGruns2023CD
        : PFJetIDSelectionFunctor::RUN2ULPUPPI;

    jetIdTight        = std::make_unique<PFJetIDSelectionFunctor>(version, PFJetIDSelectionFunctor::TIGHT);
    jetIdTightLepVeto = std::make_unique<PFJetIDSelectionFunctor>(version, PFJetIDSelectionFunctor::TIGHTLEPVETO);
}

JetAnalyzer::~JetAnalyzer() = default;

void JetAnalyzer::beginJob(TTree* outputTree)
{
    outputTree->Branch("_nJets",                    &_nJets,                    "_nJets/i");
    // Jet kinematics
    outputTree->Branch("_jetPt",                    &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetEta",                   &_jetEta,                   "_jetEta[_nJets]/D");
    outputTree->Branch("_jetPhi",                   &_jetPhi,                   "_jetPhi[_nJets]/D");
    outputTree->Branch("_jetE",                     &_jetE,                     "_jetE[_nJets]/D");
    // Jet ID
    outputTree->Branch("_jetIsTight",               &_jetIsTight,               "_jetIsTight[_nJets]/O");
    outputTree->Branch("_jetIsTightLepVeto",        &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");
    // Jet flavour
    outputTree->Branch("_jetHadronFlavor",          &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/I");
    outputTree->Branch("_jetPartonFlavor",          &_jetPartonFlavor,          "_jetPartonFlavor[_nJets]/I");
    // Flavour tagging: DeepJet
    outputTree->Branch("_jetDeepFlavor_b",          &_jetDeepFlavor_b,          "_jetDeepFlavor_b[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_bb",         &_jetDeepFlavor_bb,         "_jetDeepFlavor_bb[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_lepb",       &_jetDeepFlavor_lepb,       "_jetDeepFlavor_lepb[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_c",          &_jetDeepFlavor_c,          "_jetDeepFlavor_c[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_uds",        &_jetDeepFlavor_uds,        "_jetDeepFlavor_uds[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_g",          &_jetDeepFlavor_g,          "_jetDeepFlavor_g[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor",            &_jetDeepFlavor,            "_jetDeepFlavor[_nJets]/D");
    // Flavour tagging: ParticleNet
    outputTree->Branch("_jetParticleNet_b",         &_jetParticleNet_b,         "_jetParticleNet_b[_nJets]/D");
    outputTree->Branch("_jetParticleNet_c",         &_jetParticleNet_c,         "_jetParticleNet_c[_nJets]/D");
    outputTree->Branch("_jetParticleNet_uds",       &_jetParticleNet_uds,       "_jetParticleNet_uds[_nJets]/D");
    outputTree->Branch("_jetParticleNet_g",         &_jetParticleNet_g,         "_jetParticleNet_g[_nJets]/D");
    outputTree->Branch("_jetParticleNet_BvsAll",    &_jetParticleNet_BvsAll,    "_jetParticleNet_BvsAll[_nJets]/D");
    outputTree->Branch("_jetParticleNet_CvsL",      &_jetParticleNet_CvsL,      "_jetParticleNet_CvsL[_nJets]/D");
    outputTree->Branch("_jetParticleNet_CvsB",      &_jetParticleNet_CvsB,      "_jetParticleNet_CvsB[_nJets]/D");
    outputTree->Branch("_jetParticleNet_QvsG",      &_jetParticleNet_QvsG,      "_jetParticleNet_QvsG[_nJets]/D");
    // Flavour tagging: RobustParT
    outputTree->Branch("_jetRobustParT_BvsAll",     &_jetRobustParT_BvsAll,     "_jetRobustParT_BvsAll[_nJets]/D");
    outputTree->Branch("_jetRobustParT_CvsL",       &_jetRobustParT_CvsL,       "_jetRobustParT_CvsL[_nJets]/D");
    outputTree->Branch("_jetRobustParT_CvsB",       &_jetRobustParT_CvsB,       "_jetRobustParT_CvsB[_nJets]/D");
    outputTree->Branch("_jetRobustParT_QvsG",       &_jetRobustParT_QvsG,       "_jetRobustParT_QvsG[_nJets]/D");
    // Flavour tagging: UParT
    outputTree->Branch("_jetUParT_BvsAll",          &_jetUParT_BvsAll,          "_jetUParT_BvsAll[_nJets]/D");
    outputTree->Branch("_jetUParT_CvsL",            &_jetUParT_CvsL,            "_jetUParT_CvsL[_nJets]/D");
    outputTree->Branch("_jetUParT_CvsB",            &_jetUParT_CvsB,            "_jetUParT_CvsB[_nJets]/D");
    outputTree->Branch("_jetUParT_QvsG",            &_jetUParT_QvsG,            "_jetUParT_QvsG[_nJets]/D");
    outputTree->Branch("_jetUParT_SvsUDG",          &_jetUParT_SvsUDG,          "_jetUParT_SvsUDG[_nJets]/D");
    outputTree->Branch("_jetUParT_SvsBC",           &_jetUParT_SvsBC,           "_jetUParT_SvsBC[_nJets]/D");
    // Gen-level jet kinematics
    if( ! myAnalyzer->isData() ) {
        outputTree->Branch("_jetHasGen",            &_jetHasGen,                "_jetHasGen[_nJets]/O");
        outputTree->Branch("_jetGenPt",             &_jetGenPt,                 "_jetGenPt[_nJets]/D");
        outputTree->Branch("_jetGenEta",            &_jetGenEta,                "_jetGenEta[_nJets]/D");
        outputTree->Branch("_jetGenPhi",            &_jetGenPhi,                "_jetGenPhi[_nJets]/D");
        outputTree->Branch("_jetGenE",              &_jetGenE,                  "_jetGenE[_nJets]/D");
    }
    // MET
    outputTree->Branch("_met",                      &_met,                      "_met/D");
    outputTree->Branch("_metPhi",                   &_metPhi,                   "_metPhi/D");
    outputTree->Branch("_metPuppi",                 &_metPuppi,                 "_metPuppi/D");
    outputTree->Branch("_metPuppiPhi",              &_metPuppiPhi,              "_metPuppiPhi/D");
}

// Report once which tagger families are stored in this production.
// DeepJet and ParticleNet are present in every Run 3 MiniAOD: if they are missing the labels are wrong.
void JetAnalyzer::checkTaggers(const pat::Jet& jet)
{
    std::set<std::string> labels;
    for(const auto& pair : jet.getPairDiscri()) labels.insert(pair.first);

    const std::vector<std::pair<std::string, std::string>> families = {
        {"DeepJet",     kDeepJet     + "probb"},
        {"ParticleNet", kPNetDisc    + "BvsAll"},
        {"RobustParT",  kRobustParT  + "BvsAll"},
        {"UParT",       kUParT       + "BvsAll"},
    };

    std::ostringstream summary;
    for(const auto& [name, label] : families){
        const bool found = labels.count(label) > 0;
        summary << " " << name << (found ? "=yes" : "=no");
        if( ! found && myAnalyzer->isRun3() && (name == "DeepJet" || name == "ParticleNet") ){
            throw cms::Exception("JetAnalyzer") << "Tagger label '" << label << "' not found on slimmedJetsPuppi; "
                                                << "the tagger labels do not match this MiniAOD production.";
        }
    }
    edm::LogInfo("JetAnalyzer") << "Flavour taggers available on PUPPI jets:" << summary.str()
                                << " (missing taggers are stored as -1)";
    taggersChecked = true;
}

bool JetAnalyzer::passJetId(PFJetIDSelectionFunctor& selector, const pat::Jet& jet) const
{
    pat::strbitset ret = selector.getBitTemplate();
    return selector(jet, ret);
}

bool JetAnalyzer::analyze(const edm::Event& iEvent)
{
    edm::Handle<std::vector<pat::Jet>> jets      = getHandle(iEvent, myAnalyzer->jetPuppiToken);
    edm::Handle<std::vector<pat::MET>> mets      = getHandle(iEvent, myAnalyzer->metToken);
    edm::Handle<std::vector<pat::MET>> metsPuppi = getHandle(iEvent, myAnalyzer->metPuppiToken);

    _nJets = 0;

    for(const auto& jet : *jets){
        if(_nJets == nJets_max) break;
        if( ! taggersChecked ) checkTaggers(jet);

        // Jet kinematics
        _jetPt[_nJets]                  = jet.pt();
        _jetEta[_nJets]                 = jet.eta();
        _jetPhi[_nJets]                 = jet.phi();
        _jetE[_nJets]                   = jet.energy();

        // Jet ID
        _jetIsTight[_nJets]             = passJetId(*jetIdTight,        jet);
        _jetIsTightLepVeto[_nJets]      = passJetId(*jetIdTightLepVeto, jet);

        // Jet flavour (0 in data)
        _jetHadronFlavor[_nJets]        = jet.hadronFlavour();
        _jetPartonFlavor[_nJets]        = jet.partonFlavour();

        // Flavour tagging: DeepJet
        _jetDeepFlavor_b[_nJets]        = tagScore(jet, kDeepJet + "probb");
        _jetDeepFlavor_bb[_nJets]       = tagScore(jet, kDeepJet + "probbb");
        _jetDeepFlavor_lepb[_nJets]     = tagScore(jet, kDeepJet + "problepb");
        _jetDeepFlavor_c[_nJets]        = tagScore(jet, kDeepJet + "probc");
        _jetDeepFlavor_uds[_nJets]      = tagScore(jet, kDeepJet + "probuds");
        _jetDeepFlavor_g[_nJets]        = tagScore(jet, kDeepJet + "probg");
        const bool deepJetValid = _jetDeepFlavor_b[_nJets] >= 0 && _jetDeepFlavor_bb[_nJets] >= 0 && _jetDeepFlavor_lepb[_nJets] >= 0;
        _jetDeepFlavor[_nJets]          = deepJetValid ? _jetDeepFlavor_b[_nJets] + _jetDeepFlavor_bb[_nJets] + _jetDeepFlavor_lepb[_nJets] : -1.;

        // Flavour tagging: ParticleNet
        // The central model only runs for |eta| < 2.5 (on its own eta definition, so jets at |eta| ~ 2.4-2.5 can
        // also be skipped). For skipped jets MiniAOD stores the raw scores as 0 and BvsAll as -1:
        // use BvsAll to decide whether the tagger ran, and set the raw scores to -1 if it did not.
        _jetParticleNet_BvsAll[_nJets]  = tagScore(jet, kPNetDisc + "BvsAll");
        const bool pnetEvaluated        = _jetParticleNet_BvsAll[_nJets] >= 0.;
        const bool pnetCentral          = std::abs(jet.eta()) < 2.5;
        _jetParticleNet_b[_nJets]       = pnetEvaluated ? tagScore(jet, kPNet + "probb")   : -1.;
        _jetParticleNet_c[_nJets]       = pnetEvaluated ? tagScore(jet, kPNet + "probc")   : -1.;
        _jetParticleNet_uds[_nJets]     = pnetEvaluated ? tagScore(jet, kPNet + "probuds") : -1.;
        _jetParticleNet_g[_nJets]       = pnetEvaluated ? tagScore(jet, kPNet + "probg")   : -1.;
        _jetParticleNet_CvsL[_nJets]    = tagScore(jet, kPNetDisc + "CvsL");
        _jetParticleNet_CvsB[_nJets]    = tagScore(jet, kPNetDisc + "CvsB");
        _jetParticleNet_QvsG[_nJets]    = tagScore(jet, (pnetCentral ? kPNetDisc : kPNetFwdDisc) + "QvsG");

        // Flavour tagging: RobustParT
        _jetRobustParT_BvsAll[_nJets]   = tagScore(jet, kRobustParT + "BvsAll");
        _jetRobustParT_CvsL[_nJets]     = tagScore(jet, kRobustParT + "CvsL");
        _jetRobustParT_CvsB[_nJets]     = tagScore(jet, kRobustParT + "CvsB");
        _jetRobustParT_QvsG[_nJets]     = tagScore(jet, kRobustParT + "QvsG");

        // Flavour tagging: UParT
        _jetUParT_BvsAll[_nJets]        = tagScore(jet, kUParT + "BvsAll");
        _jetUParT_CvsL[_nJets]          = tagScore(jet, kUParT + "CvsL");
        _jetUParT_CvsB[_nJets]          = tagScore(jet, kUParT + "CvsB");
        _jetUParT_QvsG[_nJets]          = tagScore(jet, kUParT + "QvsG");
        _jetUParT_SvsUDG[_nJets]        = tagScore(jet, kUParT + "SvsUDG");
        _jetUParT_SvsBC[_nJets]         = tagScore(jet, kUParT + "SvsBC");

        // Gen-level jet kinematics
        _jetHasGen[_nJets]  = false;
        _jetGenPt[_nJets]   = 0;
        _jetGenEta[_nJets]  = 0;
        _jetGenPhi[_nJets]  = 0;
        _jetGenE[_nJets]    = 0;
        if( ! myAnalyzer->isData() && jet.genJet() != nullptr ){
            _jetHasGen[_nJets]  = true;
            _jetGenPt[_nJets]   = jet.genJet()->pt();
            _jetGenEta[_nJets]  = jet.genJet()->eta();
            _jetGenPhi[_nJets]  = jet.genJet()->phi();
            _jetGenE[_nJets]    = jet.genJet()->energy();
        }

        ++_nJets;
    }

    // MET
    const pat::MET& met      = mets->front();
    _met                     = met.pt();
    _metPhi                  = met.phi();

    const pat::MET& metPuppi = metsPuppi->front();
    _metPuppi                = metPuppi.pt();
    _metPuppiPhi             = metPuppi.phi();

    // Legacy skims
    if(myAnalyzer->skim == "singlejet" and _nJets < 1) return false;
    if(myAnalyzer->skim == "FR" and _nJets < 1)        return false;
    if(myAnalyzer->skim == "fourTopBase")  {
        int nGoodJets   = 0;
        int nBJets      = 0;
        for (unsigned jet = 0; jet < _nJets; jet++) {
            if (_jetPt[jet] < 25)                       continue;
            if (std::abs(_jetEta[jet]) > 2.4)           continue;
            if (! _jetIsTight[jet])                     continue;
            nGoodJets++;

            if (_jetDeepFlavor[jet] < 0.045)            continue;
            nBJets++;
        }
        if (nGoodJets < 2)                              return false;
        if (nBJets < 1)                                 return false;
    }

    return true;
}
