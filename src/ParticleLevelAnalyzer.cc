//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//include ROOT classes
#include "TLorentzVector.h"

//include other parts of code
#include "V0Analysis/V0Analyzer/src/ParticleLevelAnalyzer.h"

/*
   Class storing data for unfolding to particle-level in differential measurement
   Saving the products created by https://twiki.cern.ch/twiki/bin/viewauth/CMS/ParticleLevelProducer
*/

ParticleLevelAnalyzer::ParticleLevelAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* myAnalyzer):
    myAnalyzer(myAnalyzer){};

void ParticleLevelAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("_pl_met",                   &_pl_met,                   "_pl_met/D");
    outputTree->Branch("_pl_metPhi",                &_pl_metPhi,                "_pl_metPhi/D");
    outputTree->Branch("_pl_nPh",                   &_pl_nPh,                   "_pl_nPh/i");
    outputTree->Branch("_pl_phPt",                  &_pl_phPt,                  "_pl_phPt[_pl_nPh]/D");
    outputTree->Branch("_pl_phEta",                 &_pl_phEta,                 "_pl_phEta[_pl_nPh]/D");
    outputTree->Branch("_pl_phPhi",                 &_pl_phPhi,                 "_pl_phPhi[_pl_nPh]/D");
    outputTree->Branch("_pl_phE",                   &_pl_phE,                   "_pl_phE[_pl_nPh]/D");
    outputTree->Branch("_pl_nL",                    &_pl_nL,                    "_pl_nL/i");
    outputTree->Branch("_pl_lPt",                   &_pl_lPt,                   "_pl_lPt[_pl_nL]/D");
    outputTree->Branch("_pl_lEta",                  &_pl_lEta,                  "_pl_lEta[_pl_nL]/D");
    outputTree->Branch("_pl_lPhi",                  &_pl_lPhi,                  "_pl_lPhi[_pl_nL]/D");
    outputTree->Branch("_pl_lE",                    &_pl_lE,                    "_pl_lE[_pl_nL]/D");
    outputTree->Branch("_pl_lFlavor",               &_pl_lFlavor,               "_pl_lFlavor[_pl_nL]/i");
    outputTree->Branch("_pl_lCharge",               &_pl_lCharge,               "_pl_lCharge[_pl_nL]/I");
    outputTree->Branch("_pl_nJets",                 &_pl_nJets,                 "_pl_nJets/i");
    outputTree->Branch("_pl_jetPt",                 &_pl_jetPt,                 "_pl_jetPt[_pl_nJets]/D");
    outputTree->Branch("_pl_jetEta",                &_pl_jetEta,                "_pl_jetEta[_pl_nJets]/D");
    outputTree->Branch("_pl_jetPhi",                &_pl_jetPhi,                "_pl_jetPhi[_pl_nJets]/D");
    outputTree->Branch("_pl_jetE",                  &_pl_jetE,                  "_pl_jetE[_pl_nJets]/D");
    outputTree->Branch("_pl_jetHadronFlavor",       &_pl_jetHadronFlavor,       "_pl_jetHadronFlavor[_pl_nJets]/i");

}

bool ParticleLevelAnalyzer::analyze(const edm::Event& iEvent){
    edm::Handle<std::vector<reco::GenParticle>> photons = getHandle(iEvent, myAnalyzer->particleLevelPhotonsToken);
    edm::Handle<std::vector<reco::GenJet>> leptons      = getHandle(iEvent, myAnalyzer->particleLevelLeptonsToken);
    edm::Handle<std::vector<reco::GenJet>> jets         = getHandle(iEvent, myAnalyzer->particleLevelJetsToken);
    edm::Handle<std::vector<reco::MET>> mets            = getHandle(iEvent, myAnalyzer->particleLevelMetsToken);

    _pl_met    = mets->front().pt();
    _pl_metPhi = mets->front().phi();

    _pl_nPh = 0;
    for(const reco::GenParticle& p : *photons){
        if(_pl_nPh == pl_nPh_max) break;
        _pl_phPt[_pl_nPh]  = p.pt();
        _pl_phEta[_pl_nPh] = p.eta();
        _pl_phPhi[_pl_nPh] = p.phi();
        _pl_phE[_pl_nPh]   = p.energy();
        ++_pl_nPh;
    }

    _pl_nL = 0;
    for(const reco::GenJet& p : *leptons){
        if(_pl_nL == pl_nL_max) break;
        _pl_lPt[_pl_nL]     = p.pt();
        _pl_lEta[_pl_nL]    = p.eta();
        _pl_lPhi[_pl_nL]    = p.phi();
        _pl_lE[_pl_nL]      = p.energy();
        _pl_lCharge[_pl_nL] = p.charge();

        if(abs(p.pdgId()) == 11)      _pl_lFlavor[_pl_nL] = 0;
        else if(abs(p.pdgId()) == 13) _pl_lFlavor[_pl_nL] = 1;
        else                          _pl_lFlavor[_pl_nL] = 2;
        ++_pl_nL;
    }

    _pl_nJets = 0;
    for(const reco::GenJet& p : *jets){
        if(_pl_nJets == pl_nJet_max) break;
        _pl_jetPt[_pl_nJets]           = p.pt();
        _pl_jetEta[_pl_nJets]          = p.eta();
        _pl_jetPhi[_pl_nJets]          = p.phi();
        _pl_jetE[_pl_nJets]            = p.energy();
        _pl_jetHadronFlavor[_pl_nJets] = p.pdgId();
        ++_pl_nJets;
    }

    if(myAnalyzer->skim == "trilep"       and _pl_nL < 3)                    return false;
    if(myAnalyzer->skim == "dilep"        and _pl_nL < 3)                    return false;
    if(myAnalyzer->skim == "singlelep"    and _pl_nL < 1)                    return false;
    if(myAnalyzer->skim == "FR"           and (_pl_nL < 1 or _pl_nJets < 1)) return false;
    if(myAnalyzer->skim == "singlephoton" and _pl_nPh < 1)                   return false;
    if(myAnalyzer->skim == "diphoton"     and _pl_nPh < 2)                   return false;
    if(myAnalyzer->skim == "singlejet"    and _pl_nJets < 1)                 return false;
    return true;
}
