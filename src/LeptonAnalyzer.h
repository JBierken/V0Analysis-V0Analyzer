#ifndef LEPTON_ANALYZER_H
#define LEPTON_ANALYZER_H

//include CMSSW classes
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

//include ROOT classes
#include "TTree.h"

//include c++ library classes
#include <memory>                                              //for std::shared_ptr

//include other parts of the framework
#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"
#include "V0Analysis/V0Analyzer/src/GenTools.h"
#include "V0Analysis/V0Analyzer/src/TauTools.h"

class V0Analyzer;

class LeptonAnalyzer {
  //Friend classes and functions
  friend class V0Analyzer;
  
  private:
    //this has to come before the effective areas as their initialization depends on it!
    V0Analyzer* myAnalyzer;

    //maximum number of leptons to be stored 
    static const unsigned nL_max    = 20; 

    //number of leptons of each type in the event
    unsigned _nL                    = 0;
    unsigned _nMu                   = 0;
    unsigned _nEle                  = 0;
    unsigned _nLight                = 0;
    unsigned _nTau                  = 0;

    //lepton kinematics and systematic variations
    double      _lPt[nL_max];
    double      _lEta[nL_max];
    double      _lPhi[nL_max];
    double      _lE[nL_max];
    
    //lepton flavor and charge 
    unsigned    _lFlavor[nL_max];
    int         _lCharge[nL_max];

    //pointing variables
    double      _dxy[nL_max];
    double      _dz[nL_max];
    double      _3dIP[nL_max];
    double      _3dIPSig[nL_max];
    double      _tauDxyLead[nL_max];
    double      _tauDzLead[nL_max];
    
    //muon properties
    unsigned    _tauGenStatus[nL_max];                 //1: prompt ele, 2:prompt mu, 3: ele from leptonic tau, 4:mu from leptonic tau, 5: hadronically decayed tau, 6:rest 
    
    //official POG selection definitions
    bool        _lPOGVeto[nL_max];
    bool        _lPOGLoose[nL_max];
    bool        _lPOGMedium[nL_max];
    bool        _lPOGTight[nL_max];

    //MC truth information from matching 
    bool        _lIsPrompt[nL_max];
    int         _lMatchPdgId[nL_max];
    int         _lMatchCharge[nL_max];
    double      _lMatchPt[nL_max];
    bool        _lHasMatch[nL_max];
    int         _lMomPdgId[nL_max];

    template <typename Lepton> void fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles);
    void fillTauGenVars(const pat::Tau&, const std::vector<reco::GenParticle>& genParticles);
    void fillLeptonKinVars(const reco::Candidate&);
    void fillLeptonImpactParameters(const pat::Electron& );
    void fillLeptonImpactParameters(const pat::Muon& );
    void fillLeptonImpactParameters(const pat::Tau&, const reco::Vertex&);
    double tau_dz(const pat::Tau&, const reco::Vertex::Point&) const;
    bool eleMuOverlap(const pat::Electron& ele, const bool* loose) const;
    bool tauLightOverlap(const pat::Tau& tau, const bool* loose) const;
    //void fillLeptonJetVariables(const reco::Candidate&, edm::Handle<std::vector<pat::Jet>>&, const reco::Vertex&, const double rho, const bool oldMatching = false);

    // In LeptonAnalyzerId.cc
    bool  passTriggerEmulationDoubleEG(const pat::Electron*, const bool hOverE = true) const;               //For ewkino id it needs to be possible to check hOverE separately
    float slidingCut(float, float, float) const;
    bool  passingElectronMvaHZZ(const pat::Electron*, double) const;
    bool  passingElectronMvaLooseSusy(const pat::Electron*, double, double) const;
    bool  passingElectronMvaMediumSusy(const pat::Electron*, double) const;
    bool  passingElectronMvaTightSusy(const pat::Electron*, double) const;

  public:
    LeptonAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* vars);
    ~LeptonAnalyzer();

    void beginJob(TTree* outputTree);
    bool analyze(const edm::Event&, const reco::Vertex&);
};

#endif
