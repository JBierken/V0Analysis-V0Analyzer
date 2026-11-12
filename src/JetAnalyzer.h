#ifndef JET_ANALYZER_H
#define JET_ANALYZER_H
//CMSSW
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

//c++ standard library
#include <memory>
#include <string>

//ROOT
#include "TTree.h"

//include other parts of framework
#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"

class V0Analyzer;

/*
 * AK4 PUPPI jets (slimmedJetsPuppi), the default jet collection in Run 3.
 * The CHS collection (slimmedJets) is no longer stored.
 *
 * Flavour-tagging scores follow the NanoAOD convention: -1 means "not available"
 * (tagger not present in this MiniAOD production, or jet outside its acceptance).
 * Availability per production (checked on real files):
 *   DeepJet             : all Run 3
 *   ParticleNet (PNet)  : all Run 3 (central model, |eta| < 2.5)
 *   RobustParT          : 2022, 2023, 2024 PromptReco (not in 2024 MINIv6 / 2025)
 *   UParT               : 2024 MINIv6 reprocessing, 2025 (not in 2022, 2023, 2024 PromptReco)
 */
class JetAnalyzer {

    //Declare friend classes
    friend class V0Analyzer;

    private:
        // Maximum number of jets
        static const unsigned nJets_max = 100;

        // Jet counter
        unsigned _nJets = 0;

        // Jet kinematics
        double   _jetPt[nJets_max];
        double   _jetEta[nJets_max];
        double   _jetPhi[nJets_max];
        double   _jetE[nJets_max];

        // Jet ID (Run 3 PUPPI definition)
        bool     _jetIsTight[nJets_max];
        bool     _jetIsTightLepVeto[nJets_max];

        // Jet flavour (MC truth, 0 in data). Parton flavour is signed.
        int      _jetHadronFlavor[nJets_max];
        int      _jetPartonFlavor[nJets_max];

        // Flavour tagging: DeepJet
        double   _jetDeepFlavor_b[nJets_max];
        double   _jetDeepFlavor_bb[nJets_max];
        double   _jetDeepFlavor_lepb[nJets_max];
        double   _jetDeepFlavor_c[nJets_max];
        double   _jetDeepFlavor_uds[nJets_max];
        double   _jetDeepFlavor_g[nJets_max];
        double   _jetDeepFlavor[nJets_max];             // b + bb + lepb (b vs all)

        // Flavour tagging: ParticleNet (from MiniAOD, PUPPI central)
        double   _jetParticleNet_b[nJets_max];
        double   _jetParticleNet_c[nJets_max];
        double   _jetParticleNet_uds[nJets_max];
        double   _jetParticleNet_g[nJets_max];
        double   _jetParticleNet_BvsAll[nJets_max];     // use this as the b-tag score
        double   _jetParticleNet_CvsL[nJets_max];
        double   _jetParticleNet_CvsB[nJets_max];
        double   _jetParticleNet_QvsG[nJets_max];       // central or forward model, depending on |eta|

        // Flavour tagging: RobustParT
        double   _jetRobustParT_BvsAll[nJets_max];
        double   _jetRobustParT_CvsL[nJets_max];
        double   _jetRobustParT_CvsB[nJets_max];
        double   _jetRobustParT_QvsG[nJets_max];

        // Flavour tagging: UParT
        double   _jetUParT_BvsAll[nJets_max];
        double   _jetUParT_CvsL[nJets_max];
        double   _jetUParT_CvsB[nJets_max];
        double   _jetUParT_QvsG[nJets_max];
        double   _jetUParT_SvsUDG[nJets_max];
        double   _jetUParT_SvsBC[nJets_max];

        // Gen-level jet kinematics (MC only)
        bool     _jetHasGen[nJets_max];
        double   _jetGenPt[nJets_max];
        double   _jetGenEta[nJets_max];
        double   _jetGenPhi[nJets_max];
        double   _jetGenE[nJets_max];

        // MET (type-1 PF) and PUPPI MET
        double   _met;
        double   _metPhi;
        double   _metPuppi;
        double   _metPuppiPhi;

        // Jet ID selectors
        std::unique_ptr<PFJetIDSelectionFunctor> jetIdTight;
        std::unique_ptr<PFJetIDSelectionFunctor> jetIdTightLepVeto;

        // Check the available tagger labels once per job
        bool taggersChecked = false;
        void checkTaggers(const pat::Jet& jet);

        V0Analyzer* myAnalyzer;

        bool passJetId(PFJetIDSelectionFunctor& selector, const pat::Jet& jet) const;

    public:
        JetAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* vars);
        ~JetAnalyzer();

        void beginJob(TTree* outputTree);
        bool analyze(const edm::Event&);
};

#endif
