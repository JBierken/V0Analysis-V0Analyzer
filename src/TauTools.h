#ifndef TauTools_H
#define TauTools_H

//include CMSSW classes
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

//include other parts of framework
#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"

namespace TauTools{

    void getNextDaughters(const reco::GenParticle& gen, const std::vector<reco::GenParticle>& genParticles, std::set<int>& list);

    //Check whether a tau decayed hadronically
    bool decayedHadronically(const reco::GenParticle&, const std::vector<reco::GenParticle>&);

    const reco::GenParticle* findMatch(const pat::Tau& reco, const std::vector<reco::GenParticle>& genParticles);

    const reco::GenJet* findMatchedGenJet(const reco::GenParticle& genTau, const std::vector<reco::GenJet>& genJets);

    const bool considerForMatching(const pat::Tau&, const reco::GenParticle&, const std::vector<reco::GenParticle>& genParticles);

    const unsigned tauGenStatus(const reco::GenParticle*);
}

#endif
