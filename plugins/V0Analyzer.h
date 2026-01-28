// -*- C++ -*-
//
// Package:    V0Analysis/V0Analyzer
// Class:      V0Analyzer
//
/**\class V0Analyzer V0Analyzer.cc V0Analysis/V0Analyzer/plugins/V0Analyzer.cc

Description: [one line class summary]

Implementation:
Implementation of the NTuplizer for the Displaced vertexing study using V0 candidates for run 3 operation 
*/
//
// Original Author: Jas Bierkens 
//         Created:  Tue, 07 May 2024 14:22:28 GMT
//
//
#ifndef V0_ANALYZER_H
#define V0_ANALYZER_H

// system include files
#include <memory>

/********************************************************************
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
 ********************************************************************/

//include CMSSW classes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//include ROOT classes
#include "TTree.h"

//include other parts of framework
#include "V0Analysis/V0Analyzer/src/TriggerAnalyzer.h"        
#include "V0Analysis/V0Analyzer/src/LheAnalyzer.h"
#include "V0Analysis/V0Analyzer/src/GenAnalyzer.h"
#include "V0Analysis/V0Analyzer/src/ParticleLevelAnalyzer.h"
#include "V0Analysis/V0Analyzer/src/JetAnalyzer.h"
#include "V0Analysis/V0Analyzer/src/LeptonAnalyzer.h"
#include "V0Analysis/V0Analyzer/src/K0Analyzer.h"

//---------------------------------------------------------------------
//  Allow for easy way to retrieve handles
//---------------------------------------------------------------------
namespace {
    template<typename T, typename I> edm::Handle<T> getHandle(const I& iEvent,const edm::EDGetTokenT<T>& token){
        edm::Handle<T> handle;
        iEvent.getByToken(token,handle);
        return handle;
    }
}
//---------------------------------------------------------------------
//  class declarations
//---------------------------------------------------------------------

// Sub-Analyzer declaration
class TriggerAnalyzer;
class LheAnalyzer;
class GenAnalyzer;
class ParticleLevelAnalyzer;
class JetAnalyzer;
class LeptonAnalyzer;
class K0Analyzer;

using reco::TrackCollection;

//Main Analyzer class
class V0Analyzer : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks, edm::one::WatchRuns, edm::one::SharedResources> 
{
    //Define other analyzers as friends
    friend TriggerAnalyzer;
    friend LheAnalyzer;
    friend GenAnalyzer;
    friend ParticleLevelAnalyzer;
    friend JetAnalyzer;
    friend LeptonAnalyzer;
    friend K0Analyzer;

    public:

        explicit V0Analyzer(const edm::ParameterSet&);                                // default constructor
        ~V0Analyzer() override;                                                       // destructor

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        //Helper functions
        bool isData()       const{ return sampleIsData; }
        bool isMC()         const{ return !sampleIsData; }
        bool isRun2()       const{ return !(sampleIs2022 || sampleIs2022EE || sampleIs2023 || sampleIs2023BPix || sampleIs2024); }
        bool isRun3()       const{ return (sampleIs2022 || sampleIs2022EE || sampleIs2023 || sampleIs2023BPix || sampleIs2024); }
        bool is2017()       const{ return sampleIs2017; }
        bool is2018()       const{ return sampleIs2018; }
        bool is2022()       const{ return sampleIs2022; }
        bool is2022EE()     const{ return sampleIs2022EE; }
        bool is2023()       const{ return sampleIs2023; }
        bool is2023BPix()   const{ return sampleIs2023BPix; }
        bool is2024()       const{ return sampleIs2024; }

    private:

        virtual void beginJob() override;
        virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
        virtual void beginRun(const edm::Run&, edm::EventSetup const&) override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        virtual void endRun(const edm::Run&, edm::EventSetup const&) override {}
        virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {}
        virtual void endJob() override {}

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::BeamSpot>                            beamSpotToken;
        edm::EDGetTokenT<std::vector<reco::Vertex>>                 vtxToken;
        edm::EDGetTokenT<GenEventInfoProduct>                       genEventInfoToken;
        edm::EDGetTokenT<LHEEventProduct>                           lheEventInfoToken;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>>            pileUpToken;
        edm::EDGetTokenT<reco::GenParticleCollection>               genParticleToken;
        edm::EDGetTokenT<reco::GenParticleCollection>               particleLevelPhotonsToken;
        edm::EDGetTokenT<reco::GenJetCollection>                    particleLevelLeptonsToken;
        edm::EDGetTokenT<reco::GenJetCollection>                    particleLevelJetsToken;
        edm::EDGetTokenT<reco::GenJetCollection>                    genJetsToken;
        edm::EDGetTokenT<reco::METCollection>                       particleLevelMetsToken;
        edm::EDGetTokenT<std::vector<pat::Muon>>                    muonToken;
        edm::EDGetTokenT<std::vector<pat::Electron>>                eleToken;
        edm::EDGetTokenT<std::vector<pat::Tau>>                     tauToken;
        edm::EDGetTokenT<std::vector<reco::GenJet>>                 tauGenJetsToken;
        edm::EDGetTokenT<std::vector<pat::Photon>>                  photonToken;
        edm::EDGetTokenT<std::vector<pat::PackedCandidate>>         packedCandidatesToken;
        edm::EDGetTokenT<std::vector<pat::PackedCandidate>>         lostTracksToken;
        //edm::EDGetTokenT<double>                                    rhoToken;
        edm::EDGetTokenT<std::vector<pat::MET>>                     metToken;
        edm::EDGetTokenT<std::vector<pat::MET>>                     metPuppiToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>                     jetToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>                     jetPuppiToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>                     jetSmearedToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>                     jetSmearedUpToken;
        edm::EDGetTokenT<std::vector<pat::Jet>>                     jetSmearedDownToken;
        edm::EDGetTokenT<edm::TriggerResults>                       recoResultsPrimaryToken;
        edm::EDGetTokenT<edm::TriggerResults>                       recoResultsSecondaryToken;
        edm::EDGetTokenT<edm::TriggerResults>                       triggerToken;
        edm::EDGetTokenT<pat::PackedTriggerPrescales>               prescalesToken;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>    trigObjToken;
        edm::EDGetTokenT<bool>                                      ecalBadCalibFilterToken;
        //edm::EDGetTokenT<double>                                    prefireWeightToken;
        //edm::EDGetTokenT<double>                                    prefireWeightUpToken;
        //edm::EDGetTokenT<double>                                    prefireWeightDownToken;
        //edm::EDGetTokenT<double>                                    prefireWeightMuonToken;
        //edm::EDGetTokenT<double>                                    prefireWeightMuonUpToken;
        //edm::EDGetTokenT<double>                                    prefireWeightMuonDownToken;
        //edm::EDGetTokenT<double>                                    prefireWeightECALToken;
        //edm::EDGetTokenT<double>                                    prefireWeightECALUpToken;
        //edm::EDGetTokenT<double>                                    prefireWeightECALDownToken;
        edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>    bFieldToken_;

        // Helper variables
        std::string                 skim;
        bool                        sampleIsData;
        bool                        sampleIs2016preVFP;
        bool                        sampleIs2017;
        bool                        sampleIs2018;
        bool                        sampleIs2022;
        bool                        sampleIs2022EE;
        bool                        sampleIs2023;
        bool                        sampleIs2023BPix;
        bool                        sampleIs2024;
        bool                        storeLheParticles;
        bool                        storeGenParticles;
        bool                        storeParticleLevel;
        bool                        storeJecSourcesAll;
        bool                        storeJecSourcesGrouped;
        bool                        storeAllTauID;
        //bool                        storePrefireComponents;
        bool                        storeJetSubstructure;

        // TTree variables
        edm::Service<TFileService>  fs;
        TTree*                      outputTree;
        TH1D*                       nVertices;

        // Run variables
        unsigned long               _runNb;
        unsigned long               _lumiBlock;
        unsigned long               _eventNb;
        unsigned                    _nVertex;

        //float                       _prefireWeight;
        //float                       _prefireWeightUp;
        //float                       _prefireWeightDown;
        //float                       _prefireWeightMuon;
        //float                       _prefireWeightMuonUp;
        //float                       _prefireWeightMuonDown;
        //float                       _prefireWeightECAL;
        //float                       _prefireWeightECALUp;
        //float                       _prefireWeightECALDown;
        
        // Declare sub-analyzers
        TriggerAnalyzer*            triggerAnalyzer;
        LheAnalyzer*                lheAnalyzer;
        GenAnalyzer*                genAnalyzer;
        ParticleLevelAnalyzer*      particleLevelAnalyzer;
        JetAnalyzer*                jetAnalyzer;
        LeptonAnalyzer*             leptonAnalyzer;
        K0Analyzer*                 k0Analyzer;

};
#endif
