/*
   Custom analyzer class for finding secondary vertices corresponding to light neutral hadron (V0) decays.
   */

#ifndef K0_ANALYZER_H
#define K0_ANALYZER_H

// include other parts of the V0Analyzer framework
#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"

// system include files
#include <memory>
#include <unordered_map>
#include <Math/SVector.h>       // root high-performance vector class
#include <Math/SMatrix.h>       // root high-performance matrix class
#include "TMatrixDSym.h"        // for fixTrackCovariance
#include "TVectorD.h"           // for fixTrackCovariance

// main include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

// general include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "TTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"

// vertex fitter include files
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// data format include files
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

class V0Analyzer;

class K0Analyzer {
    // Define main analyzer as friend
    friend class V0Analyzer;
    
    private:
        // constants
        static constexpr double pimass      = 0.13957;
        static constexpr double pimass2     = pimass*pimass;
        static constexpr double pmass       = 0.93827;
        static constexpr double pmass2      = pmass*pmass;
        static constexpr double emass       = 0.000511;
        static constexpr double emass2      = emass*emass;
        static constexpr double mumass      = 0.1057;
        static constexpr double mumass2     = mumass*mumass;
        static constexpr double ksmass      = 0.49761;
        static constexpr double lambdamass  = 1.1157;
        static constexpr double jpsimass    = 3.0969;

        // V0 mass windows (final selection, applied after the vertex fit)
        static constexpr double ksMassWindow        = 0.07;
        static constexpr double lambdaMassWindow    = 0.05;
        // Loose pre-fit mass filter: before the (expensive) Kalman fit, the pair is skipped
        // if the K0S, Lambda and LambdaBar masses computed from the track momenta at the PCA
        // are all outside preFitMassWindowScale times the final windows above.
        // A value <= 0 disables the pre-filter. Set from the optional untracked parameter
        // v0PreFitMassWindowScale (default 5: on 1000 2024G events 1 of 40230 V0s lost, ~25% less CPU than no filter).
        double preFitMassWindowScale = 5.;

        // Per-event cache of a preselected track: everything that only depends on
        // the single track is computed once instead of once per track pair.
        struct SelTrack {
            reco::Track             track;
            reco::TransientTrack    ttrack;             // built from the original track
            bool                    validIP     = false;// impactPointTSCP() is valid
            FreeTrajectoryState     ipState;            // state at the impact point
            bool                    hasFixed    = false;
            reco::TransientTrack    ttrackFixed;        // built from fixTrackCovariance(track), lazily
            double                  relIso      = -1.;  // lazily computed (<0: not yet)
        };

        // Output of the V0 fit for one track pair
        struct V0Result {
            double type = -1., invmass = 0.;
            double vtxx = 0., vtxy = 0., vtxz = 0., vtxxunc = 0., vtxyunc = 0., vtxzunc = 0.;
            double vtxr_pv = 0., vtxr_bs = 0., vtxrunc_pv = 0., vtxrunc_bs = 0., vtxrsig_pv = 0., vtxrsig_bs = 0.;
            double px = 0., py = 0., pz = 0., pt = 0., eta = 0., phi = 0.;
            double dca = 0., pcax = 0., pcay = 0., pcaz = 0.;
            double vtxnormchi2 = 0.;
            double pxpos = 0., pypos = 0., pzpos = 0., ptpos = 0.;
            double pxneg = 0., pyneg = 0., pzneg = 0., ptneg = 0.;
            double etapos = 0., etaneg = 0., phipos = 0., phineg = 0.;
            double nhitspos = 0., nhitsneg = 0., normchi2pos = 0., normchi2neg = 0.;
            double d0pos = 0., dzpos = 0., d0neg = 0., dzneg = 0.;
            double isopos = 0., isoneg = 0.;
        };

        // vertex fitter (smoothing on, so refitted tracks are available); stateless between fits
        KalmanVertexFitter v0Fitter{true};
        
        // Root tree variable declarations
        static const unsigned nV0s_max      = 50;
        
        // beamspot position
        double      _beamSpotX;
        double      _beamSpotY;
        double      _beamSpotZ;
        
        // primary vertex position and position uncertainty
        double      _primaryVertexX;
        double      _primaryVertexY;
        double      _primaryVertexZ;
        double      _primaryVertexXUnc;
        double      _primaryVertexYUnc;
        double      _primaryVertexZUnc;
        
        // variables from customized V0 fitter
        unsigned    _nV0s = 0;    
        double      _V0InvMass[nV0s_max];
        unsigned    _V0Type[nV0s_max];
        double      _V0X[nV0s_max];
        double      _V0Y[nV0s_max];
        double      _V0Z[nV0s_max];
        double      _V0XUnc[nV0s_max];
        double      _V0YUnc[nV0s_max];
        double      _V0ZUnc[nV0s_max];
        double      _V0RPV[nV0s_max];
        double      _V0RBS[nV0s_max];
        double      _V0RPVUnc[nV0s_max];
        double      _V0RBSUnc[nV0s_max];
        double      _V0RPVSig[nV0s_max];
        double      _V0RBSSig[nV0s_max];
        double      _V0Px[nV0s_max];
        double      _V0Py[nV0s_max];
        double      _V0Pz[nV0s_max];
        double      _V0Pt[nV0s_max];
        double      _V0Eta[nV0s_max];
        double      _V0Phi[nV0s_max];
        double      _V0DCA[nV0s_max];
        double      _V0PCAX[nV0s_max];
        double      _V0PCAY[nV0s_max];
        double      _V0PCAZ[nV0s_max];
        double      _V0VtxNormChi2[nV0s_max];
        double      _V0PxPos[nV0s_max];
        double      _V0PyPos[nV0s_max];
        double      _V0PzPos[nV0s_max];
        double      _V0PtPos[nV0s_max];
        double      _V0PxNeg[nV0s_max];
        double      _V0PyNeg[nV0s_max];
        double      _V0PzNeg[nV0s_max];
        double      _V0PtNeg[nV0s_max];
        double      _V0EtaPos[nV0s_max];
        double      _V0EtaNeg[nV0s_max];
        double      _V0PhiPos[nV0s_max];
        double      _V0PhiNeg[nV0s_max];
        double      _V0NHitsPos[nV0s_max];
        double      _V0NHitsNeg[nV0s_max];
        double      _V0NormChi2Pos[nV0s_max];
        double      _V0NormChi2Neg[nV0s_max];
        double      _V0D0Pos[nV0s_max];
        double      _V0DzPos[nV0s_max];
        double      _V0D0Neg[nV0s_max];
        double      _V0DzNeg[nV0s_max];
        double      _V0IsoPos[nV0s_max];
        double      _V0IsoNeg[nV0s_max];

        unsigned debugcounter = 0;

        V0Analyzer* myAnalyzer;

    public:
        // Constructor/destructor
        K0Analyzer(const edm::ParameterSet& iConfig, V0Analyzer* vars);
        ~K0Analyzer();
        
        // template member functions
        void beginJob(TTree*);
        void analyze(const edm::Event&, const edm::EventSetup&, const reco::Vertex&);
        
        // help member functions
        bool passTrackPreselection(const reco::Track&, const reco::BeamSpot&, const MagneticField*) const;
        double getTrackRelIso(const reco::Track&,
                const std::vector<pat::PackedCandidate>&,
                const std::vector<pat::PackedCandidate>&) const;
        double getTrackRelIso(SelTrack&,
                const std::vector<pat::PackedCandidate>&,
                const std::vector<pat::PackedCandidate>&) const;
        bool VZeroFitter(SelTrack&, SelTrack&,
                const reco::BeamSpot&,
                const reco::Vertex&,
                const MagneticField*,
                const std::vector<pat::PackedCandidate>&,
                const std::vector<pat::PackedCandidate>&,
                V0Result&) const;
};

#endif
