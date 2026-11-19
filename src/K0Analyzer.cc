/***************************************************************************
 *
 *  Displaces vertexing Efficiency study Analyzer:
 *          - Custom analyzer class for finding secondary vertices corresponding to Light Neutral Hadron (V0) decays.
 *
 **************************************************************************/

// include header
// note: all other includes are defined in the header
#include "V0Analysis/V0Analyzer/src/K0Analyzer.h"

// constructor //
K0Analyzer::K0Analyzer(const edm::ParameterSet& iConfig, V0Analyzer* myAnalyzer):
    preFitMassWindowScale(iConfig.getUntrackedParameter<double>("v0PreFitMassWindowScale", 5.)),
    myAnalyzer(myAnalyzer){};

// destructor //
K0Analyzer::~K0Analyzer(){
}

// beginJob //
void K0Analyzer::beginJob(TTree* outputTree){
    // initialize branches in the output tree

    // beamspot position
    outputTree->Branch("_beamSpotX",            &_beamSpotX,            "_beamSpotX/D");
    outputTree->Branch("_beamSpotY",            &_beamSpotY,            "_beamSpotY/D");
    outputTree->Branch("_beamSpotZ",            &_beamSpotZ,            "_beamSpotZ/D");

    // primary vertex position and position uncertainty
    outputTree->Branch("_primaryVertexX",       &_primaryVertexX,       "_primaryVertexX/D");
    outputTree->Branch("_primaryVertexY",       &_primaryVertexY,       "_primaryVertexY/D");
    outputTree->Branch("_primaryVertexZ",       &_primaryVertexZ,       "_primaryVertexZ/D");
    outputTree->Branch("_primaryVertexXUnc",    &_primaryVertexXUnc,    "_primaryVertexXUnc/D");
    outputTree->Branch("_primaryVertexYUnc",    &_primaryVertexYUnc,    "_primaryVertexYUnc/D");
    outputTree->Branch("_primaryVertexZUnc",    &_primaryVertexZUnc,    "_primaryVertexZUnc/D");
    
    // V0 compound variables
    // - ID variables
    outputTree->Branch("_nV0s",                 &_nV0s,                 "_nV0s/i");
    outputTree->Branch("_V0InvMass",            &_V0InvMass,            "_V0InvMass[_nV0s]/D");
    outputTree->Branch("_V0Type",               &_V0Type,               "_V0Type[_nV0s]/i");
    // - vertex position and position uncertainty
    outputTree->Branch("_V0X",                  &_V0X,                  "_V0X[_nV0s]/D");
    outputTree->Branch("_V0Y",                  &_V0Y,                  "_V0Y[_nV0s]/D");
    outputTree->Branch("_V0Z",                  &_V0Z,                  "_V0Z[_nV0s]/D");
    outputTree->Branch("_V0XUnc",               &_V0XUnc,               "_V0XUnc[_nV0s]/D");
    outputTree->Branch("_V0YUnc",               &_V0YUnc,               "_V0YUnc[_nV0s]/D");
    outputTree->Branch("_V0ZUnc",               &_V0ZUnc,               "_V0ZUnc[_nV0s]/D");
    // - radial distance, uncertainty and significance
    outputTree->Branch("_V0RPV",                &_V0RPV,                "_V0RPV[_nV0s]/D");
    outputTree->Branch("_V0RBS",                &_V0RBS,                "_V0RBS[_nV0s]/D");
    outputTree->Branch("_V0RPVUnc",             &_V0RPVUnc,             "_V0RPVUnc[_nV0s]/D");
    outputTree->Branch("_V0RBSUnc",             &_V0RBSUnc,             "_V0RBSUnc[_nV0s]/D");
    outputTree->Branch("_V0RPVSig",             &_V0RPVSig,             "_V0RPVSig[_nV0s]/D");
    outputTree->Branch("_V0RBSSig",             &_V0RBSSig,             "_V0RBSSig[_nV0s]/D");
    // - momentum and direction
    outputTree->Branch("_V0Px",                 &_V0Px,                 "_V0Px[_nV0s]/D");
    outputTree->Branch("_V0Py",                 &_V0Py,                 "_V0Py[_nV0s]/D");
    outputTree->Branch("_V0Pz",                 &_V0Pz,                 "_V0Pz[_nV0s]/D");
    outputTree->Branch("_V0Pt",                 &_V0Pt,                 "_V0Pt[_nV0s]/D");
    outputTree->Branch("_V0Eta",                &_V0Eta,                "_V0Eta[_nV0s]/D");
    outputTree->Branch("_V0Phi",                &_V0Phi,                "_V0Phi[_nV0s]/D");
    // - point of closest approach
    outputTree->Branch("_V0DCA",                &_V0DCA,                "_V0DCA[_nV0s]/D");
    outputTree->Branch("_V0PCAX",               &_V0PCAX,               "_V0PCAX[_nV0s]/D");
    outputTree->Branch("_V0PCAY",               &_V0PCAY,               "_V0PCAY[_nV0s]/D");
    outputTree->Branch("_V0PCAZ",               &_V0PCAZ,               "_V0PCAZ[_nV0s]/D");
    // - other
    outputTree->Branch("_V0VtxNormChi2",        &_V0VtxNormChi2,        "_V0VtxNormChi2[_nV0s]/D");
    
    // V0 decay track variables
    // - momentum and direction
    outputTree->Branch("_V0PxPos",              &_V0PxPos,              "_V0PxPos[_nV0s]/D");
    outputTree->Branch("_V0PyPos",              &_V0PyPos,              "_V0PyPos[_nV0s]/D");
    outputTree->Branch("_V0PzPos",              &_V0PzPos,              "_V0PzPos[_nV0s]/D");
    outputTree->Branch("_V0PtPos",              &_V0PtPos,              "_V0PtPos[_nV0s]/D");
    outputTree->Branch("_V0PxNeg",              &_V0PxNeg,              "_V0PxNeg[_nV0s]/D");
    outputTree->Branch("_V0PyNeg",              &_V0PyNeg,              "_V0PyNeg[_nV0s]/D");
    outputTree->Branch("_V0PzNeg",              &_V0PzNeg,              "_V0PzNeg[_nV0s]/D");
    outputTree->Branch("_V0PtNeg",              &_V0PtNeg,              "_V0PtNeg[_nV0s]/D");
    outputTree->Branch("_V0EtaPos",             &_V0EtaPos,             "_V0EtaPos[_nV0s]/D");
    outputTree->Branch("_V0EtaNeg",             &_V0EtaNeg,             "_V0EtaNeg[_nV0s]/D");
    outputTree->Branch("_V0PhiPos",             &_V0PhiPos,             "_V0PhiPos[_nV0s]/D");
    outputTree->Branch("_V0PhiNeg",             &_V0PhiNeg,             "_V0PhiNeg[_nV0s]/D");
    // - track quality variables
    outputTree->Branch("_V0NHitsPos",           &_V0NHitsPos,           "_V0NHitsPos[_nV0s]/D");
    outputTree->Branch("_V0NHitsNeg",           &_V0NHitsNeg,           "_V0NHitsNeg[_nV0s]/D");
    outputTree->Branch("_V0NormChi2Pos",        &_V0NormChi2Pos,        "_V0NormChi2Pos[_nV0s]/D");
    outputTree->Branch("_V0NormChi2Neg",        &_V0NormChi2Neg,        "_V0NormChi2Neg[_nV0s]/D");
    outputTree->Branch("_V0D0Pos",              &_V0D0Pos,              "_V0D0Pos[_nV0s]/D");
    outputTree->Branch("_V0DzPos",              &_V0DzPos,              "_V0DzPos[_nV0s]/D");
    outputTree->Branch("_V0D0Neg",              &_V0D0Neg,              "_V0D0Neg[_nV0s]/D");
    outputTree->Branch("_V0DzNeg",              &_V0DzNeg,              "_V0DzNeg[_nV0s]/D");
    outputTree->Branch("_V0IsoPos",             &_V0IsoPos,             "_V0IsoPos[_nV0s]/D");
    outputTree->Branch("_V0IsoNeg",             &_V0IsoNeg,             "_V0IsoNeg[_nV0s]/D");
}

// analyze (main method) //
void K0Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::Vertex& primaryVertex){
    //get packedCandidates
    edm::Handle<std::vector<pat::PackedCandidate>> packedCandidates;
    iEvent.getByToken(myAnalyzer->packedCandidatesToken, packedCandidates);
    //get lostTracks 
    edm::Handle<std::vector<pat::PackedCandidate>> lostTracks;
    iEvent.getByToken(myAnalyzer->lostTracksToken, lostTracks);
    //get beamSpot
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(myAnalyzer->beamSpotToken, beamSpotHandle);
    const MagneticField* bfield = &iSetup.getData(myAnalyzer->bFieldToken_);

    // store primary vertex position and position uncertainty
    _primaryVertexX     = primaryVertex.position().x();
    _primaryVertexY     = primaryVertex.position().y();
    _primaryVertexZ     = primaryVertex.position().z();
    _primaryVertexXUnc  = primaryVertex.xError();
    _primaryVertexYUnc  = primaryVertex.yError();
    _primaryVertexZUnc  = primaryVertex.zError();
    
    // store beamspot position
    _beamSpotX          = beamSpotHandle->position().x();
    _beamSpotY          = beamSpotHandle->position().y();
    _beamSpotZ          = beamSpotHandle->position().z();

    // preselect PackedCandidates and lost tracks (see help function below),
    // and cache all single-track quantities needed in the pair loop
    std::vector<SelTrack> seltracks;
    seltracks.reserve(256);
    auto addTrack = [&](const pat::PackedCandidate& pc){
        if(!pc.hasTrackDetails()) return;
        const reco::Track& tr = *pc.bestTrack();
        if(!passTrackPreselection(tr, *beamSpotHandle, bfield)) return;
        SelTrack st;
        st.track    = tr;
        st.ttrack   = reco::TransientTrack(st.track, bfield);
        // if this throws, every pair with this track would throw (and be skipped) in the fit
        try{
            st.validIP = st.ttrack.impactPointTSCP().isValid();
            if(st.validIP) st.ipState = st.ttrack.impactPointTSCP().theState();
        } catch(...){ st.validIP = false; }
        seltracks.push_back(std::move(st));
    };
    for(const pat::PackedCandidate& pc : *packedCandidates) addTrack(pc);
    for(const pat::PackedCandidate& pc : *lostTracks)       addTrack(pc);

    // initialize number of V0s
    _nV0s = 0;

    // loop over pairs of tracks
    V0Result res;
    for(unsigned i=0; i < seltracks.size(); i++){
        for(unsigned j=i+1; j < seltracks.size(); j++){
            // temporary try catch for 2018D files giving a strange 
            // "BasicSingleVertexState::could not invert weight matrix" error
            // that presumably has something to do with the kalman fitter
            res = V0Result();
            bool isV0 = false;
            try{
                // (see help function below)
                isV0 = VZeroFitter(seltracks[i], seltracks[j], *beamSpotHandle, primaryVertex, bfield,
                                   *packedCandidates, *lostTracks, res);
            } catch(...){ continue; }

            // if no valid type was found, continue without writing a new V0 candidate
            if(!isV0) continue;

            // write V0 properties
            _V0InvMass[_nV0s]       = res.invmass;
            _V0Type[_nV0s]          = res.type;
            _V0X[_nV0s]             = res.vtxx;
            _V0Y[_nV0s]             = res.vtxy;
            _V0Z[_nV0s]             = res.vtxz;
            _V0XUnc[_nV0s]          = res.vtxxunc;
            _V0YUnc[_nV0s]          = res.vtxyunc;
            _V0ZUnc[_nV0s]          = res.vtxzunc;
            _V0RPV[_nV0s]           = res.vtxr_pv;
            _V0RBS[_nV0s]           = res.vtxr_bs;
            _V0RPVUnc[_nV0s]        = res.vtxrunc_pv;
            _V0RBSUnc[_nV0s]        = res.vtxrunc_bs;
            _V0RPVSig[_nV0s]        = res.vtxrsig_pv;
            _V0RBSSig[_nV0s]        = res.vtxrsig_bs;
            _V0Px[_nV0s]            = res.px;
            _V0Py[_nV0s]            = res.py;
            _V0Pz[_nV0s]            = res.pz;
            _V0Pt[_nV0s]            = res.pt;
            _V0Eta[_nV0s]           = res.eta;
            _V0Phi[_nV0s]           = res.phi;
            _V0DCA[_nV0s]           = res.dca;
            _V0PCAX[_nV0s]          = res.pcax;
            _V0PCAY[_nV0s]          = res.pcay;
            _V0PCAZ[_nV0s]          = res.pcaz;
            _V0VtxNormChi2[_nV0s]   = res.vtxnormchi2;
            _V0PxPos[_nV0s]         = res.pxpos;
            _V0PyPos[_nV0s]         = res.pypos;
            _V0PzPos[_nV0s]         = res.pzpos;
            _V0PtPos[_nV0s]         = res.ptpos;
            _V0PxNeg[_nV0s]         = res.pxneg;
            _V0PyNeg[_nV0s]         = res.pyneg;
            _V0PzNeg[_nV0s]         = res.pzneg;
            _V0PtNeg[_nV0s]         = res.ptneg;
            _V0EtaPos[_nV0s]        = res.etapos;
            _V0EtaNeg[_nV0s]        = res.etaneg;
            _V0PhiPos[_nV0s]        = res.phipos;
            _V0PhiNeg[_nV0s]        = res.phineg;
            _V0NHitsPos[_nV0s]      = res.nhitspos;
            _V0NHitsNeg[_nV0s]      = res.nhitsneg;
            _V0NormChi2Pos[_nV0s]   = res.normchi2pos;
            _V0NormChi2Neg[_nV0s]   = res.normchi2neg;
            _V0D0Pos[_nV0s]         = res.d0pos;
            _V0DzPos[_nV0s]         = res.dzpos;
            _V0D0Neg[_nV0s]         = res.d0neg;
            _V0DzNeg[_nV0s]         = res.dzneg;
            _V0IsoPos[_nV0s]        = res.isopos;
            _V0IsoNeg[_nV0s]        = res.isoneg;

            ++_nV0s; 
            if(_nV0s == nV0s_max) break;
        }
        if(_nV0s == nV0s_max) break;
    }
}


bool K0Analyzer::passTrackPreselection(
        const reco::Track& tr,
        const reco::BeamSpot& beamSpot,
        const MagneticField* bfield) const{
    // minimal selection of an individual track:
    // loose quality and transverse impact parameter significance (w.r.t. beamline) >= 2
    if(!tr.quality(reco::TrackBase::qualityByName("loose"))) return false;
    // special for impact parameter: use proper extrapolation to beamspot
    FreeTrajectoryState initialFTS = trajectoryStateTransform::initialFreeState(tr, bfield);
    TSCBLBuilderNoMaterial blsBuilder;
    TrajectoryStateClosestToBeamLine tscb(blsBuilder(initialFTS, beamSpot));
    if(!tscb.isValid()) return false;
    if(tscb.transverseImpactParameter().significance() < 2.) return false;
    return true;
}


double K0Analyzer::getTrackRelIso(const reco::Track& tr,
        const std::vector<pat::PackedCandidate>& packedCandidates,
        const std::vector<pat::PackedCandidate>& lostTracks) const{
    // calculate relative isolation of a track
    // (still to check if this works correctly or if there isn't a better way)
    // (however, not yet used in further analysis, so harmless for now)
    double relIso = 0;
    for(const pat::PackedCandidate& pc: packedCandidates){
        double dR = reco::deltaR(tr,pc);
        if((dR < 0.3) and (dR > 1e-10)){
            relIso = relIso + pc.pt();
        }
    }
    for(const pat::PackedCandidate& pc : lostTracks){
        double dR = reco::deltaR(tr,pc);
        if((dR < 0.3) and (dR > 1e-10)){
            relIso = relIso + pc.pt();
        }
    }
    relIso = relIso/tr.pt();
    return relIso;
}


double K0Analyzer::getTrackRelIso(SelTrack& st,
        const std::vector<pat::PackedCandidate>& packedCandidates,
        const std::vector<pat::PackedCandidate>& lostTracks) const{
    // cached per track: a track can enter several V0 candidates
    if(st.relIso < 0.) st.relIso = getTrackRelIso(st.track, packedCandidates, lostTracks);
    return st.relIso;
}


reco::Track fixTrackCovariance(const reco::Track& tk, double delta=1e-6){
    // see here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingPOGRecommendations
    unsigned int i, j;
    // Initialize minimum eigenvalue to a default value.
    double min_eig = 1;
    // Get the original covariance matrix.
    reco::TrackBase::CovarianceMatrix cov = tk.covariance();
    // Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues.
    TMatrixDSym new_cov(cov.kRows);
    for (i = 0; i < cov.kRows; i++) {
        for (j = 0; j < cov.kRows; j++) {
            // Need to check for nan or inf, because for some reason these cause a segfault when calling Eigenvectors()
            if (std::isnan(cov(i,j)) || std::isinf(cov(i,j))) cov(i,j) = 1e-6;
            // In all other cases, just copy the values over from cov to new_cov.
            new_cov(i,j) = cov(i,j);
        }
    }
    // Get the eigenvalues.
    TVectorD eig(cov.kRows);
    new_cov.EigenVectors(eig);
    // Find the minimum eigenvalue
    for (i = 0; i < cov.kRows; i++){
        if (eig(i) < min_eig){
            min_eig = eig(i);
        }
    }
    // If the minimum eigenvalue is less than zero, then subtract it from the diagonal and add `delta`.
    if (min_eig < 0) {
        for (i = 0; i < cov.kRows; i++){
            cov(i,i) -= min_eig - delta;
        }
    } 
    return reco::Track(tk.chi2(), tk.ndof(), tk.referencePoint(), 
            tk.momentum(), tk.charge(), cov, tk.algo(), 
            (reco::TrackBase::TrackQuality) tk.qualityMask()); 
}


bool K0Analyzer::VZeroFitter(SelTrack& st1,
        SelTrack& st2,
        const reco::BeamSpot& beamSpot,
        const reco::Vertex& primaryVertex,
        const MagneticField* bfield,
        const std::vector<pat::PackedCandidate>& packedCandidates,
        const std::vector<pat::PackedCandidate>& lostTracks,
        V0Result& out) const{
    // FIT AND SELECT V0 VERTICES FROM TWO TRACKS AND CALCULATE PROPERTIES
    // returns true if the pair is a V0 candidate (out.type = 1, 2 or 3)
    const reco::Track& tr1 = st1.track;
    const reco::Track& tr2 = st2.track;

    // candidates must have opposite charge
    // (pairs with a neutral track can never be assigned a positive and negative track, see below)
    if(tr1.charge()*tr2.charge() >= 0) return false;

    // the following fragment applies some conditions on the point of closest approach,
    // namely that the distance of closest approach must be reasonably small
    // and that the point of closest approach must be situated not too far from the center of CMS.
    if(!st1.validIP or !st2.validIP)        return false;
    ClosestApproachInRPhi capp; capp.calculate(st1.ipState, st2.ipState);
    if(!capp.status())                      return false;
    double dca = fabs(capp.distance());
    out.dca = dca;
    GlobalPoint cxpt = capp.crossingPoint();
    if(dca < 0. or dca > 1.) return false;
    if(std::sqrt(cxpt.x()*cxpt.x() + cxpt.y()*cxpt.y()) > 120.
            or std::abs(cxpt.z()) > 300.)   return false;
    out.pcax = cxpt.x();
    out.pcay = cxpt.y();
    out.pcaz = cxpt.z();

    // the following fragment calculates the trajectory states at PCA
    TrajectoryStateClosestToPoint tscp1 = 
        st1.ttrack.trajectoryStateClosestToPoint( cxpt );
    TrajectoryStateClosestToPoint tscp2 =
        st2.ttrack.trajectoryStateClosestToPoint( cxpt );
    if(!tscp1.isValid() or !tscp2.isValid()) return false;

    // loose pre-fit mass filter, using the track momenta at the PCA:
    // skip the expensive vertex fit if the pair is far from all V0 mass hypotheses.
    // (the final mass selection below is done with the refitted momenta at the vertex)
    if(preFitMassWindowScale > 0.){
        const GlobalVector pPosPCA  = (tr1.charge() > 0) ? tscp1.momentum() : tscp2.momentum();
        const GlobalVector pNegPCA  = (tr1.charge() > 0) ? tscp2.momentum() : tscp1.momentum();
        const double pPos2          = pPosPCA.mag2();
        const double pNeg2          = pNegPCA.mag2();
        const double pTot2          = (pPosPCA + pNegPCA).mag2();
        auto pairMass = [&](double m2pos, double m2neg){
            double e = std::sqrt(pPos2 + m2pos) + std::sqrt(pNeg2 + m2neg);
            return std::sqrt(std::max(0., e*e - pTot2));
        };
        const bool nearV0Mass = 
            std::abs(pairMass(pimass2, pimass2) - ksmass)       < preFitMassWindowScale*ksMassWindow     or
            std::abs(pairMass(pmass2,  pimass2) - lambdamass)   < preFitMassWindowScale*lambdaMassWindow or
            std::abs(pairMass(pimass2, pmass2)  - lambdamass)   < preFitMassWindowScale*lambdaMassWindow;
        if(!nearV0Mass) return false;
    }

    // preliminary conditions are met; now fit a vertex and perform further selections
    // note: need to apply the fix for track covariance matrices recommended here:
    //       https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingPOGRecommendations
    //       (done once per track and cached)
    for(SelTrack* st : {&st1, &st2}){
        if(!st->hasFixed){
            st->ttrackFixed = reco::TransientTrack(fixTrackCovariance(st->track), bfield);
            st->hasFixed    = true;
        }
    }
    std::vector<reco::TransientTrack> transtracks{st1.ttrackFixed, st2.ttrackFixed};
    TransientVertex v0vtx = v0Fitter.vertex(transtracks);
    // condition: vertex must be valid
    if(!v0vtx.isValid())                    return false;
    // condition: chi squared of fit must be small
    if(v0vtx.normalisedChiSquared() > 15.)  return false;
    if(v0vtx.normalisedChiSquared() < 0.)   return false;
    out.vtxnormchi2 = v0vtx.normalisedChiSquared();
    // calculate refitted tracks for further selection
    std::vector<reco::TransientTrack> refittedTracks;
    if(v0vtx.hasRefittedTracks()) refittedTracks = v0vtx.refittedTracks();
    else                                    return false; // should not happen when using smoothing in the fitter

    // convert TransientVertex to RecoVertex for later use
    reco::Vertex v0recovtx = v0vtx;
    // store vertex position
    out.vtxx       = v0recovtx.position().x();
    out.vtxy       = v0recovtx.position().y();
    out.vtxz       = v0recovtx.position().z();
    out.vtxxunc    = v0recovtx.xError();
    out.vtxyunc    = v0recovtx.yError();
    out.vtxzunc    = v0recovtx.zError();

    // convert vertex, beamspot and primary vertex position to GlobalPoints for later use
    GlobalPoint vtxXYZ(v0vtx.position().x(), v0vtx.position().y(), v0vtx.position().z());
    GlobalPoint beamSpotXYZ(beamSpot.position().x(),
            beamSpot.position().y(),
            beamSpot.position().z());
    GlobalPoint primaryVertexXYZ(primaryVertex.position().x(),
            primaryVertex.position().y(),
            primaryVertex.position().z());

    // calculate transverse displacement (+ unc and sig) with respect to the beamspot
    ROOT::Math::SVector<double, 3>  vtxrvec_bs(vtxXYZ.x()-beamSpotXYZ.x(), vtxXYZ.y()-beamSpotXYZ.y(), 0.);
    double vtxr_bs                  = ROOT::Math::Mag(vtxrvec_bs);
    ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> totalcovmat_bs = 
        beamSpot.rotatedCovariance3D() + v0recovtx.covariance();
    double vtxrunc_bs               = sqrt(ROOT::Math::Similarity(totalcovmat_bs, vtxrvec_bs))/vtxr_bs;
    double vtxrsig_bs               = vtxr_bs/vtxrunc_bs;
    out.vtxr_bs                     = vtxr_bs;
    out.vtxrunc_bs                  = vtxrunc_bs;
    out.vtxrsig_bs                  = vtxrsig_bs;

    // calculate transverse displacement (+ unc and sig) with respect to primary vertex
    ROOT::Math::SVector<double, 3> vtxrvec_pv(vtxXYZ.x()-primaryVertexXYZ.x(), 
            vtxXYZ.y()-primaryVertexXYZ.y(), 0.);
    double vtxr_pv                  = ROOT::Math::Mag(vtxrvec_pv);
    ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> totalcovmat_pv = 
        primaryVertex.covariance() + v0recovtx.covariance();
    double vtxrunc_pv               = sqrt(ROOT::Math::Similarity(totalcovmat_pv, vtxrvec_pv))/vtxr_pv;
    double vtxrsig_pv               = vtxr_pv/vtxrunc_pv;
    out.vtxr_pv                     = vtxr_pv;
    out.vtxrunc_pv                  = vtxrunc_pv;
    out.vtxrsig_pv                  = vtxrsig_pv;

    // calculate properties of tracks at the vertex
    SelTrack* posSel;
    SelTrack* negSel;
    reco::TransientTrack posrefittedtranstrack;
    reco::TransientTrack negrefittedtranstrack;

    if(tr1.charge() > 0. and tr2.charge() < 0){
        posSel                  = &st1;
        negSel                  = &st2;

        posrefittedtranstrack   = refittedTracks.at(0);
        negrefittedtranstrack   = refittedTracks.at(1);
    } 
    else if(tr1.charge() < 0. and tr2.charge() > 0){
        posSel                  = &st2;
        negSel                  = &st1;

        posrefittedtranstrack   = refittedTracks.at(1);
        negrefittedtranstrack   = refittedTracks.at(0);
    } 
    else                                        return false; // should not normally happen but just for safety
    const reco::Track& postrack = posSel->track;
    const reco::Track& negtrack = negSel->track;
    
    TrajectoryStateClosestToPoint postscp = 
        posrefittedtranstrack.trajectoryStateClosestToPoint(vtxXYZ);
    TrajectoryStateClosestToPoint negtscp = 
        negrefittedtranstrack.trajectoryStateClosestToPoint(vtxXYZ);
    
    if(!postscp.isValid() or !negtscp.isValid()) return false;
    GlobalVector positiveP(postscp.momentum());
    out.pxpos      = positiveP.x();
    out.pypos      = positiveP.y();
    out.pzpos      = positiveP.z();
    out.ptpos      = positiveP.transverse();
    out.etapos     = positiveP.eta();
    out.phipos     = positiveP.phi();
    
    GlobalVector negativeP(negtscp.momentum());
    out.pxneg      = negativeP.x();
    out.pyneg      = negativeP.y();
    out.pzneg      = negativeP.z();
    out.ptneg      = negativeP.transverse();
    out.etaneg     = negativeP.eta();
    out.phineg     = negativeP.phi();
    
    GlobalVector totalP(positiveP + negativeP);
    out.px         = totalP.x();
    out.py         = totalP.y();
    out.pz         = totalP.z();
    out.pt         = totalP.transverse();
    out.eta        = totalP.eta();
    out.phi        = totalP.phi();

    // (optional) condition: cosine of pointing angle must be close to one (not default)
    // note: not applied since it can be easily recalculated and used in selections
    //       in downstream stages of the analysis.

    // calculate energy assuming Kshort, Lambda and LambdaBar
    double piPlusE          = sqrt(positiveP.mag2() + pimass2); // pi plus
    double piMinE           = sqrt(negativeP.mag2() + pimass2); // pi minus
    double pE               = sqrt(positiveP.mag2() + pmass2); // proton
    double pBarE            = sqrt(negativeP.mag2() + pmass2); // antiproton
    double KsETot           = piPlusE + piMinE;
    double LambdaETot       = pE + piMinE;
    double LambdaBarETot    = pBarE + piPlusE;

    // calculate corresponding mass
    math::XYZTLorentzVector KsP4(totalP.x(),totalP.y(),totalP.z(),KsETot);
    double KsInvMass        = KsP4.M();

    math::XYZTLorentzVector LambdaP4(totalP.x(),totalP.y(),totalP.z(),LambdaETot);
    double LambdaInvMass    = LambdaP4.M();

    math::XYZTLorentzVector LambdaBarP4(totalP.x(),totalP.y(),totalP.z(),LambdaBarETot);
    double LambdaBarInvMass = LambdaBarP4.M();

    // mass must be close to Ks mass (default 0.07) or Lambda mass (default 0.05)
    if(std::abs(KsInvMass-ksmass) > ksMassWindow 
            and std::abs(LambdaInvMass-lambdamass) > lambdaMassWindow
            and std::abs(LambdaBarInvMass-lambdamass) > lambdaMassWindow)   return false;

    // fill track properties
    out.nhitspos       = postrack.numberOfValidHits();
    out.nhitsneg       = negtrack.numberOfValidHits();
    out.normchi2pos    = postrack.normalizedChi2();
    out.normchi2neg    = negtrack.normalizedChi2();
    reco::TrackBase::Point refPoint(primaryVertex.position().x(),
            primaryVertex.position().y(),
            primaryVertex.position().z());
    out.d0pos          = postrack.dxy(refPoint);
    out.dzpos          = postrack.dz(refPoint);
    out.d0neg          = negtrack.dxy(refPoint);
    out.dzneg          = negtrack.dz(refPoint);
    out.isopos         = getTrackRelIso(*posSel, packedCandidates, lostTracks);
    out.isoneg         = getTrackRelIso(*negSel, packedCandidates, lostTracks);

    // case 1: K0S
    if(std::abs(KsInvMass-ksmass) < ksMassWindow 
            and std::abs(KsInvMass-ksmass) < std::abs(LambdaInvMass-lambdamass)
            and std::abs(KsInvMass-ksmass) < std::abs(LambdaBarInvMass-lambdamass))
    {
        out.invmass    = KsInvMass;
        out.type       = 1.;
        return true;
    }
    // case 2a: mass close to Lambda mass (default: 0.05)
    if(std::abs(LambdaInvMass-lambdamass) < lambdaMassWindow
            and std::abs(LambdaInvMass-lambdamass) < std::abs(KsInvMass-ksmass)
            and std::abs(LambdaInvMass-lambdamass) < std::abs(LambdaBarInvMass-lambdamass))
    {
        out.invmass    = LambdaInvMass;
        out.type       = 2.;
        return true;
    }
    // case 2b: same but candidate is a LambdaBar
    if(std::abs(LambdaBarInvMass-lambdamass) < lambdaMassWindow
            and std::abs(LambdaBarInvMass-lambdamass) < std::abs(KsInvMass-ksmass)
            and std::abs(LambdaBarInvMass-lambdamass) < std::abs(LambdaInvMass-lambdamass))
    {
        out.invmass    = LambdaBarInvMass;
        out.type       = 3.;
        return true;
    }
    // if none of the above, this is not a V0 candidate
    return false;
}
