//include CMSSW classes
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

//include ROOT classes
#include "TLorentzVector.h"
#include "TRandom.h"

//include c++ library classes
#include <algorithm>

//include other parts of framework
#include "V0Analysis/V0Analyzer/src/LeptonAnalyzer.h"
#include "V0Analysis/V0Analyzer/src/GenTools.h"
#include "V0Analysis/V0Analyzer/src/TauTools.h"

LeptonAnalyzer::LeptonAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* myAnalyzer):
    myAnalyzer(myAnalyzer)
{
};

LeptonAnalyzer::~LeptonAnalyzer()
{
}

void LeptonAnalyzer::beginJob(TTree* outputTree)
{
    outputTree->Branch("_nL",                  &_nL,                "_nL/i");
    outputTree->Branch("_nMu",                 &_nMu,               "_nMu/i");
    outputTree->Branch("_nEle",                &_nEle,              "_nEle/i");
    outputTree->Branch("_nLight",              &_nLight,            "_nLight/i");
    outputTree->Branch("_nTau",                &_nTau,              "_nTau/i");
    outputTree->Branch("_lPt",                 &_lPt,               "_lPt[_nL]/D");
    outputTree->Branch("_lEta",                &_lEta,              "_lEta[_nL]/D");
    outputTree->Branch("_lPhi",                &_lPhi,              "_lPhi[_nL]/D");
    outputTree->Branch("_lE",                  &_lE,                "_lE[_nL]/D");
    outputTree->Branch("_lFlavor",             &_lFlavor,           "_lFlavor[_nL]/i");
    outputTree->Branch("_lCharge",             &_lCharge,           "_lCharge[_nL]/I");
    outputTree->Branch("_dxy",                 &_dxy,               "_dxy[_nL]/D");
    outputTree->Branch("_dz",                  &_dz,                "_dz[_nL]/D");
    outputTree->Branch("_3dIP",                &_3dIP,              "_3dIP[_nL]/D");
    outputTree->Branch("_3dIPSig",             &_3dIPSig,           "_3dIPSig[_nL]/D");
    outputTree->Branch("_lPOGVeto",            &_lPOGVeto,          "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",           &_lPOGLoose,         "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",          &_lPOGMedium,        "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",           &_lPOGTight,         "_lPOGTight[_nL]/O");
    outputTree->Branch("_tauDxyLead",          &_tauDxyLead,        "_tauDxyLead[_nL]/D");
    outputTree->Branch("_tauDzLead",           &_tauDzLead,         "_tauDzLead[_nL]/D");
    
    // initiate MC branches
    if( myAnalyzer->isMC() ){
        outputTree->Branch("_lIsPrompt",       &_lIsPrompt,         "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lMatchPdgId",     &_lMatchPdgId,       "_lMatchPdgId[_nL]/I");
        outputTree->Branch("_lMatchCharge",    &_lMatchCharge,      "_lMatchCharge[_nLight]/I");
        outputTree->Branch("_lMatchPt",        &_lMatchPt,          "_lMatchPt[_nLight]/D");
        outputTree->Branch("_lHasMatch",       &_lHasMatch,         "_lHasMatch[_nLight]/O");
        outputTree->Branch("_tauGenStatus",    &_tauGenStatus,      "_tauGenStatus[_nL]/i");
        outputTree->Branch("_lMomPdgId",       &_lMomPdgId,         "_lMomPdgId[_nL]/I");
    }
}

bool LeptonAnalyzer::analyze(const edm::Event& iEvent, const reco::Vertex& primaryVertex)
{
    edm::Handle<std::vector<pat::Electron>> electrons          = getHandle(iEvent, myAnalyzer->eleToken);
    edm::Handle<std::vector<pat::Muon>> muons                  = getHandle(iEvent, myAnalyzer->muonToken);
    edm::Handle<std::vector<pat::Tau>> taus                    = getHandle(iEvent, myAnalyzer->tauToken);
    //edm::Handle<double> rho                                    = getHandle(iEvent, myAnalyzer->rhoToken);
    //edm::Handle<std::vector<pat::Jet>> jets                    = getHandle(iEvent, myAnalyzer->jetToken);
    edm::Handle<std::vector<reco::GenParticle>> genParticles   = getHandle(iEvent, myAnalyzer->genParticleToken);

    // start lepton counters
    _nL     = 0;
    _nLight = 0;
    _nMu    = 0;
    _nEle   = 0;
    _nTau   = 0;

    // loop over muons
    // muons need to be run first, because some ID's need to calculate a muon veto for electrons
    for(const pat::Muon& mu : *muons){
        if(_nL == nL_max)                              break;
        if(mu.innerTrack().isNull())                   continue;
        if(mu.pt() < 5)                                continue;
        if(fabs(mu.eta()) > 2.4)                       continue;
        if(!mu.isPFMuon())                             continue;
        if(!(mu.isTrackerMuon() || mu.isGlobalMuon())) continue;
        
        //Fill kinematics parameters
        fillLeptonImpactParameters(mu);
        fillLeptonKinVars(mu);
        if( myAnalyzer->isMC() )            fillLeptonGenVars(mu, *genParticles);

        //store Lepton flavor and ID
        _lFlavor[_nL]      = 1;
        _lPOGVeto[_nL]     = mu.passed(reco::Muon::CutBasedIdLoose);                // no veto available, so we take loose here
        _lPOGLoose[_nL]    = mu.passed(reco::Muon::CutBasedIdLoose);
        _lPOGMedium[_nL]   = mu.passed(reco::Muon::CutBasedIdMedium);
        _lPOGTight[_nL]    = mu.passed(reco::Muon::CutBasedIdTight);

        ++_nMu;
        ++_nL;
        ++_nLight;
    }

    // Loop over electrons (note: using iterator we can easily get the ref too)
    for(auto ele = electrons->begin(); ele != electrons->end(); ++ele){
        if(_nL == nL_max)                               break;
        if(ele->gsfTrack().isNull())                    continue;
        if(ele->pt() < 7)                               continue;
        if(fabs(ele->eta()) > 2.5)                      continue;
        
        //Fill kinematics and impact parameters
        fillLeptonImpactParameters(*ele);
        fillLeptonKinVars(*ele);
        if( myAnalyzer->isMC() )            fillLeptonGenVars(*ele, *genParticles);

        //store Lepton flavor and ID
        _lFlavor[_nL]       = 0;
        //_lPOGVeto[_nL]      = ele->electronID("mvaEleID-RunIIIWinter22-V1-wp90");   // no veto available, so we take WP90 here
        //_lPOGLoose[_nL]     = ele->electronID("mvaEleID-RunIIIWinter22-V1-wp90");   // no loose available, so we take WP90 here
        //_lPOGMedium[_nL]    = ele->electronID("mvaEleID-RunIIIWinter22-V1-wp90");   // no medium available, the run-3 equivalent is WP90
        //_lPOGTight[_nL]     = ele->electronID("mvaEleID-RunIIIWinter22-V1-wp80");   // no tight available, the run-3 equivalent is WP80
        _lPOGVeto[_nL]      = ele->electronID("mvaEleID-Fall17-iso-V1-wp90");   // no veto available, so we take WP90 here
        _lPOGLoose[_nL]     = ele->electronID("mvaEleID-Fall17-iso-V1-wp90");   // no loose available, so we take WP90 here
        _lPOGMedium[_nL]    = ele->electronID("mvaEleID-Fall17-iso-V1-wp90");   // no medium available, the run-3 equivalent is WP90
        _lPOGTight[_nL]     = ele->electronID("mvaEleID-Fall17-iso-V1-wp80");   // no tight available, the run-3 equivalent is WP80

        ++_nEle;
        ++_nL;
        ++_nLight;
    }

    //loop over taus
    for(const pat::Tau& tau : *taus){
        if(_nL == nL_max)                       break;
        if(tau.pt() < 20)                       continue;                           // Minimum pt for tau reconstruction
        if(fabs(tau.eta()) > 2.3)               continue;
        
        //Fill kinematics parameters
        fillLeptonKinVars(tau);
        fillLeptonImpactParameters(tau, primaryVertex);
        if(myAnalyzer->isMC())              fillTauGenVars(tau, *genParticles);     //Still needs to be tested

        //store Lepton flavor and ID
        _lFlavor[_nL]                           = 2;
        _lPOGVeto[_nL]                          = ((tau.tauID("byUTagPUPPIDecayMode") > 0.5) 
                                                    && (tau.tauID("byUTagPUPPIVSjetraw") > 0.10) 
                                                    && (tau.tauID("byUTagPUPPIVSeraw")  > 0.50) 
                                                    && (tau.tauID("byUTagPUPPIVSmuraw") > 0.50)); 
        _lPOGLoose[_nL]                         = ((tau.tauID("byUTagPUPPIDecayMode") > 0.5) 
                                                    && (tau.tauID("byUTagPUPPIVSjetraw") > 0.10) 
                                                    && (tau.tauID("byUTagPUPPIVSeraw")  > 0.50) 
                                                    && (tau.tauID("byUTagPUPPIVSmuraw") > 0.50)); 
        _lPOGMedium[_nL]                        = ((tau.tauID("byUTagPUPPIDecayMode") > 0.5) 
                                                    && (tau.tauID("byUTagPUPPIVSjetraw") > 0.50) 
                                                    && (tau.tauID("byUTagPUPPIVSeraw")  > 0.50) 
                                                    && (tau.tauID("byUTagPUPPIVSmuraw") > 0.50)); 
        _lPOGTight[_nL]                         = ((tau.tauID("byUTagPUPPIDecayMode") > 0.5) 
                                                    && (tau.tauID("byUTagPUPPIVSjetraw") > 0.10) 
                                                    && (tau.tauID("byUTagPUPPIVSeraw")  > 0.50) 
                                                    && (tau.tauID("byUTagPUPPIVSmuraw") > 0.70)); 
        //_lPOGVeto[_nL]                          = tau.tauID("byVVVLooseDeepTau2017v2p1VSjet");
        //_lPOGLoose[_nL]                         = tau.tauID("byLooseDeepTau2017v2p1VSjet");
        //_lPOGMedium[_nL]                        = tau.tauID("byMediumDeepTau2017v2p1VSjet");
        //_lPOGTight[_nL]                         = tau.tauID("byTightDeepTau2017v2p1VSjet");
        
        ++_nTau;
        ++_nL;
    }

    //Initialize with default values for those tau-only arrays which weren't filled with electrons and muons [to allow correct comparison by the test script]
    for(auto array : {_tauDxyLead, _tauDzLead}) std::fill_n(array, _nLight, 0.);

    // Apply skim requirements for different analysis
    /*
    if(myAnalyzer->skim == "trilep"    &&  _nL     < 3)                         return false;
    if(myAnalyzer->skim == "dilep"     &&  _nL     < 2)                         return false;
    if(myAnalyzer->skim == "ttg"       &&  _nLight < 2)                         return false;
    if(myAnalyzer->skim == "singlelep" &&  _nL     < 1)                         return false;
    if(myAnalyzer->skim == "singletau" &&  _nTau   < 1)                         return false;
    if(myAnalyzer->skim == "FR"        &&  _nLight < 1)                         return false;
    if(myAnalyzer->skim == "lightdilep"   &&  _nLight < 2)                      return false;
    if(myAnalyzer->skim == "fourTopBase")  {
        // check for 2 tight leptons:
        if(_nLight < 2)                                                         return false;
        //if (_nLight == 2 && _lCharge[0] != _lCharge[1]) return false;
        int nMinimalLeptons = 0;
        int nTightLeptons   = 0;
        int sumCharge       = 0;

        for (unsigned lep=0; lep < _nL; lep++) {
            if (_lFlavor[lep] == 2)                                             continue;
            if (_lPt[lep] < 10 )                                                continue;
            if (fabs(_lEta[lep]) > 2.6 )                                        continue;
            //if (_leptonMvaTOPv2UL[lep] < 0.59 && _leptonMvaTOPUL[lep] < 0.20) continue;
            if (_3dIPSig[lep] >= 8)                                             continue;
            if (_dxy[lep] >= 0.05)                                              continue;
            if (_dz[lep] >= 0.1)                                                continue;
            if (_lFlavor[lep] == 1) {
                if( _lPOGMedium[lep])                                           continue;
            }
            nMinimalLeptons++;
            sumCharge += _lCharge[lep];

        }
        if (nMinimalLeptons-nTightLeptons > 1)                                  return false;
        if (nMinimalLeptons < 2)                                                return false;
        if (nMinimalLeptons == 2 && sumCharge == 0)                             return false;
    }*/
    return true;
}

void LeptonAnalyzer::fillLeptonKinVars(const reco::Candidate& lepton)
{
    _lPt[_nL]     = lepton.pt();
    _lEta[_nL]    = lepton.eta();
    _lPhi[_nL]    = lepton.phi();
    _lE[_nL]      = lepton.energy();
    _lCharge[_nL] = lepton.charge();
}

template <typename Lepton> void LeptonAnalyzer::fillLeptonGenVars(const Lepton& lepton, const std::vector<reco::GenParticle>& genParticles)
{
    const reco::GenParticle* match = lepton.genParticle();
    if(!match || match->pdgId() != lepton.pdgId()) match = GenTools::geometricMatch(lepton, genParticles); // if no match or pdgId is different, try the geometric match

    _tauGenStatus[_nL]          = TauTools::tauGenStatus(match); 
    _lIsPrompt[_nL]             = match && (abs(lepton.pdgId()) == abs(match->pdgId()) || match->pdgId() == 22) && GenTools::isPrompt(*match, genParticles); // only when matched to its own flavor or a photon
    _lMatchPdgId[_nL]           = match != nullptr ? match->pdgId() : 0;
    _lMatchCharge[_nL]          = match != nullptr ? match->charge() : 0;
    _lMatchPt[_nL]              = match ? match->pt() : 0.;
    _lHasMatch[_nL]             = ( match != nullptr );
    _lMomPdgId[_nL]             = match ? GenTools::getMotherPdgId(*match, genParticles) : 0;
}

void LeptonAnalyzer::fillTauGenVars(const pat::Tau& tau, const std::vector<reco::GenParticle>& genParticles)
{
    
    const reco::GenParticle* match  = TauTools::findMatch(tau, genParticles);

    
    _tauGenStatus[_nL]              = TauTools::tauGenStatus(match); 
    _lIsPrompt[_nL]                 = match && _tauGenStatus[_nL] != 6; 
    _lMatchPdgId[_nL]               = match ? match->pdgId() : 0;
    _lMomPdgId[_nL]                 = match ? GenTools::getMotherPdgId(*match, genParticles) : 0;

}

/*
 * Impact parameters:
 * Note: dB function seems to be preferred and more accurate over track->dxy and dz functions 
 * as the latter ones have some simplified extrapolation used (leading to slightly different values or opposite sign)
 * For taus: dxy is pre-computed with PV it was constructed with
 */
void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Electron& ele)
{
    _dxy[_nL]     = ele.dB(pat::Electron::PV2D);
    _dz[_nL]      = ele.dB(pat::Electron::PVDZ);
    _3dIP[_nL]    = ele.dB(pat::Electron::PV3D);
    _3dIPSig[_nL] = fabs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D));
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Muon& muon)
{
    _dxy[_nL]     = muon.dB(pat::Muon::PV2D);
    _dz[_nL]      = muon.dB(pat::Muon::PVDZ);
    _3dIP[_nL]    = muon.dB(pat::Muon::PV3D);
    _3dIPSig[_nL] = fabs(muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D));
}

void LeptonAnalyzer::fillLeptonImpactParameters(const pat::Tau& tau, const reco::Vertex& vertex)
{
   
    _dxy[_nL]     = (double) tau.dxy();                                      // warning: float while dxy of tracks are double; could also return -1000
   
   if( tau.leadChargedHadrCand().isNonnull() )
    {	
	    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
	    
	    _tauDxyLead[_nL]      = packedLeadTauCand->dxy();
	    _tauDzLead[_nL]       = packedLeadTauCand->dz();
	    _dz[_nL]              = tau_dz(tau, vertex.position());
    }
   
    _3dIP[_nL]    = tau.ip3d();
    _3dIPSig[_nL] = tau.ip3d_Sig();
}

//Function returning tau dz
double LeptonAnalyzer::tau_dz(const pat::Tau& tau, const reco::Vertex::Point& vertex) const
{
    const reco::Candidate::Point& tauVtx = tau.leadChargedHadrCand()->vertex();
    return (tauVtx.Z() - vertex.z()) - ((tauVtx.X() - vertex.x())*tau.px()+(tauVtx.Y()-vertex.y())*tau.py())/tau.pt()*tau.pz()/tau.pt();
}


//fucking disgusting lepton-jet matching based on obscure pointer defined somewhere in the abyss of CMSSW. (I feel ashamed for using this, but one has to be in sync with the POG - Willem )
template< typename T1, typename T2 > bool isSourceCandidatePtrMatch( const T1& lhs, const T2& rhs )
{
    for( size_t lhsIndex = 0; lhsIndex < lhs.numberOfSourceCandidatePtrs(); ++lhsIndex ){
        auto lhsSourcePtr = lhs.sourceCandidatePtr( lhsIndex );
        for( size_t rhsIndex = 0; rhsIndex < rhs.numberOfSourceCandidatePtrs(); ++rhsIndex ){
            auto rhsSourcePtr = rhs.sourceCandidatePtr( rhsIndex );
            if( lhsSourcePtr == rhsSourcePtr ){
                return true;
            }
        }
    }
    return false;
}


const pat::Jet* findMatchedJet( const reco::Candidate& lepton, const edm::Handle< std::vector< pat::Jet > >& jets, const bool oldMatching )
{
    
    //Look for jet that matches with lepton
    const pat::Jet* matchedJetPtr = nullptr;

    //old matching scheme looks for closest jet in terms of delta R, and required this to be within delta R 0.4 of the lepton
    if( oldMatching ){
        for( auto& jet : *jets ){
            if( jet.pt() <= 5 || fabs( jet.eta() ) >= 3 ) continue;
            if( ( matchedJetPtr == nullptr) || reco::deltaR( jet, lepton ) < reco::deltaR( *matchedJetPtr, lepton ) ){
                matchedJetPtr = &jet;
            }
        }
        if( matchedJetPtr != nullptr && reco::deltaR( lepton, *matchedJetPtr ) > 0.4 ){
            matchedJetPtr = nullptr;
        }
    } else {
        for( auto& jet : *jets ){
            if( isSourceCandidatePtrMatch( lepton, jet ) ){

                //immediately returning guarantees that the leading jet matched to the lepton is returned
                return &jet;
            }
        }
    }
    return matchedJetPtr;
}

/*
//compute closest jet variables using new matching scheme 
void LeptonAnalyzer::fillLeptonJetVariables( const reco::Candidate& lepton, edm::Handle< std::vector< pat::Jet > >& jets, const reco::Vertex& vertex, const double rho, const bool oldMatching )
{
    //find closest jet based on source candidate pointer matching
    const pat::Jet* matchedJetPtr = findMatchedJet( lepton, jets, oldMatching );

    if( matchedJetPtr == nullptr ){
        if( _lFlavor[_nL] == 1 ){
            _ptRatio[_nL]           = ( oldMatching ? 1. : 1. / ( 1. + _relIso0p4MuDeltaBeta[_nL] ) );
        } 
        else{
            _ptRatio[_nL]           = ( oldMatching ? 1. : 1. / (1. + _relIso0p4[_nL]));
            _ptRatio_Summer16[_nL]  = ( oldMatching ? 1. : 1. / (1. + _relIso0p4_Summer16[_nL]));
        }
        _ptRel[_nL]                     = 0;
        _selectedTrackMult[_nL]         = 0;
        _closestJetDeepFlavor_b[_nL]    = 0;
        _closestJetDeepFlavor_bb[_nL]   = 0;
        _closestJetDeepFlavor_lepb[_nL] = 0;
        _closestJetDeepFlavor[_nL]      = 0;
        _closestJetDeepCsv_b[_nL]       = 0;
        _closestJetDeepCsv_bb[_nL]      = 0;
        _closestJetDeepCsv[_nL]         = 0;
        _closestJetCsvV2[_nL]           = 0;
    } 
    else {
        const pat::Jet& jet     = *matchedJetPtr;

        auto rawJetP4           = jet.correctedP4("Uncorrected"); 
        auto leptonP4           = lepton.p4();

        bool leptonEqualsJet    = ( ( rawJetP4 - leptonP4 ).P() < 1e-4 );

        //if lepton and jet vector are equal set _ptRatio, _ptRel and track multipliticy to defaults 
        if( leptonEqualsJet && !oldMatching ){
            _ptRatio[_nL]           = 1;
            _ptRatio_Summer16[_nL]  = 1;
            _ptRel[_nL]             = 0;
            _selectedTrackMult[_nL] = 0;
        } 
        else {

            //remove all corrections above L1 from the lepton
            auto L1JetP4            = jet.correctedP4("L1FastJet");
            double L2L3JEC          = jet.pt()/L1JetP4.pt(); 
            auto lepAwareJetP4      = ( L1JetP4 - leptonP4 )*L2L3JEC + leptonP4;

            _ptRatio[_nL]           = lepton.pt() / lepAwareJetP4.pt();
            _ptRatio_Summer16[_nL]  = lepton.pt() / lepAwareJetP4.pt();

            //lepton momentum orthogonal to the jet axis
            //magnitude of cross-product between lepton and rest of jet 
            _ptRel[_nL ]            = leptonP4.Vect().Cross( (lepAwareJetP4 - leptonP4 ).Vect().Unit() ).R();

            _selectedTrackMult[_nL] = 0;
            for( const auto& daughterPtr : jet.daughterPtrVector() ){
                const pat::PackedCandidate& daughter = *( (const pat::PackedCandidate*) daughterPtr.get() );
            
                if( daughter.charge() == 0 )                                        continue;
                if( daughter.fromPV() < 2 )                                         continue;
                if( reco::deltaR( daughter, lepton ) > 0.4 )                        continue;
                if( !daughter.hasTrackDetails() )                                   continue;

                auto daughterTrack = daughter.pseudoTrack();
                if( daughterTrack.pt() <= 1 )                                       continue;
                if( daughterTrack.hitPattern().numberOfValidHits() < 8 )            continue;
                if( daughterTrack.hitPattern().numberOfValidPixelHits() < 2 )       continue;
                if( daughterTrack.normalizedChi2() >= 5 )                           continue;
                if( std::abs( daughterTrack.dz( vertex.position() ) ) >= 17 )       continue;
                if( std::abs( daughterTrack.dxy( vertex.position() ) ) >= 0.2 )     continue;
                ++_selectedTrackMult[_nL];
            }

        }

        //CSVv2 of closest jet
        _closestJetCsvV2[_nL]           = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        //DeepCSV of closest jet
        _closestJetDeepCsv_b[_nL]       = jet.bDiscriminator("pfDeepCSVJetTags:probb");
        _closestJetDeepCsv_bb[_nL]      = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
        _closestJetDeepCsv[_nL]         = _closestJetDeepCsv_b[_nL] + _closestJetDeepCsv_bb[_nL];
        if( std::isnan( _closestJetDeepCsv[_nL] ) ) _closestJetDeepCsv[_nL] = 0.;

        //DeepFlavor b-tag values of closest jet
        _closestJetDeepFlavor_b[_nL]    = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
        _closestJetDeepFlavor_bb[_nL]   = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
        _closestJetDeepFlavor_lepb[_nL] = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
        _closestJetDeepFlavor[_nL]      = _closestJetDeepFlavor_b[_nL] + _closestJetDeepFlavor_bb[_nL] + _closestJetDeepFlavor_lepb[_nL];
        if( std::isnan( _closestJetDeepFlavor[_nL] ) ) _closestJetDeepFlavor[_nL] = 0.;
    }
}
*/
