#ifndef JET_ANALYZER_H
#define JET_ANALYZER_H
//CMSSW
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

//c++ standard library
#include <memory>

//ROOT
#include "TTree.h"

//include other parts of framework
#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"

class V0Analyzer;

class JetAnalyzer {
   
    //Declare friend classes
    friend class V0Analyzer;
    
    private:
        JetCorrectionUncertainty* jecUnc;
        JetCorrectionUncertainty* jecUncPuppi;

        // Maxim number of Jets/PFcandidates
        static const unsigned nJets_max         = 100;
        static const unsigned nPFCandidates_max = 1000;

        std::map<std::string, std::shared_ptr< JetCorrectorParameters> > jetSourcesCorParameters;
        std::map<std::string, std::shared_ptr< JetCorrectorParameters> > jetGroupedCorParameters;

        std::shared_ptr<FactorizedJetCorrector> jetCorrector;

        // Initiate counter for number of Jets
        unsigned _nJets     = 0;

        // Jet kinematics
        double   _jetPt[nJets_max];
        double   _jetEta[nJets_max];
        double   _jetPhi[nJets_max];
        double   _jetE[nJets_max];
        
        // Flavor tagging: Deep CSV
        double   _jetCsvV2[nJets_max];
        
        // Flavor tagging: Deep CSV
        double   _jetDeepCsv_udsg[nJets_max];
        double   _jetDeepCsv_b[nJets_max];
        double   _jetDeepCsv_c[nJets_max];
        double   _jetDeepCsv_bb[nJets_max];
        double   _jetDeepCsv[nJets_max];

        // Flavor tagging: Deep CSV
        double   _jetDeepFlavor_b[nJets_max];
        double   _jetDeepFlavor_bb[nJets_max];
        double   _jetDeepFlavor_lepb[nJets_max];
        double   _jetDeepFlavor[nJets_max];
        double   _jetDeepFlavor_c[nJets_max];
        double   _jetDeepFlavor_uds[nJets_max];
        double   _jetDeepFlavor_g[nJets_max];
        
        // Flavor tagging: Puppi Jets 
        unsigned _nJetsPuppi = 0;
        double   _jetPuppiPt[nJets_max];
        double   _jetPuppiEta[nJets_max];
        double   _jetPuppiPhi[nJets_max];

        //metPuppi kinematics
        double   _metPuppi;
        double   _metPuppiPhi;
    
        //met kinematics
        double   _met;
        double   _metPhi;

        // Jet ID requirements
        unsigned _jetHadronFlavor[nJets_max];
        unsigned _jetPartonFlavor[nJets_max];
        bool     _jetIsLoose[nJets_max];
        bool     _jetIsTight[nJets_max];
        bool     _jetIsTightLepVeto[nJets_max];

        // Gen-level Jet kinemmatics
        bool     _jetHasGen[nJets_max];
        double   _jetGenPt[nJets_max];
        double   _jetGenEta[nJets_max];
        double   _jetGenPhi[nJets_max];
        double   _jetGenE[nJets_max];

        /*
        //split JEC in different sources
        std::map< std::string, double[nJets_max] >  _jetPt_groupedVariationsDown;
        std::map< std::string, double[nJets_max] >  _jetPt_groupedVariationsUp;
        std::map< std::string, double[nJets_max] >  _jetPt_allVariationsDown;
        std::map< std::string, double[nJets_max] >  _jetPt_allVariationsUp;

        std::map< std::string, double[nJets_max] >  _jetSmearedPt_groupedVariationsDown;
        std::map< std::string, double[nJets_max] >  _jetSmearedPt_groupedVariationsUp;
        std::map< std::string, double[nJets_max] >  _jetSmearedPt_allVariationsDown;
        std::map< std::string, double[nJets_max] >  _jetSmearedPt_allVariationsUp;

        // MET with propagated JEC sources and uncertainties
        std::map< std::string, double >             _corrMETx_groupedVariationsDown;
        std::map< std::string, double >             _corrMETx_groupedVariationsUp;
        std::map< std::string, double >             _corrMETy_groupedVariationsDown;
        std::map< std::string, double >             _corrMETy_groupedVariationsUp;

        std::map< std::string, double >             _corrMETx_allVariationsDown;
        std::map< std::string, double >             _corrMETx_allVariationsUp;
        std::map< std::string, double >             _corrMETy_allVariationsDown;
        std::map< std::string, double >             _corrMETy_allVariationsUp;
        */

        V0Analyzer* myAnalyzer;

        bool jetIsLoose(const pat::Jet& jet, const bool is2017) const;
        bool jetIsTight(const pat::Jet& jet, const bool is2017, const bool is2018) const;
        bool jetIsTightLepVeto(const pat::Jet& jet, const bool is2017, const bool is2018) const;

        /*
        std::vector<float> getSubCorrections(double rawPt, double eta, double rho, double area);
        std::pair<double, double> getMETCorrectionPxPy(double corrPt, double rawPt, double rawEta, double rawMuonSubtractedPt, double phi, double emf, double rho, double area, const std::string& source, unsigned jetIndex, double jecShift);
        */

    public:
        JetAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* vars);
        ~JetAnalyzer();

        void beginJob(TTree* outputTree);
        bool analyze(const edm::Event&);

        //void correctedMETAndPhi(const pat::MET& met, const std::vector< pat::Jet >& jets, const double rho);
};

#endif
