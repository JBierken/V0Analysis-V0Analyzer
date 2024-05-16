#ifndef LHE_ANALYZER_H
#define LHE_ANALYZER_H

//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

//include ROOT classes
#include "TTree.h"

//include other parts of framework
#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"

class V0Analyzer;

class LheAnalyzer {
    public:
        double getWeight() const;
    private:
        float  _nTrueInt;
        double _weight;
        double _lheHTIncoming;
        double _ctauHN;

        TH1D*  hCounter;
        TH1D*  lheCounter;
        TH1D*  psCounter;
        TH1D*  tauCounter;
        TH1D*  nTrueInteractions;

        static constexpr unsigned maxNumberOfLheWeights = 148;
        unsigned _nLheWeights = 0;
        unsigned _nTau;
        double _lheWeight[ maxNumberOfLheWeights ];

        static constexpr unsigned maxNumberOfPsWeights = 46;
        unsigned _nPsWeights = 0;
        double _psWeight[ maxNumberOfPsWeights ];

        static constexpr unsigned nLhe_max = 25;  // maximum number of LHE particles stored (the exact number of LHE particles will typically be the same for all events of a given process)
        unsigned              _nLheParticles = 0;
        int                   _lheStatus[nLhe_max];
        int                   _lhePdgId[nLhe_max];
        int                   _lheMother1[nLhe_max];
        int                   _lheMother2[nLhe_max];
        float                 _lhePt[nLhe_max];
        float                 _lheEta[nLhe_max];
        float                 _lhePhi[nLhe_max];
        float                 _lheE[nLhe_max];
        float                 _lheMass[nLhe_max];

        V0Analyzer* myAnalyzer;

    public:
        LheAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* vars);
        ~LheAnalyzer(){};

        void beginJob(TTree* outputTree, edm::Service<TFileService>& fs);
        void analyze(const edm::Event&);
};
#endif
