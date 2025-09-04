#ifndef TRIGGER_ANALYZER_H
#define TRIGGER_ANALYZER_H

//include CMSSW classes
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//include ROOT classes
#include "TTree.h"

//include other parts of framework
#include "V0Analysis/V0Analyzer/plugins/V0Analyzer.h"

class V0Analyzer;

class TriggerAnalyzer {
    private:

        std::map<TString, std::vector<TString>> allFlags;
        std::vector<TString> triggersToSave;
        std::vector<TString> filtersToSave;

        std::map<TString, bool> flag;
        //std::map<TString, int>  prescale;
        std::map<TString, double>  prescale;
        std::map<TString, int>  index;

        V0Analyzer* myAnalyzer;

        void indexFlags(const edm::Event&, edm::Handle<edm::TriggerResults>&, std::vector<TString>&);
        void getResults(const edm::Event&, edm::Handle<edm::TriggerResults>&, std::vector<TString>&, const bool);

        void initList(std::vector<TString>&, TString);
        std::vector<TString> getAllFlags();

        bool passCombinedFlagAND(TString combinedFlag);
        bool passCombinedFlagOR(TString combinedFlag);

        bool passEle32WPTight(const edm::Event&, edm::Handle<edm::TriggerResults>&);
        std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjects(const pat::Electron&, const std::vector<pat::TriggerObjectStandAlone>&,const float);

    public:
        TriggerAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* vars);
        ~TriggerAnalyzer(){};

        bool reIndex;
        void beginJob(TTree* outputTree);
        void analyze(const edm::Event&);
};

#endif
