/*****************************************************************************************
 * Triggers and MET filters:
 *          - Simply add your triggers to the list below, 
 *            the key in allFlags[key] takes the OR of the following triggers
 *          - Currently not only the combined but also the individual triggers are stored, 
 *            keeping the possibility for trigger studies
 *          - Might add a flag to switch all the storage of all those individual paths off
 *
 *  NOTE: use "pass" in the combined flag if you want to use it recursively
 *****************************************************************************************/

//include CMSSW classes
#include "FWCore/Common/interface/TriggerNames.h"

//include other parts of code
#include "V0Analysis/V0Analyzer/src/TriggerAnalyzer.h"

TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& iConfig, V0Analyzer* myAnalyzer):myAnalyzer(myAnalyzer)
{
    //--------------------------------------------------------
    // MET Filters: 
    //      - first add common ones for 2016, 2017, 2018
    //      - MET filter are taken in AND (based on the occurence of capitalized 'MET' in the allFlags key) and always start with "Flag"
    //
    // REFERENCE: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
    //--------------------------------------------------------

    if( myAnalyzer->isRun2() ) {
        allFlags["passMETFilters"] = {
            "Flag_goodVertices", 
            "Flag_HBHENoiseFilter", 
            "Flag_HBHENoiseIsoFilter",  
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter"
        };

        allFlags["passBadChargedCandidateFilter"] = {
            "Flag_BadChargedCandidateFilter"
        };

        if( myAnalyzer->isData() ){ // This one is only to be applied on data
            allFlags["passMETFilters"].push_back("Flag_eeBadScFilter");
        }

        if( myAnalyzer->is2017() || myAnalyzer->is2018() ){ // This one is only for 2017 and 2018 data
            allFlags["passMETFilters"].push_back("updated_ecalBadCalibFilter"); // improved version over the Flag_ecalBadCalibFilter, implementation manually
        }
    }
    else if(myAnalyzer->isRun3()) 
    {
        allFlags["passMETFilters"] = {
            "Flag_goodVertices",
            "Flag_globalSuperTightHalo2016Filter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter",
            "Flag_BadPFMuonDzFilter",
            "Flag_hfNoisyHitsFilter",
            "Flag_eeBadScFilter",
            "Flag_ecalBadCalibFilter"
        };

    }

    //--------------------------------------------------------
    // Triggers: 
    //      - grouped per year
    //      - Triggers are taken in OR and always start with "HLT"
    //      - To check if triggers are existing/prescaled/unused, use https://tomc.web.cern.ch/tomc/triggerPrescales
    //
    // WARNING: several triggers are off for part of the datataking, this is typically mentioned in the comments, preferably offline cuts are higher than the unprescaled ones
    // TODO: maybe we should clean up this part by storing it in some configuration file which can be analysis-specific
    //--------------------------------------------------------
    
    // ------------------ 2025 triggers ------------------
    if( myAnalyzer->is2025() ){
        allFlags["passTrigger_mm"] = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "HLT_Mu37_TkMu27"
        };
    }
    
    // ------------------ 2024 triggers ------------------
    else if( myAnalyzer->is2024() ){
        allFlags["passTrigger_mm"] = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30",
            "HLT_Mu37_TkMu27"
        };
    }
    
    // ------------------ 2023 triggers ------------------
    else if( myAnalyzer->is2023() || myAnalyzer->is2023BPix()){
        allFlags["passTrigger_mm"] = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "HLT_Mu37_TkMu27"
        };
    }
    
    // ------------------ 2022 triggers ------------------
    else if( myAnalyzer->is2022()  || myAnalyzer->is2022EE()){
        allFlags["passTrigger_mm"] = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "HLT_Mu37_TkMu27"
        };
    }
    
    // ------------------ 2018 triggers ------------------
    else if( myAnalyzer->is2018() ){
        allFlags["passTrigger_mm"] = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", 
            "HLT_Mu37_TkMu27"
        };
    } 

    // ------------------ 2017 triggers ------------------
    // WARNING: very little triggers available which were used throughout the whole 2017 dataset
    // Maybe another trigger strategy will be needed if insufficient statistics
    else if( myAnalyzer->is2017() ){
        allFlags["passTrigger_mm"] = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",                                         // first two triggers heavily prescaled at some /fb,
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", 
            "HLT_Mu37_TkMu27",                                                          // HLT_Mu37_TrkMu27 trigger missing for first ~14/fb 
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"                               // Mass3p8 trigger off at start of 2017

        };
    }

    // ------------------ 2016 triggers ------------------
    else {
        allFlags["passTrigger_mm"] = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", 
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                                      // non-DZ version heavily prescaled for a few /fb at end of 2016
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", 
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",                                    // same as above
            "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", 
            "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",                                  // non-DZ version heavily prescaled for a few /fb at end of 2016, DZ-version only introduced at end of 2016
            "HLT_Mu30_TkMu11"
        };
    }
}

void TriggerAnalyzer::beginJob(TTree* outputTree)
{
    reIndex = true;
    initList(triggersToSave, "HLT");
    initList(filtersToSave, "Flag");

    for(TString f : getAllFlags()){
        outputTree->Branch("_" + f, &flag[f], "_" + f + "/O");
        if(f.Contains("HLT")) outputTree->Branch("_" + f + "_prescale", &prescale[f], "_" + f + "_prescale/I");
    }
}


// Filters the HLT and MET flags from allFlags (maybe better to change to set...)
void TriggerAnalyzer::initList(std::vector<TString>& list, TString identifier)
{
    list.clear();
    for(auto& v : allFlags){
        for(auto& t : v.second){
            if(t.Contains(identifier)){
                if(std::find(list.begin(), list.end(), t) == list.end()){
                    list.push_back(t);
                }
            }
        }
    }
}


std::vector<TString> TriggerAnalyzer::getAllFlags()
{
    std::vector<TString> list;
    for(auto& v : allFlags){
        list.push_back(v.first);
        for(auto& t : v.second){
            if(std::find(list.begin(), list.end(), t) == list.end()){
                list.push_back(t);
            }
        }
    }
    return list;
}

bool TriggerAnalyzer::passCombinedFlagOR(TString combinedFlag)
{
    for(auto& f : allFlags[combinedFlag]){
        if(f.Contains("pass") and passCombinedFlagOR(f))        return true;
        else if(flag[f])                                        return true;
    }
    return false;
}

bool TriggerAnalyzer::passCombinedFlagAND(TString combinedFlag)
{
    for(auto& f : allFlags[combinedFlag]){
        if(f.Contains("pass") and not passCombinedFlagAND(f))   return false;
        else if(not flag[f])                                    return false;
    }
    return true;
}

// Matching trigger objects within maxDeltaR to the electron supercluster eta/phi
// It is important to match to ALL objects as there are different ways to reconstruct the same electron (e.g. L1 seeded, unseeded, jet,...)
std::vector<const pat::TriggerObjectStandAlone*> TriggerAnalyzer::getMatchedObjects(const pat::Electron& ele, const std::vector<pat::TriggerObjectStandAlone>& trigObjs, const float maxDeltaR)
{
    std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
    const float maxDR2 = maxDeltaR*maxDeltaR;
    for(auto& trigObj : trigObjs){
        const float dR2 = reco::deltaR2(ele.superCluster()->eta(), ele.superCluster()->phi(), trigObj.eta(), trigObj.phi());
        if(dR2 < maxDR2) matchedObjs.push_back(&trigObj);
    }
    return matchedObjs;
}

bool TriggerAnalyzer::passEle32WPTight(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& triggerResults)
{
    auto electrons      = getHandle(iEvent, myAnalyzer->eleToken);
    auto triggerObjects = getHandle(iEvent, myAnalyzer->trigObjToken);

    // unpack the trigger objects
    std::vector<pat::TriggerObjectStandAlone> unpackedTriggerObjects;
    for(auto& trigObj : *triggerObjects){
        unpackedTriggerObjects.push_back(trigObj);
        unpackedTriggerObjects.back().unpackFilterLabels(iEvent, *triggerResults);
    }

    // Check if there's an electron matched to a trigger object which passes the two filters
    for(auto& ele : *electrons){
        for(const auto trigObj : getMatchedObjects(ele, unpackedTriggerObjects, 0.1)){
            if(!trigObj->hasFilterLabel("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter"))  continue;
            if(!trigObj->hasFilterLabel("hltEGL1SingleEGOrFilter"))                     continue;
            return true;
        }
    }
    return false;
}

void TriggerAnalyzer::analyze(const edm::Event& iEvent)
{
    edm::Handle<edm::TriggerResults> recoResults    = getHandle(iEvent, myAnalyzer->recoResultsPrimaryToken);
    if(recoResults.failedToGet())    recoResults    = getHandle(iEvent, myAnalyzer->recoResultsSecondaryToken);
    edm::Handle<edm::TriggerResults> triggerResults = getHandle(iEvent, myAnalyzer->triggerToken);

    if( myAnalyzer->is2017() || myAnalyzer->is2018() ){         // The updated ecalBadCalibFilter
        edm::Handle<bool> passEcalBadCalibFilterUpdate = getHandle(iEvent, myAnalyzer->ecalBadCalibFilterToken);
        flag["updated_ecalBadCalibFilter"] = (*passEcalBadCalibFilterUpdate);
    }
    if( myAnalyzer->is2022() || myAnalyzer->is2022EE() || myAnalyzer->is2023() || myAnalyzer->is2023BPix() ){         // The updated ecalBadCalibFilter
        edm::Handle<bool> passEcalBadCalibFilterUpdate = getHandle(iEvent, myAnalyzer->ecalBadCalibFilterToken);
        flag["updated_ecalBadCalibFilter"] = (*passEcalBadCalibFilterUpdate);
    }

    // Get all flags
    getResults(iEvent, triggerResults, triggersToSave, true);
    getResults(iEvent, recoResults,    filtersToSave,  false);

    // In 2017: emulate the non-existing HLT_Ele32_WPTight_Gsf
    //if(myAnalyzer->is2017() && ! myAnalyzer->isFastSim() ){
    if(myAnalyzer->is2017() ){
        flag["HLT_Ele32_WPTight_Gsf"] = passEle32WPTight(iEvent, triggerResults);
    }

    reIndex = false;

    for(auto& combinedFlag : allFlags){
        if(combinedFlag.first.Contains("MET")) flag[combinedFlag.first] = passCombinedFlagAND(combinedFlag.first);
        else                                   flag[combinedFlag.first] = passCombinedFlagOR(combinedFlag.first);
    }
}

/*
 * Call this at the first event, checks for available triggers and warns for missing triggers
 * Stores indexes of wanted triggers, which minimizes string comparisons for the next events
 */
void TriggerAnalyzer::indexFlags(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& results, std::vector<TString>& toSave)
{
    for(TString t : toSave) index[t] = -1;

    std::cout << "Available triggers:" << std::endl;
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*results);
    for (unsigned i = 0; i < results->size(); ++i){
        std::cout << "  " << triggerNames.triggerName(i);
        for(TString t : toSave){
            TString tt = (t.Contains("HLT") ? t + "_v" : t);
            if(TString(triggerNames.triggerName(i)).Contains(tt)){
                index[t] = i;
                std::cout << "     --> saving to tree";
            }
        }
        std::cout << std::endl;
    }

    for(TString t : toSave){
        if(index[t] == -1) std::cout << "WARNING: " << t << " not found in triggerresult, please check!" << std::endl;
    }
}


 //Saving triggers and prescales
void TriggerAnalyzer::getResults(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& results, std::vector<TString>& toSave, const bool savePrescales)
{
    if(results.failedToGet()) return;

    if(reIndex) indexFlags(iEvent, results, toSave);

    for(TString t : toSave){
        if(index[t] == -1) flag[t] = false;
        else               flag[t] = results->accept(index[t]);
    }

    if(savePrescales){
        edm::Handle<pat::PackedTriggerPrescales> prescales = getHandle(iEvent, myAnalyzer->prescalesToken);
        for(TString t : toSave){
            if(index[t] == -1) prescale[t] = -1;
            else               prescale[t] = (double) prescales->getPrescaleForIndex<double>(index[t]);
        }
    }
}
