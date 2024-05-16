import sys, os
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

## Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

### Data global tag
process.GlobalTag.globaltag = "106X_dataRun2_v35"

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source = cms.Source("PoolSource",
        # replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
            '/store/data/Run2018A/SingleMuon/MINIAOD/UL2018_MiniAODv2-v3/2530000/002A113D-FB15-1341-A170-638E53A7261F.root'
            #'/store/data/Run2018A/EGamma/MINIAOD/UL2018_MiniAODv2-v1/230000/1DC29AF8-7091-4245-A0D8-CFDF650310CC.root'
        ),
)

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("SingleMuon_Run2018A-UL2018.root")
)
#fileName = cms.string("EGamma_Run2018A-UL2018.root"))

process.demo = cms.EDAnalyzer('ECPTreeMaker',
        l1Alg            = cms.InputTag("gtStage2Digis",                 "", "RECO"),
        l1Ext            = cms.InputTag("gtStage2Digis",                 "", "RECO"),
        triggerResults   = cms.InputTag("TriggerResults",                "", "HLT"),
        triggerPrescales = cms.InputTag("patTrigger",                    "", "PAT"),
        metFilters       = cms.InputTag("TriggerResults",                "", "PAT"),
        muons            = cms.InputTag("slimmedMuons",                  "", "PAT"),
        #electrons        = cms.InputTag("slimmedElectrons",              "", "PAT"),
        #electronLooseID  = cms.string("mvaEleID-Fall17-iso-V2-wp90"),
        #electronTightID  = cms.string("mvaEleID-Fall17-iso-V2-wp80"),
        #caloJets         = cms.InputTag("slimmedCaloJets",               "", "PAT"),
        #ecalRecHits      = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
        tracks           = cms.InputTag("isolatedTracks",                "", "PAT"),
        primaryVertices  = cms.InputTag("offlineSlimmedPrimaryVertices", "", "PAT"),
        beamSpot         = cms.InputTag("offlineBeamSpot",               "", "RECO"),
        dedx             = cms.InputTag("isolatedTracks",                "", "PAT"),
        met              = cms.InputTag("slimmedMETs",                   "", "PAT")
)


process.p = cms.Path(process.demo))
