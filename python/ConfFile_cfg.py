import sys, os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

## ---------------------------------------
### parser (for crab submission)
## ---------------------------------------
options = VarParsing.VarParsing('analysis')
 
options.register("isData",
                False, 
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.bool,
                "Are we considering Data/MC"
)
 
options.register ("campaign",
                "2022",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "Data-taking/MC Campaign"
)
 
options.register ("era",
                "Run2022B",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "Data-taking era for data"
)
 
options.register ("dataset",
                "DoubleMuon",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "which Primary-Dataset are we considering"
)
 
options.register ("globaltag",
                "130X_dataRun3_v2",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "globalTag to be used"
)

options.parseArguments()


## ---------------------------------------
# Default Arguments:
## ---------------------------------------
isData                  = options.isData
# Data taking year:
#   run-2
is2017                      = "2017"        in options.campaign 
is2018                      = "2018"        in options.campaign 
is2016preVFP                = "preVFP"      in options.campaign 
#   run-3
is2022                      = "2022"        in options.campaign 
is2022EE                    = "2022EE"      in options.campaign 
is2023                      = "2023"        in options.campaign 

## ---------------------------------------
# Start process
## ---------------------------------------
process                     = cms.Process("BlackJackAndHookers")

## Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


### Data global tag
#process.GlobalTag.globaltag= "106X_dataRun2_v35"
process.GlobalTag.globaltag = options.globaltag 

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery= 10000
process.maxEvents           = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source              = cms.Source("PoolSource",
                                fileNames           = cms.untracked.vstring(
                                                        # Run-2 (locally)
                                                        #'/store/data/Run2018A/DoubleMuon/MINIAOD/UL2018_MiniAODv2-v1/260000/00264E65-8EFD-974C-8A29-866EFA1609D3.root'
                                                        # Run-3 (locally)
                                                        #'root://cms-xrd-global.cern.ch/'    + '/store/data/Run2022B/DoubleMuon/MINIAOD/PromptReco-v1/000/355/207/00000/e6f8e5ae-70e0-4dac-9106-b81723a98c7a.root',
                                                        # for Crab job submission
                                                        options.inputFiles
                                                    ),
                                secondaryFileNames  = cms.untracked.vstring(),
                                duplicateCheckMode  = cms.untracked.string('noDuplicateCheck')
                            )

## ---------------------------------------
## Create output file
## ---------------------------------------

## Setup the service to make a ROOT TTree
process.TFileService        = cms.Service("TFileService",
                                fileName = cms.string(f"output_{options.dataset}_{options.era}.root")
                            )

#process.demo = cms.EDAnalyzer('ECPTreeMaker',
process.blackJackAndHookers = cms.EDAnalyzer('V0Analyzer',
                                # Beam spot & vertices
                                offlineBeamSpot         = cms.InputTag("offlineBeamSpot",               "", "RECO"),
                                vertices                = cms.InputTag("offlineSlimmedPrimaryVertices"),

                                # Generator information
                                genEventInfo            = cms.InputTag("generator"),
                                lheEventInfo            = cms.InputTag("externalLHEProducer"),
                                pileUpInfo              = cms.InputTag("slimmedAddPileupInfo"),
                                genParticles            = cms.InputTag("prunedGenParticles"),

                                # Particle-level objects (MINIAOD doesn't store these directly usually from NanoAOD or Rivet step)
                                particleLevelPhotons    = cms.InputTag("particleLevel",                     "photons"),     # adjust if not produced
                                particleLevelLeptons    = cms.InputTag("particleLevel",                     "leptons"),
                                particleLevelJets       = cms.InputTag("particleLevel",                     "jets"),
                                particleLevelMets       = cms.InputTag("particleLevel",                     "met"),


                                # Reco / PAT objects
                                muons                   = cms.InputTag("slimmedMuons"),
                                electrons               = cms.InputTag("slimmedElectrons"),
                                taus                    = cms.InputTag("slimmedTaus"),
                                tauGenJets              = cms.InputTag("tauGenJetsSelectorAllHadrons"),                     # from tau tools
                                photons                 = cms.InputTag("slimmedPhotons"),
                                packedCandidates        = cms.InputTag("packedPFCandidates"),
                                lostTracks              = cms.InputTag("lostTracks"),
                                fixedGridRhoFastjetAll  = cms.InputTag("fixedGridRhoFastjetAll"),
                                met                     = cms.InputTag("slimmedMETs"),
                                metPuppi                = cms.InputTag("slimmedMETsPuppi"),
                                jets                    = cms.InputTag("slimmedJets"),
                                jetsPuppi               = cms.InputTag("slimmedJetsPuppi"),

                                # JES-smeared jets (only present if you run JME modules yourself)
                                jetsSmeared             = cms.InputTag("slimmedJetsSmeared"),
                                jetsSmearedUp           = cms.InputTag("slimmedJetsSmearedUp"),
                                jetsSmearedDown         = cms.InputTag("slimmedJetsSmearedDown"),
                                
                                # Trigger info
                                recoResultsPrimary      = cms.InputTag("TriggerResults",                "", "RECO"),
                                recoResultsSecondary    = cms.InputTag("TriggerResults",                "", "RECO"),
                                triggers                = cms.InputTag("TriggerResults",                "", "HLT"),
                                prescales               = cms.InputTag("patTrigger"),
                                triggerObjects          = cms.InputTag("slimmedPatTrigger"),

                                # Additional input parameters
                                skim                    = cms.untracked.string("noskim"),
                                isData                  = cms.untracked.bool(isData),
                                is2016preVFP            = cms.untracked.bool(is2016preVFP),
                                is2017                  = cms.untracked.bool(is2017),
                                is2018                  = cms.untracked.bool(is2018),
                                is2022                  = cms.untracked.bool(is2022),
                                is2022EE                = cms.untracked.bool(is2022EE),
                                is2023                  = cms.untracked.bool(is2023),

                                # Which info should be collected
                                storeLheParticles       = cms.untracked.bool(True),
                                storeGenParticles       = cms.untracked.bool(True),
                                storeParticleLevel      = cms.untracked.bool(True),
                                storeJecSourcesAll      = cms.untracked.bool(True),
                                storeJecSourcesGrouped  = cms.untracked.bool(True),
                                storeAllTauID           = cms.untracked.bool(True),
                                storePrefireComponents  = cms.untracked.bool(True),
                                storeJetSubstructure    = cms.untracked.bool(True),
                            )

process.p                   = cms.Path(process.blackJackAndHookers)
