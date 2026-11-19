import os
import json
#import argparse
from CRABClient.UserUtilities import config
config = config()

# ---------------------------------------------------------
# HYPERPARAMETERS:
# ---------------------------------------------------------

# data configuration
year                                    = "2025" 
#era                                     = 'Run2024I_M1_v2'
era                                     = 'Run2025C_M0_PromptReco'
#era                                     = 'DYJetsTo2L_M50'
primary_dataset                         = "Muon"
process                                 = "V0Analyzer"

# user/processing configuration
version                                 = 1
user                                    = "jbierken"
cmssw                                   = "CMSSW_15_0_17"
nthreads                                = 1
cores                                   = 1
memory                                  = 2500                  # in MB
runTime                                 = 2200                  # 2750 ~45 hours (default is ~20h)

# ---------------------------------------------------------
# RUN CONFIGURATION:
# ---------------------------------------------------------

# Data or MC
isData                                  = False     if 'DYJets' in era else True
isPromptReco                            = 'PromptReco' in era
isBusyDataEra                           = isData and (era == 'Run2024F_M0' or era == 'Run2024F_M0_v2' or (year in ('2025', '2026') and isPromptReco))

dataType                                = 'data'    if isData else 'sim'

# Era-aware queue tuning:
defaultDataUnits                        = 300
# Lumis per job for busy (LumiBased) eras. Derivation: before the K0Analyzer speed-up, past tasks took
# ~0.45-0.7 h per lumi (2025C/F: 5 lumis/job -> 2.1-2.6 h median; 2025D: 10 lumis/job -> 6.7 h).
# The per-track caching + pre-fit mass filter cut CPU/event by ~3.2-4.9x (0.33 -> 0.07-0.10 s on 2024G data),
# so ~0.1-0.2 h/lumi -> 30 lumis/job gives ~4-6 h jobs and 1/6-1/3 of the old job count, far below the
# 10000-jobs-per-task cap. Check the median job runtime ('crab status --long') after the first jobs and adjust.
busyDataUnits                           = 30
nunits                                  = (busyDataUnits if isBusyDataEra else defaultDataUnits) if isData else 10

useLumiBasedForBusyDataEras             = True
# Busy eras: data is often only hosted at a few sites, so jobs stay idle/unsubmitted.
# The jobs are CPU-bound (~90% CPU efficiency), so reading the input remotely (xrootd/AAA) is cheap:
# run anywhere in the site whitelist below instead of only where the data sits.
useIgnoreLocalityForBusyDataEras        = True
useWhitelist                            = False

# Create storage location (if not already exist)
dbssavepath                             = f'/store/user/{user}/K0sAnalysis/NTuples/MINIAOD/{dataType}/v{version}'
if not os.path.exists('/pnfs/iihe/cms/' + dbssavepath):         os.makedirs('/pnfs/iihe/cms/' + dbssavepath)

# Read data from JSON:
try:
    with open(f'V0Analysis/V0Analyzer/data/Run{year}.json', 'r') as file:
        data_config                     = json.load(file)

except FileNotFoundError:
    print(f"Error: The file 'Run{year}.json' was not found.")

# ---------------------------------------------------------
# CRAB SETUP:
# ---------------------------------------------------------

## General config:
config.General.requestName              = f"{process}_Run{year}_MiniAOD_{dataType}_{era}"
config.General.workArea                 = 'crab_V0Analyzer'
config.General.transferOutputs          = True
config.General.transferLogs             = False                                                         # logs still available with 'crab getlog'

## JobType config:
config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'V0Analysis/V0Analyzer/python/ConfFile_cfg.py'                # your cmsRun config
config.JobType.allowUndistributedCMSSW  = True                                                          # useful if you're on CMSSW_14_X
config.JobType.inputFiles               = [f'V0Analysis/V0Analyzer/python/ConfFile_cfg.py']
config.JobType.pyCfgParams              = [
                                            f'isData={isData}', 
                                            f'campaign={year}', 
                                            f'era={era}', 
                                            f'dataset={primary_dataset}', 
                                            f'globaltag={data_config["eras"][era]["globalTag"]}'
                                        ]

config.JobType.numCores                 = cores
config.JobType.maxMemoryMB              = memory

## Data config:
config.Data.inputDataset                = data_config["eras"][era]["dataset"]
config.Data.inputDBS                    = 'global'
config.Data.outLFNDirBase               = dbssavepath
config.Data.outputDatasetTag            = data_config["eras"][era]["dataset"].split("/")[2]
config.Data.allowNonValidInputDataset   = True

if isData:
    # Busy eras often spend a long time idle; split into smaller lumi chunks there.
    config.Data.splitting               = 'LumiBased' if (isBusyDataEra and useLumiBasedForBusyDataEras) else 'Automatic'
    config.Data.lumiMask                = data_config["lumijson"]
    if isBusyDataEra and useIgnoreLocalityForBusyDataEras:
        config.Data.ignoreLocality      = True                                                          # requires a site whitelist (set below)
else:
    # For MC: use FileBased splitting of files 
    config.Data.splitting               = 'FileBased'
    config.JobType.maxJobRuntimeMin     = runTime
config.Data.unitsPerJob                 = nunits
config.Data.totalUnits                  = -1                                                            # -1 = process all files
config.Data.publication                 = False

## Site config:
config.Site.ignoreGlobalBlacklist       = False
config.Site.storageSite                 = 'T2_BE_IIHE'                                                  # or your site
if useWhitelist or getattr(config.Data, 'ignoreLocality', False):
    config.Site.whitelist               = [
                                            "T2_CH*", "T2_FR*", "T2_IT*", "T2_DE*", "T2_AT*", "T2_BE*", "T2_ES*",
                                            "T1_US*", "T2_US*", "T3_US*", "T1_FR*", "T1_IT*", "T1_DE*", "T2_UK*", "T3_UK*", "T2_FI*", "T2_EE*", "T1_ES*"
                                        ]


