import os
import json
#import argparse
from CRABClient.UserUtilities import config
config = config()

# ---------------------------------------------------------
# HYPERPARAMETERS:
# ---------------------------------------------------------

# data configuration
year                                    = "2024" 
#era                                     = 'Run2024G_M1'
era                                     = 'DYJetsTo2L_M50'
primary_dataset                         = "Muon"
process                                 = "V0Analyzer"

# user/processing configuration
version                                 = 1
user                                    = "jbierken"
cmssw                                   = "CMSSW_15_0_17"
nthreads                                = 1
cores                                   = 1
memory                                  = 2500                  # in MB
runTime                                 = 2750                  # ~45 hours (default is ~20h)

# ---------------------------------------------------------
# RUN CONFIGURATION:
# ---------------------------------------------------------

# Data or MC
isData                                  = False     if 'DYJets' in era else True

dataType                                = 'data'    if isData else 'sim'
nunits                                  = 50        if isData else 5

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
config.General.transferLogs             = True

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
    # For Data: Automatic splitting with lumimasking gives best configuration
    #config.Data.splitting               = 'Automatic'
    config.Data.splitting               = 'LumiBased'
    config.Data.lumiMask                = data_config["lumijson"]
else:
    # For MC: use FileBased splitting of files 
    config.Data.splitting               = 'FileBased'
    config.JobType.maxJobRuntimeMin     = runTime
config.Data.unitsPerJob                 = nunits
config.Data.totalUnits                  = -1                                                            # -1 = process all files
config.Data.publication                 = False

## Site config:
config.Site.ignoreGlobalBlacklist       = True
config.Site.storageSite                 = 'T2_BE_IIHE'                                                  # or your site
config.Site.whitelist                   = [
                                            "T2_CH*", "T2_FR*", "T2_IT*", "T2_DE*", "T2_AT*", "T2_BE*", "T2_ES*",
                                            "T1_US*", "T2_US*", "T3_US*", "T1_FR*", "T1_IT*", "T1_DE*", "T2_UK*", "T3_UK*", "T2_FI*", "T2_EE*", "T1_ES*"
                                        ]


