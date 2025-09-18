import os
import json
from CRABClient.UserUtilities import config
config = config()

# ---------------------------------------------------------
# HYPERPARAMETERS:
# ---------------------------------------------------------

# data configuration
version                                 = 0
isData                                  = False 
year                                    = "2022"
era                                     = "DYJetsToLL_M50"
primary_dataset                         = "DoubleMuon"
process                                 = "V0Analyzer"

# user/processing configuration
user                                    = 'jbierken'
cmssw                                   = "CMSSW_14_0_15"
nunits                                  = 200
nthreads                                = 1
cores                                   = 1
memory                                  = 2000                                                          # in MB

# ---------------------------------------------------------
# RUN CONFIGURATION:
# ---------------------------------------------------------

# Data or MC
if isData:                              dataType = 'data'
else:                                   dataType = 'sim'

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
config.General.requestName              = f"{process}_MiniAOD_{dataType}_{era}"
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

#config.Data.splitting                   = 'FileBased'
config.Data.splitting                   = 'Automatic'
config.Data.lumiMask                    = data_config["lumijson"]

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


