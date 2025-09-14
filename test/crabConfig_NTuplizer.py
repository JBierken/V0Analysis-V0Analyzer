import os
import datetime
from CRABClient.UserUtilities import config
config = config()

# ---------------------------------------------------------
# HYPERPARAMETERS:

user                                    = 'jbierken'
cmssw                                   = "CMSSW_14_0_15"
isData                                  = True 
year                                    = 2022
era                                     = "Run2022B"
primary_dataset                         = "DoubleMuon"
process                                 = "V0Analyzer"
version                                 = 0

nunits                                  = 200
nthreads                                = 1

# ---------------------------------------------------------
# RUN CONFIGURATION:

# Data or MC
if isData:
    dataType                            = 'data'
else:
    dataType                            = 'sim'


# get current time for folder naming
time                                    = datetime.datetime.now()
current_time                            = f'{time.year}{time.month}{time.day}_{time.hour}{time.minute}{time.second}'

dbssavepath                             = f'/store/user/{user}/K0sAnalysis/NTuples/MINIAOD/{dataType}/v{version}'
if not os.path.exists('/pnfs/iihe/cms/' + dbssavepath):     os.makedirs('/pnfs/iihe/cms/' + dbssavepath)

lumijson                                = {
                                            2022        : "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
                                            2023        : "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json",
                                            2024        : "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json",
                                        }
datasets                                = {
                                            # Run2022
                                            "Run2022A"  : f"/{primary_dataset}/Run2022A-22Sep2023-v1/MINIAOD",
                                            "Run2022B"  : f"/{primary_dataset}/Run2022B-22Sep2023-v1/MINIAOD",
                                            "Run2022C"  : f"/{primary_dataset}/Run2022C-22Sep2023-v1/MINIAOD",
                                            # Run2023
                                            # Run2024
                                            "Run2022E"  : f"/{primary_dataset}/Run2022A-PromptReco-v1/MINIAOD",     # Not correct!
                                        }


# ---------------------------------------------------------
# CRAB SETUP:

# General
config.General.requestName              = f"{process}_MiniAOD_{dataType}_{era}_v{version}"                                                                                                                                          
#config.General.requestName              = f'V0Analyzer_Run3_MINIAOD_{era}'
config.General.workArea                 = 'crab_V0Analyzer'
config.General.transferOutputs          = True
config.General.transferLogs             = True

# JobType
config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'V0Analysis/V0Analyzer/python/ConfFile_cfg.py'                # your cmsRun config
config.JobType.allowUndistributedCMSSW  = True                                                          # useful if you're on CMSSW_14_X
#config.JobType.inputFiles               = [f''V0Analysis/V0Analyzer/python/ConfFile_cfg.py']
#config.JobType.pyCfgParams              = ['isMC=False', 'campaign=2024']

config.JobType.numCores                 = 1
config.JobType.maxMemoryMB              = 2500


# Data
config.Data.inputDataset                = datasets[era]
config.Data.inputDBS                    = 'global'
config.Data.outLFNDirBase               = dbssavepath
config.Data.outputDatasetTag            = datasets[era].split("/")[2]
config.Data.allowNonValidInputDataset   = True

#config.Data.splitting                   = 'FileBased'
config.Data.splitting                   = 'Automatic'
config.Data.lumiMask                    = lumijson[year]

config.Data.unitsPerJob                 = nunits
config.Data.totalUnits                  = -1                                            # -1 = process all files
config.Data.publication                 = False

# Site
config.Site.ignoreGlobalBlacklist       = True
config.Site.storageSite                 = 'T2_BE_IIHE'                              # or your site
config.Site.whitelist                   = [
                                            "T2_CH*", "T2_FR*", "T2_IT*", "T2_DE*", "T2_AT*", "T2_BE*", "T2_ES*",
                                            "T1_US*", "T2_US*", "T3_US*", "T1_FR*", "T1_IT*", "T1_DE*", "T2_UK*", "T3_UK*", "T2_FI*", "T2_EE*", "T1_ES*"
                                        ]


