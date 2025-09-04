# Structure of the V0Analysis/V0Analyzer module

This directory follows the general CMSSW module layout, see
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBuildFile#CmsswSrcCodeDir

In particular for the V0Analysis/V0Analyzer module following directories are used:

### data
Contains JSON and .txt configuration files needed to run.
Putting those in the data directory ensures they are shipped with crab.

### interface
Header-files for the classes which are implemented in ./src
In order to structure the code, this module uses sub-analyzers, each of them aimed at calculating and storing the variables for a given object type:
* TriggerAnalyzer
* LheAnalyzer
* GenAnalyzer
* PartlcleLevelAnalyzer
* JetAnalyzer
* LeptonAnalyzer
* K0Analyzer

### src
Implementation of the sub-analyzers. Note that for the leptonAnalyzer, the identification and isolation functions are stored in separate .cc files in order to improve readability.

### plugins
Contains mulilep.h and multilep.cc, which form the main plugin of this module. Focusing on keeping track of all the tokens retrieved from the parameters the module is given, and organises the main order of how the sub-analyzers are run.
Note the LheAnalyzer should always be run before a skimming sub-analyzer, and that GenAnalyzer should be run before PhotonAnalyzer.

### python
[Placeholder]

### test
Directory with the main python config file multilep.py as well as submission scripts
