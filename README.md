# V0Analysis/V0Analyzer

NTuplizer for a displaced-vertexing study using V0 candidates
(K<sup>0</sup><sub>S</sub> → π<sup>+</sup>π<sup>−</sup>, Λ → pπ<sup>−</sup>, Λ̄ → p̄π<sup>+</sup>)
in CMS Run 3 data and simulation (2022–2026; some Run 2 flags are kept).
It reads MiniAOD and writes a flat ROOT `TTree`.

The module follows the general CMSSW package layout, see
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBuildFile#CmsswSrcCodeDir

## Directory structure

### plugins
`V0Analyzer.h` / `V0Analyzer.cc` form the main `EDAnalyzer` plugin. It:
* registers all input tokens from the configuration parameters,
* holds the era flags (`isData()`, `isRun3()`, `is2022EE()`, …, `is2026()`),
* books the output tree `blackJackAndHookersTree` and the `nVertices` histogram,
* runs the sub-analyzers in a fixed order (see [Event flow](#event-flow)).

The sub-analyzers are declared as `friend`s of `V0Analyzer` so they can access its tokens and flags.

### src
Header and implementation files of the sub-analyzers. Each one books and fills the branches for one object type:

| Sub-analyzer            | Content |
|-------------------------|---------|
| `K0Analyzer`            | V0 reconstruction and branches (`_nV0s`, `_V0*`), beamspot and primary vertex position |
| `LeptonAnalyzer`        | Muons, electrons, taus (`TauTools` contains tau helpers) |
| `JetAnalyzer`           | AK4 CHS and Puppi jets, JEC and uncertainty sources, b-tagging (DeepCSV, DeepJet, ParticleNet), MET |
| `TriggerAnalyzer`       | HLT decisions, prescales, MET filters |
| `LheAnalyzer`           | LHE and event weights, pile-up, counter histograms (MC only) |
| `GenAnalyzer`           | Generator particles (MC only; `GenTools` contains helpers) |
| `ParticleLevelAnalyzer` | Particle-level objects (MC only, if `storeParticleLevel=True`) |

`src/triggers/` holds the trigger lists for 2016–2018 (Run 2) in JSON format.

#### V0 reconstruction (`K0Analyzer`)
V0s are not taken from the standard MiniAOD V0 collections. Instead they are rebuilt from tracks:
1. Tracks from `packedPFCandidates` and `lostTracks` with track details are preselected.
2. Every pair of preselected tracks is fitted to a common vertex (`VZeroFitter`).
3. The pair's invariant mass is computed under the K<sup>0</sup><sub>S</sub>, Λ and Λ̄ hypotheses. A candidate is kept if it falls inside one of the mass windows (±70 MeV around m(K<sup>0</sup><sub>S</sub>), ±50 MeV around m(Λ)). The closest hypothesis sets `_V0Type`:
   `1` = K<sup>0</sup><sub>S</sub>, `2` = Λ, `3` = Λ̄.
4. The stored quantities include the vertex position and its uncertainty, displacement from the primary vertex and beamspot (`_V0RPV`, `_V0RBS` and their significances), DCA and point of closest approach, vertex χ², and the kinematics, hits, impact parameters and isolation of the daughter tracks.

### python
`ConfFile_cfg.py` is the `cmsRun` configuration. It configures the `V0Analyzer` module and its input tags, and takes these command-line options:

| Option      | Default            | Meaning |
|-------------|--------------------|---------|
| `isData`    | `False`            | Data or MC |
| `campaign`  | `2022`             | Data-taking/MC campaign (`2022`, `2022EE`, `2023`, `2023BPix`, `2024`, `2025`, `2026`) |
| `era`       | `Run2022B`         | Era or sample name (key in `data/Run<campaign>.json`) |
| `dataset`   | `DoubleMuon`       | Primary dataset |
| `globaltag` | `130X_dataRun3_v2` | Global tag |
| `inputFiles`| —                  | Input MiniAOD files |

The output file is named `output_<dataset>_<campaign>_<era>.root`.

Other module parameters set in the config are `skim` (default `noskim`) and the `store*` flags (`storeLheParticles`, `storeGenParticles`, `storeParticleLevel`, `storeJecSourcesAll`, `storeJecSourcesGrouped`, `storeAllTauID`, `storeJetSubstructure`).

### data
One JSON file per campaign (`Run2022.json`, `Run2022EE.json`, `Run2023.json`, `Run2024.json`, `Run2025.json`, `Run2026.json`). Each contains:
* `lumijson`: the golden JSON used as lumi mask for data,
* `eras`: for each era or MC sample, its DAS `dataset` and `globalTag`.

Files in `data/` are shipped with crab jobs.

### test
The CRAB submission script and Catch2 unit-test stubs, see [test/README.md](test/README.md).

## Event flow
For each event, `V0Analyzer::analyze` does the following:
1. `LheAnalyzer` (MC only): weights, filling `nVertices`.
2. `ParticleLevelAnalyzer` (MC, `storeParticleLevel`): if the event passes particle level, the reco-level skim below is **not** applied.
3. Events without a primary vertex are dropped.
4. `LeptonAnalyzer`: the event is dropped if it fails the lepton part of the `skim` (`dilep`, `trilep`, `singlelep`, `singletau`, `ttg`, `FR`).
5. `JetAnalyzer`: the event is dropped if it fails the jet part of the `skim` (`singlejet`, `FR`, `fourTopBase`).
6. `GenAnalyzer` (MC only), `TriggerAnalyzer`.
7. `K0Analyzer`: V0 reconstruction relative to the leading primary vertex.

With the default `skim = "noskim"`, no lepton or jet requirement is applied.

## Quick start
```bash
cd $CMSSW_BASE/src
scram b -j 8
cmsRun V0Analysis/V0Analyzer/python/ConfFile_cfg.py \
    isData=True campaign=2022EE era=Run2022G dataset=Muon \
    globaltag=<global tag from data/Run2022EE.json> \
    inputFiles=root://cms-xrd-global.cern.ch//store/data/...
```
