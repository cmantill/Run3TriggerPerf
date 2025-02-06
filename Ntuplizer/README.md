# Setup
```bash
cmsrel CMSSW 13_0_6
cd $CMSSW_BASE/src
cmsenv
git clone git@github.com:cmantill/Run3TriggerPerf.git
cd Run3TriggerPerf/Ntuplizer/
scram b -j 8

```
# Local Run
for Double Electron
```bash
cmsRun nanoanalyzercrab_cfg_23_DoubleEle_for_singleEleEff.py
```
The file to edit for Double Electron is plugins/NanoAnalyzer.cc. Recompile with `scram b -j 8` when changes are made.

The output will be a file with a TTree (by default named `JPsi_ElePlusJet_controlTrigger.root`).

# CRAB job

## Datasets
```
datasets = {
    "2023_DoubleElectron_B": [
        "/ParkingDoubleElectronLowMass/Run2023B-22Sep2023-v1/MINIAOD",
    ],
    "2023_DoubleElectron_C": [
        "/ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v1-v1/MINIAOD",
        "/ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v2-v1/MINIAOD",
        "/ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v3-v1/MINIAOD",
        "/ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v4-v1/MINIAOD",
    ],
    "2023_DoubleElectron_D": [
        "/ParkingDoubleElectronLowMass/Run2023D-22Sep2023_v1-v1/MINIAOD",
        "/ParkingDoubleElectronLowMass/Run2023D-22Sep2023_v2-v1/MINIAOD",
    ],
}
```

To run crab jobs just submit crab.py file using:
```bash
crab submit crab.py
```

Change the following:
```python
config.Site.storageSite # this needs to be changed depending on the site where it is stored.
config.Data.outLFNDirBase # location where it is stored ( directory name) CHANGE USERNAME
```

Crab jobs must be submitted for each dataset.