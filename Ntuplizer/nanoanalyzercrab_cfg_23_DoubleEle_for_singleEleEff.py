##################################################################################
# Nanoanalyzer configuration file for all years                                  #
# Use for HT Condor and VM                                                       #
# Uncomment & comment relevant lines before you run it                           #
##################################################################################

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("Nano")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

###################################### GLOBAL TAG ################################
# Change the global tag accordingly
process.GlobalTag.globaltag = '130X_mcRun3_2023_realistic_v14'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

##################################################################################

# Load jet correction services for all jet algoritms
#process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

# /ParkingDoubleElectronLowMass/Run2023B-22Sep2023-v1/MINIAOD
# /ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v1-v1/MINIAOD
# /ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v2-v1/MINIAOD
# /ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v3-v1/MINIAOD
# /ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v4-v1/MINIAOD
# /ParkingDoubleElectronLowMass/Run2023D-22Sep2023_v1-v1/MINIAOD
# /ParkingDoubleElectronLowMass/Run2023D-22Sep2023_v2-v1/MINIAOD
inputFiles= [
    '/store/data/Run2023C/ParkingDoubleElectronLowMass/MINIAOD/22Sep2023_v1-v1/2550000/0003dc3b-632b-40bc-91a1-694f8badf775.root',
]

process.source = cms.Source("PoolSource",
# for crab
	  fileNames = cms.untracked.vstring (inputFiles),
)

##################################################################################

# Process the lumi
# Only uncomment if you run in Data
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string ('JPsi_ElePlusJet_controlTrigger.root')
                                   )


# Process the analyzer
process.nano_ = cms.EDAnalyzer('NanoAnalyzerDoubleEle',
                              electrons = cms.InputTag("slimmedElectrons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              packedpfcandidates = cms.InputTag('packedPFCandidates'),
                              HLT = cms.InputTag("TriggerResults","","HLT"),
                              triggerobjects = cms.InputTag("slimmedPatTrigger"),
                              l1EG = cms.InputTag("caloStage2Digis", "EGamma"), 
)
process.p = cms.Path(process.nano_)
