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
# ParkingBPH 2021 UL
#process.GlobalTag.globaltag = '106X_dataRun2_v35'
# 2021 data
#process.GlobalTag.globaltag = '120X_dataRun3_Prompt_v2'
#process.GlobalTag.globaltag = '121X_dataRun3_v13'
# 2022 data
#process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v8'
#process.GlobalTag.globaltag = '123X_dataRun3_Express_v10'
#process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v12'
# 2023 data
process.GlobalTag.globaltag = '130X_dataRun3_Express_v2'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#####################################  JSON FILE #################################
# Change the directory and JSON file accordingly
# Only uncomment if you run in Data
# ParkingBPH 2021 UL
#goodJSON = '/nfs/dust/cms/user/yangq2/goodJson/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
# pilot 2021
#goodJSON = './BT21GOOD_withTKon.json'
# 2022A, data taken up to June 3
#goodJSON = './Collisions22AGOOD_withALLON.json'


##################################################################################

# Get the luminosity list from the JSON
# Only uncomment if you run in Data
#myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

##################################################################################

# Load jet correction services for all jet algoritms
#process.load("JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff") # 2022 data
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
#

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

inputFiles= [
    # '/store/data/Run2022D/EGamma/MINIAOD/PromptReco-v2/000/357/899/00000/990a1665-d685-4623-8ddb-803c01962243.root', # 2022
    '/store/data/Run2023B/ParkingDoubleElectronLowMass/MINIAOD/PromptReco-v1/000/366/469/00000/c917911f-2406-4cca-91fe-3deb59cf9067.root' # 2023
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
