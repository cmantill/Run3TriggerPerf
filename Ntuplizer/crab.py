from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'nanoanalyzercrab_cfg_22_DoubleEle.py'
config.JobType.psetName = 'nanoanalyzercrab_cfg_23_DoubleEle.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 1

config.Data.inputDataset = ''
config.Data.inputDBS = 'global' # Data
#config.Data.inputDBS = 'phys03' #MC
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2

#config.Data.outLFNDirBase = '/store/user/jodedra/Run3/' #% (getUsernameFromSiteDB())
#config.Data.publication = True
#config.Data.outputDatasetTag = 'winter21'
#config.Site.ignoreGlobalBlacklist = True
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_CH_*']
#config.Site.blacklist = ['T2_US_Purdue']


#config.General.workArea = 'crab_doubleele_MC/crab_doubleele_MC'

#name='/EGamma/Run2022C-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022C-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022C-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

#name='/EGamma/Run2022D-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022D-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

#name='/EGamma/Run2022D-PromptReco-v2/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)
##config.Data.outLFNDirBase = '/store/group/phys_bphys/crovelli/triggerRun3/EGamma2022/%s' % (config.General.workArea)   

#name='/EGamma/Run2022E-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

#name='/EGamma/Run2022F-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

# name='/BuTOpsi2sKEE2022/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/BuTOpsi2sKEE20220831fiftyMbettersplitting/jodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75/USER'
# config.Data.outputDatasetTag = 'BuTOpsi2sKEE2022_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/cquarant/BuTOjpsiKEE2022/%s' % (config.General.workArea)

# config.General.workArea = 'crab_doubleele_MC/crab_doubleele_MC'

# name='/BuTOjpsiKEE2022_realisticPU/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/BuTOjpsiKEE20221103FIFTYM/jodedra-SUMMER22_MINIAOD-8da3c779fd247f42d8411cf96207b146/USER'
# config.Data.outputDatasetTag = 'BuTOjpsiKEE2022_realisticPU_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/cquarant/BuTOjpsiKEE2022_realisticPU/%s' % (config.General.workArea)

####################  2023  ##############################

config.General.workArea = 'crab_doubleele_2023/'

name='/ParkingDoubleElectronLowMass/Run2023B-PromptReco-v1/MINIAOD'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/ParkingDoubleElectronLowMass/Run2023B-PromptReco-v1/MINIAOD'
config.Data.outputDatasetTag = 'ParkingDoubleEleLowMass_2023'
config.Data.outLFNDirBase = '/store/user/cquarant/%s' % (config.General.workArea)

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'
#config.Site.storageSite = 'T2_CH_CERN'

