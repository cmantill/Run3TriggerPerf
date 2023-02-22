from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'nanoanalyzercrab_cfg_22_DoubleEle.py'
config.JobType.psetName = 'nanoanalyzercrab_cfg_22_mc_DoubleEle.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 1

config.Data.inputDataset = ''
#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2

<<<<<<< HEAD
#config.Data.outLFNDirBase = '/store/user/jodedra/Run3/' #% (getUsernameFromSiteDB())
#config.Data.publication = True
#config.Data.outputDatasetTag = 'winter21'
#config.Site.ignoreGlobalBlacklist = True
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_CH_*']
#config.Site.blacklist = ['T2_US_Purdue']


#config.General.workArea = 'crab_doubleele_MC/crab_doubleele_MC'
=======
config.General.workArea = 'crab_elejet_lastversion/crab_elejet_lastversion'
>>>>>>> 681a89e5aabc6b1545159594ba94d447d850dfe5

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

<<<<<<< HEAD
# name='/EGamma/Run2022F-PromptReco-v1/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/EGamma/Run2022F-PromptReco-v1/MINIAOD'
# config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022/%s' % (config.General.workArea)

# name='/BuTOjpsiKEE2022/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/BuTOjpsiKEE20220831fiftyMbettersplitting/jodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75/USER'
# config.Data.outputDatasetTag = 'BuTOjpsiKEE2022_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/cquarant/BuTOjpsiKEE2022/%s' % (config.General.workArea)

# name='/BuTOKEE2022/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/BuTOKEE20220826bettersplitting/jodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75/USER'
# config.Data.outputDatasetTag = 'BuTOKEE2022_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/cquarant/BuTOKEE2022/%s' % (config.General.workArea)
=======
name='/EGamma/Run2022G-PromptReco-v1/MINIAOD/DoubleEleTest'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/EGamma/Run2022G-PromptReco-v1/MINIAOD'
config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleEleTest'
config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)
>>>>>>> 681a89e5aabc6b1545159594ba94d447d850dfe5

# name='/BuTOpsi2sKEE2022/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/BuTOpsi2sKEE20220831fiftyMbettersplitting/jodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75/USER'
# config.Data.outputDatasetTag = 'BuTOpsi2sKEE2022_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/cquarant/BuTOjpsiKEE2022/%s' % (config.General.workArea)

config.General.workArea = 'crab_doubleele_MC/crab_doubleele_MC'

name='/BuTOjpsiKEE2022_realisticPU/MINIAOD/DoubleEleTest'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/BuTOjpsiKEE20221103FIFTYM/jodedra-SUMMER22_MINIAOD-8da3c779fd247f42d8411cf96207b146/USER'
config.Data.outputDatasetTag = 'BuTOjpsiKEE2022_realisticPU_DoubleEleTest'
config.Data.outLFNDirBase = '/store/user/cquarant/BuTOjpsiKEE2022_realisticPU/%s' % (config.General.workArea)

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'
#config.Site.storageSite = 'T2_CH_CERN'

