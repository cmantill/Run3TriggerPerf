from CRABClient.UserUtilities import config

config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzercrab_cfg_23_DoubleEle_for_singleEleEff.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 1

config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2

config.General.workArea = 'crab_doubleele_2023/'

dataset_name = '/ParkingDoubleElectronLowMass/Run2023C-22Sep2023_v1-v1/MINIAOD' # change here dataset name
name='DoubleElectron-Run2023Cv1' # change here crab_job name
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = dataset_name
config.Data.outputDatasetTag = 'ParkingDoubleEleLowMass_2023'
config.Data.outLFNDirBase = '/store/user/cmantill/run3triggerperf/%s' % (config.General.workArea)

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

