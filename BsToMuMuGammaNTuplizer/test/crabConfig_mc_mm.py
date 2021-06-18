from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'ntuples_mm_correctGenInfov2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplize_AOD.py'
config.JobType.maxMemoryMB = 4000

config.Data.inputDataset ='/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DRPremix-102X_upgrade2018_realistic_v15-v1/AODSIM'

config.JobType.allowUndistributedCMSSW = True
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/lowPT_photonReco'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ntuples_mm_correctGenInfov2'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_CH_CERN']
