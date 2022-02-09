from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'ntuples_QCD-Pt15to30-Run2MVA'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplize_AOD.py'
config.JobType.maxMemoryMB = 4000
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.inputDataset ='/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/RunIISummer19UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM'

config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 100
#config.Data.outLFNDirBase = '/store/group/phys_diffraction/lbyl_2018/mc_lbl/ntuples'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/lowPT_photonReco'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ntuples_QCD-Pt15to30-Run2MVA'
config.Site.storageSite = 'T2_CH_CERN'
