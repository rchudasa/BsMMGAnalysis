from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'ntuples_defaultReco_flatPhoton_genInfo'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplize_AOD.py'
config.JobType.maxMemoryMB = 4000
config.Data.inputDBS = 'phys03'
config.Data.inputDataset ='/SinglePhotonFlatPt1To20_GENSIM_Run2018/rchudasa-crab_SinglePhotonFlatPt1To20_default_RecoAOD_Run2018-648c701074ba72c4e16eb15dae8f7a8f/USER'

config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#config.Data.outLFNDirBase = '/store/group/phys_diffraction/lbyl_2018/mc_lbl/ntuples'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/lowPT_photonReco'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ntuples_defaultReco_flatPhoton_genInfo'
config.Site.storageSite = 'T2_CH_CERN'
