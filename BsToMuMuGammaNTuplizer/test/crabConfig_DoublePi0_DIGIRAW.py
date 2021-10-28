from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'DoublePi0FlatPt1To20_DIGIRAW_Run2018_correctEta'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'DIGIToRaw_noPU2018_cfg.py'
#config.Data.outputPrimaryDataset = 'SinglePhotonFlatPt1To20_DIGIRAW_Run2018'
config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.inputDataset ='/DoublePi0FlatPt1To20_GENSIM_Run2018_correctEta/rchudasa-crab_DoublePi0FlatPt1To20_GENSIM_Run2018_correctEta-fd99b7ab52e8e7d9d9b1e36213873f36/USER'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 

#config.Data.outLFNDirBase = '/store/group/phys_diffraction/lbyl_2018/mc_cep/ntuples'
config.Data.outLFNDirBase = '/store/user/rchudasa/BsMMG_2018UL'
config.Data.publication = True 
config.Site.storageSite = 'T2_IN_TIFR'
