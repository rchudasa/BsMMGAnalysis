
## Edited By Raman Khurana
##
## CRAB documentation : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrab
##
## CRAB 3 parameters : https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters
##
## Once you are happy with this file, please run
## crab submit

## In CRAB3 the configuration file is in Python language. It consists of creating a Configuration object imported from the WMCore library: 

from WMCore.Configuration import Configuration
config = Configuration()

##  Once the Configuration object is created, it is possible to add new sections into it with corresponding parameters
config.section_("General")
config.General.requestName = 'charmonium2018A_ntuplizer'
config.General.workArea = 'crab_pickevents_20210803_183029'


config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PickMerge_cfg.py'
config.JobType.pyCfgParams = ['eventsToProcess_load=pickevents_runEvents.txt', 'outputFile=charmonium2018A_BMMGNtuple.root']

config.section_("Data")
config.Data.inputDataset = '/Charmonium/Run2018A-12Nov2019_UL2018_rsb-v1/AOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'pickevents.json'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/athachay/BsToMuMuGamma/Run2Studies/BsToMuMuGammaNtuples'

config.section_("Site")
## Change site name accordingly
config.Site.storageSite = "T2_IN_TIFR"

