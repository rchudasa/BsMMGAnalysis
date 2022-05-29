import FWCore.ParameterSet.Config as cms

from BsMMGAnalysis.BsToMuMuGammaNTuplizer.Ntuplizer_cff import *

process = getDefaultWorkflow()

process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
process.MessageLogger.cerr.FwkReport.reportEvery = 50

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.options.numberOfStreams =0 ;
process.options.numberOfThreads =1 ;
process.MessageLogger.cerr.FwkReport.reportEvery = 50



process.TFileService.fileName = "dev_ntuples.root"

process=customizedProcessForData(process)

process.source.fileNames =cms.untracked.vstring(
   'file:/grid_mnt/t3storage3/athachay/bs2mumug/store/data/bParking/AOD/Run2018A/ParkingBPH3/AOD/20Jun2021_UL2018-v1/2510001/17A231BA-4B20-A648-A06D-B0E2A253ED85.root',
   )
