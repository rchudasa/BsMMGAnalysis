import FWCore.ParameterSet.Config as cms

from BsMMGAnalysis.BsToMuMuGammaNTuplizer.Ntuplizer_cff import *

process = getDefaultWorkflow()

process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
process.MessageLogger.cerr.FwkReport.reportEvery = 50

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options.numberOfStreams =0 ;
process.options.numberOfThreads =1 ;
process.MessageLogger.cerr.FwkReport.reportEvery = 50



process.TFileService.fileName = "dev_ntuples.root"


process.source.fileNames =cms.untracked.vstring(
#    'file:/grid_mnt/t3storage3/athachay/bs2mumug/store/mc/2018UL/sig/RECO/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_RECO.root'
#     'file:/grid_mnt/t3storage3/athachay/bs2mumug/store/mc/2018UL/bkg/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/00000/13B83F90-3542-484C-9B6D-1C992C96FB1A.root',
    'file:/grid_mnt/t3storage3/athachay/bs2mumug/store/data/bParking/AOD/Run2018A/ParkingBPH3/AOD/20Jun2021_UL2018-v1/2510001/17A231BA-4B20-A648-A06D-B0E2A253ED85.root',
#    'root://se01.indiacms.res.in:1094//cms/store/user/athachay/BsToMuMuGamma/Data/BsToJPsiGamma/RunIISummer20UL18/106X_upgrade2018_realistic_v11_L1v1/BsToJPsiGamma_MuGFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ui2_condor_reco/bsToJPsiGamma_run2UL2018_143.root',
   )

## If the process needs to be done on MC or MC
process=customizedProcessForData(process,addJson=True)
#process=customizedProcessForData(process,addJson=False)
#process=customizedProcessForMC(process)


## Switch on Individual branches 
switchOnMuMuK(process,ptMinKaon=0.98,etaMaxKaon=2.5,minBKmmMass=4.5,maxBKmmMass=6.5)
switchOnJPsiK(process,ptMinKaon=0.98,etaMaxKaon=2.5,minBKmmMass=4.5,maxBKmmMass=6.5)
switchOnPsi2SK(process,ptMinKaon=0.98,etaMaxKaon=2.5,minBKmmMass=4.5,maxBKmmMass=6.5)
switchOnPFPhotons(process)
switchOnJPsiGamma(process)
switchOnMuMuGamma(process)
switchOnHLT(process,True)
