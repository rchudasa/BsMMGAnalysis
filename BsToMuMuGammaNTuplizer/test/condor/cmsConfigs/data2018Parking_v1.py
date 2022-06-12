import FWCore.ParameterSet.Config as cms

from BsMMGAnalysis.BsToMuMuGammaNTuplizer.Ntuplizer_cff import *

process = getDefaultWorkflow()

process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
process.MessageLogger.cerr.FwkReport.reportEvery = 50

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.options.numberOfStreams =0 ;
process.options.numberOfThreads =1 ;
process.MessageLogger.cerr.FwkReport.reportEvery = 50



process.TFileService.fileName = "dev_ntuples.root"

process=customizedProcessForData(process)

process.source.fileNames =cms.untracked.vstring(
#    'file:/grid_mnt/t3storage3/athachay/bs2mumug/store/mc/2018UL/sig/RECO/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_RECO.root'
#     'file:/grid_mnt/t3storage3/athachay/bs2mumug/store/mc/2018UL/bkg/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/00000/13B83F90-3542-484C-9B6D-1C992C96FB1A.root',
    'file:/grid_mnt/t3storage3/athachay/bs2mumug/store/data/bParking/AOD/Run2018A/ParkingBPH3/AOD/20Jun2021_UL2018-v1/2510001/17A231BA-4B20-A648-A06D-B0E2A253ED85.root',
   )


## Switch on Individual branches 
switchOnMuMuK(process,ptMinKaon=0.98,etaMaxKaon=2.5,minBKmmMass=4.5,maxBKmmMass=6.5)
switchOnJPsiK(process,ptMinKaon=0.98,etaMaxKaon=2.5,minBKmmMass=4.5,maxBKmmMass=6.5)
switchOnPsi2SK(process,ptMinKaon=0.98,etaMaxKaon=2.5,minBKmmMass=4.5,maxBKmmMass=6.5)
switchOnJPsiGamma(process)
switchOnMuMuGamma(process)
switchOnPFPhotons(process)


f=open("flist.cfg")
txt=f.readlines()
f.close()
maxEvents=int(txt.pop(0)[:-1])
fOutName=txt.pop(0)[:-1]
process.source.fileNames=cms.untracked.vstring()
for l in txt:
    process.source.fileNames.append(l[:-1])
process.maxEvents.input = maxEvents
process.TFileService.fileName = fOutName

