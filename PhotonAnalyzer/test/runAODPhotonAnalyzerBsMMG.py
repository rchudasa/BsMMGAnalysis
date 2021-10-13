import FWCore.ParameterSet.Config as cms

process = cms.Process('Demo')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/r/rchudasa/private/bsmumu/Run2_analysis/CMSSW_10_6_20/src/BsMMGAnalysis/PhotonAnalyzer/test/AEFD418A-0A8D-414C-A8AC-86EE20287BDF.root'       
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("photon_BsToMuMuGamma_Run2018_default_reco.root"))


process.demo = cms.EDAnalyzer("PhotonAnalyzer",
    doGenParticles     = cms.bool(True),
    doPhotons          = cms.bool(True),
    doPFPhotons        = cms.bool(True),
    doSuperClusters    = cms.bool(True),
    genParticleSrc     = cms.untracked.InputTag("genParticles"),
    gedPhotonSrc       = cms.untracked.InputTag("gedPhotons"),
    pfPhotonSrc        = cms.untracked.InputTag("particleFlow"),
    MustacheSCBarrelSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    MustacheSCEndcapSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower")
)

process.p = cms.Path(process.demo)
