import FWCore.ParameterSet.Config as cms

process = cms.Process('Demo')

process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet( numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),

   )



process.source = cms.Source("PoolSource",
     duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
   #'file:/eos/user/a/athachay/workarea/data/BsToMuMuGamma/RunIIAutumn18DRPremix/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-evtgen-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/606765BD-9BB4-9741-925C-A0C69B933039.root',      
 #  'file:/afs/cern.ch/work/a/athachay/public/BsToMuMuGamma/RunIIAutumn18DRPremix/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-evtgen-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/606765BD-9BB4-9741-925C-A0C69B933039.root',      
   #'file:DoublePhotonGun/DoublePhoton0To40FlatPtAODSIM_HI_Reco_1.root',      
   'file:/afs/cern.ch/work/r/rchudasa/private/bsmumu/Run2_analysis/CMSSW_10_6_20/src/BsMMGAnalysis/PhotonAnalyzer/test/AEFD418A-0A8D-414C-A8AC-86EE20287BDF.root'      
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("muonNtuplizer.root")
)


process.decayfilter = cms.EDFilter("GenDecayKineFilter",
    SimGenParticle = cms.InputTag("genParticles"),
    DaughterIDs = cms.untracked.vint32(13, -13, 22),
    MaxEta = cms.untracked.vdouble(2.5, 2.5, 9999.0),
    MinEta = cms.untracked.vdouble(-2.5, -2.5, -9999.0),
    MinPt = cms.untracked.vdouble(3.5, 3.5, -99.0),
    NumberDaughters = cms.untracked.int32(3),
    ParticleID = cms.untracked.int32(531),
    verbose = cms.untracked.int32(20)
)

process.Ntuples = cms.EDAnalyzer("BsToMuMuGammaNTuplizer",
	muons   =cms.InputTag("muons"),
	beamSpot = cms.InputTag("offlineBeamSpot"),
	vertices = cms.InputTag("offlinePrimaryVertices"),
    	MustacheSCBarrelSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    	MustacheSCEndcapSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
	muon_EtaMax      	= cms.untracked.double(1e3),	 
	muon_dcaMAX 		= cms.untracked.double(1e3),	
	muon_minPt  		= cms.untracked.double(1.0),	
	muon_zIPMax 		= cms.untracked.double(1e4),	
	muon_rIPMax 		= cms.untracked.double(1e4),	
	dimuon_minPt		= cms.untracked.double(0.0),
	dimuon_minInvMass 	= cms.untracked.double(-1e3),	
	dimuon_maxInvMass 	= cms.untracked.double(1e5),	
	dimuon_minVtxCL   	= cms.untracked.double(0.0),	
	dimuon_maxLStoBS  	= cms.untracked.double(1e5),	
	dimuon_maxDCAMuMu 	= cms.untracked.double(1e5),	
	dimuon_maxCosAlphaToBS 	= cms.untracked.double(1e5),	
    	doGenParticles     = cms.bool(True),
   	doMuons            = cms.bool(True),
   	doPhotons          = cms.bool(True),
    	doPFPhotons        = cms.bool(True),
    	doSuperClusters    = cms.bool(True),
    	genParticles       = cms.InputTag("genParticles"),
    	gedPhotonSrc       = cms.untracked.InputTag("gedPhotons"),
    	pfPhotonSrc        = cms.untracked.InputTag("particleFlow"),
	verbose  = cms.bool(False),
	doBsToMuMuGamma = cms.bool(True),
	isMC = cms.bool(True),
)


process.p = cms.Path(process.decayfilter*process.Ntuples)
