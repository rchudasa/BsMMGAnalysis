import FWCore.ParameterSet.Config as cms

process = cms.Process('Demo')

process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
process.load("Geometry.HcalEventSetup.CaloTowerTopology_cfi")
process.load("Configuration.Geometry.GeometryExtended2017_cff")
process.load("Configuration.Geometry.GeometryExtended2017Reco_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("Geometry.HcalEventSetup.hcalTopologyIdeal_cfi")

process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
#process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v15')

process.MessageLogger.cerr.FwkReport.reportEvery = 50

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

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
   'file:59B967C8-4A74-ED4F-9A8D-15A05983A34A.root',
   )
)
process.TFileService = cms.Service("TFileService",
    #fileName = cms.string("muonNtuplizer.root")
    #fileName = cms.string("flatPtPhoton_ntuple_pfIso.root")
    #fileName = cms.string("flatPtElectron_ntuple_pfIso.root")
    #fileName = cms.string("flatPtPi0_ntuple_pfIso.root")
    #fileName = cms.string("emptyBX_ntuple_pfIso.root")
    fileName = cms.string("qcd_EmEnriched_ntuples.root")
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


process.dimuonFilter = cms.EDFilter("DimuonMassFilter",
    muons = cms.InputTag("muons"),
    MinPt = cms.untracked.double(3.5),
    MaxAbsEta = cms.untracked.double(2.5),
    MinDimuonMass = cms.untracked.double(0.5),
    MaxDimuonMass = cms.untracked.double(8.5)
)

process.Ntuples = cms.EDAnalyzer("BsToMuMuGammaNTuplizer",
	muons   =cms.InputTag("muons"),
	beamSpot = cms.InputTag("offlineBeamSpot"),
    doBeamSpot = cms.bool(False),
    doPrimaryVetrices = cms.bool(True),
	vertices = cms.InputTag("offlinePrimaryVertices"),
    	MustacheSCBarrelSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    	MustacheSCEndcapSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
	GsfElectronSrc     = cms.InputTag("gedGsfElectrons"),
    Run2_2018          = cms.bool(True),
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
	isMC = cms.bool(False),
    doGenParticles     = cms.bool(False),
	doBsToMuMuGamma = cms.bool(False),
    doFlatPt           = cms.bool(False),
    genParticles       = cms.InputTag("genParticles"),
   	doMuons            = cms.bool(False),
	PFPhoton_minPt     = cms.untracked.double(0.0),	
    doSuperClusters    = cms.bool(True),
 	hbheRechitCollection            = cms.InputTag("reducedHcalRecHits","hbhereco","RECO"),
 	ebRechitCollection              = cms.InputTag("reducedEcalRecHitsEB","","RECO"),
   	eeRechitCollection              = cms.InputTag("reducedEcalRecHitsEE","","RECO"),
    pfRechitCollection              = cms.InputTag("particleFlowRecHitECAL","","RECO"),
   	doPhotons          = cms.bool(False),
    gedPhotonSrc       = cms.untracked.InputTag("gedPhotons"),
    doPFPhotons        = cms.bool(False),
    pfPhotonSrc        = cms.untracked.InputTag("particleFlow"),
    doParticleFlow     = cms.bool(False),
    particlFlowSrc     = cms.InputTag("particleFlow"),
    doGeneralTracks    = cms.bool(False),
    generalTrackSrc    = cms.InputTag("generalTracks"),
    doHCALClusters     = cms.bool(False),
    hcalClusterSrc     = cms.InputTag("particleFlowClusterHCAL"),
    doECALClusters     = cms.bool(False),
    ecalClusterSrc     = cms.InputTag("particleFlowClusterECAL"),
    doHLT              = cms.bool(False),
	TriggerNames = cms.vstring("HLT_DoubleMu4_3_Bs_v14"),
	HLTResult = cms.InputTag("TriggerResults","","HLT"),
	verbose  = cms.bool(False),
    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23)   #nbits for float compression (<=23)
)


#process.p = cms.Path(process.decayfilter*process.Ntuples)
process.p = cms.Path(process.dimuonFilter+process.Ntuples)
