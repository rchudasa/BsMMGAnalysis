import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

# add a list of strings for events to process
options.register ('eventsToProcess',
				  '',
				  VarParsing.multiplicity.list,
				  VarParsing.varType.string,
				  "Events to process")
options.register ('maxSize',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "Maximum (suggested) file size (in Kb)")
options.parseArguments()

process = cms.Process("PickEvent")
process.source = cms.Source ("PoolSource",
	  fileNames = cms.untracked.vstring(options.inputFiles),
)

if options.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange (options.eventsToProcess)

process.maxEvents = cms.untracked.PSet(
	    input = cms.untracked.int32 (options.maxEvents)
)

process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = cms.string('106X_dataRun2_v24')

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet( numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
   )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("Charmonium_Run2018A12Nov2019_UL2018_rsbv1.root")
)

process.Ntuples = cms.EDAnalyzer("BsToMuMuGammaNTuplizer",
	muons   =cms.InputTag("muons"),
	beamSpot = cms.InputTag("offlineBeamSpot"),
	vertices = cms.InputTag("offlinePrimaryVertices"),
    	MustacheSCBarrelSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    	MustacheSCEndcapSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
	GsfElectronSrc         = cms.InputTag("gedGsfElectrons"),
	muon_EtaMax      	   = cms.untracked.double(1e3),	 
	muon_dcaMAX 		   = cms.untracked.double(1e3),	
	muon_minPt  		   = cms.untracked.double(1.0),	
	muon_zIPMax 		   = cms.untracked.double(1e4),	
	muon_rIPMax 		   = cms.untracked.double(1e4),	
	dimuon_minPt		   = cms.untracked.double(0.0),
	dimuon_minInvMass 	   = cms.untracked.double(-1e3),	
	dimuon_maxInvMass 	   = cms.untracked.double(1e5),	
	dimuon_minVtxCL   	   = cms.untracked.double(0.0),	
	dimuon_maxLStoBS  	   = cms.untracked.double(1e5),	
	dimuon_maxDCAMuMu 	   = cms.untracked.double(1e5),	
	dimuon_maxCosAlphaToBS = cms.untracked.double(1e5),	
	verbose                = cms.bool(False),
	doBsToMuMuGamma        = cms.bool(False),
    doHLT                  = cms.bool(True),
	isMC                   = cms.bool(False),
    doGenParticles         = cms.bool(False),
   	doMuons                = cms.bool(False),
   	doPhotons              = cms.bool(True),
    doPFPhotons            = cms.bool(True),
	PFPhoton_minPt         = cms.untracked.double(2.0),	
    doSuperClusters        = cms.bool(True),
    genParticles           = cms.InputTag("genParticles"),
    gedPhotonSrc           = cms.untracked.InputTag("gedPhotons"),
    pfPhotonSrc            = cms.untracked.InputTag("particleFlow"),
	TriggerNames           = cms.vstring("HLT_DoubleMu4_3_Bs_v14",
	                           "HLT_DoubleMu4_3_Jpsi_v2",
	            			   "HLT_DoubleMu4_JpsiTrk_Displaced_v15",
				               "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15",
				               "HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8"),
	HLTResult              = cms.InputTag("TriggerResults","","HLT"),
 	ebRechitCollection     = cms.InputTag("reducedEcalRecHitsEB","","RECO"),
   	eeRechitCollection     = cms.InputTag("reducedEcalRecHitsEE","","RECO"),
    pfRechitCollection     = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    doCompression          = cms.bool(True),  #do the compression of floats
    nBits                  = cms.int32(23),   #nbits for float compression (<=23)
  )


process.p = cms.Path(process.Ntuples)
