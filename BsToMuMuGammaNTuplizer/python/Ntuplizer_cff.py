import trigDetails

import FWCore.ParameterSet.Config as cms

globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
options = cms.untracked.PSet( numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),

   )
source = cms.Source("PoolSource",
     duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
   '/store/mc/RunIISummer19UL18RECO/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1_ext1-v2/10000/5FE97F9B-083C-9241-A7F0-6CFF243A6956.root'
   )
)

TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root")
 )


decayfilter = cms.EDFilter("GenDecayKineFilter",
    SimGenParticle = cms.InputTag("genParticles"),
    DaughterIDs = cms.untracked.vint32(13, -13, 22),
    MaxEta = cms.untracked.vdouble(2.5, 2.5, 9999.0),
    MinEta = cms.untracked.vdouble(-2.5, -2.5, -9999.0),
    MinPt = cms.untracked.vdouble(3.5, 3.5, -99.0),
    NumberDaughters = cms.untracked.int32(3),
    ParticleID = cms.untracked.int32(531),
    verbose = cms.untracked.int32(20)
)

triggersOfInterest=trigDetails.triggersOfInterest;
triggersOfInterest.extend(trigDetails.bParkTriglist)

triggerFiltesrOfInterest=trigDetails.bParkTrigFilterModules

Ntuples = cms.EDAnalyzer("BsToMuMuGammaNTuplizer",
	isMC = cms.bool(False),
	isRECO = cms.bool(True),
	isMiniAOD = cms.bool(False),
    Run2_2018               = cms.bool(True),
    
    doBeamSpot = cms.bool(False),
	beamSpot = cms.InputTag("offlineBeamSpot"),
    
    doPrimaryVetrices = cms.bool(False),
	vertices = cms.InputTag("offlinePrimaryVertices"),
    
    doMuons                 = cms.bool(True),
	muons   =cms.InputTag("muons"),
	muon_EtaMax      	    = cms.untracked.double(2.5),	 
	muon_dcaMAX 		    = cms.untracked.double(1e3),	
	muon_minPt  		    = cms.untracked.double(3.0),	
	muon_zIPMax 		    = cms.untracked.double(1e4),	
	muon_rIPMax 		    = cms.untracked.double(1e4),	
	
    doDimuons               = cms.bool(True),  
    diMuonCharge            = cms.untracked.bool(True),
    maxTwoTrackDOCA         = cms.untracked.double(0.1),
    dimuon_minPt		    = cms.untracked.double(0.0),
	dimuon_minInvMass 	    = cms.untracked.double(-1e3),	
	dimuon_maxInvMass 	    = cms.untracked.double(1e5),	
	dimuon_minVtxCL   	    = cms.untracked.double(0.0),	
	dimuon_maxLStoBS  	    = cms.untracked.double(1e5),	
	dimuon_maxDCAMuMu 	    = cms.untracked.double(1e5),	
	dimuon_maxCosAlphaToBS 	= cms.untracked.double(1e5),	
    
    genParticles       = cms.InputTag("genParticles"),
    doGenParticles     = cms.bool(False),
	doBsToMuMuGamma = cms.bool(False),
    doFlatPt           = cms.bool(True),
   	
    doSuperClusters       = cms.bool(True),
    
    MustacheSCBarrelSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    MustacheSCEndcapSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
	GsfElectronSrc     = cms.InputTag("gedGsfElectrons"),
 	hbheRechitCollection  = cms.InputTag("reducedHcalRecHits","hbhereco","RECO"),
 	ebRechitCollection    = cms.InputTag("reducedEcalRecHitsEB","","RECO"),
   	eeRechitCollection    = cms.InputTag("reducedEcalRecHitsEE","","RECO"),
    pfRechitCollection    = cms.InputTag("particleFlowRecHitECAL","","RECO"),
   	
    doPhotons          = cms.bool(True),
    gedPhotonSrc       = cms.untracked.InputTag("gedPhotons"),
    
    doPFPhotons        = cms.bool(True),
    pfPhotonSrc        = cms.untracked.InputTag("particleFlow"),
	PFPhoton_minPt     = cms.untracked.double(2.5),	
    
    doMuMuGamma        = cms.bool(True),
    doJPsiGamma        = cms.bool(True),
    minJPsiGammaMass   = cms.double(3.00),
    maxJPsiGammaMass   = cms.double(8.00),
    minBsToMuMuGammaMass = cms.double(3.50),    
    maxBsToMuMuGammaMass = cms.double(6.00),    
    
    doMuMuK            = cms.bool(True), 
    minJPsiMass        = cms.double(2.90),
    maxJPsiMass        = cms.double(3.30),
    ptMinKaon          = cms.double(1.00),
    etaMaxKaon         = cms.double(2.4),
    minBKmmMass        = cms.double(4.5),
    maxBKmmMass        = cms.double(6.0),
    doJPsiK            = cms.bool(True), 
    doPsi2SK           = cms.bool(True), 

    doMuMuKK            = cms.bool(False), 
    
    doParticleFlow     = cms.bool(False),
    particlFlowSrc     = cms.InputTag("particleFlow"),
    
    doGeneralTracks    = cms.bool(False),
    generalTrackSrc    = cms.InputTag("generalTracks"),
    
    doHCALClusters     = cms.bool(False),
    hcalClusterSrc     = cms.InputTag("particleFlowClusterHCAL"),
    
    doECALClusters     = cms.bool(False),
    ecalClusterSrc     = cms.InputTag("particleFlowClusterECAL"),
    
    doHLT              = cms.bool(True),
	TriggerNames = triggersOfInterest,
    TriggerFilters = triggerFiltesrOfInterest,
	HLTResult = cms.InputTag("TriggerResults","","HLT"),
	TriggerEvent = cms.InputTag("hltTriggerSummaryAOD"),
    verbose  = cms.bool(False),
    
    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23)   #nbits for float compression (<=23)
)


def getDefaultWorkflow( testFileName = ''  ):
    
    if testFileName=='':
        testFileName='root://se01.indiacms.res.in//store/data/Run2018A/ParkingBPH3/AOD/20Jun2021_UL2018-v1/2510000/698BC294-E1AF-C34A-94DD-2505AE79A76D.root'
    
    process = cms.Process('BMMGNtuple')
    process.load('FWCore.MessageService.MessageLogger_cfi')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
    process.load('Configuration.Geometry.GeometryIdeal_cff')
    process.load('Configuration.StandardSequences.MagneticField_cff')
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
    # import of standard configurations
    process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
    process.load('Configuration.StandardSequences.Services_cff')
    # process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
    process.load('Configuration.StandardSequences.EndOfProcess_cff')




    process.MessageLogger.cerr.FwkReport.reportEvery = 50
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
    process.options = options
    process.source = source
    process.TFileService = TFileService
    process.Ntuples = Ntuples


    process.p = cms.Path(process.Ntuples)
    
    return process

def customizedProcessForMC(process=None, minPtOfGen=2.0):
    if process==None:
        process=getDefaultWorkflow()
    process.Ntuples.doGenParticles = False
    process.Ntuples.isMC = True
    
    return process

def customizedProcessForData(process=None, minPtOfGen=2.0):
    if process==None:
        process=getDefaultWorkflow()
    return process
    
def customizedProcessForRECO(process=None, minPtOfGen=2.0):
    if process==None:
        process=getDefaultWorkflow()
    process.isRECO= True
