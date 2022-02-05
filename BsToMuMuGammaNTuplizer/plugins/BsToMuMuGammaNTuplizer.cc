// -*- C++ -*-
//
// Package:    BsMMGAnalysis/BsToMuMuGammaNTuplizer
// Class:      BsToMuMuGammaNTuplizer
//
/**\class BsToMuMuGammaNTuplizer BsToMuMuGammaNTuplizer.cc BsMMGAnalysis/BsToMuMuGammaNTuplizer/plugins/BsToMuMuGammaNTuplizer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author: 
//         Created: 
//
//


// system include files
#include <memory>
#include <map>
#include <string>
#include "TH1.h"

#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EGHcalRecHitSelector.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/plugins/pfIsoCalculator.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "EgammaAnalysis/ElectronTools/interface/SuperClusterHelper.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/HIPhotonIsolation.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/plugins/BsToMuMuGammaNTuplizer.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/interface/Utils.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include <TTree.h>


BsToMuMuGammaNTuplizer::BsToMuMuGammaNTuplizer(const edm::ParameterSet& iConfig) :

  doGenParticles_(iConfig.getParameter<bool>("doGenParticles")),
  doMuons_(iConfig.getParameter<bool>("doMuons")),
  doPhotons_(iConfig.getParameter<bool>("doPhotons")),
  doPFPhotons_(iConfig.getParameter<bool>("doPFPhotons")),
  doSuperClusters_(iConfig.getParameter<bool>("doSuperClusters")),
  doParticleFlow(iConfig.getParameter<bool>("doParticleFlow")),
  doGeneralTracks(iConfig.getParameter<bool>("doGeneralTracks")),
  doECALClusters(iConfig.getParameter<bool>("doECALClusters")),
  doHCALClusters(iConfig.getParameter<bool>("doHCALClusters")),
  doPrimaryVetrices(iConfig.getParameter<bool>("doPrimaryVetrices")),
  doBeamSpot(iConfig.getParameter<bool>("doBeamSpot")),
  doFlatPt_(iConfig.getParameter<bool>("doFlatPt")),
  Run2_2018_(iConfig.getParameter<bool>("Run2_2018")),
  doHLT(iConfig.getParameter<bool>("doHLT")),
  energyMatrixSize_(2)
{
  
  energyMatrixSizeFull_=(2*energyMatrixSize_+1)*(2*energyMatrixSize_+1);
  Utility= new Utils();
  
  if(doMuons_) muonToken_              = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  if(doPhotons_)    gedPhotonsCollection_       = consumes<std::vector<reco::Photon>>(iConfig.getUntrackedParameter<edm::InputTag>("gedPhotonSrc"));
  
  if(doPFPhotons_ or doSuperClusters_){
    pfPhotonsCollection_        = consumes<edm::View<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("pfPhotonSrc"));
  }  
  
  if(doSuperClusters_){
    MustacheSCBarrelCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCBarrelSrc"));
    MustacheSCEndcapCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCEndcapSrc"));
    gsfElectronToken_                       = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("GsfElectronSrc"));
    hbheRechitToken_               = consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>(iConfig.getParameter<edm::InputTag>("hbheRechitCollection"));
  }
  if(doSuperClusters_ or doECALClusters) {
    ebRechitToken_                 = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
    eeRechitToken_                 = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
  }
  if(doParticleFlow)
  {
     pfCandidateCollection_        = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("particlFlowSrc"));
  }
  if(doHCALClusters)
  {
     //std::cout<<" doHCALClusters : is true \n";
     hcalClusterCollection_        = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("hcalClusterSrc"));
  }
  if(doECALClusters)
  {
     ecalClusterCollection_        = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("ecalClusterSrc"));
  }
  if(doGeneralTracks)
  {
     generalTracksCollection_        = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTrackSrc"));
  }


  if(doHLT) {
    trigTable    =iConfig.getParameter<std::vector<std::string>>("TriggerNames");
    triggerBits_ =consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTResult"));
  }


  
  beamSpotToken_  	   =consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  primaryVtxToken_       =consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  if(doGenParticles_)genParticlesCollection_          =consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  
  etaMax_muon               =  iConfig.getUntrackedParameter<double>("muon_EtaMax")        ;
  dcaMax_muon_bs            =  iConfig.getUntrackedParameter<double>("muon_dcaMAX")        ;
  pTMinMuons 		      =  iConfig.getUntrackedParameter<double>("muon_minPt");
  pTMinPFPhotons               =  iConfig.getUntrackedParameter<double>("PFPhoton_minPt");
  trackIP_zMax_muon	      =  iConfig.getUntrackedParameter<double>("muon_zIPMax")        ;
  trackIP_rMax_muon	      =  iConfig.getUntrackedParameter<double>("muon_rIPMax")        ;
  
  minDimuon_pt              =  iConfig.getUntrackedParameter<double>("dimuon_minPt")      ;
  cl_dimuon_vtx             =  iConfig.getUntrackedParameter<double>("dimuon_minVtxCL")      ;
  ls_max_dimuonBS           =  iConfig.getUntrackedParameter<double>("dimuon_maxLStoBS")       ;
  dcaMax_dimuon_mumu        =  iConfig.getUntrackedParameter<double>("dimuon_maxDCAMuMu")        ;
  cosAlphaMax_dimuonBs      =  iConfig.getUntrackedParameter<double>("dimuon_maxCosAlphaToBS") ;
  minDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_minInvMass")    ;
  maxDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_maxInvMass")    ;
  
  printMsg=iConfig.getParameter<bool>("verbose");
  isMC=iConfig.getParameter<bool>("isMC");
  doBsToMuMuGamma=iConfig.getParameter<bool>("doBsToMuMuGamma");
  nBits_                         = iConfig.getParameter<int>("nBits"); 
  doCompression_                 = iConfig.getParameter<bool>("doCompression");  
  // initialize output TTree
  
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("EventTree", "Event data");
  
  theTree->Branch("run",    &run_);
  theTree->Branch("event",  &event_);
  theTree->Branch("lumis",  &lumis_);
  theTree->Branch("isData", &isData_);
 
  if(doBeamSpot)
  {
    theTree->Branch("beamspot_x",         		 &beamspot_x_);
    theTree->Branch("beamspot_y",       	    	 &beamspot_y_);
    theTree->Branch("beamspot_z",       		     &beamspot_z_);
    theTree->Branch("beamspot_x_error",       	 &beamspot_x_error_);
    theTree->Branch("beamspot_y_error",       	 &beamspot_y_error_);
    theTree->Branch("beamspot_z_error",        	 &beamspot_z_error_);
    theTree->Branch("beamspot_covXX",         	 &beamspot_covXX);
    theTree->Branch("beamspot_covXY",         	 &beamspot_covXY);
    theTree->Branch("beamspot_covXZ",         	 &beamspot_covXZ);
    theTree->Branch("beamspot_covYY",         	 &beamspot_covYY);
    theTree->Branch("beamspot_covYZ",         	 &beamspot_covYZ);
    theTree->Branch("beamspot_covZZ",         	 &beamspot_covZZ);

    theTree->Branch("beamspot_dxdz",               &beamspot_dxdz_);
    theTree->Branch("beamspot_dydz",               &beamspot_dydz_);
    theTree->Branch("beamspot_sigmaZ",             &beamspot_sigmaZ_);
    theTree->Branch("beamspot_dxdz_error",         &beamspot_dxdz_error_);
    theTree->Branch("beamspot_dydz_error",         &beamspot_dydz_error_);
    theTree->Branch("beamspot_sigmaZError",        &beamspot_sigmaZError_);
    theTree->Branch("beamspot_beamWidthX",         &beamspot_beamWidthX_);
    theTree->Branch("beamspot_beamWidthY",         &beamspot_beamWidthY_);
    theTree->Branch("beamspot_beamWidthX_error",   &beamspot_beamWidthX_error_);
    theTree->Branch("beamspot_beamWidthY_error",   &beamspot_beamWidthY_error_);
  }

  if (doHLT) {
    // ### Trigger ###
    //theTree->Branch("trigTable",     &TrigTable);
    TrigTable_store=nullptr;
    TrigResult_store=nullptr;
    TrigPrescales_store=nullptr;
    theTree->Branch("trigResult",    &trigResult);
    theTree->Branch("trigPrescales", &trigPrescales);
    theTree->Branch("l1Table",       &l1Table);
    theTree->Branch("l1Prescales",   &l1Prescales);
  
    SetupTriggerStorageVectors();
    SetupTriggerBranches();

  }


  if (doGenParticles_) {
    theTree->Branch("gen_nBs"			,&gen_nBs_);
    theTree->Branch("gen_Bs_pt"			,&gen_Bs_pt_);
    theTree->Branch("gen_Bs_energy"		,&gen_Bs_energy_);
    theTree->Branch("gen_Bs_eta"		,&gen_Bs_eta_);
    theTree->Branch("gen_Bs_phi"		,&gen_Bs_phi_);
    theTree->Branch("gen_Bs_pz"			,&gen_Bs_pz_);
    theTree->Branch("gen_Bs_pdgId"		,&gen_Bs_pdgId_);

    theTree->Branch("gen_nBsMuonM"		,&gen_nBsMuonM_);
    theTree->Branch("gen_BsMuonM_pt"		,&gen_BsMuonM_pt_);
    theTree->Branch("gen_BsMuonM_eta"		,&gen_BsMuonM_eta_);
    theTree->Branch("gen_BsMuonM_phi"		,&gen_BsMuonM_phi_);

    theTree->Branch("gen_nBsMuonP"		,&gen_nBsMuonP_);
    theTree->Branch("gen_BsMuonP_pt"		,&gen_BsMuonP_pt_);
    theTree->Branch("gen_BsMuonP_eta"		,&gen_BsMuonP_eta_);
    theTree->Branch("gen_BsMuonP_phi"		,&gen_BsMuonP_phi_);

    theTree->Branch("gen_nBsPhoton"		,&gen_nBsPhoton_);
    theTree->Branch("gen_BsPhoton_pt"		,&gen_BsPhoton_pt_);
    theTree->Branch("gen_BsPhoton_energy"	,&gen_BsPhoton_energy_);
    theTree->Branch("gen_BsPhoton_eta"		,&gen_BsPhoton_eta_);
    theTree->Branch("gen_BsPhoton_phi"		,&gen_BsPhoton_phi_);
    
    if(doFlatPt_){
    theTree->Branch("nMC",          &nMC_);
    theTree->Branch("mcPID",        &mcPID_);
    theTree->Branch("mcStatus",     &mcStatus_);
    theTree->Branch("mcVtx_x",      &mcVtx_x_);
    theTree->Branch("mcVtx_y",      &mcVtx_y_);
    theTree->Branch("mcVtx_z",      &mcVtx_z_);
    theTree->Branch("mcPt",         &mcPt_);
    theTree->Branch("mcEta",        &mcEta_);
    theTree->Branch("mcPhi",        &mcPhi_);
    theTree->Branch("mcE",          &mcE_);
    theTree->Branch("mcEt",         &mcEt_);
    theTree->Branch("mcMass",       &mcMass_);
   }
 
  }

  if(doMuons_){

    // ### mu- ###  
    theTree->Branch("nMuM",             &nMuM_);
    theTree->Branch("mumHighPurity",    &mumHighPurity_);
    theTree->Branch("mumPt",            &mumPt_);
    theTree->Branch("mumEta",           &mumEta_);
    theTree->Branch("mumPhi",           &mumPhi_);
    theTree->Branch("mumCL",            &mumCL_);
    theTree->Branch("mumNormChi2",      &mumNormChi2_);
    theTree->Branch("mumVx",            &mumVx_);
    theTree->Branch("mumVy",            &mumVy_);
    theTree->Branch("mumVz",            &mumVz_);
    theTree->Branch("mumDCABS",         &mumDCABS_);
    theTree->Branch("mumDCABSE",        &mumDCABSE_);
    theTree->Branch("mumFracHits",      &mumFracHits_);
    theTree->Branch("mumdxyBS",         &mumdxyBS_);
    theTree->Branch("mumdzBS",          &mumdzBS_);
    theTree->Branch("mumIdx",           &mumIdx_);
    theTree->Branch("mumNPixHits",      &mumNPixHits_);
    theTree->Branch("mumNPixLayers",    &mumNPixLayers_);
    theTree->Branch("mumNTrkHits",      &mumNTrkHits_);
    theTree->Branch("mumNTrkLayers",    &mumNTrkLayers_);
    theTree->Branch("mumNMuonHits",     &mumNMuonHits_);
    theTree->Branch("mumNMatchStation", &mumNMatchStation_);
    theTree->Branch("mum_isGlobalMuon", &mum_isGlobalMuon_);
    theTree->Branch("mum_isTrackerMuon",&mum_isTrackerMuon_);
    theTree->Branch("mum_StandAloneMuon",&mum_StandAloneMuon_);
    theTree->Branch("mum_isCaloMuon",    &mum_isCaloMuon_);
    theTree->Branch("mum_isPFMuon",      &mum_isPFMuon_);

    theTree->Branch("mum_selector"          	,  &mum_selector_); 
    theTree->Branch("mum_isIsolationValid"  	,  &mum_isIsolationValid_);
    theTree->Branch("mum_isPFIsolationValid"	,  &mum_isPFIsolationValid_);
   
    theTree->Branch("mum_isolationR03_trackSumPt"			,&mum_isolationR03_trackSumPt_);
    theTree->Branch("mum_isolationR03_trackEcalSumEt"			,&mum_isolationR03_trackEcalSumEt_);
    theTree->Branch("mum_isolationR03_trackHcalSumEt"			,&mum_isolationR03_trackHcalSumEt_);
    theTree->Branch("mum_isolationR03_trackHoSumEt"			,&mum_isolationR03_trackHoSumEt_);
    theTree->Branch("mum_isolationR03_trackNTracks"			,&mum_isolationR03_trackNTracks_);
    theTree->Branch("mum_isolationR03_trackNJets"			,&mum_isolationR03_trackNJets_);
    theTree->Branch("mum_isolationR03_trackerVetoSumPt"			,&mum_isolationR03_trackerVetoSumPt_);
    theTree->Branch("mum_isolationR03_emVetoSumEt"			,&mum_isolationR03_emVetoSumEt_);
    theTree->Branch("mum_isolationR03_hadVetoSumEt"			,&mum_isolationR03_hadVetoSumEt_);
    theTree->Branch("mum_isolationR03_hoVetoEt"				,&mum_isolationR03_hoVetoEt_);

    theTree->Branch("mum_isolationR05_trackSumPt"			,&mum_isolationR05_trackSumPt_);
    theTree->Branch("mum_isolationR05_trackEcalSumEt"			,&mum_isolationR05_trackEcalSumEt_);
    theTree->Branch("mum_isolationR05_trackHcalSumEt"			,&mum_isolationR05_trackHcalSumEt_);
    theTree->Branch("mum_isolationR05_trackHoSumEt"			,&mum_isolationR05_trackHoSumEt_);
    theTree->Branch("mum_isolationR05_trackNTracks"			,&mum_isolationR05_trackNTracks_);
    theTree->Branch("mum_isolationR05_trackNJets"			,&mum_isolationR05_trackNJets_);
    theTree->Branch("mum_isolationR05_trackerVetoSumPt"			,&mum_isolationR05_trackerVetoSumPt_);
    theTree->Branch("mum_isolationR05_emVetoSumEt"			,&mum_isolationR05_emVetoSumEt_);
    theTree->Branch("mum_isolationR05_hadVetoSumEt"			,&mum_isolationR05_hadVetoSumEt_);
    theTree->Branch("mum_isolationR05_hoVetoEt"				,&mum_isolationR05_hoVetoEt_);
                                                                                                                                                                                                                                                          
    theTree->Branch("mum_PFIsolationR03_sumChargedHadronPt"		,&mum_PFIsolationR03_sumChargedHadronPt_);
    theTree->Branch("mum_PFIsolationR03_sumChargedParticlePt"		,&mum_PFIsolationR03_sumChargedParticlePt_);
    theTree->Branch("mum_PFIsolationR03_sumNeutralHadronEt"		,&mum_PFIsolationR03_sumNeutralHadronEt_);
    theTree->Branch("mum_PFIsolationR03_sumPhotonEt"			,&mum_PFIsolationR03_sumPhotonEt_);
    theTree->Branch("mum_PFIsolationR03_sumNeutralHadronEtHighThreshold",&mum_PFIsolationR03_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR03_sumPhotonEtHighThreshold"	,&mum_PFIsolationR03_sumPhotonEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR03_sumPUPt"			        ,&mum_PFIsolationR03_sumPUPt_);
                                                                                                                             
    theTree->Branch("mum_PFIsolationR04_sumChargedHadronPt"		,&mum_PFIsolationR04_sumChargedHadronPt_);
    theTree->Branch("mum_PFIsolationR04_sumChargedParticlePt"		,&mum_PFIsolationR04_sumChargedParticlePt_);
    theTree->Branch("mum_PFIsolationR04_sumNeutralHadronEt"		,&mum_PFIsolationR04_sumNeutralHadronEt_);
    theTree->Branch("mum_PFIsolationR04_sumPhotonEt"			,&mum_PFIsolationR04_sumPhotonEt_);
    theTree->Branch("mum_PFIsolationR04_sumNeutralHadronEtHighThreshold",&mum_PFIsolationR04_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR04_sumPhotonEtHighThreshold"	,&mum_PFIsolationR04_sumPhotonEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR04_sumPUPt"			        ,&mum_PFIsolationR04_sumPUPt_);

    // ### mu+ ###  
    theTree->Branch("nMuP",             &nMuP_);
    theTree->Branch("mupHighPurity",    &mupHighPurity_);
    theTree->Branch("mupPt",            &mupPt_);
    theTree->Branch("mupEta",           &mupEta_);
    theTree->Branch("mupPhi",           &mupPhi_);
    theTree->Branch("mupCL",            &mupCL_);
    theTree->Branch("mupNormChi2",      &mupNormChi2_);
    theTree->Branch("mupVx",            &mupVx_);
    theTree->Branch("mupVy",            &mupVy_);
    theTree->Branch("mupVz",            &mupVz_);
    theTree->Branch("mupDCABS",         &mupDCABS_);
    theTree->Branch("mupDCABSE",        &mupDCABSE_);
    theTree->Branch("mupFracHits",      &mupFracHits_);
    theTree->Branch("mupdxyBS",         &mupdxyBS_);
    theTree->Branch("mupdzBS",          &mupdzBS_);
    theTree->Branch("mupIdx_",          &mupIdx_);
    theTree->Branch("mupNPixHits",      &mupNPixHits_);
    theTree->Branch("mupNPixLayers",    &mupNPixLayers_);
    theTree->Branch("mupNTrkHits",      &mupNTrkHits_);
    theTree->Branch("mupNTrkLayers",    &mupNTrkLayers_);
    theTree->Branch("mupNMuonHits",     &mupNMuonHits_);
    theTree->Branch("mupNMatchStation", &mupNMatchStation_);
    theTree->Branch("mup_isGlobalMuon", &mup_isGlobalMuon_);
    theTree->Branch("mup_isTrackerMuon",&mup_isTrackerMuon_);
    theTree->Branch("mup_StandAloneMuon",&mup_StandAloneMuon_);
    theTree->Branch("mup_isCaloMuon",    &mup_isCaloMuon_);
    theTree->Branch("mup_isPFMuon",      &mup_isPFMuon_);

    theTree->Branch("mup_selector"          	,  &mup_selector_); 
    theTree->Branch("mup_isIsolationValid"  	,  &mup_isIsolationValid_);
    theTree->Branch("mup_isPFIsolationValid"	,  &mup_isPFIsolationValid_);
   
    theTree->Branch("mup_isolationR03_trackSumPt"			,&mup_isolationR03_trackSumPt_);
    theTree->Branch("mup_isolationR03_trackEcalSumEt"			,&mup_isolationR03_trackEcalSumEt_);
    theTree->Branch("mup_isolationR03_trackHcalSumEt"			,&mup_isolationR03_trackHcalSumEt_);
    theTree->Branch("mup_isolationR03_trackHoSumEt"			,&mup_isolationR03_trackHoSumEt_);
    theTree->Branch("mup_isolationR03_trackNTracks"			,&mup_isolationR03_trackNTracks_);
    theTree->Branch("mup_isolationR03_trackNJets"			,&mup_isolationR03_trackNJets_);
    theTree->Branch("mup_isolationR03_trackerVetoSumPt"			,&mup_isolationR03_trackerVetoSumPt_);
    theTree->Branch("mup_isolationR03_emVetoSumEt"			,&mup_isolationR03_emVetoSumEt_);
    theTree->Branch("mup_isolationR03_hadVetoSumEt"			,&mup_isolationR03_hadVetoSumEt_);
    theTree->Branch("mup_isolationR03_hoVetoEt"				,&mup_isolationR03_hoVetoEt_);

    theTree->Branch("mup_isolationR05_trackSumPt"			,&mup_isolationR05_trackSumPt_);
    theTree->Branch("mup_isolationR05_trackEcalSumEt"			,&mup_isolationR05_trackEcalSumEt_);
    theTree->Branch("mup_isolationR05_trackHcalSumEt"			,&mup_isolationR05_trackHcalSumEt_);
    theTree->Branch("mup_isolationR05_trackHoSumEt"			,&mup_isolationR05_trackHoSumEt_);
    theTree->Branch("mup_isolationR05_trackNTracks"			,&mup_isolationR05_trackNTracks_);
    theTree->Branch("mup_isolationR05_trackNJets"			,&mup_isolationR05_trackNJets_);
    theTree->Branch("mup_isolationR05_trackerVetoSumPt"			,&mup_isolationR05_trackerVetoSumPt_);
    theTree->Branch("mup_isolationR05_emVetoSumEt"			,&mup_isolationR05_emVetoSumEt_);
    theTree->Branch("mup_isolationR05_hadVetoSumEt"			,&mup_isolationR05_hadVetoSumEt_);
    theTree->Branch("mup_isolationR05_hoVetoEt"				,&mup_isolationR05_hoVetoEt_);
                                                                                                                                                                                                                                                          
    theTree->Branch("mup_PFIsolationR03_sumChargedHadronPt"		,&mup_PFIsolationR03_sumChargedHadronPt_);
    theTree->Branch("mup_PFIsolationR03_sumChargedParticlePt"		,&mup_PFIsolationR03_sumChargedParticlePt_);
    theTree->Branch("mup_PFIsolationR03_sumNeutralHadronEt"		,&mup_PFIsolationR03_sumNeutralHadronEt_);
    theTree->Branch("mup_PFIsolationR03_sumPhotonEt"			,&mup_PFIsolationR03_sumPhotonEt_);
    theTree->Branch("mup_PFIsolationR03_sumNeutralHadronEtHighThreshold",&mup_PFIsolationR03_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR03_sumPhotonEtHighThreshold"	,&mup_PFIsolationR03_sumPhotonEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR03_sumPUPt"			        ,&mup_PFIsolationR03_sumPUPt_);
                                                                                                                             
    theTree->Branch("mup_PFIsolationR04_sumChargedHadronPt"		,&mup_PFIsolationR04_sumChargedHadronPt_);
    theTree->Branch("mup_PFIsolationR04_sumChargedParticlePt"		,&mup_PFIsolationR04_sumChargedParticlePt_);
    theTree->Branch("mup_PFIsolationR04_sumNeutralHadronEt"		,&mup_PFIsolationR04_sumNeutralHadronEt_);
    theTree->Branch("mup_PFIsolationR04_sumPhotonEt"			,&mup_PFIsolationR04_sumPhotonEt_);
    theTree->Branch("mup_PFIsolationR04_sumNeutralHadronEtHighThreshold",&mup_PFIsolationR04_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR04_sumPhotonEtHighThreshold"	,&mup_PFIsolationR04_sumPhotonEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR04_sumPUPt"			        ,&mup_PFIsolationR04_sumPUPt_);

    // ### mu+ mu- Mass ###
    theTree->Branch("mumuPt",    &mumuPt_);
    theTree->Branch("mumuEta",   &mumuEta_);
    theTree->Branch("mumuRapidity",&mumuRapidity_);
    theTree->Branch("mumuPhi",   &mumuPhi_);
    theTree->Branch("mumuMass",  &mumuMass_);
    theTree->Branch("mumuMassE", &mumuMassE_);
    theTree->Branch("mumuPx",    &mumuPx_);
    theTree->Branch("mumuPy",    &mumuPy_);
    theTree->Branch("mumuPz",    &mumuPz_);
    theTree->Branch("mumuCovPxPx_",    &mumuCovPxPx_);
    theTree->Branch("mumuCovPxPy_",    &mumuCovPxPy_);
    theTree->Branch("mumuCovPxPz_",    &mumuCovPxPz_);
    theTree->Branch("mumuCovPyPy_",    &mumuCovPyPy_);
    theTree->Branch("mumuCovPyPz_",    &mumuCovPyPz_);
    theTree->Branch("mumuCovPzPz_",    &mumuCovPzPz_);
    theTree->Branch("mumuDR",    &mumuDR_);
    theTree->Branch("mumuParentMuP",    &mumuParentMuP_);
    theTree->Branch("mumuParentMuM",    &mumuParentMuM_);

    // ### mu+ mu- Vtx ###
    theTree->Branch("mumuVtxCL",       &mumuVtxCL_);
    theTree->Branch("mumuVtxX",        &mumuVtxX_);
    theTree->Branch("mumuVtxY",        &mumuVtxY_);
    theTree->Branch("mumuVtxZ",        &mumuVtxZ_);
    theTree->Branch("mumuVtxCovXX_",   &mumuVtxCovXX_);
    theTree->Branch("mumuVtxCovXY_",   &mumuVtxCovXY_);
    theTree->Branch("mumuVtxCovXZ_",   &mumuVtxCovYY_);
    theTree->Branch("mumuVtxCovYY_",   &mumuVtxCovYY_);
    theTree->Branch("mumuVtxCovYZ_",   &mumuVtxCovYZ_);
    theTree->Branch("mumuVtxCovZZ_",   &mumuVtxCovZZ_);
    theTree->Branch("mumuVtxChi2",     &mumuVtxChi2_);
    theTree->Branch("mumuVtxNdof",     &mumuVtxNdof_);

    theTree->Branch("mumuCosAlphaBS",  &mumuCosAlphaBS_);
    theTree->Branch("mumuCosAlphaBSE", &mumuCosAlphaBSE_);
    theTree->Branch("mumuLBS",         &mumuLBS_);
    theTree->Branch("mumuLBSE",        &mumuLBSE_);
    theTree->Branch("mumuDCA",         &mumuDCA_);


  }
  
  if (doPhotons_) {
    theTree->Branch("nPho",                  &nPho_);
    theTree->Branch("phoE",                  &phoE_);
    theTree->Branch("phoEt",                 &phoEt_);
    theTree->Branch("phoEta",                &phoEta_);
    theTree->Branch("phoPhi",                &phoPhi_);
    theTree->Branch("phoSigmaE",               &phoSigmaE_);
    theTree->Branch("phoCalibE",               &phoCalibE_);
    theTree->Branch("phoCalibEt",              &phoCalibEt_);
    theTree->Branch("phoSCE",                  &phoSCE_);
    theTree->Branch("phoSCEt",                 &phoSCEt_);
    theTree->Branch("phoSCRawE",               &phoSCRawE_);
    theTree->Branch("phoESEnP1",               &phoESEnP1_);
    theTree->Branch("phoESEnP2",               &phoESEnP2_);
    theTree->Branch("phoSCEta",                &phoSCEta_);
    theTree->Branch("phoSCPhi",                &phoSCPhi_);
    theTree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
    theTree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
    theTree->Branch("phoSCBrem",               &phoSCBrem_);
    theTree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
    theTree->Branch("phoEleVeto",              &phoEleVeto_);
    theTree->Branch("phoR9",                   &phoR9_);
    theTree->Branch("phoHoverE",               &phoHoverE_);
    theTree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
    theTree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
    theTree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
    theTree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
    theTree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
    theTree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
    theTree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
    
    theTree->Branch("phoPFChIso",              &phoPFChIso_);
    theTree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
    theTree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
    theTree->Branch("phoEcalPFClusterIso",     &phoEcalPFClusterIso_);
    theTree->Branch("phoHcalPFClusterIso",     &phoHcalPFClusterIso_);
    theTree->Branch("phoIDMVA",                &phoIDMVA_);
   
    theTree->Branch("phoSeedTime",             &phoSeedTime_);
    theTree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
    theTree->Branch("phoMIPTotEnergy",         &phoMIPTotEnergy_);
    theTree->Branch("phoMIPChi2",                      &phoMIPChi2_);
    theTree->Branch("phoMIPSlope",                     &phoMIPSlope_);
    theTree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
    theTree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
    theTree->Branch("phoMIPIsHalo",                    &phoMIPIsHalo_);

  }

  if (doPFPhotons_) {
    theTree->Branch("nPFPho",                  &nPFPho_);
    theTree->Branch("phoPFE",                  &phoPFE_);
    theTree->Branch("phoPFEt",                 &phoPFEt_);
    theTree->Branch("phoPFEta",                &phoPFEta_);
    theTree->Branch("phoPFPhi",                &phoPFPhi_);
  }

  if (doSuperClusters_) {
    theTree->Branch("nSC",                  &nSC_);
    theTree->Branch("scE",                  &scE_);
    theTree->Branch("scEt",                 &scEt_);
    theTree->Branch("scRawE",               &scRawE_);
    theTree->Branch("scEta",                &scEta_);
    theTree->Branch("scPhi",                &scPhi_);
    theTree->Branch("scX",        &scX_);
    theTree->Branch("scY",        &scY_);
    theTree->Branch("scZ",        &scZ_);
    theTree->Branch("scEtaWidth", &scEtaWidth_);
    theTree->Branch("scPhiWidth", &scPhiWidth_);
    theTree->Branch("scRawEt",    &scRawEt_);
    theTree->Branch("scMinDrWithGsfElectornSC_",  &scMinDrWithGsfElectornSC_);
    theTree->Branch("scFoundGsfMatch_" ,        &scFoundGsfMatch_);
    theTree->Branch("scE5x5",   &scE5x5_);
    theTree->Branch("scE2x2Ratio",   &scE2x2Ratio_);
    theTree->Branch("scE3x3Ratio",   &scE3x3Ratio_);
    theTree->Branch("scEMaxRatio",   &scEMaxRatio_);
    theTree->Branch("scE2ndRatio",   &scE2ndRatio_);
    theTree->Branch("scETopRatio",   &scETopRatio_);
    theTree->Branch("scERightRatio",   &scERightRatio_);
    theTree->Branch("scEBottomRatio",   &scEBottomRatio_);
    theTree->Branch("scELeftRatio",   &scELeftRatio_);
    theTree->Branch("scE2x5MaxRatio",   &scE2x5MaxRatio_);
    theTree->Branch("scE2x5TopRatio",   &scE2x5TopRatio_);
    theTree->Branch("scE2x5RightRatio",   &scE2x5RightRatio_);
    theTree->Branch("scE2x5BottomRatio",   &scE2x5BottomRatio_); 
    theTree->Branch("scE2x5LeftRatio",   &scE2x5LeftRatio_); 
    theTree->Branch("scSwissCross",   &scSwissCross_); 
    theTree->Branch("scR9",   &scR9_);
    theTree->Branch("scSigmaIetaIeta",   &scSigmaIetaIeta_);
    theTree->Branch("scSigmaIetaIphi",   &scSigmaIetaIphi_);
    theTree->Branch("scSigmaIphiIphi",   &scSigmaIphiIphi_);
    theTree->Branch("scFull5x5_e5x5",   &scFull5x5_e5x5_);
    theTree->Branch("scFull5x5_e2x2Ratio",   &scFull5x5_e2x2Ratio_);
    theTree->Branch("scFull5x5_e3x3Ratio",   &scFull5x5_e3x3Ratio_);
    theTree->Branch("scFull5x5_eMaxRatio",   &scFull5x5_eMaxRatio_);
    theTree->Branch("scFull5x5_e2ndRatio",   &scFull5x5_e2ndRatio_);
    theTree->Branch("scFull5x5_eTopRatio",   &scFull5x5_eTopRatio_);
    theTree->Branch("scFull5x5_eRightRatio",   &scFull5x5_eRightRatio_);
    theTree->Branch("scFull5x5_eBottomRatio",   &scFull5x5_eBottomRatio_);
    theTree->Branch("scFull5x5_eLeftRatio",   &scFull5x5_eLeftRatio_);
    theTree->Branch("scFull5x5_e2x5MaxRatio",   &scFull5x5_e2x5MaxRatio_);
    theTree->Branch("scFull5x5_e2x5TopRatio",   &scFull5x5_e2x5TopRatio_);
    theTree->Branch("scFull5x5_e2x5RightRatio",   &scFull5x5_e2x5RightRatio_);
    theTree->Branch("scFull5x5_e2x5BottomRatio",   &scFull5x5_e2x5BottomRatio_); 
    theTree->Branch("scFull5x5_e2x5LeftRatio",   &scFull5x5_e2x5LeftRatio_); 
    theTree->Branch("scFull5x5_swissCross",   &scFull5x5_swissCross_); 
    theTree->Branch("scFull5x5_r9",   &scFull5x5_r9_);
    theTree->Branch("scFull5x5_sigmaIetaIeta",   &scFull5x5_sigmaIetaIeta_);
    theTree->Branch("scFull5x5_sigmaIetaIphi",   &scFull5x5_sigmaIetaIphi_);
    theTree->Branch("scFull5x5_sigmaIphiIphi",   &scFull5x5_sigmaIphiIphi_);
  
    //theTree->Branch("nhcalRechit",      &nhcalRechit_);
    //theTree->Branch("hcalRechitIEta",   &hcalRechitIEta_);
    //theTree->Branch("hcalRechitIPhi",   &hcalRechitIPhi_);
    //theTree->Branch("hcalRechitEnergy", &hcalRechitEnergy_);

    theTree->Branch("scNHcalRecHitInDIEta5IPhi5",              &scNHcalRecHitInDIEta5IPhi5);
    theTree->Branch("scEFromHcalRecHitInDIEta5IPhi5",              &scEFromHcalRecHitInDIEta5IPhi5);
    theTree->Branch("scNHcalRecHitInDIEta2IPhi2",              &scNHcalRecHitInDIEta2IPhi2);
    theTree->Branch("scEFromHcalRecHitInDIEta2IPhi2",              &scEFromHcalRecHitInDIEta2IPhi2);

    theTree->Branch("scPFChIso1",              &scPFChIso1_);
    theTree->Branch("scPFChIso2",              &scPFChIso2_);
    theTree->Branch("scPFChIso3",              &scPFChIso3_);
    theTree->Branch("scPFChIso4",              &scPFChIso4_);
    theTree->Branch("scPFChIso5",              &scPFChIso5_);
    
    theTree->Branch("scPFPhoIso1",             &scPFPhoIso1_);
    theTree->Branch("scPFPhoIso2",             &scPFPhoIso2_);
    theTree->Branch("scPFPhoIso3",             &scPFPhoIso3_);
    theTree->Branch("scPFPhoIso4",             &scPFPhoIso4_);
    theTree->Branch("scPFPhoIso5",             &scPFPhoIso5_);
    
    theTree->Branch("scPFNeuIso1",             &scPFNeuIso1_);
    theTree->Branch("scPFNeuIso2",             &scPFNeuIso2_);
    theTree->Branch("scPFNeuIso3",             &scPFNeuIso3_);
    theTree->Branch("scPFNeuIso4",             &scPFNeuIso4_);
    theTree->Branch("scPFNeuIso5",             &scPFNeuIso5_);
  }
  if(doParticleFlow)
  {
     addParticleFlowBranches();
  }
  if(doHCALClusters)
  {
     addHCALCluserBranches();
  }
  if(doECALClusters)
  {
     addECALCluserBranches();
  }
  if(doGeneralTracks)
  {
     addGeneralTracksBranches();
  }
  if(doPrimaryVetrices)
  {
     addPrimaryVertexBranches();
  }

}

//BsToMuMuGammaNTuplizer::~BsToMuMuGammaNTuplizer(){
//}


// ------------ method called for each event  ------------
void 
BsToMuMuGammaNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  // ## BEAMSOPT STUFF  ## //
  beamspot_x_  = 0.0   ;
  beamspot_y_  = 0.0   ;
  beamspot_z_  = 0.0   ;
  beamspot_x_error_  = 0.0   ;
  beamspot_y_error_  = 0.0   ;
  beamspot_z_error_  = 0.0   ;
  beamspot_dxdz_  = 0.0   ;
  beamspot_dydz_  = 0.0   ;
  beamspot_sigmaZ_  = 0.0   ;
  beamspot_dxdz_error_  = 0.0   ;
  beamspot_dydz_error_  = 0.0   ;
  beamspot_sigmaZError_  = 0.0   ;
  beamspot_beamWidthX_  = 0.0   ;
  beamspot_beamWidthY_  = 0.0   ;
  beamspot_beamWidthX_error_  = 0.0   ;
  beamspot_beamWidthY_error_  = 0.0   ;


  if (doHLT) {

    // ### Trigger ###
    //TrigTable.clear();
    trigNames.clear();
    trigResult.clear();
    trigPrescales.clear();
    l1Table.clear();
    l1Prescales.clear();
    ClearTrggerStorages();
  
  }
  if (doGenParticles_) {
  
    gen_nBs_ = 0;
    gen_nBsMuonM_ = 0 ;
    gen_nBsMuonP_ = 0 ;
    gen_nBsPhoton_ = 0 ;

    gen_Bs_pt_.clear() ;
    gen_Bs_energy_.clear() ;
    gen_Bs_eta_.clear() ;
    gen_Bs_phi_.clear() ;
    gen_Bs_pz_.clear() ;
    gen_Bs_pdgId_.clear();
    gen_BsMuonM_pt_.clear() ;
    gen_BsMuonM_eta_.clear() ;
    gen_BsMuonM_phi_.clear();
    gen_BsMuonP_pt_.clear() ;
    gen_BsMuonP_eta_.clear() ;
    gen_BsMuonP_phi_.clear();
    gen_BsPhoton_pt_.clear() ;
    gen_BsPhoton_energy_.clear() ;
    gen_BsPhoton_eta_.clear() ;
    gen_BsPhoton_phi_.clear();
    if(doFlatPt_){
    nMC_ = 0;
    mcPID_                .clear();
    mcStatus_             .clear();
    mcVtx_x_              .clear();
    mcVtx_y_              .clear();
    mcVtx_z_              .clear();
    mcPt_                 .clear();
    mcEta_                .clear();
    mcPhi_                .clear();
    mcE_                  .clear();
    mcEt_                 .clear();
    mcMass_               .clear();
     }
}

  if(doMuons_){
    nMuM_ = 0;
    nMuP_ = 0; 
    mumHighPurity_.clear();
    mumPt_.clear();
    mumEta_.clear();
    mumPhi_.clear();
    mumCL_.clear(); 
    mumNormChi2_.clear();
    mumVx_.clear();
    mumVy_.clear();
    mumVz_.clear();
 
    mumDCABS_.clear();
    mumDCABSE_.clear();
    mumFracHits_.clear();
    mumdxyBS_.clear();
    mumdzBS_.clear();
    mumIdx_.clear();
    mumCharge_.clear();
    mumNPixHits_.clear();
    mumNPixLayers_.clear();
    mumNTrkHits_.clear();
    mumNTrkLayers_.clear();
    mumNMuonHits_.clear();
    mumNMatchStation_.clear();
    mum_isGlobalMuon_.clear();
    mum_isTrackerMuon_.clear();
    mum_StandAloneMuon_.clear();
    mum_isCaloMuon_.clear();
    mum_isPFMuon_.clear();

    mum_selector_.clear();
 
    mum_isIsolationValid_.clear();
    mum_isPFIsolationValid_.clear();
   
    mum_isolationR03_trackSumPt_.clear();
    mum_isolationR03_trackEcalSumEt_.clear();
    mum_isolationR03_trackHcalSumEt_.clear();
    mum_isolationR03_trackHoSumEt_.clear();
    mum_isolationR03_trackNTracks_.clear();
    mum_isolationR03_trackNJets_.clear();
    mum_isolationR03_trackerVetoSumPt_.clear();
    mum_isolationR03_emVetoSumEt_.clear();
    mum_isolationR03_hadVetoSumEt_.clear();
    mum_isolationR03_hoVetoEt_.clear();
  
    mum_isolationR05_trackSumPt_.clear();
    mum_isolationR05_trackEcalSumEt_.clear();
    mum_isolationR05_trackHcalSumEt_.clear();
    mum_isolationR05_trackHoSumEt_.clear();
    mum_isolationR05_trackNTracks_.clear();
    mum_isolationR05_trackNJets_.clear();
    mum_isolationR05_trackerVetoSumPt_.clear();
    mum_isolationR05_emVetoSumEt_.clear();
    mum_isolationR05_hadVetoSumEt_.clear();
    mum_isolationR05_hoVetoEt_.clear();
  
    mum_PFIsolationR03_sumChargedHadronPt_.clear();
    mum_PFIsolationR03_sumChargedParticlePt_.clear();
    mum_PFIsolationR03_sumNeutralHadronEt_.clear();
    mum_PFIsolationR03_sumPhotonEt_.clear();
    mum_PFIsolationR03_sumNeutralHadronEtHighThreshold_.clear();
    mum_PFIsolationR03_sumPhotonEtHighThreshold_.clear();
    mum_PFIsolationR03_sumPUPt_.clear();
    
    mum_PFIsolationR04_sumChargedHadronPt_.clear();
    mum_PFIsolationR04_sumChargedParticlePt_.clear();
    mum_PFIsolationR04_sumNeutralHadronEt_.clear();
    mum_PFIsolationR04_sumPhotonEt_.clear();
    mum_PFIsolationR04_sumNeutralHadronEtHighThreshold_.clear();
    mum_PFIsolationR04_sumPhotonEtHighThreshold_.clear();
    mum_PFIsolationR04_sumPUPt_.clear();
  
  
  
    mupHighPurity_.clear();
    mupPt_.clear();
    mupEta_.clear();
    mupPhi_.clear();
    mupCL_.clear(); 
    mupNormChi2_.clear();
    mupVx_.clear();
    mupVy_.clear();
    mupVz_.clear();
    mupDCABS_.clear();
    mupDCABSE_.clear();
  
    mupFracHits_.clear();
    mupdxyBS_.clear();
    mupdzBS_.clear();
    mupIdx_.clear();
    mupCharge_.clear();
    mupNPixHits_.clear();
    mupNPixLayers_.clear();
    mupNTrkHits_.clear();
    mupNTrkLayers_.clear();
    mupNMuonHits_.clear();
    mupNMatchStation_.clear();
    mup_isGlobalMuon_.clear();
    mup_isTrackerMuon_.clear();
    mup_StandAloneMuon_.clear();
    mup_isCaloMuon_.clear();
    mup_isPFMuon_.clear();

    mup_selector_.clear();
 
    mup_isIsolationValid_.clear();
    mup_isPFIsolationValid_.clear();
   
    mup_isolationR03_trackSumPt_.clear();
    mup_isolationR03_trackEcalSumEt_.clear();
    mup_isolationR03_trackHcalSumEt_.clear();
    mup_isolationR03_trackHoSumEt_.clear();
    mup_isolationR03_trackNTracks_.clear();
    mup_isolationR03_trackNJets_.clear();
    mup_isolationR03_trackerVetoSumPt_.clear();
    mup_isolationR03_emVetoSumEt_.clear();
    mup_isolationR03_hadVetoSumEt_.clear();
    mup_isolationR03_hoVetoEt_.clear();
  
    mup_isolationR05_trackSumPt_.clear();
    mup_isolationR05_trackEcalSumEt_.clear();
    mup_isolationR05_trackHcalSumEt_.clear();
    mup_isolationR05_trackHoSumEt_.clear();
    mup_isolationR05_trackNTracks_.clear();
    mup_isolationR05_trackNJets_.clear();
    mup_isolationR05_trackerVetoSumPt_.clear();
    mup_isolationR05_emVetoSumEt_.clear();
    mup_isolationR05_hadVetoSumEt_.clear();
    mup_isolationR05_hoVetoEt_.clear();
  
    mup_PFIsolationR03_sumChargedHadronPt_.clear();
    mup_PFIsolationR03_sumChargedParticlePt_.clear();
    mup_PFIsolationR03_sumNeutralHadronEt_.clear();
    mup_PFIsolationR03_sumPhotonEt_.clear();
    mup_PFIsolationR03_sumNeutralHadronEtHighThreshold_.clear();
    mup_PFIsolationR03_sumPhotonEtHighThreshold_.clear();
    mup_PFIsolationR03_sumPUPt_.clear();
    
    mup_PFIsolationR04_sumChargedHadronPt_.clear();
    mup_PFIsolationR04_sumChargedParticlePt_.clear();
    mup_PFIsolationR04_sumNeutralHadronEt_.clear();
    mup_PFIsolationR04_sumPhotonEt_.clear();
    mup_PFIsolationR04_sumNeutralHadronEtHighThreshold_.clear();
    mup_PFIsolationR04_sumPhotonEtHighThreshold_.clear();
    mup_PFIsolationR04_sumPUPt_.clear();

    // ### mu+ mu- variables ###
    mumuPt_.clear();
    mumuEta_.clear();
    mumuRapidity_.clear();
    mumuPhi_.clear();
    mumuMass_.clear();
    mumuMassE_.clear();
    mumuPx_.clear();
    mumuPy_.clear();
    mumuPz_.clear();
    mumuCovPxPx_.clear();
    mumuCovPxPy_.clear();
    mumuCovPxPz_.clear();
    mumuCovPyPy_.clear();
    mumuCovPyPz_.clear();
    mumuCovPzPz_.clear();
    mumuDR_.clear();
    mumuParentMuP_.clear();
    mumuParentMuM_.clear();
  
    mumuVtxCL_.clear();
    mumuVtxX_.clear();
    mumuVtxY_.clear();
    mumuVtxZ_.clear();  
    mumuVtxCovXX_.clear();  
    mumuVtxCovXY_.clear();  
    mumuVtxCovXZ_.clear();  
    mumuVtxCovYY_.clear();  
    mumuVtxCovYZ_.clear();  
    mumuVtxCovZZ_.clear();  
    mumuVtxChi2_.clear();
    mumuVtxNdof_.clear();


    mumuCosAlphaBS_.clear();
    mumuCosAlphaBSE_.clear(); 
    mumuLBS_.clear();
    mumuLBSE_.clear();
    mumuDCA_.clear();
 
    
    nMuP_=0;
    nMuM_=0;
  }
  
  if (doPhotons_) {
    nPho_ = 0;
    phoE_                 .clear();
    phoEt_                .clear();
    phoEta_               .clear();
    phoPhi_               .clear();
    phoSigmaE_              .clear();
    phoCalibE_              .clear();
    phoCalibEt_             .clear();
    phoSCE_                 .clear();
    phoSCEt_                 .clear();
    phoSCRawE_              .clear();
    phoESEnP1_              .clear();
    phoESEnP2_              .clear();
    phoSCEta_               .clear();
    phoSCPhi_               .clear();
    phoSCEtaWidth_          .clear();
    phoSCPhiWidth_          .clear();
    phoSCBrem_              .clear();
    phohasPixelSeed_        .clear();
    phoEleVeto_             .clear();
    phoR9_                  .clear();
    phoHoverE_              .clear();
    phoESEffSigmaRR_        .clear();
    phoSigmaIEtaIEtaFull5x5_.clear();
    phoSigmaIEtaIPhiFull5x5_.clear();
    phoSigmaIPhiIPhiFull5x5_.clear();
    phoE2x2Full5x5_         .clear();
    phoE5x5Full5x5_         .clear();
    phoR9Full5x5_           .clear();
    phoPFChIso_             .clear();
    phoPFPhoIso_            .clear();
    phoPFNeuIso_            .clear();
    phoEcalPFClusterIso_    .clear();
    phoHcalPFClusterIso_    .clear();
    phoIDMVA_               .clear();
    phoSeedTime_            .clear();
    phoSeedEnergy_          .clear();
    phoMIPTotEnergy_        .clear();  
    phoMIPChi2_           .clear();
    phoMIPSlope_          .clear();
    phoMIPIntercept_      .clear();
    phoMIPNhitCone_       .clear();
    phoMIPIsHalo_         .clear();
  }
  
  if (doPFPhotons_) {
    nPFPho_ = 0;
    phoPFE_                 .clear();
    phoPFEt_                .clear();
    phoPFEta_               .clear();
    phoPFPhi_               .clear();
  }

  if (doSuperClusters_) {
    nSC_ = 0;
    scE_                  .clear();
    scEt_                 .clear();
    scRawE_               .clear();
    scEta_                .clear();
    scPhi_                .clear();
    scX_          .clear();      
    scY_          .clear();      
    scZ_          .clear();      
    scEtaWidth_   .clear();         
    scPhiWidth_   .clear();         
    scRawEt_      .clear();   
    scMinDrWithGsfElectornSC_.clear();
    scFoundGsfMatch_.clear();
    
    scE5x5_.clear();
    scE2x2Ratio_.clear();
    scE3x3Ratio_.clear();
    scEMaxRatio_.clear();
    scE2ndRatio_.clear();
    scETopRatio_.clear();
    scERightRatio_.clear();
    scEBottomRatio_.clear();
    scELeftRatio_.clear();
    scE2x5MaxRatio_.clear();
    scE2x5TopRatio_.clear();
    scE2x5RightRatio_.clear();
    scE2x5BottomRatio_.clear();
    scE2x5LeftRatio_.clear();
    scSwissCross_.clear();
    scR9_.clear();
    scSigmaIetaIeta_.clear(); 
    scSigmaIetaIphi_.clear(); 
    scSigmaIphiIphi_.clear(); 
    scFull5x5_e5x5_.clear();
    scFull5x5_e2x2Ratio_.clear();
    scFull5x5_e3x3Ratio_.clear();
    scFull5x5_eMaxRatio_.clear();
    scFull5x5_e2ndRatio_.clear();
    scFull5x5_eTopRatio_.clear();
    scFull5x5_eRightRatio_.clear();
    scFull5x5_eBottomRatio_.clear();
    scFull5x5_eLeftRatio_.clear();
    scFull5x5_e2x5MaxRatio_.clear();
    scFull5x5_e2x5TopRatio_.clear();
    scFull5x5_e2x5RightRatio_.clear();
    scFull5x5_e2x5BottomRatio_.clear();
    scFull5x5_e2x5LeftRatio_.clear();
    scFull5x5_swissCross_.clear();
    scFull5x5_r9_.clear();
    scFull5x5_sigmaIetaIeta_.clear(); 
    scFull5x5_sigmaIetaIphi_.clear(); 
    scFull5x5_sigmaIphiIphi_.clear();  
  
    nhcalRechit_ =0;
    hcalRechitIEta_.clear();
    hcalRechitIPhi_.clear();
    hcalRechitEnergy_.clear();

    scNHcalRecHitInDIEta2IPhi2 .clear();
    scEFromHcalRecHitInDIEta2IPhi2             .clear();

    scNHcalRecHitInDIEta5IPhi5 .clear();
    scEFromHcalRecHitInDIEta5IPhi5             .clear();

    scPFChIso1_             .clear();
    scPFChIso2_             .clear();
    scPFChIso3_             .clear();
    scPFChIso4_             .clear();
    scPFChIso5_             .clear();
    
    scPFPhoIso1_            .clear();
    scPFPhoIso2_            .clear();
    scPFPhoIso3_            .clear();
    scPFPhoIso4_            .clear();
    scPFPhoIso5_            .clear();
    
    scPFNeuIso1_            .clear();
    scPFNeuIso2_            .clear();
    scPFNeuIso3_            .clear();
    scPFNeuIso4_            .clear();
    scPFNeuIso5_            .clear();


  }

  run_    = iEvent.id().run();
  event_  = iEvent.id().event();
  lumis_  = iEvent.luminosityBlock();
  isData_ = iEvent.isRealData();

  // Get magnetic field
    
  //  Get BeamSpot
  if(doBeamSpot)
  {
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken_, beamSpotH);
  reco::BeamSpot beamSpot = *beamSpotH;

 
 
  // adding  BEAMSOPT 
  beamspot_x_			= beamSpot.x0();  ;
  beamspot_y_			= beamSpot.y0();  ;
  beamspot_z_			= beamSpot.z0();  ;
  beamspot_x_error_		= beamSpot.x0Error();  ;
  beamspot_y_error_		= beamSpot.y0Error();  ;
  beamspot_z_error_		= beamSpot.z0Error();  ;
  beamspot_covXX        = beamSpot.covariance()(0,0);
  beamspot_covXY        = beamSpot.covariance()(0,1);
  beamspot_covXZ        = beamSpot.covariance()(0,2);
  beamspot_covYY        = beamSpot.covariance()(1,1);
  beamspot_covYZ        = beamSpot.covariance()(1,2);
  beamspot_covZZ        = beamSpot.covariance()(2,2);

  beamspot_dxdz_   		= beamSpot.dxdz();  ;
  beamspot_dydz_         	= beamSpot.dydz();  ;
  beamspot_sigmaZ_		= beamSpot.sigmaZ();  ;
  beamspot_dxdz_error_		= beamSpot.dxdzError();  ;
  beamspot_dydz_error_		= beamSpot.dydzError();  ;
  beamspot_sigmaZError_		= beamSpot.sigmaZ0Error();  ;
  beamspot_beamWidthX_		= beamSpot.BeamWidthX();  ;
  beamspot_beamWidthY_		= beamSpot.BeamWidthY();  ;
  beamspot_beamWidthX_error_	= beamSpot.BeamWidthXError();  ;
  beamspot_beamWidthY_error_	= beamSpot.BeamWidthXError();  ;
  }

  if(doPrimaryVetrices)
  {
        fillPrimaryVertexBranches(iEvent,iSetup);
  }
  // MC truth
  if (isMC and doGenParticles_) {
    fillGenParticles(iEvent);
  }

  if (doHLT) 		    fillHLT(iEvent);
  if (doMuons_)     	fillMuons(iEvent, iSetup);
  if (doPhotons_)    	fillPhotons(iEvent, iSetup);
  if (doPFPhotons_) 	fillPFPhotons(iEvent, iSetup);
  if (doSuperClusters_) {
        reco::Vertex pv(math::XYZPoint(0, 0, -999), math::Error<3>::type()); 
        fillSC(iEvent, iSetup,pv);
  }
  if (doHCALClusters) {
        fillHCALClusterCollection(iEvent,iSetup);
  }
  if (doECALClusters)  fillECALClusterCollection(iEvent,iSetup);
  if (doGeneralTracks)  fillGeneralTrackCollectionBranches(iEvent,iSetup);
  if (doParticleFlow)   fillPFCandiateCollection(iEvent,iSetup);
  theTree->Fill();
  
}



void BsToMuMuGammaNTuplizer::fillGenParticles(const edm::Event& iEvent)
{
  
  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  iEvent.getByToken(genParticlesCollection_, genParticleCollection);
  
  int phoMul(-1),muMMul(-1),muPMul(-1);
  
  if(doFlatPt_){
    for (auto p = genParticleCollection->begin(); p != genParticleCollection->end(); ++p) {
      mcPID_   .push_back(p->pdgId());
      mcStatus_.push_back(p->status());
      mcVtx_x_ .push_back(p->vx());
      mcVtx_y_ .push_back(p->vy());
      mcVtx_z_ .push_back(p->vz());
      mcPt_    .push_back(p->pt());
      mcEta_   .push_back(p->eta());
      mcPhi_   .push_back(p->phi());
      mcE_     .push_back(p->energy());
      mcEt_    .push_back(p->et());
      mcMass_  .push_back(p->mass());
      nMC_++;      
    } 
  }
  
  for(auto& aBsMeson : *genParticleCollection){
    
    if(abs(aBsMeson.pdgId())!=531) continue;
    
    for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
      
      auto& bsDaughter = *(aBsMeson.daughter(j));
      if(bsDaughter.pdgId() == -13) muPMul++;
      if(bsDaughter.pdgId() ==  13) muMMul++;
      if(bsDaughter.pdgId() ==  22) phoMul++;
    }
    if(muMMul <0 or muPMul <0 ) continue;
    if(phoMul <0 and doBsToMuMuGamma  ) continue;

    for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
      auto& bsDaughter = *(aBsMeson.daughter(j));
      if(bsDaughter.pdgId() == 13) {
	gen_BsMuonM_pt_.push_back(bsDaughter.pt());
	gen_BsMuonM_eta_.push_back(bsDaughter.eta());
	gen_BsMuonM_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonM_++;
      }
      if(bsDaughter.pdgId() == -13){
	gen_BsMuonP_pt_.push_back(bsDaughter.pt());
	gen_BsMuonP_eta_.push_back(bsDaughter.eta());
	gen_BsMuonP_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonP_++;
      }
      if(bsDaughter.pdgId() ==  22){
	gen_BsPhoton_pt_.push_back(bsDaughter.pt());
	gen_BsPhoton_energy_.push_back(bsDaughter.energy());
	gen_BsPhoton_eta_.push_back(bsDaughter.eta());
	gen_BsPhoton_phi_.push_back(bsDaughter.phi());
	gen_nBsPhoton_++;
	
      }
    }  //number of daughters

    gen_Bs_pt_.push_back(aBsMeson.pt());
    gen_Bs_energy_.push_back(aBsMeson.energy());
    gen_Bs_eta_.push_back(aBsMeson.eta());
    gen_Bs_phi_.push_back(aBsMeson.phi());
    gen_Bs_pz_.push_back(aBsMeson.pz());
    gen_Bs_pdgId_.push_back(aBsMeson.pdgId());
    gen_nBs_++;
    
  } // genparticle collection
  
} // fill gen particles


void BsToMuMuGammaNTuplizer::fillMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  
  // Get magnetic field
    
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);      
  //  Get BeamSpot
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken_, beamSpotH);
  reco::BeamSpot beamSpot = *beamSpotH;


  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);
  
  reco::TrackRef muTrackm,muTrackp, muTrack;
  reco::TransientTrack refitMupTT, refitMumTT;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  
  // variables
  float muonMass,muonMassErr;
  muonMass= Utility->muonMass;
  muonMassErr= Utility->muonMassErr;
  double chi2,ndof;

  //std::cout << " muon size:" << muons->size() << std::endl; 
  //start loop on first muon 
  for (uint32_t i=0;i<muons->size();i++){
   
    auto &mu = muons->at(i);
    if(muons->size() < 2 || mu.pt()  < pTMinMuons) continue;
    muTrack = mu.innerTrack();

    if ((muTrack.isNull() == true)) continue;
    
    const reco::TransientTrack muTrackTT( muTrack, &(*bFieldHandle));
    if (!muTrackTT.isValid()) continue;
    
    // # Compute mu- DCA to BeamSpot #
    theDCAXBS = muTrackTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
    double DCAmuBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
    double DCAmuBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();

    //// #############
    //// # Save: mu- #
    //// #############
    if(muTrack->charge() ==-1){ 
      mumHighPurity_.push_back( (int)muTrack->quality(reco::Track::highPurity));
      mumPt_.push_back(muTrack->pt());
      mumEta_.push_back(muTrack->eta());
      mumPhi_.push_back(muTrack->phi());
      mumCL_.push_back(TMath::Prob(muTrackTT.chi2(), static_cast<int>(rint(muTrackTT.ndof()))));
      mumNormChi2_.push_back(muTrack->normalizedChi2());
      mumVx_.push_back(mu.vx());
      mumVy_.push_back(mu.vy());
      mumVz_.push_back(mu.vz());
      
 
      mumDCABS_.push_back(DCAmuBS);
      mumDCABSE_.push_back(DCAmuBSErr);
      

      mumFracHits_.push_back(static_cast<double>(muTrack->hitPattern().numberOfValidHits()) / static_cast<double>(muTrack->hitPattern().numberOfValidHits() +
														  muTrack->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
														  muTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
														  muTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));

      mumdxyBS_.push_back(muTrackTT.track().dxy( (beamSpot.position() )));
      mumdzBS_.push_back(muTrackTT.track().dz(  (beamSpot.position() )));
      

      
      mumIdx_.push_back(i);
      mumNPixHits_.push_back(muTrack->hitPattern().numberOfValidPixelHits());
      mumNPixLayers_.push_back(muTrack->hitPattern().pixelLayersWithMeasurement());  
      mumNTrkHits_.push_back(muTrack->hitPattern().numberOfValidTrackerHits());
      mumNTrkLayers_.push_back(muTrack->hitPattern().trackerLayersWithMeasurement());
      if (mu.isGlobalMuon() == true) mumNMuonHits_.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
      else mumNMuonHits_.push_back(0);
      mumNMatchStation_.push_back(mu.numberOfMatchedStations());

      mum_isGlobalMuon_	      .push_back(mu.isGlobalMuon());
      mum_isTrackerMuon_	      .push_back(mu.isTrackerMuon());
      mum_StandAloneMuon_          .push_back(mu.isStandAloneMuon());
      mum_isCaloMuon_	      .push_back(mu.isCaloMuon());
      mum_isPFMuon_	              .push_back(mu.isPFMuon());

      mum_selector_.push_back(mu.selectors()); 
      mum_isIsolationValid_.push_back(mu.isIsolationValid());
      mum_isPFIsolationValid_.push_back(mu.isIsolationValid());

      auto &MuIsol03 = mu.isolationR03();
      mum_isolationR03_trackSumPt_.push_back(MuIsol03.sumPt);
      mum_isolationR03_trackEcalSumEt_.push_back(MuIsol03.emEt);
      mum_isolationR03_trackHcalSumEt_.push_back(MuIsol03.hadEt);
      mum_isolationR03_trackHoSumEt_.push_back(MuIsol03.hoEt);
      mum_isolationR03_trackNTracks_.push_back(MuIsol03.nTracks);
      mum_isolationR03_trackNJets_.push_back(MuIsol03.nJets);
      mum_isolationR03_trackerVetoSumPt_.push_back(MuIsol03.trackerVetoPt);
      mum_isolationR03_emVetoSumEt_.push_back(MuIsol03.emVetoEt);
      mum_isolationR03_hadVetoSumEt_.push_back(MuIsol03.hadVetoEt);
      mum_isolationR03_hoVetoEt_.push_back(MuIsol03.hoVetoEt);
  
      auto &MuIsol05 = mu.isolationR05();
      mum_isolationR05_trackSumPt_.push_back(MuIsol05.sumPt);
      mum_isolationR05_trackEcalSumEt_.push_back(MuIsol05.emEt);
      mum_isolationR05_trackHcalSumEt_.push_back(MuIsol05.hadEt);
      mum_isolationR05_trackHoSumEt_.push_back(MuIsol05.hoEt);
      mum_isolationR05_trackNTracks_.push_back(MuIsol05.nTracks);
      mum_isolationR05_trackNJets_.push_back(MuIsol05.nJets);
      mum_isolationR05_trackerVetoSumPt_.push_back(MuIsol05.trackerVetoPt);
      mum_isolationR05_emVetoSumEt_.push_back(MuIsol05.emVetoEt);
      mum_isolationR05_hadVetoSumEt_.push_back(MuIsol05.hadVetoEt);
      mum_isolationR05_hoVetoEt_.push_back(MuIsol05.hoVetoEt);
  
      auto &MuPFIsol = mu.pfIsolationR03();
      mum_PFIsolationR03_sumChargedHadronPt_.push_back(MuPFIsol.sumChargedHadronPt);
      mum_PFIsolationR03_sumChargedParticlePt_.push_back(MuPFIsol.sumChargedParticlePt);
      mum_PFIsolationR03_sumNeutralHadronEt_.push_back(MuPFIsol.sumNeutralHadronEt);
      mum_PFIsolationR03_sumPhotonEt_.push_back(MuPFIsol.sumPhotonEt);
      mum_PFIsolationR03_sumNeutralHadronEtHighThreshold_.push_back(MuPFIsol.sumNeutralHadronEtHighThreshold);
      mum_PFIsolationR03_sumPhotonEtHighThreshold_.push_back(MuPFIsol.sumPhotonEtHighThreshold);
      mum_PFIsolationR03_sumPUPt_.push_back(MuPFIsol.sumPUPt);
    
      auto &MuPFIsol04 = mu.pfIsolationR04();
      mum_PFIsolationR04_sumChargedHadronPt_.push_back(MuPFIsol04.sumChargedHadronPt);
      mum_PFIsolationR04_sumChargedParticlePt_.push_back(MuPFIsol04.sumChargedParticlePt);
      mum_PFIsolationR04_sumNeutralHadronEt_.push_back(MuPFIsol04.sumNeutralHadronEt);
      mum_PFIsolationR04_sumPhotonEt_.push_back(MuPFIsol04.sumPhotonEt);
      mum_PFIsolationR04_sumNeutralHadronEtHighThreshold_.push_back(MuPFIsol04.sumNeutralHadronEtHighThreshold);
      mum_PFIsolationR04_sumPhotonEtHighThreshold_.push_back(MuPFIsol04.sumPhotonEtHighThreshold);
      mum_PFIsolationR04_sumPUPt_.push_back(MuPFIsol04.sumPUPt);

      nMuM_++; 
    }
      
    //// #############
    //// # Save: mu+ #
    //// #############
    if(muTrack->charge() == 1){ 
      mupHighPurity_.push_back( (int) muTrack->quality(reco::Track::highPurity));
      mupPt_.push_back(muTrack->pt());
      mupEta_.push_back(muTrack->eta());
      mupPhi_.push_back(muTrack->phi());
      mupCL_.push_back(TMath::Prob(muTrackTT.chi2(), static_cast<int>(rint(muTrackTT.ndof()))));
      mupNormChi2_.push_back(muTrack->normalizedChi2());
      mupVx_.push_back(mu.vx());
      mupVy_.push_back(mu.vy());
      mupVz_.push_back(mu.vz());
      
      mupDCABS_.push_back(DCAmuBS);
      mupDCABSE_.push_back(DCAmuBSErr);
      mupdxyBS_.push_back(muTrackTT.track().dxy( (beamSpot.position() )));
      mupdzBS_.push_back(muTrackTT.track().dz(  (beamSpot.position() )));
      
      mupFracHits_.push_back(static_cast<double>(muTrack->hitPattern().numberOfValidHits()) / 
			     (  muTrack->hitPattern().numberOfValidHits() +
				muTrack->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
				muTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
				muTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS) ) ) ;

      
      mupIdx_.push_back(i);
      mupNPixHits_.push_back(muTrack->hitPattern().numberOfValidPixelHits());
      mupNPixLayers_.push_back(muTrack->hitPattern().pixelLayersWithMeasurement());  
      mupNTrkHits_.push_back(muTrack->hitPattern().numberOfValidTrackerHits());
      mupNTrkLayers_.push_back(muTrack->hitPattern().trackerLayersWithMeasurement());
      if (mu.isGlobalMuon() == true) mupNMuonHits_.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
      else mupNMuonHits_.push_back(0);
      mupNMatchStation_.push_back(mu.numberOfMatchedStations());

      mup_isGlobalMuon_	      .push_back(mu.isGlobalMuon());
      mup_isTrackerMuon_	      .push_back(mu.isTrackerMuon());
      mup_StandAloneMuon_          .push_back(mu.isStandAloneMuon());
      mup_isCaloMuon_	      .push_back(mu.isCaloMuon());
      mup_isPFMuon_	              .push_back(mu.isPFMuon());

      mup_selector_.push_back(mu.selectors()); 
      mup_isIsolationValid_.push_back(mu.isIsolationValid());
      mup_isPFIsolationValid_.push_back(mu.isIsolationValid());

      auto &MuIsol03 = mu.isolationR03();
      mup_isolationR03_trackSumPt_.push_back(MuIsol03.sumPt);
      mup_isolationR03_trackEcalSumEt_.push_back(MuIsol03.emEt);
      mup_isolationR03_trackHcalSumEt_.push_back(MuIsol03.hadEt);
      mup_isolationR03_trackHoSumEt_.push_back(MuIsol03.hoEt);
      mup_isolationR03_trackNTracks_.push_back(MuIsol03.nTracks);
      mup_isolationR03_trackNJets_.push_back(MuIsol03.nJets);
      mup_isolationR03_trackerVetoSumPt_.push_back(MuIsol03.trackerVetoPt);
      mup_isolationR03_emVetoSumEt_.push_back(MuIsol03.emVetoEt);
      mup_isolationR03_hadVetoSumEt_.push_back(MuIsol03.hadVetoEt);
      mup_isolationR03_hoVetoEt_.push_back(MuIsol03.hoVetoEt);
  
      auto &MuIsol05 = mu.isolationR05();
      mup_isolationR05_trackSumPt_.push_back(MuIsol05.sumPt);
      mup_isolationR05_trackEcalSumEt_.push_back(MuIsol05.emEt);
      mup_isolationR05_trackHcalSumEt_.push_back(MuIsol05.hadEt);
      mup_isolationR05_trackHoSumEt_.push_back(MuIsol05.hoEt);
      mup_isolationR05_trackNTracks_.push_back(MuIsol05.nTracks);
      mup_isolationR05_trackNJets_.push_back(MuIsol05.nJets);
      mup_isolationR05_trackerVetoSumPt_.push_back(MuIsol05.trackerVetoPt);
      mup_isolationR05_emVetoSumEt_.push_back(MuIsol05.emVetoEt);
      mup_isolationR05_hadVetoSumEt_.push_back(MuIsol05.hadVetoEt);
      mup_isolationR05_hoVetoEt_.push_back(MuIsol05.hoVetoEt);
  
      auto &MuPFIsol = mu.pfIsolationR03();
      mup_PFIsolationR03_sumChargedHadronPt_.push_back(MuPFIsol.sumChargedHadronPt);
      mup_PFIsolationR03_sumChargedParticlePt_.push_back(MuPFIsol.sumChargedParticlePt);
      mup_PFIsolationR03_sumNeutralHadronEt_.push_back(MuPFIsol.sumNeutralHadronEt);
      mup_PFIsolationR03_sumPhotonEt_.push_back(MuPFIsol.sumPhotonEt);
      mup_PFIsolationR03_sumNeutralHadronEtHighThreshold_.push_back(MuPFIsol.sumNeutralHadronEtHighThreshold);
      mup_PFIsolationR03_sumPhotonEtHighThreshold_.push_back(MuPFIsol.sumPhotonEtHighThreshold);
      mup_PFIsolationR03_sumPUPt_.push_back(MuPFIsol.sumPUPt);
    
      auto &MuPFIsol04 = mu.pfIsolationR04();
      mup_PFIsolationR04_sumChargedHadronPt_.push_back(MuPFIsol04.sumChargedHadronPt);
      mup_PFIsolationR04_sumChargedParticlePt_.push_back(MuPFIsol04.sumChargedParticlePt);
      mup_PFIsolationR04_sumNeutralHadronEt_.push_back(MuPFIsol04.sumNeutralHadronEt);
      mup_PFIsolationR04_sumPhotonEt_.push_back(MuPFIsol04.sumPhotonEt);
      mup_PFIsolationR04_sumNeutralHadronEtHighThreshold_.push_back(MuPFIsol04.sumNeutralHadronEtHighThreshold);
      mup_PFIsolationR04_sumPhotonEtHighThreshold_.push_back(MuPFIsol04.sumPhotonEtHighThreshold);
      mup_PFIsolationR04_sumPUPt_.push_back(MuPFIsol04.sumPUPt);

      nMuP_++; 
    }

    //std::cout << "1st loop:" <<  iEvent.id().event() <<  " i: " << i <<  " " << muTrack->charge() << " pt: " << muTrack->pt() << " eta:" << muTrack->eta() << std::endl;
    
    
    //start loop on second muon 
    for ( uint32_t j= i+1 ;j<muons->size();j++){

      auto &mu2 = muons->at(j);
      
      if( (mu.charge())*(mu2.charge()) == 1) continue;

      if(mu.charge() == 1){ muTrackp   =  mu.innerTrack();}
      if(mu.charge() == -1){ muTrackm  =  mu.innerTrack();}
      
      if(mu2.charge() == 1) { muTrackp  = mu2.innerTrack();}
      if(mu2.charge() == -1){ muTrackm  = mu2.innerTrack();}


      if(mu2.pt()  < pTMinMuons) continue;

      //muTrackp = mup.innerTrack();
      if ((muTrackp.isNull() == true) || (muTrackm.isNull() == true)) continue;
     
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));
      if (!muTrackpTT.isValid()) continue;

      const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle));
      if (!muTrackmTT.isValid()) continue;
      
      
      // # Check goodness of dimuons closest approach #
      ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
      XingPoint = ClosestApp.crossingPoint();
      
      double mumuDCA = ClosestApp.distance();
      
      // # dimuon inviariant mass and pT before vertex fitting #
      bsDimuon_lv.SetPxPyPzE( muTrackmTT.track().px() + muTrackpTT.track().px(), 
			      muTrackmTT.track().py() + muTrackpTT.track().py(),
			      muTrackmTT.track().pz() + muTrackpTT.track().pz(),
			      sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ) + sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) ) );
      
      
      if (true or (bsDimuon_lv.Pt() < minDimuon_pt)  || (bsDimuon_lv.M() < minDimuonInvariantMass) || (bsDimuon_lv.M() > maxDimuonInvariantMass))
	{
	  if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << bsDimuon_lv.Pt() << "\tinv. mass: " << bsDimuon_lv.M() << std::endl;
	  //        continue;
	}
      
      
      chi2 = 0.;
      ndof  = 0.;
      
      // ####################################################
      // # Try to vertex the two muons to get dimuon vertex #
      // ####################################################
      KinematicParticleFactoryFromTransientTrack partFactory;
      KinematicParticleVertexFitter PartVtxFitter;
      
      std::vector<RefCountedKinematicParticle> muonParticles;
      muonParticles.push_back(partFactory.particle(muTrackmTT, muonMass,chi2,ndof,muonMassErr));
      muonParticles.push_back(partFactory.particle(muTrackpTT, muonMass,chi2,ndof,muonMassErr));
      
      RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles); 
      if ( !mumuVertexFitTree->isValid()) continue;
      
      if (mumuVertexFitTree->isValid() == false){
	if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
	continue; 
      }
      
      mumuVertexFitTree->movePointerToTheTop();
      RefCountedKinematicVertex mumu_KV   = mumuVertexFitTree->currentDecayVertex();
      RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
      
      if ( !mumu_KV->vertexIsValid()) continue;
      
      mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
				 int(rint(mumu_KV->degreesOfFreedom())));
      
      
      // extract the re-fitted tracks
      mumuVertexFitTree->movePointerToTheTop();
      
      mumuVertexFitTree->movePointerToTheFirstChild();
      RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
      refitMumTT = refitMum->refittedTransientTrack();
      
      mumuVertexFitTree->movePointerToTheNextChild();
      RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
      refitMupTT = refitMup->refittedTransientTrack();
      
      TLorentzVector mymum, mymup, mydimu;
      
      mymum.SetXYZM(refitMumTT.track().momentum().x(),
		    refitMumTT.track().momentum().y(),
		    refitMumTT.track().momentum().z(), muonMass);
      
      mymup.SetXYZM(refitMupTT.track().momentum().x(),
		    refitMupTT.track().momentum().y(),
		    refitMupTT.track().momentum().z(), muonMass);
      
      mydimu = mymum + mymup;
      mu_mu_pt = mydimu.Perp();
      
      mu_mu_mass = mumu_KP->currentState().mass();
      mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
			    matrix()(6,6));
      
      //std::cout << mu_mu_vtx_cl << " " << mu_mu_pt << " " << mu_mu_mass << " " << mu_mu_mass_err << std::endl;
      
      // ######################################################
      // # Compute the distance between mumu vtx and BeamSpot #
      // ######################################################
      
      double MuMuLSBS;
      double MuMuLSBSErr;
      Utility->computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
			  beamSpot.position().x(),beamSpot.position().y(),0.0,
			  mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
			  mumu_KV->error().matrix()(0,1),0.0,0.0,
			  beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
			  beamSpot.covariance()(0,1),0.0,0.0,
			  &MuMuLSBS,&MuMuLSBSErr);
      
      
      
      // ###################################################################
      // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
      // ###################################################################
      double MuMuCosAlphaBS;
      double MuMuCosAlphaBSErr;
      Utility->computeCosAlpha (mumu_KP->currentState().globalMomentum().x(),
                                mumu_KP->currentState().globalMomentum().y(),
				0.0,
				mumu_KV->position().x() - beamSpot.position().x(),
				mumu_KV->position().y() - beamSpot.position().y(),
				0.0,
				mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
				mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
				0.0,
				mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
				0.0,
				0.0,
				mumu_KV->error().cxx() + beamSpot.covariance()(0,0),
				mumu_KV->error().cyy() + beamSpot.covariance()(1,1),
				0.0,
				mumu_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),
				0.0,
				0.0,
				&MuMuCosAlphaBS,&MuMuCosAlphaBSErr);
    
      
      // add loop for photons here...

       
    
      //// #################
      //// # Save: mu+ mu- #
      //// #################
      if(mu2.charge()==1) {
        mumuParentMuM_.push_back(i);
        mumuParentMuP_.push_back(j);
      }
      if(mu2.charge()==-1){
        mumuParentMuM_.push_back(j);
        mumuParentMuP_.push_back(i);

      }

      
      mumuPt_.push_back(mu_mu_pt);
      mumuEta_.push_back(mydimu.Eta());
      mumuRapidity_.push_back(mydimu.Rapidity());
      mumuPhi_.push_back(mydimu.Phi());
      mumuMass_.push_back(mu_mu_mass);
      mumuMassE_.push_back(mu_mu_mass_err);
      
      mumuPx_.push_back(mumu_KP->currentState().globalMomentum().x());
      mumuPy_.push_back(mumu_KP->currentState().globalMomentum().y());
      mumuPz_.push_back(mumu_KP->currentState().globalMomentum().z()); 
      mumuCovPxPx_.push_back(mumu_KP->currentState().kinematicParametersError().matrix()(3,3)); 
      mumuCovPxPy_.push_back(mumu_KP->currentState().kinematicParametersError().matrix()(3,4)); 
      mumuCovPxPz_.push_back(mumu_KP->currentState().kinematicParametersError().matrix()(3,5)); 
      mumuCovPyPy_.push_back(mumu_KP->currentState().kinematicParametersError().matrix()(4,4)); 
      mumuCovPyPz_.push_back(mumu_KP->currentState().kinematicParametersError().matrix()(4,5)); 
      mumuCovPzPz_.push_back(mumu_KP->currentState().kinematicParametersError().matrix()(5,5)); 
      mumuDR_.push_back(deltaR(mu,mu2));
      
      mumuVtxCL_.push_back(mu_mu_vtx_cl);
      mumuVtxChi2_.push_back(mumu_KV->chiSquared());
      mumuVtxNdof_.push_back(mumu_KV->degreesOfFreedom());
      mumuVtxX_.push_back(mumu_KV->position().x());
      mumuVtxY_.push_back(mumu_KV->position().y());
      mumuVtxZ_.push_back(mumu_KV->position().z());
      mumuVtxCovXX_.push_back(mumu_KV->error().cxx());
      mumuVtxCovXY_.push_back(mumu_KV->error().cyx());
      mumuVtxCovXZ_.push_back(mumu_KV->error().czx());
      mumuVtxCovYY_.push_back(mumu_KV->error().cyy());
      mumuVtxCovYZ_.push_back(mumu_KV->error().czy());
      mumuVtxCovZZ_.push_back(mumu_KV->error().czz());
      
      mumuCosAlphaBS_.push_back(MuMuCosAlphaBS);
      mumuCosAlphaBSE_.push_back(MuMuCosAlphaBSErr);
      mumuLBS_.push_back(MuMuLSBS);
      mumuLBSE_.push_back(MuMuLSBSErr);
      mumuDCA_.push_back(mumuDCA);
      
      // std::cout << "2nd loop after continue:" <<  iEvent.id().event() <<  " i: " << i <<  " " << mu.charge() << " pt: " << mu.pt() << " eta:" << mu.eta() << std::endl;
      // std::cout << "2nd loop cafter continue:" <<  iEvent.id().event() <<  " j: " << j <<  " " << mu2.charge() << " pt: " << mu2.pt() << " eta:" << mu2.eta() << std::endl << std::endl;
    }
  }

} //fill muons

void BsToMuMuGammaNTuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es)
{
  // Fills tree branches with photons.
  edm::Handle<std::vector<reco::Photon> > gedPhotonsHandle;
  e.getByToken(gedPhotonsCollection_, gedPhotonsHandle);

  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> recHitsEB;
  e.getByToken(ebRechitToken_, recHitsEB);
 
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> recHitsEE;
  e.getByToken(eeRechitToken_, recHitsEE);
  
  EcalClusterLazyTools       lazyTool    (e, es, ebRechitToken_, eeRechitToken_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebRechitToken_, eeRechitToken_);

  //EcalClusterLazyTools::ESGetTokens ecalClusterToolsESGetTokens_(consumesCollector());
  //EcalClusterLazyTools       lazyTool    (e,ecalClusterToolsESGetTokens_.get(es), ebRechitToken_, eeRechitToken_);
  //noZS::EcalClusterLazyTools lazyToolnoZS(e,ecalClusterToolsESGetTokens_.get(es), ebRechitToken_, eeRechitToken_);

  // loop over photons
  for (auto pho = gedPhotonsHandle->begin(); pho != gedPhotonsHandle->end(); ++pho) {
    phoE_             .push_back(pho->energy());
    phoEt_            .push_back(pho->et());
    phoEta_           .push_back(pho->eta());
    phoPhi_           .push_back(pho->phi());
    
    //phoSigmaE_          .push_back(pho->userFloat("ecalEnergyErrPostCorr"));
    //phoCalibE_          .push_back(pho->userFloat("ecalEnergyPostCorr"));
    
    //phoCalibEt_         .push_back(pho->et()*pho->userFloat("ecalEnergyPostCorr")/pho->energy());
    phoSCE_             .push_back((*pho).superCluster()->energy());
    phoSCEt_            .push_back(pho->superCluster()->energy()/cosh(pho->superCluster()->eta()));
    phoSCRawE_          .push_back((*pho).superCluster()->rawEnergy());
    phoESEnP1_          .push_back((*pho).superCluster()->preshowerEnergyPlane1());
    phoESEnP2_          .push_back((*pho).superCluster()->preshowerEnergyPlane2());
    phoSCEta_           .push_back((*pho).superCluster()->eta());
    phoSCPhi_           .push_back((*pho).superCluster()->phi());
    phoSCEtaWidth_      .push_back((*pho).superCluster()->etaWidth());
    phoSCPhiWidth_      .push_back((*pho).superCluster()->phiWidth());
    phoSCBrem_          .push_back((*pho).superCluster()->phiWidth()/(*pho).superCluster()->etaWidth());
    phohasPixelSeed_    .push_back((Int_t)pho->hasPixelSeed());
    //phoEleVeto_         .push_back((Int_t)pho->passElectronVeto());
    phoR9_              .push_back(pho->r9());
    phoHoverE_          .push_back(pho->hadTowOverEm());
    phoESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*pho).superCluster())));
    
    phoPFChIso_         .push_back(pho->chargedHadronIso()); 
    phoPFPhoIso_        .push_back(pho->photonIso());
    phoPFNeuIso_        .push_back(pho->neutralHadronIso());

    phoEcalPFClusterIso_.push_back(pho->ecalPFClusterIso());
    phoHcalPFClusterIso_.push_back(pho->hcalPFClusterIso());
    //phoIDMVA_           .push_back(pho->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));  
    
    phoSigmaIEtaIEtaFull5x5_ .push_back(pho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(pho->full5x5_showerShapeVariables().sigmaIetaIphi);
    phoSigmaIPhiIPhiFull5x5_ .push_back(pho->full5x5_showerShapeVariables().sigmaIphiIphi);
    phoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*((*pho).superCluster()->seed())));
    phoE5x5Full5x5_          .push_back(pho->full5x5_e5x5());
    phoR9Full5x5_            .push_back(pho->full5x5_r9());
    phoMIPTotEnergy_         .push_back(pho->mipTotEnergy());
    
    
    phoMIPChi2_        .push_back(pho->mipChi2());
    phoMIPSlope_       .push_back(pho->mipSlope());
    phoMIPIntercept_   .push_back(pho->mipIntercept());
    phoMIPNhitCone_    .push_back(pho->mipNhitCone());
    phoMIPIsHalo_      .push_back(pho->mipIsHalo());    
    
    ///////////////////////////////SATURATED/UNSATURATED ///from ggFlash////
    DetId seed = (pho->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
    
    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      //std::cout<<"(*theSeedHit).time()"<<(*theSeedHit).time()<<"seed energy: "<<(*theSeedHit).energy()<<std::endl;  
      
      phoSeedTime_  .push_back((*theSeedHit).time());
      phoSeedEnergy_.push_back((*theSeedHit).energy());

      
    } else{
      phoSeedTime_  .push_back(-99.);
      phoSeedEnergy_.push_back(-99.);
      
    }
 
    nPho_++;
  } // photons loop
}

void BsToMuMuGammaNTuplizer::fillPFPhotons(const edm::Event& e, const edm::EventSetup& es)
{
  // Fills tree branches with photons.
  edm::Handle<edm::View<reco::PFCandidate> > pfPhotonsHandle;
  e.getByToken(pfPhotonsCollection_, pfPhotonsHandle);

  // loop over photons
  for (auto pf = pfPhotonsHandle->begin(); pf != pfPhotonsHandle->end(); ++pf) {
    if(pf->pdgId() !=22)continue;
    if(pf->et() < pTMinPFPhotons)continue;
    phoPFE_             .push_back(pf->energy());
    phoPFEt_            .push_back(pf->et());
    phoPFEta_           .push_back(pf->eta());
    phoPFPhi_           .push_back(pf->phi());
    
    nPFPho_++;
  } // PF photons loop
}


void BsToMuMuGammaNTuplizer::fillSC(edm::Event const& e, const edm::EventSetup& es, reco::Vertex& pv) {
  edm::Handle<reco::SuperClusterCollection> barrelSCHandle;
  e.getByToken(MustacheSCBarrelCollection_, barrelSCHandle);

  edm::Handle<reco::SuperClusterCollection> endcapSCHandle;
  e.getByToken(MustacheSCEndcapCollection_, endcapSCHandle);

  edm::Handle<reco::GsfElectronCollection> gsfElectronHandle;
  e.getByToken(gsfElectronToken_, gsfElectronHandle);


  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> recHitsEB;
  //edm::Handle<EcalRecHitCollection> recHitsEB;
  e.getByToken(ebRechitToken_, recHitsEB);
  if (!recHitsEB.isValid()) {
    std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
    return;
  }

  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> recHitsEE;
  //edm::Handle<EcalRecHitCollection> recHitsEE;
  e.getByToken(eeRechitToken_, recHitsEE);
  if (!recHitsEE.isValid()) {
    std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
    return;
  } 

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> recHitsHBHE;
  e.getByToken(hbheRechitToken_, recHitsHBHE);
  if (!recHitsHBHE.isValid()) {
    std::cerr << "Analyze --> recHitsHBHE not found" << std::endl;
    return;
  }
  edm::ESHandle<CaloGeometry> theCaloGeometry;
  es.get<CaloGeometryRecord>().get(theCaloGeometry);
  
  edm::ESHandle<CaloTowerConstituentsMap> towerMap_;
  es.get<CaloGeometryRecord>().get(towerMap_);
 

  locCov_.clear();
  full5x5_locCov_.clear();
  // for(const auto& iSuperCluster : *(superClusterEB.product())){  
  for (auto const& scs : { *barrelSCHandle.product(), *endcapSCHandle.product() }) {
    //for (auto const& sc : *(barrelSCHandle.product())) {
    for (auto const& sc : scs) {
      //if(abs(sc.eta())>2.4)continue;
      //
	
      edm::ESHandle<CaloTopology> caloTopology;
      es.get<CaloTopologyRecord>().get(caloTopology);
      const CaloTopology* topology = caloTopology.product();
	
      scE_.push_back(sc.correctedEnergy());
      scEt_.push_back(sc.correctedEnergy()/cosh(sc.eta()));
      scRawE_.push_back(sc.rawEnergy());
      scRawEt_.push_back(sc.rawEnergy()/cosh(sc.eta()));
      scEta_.push_back(sc.eta());
      scPhi_.push_back(sc.phi());
      scX_.push_back(sc.x());
      scY_.push_back(sc.y());
      scZ_.push_back(sc.z());
      scEtaWidth_.push_back(sc.etaWidth());
      scPhiWidth_.push_back(sc.phiWidth());
	
      reco::CaloCluster caloBC(*sc.seed());  
      showerShapes_.clear();
      if(abs(sc.eta()) <= 1.442)showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
      if(abs(sc.eta()) >= 1.566)showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
       scE5x5_.push_back(reduceFloat(showerShapes_[0],nBits_));
       scE2x2Ratio_.push_back(reduceFloat(showerShapes_[1],nBits_));
       scE3x3Ratio_.push_back(reduceFloat(showerShapes_[2],nBits_));
       scEMaxRatio_.push_back(reduceFloat(showerShapes_[3],nBits_));
       scE2ndRatio_.push_back(reduceFloat(showerShapes_[4],nBits_));
       scETopRatio_.push_back(reduceFloat(showerShapes_[5],nBits_));
       scERightRatio_.push_back(reduceFloat(showerShapes_[6],nBits_));
       scEBottomRatio_.push_back(reduceFloat(showerShapes_[7],nBits_));
       scELeftRatio_.push_back(reduceFloat(showerShapes_[8],nBits_));
       scE2x5MaxRatio_.push_back(reduceFloat(showerShapes_[9],nBits_));
       scE2x5TopRatio_.push_back(reduceFloat(showerShapes_[10],nBits_));
       scE2x5RightRatio_.push_back(reduceFloat(showerShapes_[11],nBits_));
       scE2x5BottomRatio_.push_back(reduceFloat(showerShapes_[12],nBits_));
       scE2x5LeftRatio_.push_back(reduceFloat(showerShapes_[13],nBits_));
       scSwissCross_.push_back(reduceFloat(showerShapes_[14],nBits_));
       scR9_.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/sc.rawEnergy(),nBits_));
       scSigmaIetaIeta_.push_back(reduceFloat(showerShapes_[16],nBits_)); 
       scSigmaIetaIphi_.push_back(reduceFloat(showerShapes_[17],nBits_)); 
       scSigmaIphiIphi_.push_back(reduceFloat(showerShapes_[18],nBits_)); 
       scFull5x5_e5x5_.push_back(reduceFloat(showerShapes_[19],nBits_));
       scFull5x5_e2x2Ratio_.push_back(reduceFloat(showerShapes_[20],nBits_));
       scFull5x5_e3x3Ratio_.push_back(reduceFloat(showerShapes_[21],nBits_));
       scFull5x5_eMaxRatio_.push_back(reduceFloat(showerShapes_[22],nBits_));
       scFull5x5_e2ndRatio_.push_back(reduceFloat(showerShapes_[23],nBits_));
       scFull5x5_eTopRatio_.push_back(reduceFloat(showerShapes_[24],nBits_));
       scFull5x5_eRightRatio_.push_back(reduceFloat(showerShapes_[25],nBits_));
       scFull5x5_eBottomRatio_.push_back(reduceFloat(showerShapes_[26],nBits_));
       scFull5x5_eLeftRatio_.push_back(reduceFloat(showerShapes_[27],nBits_));
       scFull5x5_e2x5MaxRatio_.push_back(reduceFloat(showerShapes_[28],nBits_));
       scFull5x5_e2x5TopRatio_.push_back(reduceFloat(showerShapes_[29],nBits_));
       scFull5x5_e2x5RightRatio_.push_back(reduceFloat(showerShapes_[30],nBits_));
       scFull5x5_e2x5BottomRatio_.push_back(reduceFloat(showerShapes_[31],nBits_));
       scFull5x5_e2x5LeftRatio_.push_back(reduceFloat(showerShapes_[32],nBits_));
       scFull5x5_swissCross_.push_back(reduceFloat(showerShapes_[33],nBits_));
       scFull5x5_r9_.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/sc.rawEnergy(),nBits_));
       scFull5x5_sigmaIetaIeta_.push_back(reduceFloat(showerShapes_[35],nBits_)); 
       scFull5x5_sigmaIetaIphi_.push_back(reduceFloat(showerShapes_[36],nBits_)); 
       scFull5x5_sigmaIphiIphi_.push_back(reduceFloat(showerShapes_[37],nBits_)); 


      ++nSC_;
      double dRmin=1e9;
      bool foundGsfEleMatch=false;
      for(auto const& ele : *gsfElectronHandle)
	{
	  auto dr=deltaR(*(ele.superCluster()),sc);
	  dRmin = dr<dRmin ? dr : dRmin;
//	  if( &( *(ele.superCluster()) ) == &sc) 
//	    {
//	      foundGsfEleMatch=true;
//	      break;
//	    }
	}
    if(dRmin<0.01)
	{
	  foundGsfEleMatch=true;
	}

      scMinDrWithGsfElectornSC_.push_back(dRmin);
      scFoundGsfMatch_.push_back(foundGsfEleMatch);
      
      // fill H/E variable
      const reco::CaloCluster &seedCluster = *sc.seed();
      DetId seedId = seedCluster.seed() ;
     // if( seedId.det() == DetId::Forward ) return;

      CaloTowerDetId towerId(towerMap_->towerOf(seedId)); 
      int seedHcalIEta = towerId.ieta();
      int seedHcalIPhi = towerId.iphi();
    //  std::cout << e.id().event() << " Seed ID" << seedHcalIEta << std::endl;

    
     Float_t hcal2Energy(0.0);
     Int_t nhcal2Rechit_=0;
     Float_t hcalEnergy(0.0);
     nhcalRechit_=0;
     for (auto& hcalrh : e.get(hbheRechitToken_) ) {
	int dIEtaAbs = std::abs(calDIEta(seedHcalIEta, hcalrh.id().ieta()));
	int dIPhiAbs = std::abs(calDIPhi(seedHcalIPhi, hcalrh.id().iphi()));
	if ( (dIEtaAbs <= maxDIEta_) && (dIPhiAbs <= maxDIPhi_) &&  (hcalrh.energy()>getMinEnergyHCAL_(hcalrh.id()) ) ) {
      hcalEnergy+=hcalrh.energy();
      nhcalRechit_++;
      }
	if ( (dIEtaAbs <= 2) && (dIPhiAbs <= 2) &&  (hcalrh.energy()>getMinEnergyHCAL_(hcalrh.id()) ) ) {
      hcal2Energy+=hcalrh.energy();
      nhcal2Rechit_++;
      }
	} // HCAL rec hits
        
        scEFromHcalRecHitInDIEta5IPhi5.push_back(hcalEnergy);
        scNHcalRecHitInDIEta5IPhi5.push_back(float(nhcalRechit_));
        scEFromHcalRecHitInDIEta2IPhi2.push_back(hcal2Energy);
        scNHcalRecHitInDIEta2IPhi2.push_back(float(nhcal2Rechit_));
        

        pfIsoCalculator pfIso(e, pfPhotonsCollection_, pv.position());
 
        scPFChIso1_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h, 0.1, 0.02, 0.));
        scPFChIso2_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h, 0.2, 0.02, 0.));
        scPFChIso3_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h, 0.3, 0.02, 0.));
        scPFChIso4_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h, 0.4, 0.02, 0.));
        scPFChIso5_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h, 0.5, 0.02, 0.));
 
        scPFPhoIso1_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::gamma, 0.1, 0.02, 0.));
        scPFPhoIso2_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::gamma, 0.2, 0.02, 0.));
        scPFPhoIso3_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::gamma, 0.3, 0.02, 0.));
        scPFPhoIso4_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::gamma, 0.4, 0.02, 0.));
        scPFPhoIso5_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::gamma, 0.5, 0.02, 0.)); 
        
        scPFNeuIso1_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h0, 0.1, 0.0, 0.));
        scPFNeuIso2_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h0, 0.2, 0.0, 0.));
        scPFNeuIso3_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h0, 0.3, 0.0, 0.));
        scPFNeuIso4_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h0, 0.4, 0.0, 0.));
        scPFNeuIso5_          .push_back(pfIso.getPfIso(sc, reco::PFCandidate::h0, 0.5, 0.0, 0.));
    	
    } // supercluster loop
  }
}



void BsToMuMuGammaNTuplizer::SetupTriggerStorageVectors()
{
  auto numTrigs = trigTable.size();
  TrigPrescales_store = new std::vector<int> [numTrigs];
  TrigResult_store    = new std::vector<bool>[numTrigs];
}

void BsToMuMuGammaNTuplizer::ClearTrggerStorages()
{
  for(uint32_t i=0;i<trigTable.size();i++)
    {
      TrigResult_store[i].clear();
      TrigPrescales_store[i].clear();
    }
}

void BsToMuMuGammaNTuplizer::FillTrggerBranches()
{
  for(uint32_t i=0;i<trigTable.size();i++)
    {
      int foundTrig=-1;
      for(uint32_t j=0;j<trigNames.size();j++)
	{
	  if (trigNames[j].find(trigTable[i]) != std::string::npos)
	    {
	      foundTrig=j;
	    }
	}

      if(foundTrig >-1) 
	{
	  TrigResult_store[i].push_back(true);	
	  TrigPrescales_store[i].push_back(trigPrescales[foundTrig]);	
	}
      else
	{
	  TrigResult_store[i].push_back(false);	
	  TrigPrescales_store[i].push_back(-1);	
	}
    }
	
}

void BsToMuMuGammaNTuplizer::SetupTriggerBranches()
{
  std::string branchName;
  for(uint32_t i=0;i<trigTable.size();i++)
    {
      branchName=trigTable[i]+"_result";
      theTree->Branch(branchName.c_str(),&(TrigResult_store[i]));
      branchName=trigTable[i]+"_prescale";
      theTree->Branch(branchName.c_str(),&(TrigPrescales_store[i]));
    }
}

void BsToMuMuGammaNTuplizer::fillHLT(edm::Event const& iEvent)
{

  edm::Handle<edm::TriggerResults>     hltTriggerResults;
  iEvent.getByToken(triggerBits_,      hltTriggerResults);
 
  HLTConfigProvider hltConfig_;
  if (hltTriggerResults.isValid()) {
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);
    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++){
      //        std::cout<<" i = "<<itrig<<", Name  : "<<triggerNames_.triggerName(itrig)<<" \n";
      // Only consider the triggered case.                                                                                                                          
      if ((*hltTriggerResults)[itrig].accept() == 1)
	{
	  std::string triggername = triggerNames_.triggerName(itrig);
	  int triggerprescale = hltConfig_.prescaleValue(itrig, triggername);
	  // Loop over our interested HLT trigger names to find if this event contains.
	  for (unsigned int it=0; it<trigTable.size(); it++){
	    if (triggername.find(trigTable[it]) != std::string::npos) 
	      {
		// save the no versioned case
		trigNames.push_back(trigTable[it]);
		trigPrescales.push_back(triggerprescale);
	      }
	  }
	}
    }
  }
  else
    {
      std::cout<<" Trigger result Not valid !!"<<"\n";
    }
  FillTrggerBranches();
    
}

std::vector<float> BsToMuMuGammaNTuplizer::getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology)
{
  std::vector<float> shapes;
  shapes.resize(38); 
  locCov_.clear();
  full5x5_locCov_.clear();
  locCov_ = EcalClusterTools::localCovariances(*caloBC, recHits, topology);
  full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    
  float e5x5 = EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
  float e3x3 = EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
  float eMax = EcalClusterTools::eMax(*caloBC, recHits); // eMax
  float eTop = EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
  float eRight = EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
  float eBottom = EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
  float eLeft = EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
  float e4 = eTop + eRight + eBottom + eLeft;

  shapes[0] = e5x5;
  shapes[1] = EcalClusterTools::e2x2(*caloBC, recHits, topology)/e5x5; // e2x2/e5x5
  shapes[2] = EcalClusterTools::e3x3(*caloBC, recHits, topology)/e5x5; // e3x3/e5x5
  shapes[3] = EcalClusterTools::eMax(*caloBC,  recHits)/e5x5; // eMax/e5x5
  shapes[4] = EcalClusterTools::e2nd(*caloBC, recHits)/e5x5; // e2nd/e5x5
  shapes[5] = EcalClusterTools::eTop(*caloBC, recHits, topology)/e5x5; // eTop/e5x5 
  shapes[6] = EcalClusterTools::eRight(*caloBC, recHits, topology)/e5x5; // eRight/e5x5
  shapes[7] = EcalClusterTools::eBottom(*caloBC, recHits, topology)/e5x5; // eBottom/e5x5
  shapes[8] = EcalClusterTools::eLeft(*caloBC, recHits, topology)/e5x5; // eLeft/e5x5
  shapes[9] = EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/e5x5; // e2x5Max/e5x5
  shapes[10] = EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/e5x5; // e2x5Top/e5x5  
  shapes[11] = EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/e5x5; // e2x5Bottom/e5x5  
  shapes[12] = EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/e5x5; // e2x5Left/e5x5  
  shapes[13] = EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/e5x5; // e2x5Right/e5x5   
  shapes[14] = 1.-e4/eMax; // swissCross 
  shapes[15] = e3x3/caloBC->energy(); // r9     
  shapes[16] = sqrt(locCov_[0]); // sigmaIetaIeta 
  shapes[17] = locCov_[1]; // sigmaIetaIphi
  shapes[18] = !edm::isFinite(locCov_[2]) ? 0. : sqrt(locCov_[2]); // sigmaIphiIphi 

  //if(std::isnan(shapes.at(17))) std::cout << shapes[0] << " - " << shapes[1] << " - " << shapes[2] << " - " << shapes[3] << " - " << shapes[4] << " - " << shapes[5]  << " - " << shapes[6] << " - " << shapes[7] << " - " << shapes[8] << " - " << shapes[9] << " - " << shapes[10] << " - " << shapes[11] << " - " << shapes[12] << " - " << shapes[13] << " - " << shapes[14] << " - " << shapes[15]  << " - " << shapes[16] << " - " << shapes[17] << " - " << shapes[18] << std::endl;

  // full_5x5 variables
  float full5x5_e5x5 = noZS::EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
  float full5x5_e3x3 = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
  float full5x5_eMax = noZS::EcalClusterTools::eMax(*caloBC, recHits); // eMax
  float full5x5_eTop = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
  float full5x5_eRight = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
  float full5x5_eBottom = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
  float full5x5_eLeft = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
  float full5x5_e4 = full5x5_eTop + full5x5_eRight + full5x5_eBottom + full5x5_eLeft;

  shapes[19] = full5x5_e5x5;
  shapes[20] = noZS::EcalClusterTools::e2x2(*caloBC, recHits, topology)/full5x5_e5x5; // e2x2/e5x5
  shapes[21] = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology)/full5x5_e5x5; // e3x3/e5x5
  shapes[22] = noZS::EcalClusterTools::eMax(*caloBC, recHits)/full5x5_e5x5; // eMax/e5x5
  shapes[23] = noZS::EcalClusterTools::e2nd(*caloBC, recHits)/full5x5_e5x5; // e2nd/e5x5
  shapes[24] = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology)/full5x5_e5x5; // eTop/e5x5 
  shapes[25] = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology)/full5x5_e5x5; // eRight/e5x5
  shapes[26] = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology)/full5x5_e5x5; // eBottom/e5x5
  shapes[27] = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology)/full5x5_e5x5; // eLeft/e5x5
  shapes[28] = noZS::EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Max/e5x5
  shapes[29] = noZS::EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Top/e5x5  
  shapes[30] = noZS::EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Bottom/e5x5  
  shapes[31] = noZS::EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Left/e5x5  
  shapes[32] = noZS::EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Right/e5x5   
  shapes[33] = 1.-full5x5_e4/full5x5_eMax; // swissCross 
  shapes[34] = full5x5_e3x3/caloBC->energy(); // r9
  shapes[35] = sqrt(full5x5_locCov_[0]); // sigmaIetaIeta        
  shapes[36] = full5x5_locCov_[1]; // sigmaIetaIphi          
  shapes[37] = !edm::isFinite(full5x5_locCov_[2]) ? 0. : sqrt(full5x5_locCov_[2]); // sigmaIphiIphi

  for(unsigned iVar=0; iVar<shapes.size(); iVar++)
    if(std::isnan(shapes.at(iVar))) std::cout << "showerShape = " << iVar << " ---> NAN " << std::endl;  

  return shapes; 
}

void BsToMuMuGammaNTuplizer::addPrimaryVertexBranches()
{
    storageMapInt["nPrimaryVertex"]  = 0 ;
    theTree->Branch("nPrimaryVertex",   &storageMapInt["nPrimaryVertex"]);
    
    storageMapFloatArray["primaryVertex_isFake"        ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_x"             ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_y"             ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_z"             ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_t"             ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_covXX"         ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_covXY"         ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_covXZ"         ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_covYY"         ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_covYZ"         ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_covZZ"         ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_x_error"       ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_y_error"       ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_z_error"       ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_t_error"       ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_ntracks"       ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_ndof"          ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_chi2"          ] = new Float_t[N_PV_MAX];
    storageMapFloatArray["primaryVertex_normalizedChi2" ]= new Float_t[N_PV_MAX];
    
    theTree->Branch("primaryVertex_isFake",         storageMapFloatArray["primaryVertex_isFake"        ] , "primaryVertex_isFake[nPrimaryVertex]"        );
    theTree->Branch("primaryVertex_x",              storageMapFloatArray["primaryVertex_x"             ] , "primaryVertex_x[nPrimaryVertex]"             ); 
    theTree->Branch("primaryVertex_y",              storageMapFloatArray["primaryVertex_y"             ] , "primaryVertex_y[nPrimaryVertex]"             ); 
    theTree->Branch("primaryVertex_z",              storageMapFloatArray["primaryVertex_z"             ] , "primaryVertex_z[nPrimaryVertex]"             ); 
    theTree->Branch("primaryVertex_t",              storageMapFloatArray["primaryVertex_t"             ] , "primaryVertex_t[nPrimaryVertex]"             ); 
    theTree->Branch("primaryVertex_covXX",          storageMapFloatArray["primaryVertex_covXX"         ] , "primaryVertex_covXX[nPrimaryVertex]"         ); 
    theTree->Branch("primaryVertex_covXY",          storageMapFloatArray["primaryVertex_covXY"         ] , "primaryVertex_covXY[nPrimaryVertex]"         ); 
    theTree->Branch("primaryVertex_covXZ",          storageMapFloatArray["primaryVertex_covXZ"         ] , "primaryVertex_covXZ[nPrimaryVertex]"         ); 
    theTree->Branch("primaryVertex_covYY",          storageMapFloatArray["primaryVertex_covYY"         ] , "primaryVertex_covYY[nPrimaryVertex]"         ); 
    theTree->Branch("primaryVertex_covYZ",          storageMapFloatArray["primaryVertex_covYZ"         ] , "primaryVertex_covYZ[nPrimaryVertex]"         ); 
    theTree->Branch("primaryVertex_covZZ",          storageMapFloatArray["primaryVertex_covZZ"         ] , "primaryVertex_covZZ[nPrimaryVertex]"         ); 
    theTree->Branch("primaryVertex_x_error",        storageMapFloatArray["primaryVertex_x_error"       ] , "primaryVertex_x_error[nPrimaryVertex]"       ); 
    theTree->Branch("primaryVertex_y_error",        storageMapFloatArray["primaryVertex_y_error"       ] , "primaryVertex_y_error[nPrimaryVertex]"       ); 
    theTree->Branch("primaryVertex_z_error",        storageMapFloatArray["primaryVertex_z_error"       ] , "primaryVertex_z_error[nPrimaryVertex]"       ); 
    theTree->Branch("primaryVertex_t_error",        storageMapFloatArray["primaryVertex_t_error"       ] , "primaryVertex_t_error[nPrimaryVertex]"       ); 
    theTree->Branch("primaryVertex_ntracks",        storageMapFloatArray["primaryVertex_ntracks"       ] , "primaryVertex_ntracks[nPrimaryVertex]"       ); 
    theTree->Branch("primaryVertex_ndof",           storageMapFloatArray["primaryVertex_ndof"          ] , "primaryVertex_ndof[nPrimaryVertex]"          ); 
    theTree->Branch("primaryVertex_chi2",           storageMapFloatArray["primaryVertex_chi2"          ] , "primaryVertex_chi2[nPrimaryVertex]"          ); 
    theTree->Branch("primaryVertex_normalizedChi2", storageMapFloatArray["primaryVertex_normalizedChi2" ], "primaryVertex_normalizedChi2[nPrimaryVertex]");
}

void BsToMuMuGammaNTuplizer::fillPrimaryVertexBranches(const edm::Event& iEvent,  const edm::EventSetup& iSetup)
{

    edm::Handle<std::vector<reco::Vertex>> primaryVertexCollection;
    iEvent.getByToken(primaryVtxToken_, primaryVertexCollection);
    int i=0;
    storageMapInt["nPrimaryVertex"]=0;
    for(auto&  aVertex : *primaryVertexCollection){
    if( not aVertex.isValid() ) continue;
    
    // # offlinePrimaryVertices # //
    storageMapFloatArray["primaryVertex_isFake"         ][i] =   aVertex.isFake()           ;
    storageMapFloatArray["primaryVertex_x"              ][i] =   aVertex.x()                ;
    storageMapFloatArray["primaryVertex_y"              ][i] =   aVertex.y()                ;
    storageMapFloatArray["primaryVertex_z"              ][i] =   aVertex.z()                ;
    storageMapFloatArray["primaryVertex_t"              ][i] =   aVertex.t()                ;
    storageMapFloatArray["primaryVertex_covXX"          ][i] =   aVertex.covariance(0,0) ;
    storageMapFloatArray["primaryVertex_covXY"          ][i] =   aVertex.covariance(0,1) ;
    storageMapFloatArray["primaryVertex_covXZ"          ][i] =   aVertex.covariance(0,2) ;
    storageMapFloatArray["primaryVertex_covYY"          ][i] =   aVertex.covariance(1,1) ;
    storageMapFloatArray["primaryVertex_covYZ"          ][i] =   aVertex.covariance(1,2) ;
    storageMapFloatArray["primaryVertex_covZZ"          ][i] =   aVertex.covariance(2,2) ;
    storageMapFloatArray["primaryVertex_x_error"        ][i] =   aVertex.xError()        ;
    storageMapFloatArray["primaryVertex_y_error"        ][i] =   aVertex.yError()        ;
    storageMapFloatArray["primaryVertex_z_error"        ][i] =   aVertex.zError()        ;
    storageMapFloatArray["primaryVertex_t_error"        ][i] =   aVertex.tError()        ;
    storageMapFloatArray["primaryVertex_ntracks"        ][i] =   aVertex.nTracks()       ;
    storageMapFloatArray["primaryVertex_ndof"           ][i] =   aVertex.ndof() 	 	  ;
    storageMapFloatArray["primaryVertex_chi2"           ][i] =   aVertex.chi2()             ;
    storageMapFloatArray["primaryVertex_normalizedChi2" ][i] =   aVertex.normalizedChi2()  ;
    storageMapInt["nPrimaryVertex"]++;
    i++;
  } // loop over primary vertex collection
}
void BsToMuMuGammaNTuplizer::addGeneralTracksBranches()
{
    storageMapInt["nGeneralTracks"]  = 0 ;
    theTree->Branch("nGeneralTracks",   &storageMapInt["nGeneralTracks"]);
    storageMapFloatArray["generalTracks_outer_x"] = new Float_t[N_TRK_MAX];
    theTree->Branch("generalTracks_outer_x",   storageMapFloatArray["generalTracks_outer_x"],"generalTracks_outer_x[nGeneralTracks]/F");
    storageMapFloatArray["generalTracks_outer_y"] = new Float_t[N_TRK_MAX];
    theTree->Branch("generalTracks_outer_y",   storageMapFloatArray["generalTracks_outer_y"],"generalTracks_outer_y[nGeneralTracks]/F");
    storageMapFloatArray["generalTracks_outer_z"] = new Float_t[N_TRK_MAX];
    theTree->Branch("generalTracks_outer_z",   storageMapFloatArray["generalTracks_outer_z"],"generalTracks_outer_z[nGeneralTracks]/F");
    storageMapFloatArray["generalTracks_outer_px"] = new Float_t[N_TRK_MAX];
    theTree->Branch("generalTracks_outer_px",   storageMapFloatArray["generalTracks_outer_px"],"generalTracks_outer_px[nGeneralTracks]/F");
    storageMapFloatArray["generalTracks_outer_py"] = new Float_t[N_TRK_MAX];
    theTree->Branch("generalTracks_outer_py",   storageMapFloatArray["generalTracks_outer_py"],"generalTracks_outer_py[nGeneralTracks]/F");
    storageMapFloatArray["generalTracks_outer_pz"] = new Float_t[N_TRK_MAX];
    theTree->Branch("generalTracks_outer_pz",   storageMapFloatArray["generalTracks_outer_pz"],"generalTracks_outer_pz[nGeneralTracks]/F");
    storageMapFloatArray["generalTracks_charge"] = new Float_t[N_TRK_MAX];
    theTree->Branch("generalTracks_charge",   storageMapFloatArray["generalTracks_charge"],"generalTracks_charge[nGeneralTracks]/F");
}



void BsToMuMuGammaNTuplizer::fillGeneralTrackCollectionBranches( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Fills tree branches with photons.
  edm::Handle<reco::TrackCollection>  generalTracksHandle;
  iEvent.getByToken(generalTracksCollection_, generalTracksHandle);
  //std::cout<<"Size of trackColection = "<<generalTracksHandle->size()<<"\n";
  Int_t i=0;
  for (auto track = generalTracksHandle->begin();track != generalTracksHandle->end(); track++) {
    storageMapFloatArray["generalTracks_outer_x"][i]  =   track->outerX();
    storageMapFloatArray["generalTracks_outer_y"][i]  =   track->outerY();
    storageMapFloatArray["generalTracks_outer_z"][i]  =   track->outerZ();
    storageMapFloatArray["generalTracks_outer_px"][i] =   track->outerPx();
    storageMapFloatArray["generalTracks_outer_py"][i] =   track->outerPy();
    storageMapFloatArray["generalTracks_outer_pz"][i] =   track->outerPz();
    storageMapFloatArray["generalTracks_charge"][i]   =   track->charge();
    i++;
    if(i >= N_TRK_MAX) break;
  }
    
    storageMapInt["nGeneralTracks"]=i;
}

void BsToMuMuGammaNTuplizer::addParticleFlowBranches()
{
      storageMapInt["nPFCandidates"]=0;
      theTree->Branch("nPFCandidates",   &storageMapInt["nPFCandidates"]);
      storageMapFloatArray["pf_ecalEnergy"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_ecalEnergy",   storageMapFloatArray["pf_ecalEnergy"],"pf_ecalEnergy[nPFCandidates]/F");
      storageMapFloatArray["pf_ecalRawEnergy"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_ecalRawEnergy",   storageMapFloatArray["pf_ecalRawEnergy"],"pf_ecalRawEnergy[nPFCandidates]/F");
      storageMapFloatArray["pf_hcalEnergy"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_hcalEnergy",   storageMapFloatArray["pf_hcalEnergy"],"pf_hcalEnergy[nPFCandidates]/F");
      storageMapFloatArray["pf_hcalRawEnergy"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_hcalRawEnergy",   storageMapFloatArray["pf_hcalRawEnergy"],"pf_hcalRawEnergy[nPFCandidates]/F");
      storageMapFloatArray["pf_HoE"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_HoE",   storageMapFloatArray["pf_HoE"],"pf_HoE[nPFCandidates]/F");
      storageMapFloatArray["pf_mvaIso"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_mvaIso",   storageMapFloatArray["pf_mvaIso"],"pf_mvaIso[nPFCandidates]/F");
      storageMapFloatArray["pf_mvaPiMu"] = new Float_t[N_PF_MAX];
     // theTree->Branch("pf_mvaPiMu",   storageMapFloatArray["pf_mvaPiMu"],"pf_mvaPiMu[nPFCandidates]/F");
     // storageMapFloatArray["pf_mvaEMu"] = new Float_t[N_PF_MAX];
      // theTree->Branch("pf_mvaEMu",   storageMapFloatArray["pf_mvaEMu"],"pf_mvaEMu[nPFCandidates]/F");
      // storageMapFloatArray["pf_mvaEPi"] = new Float_t[N_PF_MAX];
      // theTree->Branch("pf_mvaEPi",   storageMapFloatArray["pf_mvaEPi"],"pf_mvaEPi[nPFCandidates]/F");
      // storageMapFloatArray["pf_mvaGamma"] = new Float_t[N_PF_MAX];
      // theTree->Branch("pf_mvaGamma",   storageMapFloatArray["pf_mvaGamma"],"pf_mvaGamma[nPFCandidates]/F");
      // storageMapFloatArray["pf_mvaNeuH"] = new Float_t[N_PF_MAX];
      // theTree->Branch("pf_mvaNeuH",   storageMapFloatArray["pf_mvaNeuH"],"pf_mvaNeuH[nPFCandidates]/F");
      // storageMapFloatArray["pf_mvaNeuHGamma"] = new Float_t[N_PF_MAX];
      // theTree->Branch("pf_mvaNeuHGamma",   storageMapFloatArray["pf_mvaNeuHGamma"],"pf_mvaNeuHGamma[nPFCandidates]/F");
     // storageMapFloatArray["pf_dnnGamma"] = new Float_t[N_PF_MAX];
     // theTree->Branch("pf_dnnGamma",   storageMapFloatArray["pf_dnnGamma"],"pf_dnnGamma[nPFCandidates]/F");
     // storageMapFloatArray["pf_dnnEeleBkgGamma"] = new Float_t[N_PF_MAX];
     // theTree->Branch("pf_dnnEeleBkgGamma",   storageMapFloatArray["pf_dnnEeleBkgGamma"],"pf_dnnEeleBkgGamma[nPFCandidates]/F");
     // storageMapFloatArray["pf_dnnEeleSigIso"] = new Float_t[N_PF_MAX];
     // theTree->Branch("pf_dnnEeleSigIso",   storageMapFloatArray["pf_dnnEeleSigIso"],"pf_dnnEeleSigIso[nPFCandidates]/F");
      storageMapFloatArray["pf_vertexX"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_vertexX",   storageMapFloatArray["pf_vertexX"],"pf_vertexX[nPFCandidates]/F");
      storageMapFloatArray["pf_vertexY"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_vertexY",   storageMapFloatArray["pf_vertexY"],"pf_vertexY[nPFCandidates]/F");
      storageMapFloatArray["pf_vertexZ"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_vertexZ",   storageMapFloatArray["pf_vertexZ"],"pf_vertexZ[nPFCandidates]/F");
      storageMapFloatArray["pf_ecalEntryEta"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_ecalEntryEta",   storageMapFloatArray["pf_ecalEntryEta"],"pf_ecalEntryEta[nPFCandidates]/F");
      storageMapFloatArray["pf_ecalEntryPhi"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_ecalEntryPhi",   storageMapFloatArray["pf_ecalEntryPhi"],"pf_ecalEntryPhi[nPFCandidates]/F");
      storageMapFloatArray["pf_ecalEntryX"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_ecalEntryX",   storageMapFloatArray["pf_ecalEntryX"],"pf_ecalEntryX[nPFCandidates]/F");
      storageMapFloatArray["pf_ecalEntryY"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_ecalEntryY",   storageMapFloatArray["pf_ecalEntryY"],"pf_ecalEntryY[nPFCandidates]/F");
      storageMapFloatArray["pf_ecalEntryZ"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_ecalEntryZ",   storageMapFloatArray["pf_ecalEntryZ"],"pf_ecalEntryZ[nPFCandidates]/F");
      storageMapFloatArray["pf_eta"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_eta",   storageMapFloatArray["pf_eta"],"pf_eta[nPFCandidates]/F");
      storageMapFloatArray["pf_phi"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_phi",   storageMapFloatArray["pf_phi"],"pf_phi[nPFCandidates]/F");
      storageMapFloatArray["pf_pt"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_pt",   storageMapFloatArray["pf_pt"],"pf_pt[nPFCandidates]/F");
      storageMapFloatArray["pf_id"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_id",   storageMapFloatArray["pf_id"],"pf_id[nPFCandidates]/F");
      storageMapFloatArray["pf_mass"] = new Float_t[N_PF_MAX];
      theTree->Branch("pf_mass",   storageMapFloatArray["pf_mass"],"pf_mass[nPFCandidates]/F");
}

void BsToMuMuGammaNTuplizer::fillPFCandiateCollection( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  edm::Handle<reco::PFCandidateCollection>  pfCandidateHandle;
  iEvent.getByToken(pfCandidateCollection_, pfCandidateHandle);
  //std::cout<<"Size of pfColection = "<<pfCandidateHandle->size()<<"\n";
  Int_t i=0;
  storageMapInt["nPFCandidates"]=0;
  for (auto pfCd = pfCandidateHandle->begin();pfCd != pfCandidateHandle->end(); pfCd++) 
  {
      storageMapFloatArray["pf_ecalEnergy"]     [i]= pfCd->ecalEnergy();   
      storageMapFloatArray["pf_ecalRawEnergy"]  [i]= pfCd->rawEcalEnergy(); 
      storageMapFloatArray["pf_hcalEnergy"]     [i]= pfCd->hcalEnergy(); 
      storageMapFloatArray["pf_hcalRawEnergy"]  [i]= pfCd->rawHcalEnergy(); 
      storageMapFloatArray["pf_HoE"]            [i]= pfCd->hcalEnergy()/(1e-6 +  pfCd->ecalEnergy() ); 
      storageMapFloatArray["pf_mvaIso"]         [i]= pfCd->mva_Isolated(); 
    //  storageMapFloatArray["pf_mvaPiMu"]        [i]= pfCd->mva_pi_mu(); 
    //  storageMapFloatArray["pf_mvaEMu"]         [i]= pfCd->mva_e_mu(); 
    //  storageMapFloatArray["pf_mvaEPi"]         [i]= pfCd->mva_e_pi(); 
    //  storageMapFloatArray["pf_mvaGamma"]       [i]= pfCd->mva_nothing_gamma(); 
    //  storageMapFloatArray["pf_mvaNeuH"]        [i]= pfCd->mva_nothing_nh(); 
    //  storageMapFloatArray["pf_mvaNeuHGamma"]   [i]= pfCd->mva_gamma_nh(); 
    //  storageMapFloatArray["pf_dnnGamma"]       [i]= pfCd->dnn_gamma(); 
    //  storageMapFloatArray["pf_dnnEeleBkgGamma"][i]= pfCd->dnn_e_bkgPhoton();
    //  storageMapFloatArray["pf_dnnEeleSigNonIso"]  [i]= pfCd->dnn_e_sigNonIsolated(); 
      storageMapFloatArray["pf_vertexX"]    [i] = pfCd->vx();  
      storageMapFloatArray["pf_vertexY"]    [i] = pfCd->vy();
      storageMapFloatArray["pf_vertexZ"]    [i] = pfCd->vz();
      storageMapFloatArray["pf_ecalEntryEta"] [i] = (pfCd->positionAtECALEntrance()).Eta();
      storageMapFloatArray["pf_ecalEntryPhi"] [i] = (pfCd->positionAtECALEntrance()).Phi();
      storageMapFloatArray["pf_ecalEntryX"] [i] = (pfCd->positionAtECALEntrance()).X();
      storageMapFloatArray["pf_ecalEntryY"] [i] = (pfCd->positionAtECALEntrance()).Y();
      storageMapFloatArray["pf_ecalEntryZ"] [i] = (pfCd->positionAtECALEntrance()).Z();
      storageMapFloatArray["pf_eta"]        [i] = pfCd->eta();
      storageMapFloatArray["pf_phi"]        [i] = pfCd->phi();
      storageMapFloatArray["pf_pt"]         [i] = pfCd->pt();
      storageMapFloatArray["pf_id"]         [i] = pfCd->particleId();
      storageMapFloatArray["pf_mass"]       [i] = pfCd->mass();
      i++;
    if(i >= N_PF_MAX) break;
  }
  storageMapInt["nPFCandidates"]=i;
}
void BsToMuMuGammaNTuplizer::addECALCluserBranches()
{
      storageMapInt["nECALClusters"]=0;
      theTree->Branch("nECALClusters",   &storageMapInt["nECALClusters"]);
      storageMapInt["nECALClusterEnergyMatrix"]=0;
      theTree->Branch("nECALClusterEnergyMatrix",   &storageMapInt["nECALClusterEnergyMatrix"]);
      
      storageMapFloatArray["clusterECAL_energy"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_energy",   storageMapFloatArray["clusterECAL_energy"],"clusterECAL_energy[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_correctedEnergy"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_correctedEnergy",   storageMapFloatArray["clusterECAL_correctedEnergy"],"clusterECAL_correctedEnergy[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_time"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_time",   storageMapFloatArray["clusterECAL_time"],"clusterECAL_time[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_eta"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_eta",   storageMapFloatArray["clusterECAL_eta"],"clusterECAL_eta[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_phi"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_phi",   storageMapFloatArray["clusterECAL_phi"],"clusterECAL_phi[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_x"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_x",   storageMapFloatArray["clusterECAL_x"],"clusterECAL_x[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_y"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_y",   storageMapFloatArray["clusterECAL_y"],"clusterECAL_y[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_z"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_z",   storageMapFloatArray["clusterECAL_z"],"clusterECAL_z[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_pt"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_pt",   storageMapFloatArray["clusterECAL_pt"],"clusterECAL_pt[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_size"] = new Float_t[N_ECAL_CLUSTERS];
      theTree->Branch("clusterECAL_size",   storageMapFloatArray["clusterECAL_size"],"clusterECAL_size[nECALClusters]/F");
      storageMapFloatArray["clusterECAL_energyMatrix"] = new Float_t[N_ECAL_CLUSTERS*energyMatrixSizeFull_];
      theTree->Branch("clusterECAL_energyMatrix",   storageMapFloatArray["clusterECAL_energyMatrix"],"clusterECAL_energyMatrix[nECALClusterEnergyMatrix]/F");
}


void BsToMuMuGammaNTuplizer::fillECALClusterCollection( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::PFClusterCollection>  ecalClusterHandle;
  iEvent.getByToken(ecalClusterCollection_, ecalClusterHandle);
  //std::cout<<"Size of ECAL Cluster Collection  = "<<ecalClusterHandle->size()<<"\n";
  
  noZS::EcalClusterLazyTools lazyTool(iEvent,  iSetup, ebRechitToken_, eeRechitToken_);
  
  //std::unique_ptr<noZS::EcalClusterLazyTools> lazyTools;
  //lazyTools = std::make_unique<noZS::EcalClusterLazyTools>(iEvent, ecalClusterToolsESGetTokens_.get(iSetup), ebRechitToken_, eeRechitToken_);

  //EcalClusterLazyTools       lazyTool    (e, es, ebRechitToken_, eeRechitToken_);
  //noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebRechitToken_, eeRechitToken_);
  
  Int_t i=0,energySize=0;
  storageMapInt["nECALClusters"]=0;
  for (auto ecalClus = ecalClusterHandle->begin(); ecalClus != ecalClusterHandle->end(); ecalClus++) 
  {
      storageMapFloatArray["clusterECAL_energy"]          [i]  =  ecalClus->energy() ; 
      storageMapFloatArray["clusterECAL_correctedEnergy"] [i]  =  ecalClus->correctedEnergy() ; 
      storageMapFloatArray["clusterECAL_time"]            [i]  =  ecalClus->time()   ; 
      storageMapFloatArray["clusterECAL_x"]               [i]  =  ecalClus->x()      ; 
      storageMapFloatArray["clusterECAL_y"]               [i]  =  ecalClus->y()      ; 
      storageMapFloatArray["clusterECAL_z"]               [i]  =  ecalClus->z()      ; 
      storageMapFloatArray["clusterECAL_pt"]              [i]  =  ecalClus->pt()     ; 
      storageMapFloatArray["clusterECAL_eta"]              [i]  =  ecalClus->eta()     ; 
      storageMapFloatArray["clusterECAL_phi"]              [i]  =  ecalClus->phi()     ; 
      storageMapFloatArray["clusterECAL_size"]            [i]  =  ecalClus->size()   ; 
      auto energyMatrix_ = lazyTool.energyMatrix(*ecalClus , energyMatrixSize_)  ;
      for( int j=0; j < energyMatrixSizeFull_; j++)
      {
            storageMapFloatArray["clusterECAL_energyMatrix"][ i*energyMatrixSizeFull_+ j] = energyMatrix_[j];
            energySize++;
      }
    i++;
    if(i >= N_ECAL_CLUSTERS) break;
  }
  storageMapInt["nECALClusters"]=i;
  storageMapInt["nECALClusterEnergyMatrix"]=energySize;

}

void BsToMuMuGammaNTuplizer::addHCALCluserBranches()
{
      storageMapInt["nHCALClusters"]=0;
      theTree->Branch("nHCALClusters",   &storageMapInt["nHCALClusters"]);
      storageMapInt["nHCALClusterEnergyMatrix"]=0;
      theTree->Branch("nHCALClusterEnergyMatrix",   &storageMapInt["nHCALClusterEnergyMatrix"]);
      storageMapFloatArray["clusterHCAL_energy"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_energy",   storageMapFloatArray["clusterHCAL_energy"],"clusterHCAL_energy[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_correctedEnergy"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_correctedEnergy",   storageMapFloatArray["clusterHCAL_correctedEnergy"],"clusterHCAL_correctedEnergy[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_time"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_time",   storageMapFloatArray["clusterHCAL_time"],"clusterHCAL_time[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_x"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_x",   storageMapFloatArray["clusterHCAL_x"],"clusterHCAL_x[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_y"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_y",   storageMapFloatArray["clusterHCAL_y"],"clusterHCAL_y[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_z"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_z",   storageMapFloatArray["clusterHCAL_z"],"clusterHCAL_z[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_pt"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_pt",   storageMapFloatArray["clusterHCAL_pt"],"clusterHCAL_pt[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_eta"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_eta",   storageMapFloatArray["clusterHCAL_eta"],"clusterHCAL_eta[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_phi"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_phi",   storageMapFloatArray["clusterHCAL_phi"],"clusterHCAL_phi[nHCALClusters]/F");
      storageMapFloatArray["clusterHCAL_size"] = new Float_t[N_HCAL_CLUSTERS];
      theTree->Branch("clusterHCAL_size",   storageMapFloatArray["clusterHCAL_size"],"clusterHCAL_size[nHCALClusters]/F");
}

void BsToMuMuGammaNTuplizer::fillHCALClusterCollection( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     edm::Handle<reco::PFClusterCollection>  hcalClusterHandle;
     iEvent.getByToken(hcalClusterCollection_, hcalClusterHandle);
     //std::cout<<"Size of HCAL Cluster Collection  = "<<hcalClusterHandle->size()<<"\n";
     
     //std::unique_ptr<noZS::EcalClusterLazyTools> lazyTools;
     //lazyTools = std::make_unique<noZS::EcalClusterLazyTools>(iEvent, hcalClusterToolsESGetTokens_.get(iSetup), ebRechitToken_, eeRechitToken_);
   
     //EcalClusterLazyTools       lazyTool    (e, es, ebRechitToken_, eeRechitToken_);
     //noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebRechitToken_, eeRechitToken_);
     Int_t i=0;
     for (auto hcalClus = hcalClusterHandle->begin(); hcalClus != hcalClusterHandle->end(); hcalClus++) 
     {
         storageMapFloatArray["clusterHCAL_energy"]          [i]  =  hcalClus->energy() ; 
         storageMapFloatArray["clusterHCAL_correctedEnergy"] [i]  =  hcalClus->correctedEnergy() ; 
         storageMapFloatArray["clusterHCAL_time"]            [i]  =  hcalClus->time()   ; 
         storageMapFloatArray["clusterHCAL_x"]               [i]  =  hcalClus->x()      ; 
         storageMapFloatArray["clusterHCAL_y"]               [i]  =  hcalClus->y()      ; 
         storageMapFloatArray["clusterHCAL_z"]               [i]  =  hcalClus->z()      ; 
         storageMapFloatArray["clusterHCAL_pt"]              [i]  =  hcalClus->pt()     ; 
         storageMapFloatArray["clusterHCAL_eta"]              [i]  =  hcalClus->eta()     ; 
         storageMapFloatArray["clusterHCAL_phi"]              [i]  =  hcalClus->phi()     ; 
         storageMapFloatArray["clusterHCAL_size"]            [i]  =  hcalClus->size()   ; 
         i++;
       if(i >= N_HCAL_CLUSTERS) break;
     }
     storageMapInt["nHCALClusters"]=i;
}

float BsToMuMuGammaNTuplizer::reduceFloat(float val, int bits)
{
  if(!doCompression_) return val;
  else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


int BsToMuMuGammaNTuplizer::calDIPhi(int iPhi1, int iPhi2) {
  int dPhi = iPhi1 - iPhi2;
  if (dPhi > 72 / 2)
    dPhi -= 72;
  else if (dPhi < -72 / 2)
    dPhi += 72;
  return dPhi;
}

int BsToMuMuGammaNTuplizer::calDIEta(int iEta1, int iEta2) {
  int dEta = iEta1 - iEta2;
  if (iEta1 * iEta2 < 0) {  //-ve to +ve transistion and no crystal at zero
    if (dEta < 0)
      dEta++;
    else
      dEta--;
  }
  return dEta;
}

float BsToMuMuGammaNTuplizer::getMinEnergyHCAL_(HcalDetId id) const {
  if ( (id.subdetId() == HcalBarrel)  ) {
    if ( (Run2_2018_ == 1) )
      return 0.7;
    else if ( (Run2_2018_ == 0) ) { //means Run3
      if (id.depth() == 1)
	return 0.1;
      else if (id.depth() == 2)
	return 0.2;
      else
	return 0.3;
    }
    else //neither 2018 , nor Run3, not supported
      return 9999.0;
  } 

  else if (id.subdetId() == HcalEndcap) {
    if (id.depth() == 1)
      return 0.1;
    else
      return 0.2;
  } else
    return 9999.0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(BsToMuMuGammaNTuplizer);
