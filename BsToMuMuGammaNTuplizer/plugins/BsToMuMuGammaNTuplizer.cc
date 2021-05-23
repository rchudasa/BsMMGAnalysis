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
  doHLT(iConfig.getParameter<bool>("doHLT"))

  //genParticleSrc_(iConfig.getUntrackedParameter<edm::InputTag>("genParticleSrc")),
  //gedPhotonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("gedPhotonSrc")),
  //pfPhotonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("pfPhotonSrc")),
  //MustacheSCBarrelSrc_(iConfig.getParameter<edm::InputTag>("MustacheSCBarrelSrc")),
  //MustacheSCEndcapSrc_(iConfig.getParameter<edm::InputTag>("MustacheSCEndcapSrc"))
{
  
  
  Utility= new Utils();
  
  if(doMuons_) muonToken_              = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  //caloPartToken_                 = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
  
  if(doPhotons_)    gedPhotonsCollection_       = consumes<std::vector<reco::Photon>>(iConfig.getUntrackedParameter<edm::InputTag>("gedPhotonSrc"));
  
  if(doPFPhotons_){
    pfPhotonsCollection_        = consumes<std::vector<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("pfPhotonSrc"));
   // pfRecHitToken_                 = consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag>("pfRechitCollection"));
   // pfClusterToken_                = consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection"));
  }  
  
  if(doSuperClusters_){
    MustacheSCBarrelCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCBarrelSrc"));
    MustacheSCEndcapCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCEndcapSrc"));
    gsfElectronToken_                       = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("GsfElectronSrc"));
    ebRechitToken_                 = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
    eeRechitToken_                 = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
    //ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
    //eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
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
  theTree->Branch("gen_Bs_eta"			,&gen_Bs_eta_);
  theTree->Branch("gen_Bs_phi"			,&gen_Bs_phi_);
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
  theTree->Branch("gen_BsPhoton_eta"		,&gen_BsPhoton_eta_);
  theTree->Branch("gen_BsPhoton_phi"		,&gen_BsPhoton_phi_);
 
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
  theTree->Branch("mumuDR",    &mumuDR_);

  // ### mu+ mu- Vtx ###
  theTree->Branch("mumuVtxCL",       &mumuVtxCL_);
  theTree->Branch("mumuVtxX",        &mumuVtxX_);
  theTree->Branch("mumuVtxY",        &mumuVtxY_);
  theTree->Branch("mumuVtxZ",        &mumuVtxZ_);
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
  theTree->Branch("scRawE",               &scRawE_);
  theTree->Branch("scEta",                &scEta_);
  theTree->Branch("scPhi",                &scPhi_);
  theTree->Branch("scX",        &scX_);
  theTree->Branch("scY",        &scY_);
  theTree->Branch("scZ",        &scZ_);
  theTree->Branch("scEtaWidth", &scEtaWidth_);
  theTree->Branch("scPhiWidth", &scPhiWidth_);
  theTree->Branch("scRawE",     &scRawE_);
  theTree->Branch("scRawEt",    &scRawEt_);
  theTree->Branch("scMinDrWithGsfElectornSC_",  &scMinDrWithGsfElectornSC_);
  theTree->Branch("scFoundGsfMatch_" ,        &scFoundGsfMatch_);

      theTree->Branch("superCluster_e5x5",   &superCluster_e5x5_);
      theTree->Branch("superCluster_e2x2Ratio",   &superCluster_e2x2Ratio_);
      theTree->Branch("superCluster_e3x3Ratio",   &superCluster_e3x3Ratio_);
      theTree->Branch("superCluster_eMaxRatio",   &superCluster_eMaxRatio_);
      theTree->Branch("superCluster_e2ndRatio",   &superCluster_e2ndRatio_);
      theTree->Branch("superCluster_eTopRatio",   &superCluster_eTopRatio_);
      theTree->Branch("superCluster_eRightRatio",   &superCluster_eRightRatio_);
      theTree->Branch("superCluster_eBottomRatio",   &superCluster_eBottomRatio_);
      theTree->Branch("superCluster_eLeftRatio",   &superCluster_eLeftRatio_);
      theTree->Branch("superCluster_e2x5MaxRatio",   &superCluster_e2x5MaxRatio_);
      theTree->Branch("superCluster_e2x5TopRatio",   &superCluster_e2x5TopRatio_);
      theTree->Branch("superCluster_e2x5RightRatio",   &superCluster_e2x5RightRatio_);
      theTree->Branch("superCluster_e2x5BottomRatio",   &superCluster_e2x5BottomRatio_); 
      theTree->Branch("superCluster_e2x5LeftRatio",   &superCluster_e2x5LeftRatio_); 
      theTree->Branch("superCluster_swissCross",   &superCluster_swissCross_); 
      theTree->Branch("superCluster_r9",   &superCluster_r9_);
      theTree->Branch("superCluster_sigmaIetaIeta",   &superCluster_sigmaIetaIeta_);
      theTree->Branch("superCluster_sigmaIetaIphi",   &superCluster_sigmaIetaIphi_);
      theTree->Branch("superCluster_sigmaIphiIphi",   &superCluster_sigmaIphiIphi_);
      theTree->Branch("superCluster_full5x5_e5x5",   &superCluster_full5x5_e5x5_);
      theTree->Branch("superCluster_full5x5_e2x2Ratio",   &superCluster_full5x5_e2x2Ratio_);
      theTree->Branch("superCluster_full5x5_e3x3Ratio",   &superCluster_full5x5_e3x3Ratio_);
      theTree->Branch("superCluster_full5x5_eMaxRatio",   &superCluster_full5x5_eMaxRatio_);
      theTree->Branch("superCluster_full5x5_e2ndRatio",   &superCluster_full5x5_e2ndRatio_);
      theTree->Branch("superCluster_full5x5_eTopRatio",   &superCluster_full5x5_eTopRatio_);
      theTree->Branch("superCluster_full5x5_eRightRatio",   &superCluster_full5x5_eRightRatio_);
      theTree->Branch("superCluster_full5x5_eBottomRatio",   &superCluster_full5x5_eBottomRatio_);
      theTree->Branch("superCluster_full5x5_eLeftRatio",   &superCluster_full5x5_eLeftRatio_);
      theTree->Branch("superCluster_full5x5_e2x5MaxRatio",   &superCluster_full5x5_e2x5MaxRatio_);
      theTree->Branch("superCluster_full5x5_e2x5TopRatio",   &superCluster_full5x5_e2x5TopRatio_);
      theTree->Branch("superCluster_full5x5_e2x5RightRatio",   &superCluster_full5x5_e2x5RightRatio_);
      theTree->Branch("superCluster_full5x5_e2x5BottomRatio",   &superCluster_full5x5_e2x5BottomRatio_); 
      theTree->Branch("superCluster_full5x5_e2x5LeftRatio",   &superCluster_full5x5_e2x5LeftRatio_); 
      theTree->Branch("superCluster_full5x5_swissCross",   &superCluster_full5x5_swissCross_); 
      theTree->Branch("superCluster_full5x5_r9",   &superCluster_full5x5_r9_);
      theTree->Branch("superCluster_full5x5_sigmaIetaIeta",   &superCluster_full5x5_sigmaIetaIeta_);
      theTree->Branch("superCluster_full5x5_sigmaIetaIphi",   &superCluster_full5x5_sigmaIetaIphi_);
      theTree->Branch("superCluster_full5x5_sigmaIphiIphi",   &superCluster_full5x5_sigmaIphiIphi_);
  }
}

//BsToMuMuGammaNTuplizer::~BsToMuMuGammaNTuplizer(){
//}


// ------------ method called for each event  ------------
void 
BsToMuMuGammaNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

      // ## BEAMSOPT STUFF  ## //
   beamspot_x  = 0.0   ;
   beamspot_y  = 0.0   ;
   beamspot_z  = 0.0   ;
   beamspot_x_error  = 0.0   ;
   beamspot_y_error  = 0.0   ;
   beamspot_z_error  = 0.0   ;
   beamspot_dxdz  = 0.0   ;
   beamspot_dydz  = 0.0   ;
   beamspot_sigmaZ  = 0.0   ;
   beamspot_dxdz_error  = 0.0   ;
   beamspot_dydz_error  = 0.0   ;
   beamspot_sigmaZError  = 0.0   ;
   beamspot_beamWidthX  = 0.0   ;
   beamspot_beamWidthY  = 0.0   ;
   beamspot_beamWidthX_error  = 0.0   ;
   beamspot_beamWidthY_error  = 0.0   ;

  // # offlinePrimaryVertices # //
  
  primaryVertex_isFake.clear();
  primaryVertex_x.clear();
  primaryVertex_y.clear();
  primaryVertex_z.clear();
  primaryVertex_t.clear();
  primaryVertex_x_error.clear();
  primaryVertex_y_error.clear();
  primaryVertex_z_error.clear();
  primaryVertex_t_error.clear();
  primaryVertex_ntracks.clear();
  primaryVertex_ndof.clear();
  primaryVertex_chi2.clear();
  primaryVertex_normalizedChi2.clear();

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
  gen_BsPhoton_eta_.clear() ;
  gen_BsPhoton_phi_.clear();
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
    mumuDR_.clear();
  
  
    mumuVtxCL_.clear();
    mumuVtxX_.clear();
    mumuVtxY_.clear();
    mumuVtxZ_.clear();  
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
    scRawE_               .clear();
    scEta_                .clear();
    scPhi_                .clear();
  scX_          .clear();      
  scY_          .clear();      
  scZ_          .clear();      
  scEtaWidth_   .clear();         
  scPhiWidth_   .clear();         
  scRawE_       .clear();         
  scRawEt_      .clear();   
  scMinDrWithGsfElectornSC_.clear();
  scFoundGsfMatch_.clear();

   superCluster_e5x5_.clear();
   superCluster_e2x2Ratio_.clear();
   superCluster_e3x3Ratio_.clear();
   superCluster_eMaxRatio_.clear();
   superCluster_e2ndRatio_.clear();
   superCluster_eTopRatio_.clear();
   superCluster_eRightRatio_.clear();
   superCluster_eBottomRatio_.clear();
   superCluster_eLeftRatio_.clear();
   superCluster_e2x5MaxRatio_.clear();
   superCluster_e2x5TopRatio_.clear();
   superCluster_e2x5RightRatio_.clear();
   superCluster_e2x5BottomRatio_.clear();
   superCluster_e2x5LeftRatio_.clear();
   superCluster_swissCross_.clear();
   superCluster_r9_.clear();
   superCluster_sigmaIetaIeta_.clear(); 
   superCluster_sigmaIetaIphi_.clear(); 
   superCluster_sigmaIphiIphi_.clear(); 
   superCluster_full5x5_e5x5_.clear();
   superCluster_full5x5_e2x2Ratio_.clear();
   superCluster_full5x5_e3x3Ratio_.clear();
   superCluster_full5x5_eMaxRatio_.clear();
   superCluster_full5x5_e2ndRatio_.clear();
   superCluster_full5x5_eTopRatio_.clear();
   superCluster_full5x5_eRightRatio_.clear();
   superCluster_full5x5_eBottomRatio_.clear();
   superCluster_full5x5_eLeftRatio_.clear();
   superCluster_full5x5_e2x5MaxRatio_.clear();
   superCluster_full5x5_e2x5TopRatio_.clear();
   superCluster_full5x5_e2x5RightRatio_.clear();
   superCluster_full5x5_e2x5BottomRatio_.clear();
   superCluster_full5x5_e2x5LeftRatio_.clear();
   superCluster_full5x5_swissCross_.clear();
   superCluster_full5x5_r9_.clear();
   superCluster_full5x5_sigmaIetaIeta_.clear(); 
   superCluster_full5x5_sigmaIetaIphi_.clear(); 
   superCluster_full5x5_sigmaIphiIphi_.clear();  
  }

  run_    = iEvent.id().run();
  event_  = iEvent.id().event();
  lumis_  = iEvent.luminosityBlock();
  isData_ = iEvent.isRealData();

// Get magnetic field
    
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);      
  //  Get BeamSpot
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken_, beamSpotH);
  reco::BeamSpot beamSpot = *beamSpotH;

 
  edm::Handle<std::vector<reco::Vertex>> primaryVertexCollection;
  iEvent.getByToken(primaryVtxToken_, primaryVertexCollection);
 
  // adding  BEAMSOPT 
  beamspot_x			= beamSpot.x0();  ;
  beamspot_y			= beamSpot.y0();  ;
  beamspot_z			= beamSpot.z0();  ;
  beamspot_x_error		= beamSpot.x0Error();  ;
  beamspot_y_error		= beamSpot.y0Error();  ;
  beamspot_z_error		= beamSpot.z0Error();  ;
  beamspot_dxdz   		= beamSpot.dxdz();  ;
  beamspot_dydz	         	= beamSpot.dydz();  ;
  beamspot_sigmaZ		= beamSpot.sigmaZ();  ;
  beamspot_dxdz_error		= beamSpot.dxdzError();  ;
  beamspot_dydz_error		= beamSpot.dydzError();  ;
  beamspot_sigmaZError		= beamSpot.sigmaZ0Error();  ;
  beamspot_beamWidthX		= beamSpot.BeamWidthX();  ;
  beamspot_beamWidthY		= beamSpot.BeamWidthY();  ;
  beamspot_beamWidthX_error	= beamSpot.BeamWidthXError();  ;
  beamspot_beamWidthY_error	= beamSpot.BeamWidthXError();  ;
 
  for(auto&  aVertex : *primaryVertexCollection){

    if( not aVertex.isValid() ) continue;
    
    // # offlinePrimaryVertices # //
    primaryVertex_isFake .push_back(   aVertex.isFake() );
    primaryVertex_x .push_back(   aVertex.x() );
    primaryVertex_y .push_back(   aVertex.y()  );
    primaryVertex_z .push_back(   aVertex.z()  );
    primaryVertex_t .push_back(   aVertex.t()  );
    primaryVertex_x_error .push_back(   aVertex.xError()  );
    primaryVertex_y_error .push_back(   aVertex.yError() );
    primaryVertex_z_error .push_back(   aVertex.zError()  );
    primaryVertex_t_error .push_back(   aVertex.tError()  );
    primaryVertex_ntracks .push_back(   aVertex.nTracks() );
    primaryVertex_ndof .push_back(   aVertex.ndof() 	 	  );
    primaryVertex_chi2 .push_back(   aVertex.chi2()  );
    primaryVertex_normalizedChi2 .push_back(   aVertex.normalizedChi2()  );
  } // loop over primary vertex collection

  // MC truth
  if (isMC) {
    fillGenParticles(iEvent);
  }

  if (doHLT) 		fillHLT(iEvent);
  if (doMuons_)     	fillMuons(iEvent, iSetup);
  if (doPhotons_)    	fillPhotons(iEvent, iSetup);
  if (doPFPhotons_) 	fillPFPhotons(iEvent, iSetup);
  if (doSuperClusters_) fillSC(iEvent, iSetup);
  theTree->Fill();
  
}



void BsToMuMuGammaNTuplizer::fillGenParticles(const edm::Event& iEvent)
{
  
  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  iEvent.getByToken(genParticlesCollection_, genParticleCollection);
  
  int phoMul(-1),muMMul(-1),muPMul(-1);
  
  
  for(auto& aBsMeson : *genParticleCollection){
    
    if(abs(aBsMeson.pdgId())!=531) continue;
    
    for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
      
      auto& bsDaughter = *(aBsMeson.daughter(j));
      if(bsDaughter.pdgId() == -13) muMMul++;
      if(bsDaughter.pdgId() ==  13) muPMul++;
      if(bsDaughter.pdgId() ==  22) phoMul++;
    }
    if(muMMul <0 or muPMul <0 ) continue;
    if(phoMul <0 and doBsToMuMuGamma  ) continue;

    for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
      auto& bsDaughter = *(aBsMeson.daughter(j));
      if(bsDaughter.pdgId() == -13) {
	gen_BsMuonM_pt_.push_back(bsDaughter.pt());
	gen_BsMuonM_eta_.push_back(bsDaughter.eta());
	gen_BsMuonM_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonM_++;
      }
      if(bsDaughter.pdgId() ==  13){
	gen_BsMuonP_pt_.push_back(bsDaughter.pt());
	gen_BsMuonP_eta_.push_back(bsDaughter.eta());
	gen_BsMuonP_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonP_++;
      }
      if(bsDaughter.pdgId() ==  22){
	gen_BsPhoton_pt_.push_back(bsDaughter.pt());
	gen_BsPhoton_eta_.push_back(bsDaughter.eta());
	gen_BsPhoton_phi_.push_back(bsDaughter.phi());
	gen_nBsPhoton_++;
	
      }
    }  //number of daughters

    gen_Bs_pt_.push_back(aBsMeson.pt());
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
      if(i==j) continue;
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
      mumuPt_.push_back(mu_mu_pt);
      mumuEta_.push_back(mydimu.Eta());
      mumuRapidity_.push_back(mydimu.Rapidity());
      mumuPhi_.push_back(mydimu.Phi());
      mumuMass_.push_back(mu_mu_mass);
      mumuMassE_.push_back(mu_mu_mass_err);
      
      mumuPx_.push_back(mumu_KP->currentState().globalMomentum().x());
      mumuPy_.push_back(mumu_KP->currentState().globalMomentum().y());
      mumuPz_.push_back(mumu_KP->currentState().globalMomentum().z()); 
      mumuDR_.push_back(deltaR(mu,mu2));
      
      mumuVtxCL_.push_back(mu_mu_vtx_cl);
      mumuVtxChi2_.push_back(mumu_KV->chiSquared());
      mumuVtxNdof_.push_back(mumu_KV->degreesOfFreedom());
      mumuVtxX_.push_back(mumu_KV->position().x());
      mumuVtxY_.push_back(mumu_KV->position().y());
      mumuVtxZ_.push_back(mumu_KV->position().z());
      
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

  // loop over photons
  for (auto pho = gedPhotonsHandle->begin(); pho != gedPhotonsHandle->end(); ++pho) {
    phoE_             .push_back(pho->energy());
    phoEt_            .push_back(pho->et());
    phoEta_           .push_back(pho->eta());
    phoPhi_           .push_back(pho->phi());
    
    nPho_++;
  } // photons loop
}

void BsToMuMuGammaNTuplizer::fillPFPhotons(const edm::Event& e, const edm::EventSetup& es)
{
  // Fills tree branches with photons.
  edm::Handle<std::vector<reco::PFCandidate> > pfPhotonsHandle;
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


void BsToMuMuGammaNTuplizer::fillSC(edm::Event const& e, const edm::EventSetup& es) {
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
	scRawE_.push_back(sc.rawEnergy());
	scRawEt_.push_back(sc.rawEnergy()/cosh(sc.eta()));
	scEta_.push_back(sc.eta());
	scPhi_.push_back(sc.phi());
	scEtaWidth_.push_back(sc.etaWidth());
	scPhiWidth_.push_back(sc.phiWidth());
	
	reco::CaloCluster caloBC(*sc.seed());  
	showerShapes_.clear();
	if(abs(sc.eta()) <= 1.442)showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
	if(abs(sc.eta()) >= 1.566)showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
	superCluster_e5x5_.push_back(reduceFloat(showerShapes_[0],nBits_));
	superCluster_e2x2Ratio_.push_back(reduceFloat(showerShapes_[1],nBits_));
	superCluster_e3x3Ratio_.push_back(reduceFloat(showerShapes_[2],nBits_));
	superCluster_eMaxRatio_.push_back(reduceFloat(showerShapes_[3],nBits_));
	superCluster_e2ndRatio_.push_back(reduceFloat(showerShapes_[4],nBits_));
	superCluster_eTopRatio_.push_back(reduceFloat(showerShapes_[5],nBits_));
	superCluster_eRightRatio_.push_back(reduceFloat(showerShapes_[6],nBits_));
	superCluster_eBottomRatio_.push_back(reduceFloat(showerShapes_[7],nBits_));
	superCluster_eLeftRatio_.push_back(reduceFloat(showerShapes_[8],nBits_));
	superCluster_e2x5MaxRatio_.push_back(reduceFloat(showerShapes_[9],nBits_));
	superCluster_e2x5TopRatio_.push_back(reduceFloat(showerShapes_[10],nBits_));
	superCluster_e2x5RightRatio_.push_back(reduceFloat(showerShapes_[11],nBits_));
	superCluster_e2x5BottomRatio_.push_back(reduceFloat(showerShapes_[12],nBits_));
	superCluster_e2x5LeftRatio_.push_back(reduceFloat(showerShapes_[13],nBits_));
	superCluster_swissCross_.push_back(reduceFloat(showerShapes_[14],nBits_));
	superCluster_r9_.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/sc.rawEnergy(),nBits_));
	superCluster_sigmaIetaIeta_.push_back(reduceFloat(showerShapes_[16],nBits_)); 
	superCluster_sigmaIetaIphi_.push_back(reduceFloat(showerShapes_[17],nBits_)); 
	superCluster_sigmaIphiIphi_.push_back(reduceFloat(showerShapes_[18],nBits_)); 
	superCluster_full5x5_e5x5_.push_back(reduceFloat(showerShapes_[19],nBits_));
	superCluster_full5x5_e2x2Ratio_.push_back(reduceFloat(showerShapes_[20],nBits_));
	superCluster_full5x5_e3x3Ratio_.push_back(reduceFloat(showerShapes_[21],nBits_));
	superCluster_full5x5_eMaxRatio_.push_back(reduceFloat(showerShapes_[22],nBits_));
	superCluster_full5x5_e2ndRatio_.push_back(reduceFloat(showerShapes_[23],nBits_));
	superCluster_full5x5_eTopRatio_.push_back(reduceFloat(showerShapes_[24],nBits_));
	superCluster_full5x5_eRightRatio_.push_back(reduceFloat(showerShapes_[25],nBits_));
	superCluster_full5x5_eBottomRatio_.push_back(reduceFloat(showerShapes_[26],nBits_));
	superCluster_full5x5_eLeftRatio_.push_back(reduceFloat(showerShapes_[27],nBits_));
	superCluster_full5x5_e2x5MaxRatio_.push_back(reduceFloat(showerShapes_[28],nBits_));
	superCluster_full5x5_e2x5TopRatio_.push_back(reduceFloat(showerShapes_[29],nBits_));
	superCluster_full5x5_e2x5RightRatio_.push_back(reduceFloat(showerShapes_[30],nBits_));
	superCluster_full5x5_e2x5BottomRatio_.push_back(reduceFloat(showerShapes_[31],nBits_));
	superCluster_full5x5_e2x5LeftRatio_.push_back(reduceFloat(showerShapes_[32],nBits_));
	superCluster_full5x5_swissCross_.push_back(reduceFloat(showerShapes_[33],nBits_));
	superCluster_full5x5_r9_.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/sc.rawEnergy(),nBits_));
	superCluster_full5x5_sigmaIetaIeta_.push_back(reduceFloat(showerShapes_[35],nBits_)); 
	superCluster_full5x5_sigmaIetaIphi_.push_back(reduceFloat(showerShapes_[36],nBits_)); 
	superCluster_full5x5_sigmaIphiIphi_.push_back(reduceFloat(showerShapes_[37],nBits_)); 
	
	++nSC_;
	double dRmin=1e9;
	bool foundGsfEleMatch=false;
	for(auto const& ele : *gsfElectronHandle)
	  {
	    auto dr=deltaR(*(ele.superCluster()),sc);
	    dRmin = dr<dRmin ? dr : dRmin;
	    if( &( *(ele.superCluster()) ) == &sc) 
	      {
		foundGsfEleMatch=true;
		break;
	      }
	  }
	if(gsfElectronHandle->size()<1)
	  {
	    dRmin=-0.333;
	  }
	
	scMinDrWithGsfElectornSC_.push_back(dRmin);
	scFoundGsfMatch_.push_back(foundGsfEleMatch);
	
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

float BsToMuMuGammaNTuplizer::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


//define this as a plug-in
DEFINE_FWK_MODULE(BsToMuMuGammaNTuplizer);
