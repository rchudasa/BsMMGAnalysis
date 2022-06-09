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
  isMC(iConfig.getParameter<Bool_t>("isMC")),
  isRECO(iConfig.getParameter<Bool_t>("isRECO")),
  isMiniAOD(iConfig.getParameter<Bool_t>("isMiniAOD")),
  Run2_2018_(iConfig.getParameter<Bool_t>("Run2_2018")),
  doGenParticles_(iConfig.getParameter<Bool_t>("doGenParticles")),
  doMuons_(iConfig.getParameter<Bool_t>("doMuons")),
  doDimuons_(iConfig.getParameter<Bool_t>("doDimuons")),
  doMuMuK_(iConfig.getParameter<Bool_t>("doMuMuK")),
  doJPsiK_(iConfig.getParameter<Bool_t>("doJPsiK")),
  doPsi2SK_(iConfig.getParameter<Bool_t>("doPsi2SK")),
  doMuMuKK_(iConfig.getParameter<Bool_t>("doMuMuKK")),
  doMuMuGamma(iConfig.getParameter<Bool_t>("doMuMuGamma")),
  doJPsiGamma(iConfig.getParameter<Bool_t>("doJPsiGamma")),
  doPhotons_(iConfig.getParameter<Bool_t>("doPhotons")),
  doPFPhotons_(iConfig.getParameter<Bool_t>("doPFPhotons")),
  doSuperClusters_(iConfig.getParameter<Bool_t>("doSuperClusters")),
  doParticleFlow(iConfig.getParameter<Bool_t>("doParticleFlow")),
  doGeneralTracks(iConfig.getParameter<Bool_t>("doGeneralTracks")),
  doECALClusters(iConfig.getParameter<Bool_t>("doECALClusters")),
  doHCALClusters(iConfig.getParameter<Bool_t>("doHCALClusters")),
  doPrimaryVetrices(iConfig.getParameter<Bool_t>("doPrimaryVetrices")),
  doBeamSpot(iConfig.getParameter<Bool_t>("doBeamSpot")),
  doFlatPt_(iConfig.getParameter<Bool_t>("doFlatPt")),
  doHLT(iConfig.getParameter<Bool_t>("doHLT")),
  energyMatrixSize_(2)
{  
  if(doJPsiGamma)
  {
    doMuMuGamma=true;
  }

  if(doMuMuGamma)
  {
    doMuons_=true;
    doPhotons_=true;
    doSuperClusters_=true;
    doDimuons_=true;
  }
  if(doPsi2SK_)
  {
     doMuMuK_=true;
  }
  if(doJPsiK_)
  {
     doMuMuK_=true;
  }
  if(doMuMuKK_)
  {
      doMuMuK_=true;
  }

  if(doMuMuK_)
  {
      doDimuons_=true;   
  }
  
  energyMatrixSizeFull_=(2*energyMatrixSize_+1)*(2*energyMatrixSize_+1);
  Utility= new Utils();
  
 // if(doMuons_ or do doDimuons ) 
  track_builder_token_=esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  muonToken_     = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
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

  if(doHCALClusters and isRECO)
  if(doECALClusters and isRECO)
  {
     ecalClusterCollection_        = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("ecalClusterSrc"));
  }
  if(doGeneralTracks and isRECO)
  {
     generalTracksCollection_        = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTrackSrc"));
  }


  if(doHLT) {
    triggerFilter    = iConfig.getParameter<std::vector<std::string>>("TriggerFilters");
    trigTable    = iConfig.getParameter<std::vector<std::string>>("TriggerNames");
    triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTResult"));
    triggerEvent_token_ = consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerEvent"));
  }


  
  beamSpotToken_  	   =consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  primaryVtxToken_       =consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  pfCandidateCollection_        = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("particlFlowSrc"));
  
  if(isMC)   genParticlesCollection_          =consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  
  diMuonCharge_    =  iConfig.getUntrackedParameter<bool>("diMuonCharge")    ;
  
  etaMax_muon               =  iConfig.getUntrackedParameter<double>("muon_EtaMax")        ;
  dcaMax_muon_bs            =  iConfig.getUntrackedParameter<double>("muon_dcaMAX")        ;
  pTMinMuons 		        =  iConfig.getUntrackedParameter<double>("muon_minPt");
  pTMinPFPhotons            =  iConfig.getUntrackedParameter<double>("PFPhoton_minPt");
  trackIP_zMax_muon	        =  iConfig.getUntrackedParameter<double>("muon_zIPMax")        ;
  trackIP_rMax_muon	        =  iConfig.getUntrackedParameter<double>("muon_rIPMax")        ;
  
  maxTwoTrackDOCA_          =  iConfig.getUntrackedParameter<double>("maxTwoTrackDOCA")      ;
  minDimuon_pt              =  iConfig.getUntrackedParameter<double>("dimuon_minPt")      ;
  cl_dimuon_vtx             =  iConfig.getUntrackedParameter<double>("dimuon_minVtxCL")      ;
  ls_max_dimuonBS           =  iConfig.getUntrackedParameter<double>("dimuon_maxLStoBS")       ;
  dcaMax_dimuon_mumu        =  iConfig.getUntrackedParameter<double>("dimuon_maxDCAMuMu")        ;
  cosAlphaMax_dimuonBs      =  iConfig.getUntrackedParameter<double>("dimuon_maxCosAlphaToBS") ;
  minDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_minInvMass")    ;
  maxDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_maxInvMass")    ;
  
  minJPsiMass_    =  iConfig.getParameter<double>("minJPsiMass")    ;
  maxJPsiMass_    =  iConfig.getParameter<double>("maxJPsiMass")    ;
  minJPsiGammaMass_    =  iConfig.getParameter<double>("minJPsiGammaMass")    ;
  maxJPsiGammaMass_    =  iConfig.getParameter<double>("maxJPsiGammaMass")    ;
  minBsMuMuGammaMass_  =  iConfig.getParameter<double>("minBsToMuMuGammaMass")    ;
  maxBsMuMuGammaMass_    =  iConfig.getParameter<double>("maxBsToMuMuGammaMass")    ;

  ptMinKaon_    = iConfig.getParameter<double>("ptMinKaon") ;
  etaMaxKaon_   = iConfig.getParameter<double>("etaMaxKaon") ;
  minBKmmMass_  = iConfig.getParameter<double>("minBKmmMass") ;
  maxBKmmMass_  = iConfig.getParameter<double>("maxBKmmMass") ;

  printMsg=iConfig.getParameter<Bool_t>("verbose");
  doMuMuGamma=iConfig.getParameter<Bool_t>("doMuMuGamma");
  nBits_                         = iConfig.getParameter<int>("nBits"); 
  doCompression_                 = iConfig.getParameter<Bool_t>("doCompression");  
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
    addHLTObjectBranches();

    //theTree->Branch("trigResult",    &trigResult);
    //theTree->Branch("trigPrescales", &trigPrescales);
    //theTree->Branch("l1Table",       &l1Table);
    //theTree->Branch("l1Prescales",   &l1Prescales);
  
    SetupTriggerStorageVectors();
    SetupTriggerBranches();
  }


  if (doGenParticles_) {
    addGenBranches();
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
    //theTree->Branch("phoSCRawE",               &phoSCRawE_);
    //theTree->Branch("phoESEnP1",               &phoESEnP1_);
    //theTree->Branch("phoESEnP2",               &phoESEnP2_);
    //theTree->Branch("phoSCEta",                &phoSCEta_);
    //theTree->Branch("phoSCPhi",                &phoSCPhi_);
    //theTree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
    //theTree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
    //theTree->Branch("phoSCBrem",               &phoSCBrem_);
    //theTree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
    //theTree->Branch("phoEleVeto",              &phoEleVeto_);
    //theTree->Branch("phoR9",                   &phoR9_);
    //theTree->Branch("phoHoverE",               &phoHoverE_);
    //theTree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
    //theTree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
    //theTree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
    //theTree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
    //theTree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
    //theTree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
    //theTree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
    
    theTree->Branch("phoPFChIso",              &phoPFChIso_);
    theTree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
    theTree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
    theTree->Branch("phoEcalPFClusterIso",     &phoEcalPFClusterIso_);
    theTree->Branch("phoHcalPFClusterIso",     &phoHcalPFClusterIso_);
    theTree->Branch("phoIDMVA",                &phoIDMVA_);
   
    //theTree->Branch("phoSeedTime",             &phoSeedTime_);
    //theTree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
    //theTree->Branch("phoMIPTotEnergy",         &phoMIPTotEnergy_);
    //theTree->Branch("phoMIPChi2",                      &phoMIPChi2_);
    //theTree->Branch("phoMIPSlope",                     &phoMIPSlope_);
    //theTree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
    //theTree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
    //theTree->Branch("phoMIPIsHalo",                    &phoMIPIsHalo_);

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
  if(doMuMuK_)
  {
        addMuMuKBranches();
  }
  if(doJPsiGamma)
  {
    addJPsiGammaBranches();
    if(doPFPhotons_)
    {
        addPF_JPsiGammaBranches();
    }
  }
  if(doMuMuGamma)
  {
    addMMGBranches();
    if(doPFPhotons_)
    {
        addPF_MMGBranches();
    }
  }
  if(doMuons_)
  {
        addMuonBranches();
  }
  if(doDimuons_)
  {
        addDimuonBranches();
  }
  if(doParticleFlow)
  { 
    std::cout<<" Particle Flow branches are added !! \n";
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

  run_    = iEvent.id().run();
  event_  = iEvent.id().event();
  lumis_  = iEvent.luminosityBlock();
  isData_ = iEvent.isRealData();
  
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
  theTTBuilder_ = iSetup.getHandle(track_builder_token_);
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);
  iEvent.getByToken(primaryVtxToken_, pvHandle_);
  iEvent.getByToken(pfCandidateCollection_, pfCandidateHandle);
  beamSpot = beamSpotHandle.product();
  if(isMC) 
       iEvent.getByToken(genParticlesCollection_, genParticleCollection);

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
  //  Get BeamSpot
  if(doBeamSpot)
  {
  
  //reco::BeamSpot beamSpot = *beamSpotHandle;

  // adding  BEAMSOPT 
  beamspot_x_			= beamSpot->x0();  ;
  beamspot_y_			= beamSpot->y0();  ;
  beamspot_z_			= beamSpot->z0();  ;
  beamspot_x_error_		= beamSpot->x0Error();  ;
  beamspot_y_error_		= beamSpot->y0Error();  ;
  beamspot_z_error_		= beamSpot->z0Error();  ;
  beamspot_covXX        = beamSpot->covariance()(0,0);
  beamspot_covXY        = beamSpot->covariance()(0,1);
  beamspot_covXZ        = beamSpot->covariance()(0,2);
  beamspot_covYY        = beamSpot->covariance()(1,1);
  beamspot_covYZ        = beamSpot->covariance()(1,2);
  beamspot_covZZ        = beamSpot->covariance()(2,2);

  beamspot_dxdz_   		= beamSpot->dxdz();  ;
  beamspot_dydz_       	= beamSpot->dydz();  ;
  beamspot_sigmaZ_		= beamSpot->sigmaZ();  ;
  beamspot_dxdz_error_	= beamSpot->dxdzError();  ;
  beamspot_dydz_error_	= beamSpot->dydzError();  ;
  beamspot_sigmaZError_	= beamSpot->sigmaZ0Error();  ;
  beamspot_beamWidthX_	= beamSpot->BeamWidthX();  ;
  beamspot_beamWidthY_	= beamSpot->BeamWidthY();  ;
  beamspot_beamWidthX_error_	= beamSpot->BeamWidthXError();  ;
  beamspot_beamWidthY_error_	= beamSpot->BeamWidthXError();  ;
  }

  if(doPrimaryVetrices)
  {
        fillPrimaryVertexBranches(iEvent,iSetup);
  }
  // MC truth
  if (isMC and doGenParticles_) {
    fillGenBranches(iEvent);
  }

  if (doHLT) 		    fillHLT(iEvent);
  if (doDimuons_)    	fillDimuonBranches(iEvent, iSetup);
  if (doPhotons_)    	fillPhotons(iEvent, iSetup);
  if (doPFPhotons_) 	fillPFPhotons(iEvent, iSetup);
  if (doSuperClusters_) {
        reco::Vertex pv(math::XYZPoint(0, 0, -999), math::Error<3>::type()); 
        fillSC(iEvent, iSetup,pv);
  }
  if (doMuons_)     	fillMuonBranches(iEvent, iSetup);
  if (doGeneralTracks)  fillGeneralTrackCollectionBranches(iEvent,iSetup);
  if (doParticleFlow)   fillPFCandiateCollection(iEvent,iSetup);
  
  if (isRECO and doHCALClusters) {
        fillHCALClusterCollection(iEvent,iSetup);
  }
  if (isRECO and doECALClusters)   fillECALClusterCollection(iEvent,iSetup);
  theTree->Fill();
  
}

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
    phoR9_              .push_back(pho->r9());
    phoHoverE_          .push_back(pho->hadTowOverEm());
    phoESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*pho).superCluster())));
    
    phoPFChIso_         .push_back(pho->chargedHadronIso()); 
    phoPFPhoIso_        .push_back(pho->photonIso());
    phoPFNeuIso_        .push_back(pho->neutralHadronIso());
    phoEcalPFClusterIso_.push_back(pho->ecalPFClusterIso());
    phoHcalPFClusterIso_.push_back(pho->hcalPFClusterIso());

    // phoEleVeto_         .push_back((Int_t)pho->passElectronVeto());
    // phoIDMVA_           .push_back(pho->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));  
    
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
    Bool_t isBarrel = seed.subdetId() == EcalBarrel;
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
 // edm::Handle<edm::View<reco::PFCandidate> > pfCandidateHandle;
 // e.getByToken(pfPhotonsCollection_, pfCandidateHandle);
  
  // loop over photons
  for (auto pf = pfCandidateHandle->begin(); pf != pfCandidateHandle->end(); ++pf) {
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
    for (auto const& sc : scs) {
	
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
      Bool_t foundGsfEleMatch=false;
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
        scNHcalRecHitInDIEta5IPhi5.push_back(Float_t(nhcalRechit_));
        scEFromHcalRecHitInDIEta2IPhi2.push_back(hcal2Energy);
        scNHcalRecHitInDIEta2IPhi2.push_back(Float_t(nhcal2Rechit_));
        

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
  TrigPrescales_store = new Int_t [numTrigs];
  TrigResult_store    = new Bool_t[numTrigs];
}

void BsToMuMuGammaNTuplizer::ClearTrggerStorages()
{
  for(uint32_t i=0;i<trigTable.size();i++)
    {
      TrigResult_store[i]=false;
      TrigPrescales_store[i]=1.0;
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
	  TrigResult_store[i]=true;	
	  TrigPrescales_store[i]=trigPrescales[foundTrig];	
	}
      else
	{
	  TrigResult_store[i]=false;	
	  TrigPrescales_store[i]=1.0;	
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


  // Get the trigger results
  bool validTriggerEvent = true;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  const trigger::TriggerEvent dummyTE;
  iEvent.getByToken(triggerEvent_token_, triggerEventHandle);
  if (!triggerEventHandle.isValid()) {
    std::cout<< "Error! Can't get the product: triggerEvent_" << endl;
    validTriggerEvent = false;
  }
  
  const trigger::TriggerEvent& triggerEvent(validTriggerEvent ? *(triggerEventHandle.product()) : dummyTE);

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

  vector<int> trigIdx;
  vector<int> Keys;
  for (uint filterIndex = 0; filterIndex < triggerEvent.sizeFilters(); ++filterIndex) {  //loop over all trigger filters in event (i.e. filters passed)
    string label = triggerEvent.filterTag(filterIndex).label();
    for(unsigned int  tIdx=0;tIdx < triggerFilter.size() ; tIdx++)
    {
      if (label.find(triggerFilter[tIdx]) != string::npos) {  
      for (uint filterKeyIndex = 0; filterKeyIndex < triggerEvent.filterKeys(filterIndex).size();
           ++filterKeyIndex) {  //loop over keys to objects passing this filter
        trigIdx.push_back(tIdx);  
        Keys.push_back(triggerEvent.filterKeys(filterIndex)[filterKeyIndex]);  //add keys to a vector for later reference
        }
      }
    }
   }
   fillHLTL3MuonBranches(  trigIdx, Keys, triggerEvent );

}

std::vector<Float_t> BsToMuMuGammaNTuplizer::getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology)
{
  std::vector<Float_t> shapes;
  shapes.resize(38); 
  locCov_.clear();
  full5x5_locCov_.clear();
  locCov_ = EcalClusterTools::localCovariances(*caloBC, recHits, topology);
  full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    
  Float_t e5x5 = EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
  Float_t e3x3 = EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
  Float_t eMax = EcalClusterTools::eMax(*caloBC, recHits); // eMax
  Float_t eTop = EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
  Float_t eRight = EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
  Float_t eBottom = EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
  Float_t eLeft = EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
  Float_t e4 = eTop + eRight + eBottom + eLeft;

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
  Float_t full5x5_e5x5 = noZS::EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
  Float_t full5x5_e3x3 = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
  Float_t full5x5_eMax = noZS::EcalClusterTools::eMax(*caloBC, recHits); // eMax
  Float_t full5x5_eTop = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
  Float_t full5x5_eRight = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
  Float_t full5x5_eBottom = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
  Float_t full5x5_eLeft = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
  Float_t full5x5_e4 = full5x5_eTop + full5x5_eRight + full5x5_eBottom + full5x5_eLeft;

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

    int i=0;
    storageMapInt["nPrimaryVertex"]=0;
    for(auto&  aVertex : *pvHandle_){
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
    i++;
  } // loop over primary vertex collection
    storageMapInt["nPrimaryVertex"]=i;
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
  Int_t i=0;
  storageMapInt["nPFCandidates"]=0;
  for (auto pfCd = pfCandidateHandle->begin();pfCd != pfCandidateHandle->end(); pfCd++) 
  {
    if(  pfCd->pt() < 2.0 ) continue;

      storageMapFloatArray["pf_ecalEnergy"]     [i]= pfCd->ecalEnergy();   
      storageMapFloatArray["pf_ecalRawEnergy"]  [i]= pfCd->rawEcalEnergy(); 
      storageMapFloatArray["pf_hcalEnergy"]     [i]= pfCd->hcalEnergy(); 
      storageMapFloatArray["pf_hcalRawEnergy"]  [i]= pfCd->rawHcalEnergy(); 
      storageMapFloatArray["pf_HoE"]            [i]= pfCd->hcalEnergy()/(1e-6 +  pfCd->ecalEnergy() ); 
      storageMapFloatArray["pf_mvaIso"]         [i]= pfCd->mva_Isolated(); 
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


void BsToMuMuGammaNTuplizer::addGenBranches()
{
      storageMapInt["nHCALClusters"]=0;
      theTree->Branch("nHCALClusters",   &storageMapInt["nHCALClusters"]);

      storageMapInt["nGen"]   = 0;
      theTree->Branch("nGen",   &storageMapInt["nGen"]);

      storageMapFloatArray["gen_mcPID"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcPID", storageMapFloatArray["gen_mcPID"],"gen_mcPID[nGen/F");
      storageMapFloatArray["gen_mcStatus"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcStatus", storageMapFloatArray["gen_mcStatus"],"gen_mcStatus[nGen/F");
      storageMapFloatArray["gen_mcVtx_x"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcVtx_x", storageMapFloatArray["gen_mcVtx_x"],"gen_mcVtx_x[nGen/F");
      storageMapFloatArray["gen_mcVtx_y"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcVtx_y", storageMapFloatArray["gen_mcVtx_y"],"gen_mcVtx_y[nGen/F");
      storageMapFloatArray["gen_mcVtx_z"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcVtx_z", storageMapFloatArray["gen_mcVtx_z"],"gen_mcVtx_z[nGen/F");
      storageMapFloatArray["gen_mcPt"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcPt", storageMapFloatArray["gen_mcPt"],"gen_mcPt[nGen/F");
      storageMapFloatArray["gen_mcEta"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcEta", storageMapFloatArray["gen_mcEta"],"gen_mcEta[nGen/F");
      storageMapFloatArray["gen_mcPhi"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcPhi", storageMapFloatArray["gen_mcPhi"],"gen_mcPhi[nGen/F");
      storageMapFloatArray["gen_mcE"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcE", storageMapFloatArray["gen_mcE"],"gen_mcE[nGen/F");
      storageMapFloatArray["gen_mcEt"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcEt", storageMapFloatArray["gen_mcEt"],"gen_mcEt[nGen/F");
      storageMapFloatArray["gen_mcMass"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcMass", storageMapFloatArray["gen_mcMass"],"gen_mcMass[nGen/F");
      storageMapFloatArray["gen_mcMotherPDGID"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcMotherPDGID", storageMapFloatArray["gen_mcMotherPDGID"],"gen_mcMotherPDGID[nGen/F");
      storageMapFloatArray["gen_mcMotherPt"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcMotherPt", storageMapFloatArray["gen_mcMotherPt"],"gen_mcMotherPt[nGen/F");
      storageMapFloatArray["gen_mcMotherEta"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcMotherEta", storageMapFloatArray["gen_mcMotherEta"],"gen_mcMotherEta[nGen/F");
      storageMapFloatArray["gen_mcMotherPhi"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcMotherPhi", storageMapFloatArray["gen_mcMotherPhi"],"gen_mcMotherPhi[nGen/F");
      storageMapFloatArray["gen_mcMotherMass"]   = new Float_t[N_GEN_MAX];
      theTree->Branch("gen_mcMotherMass", storageMapFloatArray["gen_mcMotherMass"],"gen_mcMotherMass[nGen/F");
}

void BsToMuMuGammaNTuplizer::fillGenBranches( const edm::Event& iEvent)
{
    Int_t nGen=0;
    for (auto p = genParticleCollection->begin(); p != genParticleCollection->end(); ++p) {
        storageMapFloatArray["gen_mcPID"][nGen]   = p->pdgId() ;
        storageMapFloatArray["gen_mcStatus"][nGen]   = p->status(); 
        storageMapFloatArray["gen_mcVtx_x"][nGen]   = p->vx() ;
        storageMapFloatArray["gen_mcVtx_y"][nGen]   = p->vy() ;
        storageMapFloatArray["gen_mcVtx_z"][nGen]   = p->vz() ;
        storageMapFloatArray["gen_mcPt"][nGen]   = p->pt() ;
        storageMapFloatArray["gen_mcEta"][nGen]   = p->eta(); 
        storageMapFloatArray["gen_mcPhi"][nGen]   = p->phi() ;
        storageMapFloatArray["gen_mcE"][nGen]   = p->energy() ;
        storageMapFloatArray["gen_mcEt"][nGen]   = p->et() ;
        storageMapFloatArray["gen_mcMass"][nGen]   = p->mass() ;
        storageMapFloatArray["gen_gen_mcMotherPDGID"][nGen]   = p->mother()->pdgId() ;
        storageMapFloatArray["gen_gen_mcMotherPt"][nGen]      = p->mother()->pt() ;
        storageMapFloatArray["gen_gen_mcMotherEta"][nGen]     = p->mother()->eta() ;
        storageMapFloatArray["gen_gen_mcMotherPhi"][nGen]     = p->mother()->phi() ;
        storageMapFloatArray["gen_gen_mcMotherMass"][nGen]    = p->mother()->mass() ;
        nGen++;
     }
        storageMapInt["nGen"]     = nGen ;
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

void BsToMuMuGammaNTuplizer::addMuonBranches()
{
      storageMapInt["nMuons"]=0;
      theTree->Branch("nMuons",   &storageMapInt["nMuons"]);
      storageMapFloatArray["muon_Charge"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Charge", storageMapFloatArray["muon_Charge"],"muon_Charge[nMuons]/F");
      storageMapFloatArray["muon_Idx"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Idx", storageMapFloatArray["muon_Idx"],"muon_Idx[nMuons]/F");
      storageMapFloatArray["muon_Pt"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Pt", storageMapFloatArray["muon_Pt"],"muon_Pt[nMuons]/F");
      storageMapFloatArray["muon_Eta"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Eta", storageMapFloatArray["muon_Eta"],"muon_Eta[nMuons]/F");
      storageMapFloatArray["muon_Phi"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Phi", storageMapFloatArray["muon_Phi"],"muon_Phi[nMuons]/F");
      storageMapFloatArray["muon_CL"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_CL", storageMapFloatArray["muon_CL"],"muon_CL[nMuons]/F");
      storageMapFloatArray["muon_NormChi2"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_NormChi2", storageMapFloatArray["muon_NormChi2"],"muon_NormChi2[nMuons]/F");
      storageMapFloatArray["muon_Vx"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Vx", storageMapFloatArray["muon_Vx"],"muon_Vx[nMuons]/F");
      storageMapFloatArray["muon_Vy"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Vy", storageMapFloatArray["muon_Vy"],"muon_Vy[nMuons]/F");
      storageMapFloatArray["muon_Vz"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_Vz", storageMapFloatArray["muon_Vz"],"muon_Vz[nMuons]/F");
      storageMapFloatArray["muon_dxyBS"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_dxyBS", storageMapFloatArray["muon_dxyBS"],"muon_dxyBS[nMuons]/F");
      storageMapFloatArray["muon_dzBS"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_dzBS", storageMapFloatArray["muon_dzBS"],"muon_dzBS[nMuons]/F");
      storageMapFloatArray["muon_MissingHits"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_MissingHits", storageMapFloatArray["muon_MissingHits"],"muon_MissingHits[nMuons]/F");
      
      storageMapUint64Array["muon_selectors"]   = new uint64_t[N_MUON_MAX];
      theTree->Branch("muon_selectors", storageMapUint64Array["muon_selectors"],"muon_selectors[nMuons]/l");
      storageMapBoolArray["muon_isGlobalMuon"]   = new Bool_t[N_MUON_MAX];
      theTree->Branch("muon_isGlobalMuon", storageMapBoolArray["muon_isGlobalMuon"],"muon_isGlobalMuon[nMuons]/O");
      storageMapBoolArray["muon_isTrackerMuon"]   = new Bool_t[N_MUON_MAX];
      theTree->Branch("muon_isTrackerMuon", storageMapBoolArray["muon_isTrackerMuon"],"muon_isTrackerMuon[nMuons]/O");
      storageMapBoolArray["muon_isStandAloneMuon"]   = new Bool_t[N_MUON_MAX];
      theTree->Branch("muon_isStandAloneMuon", storageMapBoolArray["muon_isStandAloneMuon"],"muon_isStandAloneMuon[nMuons]/O");
      storageMapBoolArray["muon_isCaloMuon"]   = new Bool_t[N_MUON_MAX];
      theTree->Branch("muon_isCaloMuon", storageMapBoolArray["muon_isCaloMuon"],"muon_isCaloMuon[nMuons]/O");
      storageMapBoolArray["muon_isPFMuon"]   = new Bool_t[N_MUON_MAX];
      theTree->Branch("muon_isPFMuon", storageMapBoolArray["muon_isPFMuon"],"muon_isPFMuon[nMuons]/O");

      storageMapFloatArray["muon_chi2LocalPosition"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_chi2LocalPosition", storageMapFloatArray["muon_chi2LocalPosition"],"muon_chi2LocalPosition[nMuons]/F");
      storageMapFloatArray["muon_glbNormChi2"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_glbNormChi2", storageMapFloatArray["muon_glbNormChi2"],"muon_glbNormChi2[nMuons]/F");
      storageMapFloatArray["muon_glbTrackProbability"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_glbTrackProbability", storageMapFloatArray["muon_glbTrackProbability"],"muon_glbTrackProbability[nMuons]/F");
      storageMapFloatArray["muon_nLostHitsInner"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_nLostHitsInner", storageMapFloatArray["muon_nLostHitsInner"],"muon_nLostHitsInner[nMuons]/F");
      storageMapFloatArray["muon_nLostHitsOn"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_nLostHitsOn", storageMapFloatArray["muon_nLostHitsOn"],"muon_nLostHitsOn[nMuons]/F");
      storageMapFloatArray["muon_nLostHitsOuter"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_nLostHitsOuter", storageMapFloatArray["muon_nLostHitsOuter"],"muon_nLostHitsOuter[nMuons]/F");
      storageMapFloatArray["muon_nPixels"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_nPixels", storageMapFloatArray["muon_nPixels"],"muon_nPixels[nMuons]/F");
      storageMapFloatArray["muon_nValidHits"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_nValidHits", storageMapFloatArray["muon_nValidHits"],"muon_nValidHits[nMuons]/F");
      storageMapFloatArray["muon_trkKink"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_trkKink", storageMapFloatArray["muon_trkKink"],"muon_trkKink[nMuons]/F");
      storageMapFloatArray["muon_trkValidFrac"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_trkValidFrac", storageMapFloatArray["muon_trkValidFrac"],"muon_trkValidFrac[nMuons]/F");
      storageMapFloatArray["muon_trkLayers"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_trkLayers", storageMapFloatArray["muon_trkLayers"],"muon_trkLayers[nMuons]/F");
      storageMapFloatArray["muon_trkLostLayersInner"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_trkLostLayersInner", storageMapFloatArray["muon_trkLostLayersInner"],"muon_trkLostLayersInner[nMuons]/F");
      storageMapFloatArray["muon_trkLostLayersOn"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_trkLostLayersOn", storageMapFloatArray["muon_trkLostLayersOn"],"muon_trkLostLayersOn[nMuons]/F");
      storageMapFloatArray["muon_trkLostLayersOuter"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_trkLostLayersOuter", storageMapFloatArray["muon_trkLostLayersOuter"],"muon_trkLostLayersOuter[nMuons]F");
      storageMapFloatArray["muon_highPurity"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_highPurity", storageMapFloatArray["muon_highPurity"],"muon_highPurity[nMuons]/F");
      storageMapFloatArray["muon_match1_dX"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match1_dX", storageMapFloatArray["muon_match1_dX"],"muon_match1_dX[nMuons]/F");
      storageMapFloatArray["muon_match1_dY"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match1_dY", storageMapFloatArray["muon_match1_dY"],"muon_match1_dY[nMuons]/F");
      storageMapFloatArray["muon_match1_pullDxDz"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match1_pullDxDz", storageMapFloatArray["muon_match1_pullDxDz"],"muon_match1_pullDxDz[nMuons]/F");
      storageMapFloatArray["muon_match1_pullDyDz"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match1_pullDyDz", storageMapFloatArray["muon_match1_pullDyDz"],"muon_match1_pullDyDz[nMuons]/F");
      storageMapFloatArray["muon_match1_pullX"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match1_pullX", storageMapFloatArray["muon_match1_pullX"],"muon_match1_pullX[nMuons]/F");
      storageMapFloatArray["muon_match1_pullY"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match1_pullY", storageMapFloatArray["muon_match1_pullY"],"muon_match1_pullY[nMuons]/F");
      storageMapFloatArray["muon_match2_dX"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match2_dX", storageMapFloatArray["muon_match2_dX"],"muon_match2_dX[nMuons]/F");
      storageMapFloatArray["muon_match2_pullDxDz"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match2_pullDxDz", storageMapFloatArray["muon_match2_pullDxDz"],"muon_match2_pullDxDz[nMuons]/F");
      storageMapFloatArray["muon_match2_pullDyDz"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match2_pullDyDz", storageMapFloatArray["muon_match2_pullDyDz"],"muon_match2_pullDyDz[nMuons]/F");
      storageMapFloatArray["muon_match2_pullX"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match2_pullX", storageMapFloatArray["muon_match2_pullX"],"muon_match2_pullX[nMuons]/F");
      storageMapFloatArray["muon_match2_pullY"]   = new Float_t[N_MUON_MAX];
      theTree->Branch("muon_match2_pullY", storageMapFloatArray["muon_match2_pullY"],"muon_match2_pullY[nMuons]/F");

}

void BsToMuMuGammaNTuplizer::fillMuonBranches( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);
   
  //std::cout << " muon size:" << muons->size() << std::endl; 
  //start loop on first muon 
  
  int nMu=0;
  for (uint32_t i=0;i<muons->size();i++){
   
    auto &mu = muons->at(i);
    auto muTrack = mu.innerTrack();

    if ((muTrack.isNull() == true)) continue;
    
    //const reco::TransientTrack muTrackTT( muTrack, &(*bFieldHandle));
    // # Compute mu- DCA to BeamSpot #
    // theDCAXBS = muTrackTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
    // double DCAmuBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
    // double DCAmuBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
    storageMapFloatArray["muon_Idx"][nMu]          = i;  
    storageMapFloatArray["muon_Charge"][nMu]       = mu.charge();  
    storageMapFloatArray["muon_Pt"][nMu]           = mu.pt() ;
    storageMapFloatArray["muon_Eta"][nMu]          = mu.eta() ;
    storageMapFloatArray["muon_Phi"][nMu]          = mu.phi() ;
    storageMapFloatArray["muon_CL"][nMu]           = TMath::Prob(muTrack->chi2(), static_cast<int>(rint(muTrack->ndof()))) ;
    storageMapFloatArray["muon_NormChi2"][nMu]     = muTrack->normalizedChi2() ;
    storageMapFloatArray["muon_Vx"][nMu]           = mu.vx() ;
    storageMapFloatArray["muon_Vy"][nMu]           = mu.vy() ;
    storageMapFloatArray["muon_Vz"][nMu]           = mu.vz() ;
    storageMapFloatArray["muon_MissingHits"][nMu]  = muTrack->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
											    muTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
											    muTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS) ;
    storageMapBoolArray["muon_isGlobalMuon"][nMu]   = mu.isGlobalMuon() ;
    storageMapBoolArray["muon_isTrackerMuon"][nMu]  = mu.isTrackerMuon() ;
    storageMapBoolArray["muon_isStandAloneMuon"][nMu] = mu.isStandAloneMuon() ;
    storageMapBoolArray["muon_isCaloMuon"][nMu]     = mu.isCaloMuon() ;
    storageMapBoolArray["muon_isPFMuon"][nMu]       = mu.isPFMuon();
    
    storageMapFloatArray["muon_trkKink"][nMu]             = mu.combinedQuality().trkKink ;
    storageMapFloatArray["muon_glbTrackProbability"][nMu] = mu.combinedQuality().glbTrackProbability ;
    storageMapFloatArray["muon_chi2LocalPosition"][nMu]   = mu.combinedQuality().chi2LocalPosition ;

      if (mu.isGlobalMuon())
    	storageMapFloatArray["muon_glbNormChi2"][nMu]=mu.globalTrack()->normalizedChi2();
      else
    	storageMapFloatArray["muon_glbNormChi2"][nMu]=9e4f;
	  
     storageMapFloatArray["muon_nPixels"][nMu]            = muTrack->hitPattern().numberOfValidPixelHits();
	 storageMapFloatArray["muon_nValidHits"][nMu]         = muTrack->hitPattern().numberOfValidTrackerHits();
	 storageMapFloatArray["muon_nLostHitsInner"][nMu]     = muTrack->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS);
	 storageMapFloatArray["muon_nLostHitsOn"][nMu]        = muTrack->hitPattern().numberOfLostTrackerHits(reco::HitPattern::TRACK_HITS);
	 storageMapFloatArray["muon_nLostHitsOuter"][nMu]     = muTrack->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS);
	 storageMapFloatArray["muon_trkLayers"][nMu]          = muTrack->hitPattern().trackerLayersWithMeasurement();
	 storageMapFloatArray["muon_trkLostLayersInner"][nMu] = muTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS);
	 storageMapFloatArray["muon_trkLostLayersOn"][nMu]    = muTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
	 storageMapFloatArray["muon_trkLostLayersOuter"][nMu] = muTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_OUTER_HITS);
	 storageMapFloatArray["muon_highPurity"][nMu]         = muTrack->quality(reco::Track::highPurity);

      if (mu.isTrackerMuon() or mu.isGlobalMuon()){
	     storageMapFloatArray["muon_trkValidFrac"][nMu]  = muTrack->validFraction() ;
      } 
      else 
      {
	     storageMapFloatArray["muon_trkValidFrac"][nMu]  =0.0f ;
      }

      storageMapUint64Array["muon_selectors"][nMu] = mu.selectors(); 
      nMu++; 
      if(nMu>= N_MUON_MAX) 
      {
            std::cout<<" nMu --> "<<nMu<<"  !! \n";
      }

    }
    storageMapInt["nMuons"]=nMu;

}
void BsToMuMuGammaNTuplizer::addDimuonBranches()
{
        std::cout<<"Initilizing Dimuon Branches !\n";
        storageMapInt["nDimuons"]=0;
        theTree->Branch("nDimuons",   &storageMapInt["nDimuons"]);
        
        storageMapIntArray["dimuon_mu1_index"]   = new Int_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu1_index", storageMapIntArray["dimuon_mu1_index"],"dimuon_mu1_index[nDimuons]/I");
        storageMapIntArray["dimuon_mu1_pdgId"]   = new Int_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu1_pdgId", storageMapIntArray["dimuon_mu1_pdgId"],"dimuon_mu1_pdgId[nDimuons]/I");
        storageMapFloatArray["dimuon_mu1_pt"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu1_pt", storageMapFloatArray["dimuon_mu1_pt"],"dimuon_mu1_pt[nDimuons]/F");
        storageMapFloatArray["dimuon_mu1_eta"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu1_eta", storageMapFloatArray["dimuon_mu1_eta"],"dimuon_mu1_eta[nDimuons]/F");
        storageMapFloatArray["dimuon_mu1_phi"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu1_phi", storageMapFloatArray["dimuon_mu1_phi"],"dimuon_mu1_phi[nDimuons]/F");
        storageMapIntArray["dimuon_mu2_index"]   = new Int_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu2_index", storageMapIntArray["dimuon_mu2_index"],"dimuon_mu2_index[nDimuons]/I");
        storageMapIntArray["dimuon_mu2_pdgId"]   = new Int_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu2_pdgId", storageMapIntArray["dimuon_mu2_pdgId"],"dimuon_mu2_pdgId[nDimuons]/I");
        storageMapFloatArray["dimuon_mu2_pt"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu2_pt", storageMapFloatArray["dimuon_mu2_pt"],"dimuon_mu2_pt[nDimuons]/F");
        storageMapFloatArray["dimuon_mu2_eta"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu2_eta", storageMapFloatArray["dimuon_mu2_eta"],"dimuon_mu2_eta[nDimuons]/F");
        storageMapFloatArray["dimuon_mu2_phi"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_mu2_phi", storageMapFloatArray["dimuon_mu2_phi"],"dimuon_mu2_phi[nDimuons]/F");
        storageMapFloatArray["dimuon_doca"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_doca", storageMapFloatArray["dimuon_doca"],"dimuon_doca[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_mu1pt"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_mu1pt", storageMapFloatArray["dimuon_kin_mu1pt"],"dimuon_kin_mu1pt[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_mu1eta"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_mu1eta", storageMapFloatArray["dimuon_kin_mu1eta"],"dimuon_kin_mu1eta[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_mu1phi"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_mu1phi", storageMapFloatArray["dimuon_kin_mu1phi"],"dimuon_kin_mu1phi[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_mu2pt"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_mu2pt", storageMapFloatArray["dimuon_kin_mu2pt"],"dimuon_kin_mu2pt[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_mu2eta"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_mu2eta", storageMapFloatArray["dimuon_kin_mu2eta"],"dimuon_kin_mu2eta[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_mu2phi"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_mu2phi", storageMapFloatArray["dimuon_kin_mu2phi"],"dimuon_kin_mu2phi[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_valid"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_valid", storageMapFloatArray["dimuon_kin_valid"],"dimuon_kin_valid[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_prob"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_prob", storageMapFloatArray["dimuon_kin_vtx_prob"],"dimuon_kin_vtx_prob[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_chi2dof"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_chi2dof", storageMapFloatArray["dimuon_kin_vtx_chi2dof"],"dimuon_kin_vtx_chi2dof[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_mass"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_mass", storageMapFloatArray["dimuon_kin_mass"],"dimuon_kin_mass[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_massErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_massErr", storageMapFloatArray["dimuon_kin_massErr"],"dimuon_kin_massErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_lxy"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_lxy", storageMapFloatArray["dimuon_kin_lxy"],"dimuon_kin_lxy[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_sigLxy"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_sigLxy", storageMapFloatArray["dimuon_kin_sigLxy"],"dimuon_kin_sigLxy[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_alphaBS"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_alphaBS", storageMapFloatArray["dimuon_kin_alphaBS"],"dimuon_kin_alphaBS[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_alphaBSErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_alphaBSErr", storageMapFloatArray["dimuon_kin_alphaBSErr"],"dimuon_kin_alphaBSErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_x"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_x", storageMapFloatArray["dimuon_kin_vtx_x"],"dimuon_kin_vtx_x[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_xErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_xErr", storageMapFloatArray["dimuon_kin_vtx_xErr"],"dimuon_kin_vtx_xErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_y"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_y", storageMapFloatArray["dimuon_kin_vtx_y"],"dimuon_kin_vtx_y[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_yErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_yErr", storageMapFloatArray["dimuon_kin_vtx_yErr"],"dimuon_kin_vtx_yErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_z"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_z", storageMapFloatArray["dimuon_kin_vtx_z"],"dimuon_kin_vtx_z[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_vtx_zErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_vtx_zErr", storageMapFloatArray["dimuon_kin_vtx_zErr"],"dimuon_kin_vtx_zErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pt"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pt", storageMapFloatArray["dimuon_kin_pt"],"dimuon_kin_pt[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_eta"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_eta", storageMapFloatArray["dimuon_kin_eta"],"dimuon_kin_eta[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_phi"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_phi", storageMapFloatArray["dimuon_kin_phi"],"dimuon_kin_phi[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_alpha"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_alpha", storageMapFloatArray["dimuon_kin_alpha"],"dimuon_kin_alpha[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_alphaErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_alphaErr", storageMapFloatArray["dimuon_kin_alphaErr"],"dimuon_kin_alphaErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_l3d"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_l3d", storageMapFloatArray["dimuon_kin_l3d"],"dimuon_kin_l3d[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_sl3d"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_sl3d", storageMapFloatArray["dimuon_kin_sl3d"],"dimuon_kin_sl3d[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pv_z"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pv_z", storageMapFloatArray["dimuon_kin_pv_z"],"dimuon_kin_pv_z[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pv_zErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pv_zErr", storageMapFloatArray["dimuon_kin_pv_zErr"],"dimuon_kin_pv_zErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pvip"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pvip", storageMapFloatArray["dimuon_kin_pvip"],"dimuon_kin_pvip[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_spvip"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_spvip", storageMapFloatArray["dimuon_kin_spvip"],"dimuon_kin_spvip[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pvipErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pvipErr", storageMapFloatArray["dimuon_kin_pvipErr"],"dimuon_kin_pvipErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pv2ip"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pv2ip", storageMapFloatArray["dimuon_kin_pv2ip"],"dimuon_kin_pv2ip[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_spv2ip"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_spv2ip", storageMapFloatArray["dimuon_kin_spv2ip"],"dimuon_kin_spv2ip[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pv2ipErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pv2ipErr", storageMapFloatArray["dimuon_kin_pv2ipErr"],"dimuon_kin_pv2ipErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pvlip"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pvlip", storageMapFloatArray["dimuon_kin_pvlip"],"dimuon_kin_pvlip[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pvlipSig"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pvlipSig", storageMapFloatArray["dimuon_kin_pvlipSig"],"dimuon_kin_pvlipSig[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pvlipErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pvlipErr", storageMapFloatArray["dimuon_kin_pvlipErr"],"dimuon_kin_pvlipErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pv2lip"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pv2lip", storageMapFloatArray["dimuon_kin_pv2lip"],"dimuon_kin_pv2lip[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pv2lipSig"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pv2lipSig", storageMapFloatArray["dimuon_kin_pv2lipSig"],"dimuon_kin_pv2lipSig[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pv2lipErr"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pv2lipErr", storageMapFloatArray["dimuon_kin_pv2lipErr"],"dimuon_kin_pv2lipErr[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_pvIndex"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_pvIndex", storageMapFloatArray["dimuon_kin_pvIndex"],"dimuon_kin_pvIndex[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_tau"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_tau", storageMapFloatArray["dimuon_kin_tau"],"dimuon_kin_tau[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_taue"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_taue", storageMapFloatArray["dimuon_kin_taue"],"dimuon_kin_taue[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_tauxy"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_tauxy", storageMapFloatArray["dimuon_kin_tauxy"],"dimuon_kin_tauxy[nDimuons]/F");
        storageMapFloatArray["dimuon_kin_tauxye"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_kin_tauxye", storageMapFloatArray["dimuon_kin_tauxye"],"dimuon_kin_tauxye[nDimuons]/F");
        storageMapFloatArray["dimuon_nTrks"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_nTrks", storageMapFloatArray["dimuon_nTrks"],"dimuon_nTrks[nDimuons]/F");
        storageMapFloatArray["dimuon_nBMTrks"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_nBMTrks", storageMapFloatArray["dimuon_nBMTrks"],"dimuon_nBMTrks[nDimuons]/F");
        storageMapFloatArray["dimuon_nDisTrks"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_nDisTrks", storageMapFloatArray["dimuon_nDisTrks"],"dimuon_nDisTrks[nDimuons]/F");
        storageMapFloatArray["dimuon_closetrk"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_closetrk", storageMapFloatArray["dimuon_closetrk"],"dimuon_closetrk[nDimuons]/F");
        storageMapFloatArray["dimuon_closetrks1"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_closetrks1", storageMapFloatArray["dimuon_closetrks1"],"dimuon_closetrks1[nDimuons]/F");
        storageMapFloatArray["dimuon_closetrks2"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_closetrks2", storageMapFloatArray["dimuon_closetrks2"],"dimuon_closetrks2[nDimuons]/F");
        storageMapFloatArray["dimuon_closetrks3"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_closetrks3", storageMapFloatArray["dimuon_closetrks3"],"dimuon_closetrks3[nDimuons]/F");
        storageMapFloatArray["dimuon_docatrk"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_docatrk", storageMapFloatArray["dimuon_docatrk"],"dimuon_docatrk[nDimuons]/F");
        storageMapFloatArray["dimuon_m1iso"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_m1iso", storageMapFloatArray["dimuon_m1iso"],"dimuon_m1iso[nDimuons]/F");
        storageMapFloatArray["dimuon_m2iso"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_m2iso", storageMapFloatArray["dimuon_m2iso"],"dimuon_m2iso[nDimuons]/F");
        storageMapFloatArray["dimuon_iso"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_iso", storageMapFloatArray["dimuon_iso"],"dimuon_iso[nDimuons]/F");
        storageMapFloatArray["dimuon_otherVtxMaxProb"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_otherVtxMaxProb", storageMapFloatArray["dimuon_otherVtxMaxProb"],"dimuon_otherVtxMaxProb[nDimuons]/F");
        storageMapFloatArray["dimuon_otherVtxMaxProb1"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_otherVtxMaxProb1", storageMapFloatArray["dimuon_otherVtxMaxProb1"],"dimuon_otherVtxMaxProb1[nDimuons]/F");
        storageMapFloatArray["dimuon_otherVtxMaxProb2"]   = new Float_t[N_DIMU_MAX];
        theTree->Branch("dimuon_otherVtxMaxProb2", storageMapFloatArray["dimuon_otherVtxMaxProb2"],"dimuon_otherVtxMaxProb2[nDimuons]/F");

        if(isMC)
        {
            storageMapIntArray["gen_mm_mu1_pdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mu1_pdgId", storageMapIntArray["gen_mm_mu1_pdgId"],"gen_mm_mu1_pdgId[nDimuons]/I");
            storageMapIntArray["gen_mm_mu1_mpdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mu1_mpdgId", storageMapIntArray["gen_mm_mu1_mpdgId"],"gen_mm_mu1_mpdgId[nDimuons]/I");
            storageMapFloatArray["gen_mm_mu1_pt"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mu1_pt", storageMapFloatArray["gen_mm_mu1_pt"],"gen_mm_mu1_pt[nDimuons]/F");
            storageMapIntArray["gen_mm_mu2_pdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mu2_pdgId", storageMapIntArray["gen_mm_mu2_pdgId"],"gen_mm_mu2_pdgId[nDimuons]/I");
            storageMapIntArray["gen_mm_mu2_mpdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mu2_mpdgId", storageMapIntArray["gen_mm_mu2_mpdgId"],"gen_mm_mu2_mpdgId[nDimuons]/I");
            storageMapFloatArray["gen_mm_mu2_pt"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mu2_pt", storageMapFloatArray["gen_mm_mu2_pt"],"gen_mm_mu2_pt[nDimuons]/F");
            storageMapFloatArray["gen_mm_mass"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mass", storageMapFloatArray["gen_mm_mass"],"gen_mm_mass[nDimuons]/F");
            storageMapFloatArray["gen_mm_pt"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_pt", storageMapFloatArray["gen_mm_pt"],"gen_mm_pt[nDimuons]/F");
            storageMapIntArray["gen_mm_pdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_pdgId", storageMapIntArray["gen_mm_pdgId"],"gen_mm_pdgId[nDimuons]/I");
            storageMapIntArray["gen_mm_mpdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_mpdgId", storageMapIntArray["gen_mm_mpdgId"],"gen_mm_mpdgId[nDimuons]/I");
            storageMapIntArray["gen_mm_cpdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_cpdgId", storageMapIntArray["gen_mm_cpdgId"],"gen_mm_cpdgId[nDimuons]/I");
            storageMapFloatArray["gen_mm_prod_x"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_prod_x", storageMapFloatArray["gen_mm_prod_x"],"gen_mm_prod_x[nDimuons]/F");
            storageMapFloatArray["gen_mm_prod_y"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_prod_y", storageMapFloatArray["gen_mm_prod_y"],"gen_mm_prod_y[nDimuons]/F");
            storageMapFloatArray["gen_mm_prod_z"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_prod_z", storageMapFloatArray["gen_mm_prod_z"],"gen_mm_prod_z[nDimuons]/F");
            storageMapFloatArray["gen_mm_vtx_x"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_vtx_x", storageMapFloatArray["gen_mm_vtx_x"],"gen_mm_vtx_x[nDimuons]/F");
            storageMapFloatArray["gen_mm_vtx_y"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_vtx_y", storageMapFloatArray["gen_mm_vtx_y"],"gen_mm_vtx_y[nDimuons]/F");
            storageMapFloatArray["gen_mm_vtx_z"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_vtx_z", storageMapFloatArray["gen_mm_vtx_z"],"gen_mm_vtx_z[nDimuons]/F");
            storageMapFloatArray["gen_mm_l3d"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_l3d", storageMapFloatArray["gen_mm_l3d"],"gen_mm_l3d[nDimuons]/F");
            storageMapFloatArray["gen_mm_lxy"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_lxy", storageMapFloatArray["gen_mm_lxy"],"gen_mm_lxy[nDimuons]/F");
            storageMapFloatArray["gen_mm_tau"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_tau", storageMapFloatArray["gen_mm_tau"],"gen_mm_tau[nDimuons]/F");
            storageMapFloatArray["gen_mm_alpha_p_phi"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_alpha_p_phi", storageMapFloatArray["gen_mm_alpha_p_phi"],"gen_mm_alpha_p_phi[nDimuons]/F");
            storageMapFloatArray["gen_mm_alpha_p_theta"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_alpha_p_theta", storageMapFloatArray["gen_mm_alpha_p_theta"],"gen_mm_alpha_p_theta[nDimuons]/F");
            storageMapFloatArray["gen_mm_alpha_ip"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_alpha_ip", storageMapFloatArray["gen_mm_alpha_ip"],"gen_mm_alpha_ip[nDimuons]/F");
            storageMapFloatArray["gen_mm_alpha_vtx"]   = new Float_t[N_COMPOSIT_PART_MAX];
            theTree->Branch("gen_mm_alpha_vtx", storageMapFloatArray["gen_mm_alpha_vtx"],"gen_mm_alpha_vtx[nDimuons]/F");
        }
}





void BsToMuMuGammaNTuplizer::fillDimuonBranches( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<std::vector<reco::Muon>> muonHandle;
    iEvent.getByToken(muonToken_, muonHandle);
    // Fills tree branches with photons.
    edm::Handle<std::vector<reco::Photon> > photonHandle;
    edm::Handle<reco::SuperClusterCollection> barrelSCHandle;
    //edm::Handle<edm::View<reco::PFCandidate> > pfCandidateHandle;
    edm::Handle<reco::SuperClusterCollection> endcapSCHandle;

    if(doSuperClusters_)
    { 
     iEvent.getByToken(MustacheSCBarrelCollection_, barrelSCHandle);
     iEvent.getByToken(MustacheSCEndcapCollection_, endcapSCHandle);
    }

    if(doPhotons_) iEvent.getByToken(gedPhotonsCollection_, photonHandle);
  


    
    AnalyticalImpactPointExtrapolator extrapolator(bFieldHandle_.product());
    
    impactPointExtrapolator_ = &extrapolator;

    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BsToMuMuGammaNTuplizer") << "No beam spot available from EventSetup" ;
    }
    
    

    
    auto nMuons   = muonHandle->size();
    auto nPhotons = photonHandle->size();
    auto nPFCands = pfCandidateHandle->size();
    
    //  setting the JPsi and Bmmg cands
    if(doJPsiGamma) { 
        storageMapInt["nJPsiGammaCands"]=0;
        if(doPFPhotons_)
        {
            storageMapInt["nJPsiGammaPFCands"]=0;
        }
    }
    if(doMuMuGamma) {
        storageMapInt["nMMGCands"]=0;

        if(doPFPhotons_)
        {
            storageMapInt["nMMGpfCands"]=0;
        }
    }
    storageMapInt["nMuMuK"]=0;
    storageMapInt["nJPsiK"]=0;
    storageMapInt["nPsi2SK"]=0;

    // Output collection
    auto dimuon  = std::make_unique<pat::CompositeCandidateCollection>();
    auto btommg  = std::make_unique<pat::CompositeCandidateCollection>();
    AddFourMomenta addP4;

    // good muon candidates
    std::vector<MuonCand> good_muon_candidates;
    for (unsigned int i = 0; i < nMuons; ++i) {
      const pat::Muon & muon = muonHandle->at(i);
      if (not isGoodMuon(muon)) continue;
      good_muon_candidates.push_back(MuonCand(muon, i));
    }
    
    // Inject B to hh candidates where hadrons are explicitely matched
    // to gen level decays 
    // if ( injectMatchedBtohh_ and isMC_ ) {
    //  injectHadronsThatMayFakeMuons(good_muon_candidates);
    //}

    // Inject reco B to hh candidates
    //if ( injectBtohh_ ) {
    //  injectBhhHadrons(good_muon_candidates);
    //}
    
    // Build dimuon candidates

    Int_t nDimu=0;
    if ( good_muon_candidates.size() > 1 ){
      for (unsigned int i = 0; i < good_muon_candidates.size(); ++i) {
	const MuonCand & muon1 = good_muon_candidates.at(i);
	for (unsigned int j = 0; j < good_muon_candidates.size(); ++j) {
	  if (i==j) continue;
	  const MuonCand & muon2 = good_muon_candidates.at(j);
	  // Ensure that muon1.pt > muon2.pt
	  if (muon2.pt() > muon1.pt()) continue;

	  auto mm_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
						   muon2.innerTrack().get());
	  if (maxTwoTrackDOCA_>0 and mm_doca > maxTwoTrackDOCA_) continue;
	  if (diMuonCharge_ && muon1.charge() * muon2.charge()>0) continue;

	  pat::CompositeCandidate dimuonCand;
	  dimuonCand.addDaughter( muon1 , "muon1");
	  dimuonCand.addDaughter( muon2 , "muon2");
	  addP4.set( dimuonCand );

	  // reco Btohh
	  //if (muon1.index() < 0 and not muon1.from_gen() and 
	  //    muon2.index() < 0 and not muon2.from_gen())
	  //  {
	  //    if (dimuonCand.mass() < minBhhMass_) continue;
	  //    if (dimuonCand.mass() > maxBhhMass_) continue;
	  //  }
	storageMapIntArray["dimuon_mu1_index"][nDimu]  =  muon1.index();
	storageMapIntArray["dimuon_mu1_pdgId"][nDimu]  =  muon1.pdgId();
	storageMapFloatArray["dimuon_mu1_pt"][nDimu]  =   muon1.pt();
	storageMapFloatArray["dimuon_mu1_eta"][nDimu]  =  muon1.eta();
	storageMapFloatArray["dimuon_mu1_phi"][nDimu]  =  muon1.phi();
	storageMapIntArray["dimuon_mu2_index"][nDimu]  =  muon2.index();
	storageMapIntArray["dimuon_mu2_pdgId"][nDimu]  =  muon2.pdgId();
	storageMapFloatArray["dimuon_mu2_pt"][nDimu]  =   muon2.pt();
	storageMapFloatArray["dimuon_mu2_eta"][nDimu]  =  muon2.eta();
	storageMapFloatArray["dimuon_mu2_phi"][nDimu]  =  muon2.phi();
	storageMapFloatArray["dimuon_doca"][nDimu]  =  mm_doca;

	  // Kinematic Fits
	  // auto kinematicMuMuVertexFit = fillMuMuInfo(dimuonCand,iEvent,muon1,muon2);

	  // reco Btohh
	  // if (muon1.index() < 0 and not muon1.from_gen() and 
	  //     muon2.index() < 0 and not muon2.from_gen())
	  //   {
	  //     if (not kinematicMuMuVertexFit.valid()) continue;
	  //     if (kinematicMuMuVertexFit.vtxProb() < minBhhVtxProb_) continue;
	  //     if (kinematicMuMuVertexFit.sigLxy < minBhhSigLxy_) continue;
	  //   }

      auto dimuFit = vertexMuonsWithKinematicFitter(muon1, muon2);
      dimuFit.postprocess(*beamSpot);
      auto displacement3d = compute3dDisplacement(dimuFit, *pvHandle_.product(),true);
    
    //addFitInfo(dimuonCand, kinematicMuMuVertexFit, "kin", displacement3D,0,1);
    storageMapFloatArray["dimuon_kin_valid"][nDimu]  =        dimuFit.valid() ;
    storageMapFloatArray["dimuon_kin_vtx_prob"][nDimu]  =     dimuFit.vtxProb() ;
    storageMapFloatArray["dimuon_kin_vtx_chi2dof"][nDimu]  =  dimuFit.chi2()>0?dimuFit.chi2()/dimuFit.ndof():-1;
    storageMapFloatArray["dimuon_kin_mass"][nDimu]  =         dimuFit.mass() ;
    storageMapFloatArray["dimuon_kin_massErr"][nDimu]  =      dimuFit.massErr() ;
    storageMapFloatArray["dimuon_kin_lxy"][nDimu]  =          dimuFit.lxy ;
    storageMapFloatArray["dimuon_kin_sigLxy"][nDimu]  =       dimuFit.sigLxy ;
    storageMapFloatArray["dimuon_kin_alphaBS"][nDimu]  =      dimuFit.alphaBS;
    storageMapFloatArray["dimuon_kin_alphaBSErr"][nDimu]  =   dimuFit.alphaBSErr;
    storageMapFloatArray["dimuon_kin_vtx_x"][nDimu]  =        dimuFit.valid()?dimuFit.refitVertex->position().x():0 ;
    storageMapFloatArray["dimuon_kin_vtx_xErr"][nDimu]  =     dimuFit.valid()?sqrt(dimuFit.refitVertex->error().cxx()):0 ;
    storageMapFloatArray["dimuon_kin_vtx_y"][nDimu]  =        dimuFit.valid()?dimuFit.refitVertex->position().y():0 ;
    storageMapFloatArray["dimuon_kin_vtx_yErr"][nDimu]  =     dimuFit.valid()?sqrt(dimuFit.refitVertex->error().cyy()):0 ;
    storageMapFloatArray["dimuon_kin_vtx_z"][nDimu]  =        dimuFit.valid()?dimuFit.refitVertex->position().z():0 ;
    storageMapFloatArray["dimuon_kin_vtx_zErr"][nDimu]  =     dimuFit.valid()?sqrt(dimuFit.refitVertex->error().czz()):0 ;
    storageMapFloatArray["dimuon_kin_pt"][nDimu]  =           dimuFit.p3().perp() ;
    storageMapFloatArray["dimuon_kin_eta"][nDimu]  =          dimuFit.p3().eta() ;
    storageMapFloatArray["dimuon_kin_phi"][nDimu]  =          dimuFit.p3().phi() ;
    storageMapFloatArray["dimuon_kin_mu1pt"][nDimu]   =       dimuFit.dau_p3(0).perp();
    storageMapFloatArray["dimuon_kin_mu1eta"][nDimu]  =       dimuFit.dau_p3(0).eta() ;
    storageMapFloatArray["dimuon_kin_mu1phi"][nDimu]  =       dimuFit.dau_p3(0).phi() ;
    storageMapFloatArray["dimuon_kin_mu1pt"][nDimu]   =       dimuFit.dau_p3(1).perp();
    storageMapFloatArray["dimuon_kin_mu1eta"][nDimu]  =       dimuFit.dau_p3(1).eta() ;
    storageMapFloatArray["dimuon_kin_mu1phi"][nDimu]  =       dimuFit.dau_p3(1).phi() ;
   
    // IP info
    storageMapFloatArray["dimuon_kin_alpha"][nDimu]  =        displacement3d.alpha;
    storageMapFloatArray["dimuon_kin_alphaErr"][nDimu]  =     displacement3d.alphaErr;
    storageMapFloatArray["dimuon_kin_l3d"][nDimu]  =          displacement3d.decayLength;
    storageMapFloatArray["dimuon_kin_sl3d"][nDimu]  =         displacement3d.decayLengthErr>0?displacement3d.decayLength/displacement3d.decayLengthErr:0;
    storageMapFloatArray["dimuon_kin_pv_z"][nDimu]  =         displacement3d.pv?displacement3d.pv->position().z():0;
    storageMapFloatArray["dimuon_kin_pv_zErr"][nDimu]  =      displacement3d.pv?displacement3d.pv->zError():0;
    storageMapFloatArray["dimuon_kin_pvip"][nDimu]  =         displacement3d.distaceOfClosestApproach;
    storageMapFloatArray["dimuon_kin_spvip"][nDimu]  =        displacement3d.distaceOfClosestApproachSig;
    storageMapFloatArray["dimuon_kin_pvipErr"][nDimu]  =      displacement3d.distaceOfClosestApproachErr;
    storageMapFloatArray["dimuon_kin_pv2ip"][nDimu]  =        displacement3d.distaceOfClosestApproach2;
    storageMapFloatArray["dimuon_kin_spv2ip"][nDimu]  =       displacement3d.distaceOfClosestApproach2Sig;
    storageMapFloatArray["dimuon_kin_pv2ipErr"][nDimu]  =     displacement3d.distaceOfClosestApproach2Err;
    storageMapFloatArray["dimuon_kin_pvlip"][nDimu]  =        displacement3d.longitudinalImpactParameter;
    storageMapFloatArray["dimuon_kin_pvlipSig"][nDimu]  =     displacement3d.longitudinalImpactParameterSig;
    storageMapFloatArray["dimuon_kin_pvlipErr"][nDimu]  =     displacement3d.longitudinalImpactParameterErr;
    storageMapFloatArray["dimuon_kin_pv2lip"][nDimu]  =       displacement3d.longitudinalImpactParameter2;
    storageMapFloatArray["dimuon_kin_pv2lipSig"][nDimu]  =    displacement3d.longitudinalImpactParameter2Sig;
    storageMapFloatArray["dimuon_kin_pv2lipErr"][nDimu]  =    displacement3d.longitudinalImpactParameter2Err;
    storageMapFloatArray["dimuon_kin_pvIndex"][nDimu]  =      displacement3d.pvIndex;

    // DecayTime
    storageMapFloatArray["dimuon_kin_tau"][nDimu]  =          displacement3d.decayTime;
    storageMapFloatArray["dimuon_kin_taue"][nDimu]  =         displacement3d.decayTimeError;
    storageMapFloatArray["dimuon_kin_tauxy"][nDimu]  =        displacement3d.decayTimeXY;
    storageMapFloatArray["dimuon_kin_tauxye"][nDimu]  =       displacement3d.decayTimeXYError;
 
    int pvIndex = displacement3d.pvIndex;
    const reco::Vertex *vertex(nullptr) ;
    if(pvIndex >=0 ) vertex=&(pvHandle_->at(pvIndex));
    // Look for additional tracks compatible with the dimuon vertex
    auto closeTracks = findTracksCompatibleWithTheVertex(muon1,muon2,dimuFit);
    //closeTracks.fillCandInfo(dimuonCand, pvIndex, "");
    storageMapFloatArray["dimuon_nTrks"][nDimu]  =        closeTracks.nTracksByVertexProbability(0.1, -1.0, vertex);
    storageMapFloatArray["dimuon_nBMTrks"][nDimu]  =      closeTracks.nTracksByBetterMatch();
    storageMapFloatArray["dimuon_nDisTrks"][nDimu]  =     closeTracks.nTracksByVertexProbability(0.1,  2.0 ,  vertex);
    storageMapFloatArray["dimuon_closetrk"][nDimu]  =     closeTracks.nTracksByDisplacementSignificance(0.03 ,-1 ,  vertex);
    storageMapFloatArray["dimuon_closetrks1"][nDimu]  =   closeTracks.nTracksByDisplacementSignificance(0.03 , 1 ,  vertex);
    storageMapFloatArray["dimuon_closetrks2"][nDimu]  =   closeTracks.nTracksByDisplacementSignificance(0.03 , 2 ,  vertex);
    storageMapFloatArray["dimuon_closetrks3"][nDimu]  =   closeTracks.nTracksByDisplacementSignificance(0.03 , 3 ,  vertex);
    storageMapFloatArray["dimuon_docatrk"][nDimu]  =      closeTracks.minDoca(0.03,vertex);
    storageMapFloatArray["dimuon_m1iso"][nDimu]             =   computeTrkMuonIsolation(muon1,muon2,pvIndex,0.5,0.5);
    storageMapFloatArray["dimuon_m2iso"][nDimu]             =   computeTrkMuonIsolation(muon2,muon1,pvIndex,0.5,0.5);
    storageMapFloatArray["dimuon_iso"][nDimu]               =   computeTrkMuMuIsolation(muon2,muon1,pvIndex,0.9,0.7);
    storageMapFloatArray["dimuon_otherVtxMaxProb"][nDimu]   =   otherVertexMaxProb(muon1,muon2,0.5);
    storageMapFloatArray["dimuon_otherVtxMaxProb1"][nDimu]  =   otherVertexMaxProb(muon1,muon2,1.0);
    storageMapFloatArray["dimuon_otherVtxMaxProb2"][nDimu]  =   otherVertexMaxProb(muon1,muon2,2.0);

    if(isMC)
    {
        auto gen_mm = getGenMatchInfo(muon1,muon2);
        storageMapIntArray["gen_mm_mu1_pdgId"][nDimu]    = gen_mm.mu1_pdgId                                    		;
        storageMapIntArray["gen_mm_mu1_mpdgId"][nDimu]   = gen_mm.mu1_motherPdgId                               		;
        storageMapFloatArray["gen_mm_mu1_pt"][nDimu]   	  = gen_mm.mu1_pt                                          	;
        storageMapIntArray["gen_mm_mu2_pdgId"][nDimu]    = gen_mm.mu2_pdgId                                    		;
        storageMapIntArray["gen_mm_mu2_mpdgId"][nDimu]   = gen_mm.mu2_motherPdgId                               		;
        storageMapFloatArray["gen_mm_mu2_pt"][nDimu]   	  = gen_mm.mu2_pt                                          ;
        storageMapFloatArray["gen_mm_mass"][nDimu]   	  = gen_mm.mm_mass                                         ;
        storageMapFloatArray["gen_mm_pt"][nDimu]   	      = gen_mm.mm_pt                                           ;
        storageMapIntArray["gen_mm_pdgId"][nDimu]   	  = gen_mm.mm_pdgId                                        ;
        storageMapIntArray["gen_mm_mpdgId"][nDimu]   	  = gen_mm.mm_motherPdgId                                  ;
        storageMapIntArray["gen_mm_cpdgId"][nDimu]   	  = gen_mm.common_mother?gen_mm.common_mother->pdgId():0   ;
        storageMapFloatArray["gen_mm_prod_x"][nDimu]   	  = gen_mm.mm_prod_vtx.x()                                 ;
        storageMapFloatArray["gen_mm_prod_y"][nDimu]   	  = gen_mm.mm_prod_vtx.y()                                 ;
        storageMapFloatArray["gen_mm_prod_z"][nDimu]   	  = gen_mm.mm_prod_vtx.z()                                 ;
        storageMapFloatArray["gen_mm_vtx_x"][nDimu]   	  = gen_mm.mm_vtx.x()                                      ;
        storageMapFloatArray["gen_mm_vtx_y"][nDimu]   	  = gen_mm.mm_vtx.y()                                      ;
        storageMapFloatArray["gen_mm_vtx_z"][nDimu]   	  = gen_mm.mm_vtx.z()                                      ;
        storageMapFloatArray["gen_mm_l3d"][nDimu]   		  = (gen_mm.mm_prod_vtx-gen_mm.mm_vtx).r()                 ;
        storageMapFloatArray["gen_mm_lxy"][nDimu]   		  = (gen_mm.mm_prod_vtx-gen_mm.mm_vtx).rho()               ;
        storageMapFloatArray["gen_mm_tau"][nDimu]   		  = computeDecayTime(gen_mm)                               ;
        if (gen_mm.match and dimuFit.valid()){
          storageMapFloatArray["gen_mm_alpha_p_phi"][nDimu]  = dimuFit.refitMother->currentState().globalMomentum().phi() - gen_mm.match->phi();
          storageMapFloatArray["gen_mm_alpha_p_theta"][nDimu]= dimuFit.refitMother->currentState().globalMomentum().phi() - gen_mm.match->theta();
          TVector3 p_gen(gen_mm.match->px(), gen_mm.match->py(),   gen_mm.match->pz());
          TVector3 ip_reco(displacement3d.pv->x(),  displacement3d.pv->y(),    displacement3d.pv->z());
          TVector3 ip_gen(gen_mm.mm_prod_vtx.x(),   gen_mm.mm_prod_vtx.y(),    gen_mm.mm_prod_vtx.z());
          TVector3 vtx_reco(dimuFit.refitVertex->vertexState().position().x(), 
    			dimuFit.refitVertex->vertexState().position().y(), 
    			dimuFit.refitVertex->vertexState().position().z());
          TVector3 vtx_gen(gen_mm.mm_vtx.x(),
    		       gen_mm.mm_vtx.y(),
    		       gen_mm.mm_vtx.z());
          float cosAlpha_ip  = p_gen.Dot(vtx_gen - ip_reco) / (p_gen.Mag() * (vtx_gen - ip_reco).Mag());
          float cosAlpha_vtx = p_gen.Dot(vtx_reco - ip_gen) / (p_gen.Mag() * (vtx_reco - ip_gen).Mag());
          
          storageMapFloatArray["gen_mm_alpha_ip"][nDimu]  = acos(cosAlpha_ip);
          storageMapFloatArray["gen_mm_alpha_vtx"][nDimu] = acos(cosAlpha_vtx);
        } else {
          storageMapFloatArray["gen_mm_alpha_p_phi"][nDimu]  = 999;
          storageMapFloatArray["gen_mm_alpha_p_theta"][nDimu]= 999;
          storageMapFloatArray["gen_mm_alpha_ip"][nDimu]     = 999;
          storageMapFloatArray["gen_mm_alpha_vtx"][nDimu]    = 999;
        }
    }
    
	  if (muon1.index() >= 0 and muon2.index() >= 0){
   		        auto dimuon_p4(makeLorentzVectorFromPxPyPzM(dimuFit.p3().x(),
		    					    dimuFit.p3().y(),
		    					    dimuFit.p3().z(),
		    					    dimuFit.mass()));
		    const auto & vtx_point = dimuFit.refitVertex->vertexState().position();
            math::XYZVectorF  secVtx(vtx_point.x(),vtx_point.y(),vtx_point.z());
            auto isJpsiCand = ( dimuFit.mass() >= minJPsiMass_ and dimuFit.mass() <= maxJPsiMass_ );

	        // MuMuGamma
	        if (doMuMuGamma && dimuFit.valid()){
	        for (unsigned int k=0; k < nPhotons; ++k){
		        auto photon(photonHandle->at(k));
		        photon.setVertex(reco::Photon::Point(vtx_point.x(), vtx_point.y(), vtx_point.z()));
		        
                double mmg_mass = (dimuon_p4 + photon.p4()).mass();
                
                 GenMatchInfo *gen_mmg(nullptr) ;
                 if(isMC) {
                            reco::PFCandidate pf;
                            pf.setP4(photon.p4());
                           auto gen_mmg_ = getGenMatchInfo(muon1,muon2,nullptr,nullptr,&pf);
                           gen_mmg=&gen_mmg_;
                         }
		        
                if (mmg_mass >= minBsMuMuGammaMass_ and mmg_mass <= maxBsMuMuGammaMass_){
                     fillBmmgBranchs(nDimu,k,-1,mmg_mass,gen_mmg);
	  	         }
                if( doJPsiGamma and isJpsiCand)
                {
                    if( mmg_mass >= minJPsiGammaMass_ and mmg_mass <= maxJPsiGammaMass_)
                    {
                        fillJPsiGammaBranches(nDimu,k,-1,mmg_mass,gen_mmg);
                    }
                }
	          }
            //   make the mmg with SC
            if(doSuperClusters_)
            {
             Int_t k=0;
             for (auto const& scs : { *barrelSCHandle.product(), *endcapSCHandle.product() }) {
                for (auto const& sc : scs) {
                  
                 math::XYZVectorF  scPos(sc.position().x(),sc.position().y(),sc.position().z());
                 auto scP3 = (scPos - secVtx).unit()*sc.energy();
                 auto scP4(makeLorentzVectorFromPxPyPzM(scP3.x(),scP3.y(),scP3.z(),0.0));
		         double mmg_mass = (dimuon_p4 + scP4).mass();
                 
                 GenMatchInfo *gen_mmg(nullptr) ;
                 if(isMC) {
                            reco::PFCandidate pf;
                            pf.setP4(scP4);
                           auto gen_mmg_ = getGenMatchInfo(muon1,muon2,nullptr,nullptr,&pf);
                           gen_mmg=&gen_mmg_;
                         }
                 
                if (mmg_mass >= minBsMuMuGammaMass_ and mmg_mass <= maxBsMuMuGammaMass_){
                    fillBmmgBranchs(nDimu,-1,k,mmg_mass,gen_mmg);
	  	         }
                if( doJPsiGamma and isJpsiCand)
                  {
                    if( mmg_mass >= minJPsiGammaMass_ and mmg_mass <= maxJPsiGammaMass_)
                    {
                       fillJPsiGammaBranches(nDimu,-1,k,mmg_mass,gen_mmg);
                    }
                  }
                  k++;
                } // supercluster loop
                }
             }

            //   make the mmg with PF
            if(doPFPhotons_)
            {
             Int_t k=0;

               for (auto pf = pfCandidateHandle->begin(); pf != pfCandidateHandle->end(); ++pf) {
                  
                  if(pf->pdgId() !=22)continue;
                  if(pf->et() < pTMinPFPhotons)continue;
                  
                  math::XYZVectorF  candPos(pf->vx(),pf->vy(),pf->vz());
                  auto candP3 = (candPos - secVtx).unit()*pf->energy();
                  auto candP4(makeLorentzVectorFromPxPyPzM(candP3.x(),candP3.y(),candP3.z(),0.0));
		          double mmg_mass = (dimuon_p4 + candP4).mass();
                  
               GenMatchInfo *gen_mmg(nullptr) ;
               if(isMC) 
               {
                   auto gen_mmg_ = getGenMatchInfo(muon1,muon2,nullptr,nullptr,&(*pf));
                   gen_mmg=&gen_mmg_;
               }

                if (mmg_mass >= minBsMuMuGammaMass_ and mmg_mass <= maxBsMuMuGammaMass_){
                       fillPF_BmmgBranchs(nDimu,k,mmg_mass,gen_mmg);
	  	         }
                if( doJPsiGamma and isJpsiCand)
                  {
                    if( mmg_mass >= minJPsiGammaMass_ and mmg_mass <= maxJPsiGammaMass_)
                    {
                        fillPF_JPsiGammaBranches(nDimu,k,mmg_mass,gen_mmg);
                    }
                  }
                  k++;
                }              
               } // PF photons loop

             }
            
	  
              // dimuon
         
        	    // mmK and mmKK
                
                if(doMuMuK_)
                {
        	        for (unsigned int k = 0; k < nPFCands; ++k) {
        	          reco::PFCandidate kaonCand1((*pfCandidateHandle)[k]);
        	          kaonCand1.setMass(KaonMass_);
        	          if (deltaR(muon1, kaonCand1) < 0.01 || deltaR(muon2, kaonCand1) < 0.01) continue;
        	          if (kaonCand1.charge() == 0 ) continue;
        	          //if (!kaonCand1.hasTrackDetails()) continue;
                      if (not kaonCand1.trackRef()) continue;
        	          double mu1_kaon_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
        		    					       kaonCand1.bestTrack());
        	          double mu2_kaon_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
        		    					       kaonCand1.bestTrack());
        	          // BtoMuMuK
        	          bool goodBtoMuMuK = true;
        	          if (dimuFit.mass() < 2.9) goodBtoMuMuK = false;
        	          if (abs(kaonCand1.pdgId()) != 211) goodBtoMuMuK = false; //Charged hadrons
        	          if (kaonCand1.pt() < ptMinKaon_ || abs(kaonCand1.eta()) > etaMaxKaon_) goodBtoMuMuK = false;
        	          if (maxTwoTrackDOCA_ > 0 and mu1_kaon_doca > maxTwoTrackDOCA_) goodBtoMuMuK = false;
        	          if (maxTwoTrackDOCA_ > 0 and mu2_kaon_doca > maxTwoTrackDOCA_) goodBtoMuMuK = false;
        	          
        
                      double kmm_mass = (muon1.p4() + muon2.p4() + kaonCand1.p4()).mass();
        	          if (kmm_mass < minBKmmMass_ || kmm_mass > maxBKmmMass_) goodBtoMuMuK = false;
        	        
        	          if (goodBtoMuMuK){
        		          // fill BtoMuMuK candidate info
                          fillBtoMuMuKInfo(nDimu,iEvent,dimuFit,muon1,muon2,kaonCand1);
        	          }
        
        	          // BToJpsiKK
        	          // bool goodBtoJpsiKK = goodBtoMuMuK;
        	          // if (fabs(dimuFit.mass()-3.1) > 0.2) goodBtoJpsiKK = false;
        
        	          // if (goodBtoJpsiKK){ // good candidate to consider for JpsiKK
        		      //    for (unsigned int k2 = k+1; k2 < nPFCands; ++k2) { // only works if selection requirements for both kaons are identical
        		      //    reco::PFCandidate kaonCand2((*pfCandidateHandle)[k2]);
        		      //    kaonCand2.setMass(KaonMass_);
        		      // if (deltaR(muon1, kaonCand2) < 0.01 || deltaR(muon2, kaonCand2) < 0.01) continue;
        		      // if (kaonCand2.charge() == 0 ) continue;
        		      // if (!kaonCand2.hasTrackDetails()) continue;
        		      // double mu1_kaon2_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
        		      //   					    kaonCand2.bestTrack());
        		      // double mu2_kaon2_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
        		      //   					    kaonCand2.bestTrack());
        	          // 
        		      // if (abs(kaonCand2.pdgId())!=211) goodBtoJpsiKK = false; //Charged hadrons
        		      // if (kaonCand2.pt()<ptMinKaon_ || abs(kaonCand2.eta())>etaMaxKaon_) goodBtoJpsiKK = false;
        		      // if (maxTwoTrackDOCA_>0 and mu1_kaon2_doca> maxTwoTrackDOCA_) goodBtoJpsiKK = false;	      
        		      // if (maxTwoTrackDOCA_>0 and mu2_kaon2_doca> maxTwoTrackDOCA_) goodBtoJpsiKK = false;	      
        		      // 
        		      // double kkmm_mass = (muon1.p4()+muon2.p4()+kaonCand1.p4()+kaonCand2.p4()).mass();
        		      // if ( kkmm_mass<minBKKmmMass_ || kkmm_mass>maxBKKmmMass_ ) goodBtoJpsiKK = false;
        		      // 
        		      // if (goodBtoJpsiKK){
        		      //   // fill BtoJpsiKK candidate info
        		      //   
        		      //   pat::CompositeCandidate btokkmmCand;
        		      //   btokkmmCand.addUserInt("mm_index", imm);
        		      //   int ikmm = -1;
        		      //   if (goodBtoMuMuK) ikmm = btokmm->size()-1;
        		      //   btokkmmCand.addUserInt("kmm_index", ikmm);
        		      //   btokkmmCand.addUserFloat("kaon_mu1_doca", mu1_kaon2_doca);
        		      //   btokkmmCand.addUserFloat("kaon_mu2_doca", mu2_kaon2_doca);
        		      //   
        		      //   fillBstoJpsiKKInfo(btokkmmCand,iEvent,kinematicMuMuVertexFit,muon1,muon2,kaonCand1,kaonCand2);
                      //   fillBtoMuMuKInfo(nDimu,btokkmmCand,iEvent,kinematicMuMuVertexFit,muon1,muon2,kaonCand1,kaonCand2);
        		      //   // FIXME
        		      //   // fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(btokkmmCand,dimuonCand,iEvent,kinematicMuMuVertexFit,muon1,muon2,kaonCand1);
        
        		      //   btokkmm->push_back(btokkmmCand);
        		      //      }
        		      //    }
        	          // }
        	        }
                }
        
            // printf("kinematicMuMuVertexFit (x,y,z): (%7.3f,%7.3f,%7.3f)\n", 
        
            }

          nDimu++;
	   }
      }
     }
    
    storageMapInt["nDimuons"]=nDimu;

}


void BsToMuMuGammaNTuplizer::addGenBranches(TString sTag,TString nTag)
{
       std::string  tag(sTag);

       storageMapIntArray["gen_"+tag+"_cpdgId"]   = new Int_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_cpdgId", storageMapIntArray["gen_"+tag+"_cpdgId"],"gen_"+sTag+"_cpdgId[nMuMuK]/I");
       storageMapFloatArray["gen_"+tag+"_l3d"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_l3d", storageMapFloatArray["gen_"+tag+"_l3d"],"gen_"+sTag+"_l3d["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_lxy"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_lxy", storageMapFloatArray["gen_"+tag+"_lxy"],"gen_"+sTag+"_lxy["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_mass"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_mass", storageMapFloatArray["gen_"+tag+"_mass"],"gen_"+sTag+"_mass["+nTag+"]/F");
       storageMapIntArray["gen_"+tag+"_pdgId"]   = new Int_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_pdgId", storageMapIntArray["gen_"+tag+"_pdgId"],"gen_"+sTag+"_pdgId["+nTag+"]/I");
       storageMapFloatArray["gen_"+tag+"_prod_x"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_prod_x", storageMapFloatArray["gen_"+tag+"_prod_x"],"gen_"+sTag+"_prod_x["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_prod_y"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_prod_y", storageMapFloatArray["gen_"+tag+"_prod_y"],"gen_"+sTag+"_prod_y["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_prod_z"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_prod_z", storageMapFloatArray["gen_"+tag+"_prod_z"],"gen_"+sTag+"_prod_z["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_pt"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_pt", storageMapFloatArray["gen_"+tag+"_pt"],"gen_"+sTag+"_pt["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_tau"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_tau", storageMapFloatArray["gen_"+tag+"_tau"],"gen_"+sTag+"_tau["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_alpha_p_phi"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_alpha_p_phi", storageMapFloatArray["gen_"+tag+"_alpha_p_phi"],"gen_"+sTag+"_alpha_p_phi["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_alpha_p_theta"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_alpha_p_theta", storageMapFloatArray["gen_"+tag+"_alpha_p_theta"],"gen_"+sTag+"_alpha_p_theta["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_alpha_ip"]   = new Float_t[N_GEN_PARTS_ALL];
       theTree->Branch("gen_"+sTag+"_alpha_ip", storageMapFloatArray["gen_"+tag+"_alpha_ip"],"gen_"+sTag+"_alpha_ip["+nTag+"]/F");
       storageMapFloatArray["gen_"+tag+"_alpha_vtx"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_"+sTag+"_alpha_vtx", storageMapFloatArray["gen_"+tag+"_alpha_vtx"],"gen_"+sTag+"_alpha_vtx["+nTag+"]/F");
}


void BsToMuMuGammaNTuplizer::addMMGBranches()
{
      storageMapInt["nMMGCands"]=0;
      theTree->Branch("nMMGCands",   &storageMapInt["nMMGCands"]);
      
      storageMapIntArray["mmg_phoIdx"] = new Int_t[NMAX_MMG] ;
      theTree->Branch("mmg_phoIdx",   storageMapIntArray["mmg_phoIdx"],"mmg_phoIdx[nMMGCands]/I");
      storageMapIntArray["mmg_scIdx"] = new Int_t[NMAX_MMG] ;
      theTree->Branch("mmg_scIdx",   storageMapIntArray["mmg_scIdx"],"mmg_scIdx[nMMGCands]/I");
      storageMapFloatArray["mmg_mass"] = new Float_t[NMAX_MMG] ;
      theTree->Branch("mmg_mass",   storageMapFloatArray["mmg_mass"],"mmg_mass[nMMGCands]/F");
      storageMapFloatArray["mmg_diMuMass"] = new Float_t[NMAX_MMG] ;
      theTree->Branch("mmg_diMuMass",   storageMapFloatArray["mmg_diMuMass"],"mmg_diMuMass[nMMGCands]/F");

      if(isMC)
      {
        addGenBranches("mmg","nMMGCands");
        storageMapIntArray["gen_mmg_photon_mpdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_mmg_photon_mpdgId", storageMapIntArray["gen_mmg_photon_mpdgId"],"gen_mmg_photon_mpdgId[nMMGCands]/I");
        storageMapIntArray["gen_mmg_photon_pdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_mmg_photon_pdgId", storageMapIntArray["gen_mmg_photon_pdgId"],"gen_mmg_photon_pdgId[nMMGCands]/I");
        storageMapFloatArray["gen_mmg_photon_pt"]   = new Float_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_mmg_photon_pt", storageMapFloatArray["gen_mmg_photon_pt"],"gen_mmg_photon_pt[nMMGCands]/F");
      }

}



void BsToMuMuGammaNTuplizer:: fillBmmgBranchs(Int_t nDimu,Int_t phoIdx,Int_t scIdx , Float_t mmg_mass, GenMatchInfo *gen_mmg_)
{
      auto idx= storageMapInt["nMMGCands"];
      storageMapIntArray["mmg_phoIdx"][idx]=phoIdx;
      storageMapIntArray["mmg_scIdx"][idx] =scIdx;
      storageMapFloatArray["mmg_mass"][idx]  =mmg_mass;
      storageMapFloatArray["mmg_diMuMass"][idx]  = storageMapFloatArray["dimuon_kin_mass"][nDimu];
      if(isMC)
      {
        //auto gen_mmg = getGenMatchInfo(muon1,muon2,nullptr,nullptr,&photon);
        auto &gen_mmg = *gen_mmg_;
        storageMapIntArray[  "gen_mmg_photon_pdgId" ][idx] =     gen_mmg.photon_pdgId;
        storageMapIntArray[  "gen_mmg_photon_mpdgId"][idx] =     gen_mmg.photon_motherPdgId;
        storageMapFloatArray["gen_mmg_photon_pt"    ][idx] =     gen_mmg.photon_pt;
        storageMapFloatArray["gen_mmg_mass"       ][idx] =     gen_mmg.mmg_mass;
        storageMapFloatArray["gen_mmg_pt"         ][idx] =     gen_mmg.mmg_pt;
        storageMapIntArray[  "gen_mmg_pdgId"      ][idx] =     gen_mmg.mmg_pdgId;
        storageMapFloatArray["gen_mmg_prod_x"     ][idx] =     gen_mmg.mmg_prod_vtx.x();
        storageMapFloatArray["gen_mmg_prod_y"     ][idx] =     gen_mmg.mmg_prod_vtx.y();
        storageMapFloatArray["gen_mmg_prod_z"     ][idx] =     gen_mmg.mmg_prod_vtx.z();
        storageMapFloatArray["gen_mmg_l3d"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).r();
        storageMapFloatArray["gen_mmg_lxy"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).rho();
        storageMapFloatArray["gen_mmg_tau"        ][idx] =     computeDecayTime(gen_mmg);
        storageMapIntArray["gen_mmg_cpdgId"     ][idx] =     gen_mmg.common_mother?gen_mmg.common_mother->pdgId():0;
      }
      storageMapInt["nMMGCands"]+=1;
}

void BsToMuMuGammaNTuplizer::addJPsiGammaBranches()
{
      storageMapInt["nJPsiGammaCands"]=0;
      theTree->Branch("nJPsiGammaCands",   &storageMapInt["nJPsiGammaCands"]);
      storageMapIntArray["jPsiGamma_phoIdx"] = new Int_t[NMAX_MMG] ;
      theTree->Branch("jPsiGamma_phoIdx",   storageMapIntArray["jPsiGamma_phoIdx"],"jPsiGamma_phoIdx[nJPsiGammaCands]/I");
      storageMapIntArray["jPsiGamma_scIdx"] = new Int_t[NMAX_MMG] ;
      theTree->Branch("jPsiGamma_scIdx",   storageMapIntArray["jPsiGamma_scIdx"],"jPsiGamma_scIdx[nJPsiGammaCands]/I");
      storageMapFloatArray["jPsiGamma_mass"] = new Float_t[NMAX_MMG] ;
      theTree->Branch("jPsiGamma_mass",   storageMapFloatArray["jPsiGamma_mass"],"jPsiGamma_phoIdx[nJPsiGammaCands]/f");
      storageMapFloatArray["jPsiGamma_diMuMass"] = new Float_t[NMAX_MMG] ;
      if(isMC)
      {
        addGenBranches("jpsiGamma","nJPsiGammaCands");
        storageMapIntArray["gen_jpsiGamma_photon_mpdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_jpsiGamma_photon_mpdgId", storageMapIntArray["gen_jpsiGamma_photon_mpdgId"],"gen_jpsiGamma_photon_mpdgId[nJPsiGammaCands]/I");
        storageMapIntArray["gen_jpsiGamma_photon_pdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_jpsiGamma_photon_pdgId", storageMapIntArray["gen_jpsiGamma_photon_pdgId"],"gen_jpsiGamma_photon_pdgId[nJPsiGammaCands]/I");
        storageMapFloatArray["gen_jpsiGamma_photon_pt"]   = new Float_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_jpsiGamma_photon_pt", storageMapFloatArray["gen_jpsiGamma_photon_pt"],"gen_jpsiGamma_photon_pt[nJPsiGammaCands]/F");
      }

     theTree->Branch("jPsiGamma_diMuMass",   storageMapFloatArray["jPsiGamma_diMuMass"],"jPsiGamma_diMuMass[nJPsiGammaCands]/f");
}

void BsToMuMuGammaNTuplizer::fillJPsiGammaBranches(Int_t nDimu,Int_t phoIdx,Int_t scIdx , Float_t mmg_mass, GenMatchInfo *gen_mmg_)
{
      storageMapIntArray["jPsiGamma_phoIdx"][storageMapInt["nJPsiGammaCands"]]=phoIdx;
      storageMapIntArray["jPsiGamma_scIdx"][storageMapInt["nJPsiGammaCands"]] =scIdx;
      storageMapFloatArray["jPsiGamma_mass"][storageMapInt["nJPsiGammaCands"]]  =mmg_mass;
      storageMapFloatArray["jPsiGamma_diMuMass"][storageMapInt["nJPsiGammaCands"]]  = storageMapFloatArray["dimuon_kin_mass"][nDimu];
      auto idx=storageMapInt["nJPsiGammaCands"];
      if(isMC)
      {
        auto &gen_mmg = *gen_mmg_;
        storageMapIntArray[  "gen_jpsiGamma_photon_pdgId" ][idx] =     gen_mmg.photon_pdgId;
        storageMapIntArray[  "gen_jpsiGamma_photon_mpdgId"][idx] =     gen_mmg.photon_motherPdgId;
        storageMapFloatArray["gen_jpsiGamma_photon_pt"    ][idx] =     gen_mmg.photon_pt;
        storageMapFloatArray["gen_jpsiGamma_mass"       ][idx] =     gen_mmg.mmg_mass;
        storageMapFloatArray["gen_jpsiGamma_pt"         ][idx] =     gen_mmg.mmg_pt;
        storageMapIntArray[  "gen_jpsiGamma_pdgId"      ][idx] =     gen_mmg.mmg_pdgId;
        storageMapFloatArray["gen_jpsiGamma_prod_x"     ][idx] =     gen_mmg.mmg_prod_vtx.x();
        storageMapFloatArray["gen_jpsiGamma_prod_y"     ][idx] =     gen_mmg.mmg_prod_vtx.y();
        storageMapFloatArray["gen_jpsiGamma_prod_z"     ][idx] =     gen_mmg.mmg_prod_vtx.z();
        storageMapFloatArray["gen_jpsiGamma_l3d"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).r();
        storageMapFloatArray["gen_jpsiGamma_lxy"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).rho();
        storageMapFloatArray["gen_jpsiGamma_tau"        ][idx] =     computeDecayTime(gen_mmg);
        storageMapIntArray["gen_jpsiGamma_cpdgId"     ][idx] =     gen_mmg.common_mother?gen_mmg.common_mother->pdgId():0;
      }
      storageMapInt["nJPsiGammaCands"]+=1;
}


void BsToMuMuGammaNTuplizer::addPF_MMGBranches()
{
      storageMapInt["nMMGpfCands"]=0;
      theTree->Branch("nMMGpfCands",   &storageMapInt["nMMGpfCands"]);
      
      storageMapIntArray["mmgPF_pfPhoIdx"] = new Int_t[NMAX_MMG] ;
      theTree->Branch("mmgPF_pfPhoIdx",   storageMapIntArray["mmgPF_pfPhoIdx"],"mmgPF_pfPhoIdx[nMMGpfCands]/I");
      storageMapFloatArray["mmgPF_mass"] = new Float_t[NMAX_MMG] ;
      theTree->Branch("mmgPF_mass",   storageMapFloatArray["mmgPF_mass"],"mmgPF_mass[nMMGpfCands]/F");
      storageMapFloatArray["mmgPF_diMuMass"] = new Float_t[NMAX_MMG] ;
      theTree->Branch("mmgPF_diMuMass",   storageMapFloatArray["mmgPF_diMuMass"],"mmgPF_diMuMass[nMMGpfCands]/F");

      if(isMC)
      {
        addGenBranches("mmgPF","nMMGCands");
        storageMapIntArray["gen_mmgPF_photon_mpdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_mmgPF_photon_mpdgId", storageMapIntArray["gen_mmgPF_photon_mpdgId"],"gen_mmgPF_photon_mpdgId[nMMGCands]/I");
        storageMapIntArray["gen_mmgPF_photon_pdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_mmgPF_photon_pdgId", storageMapIntArray["gen_mmgPF_photon_pdgId"],"gen_mmgPF_photon_pdgId[nMMGCands]/I");
        storageMapFloatArray["gen_mmgPF_photon_pt"]   = new Float_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_mmgPF_photon_pt", storageMapFloatArray["gen_mmgPF_photon_pt"],"gen_mmgPF_photon_pt[nMMGCands]/F");
      }


}

void BsToMuMuGammaNTuplizer:: fillPF_BmmgBranchs(Int_t nDimu,Int_t phoIdx, Float_t mmg_mass, GenMatchInfo* gen_mmg_)
{

      auto idx=storageMapInt["nMMGpfCands"];
      storageMapIntArray["mmgPF_pfPhoIdx"][idx]=phoIdx;
      storageMapFloatArray["mmgPF_mass"][idx]  =mmg_mass;
      storageMapFloatArray["mmgPF_diMuMass"][idx]  = storageMapFloatArray["dimuon_kin_mass"][nDimu];

      if(isMC)
      {
        auto& gen_mmg= *gen_mmg_;
        storageMapIntArray[  "gen_mmgPF_photon_pdgId" ][idx] =     gen_mmg.photon_pdgId;
        storageMapIntArray[  "gen_mmgPF_photon_mpdgId"][idx] =     gen_mmg.photon_motherPdgId;
        storageMapFloatArray["gen_mmgPF_photon_pt"    ][idx] =     gen_mmg.photon_pt;
        storageMapFloatArray["gen_mmgPF_mass"       ][idx] =     gen_mmg.mmg_mass;
        storageMapFloatArray["gen_mmgPF_pt"         ][idx] =     gen_mmg.mmg_pt;
        storageMapIntArray[  "gen_mmgPF_pdgId"      ][idx] =     gen_mmg.mmg_pdgId;
        storageMapFloatArray["gen_mmgPF_prod_x"     ][idx] =     gen_mmg.mmg_prod_vtx.x();
        storageMapFloatArray["gen_mmgPF_prod_y"     ][idx] =     gen_mmg.mmg_prod_vtx.y();
        storageMapFloatArray["gen_mmgPF_prod_z"     ][idx] =     gen_mmg.mmg_prod_vtx.z();
        storageMapFloatArray["gen_mmgPF_l3d"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).r();
        storageMapFloatArray["gen_mmgPF_lxy"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).rho();
        storageMapFloatArray["gen_mmgPF_tau"        ][idx] =     computeDecayTime(gen_mmg);
        storageMapIntArray["gen_mmgPF_cpdgId"     ][idx] =     gen_mmg.common_mother?gen_mmg.common_mother->pdgId():0;
      }

      storageMapInt["nMMGpfCands"]+=1;

}


void BsToMuMuGammaNTuplizer::addPF_JPsiGammaBranches()
{


      storageMapInt["nJPsiGammaPFCands"]=0;
      theTree->Branch("nJPsiGammaPFCands",   &storageMapInt["nJPsiGammaPFCands"]);
      
      storageMapIntArray["jPsiGammaPF_mumuIdx"] = new Int_t[NMAX_MMG] ;
      theTree->Branch("jPsiGammaPF_mumuIdx",   storageMapIntArray["jPsiGammaPF_mumuIdx"],"jPsiGammaPF_mumuIdx[nJPsiGammaPFCands]/I");
      storageMapIntArray["jPsiGammaPF_pfPhoIdx"] = new Int_t[NMAX_MMG] ;
      theTree->Branch("jPsiGammaPF_pfPhoIdx",   storageMapIntArray["jPsiGammaPF_pfPhoIdx"],"jPsiGammaPF_pfPhoIdx[nJPsiGammaPFCands]/I");
      storageMapFloatArray["jPsiGammaPF_mass"] = new Float_t[NMAX_MMG] ;
      theTree->Branch("jPsiGammaPF_mass",   storageMapFloatArray["jPsiGammaPF_mass"],"jPsiGammaPF_mass[nJPsiGammaPFCands]/f");
      storageMapFloatArray["jPsiGammaPF_diMuMass"] = new Float_t[NMAX_MMG] ;
      theTree->Branch("jPsiGammaPF_diMuMass",   storageMapFloatArray["jPsiGammaPF_diMuMass"],"jPsiGammaPF_diMuMass[nJPsiGammaPFCands]/f");
      
      if(isMC)
      {
        addGenBranches("jpsiGammaPF","nJPsiGammaCands");
        storageMapIntArray["gen_jpsiGammaPF_photon_mpdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_jpsiGammaPF_photon_mpdgId", storageMapIntArray["gen_jpsiGammaPF_photon_mpdgId"],"gen_jpsiGammaPF_photon_mpdgId[nJPsiGammaCands]/I");
        storageMapIntArray["gen_jpsiGammaPF_photon_pdgId"]   = new Int_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_jpsiGammaPF_photon_pdgId", storageMapIntArray["gen_jpsiGammaPF_photon_pdgId"],"gen_jpsiGammaPF_photon_pdgId[nJPsiGammaCands]/I");
        storageMapFloatArray["gen_jpsiGammaPF_photon_pt"]   = new Float_t[N_GEN_PARTS_ALL];
        theTree->Branch("gen_jpsiGammaPF_photon_pt", storageMapFloatArray["gen_jpsiGammaPF_photon_pt"],"gen_jpsiGammaPF_photon_pt[nJPsiGammaCands]/F");
      }


}

void BsToMuMuGammaNTuplizer::fillPF_JPsiGammaBranches(Int_t nDimu,Int_t phoIdx , Float_t mmg_mass, GenMatchInfo * gen_mmg_)
{
      auto idx=storageMapInt["nJPsiGammaPFCands"];
      storageMapIntArray["jPsiGammaPF_mumuIdx"][idx]= nDimu;
      storageMapIntArray["jPsiGammaPF_pfPhoIdx"][idx]=phoIdx;
      storageMapFloatArray["jPsiGammaPF_mass"][idx]  =mmg_mass;
      storageMapFloatArray["jPsiGammaPF_diMuMass"][idx]  = storageMapFloatArray["dimuon_kin_mass"][nDimu];
      
      if(isMC)
      {
        auto& gen_mmg= *gen_mmg_;
        storageMapIntArray[  "gen_jpsiGammaPF_photon_pdgId" ][idx] =     gen_mmg.photon_pdgId;
        storageMapIntArray[  "gen_jpsiGammaPF_photon_mpdgId"][idx] =     gen_mmg.photon_motherPdgId;
        storageMapFloatArray["gen_jpsiGammaPF_photon_pt"    ][idx] =     gen_mmg.photon_pt;
        storageMapFloatArray["gen_jpsiGammaPF_mass"       ][idx] =     gen_mmg.mmg_mass;
        storageMapFloatArray["gen_jpsiGammaPF_pt"         ][idx] =     gen_mmg.mmg_pt;
        storageMapIntArray[  "gen_jpsiGammaPF_pdgId"      ][idx] =     gen_mmg.mmg_pdgId;
        storageMapFloatArray["gen_jpsiGammaPF_prod_x"     ][idx] =     gen_mmg.mmg_prod_vtx.x();
        storageMapFloatArray["gen_jpsiGammaPF_prod_y"     ][idx] =     gen_mmg.mmg_prod_vtx.y();
        storageMapFloatArray["gen_jpsiGammaPF_prod_z"     ][idx] =     gen_mmg.mmg_prod_vtx.z();
        storageMapFloatArray["gen_jpsiGammaPF_l3d"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).r();
        storageMapFloatArray["gen_jpsiGammaPF_lxy"        ][idx] =     (gen_mmg.mmg_prod_vtx-gen_mmg.mm_vtx).rho();
        storageMapFloatArray["gen_jpsiGammaPF_tau"        ][idx] =     computeDecayTime(gen_mmg);
        storageMapIntArray["gen_jpsiGammaPF_cpdgId"     ][idx] =     gen_mmg.common_mother?gen_mmg.common_mother->pdgId():0;
      }
      storageMapInt["nJPsiGammaPFCands"]+=1;

}

void BsToMuMuGammaNTuplizer::fillBtoMuMuKInfo(
                      Int_t mumuIdx ,
					  const edm::Event& iEvent,
					  const KinematicFitResult& kinematicMuMuVertexFit,
					  const pat::Muon& muon1,
					  const pat::Muon& muon2,
					  const reco::PFCandidate & kaon
					) 
{
    std::map<std::string,KinematicFitResult> allFits;
    std::map<std::string,DisplacementInformationIn3D> allDisplacements;
    std::map<std::string,CloseTrackInfo> allCloseTracks;
    
    // MuMuK
    // auto bToKJPsiMuMu_NoMassConstraint = fitBToKJPsiMuMu(kinematicMuMuVertexFit.refitMother, kaon, false);
    allFits["muMuK"] = fitBToKJPsiMuMu(kinematicMuMuVertexFit.refitMother, kaon, false);
    allFits["muMuK"].postprocess(*beamSpot);
    allDisplacements["muMuK"]=compute3dDisplacement(allFits["muMuK"], *pvHandle_.product(),true);
    
    // Look for additional tracks compatible with the dimuon vertex
    // int pvIndex = allDisplacements["muMuK"].pvIndex;
    std::vector<const reco::PFCandidate*>  ignoreTracks;
    ignoreTracks.push_back(&kaon);
    allCloseTracks["muMuK"] = findTracksCompatibleWithTheVertex(muon1,muon2,kinematicMuMuVertexFit,0.03,ignoreTracks);
    Int_t idx=storageMapInt["nMuMuK"];

    storageMapIntArray["muMuK_dimuonIdx"][idx]   = mumuIdx;
    storageMapFloatArray["muMuK_kaonPt"][idx]      = kaon.pt();
    storageMapFloatArray["muMuK_kaonEta"][idx]     = kaon.eta();
    storageMapFloatArray["muMuK_kaonPhi"][idx]     = kaon.phi();
    storageMapFloatArray["muMuK_kaonMass"][idx]    = kaon.mass();
    storageMapFloatArray["muMuK_kaonCharge"][idx]  = kaon.charge();
    storageMapFloatArray["muMuK_kaon_sdxy_bs"][idx]=  kaon.bestTrack()->dxyError()>0 ? fabs(kaon.bestTrack()->dxy(*beamSpot))/kaon.bestTrack()->dxyError():0.0;
    storageMapFloatArray["muMuK_otherVtxMaxProb"][idx]   =   otherVertexMaxProb(muon1,muon2,0.5,0.1,ignoreTracks);
    storageMapFloatArray["muMuK_otherVtxMaxProb1"][idx]  =   otherVertexMaxProb(muon1,muon2,1.0,0.1,ignoreTracks);
    storageMapFloatArray["muMuK_otherVtxMaxProb2"][idx]  =   otherVertexMaxProb(muon1,muon2,2.0,0.1,ignoreTracks);
    fillCompositParticleBranches("muMuK",idx,allFits["muMuK"],allDisplacements["muMuK"],allCloseTracks["muMuK"]);
   if(isMC) {   
        auto gen_kmm = getGenMatchInfo(muon1,muon2,&kaon);
        storageMapIntArray[  "gen_mmk_kaon_pdgId" ][idx] =  gen_kmm.kaon1_pdgId;
        storageMapIntArray[  "gen_mmk_kaon_mpdgId"][idx] =  gen_kmm.kaon1_motherPdgId;
        storageMapFloatArray["gen_mmk_kaon_pt"    ][idx] =  gen_kmm.kaon1_pt;
        storageMapFloatArray["gen_mmk_mass"       ][idx] =  gen_kmm.kmm_mass;
        storageMapFloatArray["gen_mmk_pt"         ][idx] =  gen_kmm.kmm_pt;
        storageMapIntArray[  "gen_mmk_pdgId"      ][idx] =  gen_kmm.kmm_pdgId;
        storageMapFloatArray["gen_mmk_prod_x"     ][idx] =  gen_kmm.kmm_prod_vtx.x();
        storageMapFloatArray["gen_mmk_prod_y"     ][idx] =  gen_kmm.kmm_prod_vtx.y();
        storageMapFloatArray["gen_mmk_prod_z"     ][idx] =  gen_kmm.kmm_prod_vtx.z();
        storageMapFloatArray["gen_mmk_l3d"        ][idx] =  (gen_kmm.kmm_prod_vtx-gen_kmm.mm_vtx).r();
        storageMapFloatArray["gen_mmk_lxy"        ][idx] =  (gen_kmm.kmm_prod_vtx-gen_kmm.mm_vtx).rho();
        storageMapFloatArray["gen_mmk_tau"        ][idx] =  computeDecayTime(gen_kmm);
        storageMapIntArray["gen_mmk_cpdgId"     ][idx] =  gen_kmm.common_mother?gen_kmm.common_mother->pdgId():0;

        if (gen_kmm.match and kinematicMuMuVertexFit.valid()){
          storageMapFloatArray["gen_mmk_alpha_p_phi"][idx]= allFits["muMuK"].refitMother->currentState().globalMomentum().phi() - gen_kmm.match->phi();
          storageMapFloatArray["gen_mmk_alpha_p_theta"][idx]= allFits["muMuK"].refitMother->currentState().globalMomentum().phi() - gen_kmm.match->theta();
          TVector3 p_gen(gen_kmm.match->px(), gen_kmm.match->py(),   gen_kmm.match->pz());
          TVector3 ip_reco(allDisplacements["muMuK"].pv->x(),  allDisplacements["muMuK"].pv->y(),    allDisplacements["muMuK"].pv->z());
          TVector3 ip_gen(gen_kmm.mm_prod_vtx.x(),   gen_kmm.mm_prod_vtx.y(),    gen_kmm.mm_prod_vtx.z());
          TVector3 vtx_reco(allFits["muMuK"].refitVertex->vertexState().position().x(), 
    			allFits["muMuK"].refitVertex->vertexState().position().y(), 
    			allFits["muMuK"].refitVertex->vertexState().position().z());
          TVector3 vtx_gen(gen_kmm.mm_vtx.x(),
    		       gen_kmm.mm_vtx.y(),
    		       gen_kmm.mm_vtx.z());
          float cosAlpha_ip  = p_gen.Dot(vtx_gen - ip_reco) / (p_gen.Mag() * (vtx_gen - ip_reco).Mag());
          float cosAlpha_vtx = p_gen.Dot(vtx_reco - ip_gen) / (p_gen.Mag() * (vtx_reco - ip_gen).Mag());
          
          storageMapFloatArray["gen_mmk_alpha_ip"][idx]  = acos(cosAlpha_ip);
          storageMapFloatArray["gen_mmk_alpha_vtx"][idx] = acos(cosAlpha_vtx);
        } else {
          storageMapFloatArray["gen_mmk_alpha_p_phi"][idx]  = 999;
          storageMapFloatArray["gen_mmk_alpha_p_theta"][idx]= 999;
          storageMapFloatArray["gen_mmk_alpha_ip"][idx]     = 999;
          storageMapFloatArray["gen_mmk_alpha_vtx"][idx]    = 999;
        }
    }
 
    // JPsiK
    auto goodBtoJpsiK = true;
    if (fabs(kinematicMuMuVertexFit.mass()-3.1) > 0.2) goodBtoJpsiK = false;
    if(doJPsiK_ and goodBtoJpsiK) {
      allFits["jpsiK"] = fitBToKMuMu(kinematicMuMuVertexFit.refitTree, kaon, JPsiMass_);
      allFits["jpsiK"].postprocess(*beamSpot);
      allDisplacements["jpsiK"]=compute3dDisplacement(allFits["jpsiK"], *pvHandle_.product(),true);
      allCloseTracks["jpsiK"] = findTracksCompatibleWithTheVertex(muon1,muon2,kinematicMuMuVertexFit,0.03,ignoreTracks);
      idx=storageMapInt["nJPsiK"];

        storageMapIntArray["jpsiK_dimuonIdx"][idx]   = mumuIdx;
        storageMapIntArray["jpsiK_muMuKIdx"][idx]    = storageMapInt["nMuMuK"];
      storageMapFloatArray["jpsiK_kaonPt"][idx]      = kaon.pt();
      storageMapFloatArray["jpsiK_kaonEta"][idx]     = kaon.eta();
      storageMapFloatArray["jpsiK_kaonPhi"][idx]     = kaon.phi();
      storageMapFloatArray["jpsiK_kaonMass"][idx]    = kaon.mass();
      storageMapFloatArray["jpsiK_kaonCharge"][idx]  = kaon.charge();
      storageMapFloatArray["jpsiK_kaon_sdxy_bs"][idx]=  kaon.bestTrack()->dxyError()>0 ? fabs(kaon.bestTrack()->dxy(*beamSpot))/kaon.bestTrack()->dxyError():0.0;
      storageMapFloatArray["jpsiK_otherVtxMaxProb"][idx]   =   otherVertexMaxProb(muon1,muon2,0.5,0.1,ignoreTracks);
      storageMapFloatArray["jpsiK_otherVtxMaxProb1"][idx]  =   otherVertexMaxProb(muon1,muon2,1.0,0.1,ignoreTracks);
      storageMapFloatArray["jpsiK_otherVtxMaxProb2"][idx]  =   otherVertexMaxProb(muon1,muon2,2.0,0.1,ignoreTracks);
      fillCompositParticleBranches("jpsiK",idx,allFits["jpsiK"],allDisplacements["jpsiK"],allCloseTracks["jpsiK"]);
      storageMapInt["nJPsiK"]++;
    }
    
    // Psi2s
    auto goodBtoPsi2sK = true;
    if (fabs(kinematicMuMuVertexFit.mass()-3.68) > 0.2) goodBtoPsi2sK = false;
    if(doPsi2SK_ and goodBtoPsi2sK) {
      allFits["psi2sK"] = fitBToKMuMu(kinematicMuMuVertexFit.refitTree, kaon, JPsiMass_);
      allFits["psi2sK"].postprocess(*beamSpot);
      allDisplacements["psi2sK"]=compute3dDisplacement(allFits["psi2sK"], *pvHandle_.product(),true);
      allCloseTracks["psi2sK"] = findTracksCompatibleWithTheVertex(muon1,muon2,kinematicMuMuVertexFit,0.03,ignoreTracks);
      idx=storageMapInt["nPsi2SK"];

        storageMapIntArray["psi2sK_dimuonIdx"][idx]   = mumuIdx;
        storageMapIntArray["psi2sK_muMuKIdx"][idx]    = storageMapInt["nMuMuK"];
      storageMapFloatArray["psi2sK_kaonPt"][idx]      = kaon.pt();
      storageMapFloatArray["psi2sK_kaonEta"][idx]     = kaon.eta();
      storageMapFloatArray["psi2sK_kaonPhi"][idx]     = kaon.phi();
      storageMapFloatArray["psi2sK_kaonMass"][idx]    = kaon.mass();
      storageMapFloatArray["psi2sK_kaonCharge"][idx]  = kaon.charge();
      storageMapFloatArray["psi2sK_kaon_sdxy_bs"][idx]=  kaon.bestTrack()->dxyError()>0 ? fabs(kaon.bestTrack()->dxy(*beamSpot))/kaon.bestTrack()->dxyError():0.0;
      storageMapFloatArray["psi2sK_otherVtxMaxProb"][idx]   =   otherVertexMaxProb(muon1,muon2,0.5,0.1,ignoreTracks);
      storageMapFloatArray["psi2sK_otherVtxMaxProb1"][idx]  =   otherVertexMaxProb(muon1,muon2,1.0,0.1,ignoreTracks);
      storageMapFloatArray["psi2sK_otherVtxMaxProb2"][idx]  =   otherVertexMaxProb(muon1,muon2,2.0,0.1,ignoreTracks);
      fillCompositParticleBranches("psi2sK",idx,allFits["psi2sK"],allDisplacements["psi2sK"],allCloseTracks["psi2sK"]);
      storageMapInt["nPsi2SK"]++;
    }


    storageMapInt["nMuMuK"]++;

    // // JpsiK
    // KinematicFitResult bToKJPsiMuMu_MassConstraint;
    // DisplacementInformationIn3D bToKJPsiMuMu_MassConstraint_displacement;
    // if (fabs(kinematicMuMuVertexFit.mass()-3.1) < 0.2) {
    //   bToKJPsiMuMu_MassConstraint = fitBToKMuMu(kinematicMuMuVertexFit.refitTree, kaon, JPsiMass_);
    //   bToKJPsiMuMu_MassConstraint.postprocess(*beamSpot_);
    //   bToKJPsiMuMu_MassConstraint_displacement = compute3dDisplacement(bToKJPsiMuMu_MassConstraint, *pvHandle_.product(),true);
    // }
    // addFitInfo(btokmmCand, bToKJPsiMuMu_MassConstraint, "jpsimc", bToKJPsiMuMu_MassConstraint_displacement,-1,-1,1);
  
    // // Psi(2S)K
    // KinematicFitResult bToKPsi2SMuMu_MassConstraint;
    // DisplacementInformationIn3D bToKPsi2SMuMu_MassConstraint_displacement;
    // if (fabs(kinematicMuMuVertexFit.mass()-3.7) < 0.2) {
    //   bToKPsi2SMuMu_MassConstraint = fitBToKMuMu(kinematicMuMuVertexFit.refitTree, kaon, Psi2SMass_);
    //   bToKPsi2SMuMu_MassConstraint.postprocess(*beamSpot_);
    //   bToKPsi2SMuMu_MassConstraint_displacement = compute3dDisplacement(bToKPsi2SMuMu_MassConstraint, *pvHandle_.product(),true);
    // }
    // addFitInfo(btokmmCand, bToKPsi2SMuMu_MassConstraint, "psimc", bToKPsi2SMuMu_MassConstraint_displacement,-1,-1,1);



}

void BsToMuMuGammaNTuplizer::addMuMuKBranches()
{
    addCompositParticleBranches("muMuK","nMuMuK");
    storageMapIntArray["muMuK_dimuonIdx"]   =  new Int_t[N_COMPOSIT_PART_MAX];
    theTree->Branch("muMuK_dimuonIdx", storageMapIntArray["muMuK_dimuonIdx"],"muMuK_dimuonIdx[nMuMuK]/I");
    storageMapFloatArray["muMuK_kaonPt"]      =  new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch("muMuK_kaonPt", storageMapFloatArray["muMuK_kaonPt"],"muMuK_kaonPt[nMuMuK]/F");
    storageMapFloatArray["muMuK_kaonEta"]      =  new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch("muMuK_kaonEta", storageMapFloatArray["muMuK_kaonEta"],"muMuK_kaonEta[nMuMuK]/F");
    storageMapFloatArray["muMuK_kaonPhi"]      =  new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch("muMuK_kaonPhi", storageMapFloatArray["muMuK_kaonPhi"],"muMuK_kaonPhi[nMuMuK]/F");
    storageMapFloatArray["muMuK_kaonMass"]      =  new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch("muMuK_kaonMass", storageMapFloatArray["muMuK_kaonMass"],"muMuK_kaonMass[nMuMuK]/F");
    storageMapFloatArray["muMuK_kaonCharge"]      =  new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch("muMuK_kaonCharge", storageMapFloatArray["muMuK_kaonCharge"],"muMuK_kaonCharge[nMuMuK]/F");
    storageMapFloatArray["muMuK_kaon_sdxy_bs"]      =  new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch("muMuK_kaon_sdxy_bs", storageMapFloatArray["muMuK_kaon_sdxy_bs"],"muMuK_kaon_sdxy_bs[nMuMuK]/F");
    
    if(isMC)
    {
       storageMapIntArray["gen_mmk_cpdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_cpdgId", storageMapIntArray["gen_mmk_cpdgId"],"gen_mmk_cpdgId[nMuMuK]/I");
       storageMapIntArray["gen_mmk_kaon_mpdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_kaon_mpdgId", storageMapIntArray["gen_mmk_kaon_mpdgId"],"gen_mmk_kaon_mpdgId[nMuMuK]/I");
       storageMapIntArray["gen_mmk_kaon_pdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_kaon_pdgId", storageMapIntArray["gen_mmk_kaon_pdgId"],"gen_mmk_kaon_pdgId[nMuMuK]/I");
       storageMapFloatArray["gen_mmk_kaon_pt"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_kaon_pt", storageMapFloatArray["gen_mmk_kaon_pt"],"gen_mmk_kaon_pt[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_l3d"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_l3d", storageMapFloatArray["gen_mmk_l3d"],"gen_mmk_l3d[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_lxy"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_lxy", storageMapFloatArray["gen_mmk_lxy"],"gen_mmk_lxy[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_mass"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_mass", storageMapFloatArray["gen_mmk_mass"],"gen_mmk_mass[nMuMuK]/F");
       storageMapIntArray["gen_mmk_pdgId"]   = new Int_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_pdgId", storageMapIntArray["gen_mmk_pdgId"],"gen_mmk_pdgId[nMuMuK]/I");
       storageMapFloatArray["gen_mmk_prod_x"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_prod_x", storageMapFloatArray["gen_mmk_prod_x"],"gen_mmk_prod_x[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_prod_y"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_prod_y", storageMapFloatArray["gen_mmk_prod_y"],"gen_mmk_prod_y[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_prod_z"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_prod_z", storageMapFloatArray["gen_mmk_prod_z"],"gen_mmk_prod_z[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_pt"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_pt", storageMapFloatArray["gen_mmk_pt"],"gen_mmk_pt[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_tau"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_tau", storageMapFloatArray["gen_mmk_tau"],"gen_mmk_tau[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_alpha_p_phi"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_alpha_p_phi", storageMapFloatArray["gen_mmk_alpha_p_phi"],"gen_mmk_alpha_p_phi[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_alpha_p_theta"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_alpha_p_theta", storageMapFloatArray["gen_mmk_alpha_p_theta"],"gen_mmk_alpha_p_theta[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_alpha_ip"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_alpha_ip", storageMapFloatArray["gen_mmk_alpha_ip"],"gen_mmk_alpha_ip[nMuMuK]/F");
       storageMapFloatArray["gen_mmk_alpha_vtx"]   = new Float_t[N_COMPOSIT_PART_MAX];
       theTree->Branch("gen_mmk_alpha_vtx", storageMapFloatArray["gen_mmk_alpha_vtx"],"gen_mmk_alpha_vtx[nMuMuK]/F");
    }
    // JPsiK
    if(doJPsiK_)
    {
        addCompositParticleBranches("jpsiK","nJPsiK");
        storageMapIntArray["jpsiK_muMuKIdx"]   =  new Int_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_muMuKIdx", storageMapIntArray["jpsiK_muMuKIdx"],"jpsiK_muMuKIdx[nJPsiK]/I");
        storageMapIntArray["jpsiK_dimuonIdx"]   =  new Int_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_dimuonIdx", storageMapIntArray["jpsiK_dimuonIdx"],"jpsiK_dimuonIdx[nJPsiK]/I");
        storageMapFloatArray["jpsiK_kaonPt"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_kaonPt", storageMapFloatArray["jpsiK_kaonPt"],"jpsiK_kaonPt[nJPsiK]/F");
        storageMapFloatArray["jpsiK_kaonEta"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_kaonEta", storageMapFloatArray["jpsiK_kaonEta"],"jpsiK_kaonEta[nJPsiK]/F");
        storageMapFloatArray["jpsiK_kaonPhi"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_kaonPhi", storageMapFloatArray["jpsiK_kaonPhi"],"jpsiK_kaonPhi[nJPsiK]/F");
        storageMapFloatArray["jpsiK_kaonMass"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_kaonMass", storageMapFloatArray["jpsiK_kaonMass"],"jpsiK_kaonMass[nJPsiK]/F");
        storageMapFloatArray["jpsiK_kaonCharge"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_kaonCharge", storageMapFloatArray["jpsiK_kaonCharge"],"jpsiK_kaonCharge[nJPsiK]/F");
        storageMapFloatArray["jpsiK_kaon_sdxy_bs"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("jpsiK_kaon_sdxy_bs", storageMapFloatArray["jpsiK_kaon_sdxy_bs"],"jpsiK_kaon_sdxy_bs[nJPsiK]/F");

    }

    // Psi2sK
    if(doPsi2SK_)
    {
        addCompositParticleBranches("psi2sK","nPsi2SK");
        storageMapIntArray["psi2sK_muMuKIdx"]   =  new Int_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_muMuKIdx", storageMapIntArray["psi2sK_muMuKIdx"],"psi2sK_muMuKIdx[nPsi2SK]/I");
        storageMapIntArray["psi2sK_dimuonIdx"]   =  new Int_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_dimuonIdx", storageMapIntArray["psi2sK_dimuonIdx"],"psi2sK_dimuonIdx[nPsi2SK]/I");
        storageMapFloatArray["psi2sK_kaonPt"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_kaonPt", storageMapFloatArray["psi2sK_kaonPt"],"psi2sK_kaonPt[nPsi2SK]/F");
        storageMapFloatArray["psi2sK_kaonEta"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_kaonEta", storageMapFloatArray["psi2sK_kaonEta"],"psi2sK_kaonEta[nPsi2SK]/F");
        storageMapFloatArray["psi2sK_kaonPhi"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_kaonPhi", storageMapFloatArray["psi2sK_kaonPhi"],"psi2sK_kaonPhi[nPsi2SK]/F");
        storageMapFloatArray["psi2sK_kaonMass"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_kaonMass", storageMapFloatArray["psi2sK_kaonMass"],"psi2sK_kaonMass[nPsi2SK]/F");
        storageMapFloatArray["psi2sK_kaonCharge"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_kaonCharge", storageMapFloatArray["psi2sK_kaonCharge"],"psi2sK_kaonCharge[nPsi2SK]/F");
        storageMapFloatArray["psi2sK_kaon_sdxy_bs"]      =  new Float_t[N_COMPOSIT_PART_MAX];
        theTree->Branch("psi2sK_kaon_sdxy_bs", storageMapFloatArray["psi2sK_kaon_sdxy_bs"],"psi2sK_kaon_sdxy_bs[nPsi2SK]/F");

    }

}

void BsToMuMuGammaNTuplizer::addCompositParticleBranches(TString tag,TString nTag)
{
    std::string sTag(tag) ;
    std::string sNtag(nTag) ;
    storageMapInt[sNtag]=0;
    theTree->Branch(nTag,   &storageMapInt[sNtag]);
    storageMapFloatArray[sTag+"_valid"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_valid", storageMapFloatArray[sTag+"_valid"],tag+"_valid["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_prob"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_prob", storageMapFloatArray[sTag+"_vtx_prob"],tag+"_vtx_prob["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_chi2dof"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_chi2dof", storageMapFloatArray[sTag+"_vtx_chi2dof"],tag+"_vtx_chi2dof["+nTag+"]/F");
    storageMapFloatArray[sTag+"_mass"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_mass", storageMapFloatArray[sTag+"_mass"],tag+"_mass["+nTag+"]/F");
    storageMapFloatArray[sTag+"_massErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_massErr", storageMapFloatArray[sTag+"_massErr"],tag+"_massErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_lxy"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_lxy", storageMapFloatArray[sTag+"_lxy"],tag+"_lxy["+nTag+"]/F");
    storageMapFloatArray[sTag+"_sigLxy"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_sigLxy", storageMapFloatArray[sTag+"_sigLxy"],tag+"_sigLxy["+nTag+"]/F");
    storageMapFloatArray[sTag+"_alphaBS"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_alphaBS", storageMapFloatArray[sTag+"_alphaBS"],tag+"_alphaBS["+nTag+"]/F");
    storageMapFloatArray[sTag+"_alphaBSErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_alphaBSErr", storageMapFloatArray[sTag+"_alphaBSErr"],tag+"_alphaBSErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_x"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_x", storageMapFloatArray[sTag+"_vtx_x"],tag+"_vtx_x["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_xErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_xErr", storageMapFloatArray[sTag+"_vtx_xErr"],tag+"_vtx_xErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_y"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_y", storageMapFloatArray[sTag+"_vtx_y"],tag+"_vtx_y["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_yErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_yErr", storageMapFloatArray[sTag+"_vtx_yErr"],tag+"_vtx_yErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_z"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_z", storageMapFloatArray[sTag+"_vtx_z"],tag+"_vtx_z["+nTag+"]/F");
    storageMapFloatArray[sTag+"_vtx_zErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_vtx_zErr", storageMapFloatArray[sTag+"_vtx_zErr"],tag+"_vtx_zErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pt"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pt", storageMapFloatArray[sTag+"_pt"],tag+"_pt["+nTag+"]/F");
    storageMapFloatArray[sTag+"_eta"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_eta", storageMapFloatArray[sTag+"_eta"],tag+"_eta["+nTag+"]/F");
    storageMapFloatArray[sTag+"_phi"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_phi", storageMapFloatArray[sTag+"_phi"],tag+"_phi["+nTag+"]/F");
    storageMapFloatArray[sTag+"_alpha"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_alpha", storageMapFloatArray[sTag+"_alpha"],tag+"_alpha["+nTag+"]/F");
    storageMapFloatArray[sTag+"_alphaErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_alphaErr", storageMapFloatArray[sTag+"_alphaErr"],tag+"_alphaErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_l3d"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_l3d", storageMapFloatArray[sTag+"_l3d"],tag+"_l3d["+nTag+"]/F");
    storageMapFloatArray[sTag+"_sl3d"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_sl3d", storageMapFloatArray[sTag+"_sl3d"],tag+"_sl3d["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pv_z"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pv_z", storageMapFloatArray[sTag+"_pv_z"],tag+"_pv_z["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pv_zErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pv_zErr", storageMapFloatArray[sTag+"_pv_zErr"],tag+"_pv_zErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pvip"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pvip", storageMapFloatArray[sTag+"_pvip"],tag+"_pvip["+nTag+"]/F");
    storageMapFloatArray[sTag+"_spvip"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_spvip", storageMapFloatArray[sTag+"_spvip"],tag+"_spvip["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pvipErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pvipErr", storageMapFloatArray[sTag+"_pvipErr"],tag+"_pvipErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pv2ip"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pv2ip", storageMapFloatArray[sTag+"_pv2ip"],tag+"_pv2ip["+nTag+"]/F");
    storageMapFloatArray[sTag+"_spv2ip"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_spv2ip", storageMapFloatArray[sTag+"_spv2ip"],tag+"_spv2ip["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pv2ipErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pv2ipErr", storageMapFloatArray[sTag+"_pv2ipErr"],tag+"_pv2ipErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pvlip"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pvlip", storageMapFloatArray[sTag+"_pvlip"],tag+"_pvlip["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pvlipSig"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pvlipSig", storageMapFloatArray[sTag+"_pvlipSig"],tag+"_pvlipSig["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pvlipErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pvlipErr", storageMapFloatArray[sTag+"_pvlipErr"],tag+"_pvlipErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pv2lip"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pv2lip", storageMapFloatArray[sTag+"_pv2lip"],tag+"_pv2lip["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pv2lipSig"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pv2lipSig", storageMapFloatArray[sTag+"_pv2lipSig"],tag+"_pv2lipSig["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pv2lipErr"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pv2lipErr", storageMapFloatArray[sTag+"_pv2lipErr"],tag+"_pv2lipErr["+nTag+"]/F");
    storageMapFloatArray[sTag+"_pvIndex"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_pvIndex", storageMapFloatArray[sTag+"_pvIndex"],tag+"_pvIndex["+nTag+"]/F");
    storageMapFloatArray[sTag+"_tau"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_tau", storageMapFloatArray[sTag+"_tau"],tag+"_tau["+nTag+"]/F");
    storageMapFloatArray[sTag+"_taue"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_taue", storageMapFloatArray[sTag+"_taue"],tag+"_taue["+nTag+"]/F");
    storageMapFloatArray[sTag+"_tauxy"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_tauxy", storageMapFloatArray[sTag+"_tauxy"],tag+"_tauxy["+nTag+"]/F");
    storageMapFloatArray[sTag+"_tauxye"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_tauxye", storageMapFloatArray[sTag+"_tauxye"],tag+"_tauxye["+nTag+"]/F");
    storageMapFloatArray[sTag+"_nTrks"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_nTrks", storageMapFloatArray[sTag+"_nTrks"],tag+"_nTrks["+nTag+"]/F");
    storageMapFloatArray[sTag+"_nBMTrks"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_nBMTrks", storageMapFloatArray[sTag+"_nBMTrks"],tag+"_nBMTrks["+nTag+"]/F");
    storageMapFloatArray[sTag+"_nDisTrks"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_nDisTrks", storageMapFloatArray[sTag+"_nDisTrks"],tag+"_nDisTrks["+nTag+"]/F");
    storageMapFloatArray[sTag+"_closetrk"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_closetrk", storageMapFloatArray[sTag+"_closetrk"],tag+"_closetrk["+nTag+"]/F");
    storageMapFloatArray[sTag+"_closetrks1"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_closetrks1", storageMapFloatArray[sTag+"_closetrks1"],tag+"_closetrks1["+nTag+"]/F");
    storageMapFloatArray[sTag+"_closetrks2"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_closetrks2", storageMapFloatArray[sTag+"_closetrks2"],tag+"_closetrks2["+nTag+"]/F");
    storageMapFloatArray[sTag+"_closetrks3"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_closetrks3", storageMapFloatArray[sTag+"_closetrks3"],tag+"_closetrks3["+nTag+"]/F");
    storageMapFloatArray[sTag+"_docatrk"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_docatrk", storageMapFloatArray[sTag+"_docatrk"],tag+"_docatrk["+nTag+"]/F");
    storageMapFloatArray[sTag+"_m1iso"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_m1iso", storageMapFloatArray[sTag+"_m1iso"],tag+"_m1iso["+nTag+"]/F");
    storageMapFloatArray[sTag+"_m2iso"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_m2iso", storageMapFloatArray[sTag+"_m2iso"],tag+"_m2iso["+nTag+"]/F");
    storageMapFloatArray[sTag+"_iso"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_iso", storageMapFloatArray[sTag+"_iso"],tag+"_iso["+nTag+"]/F");
    storageMapFloatArray[sTag+"_otherVtxMaxProb"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_otherVtxMaxProb", storageMapFloatArray[sTag+"_otherVtxMaxProb"],tag+"_otherVtxMaxProb["+nTag+"]/F");
    storageMapFloatArray[sTag+"_otherVtxMaxProb1"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_otherVtxMaxProb1", storageMapFloatArray[sTag+"_otherVtxMaxProb1"],tag+"_otherVtxMaxProb1["+nTag+"]/F");
    storageMapFloatArray[sTag+"_otherVtxMaxProb2"]   = new Float_t[N_COMPOSIT_PART_MAX];
    theTree->Branch(tag+"_otherVtxMaxProb2", storageMapFloatArray[sTag+"_otherVtxMaxProb2"],tag+"_otherVtxMaxProb2["+nTag+"]/F");
}

void BsToMuMuGammaNTuplizer::fillCompositParticleBranches( std::string tag, Int_t idx ,const KinematicFitResult &fit ,
                const DisplacementInformationIn3D &displacement ,CloseTrackInfo &closeTracks)
{
    storageMapFloatArray[tag+"_valid"][idx]  =        fit.valid() ;
    storageMapFloatArray[tag+"_vtx_prob"][idx]  =     fit.vtxProb() ;
    storageMapFloatArray[tag+"_vtx_chi2dof"][idx]  =  fit.chi2()>0?fit.chi2()/fit.ndof():-1;
    storageMapFloatArray[tag+"_mass"][idx]  =         fit.mass() ;
    storageMapFloatArray[tag+"_massErr"][idx]  =      fit.massErr() ;
    storageMapFloatArray[tag+"_lxy"][idx]  =          fit.lxy ;
    storageMapFloatArray[tag+"_sigLxy"][idx]  =       fit.sigLxy ;
    storageMapFloatArray[tag+"_alphaBS"][idx]  =      fit.alphaBS;
    storageMapFloatArray[tag+"_alphaBSErr"][idx]  =   fit.alphaBSErr;
    storageMapFloatArray[tag+"_vtx_x"][idx]  =        fit.valid()?fit.refitVertex->position().x():0 ;
    storageMapFloatArray[tag+"_vtx_xErr"][idx]  =     fit.valid()?sqrt(fit.refitVertex->error().cxx()):0 ;
    storageMapFloatArray[tag+"_vtx_y"][idx]  =        fit.valid()?fit.refitVertex->position().y():0 ;
    storageMapFloatArray[tag+"_vtx_yErr"][idx]  =     fit.valid()?sqrt(fit.refitVertex->error().cyy()):0 ;
    storageMapFloatArray[tag+"_vtx_z"][idx]  =        fit.valid()?fit.refitVertex->position().z():0 ;
    storageMapFloatArray[tag+"_vtx_zErr"][idx]  =     fit.valid()?sqrt(fit.refitVertex->error().czz()):0 ;
    storageMapFloatArray[tag+"_pt"][idx]  =           fit.p3().perp() ;
    storageMapFloatArray[tag+"_eta"][idx]  =          fit.p3().eta() ;
    storageMapFloatArray[tag+"_phi"][idx]  =          fit.p3().phi() ;
    
    // IP info
    storageMapFloatArray[tag+"_alpha"][idx]  =        displacement.alpha;
    storageMapFloatArray[tag+"_alphaErr"][idx]  =     displacement.alphaErr;
    storageMapFloatArray[tag+"_l3d"][idx]  =          displacement.decayLength;
    storageMapFloatArray[tag+"_sl3d"][idx]  =         displacement.decayLengthErr>0?displacement.decayLength/displacement.decayLengthErr:0;
    storageMapFloatArray[tag+"_pv_z"][idx]  =         displacement.pv?displacement.pv->position().z():0;
    storageMapFloatArray[tag+"_pv_zErr"][idx]  =      displacement.pv?displacement.pv->zError():0;
    storageMapFloatArray[tag+"_pvip"][idx]  =         displacement.distaceOfClosestApproach;
    storageMapFloatArray[tag+"_spvip"][idx]  =        displacement.distaceOfClosestApproachSig;
    storageMapFloatArray[tag+"_pvipErr"][idx]  =      displacement.distaceOfClosestApproachErr;
    storageMapFloatArray[tag+"_pv2ip"][idx]  =        displacement.distaceOfClosestApproach2;
    storageMapFloatArray[tag+"_spv2ip"][idx]  =       displacement.distaceOfClosestApproach2Sig;
    storageMapFloatArray[tag+"_pv2ipErr"][idx]  =     displacement.distaceOfClosestApproach2Err;
    storageMapFloatArray[tag+"_pvlip"][idx]  =        displacement.longitudinalImpactParameter;
    storageMapFloatArray[tag+"_pvlipSig"][idx]  =     displacement.longitudinalImpactParameterSig;
    storageMapFloatArray[tag+"_pvlipErr"][idx]  =     displacement.longitudinalImpactParameterErr;
    storageMapFloatArray[tag+"_pv2lip"][idx]  =       displacement.longitudinalImpactParameter2;
    storageMapFloatArray[tag+"_pv2lipSig"][idx]  =    displacement.longitudinalImpactParameter2Sig;
    storageMapFloatArray[tag+"_pv2lipErr"][idx]  =    displacement.longitudinalImpactParameter2Err;
    storageMapFloatArray[tag+"_pvIndex"][idx]  =      displacement.pvIndex;

    // DecayTime
    storageMapFloatArray[tag+"_tau"][idx]  =          displacement.decayTime;
    storageMapFloatArray[tag+"_taue"][idx]  =         displacement.decayTimeError;
    storageMapFloatArray[tag+"_tauxy"][idx]  =        displacement.decayTimeXY;
    storageMapFloatArray[tag+"_tauxye"][idx]  =       displacement.decayTimeXYError;
 
    // Close Tracks
    auto pvIndex = displacement.pvIndex;
    const reco::Vertex *vertex(nullptr) ;
    if(pvIndex >=0 ){ 
        vertex=&(pvHandle_->at(pvIndex)); 
       }
    storageMapFloatArray[tag+"_nTrks"][idx]       =   closeTracks.nTracksByVertexProbability(0.1, -1.0, vertex);
    storageMapFloatArray[tag+"_nBMTrks"][idx]     =   closeTracks.nTracksByBetterMatch();
    storageMapFloatArray[tag+"_nDisTrks"][idx]    =   closeTracks.nTracksByVertexProbability(0.1,  2.0 ,  vertex);
    storageMapFloatArray[tag+"_closetrk"][idx]    =   closeTracks.nTracksByDisplacementSignificance(0.03 ,-1 ,   vertex);
    storageMapFloatArray[tag+"_closetrks1"][idx]  =   closeTracks.nTracksByDisplacementSignificance(0.03 , 1 ,   vertex);
    storageMapFloatArray[tag+"_closetrks2"][idx]  =   closeTracks.nTracksByDisplacementSignificance(0.03 , 2 ,   vertex);
    storageMapFloatArray[tag+"_closetrks3"][idx]  =   closeTracks.nTracksByDisplacementSignificance(0.03 , 3 ,   vertex);
    storageMapFloatArray[tag+"_docatrk"][idx]     =   closeTracks.minDoca(0.03,vertex);

   
}


void BsToMuMuGammaNTuplizer::addHLTObjectBranches()
{
      storageMapInt["nHLTObj"]=0;
      theTree->Branch("nHLTObj", &storageMapInt["nHLTObj"]);
      storageMapIntArray["hltObj_trigIdx"] = new Int_t[N_L3MUON];
      theTree->Branch("hltObj_trigIdx",   storageMapIntArray["hltObj_trigIdx"],"hltObj_trigIdx[nHLTObj]/I");
      storageMapIntArray["hltObj_idx"] = new Int_t[N_L3MUON];
      theTree->Branch("hltObj_idx",   storageMapIntArray["hltObj_idx"],"hltObj_idx[nHLTObj]/I");
      storageMapFloatArray["hltObj_pt"] = new Float_t[N_L3MUON];
      theTree->Branch("hltObj_pt",   storageMapFloatArray["hltObj_pt"],"hltObj_pt[nHLTObj]/F");
      storageMapFloatArray["hltObj_eta"] = new Float_t[N_L3MUON];
      theTree->Branch("hltObj_eta",   storageMapFloatArray["hltObj_eta"],"hltObj_eta[nHLTObj]/F");
      storageMapFloatArray["hltObj_phi"] = new Float_t[N_L3MUON];
      theTree->Branch("hltObj_phi",   storageMapFloatArray["hltObj_phi"],"hltObj_phi[nHLTObj]/F");
      storageMapFloatArray["hltObj_mass"] = new Float_t[N_L3MUON];
      theTree->Branch("hltObj_mass",   storageMapFloatArray["hltObj_mass"],"hltObj_mass[nHLTObj]/F");
}

void BsToMuMuGammaNTuplizer::fillHLTL3MuonBranches(  std::vector<Int_t> trigIdx, std::vector<Int_t> &Keys,const trigger::TriggerEvent & triggerEvent )
{
    int idx=0;
    storageMapInt["nHLTObj"]=0;
    for( ; idx <int(trigIdx.size()) ; )
    {  
         auto objectKey = Keys[idx];
         storageMapIntArray["hltObj_trigIdx"][idx] = trigIdx[idx];
         storageMapIntArray["hltObj_idx"][idx]     = triggerEvent.getObjects()[objectKey].id(); 
         storageMapFloatArray["hltObj_pt"][idx]    = triggerEvent.getObjects()[objectKey].pt(); 
         storageMapFloatArray["hltObj_eta"][idx]   = triggerEvent.getObjects()[objectKey].eta(); 
         storageMapFloatArray["hltObj_phi"][idx]   = triggerEvent.getObjects()[objectKey].phi(); 
         storageMapFloatArray["hltObj_mass"][idx]  = triggerEvent.getObjects()[objectKey].mass(); 
         idx++;
    }
     storageMapInt["nHLTObj"]=idx;
}


///  GEN MATCH Stuff

const reco::Candidate* BsToMuMuGammaNTuplizer::getGenParticle(const reco::Candidate* cand)
{
  if (not cand) return nullptr;
  auto muon = dynamic_cast<const pat::Muon*>(cand);
  if (muon and muon->genParticle()) return muon->genParticle();
  for (auto const & genParticle: *genParticleCollection){
      if (dr_match(cand->p4(), genParticle.p4()))
	return &genParticle;
  }
  return nullptr;
}

GenMatchInfo BsToMuMuGammaNTuplizer::getGenMatchInfo( const reco::Muon& muon1,
						const reco::Muon& muon2,
						const reco::PFCandidate* kaon1,
						const reco::PFCandidate* kaon2,
						const reco::PFCandidate* photon)
{
  auto result = GenMatchInfo();
  const reco::Candidate*   mm_mother(0);
  //assert(prunedGenParticles_);
  //assert(genParticleCollection);
  std::vector<const reco::Candidate*> daughters;

  result.mc_mu1 = getGenParticle(&muon1);
  if (result.mc_mu1){
    result.mu1_pdgId = result.mc_mu1->pdgId();
    result.mu1_pt    = result.mc_mu1->pt();
    if (result.mc_mu1->mother()){
      result.mu1_motherPdgId = result.mc_mu1->mother()->pdgId();
    }
    daughters.push_back(result.mc_mu1);
  }

  result.mc_mu2 = getGenParticle(&muon2);
  if (result.mc_mu2){
    result.mu2_pdgId = result.mc_mu2->pdgId();
    result.mu2_pt    = result.mc_mu2->pt();
    if (result.mc_mu2->mother()){
      result.mu2_motherPdgId = result.mc_mu2->mother()->pdgId();
    }
    daughters.push_back(result.mc_mu2);
  }

  if ( result.mc_mu1 and result.mc_mu2 ){
    if ( (result.mc_mu1->vertex()-result.mc_mu2->vertex()).r() < 1e-4)
      result.mm_vtx    = result.mc_mu1->vertex();
    if ( result.mc_mu1->mother() and result.mc_mu1->mother() == result.mc_mu2->mother() ){
      mm_mother = result.mc_mu1->mother();
      result.match = result.mc_mu1->mother();
      result.mm_mass      = mm_mother->mass();
      result.mm_pt        = mm_mother->pt();
      result.mm_pdgId     = mm_mother->pdgId();
      if (mm_mother->mother()) result.mm_motherPdgId = mm_mother->mother()->pdgId();
      result.mm_prod_vtx = getProductionVertex(mm_mother);
    }
  }
  
  if (kaon1){
    for (auto const & genParticle: *genParticleCollection){
      if (dr_match(kaon1->p4(),genParticle.p4())){
	result.mc_kaon1 = &genParticle;
	daughters.push_back(result.mc_kaon1);
	result.kaon1_pdgId = genParticle.pdgId();
	result.kaon1_pt    = genParticle.pt();
	if (genParticle.mother(0)){
	  result.kaon1_motherPdgId = genParticle.mother(0)->pdgId();
	}
	break;
      }
    }
    if (daughters.size()==3){
      const auto* mother = find_common_ancestor(daughters);
      if (mother){
	result.match        = mother;
	result.kmm_pdgId    = mother->pdgId();
	result.kmm_mass     = mother->mass();
	result.kmm_pt       = mother->pt();
	result.kmm_prod_vtx = getProductionVertex(mother);
      }
    }
  }
  if (kaon2){
    for (auto const & genParticle: *genParticleCollection){
      if (dr_match(kaon2->p4(),genParticle.p4())){
	result.mc_kaon2 = &genParticle;
	daughters.push_back(result.mc_kaon2);
	result.kaon2_pdgId = genParticle.pdgId();
	result.kaon2_pt    = genParticle.pt();
	if (genParticle.mother(0)){
	  result.kaon2_motherPdgId = genParticle.mother(0)->pdgId();
	}
	break;
      }
    }
    if (daughters.size()==4){
      const auto* mother = find_common_ancestor(daughters);
      if (mother){
	result.match         = mother;
	result.kkmm_pdgId    = mother->pdgId();
	result.kkmm_mass     = mother->mass();
	result.kkmm_pt       = mother->pt();
	result.kkmm_prod_vtx = getProductionVertex(mother);
      }
    }
  }
  
  if (photon){
    for (auto const & genParticle: *genParticleCollection){
      if (dr_match(photon->p4(),genParticle.p4())){
	result.mc_photon = &genParticle;
	daughters.push_back(result.mc_photon);
	result.photon_pdgId = genParticle.pdgId();
	result.photon_pt    = genParticle.pt();
	if (genParticle.mother(0)){
	  result.photon_motherPdgId = genParticle.mother(0)->pdgId();
	}
	break;
      }
    }
    if (daughters.size()==3){
      const auto* mother = find_common_ancestor(daughters);
      if (mother){
	result.match        = mother;
	result.mmg_pdgId    = mother->pdgId();
	result.mmg_mass     = mother->mass();
	result.mmg_pt       = mother->pt();
	result.mmg_prod_vtx = getProductionVertex(mother);
      }
    }
  }

  if (daughters.size() > 1){
    const auto* mother = find_common_ancestor(daughters); 
    if (mother){ 
      result.common_mother = mother;
    }
  }

  return result;
}

float BsToMuMuGammaNTuplizer::distanceOfClosestApproach( const reco::GenParticle* track1,
						   const reco::GenParticle* track2)
{
  TwoTrackMinimumDistance md;
  GlobalPoint trk1_pos(track1->vertex().x(), 
		       track1->vertex().y(), 
		       track1->vertex().z());
  GlobalVector trk1_mom(track1->px(),track1->py(),track1->pz());

  GlobalTrajectoryParameters trk1(trk1_pos,trk1_mom,track1->charge(),bFieldHandle_.product());
  GlobalPoint trk2_pos(track2->vertex().x(), 
		       track2->vertex().y(), 
		       track2->vertex().z());
  GlobalVector trk2_mom(track2->px(),track2->py(),track2->pz());

  GlobalTrajectoryParameters trk2(trk2_pos,trk2_mom,track2->charge(),bFieldHandle_.product());
  if ( not md.calculate( trk1, trk2 ) ) return -1.0;
  return md.distance();
}




#include "BMMGUtil.h"


//define this as a plug-in
DEFINE_FWK_MODULE(BsToMuMuGammaNTuplizer);
