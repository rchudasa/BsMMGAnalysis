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
#include <TTree.h>


BsToMuMuGammaNTuplizer::BsToMuMuGammaNTuplizer(const edm::ParameterSet& iConfig) :

  doGenParticles_(iConfig.getParameter<bool>("doGenParticles")),
  doMuons_(iConfig.getParameter<bool>("doMuons")),
  doPhotons_(iConfig.getParameter<bool>("doPhotons")),
  doPFPhotons_(iConfig.getParameter<bool>("doPFPhotons")),
  doSuperClusters_(iConfig.getParameter<bool>("doSuperClusters"))

  //genParticleSrc_(iConfig.getUntrackedParameter<edm::InputTag>("genParticleSrc")),
  //gedPhotonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("gedPhotonSrc")),
  //pfPhotonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("pfPhotonSrc")),
  //MustacheSCBarrelSrc_(iConfig.getParameter<edm::InputTag>("MustacheSCBarrelSrc")),
  //MustacheSCEndcapSrc_(iConfig.getParameter<edm::InputTag>("MustacheSCEndcapSrc"))
{
   if(doMuons_) muonToken_              = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  //now do what ever initialization is needed
  // if(doGenParticles_){
//	 genParticlesCollection_   = consumes<edm::View<reco::GenParticle>>(iConfig.getUntrackedParameter<edm::InputTag>("genParticleSrc"));
 //  }
   if(doPhotons_)    gedPhotonsCollection_       = consumes<std::vector<reco::Photon>>(iConfig.getUntrackedParameter<edm::InputTag>("gedPhotonSrc"));
   if(doPFPhotons_)  pfPhotonsCollection_        = consumes<std::vector<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("pfPhotonSrc"));
   if(doSuperClusters_){
	MustacheSCBarrelCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCBarrelSrc"));
   	MustacheSCEndcapCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCEndcapSrc"));
   }
  
  beamSpotToken_  	 =consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  primaryVtxToken_       =consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  if(doGenParticles_)simGenTocken_          =consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));

  etaMax_muon               =  iConfig.getUntrackedParameter<double>("muon_EtaMax")        ;
  dcaMax_muon_bs            =  iConfig.getUntrackedParameter<double>("muon_dcaMAX")        ;
  pTMinMuons 		      =  iConfig.getUntrackedParameter<double>("muon_minPt");
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

  
  // initialize output TTree
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("EventTree", "Event data");
  
  theTree->Branch("run",    &run_);
  theTree->Branch("event",  &event_);
  theTree->Branch("lumis",  &lumis_);
  theTree->Branch("isData", &isData_);
  
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
 
  theTree->Branch("gen_hasAValid_candidate"	,&gen_hasAValid_candidate_);
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
  theTree->Branch("mumPx",            &mumPx_);
  theTree->Branch("mumPy",            &mumPy_);
  theTree->Branch("mumPz",            &mumPz_);
  theTree->Branch("mumDCAVtx",        &mumDCAVtx_);
  theTree->Branch("mumDCAVtxE",       &mumDCAVtxE_);
  theTree->Branch("mumDCABS",         &mumDCABS_);
  theTree->Branch("mumDCABSE",        &mumDCABSE_);
  theTree->Branch("mumKinkChi2",      &mumKinkChi2_);
  theTree->Branch("mumFracHits",      &mumFracHits_);
  theTree->Branch("mumdxyBS",         &mumdxyBS_);
  theTree->Branch("mumdzBS",          &mumdzBS_);
  theTree->Branch("mumMinIP2D",       &mumMinIP2D_);
  theTree->Branch("mumMinIP2DE",      &mumMinIP2DE_);
  theTree->Branch("mumMinIP",         &mumMinIP_);
  theTree->Branch("mumMinIPS",        &mumMinIPS_);
  theTree->Branch("mumDeltaRwithMC",  &mumDeltaRwithMC_);
  theTree->Branch("mumCat",           &mumCat_);
  theTree->Branch("mumNPixHits",      &mumNPixHits_);
  theTree->Branch("mumNPixLayers",    &mumNPixLayers_);
  theTree->Branch("mumNTrkHits",      &mumNTrkHits_);
  theTree->Branch("mumNTrkLayers",    &mumNTrkLayers_);
  theTree->Branch("mumNMuonHits",     &mumNMuonHits_);
  theTree->Branch("mumNMatchStation", &mumNMatchStation_);
  theTree->Branch("mumIso",           &mumIso_);
  theTree->Branch("mumIsoPt",         &mumIsoPt_);
  theTree->Branch("mumIsodR",         &mumIsodR_);
  theTree->Branch("mum_isGlobalMuon", &mum_isGlobalMuon_);
  theTree->Branch("mum_isTrackerMuon",&mum_isTrackerMuon_);
  theTree->Branch("mum_StandAloneMuon",&mum_StandAloneMuon_);
  theTree->Branch("mum_isCaloMuon",    &mum_isCaloMuon_);
  theTree->Branch("mum_isPFMuon",      &mum_isPFMuon_);

  // ### mu+ ###  
  theTree->Branch("nMuP",             &nMuP_);
  theTree->Branch("mupHighPurity",    &mupHighPurity_);
  theTree->Branch("mupPt",            &mupPt_);
  theTree->Branch("mupEta",           &mupEta_);
  theTree->Branch("mupPhi",           &mupPhi_);
  theTree->Branch("mupCL",            &mupCL_);
  theTree->Branch("mupNormChi2",      &mupNormChi2_);
  theTree->Branch("mupPx",            &mupPx_);
  theTree->Branch("mupPy",            &mupPy_);
  theTree->Branch("mupPz",            &mupPz_);
  theTree->Branch("mupDCAVtx",        &mupDCAVtx_);
  theTree->Branch("mupDCAVtxE",       &mupDCAVtxE_);
  theTree->Branch("mupDCABS",         &mupDCABS_);
  theTree->Branch("mupDCABSE",        &mupDCABSE_);
  theTree->Branch("mupKinkChi2",      &mupKinkChi2_);
  theTree->Branch("mupFracHits",      &mupFracHits_);
  theTree->Branch("mupdxyBS",         &mupdxyBS_);
  theTree->Branch("mupdzBS",          &mupdzBS_);
  theTree->Branch("mupMinIP2D",       &mupMinIP2D_);
  theTree->Branch("mupMinIP2DE",      &mupMinIP2DE_);
  theTree->Branch("mupMinIP",         &mupMinIP_);
  theTree->Branch("mupMinIPS",        &mupMinIPS_);
  theTree->Branch("mupDeltaRwithMC",  &mupDeltaRwithMC_);
  theTree->Branch("mupCat",           &mupCat_);
  theTree->Branch("mupNPixHits",      &mupNPixHits_);
  theTree->Branch("mupNPixLayers",    &mupNPixLayers_);
  theTree->Branch("mupNTrkHits",      &mupNTrkHits_);
  theTree->Branch("mupNTrkLayers",    &mupNTrkLayers_);
  theTree->Branch("mupNMuonHits",     &mupNMuonHits_);
  theTree->Branch("mupNMatchStation", &mupNMatchStation_);
  theTree->Branch("mupIso",           &mupIso_);
  theTree->Branch("mupIsoPt",         &mupIsoPt_);
  theTree->Branch("mupIsodR",         &mupIsodR_);
  theTree->Branch("mup_isGlobalMuon", &mup_isGlobalMuon_);
  theTree->Branch("mup_isTrackerMuon",&mup_isTrackerMuon_);
  theTree->Branch("mup_StandAloneMuon",&mup_StandAloneMuon_);
  theTree->Branch("mup_isCaloMuon",    &mup_isCaloMuon_);
  theTree->Branch("mup_isPFMuon",      &mup_isPFMuon_);

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
  theTree->Branch("mumuVtxProb",     &mumuVtxProb_);
  theTree->Branch("mumuVtxIsGoodFit",&mumuVtxIsGoodFit_);
  theTree->Branch("mumuCosAlphaBS",  &mumuCosAlphaBS_);
  theTree->Branch("mumuCosAlphaBSE", &mumuCosAlphaBSE_);
  theTree->Branch("mumuLBS",         &mumuLBS_);
  theTree->Branch("mumuLBSE",        &mumuLBSE_);
  theTree->Branch("mumuDCA",         &mumuDCA_);
  theTree->Branch("mumuLS",          &mumuLS_);
  theTree->Branch("mumuLSErr",       &mumuLSErr_);

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
  gen_hasAValid_candidate_.clear();
  }

  if(doMuons_){
  
   mumHighPurity_.clear();
   mumPt_.clear();
   mumEta_.clear();
   mumPhi_.clear();
   mumCL_.clear(); 
   mumNormChi2_.clear();
   mumPx_.clear();
   mumPy_.clear();
   mumPz_.clear();
   mumDCAVtx_.clear();
   mumDCAVtxE_.clear();
   mumDCABS_.clear();
   mumDCABSE_.clear();
   mumKinkChi2_.clear();
   mumFracHits_.clear();
   mumdxyBS_.clear();
   mumdzBS_.clear();
   mumMinIP2D_.clear();
   mumMinIP2DE_.clear();
   mumMinIP_.clear();
   mumMinIPS_.clear();
   mumDeltaRwithMC_.clear();
   mumCat_.clear();
   mumCharge_.clear();
   mumNPixHits_.clear();
   mumNPixLayers_.clear();
   mumNTrkHits_.clear();
   mumNTrkLayers_.clear();
   mumNMuonHits_.clear();
   mumNMatchStation_.clear();
   mumIso_.clear();
   mumIsoPt_.clear();
   mumIsodR_.clear();
   mum_isGlobalMuon_.clear();
   mum_isTrackerMuon_.clear();
   mum_StandAloneMuon_.clear();
   mum_isCaloMuon_.clear();
   mum_isPFMuon_.clear();
  
  
  
   mupHighPurity_.clear();
   mupPt_.clear();
   mupEta_.clear();
   mupPhi_.clear();
   mupCL_.clear(); 
   mupNormChi2_.clear();
   mupPx_.clear();
   mupPy_.clear();
   mupPz_.clear();
   mupDCAVtx_.clear();
   mupDCAVtxE_.clear();
   mupDCABS_.clear();
   mupDCABSE_.clear();
   mupKinkChi2_.clear();
   mupFracHits_.clear();
   mupdxyBS_.clear();
   mupdzBS_.clear();
   mupMinIP2D_.clear();
   mupMinIP2DE_.clear();
   mupMinIP_.clear();
   mupMinIPS_.clear();
   mupDeltaRwithMC_.clear();
   mupCat_.clear();
   mupCharge_.clear();
   mupNPixHits_.clear();
   mupNPixLayers_.clear();
   mupNTrkHits_.clear();
   mupNTrkLayers_.clear();
   mupNMuonHits_.clear();
   mupNMatchStation_.clear();
   mupIso_.clear();
   mupIsoPt_.clear();
   mupIsodR_.clear();
   mup_isGlobalMuon_.clear();
   mup_isTrackerMuon_.clear();
   mup_StandAloneMuon_.clear();
   mup_isCaloMuon_.clear();
   mup_isPFMuon_.clear();

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
    mumuVtxProb_.clear();
    mumuVtxIsGoodFit_.clear();
    mumuCosAlphaBS_.clear();
    mumuCosAlphaBSE_.clear(); 
    mumuLBS_.clear();
    mumuLBSE_.clear();
    mumuDCA_.clear();
    mumuLS_.clear();
    mumuLSErr_.clear();
    
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

  if (doMuons_)     fillMuons(iEvent, iSetup);
  if (doPhotons_)    fillPhotons(iEvent, iSetup);
  if (doPFPhotons_) fillPFPhotons(iEvent, iSetup);
  if (doSuperClusters_) fillSC(iEvent);
  theTree->Fill();
  
}



void BsToMuMuGammaNTuplizer::fillGenParticles(const edm::Event& iEvent)
{
  
  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  iEvent.getByToken(simGenTocken_, genParticleCollection);
  
  int phoMul(-1),muMMul(-1),muPMul(-1);
  
  
  for(auto& aBsMeson : *genParticleCollection){
    
    if(abs(aBsMeson.pdgId())!=531) continue;
       std::cout << "event:" << iEvent.id().event() << "  Bs found:" << aBsMeson.pdgId() << std::endl;
    /* for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++)
       {
       auto& bsDaughter = *(aBsMeson.daughter(j));
       
       if(bsDaughter.pdgId() == -13) muMMul++;
       if(bsDaughter.pdgId() ==  13) muPMul++;
       if(bsDaughter.pdgId() ==  22) phoMul++;
       }
       if(muMMul!=1 or muPMul!=1)
       {
       muMMul=0;
       muPMul=0;
       phoMul=0;
       continue;
       }
       
       if(doBsToMuMuGamma and phoMul!=1)
       {
       muMMul=0;
       muPMul=0;
       phoMul=0;
       continue;
       }*/
    
    for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
      
      auto& bsDaughter = *(aBsMeson.daughter(j));
      if(bsDaughter.pdgId() == -13) {
	muMMul = j;
	gen_BsMuonM_pt_.push_back(bsDaughter.pt());
	gen_BsMuonM_eta_.push_back(bsDaughter.eta());
	gen_BsMuonM_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonM_++;
      }
      if(bsDaughter.pdgId() ==  13){
	muPMul = j;
	gen_BsMuonP_pt_.push_back(bsDaughter.pt());
	gen_BsMuonP_eta_.push_back(bsDaughter.eta());
	gen_BsMuonP_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonP_++;
      }
      if(bsDaughter.pdgId() ==  22){
	phoMul = j;
	gen_BsPhoton_pt_.push_back(bsDaughter.pt());
	gen_BsPhoton_eta_.push_back(bsDaughter.eta());
	gen_BsPhoton_phi_.push_back(bsDaughter.phi());
	gen_nBsPhoton_++;
	
      }
    }  //number of daughters
    
    if(phoMul == -1) continue;
    const reco::Candidate *mum = NULL;
    const reco::Candidate *mup = NULL;
    
    if (muMMul != -1 && muPMul != -1) {
      mum = aBsMeson.daughter(muMMul);
      mup = aBsMeson.daughter(muPMul);
    }
    
    if ( mum == NULL || mup == NULL) continue;

    std::cout << "event:" << iEvent.id().event() << "  Bs ID:" << aBsMeson.pdgId() << "  photon multiplicity:" << phoMul <<  "  negative muon number:" << muMMul << "   positive: " << muPMul << std::endl;

    gen_Bs_pt_.push_back(aBsMeson.pt());
    gen_Bs_eta_.push_back(aBsMeson.eta());
    gen_Bs_phi_.push_back(aBsMeson.phi());
    gen_Bs_pz_.push_back(aBsMeson.pz());
    gen_Bs_pdgId_.push_back(aBsMeson.pdgId());
    gen_nBs_++;
    
    // break;
    
  } // genparticle collection
  // gen_BsMuonMMultiplicity_.push_back(muMMul);
  // gen_BsMuonPMultiplicity_.push_back(muPMul);
 // gen_BsPhotonMultiplicity_.push_back(phoMul);
  
/*  bool hasAValidMCCandidate= true;
  if(muMMul!=1 or muPMul!=1 ) 
  {
  if(printMsg) std::cout<<" Ghost event found !! Mu+ Mu- from any of the Bs not found to == 1 "<<std::endl;
  hasAValidMCCandidate = false;
  }
  else if(doBsToMuMuGamma and phoMul!=1)
  {
  if(printMsg) std::cout<<" Ghost event found !! gamma multiplicity from any of the Bs not found to == 1 "<<std::endl;
  hasAValidMCCandidate = false;
  }
  
  gen_hasAValid_candidate_.push_back(hasAValidMCCandidate);*/
  
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
  
  reco::TrackRef muTrackm,muTrackp;
  reco::TransientTrack refitMupTT, refitMumTT;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  int num_SC = 0; int nMuon_pos = 0; int nMuon_neg  = 0; 
  
  // variables
  Utils* Utility;
  Utility= new Utils();
  float muonMass,muonMassErr;
  muonMass= Utility->muonMass;
  muonMassErr= Utility->muonMassErr;
  double chi2,ndof;

  
  //start loop for positive muons 
  for (uint32_t i=0;i<muons->size();i++){
    
    auto &mum=muons->at(i);
    if(mum.pt()  < pTMinMuons) continue;
    muTrackm= mum.innerTrack();
    if ((muTrackm.isNull() == true) || (muTrackm->charge() != -1))   continue;
    
    const reco::TransientTrack muTrackmTT( muTrackm, &(*bFieldHandle));
    if (!muTrackmTT.isValid()) continue;
    
    // # Compute mu- DCA to BeamSpot #
    theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
    double DCAmumBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
    double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
    
    //start loop for positive muons 
    for ( uint32_t j=0;j<muons->size();j++){
      auto &mup=muons->at(j);
      muTrackp = mup.innerTrack();
      if ((muTrackp.isNull() == true) || (muTrackp->charge() != 1)) continue;
      
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));
      if (!muTrackpTT.isValid()) continue;
      
      // # Compute mu+ DCA to BeamSpot #
      theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
      double DCAmupBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
      double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
      
      
      // # Check goodness of muons closest approach #
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
      mumuMass_.push_back(mu_mu_mass);
      mumuMassE_.push_back(mu_mu_mass_err);
      
      mumuPx_.push_back(mumu_KP->currentState().globalMomentum().x());
      mumuPy_.push_back(mumu_KP->currentState().globalMomentum().y());
      mumuPz_.push_back(mumu_KP->currentState().globalMomentum().z());
      
      mumuVtxCL_.push_back(mu_mu_vtx_cl);
      mumuVtxX_.push_back(mumu_KV->position().x());
      mumuVtxY_.push_back(mumu_KV->position().y());
      mumuVtxZ_.push_back(mumu_KV->position().z());
      
      mumuCosAlphaBS_.push_back(MuMuCosAlphaBS);
      mumuCosAlphaBSE_.push_back(MuMuCosAlphaBSErr);
      mumuLBS_.push_back(MuMuLSBS);
      mumuLBSE_.push_back(MuMuLSBSErr);
      mumuDCA_.push_back(mumuDCA);
      
      std::cout << "Entry:" <<  iEvent.id().event() <<  " Neg muon " << i <<  " " << muTrackm->charge() << " pt: " << muTrackm->pt() << " eta:" << muTrackm->eta() << std::endl;
      std::cout << "Entry:" <<  iEvent.id().event() <<  " Pos muon " << j <<  " " < muTrackp->charge() << " pt: " << muTrackp->pt() << " eta:" << muTrackp->eta() << std::endl;

      //// #############
      //// # Save: mu- #
      //// #############
      mumHighPurity_.push_back( (int)muTrackm->quality(reco::Track::highPurity));
      mumPt_.push_back(muTrackm->pt());
      mumEta_.push_back(muTrackm->eta());
      mumPhi_.push_back(muTrackm->phi());
      mumCL_.push_back(TMath::Prob(muTrackmTT.chi2(), static_cast<int>(rint(muTrackmTT.ndof()))));
      mumNormChi2_.push_back(muTrackm->normalizedChi2());
      mumPx_.push_back(refitMumTT.track().momentum().x());
      mumPy_.push_back(refitMumTT.track().momentum().y());
      mumPz_.push_back(refitMumTT.track().momentum().z());
      
      //                     mumDCAVtx_.push_back(DCAmumVtx);
      //                     mumDCAVtxE_.push_back(DCAmumVtxErr);
      mumDCABS_.push_back(DCAmumBS);
      mumDCABSE_.push_back(DCAmumBSErr);
      
      //                     mumKinkChi2_.push_back(iMuonM->combinedQuality().trkKink);
      mumFracHits_.push_back(static_cast<double>(muTrackm->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackm->hitPattern().numberOfValidHits() +
															   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
															   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
															   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
      //                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackmTT, bestVtxReFit);
      //                     mumdxyVtx_.push_back(theDCAXVtx.second.value());
      //                     mumdzVtx_.push_back(muTrackmTT.track().dz( ));
      mumdxyBS_.push_back(muTrackmTT.track().dxy( (beamSpot.position() )));
      mumdzBS_.push_back(muTrackmTT.track().dz(  (beamSpot.position() )));
      
      //mumCat_.push_back(getMuCat(mum));
      
      mumNPixHits_.push_back(muTrackm->hitPattern().numberOfValidPixelHits());
      mumNPixLayers_.push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());  
      mumNTrkHits_.push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
      mumNTrkLayers_.push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
      if (mum.isGlobalMuon() == true) mumNMuonHits_.push_back(mum.globalTrack()->hitPattern().numberOfValidMuonHits());
      else mumNMuonHits_.push_back(0);
      mumNMatchStation_.push_back(mum.numberOfMatchedStations());
      
      
      //// #############
      //// # Save: mu+ #
      //// #############
      mupHighPurity_.push_back( (int) muTrackp->quality(reco::Track::highPurity));
      mupPt_.push_back(muTrackp->pt());
      mupEta_.push_back(muTrackp->eta());
      mupPhi_.push_back(muTrackp->phi());
      //mupCL_.push_back(TMath::Prob(muTrackpTT.chi2(), static_cast<int>(rint(muTrackpTT.ndof()))));
      mupNormChi2_.push_back(muTrackp->normalizedChi2());
      mupPx_.push_back(refitMupTT.track().momentum().x());
      mupPy_.push_back(refitMupTT.track().momentum().y());
      mupPz_.push_back(refitMupTT.track().momentum().z());
      
      mupDCABS_.push_back(DCAmupBS);
      mupDCABSE_.push_back(DCAmupBSErr);
      mupdxyBS_.push_back(muTrackpTT.track().dxy( (beamSpot.position() )));
      mupdzBS_.push_back(muTrackpTT.track().dz(  (beamSpot.position() )));
      
      //                     mupKinkChi2_.push_back(iMuonP->combinedQuality().trkKink);
      mupFracHits_.push_back(static_cast<double>(muTrackp->hitPattern().numberOfValidHits()) / 
                                                   (  muTrackp->hitPattern().numberOfValidHits() +
						      muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
						      muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
						      muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS) ) ) ;
      //                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackpTT, bestVtxReFit);
      //                     mupdxyVtx_.push_back(theDCAXVtx.second.value());
      //                     mupdzVtx_.push_back(muTrackpTT.track().dz(bestVtxReFit.position()));
      
      //mupCat_.push_back(getMuCat(mup));
      
      mupNPixHits_.push_back(muTrackp->hitPattern().numberOfValidPixelHits());
      mupNPixLayers_.push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());  
      mupNTrkHits_.push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
      mupNTrkLayers_.push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
      if (mup.isGlobalMuon() == true) mupNMuonHits_.push_back(mup.globalTrack()->hitPattern().numberOfValidMuonHits());
      else mupNMuonHits_.push_back(0);
      mupNMatchStation_.push_back(mup.numberOfMatchedStations());
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
    phoPFE_             .push_back(pf->energy());
    phoPFEt_            .push_back(pf->et());
    phoPFEta_           .push_back(pf->eta());
    phoPFPhi_           .push_back(pf->phi());
    
    nPFPho_++;
  } // PF photons loop
}


void BsToMuMuGammaNTuplizer::fillSC(edm::Event const& e) {
  edm::Handle<reco::SuperClusterCollection> barrelSCHandle;
  e.getByToken(MustacheSCBarrelCollection_, barrelSCHandle);

  edm::Handle<reco::SuperClusterCollection> endcapSCHandle;
  e.getByToken(MustacheSCEndcapCollection_, endcapSCHandle);

  for (auto const& scs : { *barrelSCHandle, *endcapSCHandle }) {
    for (auto const& sc : scs) {
      //if(abs(sc.eta())>2.4)continue;
      scE_.push_back(sc.energy());
      scRawE_.push_back(sc.rawEnergy());
      scEta_.push_back(sc.eta());
      scPhi_.push_back(sc.phi());

      ++nSC_;
    }
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(BsToMuMuGammaNTuplizer);
