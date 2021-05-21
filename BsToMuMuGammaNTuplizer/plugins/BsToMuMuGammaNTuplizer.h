#ifndef BsToMuMuGammaNTuplizer_h
#define BsToMuMuGammaNTuplizer_h

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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"



// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

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

#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/interface/Utils.h"


#include <TTree.h>

class BsToMuMuGammaNTuplizer : public edm::EDAnalyzer {

public:

   BsToMuMuGammaNTuplizer(const edm::ParameterSet&);
   virtual ~BsToMuMuGammaNTuplizer() {};

  bool printMsg;
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;

   void fillGenParticles (const edm::Event&);
   void fillMuons        (const edm::Event&, const edm::EventSetup&);
   void fillPhotons      (const edm::Event&, const edm::EventSetup&);
   void fillPFPhotons    (const edm::Event&, const edm::EventSetup&);
   void fillSC           (const edm::Event&);

     // switches
   bool doGenParticles_;
   bool doMuons_;
   bool doPhotons_;
   bool doPFPhotons_;
   bool doSuperClusters_;
  bool isMC;
  bool doBsToMuMuGamma;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::BeamSpot>                  beamSpotToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>   simGenTocken_;
  edm::EDGetTokenT<std::vector<reco::Photon>>  gedPhotonsCollection_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>>    pfPhotonsCollection_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>>   MustacheSCBarrelCollection_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>>   MustacheSCEndcapCollection_;
  edm::EDGetTokenT<reco::VertexCollection>          primaryVtxToken_;
  edm::EDGetTokenT<std::vector<reco::Muon>>         muonToken_;

  // input tags 
  //edm::InputTag genParticleSrc_;
  edm::InputTag gedPhotonSrc_;
  edm::InputTag pfPhotonSrc_;
  edm::InputTag MustacheSCBarrelSrc_;
  edm::InputTag MustacheSCEndcapSrc_;


  // selection cuts;
  double pTMinMuons;
  double cl_dimuon_vtx;
  double ls_max_dimuonBS          ;
  double dcaMax_dimuon_mumu           ;
  double dcaMax_muon_bs           ;
  double cosAlphaMax_dimuonBs    ;
  double MUMINPT           ;
  double etaMax_muon          ;
  double minDimuon_pt         ;
  double minDimuonInvariantMass    ;
  double maxDimuonInvariantMass    ;
  double trackIP_zMax_muon 		;
  double trackIP_rMax_muon	        ;
    
  // Dimuon Reco vars
  TrajectoryStateClosestToPoint theDCAXBS;
  ClosestApproachInRPhi ClosestApp;
  GlobalPoint XingPoint;


  RefCountedKinematicTree mumuVertexFitTree;
	 
  TLorentzVector bsDimuon_lv,sc_p4_;
  //LorentzVector scp4;
  reco::LeafCandidate scPhoton;
  reco::CompositeCandidate mmg;



   TTree*         theTree;

   // variables associated with tree branches
   UInt_t         run_;
   ULong64_t      event_;
   UInt_t         lumis_;
   Bool_t         isData_;
   Float_t        rho_;


  // # BeamSpot # //
  double beamspot_x,beamspot_y,beamspot_z,beamspot_x_error,beamspot_y_error,beamspot_z_error;
  double beamspot_dxdz,beamspot_dydz,beamspot_sigmaZ,beamspot_dxdz_error,beamspot_dydz_error,beamspot_sigmaZError;
  double beamspot_beamWidthX,beamspot_beamWidthY,beamspot_beamWidthX_error,beamspot_beamWidthY_error;
 
  // # offlinePrimaryVertices # //
  std::vector<bool> primaryVertex_isFake;
  std::vector<double> primaryVertex_x, primaryVertex_y,primaryVertex_z,primaryVertex_t;
  std::vector<double> primaryVertex_x_error, primaryVertex_y_error,primaryVertex_z_error,primaryVertex_t_error;
  std::vector<double> primaryVertex_ntracks,primaryVertex_ndof,primaryVertex_chi2,primaryVertex_normalizedChi2;

   // reco::GenParticle
  int  gen_nBs_, gen_nBsMuonM_, gen_nBsMuonP_ , gen_nBsPhoton_ ;
  std::vector<bool>   gen_hasAValid_candidate_;
  std::vector<double> gen_Bs_pt_,      gen_Bs_eta_,      gen_Bs_phi_,   gen_Bs_pz_,  gen_Bs_pdgId_;
  std::vector<double> gen_BsMuonM_pt_, gen_BsMuonM_eta_, gen_BsMuonM_phi_;
  std::vector<double> gen_BsMuonP_pt_, gen_BsMuonP_eta_, gen_BsMuonP_phi_;
  std::vector<double> gen_BsPhoton_pt_, gen_BsPhoton_eta_, gen_BsPhoton_phi_;


  // ### mu+ mu- variables ###
  std::vector<double>   mumuPt_;
  std::vector<double>   mumuEta_;
  std::vector<double>   mumuRapidity_;
  std::vector<double>   mumuPhi_;
  std::vector<double>   mumuMass_;
  std::vector<double>   mumuMassE_;
  std::vector<double>   mumuPx_;
  std::vector<double>   mumuPy_;
  std::vector<double>   mumuPz_;
  std::vector<double>   mumuDR_;

  // ### mu+ mu- Vtx ###
  std::vector<double>  mumuVtxCL_;
  std::vector<double>  mumuVtxX_;
  std::vector<double>  mumuVtxY_;
  std::vector<double>  mumuVtxZ_;
  std::vector<double>  mumuVtxChi2_;
  std::vector<double>  mumuVtxNdof_;
  std::vector<double>  mumuVtxProb_;
  std::vector<bool>    mumuVtxIsGoodFit_;
  std::vector<double>  mumuCosAlphaBS_;
  std::vector<double>  mumuCosAlphaBSE_; 
  std::vector<double>  mumuLBS_;
  std::vector<double>  mumuLBSE_;
  std::vector<double>  mumuDCA_;
  std::vector<double>   mumuLS_;
  std::vector<double>   mumuLSErr_;


  
  // ### mu- ###
  int 	               nMuM_; 
  std::vector<bool>    mumHighPurity_;
  std::vector<double>  mumPt_;
  std::vector<double>  mumEta_;
  std::vector<double>  mumPhi_;
  std::vector<double>  mumCL_; 
  std::vector<double>  mumNormChi2_;
  std::vector<double>  mumPx_;
  std::vector<double>  mumPy_;
  std::vector<double>  mumPz_;
  std::vector<double>  mumDCAVtx_;
  std::vector<double>  mumDCAVtxE_;
  std::vector<double>  mumDCABS_;
  std::vector<double>  mumDCABSE_;
  std::vector<double>  mumKinkChi2_;
  std::vector<double>  mumFracHits_;
  std::vector<double>  mumdxyBS_;
  std::vector<double>  mumdzBS_;
  std::vector<double>  mumMinIP2D_;
  std::vector<double>  mumMinIP2DE_;
  std::vector<double>  mumMinIP_;
  std::vector<double>  mumMinIPS_;
  std::vector<double>  mumDeltaRwithMC_;
  std::vector<std::string> mumCat_;
  std::vector<int>     mumCharge_;
  std::vector<int>     mumNPixHits_;
  std::vector<int>     mumNPixLayers_;
  std::vector<int>     mumNTrkHits_;
  std::vector<int>     mumNTrkLayers_;
  std::vector<int>     mumNMuonHits_;
  std::vector<int>     mumNMatchStation_;
  std::vector<float>   mumIso_;
  std::vector<float>   mumIsoPt_;
  std::vector<float>   mumIsodR_;
  std::vector<bool>    mum_isGlobalMuon_;
  std::vector<bool>    mum_isTrackerMuon_;
  std::vector<bool>    mum_StandAloneMuon_;
  std::vector<bool>    mum_isCaloMuon_;
  std::vector<bool>    mum_isPFMuon_;


  // ### mu+ ###
  int 	               nMuP_; 
  std::vector<bool>    mupHighPurity_;
  std::vector<double>  mupPt_;
  std::vector<double>  mupEta_;
  std::vector<double>  mupPhi_;
  std::vector<double>  mupCL_; 
  std::vector<double>  mupNormChi2_;
  std::vector<double>  mupPx_;
  std::vector<double>  mupPy_;
  std::vector<double>  mupPz_;
  std::vector<double>  mupDCAVtx_;
  std::vector<double>  mupDCAVtxE_;
  std::vector<double>  mupDCABS_;
  std::vector<double>  mupDCABSE_;
  std::vector<double>  mupKinkChi2_;
  std::vector<double>  mupFracHits_;
  std::vector<double>  mupdxyBS_;
  std::vector<double>  mupdzBS_;
  std::vector<double>  mupMinIP2D_;
  std::vector<double>  mupMinIP2DE_;
  std::vector<double>  mupMinIP_;
  std::vector<double>  mupMinIPS_;
  std::vector<double>  mupDeltaRwithMC_;
  std::vector<std::string> mupCat_;
  std::vector<int>     mupCharge_;
  std::vector<int>     mupNPixHits_;
  std::vector<int>     mupNPixLayers_;
  std::vector<int>     mupNTrkHits_;
  std::vector<int>     mupNTrkLayers_;
  std::vector<int>     mupNMuonHits_;
  std::vector<int>     mupNMatchStation_;
  std::vector<float>   mupIso_;
  std::vector<float>   mupIsoPt_;
  std::vector<float>   mupIsodR_;
  std::vector<bool>    mup_isGlobalMuon_;
  std::vector<bool>    mup_isTrackerMuon_;
  std::vector<bool>    mup_StandAloneMuon_;
  std::vector<bool>    mup_isCaloMuon_;
  std::vector<bool>    mup_isPFMuon_;


   // reco::Photon
   Int_t          nPho_;
   std::vector<float>  phoE_;
   std::vector<float>  phoEt_;
   std::vector<float>  phoEta_;
   std::vector<float>  phoPhi_;

  // reco::PFPhoton
   Int_t          nPFPho_;
   std::vector<float>  phoPFE_;
   std::vector<float>  phoPFEt_;
   std::vector<float>  phoPFEta_;
   std::vector<float>  phoPFPhi_;

   /* supercluster info */
   int nSC_;
   std::vector<float> scE_;
   std::vector<float> scEta_;
   std::vector<float> scPhi_;
   std::vector<float>  scX_;
   std::vector<float>  scY_;
   std::vector<float>  scZ_;
   std::vector<float>  scEtaWidth_;
   std::vector<float>  scPhiWidth_;
   std::vector<float>  scRawE_;
   std::vector<float>  scRawEt_;
};

#endif
