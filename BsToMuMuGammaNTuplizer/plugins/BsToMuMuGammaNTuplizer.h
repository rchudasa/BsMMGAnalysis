#ifndef BsToMuMuGammaNTuplizer_h
#define BsToMuMuGammaNTuplizer_h

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
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
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "FWCore/Common/interface/TriggerNames.h"


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
#include "FWCore/Utilities/interface/isFinite.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/interface/Utils.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/plugins/CaloTopologyBuilder.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include <TTree.h>
#include "TMath.h"

#define N_PV_MAX 500
#define N_TRK_MAX 5000
#define N_PF_MAX 10000
#define N_ECAL_CLUSTERS 1000
#define N_HCAL_CLUSTERS 1000
#define N_MUON_MAX 1000
#define N_DIMU_MAX 1000
#define N_GEN_MAX 4000
#define NMAX_MMG 1000

#include "BmmgInterface.h"

using namespace std;

typedef std::pair<const reco::MuonChamberMatch*, const reco::MuonSegmentMatch*> MatchPair;

class BsToMuMuGammaNTuplizer : public edm::EDAnalyzer {

 public:
  
  BsToMuMuGammaNTuplizer(const edm::ParameterSet&);
  virtual ~BsToMuMuGammaNTuplizer() {};
  
  Bool_t printMsg;
  
private:
  
  /////////////////////////////////////

  bool isGoodMuon(const pat::Muon& muon);
    
  KinematicFitResult vertexWithKinematicFitter(std::vector<const reco::Track*> trks, std::vector<float> masses);

  KinematicFitResult 
  vertexMuonsWithKinematicFitter(const pat::Muon& muon1,
				 const pat::Muon& muon2);

  KinematicFitResult 
  vertexWithKinematicFitter(const pat::Muon& muon1,
			    const pat::Muon& muon2,
			    const reco::PFCandidate& pfCand);

  pair<double,double> computeDCA(const reco::PFCandidate &kaon,
   				 reco::BeamSpot beamSpot);

  // Two track DOCA
  float 
  distanceOfClosestApproach( const reco::Track* track1,
			     const reco::Track* track2 );
  float 
  distanceOfClosestApproach( const reco::GenParticle* track1,
			     const reco::GenParticle* track2);
  
  // Track to vertex DOCA
  Measurement1D
  distanceOfClosestApproach( const reco::Track* track,
			     RefCountedKinematicVertex vertex);
  Measurement1D 
  distanceOfClosestApproach( const reco::Track* track,
			     const reco::Vertex& vertex);

  DisplacementInformationIn3D 
  compute3dDisplacement(const KinematicFitResult& fit,
			const reco::VertexCollection& vertices,
			bool closestIn3D = true);

  CloseTrackInfo 
  findTracksCompatibleWithTheVertex(const pat::Muon& muon1,
				    const pat::Muon& muon2,
				    const KinematicFitResult& fit,
				    double maxDoca=0.03,
				    std::vector<const reco::PFCandidate*> ignoreTracks = 
				    std::vector<const reco::PFCandidate*>());
  float
  computeTrkMuonIsolation(const pat::Muon& muon, 
			  const pat::Muon& the_other_muon,
			  unsigned int primaryVertexIndex,
			  float minPt=0.5, float dR=0.5,
			  std::vector<const reco::PFCandidate*> ignoreTracks = 
			  std::vector<const reco::PFCandidate*>());
  float
  computeTrkMuMuIsolation(const pat::Muon& muon1, 
			  const pat::Muon& muon2,
			  unsigned int primaryVertexIndex,
			  float minPt=0.9, float dR=0.7,
			  std::vector<const reco::PFCandidate*> ignoreTracks = 
			  std::vector<const reco::PFCandidate*>());

  float
  otherVertexMaxProb(const pat::Muon& muon1, 
		     const pat::Muon& muon2,
		     float min_pt = 0.5,
		     float max_doca = 0.1,
		     std::vector<const reco::PFCandidate*> ignoreTracks = 
		     std::vector<const reco::PFCandidate*>());
  void 
  injectHadronsThatMayFakeMuons(std::vector<MuonCand>& good_muon_candidates);

  void 
  injectBhhHadrons(std::vector<MuonCand>& good_muon_candidates);
//  float computeTrkMuonIsolation(const pat::Muon& the_muon, 
//                      const pat::Muon& the_other_muon, 
//					  unsigned int primaryVertexIndex,
//					  float minPt, float dR,
//					  std::vector<const reco::PFCandidate*> ignoreTracks);
//  float otherVertexMaxProb(const pat::Muon& muon1, 
//				     const pat::Muon& muon2,
//				     float minPt,
//				     float max_doca,
//				     std::vector<const reco::PFCandidate*> ignoreTracks);
//
  
  void addDimuonBranches();
  void fillDimuonBranches      (const edm::Event&, const edm::EventSetup&);
  /////////////////////////////////////



  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  
  void fillGenParticles (const edm::Event&);
  void fillMuons        (const edm::Event&, const edm::EventSetup&);
  void fillPhotons      (const edm::Event&, const edm::EventSetup&);
  void fillPFPhotons    (const edm::Event&, const edm::EventSetup&);
  void fillSC           (const edm::Event&, const edm::EventSetup&, reco::Vertex& pv);
  void fillHLT          (edm::Event const& );
  std::vector<Float_t> getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology);
  
  void fillMatchInfo(Int_t idxToFill, const reco::Muon &muon);
  void fillMatchInfoForStation(std::string prefix,Int_t idxToFill, const MatchPair& match);
  
  void addGenBranches();
  void fillGenBranches( const edm::Event&);

  void addMuonBranches();
  void fillMuonBranches( const edm::Event&, const edm::EventSetup& );
  
  void addMMGBranches();
  void fillBmmgBranchs(Int_t nDimu,Int_t phoIdx,Int_t scIdx , Float_t mmg_mass);
  void addJPsiGammaBranches();
  void fillJPsiGammaBranches(Int_t nDimu,Int_t phoIdx,Int_t scIdx , Float_t mmg_mass);

  void addPF_MMGBranches();
  void fillPF_BmmgBranchs(Int_t nDimu,Int_t phoIdx,Float_t mmg_mass);
  void addPF_JPsiGammaBranches();
  void fillPF_JPsiGammaBranches(Int_t nDimu,Int_t phoIdx, Float_t mmg_mass);

  void fillPFCandiateCollection( const edm::Event&, const edm::EventSetup& );
  void addParticleFlowBranches();
  void fillGeneralTrackCollectionBranches( const edm::Event&, const edm::EventSetup& );
  void addGeneralTracksBranches();
  void fillHCALClusterCollection( const edm::Event&, const edm::EventSetup& );
  void addHCALCluserBranches();
  void fillECALClusterCollection( const edm::Event&, const edm::EventSetup& );
  void addECALCluserBranches();
  void addPrimaryVertexBranches();
  void fillPrimaryVertexBranches(const edm::Event&, const edm::EventSetup&);
  Float_t reduceFloat(Float_t val, int bits);
  static int calDIEta(int iEta1, int iEta2);
  static int calDIPhi(int iPhi1, int iPhi2);
  Float_t getMinEnergyHCAL_(HcalDetId id) const;

  // master switches
  Bool_t         isMC;
  Bool_t         isRECO;
  Bool_t         isMiniAOD;
  Bool_t         Run2_2018_;
  
  // switches
  Bool_t doGenParticles_;
  Bool_t doMuons_;
  Bool_t doDimuons_;
  Bool_t doMuMuGamma;
  Bool_t doJPsiGamma;
  Bool_t doPhotons_;
  Bool_t doPFPhotons_;
  Bool_t doSuperClusters_;
  Bool_t doBsToMuMuGamma;
  Bool_t doCompression_;
  Bool_t doParticleFlow;
  Bool_t doGeneralTracks;
  Bool_t doECALClusters;
  Bool_t doHCALClusters;
  Bool_t doPrimaryVetrices;
  Bool_t doBeamSpot;
  Bool_t doFlatPt_;
  Bool_t doHLT;
  int maxDIEta_=5;
  int maxDIPhi_=5;
  
  // ----------

  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::GenParticleCollection>     genParticlesCollection_;
  edm::EDGetTokenT<std::vector<reco::Photon>>       gedPhotonsCollection_;
  edm::EDGetTokenT<edm::View<reco::PFCandidate>>    pfPhotonsCollection_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> MustacheSCBarrelCollection_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> MustacheSCEndcapCollection_;
  edm::EDGetTokenT<reco::GsfElectronCollection>     gsfElectronToken_;
  edm::EDGetTokenT<edm::TriggerResults>             triggerBits_;
  edm::EDGetTokenT<std::vector<reco::Muon>>         muonToken_;
  edm::EDGetTokenT<reco::TrackCollection>           generalTracksCollection_;
  edm::EDGetTokenT<reco::PFClusterCollection>       ecalClusterCollection_;
  edm::EDGetTokenT<reco::PFClusterCollection>       hcalClusterCollection_;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> hbheRechitToken_; 
  
  edm::EDGetTokenT<reco::PFCandidateCollection>     pfCandidateCollection_;
  edm::Handle<reco::PFCandidateCollection> pfCandidateHandle;
  //edm::Handle<edm::View<reco::PFCandidateCollection>> pfCandHandle_;
  //edm::Handle<reco::PFCandidateCollection> pfCandHandle_;
  //edm::Handle<edm::View<reco::PFCandidate> > pfCandidateHandle;
  
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;
  edm::ESHandle<TransientTrackBuilder> theTTBuilder_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  
  edm::EDGetTokenT<reco::BeamSpot>                  beamSpotToken_;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  
  edm::EDGetTokenT<reco::VertexCollection>          primaryVtxToken_;
  edm::Handle<reco::VertexCollection> pvHandle_;
  
  const Int_t energyMatrixSize_;
  Int_t energyMatrixSizeFull_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRechitToken_; 
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> eeRechitToken_; 

  // input tags 
  edm::InputTag gedPhotonSrc_;
  edm::InputTag pfPhotonSrc_;
  edm::InputTag MustacheSCBarrelSrc_;
  edm::InputTag MustacheSCEndcapSrc_;


  // selection cuts;
  Float_t pTMinMuons;
  Float_t pTMinPFPhotons;
  Float_t cl_dimuon_vtx;
  Float_t ls_max_dimuonBS          ;
  Float_t dcaMax_dimuon_mumu           ;
  Float_t dcaMax_muon_bs           ;
  Float_t cosAlphaMax_dimuonBs    ;
  Float_t muMinpt           ;
  Float_t etaMax_muon          ;
  Float_t minDimuon_pt         ;
  Float_t minDimuonInvariantMass    ;
  Float_t maxDimuonInvariantMass    ;
  Float_t trackIP_zMax_muon 		;
  Float_t trackIP_rMax_muon	        ;
  Float_t maxTwoTrackDOCA_          ;
  Float_t minJPsiMass_;
  Float_t maxJPsiMass_;
  Float_t minJPsiGammaMass_;
  Float_t maxJPsiGammaMass_;
  Float_t maxBsMuMuGammaMass_;
  Float_t minBsMuMuGammaMass_;
  bool diMuonCharge_ ;

  // Utility 
  Utils* Utility;
    
  // Dimuon Reco vars
  TrajectoryStateClosestToPoint theDCAXBS;
  ClosestApproachInRPhi ClosestApp;
  GlobalPoint XingPoint;


  RefCountedKinematicTree mumuVertexFitTree;
	 
  TLorentzVector bsDimuon_lv,sc_p4_;
  reco::LeafCandidate scPhoton;
  reco::CompositeCandidate mmg;


  // for showershape
  int nBits_;
  std::vector<Float_t> locCov_;
  std::vector<Float_t> full5x5_locCov_;
  std::vector<Float_t> showerShapes_;
 
  TTree*         theTree;

  // variables associated with tree branches
  UInt_t         run_;
  ULong64_t      event_;
  UInt_t         lumis_;
  Bool_t         isData_;
  Float_t        rho_;


  // # BeamSpot # //
  Float_t beamspot_x_,beamspot_y_,beamspot_z_,beamspot_x_error_,beamspot_y_error_,beamspot_z_error_;
  Float_t beamspot_dxdz_,beamspot_dydz_,beamspot_sigmaZ_,beamspot_dxdz_error_,beamspot_dydz_error_,beamspot_sigmaZError_;
  Float_t beamspot_beamWidthX_,beamspot_beamWidthY_,beamspot_beamWidthX_error_,beamspot_beamWidthY_error_;
  Float_t beamspot_covXX,beamspot_covXY,beamspot_covXZ,beamspot_covYY,beamspot_covYZ,beamspot_covZZ;  
  // # offlinePrimaryVertices # //
  int nPrimaryVertex_;
  std::vector<Bool_t> primaryVertex_isFake_;
  std::vector<Float_t> primaryVertex_x_, primaryVertex_y_,primaryVertex_z_,primaryVertex_t_;
  std::vector<Float_t> primaryVertex_x_error_, primaryVertex_y_error_,primaryVertex_z_error_,primaryVertex_t_error_;
  std::vector<Float_t> primaryVertex_ntracks_,primaryVertex_ndof_,primaryVertex_chi2_,primaryVertex_normalizedChi2_;
  std::vector<Float_t> primaryVertex_covXX,primaryVertex_covXY,primaryVertex_covXZ,primaryVertex_covYY,primaryVertex_covYZ,primaryVertex_covZZ;  

  // # Trigger #
  std::vector<std::string>  trigTable;
  std::vector<std::string>  trigNames;
  std::vector<Bool_t>         trigResult;
  std::vector<int>          trigPrescales;
  
  std::vector<std::string>  *TrigTable_store;
  std::vector<Bool_t>         *TrigResult_store;
  std::vector<int>          *TrigPrescales_store;
  
  std::vector<std::string>  l1Table;
  std::vector<int>          l1Prescales;

  void SetupTriggerStorageVectors();
  void SetupTriggerBranches();
  void FillTrggerBranches();
  void ClearTrggerStorages();

  // generic storage dictionary

  std::map<std::string,Int_t > storageMapInt;
  std::map<std::string,Int_t *> storageMapIntArray;
  std::map<std::string,Float_t * > storageMapFloatArray;
  std::map<std::string,Bool_t * > storageMapBoolArray;
  std::map<std::string,uint64_t * > storageMapUint64Array;


  // reco::GenParticle
  int  gen_nBs_, gen_nBsMuonM_, gen_nBsMuonP_ , gen_nBsPhoton_ ;
  std::vector<Float_t> gen_Bs_pt_,  gen_Bs_energy_,    gen_Bs_eta_,      gen_Bs_phi_,   gen_Bs_pz_,  gen_Bs_pdgId_;
  std::vector<Float_t> gen_BsMuonM_pt_, gen_BsMuonM_eta_, gen_BsMuonM_phi_;
  std::vector<Float_t> gen_BsMuonP_pt_, gen_BsMuonP_eta_, gen_BsMuonP_phi_;
  std::vector<Float_t> gen_BsPhoton_pt_, gen_BsPhoton_energy_, gen_BsPhoton_eta_, gen_BsPhoton_phi_;


  // ### mu+ mu- variables ###
  std::vector<Float_t>   mumuPt_;
  std::vector<Float_t>   mumuEta_;
  std::vector<Float_t>   mumuRapidity_;
  std::vector<Float_t>   mumuPhi_;
  std::vector<Float_t>   mumuMass_;
  std::vector<Float_t>   mumuMassE_;
  std::vector<Float_t>   mumuPx_;
  std::vector<Float_t>   mumuPy_;
  std::vector<Float_t>   mumuPz_;
  std::vector<Float_t>   mumuDR_;
  std::vector<Float_t>   mumuCovPxPx_;
  std::vector<Float_t>   mumuCovPxPy_;
  std::vector<Float_t>   mumuCovPxPz_;
  std::vector<Float_t>   mumuCovPyPy_;
  std::vector<Float_t>   mumuCovPyPz_;
  std::vector<Float_t>   mumuCovPzPz_;
  

  std::vector<int>   mumuParentMuM_;
  std::vector<int>   mumuParentMuP_;
  // ### mu+ mu- Vtx ###
  std::vector<Float_t>  mumuVtxCL_;
  std::vector<Float_t>  mumuVtxX_;
  std::vector<Float_t>  mumuVtxY_;
  std::vector<Float_t>  mumuVtxZ_;
  std::vector<Float_t>  mumuVtxCovXX_;
  std::vector<Float_t>  mumuVtxCovXY_;
  std::vector<Float_t>  mumuVtxCovXZ_;
  std::vector<Float_t>  mumuVtxCovYY_;
  std::vector<Float_t>  mumuVtxCovYZ_;
  std::vector<Float_t>  mumuVtxCovZZ_;
  
  std::vector<Float_t>  mumuVtxChi2_;
  std::vector<Float_t>  mumuVtxNdof_;

  std::vector<Float_t>  mumuCosAlphaBS_;
  std::vector<Float_t>  mumuCosAlphaBSE_; 
  std::vector<Float_t>  mumuLBS_;
  std::vector<Float_t>  mumuLBSE_;
  std::vector<Float_t>  mumuDCA_;
  


  
  // ### mu- ###
  int 	               nMuM_; 
  std::vector<Bool_t>    mumHighPurity_;
  std::vector<Float_t>  mumPt_;
  std::vector<Float_t>  mumEta_;
  std::vector<Float_t>  mumPhi_;
  std::vector<Float_t>  mumCL_; 
  std::vector<Float_t>  mumNormChi2_;
  std::vector<Float_t>  mumVx_;
  std::vector<Float_t>  mumVy_;
  std::vector<Float_t>  mumVz_;

  std::vector<Float_t>  mumDCABS_;
  std::vector<Float_t>  mumDCABSE_;

  std::vector<Float_t>  mumFracHits_;
  std::vector<Float_t>  mumdxyBS_;
  std::vector<Float_t>  mumdzBS_;
 

  std::vector<int>     mumIdx_;
  std::vector<int>     mumCharge_;
  std::vector<int>     mumNPixHits_;
  std::vector<int>     mumNPixLayers_;
  std::vector<int>     mumNTrkHits_;
  std::vector<int>     mumNTrkLayers_;
  std::vector<int>     mumNMuonHits_;
  std::vector<int>     mumNMatchStation_;
  std::vector<Bool_t>    mum_isGlobalMuon_;
  std::vector<Bool_t>    mum_isTrackerMuon_;
  std::vector<Bool_t>    mum_StandAloneMuon_;
  std::vector<Bool_t>    mum_isCaloMuon_;
  std::vector<Bool_t>    mum_isPFMuon_;

  std::vector<uint64_t> mum_selector_; 
  std::vector<Bool_t>	mum_isIsolationValid_;
  std::vector<Bool_t>	mum_isPFIsolationValid_;
  
  std::vector<Float_t>  mum_isolationR03_trackSumPt_;
  std::vector<Float_t>  mum_isolationR03_trackEcalSumEt_;
  std::vector<Float_t>  mum_isolationR03_trackHcalSumEt_;
  std::vector<Float_t>  mum_isolationR03_trackHoSumEt_;
  std::vector<int>     mum_isolationR03_trackNTracks_;
  std::vector<int>     mum_isolationR03_trackNJets_;
  std::vector<Float_t>  mum_isolationR03_trackerVetoSumPt_;
  std::vector<Float_t>  mum_isolationR03_emVetoSumEt_;
  std::vector<Float_t>  mum_isolationR03_hadVetoSumEt_;
  std::vector<Float_t>  mum_isolationR03_hoVetoEt_;
  
  std::vector<Float_t>  mum_isolationR05_trackSumPt_;
  std::vector<Float_t>  mum_isolationR05_trackEcalSumEt_;
  std::vector<Float_t>  mum_isolationR05_trackHcalSumEt_;
  std::vector<Float_t>  mum_isolationR05_trackHoSumEt_;
  std::vector<int>     mum_isolationR05_trackNTracks_;
  std::vector<int>     mum_isolationR05_trackNJets_;
  std::vector<Float_t>  mum_isolationR05_trackerVetoSumPt_;
  std::vector<Float_t>  mum_isolationR05_emVetoSumEt_;
  std::vector<Float_t>  mum_isolationR05_hadVetoSumEt_;
  std::vector<Float_t>  mum_isolationR05_hoVetoEt_;
  
  std::vector<Float_t>  mum_PFIsolationR03_sumChargedHadronPt_;
  std::vector<Float_t>  mum_PFIsolationR03_sumChargedParticlePt_;
  std::vector<Float_t>  mum_PFIsolationR03_sumNeutralHadronEt_;
  std::vector<Float_t>  mum_PFIsolationR03_sumPhotonEt_;
  std::vector<Float_t>  mum_PFIsolationR03_sumNeutralHadronEtHighThreshold_;
  std::vector<Float_t>  mum_PFIsolationR03_sumPhotonEtHighThreshold_;
  std::vector<Float_t>  mum_PFIsolationR03_sumPUPt_;
  
  std::vector<Float_t>  mum_PFIsolationR04_sumChargedHadronPt_;
  std::vector<Float_t>  mum_PFIsolationR04_sumChargedParticlePt_;
  std::vector<Float_t>  mum_PFIsolationR04_sumNeutralHadronEt_;
  std::vector<Float_t>  mum_PFIsolationR04_sumPhotonEt_;
  std::vector<Float_t>  mum_PFIsolationR04_sumNeutralHadronEtHighThreshold_;
  std::vector<Float_t>  mum_PFIsolationR04_sumPhotonEtHighThreshold_;
  std::vector<Float_t>  mum_PFIsolationR04_sumPUPt_;

  // ### mu+ ###
  int 	               nMuP_; 
  std::vector<Bool_t>    mupHighPurity_;
  std::vector<Float_t>  mupPt_;
  std::vector<Float_t>  mupEta_;
  std::vector<Float_t>  mupPhi_;
  std::vector<Float_t>  mupCL_; 
  std::vector<Float_t>  mupNormChi2_;
  std::vector<Float_t>  mupVx_;
  std::vector<Float_t>  mupVy_;
  std::vector<Float_t>  mupVz_;
 
  std::vector<Float_t>  mupDCABS_;
  std::vector<Float_t>  mupDCABSE_;

  std::vector<Float_t>  mupFracHits_;
  std::vector<Float_t>  mupdxyBS_;
  std::vector<Float_t>  mupdzBS_;
 
  std::vector<int>     mupIdx_;
  std::vector<int>     mupCharge_;
  std::vector<int>     mupNPixHits_;
  std::vector<int>     mupNPixLayers_;
  std::vector<int>     mupNTrkHits_;
  std::vector<int>     mupNTrkLayers_;
  std::vector<int>     mupNMuonHits_;
  std::vector<int>     mupNMatchStation_;
  std::vector<Bool_t>    mup_isGlobalMuon_;
  std::vector<Bool_t>    mup_isTrackerMuon_;
  std::vector<Bool_t>    mup_StandAloneMuon_;
  std::vector<Bool_t>    mup_isCaloMuon_;
  std::vector<Bool_t>    mup_isPFMuon_;

  std::vector<uint64_t> mup_selector_; 
  std::vector<Bool_t>	mup_isIsolationValid_;
  std::vector<Bool_t>	mup_isPFIsolationValid_;
  
  std::vector<Float_t>  mup_isolationR03_trackSumPt_;
  std::vector<Float_t>  mup_isolationR03_trackEcalSumEt_;
  std::vector<Float_t>  mup_isolationR03_trackHcalSumEt_;
  std::vector<Float_t>  mup_isolationR03_trackHoSumEt_;
  std::vector<int>     mup_isolationR03_trackNTracks_;
  std::vector<int>     mup_isolationR03_trackNJets_;
  std::vector<Float_t>  mup_isolationR03_trackerVetoSumPt_;
  std::vector<Float_t>  mup_isolationR03_emVetoSumEt_;
  std::vector<Float_t>  mup_isolationR03_hadVetoSumEt_;
  std::vector<Float_t>  mup_isolationR03_hoVetoEt_;
  
  std::vector<Float_t>  mup_isolationR05_trackSumPt_;
  std::vector<Float_t>  mup_isolationR05_trackEcalSumEt_;
  std::vector<Float_t>  mup_isolationR05_trackHcalSumEt_;
  std::vector<Float_t>  mup_isolationR05_trackHoSumEt_;
  std::vector<int>     mup_isolationR05_trackNTracks_;
  std::vector<int>     mup_isolationR05_trackNJets_;
  std::vector<Float_t>  mup_isolationR05_trackerVetoSumPt_;
  std::vector<Float_t>  mup_isolationR05_emVetoSumEt_;
  std::vector<Float_t>  mup_isolationR05_hadVetoSumEt_;
  std::vector<Float_t>  mup_isolationR05_hoVetoEt_;
  
  std::vector<Float_t>  mup_PFIsolationR03_sumChargedHadronPt_;
  std::vector<Float_t>  mup_PFIsolationR03_sumChargedParticlePt_;
  std::vector<Float_t>  mup_PFIsolationR03_sumNeutralHadronEt_;
  std::vector<Float_t>  mup_PFIsolationR03_sumPhotonEt_;
  std::vector<Float_t>  mup_PFIsolationR03_sumNeutralHadronEtHighThreshold_;
  std::vector<Float_t>  mup_PFIsolationR03_sumPhotonEtHighThreshold_;
  std::vector<Float_t>  mup_PFIsolationR03_sumPUPt_;
  
  std::vector<Float_t>  mup_PFIsolationR04_sumChargedHadronPt_;
  std::vector<Float_t>  mup_PFIsolationR04_sumChargedParticlePt_;
  std::vector<Float_t>  mup_PFIsolationR04_sumNeutralHadronEt_;
  std::vector<Float_t>  mup_PFIsolationR04_sumPhotonEt_;
  std::vector<Float_t>  mup_PFIsolationR04_sumNeutralHadronEtHighThreshold_;
  std::vector<Float_t>  mup_PFIsolationR04_sumPhotonEtHighThreshold_;
  std::vector<Float_t>  mup_PFIsolationR04_sumPUPt_;


  // reco::Photon
  Int_t          nPho_;
  std::vector<Float_t>  phoE_;
  std::vector<Float_t>  phoEt_;
  std::vector<Float_t>  phoEta_;
  std::vector<Float_t>  phoPhi_;

  std::vector<Float_t>  phoSigmaE_;
  std::vector<Float_t>  phoCalibE_;
  std::vector<Float_t>  phoCalibEt_;
  std::vector<Float_t>  phoSCE_;
  std::vector<Float_t>  phoSCEt_;
  std::vector<Float_t>  phoSCRawE_;
  std::vector<Float_t>  phoESEnP1_;
  std::vector<Float_t>  phoESEnP2_;
  std::vector<Float_t>  phoSCEta_;
  std::vector<Float_t>  phoSCPhi_;
  std::vector<Float_t>  phoSCEtaWidth_;
  std::vector<Float_t>  phoSCPhiWidth_;
  std::vector<Float_t>  phoSCBrem_;
  std::vector<int>    phohasPixelSeed_;
  std::vector<int>    phoEleVeto_;
  std::vector<Float_t>  phoR9_;
  std::vector<Float_t>  phoHoverE_;
  std::vector<Float_t>  phoESEffSigmaRR_;

  std::vector<Float_t>  phoSigmaIEtaIEtaFull5x5_;
  std::vector<Float_t>  phoSigmaIEtaIPhiFull5x5_;
  std::vector<Float_t>  phoSigmaIPhiIPhiFull5x5_;
  std::vector<Float_t>  phoE2x2Full5x5_;
  std::vector<Float_t>  phoE5x5Full5x5_;
  std::vector<Float_t>  phoR9Full5x5_;

  std::vector<Float_t>  phoPFChIso_;
  std::vector<Float_t>  phoPFPhoIso_;
  std::vector<Float_t>  phoPFNeuIso_;
  std::vector<Float_t>  phoEcalPFClusterIso_;
  std::vector<Float_t>  phoHcalPFClusterIso_;
  std::vector<Float_t>  phoIDMVA_;

  std::vector<Float_t>  phoSeedTime_;
  std::vector<Float_t>  phoSeedEnergy_;
  std::vector<Float_t>  phoMIPTotEnergy_;
  std::vector<Float_t>  phoMIPChi2_;
  std::vector<Float_t>  phoMIPSlope_;
  std::vector<Float_t>  phoMIPIntercept_;
  std::vector<Float_t>  phoMIPNhitCone_;
  std::vector<Float_t>  phoMIPIsHalo_;


  // reco::PFPhoton
  Int_t          nPFPho_;
  std::vector<Float_t>  phoPFE_;
  std::vector<Float_t>  phoPFEt_;
  std::vector<Float_t>  phoPFEta_;
  std::vector<Float_t>  phoPFPhi_;

  /* supercluster info */
  int nSC_;
  std::vector<Float_t> scE_;
  std::vector<Float_t> scEta_;
  std::vector<Float_t> scPhi_;
  std::vector<Float_t>  scX_;
  std::vector<Float_t>  scY_;
  std::vector<Float_t>  scZ_;
  std::vector<Float_t>  scEtaWidth_;
  std::vector<Float_t>  scPhiWidth_;
  std::vector<Float_t>  scRawE_;
  std::vector<Float_t>  scRawEt_;
  std::vector<Float_t>  scMinDrWithGsfElectornSC_;
  std::vector< Bool_t>  scFoundGsfMatch_;
  std::vector<Float_t> superCluster_e5x5_;
  std::vector<Float_t> superCluster_e2x2Ratio_;
  std::vector<Float_t> superCluster_e3x3Ratio_;
  std::vector<Float_t> superCluster_eMaxRatio_;
  std::vector<Float_t> superCluster_e2ndRatio_;
  std::vector<Float_t> superCluster_eTopRatio_;
  std::vector<Float_t> superCluster_eRightRatio_;
  std::vector<Float_t> superCluster_eBottomRatio_;
  std::vector<Float_t> superCluster_eLeftRatio_;
  std::vector<Float_t> superCluster_e2x5MaxRatio_;
  std::vector<Float_t> superCluster_e2x5TopRatio_;
  std::vector<Float_t> superCluster_e2x5RightRatio_;
  std::vector<Float_t> superCluster_e2x5BottomRatio_;
  std::vector<Float_t> superCluster_e2x5LeftRatio_;
  std::vector<Float_t> superCluster_swissCross_;
  std::vector<Float_t> superCluster_r9_;
  std::vector<Float_t> superCluster_sigmaIetaIeta_; 
  std::vector<Float_t> superCluster_sigmaIetaIphi_; 
  std::vector<Float_t> superCluster_sigmaIphiIphi_; 
  std::vector<Float_t> superCluster_full5x5_e5x5_;
  std::vector<Float_t> superCluster_full5x5_e2x2Ratio_;
  std::vector<Float_t> superCluster_full5x5_e3x3Ratio_;
  std::vector<Float_t> superCluster_full5x5_eMaxRatio_;
  std::vector<Float_t> superCluster_full5x5_e2ndRatio_;
  std::vector<Float_t> superCluster_full5x5_eTopRatio_;
  std::vector<Float_t> superCluster_full5x5_eRightRatio_;
  std::vector<Float_t> superCluster_full5x5_eBottomRatio_;
  std::vector<Float_t> superCluster_full5x5_eLeftRatio_;
  std::vector<Float_t> superCluster_full5x5_e2x5MaxRatio_;
  std::vector<Float_t> superCluster_full5x5_e2x5TopRatio_;
  std::vector<Float_t> superCluster_full5x5_e2x5RightRatio_;
  std::vector<Float_t> superCluster_full5x5_e2x5BottomRatio_;
  std::vector<Float_t> superCluster_full5x5_e2x5LeftRatio_;
  std::vector<Float_t> superCluster_full5x5_swissCross_;
  std::vector<Float_t> superCluster_full5x5_r9_;
  std::vector<Float_t> superCluster_full5x5_sigmaIetaIeta_; 
  std::vector<Float_t> superCluster_full5x5_sigmaIetaIphi_; 
  std::vector<Float_t> superCluster_full5x5_sigmaIphiIphi_;   

   std::vector<Float_t> scE5x5_;
   std::vector<Float_t> scE2x2Ratio_;
   std::vector<Float_t> scE3x3Ratio_;
   std::vector<Float_t> scEMaxRatio_;
   std::vector<Float_t> scE2ndRatio_;
   std::vector<Float_t> scETopRatio_;
   std::vector<Float_t> scERightRatio_;
   std::vector<Float_t> scEBottomRatio_;
   std::vector<Float_t> scELeftRatio_;
   std::vector<Float_t> scE2x5MaxRatio_;
   std::vector<Float_t> scE2x5TopRatio_;
   std::vector<Float_t> scE2x5RightRatio_;
   std::vector<Float_t> scE2x5BottomRatio_;
   std::vector<Float_t> scE2x5LeftRatio_;
   std::vector<Float_t> scSwissCross_;
   std::vector<Float_t> scR9_;
   std::vector<Float_t> scSigmaIetaIeta_; 
   std::vector<Float_t> scSigmaIetaIphi_; 
   std::vector<Float_t> scSigmaIphiIphi_; 
   std::vector<Float_t> scFull5x5_e5x5_;
   std::vector<Float_t> scFull5x5_e2x2Ratio_;
   std::vector<Float_t> scFull5x5_e3x3Ratio_;
   std::vector<Float_t> scFull5x5_eMaxRatio_;
   std::vector<Float_t> scFull5x5_e2ndRatio_;
   std::vector<Float_t> scFull5x5_eTopRatio_;
   std::vector<Float_t> scFull5x5_eRightRatio_;
   std::vector<Float_t> scFull5x5_eBottomRatio_;
   std::vector<Float_t> scFull5x5_eLeftRatio_;
   std::vector<Float_t> scFull5x5_e2x5MaxRatio_;
   std::vector<Float_t> scFull5x5_e2x5TopRatio_;
   std::vector<Float_t> scFull5x5_e2x5RightRatio_;
   std::vector<Float_t> scFull5x5_e2x5BottomRatio_;
   std::vector<Float_t> scFull5x5_e2x5LeftRatio_;
   std::vector<Float_t> scFull5x5_swissCross_;
   std::vector<Float_t> scFull5x5_r9_;
   std::vector<Float_t> scFull5x5_sigmaIetaIeta_; 
   std::vector<Float_t> scFull5x5_sigmaIetaIphi_; 
   std::vector<Float_t> scFull5x5_sigmaIphiIphi_;   
 
 
   int nhcalRechit_;
   std::vector<Float_t> hcalRechitIEta_;
   std::vector<Float_t> hcalRechitIPhi_;
   std::vector<Float_t> hcalRechitEnergy_;
 
   std::vector<Float_t>  scNHcalRecHitInDIEta2IPhi2;
   std::vector<Float_t>  scEFromHcalRecHitInDIEta2IPhi2;
   
   std::vector<Float_t>  scNHcalRecHitInDIEta5IPhi5;
   std::vector<Float_t>  scEFromHcalRecHitInDIEta5IPhi5;
   
   std::vector<Float_t>  scPFChIso1_;
   std::vector<Float_t>  scPFChIso2_;
   std::vector<Float_t>  scPFChIso3_;
   std::vector<Float_t>  scPFChIso4_;
   std::vector<Float_t>  scPFChIso5_;
   
   std::vector<Float_t>  scPFPhoIso1_;
   std::vector<Float_t>  scPFPhoIso2_;
   std::vector<Float_t>  scPFPhoIso3_;
   std::vector<Float_t>  scPFPhoIso4_;
   std::vector<Float_t>  scPFPhoIso5_;
   
   std::vector<Float_t>  scPFNeuIso1_;
   std::vector<Float_t>  scPFNeuIso2_;
   std::vector<Float_t>  scPFNeuIso3_;
   std::vector<Float_t>  scPFNeuIso4_;
   std::vector<Float_t>  scPFNeuIso5_;



    Int_t               nMC_;
    std::vector<int>    mcPID_;
    std::vector<int>    mcStatus_;
    std::vector<Float_t>  mcVtx_x_;
    std::vector<Float_t>  mcVtx_y_;
    std::vector<Float_t>  mcVtx_z_;
    std::vector<Float_t>  mcPt_;
    std::vector<Float_t>  mcEta_;
    std::vector<Float_t>  mcPhi_;
    std::vector<Float_t>  mcE_;
    std::vector<Float_t>  mcEt_;
    std::vector<Float_t>  mcMass_;
    std::vector<Float_t> scEt_;

};

#endif
