#ifndef BsToMuMuGammaNTuplizer_h
#define BsToMuMuGammaNTuplizer_h

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
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
#define N_PF_PHOTONS 300
#define N_PHOTONS_MAX 300
#define N_SC_MAX 400
#define N_GEN_MAX 4000
#define N_GEN_PARTS_ALL 500        
#define NMAX_MMG 1000
#define N_COMPOSIT_PART_MAX 200
#define N_L3MUON 500

#define N_COMPRESS_MAX 100  

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
  findTracksCompatibleWithTheVertex(const reco::Muon& muon1,
				    const reco::Muon& muon2,
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

  const reco::Candidate* getGenParticle(const reco::Candidate* cand);
  GenMatchInfo getGenMatchInfo( const reco::Muon& muon1,
						const reco::Muon& muon2,
						const reco::PFCandidate* kaon1=nullptr,
						const reco::PFCandidate* kaon2=nullptr,
						const reco::PFCandidate* photon=nullptr);

  void addGenBranches(TString sTag,TString nTag);
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
    KinematicFitResult fitBToKJPsiMuMu( RefCountedKinematicParticle refitMuMu,
				   const reco::PFCandidate &kaon,
				   bool applyJpsiMassConstraint);
    KinematicFitResult fitBToKMuMu( RefCountedKinematicTree tree,
				  const reco::PFCandidate& kaon,
				  float mass_constraint);
  
    void fillBtoMuMuKInfo(Int_t mumuIdx ,
					  const edm::Event& iEvent,
					  const KinematicFitResult& kinematicMuMuVertexFit,
					  const pat::Muon& muon1,
					  const pat::Muon& muon2,
					  const reco::PFCandidate & kaon
					) ;
  void addPhotonBranches();
  void fillPhotons      (const edm::Event&, const edm::EventSetup&);
  
  void addPFPhotonBranches();
  void fillPFPhotons      (const edm::Event&, const edm::EventSetup&);
  
  void addSCBranches();
  void fillSC           (const edm::Event&, const edm::EventSetup&, reco::Vertex& pv);

  void addDimuonBranches();
  void fillDimuonBranches      (const edm::Event&, const edm::EventSetup&);
  
  void  addMuMuKBranches();
  void addCompositParticleBranches(TString tag,TString nTag);
  void fillCompositParticleBranches( std::string tag, Int_t idx ,const KinematicFitResult &fit ,
                        const DisplacementInformationIn3D &displacement , CloseTrackInfo &closeTracks );
  /////////////////////////////////////



  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  
  void fillGenParticles (const edm::Event&);
  void fillMuons        (const edm::Event&, const edm::EventSetup&);
  std::vector<Float_t> getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology);
  void fillHLT          (edm::Event const& );
  
  void addHLTObjectBranches();
  void fillHLTL3MuonBranches(  std::vector<Int_t> trigIdx, std::vector<Int_t> &Keys,const trigger::TriggerEvent & triggerObj );
  
  void fillMatchInfo(Int_t idxToFill, const reco::Muon &muon);
  void fillMatchInfoForStation(std::string prefix,Int_t idxToFill, const MatchPair& match);
  
  void addGenBranches();
  void fillGenBranches( const edm::Event&);

  void addMuonBranches();
  void fillMuonBranches( const edm::Event&, const edm::EventSetup& );
  
  void addMMGBranches();
  void fillBmmgBranchs(Int_t nDimu,Int_t phoIdx,Int_t scIdx , Float_t mmg_mass, GenMatchInfo *pa=nullptr);
  void addJPsiGammaBranches();
  void fillJPsiGammaBranches(Int_t nDimu,Int_t phoIdx,Int_t scIdx , Float_t mmg_mass,GenMatchInfo *d=nullptr);

  void addPF_MMGBranches();
  void fillPF_BmmgBranchs(Int_t nDimu,Int_t phoIdx,Float_t mmg_mass, GenMatchInfo *genm=nullptr);
  void addPF_JPsiGammaBranches();
  void fillPF_JPsiGammaBranches(Int_t nDimu,Int_t phoIdx, Float_t mmg_mass, GenMatchInfo *grnm=nullptr);

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
 
  void  compressAllFloatStorage();
  // master switches
  Bool_t         isMC;
  Bool_t         isRECO;
  Bool_t         isMiniAOD;
  Bool_t         Run2_2018_;
  
  // switches
  Bool_t doGenParticles_;
  Bool_t doMuons_;
  Bool_t doDimuons_;
  Bool_t doMuMuK_;
  Bool_t doJPsiK_;
  Bool_t doPsi2SK_;
  Bool_t doMuMuKK_;
  Bool_t doMuMuGamma;
  Bool_t doJPsiGamma;
  Bool_t doPhotons_;
  Bool_t doPhotonID_;
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
  Int_t nCountErr;
  // ----------

  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;
  const reco::BeamSpot* beamSpot;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::GenParticleCollection>     genParticlesCollection_;
  edm::EDGetTokenT<edm::View<reco::Photon>>       gedPhotonsCollection_;
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
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent_token_;
  edm::EDGetTokenT<reco::PFCandidateCollection>     pfCandidateCollection_;
  edm::Handle<reco::PFCandidateCollection> pfCandidateHandle;
  
  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  

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

  // Photon MVA

  std::string mvaValMapTag_;
  edm::EDGetTokenT<edm::ValueMap<float>> valMapToken_;


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
  Float_t ptMinKaon_;
  Float_t etaMaxKaon_;
  Float_t minBKmmMass_;
  Float_t maxBKmmMass_;

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
  ULong64_t      nPU;
  UInt_t         lumis_;
  Bool_t         isData_;
  Float_t        rho_;


  // # BeamSpot # //
  Float_t beamspot_x_,beamspot_y_,beamspot_z_,beamspot_x_error_,beamspot_y_error_,beamspot_z_error_;
  Float_t beamspot_dxdz_,beamspot_dydz_,beamspot_sigmaZ_,beamspot_dxdz_error_,beamspot_dydz_error_,beamspot_sigmaZError_;
  Float_t beamspot_beamWidthX_,beamspot_beamWidthY_,beamspot_beamWidthX_error_,beamspot_beamWidthY_error_;
  Float_t beamspot_covXX,beamspot_covXY,beamspot_covXZ,beamspot_covYY,beamspot_covYZ,beamspot_covZZ;  
  
  // # Trigger #
  std::vector<std::string>  triggerFilter;
  std::vector<std::string>  trigTable;
  std::vector<std::string>  trigNames;
  std::vector<Bool_t>         trigResult;
  std::vector<int>          trigPrescales;
  
  std::string  *TrigTable_store;
  Bool_t         *TrigResult_store;
  Int_t          *TrigPrescales_store;
  
  std::vector<std::string>  l1Table;
  std::vector<int>          l1Prescales;

  void SetupTriggerStorageVectors();
  void SetupTriggerBranches();
  void FillTrggerBranches();
  void ClearTrggerStorages();

  // generic storage dictionary

  std::map<std::string,Int_t > storageMapInt;
  std::map<std::string,Int_t *> storageMapIntArray;
  std::map<std::string,std::string > storageMapFloatArrayLengthTags;
  std::map<std::string,Float_t * > storageMapFloatArray;
  std::map<std::string,Bool_t * > storageMapBoolArray;
  std::map<std::string,uint64_t * > storageMapUint64Array;


  // reco::GenParticle
  int  gen_nBs_, gen_nBsMuonM_, gen_nBsMuonP_ , gen_nBsPhoton_ ;
  std::vector<Float_t> gen_Bs_pt_,  gen_Bs_energy_,    gen_Bs_eta_,      gen_Bs_phi_,   gen_Bs_pz_,  gen_Bs_pdgId_;
  std::vector<Float_t> gen_BsMuonM_pt_, gen_BsMuonM_eta_, gen_BsMuonM_phi_;
  std::vector<Float_t> gen_BsMuonP_pt_, gen_BsMuonP_eta_, gen_BsMuonP_phi_;
  std::vector<Float_t> gen_BsPhoton_pt_, gen_BsPhoton_energy_, gen_BsPhoton_eta_, gen_BsPhoton_phi_;

};

#endif
