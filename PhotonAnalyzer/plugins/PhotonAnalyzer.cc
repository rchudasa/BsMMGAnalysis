// -*- C++ -*-
//
// Package:    BsMMGAnalysis/PhotonAnalyzer
// Class:      PhotonAnalyzer
//
/**\class PhotonAnalyzer PhotonAnalyzer.cc BsMMGAnalysis/PhotonAnalyzer/plugins/PhotonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ruchi Chudasama
//         Created:  Tue, 04 May 2021 06:26:56 GMT
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

#include <TTree.h>


class PhotonAnalyzer : public edm::EDAnalyzer {

public:
   explicit PhotonAnalyzer(const edm::ParameterSet&);
   ~PhotonAnalyzer();
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;

   void fillGenParticles (const edm::Event&);
   void fillPhotons      (const edm::Event&, const edm::EventSetup&);
   void fillPFPhotons    (const edm::Event&, const edm::EventSetup&);
   void fillSC           (const edm::Event&);

     // switches
   bool doGenParticles_;
   bool doPhotons_;
   bool doPFPhotons_;
   bool doSuperClusters_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::GenParticle>>      genParticlesCollection_;
  edm::EDGetTokenT<std::vector<reco::Photon>>         gedPhotonsCollection_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>>    pfPhotonsCollection_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>>   MustacheSCBarrelCollection_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>>   MustacheSCEndcapCollection_;

  // input tags 
  edm::InputTag genParticleSrc_;
  edm::InputTag gedPhotonSrc_;
  edm::InputTag pfPhotonSrc_;
  edm::InputTag MustacheSCBarrelSrc_;
  edm::InputTag MustacheSCEndcapSrc_;

   TTree*         tree_;

   // variables associated with tree branches
   UInt_t         run_;
   ULong64_t      event_;
   UInt_t         lumis_;
   Bool_t         isData_;
   Float_t        rho_;

   // reco::GenParticle
   Int_t          nMC_;
   std::vector<int>    mcPID_;
   std::vector<int>    mcStatus_;
   std::vector<float>  mcPt_;
   std::vector<float>  mcEta_;
   std::vector<float>  mcPhi_;
   std::vector<float>  mcE_;
   std::vector<float>  mcEt_;
   std::vector<float>  mcMass_;

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
   std::vector<float> scRawE_;
   std::vector<float> scEta_;
   std::vector<float> scPhi_;

};

PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& ps) :

  doGenParticles_(ps.getParameter<bool>("doGenParticles")),
  doPFPhotons_(ps.getParameter<bool>("doPFPhotons")),
  doSuperClusters_(ps.getParameter<bool>("doSuperClusters")),

  genParticleSrc_(ps.getUntrackedParameter<edm::InputTag>("genParticleSrc")),
  gedPhotonSrc_(ps.getUntrackedParameter<edm::InputTag>("gedPhotonSrc")),
  pfPhotonSrc_(ps.getUntrackedParameter<edm::InputTag>("pfPhotonSrc")),
  MustacheSCBarrelSrc_(ps.getParameter<edm::InputTag>("MustacheSCBarrelSrc")),
  MustacheSCEndcapSrc_(ps.getParameter<edm::InputTag>("MustacheSCEndcapSrc"))
{
  
  //now do what ever initialization is needed
   if(doGenParticles_) genParticlesCollection_   = consumes<edm::View<reco::GenParticle>>(genParticleSrc_);
   if(doPhotons_)    gedPhotonsCollection_       = consumes<std::vector<reco::Photon>>(gedPhotonSrc_);
   if(doPFPhotons_)  pfPhotonsCollection_        = consumes<std::vector<reco::PFCandidate>>(pfPhotonSrc_);
   if(doSuperClusters_){
	MustacheSCBarrelCollection_             = consumes<std::vector<reco::SuperCluster>>(MustacheSCBarrelSrc_);
   	MustacheSCEndcapCollection_             = consumes<std::vector<reco::SuperCluster>>(MustacheSCEndcapSrc_);
   }
  
  
  // initialize output TTree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("EventTree", "Event data");
  
  tree_->Branch("run",    &run_);
  tree_->Branch("event",  &event_);
  tree_->Branch("lumis",  &lumis_);
  tree_->Branch("isData", &isData_);
  
  if (doGenParticles_) {
    tree_->Branch("nMC",          &nMC_);
    tree_->Branch("mcPID",        &mcPID_);
    tree_->Branch("mcStatus",     &mcStatus_);
    tree_->Branch("mcPt",         &mcPt_);
    tree_->Branch("mcEta",        &mcEta_);
    tree_->Branch("mcPhi",        &mcPhi_);
    tree_->Branch("mcE",          &mcE_);
    tree_->Branch("mcEt",         &mcEt_);
    tree_->Branch("mcMass",       &mcMass_);  
  }
  
  if (doPhotons_) {
    tree_->Branch("nPho",                  &nPho_);
    tree_->Branch("phoE",                  &phoE_);
    tree_->Branch("phoEt",                 &phoEt_);
    tree_->Branch("phoEta",                &phoEta_);
    tree_->Branch("phoPhi",                &phoPhi_);
  }

  if (doPFPhotons_) {
    tree_->Branch("nPFPho",                  &nPFPho_);
    tree_->Branch("phoPFE",                  &phoPFE_);
    tree_->Branch("phoPFEt",                 &phoPFEt_);
    tree_->Branch("phoPFEta",                &phoPFEta_);
    tree_->Branch("phoPFPhi",                &phoPFPhi_);
  }

  if (doSuperClusters_) {
    tree_->Branch("nSC",                  &nSC_);
    tree_->Branch("scE",                  &scE_);
    tree_->Branch("scRawE",               &scRawE_);
    tree_->Branch("scEta",                &scEta_);
    tree_->Branch("scPhi",                &scPhi_);
  }
}

PhotonAnalyzer::~PhotonAnalyzer(){
}


// ------------ method called for each event  ------------
void 
PhotonAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es) {
  using namespace edm;
  
  if (doGenParticles_) {
    nMC_ = 0;
    mcPID_                .clear();
    mcStatus_             .clear();
    mcPt_                 .clear();
    mcEta_                .clear();
    mcPhi_                .clear();
    mcE_                  .clear();
    mcEt_                 .clear();
    mcMass_               .clear();
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
  }

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();

  // MC truth
  if (doGenParticles_ && !isData_) {
    fillGenParticles(e);
  }
  
  if (doPhotons_)    fillPhotons(e, es);
  if (doPFPhotons_) fillPFPhotons(e, es);
  if (doSuperClusters_) fillSC(e);
  tree_->Fill();
  
}



void PhotonAnalyzer::fillGenParticles(const edm::Event& e)
{
  // Fills tree branches with generated particle info.
  
  //std::cout << " fill gen loop "<< std::endl;
  edm::Handle<edm::View<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);
  
  // loop over MC particles
  for (auto p = genParticlesHandle->begin(); p != genParticlesHandle->end(); ++p) {
    // skip all primary particles if not particle gun MC
    //if (!runOnParticleGun_ && !p->mother()) continue;
    
    //std::cout << " in gen 1st loop "<< std::endl;
    // stable particles with pT > 5 GeV
    bool isStableFast = (p->status() == 1 && p->pdgId()==22 && p->pt() > 0.0);
    
    // stable leptons
    bool isStableLepton = (p->status() == 1 && abs(p->pdgId()) >= 11 && abs(p->pdgId()) <= 16);
    
    // (unstable) Z, W, H, top, bottom
    bool isHeavy = (p->pdgId() == 23 || abs(p->pdgId()) == 24 || p->pdgId() == 25 ||
		    abs(p->pdgId()) == 6 || abs(p->pdgId()) == 5);
    
    // reduce size of output root file
    //if (!isStableFast && !isStableLepton && !isHeavy) continue;
    //std::cout << " in gen loop "<< std::endl;
    mcPID_   .push_back(p->pdgId());
    mcStatus_.push_back(p->status());
    mcPt_    .push_back(p->pt());
    mcEta_   .push_back(p->eta());
    mcPhi_   .push_back(p->phi());
    mcE_     .push_back(p->energy());
    mcEt_    .push_back(p->et());
    mcMass_  .push_back(p->mass());
    nMC_++;
  } // gen-level particles loop
  
}

void PhotonAnalyzer::fillPhotons(const edm::Event& e, const edm::EventSetup& es)
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

void PhotonAnalyzer::fillPFPhotons(const edm::Event& e, const edm::EventSetup& es)
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


void PhotonAnalyzer::fillSC(edm::Event const& e) {
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
DEFINE_FWK_MODULE(PhotonAnalyzer);
