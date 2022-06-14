#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>
#include <algorithm>

#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/interface/XGBooster.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

// 
// BmmMuonIdProducer is designed for Bs/d->mumu analysis
//

using namespace std;
typedef reco::Candidate::LorentzVector LorentzVector;
typedef pair<const reco::MuonChamberMatch*, const reco::MuonSegmentMatch*> MatchPair;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class BmmMuonIdProducer : public edm::EDProducer {
    
public:
    
  explicit BmmMuonIdProducer(const edm::ParameterSet &iConfig);
    
  ~BmmMuonIdProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);
  void fillMatchInfo(pat::CompositeCandidate& cand, const reco::Muon& muon);
  void fillSoftMva(pat::CompositeCandidate& mu_cand);
  const l1t::Muon* getL1Muon( const reco::Candidate& cand);
  
  // ----------member data ---------------------------
    
  edm::EDGetTokenT<vector<reco::Muon> > muonToken_;
  bool isMC_;
  string triggerCollection_;
  vector<string> triggers_;
  XGBooster softMuonMva_;
  edm::Handle<BXVector<l1t::Muon> >   l1Handle_;
};

BmmMuonIdProducer::BmmMuonIdProducer(const edm::ParameterSet &iConfig):
muonToken_( consumes<vector<reco::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
softMuonMva_(iConfig.getParameter<edm::FileInPath>("softMuonMva").fullPath())
{
    produces<std::vector<float>>();
    
    vector<string> features = {"trkValidFrac", "glbTrackProbability", "nLostHitsInner",
      "nLostHitsOuter", "trkKink", "chi2LocalPosition", "match2_dX", "match2_pullX", "match1_dX", "match1_pullX",
      "nPixels", "nValidHits", "nLostHitsOn", "match2_dY", "match1_dY", "match2_pullY", "match1_pullY",
      "match2_pullDyDz", "match1_pullDyDz", "match2_pullDxDz", "match1_pullDxDz"};
    for (const auto& feature: features)
      softMuonMva_.addFeature(feature);
}

void BmmMuonIdProducer::fillSoftMva(pat::CompositeCandidate& mu_cand){
  
  softMuonMva_.set("trkValidFrac",        mu_cand.userFloat("trkValidFrac"));
  softMuonMva_.set("glbTrackProbability", mu_cand.userFloat("glbTrackProbability"));
  softMuonMva_.set("nLostHitsInner",      mu_cand.userInt(  "nLostHitsInner"));
  softMuonMva_.set("nLostHitsOuter",      mu_cand.userInt(  "nLostHitsOuter"));
  softMuonMva_.set("trkKink",             mu_cand.userFloat("trkKink"));
  softMuonMva_.set("chi2LocalPosition",   mu_cand.userFloat("chi2LocalPosition"));
  softMuonMva_.set("match2_dX",           mu_cand.userFloat("match2_dX"));
  softMuonMva_.set("match2_pullX",        mu_cand.userFloat("match2_pullX"));
  softMuonMva_.set("match1_dX",           mu_cand.userFloat("match1_dX"));
  softMuonMva_.set("match1_pullX",        mu_cand.userFloat("match1_pullX"));
  softMuonMva_.set("nPixels",             mu_cand.userInt(  "nPixels"));
  softMuonMva_.set("nValidHits",          mu_cand.userInt(  "nValidHits"));
  softMuonMva_.set("nLostHitsOn",         mu_cand.userInt(  "nLostHitsOn"));
  softMuonMva_.set("match2_dY",           mu_cand.userFloat("match2_dY"));
  softMuonMva_.set("match2_pullY",        mu_cand.userFloat("match2_pullY"));
  softMuonMva_.set("match1_dY",           mu_cand.userFloat("match1_dY"));
  softMuonMva_.set("match1_pullY",        mu_cand.userFloat("match1_pullY"));
  softMuonMva_.set("match2_pullDyDz",     mu_cand.userFloat("match2_pullDyDz"));
  softMuonMva_.set("match1_pullDyDz",     mu_cand.userFloat("match1_pullDyDz"));
  softMuonMva_.set("match2_pullDxDz",     mu_cand.userFloat("match2_pullDxDz"));
  softMuonMva_.set("match1_pullDxDz",     mu_cand.userFloat("match1_pullDxDz"));

  mu_cand.addUserFloat("newSoftMuonMva", softMuonMva_.predict());
}

const l1t::Muon* BmmMuonIdProducer::getL1Muon( const reco::Candidate& cand ){
  const l1t::Muon* match = nullptr;
  double best_dr = 999.;
  // Loop over L1 candidates from BX 0 only
  for (auto it = l1Handle_->begin(0); it != l1Handle_->end(0); it++){
    double dr = deltaR(*it, cand);
    if (match == nullptr or dr < best_dr){
      best_dr = dr;
      match = &*it;
    }
  }
  return match;
}


void BmmMuonIdProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<vector<reco::Muon> > muonHandle;
    iEvent.getByToken(muonToken_, muonHandle);

    // Output collection
    auto muonsMVA = make_unique<std::vector<float>>();
    muonsMVA->reserve(128);
    for ( const auto& muon: *muonHandle.product()){
      pat::CompositeCandidate mu_cand;
      mu_cand.addUserFloat("trkKink", muon.combinedQuality().trkKink);
      mu_cand.addUserFloat("glbTrackProbability", muon.combinedQuality().glbTrackProbability);
      mu_cand.addUserFloat("chi2LocalPosition", muon.combinedQuality().chi2LocalPosition);
      if (muon.isGlobalMuon())
	mu_cand.addUserFloat("glbNormChi2", muon.globalTrack()->normalizedChi2());
      else
	mu_cand.addUserFloat("glbNormChi2", 9999.);

      if (muon.isTrackerMuon() or muon.isGlobalMuon()){
	mu_cand.addUserFloat("trkValidFrac",  muon.innerTrack()->validFraction());
	
	mu_cand.addUserInt("nPixels",         muon.innerTrack()->hitPattern().numberOfValidPixelHits());
	mu_cand.addUserInt("nValidHits",      muon.innerTrack()->hitPattern().numberOfValidTrackerHits());
	mu_cand.addUserInt("nLostHitsInner",  muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS));
	mu_cand.addUserInt("nLostHitsOn",     muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::TRACK_HITS));
	mu_cand.addUserInt("nLostHitsOuter",  muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS));
	
	mu_cand.addUserInt("trkLayers",           muon.innerTrack()->hitPattern().trackerLayersWithMeasurement());
	mu_cand.addUserInt("trkLostLayersInner",  muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS));
	mu_cand.addUserInt("trkLostLayersOn",     muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));
	mu_cand.addUserInt("trkLostLayersOuter",  muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_OUTER_HITS));

	mu_cand.addUserInt("highPurity",   muon.innerTrack()->quality(reco::Track::highPurity));

      } else {
	mu_cand.addUserFloat("trkValidFrac",  0);
	
	mu_cand.addUserInt("nPixels",         0);
	mu_cand.addUserInt("nValidHits",      0);
	mu_cand.addUserInt("nLostHitsInner",  0);
	mu_cand.addUserInt("nLostHitsOn",     0);
	mu_cand.addUserInt("nLostHitsOuter", 0);
	
	mu_cand.addUserInt("trkLayers",           0);
	mu_cand.addUserInt("trkLostLayersInner",  0);
	mu_cand.addUserInt("trkLostLayersOn",     0);
	mu_cand.addUserInt("trkLostLayersOuter",  0);

	mu_cand.addUserInt("highPurity",   0);

      }
	
      fillMatchInfo(mu_cand, muon);
      fillSoftMva(mu_cand);
      
     // if (isMC_){
	 //    mu_cand.addUserInt("simType", muon.simType());
	 //    mu_cand.addUserInt("simExtType", muon.simExtType());
	 //    mu_cand.addUserInt("simPdgId", muon.simPdgId());
	 //    mu_cand.addUserInt("simMotherPdgId", muon.simMotherPdgId());
	 //    mu_cand.addUserFloat("simProdRho", muon.simProdRho());
	 //    mu_cand.addUserFloat("simProdZ", muon.simProdZ());
     // }

      muonsMVA->push_back(mu_cand.userFloat("newSoftMuonMva"));
    }
    
    iEvent.put(move(muonsMVA));
}

namespace {

        const MatchPair&
        getBetterMatch(const MatchPair& match1, const MatchPair& match2){
        
          // Prefer DT over CSC simply because it's closer to IP
          // and will have less multiple scattering (at least for
          // RB1 vs ME1/3 case). RB1 & ME1/2 overlap is tiny
          if (match2.first->detector() == MuonSubdetId::DT and
              match1.first->detector() != MuonSubdetId::DT)
            return match2;
        
          // For the rest compare local x match. We expect that
          // segments belong to the muon, so the difference in
          // local x is a reflection on how well we can measure it
          if ( abs(match1.first->x - match1.second->x) >
               abs(match2.first->x - match2.second->x) )
            return match2;
            
          return match1;
        }
        
        float dX(const MatchPair& match){
          if (match.first and match.second->hasPhi())
            return (match.first->x - match.second->x);
          else
            return 9999.;
        }
        
        float pullX(const MatchPair& match){
          if (match.first and match.second->hasPhi())
            return dX(match) /
              sqrt(pow(match.first->xErr, 2) + pow(match.second->xErr, 2));
          else
            return 9999.;
        }
        
        float pullDxDz(const MatchPair& match){
          if (match.first and match.second->hasPhi())
            return (match.first->dXdZ - match.second->dXdZ) /
                   sqrt(pow(match.first->dXdZErr, 2) + pow(match.second->dXdZErr, 2));
          else
            return 9999.;
        }
        
        float dY(const MatchPair& match){
          if (match.first and match.second->hasZed())
            return (match.first->y - match.second->y);
          else
            return 9999.;
        }
        
        float pullY(const MatchPair& match){
          if (match.first and match.second->hasZed())
            return dY(match) /
              sqrt(pow(match.first->yErr, 2) + pow(match.second->yErr, 2));
          else
            return 9999.;
        }
        
        float pullDyDz(const MatchPair& match){
          if (match.first and match.second->hasZed())
            return (match.first->dYdZ - match.second->dYdZ) /
                   sqrt(pow(match.first->dYdZErr, 2) + pow(match.second->dYdZErr, 2));
          else
            return 9999.;
        }
        
        void fillMatchInfoForStation(string prefix,
        			     pat::CompositeCandidate& cand,
        			     const MatchPair& match){
          cand.addUserFloat(prefix + "_dX",       dX(match));
          cand.addUserFloat(prefix + "_pullX",    pullX(match));
          cand.addUserFloat(prefix + "_pullDxDz", pullDxDz(match));
          cand.addUserFloat(prefix + "_dY",       dY(match));
          cand.addUserFloat(prefix + "_pullY",    pullY(match));
          cand.addUserFloat(prefix + "_pullDyDz", pullDyDz(match));
        }
        
}

void BmmMuonIdProducer::fillMatchInfo(pat::CompositeCandidate& cand,
        				      const reco::Muon& muon){
          // Initiate containter for results
          const int n_stations = 2;
          vector<MatchPair> matches;
          for (unsigned int i=0; i < n_stations; ++i)
            matches.push_back(pair(nullptr, nullptr));
        
          // Find best matches
          for (auto& chamberMatch : muon.matches()){
            unsigned int station = chamberMatch.station() - 1;
            if (station >= n_stations) continue;
        
            // Find best segment match.
            // We could consider all segments, but we will restrict to segments
            // that match to this candidate better than to other muon candidates
            for (auto& segmentMatch : chamberMatch.segmentMatches){
              if ( not segmentMatch.isMask(reco::MuonSegmentMatch::BestInStationByDR) ||
        	   not segmentMatch.isMask(reco::MuonSegmentMatch::BelongsToTrackByDR) )
        	continue;
        
              // Multiple segment matches are possible in different
              // chambers that are either overlapping or belong to
              // different detectors. We need to select one.
              auto match_pair = MatchPair(&chamberMatch, &segmentMatch);
              
              if (matches[station].first)
        	matches[station] = getBetterMatch(matches[station], match_pair);
              else
        	matches[station] = match_pair;
            }
          }
        
          // Fill matching information
          fillMatchInfoForStation("match1", cand, matches[0]);
          fillMatchInfoForStation("match2", cand, matches[1]);
        }



DEFINE_FWK_MODULE(BmmMuonIdProducer);

//  LocalWords:  vertices
