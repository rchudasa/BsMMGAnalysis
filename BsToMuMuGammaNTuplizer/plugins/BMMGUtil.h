#ifndef __BMMGUtil_H__
#define __BMMGUtil_H__


bool BsToMuMuGammaNTuplizer::isGoodMuon(const pat::Muon& muon){
  if ( not muon.isLooseMuon() ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  if ( muon.pt() < pTMinMuons || fabs(muon.eta()) > etaMax_muon ) return false;
  return true;
}


float
BsToMuMuGammaNTuplizer::computeTrkMuonIsolation(const pat::Muon& the_muon, const pat::Muon& the_other_muon, 
					  unsigned int primaryVertexIndex,
					  float minPt, float dR,
					  std::vector<const reco::PFCandidate*> ignoreTracks)
{
  float sumPt(0);
  const reco::VertexCollection& vertices = *pvHandle_.product();
  if( primaryVertexIndex >= vertices.size()) return 0.0;
  const auto & vertex = vertices.at(primaryVertexIndex);
  Float_t vx(vertex.x());
  Float_t vy(vertex.y());
  Float_t vz(vertex.z());
  for (const auto& pfCand: *pfCandidateHandle.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(the_muon, pfCand) < 0.01 || deltaR(the_other_muon, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (not pfCand.trackRef()) continue;
    if (pfCand.pt()<minPt) continue;
    //if (pfCand.vertexRef().key()!=primaryVertexIndex) continue;
    if ((abs(pfCand.vz() - vz) < 0.01 ) and  (abs(pfCand.vy() - vy) < 0.01 ) and ( abs(pfCand.vx() - vx) < 0.01 ) ) continue;
    if (deltaR(the_muon, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }

  return the_muon.pt()/(the_muon.pt()+sumPt);
}


float
BsToMuMuGammaNTuplizer::otherVertexMaxProb(const pat::Muon& muon1, 
				     const pat::Muon& muon2,
				     float minPt,
				     float max_doca,
				     std::vector<const reco::PFCandidate*> ignoreTracks){
  float bestMu1Vtx = 0;
  float bestMu2Vtx = 0;
  KalmanVertexFitter kvf;
  std::vector<reco::TransientTrack> transTrksForMu1Vertex;
  transTrksForMu1Vertex.push_back((*theTTBuilder_).build(muon1.innerTrack().get()));
  std::vector<reco::TransientTrack> transTrksForMu2Vertex;
  transTrksForMu2Vertex.push_back((*theTTBuilder_).build(muon2.innerTrack().get()));


  for (const auto& pfCand: *pfCandidateHandle.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    //if (!pfCand.hasTrackDetails()) continue;
    if ( not pfCand.trackRef() ) continue;
    if (pfCand.pt()<minPt) continue;
    double mu1_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
						pfCand.bestTrack());
    double mu2_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
						pfCand.bestTrack());
    if (mu1_doca < max_doca and mu1_doca < mu2_doca){
      // first  muon is closer - check vertex probability
      transTrksForMu1Vertex.push_back((*theTTBuilder_).build(pfCand.bestTrack()));
      TransientVertex tv = kvf.vertex(transTrksForMu1Vertex);
      if ( tv.isValid() ){
	float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
	if (vtxProb > bestMu1Vtx) bestMu1Vtx = vtxProb;
      }
      transTrksForMu1Vertex.pop_back();
    }
    if (mu2_doca < max_doca and mu2_doca < mu1_doca){
      // second  muon is closer - check vertex probability
      transTrksForMu2Vertex.push_back((*theTTBuilder_).build(pfCand.bestTrack()));
      TransientVertex tv = kvf.vertex(transTrksForMu2Vertex);
      if ( tv.isValid() ){
	float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
	if (vtxProb > bestMu2Vtx) bestMu2Vtx = vtxProb;
      }
      transTrksForMu2Vertex.pop_back();
    }
  }
  return max(bestMu1Vtx,bestMu2Vtx);
}

float BsToMuMuGammaNTuplizer::distanceOfClosestApproach( const reco::Track* track1,
						   const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder_->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder_->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

Measurement1D 
BsToMuMuGammaNTuplizer::distanceOfClosestApproach( const reco::Track* track,
					     RefCountedKinematicVertex vertex)
{
  if (not vertex->vertexIsValid()) return Measurement1D(-1.0,-1.0);
  VertexDistance3D distance3D;
  const reco::TransientTrack tt = theTTBuilder_->build(track);
  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), vertex->position());
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex->vertexState());
  return doca;
}

Measurement1D 
BsToMuMuGammaNTuplizer::distanceOfClosestApproach( const reco::Track* track,
					     const reco::Vertex& vertex)
{
  VertexDistance3D distance3D;
  const reco::TransientTrack tt = theTTBuilder_->build(track);
  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), GlobalPoint(Basic3DVector<float>(vertex.position())));
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex);
  return doca;
}

float
BsToMuMuGammaNTuplizer::computeTrkMuMuIsolation(const pat::Muon& muon1, const pat::Muon& muon2, 
					  unsigned int primaryVertexIndex,
					  float minPt, float dR,
					  std::vector<const reco::PFCandidate*> ignoreTracks)
{
  float sumPt(0);
  auto b_p4 = muon1.p4()+muon2.p4();
  const reco::VertexCollection& vertices = *pvHandle_.product();
  if( primaryVertexIndex >= vertices.size()) return 0.0;
  const auto & vertex = vertices.at(primaryVertexIndex);
  Float_t vx(vertex.x());
  Float_t vy(vertex.y());
  Float_t vz(vertex.z());
  for (const auto& pfCand: *pfCandidateHandle.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if ( not pfCand.trackRef()) continue;
    if (pfCand.pt()<minPt) continue;
    if (deltaR(b_p4, pfCand) > dR) continue;
    if ((abs(pfCand.vz() - vz) < 0.01 ) and  (abs(pfCand.vy() - vy) < 0.01 ) and ( abs(pfCand.vx() - vx) < 0.01 ) ) continue;
    sumPt += pfCand.pt();
  }

  return b_p4.pt()/(b_p4.pt()+sumPt);
}

DisplacementInformationIn3D BsToMuMuGammaNTuplizer::compute3dDisplacement(const KinematicFitResult& fit,
								    const reco::VertexCollection& vertices,
								    bool closestIn3D)
{
  DisplacementInformationIn3D result;
  if (not fit.valid()) return result;

  // Potential issue: tracks used to build the candidate could 
  // also be used in the primary vertex fit. One can refit the vertices
  // excluding tracks from the cadndidate. It's not done at the moment 
  // due to non-trivial linkig between primary vertex and its tracks 
  // in MiniAOD. Also not all muons associated to a vertex are really 
  // used in the fit, so the potential bias most likely small.
  
  auto candTransientTrack = fit.refitMother->refittedTransientTrack();

  const reco::Vertex* bestVertex(0);
  int bestVertexIndex(-1);
  double minDistance(999.);

  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    if (closestIn3D){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
      if (impactParameter3D.first and impactParameter3D.second.value() < minDistance){
	minDistance = impactParameter3D.second.value();
	bestVertex = &vertex;
	bestVertexIndex = i;
      }
    } else{
      auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), vertex);
      double distance = fabs(impactParameterZ.second.value());
      if (impactParameterZ.first and distance < minDistance){
	minDistance = distance;
	bestVertex = &vertex;
	bestVertexIndex = i;
      }
    }
  }

  // find second best vertex
  const reco::Vertex* bestVertex2(0);
  int bestVertexIndex2(-1);
  double minDistance2(999.);
  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    if (closestIn3D){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
      if (impactParameter3D.first and impactParameter3D.second.value() < minDistance2 and impactParameter3D.second.value() > minDistance){
	minDistance2 = impactParameter3D.second.value();
	bestVertex2 = &vertex;
	bestVertexIndex2 = i;
      }
    } else{
      auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), vertex);
      double distance = fabs(impactParameterZ.second.value());
      if (impactParameterZ.first and distance < minDistance2 and distance > minDistance){
	minDistance2 = distance;
	bestVertex2 = &vertex;
	bestVertexIndex2 = i;
      }
    }
  }

  if (! bestVertex) return result;

  auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex);
  auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex);
  result.pv = bestVertex;
  result.pvIndex = bestVertexIndex;
  if (impactParameterZ.first) {
    result.longitudinalImpactParameter    = impactParameterZ.second.value();
    result.longitudinalImpactParameterSig = impactParameterZ.second.significance();
    result.longitudinalImpactParameterErr = impactParameterZ.second.error();
  }
  if (impactParameter3D.first and not isnan(impactParameter3D.second.error())) {
    result.distaceOfClosestApproach       = impactParameter3D.second.value();
    result.distaceOfClosestApproachSig    = impactParameter3D.second.significance();
    result.distaceOfClosestApproachErr    = impactParameter3D.second.error();
  }

  // compute decay length
  VertexDistance3D distance3D;
  auto dist = distance3D.distance(*bestVertex, fit.refitVertex->vertexState() );
  result.decayLength    = dist.value();
  result.decayLengthErr = dist.error();
  
  VertexDistanceXY distanceXY;
  auto distXY = distanceXY.distance(*bestVertex, fit.refitVertex->vertexState() );

  if (bestVertex2){
    auto impactParameter3D2 = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex2);
    auto impactParameterZ2  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex2);
    result.pv2 = bestVertex2;
    result.pv2Index = bestVertexIndex2;
    if (impactParameterZ2.first) {
      result.longitudinalImpactParameter2    = impactParameterZ2.second.value();
      result.longitudinalImpactParameter2Sig = impactParameterZ2.second.significance();
      result.longitudinalImpactParameter2Err = impactParameterZ2.second.error();
    }
    if (impactParameter3D2.first) {
      result.distaceOfClosestApproach2       = impactParameter3D2.second.value();
      result.distaceOfClosestApproach2Sig    = impactParameter3D2.second.value();
      result.distaceOfClosestApproach2Err    = impactParameter3D2.second.error();
    }

    // compute decay length
    VertexDistance3D distance3D;
    auto dist = distance3D.distance(*bestVertex2, fit.refitVertex->vertexState() );
    result.decayLength2    = dist.value();
    result.decayLength2Err = dist.error();

  }

  //
  // Pointing angle
  //
  auto alpha = getAlpha(fit.refitVertex->vertexState().position(),
			fit.refitVertex->vertexState().error(),
			GlobalPoint(Basic3DVector<float>(bestVertex->position())),
			GlobalError(bestVertex->covariance()),
			fit.refitMother->currentState().globalMomentum());

  auto alphaXY = getAlpha(fit.refitVertex->vertexState().position(),
			  fit.refitVertex->vertexState().error(),
			  GlobalPoint(Basic3DVector<float>(bestVertex->position())),
			  GlobalError(bestVertex->covariance()),
			  fit.refitMother->currentState().globalMomentum(),
			  true);

  result.alpha    = alpha.first;
  result.alphaErr = alpha.second;

  result.alphaXY    = alphaXY.first;
  result.alphaXYErr = alphaXY.second;

  
  //
  // Decay time information
  //
  TVector3 plab(fit.refitMother->currentState().globalMomentum().x(),
		fit.refitMother->currentState().globalMomentum().y(),
                fit.refitMother->currentState().globalMomentum().z());
  const double massOverC = fit.mass()/TMath::Ccgs();

  // get covariance matrix for error propagation in decayTime calculation
  auto vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(bestVertex->error()),
					     fit.refitMother->currentState().kinematicParametersError().matrix());
  auto vtxDistanceJac3d = makeJacobianVector3d(bestVertex->position(), fit.refitVertex->vertexState().position(), plab);
  auto vtxDistanceJac2d = makeJacobianVector2d(bestVertex->position(), fit.refitVertex->vertexState().position(), plab);

  result.decayTime = dist.value() / plab.Mag() * cos(result.alpha) * massOverC;
  result.decayTimeError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac3d)) * massOverC;

  result.decayTimeXY = distXY.value() / plab.Perp() * cos(result.alphaXY) * massOverC;
  result.decayTimeXYError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac2d)) * massOverC;
    
  return result;
}

KinematicFitResult BsToMuMuGammaNTuplizer::vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
					    std::vector<float> masses)
{
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit
  if ( trks.size() != masses.size() ) 
    throw cms::Exception("Error") << "number of tracks and number of masses should match";

  std::vector<reco::TransientTrack> transTrks;

  KinematicParticleFactoryFromTransientTrack factory;
  KinematicParticleVertexFitter fitter;
    
  std::vector<RefCountedKinematicParticle> particles;

  double chi = 0.;
  double ndf = 0.;
  float muonMassErr(MuonMassErr_);
  for (unsigned int i=0; i<trks.size(); ++i){
    transTrks.push_back((*theTTBuilder_).build(trks[i]));
    particles.push_back(factory.particle(transTrks.back(),masses[i],chi,ndf,muonMassErr));
  }

  RefCountedKinematicTree vertexFitTree;
  KinematicFitResult result;

  try {
    vertexFitTree = fitter.fit(particles);
  } catch (const std::exception& e) {
    return result;
  }

  if ( !vertexFitTree->isValid()) return result;
  
  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;
  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();
  
  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}


KinematicFitResult
BsToMuMuGammaNTuplizer::vertexMuonsWithKinematicFitter(const pat::Muon& muon1,
						 const pat::Muon& muon2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( muon1.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( muon2.innerTrack().get() );
  masses.push_back(MuonMass_);
  return vertexWithKinematicFitter(trks,masses);
}

KinematicFitResult
BsToMuMuGammaNTuplizer::vertexWithKinematicFitter(const pat::Muon& muon1,
					    const pat::Muon& muon2,
					    const reco::PFCandidate& pion)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( muon1.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( muon2.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( pion.bestTrack() );
  masses.push_back(pionMass_);
  return vertexWithKinematicFitter(trks,masses);
}

CloseTrackInfo 
BsToMuMuGammaNTuplizer::findTracksCompatibleWithTheVertex(const reco::Muon& muon1,
						    const reco::Muon& muon2,
						    const KinematicFitResult& fit, 
						    double maxDoca,
						    std::vector<const reco::PFCandidate*> ignoreTracks)
{
  CloseTrackInfo result;
  result.pvHandle_ = pvHandle_;
  if (not fit.valid()) { 
    return result;
    }

  for (const auto& pfCand: *pfCandidateHandle.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (deltaR(*trk, pfCand) < 0.01){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    
    if (deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (not pfCand.trackRef()) continue;

    double mu1_kaon_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
						     pfCand.bestTrack());
    double mu2_kaon_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
						     pfCand.bestTrack());
    if (mu1_kaon_doca>maxDoca or mu2_kaon_doca>maxDoca) continue;
    
    CloseTrack track;
    track.pfCand = &pfCand;
    auto doca = distanceOfClosestApproach(pfCand.bestTrack(),fit.refitVertex);
    track.svDoca = doca.value();
    track.svDocaErr = doca.error();

    // add PV doca
    // if (not  pfCand.v0Ref()){
      //doca = distanceOfClosestApproach(pfCand.bestTrack(),pvHandle_->at(pfCand.vertexRef().key()) );a
      auto DrMin=1e9;
      auto pv_idx=-1;
      auto idx=-1;
        for(auto&  aVertex : *pvHandle_){
            idx++;
            if( not aVertex.isValid() ) continue;
            auto dr= sqrt( (pfCand.vx()-aVertex.x())*(pfCand.vx()-aVertex.x())  
                       + (pfCand.vy()-aVertex.y())*(pfCand.vy()-aVertex.y())
                       + (pfCand.vz()-aVertex.z())*(pfCand.vz()-aVertex.z())
                       );
            if(dr < DrMin)
            {
                DrMin=dr;
                pv_idx=idx;
            }
           }
      if( DrMin > 0.1){
         return result;
      }
      
      doca = distanceOfClosestApproach(pfCand.bestTrack(),pvHandle_->at(pv_idx) );
      track.pvDoca = doca.value();
      track.pvDocaErr = doca.error();
   // }
    
    auto fit_result = vertexWithKinematicFitter(muon1, muon2, pfCand);
    if (fit_result.valid()){
      track.svProb = fit_result.vtxProb();
      track.impactParameterSignificanceBS = pfCand.bestTrack()->dxyError()>0 ? fabs(pfCand.bestTrack()->dxy(*beamSpotHandle))/pfCand.bestTrack()->dxyError():0.0;
    }
    result.tracks.push_back(track);
  }

  return result;
}



KinematicFitResult
BsToMuMuGammaNTuplizer::fitBToKMuMu( RefCountedKinematicTree tree,
				  const reco::PFCandidate& kaon,
				  float mass_constraint)
{
  KinematicFitResult result; 
  if ( !tree->isValid()) return result;

  KinematicConstraint* mc(0);
  if (mass_constraint > 0){
    ParticleMass mass = mass_constraint;
    // mass constraint fit
    KinematicParticleFitter csFitter;
    float mass_sigma = JPsiMassErr_;
    // FIXME: memory leak
    mc = new MassKinematicConstraint(mass, mass_sigma);
    try {
      tree = csFitter.fit(mc, tree);
    } catch (const std::exception& e) {
      return result;
    }
  }

  const reco::TransientTrack kaonTT = theTTBuilder_->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  tree->movePointerToTheTop();
  BToKMuMuParticles.push_back(tree->currentParticle());
  float kaonMassErr(KaonMassErr_);
  BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,kaonMassErr));

  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(BToKMuMuParticles);
  } catch (const std::exception& e) {
    return result;
  }

  if ( !vertexFitTree->isValid()) return result;

  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;

  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();

  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}






KinematicFitResult BsToMuMuGammaNTuplizer::fitBToKJPsiMuMu( RefCountedKinematicParticle refitMuMu,
				   const reco::PFCandidate &kaon,
				   bool applyJpsiMassConstraint)
{
  const reco::TransientTrack mmTT = refitMuMu->refittedTransientTrack();
  const reco::TransientTrack kaonTT = theTTBuilder_->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  float MuMu_mass = refitMuMu->currentState().mass();
  float MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));

  if ( applyJpsiMassConstraint ){
    MuMu_mass = JPsiMass_;
    MuMu_mass_err = JPsiMassErr_;
  }

  BToKMuMuParticles.push_back(partFactory.particle(mmTT,MuMu_mass,chi,ndf,MuMu_mass_err));
  float kaonMassErr(KaonMassErr_);
  BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,kaonMassErr));

  KinematicFitResult result; 
  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(BToKMuMuParticles);
  } catch (const std::exception& e) {
    return result;
  }

  if ( !vertexFitTree->isValid()) return result;

  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;

  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();

  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}






Float_t dX(const MatchPair& match){
  if (match.first and match.second->hasPhi())
    return (match.first->x - match.second->x);
  else
    return 9999.;
}

Float_t pullX(const MatchPair& match){
  if (match.first and match.second->hasPhi())
    return dX(match) /
      sqrt(pow(match.first->xErr, 2) + pow(match.second->xErr, 2));
  else
    return 9999.;
}

Float_t pullDxDz(const MatchPair& match){
  if (match.first and match.second->hasPhi())
    return (match.first->dXdZ - match.second->dXdZ) /
           sqrt(pow(match.first->dXdZErr, 2) + pow(match.second->dXdZErr, 2));
  else
    return 9999.;
}

Float_t dY(const MatchPair& match){
  if (match.first and match.second->hasZed())
    return (match.first->y - match.second->y);
  else
    return 9999.;
}

Float_t pullY(const MatchPair& match){
  if (match.first and match.second->hasZed())
    return dY(match) /
      sqrt(pow(match.first->yErr, 2) + pow(match.second->yErr, 2));
  else
    return 9999.;
}

Float_t pullDyDz(const MatchPair& match){
  if (match.first and match.second->hasZed())
    return (match.first->dYdZ - match.second->dYdZ) /
           sqrt(pow(match.first->dYdZErr, 2) + pow(match.second->dYdZErr, 2));
  else
    return 9999.;
}

void BsToMuMuGammaNTuplizer::fillMatchInfoForStation(std::string prefix,
			     Int_t idxToFill,
			     const MatchPair& match){
  storageMapFloatArray["muon_"+prefix + "_dX"      ][idxToFill] =  dX(match);
  storageMapFloatArray["muon_"+prefix + "_pullX"   ][idxToFill] =  pullX(match);
  storageMapFloatArray["muon_"+prefix + "_pullDxDz"][idxToFill] =  pullDxDz(match);
  storageMapFloatArray["muon_"+prefix + "_dY"      ][idxToFill] =  dY(match);
  storageMapFloatArray["muon_"+prefix + "_pullY"   ][idxToFill] =  pullY(match);
  storageMapFloatArray["muon_"+prefix + "_pullDyDz"][idxToFill] =  pullDyDz(match);
}

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


void BsToMuMuGammaNTuplizer::fillMatchInfo(Int_t idxToFill, const reco::Muon &muon){
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
  fillMatchInfoForStation("match1", idxToFill, matches[0]);
  fillMatchInfoForStation("match2", idxToFill, matches[1]);
}


Float_t BsToMuMuGammaNTuplizer::reduceFloat(Float_t val, int bits)
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

Float_t BsToMuMuGammaNTuplizer::getMinEnergyHCAL_(HcalDetId id) const {
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

void BsToMuMuGammaNTuplizer::compressAllFloatStorage()
{

    for(std::map<string,Float_t*>::iterator it=storageMapFloatArray.begin() ; it!=storageMapFloatArray.end(); ++it)
    {
        auto Name=it->first;
        auto array = it->second; 
        for(int i=0;i<N_COMPRESS_MAX;i++)
            array[i]=reduceFloat(array[i],nBits_);
   }     
}

#endif
