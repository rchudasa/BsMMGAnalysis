#ifndef __BMMGINTERFACE__H__
#define __BMMGINTERFACE__H__
//  Many of these function are obtained from the Bmm groups ntuplizer codebase
//  see here https://github.com/drkovalskyi/Bmm5
namespace {
  const float MuonMass_    = 0.10565837;
  const float MuonMassErr_ = 3.5*1e-9;
  const float KaonMass_    = 0.493677;
  const float KaonMassErr_ = 1.6e-5;
  const float pionMass_    = 0.139570;
  const float pionMassErr_ = 3.5e-7;
  const float JPsiMass_    = 3.0969;
  const float Psi2SMass_   = 3.6861;
  const float JPsiMassErr_ = 92.9e-6;
};



#include <TVector.h>
#include <TMatrix.h>
#include <algorithm>
#include <math.h>

typedef reco::Candidate::LorentzVector LorentzVector;

std::pair<float, float> getAlpha(const GlobalPoint& vtx_position, const GlobalError& vtx_error,
	 const GlobalPoint& ip_position,  const GlobalError& ip_error,
			     const GlobalVector &momentum,
			     bool transverse = true){
  AlgebraicSymMatrix33 error_matrix(vtx_error.matrix() + ip_error.matrix());
  GlobalVector dir(vtx_position - ip_position);
  if (dir.mag() == 0)
    return std::pair<float, float>(999., 999.);

  GlobalVector p(momentum);
  if (transverse){
    dir = GlobalVector(dir.x(), dir.y(), 0);
    p = GlobalVector(p.x(), p.y(), 0);
  }
  
  double dot_product = dir.dot(p);
  double cosAlpha = dot_product / p.mag() / dir.mag();
  if (cosAlpha > 1) cosAlpha = 1;
  if (cosAlpha < -1) cosAlpha = -1;

  // Error propagation

  double c1 = 1 / dir.mag() / p.mag();
  double c2 = dot_product / pow(dir.mag(), 3) / p.mag();
  
  double dfdx = p.x() * c1 - dir.x() * c2;
  double dfdy = p.y() * c1 - dir.y() * c2;
  double dfdz = p.z() * c1 - dir.z() * c2;

  double err2_cosAlpha =
    pow(dfdx, 2) * error_matrix(0, 0) +
    pow(dfdy, 2) * error_matrix(1, 1) +
    pow(dfdz, 2) * error_matrix(2, 2) +
    2 * dfdx * dfdy * error_matrix(0, 1) +
    2 * dfdx * dfdz * error_matrix(0, 2) +
    2 * dfdy * dfdz * error_matrix(1, 2);

  float err_alpha = fabs(cosAlpha) <= 1 and err2_cosAlpha >=0 ? sqrt(err2_cosAlpha) / sqrt(1-pow(cosAlpha, 2)) : 999;
  float alpha = acos(cosAlpha);
  if (isnan(alpha) or isnan(err_alpha))
    return std::pair<float, float>(999., 999.);
  else
    return std::pair<float, float>(alpha, err_alpha);
}


struct KinematicFitResult{
  bool treeIsValid;
  bool vertexIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
  float lxy, lxyErr, sigLxy, alphaBS, alphaBSErr;
  KinematicFitResult():treeIsValid(false),vertexIsValid(false),
		       lxy(-1.0), lxyErr(-1.0), sigLxy(-1.0),
		       alphaBS(-999.), alphaBSErr(-999.)
  {}

  bool valid() const {
    return treeIsValid and vertexIsValid;
  }

  void postprocess(const reco::BeamSpot& beamSpot)
  {
    if ( not valid() ) return;
    // displacement information
    TVector v(2);
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();

    TMatrix errVtx(2,2);
    errVtx(0,0) = refitVertex->error().cxx();
    errVtx(0,1) = refitVertex->error().matrix()(0,1);
    errVtx(1,0) = errVtx(0,1);
    errVtx(1,1) = refitVertex->error().cyy();

    TMatrix errBS(2,2);
    errBS(0,0) = beamSpot.covariance()(0,0);
    errBS(0,1) = beamSpot.covariance()(0,1);
    errBS(1,0) = beamSpot.covariance()(1,0);
    errBS(1,1) = beamSpot.covariance()(1,1);
    
    lxy = sqrt(v.Norm2Sqr());
    lxyErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / lxy;
    if (lxyErr > 0) sigLxy = lxy/lxyErr;
    
    // compute pointing angle wrt BeamSpot (2D)

    // rotatedCovariance3D - is a proper covariance matrix for Beam Spot,
    // which includes the beam spot width, not just uncertainty on the
    // absolute beamspot position
    auto alphaXY = getAlpha(refitVertex->vertexState().position(),
			    refitVertex->vertexState().error(),
			    GlobalPoint(Basic3DVector<float>(beamSpot.position())),
			    GlobalError(beamSpot.rotatedCovariance3D()),
			    refitMother->currentState().globalMomentum(),
			    true);
    alphaBS    = alphaXY.first;
    alphaBSErr = alphaXY.second;
  }
  
  float mass() const
  {
    if ( not valid() ) return -1.0;
    return refitMother->currentState().mass();
  }

  float refit_mass(unsigned int i, unsigned int j) const
  {
    if ( not valid() ) return -1.0;
    if (i >= refitDaughters.size()) return -2.0;
    if (j >= refitDaughters.size()) return -3.0;
    if (refitDaughters.at(i)->currentState().globalMomentum().mag2()<0) return -4.0;
    if (refitDaughters.at(j)->currentState().globalMomentum().mag2()<0) return -5.0;
    auto momentum = refitDaughters.at(i)->currentState().globalMomentum() + 
      refitDaughters.at(j)->currentState().globalMomentum();
    auto energy1 = sqrt(refitDaughters.at(i)->currentState().globalMomentum().mag2() + 
			pow(refitDaughters.at(i)->currentState().mass(),2));
    auto energy2 = sqrt(refitDaughters.at(j)->currentState().globalMomentum().mag2() + 
			pow(refitDaughters.at(j)->currentState().mass(),2));
    return sqrt(pow(energy1+energy2,2)-momentum.mag2());
  }

  GlobalVector p3() const
  {
    if ( not valid() ) return GlobalVector();
    return refitMother->currentState().globalMomentum();
  }

  GlobalVector dau_p3(unsigned int i) const
  {
    if ( not valid() or i>=refitDaughters.size() ) return GlobalVector();
    return refitDaughters.at(i)->currentState().globalMomentum();
  }

  float massErr() const
  {
    if ( not valid() ) return -1.0;
    return sqrt(refitMother->currentState().kinematicParametersError().matrix()(6,6));
  }

  float chi2() const
  {
    if ( not valid() ) return -1.0;
    return refitVertex->chiSquared();
  }

  float ndof() const
  {
    return refitVertex->degreesOfFreedom();
  }

  float vtxProb() const
  {
    if ( not valid() ) return -1.0;
    return TMath::Prob((double)refitVertex->chiSquared(), int(rint(refitVertex->degreesOfFreedom())));
  }
  
};

struct KalmanVertexFitResult{
  float vtxProb;
  bool  valid;
  std::vector<LorentzVector> refitVectors;
  GlobalPoint position;
  GlobalError err;
  float lxy, lxyErr, sigLxy;

  KalmanVertexFitResult():vtxProb(-1.0),valid(false),lxy(-1.0),lxyErr(-1.0),sigLxy(-1.0){}

  float mass() const
  {
    if (not valid) return -1.0;
    LorentzVector p4;
    for (auto v: refitVectors)
      p4 += v;
    return p4.mass();
  }
  
  void postprocess(const reco::BeamSpot& bs)
  {
    if (not valid) return;
    // position of the beam spot at a given z value (it takes into account the dxdz and dydz slopes)
    reco::BeamSpot::Point bs_at_z(bs.position(position.z()));
    GlobalPoint xy_displacement(position.x() - bs_at_z.x(),
				position.y() - bs_at_z.y(),
				0);
    lxy = xy_displacement.perp();
    lxyErr = sqrt(err.rerr(xy_displacement));
    if (lxyErr > 0) sigLxy = lxy/lxyErr;
  }
};

struct DisplacementInformationIn3D{
  double decayLength, decayLengthErr, decayLength2, decayLength2Err, 
    distaceOfClosestApproach, distaceOfClosestApproachErr, distaceOfClosestApproachSig,
    distaceOfClosestApproach2, distaceOfClosestApproach2Err, distaceOfClosestApproach2Sig,
    longitudinalImpactParameter, longitudinalImpactParameterErr, longitudinalImpactParameterSig,
    longitudinalImpactParameter2, longitudinalImpactParameter2Err,longitudinalImpactParameter2Sig,
    decayTime, decayTimeError, decayTimeXY, decayTimeXYError,
    alpha, alphaErr, alphaXY, alphaXYErr;
  const reco::Vertex *pv,*pv2;
  int pvIndex,pv2Index;
  DisplacementInformationIn3D():decayLength(-1.0), decayLengthErr(0.), decayLength2(-1.0), decayLength2Err(0.),
				distaceOfClosestApproach(-1.0), distaceOfClosestApproachErr(0.0), distaceOfClosestApproachSig(0.0),
				distaceOfClosestApproach2(-1.0), distaceOfClosestApproach2Err(0.0), distaceOfClosestApproach2Sig(0.0),
				longitudinalImpactParameter(0.0), longitudinalImpactParameterErr(0.), longitudinalImpactParameterSig(0.),
				longitudinalImpactParameter2(0.0), longitudinalImpactParameter2Err(0.), longitudinalImpactParameter2Sig(0.),
				decayTime(-999.), decayTimeError(-999.),
				decayTimeXY(-999.), decayTimeXYError(-999.),
				alpha(-999.), alphaErr(-999.),
				alphaXY(-999.), alphaXYErr(-999.),
				pv(0), pv2(0),
				pvIndex(-1), pv2Index(-1)
  {};
};

LorentzVector makeLorentzVectorFromPxPyPzM(double px, double py, double pz, double m){
  double p2 = px*px+py*py+pz*pz;
  return LorentzVector(px,py,pz,sqrt(p2+m*m));
}

struct GenMatchInfo{
  int mu1_pdgId, mu1_motherPdgId, mu2_pdgId, mu2_motherPdgId,
    kaon1_pdgId, kaon1_motherPdgId, kaon2_pdgId, kaon2_motherPdgId, photon_pdgId, photon_motherPdgId,
    mm_pdgId, mm_motherPdgId, kmm_pdgId, kkmm_pdgId, mmg_pdgId;
  float mu1_pt, mu2_pt, kaon1_pt, kaon2_pt, photon_pt, mm_mass, mm_pt, kmm_mass, kkmm_mass, mmg_mass,
    kmm_pt, kkmm_pt, mmg_pt;
  math::XYZPoint mm_prod_vtx, mm_vtx, kmm_prod_vtx, kkmm_prod_vtx, mmg_prod_vtx;
  const reco::Candidate* mc_mu1;
  const reco::Candidate* mc_mu2;
  const reco::Candidate* mc_kaon1;
  const reco::Candidate* mc_kaon2;
  const reco::Candidate* mc_photon;
  const reco::Candidate* match;
  const reco::Candidate* common_mother;
  GenMatchInfo():mu1_pdgId(0), mu1_motherPdgId(0), mu2_pdgId(0), mu2_motherPdgId(0), 
		 kaon1_pdgId(0), kaon1_motherPdgId(0), kaon2_pdgId(0), kaon2_motherPdgId(0),
		 photon_pdgId(0), photon_motherPdgId(0), mm_pdgId(0), mm_motherPdgId(0), 
		 kmm_pdgId(0), kkmm_pdgId(0), mmg_pdgId(0), mu1_pt(0), mu2_pt(0), 
		 kaon1_pt(0), kaon2_pt(0), photon_pt(0), mm_mass(0), mm_pt(0), 
		 kmm_mass(0), kkmm_mass(0), mmg_mass(0), kmm_pt(0), kkmm_pt(0), mmg_pt(0),
		 mc_mu1(0), mc_mu2(0), mc_kaon1(0), mc_kaon2(0), mc_photon(0),
		 match(0), common_mother(0)
  {}
  const reco::GenParticle* gen_mu1(){
    return dynamic_cast<const reco::GenParticle*>(mc_mu1);
  }
  const reco::GenParticle* gen_mu2(){
    return dynamic_cast<const reco::GenParticle*>(mc_mu2);
  }
    

};

struct GenEventInfo{};

struct CloseTrack{
  float svDoca, svDocaErr, svProb,
    pvDoca, pvDocaErr,
    impactParameterSignificanceBS;
  //const reco::PFCandidate* pfCand;
  const reco::PFCandidate* pfCand;
  CloseTrack(): svDoca(-1), svDocaErr(-1), svProb(-1),
		pvDoca(-1), pvDocaErr(-1),
		impactParameterSignificanceBS(-1),
		pfCand(0)
  {};
};


struct CloseTrackInfo{
  edm::Handle<reco::VertexCollection> pvHandle_;
  std::vector<CloseTrack> tracks;
  unsigned int nTracksByVertexProbability(double minProb = 0.1, 
					  double minIpSignificance = -1,
                      const reco::Vertex* vertex=nullptr,
					  const reco::PFCandidate* ignoreTrack1 = 0)
  {
  //std::cout<<" pIndex = "<<pvIndex<<"\n";
  //const reco::VertexCollection& vertices = *pvHandle_.product();
  //auto pvIndex_=-1;
  //if( pvIndex >=0 ) pvIndex_=pvIndex;
  //const auto & vertex = vertices.at(pvIndex_);
  Float_t vx(0.0);
  Float_t vy(0.0);
  Float_t vz(0.0);

  if(vertex)
  {
     vx =    vertex->x()   ;
     vy =    vertex->y()   ;
     vz =    vertex->z()   ;
  }
    unsigned int n = 0;
    for (auto track: tracks){
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if (minIpSignificance>0 and track.impactParameterSignificanceBS<minIpSignificance) continue;
      if (track.svProb<minProb) continue;
      if (vertex )
          if ((abs(track.pfCand->vz() - vz) < 0.01 ) and  (abs(track.pfCand->vy() - vy) < 0.01 ) and ( abs(track.pfCand->vx() - vx) < 0.01 ) ) continue;
      n++;
    }
    return n;
  }
  unsigned int nTracksByDisplacementSignificance(double max_svDoca = 0.03, 
						 double maxSignificance = -1,
                         const reco::Vertex* vertex=nullptr,
						 const reco::PFCandidate* ignoreTrack1 = 0)
  {
  // const reco::VertexCollection& vertices = *pvHandle_.product();
  //auto pvIndex_=-1;
  //if( pvIndex >=0 ) pvIndex_=pvIndex;
  //const auto & vertex = vertices.at(pvIndex_);

  Float_t vx(0.0);
  Float_t vy(0.0);
  Float_t vz(0.0);

  if(vertex)
  {
     vx =    vertex->x()   ;
     vy =    vertex->y()   ;
     vz =    vertex->z()   ;
  }
//   Float_t vx(vertex->x());
//   Float_t vy(vertex->y());
//   Float_t vz(vertex->z());
    unsigned int n = 0;
    for (auto track: tracks){
      if (track.svDoca>max_svDoca) continue;
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if (maxSignificance>0 and (track.svDocaErr<=0 or 
				 track.svDoca/track.svDocaErr > maxSignificance) ) continue;
      if ( vertex )
          if ((abs(track.pfCand->vz() - vz) < 0.01 ) and  (abs(track.pfCand->vy() - vy) < 0.01 ) and ( abs(track.pfCand->vx() - vx) < 0.01 ) ) continue;
      n++;
    }
    return n;
  }
  unsigned int nTracksByBetterMatch(double max_svDoca = 0.03, 
				    double maxSignificance = 2,
                    const reco::Vertex* vertex=nullptr,
				    const reco::PFCandidate* ignoreTrack1 = 0)
  {
  //const reco::VertexCollection& vertices = *pvHandle_.product();
  //auto pvIndex_=0;
  //if( pvIndex >=0 ) pvIndex_=pvIndex;
  //const auto & vertex = vertices.at(pvIndex_);
//  Float_t vx(vertex->x());
//  Float_t vy(vertex->y());
//  Float_t vz(vertex->z());
 
  Float_t vx(0.0);
  Float_t vy(0.0);
  Float_t vz(0.0);

  if(vertex)
  {
     vx =    vertex->x()  ; 
     vy =    vertex->y()  ;
     vz =    vertex->z()  ;
  } 

  unsigned int n = 0;
    for (auto track: tracks){
      if (track.svDoca>max_svDoca) continue;
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if (maxSignificance>0 and (track.svDocaErr<=0 or 
				 track.svDoca/track.svDocaErr > maxSignificance) ) continue;
      if (track.svDocaErr<=0 or (track.pvDocaErr>0 and track.svDoca/track.svDocaErr > track.pvDoca/track.pvDocaErr) ) continue;
      if (vertex )
          if ((abs(track.pfCand->vz() - vz) < 0.01 ) and  (abs(track.pfCand->vy() - vy) < 0.01 ) and ( abs(track.pfCand->vx() - vx) < 0.01 ) ) continue;
      n++;
    }
    return n;
  }
  float minDoca(double max_svDoca = 0.03, 
        const reco::Vertex* vertex=nullptr,
		const reco::PFCandidate* ignoreTrack1 = 0)
  {
  //const reco::VertexCollection& vertices = *pvHandle_.product();
  //auto pvIndex_=0;
  //if( pvIndex >=0 ) pvIndex_=pvIndex;
  //const auto & vertex = vertices.at(pvIndex_);
//  Float_t vx(vertex->x());
//  Float_t vy(vertex->y());
//  Float_t vz(vertex->z());

  Float_t vx(0.0);
  Float_t vy(0.0);
  Float_t vz(0.0);

  if(vertex)
  {
     vx =    vertex->x() ;   
     vy =    vertex->y() ;
     vz =    vertex->z() ;
  }

    float doca = 99.;
    for (auto track: tracks){
      if (track.svDoca>max_svDoca) continue;
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if ( vertex )
          if ((abs(track.pfCand->vz() - vz) < 0.01 ) and  (abs(track.pfCand->vy() - vy) < 0.01 ) and ( abs(track.pfCand->vx() - vx) < 0.01 ) ) continue;
      if (doca>track.svDoca) doca = track.svDoca;
    }
    return doca;
  }

  void fillCandInfo(pat::CompositeCandidate& cand, int pvIndex, std::string name)
  {
   // if (name!="") name += "_";
   // cand.addUserInt(   name + "nTrks",       nTracksByVertexProbability(0.1,-1.0,pvIndex) );
   // cand.addUserInt(   name + "nBMTrks",     nTracksByBetterMatch() );
   // cand.addUserInt(   name + "nDisTrks",    nTracksByVertexProbability(0.1, 2.0,pvIndex) );
   // cand.addUserInt(   name + "closetrk",    nTracksByDisplacementSignificance(0.03, -1, pvIndex) );
   // cand.addUserInt(   name + "closetrks1",  nTracksByDisplacementSignificance(0.03, 1, pvIndex) );
   // cand.addUserInt(   name + "closetrks2",  nTracksByDisplacementSignificance(0.03, 2, pvIndex) );
   // cand.addUserInt(   name + "closetrks3",  nTracksByDisplacementSignificance(0.03, 3, pvIndex) );
   // cand.addUserFloat( name + "docatrk",     minDoca(0.03, pvIndex) );
  }
};

struct BdtReaderData {
  float fls3d, alpha, pvips, iso, chi2dof, docatrk, closetrk, m1iso, m2iso, eta, m;
};

namespace {
  // Muon container to hold muons and hadrons that may decays to muons
  // Index is the position of the muon in the original muon collection
  // For hadrons Index is -1

  class MuonCand: public pat::Muon{
  public:
    MuonCand(const pat::Muon& muon, int index):
      pat::Muon(muon), index_(index), gen_(false)
    {
    }
    MuonCand(const reco::PFCandidate& hadron, bool from_gen):
      pat::Muon(reco::Muon(hadron.charge(), hadron.p4())),
      index_(-1), gen_(from_gen)
    {
      std::vector<reco::Track> tracks;
      assert(not hadron.trackRef());
      tracks.push_back(*hadron.bestTrack());
      setInnerTrack(reco::TrackRef(&tracks,0));
      embedTrack();
      setPdgId(hadron.pdgId());
    }
    int index() const { return index_; }
    bool from_gen() const { return gen_; }
  private:
    int index_;
    bool gen_;
  };
}

namespace {
  typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > cov33_t;
  typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6> > cov66_t;
  typedef ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > cov77_t;
  typedef ROOT::Math::SMatrix<double,9,9,ROOT::Math::MatRepSym<double,9> > cov99_t;
  typedef ROOT::Math::SVector<double,9> jac9_t;
  
  cov33_t GlobalError2SMatrix_33(GlobalError m_in) 
  {
    cov33_t m_out;
    for (int i=0; i<3; i++) {
      for (int j=i; j<3; j++)  {
	m_out(i,j) = m_in.matrix()(i,j);
      }
    }
    return m_out;
  }
  
  cov99_t makeCovarianceMatrix(const cov33_t cov_vtx1,
			       const cov77_t cov_vtx2) 
  {
    cov99_t cov;
    cov.Place_at(cov_vtx1,0,0);
    cov.Place_at(cov_vtx2.Sub<cov66_t>(0,0),3,3);
    return cov;
  }

  jac9_t makeJacobianVector3d(const AlgebraicVector3 &vtx1, 
			      const AlgebraicVector3 &vtx2, 
			      const AlgebraicVector3 &momentum) 
  {
    jac9_t jac;
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double factor2 = 1. / ROOT::Math::Mag2(momentum);
    const double lifetime = ROOT::Math::Dot(dist, momentum) * factor2;
    jac.Place_at(-momentum*factor2,0);
    jac.Place_at( momentum*factor2,3);
    jac.Place_at( factor2*(dist-2*lifetime*momentum*factor2),6);
    return jac;
  }
  
  jac9_t makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			      ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum) 
  {
    return makeJacobianVector3d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
				AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
				AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
  }

  jac9_t makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
			      const AlgebraicVector3 &momentum) {
    jac9_t jac;
    const double momentumMag = ROOT::Math::Mag(momentum);
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double distMag = ROOT::Math::Mag(dist);
    const double factorPositionComponent = 1./(distMag*momentumMag);
    const double factorMomentumComponent = 1./pow(momentumMag,3);
    jac(0)=-dist(0)*factorPositionComponent;
    jac(1)=-dist(1)*factorPositionComponent;
    jac(3)= dist(0)*factorPositionComponent;
    jac(4)= dist(1)*factorPositionComponent;
    jac(6)= momentum(0)*factorMomentumComponent;
    jac(7)= momentum(1)*factorMomentumComponent;
    return jac;
  }
  
  jac9_t makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			      ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
    return makeJacobianVector2d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
				AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
				AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
  }
}

void crossProduct(Double_t vect_A[], Double_t vect_B[], Double_t cross_P[])

{
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

Double_t getMag(Double_t x, Double_t y,Double_t z)
{
        return sqrt(x*x+y*y+z*z);
}


Double_t getDCALineAndPoint(
                                Double_t x1[],
                                Double_t x2[],
                                Double_t p[]
                            )
{

    Double_t v1[3],v2[3],v3[3];
    v1[0]= x1[0]-p[0]; v1[1] = x1[1]- p[1] ; v1[2] = x1[1] - p[2];
    v2[0]= x2[0]-p[0]; v2[1] = x2[1]- p[1] ; v2[2] = x2[2] - p[2];

    crossProduct(v1,v2,v3);

//    std::cout<<"\nV1  : "<<v1[0]<<","<<v1[1]<<" , "<<v1[2]<<"\n";
//    std::cout<<"\nV2  : "<<v2[0]<<","<<v2[1]<<" , "<<v2[2]<<"\n";
//    std::cout<<"\nc p  : "<<v3[0]<<","<<v3[1]<<" , "<<v3[2]<<"\n";
//    std::cout<<" Mag Num : "<<getMag(v3[0],v3[1],v3[2])<<" , Mag Den : "<<getMag(x1[0]-x2[0],x1[1]-x2[1],x1[2]-x2[2])<<"\n";
    return getMag(v3[0],v3[1],v3[2]) / getMag(x1[0]-x2[0],x1[1]-x2[1],x1[2]-x2[2]);

}


namespace{

  bool dr_match(const LorentzVector& reco , const LorentzVector& gen){
    if (fabs(reco.pt()-gen.pt())/gen.pt()<0.1 and deltaR(reco,gen)<0.02)
      return true;
    return false;
  }

  std::vector<unsigned int> 
    get_depth_from_permutation(const std::vector<unsigned int>& elements){
    std::vector<unsigned int> result;
    unsigned int counter(0);
    for (auto element: elements){
      if (element==0){
	counter++;
      } else {
	result.push_back(counter);
	counter = 0;
      }
    }
    result.push_back(counter);
    return result;
  }

  bool is_acceptable(const reco::Candidate* cand){
    if ( not cand) return false; 
    // skip quarks
    if ( abs(cand->pdgId())<10 ) return false;
    // skip protons
    if ( abs(cand->pdgId())==2212 ) return false;
    // skip gluons
    if ( abs(cand->pdgId())==21 ) return false;
    return true;
  }

  // depth 0 - first mother

  const reco::Candidate* get_mother(const reco::Candidate* cand, unsigned int depth){
    if (not cand) return 0;
    const reco::Candidate* mother = cand->mother();
    unsigned int i = 0;
    while ( is_acceptable(mother) and i<depth ){
      i++;
      mother = mother->mother();
    }
    if (is_acceptable(mother))
      return mother;
    else
      return 0;
  }

  const reco::Candidate* 
    find_common_ancestor(const std::vector<const reco::Candidate*>& particles, 
			 unsigned int max_depth=10){
    auto n = particles.size();
    for (unsigned int depth=0; depth<max_depth; ++depth){
      // make a list of elements (0) and separators (1) and
      // find all possible permutations of the elements
      std::vector<unsigned int> elements;
      for (unsigned int i=0; i<depth; ++i)
	elements.push_back(0);
      for (unsigned int i=0; i<n-1; ++i)
	elements.push_back(1);
      do {
	auto depth_vector = get_depth_from_permutation(elements);
	const reco::Candidate* common_mother(0);
	for (unsigned int i=0; i<n; ++i){
	  auto mother = get_mother(particles[i],depth_vector[i]);
	  if (not mother) {
	    common_mother = 0;
	    break;
	  }
	  if (not common_mother) common_mother = mother;
	  if (common_mother != mother) {
	    common_mother = 0;
	    break;
	  }	  
	}
	if (common_mother) return common_mother;
      } while(std::next_permutation(elements.begin(), elements.end()));
    }
    return 0;
  }

}


namespace {
  math::XYZPoint getProductionVertex( const reco::Candidate* cand){
    if (not cand) return math::XYZPoint();
    const reco::Candidate* primary = cand;
    // handle oscillation and radiation
    while (primary->mother() and abs(primary->pdgId())==abs(primary->mother()->pdgId()))
      primary = primary->mother();
    return primary->vertex();
  }

  double computeDecayTime( const GenMatchInfo& info ){
    if (not info.match) return -1.0;
    auto prod_vtx = getProductionVertex(info.match);
    if (prod_vtx.r()<1e-12) return -2.0;
    return (prod_vtx-info.mm_vtx).r()/TMath::Ccgs()*info.match->mass()/info.match->p();
  }
}


#endif
