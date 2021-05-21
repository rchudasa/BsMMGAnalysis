#ifndef TREECONTENT_H
#define TREECONTENT_H

#include <vector>
#include <string>
#include "TTree.h"

#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"


class TreeContent
{
 public:
  
   TreeContent ();
  ~TreeContent ();

  void Init ();
  void ClearNTuple ();
  void ClearScalars();
  void ClearScalarsMonteCarlo();
  void ClearVectors();
  void ClearVectorsMonteCarlo();
  void ClearMonteCarlo ();
  void MakeTreeBranches (TTree* theTree);


  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  unsigned int              runN;
  unsigned int              eventN;
  unsigned int              recoVtxN;
  double                    evWeight;
  double                    evWeightE2;
  unsigned int              numEventsTried;
  unsigned int              numEventsPassed;


  // ###########
  // # Trigger #
  // ###########
  std::vector<std::string>  TrigTable;
  std::vector<int>          TrigPrescales;
  std::vector<std::string>  L1Table;
  std::vector<int>          L1Prescales;

  
  
  // ###########################
  // # Number of B0 candidates #
  // ###########################
  unsigned int              nB;
  
  // ############################
  // # Pileup information in MC #
  // ############################
  std::vector<double>       bunchXingMC, numInteractionsMC, trueNumInteractionsMC;
  // Comment:
  // - PileupSummaryInfo::getTrueNumInteractions() gives the distribution of the mean number of interactions per crossing.
  // Since this is the mean value of the poisson distribution from which the number of interactions in- and out-of-time are
  // generated, no additional information should be required for reweighting if these values are matched in data and Monte Carlo.
  // - PileupSummaryInfo::getPU_NumInteractions() gives the expected mean number of interactions per crossing for each LumiSection.
  // Therefore the pileup histogram will contain the distribution of the number of interactions one would actually observe given
  // a poisson of that mean. So, this distribution is what one would see if one counted the number of events seen in a given beam
  // crossing (by looking at the number of vertices in data, for example. This would be appropriate for pileup reweighting based
  // on in-time-only distributions.

  // ################################
  // # Primary Vertex and Beam Spot #
  // ################################
  double                    bsX, bsY;

  // ### mu+ mu- variables ###
  std::vector<double>   mumuPt;
  std::vector<double>   mumuEta;
  std::vector<double>   mumuRapidity;
  std::vector<double>   mumuPhi;
  std::vector<double>   mumuMass;
  std::vector<double>   mumuMassE;
  std::vector<double>   mumuPx;
  std::vector<double>   mumuPy;
  std::vector<double>   mumuPz;
  std::vector<double>   mumuDR;

  // ### mu+ mu- Vtx ###
  std::vector<double>  mumuVtxCL       ;
  std::vector<double>  mumuVtxX        ;
  std::vector<double>  mumuVtxY        ;
  std::vector<double>  mumuVtxZ        ;
  std::vector<double>  mumuVtxChi2     ;
  std::vector<double>  mumuVtxNdof     ;
  std::vector<double>  mumuVtxProb     ;
  std::vector<bool>    mumuVtxIsGoodFit;
  std::vector<double>  mumuCosAlphaBS  ;
  std::vector<double>  mumuCosAlphaBSE ; 
  std::vector<double>  mumuLBS         ;
  std::vector<double>  mumuLBSE        ;
  std::vector<double>  mumuDCA         ;
  std::vector<double>   mumuLS;
  std::vector<double>   mumuLSErr;


  
  // ### mu- ###
  int 	               nMuM; 
  std::vector<bool>    mumHighPurity;
  std::vector<double>  mumPt;
  std::vector<double>  mumEta;
  std::vector<double>  mumPhi;
  std::vector<double>  mumCL; 
  std::vector<double>  mumNormChi2;
  std::vector<double>  mumPx;
  std::vector<double>  mumPy;
  std::vector<double>  mumPz;
  std::vector<double>  mumDCAVtx;
  std::vector<double>  mumDCAVtxE;
  std::vector<double>  mumDCABS;
  std::vector<double>  mumDCABSE;
  std::vector<double>  mumKinkChi2;
  std::vector<double>  mumFracHits;
  std::vector<double>  mumdxyBS;
  std::vector<double>  mumdzBS;
  std::vector<double>  mumMinIP2D;
  std::vector<double>  mumMinIP2DE;
  std::vector<double>  mumMinIP;
  std::vector<double>  mumMinIPS;
  std::vector<double>  mumDeltaRwithMC;
  std::vector<std::string> mumCat;
  std::vector<int>     mumCharge;
  std::vector<int>     mumNPixHits;
  std::vector<int>     mumNPixLayers;
  std::vector<int>     mumNTrkHits;
  std::vector<int>     mumNTrkLayers;
  std::vector<int>     mumNMuonHits;
  std::vector<int>     mumNMatchStation;
  std::vector<float>   mumIso;
  std::vector<float>   mumIsoPt;
  std::vector<float>   mumIsodR;
  std::vector<bool>    mum_isGlobalMuon;
  std::vector<bool>    mum_isTrackerMuon;
  std::vector<bool>    mum_StandAloneMuon;
  std::vector<bool>    mum_isCaloMuon;
  std::vector<bool>    mum_isPFMuon;


  // ### mu+ ###
  int 	               nMuP; 
  std::vector<bool>    mupHighPurity;
  std::vector<double>  mupPt;
  std::vector<double>  mupEta;
  std::vector<double>  mupPhi;
  std::vector<double>  mupCL; 
  std::vector<double>  mupNormChi2;
  std::vector<double>  mupPx;
  std::vector<double>  mupPy;
  std::vector<double>  mupPz;
  std::vector<double>  mupDCAVtx;
  std::vector<double>  mupDCAVtxE;
  std::vector<double>  mupDCABS;
  std::vector<double>  mupDCABSE;
  std::vector<double>  mupKinkChi2;
  std::vector<double>  mupFracHits;
  std::vector<double>  mupdxyBS;
  std::vector<double>  mupdzBS;
  std::vector<double>  mupMinIP2D;
  std::vector<double>  mupMinIP2DE;
  std::vector<double>  mupMinIP;
  std::vector<double>  mupMinIPS;
  std::vector<double>  mupDeltaRwithMC;
  std::vector<std::string> mupCat;
  std::vector<int>     mupCharge;
  std::vector<int>     mupNPixHits;
  std::vector<int>     mupNPixLayers;
  std::vector<int>     mupNTrkHits;
  std::vector<int>     mupNTrkLayers;
  std::vector<int>     mupNMuonHits;
  std::vector<int>     mupNMatchStation;
  std::vector<float>   mupIso;
  std::vector<float>   mupIsoPt;
  std::vector<float>   mupIsodR;
  std::vector<bool>    mup_isGlobalMuon;
  std::vector<bool>    mup_isTrackerMuon;
  std::vector<bool>    mup_StandAloneMuon;
  std::vector<bool>    mup_isCaloMuon;
  std::vector<bool>    mup_isPFMuon;

  int         nSC;
  std::vector<float>  scE;
  std::vector<float>  scEt;
  std::vector<float>  scEta;
  std::vector<float>  scPhi;
  std::vector<float>  scX;
  std::vector<float>  scY;
  std::vector<float>  scZ;
  std::vector<float>  scEtaWidth;
  std::vector<float>  scPhiWidth;
  std::vector<float>  scRawE;
  std::vector<float>  scRawEt;
  
  
 /* std::vector<double>   muon_pt, muon_eta, muon_phi, mum_dz, muon_dxy;
  std::vector<double>   mum_dz_error, muon_dxy_error;
  std::vector<double>   muon_vx,muon_vy,muon_vz,muon_vertexChi2,muon_vertexNDoF;

  std::vector<int>	muon_charge;
  std::vector<bool>	muon_isGlobalMuon,muon_isTrackerMuon,muon_StandAloneMuon,muon_isCaloMuon,muon_isPFMuon;
  int 	            	nMuons;  // total number of muons in the reco collction before basic cuts are applied reco::MuonCollection->size()

  std::vector<uint64_t> muon_selector; 
  std::vector<bool>	muon_isIsolationValid;
  std::vector<bool>	muon_isPFIsolationValid;
  
  std::vector<double>  muon_isolationR03_trackSumPt;
  std::vector<double>  muon_isolationR03_trackEcalSumEt;
  std::vector<double>  muon_isolationR03_trackHcalSumEt;
  std::vector<double>  muon_isolationR03_trackHoSumEt;
  std::vector<int>     muon_isolationR03_trackNTracks;
  std::vector<int>     muon_isolationR03_trackNJets;
  std::vector<double>  muon_isolationR03_trackerVetoSumPt;
  std::vector<double>  muon_isolationR03_emVetoSumEt;
  std::vector<double>  muon_isolationR03_hadVetoSumEt;
  std::vector<double>  muon_isolationR03_hoVetoEt;
  
  std::vector<double>  muon_isolationR05_trackSumPt;
  std::vector<double>  muon_isolationR05_trackEcalSumEt;
  std::vector<double>  muon_isolationR05_trackHcalSumEt;
  std::vector<double>  muon_isolationR05_trackHoSumEt;
  std::vector<int>     muon_isolationR05_trackNTracks;
  std::vector<int>     muon_isolationR05_trackNJets;
  std::vector<double>  muon_isolationR05_trackerVetoSumPt;
  std::vector<double>  muon_isolationR05_emVetoSumEt;
  std::vector<double>  muon_isolationR05_hadVetoSumEt;
  std::vector<double>  muon_isolationR05_hoVetoEt;
  
  std::vector<double>  muon_PFIsolationR03_sumChargedHadronPt;
  std::vector<double>  muon_PFIsolationR03_sumChargedParticlePt;
  std::vector<double>  muon_PFIsolationR03_sumNeutralHadronEt;
  std::vector<double>  muon_PFIsolationR03_sumPhotonEt;
  std::vector<double>  muon_PFIsolationR03_sumNeutralHadronEtHighThreshold;
  std::vector<double>  muon_PFIsolationR03_sumPhotonEtHighThreshold;
  std::vector<double>  muon_PFIsolationR03_sumPUPt;
  
  std::vector<double>  muon_PFIsolationR04_sumChargedHadronPt;
  std::vector<double>  muon_PFIsolationR04_sumChargedParticlePt;
  std::vector<double>  muon_PFIsolationR04_sumNeutralHadronEt;
  std::vector<double>  muon_PFIsolationR04_sumPhotonEt;
  std::vector<double>  muon_PFIsolationR04_sumNeutralHadronEtHighThreshold;
  std::vector<double>  muon_PFIsolationR04_sumPhotonEtHighThreshold;
  std::vector<double>  muon_PFIsolationR04_sumPUPt;
  
  std::vector<double> muon_dcaToBS;
  std::vector<double> muon_dcaToBS_error;*/
  

  
  // # BeamSpot # //
  double beamspot_x,beamspot_y,beamspot_z,beamspot_x_error,beamspot_y_error,beamspot_z_error;
  double beamspot_dxdz,beamspot_dydz,beamspot_sigmaZ,beamspot_dxdz_error,beamspot_dydz_error,beamspot_sigmaZError;
  double beamspot_beamWidthX,beamspot_beamWidthY,beamspot_beamWidthX_error,beamspot_beamWidthY_error;
 
  // # offlinePrimaryVertices # //
  std::vector<bool> primaryVertex_isFake;
  std::vector<double> primaryVertex_x, primaryVertex_y,primaryVertex_z,primaryVertex_t;
  std::vector<double> primaryVertex_x_error, primaryVertex_y_error,primaryVertex_z_error,primaryVertex_t_error;
  std::vector<double> primaryVertex_ntracks,primaryVertex_ndof,primaryVertex_chi2,primaryVertex_normalizedChi2;

 // # Gen Data # //
 std::vector<bool>  gen_hasAValid_candidate;
 std::vector<double> gen_Bs_pt,gen_Bs_eta,gen_Bs_phi,gen_Bs_pz,gen_Bs_pdgId;
 std::vector<double> gen_BsMuonM_pt,gen_BsMuonM_eta,gen_BsMuonM_phi;
 std::vector<double> gen_BsMuonP_pt,gen_BsMuonP_eta,gen_BsMuonP_phi;
 std::vector<double> gen_BsPhoton_pt,gen_BsPhoton_eta,gen_BsPhoton_phi;
 std::vector<int> gen_BsPhotonMultiplicity,gen_BsMuonMMultiplicity,gen_BsMuonMPultiplicity;

};

#endif
