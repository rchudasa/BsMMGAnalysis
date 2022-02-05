//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  2 16:01:34 2022 by ROOT version 6.14/09
// from TTree mergedTree/bmm5 and bmmg merged tree
// found on file: /grid_mnt/t3storage3/athachay/bs2mumug/run2studies/bParkingAnalysis/analysis/CMSSW_10_6_29/src/BsMMGAnalysis/MergeWithBMMNtuples/MergeFiles/results/bph6A.root
//////////////////////////////////////////////////////////

#ifndef MergedBMMX2018Data_h
#define MergedBMMX2018Data_h
#define MergedBMMX2018Data_cxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MergedBMMX2018Data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxbG_nscMinDrWithGsfElectornSC = 37;
   static constexpr Int_t kMaxbG_scMinDrWithGsfElectornSC = 1;
   static constexpr Int_t kMaxbG_nscFoundGsfMatch = 37;
   static constexpr Int_t kMaxbG_scFoundGsfMatch = 1;

   // Declaration of leaf types
   UInt_t          b5_run;
   UInt_t          b5_luminosityBlock;
   ULong64_t       b5_event;
   UInt_t          b5_nMuonId;
   Float_t         b5_MuonId_chi2LocalPosition[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_glbNormChi2[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_glbTrackProbability[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match1_dX[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match1_dY[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match1_pullDxDz[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match1_pullDyDz[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match1_pullX[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match1_pullY[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match2_dX[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match2_dY[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match2_pullDxDz[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match2_pullDyDz[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match2_pullX[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_match2_pullY[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_newSoftMuonMva[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_trkKink[8];   //[b5_nMuonId]
   Float_t         b5_MuonId_trkValidFrac[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_highPurity[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_nLostHitsInner[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_nLostHitsOn[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_nLostHitsOuter[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_nPixels[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_nValidHits[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_trkLayers[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_trkLostLayersInner[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_trkLostLayersOn[8];   //[b5_nMuonId]
   Int_t           b5_MuonId_trkLostLayersOuter[8];   //[b5_nMuonId]
   UInt_t          b5_nmm;
   Float_t         b5_mm_bdt[10];   //[b5_nmm]
   Float_t         b5_mm_doca[10];   //[b5_nmm]
   Float_t         b5_mm_docatrk[10];   //[b5_nmm]
   Float_t         b5_mm_iso[10];   //[b5_nmm]
   Float_t         b5_mm_kal_lxy[10];   //[b5_nmm]
   Float_t         b5_mm_kal_mass[10];   //[b5_nmm]
   Float_t         b5_mm_kal_slxy[10];   //[b5_nmm]
   Float_t         b5_mm_kal_vtx_prob[10];   //[b5_nmm]
   Float_t         b5_mm_kin_alpha[10];   //[b5_nmm]
   Float_t         b5_mm_kin_alphaBS[10];   //[b5_nmm]
   Float_t         b5_mm_kin_alphaBSErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_alphaErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_eta[10];   //[b5_nmm]
   Float_t         b5_mm_kin_l3d[10];   //[b5_nmm]
   Float_t         b5_mm_kin_lxy[10];   //[b5_nmm]
   Float_t         b5_mm_kin_mass[10];   //[b5_nmm]
   Float_t         b5_mm_kin_massErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_mu1eta[10];   //[b5_nmm]
   Float_t         b5_mm_kin_mu1phi[10];   //[b5_nmm]
   Float_t         b5_mm_kin_mu1pt[10];   //[b5_nmm]
   Float_t         b5_mm_kin_mu2eta[10];   //[b5_nmm]
   Float_t         b5_mm_kin_mu2phi[10];   //[b5_nmm]
   Float_t         b5_mm_kin_mu2pt[10];   //[b5_nmm]
   Float_t         b5_mm_kin_phi[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pt[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pv2ip[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pv2ipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pv2lip[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pv2lipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pv2lipSig[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pv_z[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pv_zErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pvip[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pvipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pvlip[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pvlipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_pvlipSig[10];   //[b5_nmm]
   Float_t         b5_mm_kin_sl3d[10];   //[b5_nmm]
   Float_t         b5_mm_kin_slxy[10];   //[b5_nmm]
   Float_t         b5_mm_kin_spv2ip[10];   //[b5_nmm]
   Float_t         b5_mm_kin_spvip[10];   //[b5_nmm]
   Float_t         b5_mm_kin_tau[10];   //[b5_nmm]
   Float_t         b5_mm_kin_taue[10];   //[b5_nmm]
   Float_t         b5_mm_kin_tauxy[10];   //[b5_nmm]
   Float_t         b5_mm_kin_tauxye[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_chi2dof[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_prob[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_x[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_xErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_y[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_yErr[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_z[10];   //[b5_nmm]
   Float_t         b5_mm_kin_vtx_zErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_alpha[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_alphaBS[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_alphaBSErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_alphaErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_eta[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_l3d[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_lxy[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_mass[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_massErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_phi[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pt[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pv2ip[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pv2ipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pv2lip[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pv2lipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pv2lipSig[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pv_z[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pv_zErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pvip[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pvipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pvlip[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pvlipErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_pvlipSig[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_sl3d[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_slxy[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_spv2ip[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_spvip[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_tau[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_taue[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_tauxy[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_tauxye[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_chi2dof[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_prob[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_x[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_xErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_y[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_yErr[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_z[10];   //[b5_nmm]
   Float_t         b5_mm_kinpc_vtx_zErr[10];   //[b5_nmm]
   Float_t         b5_mm_m1iso[10];   //[b5_nmm]
   Float_t         b5_mm_m2iso[10];   //[b5_nmm]
   Float_t         b5_mm_mass[10];   //[b5_nmm]
   Float_t         b5_mm_mu1_eta[10];   //[b5_nmm]
   Float_t         b5_mm_mu1_phi[10];   //[b5_nmm]
   Float_t         b5_mm_mu1_pt[10];   //[b5_nmm]
   Float_t         b5_mm_mu2_eta[10];   //[b5_nmm]
   Float_t         b5_mm_mu2_phi[10];   //[b5_nmm]
   Float_t         b5_mm_mu2_pt[10];   //[b5_nmm]
   Float_t         b5_mm_mva[10];   //[b5_nmm]
   Float_t         b5_mm_otherVtxMaxProb[10];   //[b5_nmm]
   Float_t         b5_mm_otherVtxMaxProb1[10];   //[b5_nmm]
   Float_t         b5_mm_otherVtxMaxProb2[10];   //[b5_nmm]
   Int_t           b5_mm_closetrk[10];   //[b5_nmm]
   Int_t           b5_mm_closetrks1[10];   //[b5_nmm]
   Int_t           b5_mm_closetrks2[10];   //[b5_nmm]
   Int_t           b5_mm_closetrks3[10];   //[b5_nmm]
   Int_t           b5_mm_kal_valid[10];   //[b5_nmm]
   Int_t           b5_mm_kin_valid[10];   //[b5_nmm]
   Int_t           b5_mm_kinpc_valid[10];   //[b5_nmm]
   Int_t           b5_mm_mu1_index[10];   //[b5_nmm]
   Int_t           b5_mm_mu1_pdgId[10];   //[b5_nmm]
   Int_t           b5_mm_mu2_index[10];   //[b5_nmm]
   Int_t           b5_mm_mu2_pdgId[10];   //[b5_nmm]
   Int_t           b5_mm_nBMTrks[10];   //[b5_nmm]
   Int_t           b5_mm_nDisTrks[10];   //[b5_nmm]
   Int_t           b5_mm_nTrks[10];   //[b5_nmm]
   UInt_t          b5_nd0;
   Float_t         b5_d0_doca[7];   //[b5_nd0]
   Float_t         b5_d0_kaon_eta[7];   //[b5_nd0]
   Float_t         b5_d0_kaon_phi[7];   //[b5_nd0]
   Float_t         b5_d0_kaon_pt[7];   //[b5_nd0]
   Float_t         b5_d0_kaon_sip[7];   //[b5_nd0]
   Float_t         b5_d0_kin_cosAlphaXY[7];   //[b5_nd0]
   Float_t         b5_d0_kin_eta[7];   //[b5_nd0]
   Float_t         b5_d0_kin_lxy[7];   //[b5_nd0]
   Float_t         b5_d0_kin_mass[7];   //[b5_nd0]
   Float_t         b5_d0_kin_massErr[7];   //[b5_nd0]
   Float_t         b5_d0_kin_phi[7];   //[b5_nd0]
   Float_t         b5_d0_kin_pt[7];   //[b5_nd0]
   Float_t         b5_d0_kin_sipBS[7];   //[b5_nd0]
   Float_t         b5_d0_kin_sipPV[7];   //[b5_nd0]
   Float_t         b5_d0_kin_slxy[7];   //[b5_nd0]
   Float_t         b5_d0_kin_vtx_chi2dof[7];   //[b5_nd0]
   Float_t         b5_d0_kin_vtx_prob[7];   //[b5_nd0]
   Float_t         b5_d0_mass[7];   //[b5_nd0]
   Float_t         b5_d0_pion_eta[7];   //[b5_nd0]
   Float_t         b5_d0_pion_phi[7];   //[b5_nd0]
   Float_t         b5_d0_pion_pt[7];   //[b5_nd0]
   Float_t         b5_d0_pion_sip[7];   //[b5_nd0]
   Int_t           b5_d0_kaon_mu_index[7];   //[b5_nd0]
   Int_t           b5_d0_kin_valid[7];   //[b5_nd0]
   Int_t           b5_d0_pion_mu_index[7];   //[b5_nd0]
   UInt_t          b5_nks;
   Float_t         b5_ks_doca[9];   //[b5_nks]
   Float_t         b5_ks_kin_cosAlphaXY[9];   //[b5_nks]
   Float_t         b5_ks_kin_eta[9];   //[b5_nks]
   Float_t         b5_ks_kin_lxy[9];   //[b5_nks]
   Float_t         b5_ks_kin_mass[9];   //[b5_nks]
   Float_t         b5_ks_kin_massErr[9];   //[b5_nks]
   Float_t         b5_ks_kin_phi[9];   //[b5_nks]
   Float_t         b5_ks_kin_pt[9];   //[b5_nks]
   Float_t         b5_ks_kin_sipBS[9];   //[b5_nks]
   Float_t         b5_ks_kin_sipPV[9];   //[b5_nks]
   Float_t         b5_ks_kin_slxy[9];   //[b5_nks]
   Float_t         b5_ks_kin_vtx_chi2dof[9];   //[b5_nks]
   Float_t         b5_ks_kin_vtx_prob[9];   //[b5_nks]
   Float_t         b5_ks_mass[9];   //[b5_nks]
   Float_t         b5_ks_trk1_eta[9];   //[b5_nks]
   Float_t         b5_ks_trk1_phi[9];   //[b5_nks]
   Float_t         b5_ks_trk1_pt[9];   //[b5_nks]
   Float_t         b5_ks_trk1_sip[9];   //[b5_nks]
   Float_t         b5_ks_trk2_eta[9];   //[b5_nks]
   Float_t         b5_ks_trk2_phi[9];   //[b5_nks]
   Float_t         b5_ks_trk2_pt[9];   //[b5_nks]
   Float_t         b5_ks_trk2_sip[9];   //[b5_nks]
   Int_t           b5_ks_kin_valid[9];   //[b5_nks]
   Int_t           b5_ks_trk1_mu_index[9];   //[b5_nks]
   Int_t           b5_ks_trk2_mu_index[9];   //[b5_nks]
   UInt_t          b5_nlambda;
   Float_t         b5_lambda_doca[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_cosAlphaXY[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_eta[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_lxy[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_mass[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_massErr[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_phi[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_pt[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_sipBS[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_sipPV[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_slxy[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_vtx_chi2dof[9];   //[b5_nlambda]
   Float_t         b5_lambda_kin_vtx_prob[9];   //[b5_nlambda]
   Float_t         b5_lambda_mass[9];   //[b5_nlambda]
   Float_t         b5_lambda_pion_eta[9];   //[b5_nlambda]
   Float_t         b5_lambda_pion_phi[9];   //[b5_nlambda]
   Float_t         b5_lambda_pion_pt[9];   //[b5_nlambda]
   Float_t         b5_lambda_pion_sip[9];   //[b5_nlambda]
   Float_t         b5_lambda_proton_eta[9];   //[b5_nlambda]
   Float_t         b5_lambda_proton_phi[9];   //[b5_nlambda]
   Float_t         b5_lambda_proton_pt[9];   //[b5_nlambda]
   Float_t         b5_lambda_proton_sip[9];   //[b5_nlambda]
   Int_t           b5_lambda_kin_valid[9];   //[b5_nlambda]
   Int_t           b5_lambda_pion_mu_index[9];   //[b5_nlambda]
   Int_t           b5_lambda_proton_mu_index[9];   //[b5_nlambda]
   UInt_t          b5_nphi;
   Float_t         b5_phi_doca[25];   //[b5_nphi]
   Float_t         b5_phi_ds_cosAlphaXY[25];   //[b5_nphi]
   Float_t         b5_phi_ds_eta[25];   //[b5_nphi]
   Float_t         b5_phi_ds_lxy[25];   //[b5_nphi]
   Float_t         b5_phi_ds_mass[25];   //[b5_nphi]
   Float_t         b5_phi_ds_massErr[25];   //[b5_nphi]
   Float_t         b5_phi_ds_phi[25];   //[b5_nphi]
   Float_t         b5_phi_ds_pion_eta[25];   //[b5_nphi]
   Float_t         b5_phi_ds_pion_mu_index[25];   //[b5_nphi]
   Float_t         b5_phi_ds_pion_phi[25];   //[b5_nphi]
   Float_t         b5_phi_ds_pion_pt[25];   //[b5_nphi]
   Float_t         b5_phi_ds_pt[25];   //[b5_nphi]
   Float_t         b5_phi_ds_sipBS[25];   //[b5_nphi]
   Float_t         b5_phi_ds_sipPV[25];   //[b5_nphi]
   Float_t         b5_phi_ds_slxy[25];   //[b5_nphi]
   Float_t         b5_phi_ds_vtx_chi2dof[25];   //[b5_nphi]
   Float_t         b5_phi_ds_vtx_prob[25];   //[b5_nphi]
   Float_t         b5_phi_kin_cosAlphaXY[25];   //[b5_nphi]
   Float_t         b5_phi_kin_eta[25];   //[b5_nphi]
   Float_t         b5_phi_kin_lxy[25];   //[b5_nphi]
   Float_t         b5_phi_kin_mass[25];   //[b5_nphi]
   Float_t         b5_phi_kin_massErr[25];   //[b5_nphi]
   Float_t         b5_phi_kin_phi[25];   //[b5_nphi]
   Float_t         b5_phi_kin_pt[25];   //[b5_nphi]
   Float_t         b5_phi_kin_sipBS[25];   //[b5_nphi]
   Float_t         b5_phi_kin_sipPV[25];   //[b5_nphi]
   Float_t         b5_phi_kin_slxy[25];   //[b5_nphi]
   Float_t         b5_phi_kin_vtx_chi2dof[25];   //[b5_nphi]
   Float_t         b5_phi_kin_vtx_prob[25];   //[b5_nphi]
   Float_t         b5_phi_mass[25];   //[b5_nphi]
   Float_t         b5_phi_trk1_eta[25];   //[b5_nphi]
   Float_t         b5_phi_trk1_phi[25];   //[b5_nphi]
   Float_t         b5_phi_trk1_pt[25];   //[b5_nphi]
   Float_t         b5_phi_trk1_sip[25];   //[b5_nphi]
   Float_t         b5_phi_trk2_eta[25];   //[b5_nphi]
   Float_t         b5_phi_trk2_phi[25];   //[b5_nphi]
   Float_t         b5_phi_trk2_pt[25];   //[b5_nphi]
   Float_t         b5_phi_trk2_sip[25];   //[b5_nphi]
   Int_t           b5_phi_kin_valid[25];   //[b5_nphi]
   Int_t           b5_phi_trk1_mu_index[25];   //[b5_nphi]
   Int_t           b5_phi_trk2_mu_index[25];   //[b5_nphi]
   Float_t         b5_CaloMET_phi;
   Float_t         b5_CaloMET_pt;
   Float_t         b5_CaloMET_sumEt;
   Float_t         b5_ChsMET_phi;
   Float_t         b5_ChsMET_pt;
   Float_t         b5_ChsMET_sumEt;
   UInt_t          b5_nCorrT1METJet;
   Float_t         b5_CorrT1METJet_area[22];   //[b5_nCorrT1METJet]
   Float_t         b5_CorrT1METJet_eta[22];   //[b5_nCorrT1METJet]
   Float_t         b5_CorrT1METJet_muonSubtrFactor[22];   //[b5_nCorrT1METJet]
   Float_t         b5_CorrT1METJet_phi[22];   //[b5_nCorrT1METJet]
   Float_t         b5_CorrT1METJet_rawPt[22];   //[b5_nCorrT1METJet]
   Float_t         b5_DeepMETResolutionTune_phi;
   Float_t         b5_DeepMETResolutionTune_pt;
   Float_t         b5_DeepMETResponseTune_phi;
   Float_t         b5_DeepMETResponseTune_pt;
   UInt_t          b5_nElectron;
   Float_t         b5_Electron_deltaEtaSC[6];   //[b5_nElectron]
   Float_t         b5_Electron_dr03EcalRecHitSumEt[6];   //[b5_nElectron]
   Float_t         b5_Electron_dr03HcalDepth1TowerSumEt[6];   //[b5_nElectron]
   Float_t         b5_Electron_dr03TkSumPt[6];   //[b5_nElectron]
   Float_t         b5_Electron_dr03TkSumPtHEEP[6];   //[b5_nElectron]
   Float_t         b5_Electron_dxy[6];   //[b5_nElectron]
   Float_t         b5_Electron_dxyErr[6];   //[b5_nElectron]
   Float_t         b5_Electron_dz[6];   //[b5_nElectron]
   Float_t         b5_Electron_dzErr[6];   //[b5_nElectron]
   Float_t         b5_Electron_eCorr[6];   //[b5_nElectron]
   Float_t         b5_Electron_eInvMinusPInv[6];   //[b5_nElectron]
   Float_t         b5_Electron_energyErr[6];   //[b5_nElectron]
   Float_t         b5_Electron_eta[6];   //[b5_nElectron]
   Float_t         b5_Electron_hoe[6];   //[b5_nElectron]
   Float_t         b5_Electron_ip3d[6];   //[b5_nElectron]
   Float_t         b5_Electron_jetPtRelv2[6];   //[b5_nElectron]
   Float_t         b5_Electron_jetRelIso[6];   //[b5_nElectron]
   Float_t         b5_Electron_mass[6];   //[b5_nElectron]
   Float_t         b5_Electron_miniPFRelIso_all[6];   //[b5_nElectron]
   Float_t         b5_Electron_miniPFRelIso_chg[6];   //[b5_nElectron]
   Float_t         b5_Electron_mvaFall17V1Iso[6];   //[b5_nElectron]
   Float_t         b5_Electron_mvaFall17V1noIso[6];   //[b5_nElectron]
   Float_t         b5_Electron_mvaFall17V2Iso[6];   //[b5_nElectron]
   Float_t         b5_Electron_mvaFall17V2noIso[6];   //[b5_nElectron]
   Float_t         b5_Electron_pfRelIso03_all[6];   //[b5_nElectron]
   Float_t         b5_Electron_pfRelIso03_chg[6];   //[b5_nElectron]
   Float_t         b5_Electron_phi[6];   //[b5_nElectron]
   Float_t         b5_Electron_pt[6];   //[b5_nElectron]
   Float_t         b5_Electron_r9[6];   //[b5_nElectron]
   Float_t         b5_Electron_scEtOverPt[6];   //[b5_nElectron]
   Float_t         b5_Electron_sieie[6];   //[b5_nElectron]
   Float_t         b5_Electron_sip3d[6];   //[b5_nElectron]
   Float_t         b5_Electron_mvaTTH[6];   //[b5_nElectron]
   Int_t           b5_Electron_charge[6];   //[b5_nElectron]
   Int_t           b5_Electron_cutBased[6];   //[b5_nElectron]
   Int_t           b5_Electron_cutBased_Fall17_V1[6];   //[b5_nElectron]
   Int_t           b5_Electron_jetIdx[6];   //[b5_nElectron]
   Int_t           b5_Electron_pdgId[6];   //[b5_nElectron]
   Int_t           b5_Electron_photonIdx[6];   //[b5_nElectron]
   Int_t           b5_Electron_tightCharge[6];   //[b5_nElectron]
   Int_t           b5_Electron_vidNestedWPBitmap[6];   //[b5_nElectron]
   Int_t           b5_Electron_vidNestedWPBitmapHEEP[6];   //[b5_nElectron]
   Bool_t          b5_Electron_convVeto[6];   //[b5_nElectron]
   Bool_t          b5_Electron_cutBased_HEEP[6];   //[b5_nElectron]
   Bool_t          b5_Electron_isPFcand[6];   //[b5_nElectron]
   UChar_t         b5_Electron_jetNDauCharged[6];   //[b5_nElectron]
   UChar_t         b5_Electron_lostHits[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V1Iso_WP80[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V1Iso_WP90[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V1Iso_WPL[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V1noIso_WP80[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V1noIso_WP90[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V1noIso_WPL[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V2Iso_WP80[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V2Iso_WP90[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V2Iso_WPL[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V2noIso_WP80[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V2noIso_WP90[6];   //[b5_nElectron]
   Bool_t          b5_Electron_mvaFall17V2noIso_WPL[6];   //[b5_nElectron]
   UChar_t         b5_Electron_seedGain[6];   //[b5_nElectron]
   UInt_t          b5_nFatJet;
   Float_t         b5_FatJet_area[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagCMVA[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagCSVV2[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDBvL[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDBvLV2[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDBvL_noMD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDCvB[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDCvBV2[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDCvB_noMD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDCvL[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDCvLV2[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDDCvL_noMD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagDeepB[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_btagHbb[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_H4qvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_HbbvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_TvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_WvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_ZHbbvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_ZHccvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_ZbbvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_ZvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_bbvsLight[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTagMD_ccvsLight[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTag_H[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTag_QCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTag_QCDothers[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTag_TvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTag_WvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_deepTag_ZvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_eta[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_mass[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_msoftdrop[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_n2b1[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_n3b1[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNetMD_QCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNetMD_Xbb[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNetMD_Xcc[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNetMD_Xqq[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNet_H4qvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNet_HbbvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNet_HccvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNet_QCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNet_TvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNet_WvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_particleNet_ZvsQCD[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_phi[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_pt[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_rawFactor[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_tau1[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_tau2[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_tau3[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_tau4[5];   //[b5_nFatJet]
   Float_t         b5_FatJet_lsf3[5];   //[b5_nFatJet]
   Int_t           b5_FatJet_jetId[5];   //[b5_nFatJet]
   Int_t           b5_FatJet_subJetIdx1[5];   //[b5_nFatJet]
   Int_t           b5_FatJet_subJetIdx2[5];   //[b5_nFatJet]
   Int_t           b5_FatJet_electronIdx3SJ[5];   //[b5_nFatJet]
   Int_t           b5_FatJet_muonIdx3SJ[5];   //[b5_nFatJet]
   UInt_t          b5_nFsrPhoton;
   Float_t         b5_FsrPhoton_dROverEt2[2];   //[b5_nFsrPhoton]
   Float_t         b5_FsrPhoton_eta[2];   //[b5_nFsrPhoton]
   Float_t         b5_FsrPhoton_phi[2];   //[b5_nFsrPhoton]
   Float_t         b5_FsrPhoton_pt[2];   //[b5_nFsrPhoton]
   Float_t         b5_FsrPhoton_relIso03[2];   //[b5_nFsrPhoton]
   Int_t           b5_FsrPhoton_muonIdx[2];   //[b5_nFsrPhoton]
   UInt_t          b5_nIsoTrack;
   Float_t         b5_IsoTrack_dxy[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_dz[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_eta[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_pfRelIso03_all[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_pfRelIso03_chg[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_phi[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_pt[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_miniPFRelIso_all[5];   //[b5_nIsoTrack]
   Float_t         b5_IsoTrack_miniPFRelIso_chg[5];   //[b5_nIsoTrack]
   Int_t           b5_IsoTrack_fromPV[5];   //[b5_nIsoTrack]
   Int_t           b5_IsoTrack_pdgId[5];   //[b5_nIsoTrack]
   Bool_t          b5_IsoTrack_isHighPurityTrack[5];   //[b5_nIsoTrack]
   Bool_t          b5_IsoTrack_isPFcand[5];   //[b5_nIsoTrack]
   Bool_t          b5_IsoTrack_isFromLostTrack[5];   //[b5_nIsoTrack]
   UInt_t          b5_nJet;
   Float_t         b5_Jet_area[21];   //[b5_nJet]
   Float_t         b5_Jet_btagCMVA[21];   //[b5_nJet]
   Float_t         b5_Jet_btagCSVV2[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepB[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepC[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepCvB[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepCvL[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepFlavB[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepFlavC[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepFlavCvB[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepFlavCvL[21];   //[b5_nJet]
   Float_t         b5_Jet_btagDeepFlavQG[21];   //[b5_nJet]
   Float_t         b5_Jet_chEmEF[21];   //[b5_nJet]
   Float_t         b5_Jet_chFPV0EF[21];   //[b5_nJet]
   Float_t         b5_Jet_chFPV1EF[21];   //[b5_nJet]
   Float_t         b5_Jet_chFPV2EF[21];   //[b5_nJet]
   Float_t         b5_Jet_chFPV3EF[21];   //[b5_nJet]
   Float_t         b5_Jet_chHEF[21];   //[b5_nJet]
   Float_t         b5_Jet_eta[21];   //[b5_nJet]
   Float_t         b5_Jet_hfsigmaEtaEta[21];   //[b5_nJet]
   Float_t         b5_Jet_hfsigmaPhiPhi[21];   //[b5_nJet]
   Float_t         b5_Jet_mass[21];   //[b5_nJet]
   Float_t         b5_Jet_muEF[21];   //[b5_nJet]
   Float_t         b5_Jet_muonSubtrFactor[21];   //[b5_nJet]
   Float_t         b5_Jet_neEmEF[21];   //[b5_nJet]
   Float_t         b5_Jet_neHEF[21];   //[b5_nJet]
   Float_t         b5_Jet_phi[21];   //[b5_nJet]
   Float_t         b5_Jet_pt[21];   //[b5_nJet]
   Float_t         b5_Jet_puIdDisc[21];   //[b5_nJet]
   Float_t         b5_Jet_qgl[21];   //[b5_nJet]
   Float_t         b5_Jet_rawFactor[21];   //[b5_nJet]
   Float_t         b5_Jet_bRegCorr[21];   //[b5_nJet]
   Float_t         b5_Jet_bRegRes[21];   //[b5_nJet]
   Float_t         b5_Jet_cRegCorr[21];   //[b5_nJet]
   Float_t         b5_Jet_cRegRes[21];   //[b5_nJet]
   Int_t           b5_Jet_electronIdx1[21];   //[b5_nJet]
   Int_t           b5_Jet_electronIdx2[21];   //[b5_nJet]
   Int_t           b5_Jet_hfadjacentEtaStripsSize[21];   //[b5_nJet]
   Int_t           b5_Jet_hfcentralEtaStripSize[21];   //[b5_nJet]
   Int_t           b5_Jet_jetId[21];   //[b5_nJet]
   Int_t           b5_Jet_muonIdx1[21];   //[b5_nJet]
   Int_t           b5_Jet_muonIdx2[21];   //[b5_nJet]
   Int_t           b5_Jet_nElectrons[21];   //[b5_nJet]
   Int_t           b5_Jet_nMuons[21];   //[b5_nJet]
   Int_t           b5_Jet_puId[21];   //[b5_nJet]
   UChar_t         b5_Jet_nConstituents[21];   //[b5_nJet]
   Float_t         b5_L1PreFiringWeight_Dn;
   Float_t         b5_L1PreFiringWeight_Nom;
   Float_t         b5_L1PreFiringWeight_Up;
   Float_t         b5_MET_MetUnclustEnUpDeltaX;
   Float_t         b5_MET_MetUnclustEnUpDeltaY;
   Float_t         b5_MET_covXX;
   Float_t         b5_MET_covXY;
   Float_t         b5_MET_covYY;
   Float_t         b5_MET_phi;
   Float_t         b5_MET_pt;
   Float_t         b5_MET_significance;
   Float_t         b5_MET_sumEt;
   Float_t         b5_MET_sumPtUnclustered;
   UInt_t          b5_nMuon;
   Float_t         b5_Muon_dxy[8];   //[b5_nMuon]
   Float_t         b5_Muon_dxyErr[8];   //[b5_nMuon]
   Float_t         b5_Muon_dxybs[8];   //[b5_nMuon]
   Float_t         b5_Muon_dz[8];   //[b5_nMuon]
   Float_t         b5_Muon_dzErr[8];   //[b5_nMuon]
   Float_t         b5_Muon_eta[8];   //[b5_nMuon]
   Float_t         b5_Muon_ip3d[8];   //[b5_nMuon]
   Float_t         b5_Muon_jetPtRelv2[8];   //[b5_nMuon]
   Float_t         b5_Muon_jetRelIso[8];   //[b5_nMuon]
   Float_t         b5_Muon_mass[8];   //[b5_nMuon]
   Float_t         b5_Muon_miniPFRelIso_all[8];   //[b5_nMuon]
   Float_t         b5_Muon_miniPFRelIso_chg[8];   //[b5_nMuon]
   Float_t         b5_Muon_pfRelIso03_all[8];   //[b5_nMuon]
   Float_t         b5_Muon_pfRelIso03_chg[8];   //[b5_nMuon]
   Float_t         b5_Muon_pfRelIso04_all[8];   //[b5_nMuon]
   Float_t         b5_Muon_phi[8];   //[b5_nMuon]
   Float_t         b5_Muon_pt[8];   //[b5_nMuon]
   Float_t         b5_Muon_ptErr[8];   //[b5_nMuon]
   Float_t         b5_Muon_segmentComp[8];   //[b5_nMuon]
   Float_t         b5_Muon_sip3d[8];   //[b5_nMuon]
   Float_t         b5_Muon_softMva[8];   //[b5_nMuon]
   Float_t         b5_Muon_tkRelIso[8];   //[b5_nMuon]
   Float_t         b5_Muon_tunepRelPt[8];   //[b5_nMuon]
   Float_t         b5_Muon_mvaLowPt[8];   //[b5_nMuon]
   Float_t         b5_Muon_mvaTTH[8];   //[b5_nMuon]
   Int_t           b5_Muon_charge[8];   //[b5_nMuon]
   Int_t           b5_Muon_jetIdx[8];   //[b5_nMuon]
   Int_t           b5_Muon_nStations[8];   //[b5_nMuon]
   Int_t           b5_Muon_nTrackerLayers[8];   //[b5_nMuon]
   Int_t           b5_Muon_pdgId[8];   //[b5_nMuon]
   Int_t           b5_Muon_tightCharge[8];   //[b5_nMuon]
   Int_t           b5_Muon_fsrPhotonIdx[8];   //[b5_nMuon]
   UChar_t         b5_Muon_highPtId[8];   //[b5_nMuon]
   Bool_t          b5_Muon_highPurity[8];   //[b5_nMuon]
   Bool_t          b5_Muon_inTimeMuon[8];   //[b5_nMuon]
   Bool_t          b5_Muon_isGlobal[8];   //[b5_nMuon]
   Bool_t          b5_Muon_isPFcand[8];   //[b5_nMuon]
   Bool_t          b5_Muon_isTracker[8];   //[b5_nMuon]
   UChar_t         b5_Muon_jetNDauCharged[8];   //[b5_nMuon]
   Bool_t          b5_Muon_looseId[8];   //[b5_nMuon]
   Bool_t          b5_Muon_mediumId[8];   //[b5_nMuon]
   Bool_t          b5_Muon_mediumPromptId[8];   //[b5_nMuon]
   UChar_t         b5_Muon_miniIsoId[8];   //[b5_nMuon]
   UChar_t         b5_Muon_multiIsoId[8];   //[b5_nMuon]
   UChar_t         b5_Muon_mvaId[8];   //[b5_nMuon]
   UChar_t         b5_Muon_mvaLowPtId[8];   //[b5_nMuon]
   UChar_t         b5_Muon_pfIsoId[8];   //[b5_nMuon]
   UChar_t         b5_Muon_puppiIsoId[8];   //[b5_nMuon]
   Bool_t          b5_Muon_softId[8];   //[b5_nMuon]
   Bool_t          b5_Muon_softMvaId[8];   //[b5_nMuon]
   Bool_t          b5_Muon_tightId[8];   //[b5_nMuon]
   UChar_t         b5_Muon_tkIsoId[8];   //[b5_nMuon]
   Bool_t          b5_Muon_triggerIdLoose[8];   //[b5_nMuon]
   UInt_t          b5_nPhoton;
   Float_t         b5_Photon_eCorr[7];   //[b5_nPhoton]
   Float_t         b5_Photon_energyErr[7];   //[b5_nPhoton]
   Float_t         b5_Photon_eta[7];   //[b5_nPhoton]
   Float_t         b5_Photon_hoe[7];   //[b5_nPhoton]
   Float_t         b5_Photon_mass[7];   //[b5_nPhoton]
   Float_t         b5_Photon_mvaID[7];   //[b5_nPhoton]
   Float_t         b5_Photon_mvaID_Fall17V1p1[7];   //[b5_nPhoton]
   Float_t         b5_Photon_pfRelIso03_all[7];   //[b5_nPhoton]
   Float_t         b5_Photon_pfRelIso03_chg[7];   //[b5_nPhoton]
   Float_t         b5_Photon_phi[7];   //[b5_nPhoton]
   Float_t         b5_Photon_pt[7];   //[b5_nPhoton]
   Float_t         b5_Photon_r9[7];   //[b5_nPhoton]
   Float_t         b5_Photon_sieie[7];   //[b5_nPhoton]
   Int_t           b5_Photon_charge[7];   //[b5_nPhoton]
   Int_t           b5_Photon_cutBased[7];   //[b5_nPhoton]
   Int_t           b5_Photon_cutBased_Fall17V1Bitmap[7];   //[b5_nPhoton]
   Int_t           b5_Photon_electronIdx[7];   //[b5_nPhoton]
   Int_t           b5_Photon_jetIdx[7];   //[b5_nPhoton]
   Int_t           b5_Photon_pdgId[7];   //[b5_nPhoton]
   Int_t           b5_Photon_vidNestedWPBitmap[7];   //[b5_nPhoton]
   Bool_t          b5_Photon_electronVeto[7];   //[b5_nPhoton]
   Bool_t          b5_Photon_isScEtaEB[7];   //[b5_nPhoton]
   Bool_t          b5_Photon_isScEtaEE[7];   //[b5_nPhoton]
   Bool_t          b5_Photon_mvaID_WP80[7];   //[b5_nPhoton]
   Bool_t          b5_Photon_mvaID_WP90[7];   //[b5_nPhoton]
   Bool_t          b5_Photon_pixelSeed[7];   //[b5_nPhoton]
   UChar_t         b5_Photon_seedGain[7];   //[b5_nPhoton]
   Float_t         b5_PuppiMET_phi;
   Float_t         b5_PuppiMET_phiJERDown;
   Float_t         b5_PuppiMET_phiJERUp;
   Float_t         b5_PuppiMET_phiJESDown;
   Float_t         b5_PuppiMET_phiJESUp;
   Float_t         b5_PuppiMET_phiUnclusteredDown;
   Float_t         b5_PuppiMET_phiUnclusteredUp;
   Float_t         b5_PuppiMET_pt;
   Float_t         b5_PuppiMET_ptJERDown;
   Float_t         b5_PuppiMET_ptJERUp;
   Float_t         b5_PuppiMET_ptJESDown;
   Float_t         b5_PuppiMET_ptJESUp;
   Float_t         b5_PuppiMET_ptUnclusteredDown;
   Float_t         b5_PuppiMET_ptUnclusteredUp;
   Float_t         b5_PuppiMET_sumEt;
   Float_t         b5_RawMET_phi;
   Float_t         b5_RawMET_pt;
   Float_t         b5_RawMET_sumEt;
   Float_t         b5_RawPuppiMET_phi;
   Float_t         b5_RawPuppiMET_pt;
   Float_t         b5_RawPuppiMET_sumEt;
   Float_t         b5_fixedGridRhoFastjetAll;
   Float_t         b5_fixedGridRhoFastjetCentral;
   Float_t         b5_fixedGridRhoFastjetCentralCalo;
   Float_t         b5_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         b5_fixedGridRhoFastjetCentralNeutral;
   UInt_t          b5_nSoftActivityJet;
   Float_t         b5_SoftActivityJet_eta[6];   //[b5_nSoftActivityJet]
   Float_t         b5_SoftActivityJet_phi[6];   //[b5_nSoftActivityJet]
   Float_t         b5_SoftActivityJet_pt[6];   //[b5_nSoftActivityJet]
   Float_t         b5_SoftActivityJetHT;
   Float_t         b5_SoftActivityJetHT10;
   Float_t         b5_SoftActivityJetHT2;
   Float_t         b5_SoftActivityJetHT5;
   Int_t           b5_SoftActivityJetNjets10;
   Int_t           b5_SoftActivityJetNjets2;
   Int_t           b5_SoftActivityJetNjets5;
   UInt_t          b5_nSubJet;
   Float_t         b5_SubJet_btagCMVA[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_btagCSVV2[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_btagDeepB[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_eta[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_mass[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_n2b1[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_n3b1[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_phi[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_pt[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_rawFactor[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_tau1[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_tau2[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_tau3[10];   //[b5_nSubJet]
   Float_t         b5_SubJet_tau4[10];   //[b5_nSubJet]
   UInt_t          b5_nTau;
   Float_t         b5_Tau_chargedIso[5];   //[b5_nTau]
   Float_t         b5_Tau_dxy[5];   //[b5_nTau]
   Float_t         b5_Tau_dz[5];   //[b5_nTau]
   Float_t         b5_Tau_eta[5];   //[b5_nTau]
   Float_t         b5_Tau_leadTkDeltaEta[5];   //[b5_nTau]
   Float_t         b5_Tau_leadTkDeltaPhi[5];   //[b5_nTau]
   Float_t         b5_Tau_leadTkPtOverTauPt[5];   //[b5_nTau]
   Float_t         b5_Tau_mass[5];   //[b5_nTau]
   Float_t         b5_Tau_neutralIso[5];   //[b5_nTau]
   Float_t         b5_Tau_phi[5];   //[b5_nTau]
   Float_t         b5_Tau_photonsOutsideSignalCone[5];   //[b5_nTau]
   Float_t         b5_Tau_pt[5];   //[b5_nTau]
   Float_t         b5_Tau_puCorr[5];   //[b5_nTau]
   Float_t         b5_Tau_rawAntiEle[5];   //[b5_nTau]
   Float_t         b5_Tau_rawAntiEle2018[5];   //[b5_nTau]
   Float_t         b5_Tau_rawDeepTau2017v2p1VSe[5];   //[b5_nTau]
   Float_t         b5_Tau_rawDeepTau2017v2p1VSjet[5];   //[b5_nTau]
   Float_t         b5_Tau_rawDeepTau2017v2p1VSmu[5];   //[b5_nTau]
   Float_t         b5_Tau_rawIso[5];   //[b5_nTau]
   Float_t         b5_Tau_rawIsodR03[5];   //[b5_nTau]
   Float_t         b5_Tau_rawMVAnewDM2017v2[5];   //[b5_nTau]
   Float_t         b5_Tau_rawMVAoldDM[5];   //[b5_nTau]
   Float_t         b5_Tau_rawMVAoldDM2017v1[5];   //[b5_nTau]
   Float_t         b5_Tau_rawMVAoldDM2017v2[5];   //[b5_nTau]
   Float_t         b5_Tau_rawMVAoldDMdR032017v2[5];   //[b5_nTau]
   Int_t           b5_Tau_charge[5];   //[b5_nTau]
   Int_t           b5_Tau_decayMode[5];   //[b5_nTau]
   Int_t           b5_Tau_jetIdx[5];   //[b5_nTau]
   Int_t           b5_Tau_rawAntiEleCat[5];   //[b5_nTau]
   Int_t           b5_Tau_rawAntiEleCat2018[5];   //[b5_nTau]
   UChar_t         b5_Tau_idAntiEle[5];   //[b5_nTau]
   UChar_t         b5_Tau_idAntiEle2018[5];   //[b5_nTau]
   Bool_t          b5_Tau_idAntiEleDeadECal[5];   //[b5_nTau]
   UChar_t         b5_Tau_idAntiMu[5];   //[b5_nTau]
   Bool_t          b5_Tau_idDecayMode[5];   //[b5_nTau]
   Bool_t          b5_Tau_idDecayModeNewDMs[5];   //[b5_nTau]
   UChar_t         b5_Tau_idDeepTau2017v2p1VSe[5];   //[b5_nTau]
   UChar_t         b5_Tau_idDeepTau2017v2p1VSjet[5];   //[b5_nTau]
   UChar_t         b5_Tau_idDeepTau2017v2p1VSmu[5];   //[b5_nTau]
   UChar_t         b5_Tau_idMVAnewDM2017v2[5];   //[b5_nTau]
   UChar_t         b5_Tau_idMVAoldDM[5];   //[b5_nTau]
   UChar_t         b5_Tau_idMVAoldDM2017v1[5];   //[b5_nTau]
   UChar_t         b5_Tau_idMVAoldDM2017v2[5];   //[b5_nTau]
   UChar_t         b5_Tau_idMVAoldDMdR032017v2[5];   //[b5_nTau]
   Float_t         b5_TkMET_phi;
   Float_t         b5_TkMET_pt;
   Float_t         b5_TkMET_sumEt;
   UInt_t          b5_nTrigObj;
   Float_t         b5_TrigObj_pt[34];   //[b5_nTrigObj]
   Float_t         b5_TrigObj_eta[34];   //[b5_nTrigObj]
   Float_t         b5_TrigObj_phi[34];   //[b5_nTrigObj]
   Float_t         b5_TrigObj_l1pt[34];   //[b5_nTrigObj]
   Float_t         b5_TrigObj_l1pt_2[34];   //[b5_nTrigObj]
   Float_t         b5_TrigObj_l2pt[34];   //[b5_nTrigObj]
   Int_t           b5_TrigObj_id[34];   //[b5_nTrigObj]
   Int_t           b5_TrigObj_l1iso[34];   //[b5_nTrigObj]
   Int_t           b5_TrigObj_l1charge[34];   //[b5_nTrigObj]
   Int_t           b5_TrigObj_filterBits[34];   //[b5_nTrigObj]
   UInt_t          b5_nOtherPV;
   Float_t         b5_OtherPV_z[3];   //[b5_nOtherPV]
   Float_t         b5_PV_ndof;
   Float_t         b5_PV_x;
   Float_t         b5_PV_y;
   Float_t         b5_PV_z;
   Float_t         b5_PV_chi2;
   Float_t         b5_PV_score;
   Int_t           b5_PV_npvs;
   Int_t           b5_PV_npvsGood;
   UInt_t          b5_nSV;
   Float_t         b5_SV_dlen[13];   //[b5_nSV]
   Float_t         b5_SV_dlenSig[13];   //[b5_nSV]
   Float_t         b5_SV_dxy[13];   //[b5_nSV]
   Float_t         b5_SV_dxySig[13];   //[b5_nSV]
   Float_t         b5_SV_pAngle[13];   //[b5_nSV]
   UChar_t         b5_Electron_cleanmask[6];   //[b5_nElectron]
   UChar_t         b5_Jet_cleanmask[21];   //[b5_nJet]
   UChar_t         b5_Muon_cleanmask[8];   //[b5_nMuon]
   UChar_t         b5_Photon_cleanmask[7];   //[b5_nPhoton]
   UChar_t         b5_Tau_cleanmask[5];   //[b5_nTau]
   Float_t         b5_SV_chi2[13];   //[b5_nSV]
   Float_t         b5_SV_eta[13];   //[b5_nSV]
   Float_t         b5_SV_mass[13];   //[b5_nSV]
   Float_t         b5_SV_ndof[13];   //[b5_nSV]
   Float_t         b5_SV_phi[13];   //[b5_nSV]
   Float_t         b5_SV_pt[13];   //[b5_nSV]
   Float_t         b5_SV_x[13];   //[b5_nSV]
   Float_t         b5_SV_y[13];   //[b5_nSV]
   Float_t         b5_SV_z[13];   //[b5_nSV]
   UChar_t         b5_SV_ntracks[13];   //[b5_nSV]
   Bool_t          b5_L1_AlwaysTrue;
   Bool_t          b5_L1_BPTX_AND_Ref1_VME;
   Bool_t          b5_L1_BPTX_AND_Ref3_VME;
   Bool_t          b5_L1_BPTX_AND_Ref4_VME;
   Bool_t          b5_L1_BPTX_BeamGas_B1_VME;
   Bool_t          b5_L1_BPTX_BeamGas_B2_VME;
   Bool_t          b5_L1_BPTX_BeamGas_Ref1_VME;
   Bool_t          b5_L1_BPTX_BeamGas_Ref2_VME;
   Bool_t          b5_L1_BPTX_NotOR_VME;
   Bool_t          b5_L1_BPTX_OR_Ref3_VME;
   Bool_t          b5_L1_BPTX_OR_Ref4_VME;
   Bool_t          b5_L1_BPTX_RefAND_VME;
   Bool_t          b5_L1_BptxMinus;
   Bool_t          b5_L1_BptxOR;
   Bool_t          b5_L1_BptxPlus;
   Bool_t          b5_L1_BptxXOR;
   Bool_t          b5_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          b5_L1_DoubleEG8er2p5_HTT260er;
   Bool_t          b5_L1_DoubleEG8er2p5_HTT280er;
   Bool_t          b5_L1_DoubleEG8er2p5_HTT300er;
   Bool_t          b5_L1_DoubleEG8er2p5_HTT320er;
   Bool_t          b5_L1_DoubleEG8er2p5_HTT340er;
   Bool_t          b5_L1_DoubleEG_15_10_er2p5;
   Bool_t          b5_L1_DoubleEG_20_10_er2p5;
   Bool_t          b5_L1_DoubleEG_22_10_er2p5;
   Bool_t          b5_L1_DoubleEG_25_12_er2p5;
   Bool_t          b5_L1_DoubleEG_25_14_er2p5;
   Bool_t          b5_L1_DoubleEG_27_14_er2p5;
   Bool_t          b5_L1_DoubleEG_LooseIso20_10_er2p5;
   Bool_t          b5_L1_DoubleEG_LooseIso22_10_er2p5;
   Bool_t          b5_L1_DoubleEG_LooseIso22_12_er2p5;
   Bool_t          b5_L1_DoubleEG_LooseIso25_12_er2p5;
   Bool_t          b5_L1_DoubleIsoTau32er2p1;
   Bool_t          b5_L1_DoubleIsoTau34er2p1;
   Bool_t          b5_L1_DoubleIsoTau36er2p1;
   Bool_t          b5_L1_DoubleJet100er2p3_dEta_Max1p6;
   Bool_t          b5_L1_DoubleJet100er2p5;
   Bool_t          b5_L1_DoubleJet112er2p3_dEta_Max1p6;
   Bool_t          b5_L1_DoubleJet120er2p5;
   Bool_t          b5_L1_DoubleJet150er2p5;
   Bool_t          b5_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5;
   Bool_t          b5_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5;
   Bool_t          b5_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5;
   Bool_t          b5_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5;
   Bool_t          b5_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5;
   Bool_t          b5_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5;
   Bool_t          b5_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp;
   Bool_t          b5_L1_DoubleJet40er2p5;
   Bool_t          b5_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;
   Bool_t          b5_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;
   Bool_t          b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;
   Bool_t          b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28;
   Bool_t          b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620;
   Bool_t          b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28;
   Bool_t          b5_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ;
   Bool_t          b5_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp;
   Bool_t          b5_L1_DoubleJet_80_30_Mass_Min420_Mu8;
   Bool_t          b5_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;
   Bool_t          b5_L1_DoubleLooseIsoEG22er2p1;
   Bool_t          b5_L1_DoubleLooseIsoEG24er2p1;
   Bool_t          b5_L1_DoubleMu0;
   Bool_t          b5_L1_DoubleMu0_Mass_Min1;
   Bool_t          b5_L1_DoubleMu0_OQ;
   Bool_t          b5_L1_DoubleMu0_SQ;
   Bool_t          b5_L1_DoubleMu0_SQ_OS;
   Bool_t          b5_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8;
   Bool_t          b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
   Bool_t          b5_L1_DoubleMu0er1p5_SQ;
   Bool_t          b5_L1_DoubleMu0er1p5_SQ_OS;
   Bool_t          b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
   Bool_t          b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4;
   Bool_t          b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4;
   Bool_t          b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4;
   Bool_t          b5_L1_DoubleMu10_SQ;
   Bool_t          b5_L1_DoubleMu18er2p1;
   Bool_t          b5_L1_DoubleMu3_OS_DoubleEG7p5Upsilon;
   Bool_t          b5_L1_DoubleMu3_SQ_ETMHF50_HTT60er;
   Bool_t          b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5;
   Bool_t          b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t          b5_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5;
   Bool_t          b5_L1_DoubleMu3_SQ_HTT220er;
   Bool_t          b5_L1_DoubleMu3_SQ_HTT240er;
   Bool_t          b5_L1_DoubleMu3_SQ_HTT260er;
   Bool_t          b5_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8;
   Bool_t          b5_L1_DoubleMu4_SQ_EG9er2p5;
   Bool_t          b5_L1_DoubleMu4_SQ_OS;
   Bool_t          b5_L1_DoubleMu4_SQ_OS_dR_Max1p2;
   Bool_t          b5_L1_DoubleMu4p5_SQ_OS;
   Bool_t          b5_L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
   Bool_t          b5_L1_DoubleMu4p5er2p0_SQ_OS;
   Bool_t          b5_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;
   Bool_t          b5_L1_DoubleMu5Upsilon_OS_DoubleEG3;
   Bool_t          b5_L1_DoubleMu5_SQ_EG9er2p5;
   Bool_t          b5_L1_DoubleMu9_SQ;
   Bool_t          b5_L1_DoubleMu_12_5;
   Bool_t          b5_L1_DoubleMu_15_5_SQ;
   Bool_t          b5_L1_DoubleMu_15_7;
   Bool_t          b5_L1_DoubleMu_15_7_Mass_Min1;
   Bool_t          b5_L1_DoubleMu_15_7_SQ;
   Bool_t          b5_L1_DoubleTau70er2p1;
   Bool_t          b5_L1_ETM120;
   Bool_t          b5_L1_ETM150;
   Bool_t          b5_L1_ETMHF100;
   Bool_t          b5_L1_ETMHF100_HTT60er;
   Bool_t          b5_L1_ETMHF110;
   Bool_t          b5_L1_ETMHF110_HTT60er;
   Bool_t          b5_L1_ETMHF110_HTT60er_NotSecondBunchInTrain;
   Bool_t          b5_L1_ETMHF120;
   Bool_t          b5_L1_ETMHF120_HTT60er;
   Bool_t          b5_L1_ETMHF120_NotSecondBunchInTrain;
   Bool_t          b5_L1_ETMHF130;
   Bool_t          b5_L1_ETMHF130_HTT60er;
   Bool_t          b5_L1_ETMHF140;
   Bool_t          b5_L1_ETMHF150;
   Bool_t          b5_L1_ETMHF90_HTT60er;
   Bool_t          b5_L1_ETT1200;
   Bool_t          b5_L1_ETT1600;
   Bool_t          b5_L1_ETT2000;
   Bool_t          b5_L1_FirstBunchAfterTrain;
   Bool_t          b5_L1_FirstBunchBeforeTrain;
   Bool_t          b5_L1_FirstBunchInTrain;
   Bool_t          b5_L1_FirstCollisionInOrbit;
   Bool_t          b5_L1_FirstCollisionInTrain;
   Bool_t          b5_L1_HCAL_LaserMon_Trig;
   Bool_t          b5_L1_HCAL_LaserMon_Veto;
   Bool_t          b5_L1_HTT120er;
   Bool_t          b5_L1_HTT160er;
   Bool_t          b5_L1_HTT200er;
   Bool_t          b5_L1_HTT255er;
   Bool_t          b5_L1_HTT280er;
   Bool_t          b5_L1_HTT280er_QuadJet_70_55_40_35_er2p4;
   Bool_t          b5_L1_HTT320er;
   Bool_t          b5_L1_HTT320er_QuadJet_70_55_40_40_er2p4;
   Bool_t          b5_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;
   Bool_t          b5_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;
   Bool_t          b5_L1_HTT360er;
   Bool_t          b5_L1_HTT400er;
   Bool_t          b5_L1_HTT450er;
   Bool_t          b5_L1_IsoEG32er2p5_Mt40;
   Bool_t          b5_L1_IsoEG32er2p5_Mt44;
   Bool_t          b5_L1_IsoEG32er2p5_Mt48;
   Bool_t          b5_L1_IsoTau40er2p1_ETMHF100;
   Bool_t          b5_L1_IsoTau40er2p1_ETMHF110;
   Bool_t          b5_L1_IsoTau40er2p1_ETMHF120;
   Bool_t          b5_L1_IsoTau40er2p1_ETMHF90;
   Bool_t          b5_L1_IsolatedBunch;
   Bool_t          b5_L1_LastBunchInTrain;
   Bool_t          b5_L1_LastCollisionInTrain;
   Bool_t          b5_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;
   Bool_t          b5_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3;
   Bool_t          b5_L1_LooseIsoEG24er2p1_HTT100er;
   Bool_t          b5_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;
   Bool_t          b5_L1_LooseIsoEG26er2p1_HTT100er;
   Bool_t          b5_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          b5_L1_LooseIsoEG28er2p1_HTT100er;
   Bool_t          b5_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          b5_L1_LooseIsoEG30er2p1_HTT100er;
   Bool_t          b5_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          b5_L1_MinimumBiasHF0_AND_BptxAND;
   Bool_t          b5_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;
   Bool_t          b5_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6;
   Bool_t          b5_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;
   Bool_t          b5_L1_Mu18er2p1_Tau24er2p1;
   Bool_t          b5_L1_Mu18er2p1_Tau26er2p1;
   Bool_t          b5_L1_Mu20_EG10er2p5;
   Bool_t          b5_L1_Mu22er2p1_IsoTau32er2p1;
   Bool_t          b5_L1_Mu22er2p1_IsoTau34er2p1;
   Bool_t          b5_L1_Mu22er2p1_IsoTau36er2p1;
   Bool_t          b5_L1_Mu22er2p1_IsoTau40er2p1;
   Bool_t          b5_L1_Mu22er2p1_Tau70er2p1;
   Bool_t          b5_L1_Mu3_Jet120er2p5_dR_Max0p4;
   Bool_t          b5_L1_Mu3_Jet120er2p5_dR_Max0p8;
   Bool_t          b5_L1_Mu3_Jet16er2p5_dR_Max0p4;
   Bool_t          b5_L1_Mu3_Jet30er2p5;
   Bool_t          b5_L1_Mu3_Jet35er2p5_dR_Max0p4;
   Bool_t          b5_L1_Mu3_Jet60er2p5_dR_Max0p4;
   Bool_t          b5_L1_Mu3_Jet80er2p5_dR_Max0p4;
   Bool_t          b5_L1_Mu3er1p5_Jet100er2p5_ETMHF40;
   Bool_t          b5_L1_Mu3er1p5_Jet100er2p5_ETMHF50;
   Bool_t          b5_L1_Mu5_EG23er2p5;
   Bool_t          b5_L1_Mu5_LooseIsoEG20er2p5;
   Bool_t          b5_L1_Mu6_DoubleEG10er2p5;
   Bool_t          b5_L1_Mu6_DoubleEG12er2p5;
   Bool_t          b5_L1_Mu6_DoubleEG15er2p5;
   Bool_t          b5_L1_Mu6_DoubleEG17er2p5;
   Bool_t          b5_L1_Mu6_HTT240er;
   Bool_t          b5_L1_Mu6_HTT250er;
   Bool_t          b5_L1_Mu7_EG23er2p5;
   Bool_t          b5_L1_Mu7_LooseIsoEG20er2p5;
   Bool_t          b5_L1_Mu7_LooseIsoEG23er2p5;
   Bool_t          b5_L1_NotBptxOR;
   Bool_t          b5_L1_QuadJet36er2p5_IsoTau52er2p1;
   Bool_t          b5_L1_QuadJet60er2p5;
   Bool_t          b5_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0;
   Bool_t          b5_L1_QuadMu0;
   Bool_t          b5_L1_QuadMu0_OQ;
   Bool_t          b5_L1_QuadMu0_SQ;
   Bool_t          b5_L1_SecondBunchInTrain;
   Bool_t          b5_L1_SecondLastBunchInTrain;
   Bool_t          b5_L1_SingleEG10er2p5;
   Bool_t          b5_L1_SingleEG15er2p5;
   Bool_t          b5_L1_SingleEG26er2p5;
   Bool_t          b5_L1_SingleEG34er2p5;
   Bool_t          b5_L1_SingleEG36er2p5;
   Bool_t          b5_L1_SingleEG38er2p5;
   Bool_t          b5_L1_SingleEG40er2p5;
   Bool_t          b5_L1_SingleEG42er2p5;
   Bool_t          b5_L1_SingleEG45er2p5;
   Bool_t          b5_L1_SingleEG50;
   Bool_t          b5_L1_SingleEG60;
   Bool_t          b5_L1_SingleEG8er2p5;
   Bool_t          b5_L1_SingleIsoEG24er1p5;
   Bool_t          b5_L1_SingleIsoEG24er2p1;
   Bool_t          b5_L1_SingleIsoEG26er1p5;
   Bool_t          b5_L1_SingleIsoEG26er2p1;
   Bool_t          b5_L1_SingleIsoEG26er2p5;
   Bool_t          b5_L1_SingleIsoEG28er1p5;
   Bool_t          b5_L1_SingleIsoEG28er2p1;
   Bool_t          b5_L1_SingleIsoEG28er2p5;
   Bool_t          b5_L1_SingleIsoEG30er2p1;
   Bool_t          b5_L1_SingleIsoEG30er2p5;
   Bool_t          b5_L1_SingleIsoEG32er2p1;
   Bool_t          b5_L1_SingleIsoEG32er2p5;
   Bool_t          b5_L1_SingleIsoEG34er2p5;
   Bool_t          b5_L1_SingleJet10erHE;
   Bool_t          b5_L1_SingleJet120;
   Bool_t          b5_L1_SingleJet120_FWD3p0;
   Bool_t          b5_L1_SingleJet120er2p5;
   Bool_t          b5_L1_SingleJet12erHE;
   Bool_t          b5_L1_SingleJet140er2p5;
   Bool_t          b5_L1_SingleJet140er2p5_ETMHF80;
   Bool_t          b5_L1_SingleJet140er2p5_ETMHF90;
   Bool_t          b5_L1_SingleJet160er2p5;
   Bool_t          b5_L1_SingleJet180;
   Bool_t          b5_L1_SingleJet180er2p5;
   Bool_t          b5_L1_SingleJet200;
   Bool_t          b5_L1_SingleJet20er2p5_NotBptxOR;
   Bool_t          b5_L1_SingleJet20er2p5_NotBptxOR_3BX;
   Bool_t          b5_L1_SingleJet35;
   Bool_t          b5_L1_SingleJet35_FWD3p0;
   Bool_t          b5_L1_SingleJet35er2p5;
   Bool_t          b5_L1_SingleJet43er2p5_NotBptxOR_3BX;
   Bool_t          b5_L1_SingleJet46er2p5_NotBptxOR_3BX;
   Bool_t          b5_L1_SingleJet60;
   Bool_t          b5_L1_SingleJet60_FWD3p0;
   Bool_t          b5_L1_SingleJet60er2p5;
   Bool_t          b5_L1_SingleJet8erHE;
   Bool_t          b5_L1_SingleJet90;
   Bool_t          b5_L1_SingleJet90_FWD3p0;
   Bool_t          b5_L1_SingleJet90er2p5;
   Bool_t          b5_L1_SingleLooseIsoEG28er1p5;
   Bool_t          b5_L1_SingleLooseIsoEG30er1p5;
   Bool_t          b5_L1_SingleMu0_BMTF;
   Bool_t          b5_L1_SingleMu0_DQ;
   Bool_t          b5_L1_SingleMu0_EMTF;
   Bool_t          b5_L1_SingleMu0_OMTF;
   Bool_t          b5_L1_SingleMu10er1p5;
   Bool_t          b5_L1_SingleMu12_DQ_BMTF;
   Bool_t          b5_L1_SingleMu12_DQ_EMTF;
   Bool_t          b5_L1_SingleMu12_DQ_OMTF;
   Bool_t          b5_L1_SingleMu12er1p5;
   Bool_t          b5_L1_SingleMu14er1p5;
   Bool_t          b5_L1_SingleMu15_DQ;
   Bool_t          b5_L1_SingleMu16er1p5;
   Bool_t          b5_L1_SingleMu18;
   Bool_t          b5_L1_SingleMu18er1p5;
   Bool_t          b5_L1_SingleMu20;
   Bool_t          b5_L1_SingleMu22;
   Bool_t          b5_L1_SingleMu22_BMTF;
   Bool_t          b5_L1_SingleMu22_EMTF;
   Bool_t          b5_L1_SingleMu22_OMTF;
   Bool_t          b5_L1_SingleMu25;
   Bool_t          b5_L1_SingleMu3;
   Bool_t          b5_L1_SingleMu5;
   Bool_t          b5_L1_SingleMu6er1p5;
   Bool_t          b5_L1_SingleMu7;
   Bool_t          b5_L1_SingleMu7_DQ;
   Bool_t          b5_L1_SingleMu7er1p5;
   Bool_t          b5_L1_SingleMu8er1p5;
   Bool_t          b5_L1_SingleMu9er1p5;
   Bool_t          b5_L1_SingleMuCosmics;
   Bool_t          b5_L1_SingleMuCosmics_BMTF;
   Bool_t          b5_L1_SingleMuCosmics_EMTF;
   Bool_t          b5_L1_SingleMuCosmics_OMTF;
   Bool_t          b5_L1_SingleMuOpen;
   Bool_t          b5_L1_SingleMuOpen_NotBptxOR;
   Bool_t          b5_L1_SingleMuOpen_er1p1_NotBptxOR_3BX;
   Bool_t          b5_L1_SingleMuOpen_er1p4_NotBptxOR_3BX;
   Bool_t          b5_L1_SingleTau120er2p1;
   Bool_t          b5_L1_SingleTau130er2p1;
   Bool_t          b5_L1_TOTEM_1;
   Bool_t          b5_L1_TOTEM_2;
   Bool_t          b5_L1_TOTEM_3;
   Bool_t          b5_L1_TOTEM_4;
   Bool_t          b5_L1_TripleEG16er2p5;
   Bool_t          b5_L1_TripleEG_16_12_8_er2p5;
   Bool_t          b5_L1_TripleEG_16_15_8_er2p5;
   Bool_t          b5_L1_TripleEG_18_17_8_er2p5;
   Bool_t          b5_L1_TripleEG_18_18_12_er2p5;
   Bool_t          b5_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5;
   Bool_t          b5_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5;
   Bool_t          b5_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5;
   Bool_t          b5_L1_TripleMu0;
   Bool_t          b5_L1_TripleMu0_OQ;
   Bool_t          b5_L1_TripleMu0_SQ;
   Bool_t          b5_L1_TripleMu3;
   Bool_t          b5_L1_TripleMu3_SQ;
   Bool_t          b5_L1_TripleMu_5SQ_3SQ_0OQ;
   Bool_t          b5_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          b5_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          b5_L1_TripleMu_5_3_3;
   Bool_t          b5_L1_TripleMu_5_3_3_SQ;
   Bool_t          b5_L1_TripleMu_5_3p5_2p5;
   Bool_t          b5_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          b5_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;
   Bool_t          b5_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          b5_L1_TripleMu_5_5_3;
   Bool_t          b5_L1_UnpairedBunchBptxMinus;
   Bool_t          b5_L1_UnpairedBunchBptxPlus;
   Bool_t          b5_L1_ZeroBias;
   Bool_t          b5_L1_ZeroBias_copy;
   Bool_t          b5_L1_UnprefireableEvent;
   Bool_t          b5_Flag_HBHENoiseFilter;
   Bool_t          b5_Flag_HBHENoiseIsoFilter;
   Bool_t          b5_Flag_CSCTightHaloFilter;
   Bool_t          b5_Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          b5_Flag_CSCTightHalo2015Filter;
   Bool_t          b5_Flag_globalTightHalo2016Filter;
   Bool_t          b5_Flag_globalSuperTightHalo2016Filter;
   Bool_t          b5_Flag_HcalStripHaloFilter;
   Bool_t          b5_Flag_hcalLaserEventFilter;
   Bool_t          b5_Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          b5_Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          b5_Flag_ecalBadCalibFilter;
   Bool_t          b5_Flag_goodVertices;
   Bool_t          b5_Flag_eeBadScFilter;
   Bool_t          b5_Flag_ecalLaserCorrFilter;
   Bool_t          b5_Flag_trkPOGFilters;
   Bool_t          b5_Flag_chargedHadronTrackResolutionFilter;
   Bool_t          b5_Flag_muonBadTrackFilter;
   Bool_t          b5_Flag_BadChargedCandidateFilter;
   Bool_t          b5_Flag_BadPFMuonFilter;
   Bool_t          b5_Flag_BadPFMuonDzFilter;
   Bool_t          b5_Flag_hfNoisyHitsFilter;
   Bool_t          b5_Flag_BadChargedCandidateSummer16Filter;
   Bool_t          b5_Flag_BadPFMuonSummer16Filter;
   Bool_t          b5_Flag_trkPOG_manystripclus53X;
   Bool_t          b5_Flag_trkPOG_toomanystripclus53X;
   Bool_t          b5_Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          b5_Flag_METFilters;
   Bool_t          b5_L1Reco_step;
   Bool_t          b5_HLTriggerFirstPath;
   Bool_t          b5_HLT_AK8PFJet360_TrimMass30;
   Bool_t          b5_HLT_AK8PFJet380_TrimMass30;
   Bool_t          b5_HLT_AK8PFJet400_TrimMass30;
   Bool_t          b5_HLT_AK8PFJet420_TrimMass30;
   Bool_t          b5_HLT_AK8PFHT750_TrimMass50;
   Bool_t          b5_HLT_AK8PFHT800_TrimMass50;
   Bool_t          b5_HLT_AK8PFHT850_TrimMass50;
   Bool_t          b5_HLT_AK8PFHT900_TrimMass50;
   Bool_t          b5_HLT_CaloJet500_NoJetID;
   Bool_t          b5_HLT_CaloJet550_NoJetID;
   Bool_t          b5_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;
   Bool_t          b5_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;
   Bool_t          b5_HLT_Trimuon5_3p5_2_Upsilon_Muon;
   Bool_t          b5_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;
   Bool_t          b5_HLT_DoubleEle25_CaloIdL_MW;
   Bool_t          b5_HLT_DoubleEle27_CaloIdL_MW;
   Bool_t          b5_HLT_DoubleEle33_CaloIdL_MW;
   Bool_t          b5_HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   Bool_t          b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   Bool_t          b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          b5_HLT_Ele27_Ele37_CaloIdL_MW;
   Bool_t          b5_HLT_Mu27_Ele37_CaloIdL_MW;
   Bool_t          b5_HLT_Mu37_Ele27_CaloIdL_MW;
   Bool_t          b5_HLT_Mu37_TkMu27;
   Bool_t          b5_HLT_DoubleMu4_3_Bs;
   Bool_t          b5_HLT_DoubleMu4_3_Jpsi;
   Bool_t          b5_HLT_DoubleMu4_JpsiTrk_Displaced;
   Bool_t          b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;
   Bool_t          b5_HLT_DoubleMu3_Trk_Tau3mu;
   Bool_t          b5_HLT_DoubleMu3_TkMu_DsTau3Mu;
   Bool_t          b5_HLT_DoubleMu4_PsiPrimeTrk_Displaced;
   Bool_t          b5_HLT_DoubleMu4_Mass3p8_DZ_PFHT350;
   Bool_t          b5_HLT_Mu3_PFJet40;
   Bool_t          b5_HLT_Mu7p5_L2Mu2_Jpsi;
   Bool_t          b5_HLT_Mu7p5_L2Mu2_Upsilon;
   Bool_t          b5_HLT_Mu7p5_Track2_Jpsi;
   Bool_t          b5_HLT_Mu7p5_Track3p5_Jpsi;
   Bool_t          b5_HLT_Mu7p5_Track7_Jpsi;
   Bool_t          b5_HLT_Mu7p5_Track2_Upsilon;
   Bool_t          b5_HLT_Mu7p5_Track3p5_Upsilon;
   Bool_t          b5_HLT_Mu7p5_Track7_Upsilon;
   Bool_t          b5_HLT_DoublePhoton33_CaloIdL;
   Bool_t          b5_HLT_DoublePhoton70;
   Bool_t          b5_HLT_DoublePhoton85;
   Bool_t          b5_HLT_Ele20_WPTight_Gsf;
   Bool_t          b5_HLT_Ele15_WPLoose_Gsf;
   Bool_t          b5_HLT_Ele17_WPLoose_Gsf;
   Bool_t          b5_HLT_Ele20_WPLoose_Gsf;
   Bool_t          b5_HLT_Ele20_eta2p1_WPLoose_Gsf;
   Bool_t          b5_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   Bool_t          b5_HLT_Ele27_WPTight_Gsf;
   Bool_t          b5_HLT_Ele32_WPTight_Gsf;
   Bool_t          b5_HLT_Ele35_WPTight_Gsf;
   Bool_t          b5_HLT_Ele35_WPTight_Gsf_L1EGMT;
   Bool_t          b5_HLT_Ele38_WPTight_Gsf;
   Bool_t          b5_HLT_Ele40_WPTight_Gsf;
   Bool_t          b5_HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t          b5_HLT_HT450_Beamspot;
   Bool_t          b5_HLT_HT300_Beamspot;
   Bool_t          b5_HLT_ZeroBias_Beamspot;
   Bool_t          b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1;
   Bool_t          b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1;
   Bool_t          b5_HLT_IsoMu20;
   Bool_t          b5_HLT_IsoMu24;
   Bool_t          b5_HLT_IsoMu24_eta2p1;
   Bool_t          b5_HLT_IsoMu27;
   Bool_t          b5_HLT_IsoMu30;
   Bool_t          b5_HLT_UncorrectedJetE30_NoBPTX;
   Bool_t          b5_HLT_UncorrectedJetE30_NoBPTX3BX;
   Bool_t          b5_HLT_UncorrectedJetE60_NoBPTX3BX;
   Bool_t          b5_HLT_UncorrectedJetE70_NoBPTX3BX;
   Bool_t          b5_HLT_L1SingleMu18;
   Bool_t          b5_HLT_L1SingleMu25;
   Bool_t          b5_HLT_L2Mu10;
   Bool_t          b5_HLT_L2Mu10_NoVertex_NoBPTX3BX;
   Bool_t          b5_HLT_L2Mu10_NoVertex_NoBPTX;
   Bool_t          b5_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          b5_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          b5_HLT_L2Mu50;
   Bool_t          b5_HLT_L2Mu23NoVtx_2Cha;
   Bool_t          b5_HLT_L2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t          b5_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t          b5_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;
   Bool_t          b5_HLT_DoubleL2Mu50;
   Bool_t          b5_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t          b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;
   Bool_t          b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t          b5_HLT_DoubleL2Mu23NoVtx_2Cha;
   Bool_t          b5_HLT_DoubleL2Mu25NoVtx_2Cha;
   Bool_t          b5_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;
   Bool_t          b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          b5_HLT_Mu25_TkMu0_Onia;
   Bool_t          b5_HLT_Mu30_TkMu0_Onia;
   Bool_t          b5_HLT_Mu20_TkMu0_Phi;
   Bool_t          b5_HLT_Mu25_TkMu0_Phi;
   Bool_t          b5_HLT_Mu12;
   Bool_t          b5_HLT_Mu15;
   Bool_t          b5_HLT_Mu20;
   Bool_t          b5_HLT_Mu27;
   Bool_t          b5_HLT_Mu50;
   Bool_t          b5_HLT_Mu55;
   Bool_t          b5_HLT_OldMu100;
   Bool_t          b5_HLT_TkMu100;
   Bool_t          b5_HLT_DiPFJetAve40;
   Bool_t          b5_HLT_DiPFJetAve60;
   Bool_t          b5_HLT_DiPFJetAve80;
   Bool_t          b5_HLT_DiPFJetAve140;
   Bool_t          b5_HLT_DiPFJetAve200;
   Bool_t          b5_HLT_DiPFJetAve260;
   Bool_t          b5_HLT_DiPFJetAve320;
   Bool_t          b5_HLT_DiPFJetAve400;
   Bool_t          b5_HLT_DiPFJetAve500;
   Bool_t          b5_HLT_DiPFJetAve15_HFJEC;
   Bool_t          b5_HLT_DiPFJetAve25_HFJEC;
   Bool_t          b5_HLT_DiPFJetAve60_HFJEC;
   Bool_t          b5_HLT_DiPFJetAve80_HFJEC;
   Bool_t          b5_HLT_DiPFJetAve100_HFJEC;
   Bool_t          b5_HLT_DiPFJetAve160_HFJEC;
   Bool_t          b5_HLT_DiPFJetAve220_HFJEC;
   Bool_t          b5_HLT_DiPFJetAve300_HFJEC;
   Bool_t          b5_HLT_AK8PFJet15;
   Bool_t          b5_HLT_AK8PFJet25;
   Bool_t          b5_HLT_AK8PFJet40;
   Bool_t          b5_HLT_AK8PFJet60;
   Bool_t          b5_HLT_AK8PFJet80;
   Bool_t          b5_HLT_AK8PFJet140;
   Bool_t          b5_HLT_AK8PFJet200;
   Bool_t          b5_HLT_AK8PFJet260;
   Bool_t          b5_HLT_AK8PFJet320;
   Bool_t          b5_HLT_AK8PFJet400;
   Bool_t          b5_HLT_AK8PFJet450;
   Bool_t          b5_HLT_AK8PFJet500;
   Bool_t          b5_HLT_AK8PFJet550;
   Bool_t          b5_HLT_PFJet15;
   Bool_t          b5_HLT_PFJet25;
   Bool_t          b5_HLT_PFJet40;
   Bool_t          b5_HLT_PFJet60;
   Bool_t          b5_HLT_PFJet80;
   Bool_t          b5_HLT_PFJet140;
   Bool_t          b5_HLT_PFJet200;
   Bool_t          b5_HLT_PFJet260;
   Bool_t          b5_HLT_PFJet320;
   Bool_t          b5_HLT_PFJet400;
   Bool_t          b5_HLT_PFJet450;
   Bool_t          b5_HLT_PFJet500;
   Bool_t          b5_HLT_PFJet550;
   Bool_t          b5_HLT_PFJetFwd15;
   Bool_t          b5_HLT_PFJetFwd25;
   Bool_t          b5_HLT_PFJetFwd40;
   Bool_t          b5_HLT_PFJetFwd60;
   Bool_t          b5_HLT_PFJetFwd80;
   Bool_t          b5_HLT_PFJetFwd140;
   Bool_t          b5_HLT_PFJetFwd200;
   Bool_t          b5_HLT_PFJetFwd260;
   Bool_t          b5_HLT_PFJetFwd320;
   Bool_t          b5_HLT_PFJetFwd400;
   Bool_t          b5_HLT_PFJetFwd450;
   Bool_t          b5_HLT_PFJetFwd500;
   Bool_t          b5_HLT_AK8PFJetFwd15;
   Bool_t          b5_HLT_AK8PFJetFwd25;
   Bool_t          b5_HLT_AK8PFJetFwd40;
   Bool_t          b5_HLT_AK8PFJetFwd60;
   Bool_t          b5_HLT_AK8PFJetFwd80;
   Bool_t          b5_HLT_AK8PFJetFwd140;
   Bool_t          b5_HLT_AK8PFJetFwd200;
   Bool_t          b5_HLT_AK8PFJetFwd260;
   Bool_t          b5_HLT_AK8PFJetFwd320;
   Bool_t          b5_HLT_AK8PFJetFwd400;
   Bool_t          b5_HLT_AK8PFJetFwd450;
   Bool_t          b5_HLT_AK8PFJetFwd500;
   Bool_t          b5_HLT_PFHT180;
   Bool_t          b5_HLT_PFHT250;
   Bool_t          b5_HLT_PFHT370;
   Bool_t          b5_HLT_PFHT430;
   Bool_t          b5_HLT_PFHT510;
   Bool_t          b5_HLT_PFHT590;
   Bool_t          b5_HLT_PFHT680;
   Bool_t          b5_HLT_PFHT780;
   Bool_t          b5_HLT_PFHT890;
   Bool_t          b5_HLT_PFHT1050;
   Bool_t          b5_HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t          b5_HLT_PFHT500_PFMET110_PFMHT110_IDTight;
   Bool_t          b5_HLT_PFHT700_PFMET85_PFMHT85_IDTight;
   Bool_t          b5_HLT_PFHT700_PFMET95_PFMHT95_IDTight;
   Bool_t          b5_HLT_PFHT800_PFMET75_PFMHT75_IDTight;
   Bool_t          b5_HLT_PFHT800_PFMET85_PFMHT85_IDTight;
   Bool_t          b5_HLT_PFMET110_PFMHT110_IDTight;
   Bool_t          b5_HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          b5_HLT_PFMET130_PFMHT130_IDTight;
   Bool_t          b5_HLT_PFMET140_PFMHT140_IDTight;
   Bool_t          b5_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          b5_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          b5_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          b5_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          b5_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1;
   Bool_t          b5_HLT_PFMET120_PFMHT120_IDTight_PFHT60;
   Bool_t          b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          b5_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;
   Bool_t          b5_HLT_PFMETTypeOne110_PFMHT110_IDTight;
   Bool_t          b5_HLT_PFMETTypeOne120_PFMHT120_IDTight;
   Bool_t          b5_HLT_PFMETTypeOne130_PFMHT130_IDTight;
   Bool_t          b5_HLT_PFMETTypeOne140_PFMHT140_IDTight;
   Bool_t          b5_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          b5_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          b5_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          b5_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;
   Bool_t          b5_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          b5_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          b5_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          b5_HLT_L1ETMHadSeeds;
   Bool_t          b5_HLT_CaloMHT90;
   Bool_t          b5_HLT_CaloMET80_NotCleaned;
   Bool_t          b5_HLT_CaloMET90_NotCleaned;
   Bool_t          b5_HLT_CaloMET100_NotCleaned;
   Bool_t          b5_HLT_CaloMET110_NotCleaned;
   Bool_t          b5_HLT_CaloMET250_NotCleaned;
   Bool_t          b5_HLT_CaloMET70_HBHECleaned;
   Bool_t          b5_HLT_CaloMET80_HBHECleaned;
   Bool_t          b5_HLT_CaloMET90_HBHECleaned;
   Bool_t          b5_HLT_CaloMET100_HBHECleaned;
   Bool_t          b5_HLT_CaloMET250_HBHECleaned;
   Bool_t          b5_HLT_CaloMET300_HBHECleaned;
   Bool_t          b5_HLT_CaloMET350_HBHECleaned;
   Bool_t          b5_HLT_PFMET200_NotCleaned;
   Bool_t          b5_HLT_PFMET200_HBHECleaned;
   Bool_t          b5_HLT_PFMET250_HBHECleaned;
   Bool_t          b5_HLT_PFMET300_HBHECleaned;
   Bool_t          b5_HLT_PFMET200_HBHE_BeamHaloCleaned;
   Bool_t          b5_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;
   Bool_t          b5_HLT_MET105_IsoTrk50;
   Bool_t          b5_HLT_MET120_IsoTrk50;
   Bool_t          b5_HLT_SingleJet30_Mu12_SinglePFJet40;
   Bool_t          b5_HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_DoublePFJets40_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_DoublePFJets100_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_DoublePFJets200_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_DoublePFJets350_CaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;
   Bool_t          b5_HLT_Photon300_NoHE;
   Bool_t          b5_HLT_Mu8_TrkIsoVVL;
   Bool_t          b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
   Bool_t          b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
   Bool_t          b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
   Bool_t          b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          b5_HLT_Mu17_TrkIsoVVL;
   Bool_t          b5_HLT_Mu19_TrkIsoVVL;
   Bool_t          b5_HLT_BTagMu_AK4DiJet20_Mu5;
   Bool_t          b5_HLT_BTagMu_AK4DiJet40_Mu5;
   Bool_t          b5_HLT_BTagMu_AK4DiJet70_Mu5;
   Bool_t          b5_HLT_BTagMu_AK4DiJet110_Mu5;
   Bool_t          b5_HLT_BTagMu_AK4DiJet170_Mu5;
   Bool_t          b5_HLT_BTagMu_AK4Jet300_Mu5;
   Bool_t          b5_HLT_BTagMu_AK8DiJet170_Mu5;
   Bool_t          b5_HLT_BTagMu_AK8Jet170_DoubleMu5;
   Bool_t          b5_HLT_BTagMu_AK8Jet300_Mu5;
   Bool_t          b5_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL;
   Bool_t          b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          b5_HLT_Mu12_DoublePhoton20;
   Bool_t          b5_HLT_TriplePhoton_20_20_20_CaloIdLV2;
   Bool_t          b5_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL;
   Bool_t          b5_HLT_TriplePhoton_30_30_10_CaloIdLV2;
   Bool_t          b5_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL;
   Bool_t          b5_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL;
   Bool_t          b5_HLT_Photon20;
   Bool_t          b5_HLT_Photon33;
   Bool_t          b5_HLT_Photon50;
   Bool_t          b5_HLT_Photon75;
   Bool_t          b5_HLT_Photon90;
   Bool_t          b5_HLT_Photon120;
   Bool_t          b5_HLT_Photon150;
   Bool_t          b5_HLT_Photon175;
   Bool_t          b5_HLT_Photon200;
   Bool_t          b5_HLT_Photon100EB_TightID_TightIso;
   Bool_t          b5_HLT_Photon110EB_TightID_TightIso;
   Bool_t          b5_HLT_Photon120EB_TightID_TightIso;
   Bool_t          b5_HLT_Photon100EBHE10;
   Bool_t          b5_HLT_Photon100EEHE10;
   Bool_t          b5_HLT_Photon100EE_TightID_TightIso;
   Bool_t          b5_HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t          b5_HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t          b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3;
   Bool_t          b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3;
   Bool_t          b5_HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t          b5_HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t          b5_HLT_Photon165_R9Id90_HE10_IsoM;
   Bool_t          b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
   Bool_t          b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
   Bool_t          b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;
   Bool_t          b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;
   Bool_t          b5_HLT_Dimuon0_Jpsi_L1_NoOS;
   Bool_t          b5_HLT_Dimuon0_Jpsi_NoVertexing_NoOS;
   Bool_t          b5_HLT_Dimuon0_Jpsi;
   Bool_t          b5_HLT_Dimuon0_Jpsi_NoVertexing;
   Bool_t          b5_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;
   Bool_t          b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;
   Bool_t          b5_HLT_Dimuon0_Jpsi3p5_Muon2;
   Bool_t          b5_HLT_Dimuon0_Upsilon_L1_4p5;
   Bool_t          b5_HLT_Dimuon0_Upsilon_L1_5;
   Bool_t          b5_HLT_Dimuon0_Upsilon_L1_4p5NoOS;
   Bool_t          b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0;
   Bool_t          b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0M;
   Bool_t          b5_HLT_Dimuon0_Upsilon_NoVertexing;
   Bool_t          b5_HLT_Dimuon0_Upsilon_L1_5M;
   Bool_t          b5_HLT_Dimuon0_LowMass_L1_0er1p5R;
   Bool_t          b5_HLT_Dimuon0_LowMass_L1_0er1p5;
   Bool_t          b5_HLT_Dimuon0_LowMass;
   Bool_t          b5_HLT_Dimuon0_LowMass_L1_4;
   Bool_t          b5_HLT_Dimuon0_LowMass_L1_4R;
   Bool_t          b5_HLT_Dimuon0_LowMass_L1_TM530;
   Bool_t          b5_HLT_Dimuon0_Upsilon_Muon_L1_TM0;
   Bool_t          b5_HLT_Dimuon0_Upsilon_Muon_NoL1Mass;
   Bool_t          b5_HLT_TripleMu_5_3_3_Mass3p8_DZ;
   Bool_t          b5_HLT_TripleMu_10_5_5_DZ;
   Bool_t          b5_HLT_TripleMu_12_10_5;
   Bool_t          b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;
   Bool_t          b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;
   Bool_t          b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;
   Bool_t          b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;
   Bool_t          b5_HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
   Bool_t          b5_HLT_DoubleMu3_DZ_PFMET70_PFMHT70;
   Bool_t          b5_HLT_DoubleMu3_DZ_PFMET90_PFMHT90;
   Bool_t          b5_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;
   Bool_t          b5_HLT_DoubleMu4_Jpsi_Displaced;
   Bool_t          b5_HLT_DoubleMu4_Jpsi_NoVertexing;
   Bool_t          b5_HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   Bool_t          b5_HLT_DoubleMu43NoFiltersNoVtx;
   Bool_t          b5_HLT_DoubleMu48NoFiltersNoVtx;
   Bool_t          b5_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;
   Bool_t          b5_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;
   Bool_t          b5_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;
   Bool_t          b5_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;
   Bool_t          b5_HLT_DoubleMu33NoFiltersNoVtxDisplaced;
   Bool_t          b5_HLT_DoubleMu40NoFiltersNoVtxDisplaced;
   Bool_t          b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4;
   Bool_t          b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG;
   Bool_t          b5_HLT_HT425;
   Bool_t          b5_HLT_HT430_DisplacedDijet40_DisplacedTrack;
   Bool_t          b5_HLT_HT500_DisplacedDijet40_DisplacedTrack;
   Bool_t          b5_HLT_HT430_DisplacedDijet60_DisplacedTrack;
   Bool_t          b5_HLT_HT400_DisplacedDijet40_DisplacedTrack;
   Bool_t          b5_HLT_HT650_DisplacedDijet60_Inclusive;
   Bool_t          b5_HLT_HT550_DisplacedDijet60_Inclusive;
   Bool_t          b5_HLT_DiJet110_35_Mjj650_PFMET110;
   Bool_t          b5_HLT_DiJet110_35_Mjj650_PFMET120;
   Bool_t          b5_HLT_DiJet110_35_Mjj650_PFMET130;
   Bool_t          b5_HLT_TripleJet110_35_35_Mjj650_PFMET110;
   Bool_t          b5_HLT_TripleJet110_35_35_Mjj650_PFMET120;
   Bool_t          b5_HLT_TripleJet110_35_35_Mjj650_PFMET130;
   Bool_t          b5_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1;
   Bool_t          b5_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1;
   Bool_t          b5_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1;
   Bool_t          b5_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
   Bool_t          b5_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
   Bool_t          b5_HLT_Ele28_HighEta_SC20_Mass55;
   Bool_t          b5_HLT_DoubleMu20_7_Mass0to30_Photon23;
   Bool_t          b5_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
   Bool_t          b5_HLT_Ele15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          b5_HLT_Ele15_IsoVVVL_PFHT450;
   Bool_t          b5_HLT_Ele50_IsoVVVL_PFHT450;
   Bool_t          b5_HLT_Ele15_IsoVVVL_PFHT600;
   Bool_t          b5_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t          b5_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;
   Bool_t          b5_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;
   Bool_t          b5_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;
   Bool_t          b5_HLT_Mu15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          b5_HLT_Mu15_IsoVVVL_PFHT450;
   Bool_t          b5_HLT_Mu50_IsoVVVL_PFHT450;
   Bool_t          b5_HLT_Mu15_IsoVVVL_PFHT600;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;
   Bool_t          b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;
   Bool_t          b5_HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
   Bool_t          b5_HLT_Dimuon20_Jpsi_Barrel_Seagulls;
   Bool_t          b5_HLT_Dimuon12_Upsilon_y1p4;
   Bool_t          b5_HLT_Dimuon14_Phi_Barrel_Seagulls;
   Bool_t          b5_HLT_Dimuon18_PsiPrime;
   Bool_t          b5_HLT_Dimuon25_Jpsi;
   Bool_t          b5_HLT_Dimuon18_PsiPrime_noCorrL1;
   Bool_t          b5_HLT_Dimuon24_Upsilon_noCorrL1;
   Bool_t          b5_HLT_Dimuon24_Phi_noCorrL1;
   Bool_t          b5_HLT_Dimuon25_Jpsi_noCorrL1;
   Bool_t          b5_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
   Bool_t          b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Bool_t          b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
   Bool_t          b5_HLT_DoubleIsoMu20_eta2p1;
   Bool_t          b5_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;
   Bool_t          b5_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx;
   Bool_t          b5_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;
   Bool_t          b5_HLT_Mu8;
   Bool_t          b5_HLT_Mu17;
   Bool_t          b5_HLT_Mu19;
   Bool_t          b5_HLT_Mu17_Photon30_IsoCaloId;
   Bool_t          b5_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          b5_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          b5_HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          b5_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          b5_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Bool_t          b5_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t          b5_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Bool_t          b5_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
   Bool_t          b5_HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Bool_t          b5_HLT_Ele135_CaloIdVT_GsfTrkIdT;
   Bool_t          b5_HLT_Ele145_CaloIdVT_GsfTrkIdT;
   Bool_t          b5_HLT_Ele200_CaloIdVT_GsfTrkIdT;
   Bool_t          b5_HLT_Ele250_CaloIdVT_GsfTrkIdT;
   Bool_t          b5_HLT_Ele300_CaloIdVT_GsfTrkIdT;
   Bool_t          b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5;
   Bool_t          b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40;
   Bool_t          b5_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2;
   Bool_t          b5_HLT_PFHT380_SixPFJet32;
   Bool_t          b5_HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5;
   Bool_t          b5_HLT_PFHT430_SixPFJet40;
   Bool_t          b5_HLT_PFHT350;
   Bool_t          b5_HLT_PFHT350MinPFJet15;
   Bool_t          b5_HLT_Photon60_R9Id90_CaloIdL_IsoL;
   Bool_t          b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL;
   Bool_t          b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;
   Bool_t          b5_HLT_ECALHT800;
   Bool_t          b5_HLT_DiSC30_18_EIso_AND_HE_Mass70;
   Bool_t          b5_HLT_Physics;
   Bool_t          b5_HLT_Physics_part0;
   Bool_t          b5_HLT_Physics_part1;
   Bool_t          b5_HLT_Physics_part2;
   Bool_t          b5_HLT_Physics_part3;
   Bool_t          b5_HLT_Physics_part4;
   Bool_t          b5_HLT_Physics_part5;
   Bool_t          b5_HLT_Physics_part6;
   Bool_t          b5_HLT_Physics_part7;
   Bool_t          b5_HLT_Random;
   Bool_t          b5_HLT_ZeroBias;
   Bool_t          b5_HLT_ZeroBias_part0;
   Bool_t          b5_HLT_ZeroBias_part1;
   Bool_t          b5_HLT_ZeroBias_part2;
   Bool_t          b5_HLT_ZeroBias_part3;
   Bool_t          b5_HLT_ZeroBias_part4;
   Bool_t          b5_HLT_ZeroBias_part5;
   Bool_t          b5_HLT_ZeroBias_part6;
   Bool_t          b5_HLT_ZeroBias_part7;
   Bool_t          b5_HLT_AK4CaloJet30;
   Bool_t          b5_HLT_AK4CaloJet40;
   Bool_t          b5_HLT_AK4CaloJet50;
   Bool_t          b5_HLT_AK4CaloJet80;
   Bool_t          b5_HLT_AK4CaloJet100;
   Bool_t          b5_HLT_AK4CaloJet120;
   Bool_t          b5_HLT_AK4PFJet30;
   Bool_t          b5_HLT_AK4PFJet50;
   Bool_t          b5_HLT_AK4PFJet80;
   Bool_t          b5_HLT_AK4PFJet100;
   Bool_t          b5_HLT_AK4PFJet120;
   Bool_t          b5_HLT_SinglePhoton10_Eta3p1ForPPRef;
   Bool_t          b5_HLT_SinglePhoton20_Eta3p1ForPPRef;
   Bool_t          b5_HLT_SinglePhoton30_Eta3p1ForPPRef;
   Bool_t          b5_HLT_Photon20_HoverELoose;
   Bool_t          b5_HLT_Photon30_HoverELoose;
   Bool_t          b5_HLT_EcalCalibration;
   Bool_t          b5_HLT_HcalCalibration;
   Bool_t          b5_HLT_L1UnpairedBunchBptxMinus;
   Bool_t          b5_HLT_L1UnpairedBunchBptxPlus;
   Bool_t          b5_HLT_L1NotBptxOR;
   Bool_t          b5_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          b5_HLT_HcalNZS;
   Bool_t          b5_HLT_HcalPhiSym;
   Bool_t          b5_HLT_HcalIsolatedbunch;
   Bool_t          b5_HLT_IsoTrackHB;
   Bool_t          b5_HLT_IsoTrackHE;
   Bool_t          b5_HLT_ZeroBias_FirstCollisionAfterAbortGap;
   Bool_t          b5_HLT_ZeroBias_IsolatedBunches;
   Bool_t          b5_HLT_ZeroBias_FirstCollisionInTrain;
   Bool_t          b5_HLT_ZeroBias_LastCollisionInTrain;
   Bool_t          b5_HLT_ZeroBias_FirstBXAfterTrain;
   Bool_t          b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1;
   Bool_t          b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1;
   Bool_t          b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg;
   Bool_t          b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg;
   Bool_t          b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;
   Bool_t          b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;
   Bool_t          b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;
   Bool_t          b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;
   Bool_t          b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110;
   Bool_t          b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120;
   Bool_t          b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130;
   Bool_t          b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140;
   Bool_t          b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;
   Bool_t          b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;
   Bool_t          b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          b5_HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          b5_HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1;
   Bool_t          b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;
   Bool_t          b5_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
   Bool_t          b5_HLT_Rsq0p35;
   Bool_t          b5_HLT_Rsq0p40;
   Bool_t          b5_HLT_RsqMR300_Rsq0p09_MR200;
   Bool_t          b5_HLT_RsqMR320_Rsq0p09_MR200;
   Bool_t          b5_HLT_RsqMR300_Rsq0p09_MR200_4jet;
   Bool_t          b5_HLT_RsqMR320_Rsq0p09_MR200_4jet;
   Bool_t          b5_HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1;
   Bool_t          b5_HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1;
   Bool_t          b5_HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1;
   Bool_t          b5_HLT_IsoMu27_MET90;
   Bool_t          b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1;
   Bool_t          b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1;
   Bool_t          b5_HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg;
   Bool_t          b5_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50;
   Bool_t          b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t          b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3;
   Bool_t          b5_HLT_PFMET100_PFMHT100_IDTight_PFHT60;
   Bool_t          b5_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
   Bool_t          b5_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;
   Bool_t          b5_HLT_Mu18_Mu9_SameSign;
   Bool_t          b5_HLT_Mu18_Mu9_SameSign_DZ;
   Bool_t          b5_HLT_Mu18_Mu9;
   Bool_t          b5_HLT_Mu18_Mu9_DZ;
   Bool_t          b5_HLT_Mu20_Mu10_SameSign;
   Bool_t          b5_HLT_Mu20_Mu10_SameSign_DZ;
   Bool_t          b5_HLT_Mu20_Mu10;
   Bool_t          b5_HLT_Mu20_Mu10_DZ;
   Bool_t          b5_HLT_Mu23_Mu12_SameSign;
   Bool_t          b5_HLT_Mu23_Mu12_SameSign_DZ;
   Bool_t          b5_HLT_Mu23_Mu12;
   Bool_t          b5_HLT_Mu23_Mu12_DZ;
   Bool_t          b5_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
   Bool_t          b5_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;
   Bool_t          b5_HLT_DoubleMu3_DCA_PFMET50_PFMHT60;
   Bool_t          b5_HLT_TripleMu_5_3_3_Mass3p8_DCA;
   Bool_t          b5_HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          b5_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          b5_HLT_QuadPFJet105_90_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          b5_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;
   Bool_t          b5_HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          b5_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          b5_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          b5_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2;
   Bool_t          b5_HLT_QuadPFJet98_83_71_15;
   Bool_t          b5_HLT_QuadPFJet103_88_75_15;
   Bool_t          b5_HLT_QuadPFJet105_88_76_15;
   Bool_t          b5_HLT_QuadPFJet111_90_80_15;
   Bool_t          b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17;
   Bool_t          b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1;
   Bool_t          b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;
   Bool_t          b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;
   Bool_t          b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;
   Bool_t          b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55;
   Bool_t          b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto;
   Bool_t          b5_HLT_Mu8p5_IP3p5_part0;
   Bool_t          b5_HLT_Mu8p5_IP3p5_part1;
   Bool_t          b5_HLT_Mu8p5_IP3p5_part2;
   Bool_t          b5_HLT_Mu8p5_IP3p5_part3;
   Bool_t          b5_HLT_Mu8p5_IP3p5_part4;
   Bool_t          b5_HLT_Mu8p5_IP3p5_part5;
   Bool_t          b5_HLT_Mu10p5_IP3p5_part0;
   Bool_t          b5_HLT_Mu10p5_IP3p5_part1;
   Bool_t          b5_HLT_Mu10p5_IP3p5_part2;
   Bool_t          b5_HLT_Mu10p5_IP3p5_part3;
   Bool_t          b5_HLT_Mu10p5_IP3p5_part4;
   Bool_t          b5_HLT_Mu10p5_IP3p5_part5;
   Bool_t          b5_HLT_Mu9_IP6_part0;
   Bool_t          b5_HLT_Mu9_IP6_part1;
   Bool_t          b5_HLT_Mu9_IP6_part2;
   Bool_t          b5_HLT_Mu9_IP6_part3;
   Bool_t          b5_HLT_Mu9_IP6_part4;
   Bool_t          b5_HLT_Mu9_IP6_part5;
   Bool_t          b5_HLT_Mu8_IP3_part0;
   Bool_t          b5_HLT_Mu8_IP3_part1;
   Bool_t          b5_HLT_Mu8_IP3_part2;
   Bool_t          b5_HLT_Mu8_IP3_part3;
   Bool_t          b5_HLT_Mu8_IP3_part4;
   Bool_t          b5_HLT_Mu8_IP3_part5;
   Bool_t          b5_HLTriggerFinalPath;
   UInt_t          bG_run;
   ULong64_t       bG_event;
   UInt_t          bG_lumis;
   Bool_t          bG_isData;
   Int_t           bG_nSC;
   UInt_t          bG_nscE;
   Float_t         bG_scE[37];   //[bG_nscE]
   UInt_t          bG_nscEt;
   Float_t         bG_scEt[37];   //[bG_nscEt]
   UInt_t          bG_nscRawE;
   Float_t         bG_scRawE[37];   //[bG_nscRawE]
   UInt_t          bG_nscEta;
   Float_t         bG_scEta[37];   //[bG_nscEta]
   UInt_t          bG_nscPhi;
   Float_t         bG_scPhi[37];   //[bG_nscPhi]
   UInt_t          bG_nscX;
   Float_t         bG_scX[37];   //[bG_nscX]
   UInt_t          bG_nscY;
   Float_t         bG_scY[37];   //[bG_nscY]
   UInt_t          bG_nscZ;
   Float_t         bG_scZ[37];   //[bG_nscZ]
   UInt_t          bG_nscEtaWidth;
   Float_t         bG_scEtaWidth[37];   //[bG_nscEtaWidth]
   UInt_t          bG_nscPhiWidth;
   Float_t         bG_scPhiWidth[37];   //[bG_nscPhiWidth]
   UInt_t          bG_nscRawEt;
   Float_t         bG_scRawEt[37];   //[bG_nscRawEt]
   UInt_t          bG_nscMinDrWithGsfElectornSC_;
   Float_t         bG_scMinDrWithGsfElectornSC_[kMaxbG_nscMinDrWithGsfElectornSC];   //[bG_nscMinDrWithGsfElectornSC_]
   UInt_t          bG_nscFoundGsfMatch_;
   Bool_t          bG_scFoundGsfMatch_[kMaxbG_nscFoundGsfMatch];   //[bG_nscFoundGsfMatch_]
   UInt_t          bG_nscE5x5;
   Float_t         bG_scE5x5[37];   //[bG_nscE5x5]
   UInt_t          bG_nscE2x2Ratio;
   Float_t         bG_scE2x2Ratio[37];   //[bG_nscE2x2Ratio]
   UInt_t          bG_nscE3x3Ratio;
   Float_t         bG_scE3x3Ratio[37];   //[bG_nscE3x3Ratio]
   UInt_t          bG_nscEMaxRatio;
   Float_t         bG_scEMaxRatio[37];   //[bG_nscEMaxRatio]
   UInt_t          bG_nscE2ndRatio;
   Float_t         bG_scE2ndRatio[37];   //[bG_nscE2ndRatio]
   UInt_t          bG_nscETopRatio;
   Float_t         bG_scETopRatio[37];   //[bG_nscETopRatio]
   UInt_t          bG_nscERightRatio;
   Float_t         bG_scERightRatio[37];   //[bG_nscERightRatio]
   UInt_t          bG_nscEBottomRatio;
   Float_t         bG_scEBottomRatio[37];   //[bG_nscEBottomRatio]
   UInt_t          bG_nscELeftRatio;
   Float_t         bG_scELeftRatio[37];   //[bG_nscELeftRatio]
   UInt_t          bG_nscE2x5MaxRatio;
   Float_t         bG_scE2x5MaxRatio[37];   //[bG_nscE2x5MaxRatio]
   UInt_t          bG_nscE2x5TopRatio;
   Float_t         bG_scE2x5TopRatio[37];   //[bG_nscE2x5TopRatio]
   UInt_t          bG_nscE2x5RightRatio;
   Float_t         bG_scE2x5RightRatio[37];   //[bG_nscE2x5RightRatio]
   UInt_t          bG_nscE2x5BottomRatio;
   Float_t         bG_scE2x5BottomRatio[37];   //[bG_nscE2x5BottomRatio]
   UInt_t          bG_nscE2x5LeftRatio;
   Float_t         bG_scE2x5LeftRatio[37];   //[bG_nscE2x5LeftRatio]
   UInt_t          bG_nscSwissCross;
   Float_t         bG_scSwissCross[37];   //[bG_nscSwissCross]
   UInt_t          bG_nscR9;
   Float_t         bG_scR9[37];   //[bG_nscR9]
   UInt_t          bG_nscSigmaIetaIeta;
   Float_t         bG_scSigmaIetaIeta[37];   //[bG_nscSigmaIetaIeta]
   UInt_t          bG_nscSigmaIetaIphi;
   Float_t         bG_scSigmaIetaIphi[37];   //[bG_nscSigmaIetaIphi]
   UInt_t          bG_nscSigmaIphiIphi;
   Float_t         bG_scSigmaIphiIphi[37];   //[bG_nscSigmaIphiIphi]
   UInt_t          bG_nscFull5x5_e5x5;
   Float_t         bG_scFull5x5_e5x5[37];   //[bG_nscFull5x5_e5x5]
   UInt_t          bG_nscFull5x5_e2x2Ratio;
   Float_t         bG_scFull5x5_e2x2Ratio[37];   //[bG_nscFull5x5_e2x2Ratio]
   UInt_t          bG_nscFull5x5_e3x3Ratio;
   Float_t         bG_scFull5x5_e3x3Ratio[37];   //[bG_nscFull5x5_e3x3Ratio]
   UInt_t          bG_nscFull5x5_eMaxRatio;
   Float_t         bG_scFull5x5_eMaxRatio[37];   //[bG_nscFull5x5_eMaxRatio]
   UInt_t          bG_nscFull5x5_e2ndRatio;
   Float_t         bG_scFull5x5_e2ndRatio[37];   //[bG_nscFull5x5_e2ndRatio]
   UInt_t          bG_nscFull5x5_eTopRatio;
   Float_t         bG_scFull5x5_eTopRatio[37];   //[bG_nscFull5x5_eTopRatio]
   UInt_t          bG_nscFull5x5_eRightRatio;
   Float_t         bG_scFull5x5_eRightRatio[37];   //[bG_nscFull5x5_eRightRatio]
   UInt_t          bG_nscFull5x5_eBottomRatio;
   Float_t         bG_scFull5x5_eBottomRatio[37];   //[bG_nscFull5x5_eBottomRatio]
   UInt_t          bG_nscFull5x5_eLeftRatio;
   Float_t         bG_scFull5x5_eLeftRatio[37];   //[bG_nscFull5x5_eLeftRatio]
   UInt_t          bG_nscFull5x5_e2x5MaxRatio;
   Float_t         bG_scFull5x5_e2x5MaxRatio[37];   //[bG_nscFull5x5_e2x5MaxRatio]
   UInt_t          bG_nscFull5x5_e2x5TopRatio;
   Float_t         bG_scFull5x5_e2x5TopRatio[37];   //[bG_nscFull5x5_e2x5TopRatio]
   UInt_t          bG_nscFull5x5_e2x5RightRatio;
   Float_t         bG_scFull5x5_e2x5RightRatio[37];   //[bG_nscFull5x5_e2x5RightRatio]
   UInt_t          bG_nscFull5x5_e2x5BottomRatio;
   Float_t         bG_scFull5x5_e2x5BottomRatio[37];   //[bG_nscFull5x5_e2x5BottomRatio]
   UInt_t          bG_nscFull5x5_e2x5LeftRatio;
   Float_t         bG_scFull5x5_e2x5LeftRatio[37];   //[bG_nscFull5x5_e2x5LeftRatio]
   UInt_t          bG_nscFull5x5_swissCross;
   Float_t         bG_scFull5x5_swissCross[37];   //[bG_nscFull5x5_swissCross]
   UInt_t          bG_nscFull5x5_r9;
   Float_t         bG_scFull5x5_r9[37];   //[bG_nscFull5x5_r9]
   UInt_t          bG_nscFull5x5_sigmaIetaIeta;
   Float_t         bG_scFull5x5_sigmaIetaIeta[37];   //[bG_nscFull5x5_sigmaIetaIeta]
   UInt_t          bG_nscFull5x5_sigmaIetaIphi;
   Float_t         bG_scFull5x5_sigmaIetaIphi[37];   //[bG_nscFull5x5_sigmaIetaIphi]
   UInt_t          bG_nscFull5x5_sigmaIphiIphi;
   Float_t         bG_scFull5x5_sigmaIphiIphi[37];   //[bG_nscFull5x5_sigmaIphiIphi]
   UInt_t          bG_nscNHcalRecHitInDIEta5IPhi5;
   Float_t         bG_scNHcalRecHitInDIEta5IPhi5[37];   //[bG_nscNHcalRecHitInDIEta5IPhi5]
   UInt_t          bG_nscEFromHcalRecHitInDIEta5IPhi5;
   Float_t         bG_scEFromHcalRecHitInDIEta5IPhi5[37];   //[bG_nscEFromHcalRecHitInDIEta5IPhi5]
   UInt_t          bG_nscNHcalRecHitInDIEta2IPhi2;
   Float_t         bG_scNHcalRecHitInDIEta2IPhi2[37];   //[bG_nscNHcalRecHitInDIEta2IPhi2]
   UInt_t          bG_nscEFromHcalRecHitInDIEta2IPhi2;
   Float_t         bG_scEFromHcalRecHitInDIEta2IPhi2[37];   //[bG_nscEFromHcalRecHitInDIEta2IPhi2]
   UInt_t          bG_nscPFChIso1;
   Float_t         bG_scPFChIso1[37];   //[bG_nscPFChIso1]
   UInt_t          bG_nscPFChIso2;
   Float_t         bG_scPFChIso2[37];   //[bG_nscPFChIso2]
   UInt_t          bG_nscPFChIso3;
   Float_t         bG_scPFChIso3[37];   //[bG_nscPFChIso3]
   UInt_t          bG_nscPFChIso4;
   Float_t         bG_scPFChIso4[37];   //[bG_nscPFChIso4]
   UInt_t          bG_nscPFChIso5;
   Float_t         bG_scPFChIso5[37];   //[bG_nscPFChIso5]
   UInt_t          bG_nscPFPhoIso1;
   Float_t         bG_scPFPhoIso1[37];   //[bG_nscPFPhoIso1]
   UInt_t          bG_nscPFPhoIso2;
   Float_t         bG_scPFPhoIso2[37];   //[bG_nscPFPhoIso2]
   UInt_t          bG_nscPFPhoIso3;
   Float_t         bG_scPFPhoIso3[37];   //[bG_nscPFPhoIso3]
   UInt_t          bG_nscPFPhoIso4;
   Float_t         bG_scPFPhoIso4[37];   //[bG_nscPFPhoIso4]
   UInt_t          bG_nscPFPhoIso5;
   Float_t         bG_scPFPhoIso5[37];   //[bG_nscPFPhoIso5]
   UInt_t          bG_nscPFNeuIso1;
   Float_t         bG_scPFNeuIso1[37];   //[bG_nscPFNeuIso1]
   UInt_t          bG_nscPFNeuIso2;
   Float_t         bG_scPFNeuIso2[37];   //[bG_nscPFNeuIso2]
   UInt_t          bG_nscPFNeuIso3;
   Float_t         bG_scPFNeuIso3[37];   //[bG_nscPFNeuIso3]
   UInt_t          bG_nscPFNeuIso4;
   Float_t         bG_scPFNeuIso4[37];   //[bG_nscPFNeuIso4]
   UInt_t          bG_nscPFNeuIso5;
   Float_t         bG_scPFNeuIso5[37];   //[bG_nscPFNeuIso5]
   Int_t           bG_nPrimaryVertex;
   Float_t         bG_primaryVertex_isFake[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_x[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_y[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_z[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_t[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_covXX[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_covXY[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_covXZ[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_covYY[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_covYZ[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_covZZ[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_x_error[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_y_error[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_z_error[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_t_error[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_ntracks[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_ndof[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_chi2[84];   //[bG_nPrimaryVertex]
   Float_t         bG_primaryVertex_normalizedChi2[84];   //[bG_nPrimaryVertex]

   // List of branches
   TBranch        *b_b5_run;   //!
   TBranch        *b_b5_luminosityBlock;   //!
   TBranch        *b_b5_event;   //!
   TBranch        *b_b5_nMuonId;   //!
   TBranch        *b_b5_MuonId_chi2LocalPosition;   //!
   TBranch        *b_b5_MuonId_glbNormChi2;   //!
   TBranch        *b_b5_MuonId_glbTrackProbability;   //!
   TBranch        *b_b5_MuonId_match1_dX;   //!
   TBranch        *b_b5_MuonId_match1_dY;   //!
   TBranch        *b_b5_MuonId_match1_pullDxDz;   //!
   TBranch        *b_b5_MuonId_match1_pullDyDz;   //!
   TBranch        *b_b5_MuonId_match1_pullX;   //!
   TBranch        *b_b5_MuonId_match1_pullY;   //!
   TBranch        *b_b5_MuonId_match2_dX;   //!
   TBranch        *b_b5_MuonId_match2_dY;   //!
   TBranch        *b_b5_MuonId_match2_pullDxDz;   //!
   TBranch        *b_b5_MuonId_match2_pullDyDz;   //!
   TBranch        *b_b5_MuonId_match2_pullX;   //!
   TBranch        *b_b5_MuonId_match2_pullY;   //!
   TBranch        *b_b5_MuonId_newSoftMuonMva;   //!
   TBranch        *b_b5_MuonId_trkKink;   //!
   TBranch        *b_b5_MuonId_trkValidFrac;   //!
   TBranch        *b_b5_MuonId_highPurity;   //!
   TBranch        *b_b5_MuonId_nLostHitsInner;   //!
   TBranch        *b_b5_MuonId_nLostHitsOn;   //!
   TBranch        *b_b5_MuonId_nLostHitsOuter;   //!
   TBranch        *b_b5_MuonId_nPixels;   //!
   TBranch        *b_b5_MuonId_nValidHits;   //!
   TBranch        *b_b5_MuonId_trkLayers;   //!
   TBranch        *b_b5_MuonId_trkLostLayersInner;   //!
   TBranch        *b_b5_MuonId_trkLostLayersOn;   //!
   TBranch        *b_b5_MuonId_trkLostLayersOuter;   //!
   TBranch        *b_b5_nmm;   //!
   TBranch        *b_b5_mm_bdt;   //!
   TBranch        *b_b5_mm_doca;   //!
   TBranch        *b_b5_mm_docatrk;   //!
   TBranch        *b_b5_mm_iso;   //!
   TBranch        *b_b5_mm_kal_lxy;   //!
   TBranch        *b_b5_mm_kal_mass;   //!
   TBranch        *b_b5_mm_kal_slxy;   //!
   TBranch        *b_b5_mm_kal_vtx_prob;   //!
   TBranch        *b_b5_mm_kin_alpha;   //!
   TBranch        *b_b5_mm_kin_alphaBS;   //!
   TBranch        *b_b5_mm_kin_alphaBSErr;   //!
   TBranch        *b_b5_mm_kin_alphaErr;   //!
   TBranch        *b_b5_mm_kin_eta;   //!
   TBranch        *b_b5_mm_kin_l3d;   //!
   TBranch        *b_b5_mm_kin_lxy;   //!
   TBranch        *b_b5_mm_kin_mass;   //!
   TBranch        *b_b5_mm_kin_massErr;   //!
   TBranch        *b_b5_mm_kin_mu1eta;   //!
   TBranch        *b_b5_mm_kin_mu1phi;   //!
   TBranch        *b_b5_mm_kin_mu1pt;   //!
   TBranch        *b_b5_mm_kin_mu2eta;   //!
   TBranch        *b_b5_mm_kin_mu2phi;   //!
   TBranch        *b_b5_mm_kin_mu2pt;   //!
   TBranch        *b_b5_mm_kin_phi;   //!
   TBranch        *b_b5_mm_kin_pt;   //!
   TBranch        *b_b5_mm_kin_pv2ip;   //!
   TBranch        *b_b5_mm_kin_pv2ipErr;   //!
   TBranch        *b_b5_mm_kin_pv2lip;   //!
   TBranch        *b_b5_mm_kin_pv2lipErr;   //!
   TBranch        *b_b5_mm_kin_pv2lipSig;   //!
   TBranch        *b_b5_mm_kin_pv_z;   //!
   TBranch        *b_b5_mm_kin_pv_zErr;   //!
   TBranch        *b_b5_mm_kin_pvip;   //!
   TBranch        *b_b5_mm_kin_pvipErr;   //!
   TBranch        *b_b5_mm_kin_pvlip;   //!
   TBranch        *b_b5_mm_kin_pvlipErr;   //!
   TBranch        *b_b5_mm_kin_pvlipSig;   //!
   TBranch        *b_b5_mm_kin_sl3d;   //!
   TBranch        *b_b5_mm_kin_slxy;   //!
   TBranch        *b_b5_mm_kin_spv2ip;   //!
   TBranch        *b_b5_mm_kin_spvip;   //!
   TBranch        *b_b5_mm_kin_tau;   //!
   TBranch        *b_b5_mm_kin_taue;   //!
   TBranch        *b_b5_mm_kin_tauxy;   //!
   TBranch        *b_b5_mm_kin_tauxye;   //!
   TBranch        *b_b5_mm_kin_vtx_chi2dof;   //!
   TBranch        *b_b5_mm_kin_vtx_prob;   //!
   TBranch        *b_b5_mm_kin_vtx_x;   //!
   TBranch        *b_b5_mm_kin_vtx_xErr;   //!
   TBranch        *b_b5_mm_kin_vtx_y;   //!
   TBranch        *b_b5_mm_kin_vtx_yErr;   //!
   TBranch        *b_b5_mm_kin_vtx_z;   //!
   TBranch        *b_b5_mm_kin_vtx_zErr;   //!
   TBranch        *b_b5_mm_kinpc_alpha;   //!
   TBranch        *b_b5_mm_kinpc_alphaBS;   //!
   TBranch        *b_b5_mm_kinpc_alphaBSErr;   //!
   TBranch        *b_b5_mm_kinpc_alphaErr;   //!
   TBranch        *b_b5_mm_kinpc_eta;   //!
   TBranch        *b_b5_mm_kinpc_l3d;   //!
   TBranch        *b_b5_mm_kinpc_lxy;   //!
   TBranch        *b_b5_mm_kinpc_mass;   //!
   TBranch        *b_b5_mm_kinpc_massErr;   //!
   TBranch        *b_b5_mm_kinpc_phi;   //!
   TBranch        *b_b5_mm_kinpc_pt;   //!
   TBranch        *b_b5_mm_kinpc_pv2ip;   //!
   TBranch        *b_b5_mm_kinpc_pv2ipErr;   //!
   TBranch        *b_b5_mm_kinpc_pv2lip;   //!
   TBranch        *b_b5_mm_kinpc_pv2lipErr;   //!
   TBranch        *b_b5_mm_kinpc_pv2lipSig;   //!
   TBranch        *b_b5_mm_kinpc_pv_z;   //!
   TBranch        *b_b5_mm_kinpc_pv_zErr;   //!
   TBranch        *b_b5_mm_kinpc_pvip;   //!
   TBranch        *b_b5_mm_kinpc_pvipErr;   //!
   TBranch        *b_b5_mm_kinpc_pvlip;   //!
   TBranch        *b_b5_mm_kinpc_pvlipErr;   //!
   TBranch        *b_b5_mm_kinpc_pvlipSig;   //!
   TBranch        *b_b5_mm_kinpc_sl3d;   //!
   TBranch        *b_b5_mm_kinpc_slxy;   //!
   TBranch        *b_b5_mm_kinpc_spv2ip;   //!
   TBranch        *b_b5_mm_kinpc_spvip;   //!
   TBranch        *b_b5_mm_kinpc_tau;   //!
   TBranch        *b_b5_mm_kinpc_taue;   //!
   TBranch        *b_b5_mm_kinpc_tauxy;   //!
   TBranch        *b_b5_mm_kinpc_tauxye;   //!
   TBranch        *b_b5_mm_kinpc_vtx_chi2dof;   //!
   TBranch        *b_b5_mm_kinpc_vtx_prob;   //!
   TBranch        *b_b5_mm_kinpc_vtx_x;   //!
   TBranch        *b_b5_mm_kinpc_vtx_xErr;   //!
   TBranch        *b_b5_mm_kinpc_vtx_y;   //!
   TBranch        *b_b5_mm_kinpc_vtx_yErr;   //!
   TBranch        *b_b5_mm_kinpc_vtx_z;   //!
   TBranch        *b_b5_mm_kinpc_vtx_zErr;   //!
   TBranch        *b_b5_mm_m1iso;   //!
   TBranch        *b_b5_mm_m2iso;   //!
   TBranch        *b_b5_mm_mass;   //!
   TBranch        *b_b5_mm_mu1_eta;   //!
   TBranch        *b_b5_mm_mu1_phi;   //!
   TBranch        *b_b5_mm_mu1_pt;   //!
   TBranch        *b_b5_mm_mu2_eta;   //!
   TBranch        *b_b5_mm_mu2_phi;   //!
   TBranch        *b_b5_mm_mu2_pt;   //!
   TBranch        *b_b5_mm_mva;   //!
   TBranch        *b_b5_mm_otherVtxMaxProb;   //!
   TBranch        *b_b5_mm_otherVtxMaxProb1;   //!
   TBranch        *b_b5_mm_otherVtxMaxProb2;   //!
   TBranch        *b_b5_mm_closetrk;   //!
   TBranch        *b_b5_mm_closetrks1;   //!
   TBranch        *b_b5_mm_closetrks2;   //!
   TBranch        *b_b5_mm_closetrks3;   //!
   TBranch        *b_b5_mm_kal_valid;   //!
   TBranch        *b_b5_mm_kin_valid;   //!
   TBranch        *b_b5_mm_kinpc_valid;   //!
   TBranch        *b_b5_mm_mu1_index;   //!
   TBranch        *b_b5_mm_mu1_pdgId;   //!
   TBranch        *b_b5_mm_mu2_index;   //!
   TBranch        *b_b5_mm_mu2_pdgId;   //!
   TBranch        *b_b5_mm_nBMTrks;   //!
   TBranch        *b_b5_mm_nDisTrks;   //!
   TBranch        *b_b5_mm_nTrks;   //!
   TBranch        *b_b5_nd0;   //!
   TBranch        *b_b5_d0_doca;   //!
   TBranch        *b_b5_d0_kaon_eta;   //!
   TBranch        *b_b5_d0_kaon_phi;   //!
   TBranch        *b_b5_d0_kaon_pt;   //!
   TBranch        *b_b5_d0_kaon_sip;   //!
   TBranch        *b_b5_d0_kin_cosAlphaXY;   //!
   TBranch        *b_b5_d0_kin_eta;   //!
   TBranch        *b_b5_d0_kin_lxy;   //!
   TBranch        *b_b5_d0_kin_mass;   //!
   TBranch        *b_b5_d0_kin_massErr;   //!
   TBranch        *b_b5_d0_kin_phi;   //!
   TBranch        *b_b5_d0_kin_pt;   //!
   TBranch        *b_b5_d0_kin_sipBS;   //!
   TBranch        *b_b5_d0_kin_sipPV;   //!
   TBranch        *b_b5_d0_kin_slxy;   //!
   TBranch        *b_b5_d0_kin_vtx_chi2dof;   //!
   TBranch        *b_b5_d0_kin_vtx_prob;   //!
   TBranch        *b_b5_d0_mass;   //!
   TBranch        *b_b5_d0_pion_eta;   //!
   TBranch        *b_b5_d0_pion_phi;   //!
   TBranch        *b_b5_d0_pion_pt;   //!
   TBranch        *b_b5_d0_pion_sip;   //!
   TBranch        *b_b5_d0_kaon_mu_index;   //!
   TBranch        *b_b5_d0_kin_valid;   //!
   TBranch        *b_b5_d0_pion_mu_index;   //!
   TBranch        *b_b5_nks;   //!
   TBranch        *b_b5_ks_doca;   //!
   TBranch        *b_b5_ks_kin_cosAlphaXY;   //!
   TBranch        *b_b5_ks_kin_eta;   //!
   TBranch        *b_b5_ks_kin_lxy;   //!
   TBranch        *b_b5_ks_kin_mass;   //!
   TBranch        *b_b5_ks_kin_massErr;   //!
   TBranch        *b_b5_ks_kin_phi;   //!
   TBranch        *b_b5_ks_kin_pt;   //!
   TBranch        *b_b5_ks_kin_sipBS;   //!
   TBranch        *b_b5_ks_kin_sipPV;   //!
   TBranch        *b_b5_ks_kin_slxy;   //!
   TBranch        *b_b5_ks_kin_vtx_chi2dof;   //!
   TBranch        *b_b5_ks_kin_vtx_prob;   //!
   TBranch        *b_b5_ks_mass;   //!
   TBranch        *b_b5_ks_trk1_eta;   //!
   TBranch        *b_b5_ks_trk1_phi;   //!
   TBranch        *b_b5_ks_trk1_pt;   //!
   TBranch        *b_b5_ks_trk1_sip;   //!
   TBranch        *b_b5_ks_trk2_eta;   //!
   TBranch        *b_b5_ks_trk2_phi;   //!
   TBranch        *b_b5_ks_trk2_pt;   //!
   TBranch        *b_b5_ks_trk2_sip;   //!
   TBranch        *b_b5_ks_kin_valid;   //!
   TBranch        *b_b5_ks_trk1_mu_index;   //!
   TBranch        *b_b5_ks_trk2_mu_index;   //!
   TBranch        *b_b5_nlambda;   //!
   TBranch        *b_b5_lambda_doca;   //!
   TBranch        *b_b5_lambda_kin_cosAlphaXY;   //!
   TBranch        *b_b5_lambda_kin_eta;   //!
   TBranch        *b_b5_lambda_kin_lxy;   //!
   TBranch        *b_b5_lambda_kin_mass;   //!
   TBranch        *b_b5_lambda_kin_massErr;   //!
   TBranch        *b_b5_lambda_kin_phi;   //!
   TBranch        *b_b5_lambda_kin_pt;   //!
   TBranch        *b_b5_lambda_kin_sipBS;   //!
   TBranch        *b_b5_lambda_kin_sipPV;   //!
   TBranch        *b_b5_lambda_kin_slxy;   //!
   TBranch        *b_b5_lambda_kin_vtx_chi2dof;   //!
   TBranch        *b_b5_lambda_kin_vtx_prob;   //!
   TBranch        *b_b5_lambda_mass;   //!
   TBranch        *b_b5_lambda_pion_eta;   //!
   TBranch        *b_b5_lambda_pion_phi;   //!
   TBranch        *b_b5_lambda_pion_pt;   //!
   TBranch        *b_b5_lambda_pion_sip;   //!
   TBranch        *b_b5_lambda_proton_eta;   //!
   TBranch        *b_b5_lambda_proton_phi;   //!
   TBranch        *b_b5_lambda_proton_pt;   //!
   TBranch        *b_b5_lambda_proton_sip;   //!
   TBranch        *b_b5_lambda_kin_valid;   //!
   TBranch        *b_b5_lambda_pion_mu_index;   //!
   TBranch        *b_b5_lambda_proton_mu_index;   //!
   TBranch        *b_b5_nphi;   //!
   TBranch        *b_b5_phi_doca;   //!
   TBranch        *b_b5_phi_ds_cosAlphaXY;   //!
   TBranch        *b_b5_phi_ds_eta;   //!
   TBranch        *b_b5_phi_ds_lxy;   //!
   TBranch        *b_b5_phi_ds_mass;   //!
   TBranch        *b_b5_phi_ds_massErr;   //!
   TBranch        *b_b5_phi_ds_phi;   //!
   TBranch        *b_b5_phi_ds_pion_eta;   //!
   TBranch        *b_b5_phi_ds_pion_mu_index;   //!
   TBranch        *b_b5_phi_ds_pion_phi;   //!
   TBranch        *b_b5_phi_ds_pion_pt;   //!
   TBranch        *b_b5_phi_ds_pt;   //!
   TBranch        *b_b5_phi_ds_sipBS;   //!
   TBranch        *b_b5_phi_ds_sipPV;   //!
   TBranch        *b_b5_phi_ds_slxy;   //!
   TBranch        *b_b5_phi_ds_vtx_chi2dof;   //!
   TBranch        *b_b5_phi_ds_vtx_prob;   //!
   TBranch        *b_b5_phi_kin_cosAlphaXY;   //!
   TBranch        *b_b5_phi_kin_eta;   //!
   TBranch        *b_b5_phi_kin_lxy;   //!
   TBranch        *b_b5_phi_kin_mass;   //!
   TBranch        *b_b5_phi_kin_massErr;   //!
   TBranch        *b_b5_phi_kin_phi;   //!
   TBranch        *b_b5_phi_kin_pt;   //!
   TBranch        *b_b5_phi_kin_sipBS;   //!
   TBranch        *b_b5_phi_kin_sipPV;   //!
   TBranch        *b_b5_phi_kin_slxy;   //!
   TBranch        *b_b5_phi_kin_vtx_chi2dof;   //!
   TBranch        *b_b5_phi_kin_vtx_prob;   //!
   TBranch        *b_b5_phi_mass;   //!
   TBranch        *b_b5_phi_trk1_eta;   //!
   TBranch        *b_b5_phi_trk1_phi;   //!
   TBranch        *b_b5_phi_trk1_pt;   //!
   TBranch        *b_b5_phi_trk1_sip;   //!
   TBranch        *b_b5_phi_trk2_eta;   //!
   TBranch        *b_b5_phi_trk2_phi;   //!
   TBranch        *b_b5_phi_trk2_pt;   //!
   TBranch        *b_b5_phi_trk2_sip;   //!
   TBranch        *b_b5_phi_kin_valid;   //!
   TBranch        *b_b5_phi_trk1_mu_index;   //!
   TBranch        *b_b5_phi_trk2_mu_index;   //!
   TBranch        *b_b5_CaloMET_phi;   //!
   TBranch        *b_b5_CaloMET_pt;   //!
   TBranch        *b_b5_CaloMET_sumEt;   //!
   TBranch        *b_b5_ChsMET_phi;   //!
   TBranch        *b_b5_ChsMET_pt;   //!
   TBranch        *b_b5_ChsMET_sumEt;   //!
   TBranch        *b_b5_nCorrT1METJet;   //!
   TBranch        *b_b5_CorrT1METJet_area;   //!
   TBranch        *b_b5_CorrT1METJet_eta;   //!
   TBranch        *b_b5_CorrT1METJet_muonSubtrFactor;   //!
   TBranch        *b_b5_CorrT1METJet_phi;   //!
   TBranch        *b_b5_CorrT1METJet_rawPt;   //!
   TBranch        *b_b5_DeepMETResolutionTune_phi;   //!
   TBranch        *b_b5_DeepMETResolutionTune_pt;   //!
   TBranch        *b_b5_DeepMETResponseTune_phi;   //!
   TBranch        *b_b5_DeepMETResponseTune_pt;   //!
   TBranch        *b_b5_nElectron;   //!
   TBranch        *b_b5_Electron_deltaEtaSC;   //!
   TBranch        *b_b5_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_b5_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_b5_Electron_dr03TkSumPt;   //!
   TBranch        *b_b5_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_b5_Electron_dxy;   //!
   TBranch        *b_b5_Electron_dxyErr;   //!
   TBranch        *b_b5_Electron_dz;   //!
   TBranch        *b_b5_Electron_dzErr;   //!
   TBranch        *b_b5_Electron_eCorr;   //!
   TBranch        *b_b5_Electron_eInvMinusPInv;   //!
   TBranch        *b_b5_Electron_energyErr;   //!
   TBranch        *b_b5_Electron_eta;   //!
   TBranch        *b_b5_Electron_hoe;   //!
   TBranch        *b_b5_Electron_ip3d;   //!
   TBranch        *b_b5_Electron_jetPtRelv2;   //!
   TBranch        *b_b5_Electron_jetRelIso;   //!
   TBranch        *b_b5_Electron_mass;   //!
   TBranch        *b_b5_Electron_miniPFRelIso_all;   //!
   TBranch        *b_b5_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_b5_Electron_mvaFall17V1Iso;   //!
   TBranch        *b_b5_Electron_mvaFall17V1noIso;   //!
   TBranch        *b_b5_Electron_mvaFall17V2Iso;   //!
   TBranch        *b_b5_Electron_mvaFall17V2noIso;   //!
   TBranch        *b_b5_Electron_pfRelIso03_all;   //!
   TBranch        *b_b5_Electron_pfRelIso03_chg;   //!
   TBranch        *b_b5_Electron_phi;   //!
   TBranch        *b_b5_Electron_pt;   //!
   TBranch        *b_b5_Electron_r9;   //!
   TBranch        *b_b5_Electron_scEtOverPt;   //!
   TBranch        *b_b5_Electron_sieie;   //!
   TBranch        *b_b5_Electron_sip3d;   //!
   TBranch        *b_b5_Electron_mvaTTH;   //!
   TBranch        *b_b5_Electron_charge;   //!
   TBranch        *b_b5_Electron_cutBased;   //!
   TBranch        *b_b5_Electron_cutBased_Fall17_V1;   //!
   TBranch        *b_b5_Electron_jetIdx;   //!
   TBranch        *b_b5_Electron_pdgId;   //!
   TBranch        *b_b5_Electron_photonIdx;   //!
   TBranch        *b_b5_Electron_tightCharge;   //!
   TBranch        *b_b5_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_b5_Electron_vidNestedWPBitmapHEEP;   //!
   TBranch        *b_b5_Electron_convVeto;   //!
   TBranch        *b_b5_Electron_cutBased_HEEP;   //!
   TBranch        *b_b5_Electron_isPFcand;   //!
   TBranch        *b_b5_Electron_jetNDauCharged;   //!
   TBranch        *b_b5_Electron_lostHits;   //!
   TBranch        *b_b5_Electron_mvaFall17V1Iso_WP80;   //!
   TBranch        *b_b5_Electron_mvaFall17V1Iso_WP90;   //!
   TBranch        *b_b5_Electron_mvaFall17V1Iso_WPL;   //!
   TBranch        *b_b5_Electron_mvaFall17V1noIso_WP80;   //!
   TBranch        *b_b5_Electron_mvaFall17V1noIso_WP90;   //!
   TBranch        *b_b5_Electron_mvaFall17V1noIso_WPL;   //!
   TBranch        *b_b5_Electron_mvaFall17V2Iso_WP80;   //!
   TBranch        *b_b5_Electron_mvaFall17V2Iso_WP90;   //!
   TBranch        *b_b5_Electron_mvaFall17V2Iso_WPL;   //!
   TBranch        *b_b5_Electron_mvaFall17V2noIso_WP80;   //!
   TBranch        *b_b5_Electron_mvaFall17V2noIso_WP90;   //!
   TBranch        *b_b5_Electron_mvaFall17V2noIso_WPL;   //!
   TBranch        *b_b5_Electron_seedGain;   //!
   TBranch        *b_b5_nFatJet;   //!
   TBranch        *b_b5_FatJet_area;   //!
   TBranch        *b_b5_FatJet_btagCMVA;   //!
   TBranch        *b_b5_FatJet_btagCSVV2;   //!
   TBranch        *b_b5_FatJet_btagDDBvL;   //!
   TBranch        *b_b5_FatJet_btagDDBvLV2;   //!
   TBranch        *b_b5_FatJet_btagDDBvL_noMD;   //!
   TBranch        *b_b5_FatJet_btagDDCvB;   //!
   TBranch        *b_b5_FatJet_btagDDCvBV2;   //!
   TBranch        *b_b5_FatJet_btagDDCvB_noMD;   //!
   TBranch        *b_b5_FatJet_btagDDCvL;   //!
   TBranch        *b_b5_FatJet_btagDDCvLV2;   //!
   TBranch        *b_b5_FatJet_btagDDCvL_noMD;   //!
   TBranch        *b_b5_FatJet_btagDeepB;   //!
   TBranch        *b_b5_FatJet_btagHbb;   //!
   TBranch        *b_b5_FatJet_deepTagMD_H4qvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_HbbvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_TvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_WvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_ZHbbvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_ZHccvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_ZbbvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_ZvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTagMD_bbvsLight;   //!
   TBranch        *b_b5_FatJet_deepTagMD_ccvsLight;   //!
   TBranch        *b_b5_FatJet_deepTag_H;   //!
   TBranch        *b_b5_FatJet_deepTag_QCD;   //!
   TBranch        *b_b5_FatJet_deepTag_QCDothers;   //!
   TBranch        *b_b5_FatJet_deepTag_TvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTag_WvsQCD;   //!
   TBranch        *b_b5_FatJet_deepTag_ZvsQCD;   //!
   TBranch        *b_b5_FatJet_eta;   //!
   TBranch        *b_b5_FatJet_mass;   //!
   TBranch        *b_b5_FatJet_msoftdrop;   //!
   TBranch        *b_b5_FatJet_n2b1;   //!
   TBranch        *b_b5_FatJet_n3b1;   //!
   TBranch        *b_b5_FatJet_particleNetMD_QCD;   //!
   TBranch        *b_b5_FatJet_particleNetMD_Xbb;   //!
   TBranch        *b_b5_FatJet_particleNetMD_Xcc;   //!
   TBranch        *b_b5_FatJet_particleNetMD_Xqq;   //!
   TBranch        *b_b5_FatJet_particleNet_H4qvsQCD;   //!
   TBranch        *b_b5_FatJet_particleNet_HbbvsQCD;   //!
   TBranch        *b_b5_FatJet_particleNet_HccvsQCD;   //!
   TBranch        *b_b5_FatJet_particleNet_QCD;   //!
   TBranch        *b_b5_FatJet_particleNet_TvsQCD;   //!
   TBranch        *b_b5_FatJet_particleNet_WvsQCD;   //!
   TBranch        *b_b5_FatJet_particleNet_ZvsQCD;   //!
   TBranch        *b_b5_FatJet_phi;   //!
   TBranch        *b_b5_FatJet_pt;   //!
   TBranch        *b_b5_FatJet_rawFactor;   //!
   TBranch        *b_b5_FatJet_tau1;   //!
   TBranch        *b_b5_FatJet_tau2;   //!
   TBranch        *b_b5_FatJet_tau3;   //!
   TBranch        *b_b5_FatJet_tau4;   //!
   TBranch        *b_b5_FatJet_lsf3;   //!
   TBranch        *b_b5_FatJet_jetId;   //!
   TBranch        *b_b5_FatJet_subJetIdx1;   //!
   TBranch        *b_b5_FatJet_subJetIdx2;   //!
   TBranch        *b_b5_FatJet_electronIdx3SJ;   //!
   TBranch        *b_b5_FatJet_muonIdx3SJ;   //!
   TBranch        *b_b5_nFsrPhoton;   //!
   TBranch        *b_b5_FsrPhoton_dROverEt2;   //!
   TBranch        *b_b5_FsrPhoton_eta;   //!
   TBranch        *b_b5_FsrPhoton_phi;   //!
   TBranch        *b_b5_FsrPhoton_pt;   //!
   TBranch        *b_b5_FsrPhoton_relIso03;   //!
   TBranch        *b_b5_FsrPhoton_muonIdx;   //!
   TBranch        *b_b5_nIsoTrack;   //!
   TBranch        *b_b5_IsoTrack_dxy;   //!
   TBranch        *b_b5_IsoTrack_dz;   //!
   TBranch        *b_b5_IsoTrack_eta;   //!
   TBranch        *b_b5_IsoTrack_pfRelIso03_all;   //!
   TBranch        *b_b5_IsoTrack_pfRelIso03_chg;   //!
   TBranch        *b_b5_IsoTrack_phi;   //!
   TBranch        *b_b5_IsoTrack_pt;   //!
   TBranch        *b_b5_IsoTrack_miniPFRelIso_all;   //!
   TBranch        *b_b5_IsoTrack_miniPFRelIso_chg;   //!
   TBranch        *b_b5_IsoTrack_fromPV;   //!
   TBranch        *b_b5_IsoTrack_pdgId;   //!
   TBranch        *b_b5_IsoTrack_isHighPurityTrack;   //!
   TBranch        *b_b5_IsoTrack_isPFcand;   //!
   TBranch        *b_b5_IsoTrack_isFromLostTrack;   //!
   TBranch        *b_b5_nJet;   //!
   TBranch        *b_b5_Jet_area;   //!
   TBranch        *b_b5_Jet_btagCMVA;   //!
   TBranch        *b_b5_Jet_btagCSVV2;   //!
   TBranch        *b_b5_Jet_btagDeepB;   //!
   TBranch        *b_b5_Jet_btagDeepC;   //!
   TBranch        *b_b5_Jet_btagDeepCvB;   //!
   TBranch        *b_b5_Jet_btagDeepCvL;   //!
   TBranch        *b_b5_Jet_btagDeepFlavB;   //!
   TBranch        *b_b5_Jet_btagDeepFlavC;   //!
   TBranch        *b_b5_Jet_btagDeepFlavCvB;   //!
   TBranch        *b_b5_Jet_btagDeepFlavCvL;   //!
   TBranch        *b_b5_Jet_btagDeepFlavQG;   //!
   TBranch        *b_b5_Jet_chEmEF;   //!
   TBranch        *b_b5_Jet_chFPV0EF;   //!
   TBranch        *b_b5_Jet_chFPV1EF;   //!
   TBranch        *b_b5_Jet_chFPV2EF;   //!
   TBranch        *b_b5_Jet_chFPV3EF;   //!
   TBranch        *b_b5_Jet_chHEF;   //!
   TBranch        *b_b5_Jet_eta;   //!
   TBranch        *b_b5_Jet_hfsigmaEtaEta;   //!
   TBranch        *b_b5_Jet_hfsigmaPhiPhi;   //!
   TBranch        *b_b5_Jet_mass;   //!
   TBranch        *b_b5_Jet_muEF;   //!
   TBranch        *b_b5_Jet_muonSubtrFactor;   //!
   TBranch        *b_b5_Jet_neEmEF;   //!
   TBranch        *b_b5_Jet_neHEF;   //!
   TBranch        *b_b5_Jet_phi;   //!
   TBranch        *b_b5_Jet_pt;   //!
   TBranch        *b_b5_Jet_puIdDisc;   //!
   TBranch        *b_b5_Jet_qgl;   //!
   TBranch        *b_b5_Jet_rawFactor;   //!
   TBranch        *b_b5_Jet_bRegCorr;   //!
   TBranch        *b_b5_Jet_bRegRes;   //!
   TBranch        *b_b5_Jet_cRegCorr;   //!
   TBranch        *b_b5_Jet_cRegRes;   //!
   TBranch        *b_b5_Jet_electronIdx1;   //!
   TBranch        *b_b5_Jet_electronIdx2;   //!
   TBranch        *b_b5_Jet_hfadjacentEtaStripsSize;   //!
   TBranch        *b_b5_Jet_hfcentralEtaStripSize;   //!
   TBranch        *b_b5_Jet_jetId;   //!
   TBranch        *b_b5_Jet_muonIdx1;   //!
   TBranch        *b_b5_Jet_muonIdx2;   //!
   TBranch        *b_b5_Jet_nElectrons;   //!
   TBranch        *b_b5_Jet_nMuons;   //!
   TBranch        *b_b5_Jet_puId;   //!
   TBranch        *b_b5_Jet_nConstituents;   //!
   TBranch        *b_b5_L1PreFiringWeight_Dn;   //!
   TBranch        *b_b5_L1PreFiringWeight_Nom;   //!
   TBranch        *b_b5_L1PreFiringWeight_Up;   //!
   TBranch        *b_b5_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_b5_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_b5_MET_covXX;   //!
   TBranch        *b_b5_MET_covXY;   //!
   TBranch        *b_b5_MET_covYY;   //!
   TBranch        *b_b5_MET_phi;   //!
   TBranch        *b_b5_MET_pt;   //!
   TBranch        *b_b5_MET_significance;   //!
   TBranch        *b_b5_MET_sumEt;   //!
   TBranch        *b_b5_MET_sumPtUnclustered;   //!
   TBranch        *b_b5_nMuon;   //!
   TBranch        *b_b5_Muon_dxy;   //!
   TBranch        *b_b5_Muon_dxyErr;   //!
   TBranch        *b_b5_Muon_dxybs;   //!
   TBranch        *b_b5_Muon_dz;   //!
   TBranch        *b_b5_Muon_dzErr;   //!
   TBranch        *b_b5_Muon_eta;   //!
   TBranch        *b_b5_Muon_ip3d;   //!
   TBranch        *b_b5_Muon_jetPtRelv2;   //!
   TBranch        *b_b5_Muon_jetRelIso;   //!
   TBranch        *b_b5_Muon_mass;   //!
   TBranch        *b_b5_Muon_miniPFRelIso_all;   //!
   TBranch        *b_b5_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_b5_Muon_pfRelIso03_all;   //!
   TBranch        *b_b5_Muon_pfRelIso03_chg;   //!
   TBranch        *b_b5_Muon_pfRelIso04_all;   //!
   TBranch        *b_b5_Muon_phi;   //!
   TBranch        *b_b5_Muon_pt;   //!
   TBranch        *b_b5_Muon_ptErr;   //!
   TBranch        *b_b5_Muon_segmentComp;   //!
   TBranch        *b_b5_Muon_sip3d;   //!
   TBranch        *b_b5_Muon_softMva;   //!
   TBranch        *b_b5_Muon_tkRelIso;   //!
   TBranch        *b_b5_Muon_tunepRelPt;   //!
   TBranch        *b_b5_Muon_mvaLowPt;   //!
   TBranch        *b_b5_Muon_mvaTTH;   //!
   TBranch        *b_b5_Muon_charge;   //!
   TBranch        *b_b5_Muon_jetIdx;   //!
   TBranch        *b_b5_Muon_nStations;   //!
   TBranch        *b_b5_Muon_nTrackerLayers;   //!
   TBranch        *b_b5_Muon_pdgId;   //!
   TBranch        *b_b5_Muon_tightCharge;   //!
   TBranch        *b_b5_Muon_fsrPhotonIdx;   //!
   TBranch        *b_b5_Muon_highPtId;   //!
   TBranch        *b_b5_Muon_highPurity;   //!
   TBranch        *b_b5_Muon_inTimeMuon;   //!
   TBranch        *b_b5_Muon_isGlobal;   //!
   TBranch        *b_b5_Muon_isPFcand;   //!
   TBranch        *b_b5_Muon_isTracker;   //!
   TBranch        *b_b5_Muon_jetNDauCharged;   //!
   TBranch        *b_b5_Muon_looseId;   //!
   TBranch        *b_b5_Muon_mediumId;   //!
   TBranch        *b_b5_Muon_mediumPromptId;   //!
   TBranch        *b_b5_Muon_miniIsoId;   //!
   TBranch        *b_b5_Muon_multiIsoId;   //!
   TBranch        *b_b5_Muon_mvaId;   //!
   TBranch        *b_b5_Muon_mvaLowPtId;   //!
   TBranch        *b_b5_Muon_pfIsoId;   //!
   TBranch        *b_b5_Muon_puppiIsoId;   //!
   TBranch        *b_b5_Muon_softId;   //!
   TBranch        *b_b5_Muon_softMvaId;   //!
   TBranch        *b_b5_Muon_tightId;   //!
   TBranch        *b_b5_Muon_tkIsoId;   //!
   TBranch        *b_b5_Muon_triggerIdLoose;   //!
   TBranch        *b_b5_nPhoton;   //!
   TBranch        *b_b5_Photon_eCorr;   //!
   TBranch        *b_b5_Photon_energyErr;   //!
   TBranch        *b_b5_Photon_eta;   //!
   TBranch        *b_b5_Photon_hoe;   //!
   TBranch        *b_b5_Photon_mass;   //!
   TBranch        *b_b5_Photon_mvaID;   //!
   TBranch        *b_b5_Photon_mvaID_Fall17V1p1;   //!
   TBranch        *b_b5_Photon_pfRelIso03_all;   //!
   TBranch        *b_b5_Photon_pfRelIso03_chg;   //!
   TBranch        *b_b5_Photon_phi;   //!
   TBranch        *b_b5_Photon_pt;   //!
   TBranch        *b_b5_Photon_r9;   //!
   TBranch        *b_b5_Photon_sieie;   //!
   TBranch        *b_b5_Photon_charge;   //!
   TBranch        *b_b5_Photon_cutBased;   //!
   TBranch        *b_b5_Photon_cutBased_Fall17V1Bitmap;   //!
   TBranch        *b_b5_Photon_electronIdx;   //!
   TBranch        *b_b5_Photon_jetIdx;   //!
   TBranch        *b_b5_Photon_pdgId;   //!
   TBranch        *b_b5_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_b5_Photon_electronVeto;   //!
   TBranch        *b_b5_Photon_isScEtaEB;   //!
   TBranch        *b_b5_Photon_isScEtaEE;   //!
   TBranch        *b_b5_Photon_mvaID_WP80;   //!
   TBranch        *b_b5_Photon_mvaID_WP90;   //!
   TBranch        *b_b5_Photon_pixelSeed;   //!
   TBranch        *b_b5_Photon_seedGain;   //!
   TBranch        *b_b5_PuppiMET_phi;   //!
   TBranch        *b_b5_PuppiMET_phiJERDown;   //!
   TBranch        *b_b5_PuppiMET_phiJERUp;   //!
   TBranch        *b_b5_PuppiMET_phiJESDown;   //!
   TBranch        *b_b5_PuppiMET_phiJESUp;   //!
   TBranch        *b_b5_PuppiMET_phiUnclusteredDown;   //!
   TBranch        *b_b5_PuppiMET_phiUnclusteredUp;   //!
   TBranch        *b_b5_PuppiMET_pt;   //!
   TBranch        *b_b5_PuppiMET_ptJERDown;   //!
   TBranch        *b_b5_PuppiMET_ptJERUp;   //!
   TBranch        *b_b5_PuppiMET_ptJESDown;   //!
   TBranch        *b_b5_PuppiMET_ptJESUp;   //!
   TBranch        *b_b5_PuppiMET_ptUnclusteredDown;   //!
   TBranch        *b_b5_PuppiMET_ptUnclusteredUp;   //!
   TBranch        *b_b5_PuppiMET_sumEt;   //!
   TBranch        *b_b5_RawMET_phi;   //!
   TBranch        *b_b5_RawMET_pt;   //!
   TBranch        *b_b5_RawMET_sumEt;   //!
   TBranch        *b_b5_RawPuppiMET_phi;   //!
   TBranch        *b_b5_RawPuppiMET_pt;   //!
   TBranch        *b_b5_RawPuppiMET_sumEt;   //!
   TBranch        *b_b5_fixedGridRhoFastjetAll;   //!
   TBranch        *b_b5_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_b5_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_b5_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_b5_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_b5_nSoftActivityJet;   //!
   TBranch        *b_b5_SoftActivityJet_eta;   //!
   TBranch        *b_b5_SoftActivityJet_phi;   //!
   TBranch        *b_b5_SoftActivityJet_pt;   //!
   TBranch        *b_b5_SoftActivityJetHT;   //!
   TBranch        *b_b5_SoftActivityJetHT10;   //!
   TBranch        *b_b5_SoftActivityJetHT2;   //!
   TBranch        *b_b5_SoftActivityJetHT5;   //!
   TBranch        *b_b5_SoftActivityJetNjets10;   //!
   TBranch        *b_b5_SoftActivityJetNjets2;   //!
   TBranch        *b_b5_SoftActivityJetNjets5;   //!
   TBranch        *b_b5_nSubJet;   //!
   TBranch        *b_b5_SubJet_btagCMVA;   //!
   TBranch        *b_b5_SubJet_btagCSVV2;   //!
   TBranch        *b_b5_SubJet_btagDeepB;   //!
   TBranch        *b_b5_SubJet_eta;   //!
   TBranch        *b_b5_SubJet_mass;   //!
   TBranch        *b_b5_SubJet_n2b1;   //!
   TBranch        *b_b5_SubJet_n3b1;   //!
   TBranch        *b_b5_SubJet_phi;   //!
   TBranch        *b_b5_SubJet_pt;   //!
   TBranch        *b_b5_SubJet_rawFactor;   //!
   TBranch        *b_b5_SubJet_tau1;   //!
   TBranch        *b_b5_SubJet_tau2;   //!
   TBranch        *b_b5_SubJet_tau3;   //!
   TBranch        *b_b5_SubJet_tau4;   //!
   TBranch        *b_b5_nTau;   //!
   TBranch        *b_b5_Tau_chargedIso;   //!
   TBranch        *b_b5_Tau_dxy;   //!
   TBranch        *b_b5_Tau_dz;   //!
   TBranch        *b_b5_Tau_eta;   //!
   TBranch        *b_b5_Tau_leadTkDeltaEta;   //!
   TBranch        *b_b5_Tau_leadTkDeltaPhi;   //!
   TBranch        *b_b5_Tau_leadTkPtOverTauPt;   //!
   TBranch        *b_b5_Tau_mass;   //!
   TBranch        *b_b5_Tau_neutralIso;   //!
   TBranch        *b_b5_Tau_phi;   //!
   TBranch        *b_b5_Tau_photonsOutsideSignalCone;   //!
   TBranch        *b_b5_Tau_pt;   //!
   TBranch        *b_b5_Tau_puCorr;   //!
   TBranch        *b_b5_Tau_rawAntiEle;   //!
   TBranch        *b_b5_Tau_rawAntiEle2018;   //!
   TBranch        *b_b5_Tau_rawDeepTau2017v2p1VSe;   //!
   TBranch        *b_b5_Tau_rawDeepTau2017v2p1VSjet;   //!
   TBranch        *b_b5_Tau_rawDeepTau2017v2p1VSmu;   //!
   TBranch        *b_b5_Tau_rawIso;   //!
   TBranch        *b_b5_Tau_rawIsodR03;   //!
   TBranch        *b_b5_Tau_rawMVAnewDM2017v2;   //!
   TBranch        *b_b5_Tau_rawMVAoldDM;   //!
   TBranch        *b_b5_Tau_rawMVAoldDM2017v1;   //!
   TBranch        *b_b5_Tau_rawMVAoldDM2017v2;   //!
   TBranch        *b_b5_Tau_rawMVAoldDMdR032017v2;   //!
   TBranch        *b_b5_Tau_charge;   //!
   TBranch        *b_b5_Tau_decayMode;   //!
   TBranch        *b_b5_Tau_jetIdx;   //!
   TBranch        *b_b5_Tau_rawAntiEleCat;   //!
   TBranch        *b_b5_Tau_rawAntiEleCat2018;   //!
   TBranch        *b_b5_Tau_idAntiEle;   //!
   TBranch        *b_b5_Tau_idAntiEle2018;   //!
   TBranch        *b_b5_Tau_idAntiEleDeadECal;   //!
   TBranch        *b_b5_Tau_idAntiMu;   //!
   TBranch        *b_b5_Tau_idDecayMode;   //!
   TBranch        *b_b5_Tau_idDecayModeNewDMs;   //!
   TBranch        *b_b5_Tau_idDeepTau2017v2p1VSe;   //!
   TBranch        *b_b5_Tau_idDeepTau2017v2p1VSjet;   //!
   TBranch        *b_b5_Tau_idDeepTau2017v2p1VSmu;   //!
   TBranch        *b_b5_Tau_idMVAnewDM2017v2;   //!
   TBranch        *b_b5_Tau_idMVAoldDM;   //!
   TBranch        *b_b5_Tau_idMVAoldDM2017v1;   //!
   TBranch        *b_b5_Tau_idMVAoldDM2017v2;   //!
   TBranch        *b_b5_Tau_idMVAoldDMdR032017v2;   //!
   TBranch        *b_b5_TkMET_phi;   //!
   TBranch        *b_b5_TkMET_pt;   //!
   TBranch        *b_b5_TkMET_sumEt;   //!
   TBranch        *b_b5_nTrigObj;   //!
   TBranch        *b_b5_TrigObj_pt;   //!
   TBranch        *b_b5_TrigObj_eta;   //!
   TBranch        *b_b5_TrigObj_phi;   //!
   TBranch        *b_b5_TrigObj_l1pt;   //!
   TBranch        *b_b5_TrigObj_l1pt_2;   //!
   TBranch        *b_b5_TrigObj_l2pt;   //!
   TBranch        *b_b5_TrigObj_id;   //!
   TBranch        *b_b5_TrigObj_l1iso;   //!
   TBranch        *b_b5_TrigObj_l1charge;   //!
   TBranch        *b_b5_TrigObj_filterBits;   //!
   TBranch        *b_b5_nOtherPV;   //!
   TBranch        *b_b5_OtherPV_z;   //!
   TBranch        *b_b5_PV_ndof;   //!
   TBranch        *b_b5_PV_x;   //!
   TBranch        *b_b5_PV_y;   //!
   TBranch        *b_b5_PV_z;   //!
   TBranch        *b_b5_PV_chi2;   //!
   TBranch        *b_b5_PV_score;   //!
   TBranch        *b_b5_PV_npvs;   //!
   TBranch        *b_b5_PV_npvsGood;   //!
   TBranch        *b_b5_nSV;   //!
   TBranch        *b_b5_SV_dlen;   //!
   TBranch        *b_b5_SV_dlenSig;   //!
   TBranch        *b_b5_SV_dxy;   //!
   TBranch        *b_b5_SV_dxySig;   //!
   TBranch        *b_b5_SV_pAngle;   //!
   TBranch        *b_b5_Electron_cleanmask;   //!
   TBranch        *b_b5_Jet_cleanmask;   //!
   TBranch        *b_b5_Muon_cleanmask;   //!
   TBranch        *b_b5_Photon_cleanmask;   //!
   TBranch        *b_b5_Tau_cleanmask;   //!
   TBranch        *b_b5_SV_chi2;   //!
   TBranch        *b_b5_SV_eta;   //!
   TBranch        *b_b5_SV_mass;   //!
   TBranch        *b_b5_SV_ndof;   //!
   TBranch        *b_b5_SV_phi;   //!
   TBranch        *b_b5_SV_pt;   //!
   TBranch        *b_b5_SV_x;   //!
   TBranch        *b_b5_SV_y;   //!
   TBranch        *b_b5_SV_z;   //!
   TBranch        *b_b5_SV_ntracks;   //!
   TBranch        *b_b5_L1_AlwaysTrue;   //!
   TBranch        *b_b5_L1_BPTX_AND_Ref1_VME;   //!
   TBranch        *b_b5_L1_BPTX_AND_Ref3_VME;   //!
   TBranch        *b_b5_L1_BPTX_AND_Ref4_VME;   //!
   TBranch        *b_b5_L1_BPTX_BeamGas_B1_VME;   //!
   TBranch        *b_b5_L1_BPTX_BeamGas_B2_VME;   //!
   TBranch        *b_b5_L1_BPTX_BeamGas_Ref1_VME;   //!
   TBranch        *b_b5_L1_BPTX_BeamGas_Ref2_VME;   //!
   TBranch        *b_b5_L1_BPTX_NotOR_VME;   //!
   TBranch        *b_b5_L1_BPTX_OR_Ref3_VME;   //!
   TBranch        *b_b5_L1_BPTX_OR_Ref4_VME;   //!
   TBranch        *b_b5_L1_BPTX_RefAND_VME;   //!
   TBranch        *b_b5_L1_BptxMinus;   //!
   TBranch        *b_b5_L1_BptxOR;   //!
   TBranch        *b_b5_L1_BptxPlus;   //!
   TBranch        *b_b5_L1_BptxXOR;   //!
   TBranch        *b_b5_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_b5_L1_DoubleEG8er2p5_HTT260er;   //!
   TBranch        *b_b5_L1_DoubleEG8er2p5_HTT280er;   //!
   TBranch        *b_b5_L1_DoubleEG8er2p5_HTT300er;   //!
   TBranch        *b_b5_L1_DoubleEG8er2p5_HTT320er;   //!
   TBranch        *b_b5_L1_DoubleEG8er2p5_HTT340er;   //!
   TBranch        *b_b5_L1_DoubleEG_15_10_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_20_10_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_22_10_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_25_12_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_25_14_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_27_14_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_LooseIso20_10_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_LooseIso22_10_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_LooseIso22_12_er2p5;   //!
   TBranch        *b_b5_L1_DoubleEG_LooseIso25_12_er2p5;   //!
   TBranch        *b_b5_L1_DoubleIsoTau32er2p1;   //!
   TBranch        *b_b5_L1_DoubleIsoTau34er2p1;   //!
   TBranch        *b_b5_L1_DoubleIsoTau36er2p1;   //!
   TBranch        *b_b5_L1_DoubleJet100er2p3_dEta_Max1p6;   //!
   TBranch        *b_b5_L1_DoubleJet100er2p5;   //!
   TBranch        *b_b5_L1_DoubleJet112er2p3_dEta_Max1p6;   //!
   TBranch        *b_b5_L1_DoubleJet120er2p5;   //!
   TBranch        *b_b5_L1_DoubleJet150er2p5;   //!
   TBranch        *b_b5_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5;   //!
   TBranch        *b_b5_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5;   //!
   TBranch        *b_b5_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5;   //!
   TBranch        *b_b5_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5;   //!
   TBranch        *b_b5_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5;   //!
   TBranch        *b_b5_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5;   //!
   TBranch        *b_b5_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp;   //!
   TBranch        *b_b5_L1_DoubleJet40er2p5;   //!
   TBranch        *b_b5_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_b5_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;   //!
   TBranch        *b_b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28;   //!
   TBranch        *b_b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620;   //!
   TBranch        *b_b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28;   //!
   TBranch        *b_b5_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ;   //!
   TBranch        *b_b5_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp;   //!
   TBranch        *b_b5_L1_DoubleJet_80_30_Mass_Min420_Mu8;   //!
   TBranch        *b_b5_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_b5_L1_DoubleLooseIsoEG22er2p1;   //!
   TBranch        *b_b5_L1_DoubleLooseIsoEG24er2p1;   //!
   TBranch        *b_b5_L1_DoubleMu0;   //!
   TBranch        *b_b5_L1_DoubleMu0_Mass_Min1;   //!
   TBranch        *b_b5_L1_DoubleMu0_OQ;   //!
   TBranch        *b_b5_L1_DoubleMu0_SQ;   //!
   TBranch        *b_b5_L1_DoubleMu0_SQ_OS;   //!
   TBranch        *b_b5_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8;   //!
   TBranch        *b_b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_b5_L1_DoubleMu0er1p5_SQ;   //!
   TBranch        *b_b5_L1_DoubleMu0er1p5_SQ_OS;   //!
   TBranch        *b_b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4;   //!
   TBranch        *b_b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4;   //!
   TBranch        *b_b5_L1_DoubleMu10_SQ;   //!
   TBranch        *b_b5_L1_DoubleMu18er2p1;   //!
   TBranch        *b_b5_L1_DoubleMu3_OS_DoubleEG7p5Upsilon;   //!
   TBranch        *b_b5_L1_DoubleMu3_SQ_ETMHF50_HTT60er;   //!
   TBranch        *b_b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5;   //!
   TBranch        *b_b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5;   //!
   TBranch        *b_b5_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5;   //!
   TBranch        *b_b5_L1_DoubleMu3_SQ_HTT220er;   //!
   TBranch        *b_b5_L1_DoubleMu3_SQ_HTT240er;   //!
   TBranch        *b_b5_L1_DoubleMu3_SQ_HTT260er;   //!
   TBranch        *b_b5_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8;   //!
   TBranch        *b_b5_L1_DoubleMu4_SQ_EG9er2p5;   //!
   TBranch        *b_b5_L1_DoubleMu4_SQ_OS;   //!
   TBranch        *b_b5_L1_DoubleMu4_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_b5_L1_DoubleMu4p5_SQ_OS;   //!
   TBranch        *b_b5_L1_DoubleMu4p5_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_b5_L1_DoubleMu4p5er2p0_SQ_OS;   //!
   TBranch        *b_b5_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;   //!
   TBranch        *b_b5_L1_DoubleMu5Upsilon_OS_DoubleEG3;   //!
   TBranch        *b_b5_L1_DoubleMu5_SQ_EG9er2p5;   //!
   TBranch        *b_b5_L1_DoubleMu9_SQ;   //!
   TBranch        *b_b5_L1_DoubleMu_12_5;   //!
   TBranch        *b_b5_L1_DoubleMu_15_5_SQ;   //!
   TBranch        *b_b5_L1_DoubleMu_15_7;   //!
   TBranch        *b_b5_L1_DoubleMu_15_7_Mass_Min1;   //!
   TBranch        *b_b5_L1_DoubleMu_15_7_SQ;   //!
   TBranch        *b_b5_L1_DoubleTau70er2p1;   //!
   TBranch        *b_b5_L1_ETM120;   //!
   TBranch        *b_b5_L1_ETM150;   //!
   TBranch        *b_b5_L1_ETMHF100;   //!
   TBranch        *b_b5_L1_ETMHF100_HTT60er;   //!
   TBranch        *b_b5_L1_ETMHF110;   //!
   TBranch        *b_b5_L1_ETMHF110_HTT60er;   //!
   TBranch        *b_b5_L1_ETMHF110_HTT60er_NotSecondBunchInTrain;   //!
   TBranch        *b_b5_L1_ETMHF120;   //!
   TBranch        *b_b5_L1_ETMHF120_HTT60er;   //!
   TBranch        *b_b5_L1_ETMHF120_NotSecondBunchInTrain;   //!
   TBranch        *b_b5_L1_ETMHF130;   //!
   TBranch        *b_b5_L1_ETMHF130_HTT60er;   //!
   TBranch        *b_b5_L1_ETMHF140;   //!
   TBranch        *b_b5_L1_ETMHF150;   //!
   TBranch        *b_b5_L1_ETMHF90_HTT60er;   //!
   TBranch        *b_b5_L1_ETT1200;   //!
   TBranch        *b_b5_L1_ETT1600;   //!
   TBranch        *b_b5_L1_ETT2000;   //!
   TBranch        *b_b5_L1_FirstBunchAfterTrain;   //!
   TBranch        *b_b5_L1_FirstBunchBeforeTrain;   //!
   TBranch        *b_b5_L1_FirstBunchInTrain;   //!
   TBranch        *b_b5_L1_FirstCollisionInOrbit;   //!
   TBranch        *b_b5_L1_FirstCollisionInTrain;   //!
   TBranch        *b_b5_L1_HCAL_LaserMon_Trig;   //!
   TBranch        *b_b5_L1_HCAL_LaserMon_Veto;   //!
   TBranch        *b_b5_L1_HTT120er;   //!
   TBranch        *b_b5_L1_HTT160er;   //!
   TBranch        *b_b5_L1_HTT200er;   //!
   TBranch        *b_b5_L1_HTT255er;   //!
   TBranch        *b_b5_L1_HTT280er;   //!
   TBranch        *b_b5_L1_HTT280er_QuadJet_70_55_40_35_er2p4;   //!
   TBranch        *b_b5_L1_HTT320er;   //!
   TBranch        *b_b5_L1_HTT320er_QuadJet_70_55_40_40_er2p4;   //!
   TBranch        *b_b5_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;   //!
   TBranch        *b_b5_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;   //!
   TBranch        *b_b5_L1_HTT360er;   //!
   TBranch        *b_b5_L1_HTT400er;   //!
   TBranch        *b_b5_L1_HTT450er;   //!
   TBranch        *b_b5_L1_IsoEG32er2p5_Mt40;   //!
   TBranch        *b_b5_L1_IsoEG32er2p5_Mt44;   //!
   TBranch        *b_b5_L1_IsoEG32er2p5_Mt48;   //!
   TBranch        *b_b5_L1_IsoTau40er2p1_ETMHF100;   //!
   TBranch        *b_b5_L1_IsoTau40er2p1_ETMHF110;   //!
   TBranch        *b_b5_L1_IsoTau40er2p1_ETMHF120;   //!
   TBranch        *b_b5_L1_IsoTau40er2p1_ETMHF90;   //!
   TBranch        *b_b5_L1_IsolatedBunch;   //!
   TBranch        *b_b5_L1_LastBunchInTrain;   //!
   TBranch        *b_b5_L1_LastCollisionInTrain;   //!
   TBranch        *b_b5_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;   //!
   TBranch        *b_b5_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3;   //!
   TBranch        *b_b5_L1_LooseIsoEG24er2p1_HTT100er;   //!
   TBranch        *b_b5_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;   //!
   TBranch        *b_b5_L1_LooseIsoEG26er2p1_HTT100er;   //!
   TBranch        *b_b5_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;   //!
   TBranch        *b_b5_L1_LooseIsoEG28er2p1_HTT100er;   //!
   TBranch        *b_b5_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;   //!
   TBranch        *b_b5_L1_LooseIsoEG30er2p1_HTT100er;   //!
   TBranch        *b_b5_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;   //!
   TBranch        *b_b5_L1_MinimumBiasHF0_AND_BptxAND;   //!
   TBranch        *b_b5_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;   //!
   TBranch        *b_b5_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6;   //!
   TBranch        *b_b5_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;   //!
   TBranch        *b_b5_L1_Mu18er2p1_Tau24er2p1;   //!
   TBranch        *b_b5_L1_Mu18er2p1_Tau26er2p1;   //!
   TBranch        *b_b5_L1_Mu20_EG10er2p5;   //!
   TBranch        *b_b5_L1_Mu22er2p1_IsoTau32er2p1;   //!
   TBranch        *b_b5_L1_Mu22er2p1_IsoTau34er2p1;   //!
   TBranch        *b_b5_L1_Mu22er2p1_IsoTau36er2p1;   //!
   TBranch        *b_b5_L1_Mu22er2p1_IsoTau40er2p1;   //!
   TBranch        *b_b5_L1_Mu22er2p1_Tau70er2p1;   //!
   TBranch        *b_b5_L1_Mu3_Jet120er2p5_dR_Max0p4;   //!
   TBranch        *b_b5_L1_Mu3_Jet120er2p5_dR_Max0p8;   //!
   TBranch        *b_b5_L1_Mu3_Jet16er2p5_dR_Max0p4;   //!
   TBranch        *b_b5_L1_Mu3_Jet30er2p5;   //!
   TBranch        *b_b5_L1_Mu3_Jet35er2p5_dR_Max0p4;   //!
   TBranch        *b_b5_L1_Mu3_Jet60er2p5_dR_Max0p4;   //!
   TBranch        *b_b5_L1_Mu3_Jet80er2p5_dR_Max0p4;   //!
   TBranch        *b_b5_L1_Mu3er1p5_Jet100er2p5_ETMHF40;   //!
   TBranch        *b_b5_L1_Mu3er1p5_Jet100er2p5_ETMHF50;   //!
   TBranch        *b_b5_L1_Mu5_EG23er2p5;   //!
   TBranch        *b_b5_L1_Mu5_LooseIsoEG20er2p5;   //!
   TBranch        *b_b5_L1_Mu6_DoubleEG10er2p5;   //!
   TBranch        *b_b5_L1_Mu6_DoubleEG12er2p5;   //!
   TBranch        *b_b5_L1_Mu6_DoubleEG15er2p5;   //!
   TBranch        *b_b5_L1_Mu6_DoubleEG17er2p5;   //!
   TBranch        *b_b5_L1_Mu6_HTT240er;   //!
   TBranch        *b_b5_L1_Mu6_HTT250er;   //!
   TBranch        *b_b5_L1_Mu7_EG23er2p5;   //!
   TBranch        *b_b5_L1_Mu7_LooseIsoEG20er2p5;   //!
   TBranch        *b_b5_L1_Mu7_LooseIsoEG23er2p5;   //!
   TBranch        *b_b5_L1_NotBptxOR;   //!
   TBranch        *b_b5_L1_QuadJet36er2p5_IsoTau52er2p1;   //!
   TBranch        *b_b5_L1_QuadJet60er2p5;   //!
   TBranch        *b_b5_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0;   //!
   TBranch        *b_b5_L1_QuadMu0;   //!
   TBranch        *b_b5_L1_QuadMu0_OQ;   //!
   TBranch        *b_b5_L1_QuadMu0_SQ;   //!
   TBranch        *b_b5_L1_SecondBunchInTrain;   //!
   TBranch        *b_b5_L1_SecondLastBunchInTrain;   //!
   TBranch        *b_b5_L1_SingleEG10er2p5;   //!
   TBranch        *b_b5_L1_SingleEG15er2p5;   //!
   TBranch        *b_b5_L1_SingleEG26er2p5;   //!
   TBranch        *b_b5_L1_SingleEG34er2p5;   //!
   TBranch        *b_b5_L1_SingleEG36er2p5;   //!
   TBranch        *b_b5_L1_SingleEG38er2p5;   //!
   TBranch        *b_b5_L1_SingleEG40er2p5;   //!
   TBranch        *b_b5_L1_SingleEG42er2p5;   //!
   TBranch        *b_b5_L1_SingleEG45er2p5;   //!
   TBranch        *b_b5_L1_SingleEG50;   //!
   TBranch        *b_b5_L1_SingleEG60;   //!
   TBranch        *b_b5_L1_SingleEG8er2p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG24er1p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG24er2p1;   //!
   TBranch        *b_b5_L1_SingleIsoEG26er1p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG26er2p1;   //!
   TBranch        *b_b5_L1_SingleIsoEG26er2p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG28er1p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG28er2p1;   //!
   TBranch        *b_b5_L1_SingleIsoEG28er2p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG30er2p1;   //!
   TBranch        *b_b5_L1_SingleIsoEG30er2p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG32er2p1;   //!
   TBranch        *b_b5_L1_SingleIsoEG32er2p5;   //!
   TBranch        *b_b5_L1_SingleIsoEG34er2p5;   //!
   TBranch        *b_b5_L1_SingleJet10erHE;   //!
   TBranch        *b_b5_L1_SingleJet120;   //!
   TBranch        *b_b5_L1_SingleJet120_FWD3p0;   //!
   TBranch        *b_b5_L1_SingleJet120er2p5;   //!
   TBranch        *b_b5_L1_SingleJet12erHE;   //!
   TBranch        *b_b5_L1_SingleJet140er2p5;   //!
   TBranch        *b_b5_L1_SingleJet140er2p5_ETMHF80;   //!
   TBranch        *b_b5_L1_SingleJet140er2p5_ETMHF90;   //!
   TBranch        *b_b5_L1_SingleJet160er2p5;   //!
   TBranch        *b_b5_L1_SingleJet180;   //!
   TBranch        *b_b5_L1_SingleJet180er2p5;   //!
   TBranch        *b_b5_L1_SingleJet200;   //!
   TBranch        *b_b5_L1_SingleJet20er2p5_NotBptxOR;   //!
   TBranch        *b_b5_L1_SingleJet20er2p5_NotBptxOR_3BX;   //!
   TBranch        *b_b5_L1_SingleJet35;   //!
   TBranch        *b_b5_L1_SingleJet35_FWD3p0;   //!
   TBranch        *b_b5_L1_SingleJet35er2p5;   //!
   TBranch        *b_b5_L1_SingleJet43er2p5_NotBptxOR_3BX;   //!
   TBranch        *b_b5_L1_SingleJet46er2p5_NotBptxOR_3BX;   //!
   TBranch        *b_b5_L1_SingleJet60;   //!
   TBranch        *b_b5_L1_SingleJet60_FWD3p0;   //!
   TBranch        *b_b5_L1_SingleJet60er2p5;   //!
   TBranch        *b_b5_L1_SingleJet8erHE;   //!
   TBranch        *b_b5_L1_SingleJet90;   //!
   TBranch        *b_b5_L1_SingleJet90_FWD3p0;   //!
   TBranch        *b_b5_L1_SingleJet90er2p5;   //!
   TBranch        *b_b5_L1_SingleLooseIsoEG28er1p5;   //!
   TBranch        *b_b5_L1_SingleLooseIsoEG30er1p5;   //!
   TBranch        *b_b5_L1_SingleMu0_BMTF;   //!
   TBranch        *b_b5_L1_SingleMu0_DQ;   //!
   TBranch        *b_b5_L1_SingleMu0_EMTF;   //!
   TBranch        *b_b5_L1_SingleMu0_OMTF;   //!
   TBranch        *b_b5_L1_SingleMu10er1p5;   //!
   TBranch        *b_b5_L1_SingleMu12_DQ_BMTF;   //!
   TBranch        *b_b5_L1_SingleMu12_DQ_EMTF;   //!
   TBranch        *b_b5_L1_SingleMu12_DQ_OMTF;   //!
   TBranch        *b_b5_L1_SingleMu12er1p5;   //!
   TBranch        *b_b5_L1_SingleMu14er1p5;   //!
   TBranch        *b_b5_L1_SingleMu15_DQ;   //!
   TBranch        *b_b5_L1_SingleMu16er1p5;   //!
   TBranch        *b_b5_L1_SingleMu18;   //!
   TBranch        *b_b5_L1_SingleMu18er1p5;   //!
   TBranch        *b_b5_L1_SingleMu20;   //!
   TBranch        *b_b5_L1_SingleMu22;   //!
   TBranch        *b_b5_L1_SingleMu22_BMTF;   //!
   TBranch        *b_b5_L1_SingleMu22_EMTF;   //!
   TBranch        *b_b5_L1_SingleMu22_OMTF;   //!
   TBranch        *b_b5_L1_SingleMu25;   //!
   TBranch        *b_b5_L1_SingleMu3;   //!
   TBranch        *b_b5_L1_SingleMu5;   //!
   TBranch        *b_b5_L1_SingleMu6er1p5;   //!
   TBranch        *b_b5_L1_SingleMu7;   //!
   TBranch        *b_b5_L1_SingleMu7_DQ;   //!
   TBranch        *b_b5_L1_SingleMu7er1p5;   //!
   TBranch        *b_b5_L1_SingleMu8er1p5;   //!
   TBranch        *b_b5_L1_SingleMu9er1p5;   //!
   TBranch        *b_b5_L1_SingleMuCosmics;   //!
   TBranch        *b_b5_L1_SingleMuCosmics_BMTF;   //!
   TBranch        *b_b5_L1_SingleMuCosmics_EMTF;   //!
   TBranch        *b_b5_L1_SingleMuCosmics_OMTF;   //!
   TBranch        *b_b5_L1_SingleMuOpen;   //!
   TBranch        *b_b5_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_b5_L1_SingleMuOpen_er1p1_NotBptxOR_3BX;   //!
   TBranch        *b_b5_L1_SingleMuOpen_er1p4_NotBptxOR_3BX;   //!
   TBranch        *b_b5_L1_SingleTau120er2p1;   //!
   TBranch        *b_b5_L1_SingleTau130er2p1;   //!
   TBranch        *b_b5_L1_TOTEM_1;   //!
   TBranch        *b_b5_L1_TOTEM_2;   //!
   TBranch        *b_b5_L1_TOTEM_3;   //!
   TBranch        *b_b5_L1_TOTEM_4;   //!
   TBranch        *b_b5_L1_TripleEG16er2p5;   //!
   TBranch        *b_b5_L1_TripleEG_16_12_8_er2p5;   //!
   TBranch        *b_b5_L1_TripleEG_16_15_8_er2p5;   //!
   TBranch        *b_b5_L1_TripleEG_18_17_8_er2p5;   //!
   TBranch        *b_b5_L1_TripleEG_18_18_12_er2p5;   //!
   TBranch        *b_b5_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5;   //!
   TBranch        *b_b5_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5;   //!
   TBranch        *b_b5_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5;   //!
   TBranch        *b_b5_L1_TripleMu0;   //!
   TBranch        *b_b5_L1_TripleMu0_OQ;   //!
   TBranch        *b_b5_L1_TripleMu0_SQ;   //!
   TBranch        *b_b5_L1_TripleMu3;   //!
   TBranch        *b_b5_L1_TripleMu3_SQ;   //!
   TBranch        *b_b5_L1_TripleMu_5SQ_3SQ_0OQ;   //!
   TBranch        *b_b5_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_b5_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_b5_L1_TripleMu_5_3_3;   //!
   TBranch        *b_b5_L1_TripleMu_5_3_3_SQ;   //!
   TBranch        *b_b5_L1_TripleMu_5_3p5_2p5;   //!
   TBranch        *b_b5_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_b5_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;   //!
   TBranch        *b_b5_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_b5_L1_TripleMu_5_5_3;   //!
   TBranch        *b_b5_L1_UnpairedBunchBptxMinus;   //!
   TBranch        *b_b5_L1_UnpairedBunchBptxPlus;   //!
   TBranch        *b_b5_L1_ZeroBias;   //!
   TBranch        *b_b5_L1_ZeroBias_copy;   //!
   TBranch        *b_b5_L1_UnprefireableEvent;   //!
   TBranch        *b_b5_Flag_HBHENoiseFilter;   //!
   TBranch        *b_b5_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_b5_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_b5_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_b5_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_b5_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_b5_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_b5_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_b5_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_b5_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_b5_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_b5_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_b5_Flag_goodVertices;   //!
   TBranch        *b_b5_Flag_eeBadScFilter;   //!
   TBranch        *b_b5_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_b5_Flag_trkPOGFilters;   //!
   TBranch        *b_b5_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_b5_Flag_muonBadTrackFilter;   //!
   TBranch        *b_b5_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_b5_Flag_BadPFMuonFilter;   //!
   TBranch        *b_b5_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_b5_Flag_hfNoisyHitsFilter;   //!
   TBranch        *b_b5_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_b5_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_b5_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_b5_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_b5_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_b5_Flag_METFilters;   //!
   TBranch        *b_b5_L1Reco_step;   //!
   TBranch        *b_b5_HLTriggerFirstPath;   //!
   TBranch        *b_b5_HLT_AK8PFJet360_TrimMass30;   //!
   TBranch        *b_b5_HLT_AK8PFJet380_TrimMass30;   //!
   TBranch        *b_b5_HLT_AK8PFJet400_TrimMass30;   //!
   TBranch        *b_b5_HLT_AK8PFJet420_TrimMass30;   //!
   TBranch        *b_b5_HLT_AK8PFHT750_TrimMass50;   //!
   TBranch        *b_b5_HLT_AK8PFHT800_TrimMass50;   //!
   TBranch        *b_b5_HLT_AK8PFHT850_TrimMass50;   //!
   TBranch        *b_b5_HLT_AK8PFHT900_TrimMass50;   //!
   TBranch        *b_b5_HLT_CaloJet500_NoJetID;   //!
   TBranch        *b_b5_HLT_CaloJet550_NoJetID;   //!
   TBranch        *b_b5_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;   //!
   TBranch        *b_b5_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;   //!
   TBranch        *b_b5_HLT_Trimuon5_3p5_2_Upsilon_Muon;   //!
   TBranch        *b_b5_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;   //!
   TBranch        *b_b5_HLT_DoubleEle25_CaloIdL_MW;   //!
   TBranch        *b_b5_HLT_DoubleEle27_CaloIdL_MW;   //!
   TBranch        *b_b5_HLT_DoubleEle33_CaloIdL_MW;   //!
   TBranch        *b_b5_HLT_DoubleEle24_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;   //!
   TBranch        *b_b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_b5_HLT_Ele27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_b5_HLT_Mu27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_b5_HLT_Mu37_Ele27_CaloIdL_MW;   //!
   TBranch        *b_b5_HLT_Mu37_TkMu27;   //!
   TBranch        *b_b5_HLT_DoubleMu4_3_Bs;   //!
   TBranch        *b_b5_HLT_DoubleMu4_3_Jpsi;   //!
   TBranch        *b_b5_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced;   //!
   TBranch        *b_b5_HLT_DoubleMu3_Trk_Tau3mu;   //!
   TBranch        *b_b5_HLT_DoubleMu3_TkMu_DsTau3Mu;   //!
   TBranch        *b_b5_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_b5_HLT_DoubleMu4_Mass3p8_DZ_PFHT350;   //!
   TBranch        *b_b5_HLT_Mu3_PFJet40;   //!
   TBranch        *b_b5_HLT_Mu7p5_L2Mu2_Jpsi;   //!
   TBranch        *b_b5_HLT_Mu7p5_L2Mu2_Upsilon;   //!
   TBranch        *b_b5_HLT_Mu7p5_Track2_Jpsi;   //!
   TBranch        *b_b5_HLT_Mu7p5_Track3p5_Jpsi;   //!
   TBranch        *b_b5_HLT_Mu7p5_Track7_Jpsi;   //!
   TBranch        *b_b5_HLT_Mu7p5_Track2_Upsilon;   //!
   TBranch        *b_b5_HLT_Mu7p5_Track3p5_Upsilon;   //!
   TBranch        *b_b5_HLT_Mu7p5_Track7_Upsilon;   //!
   TBranch        *b_b5_HLT_DoublePhoton33_CaloIdL;   //!
   TBranch        *b_b5_HLT_DoublePhoton70;   //!
   TBranch        *b_b5_HLT_DoublePhoton85;   //!
   TBranch        *b_b5_HLT_Ele20_WPTight_Gsf;   //!
   TBranch        *b_b5_HLT_Ele15_WPLoose_Gsf;   //!
   TBranch        *b_b5_HLT_Ele17_WPLoose_Gsf;   //!
   TBranch        *b_b5_HLT_Ele20_WPLoose_Gsf;   //!
   TBranch        *b_b5_HLT_Ele20_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_b5_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;   //!
   TBranch        *b_b5_HLT_Ele27_WPTight_Gsf;   //!
   TBranch        *b_b5_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_b5_HLT_Ele35_WPTight_Gsf;   //!
   TBranch        *b_b5_HLT_Ele35_WPTight_Gsf_L1EGMT;   //!
   TBranch        *b_b5_HLT_Ele38_WPTight_Gsf;   //!
   TBranch        *b_b5_HLT_Ele40_WPTight_Gsf;   //!
   TBranch        *b_b5_HLT_Ele32_WPTight_Gsf_L1DoubleEG;   //!
   TBranch        *b_b5_HLT_HT450_Beamspot;   //!
   TBranch        *b_b5_HLT_HT300_Beamspot;   //!
   TBranch        *b_b5_HLT_ZeroBias_Beamspot;   //!
   TBranch        *b_b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu20;   //!
   TBranch        *b_b5_HLT_IsoMu24;   //!
   TBranch        *b_b5_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_b5_HLT_IsoMu27;   //!
   TBranch        *b_b5_HLT_IsoMu30;   //!
   TBranch        *b_b5_HLT_UncorrectedJetE30_NoBPTX;   //!
   TBranch        *b_b5_HLT_UncorrectedJetE30_NoBPTX3BX;   //!
   TBranch        *b_b5_HLT_UncorrectedJetE60_NoBPTX3BX;   //!
   TBranch        *b_b5_HLT_UncorrectedJetE70_NoBPTX3BX;   //!
   TBranch        *b_b5_HLT_L1SingleMu18;   //!
   TBranch        *b_b5_HLT_L1SingleMu25;   //!
   TBranch        *b_b5_HLT_L2Mu10;   //!
   TBranch        *b_b5_HLT_L2Mu10_NoVertex_NoBPTX3BX;   //!
   TBranch        *b_b5_HLT_L2Mu10_NoVertex_NoBPTX;   //!
   TBranch        *b_b5_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_b5_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_b5_HLT_L2Mu50;   //!
   TBranch        *b_b5_HLT_L2Mu23NoVtx_2Cha;   //!
   TBranch        *b_b5_HLT_L2Mu23NoVtx_2Cha_CosmicSeed;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu50;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu23NoVtx_2Cha;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu25NoVtx_2Cha;   //!
   TBranch        *b_b5_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;   //!
   TBranch        *b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;   //!
   TBranch        *b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;   //!
   TBranch        *b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_b5_HLT_Mu25_TkMu0_Onia;   //!
   TBranch        *b_b5_HLT_Mu30_TkMu0_Onia;   //!
   TBranch        *b_b5_HLT_Mu20_TkMu0_Phi;   //!
   TBranch        *b_b5_HLT_Mu25_TkMu0_Phi;   //!
   TBranch        *b_b5_HLT_Mu12;   //!
   TBranch        *b_b5_HLT_Mu15;   //!
   TBranch        *b_b5_HLT_Mu20;   //!
   TBranch        *b_b5_HLT_Mu27;   //!
   TBranch        *b_b5_HLT_Mu50;   //!
   TBranch        *b_b5_HLT_Mu55;   //!
   TBranch        *b_b5_HLT_OldMu100;   //!
   TBranch        *b_b5_HLT_TkMu100;   //!
   TBranch        *b_b5_HLT_DiPFJetAve40;   //!
   TBranch        *b_b5_HLT_DiPFJetAve60;   //!
   TBranch        *b_b5_HLT_DiPFJetAve80;   //!
   TBranch        *b_b5_HLT_DiPFJetAve140;   //!
   TBranch        *b_b5_HLT_DiPFJetAve200;   //!
   TBranch        *b_b5_HLT_DiPFJetAve260;   //!
   TBranch        *b_b5_HLT_DiPFJetAve320;   //!
   TBranch        *b_b5_HLT_DiPFJetAve400;   //!
   TBranch        *b_b5_HLT_DiPFJetAve500;   //!
   TBranch        *b_b5_HLT_DiPFJetAve15_HFJEC;   //!
   TBranch        *b_b5_HLT_DiPFJetAve25_HFJEC;   //!
   TBranch        *b_b5_HLT_DiPFJetAve60_HFJEC;   //!
   TBranch        *b_b5_HLT_DiPFJetAve80_HFJEC;   //!
   TBranch        *b_b5_HLT_DiPFJetAve100_HFJEC;   //!
   TBranch        *b_b5_HLT_DiPFJetAve160_HFJEC;   //!
   TBranch        *b_b5_HLT_DiPFJetAve220_HFJEC;   //!
   TBranch        *b_b5_HLT_DiPFJetAve300_HFJEC;   //!
   TBranch        *b_b5_HLT_AK8PFJet15;   //!
   TBranch        *b_b5_HLT_AK8PFJet25;   //!
   TBranch        *b_b5_HLT_AK8PFJet40;   //!
   TBranch        *b_b5_HLT_AK8PFJet60;   //!
   TBranch        *b_b5_HLT_AK8PFJet80;   //!
   TBranch        *b_b5_HLT_AK8PFJet140;   //!
   TBranch        *b_b5_HLT_AK8PFJet200;   //!
   TBranch        *b_b5_HLT_AK8PFJet260;   //!
   TBranch        *b_b5_HLT_AK8PFJet320;   //!
   TBranch        *b_b5_HLT_AK8PFJet400;   //!
   TBranch        *b_b5_HLT_AK8PFJet450;   //!
   TBranch        *b_b5_HLT_AK8PFJet500;   //!
   TBranch        *b_b5_HLT_AK8PFJet550;   //!
   TBranch        *b_b5_HLT_PFJet15;   //!
   TBranch        *b_b5_HLT_PFJet25;   //!
   TBranch        *b_b5_HLT_PFJet40;   //!
   TBranch        *b_b5_HLT_PFJet60;   //!
   TBranch        *b_b5_HLT_PFJet80;   //!
   TBranch        *b_b5_HLT_PFJet140;   //!
   TBranch        *b_b5_HLT_PFJet200;   //!
   TBranch        *b_b5_HLT_PFJet260;   //!
   TBranch        *b_b5_HLT_PFJet320;   //!
   TBranch        *b_b5_HLT_PFJet400;   //!
   TBranch        *b_b5_HLT_PFJet450;   //!
   TBranch        *b_b5_HLT_PFJet500;   //!
   TBranch        *b_b5_HLT_PFJet550;   //!
   TBranch        *b_b5_HLT_PFJetFwd15;   //!
   TBranch        *b_b5_HLT_PFJetFwd25;   //!
   TBranch        *b_b5_HLT_PFJetFwd40;   //!
   TBranch        *b_b5_HLT_PFJetFwd60;   //!
   TBranch        *b_b5_HLT_PFJetFwd80;   //!
   TBranch        *b_b5_HLT_PFJetFwd140;   //!
   TBranch        *b_b5_HLT_PFJetFwd200;   //!
   TBranch        *b_b5_HLT_PFJetFwd260;   //!
   TBranch        *b_b5_HLT_PFJetFwd320;   //!
   TBranch        *b_b5_HLT_PFJetFwd400;   //!
   TBranch        *b_b5_HLT_PFJetFwd450;   //!
   TBranch        *b_b5_HLT_PFJetFwd500;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd15;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd25;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd40;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd60;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd80;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd140;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd200;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd260;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd320;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd400;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd450;   //!
   TBranch        *b_b5_HLT_AK8PFJetFwd500;   //!
   TBranch        *b_b5_HLT_PFHT180;   //!
   TBranch        *b_b5_HLT_PFHT250;   //!
   TBranch        *b_b5_HLT_PFHT370;   //!
   TBranch        *b_b5_HLT_PFHT430;   //!
   TBranch        *b_b5_HLT_PFHT510;   //!
   TBranch        *b_b5_HLT_PFHT590;   //!
   TBranch        *b_b5_HLT_PFHT680;   //!
   TBranch        *b_b5_HLT_PFHT780;   //!
   TBranch        *b_b5_HLT_PFHT890;   //!
   TBranch        *b_b5_HLT_PFHT1050;   //!
   TBranch        *b_b5_HLT_PFHT500_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_b5_HLT_PFHT500_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_b5_HLT_PFHT700_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_b5_HLT_PFHT700_PFMET95_PFMHT95_IDTight;   //!
   TBranch        *b_b5_HLT_PFHT800_PFMET75_PFMHT75_IDTight;   //!
   TBranch        *b_b5_HLT_PFHT800_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_b5_HLT_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_b5_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_b5_HLT_PFMET130_PFMHT130_IDTight;   //!
   TBranch        *b_b5_HLT_PFMET140_PFMHT140_IDTight;   //!
   TBranch        *b_b5_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1;   //!
   TBranch        *b_b5_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1;   //!
   TBranch        *b_b5_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1;   //!
   TBranch        *b_b5_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1;   //!
   TBranch        *b_b5_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1;   //!
   TBranch        *b_b5_HLT_PFMET120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_b5_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_b5_HLT_PFMETTypeOne110_PFMHT110_IDTight;   //!
   TBranch        *b_b5_HLT_PFMETTypeOne120_PFMHT120_IDTight;   //!
   TBranch        *b_b5_HLT_PFMETTypeOne130_PFMHT130_IDTight;   //!
   TBranch        *b_b5_HLT_PFMETTypeOne140_PFMHT140_IDTight;   //!
   TBranch        *b_b5_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
   TBranch        *b_b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_b5_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_b5_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_b5_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight;   //!
   TBranch        *b_b5_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_b5_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_b5_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_b5_HLT_L1ETMHadSeeds;   //!
   TBranch        *b_b5_HLT_CaloMHT90;   //!
   TBranch        *b_b5_HLT_CaloMET80_NotCleaned;   //!
   TBranch        *b_b5_HLT_CaloMET90_NotCleaned;   //!
   TBranch        *b_b5_HLT_CaloMET100_NotCleaned;   //!
   TBranch        *b_b5_HLT_CaloMET110_NotCleaned;   //!
   TBranch        *b_b5_HLT_CaloMET250_NotCleaned;   //!
   TBranch        *b_b5_HLT_CaloMET70_HBHECleaned;   //!
   TBranch        *b_b5_HLT_CaloMET80_HBHECleaned;   //!
   TBranch        *b_b5_HLT_CaloMET90_HBHECleaned;   //!
   TBranch        *b_b5_HLT_CaloMET100_HBHECleaned;   //!
   TBranch        *b_b5_HLT_CaloMET250_HBHECleaned;   //!
   TBranch        *b_b5_HLT_CaloMET300_HBHECleaned;   //!
   TBranch        *b_b5_HLT_CaloMET350_HBHECleaned;   //!
   TBranch        *b_b5_HLT_PFMET200_NotCleaned;   //!
   TBranch        *b_b5_HLT_PFMET200_HBHECleaned;   //!
   TBranch        *b_b5_HLT_PFMET250_HBHECleaned;   //!
   TBranch        *b_b5_HLT_PFMET300_HBHECleaned;   //!
   TBranch        *b_b5_HLT_PFMET200_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_b5_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_b5_HLT_MET105_IsoTrk50;   //!
   TBranch        *b_b5_HLT_MET120_IsoTrk50;   //!
   TBranch        *b_b5_HLT_SingleJet30_Mu12_SinglePFJet40;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_DoublePFJets40_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_DoublePFJets100_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_DoublePFJets200_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_DoublePFJets350_CaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71;   //!
   TBranch        *b_b5_HLT_Photon300_NoHE;   //!
   TBranch        *b_b5_HLT_Mu8_TrkIsoVVL;   //!
   TBranch        *b_b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;   //!
   TBranch        *b_b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_b5_HLT_Mu17_TrkIsoVVL;   //!
   TBranch        *b_b5_HLT_Mu19_TrkIsoVVL;   //!
   TBranch        *b_b5_HLT_BTagMu_AK4DiJet20_Mu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK4DiJet40_Mu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK4DiJet70_Mu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK4DiJet110_Mu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK4DiJet170_Mu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK4Jet300_Mu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK8DiJet170_Mu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK8Jet170_DoubleMu5;   //!
   TBranch        *b_b5_HLT_BTagMu_AK8Jet300_Mu5;   //!
   TBranch        *b_b5_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_b5_HLT_Mu12_DoublePhoton20;   //!
   TBranch        *b_b5_HLT_TriplePhoton_20_20_20_CaloIdLV2;   //!
   TBranch        *b_b5_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_b5_HLT_TriplePhoton_30_30_10_CaloIdLV2;   //!
   TBranch        *b_b5_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_b5_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL;   //!
   TBranch        *b_b5_HLT_Photon20;   //!
   TBranch        *b_b5_HLT_Photon33;   //!
   TBranch        *b_b5_HLT_Photon50;   //!
   TBranch        *b_b5_HLT_Photon75;   //!
   TBranch        *b_b5_HLT_Photon90;   //!
   TBranch        *b_b5_HLT_Photon120;   //!
   TBranch        *b_b5_HLT_Photon150;   //!
   TBranch        *b_b5_HLT_Photon175;   //!
   TBranch        *b_b5_HLT_Photon200;   //!
   TBranch        *b_b5_HLT_Photon100EB_TightID_TightIso;   //!
   TBranch        *b_b5_HLT_Photon110EB_TightID_TightIso;   //!
   TBranch        *b_b5_HLT_Photon120EB_TightID_TightIso;   //!
   TBranch        *b_b5_HLT_Photon100EBHE10;   //!
   TBranch        *b_b5_HLT_Photon100EEHE10;   //!
   TBranch        *b_b5_HLT_Photon100EE_TightID_TightIso;   //!
   TBranch        *b_b5_HLT_Photon50_R9Id90_HE10_IsoM;   //!
   TBranch        *b_b5_HLT_Photon75_R9Id90_HE10_IsoM;   //!
   TBranch        *b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3;   //!
   TBranch        *b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3;   //!
   TBranch        *b_b5_HLT_Photon90_R9Id90_HE10_IsoM;   //!
   TBranch        *b_b5_HLT_Photon120_R9Id90_HE10_IsoM;   //!
   TBranch        *b_b5_HLT_Photon165_R9Id90_HE10_IsoM;   //!
   TBranch        *b_b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;   //!
   TBranch        *b_b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;   //!
   TBranch        *b_b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55;   //!
   TBranch        *b_b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55;   //!
   TBranch        *b_b5_HLT_Dimuon0_Jpsi_L1_NoOS;   //!
   TBranch        *b_b5_HLT_Dimuon0_Jpsi_NoVertexing_NoOS;   //!
   TBranch        *b_b5_HLT_Dimuon0_Jpsi;   //!
   TBranch        *b_b5_HLT_Dimuon0_Jpsi_NoVertexing;   //!
   TBranch        *b_b5_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;   //!
   TBranch        *b_b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;   //!
   TBranch        *b_b5_HLT_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_L1_4p5;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_L1_5;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_L1_4p5NoOS;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0M;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_NoVertexing;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_L1_5M;   //!
   TBranch        *b_b5_HLT_Dimuon0_LowMass_L1_0er1p5R;   //!
   TBranch        *b_b5_HLT_Dimuon0_LowMass_L1_0er1p5;   //!
   TBranch        *b_b5_HLT_Dimuon0_LowMass;   //!
   TBranch        *b_b5_HLT_Dimuon0_LowMass_L1_4;   //!
   TBranch        *b_b5_HLT_Dimuon0_LowMass_L1_4R;   //!
   TBranch        *b_b5_HLT_Dimuon0_LowMass_L1_TM530;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_Muon_L1_TM0;   //!
   TBranch        *b_b5_HLT_Dimuon0_Upsilon_Muon_NoL1Mass;   //!
   TBranch        *b_b5_HLT_TripleMu_5_3_3_Mass3p8_DZ;   //!
   TBranch        *b_b5_HLT_TripleMu_10_5_5_DZ;   //!
   TBranch        *b_b5_HLT_TripleMu_12_10_5;   //!
   TBranch        *b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;   //!
   TBranch        *b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;   //!
   TBranch        *b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;   //!
   TBranch        *b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;   //!
   TBranch        *b_b5_HLT_DoubleMu3_DZ_PFMET50_PFMHT60;   //!
   TBranch        *b_b5_HLT_DoubleMu3_DZ_PFMET70_PFMHT70;   //!
   TBranch        *b_b5_HLT_DoubleMu3_DZ_PFMET90_PFMHT90;   //!
   TBranch        *b_b5_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;   //!
   TBranch        *b_b5_HLT_DoubleMu4_Jpsi_Displaced;   //!
   TBranch        *b_b5_HLT_DoubleMu4_Jpsi_NoVertexing;   //!
   TBranch        *b_b5_HLT_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_b5_HLT_DoubleMu43NoFiltersNoVtx;   //!
   TBranch        *b_b5_HLT_DoubleMu48NoFiltersNoVtx;   //!
   TBranch        *b_b5_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;   //!
   TBranch        *b_b5_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;   //!
   TBranch        *b_b5_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;   //!
   TBranch        *b_b5_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;   //!
   TBranch        *b_b5_HLT_DoubleMu33NoFiltersNoVtxDisplaced;   //!
   TBranch        *b_b5_HLT_DoubleMu40NoFiltersNoVtxDisplaced;   //!
   TBranch        *b_b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4;   //!
   TBranch        *b_b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG;   //!
   TBranch        *b_b5_HLT_HT425;   //!
   TBranch        *b_b5_HLT_HT430_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_b5_HLT_HT500_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_b5_HLT_HT430_DisplacedDijet60_DisplacedTrack;   //!
   TBranch        *b_b5_HLT_HT400_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_b5_HLT_HT650_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_b5_HLT_HT550_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_b5_HLT_DiJet110_35_Mjj650_PFMET110;   //!
   TBranch        *b_b5_HLT_DiJet110_35_Mjj650_PFMET120;   //!
   TBranch        *b_b5_HLT_DiJet110_35_Mjj650_PFMET130;   //!
   TBranch        *b_b5_HLT_TripleJet110_35_35_Mjj650_PFMET110;   //!
   TBranch        *b_b5_HLT_TripleJet110_35_35_Mjj650_PFMET120;   //!
   TBranch        *b_b5_HLT_TripleJet110_35_35_Mjj650_PFMET130;   //!
   TBranch        *b_b5_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1;   //!
   TBranch        *b_b5_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1;   //!
   TBranch        *b_b5_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1;   //!
   TBranch        *b_b5_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;   //!
   TBranch        *b_b5_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;   //!
   TBranch        *b_b5_HLT_Ele28_HighEta_SC20_Mass55;   //!
   TBranch        *b_b5_HLT_DoubleMu20_7_Mass0to30_Photon23;   //!
   TBranch        *b_b5_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;   //!
   TBranch        *b_b5_HLT_Ele15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_b5_HLT_Ele15_IsoVVVL_PFHT450;   //!
   TBranch        *b_b5_HLT_Ele50_IsoVVVL_PFHT450;   //!
   TBranch        *b_b5_HLT_Ele15_IsoVVVL_PFHT600;   //!
   TBranch        *b_b5_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;   //!
   TBranch        *b_b5_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60;   //!
   TBranch        *b_b5_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60;   //!
   TBranch        *b_b5_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5;   //!
   TBranch        *b_b5_HLT_Mu15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_b5_HLT_Mu15_IsoVVVL_PFHT450;   //!
   TBranch        *b_b5_HLT_Mu50_IsoVVVL_PFHT450;   //!
   TBranch        *b_b5_HLT_Mu15_IsoVVVL_PFHT600;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;   //!
   TBranch        *b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;   //!
   TBranch        *b_b5_HLT_Dimuon10_PsiPrime_Barrel_Seagulls;   //!
   TBranch        *b_b5_HLT_Dimuon20_Jpsi_Barrel_Seagulls;   //!
   TBranch        *b_b5_HLT_Dimuon12_Upsilon_y1p4;   //!
   TBranch        *b_b5_HLT_Dimuon14_Phi_Barrel_Seagulls;   //!
   TBranch        *b_b5_HLT_Dimuon18_PsiPrime;   //!
   TBranch        *b_b5_HLT_Dimuon25_Jpsi;   //!
   TBranch        *b_b5_HLT_Dimuon18_PsiPrime_noCorrL1;   //!
   TBranch        *b_b5_HLT_Dimuon24_Upsilon_noCorrL1;   //!
   TBranch        *b_b5_HLT_Dimuon24_Phi_noCorrL1;   //!
   TBranch        *b_b5_HLT_Dimuon25_Jpsi_noCorrL1;   //!
   TBranch        *b_b5_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;   //!
   TBranch        *b_b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   //!
   TBranch        *b_b5_HLT_DoubleIsoMu20_eta2p1;   //!
   TBranch        *b_b5_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;   //!
   TBranch        *b_b5_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx;   //!
   TBranch        *b_b5_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx;   //!
   TBranch        *b_b5_HLT_Mu8;   //!
   TBranch        *b_b5_HLT_Mu17;   //!
   TBranch        *b_b5_HLT_Mu19;   //!
   TBranch        *b_b5_HLT_Mu17_Photon30_IsoCaloId;   //!
   TBranch        *b_b5_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_b5_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_b5_HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_b5_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_b5_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_b5_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_b5_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_b5_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;   //!
   TBranch        *b_b5_HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_b5_HLT_Ele135_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_b5_HLT_Ele145_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_b5_HLT_Ele200_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_b5_HLT_Ele250_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_b5_HLT_Ele300_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5;   //!
   TBranch        *b_b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40;   //!
   TBranch        *b_b5_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2;   //!
   TBranch        *b_b5_HLT_PFHT380_SixPFJet32;   //!
   TBranch        *b_b5_HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5;   //!
   TBranch        *b_b5_HLT_PFHT430_SixPFJet40;   //!
   TBranch        *b_b5_HLT_PFHT350;   //!
   TBranch        *b_b5_HLT_PFHT350MinPFJet15;   //!
   TBranch        *b_b5_HLT_Photon60_R9Id90_CaloIdL_IsoL;   //!
   TBranch        *b_b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL;   //!
   TBranch        *b_b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15;   //!
   TBranch        *b_b5_HLT_ECALHT800;   //!
   TBranch        *b_b5_HLT_DiSC30_18_EIso_AND_HE_Mass70;   //!
   TBranch        *b_b5_HLT_Physics;   //!
   TBranch        *b_b5_HLT_Physics_part0;   //!
   TBranch        *b_b5_HLT_Physics_part1;   //!
   TBranch        *b_b5_HLT_Physics_part2;   //!
   TBranch        *b_b5_HLT_Physics_part3;   //!
   TBranch        *b_b5_HLT_Physics_part4;   //!
   TBranch        *b_b5_HLT_Physics_part5;   //!
   TBranch        *b_b5_HLT_Physics_part6;   //!
   TBranch        *b_b5_HLT_Physics_part7;   //!
   TBranch        *b_b5_HLT_Random;   //!
   TBranch        *b_b5_HLT_ZeroBias;   //!
   TBranch        *b_b5_HLT_ZeroBias_part0;   //!
   TBranch        *b_b5_HLT_ZeroBias_part1;   //!
   TBranch        *b_b5_HLT_ZeroBias_part2;   //!
   TBranch        *b_b5_HLT_ZeroBias_part3;   //!
   TBranch        *b_b5_HLT_ZeroBias_part4;   //!
   TBranch        *b_b5_HLT_ZeroBias_part5;   //!
   TBranch        *b_b5_HLT_ZeroBias_part6;   //!
   TBranch        *b_b5_HLT_ZeroBias_part7;   //!
   TBranch        *b_b5_HLT_AK4CaloJet30;   //!
   TBranch        *b_b5_HLT_AK4CaloJet40;   //!
   TBranch        *b_b5_HLT_AK4CaloJet50;   //!
   TBranch        *b_b5_HLT_AK4CaloJet80;   //!
   TBranch        *b_b5_HLT_AK4CaloJet100;   //!
   TBranch        *b_b5_HLT_AK4CaloJet120;   //!
   TBranch        *b_b5_HLT_AK4PFJet30;   //!
   TBranch        *b_b5_HLT_AK4PFJet50;   //!
   TBranch        *b_b5_HLT_AK4PFJet80;   //!
   TBranch        *b_b5_HLT_AK4PFJet100;   //!
   TBranch        *b_b5_HLT_AK4PFJet120;   //!
   TBranch        *b_b5_HLT_SinglePhoton10_Eta3p1ForPPRef;   //!
   TBranch        *b_b5_HLT_SinglePhoton20_Eta3p1ForPPRef;   //!
   TBranch        *b_b5_HLT_SinglePhoton30_Eta3p1ForPPRef;   //!
   TBranch        *b_b5_HLT_Photon20_HoverELoose;   //!
   TBranch        *b_b5_HLT_Photon30_HoverELoose;   //!
   TBranch        *b_b5_HLT_EcalCalibration;   //!
   TBranch        *b_b5_HLT_HcalCalibration;   //!
   TBranch        *b_b5_HLT_L1UnpairedBunchBptxMinus;   //!
   TBranch        *b_b5_HLT_L1UnpairedBunchBptxPlus;   //!
   TBranch        *b_b5_HLT_L1NotBptxOR;   //!
   TBranch        *b_b5_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_b5_HLT_HcalNZS;   //!
   TBranch        *b_b5_HLT_HcalPhiSym;   //!
   TBranch        *b_b5_HLT_HcalIsolatedbunch;   //!
   TBranch        *b_b5_HLT_IsoTrackHB;   //!
   TBranch        *b_b5_HLT_IsoTrackHE;   //!
   TBranch        *b_b5_HLT_ZeroBias_FirstCollisionAfterAbortGap;   //!
   TBranch        *b_b5_HLT_ZeroBias_IsolatedBunches;   //!
   TBranch        *b_b5_HLT_ZeroBias_FirstCollisionInTrain;   //!
   TBranch        *b_b5_HLT_ZeroBias_LastCollisionInTrain;   //!
   TBranch        *b_b5_HLT_ZeroBias_FirstBXAfterTrain;   //!
   TBranch        *b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1;   //!
   TBranch        *b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1;   //!
   TBranch        *b_b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1;   //!
   TBranch        *b_b5_HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1;   //!
   TBranch        *b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_b5_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   //!
   TBranch        *b_b5_HLT_Rsq0p35;   //!
   TBranch        *b_b5_HLT_Rsq0p40;   //!
   TBranch        *b_b5_HLT_RsqMR300_Rsq0p09_MR200;   //!
   TBranch        *b_b5_HLT_RsqMR320_Rsq0p09_MR200;   //!
   TBranch        *b_b5_HLT_RsqMR300_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_b5_HLT_RsqMR320_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_b5_HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1;   //!
   TBranch        *b_b5_HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1;   //!
   TBranch        *b_b5_HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1;   //!
   TBranch        *b_b5_HLT_IsoMu27_MET90;   //!
   TBranch        *b_b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1;   //!
   TBranch        *b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1;   //!
   TBranch        *b_b5_HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg;   //!
   TBranch        *b_b5_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50;   //!
   TBranch        *b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;   //!
   TBranch        *b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3;   //!
   TBranch        *b_b5_HLT_PFMET100_PFMHT100_IDTight_PFHT60;   //!
   TBranch        *b_b5_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;   //!
   TBranch        *b_b5_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60;   //!
   TBranch        *b_b5_HLT_Mu18_Mu9_SameSign;   //!
   TBranch        *b_b5_HLT_Mu18_Mu9_SameSign_DZ;   //!
   TBranch        *b_b5_HLT_Mu18_Mu9;   //!
   TBranch        *b_b5_HLT_Mu18_Mu9_DZ;   //!
   TBranch        *b_b5_HLT_Mu20_Mu10_SameSign;   //!
   TBranch        *b_b5_HLT_Mu20_Mu10_SameSign_DZ;   //!
   TBranch        *b_b5_HLT_Mu20_Mu10;   //!
   TBranch        *b_b5_HLT_Mu20_Mu10_DZ;   //!
   TBranch        *b_b5_HLT_Mu23_Mu12_SameSign;   //!
   TBranch        *b_b5_HLT_Mu23_Mu12_SameSign_DZ;   //!
   TBranch        *b_b5_HLT_Mu23_Mu12;   //!
   TBranch        *b_b5_HLT_Mu23_Mu12_DZ;   //!
   TBranch        *b_b5_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;   //!
   TBranch        *b_b5_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;   //!
   TBranch        *b_b5_HLT_DoubleMu3_DCA_PFMET50_PFMHT60;   //!
   TBranch        *b_b5_HLT_TripleMu_5_3_3_Mass3p8_DCA;   //!
   TBranch        *b_b5_HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;   //!
   TBranch        *b_b5_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;   //!
   TBranch        *b_b5_HLT_QuadPFJet105_90_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;   //!
   TBranch        *b_b5_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1;   //!
   TBranch        *b_b5_HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2;   //!
   TBranch        *b_b5_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2;   //!
   TBranch        *b_b5_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2;   //!
   TBranch        *b_b5_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2;   //!
   TBranch        *b_b5_HLT_QuadPFJet98_83_71_15;   //!
   TBranch        *b_b5_HLT_QuadPFJet103_88_75_15;   //!
   TBranch        *b_b5_HLT_QuadPFJet105_88_76_15;   //!
   TBranch        *b_b5_HLT_QuadPFJet111_90_80_15;   //!
   TBranch        *b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17;   //!
   TBranch        *b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1;   //!
   TBranch        *b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02;   //!
   TBranch        *b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2;   //!
   TBranch        *b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4;   //!
   TBranch        *b_b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55;   //!
   TBranch        *b_b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto;   //!
   TBranch        *b_b5_HLT_Mu8p5_IP3p5_part0;   //!
   TBranch        *b_b5_HLT_Mu8p5_IP3p5_part1;   //!
   TBranch        *b_b5_HLT_Mu8p5_IP3p5_part2;   //!
   TBranch        *b_b5_HLT_Mu8p5_IP3p5_part3;   //!
   TBranch        *b_b5_HLT_Mu8p5_IP3p5_part4;   //!
   TBranch        *b_b5_HLT_Mu8p5_IP3p5_part5;   //!
   TBranch        *b_b5_HLT_Mu10p5_IP3p5_part0;   //!
   TBranch        *b_b5_HLT_Mu10p5_IP3p5_part1;   //!
   TBranch        *b_b5_HLT_Mu10p5_IP3p5_part2;   //!
   TBranch        *b_b5_HLT_Mu10p5_IP3p5_part3;   //!
   TBranch        *b_b5_HLT_Mu10p5_IP3p5_part4;   //!
   TBranch        *b_b5_HLT_Mu10p5_IP3p5_part5;   //!
   TBranch        *b_b5_HLT_Mu9_IP6_part0;   //!
   TBranch        *b_b5_HLT_Mu9_IP6_part1;   //!
   TBranch        *b_b5_HLT_Mu9_IP6_part2;   //!
   TBranch        *b_b5_HLT_Mu9_IP6_part3;   //!
   TBranch        *b_b5_HLT_Mu9_IP6_part4;   //!
   TBranch        *b_b5_HLT_Mu9_IP6_part5;   //!
   TBranch        *b_b5_HLT_Mu8_IP3_part0;   //!
   TBranch        *b_b5_HLT_Mu8_IP3_part1;   //!
   TBranch        *b_b5_HLT_Mu8_IP3_part2;   //!
   TBranch        *b_b5_HLT_Mu8_IP3_part3;   //!
   TBranch        *b_b5_HLT_Mu8_IP3_part4;   //!
   TBranch        *b_b5_HLT_Mu8_IP3_part5;   //!
   TBranch        *b_b5_HLTriggerFinalPath;   //!
   TBranch        *b_bG_run;   //!
   TBranch        *b_bG_event;   //!
   TBranch        *b_bG_lumis;   //!
   TBranch        *b_bG_isData;   //!
   TBranch        *b_bG_nSC;   //!
   TBranch        *b_bG_nscE;   //!
   TBranch        *b_bG_scE;   //!
   TBranch        *b_bG_nscEt;   //!
   TBranch        *b_bG_scEt;   //!
   TBranch        *b_bG_nscRawE;   //!
   TBranch        *b_bG_scRawE;   //!
   TBranch        *b_bG_nscEta;   //!
   TBranch        *b_bG_scEta;   //!
   TBranch        *b_bG_nscPhi;   //!
   TBranch        *b_bG_scPhi;   //!
   TBranch        *b_bG_nscX;   //!
   TBranch        *b_bG_scX;   //!
   TBranch        *b_bG_nscY;   //!
   TBranch        *b_bG_scY;   //!
   TBranch        *b_bG_nscZ;   //!
   TBranch        *b_bG_scZ;   //!
   TBranch        *b_bG_nscEtaWidth;   //!
   TBranch        *b_bG_scEtaWidth;   //!
   TBranch        *b_bG_nscPhiWidth;   //!
   TBranch        *b_bG_scPhiWidth;   //!
   TBranch        *b_bG_nscRawEt;   //!
   TBranch        *b_bG_scRawEt;   //!
   TBranch        *b_bG_nscMinDrWithGsfElectornSC_;   //!
   TBranch        *b_bG_scMinDrWithGsfElectornSC_;   //!
   TBranch        *b_bG_nscFoundGsfMatch_;   //!
   TBranch        *b_bG_scFoundGsfMatch_;   //!
   TBranch        *b_bG_nscE5x5;   //!
   TBranch        *b_bG_scE5x5;   //!
   TBranch        *b_bG_nscE2x2Ratio;   //!
   TBranch        *b_bG_scE2x2Ratio;   //!
   TBranch        *b_bG_nscE3x3Ratio;   //!
   TBranch        *b_bG_scE3x3Ratio;   //!
   TBranch        *b_bG_nscEMaxRatio;   //!
   TBranch        *b_bG_scEMaxRatio;   //!
   TBranch        *b_bG_nscE2ndRatio;   //!
   TBranch        *b_bG_scE2ndRatio;   //!
   TBranch        *b_bG_nscETopRatio;   //!
   TBranch        *b_bG_scETopRatio;   //!
   TBranch        *b_bG_nscERightRatio;   //!
   TBranch        *b_bG_scERightRatio;   //!
   TBranch        *b_bG_nscEBottomRatio;   //!
   TBranch        *b_bG_scEBottomRatio;   //!
   TBranch        *b_bG_nscELeftRatio;   //!
   TBranch        *b_bG_scELeftRatio;   //!
   TBranch        *b_bG_nscE2x5MaxRatio;   //!
   TBranch        *b_bG_scE2x5MaxRatio;   //!
   TBranch        *b_bG_nscE2x5TopRatio;   //!
   TBranch        *b_bG_scE2x5TopRatio;   //!
   TBranch        *b_bG_nscE2x5RightRatio;   //!
   TBranch        *b_bG_scE2x5RightRatio;   //!
   TBranch        *b_bG_nscE2x5BottomRatio;   //!
   TBranch        *b_bG_scE2x5BottomRatio;   //!
   TBranch        *b_bG_nscE2x5LeftRatio;   //!
   TBranch        *b_bG_scE2x5LeftRatio;   //!
   TBranch        *b_bG_nscSwissCross;   //!
   TBranch        *b_bG_scSwissCross;   //!
   TBranch        *b_bG_nscR9;   //!
   TBranch        *b_bG_scR9;   //!
   TBranch        *b_bG_nscSigmaIetaIeta;   //!
   TBranch        *b_bG_scSigmaIetaIeta;   //!
   TBranch        *b_bG_nscSigmaIetaIphi;   //!
   TBranch        *b_bG_scSigmaIetaIphi;   //!
   TBranch        *b_bG_nscSigmaIphiIphi;   //!
   TBranch        *b_bG_scSigmaIphiIphi;   //!
   TBranch        *b_bG_nscFull5x5_e5x5;   //!
   TBranch        *b_bG_scFull5x5_e5x5;   //!
   TBranch        *b_bG_nscFull5x5_e2x2Ratio;   //!
   TBranch        *b_bG_scFull5x5_e2x2Ratio;   //!
   TBranch        *b_bG_nscFull5x5_e3x3Ratio;   //!
   TBranch        *b_bG_scFull5x5_e3x3Ratio;   //!
   TBranch        *b_bG_nscFull5x5_eMaxRatio;   //!
   TBranch        *b_bG_scFull5x5_eMaxRatio;   //!
   TBranch        *b_bG_nscFull5x5_e2ndRatio;   //!
   TBranch        *b_bG_scFull5x5_e2ndRatio;   //!
   TBranch        *b_bG_nscFull5x5_eTopRatio;   //!
   TBranch        *b_bG_scFull5x5_eTopRatio;   //!
   TBranch        *b_bG_nscFull5x5_eRightRatio;   //!
   TBranch        *b_bG_scFull5x5_eRightRatio;   //!
   TBranch        *b_bG_nscFull5x5_eBottomRatio;   //!
   TBranch        *b_bG_scFull5x5_eBottomRatio;   //!
   TBranch        *b_bG_nscFull5x5_eLeftRatio;   //!
   TBranch        *b_bG_scFull5x5_eLeftRatio;   //!
   TBranch        *b_bG_nscFull5x5_e2x5MaxRatio;   //!
   TBranch        *b_bG_scFull5x5_e2x5MaxRatio;   //!
   TBranch        *b_bG_nscFull5x5_e2x5TopRatio;   //!
   TBranch        *b_bG_scFull5x5_e2x5TopRatio;   //!
   TBranch        *b_bG_nscFull5x5_e2x5RightRatio;   //!
   TBranch        *b_bG_scFull5x5_e2x5RightRatio;   //!
   TBranch        *b_bG_nscFull5x5_e2x5BottomRatio;   //!
   TBranch        *b_bG_scFull5x5_e2x5BottomRatio;   //!
   TBranch        *b_bG_nscFull5x5_e2x5LeftRatio;   //!
   TBranch        *b_bG_scFull5x5_e2x5LeftRatio;   //!
   TBranch        *b_bG_nscFull5x5_swissCross;   //!
   TBranch        *b_bG_scFull5x5_swissCross;   //!
   TBranch        *b_bG_nscFull5x5_r9;   //!
   TBranch        *b_bG_scFull5x5_r9;   //!
   TBranch        *b_bG_nscFull5x5_sigmaIetaIeta;   //!
   TBranch        *b_bG_scFull5x5_sigmaIetaIeta;   //!
   TBranch        *b_bG_nscFull5x5_sigmaIetaIphi;   //!
   TBranch        *b_bG_scFull5x5_sigmaIetaIphi;   //!
   TBranch        *b_bG_nscFull5x5_sigmaIphiIphi;   //!
   TBranch        *b_bG_scFull5x5_sigmaIphiIphi;   //!
   TBranch        *b_bG_nscNHcalRecHitInDIEta5IPhi5;   //!
   TBranch        *b_bG_scNHcalRecHitInDIEta5IPhi5;   //!
   TBranch        *b_bG_nscEFromHcalRecHitInDIEta5IPhi5;   //!
   TBranch        *b_bG_scEFromHcalRecHitInDIEta5IPhi5;   //!
   TBranch        *b_bG_nscNHcalRecHitInDIEta2IPhi2;   //!
   TBranch        *b_bG_scNHcalRecHitInDIEta2IPhi2;   //!
   TBranch        *b_bG_nscEFromHcalRecHitInDIEta2IPhi2;   //!
   TBranch        *b_bG_scEFromHcalRecHitInDIEta2IPhi2;   //!
   TBranch        *b_bG_nscPFChIso1;   //!
   TBranch        *b_bG_scPFChIso1;   //!
   TBranch        *b_bG_nscPFChIso2;   //!
   TBranch        *b_bG_scPFChIso2;   //!
   TBranch        *b_bG_nscPFChIso3;   //!
   TBranch        *b_bG_scPFChIso3;   //!
   TBranch        *b_bG_nscPFChIso4;   //!
   TBranch        *b_bG_scPFChIso4;   //!
   TBranch        *b_bG_nscPFChIso5;   //!
   TBranch        *b_bG_scPFChIso5;   //!
   TBranch        *b_bG_nscPFPhoIso1;   //!
   TBranch        *b_bG_scPFPhoIso1;   //!
   TBranch        *b_bG_nscPFPhoIso2;   //!
   TBranch        *b_bG_scPFPhoIso2;   //!
   TBranch        *b_bG_nscPFPhoIso3;   //!
   TBranch        *b_bG_scPFPhoIso3;   //!
   TBranch        *b_bG_nscPFPhoIso4;   //!
   TBranch        *b_bG_scPFPhoIso4;   //!
   TBranch        *b_bG_nscPFPhoIso5;   //!
   TBranch        *b_bG_scPFPhoIso5;   //!
   TBranch        *b_bG_nscPFNeuIso1;   //!
   TBranch        *b_bG_scPFNeuIso1;   //!
   TBranch        *b_bG_nscPFNeuIso2;   //!
   TBranch        *b_bG_scPFNeuIso2;   //!
   TBranch        *b_bG_nscPFNeuIso3;   //!
   TBranch        *b_bG_scPFNeuIso3;   //!
   TBranch        *b_bG_nscPFNeuIso4;   //!
   TBranch        *b_bG_scPFNeuIso4;   //!
   TBranch        *b_bG_nscPFNeuIso5;   //!
   TBranch        *b_bG_scPFNeuIso5;   //!
   TBranch        *b_bG_nPrimaryVertex;   //!
   TBranch        *b_bG_primaryVertex_isFake;   //!
   TBranch        *b_bG_primaryVertex_x;   //!
   TBranch        *b_bG_primaryVertex_y;   //!
   TBranch        *b_bG_primaryVertex_z;   //!
   TBranch        *b_bG_primaryVertex_t;   //!
   TBranch        *b_bG_primaryVertex_covXX;   //!
   TBranch        *b_bG_primaryVertex_covXY;   //!
   TBranch        *b_bG_primaryVertex_covXZ;   //!
   TBranch        *b_bG_primaryVertex_covYY;   //!
   TBranch        *b_bG_primaryVertex_covYZ;   //!
   TBranch        *b_bG_primaryVertex_covZZ;   //!
   TBranch        *b_bG_primaryVertex_x_error;   //!
   TBranch        *b_bG_primaryVertex_y_error;   //!
   TBranch        *b_bG_primaryVertex_z_error;   //!
   TBranch        *b_bG_primaryVertex_t_error;   //!
   TBranch        *b_bG_primaryVertex_ntracks;   //!
   TBranch        *b_bG_primaryVertex_ndof;   //!
   TBranch        *b_bG_primaryVertex_chi2;   //!
   TBranch        *b_bG_primaryVertex_normalizedChi2;   //!

   MergedBMMX2018Data(TTree *tree=0);
   virtual ~MergedBMMX2018Data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
 //  virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MergedBMMX2018Data_cxx
MergedBMMX2018Data::MergedBMMX2018Data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/bParkingAnalysis/analysis/CMSSW_10_6_29/src/BsMMGAnalysis/MergeWithBMMNtuples/MergeFiles/results/bph6A.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/bParkingAnalysis/analysis/CMSSW_10_6_29/src/BsMMGAnalysis/MergeWithBMMNtuples/MergeFiles/results/bph6A.root");
      }
      f->GetObject("mergedTree",tree);

   }
   Init(tree);
}

MergedBMMX2018Data::~MergedBMMX2018Data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MergedBMMX2018Data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MergedBMMX2018Data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MergedBMMX2018Data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("b5_run", &b5_run, &b_b5_run);
   fChain->SetBranchAddress("b5_luminosityBlock", &b5_luminosityBlock, &b_b5_luminosityBlock);
   fChain->SetBranchAddress("b5_event", &b5_event, &b_b5_event);
   fChain->SetBranchAddress("b5_nMuonId", &b5_nMuonId, &b_b5_nMuonId);
   fChain->SetBranchAddress("b5_MuonId_chi2LocalPosition", b5_MuonId_chi2LocalPosition, &b_b5_MuonId_chi2LocalPosition);
   fChain->SetBranchAddress("b5_MuonId_glbNormChi2", b5_MuonId_glbNormChi2, &b_b5_MuonId_glbNormChi2);
   fChain->SetBranchAddress("b5_MuonId_glbTrackProbability", b5_MuonId_glbTrackProbability, &b_b5_MuonId_glbTrackProbability);
   fChain->SetBranchAddress("b5_MuonId_match1_dX", b5_MuonId_match1_dX, &b_b5_MuonId_match1_dX);
   fChain->SetBranchAddress("b5_MuonId_match1_dY", b5_MuonId_match1_dY, &b_b5_MuonId_match1_dY);
   fChain->SetBranchAddress("b5_MuonId_match1_pullDxDz", b5_MuonId_match1_pullDxDz, &b_b5_MuonId_match1_pullDxDz);
   fChain->SetBranchAddress("b5_MuonId_match1_pullDyDz", b5_MuonId_match1_pullDyDz, &b_b5_MuonId_match1_pullDyDz);
   fChain->SetBranchAddress("b5_MuonId_match1_pullX", b5_MuonId_match1_pullX, &b_b5_MuonId_match1_pullX);
   fChain->SetBranchAddress("b5_MuonId_match1_pullY", b5_MuonId_match1_pullY, &b_b5_MuonId_match1_pullY);
   fChain->SetBranchAddress("b5_MuonId_match2_dX", b5_MuonId_match2_dX, &b_b5_MuonId_match2_dX);
   fChain->SetBranchAddress("b5_MuonId_match2_dY", b5_MuonId_match2_dY, &b_b5_MuonId_match2_dY);
   fChain->SetBranchAddress("b5_MuonId_match2_pullDxDz", b5_MuonId_match2_pullDxDz, &b_b5_MuonId_match2_pullDxDz);
   fChain->SetBranchAddress("b5_MuonId_match2_pullDyDz", b5_MuonId_match2_pullDyDz, &b_b5_MuonId_match2_pullDyDz);
   fChain->SetBranchAddress("b5_MuonId_match2_pullX", b5_MuonId_match2_pullX, &b_b5_MuonId_match2_pullX);
   fChain->SetBranchAddress("b5_MuonId_match2_pullY", b5_MuonId_match2_pullY, &b_b5_MuonId_match2_pullY);
   fChain->SetBranchAddress("b5_MuonId_newSoftMuonMva", b5_MuonId_newSoftMuonMva, &b_b5_MuonId_newSoftMuonMva);
   fChain->SetBranchAddress("b5_MuonId_trkKink", b5_MuonId_trkKink, &b_b5_MuonId_trkKink);
   fChain->SetBranchAddress("b5_MuonId_trkValidFrac", b5_MuonId_trkValidFrac, &b_b5_MuonId_trkValidFrac);
   fChain->SetBranchAddress("b5_MuonId_highPurity", b5_MuonId_highPurity, &b_b5_MuonId_highPurity);
   fChain->SetBranchAddress("b5_MuonId_nLostHitsInner", b5_MuonId_nLostHitsInner, &b_b5_MuonId_nLostHitsInner);
   fChain->SetBranchAddress("b5_MuonId_nLostHitsOn", b5_MuonId_nLostHitsOn, &b_b5_MuonId_nLostHitsOn);
   fChain->SetBranchAddress("b5_MuonId_nLostHitsOuter", b5_MuonId_nLostHitsOuter, &b_b5_MuonId_nLostHitsOuter);
   fChain->SetBranchAddress("b5_MuonId_nPixels", b5_MuonId_nPixels, &b_b5_MuonId_nPixels);
   fChain->SetBranchAddress("b5_MuonId_nValidHits", b5_MuonId_nValidHits, &b_b5_MuonId_nValidHits);
   fChain->SetBranchAddress("b5_MuonId_trkLayers", b5_MuonId_trkLayers, &b_b5_MuonId_trkLayers);
   fChain->SetBranchAddress("b5_MuonId_trkLostLayersInner", b5_MuonId_trkLostLayersInner, &b_b5_MuonId_trkLostLayersInner);
   fChain->SetBranchAddress("b5_MuonId_trkLostLayersOn", b5_MuonId_trkLostLayersOn, &b_b5_MuonId_trkLostLayersOn);
   fChain->SetBranchAddress("b5_MuonId_trkLostLayersOuter", b5_MuonId_trkLostLayersOuter, &b_b5_MuonId_trkLostLayersOuter);
   fChain->SetBranchAddress("b5_nmm", &b5_nmm, &b_b5_nmm);
   fChain->SetBranchAddress("b5_mm_bdt", b5_mm_bdt, &b_b5_mm_bdt);
   fChain->SetBranchAddress("b5_mm_doca", b5_mm_doca, &b_b5_mm_doca);
   fChain->SetBranchAddress("b5_mm_docatrk", b5_mm_docatrk, &b_b5_mm_docatrk);
   fChain->SetBranchAddress("b5_mm_iso", b5_mm_iso, &b_b5_mm_iso);
   fChain->SetBranchAddress("b5_mm_kal_lxy", b5_mm_kal_lxy, &b_b5_mm_kal_lxy);
   fChain->SetBranchAddress("b5_mm_kal_mass", b5_mm_kal_mass, &b_b5_mm_kal_mass);
   fChain->SetBranchAddress("b5_mm_kal_slxy", b5_mm_kal_slxy, &b_b5_mm_kal_slxy);
   fChain->SetBranchAddress("b5_mm_kal_vtx_prob", b5_mm_kal_vtx_prob, &b_b5_mm_kal_vtx_prob);
   fChain->SetBranchAddress("b5_mm_kin_alpha", b5_mm_kin_alpha, &b_b5_mm_kin_alpha);
   fChain->SetBranchAddress("b5_mm_kin_alphaBS", b5_mm_kin_alphaBS, &b_b5_mm_kin_alphaBS);
   fChain->SetBranchAddress("b5_mm_kin_alphaBSErr", b5_mm_kin_alphaBSErr, &b_b5_mm_kin_alphaBSErr);
   fChain->SetBranchAddress("b5_mm_kin_alphaErr", b5_mm_kin_alphaErr, &b_b5_mm_kin_alphaErr);
   fChain->SetBranchAddress("b5_mm_kin_eta", b5_mm_kin_eta, &b_b5_mm_kin_eta);
   fChain->SetBranchAddress("b5_mm_kin_l3d", b5_mm_kin_l3d, &b_b5_mm_kin_l3d);
   fChain->SetBranchAddress("b5_mm_kin_lxy", b5_mm_kin_lxy, &b_b5_mm_kin_lxy);
   fChain->SetBranchAddress("b5_mm_kin_mass", b5_mm_kin_mass, &b_b5_mm_kin_mass);
   fChain->SetBranchAddress("b5_mm_kin_massErr", b5_mm_kin_massErr, &b_b5_mm_kin_massErr);
   fChain->SetBranchAddress("b5_mm_kin_mu1eta", b5_mm_kin_mu1eta, &b_b5_mm_kin_mu1eta);
   fChain->SetBranchAddress("b5_mm_kin_mu1phi", b5_mm_kin_mu1phi, &b_b5_mm_kin_mu1phi);
   fChain->SetBranchAddress("b5_mm_kin_mu1pt", b5_mm_kin_mu1pt, &b_b5_mm_kin_mu1pt);
   fChain->SetBranchAddress("b5_mm_kin_mu2eta", b5_mm_kin_mu2eta, &b_b5_mm_kin_mu2eta);
   fChain->SetBranchAddress("b5_mm_kin_mu2phi", b5_mm_kin_mu2phi, &b_b5_mm_kin_mu2phi);
   fChain->SetBranchAddress("b5_mm_kin_mu2pt", b5_mm_kin_mu2pt, &b_b5_mm_kin_mu2pt);
   fChain->SetBranchAddress("b5_mm_kin_phi", b5_mm_kin_phi, &b_b5_mm_kin_phi);
   fChain->SetBranchAddress("b5_mm_kin_pt", b5_mm_kin_pt, &b_b5_mm_kin_pt);
   fChain->SetBranchAddress("b5_mm_kin_pv2ip", b5_mm_kin_pv2ip, &b_b5_mm_kin_pv2ip);
   fChain->SetBranchAddress("b5_mm_kin_pv2ipErr", b5_mm_kin_pv2ipErr, &b_b5_mm_kin_pv2ipErr);
   fChain->SetBranchAddress("b5_mm_kin_pv2lip", b5_mm_kin_pv2lip, &b_b5_mm_kin_pv2lip);
   fChain->SetBranchAddress("b5_mm_kin_pv2lipErr", b5_mm_kin_pv2lipErr, &b_b5_mm_kin_pv2lipErr);
   fChain->SetBranchAddress("b5_mm_kin_pv2lipSig", b5_mm_kin_pv2lipSig, &b_b5_mm_kin_pv2lipSig);
   fChain->SetBranchAddress("b5_mm_kin_pv_z", b5_mm_kin_pv_z, &b_b5_mm_kin_pv_z);
   fChain->SetBranchAddress("b5_mm_kin_pv_zErr", b5_mm_kin_pv_zErr, &b_b5_mm_kin_pv_zErr);
   fChain->SetBranchAddress("b5_mm_kin_pvip", b5_mm_kin_pvip, &b_b5_mm_kin_pvip);
   fChain->SetBranchAddress("b5_mm_kin_pvipErr", b5_mm_kin_pvipErr, &b_b5_mm_kin_pvipErr);
   fChain->SetBranchAddress("b5_mm_kin_pvlip", b5_mm_kin_pvlip, &b_b5_mm_kin_pvlip);
   fChain->SetBranchAddress("b5_mm_kin_pvlipErr", b5_mm_kin_pvlipErr, &b_b5_mm_kin_pvlipErr);
   fChain->SetBranchAddress("b5_mm_kin_pvlipSig", b5_mm_kin_pvlipSig, &b_b5_mm_kin_pvlipSig);
   fChain->SetBranchAddress("b5_mm_kin_sl3d", b5_mm_kin_sl3d, &b_b5_mm_kin_sl3d);
   fChain->SetBranchAddress("b5_mm_kin_slxy", b5_mm_kin_slxy, &b_b5_mm_kin_slxy);
   fChain->SetBranchAddress("b5_mm_kin_spv2ip", b5_mm_kin_spv2ip, &b_b5_mm_kin_spv2ip);
   fChain->SetBranchAddress("b5_mm_kin_spvip", b5_mm_kin_spvip, &b_b5_mm_kin_spvip);
   fChain->SetBranchAddress("b5_mm_kin_tau", b5_mm_kin_tau, &b_b5_mm_kin_tau);
   fChain->SetBranchAddress("b5_mm_kin_taue", b5_mm_kin_taue, &b_b5_mm_kin_taue);
   fChain->SetBranchAddress("b5_mm_kin_tauxy", b5_mm_kin_tauxy, &b_b5_mm_kin_tauxy);
   fChain->SetBranchAddress("b5_mm_kin_tauxye", b5_mm_kin_tauxye, &b_b5_mm_kin_tauxye);
   fChain->SetBranchAddress("b5_mm_kin_vtx_chi2dof", b5_mm_kin_vtx_chi2dof, &b_b5_mm_kin_vtx_chi2dof);
   fChain->SetBranchAddress("b5_mm_kin_vtx_prob", b5_mm_kin_vtx_prob, &b_b5_mm_kin_vtx_prob);
   fChain->SetBranchAddress("b5_mm_kin_vtx_x", b5_mm_kin_vtx_x, &b_b5_mm_kin_vtx_x);
   fChain->SetBranchAddress("b5_mm_kin_vtx_xErr", b5_mm_kin_vtx_xErr, &b_b5_mm_kin_vtx_xErr);
   fChain->SetBranchAddress("b5_mm_kin_vtx_y", b5_mm_kin_vtx_y, &b_b5_mm_kin_vtx_y);
   fChain->SetBranchAddress("b5_mm_kin_vtx_yErr", b5_mm_kin_vtx_yErr, &b_b5_mm_kin_vtx_yErr);
   fChain->SetBranchAddress("b5_mm_kin_vtx_z", b5_mm_kin_vtx_z, &b_b5_mm_kin_vtx_z);
   fChain->SetBranchAddress("b5_mm_kin_vtx_zErr", b5_mm_kin_vtx_zErr, &b_b5_mm_kin_vtx_zErr);
   fChain->SetBranchAddress("b5_mm_kinpc_alpha", b5_mm_kinpc_alpha, &b_b5_mm_kinpc_alpha);
   fChain->SetBranchAddress("b5_mm_kinpc_alphaBS", b5_mm_kinpc_alphaBS, &b_b5_mm_kinpc_alphaBS);
   fChain->SetBranchAddress("b5_mm_kinpc_alphaBSErr", b5_mm_kinpc_alphaBSErr, &b_b5_mm_kinpc_alphaBSErr);
   fChain->SetBranchAddress("b5_mm_kinpc_alphaErr", b5_mm_kinpc_alphaErr, &b_b5_mm_kinpc_alphaErr);
   fChain->SetBranchAddress("b5_mm_kinpc_eta", b5_mm_kinpc_eta, &b_b5_mm_kinpc_eta);
   fChain->SetBranchAddress("b5_mm_kinpc_l3d", b5_mm_kinpc_l3d, &b_b5_mm_kinpc_l3d);
   fChain->SetBranchAddress("b5_mm_kinpc_lxy", b5_mm_kinpc_lxy, &b_b5_mm_kinpc_lxy);
   fChain->SetBranchAddress("b5_mm_kinpc_mass", b5_mm_kinpc_mass, &b_b5_mm_kinpc_mass);
   fChain->SetBranchAddress("b5_mm_kinpc_massErr", b5_mm_kinpc_massErr, &b_b5_mm_kinpc_massErr);
   fChain->SetBranchAddress("b5_mm_kinpc_phi", b5_mm_kinpc_phi, &b_b5_mm_kinpc_phi);
   fChain->SetBranchAddress("b5_mm_kinpc_pt", b5_mm_kinpc_pt, &b_b5_mm_kinpc_pt);
   fChain->SetBranchAddress("b5_mm_kinpc_pv2ip", b5_mm_kinpc_pv2ip, &b_b5_mm_kinpc_pv2ip);
   fChain->SetBranchAddress("b5_mm_kinpc_pv2ipErr", b5_mm_kinpc_pv2ipErr, &b_b5_mm_kinpc_pv2ipErr);
   fChain->SetBranchAddress("b5_mm_kinpc_pv2lip", b5_mm_kinpc_pv2lip, &b_b5_mm_kinpc_pv2lip);
   fChain->SetBranchAddress("b5_mm_kinpc_pv2lipErr", b5_mm_kinpc_pv2lipErr, &b_b5_mm_kinpc_pv2lipErr);
   fChain->SetBranchAddress("b5_mm_kinpc_pv2lipSig", b5_mm_kinpc_pv2lipSig, &b_b5_mm_kinpc_pv2lipSig);
   fChain->SetBranchAddress("b5_mm_kinpc_pv_z", b5_mm_kinpc_pv_z, &b_b5_mm_kinpc_pv_z);
   fChain->SetBranchAddress("b5_mm_kinpc_pv_zErr", b5_mm_kinpc_pv_zErr, &b_b5_mm_kinpc_pv_zErr);
   fChain->SetBranchAddress("b5_mm_kinpc_pvip", b5_mm_kinpc_pvip, &b_b5_mm_kinpc_pvip);
   fChain->SetBranchAddress("b5_mm_kinpc_pvipErr", b5_mm_kinpc_pvipErr, &b_b5_mm_kinpc_pvipErr);
   fChain->SetBranchAddress("b5_mm_kinpc_pvlip", b5_mm_kinpc_pvlip, &b_b5_mm_kinpc_pvlip);
   fChain->SetBranchAddress("b5_mm_kinpc_pvlipErr", b5_mm_kinpc_pvlipErr, &b_b5_mm_kinpc_pvlipErr);
   fChain->SetBranchAddress("b5_mm_kinpc_pvlipSig", b5_mm_kinpc_pvlipSig, &b_b5_mm_kinpc_pvlipSig);
   fChain->SetBranchAddress("b5_mm_kinpc_sl3d", b5_mm_kinpc_sl3d, &b_b5_mm_kinpc_sl3d);
   fChain->SetBranchAddress("b5_mm_kinpc_slxy", b5_mm_kinpc_slxy, &b_b5_mm_kinpc_slxy);
   fChain->SetBranchAddress("b5_mm_kinpc_spv2ip", b5_mm_kinpc_spv2ip, &b_b5_mm_kinpc_spv2ip);
   fChain->SetBranchAddress("b5_mm_kinpc_spvip", b5_mm_kinpc_spvip, &b_b5_mm_kinpc_spvip);
   fChain->SetBranchAddress("b5_mm_kinpc_tau", b5_mm_kinpc_tau, &b_b5_mm_kinpc_tau);
   fChain->SetBranchAddress("b5_mm_kinpc_taue", b5_mm_kinpc_taue, &b_b5_mm_kinpc_taue);
   fChain->SetBranchAddress("b5_mm_kinpc_tauxy", b5_mm_kinpc_tauxy, &b_b5_mm_kinpc_tauxy);
   fChain->SetBranchAddress("b5_mm_kinpc_tauxye", b5_mm_kinpc_tauxye, &b_b5_mm_kinpc_tauxye);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_chi2dof", b5_mm_kinpc_vtx_chi2dof, &b_b5_mm_kinpc_vtx_chi2dof);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_prob", b5_mm_kinpc_vtx_prob, &b_b5_mm_kinpc_vtx_prob);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_x", b5_mm_kinpc_vtx_x, &b_b5_mm_kinpc_vtx_x);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_xErr", b5_mm_kinpc_vtx_xErr, &b_b5_mm_kinpc_vtx_xErr);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_y", b5_mm_kinpc_vtx_y, &b_b5_mm_kinpc_vtx_y);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_yErr", b5_mm_kinpc_vtx_yErr, &b_b5_mm_kinpc_vtx_yErr);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_z", b5_mm_kinpc_vtx_z, &b_b5_mm_kinpc_vtx_z);
   fChain->SetBranchAddress("b5_mm_kinpc_vtx_zErr", b5_mm_kinpc_vtx_zErr, &b_b5_mm_kinpc_vtx_zErr);
   fChain->SetBranchAddress("b5_mm_m1iso", b5_mm_m1iso, &b_b5_mm_m1iso);
   fChain->SetBranchAddress("b5_mm_m2iso", b5_mm_m2iso, &b_b5_mm_m2iso);
   fChain->SetBranchAddress("b5_mm_mass", b5_mm_mass, &b_b5_mm_mass);
   fChain->SetBranchAddress("b5_mm_mu1_eta", b5_mm_mu1_eta, &b_b5_mm_mu1_eta);
   fChain->SetBranchAddress("b5_mm_mu1_phi", b5_mm_mu1_phi, &b_b5_mm_mu1_phi);
   fChain->SetBranchAddress("b5_mm_mu1_pt", b5_mm_mu1_pt, &b_b5_mm_mu1_pt);
   fChain->SetBranchAddress("b5_mm_mu2_eta", b5_mm_mu2_eta, &b_b5_mm_mu2_eta);
   fChain->SetBranchAddress("b5_mm_mu2_phi", b5_mm_mu2_phi, &b_b5_mm_mu2_phi);
   fChain->SetBranchAddress("b5_mm_mu2_pt", b5_mm_mu2_pt, &b_b5_mm_mu2_pt);
   fChain->SetBranchAddress("b5_mm_mva", b5_mm_mva, &b_b5_mm_mva);
   fChain->SetBranchAddress("b5_mm_otherVtxMaxProb", b5_mm_otherVtxMaxProb, &b_b5_mm_otherVtxMaxProb);
   fChain->SetBranchAddress("b5_mm_otherVtxMaxProb1", b5_mm_otherVtxMaxProb1, &b_b5_mm_otherVtxMaxProb1);
   fChain->SetBranchAddress("b5_mm_otherVtxMaxProb2", b5_mm_otherVtxMaxProb2, &b_b5_mm_otherVtxMaxProb2);
   fChain->SetBranchAddress("b5_mm_closetrk", b5_mm_closetrk, &b_b5_mm_closetrk);
   fChain->SetBranchAddress("b5_mm_closetrks1", b5_mm_closetrks1, &b_b5_mm_closetrks1);
   fChain->SetBranchAddress("b5_mm_closetrks2", b5_mm_closetrks2, &b_b5_mm_closetrks2);
   fChain->SetBranchAddress("b5_mm_closetrks3", b5_mm_closetrks3, &b_b5_mm_closetrks3);
   fChain->SetBranchAddress("b5_mm_kal_valid", b5_mm_kal_valid, &b_b5_mm_kal_valid);
   fChain->SetBranchAddress("b5_mm_kin_valid", b5_mm_kin_valid, &b_b5_mm_kin_valid);
   fChain->SetBranchAddress("b5_mm_kinpc_valid", b5_mm_kinpc_valid, &b_b5_mm_kinpc_valid);
   fChain->SetBranchAddress("b5_mm_mu1_index", b5_mm_mu1_index, &b_b5_mm_mu1_index);
   fChain->SetBranchAddress("b5_mm_mu1_pdgId", b5_mm_mu1_pdgId, &b_b5_mm_mu1_pdgId);
   fChain->SetBranchAddress("b5_mm_mu2_index", b5_mm_mu2_index, &b_b5_mm_mu2_index);
   fChain->SetBranchAddress("b5_mm_mu2_pdgId", b5_mm_mu2_pdgId, &b_b5_mm_mu2_pdgId);
   fChain->SetBranchAddress("b5_mm_nBMTrks", b5_mm_nBMTrks, &b_b5_mm_nBMTrks);
   fChain->SetBranchAddress("b5_mm_nDisTrks", b5_mm_nDisTrks, &b_b5_mm_nDisTrks);
   fChain->SetBranchAddress("b5_mm_nTrks", b5_mm_nTrks, &b_b5_mm_nTrks);
   fChain->SetBranchAddress("b5_nd0", &b5_nd0, &b_b5_nd0);
   fChain->SetBranchAddress("b5_d0_doca", b5_d0_doca, &b_b5_d0_doca);
   fChain->SetBranchAddress("b5_d0_kaon_eta", b5_d0_kaon_eta, &b_b5_d0_kaon_eta);
   fChain->SetBranchAddress("b5_d0_kaon_phi", b5_d0_kaon_phi, &b_b5_d0_kaon_phi);
   fChain->SetBranchAddress("b5_d0_kaon_pt", b5_d0_kaon_pt, &b_b5_d0_kaon_pt);
   fChain->SetBranchAddress("b5_d0_kaon_sip", b5_d0_kaon_sip, &b_b5_d0_kaon_sip);
   fChain->SetBranchAddress("b5_d0_kin_cosAlphaXY", b5_d0_kin_cosAlphaXY, &b_b5_d0_kin_cosAlphaXY);
   fChain->SetBranchAddress("b5_d0_kin_eta", b5_d0_kin_eta, &b_b5_d0_kin_eta);
   fChain->SetBranchAddress("b5_d0_kin_lxy", b5_d0_kin_lxy, &b_b5_d0_kin_lxy);
   fChain->SetBranchAddress("b5_d0_kin_mass", b5_d0_kin_mass, &b_b5_d0_kin_mass);
   fChain->SetBranchAddress("b5_d0_kin_massErr", b5_d0_kin_massErr, &b_b5_d0_kin_massErr);
   fChain->SetBranchAddress("b5_d0_kin_phi", b5_d0_kin_phi, &b_b5_d0_kin_phi);
   fChain->SetBranchAddress("b5_d0_kin_pt", b5_d0_kin_pt, &b_b5_d0_kin_pt);
   fChain->SetBranchAddress("b5_d0_kin_sipBS", b5_d0_kin_sipBS, &b_b5_d0_kin_sipBS);
   fChain->SetBranchAddress("b5_d0_kin_sipPV", b5_d0_kin_sipPV, &b_b5_d0_kin_sipPV);
   fChain->SetBranchAddress("b5_d0_kin_slxy", b5_d0_kin_slxy, &b_b5_d0_kin_slxy);
   fChain->SetBranchAddress("b5_d0_kin_vtx_chi2dof", b5_d0_kin_vtx_chi2dof, &b_b5_d0_kin_vtx_chi2dof);
   fChain->SetBranchAddress("b5_d0_kin_vtx_prob", b5_d0_kin_vtx_prob, &b_b5_d0_kin_vtx_prob);
   fChain->SetBranchAddress("b5_d0_mass", b5_d0_mass, &b_b5_d0_mass);
   fChain->SetBranchAddress("b5_d0_pion_eta", b5_d0_pion_eta, &b_b5_d0_pion_eta);
   fChain->SetBranchAddress("b5_d0_pion_phi", b5_d0_pion_phi, &b_b5_d0_pion_phi);
   fChain->SetBranchAddress("b5_d0_pion_pt", b5_d0_pion_pt, &b_b5_d0_pion_pt);
   fChain->SetBranchAddress("b5_d0_pion_sip", b5_d0_pion_sip, &b_b5_d0_pion_sip);
   fChain->SetBranchAddress("b5_d0_kaon_mu_index", b5_d0_kaon_mu_index, &b_b5_d0_kaon_mu_index);
   fChain->SetBranchAddress("b5_d0_kin_valid", b5_d0_kin_valid, &b_b5_d0_kin_valid);
   fChain->SetBranchAddress("b5_d0_pion_mu_index", b5_d0_pion_mu_index, &b_b5_d0_pion_mu_index);
   fChain->SetBranchAddress("b5_nks", &b5_nks, &b_b5_nks);
   fChain->SetBranchAddress("b5_ks_doca", b5_ks_doca, &b_b5_ks_doca);
   fChain->SetBranchAddress("b5_ks_kin_cosAlphaXY", b5_ks_kin_cosAlphaXY, &b_b5_ks_kin_cosAlphaXY);
   fChain->SetBranchAddress("b5_ks_kin_eta", b5_ks_kin_eta, &b_b5_ks_kin_eta);
   fChain->SetBranchAddress("b5_ks_kin_lxy", b5_ks_kin_lxy, &b_b5_ks_kin_lxy);
   fChain->SetBranchAddress("b5_ks_kin_mass", b5_ks_kin_mass, &b_b5_ks_kin_mass);
   fChain->SetBranchAddress("b5_ks_kin_massErr", b5_ks_kin_massErr, &b_b5_ks_kin_massErr);
   fChain->SetBranchAddress("b5_ks_kin_phi", b5_ks_kin_phi, &b_b5_ks_kin_phi);
   fChain->SetBranchAddress("b5_ks_kin_pt", b5_ks_kin_pt, &b_b5_ks_kin_pt);
   fChain->SetBranchAddress("b5_ks_kin_sipBS", b5_ks_kin_sipBS, &b_b5_ks_kin_sipBS);
   fChain->SetBranchAddress("b5_ks_kin_sipPV", b5_ks_kin_sipPV, &b_b5_ks_kin_sipPV);
   fChain->SetBranchAddress("b5_ks_kin_slxy", b5_ks_kin_slxy, &b_b5_ks_kin_slxy);
   fChain->SetBranchAddress("b5_ks_kin_vtx_chi2dof", b5_ks_kin_vtx_chi2dof, &b_b5_ks_kin_vtx_chi2dof);
   fChain->SetBranchAddress("b5_ks_kin_vtx_prob", b5_ks_kin_vtx_prob, &b_b5_ks_kin_vtx_prob);
   fChain->SetBranchAddress("b5_ks_mass", b5_ks_mass, &b_b5_ks_mass);
   fChain->SetBranchAddress("b5_ks_trk1_eta", b5_ks_trk1_eta, &b_b5_ks_trk1_eta);
   fChain->SetBranchAddress("b5_ks_trk1_phi", b5_ks_trk1_phi, &b_b5_ks_trk1_phi);
   fChain->SetBranchAddress("b5_ks_trk1_pt", b5_ks_trk1_pt, &b_b5_ks_trk1_pt);
   fChain->SetBranchAddress("b5_ks_trk1_sip", b5_ks_trk1_sip, &b_b5_ks_trk1_sip);
   fChain->SetBranchAddress("b5_ks_trk2_eta", b5_ks_trk2_eta, &b_b5_ks_trk2_eta);
   fChain->SetBranchAddress("b5_ks_trk2_phi", b5_ks_trk2_phi, &b_b5_ks_trk2_phi);
   fChain->SetBranchAddress("b5_ks_trk2_pt", b5_ks_trk2_pt, &b_b5_ks_trk2_pt);
   fChain->SetBranchAddress("b5_ks_trk2_sip", b5_ks_trk2_sip, &b_b5_ks_trk2_sip);
   fChain->SetBranchAddress("b5_ks_kin_valid", b5_ks_kin_valid, &b_b5_ks_kin_valid);
   fChain->SetBranchAddress("b5_ks_trk1_mu_index", b5_ks_trk1_mu_index, &b_b5_ks_trk1_mu_index);
   fChain->SetBranchAddress("b5_ks_trk2_mu_index", b5_ks_trk2_mu_index, &b_b5_ks_trk2_mu_index);
   fChain->SetBranchAddress("b5_nlambda", &b5_nlambda, &b_b5_nlambda);
   fChain->SetBranchAddress("b5_lambda_doca", b5_lambda_doca, &b_b5_lambda_doca);
   fChain->SetBranchAddress("b5_lambda_kin_cosAlphaXY", b5_lambda_kin_cosAlphaXY, &b_b5_lambda_kin_cosAlphaXY);
   fChain->SetBranchAddress("b5_lambda_kin_eta", b5_lambda_kin_eta, &b_b5_lambda_kin_eta);
   fChain->SetBranchAddress("b5_lambda_kin_lxy", b5_lambda_kin_lxy, &b_b5_lambda_kin_lxy);
   fChain->SetBranchAddress("b5_lambda_kin_mass", b5_lambda_kin_mass, &b_b5_lambda_kin_mass);
   fChain->SetBranchAddress("b5_lambda_kin_massErr", b5_lambda_kin_massErr, &b_b5_lambda_kin_massErr);
   fChain->SetBranchAddress("b5_lambda_kin_phi", b5_lambda_kin_phi, &b_b5_lambda_kin_phi);
   fChain->SetBranchAddress("b5_lambda_kin_pt", b5_lambda_kin_pt, &b_b5_lambda_kin_pt);
   fChain->SetBranchAddress("b5_lambda_kin_sipBS", b5_lambda_kin_sipBS, &b_b5_lambda_kin_sipBS);
   fChain->SetBranchAddress("b5_lambda_kin_sipPV", b5_lambda_kin_sipPV, &b_b5_lambda_kin_sipPV);
   fChain->SetBranchAddress("b5_lambda_kin_slxy", b5_lambda_kin_slxy, &b_b5_lambda_kin_slxy);
   fChain->SetBranchAddress("b5_lambda_kin_vtx_chi2dof", b5_lambda_kin_vtx_chi2dof, &b_b5_lambda_kin_vtx_chi2dof);
   fChain->SetBranchAddress("b5_lambda_kin_vtx_prob", b5_lambda_kin_vtx_prob, &b_b5_lambda_kin_vtx_prob);
   fChain->SetBranchAddress("b5_lambda_mass", b5_lambda_mass, &b_b5_lambda_mass);
   fChain->SetBranchAddress("b5_lambda_pion_eta", b5_lambda_pion_eta, &b_b5_lambda_pion_eta);
   fChain->SetBranchAddress("b5_lambda_pion_phi", b5_lambda_pion_phi, &b_b5_lambda_pion_phi);
   fChain->SetBranchAddress("b5_lambda_pion_pt", b5_lambda_pion_pt, &b_b5_lambda_pion_pt);
   fChain->SetBranchAddress("b5_lambda_pion_sip", b5_lambda_pion_sip, &b_b5_lambda_pion_sip);
   fChain->SetBranchAddress("b5_lambda_proton_eta", b5_lambda_proton_eta, &b_b5_lambda_proton_eta);
   fChain->SetBranchAddress("b5_lambda_proton_phi", b5_lambda_proton_phi, &b_b5_lambda_proton_phi);
   fChain->SetBranchAddress("b5_lambda_proton_pt", b5_lambda_proton_pt, &b_b5_lambda_proton_pt);
   fChain->SetBranchAddress("b5_lambda_proton_sip", b5_lambda_proton_sip, &b_b5_lambda_proton_sip);
   fChain->SetBranchAddress("b5_lambda_kin_valid", b5_lambda_kin_valid, &b_b5_lambda_kin_valid);
   fChain->SetBranchAddress("b5_lambda_pion_mu_index", b5_lambda_pion_mu_index, &b_b5_lambda_pion_mu_index);
   fChain->SetBranchAddress("b5_lambda_proton_mu_index", b5_lambda_proton_mu_index, &b_b5_lambda_proton_mu_index);
   fChain->SetBranchAddress("b5_nphi", &b5_nphi, &b_b5_nphi);
   fChain->SetBranchAddress("b5_phi_doca", b5_phi_doca, &b_b5_phi_doca);
   fChain->SetBranchAddress("b5_phi_ds_cosAlphaXY", b5_phi_ds_cosAlphaXY, &b_b5_phi_ds_cosAlphaXY);
   fChain->SetBranchAddress("b5_phi_ds_eta", b5_phi_ds_eta, &b_b5_phi_ds_eta);
   fChain->SetBranchAddress("b5_phi_ds_lxy", b5_phi_ds_lxy, &b_b5_phi_ds_lxy);
   fChain->SetBranchAddress("b5_phi_ds_mass", b5_phi_ds_mass, &b_b5_phi_ds_mass);
   fChain->SetBranchAddress("b5_phi_ds_massErr", b5_phi_ds_massErr, &b_b5_phi_ds_massErr);
   fChain->SetBranchAddress("b5_phi_ds_phi", b5_phi_ds_phi, &b_b5_phi_ds_phi);
   fChain->SetBranchAddress("b5_phi_ds_pion_eta", b5_phi_ds_pion_eta, &b_b5_phi_ds_pion_eta);
   fChain->SetBranchAddress("b5_phi_ds_pion_mu_index", b5_phi_ds_pion_mu_index, &b_b5_phi_ds_pion_mu_index);
   fChain->SetBranchAddress("b5_phi_ds_pion_phi", b5_phi_ds_pion_phi, &b_b5_phi_ds_pion_phi);
   fChain->SetBranchAddress("b5_phi_ds_pion_pt", b5_phi_ds_pion_pt, &b_b5_phi_ds_pion_pt);
   fChain->SetBranchAddress("b5_phi_ds_pt", b5_phi_ds_pt, &b_b5_phi_ds_pt);
   fChain->SetBranchAddress("b5_phi_ds_sipBS", b5_phi_ds_sipBS, &b_b5_phi_ds_sipBS);
   fChain->SetBranchAddress("b5_phi_ds_sipPV", b5_phi_ds_sipPV, &b_b5_phi_ds_sipPV);
   fChain->SetBranchAddress("b5_phi_ds_slxy", b5_phi_ds_slxy, &b_b5_phi_ds_slxy);
   fChain->SetBranchAddress("b5_phi_ds_vtx_chi2dof", b5_phi_ds_vtx_chi2dof, &b_b5_phi_ds_vtx_chi2dof);
   fChain->SetBranchAddress("b5_phi_ds_vtx_prob", b5_phi_ds_vtx_prob, &b_b5_phi_ds_vtx_prob);
   fChain->SetBranchAddress("b5_phi_kin_cosAlphaXY", b5_phi_kin_cosAlphaXY, &b_b5_phi_kin_cosAlphaXY);
   fChain->SetBranchAddress("b5_phi_kin_eta", b5_phi_kin_eta, &b_b5_phi_kin_eta);
   fChain->SetBranchAddress("b5_phi_kin_lxy", b5_phi_kin_lxy, &b_b5_phi_kin_lxy);
   fChain->SetBranchAddress("b5_phi_kin_mass", b5_phi_kin_mass, &b_b5_phi_kin_mass);
   fChain->SetBranchAddress("b5_phi_kin_massErr", b5_phi_kin_massErr, &b_b5_phi_kin_massErr);
   fChain->SetBranchAddress("b5_phi_kin_phi", b5_phi_kin_phi, &b_b5_phi_kin_phi);
   fChain->SetBranchAddress("b5_phi_kin_pt", b5_phi_kin_pt, &b_b5_phi_kin_pt);
   fChain->SetBranchAddress("b5_phi_kin_sipBS", b5_phi_kin_sipBS, &b_b5_phi_kin_sipBS);
   fChain->SetBranchAddress("b5_phi_kin_sipPV", b5_phi_kin_sipPV, &b_b5_phi_kin_sipPV);
   fChain->SetBranchAddress("b5_phi_kin_slxy", b5_phi_kin_slxy, &b_b5_phi_kin_slxy);
   fChain->SetBranchAddress("b5_phi_kin_vtx_chi2dof", b5_phi_kin_vtx_chi2dof, &b_b5_phi_kin_vtx_chi2dof);
   fChain->SetBranchAddress("b5_phi_kin_vtx_prob", b5_phi_kin_vtx_prob, &b_b5_phi_kin_vtx_prob);
   fChain->SetBranchAddress("b5_phi_mass", b5_phi_mass, &b_b5_phi_mass);
   fChain->SetBranchAddress("b5_phi_trk1_eta", b5_phi_trk1_eta, &b_b5_phi_trk1_eta);
   fChain->SetBranchAddress("b5_phi_trk1_phi", b5_phi_trk1_phi, &b_b5_phi_trk1_phi);
   fChain->SetBranchAddress("b5_phi_trk1_pt", b5_phi_trk1_pt, &b_b5_phi_trk1_pt);
   fChain->SetBranchAddress("b5_phi_trk1_sip", b5_phi_trk1_sip, &b_b5_phi_trk1_sip);
   fChain->SetBranchAddress("b5_phi_trk2_eta", b5_phi_trk2_eta, &b_b5_phi_trk2_eta);
   fChain->SetBranchAddress("b5_phi_trk2_phi", b5_phi_trk2_phi, &b_b5_phi_trk2_phi);
   fChain->SetBranchAddress("b5_phi_trk2_pt", b5_phi_trk2_pt, &b_b5_phi_trk2_pt);
   fChain->SetBranchAddress("b5_phi_trk2_sip", b5_phi_trk2_sip, &b_b5_phi_trk2_sip);
   fChain->SetBranchAddress("b5_phi_kin_valid", b5_phi_kin_valid, &b_b5_phi_kin_valid);
   fChain->SetBranchAddress("b5_phi_trk1_mu_index", b5_phi_trk1_mu_index, &b_b5_phi_trk1_mu_index);
   fChain->SetBranchAddress("b5_phi_trk2_mu_index", b5_phi_trk2_mu_index, &b_b5_phi_trk2_mu_index);
   fChain->SetBranchAddress("b5_CaloMET_phi", &b5_CaloMET_phi, &b_b5_CaloMET_phi);
   fChain->SetBranchAddress("b5_CaloMET_pt", &b5_CaloMET_pt, &b_b5_CaloMET_pt);
   fChain->SetBranchAddress("b5_CaloMET_sumEt", &b5_CaloMET_sumEt, &b_b5_CaloMET_sumEt);
   fChain->SetBranchAddress("b5_ChsMET_phi", &b5_ChsMET_phi, &b_b5_ChsMET_phi);
   fChain->SetBranchAddress("b5_ChsMET_pt", &b5_ChsMET_pt, &b_b5_ChsMET_pt);
   fChain->SetBranchAddress("b5_ChsMET_sumEt", &b5_ChsMET_sumEt, &b_b5_ChsMET_sumEt);
   fChain->SetBranchAddress("b5_nCorrT1METJet", &b5_nCorrT1METJet, &b_b5_nCorrT1METJet);
   fChain->SetBranchAddress("b5_CorrT1METJet_area", b5_CorrT1METJet_area, &b_b5_CorrT1METJet_area);
   fChain->SetBranchAddress("b5_CorrT1METJet_eta", b5_CorrT1METJet_eta, &b_b5_CorrT1METJet_eta);
   fChain->SetBranchAddress("b5_CorrT1METJet_muonSubtrFactor", b5_CorrT1METJet_muonSubtrFactor, &b_b5_CorrT1METJet_muonSubtrFactor);
   fChain->SetBranchAddress("b5_CorrT1METJet_phi", b5_CorrT1METJet_phi, &b_b5_CorrT1METJet_phi);
   fChain->SetBranchAddress("b5_CorrT1METJet_rawPt", b5_CorrT1METJet_rawPt, &b_b5_CorrT1METJet_rawPt);
   fChain->SetBranchAddress("b5_DeepMETResolutionTune_phi", &b5_DeepMETResolutionTune_phi, &b_b5_DeepMETResolutionTune_phi);
   fChain->SetBranchAddress("b5_DeepMETResolutionTune_pt", &b5_DeepMETResolutionTune_pt, &b_b5_DeepMETResolutionTune_pt);
   fChain->SetBranchAddress("b5_DeepMETResponseTune_phi", &b5_DeepMETResponseTune_phi, &b_b5_DeepMETResponseTune_phi);
   fChain->SetBranchAddress("b5_DeepMETResponseTune_pt", &b5_DeepMETResponseTune_pt, &b_b5_DeepMETResponseTune_pt);
   fChain->SetBranchAddress("b5_nElectron", &b5_nElectron, &b_b5_nElectron);
   fChain->SetBranchAddress("b5_Electron_deltaEtaSC", b5_Electron_deltaEtaSC, &b_b5_Electron_deltaEtaSC);
   fChain->SetBranchAddress("b5_Electron_dr03EcalRecHitSumEt", b5_Electron_dr03EcalRecHitSumEt, &b_b5_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("b5_Electron_dr03HcalDepth1TowerSumEt", b5_Electron_dr03HcalDepth1TowerSumEt, &b_b5_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("b5_Electron_dr03TkSumPt", b5_Electron_dr03TkSumPt, &b_b5_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("b5_Electron_dr03TkSumPtHEEP", b5_Electron_dr03TkSumPtHEEP, &b_b5_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("b5_Electron_dxy", b5_Electron_dxy, &b_b5_Electron_dxy);
   fChain->SetBranchAddress("b5_Electron_dxyErr", b5_Electron_dxyErr, &b_b5_Electron_dxyErr);
   fChain->SetBranchAddress("b5_Electron_dz", b5_Electron_dz, &b_b5_Electron_dz);
   fChain->SetBranchAddress("b5_Electron_dzErr", b5_Electron_dzErr, &b_b5_Electron_dzErr);
   fChain->SetBranchAddress("b5_Electron_eCorr", b5_Electron_eCorr, &b_b5_Electron_eCorr);
   fChain->SetBranchAddress("b5_Electron_eInvMinusPInv", b5_Electron_eInvMinusPInv, &b_b5_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("b5_Electron_energyErr", b5_Electron_energyErr, &b_b5_Electron_energyErr);
   fChain->SetBranchAddress("b5_Electron_eta", b5_Electron_eta, &b_b5_Electron_eta);
   fChain->SetBranchAddress("b5_Electron_hoe", b5_Electron_hoe, &b_b5_Electron_hoe);
   fChain->SetBranchAddress("b5_Electron_ip3d", b5_Electron_ip3d, &b_b5_Electron_ip3d);
   fChain->SetBranchAddress("b5_Electron_jetPtRelv2", b5_Electron_jetPtRelv2, &b_b5_Electron_jetPtRelv2);
   fChain->SetBranchAddress("b5_Electron_jetRelIso", b5_Electron_jetRelIso, &b_b5_Electron_jetRelIso);
   fChain->SetBranchAddress("b5_Electron_mass", b5_Electron_mass, &b_b5_Electron_mass);
   fChain->SetBranchAddress("b5_Electron_miniPFRelIso_all", b5_Electron_miniPFRelIso_all, &b_b5_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("b5_Electron_miniPFRelIso_chg", b5_Electron_miniPFRelIso_chg, &b_b5_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1Iso", b5_Electron_mvaFall17V1Iso, &b_b5_Electron_mvaFall17V1Iso);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1noIso", b5_Electron_mvaFall17V1noIso, &b_b5_Electron_mvaFall17V1noIso);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2Iso", b5_Electron_mvaFall17V2Iso, &b_b5_Electron_mvaFall17V2Iso);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2noIso", b5_Electron_mvaFall17V2noIso, &b_b5_Electron_mvaFall17V2noIso);
   fChain->SetBranchAddress("b5_Electron_pfRelIso03_all", b5_Electron_pfRelIso03_all, &b_b5_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("b5_Electron_pfRelIso03_chg", b5_Electron_pfRelIso03_chg, &b_b5_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("b5_Electron_phi", b5_Electron_phi, &b_b5_Electron_phi);
   fChain->SetBranchAddress("b5_Electron_pt", b5_Electron_pt, &b_b5_Electron_pt);
   fChain->SetBranchAddress("b5_Electron_r9", b5_Electron_r9, &b_b5_Electron_r9);
   fChain->SetBranchAddress("b5_Electron_scEtOverPt", b5_Electron_scEtOverPt, &b_b5_Electron_scEtOverPt);
   fChain->SetBranchAddress("b5_Electron_sieie", b5_Electron_sieie, &b_b5_Electron_sieie);
   fChain->SetBranchAddress("b5_Electron_sip3d", b5_Electron_sip3d, &b_b5_Electron_sip3d);
   fChain->SetBranchAddress("b5_Electron_mvaTTH", b5_Electron_mvaTTH, &b_b5_Electron_mvaTTH);
   fChain->SetBranchAddress("b5_Electron_charge", b5_Electron_charge, &b_b5_Electron_charge);
   fChain->SetBranchAddress("b5_Electron_cutBased", b5_Electron_cutBased, &b_b5_Electron_cutBased);
   fChain->SetBranchAddress("b5_Electron_cutBased_Fall17_V1", b5_Electron_cutBased_Fall17_V1, &b_b5_Electron_cutBased_Fall17_V1);
   fChain->SetBranchAddress("b5_Electron_jetIdx", b5_Electron_jetIdx, &b_b5_Electron_jetIdx);
   fChain->SetBranchAddress("b5_Electron_pdgId", b5_Electron_pdgId, &b_b5_Electron_pdgId);
   fChain->SetBranchAddress("b5_Electron_photonIdx", b5_Electron_photonIdx, &b_b5_Electron_photonIdx);
   fChain->SetBranchAddress("b5_Electron_tightCharge", b5_Electron_tightCharge, &b_b5_Electron_tightCharge);
   fChain->SetBranchAddress("b5_Electron_vidNestedWPBitmap", b5_Electron_vidNestedWPBitmap, &b_b5_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("b5_Electron_vidNestedWPBitmapHEEP", b5_Electron_vidNestedWPBitmapHEEP, &b_b5_Electron_vidNestedWPBitmapHEEP);
   fChain->SetBranchAddress("b5_Electron_convVeto", b5_Electron_convVeto, &b_b5_Electron_convVeto);
   fChain->SetBranchAddress("b5_Electron_cutBased_HEEP", b5_Electron_cutBased_HEEP, &b_b5_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("b5_Electron_isPFcand", b5_Electron_isPFcand, &b_b5_Electron_isPFcand);
   fChain->SetBranchAddress("b5_Electron_jetNDauCharged", b5_Electron_jetNDauCharged, &b_b5_Electron_jetNDauCharged);
   fChain->SetBranchAddress("b5_Electron_lostHits", b5_Electron_lostHits, &b_b5_Electron_lostHits);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1Iso_WP80", b5_Electron_mvaFall17V1Iso_WP80, &b_b5_Electron_mvaFall17V1Iso_WP80);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1Iso_WP90", b5_Electron_mvaFall17V1Iso_WP90, &b_b5_Electron_mvaFall17V1Iso_WP90);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1Iso_WPL", b5_Electron_mvaFall17V1Iso_WPL, &b_b5_Electron_mvaFall17V1Iso_WPL);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1noIso_WP80", b5_Electron_mvaFall17V1noIso_WP80, &b_b5_Electron_mvaFall17V1noIso_WP80);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1noIso_WP90", b5_Electron_mvaFall17V1noIso_WP90, &b_b5_Electron_mvaFall17V1noIso_WP90);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V1noIso_WPL", b5_Electron_mvaFall17V1noIso_WPL, &b_b5_Electron_mvaFall17V1noIso_WPL);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2Iso_WP80", b5_Electron_mvaFall17V2Iso_WP80, &b_b5_Electron_mvaFall17V2Iso_WP80);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2Iso_WP90", b5_Electron_mvaFall17V2Iso_WP90, &b_b5_Electron_mvaFall17V2Iso_WP90);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2Iso_WPL", b5_Electron_mvaFall17V2Iso_WPL, &b_b5_Electron_mvaFall17V2Iso_WPL);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2noIso_WP80", b5_Electron_mvaFall17V2noIso_WP80, &b_b5_Electron_mvaFall17V2noIso_WP80);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2noIso_WP90", b5_Electron_mvaFall17V2noIso_WP90, &b_b5_Electron_mvaFall17V2noIso_WP90);
   fChain->SetBranchAddress("b5_Electron_mvaFall17V2noIso_WPL", b5_Electron_mvaFall17V2noIso_WPL, &b_b5_Electron_mvaFall17V2noIso_WPL);
   fChain->SetBranchAddress("b5_Electron_seedGain", b5_Electron_seedGain, &b_b5_Electron_seedGain);
   fChain->SetBranchAddress("b5_nFatJet", &b5_nFatJet, &b_b5_nFatJet);
   fChain->SetBranchAddress("b5_FatJet_area", b5_FatJet_area, &b_b5_FatJet_area);
   fChain->SetBranchAddress("b5_FatJet_btagCMVA", b5_FatJet_btagCMVA, &b_b5_FatJet_btagCMVA);
   fChain->SetBranchAddress("b5_FatJet_btagCSVV2", b5_FatJet_btagCSVV2, &b_b5_FatJet_btagCSVV2);
   fChain->SetBranchAddress("b5_FatJet_btagDDBvL", b5_FatJet_btagDDBvL, &b_b5_FatJet_btagDDBvL);
   fChain->SetBranchAddress("b5_FatJet_btagDDBvLV2", b5_FatJet_btagDDBvLV2, &b_b5_FatJet_btagDDBvLV2);
   fChain->SetBranchAddress("b5_FatJet_btagDDBvL_noMD", b5_FatJet_btagDDBvL_noMD, &b_b5_FatJet_btagDDBvL_noMD);
   fChain->SetBranchAddress("b5_FatJet_btagDDCvB", b5_FatJet_btagDDCvB, &b_b5_FatJet_btagDDCvB);
   fChain->SetBranchAddress("b5_FatJet_btagDDCvBV2", b5_FatJet_btagDDCvBV2, &b_b5_FatJet_btagDDCvBV2);
   fChain->SetBranchAddress("b5_FatJet_btagDDCvB_noMD", b5_FatJet_btagDDCvB_noMD, &b_b5_FatJet_btagDDCvB_noMD);
   fChain->SetBranchAddress("b5_FatJet_btagDDCvL", b5_FatJet_btagDDCvL, &b_b5_FatJet_btagDDCvL);
   fChain->SetBranchAddress("b5_FatJet_btagDDCvLV2", b5_FatJet_btagDDCvLV2, &b_b5_FatJet_btagDDCvLV2);
   fChain->SetBranchAddress("b5_FatJet_btagDDCvL_noMD", b5_FatJet_btagDDCvL_noMD, &b_b5_FatJet_btagDDCvL_noMD);
   fChain->SetBranchAddress("b5_FatJet_btagDeepB", b5_FatJet_btagDeepB, &b_b5_FatJet_btagDeepB);
   fChain->SetBranchAddress("b5_FatJet_btagHbb", b5_FatJet_btagHbb, &b_b5_FatJet_btagHbb);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_H4qvsQCD", b5_FatJet_deepTagMD_H4qvsQCD, &b_b5_FatJet_deepTagMD_H4qvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_HbbvsQCD", b5_FatJet_deepTagMD_HbbvsQCD, &b_b5_FatJet_deepTagMD_HbbvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_TvsQCD", b5_FatJet_deepTagMD_TvsQCD, &b_b5_FatJet_deepTagMD_TvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_WvsQCD", b5_FatJet_deepTagMD_WvsQCD, &b_b5_FatJet_deepTagMD_WvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_ZHbbvsQCD", b5_FatJet_deepTagMD_ZHbbvsQCD, &b_b5_FatJet_deepTagMD_ZHbbvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_ZHccvsQCD", b5_FatJet_deepTagMD_ZHccvsQCD, &b_b5_FatJet_deepTagMD_ZHccvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_ZbbvsQCD", b5_FatJet_deepTagMD_ZbbvsQCD, &b_b5_FatJet_deepTagMD_ZbbvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_ZvsQCD", b5_FatJet_deepTagMD_ZvsQCD, &b_b5_FatJet_deepTagMD_ZvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_bbvsLight", b5_FatJet_deepTagMD_bbvsLight, &b_b5_FatJet_deepTagMD_bbvsLight);
   fChain->SetBranchAddress("b5_FatJet_deepTagMD_ccvsLight", b5_FatJet_deepTagMD_ccvsLight, &b_b5_FatJet_deepTagMD_ccvsLight);
   fChain->SetBranchAddress("b5_FatJet_deepTag_H", b5_FatJet_deepTag_H, &b_b5_FatJet_deepTag_H);
   fChain->SetBranchAddress("b5_FatJet_deepTag_QCD", b5_FatJet_deepTag_QCD, &b_b5_FatJet_deepTag_QCD);
   fChain->SetBranchAddress("b5_FatJet_deepTag_QCDothers", b5_FatJet_deepTag_QCDothers, &b_b5_FatJet_deepTag_QCDothers);
   fChain->SetBranchAddress("b5_FatJet_deepTag_TvsQCD", b5_FatJet_deepTag_TvsQCD, &b_b5_FatJet_deepTag_TvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTag_WvsQCD", b5_FatJet_deepTag_WvsQCD, &b_b5_FatJet_deepTag_WvsQCD);
   fChain->SetBranchAddress("b5_FatJet_deepTag_ZvsQCD", b5_FatJet_deepTag_ZvsQCD, &b_b5_FatJet_deepTag_ZvsQCD);
   fChain->SetBranchAddress("b5_FatJet_eta", b5_FatJet_eta, &b_b5_FatJet_eta);
   fChain->SetBranchAddress("b5_FatJet_mass", b5_FatJet_mass, &b_b5_FatJet_mass);
   fChain->SetBranchAddress("b5_FatJet_msoftdrop", b5_FatJet_msoftdrop, &b_b5_FatJet_msoftdrop);
   fChain->SetBranchAddress("b5_FatJet_n2b1", b5_FatJet_n2b1, &b_b5_FatJet_n2b1);
   fChain->SetBranchAddress("b5_FatJet_n3b1", b5_FatJet_n3b1, &b_b5_FatJet_n3b1);
   fChain->SetBranchAddress("b5_FatJet_particleNetMD_QCD", b5_FatJet_particleNetMD_QCD, &b_b5_FatJet_particleNetMD_QCD);
   fChain->SetBranchAddress("b5_FatJet_particleNetMD_Xbb", b5_FatJet_particleNetMD_Xbb, &b_b5_FatJet_particleNetMD_Xbb);
   fChain->SetBranchAddress("b5_FatJet_particleNetMD_Xcc", b5_FatJet_particleNetMD_Xcc, &b_b5_FatJet_particleNetMD_Xcc);
   fChain->SetBranchAddress("b5_FatJet_particleNetMD_Xqq", b5_FatJet_particleNetMD_Xqq, &b_b5_FatJet_particleNetMD_Xqq);
   fChain->SetBranchAddress("b5_FatJet_particleNet_H4qvsQCD", b5_FatJet_particleNet_H4qvsQCD, &b_b5_FatJet_particleNet_H4qvsQCD);
   fChain->SetBranchAddress("b5_FatJet_particleNet_HbbvsQCD", b5_FatJet_particleNet_HbbvsQCD, &b_b5_FatJet_particleNet_HbbvsQCD);
   fChain->SetBranchAddress("b5_FatJet_particleNet_HccvsQCD", b5_FatJet_particleNet_HccvsQCD, &b_b5_FatJet_particleNet_HccvsQCD);
   fChain->SetBranchAddress("b5_FatJet_particleNet_QCD", b5_FatJet_particleNet_QCD, &b_b5_FatJet_particleNet_QCD);
   fChain->SetBranchAddress("b5_FatJet_particleNet_TvsQCD", b5_FatJet_particleNet_TvsQCD, &b_b5_FatJet_particleNet_TvsQCD);
   fChain->SetBranchAddress("b5_FatJet_particleNet_WvsQCD", b5_FatJet_particleNet_WvsQCD, &b_b5_FatJet_particleNet_WvsQCD);
   fChain->SetBranchAddress("b5_FatJet_particleNet_ZvsQCD", b5_FatJet_particleNet_ZvsQCD, &b_b5_FatJet_particleNet_ZvsQCD);
   fChain->SetBranchAddress("b5_FatJet_phi", b5_FatJet_phi, &b_b5_FatJet_phi);
   fChain->SetBranchAddress("b5_FatJet_pt", b5_FatJet_pt, &b_b5_FatJet_pt);
   fChain->SetBranchAddress("b5_FatJet_rawFactor", b5_FatJet_rawFactor, &b_b5_FatJet_rawFactor);
   fChain->SetBranchAddress("b5_FatJet_tau1", b5_FatJet_tau1, &b_b5_FatJet_tau1);
   fChain->SetBranchAddress("b5_FatJet_tau2", b5_FatJet_tau2, &b_b5_FatJet_tau2);
   fChain->SetBranchAddress("b5_FatJet_tau3", b5_FatJet_tau3, &b_b5_FatJet_tau3);
   fChain->SetBranchAddress("b5_FatJet_tau4", b5_FatJet_tau4, &b_b5_FatJet_tau4);
   fChain->SetBranchAddress("b5_FatJet_lsf3", b5_FatJet_lsf3, &b_b5_FatJet_lsf3);
   fChain->SetBranchAddress("b5_FatJet_jetId", b5_FatJet_jetId, &b_b5_FatJet_jetId);
   fChain->SetBranchAddress("b5_FatJet_subJetIdx1", b5_FatJet_subJetIdx1, &b_b5_FatJet_subJetIdx1);
   fChain->SetBranchAddress("b5_FatJet_subJetIdx2", b5_FatJet_subJetIdx2, &b_b5_FatJet_subJetIdx2);
   fChain->SetBranchAddress("b5_FatJet_electronIdx3SJ", b5_FatJet_electronIdx3SJ, &b_b5_FatJet_electronIdx3SJ);
   fChain->SetBranchAddress("b5_FatJet_muonIdx3SJ", b5_FatJet_muonIdx3SJ, &b_b5_FatJet_muonIdx3SJ);
   fChain->SetBranchAddress("b5_nFsrPhoton", &b5_nFsrPhoton, &b_b5_nFsrPhoton);
   fChain->SetBranchAddress("b5_FsrPhoton_dROverEt2", b5_FsrPhoton_dROverEt2, &b_b5_FsrPhoton_dROverEt2);
   fChain->SetBranchAddress("b5_FsrPhoton_eta", b5_FsrPhoton_eta, &b_b5_FsrPhoton_eta);
   fChain->SetBranchAddress("b5_FsrPhoton_phi", b5_FsrPhoton_phi, &b_b5_FsrPhoton_phi);
   fChain->SetBranchAddress("b5_FsrPhoton_pt", b5_FsrPhoton_pt, &b_b5_FsrPhoton_pt);
   fChain->SetBranchAddress("b5_FsrPhoton_relIso03", b5_FsrPhoton_relIso03, &b_b5_FsrPhoton_relIso03);
   fChain->SetBranchAddress("b5_FsrPhoton_muonIdx", b5_FsrPhoton_muonIdx, &b_b5_FsrPhoton_muonIdx);
   fChain->SetBranchAddress("b5_nIsoTrack", &b5_nIsoTrack, &b_b5_nIsoTrack);
   fChain->SetBranchAddress("b5_IsoTrack_dxy", b5_IsoTrack_dxy, &b_b5_IsoTrack_dxy);
   fChain->SetBranchAddress("b5_IsoTrack_dz", b5_IsoTrack_dz, &b_b5_IsoTrack_dz);
   fChain->SetBranchAddress("b5_IsoTrack_eta", b5_IsoTrack_eta, &b_b5_IsoTrack_eta);
   fChain->SetBranchAddress("b5_IsoTrack_pfRelIso03_all", b5_IsoTrack_pfRelIso03_all, &b_b5_IsoTrack_pfRelIso03_all);
   fChain->SetBranchAddress("b5_IsoTrack_pfRelIso03_chg", b5_IsoTrack_pfRelIso03_chg, &b_b5_IsoTrack_pfRelIso03_chg);
   fChain->SetBranchAddress("b5_IsoTrack_phi", b5_IsoTrack_phi, &b_b5_IsoTrack_phi);
   fChain->SetBranchAddress("b5_IsoTrack_pt", b5_IsoTrack_pt, &b_b5_IsoTrack_pt);
   fChain->SetBranchAddress("b5_IsoTrack_miniPFRelIso_all", b5_IsoTrack_miniPFRelIso_all, &b_b5_IsoTrack_miniPFRelIso_all);
   fChain->SetBranchAddress("b5_IsoTrack_miniPFRelIso_chg", b5_IsoTrack_miniPFRelIso_chg, &b_b5_IsoTrack_miniPFRelIso_chg);
   fChain->SetBranchAddress("b5_IsoTrack_fromPV", b5_IsoTrack_fromPV, &b_b5_IsoTrack_fromPV);
   fChain->SetBranchAddress("b5_IsoTrack_pdgId", b5_IsoTrack_pdgId, &b_b5_IsoTrack_pdgId);
   fChain->SetBranchAddress("b5_IsoTrack_isHighPurityTrack", b5_IsoTrack_isHighPurityTrack, &b_b5_IsoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("b5_IsoTrack_isPFcand", b5_IsoTrack_isPFcand, &b_b5_IsoTrack_isPFcand);
   fChain->SetBranchAddress("b5_IsoTrack_isFromLostTrack", b5_IsoTrack_isFromLostTrack, &b_b5_IsoTrack_isFromLostTrack);
   fChain->SetBranchAddress("b5_nJet", &b5_nJet, &b_b5_nJet);
   fChain->SetBranchAddress("b5_Jet_area", b5_Jet_area, &b_b5_Jet_area);
   fChain->SetBranchAddress("b5_Jet_btagCMVA", b5_Jet_btagCMVA, &b_b5_Jet_btagCMVA);
   fChain->SetBranchAddress("b5_Jet_btagCSVV2", b5_Jet_btagCSVV2, &b_b5_Jet_btagCSVV2);
   fChain->SetBranchAddress("b5_Jet_btagDeepB", b5_Jet_btagDeepB, &b_b5_Jet_btagDeepB);
   fChain->SetBranchAddress("b5_Jet_btagDeepC", b5_Jet_btagDeepC, &b_b5_Jet_btagDeepC);
   fChain->SetBranchAddress("b5_Jet_btagDeepCvB", b5_Jet_btagDeepCvB, &b_b5_Jet_btagDeepCvB);
   fChain->SetBranchAddress("b5_Jet_btagDeepCvL", b5_Jet_btagDeepCvL, &b_b5_Jet_btagDeepCvL);
   fChain->SetBranchAddress("b5_Jet_btagDeepFlavB", b5_Jet_btagDeepFlavB, &b_b5_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("b5_Jet_btagDeepFlavC", b5_Jet_btagDeepFlavC, &b_b5_Jet_btagDeepFlavC);
   fChain->SetBranchAddress("b5_Jet_btagDeepFlavCvB", b5_Jet_btagDeepFlavCvB, &b_b5_Jet_btagDeepFlavCvB);
   fChain->SetBranchAddress("b5_Jet_btagDeepFlavCvL", b5_Jet_btagDeepFlavCvL, &b_b5_Jet_btagDeepFlavCvL);
   fChain->SetBranchAddress("b5_Jet_btagDeepFlavQG", b5_Jet_btagDeepFlavQG, &b_b5_Jet_btagDeepFlavQG);
   fChain->SetBranchAddress("b5_Jet_chEmEF", b5_Jet_chEmEF, &b_b5_Jet_chEmEF);
   fChain->SetBranchAddress("b5_Jet_chFPV0EF", b5_Jet_chFPV0EF, &b_b5_Jet_chFPV0EF);
   fChain->SetBranchAddress("b5_Jet_chFPV1EF", b5_Jet_chFPV1EF, &b_b5_Jet_chFPV1EF);
   fChain->SetBranchAddress("b5_Jet_chFPV2EF", b5_Jet_chFPV2EF, &b_b5_Jet_chFPV2EF);
   fChain->SetBranchAddress("b5_Jet_chFPV3EF", b5_Jet_chFPV3EF, &b_b5_Jet_chFPV3EF);
   fChain->SetBranchAddress("b5_Jet_chHEF", b5_Jet_chHEF, &b_b5_Jet_chHEF);
   fChain->SetBranchAddress("b5_Jet_eta", b5_Jet_eta, &b_b5_Jet_eta);
   fChain->SetBranchAddress("b5_Jet_hfsigmaEtaEta", b5_Jet_hfsigmaEtaEta, &b_b5_Jet_hfsigmaEtaEta);
   fChain->SetBranchAddress("b5_Jet_hfsigmaPhiPhi", b5_Jet_hfsigmaPhiPhi, &b_b5_Jet_hfsigmaPhiPhi);
   fChain->SetBranchAddress("b5_Jet_mass", b5_Jet_mass, &b_b5_Jet_mass);
   fChain->SetBranchAddress("b5_Jet_muEF", b5_Jet_muEF, &b_b5_Jet_muEF);
   fChain->SetBranchAddress("b5_Jet_muonSubtrFactor", b5_Jet_muonSubtrFactor, &b_b5_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("b5_Jet_neEmEF", b5_Jet_neEmEF, &b_b5_Jet_neEmEF);
   fChain->SetBranchAddress("b5_Jet_neHEF", b5_Jet_neHEF, &b_b5_Jet_neHEF);
   fChain->SetBranchAddress("b5_Jet_phi", b5_Jet_phi, &b_b5_Jet_phi);
   fChain->SetBranchAddress("b5_Jet_pt", b5_Jet_pt, &b_b5_Jet_pt);
   fChain->SetBranchAddress("b5_Jet_puIdDisc", b5_Jet_puIdDisc, &b_b5_Jet_puIdDisc);
   fChain->SetBranchAddress("b5_Jet_qgl", b5_Jet_qgl, &b_b5_Jet_qgl);
   fChain->SetBranchAddress("b5_Jet_rawFactor", b5_Jet_rawFactor, &b_b5_Jet_rawFactor);
   fChain->SetBranchAddress("b5_Jet_bRegCorr", b5_Jet_bRegCorr, &b_b5_Jet_bRegCorr);
   fChain->SetBranchAddress("b5_Jet_bRegRes", b5_Jet_bRegRes, &b_b5_Jet_bRegRes);
   fChain->SetBranchAddress("b5_Jet_cRegCorr", b5_Jet_cRegCorr, &b_b5_Jet_cRegCorr);
   fChain->SetBranchAddress("b5_Jet_cRegRes", b5_Jet_cRegRes, &b_b5_Jet_cRegRes);
   fChain->SetBranchAddress("b5_Jet_electronIdx1", b5_Jet_electronIdx1, &b_b5_Jet_electronIdx1);
   fChain->SetBranchAddress("b5_Jet_electronIdx2", b5_Jet_electronIdx2, &b_b5_Jet_electronIdx2);
   fChain->SetBranchAddress("b5_Jet_hfadjacentEtaStripsSize", b5_Jet_hfadjacentEtaStripsSize, &b_b5_Jet_hfadjacentEtaStripsSize);
   fChain->SetBranchAddress("b5_Jet_hfcentralEtaStripSize", b5_Jet_hfcentralEtaStripSize, &b_b5_Jet_hfcentralEtaStripSize);
   fChain->SetBranchAddress("b5_Jet_jetId", b5_Jet_jetId, &b_b5_Jet_jetId);
   fChain->SetBranchAddress("b5_Jet_muonIdx1", b5_Jet_muonIdx1, &b_b5_Jet_muonIdx1);
   fChain->SetBranchAddress("b5_Jet_muonIdx2", b5_Jet_muonIdx2, &b_b5_Jet_muonIdx2);
   fChain->SetBranchAddress("b5_Jet_nElectrons", b5_Jet_nElectrons, &b_b5_Jet_nElectrons);
   fChain->SetBranchAddress("b5_Jet_nMuons", b5_Jet_nMuons, &b_b5_Jet_nMuons);
   fChain->SetBranchAddress("b5_Jet_puId", b5_Jet_puId, &b_b5_Jet_puId);
   fChain->SetBranchAddress("b5_Jet_nConstituents", b5_Jet_nConstituents, &b_b5_Jet_nConstituents);
   fChain->SetBranchAddress("b5_L1PreFiringWeight_Dn", &b5_L1PreFiringWeight_Dn, &b_b5_L1PreFiringWeight_Dn);
   fChain->SetBranchAddress("b5_L1PreFiringWeight_Nom", &b5_L1PreFiringWeight_Nom, &b_b5_L1PreFiringWeight_Nom);
   fChain->SetBranchAddress("b5_L1PreFiringWeight_Up", &b5_L1PreFiringWeight_Up, &b_b5_L1PreFiringWeight_Up);
   fChain->SetBranchAddress("b5_MET_MetUnclustEnUpDeltaX", &b5_MET_MetUnclustEnUpDeltaX, &b_b5_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("b5_MET_MetUnclustEnUpDeltaY", &b5_MET_MetUnclustEnUpDeltaY, &b_b5_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("b5_MET_covXX", &b5_MET_covXX, &b_b5_MET_covXX);
   fChain->SetBranchAddress("b5_MET_covXY", &b5_MET_covXY, &b_b5_MET_covXY);
   fChain->SetBranchAddress("b5_MET_covYY", &b5_MET_covYY, &b_b5_MET_covYY);
   fChain->SetBranchAddress("b5_MET_phi", &b5_MET_phi, &b_b5_MET_phi);
   fChain->SetBranchAddress("b5_MET_pt", &b5_MET_pt, &b_b5_MET_pt);
   fChain->SetBranchAddress("b5_MET_significance", &b5_MET_significance, &b_b5_MET_significance);
   fChain->SetBranchAddress("b5_MET_sumEt", &b5_MET_sumEt, &b_b5_MET_sumEt);
   fChain->SetBranchAddress("b5_MET_sumPtUnclustered", &b5_MET_sumPtUnclustered, &b_b5_MET_sumPtUnclustered);
   fChain->SetBranchAddress("b5_nMuon", &b5_nMuon, &b_b5_nMuon);
   fChain->SetBranchAddress("b5_Muon_dxy", b5_Muon_dxy, &b_b5_Muon_dxy);
   fChain->SetBranchAddress("b5_Muon_dxyErr", b5_Muon_dxyErr, &b_b5_Muon_dxyErr);
   fChain->SetBranchAddress("b5_Muon_dxybs", b5_Muon_dxybs, &b_b5_Muon_dxybs);
   fChain->SetBranchAddress("b5_Muon_dz", b5_Muon_dz, &b_b5_Muon_dz);
   fChain->SetBranchAddress("b5_Muon_dzErr", b5_Muon_dzErr, &b_b5_Muon_dzErr);
   fChain->SetBranchAddress("b5_Muon_eta", b5_Muon_eta, &b_b5_Muon_eta);
   fChain->SetBranchAddress("b5_Muon_ip3d", b5_Muon_ip3d, &b_b5_Muon_ip3d);
   fChain->SetBranchAddress("b5_Muon_jetPtRelv2", b5_Muon_jetPtRelv2, &b_b5_Muon_jetPtRelv2);
   fChain->SetBranchAddress("b5_Muon_jetRelIso", b5_Muon_jetRelIso, &b_b5_Muon_jetRelIso);
   fChain->SetBranchAddress("b5_Muon_mass", b5_Muon_mass, &b_b5_Muon_mass);
   fChain->SetBranchAddress("b5_Muon_miniPFRelIso_all", b5_Muon_miniPFRelIso_all, &b_b5_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("b5_Muon_miniPFRelIso_chg", b5_Muon_miniPFRelIso_chg, &b_b5_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("b5_Muon_pfRelIso03_all", b5_Muon_pfRelIso03_all, &b_b5_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("b5_Muon_pfRelIso03_chg", b5_Muon_pfRelIso03_chg, &b_b5_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("b5_Muon_pfRelIso04_all", b5_Muon_pfRelIso04_all, &b_b5_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("b5_Muon_phi", b5_Muon_phi, &b_b5_Muon_phi);
   fChain->SetBranchAddress("b5_Muon_pt", b5_Muon_pt, &b_b5_Muon_pt);
   fChain->SetBranchAddress("b5_Muon_ptErr", b5_Muon_ptErr, &b_b5_Muon_ptErr);
   fChain->SetBranchAddress("b5_Muon_segmentComp", b5_Muon_segmentComp, &b_b5_Muon_segmentComp);
   fChain->SetBranchAddress("b5_Muon_sip3d", b5_Muon_sip3d, &b_b5_Muon_sip3d);
   fChain->SetBranchAddress("b5_Muon_softMva", b5_Muon_softMva, &b_b5_Muon_softMva);
   fChain->SetBranchAddress("b5_Muon_tkRelIso", b5_Muon_tkRelIso, &b_b5_Muon_tkRelIso);
   fChain->SetBranchAddress("b5_Muon_tunepRelPt", b5_Muon_tunepRelPt, &b_b5_Muon_tunepRelPt);
   fChain->SetBranchAddress("b5_Muon_mvaLowPt", b5_Muon_mvaLowPt, &b_b5_Muon_mvaLowPt);
   fChain->SetBranchAddress("b5_Muon_mvaTTH", b5_Muon_mvaTTH, &b_b5_Muon_mvaTTH);
   fChain->SetBranchAddress("b5_Muon_charge", b5_Muon_charge, &b_b5_Muon_charge);
   fChain->SetBranchAddress("b5_Muon_jetIdx", b5_Muon_jetIdx, &b_b5_Muon_jetIdx);
   fChain->SetBranchAddress("b5_Muon_nStations", b5_Muon_nStations, &b_b5_Muon_nStations);
   fChain->SetBranchAddress("b5_Muon_nTrackerLayers", b5_Muon_nTrackerLayers, &b_b5_Muon_nTrackerLayers);
   fChain->SetBranchAddress("b5_Muon_pdgId", b5_Muon_pdgId, &b_b5_Muon_pdgId);
   fChain->SetBranchAddress("b5_Muon_tightCharge", b5_Muon_tightCharge, &b_b5_Muon_tightCharge);
   fChain->SetBranchAddress("b5_Muon_fsrPhotonIdx", b5_Muon_fsrPhotonIdx, &b_b5_Muon_fsrPhotonIdx);
   fChain->SetBranchAddress("b5_Muon_highPtId", b5_Muon_highPtId, &b_b5_Muon_highPtId);
   fChain->SetBranchAddress("b5_Muon_highPurity", b5_Muon_highPurity, &b_b5_Muon_highPurity);
   fChain->SetBranchAddress("b5_Muon_inTimeMuon", b5_Muon_inTimeMuon, &b_b5_Muon_inTimeMuon);
   fChain->SetBranchAddress("b5_Muon_isGlobal", b5_Muon_isGlobal, &b_b5_Muon_isGlobal);
   fChain->SetBranchAddress("b5_Muon_isPFcand", b5_Muon_isPFcand, &b_b5_Muon_isPFcand);
   fChain->SetBranchAddress("b5_Muon_isTracker", b5_Muon_isTracker, &b_b5_Muon_isTracker);
   fChain->SetBranchAddress("b5_Muon_jetNDauCharged", b5_Muon_jetNDauCharged, &b_b5_Muon_jetNDauCharged);
   fChain->SetBranchAddress("b5_Muon_looseId", b5_Muon_looseId, &b_b5_Muon_looseId);
   fChain->SetBranchAddress("b5_Muon_mediumId", b5_Muon_mediumId, &b_b5_Muon_mediumId);
   fChain->SetBranchAddress("b5_Muon_mediumPromptId", b5_Muon_mediumPromptId, &b_b5_Muon_mediumPromptId);
   fChain->SetBranchAddress("b5_Muon_miniIsoId", b5_Muon_miniIsoId, &b_b5_Muon_miniIsoId);
   fChain->SetBranchAddress("b5_Muon_multiIsoId", b5_Muon_multiIsoId, &b_b5_Muon_multiIsoId);
   fChain->SetBranchAddress("b5_Muon_mvaId", b5_Muon_mvaId, &b_b5_Muon_mvaId);
   fChain->SetBranchAddress("b5_Muon_mvaLowPtId", b5_Muon_mvaLowPtId, &b_b5_Muon_mvaLowPtId);
   fChain->SetBranchAddress("b5_Muon_pfIsoId", b5_Muon_pfIsoId, &b_b5_Muon_pfIsoId);
   fChain->SetBranchAddress("b5_Muon_puppiIsoId", b5_Muon_puppiIsoId, &b_b5_Muon_puppiIsoId);
   fChain->SetBranchAddress("b5_Muon_softId", b5_Muon_softId, &b_b5_Muon_softId);
   fChain->SetBranchAddress("b5_Muon_softMvaId", b5_Muon_softMvaId, &b_b5_Muon_softMvaId);
   fChain->SetBranchAddress("b5_Muon_tightId", b5_Muon_tightId, &b_b5_Muon_tightId);
   fChain->SetBranchAddress("b5_Muon_tkIsoId", b5_Muon_tkIsoId, &b_b5_Muon_tkIsoId);
   fChain->SetBranchAddress("b5_Muon_triggerIdLoose", b5_Muon_triggerIdLoose, &b_b5_Muon_triggerIdLoose);
   fChain->SetBranchAddress("b5_nPhoton", &b5_nPhoton, &b_b5_nPhoton);
   fChain->SetBranchAddress("b5_Photon_eCorr", b5_Photon_eCorr, &b_b5_Photon_eCorr);
   fChain->SetBranchAddress("b5_Photon_energyErr", b5_Photon_energyErr, &b_b5_Photon_energyErr);
   fChain->SetBranchAddress("b5_Photon_eta", b5_Photon_eta, &b_b5_Photon_eta);
   fChain->SetBranchAddress("b5_Photon_hoe", b5_Photon_hoe, &b_b5_Photon_hoe);
   fChain->SetBranchAddress("b5_Photon_mass", b5_Photon_mass, &b_b5_Photon_mass);
   fChain->SetBranchAddress("b5_Photon_mvaID", b5_Photon_mvaID, &b_b5_Photon_mvaID);
   fChain->SetBranchAddress("b5_Photon_mvaID_Fall17V1p1", b5_Photon_mvaID_Fall17V1p1, &b_b5_Photon_mvaID_Fall17V1p1);
   fChain->SetBranchAddress("b5_Photon_pfRelIso03_all", b5_Photon_pfRelIso03_all, &b_b5_Photon_pfRelIso03_all);
   fChain->SetBranchAddress("b5_Photon_pfRelIso03_chg", b5_Photon_pfRelIso03_chg, &b_b5_Photon_pfRelIso03_chg);
   fChain->SetBranchAddress("b5_Photon_phi", b5_Photon_phi, &b_b5_Photon_phi);
   fChain->SetBranchAddress("b5_Photon_pt", b5_Photon_pt, &b_b5_Photon_pt);
   fChain->SetBranchAddress("b5_Photon_r9", b5_Photon_r9, &b_b5_Photon_r9);
   fChain->SetBranchAddress("b5_Photon_sieie", b5_Photon_sieie, &b_b5_Photon_sieie);
   fChain->SetBranchAddress("b5_Photon_charge", b5_Photon_charge, &b_b5_Photon_charge);
   fChain->SetBranchAddress("b5_Photon_cutBased", b5_Photon_cutBased, &b_b5_Photon_cutBased);
   fChain->SetBranchAddress("b5_Photon_cutBased_Fall17V1Bitmap", b5_Photon_cutBased_Fall17V1Bitmap, &b_b5_Photon_cutBased_Fall17V1Bitmap);
   fChain->SetBranchAddress("b5_Photon_electronIdx", b5_Photon_electronIdx, &b_b5_Photon_electronIdx);
   fChain->SetBranchAddress("b5_Photon_jetIdx", b5_Photon_jetIdx, &b_b5_Photon_jetIdx);
   fChain->SetBranchAddress("b5_Photon_pdgId", b5_Photon_pdgId, &b_b5_Photon_pdgId);
   fChain->SetBranchAddress("b5_Photon_vidNestedWPBitmap", b5_Photon_vidNestedWPBitmap, &b_b5_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("b5_Photon_electronVeto", b5_Photon_electronVeto, &b_b5_Photon_electronVeto);
   fChain->SetBranchAddress("b5_Photon_isScEtaEB", b5_Photon_isScEtaEB, &b_b5_Photon_isScEtaEB);
   fChain->SetBranchAddress("b5_Photon_isScEtaEE", b5_Photon_isScEtaEE, &b_b5_Photon_isScEtaEE);
   fChain->SetBranchAddress("b5_Photon_mvaID_WP80", b5_Photon_mvaID_WP80, &b_b5_Photon_mvaID_WP80);
   fChain->SetBranchAddress("b5_Photon_mvaID_WP90", b5_Photon_mvaID_WP90, &b_b5_Photon_mvaID_WP90);
   fChain->SetBranchAddress("b5_Photon_pixelSeed", b5_Photon_pixelSeed, &b_b5_Photon_pixelSeed);
   fChain->SetBranchAddress("b5_Photon_seedGain", b5_Photon_seedGain, &b_b5_Photon_seedGain);
   fChain->SetBranchAddress("b5_PuppiMET_phi", &b5_PuppiMET_phi, &b_b5_PuppiMET_phi);
   fChain->SetBranchAddress("b5_PuppiMET_phiJERDown", &b5_PuppiMET_phiJERDown, &b_b5_PuppiMET_phiJERDown);
   fChain->SetBranchAddress("b5_PuppiMET_phiJERUp", &b5_PuppiMET_phiJERUp, &b_b5_PuppiMET_phiJERUp);
   fChain->SetBranchAddress("b5_PuppiMET_phiJESDown", &b5_PuppiMET_phiJESDown, &b_b5_PuppiMET_phiJESDown);
   fChain->SetBranchAddress("b5_PuppiMET_phiJESUp", &b5_PuppiMET_phiJESUp, &b_b5_PuppiMET_phiJESUp);
   fChain->SetBranchAddress("b5_PuppiMET_phiUnclusteredDown", &b5_PuppiMET_phiUnclusteredDown, &b_b5_PuppiMET_phiUnclusteredDown);
   fChain->SetBranchAddress("b5_PuppiMET_phiUnclusteredUp", &b5_PuppiMET_phiUnclusteredUp, &b_b5_PuppiMET_phiUnclusteredUp);
   fChain->SetBranchAddress("b5_PuppiMET_pt", &b5_PuppiMET_pt, &b_b5_PuppiMET_pt);
   fChain->SetBranchAddress("b5_PuppiMET_ptJERDown", &b5_PuppiMET_ptJERDown, &b_b5_PuppiMET_ptJERDown);
   fChain->SetBranchAddress("b5_PuppiMET_ptJERUp", &b5_PuppiMET_ptJERUp, &b_b5_PuppiMET_ptJERUp);
   fChain->SetBranchAddress("b5_PuppiMET_ptJESDown", &b5_PuppiMET_ptJESDown, &b_b5_PuppiMET_ptJESDown);
   fChain->SetBranchAddress("b5_PuppiMET_ptJESUp", &b5_PuppiMET_ptJESUp, &b_b5_PuppiMET_ptJESUp);
   fChain->SetBranchAddress("b5_PuppiMET_ptUnclusteredDown", &b5_PuppiMET_ptUnclusteredDown, &b_b5_PuppiMET_ptUnclusteredDown);
   fChain->SetBranchAddress("b5_PuppiMET_ptUnclusteredUp", &b5_PuppiMET_ptUnclusteredUp, &b_b5_PuppiMET_ptUnclusteredUp);
   fChain->SetBranchAddress("b5_PuppiMET_sumEt", &b5_PuppiMET_sumEt, &b_b5_PuppiMET_sumEt);
   fChain->SetBranchAddress("b5_RawMET_phi", &b5_RawMET_phi, &b_b5_RawMET_phi);
   fChain->SetBranchAddress("b5_RawMET_pt", &b5_RawMET_pt, &b_b5_RawMET_pt);
   fChain->SetBranchAddress("b5_RawMET_sumEt", &b5_RawMET_sumEt, &b_b5_RawMET_sumEt);
   fChain->SetBranchAddress("b5_RawPuppiMET_phi", &b5_RawPuppiMET_phi, &b_b5_RawPuppiMET_phi);
   fChain->SetBranchAddress("b5_RawPuppiMET_pt", &b5_RawPuppiMET_pt, &b_b5_RawPuppiMET_pt);
   fChain->SetBranchAddress("b5_RawPuppiMET_sumEt", &b5_RawPuppiMET_sumEt, &b_b5_RawPuppiMET_sumEt);
   fChain->SetBranchAddress("b5_fixedGridRhoFastjetAll", &b5_fixedGridRhoFastjetAll, &b_b5_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("b5_fixedGridRhoFastjetCentral", &b5_fixedGridRhoFastjetCentral, &b_b5_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("b5_fixedGridRhoFastjetCentralCalo", &b5_fixedGridRhoFastjetCentralCalo, &b_b5_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("b5_fixedGridRhoFastjetCentralChargedPileUp", &b5_fixedGridRhoFastjetCentralChargedPileUp, &b_b5_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("b5_fixedGridRhoFastjetCentralNeutral", &b5_fixedGridRhoFastjetCentralNeutral, &b_b5_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("b5_nSoftActivityJet", &b5_nSoftActivityJet, &b_b5_nSoftActivityJet);
   fChain->SetBranchAddress("b5_SoftActivityJet_eta", b5_SoftActivityJet_eta, &b_b5_SoftActivityJet_eta);
   fChain->SetBranchAddress("b5_SoftActivityJet_phi", b5_SoftActivityJet_phi, &b_b5_SoftActivityJet_phi);
   fChain->SetBranchAddress("b5_SoftActivityJet_pt", b5_SoftActivityJet_pt, &b_b5_SoftActivityJet_pt);
   fChain->SetBranchAddress("b5_SoftActivityJetHT", &b5_SoftActivityJetHT, &b_b5_SoftActivityJetHT);
   fChain->SetBranchAddress("b5_SoftActivityJetHT10", &b5_SoftActivityJetHT10, &b_b5_SoftActivityJetHT10);
   fChain->SetBranchAddress("b5_SoftActivityJetHT2", &b5_SoftActivityJetHT2, &b_b5_SoftActivityJetHT2);
   fChain->SetBranchAddress("b5_SoftActivityJetHT5", &b5_SoftActivityJetHT5, &b_b5_SoftActivityJetHT5);
   fChain->SetBranchAddress("b5_SoftActivityJetNjets10", &b5_SoftActivityJetNjets10, &b_b5_SoftActivityJetNjets10);
   fChain->SetBranchAddress("b5_SoftActivityJetNjets2", &b5_SoftActivityJetNjets2, &b_b5_SoftActivityJetNjets2);
   fChain->SetBranchAddress("b5_SoftActivityJetNjets5", &b5_SoftActivityJetNjets5, &b_b5_SoftActivityJetNjets5);
   fChain->SetBranchAddress("b5_nSubJet", &b5_nSubJet, &b_b5_nSubJet);
   fChain->SetBranchAddress("b5_SubJet_btagCMVA", b5_SubJet_btagCMVA, &b_b5_SubJet_btagCMVA);
   fChain->SetBranchAddress("b5_SubJet_btagCSVV2", b5_SubJet_btagCSVV2, &b_b5_SubJet_btagCSVV2);
   fChain->SetBranchAddress("b5_SubJet_btagDeepB", b5_SubJet_btagDeepB, &b_b5_SubJet_btagDeepB);
   fChain->SetBranchAddress("b5_SubJet_eta", b5_SubJet_eta, &b_b5_SubJet_eta);
   fChain->SetBranchAddress("b5_SubJet_mass", b5_SubJet_mass, &b_b5_SubJet_mass);
   fChain->SetBranchAddress("b5_SubJet_n2b1", b5_SubJet_n2b1, &b_b5_SubJet_n2b1);
   fChain->SetBranchAddress("b5_SubJet_n3b1", b5_SubJet_n3b1, &b_b5_SubJet_n3b1);
   fChain->SetBranchAddress("b5_SubJet_phi", b5_SubJet_phi, &b_b5_SubJet_phi);
   fChain->SetBranchAddress("b5_SubJet_pt", b5_SubJet_pt, &b_b5_SubJet_pt);
   fChain->SetBranchAddress("b5_SubJet_rawFactor", b5_SubJet_rawFactor, &b_b5_SubJet_rawFactor);
   fChain->SetBranchAddress("b5_SubJet_tau1", b5_SubJet_tau1, &b_b5_SubJet_tau1);
   fChain->SetBranchAddress("b5_SubJet_tau2", b5_SubJet_tau2, &b_b5_SubJet_tau2);
   fChain->SetBranchAddress("b5_SubJet_tau3", b5_SubJet_tau3, &b_b5_SubJet_tau3);
   fChain->SetBranchAddress("b5_SubJet_tau4", b5_SubJet_tau4, &b_b5_SubJet_tau4);
   fChain->SetBranchAddress("b5_nTau", &b5_nTau, &b_b5_nTau);
   fChain->SetBranchAddress("b5_Tau_chargedIso", b5_Tau_chargedIso, &b_b5_Tau_chargedIso);
   fChain->SetBranchAddress("b5_Tau_dxy", b5_Tau_dxy, &b_b5_Tau_dxy);
   fChain->SetBranchAddress("b5_Tau_dz", b5_Tau_dz, &b_b5_Tau_dz);
   fChain->SetBranchAddress("b5_Tau_eta", b5_Tau_eta, &b_b5_Tau_eta);
   fChain->SetBranchAddress("b5_Tau_leadTkDeltaEta", b5_Tau_leadTkDeltaEta, &b_b5_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("b5_Tau_leadTkDeltaPhi", b5_Tau_leadTkDeltaPhi, &b_b5_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("b5_Tau_leadTkPtOverTauPt", b5_Tau_leadTkPtOverTauPt, &b_b5_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("b5_Tau_mass", b5_Tau_mass, &b_b5_Tau_mass);
   fChain->SetBranchAddress("b5_Tau_neutralIso", b5_Tau_neutralIso, &b_b5_Tau_neutralIso);
   fChain->SetBranchAddress("b5_Tau_phi", b5_Tau_phi, &b_b5_Tau_phi);
   fChain->SetBranchAddress("b5_Tau_photonsOutsideSignalCone", b5_Tau_photonsOutsideSignalCone, &b_b5_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("b5_Tau_pt", b5_Tau_pt, &b_b5_Tau_pt);
   fChain->SetBranchAddress("b5_Tau_puCorr", b5_Tau_puCorr, &b_b5_Tau_puCorr);
   fChain->SetBranchAddress("b5_Tau_rawAntiEle", b5_Tau_rawAntiEle, &b_b5_Tau_rawAntiEle);
   fChain->SetBranchAddress("b5_Tau_rawAntiEle2018", b5_Tau_rawAntiEle2018, &b_b5_Tau_rawAntiEle2018);
   fChain->SetBranchAddress("b5_Tau_rawDeepTau2017v2p1VSe", b5_Tau_rawDeepTau2017v2p1VSe, &b_b5_Tau_rawDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("b5_Tau_rawDeepTau2017v2p1VSjet", b5_Tau_rawDeepTau2017v2p1VSjet, &b_b5_Tau_rawDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("b5_Tau_rawDeepTau2017v2p1VSmu", b5_Tau_rawDeepTau2017v2p1VSmu, &b_b5_Tau_rawDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("b5_Tau_rawIso", b5_Tau_rawIso, &b_b5_Tau_rawIso);
   fChain->SetBranchAddress("b5_Tau_rawIsodR03", b5_Tau_rawIsodR03, &b_b5_Tau_rawIsodR03);
   fChain->SetBranchAddress("b5_Tau_rawMVAnewDM2017v2", b5_Tau_rawMVAnewDM2017v2, &b_b5_Tau_rawMVAnewDM2017v2);
   fChain->SetBranchAddress("b5_Tau_rawMVAoldDM", b5_Tau_rawMVAoldDM, &b_b5_Tau_rawMVAoldDM);
   fChain->SetBranchAddress("b5_Tau_rawMVAoldDM2017v1", b5_Tau_rawMVAoldDM2017v1, &b_b5_Tau_rawMVAoldDM2017v1);
   fChain->SetBranchAddress("b5_Tau_rawMVAoldDM2017v2", b5_Tau_rawMVAoldDM2017v2, &b_b5_Tau_rawMVAoldDM2017v2);
   fChain->SetBranchAddress("b5_Tau_rawMVAoldDMdR032017v2", b5_Tau_rawMVAoldDMdR032017v2, &b_b5_Tau_rawMVAoldDMdR032017v2);
   fChain->SetBranchAddress("b5_Tau_charge", b5_Tau_charge, &b_b5_Tau_charge);
   fChain->SetBranchAddress("b5_Tau_decayMode", b5_Tau_decayMode, &b_b5_Tau_decayMode);
   fChain->SetBranchAddress("b5_Tau_jetIdx", b5_Tau_jetIdx, &b_b5_Tau_jetIdx);
   fChain->SetBranchAddress("b5_Tau_rawAntiEleCat", b5_Tau_rawAntiEleCat, &b_b5_Tau_rawAntiEleCat);
   fChain->SetBranchAddress("b5_Tau_rawAntiEleCat2018", b5_Tau_rawAntiEleCat2018, &b_b5_Tau_rawAntiEleCat2018);
   fChain->SetBranchAddress("b5_Tau_idAntiEle", b5_Tau_idAntiEle, &b_b5_Tau_idAntiEle);
   fChain->SetBranchAddress("b5_Tau_idAntiEle2018", b5_Tau_idAntiEle2018, &b_b5_Tau_idAntiEle2018);
   fChain->SetBranchAddress("b5_Tau_idAntiEleDeadECal", b5_Tau_idAntiEleDeadECal, &b_b5_Tau_idAntiEleDeadECal);
   fChain->SetBranchAddress("b5_Tau_idAntiMu", b5_Tau_idAntiMu, &b_b5_Tau_idAntiMu);
   fChain->SetBranchAddress("b5_Tau_idDecayMode", b5_Tau_idDecayMode, &b_b5_Tau_idDecayMode);
   fChain->SetBranchAddress("b5_Tau_idDecayModeNewDMs", b5_Tau_idDecayModeNewDMs, &b_b5_Tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("b5_Tau_idDeepTau2017v2p1VSe", b5_Tau_idDeepTau2017v2p1VSe, &b_b5_Tau_idDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("b5_Tau_idDeepTau2017v2p1VSjet", b5_Tau_idDeepTau2017v2p1VSjet, &b_b5_Tau_idDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("b5_Tau_idDeepTau2017v2p1VSmu", b5_Tau_idDeepTau2017v2p1VSmu, &b_b5_Tau_idDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("b5_Tau_idMVAnewDM2017v2", b5_Tau_idMVAnewDM2017v2, &b_b5_Tau_idMVAnewDM2017v2);
   fChain->SetBranchAddress("b5_Tau_idMVAoldDM", b5_Tau_idMVAoldDM, &b_b5_Tau_idMVAoldDM);
   fChain->SetBranchAddress("b5_Tau_idMVAoldDM2017v1", b5_Tau_idMVAoldDM2017v1, &b_b5_Tau_idMVAoldDM2017v1);
   fChain->SetBranchAddress("b5_Tau_idMVAoldDM2017v2", b5_Tau_idMVAoldDM2017v2, &b_b5_Tau_idMVAoldDM2017v2);
   fChain->SetBranchAddress("b5_Tau_idMVAoldDMdR032017v2", b5_Tau_idMVAoldDMdR032017v2, &b_b5_Tau_idMVAoldDMdR032017v2);
   fChain->SetBranchAddress("b5_TkMET_phi", &b5_TkMET_phi, &b_b5_TkMET_phi);
   fChain->SetBranchAddress("b5_TkMET_pt", &b5_TkMET_pt, &b_b5_TkMET_pt);
   fChain->SetBranchAddress("b5_TkMET_sumEt", &b5_TkMET_sumEt, &b_b5_TkMET_sumEt);
   fChain->SetBranchAddress("b5_nTrigObj", &b5_nTrigObj, &b_b5_nTrigObj);
   fChain->SetBranchAddress("b5_TrigObj_pt", b5_TrigObj_pt, &b_b5_TrigObj_pt);
   fChain->SetBranchAddress("b5_TrigObj_eta", b5_TrigObj_eta, &b_b5_TrigObj_eta);
   fChain->SetBranchAddress("b5_TrigObj_phi", b5_TrigObj_phi, &b_b5_TrigObj_phi);
   fChain->SetBranchAddress("b5_TrigObj_l1pt", b5_TrigObj_l1pt, &b_b5_TrigObj_l1pt);
   fChain->SetBranchAddress("b5_TrigObj_l1pt_2", b5_TrigObj_l1pt_2, &b_b5_TrigObj_l1pt_2);
   fChain->SetBranchAddress("b5_TrigObj_l2pt", b5_TrigObj_l2pt, &b_b5_TrigObj_l2pt);
   fChain->SetBranchAddress("b5_TrigObj_id", b5_TrigObj_id, &b_b5_TrigObj_id);
   fChain->SetBranchAddress("b5_TrigObj_l1iso", b5_TrigObj_l1iso, &b_b5_TrigObj_l1iso);
   fChain->SetBranchAddress("b5_TrigObj_l1charge", b5_TrigObj_l1charge, &b_b5_TrigObj_l1charge);
   fChain->SetBranchAddress("b5_TrigObj_filterBits", b5_TrigObj_filterBits, &b_b5_TrigObj_filterBits);
   fChain->SetBranchAddress("b5_nOtherPV", &b5_nOtherPV, &b_b5_nOtherPV);
   fChain->SetBranchAddress("b5_OtherPV_z", b5_OtherPV_z, &b_b5_OtherPV_z);
   fChain->SetBranchAddress("b5_PV_ndof", &b5_PV_ndof, &b_b5_PV_ndof);
   fChain->SetBranchAddress("b5_PV_x", &b5_PV_x, &b_b5_PV_x);
   fChain->SetBranchAddress("b5_PV_y", &b5_PV_y, &b_b5_PV_y);
   fChain->SetBranchAddress("b5_PV_z", &b5_PV_z, &b_b5_PV_z);
   fChain->SetBranchAddress("b5_PV_chi2", &b5_PV_chi2, &b_b5_PV_chi2);
   fChain->SetBranchAddress("b5_PV_score", &b5_PV_score, &b_b5_PV_score);
   fChain->SetBranchAddress("b5_PV_npvs", &b5_PV_npvs, &b_b5_PV_npvs);
   fChain->SetBranchAddress("b5_PV_npvsGood", &b5_PV_npvsGood, &b_b5_PV_npvsGood);
   fChain->SetBranchAddress("b5_nSV", &b5_nSV, &b_b5_nSV);
   fChain->SetBranchAddress("b5_SV_dlen", b5_SV_dlen, &b_b5_SV_dlen);
   fChain->SetBranchAddress("b5_SV_dlenSig", b5_SV_dlenSig, &b_b5_SV_dlenSig);
   fChain->SetBranchAddress("b5_SV_dxy", b5_SV_dxy, &b_b5_SV_dxy);
   fChain->SetBranchAddress("b5_SV_dxySig", b5_SV_dxySig, &b_b5_SV_dxySig);
   fChain->SetBranchAddress("b5_SV_pAngle", b5_SV_pAngle, &b_b5_SV_pAngle);
   fChain->SetBranchAddress("b5_Electron_cleanmask", b5_Electron_cleanmask, &b_b5_Electron_cleanmask);
   fChain->SetBranchAddress("b5_Jet_cleanmask", b5_Jet_cleanmask, &b_b5_Jet_cleanmask);
   fChain->SetBranchAddress("b5_Muon_cleanmask", b5_Muon_cleanmask, &b_b5_Muon_cleanmask);
   fChain->SetBranchAddress("b5_Photon_cleanmask", b5_Photon_cleanmask, &b_b5_Photon_cleanmask);
   fChain->SetBranchAddress("b5_Tau_cleanmask", b5_Tau_cleanmask, &b_b5_Tau_cleanmask);
   fChain->SetBranchAddress("b5_SV_chi2", b5_SV_chi2, &b_b5_SV_chi2);
   fChain->SetBranchAddress("b5_SV_eta", b5_SV_eta, &b_b5_SV_eta);
   fChain->SetBranchAddress("b5_SV_mass", b5_SV_mass, &b_b5_SV_mass);
   fChain->SetBranchAddress("b5_SV_ndof", b5_SV_ndof, &b_b5_SV_ndof);
   fChain->SetBranchAddress("b5_SV_phi", b5_SV_phi, &b_b5_SV_phi);
   fChain->SetBranchAddress("b5_SV_pt", b5_SV_pt, &b_b5_SV_pt);
   fChain->SetBranchAddress("b5_SV_x", b5_SV_x, &b_b5_SV_x);
   fChain->SetBranchAddress("b5_SV_y", b5_SV_y, &b_b5_SV_y);
   fChain->SetBranchAddress("b5_SV_z", b5_SV_z, &b_b5_SV_z);
   fChain->SetBranchAddress("b5_SV_ntracks", b5_SV_ntracks, &b_b5_SV_ntracks);
   fChain->SetBranchAddress("b5_L1_AlwaysTrue", &b5_L1_AlwaysTrue, &b_b5_L1_AlwaysTrue);
   fChain->SetBranchAddress("b5_L1_BPTX_AND_Ref1_VME", &b5_L1_BPTX_AND_Ref1_VME, &b_b5_L1_BPTX_AND_Ref1_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_AND_Ref3_VME", &b5_L1_BPTX_AND_Ref3_VME, &b_b5_L1_BPTX_AND_Ref3_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_AND_Ref4_VME", &b5_L1_BPTX_AND_Ref4_VME, &b_b5_L1_BPTX_AND_Ref4_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_BeamGas_B1_VME", &b5_L1_BPTX_BeamGas_B1_VME, &b_b5_L1_BPTX_BeamGas_B1_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_BeamGas_B2_VME", &b5_L1_BPTX_BeamGas_B2_VME, &b_b5_L1_BPTX_BeamGas_B2_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_BeamGas_Ref1_VME", &b5_L1_BPTX_BeamGas_Ref1_VME, &b_b5_L1_BPTX_BeamGas_Ref1_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_BeamGas_Ref2_VME", &b5_L1_BPTX_BeamGas_Ref2_VME, &b_b5_L1_BPTX_BeamGas_Ref2_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_NotOR_VME", &b5_L1_BPTX_NotOR_VME, &b_b5_L1_BPTX_NotOR_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_OR_Ref3_VME", &b5_L1_BPTX_OR_Ref3_VME, &b_b5_L1_BPTX_OR_Ref3_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_OR_Ref4_VME", &b5_L1_BPTX_OR_Ref4_VME, &b_b5_L1_BPTX_OR_Ref4_VME);
   fChain->SetBranchAddress("b5_L1_BPTX_RefAND_VME", &b5_L1_BPTX_RefAND_VME, &b_b5_L1_BPTX_RefAND_VME);
   fChain->SetBranchAddress("b5_L1_BptxMinus", &b5_L1_BptxMinus, &b_b5_L1_BptxMinus);
   fChain->SetBranchAddress("b5_L1_BptxOR", &b5_L1_BptxOR, &b_b5_L1_BptxOR);
   fChain->SetBranchAddress("b5_L1_BptxPlus", &b5_L1_BptxPlus, &b_b5_L1_BptxPlus);
   fChain->SetBranchAddress("b5_L1_BptxXOR", &b5_L1_BptxXOR, &b_b5_L1_BptxXOR);
   fChain->SetBranchAddress("b5_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &b5_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_b5_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("b5_L1_DoubleEG8er2p5_HTT260er", &b5_L1_DoubleEG8er2p5_HTT260er, &b_b5_L1_DoubleEG8er2p5_HTT260er);
   fChain->SetBranchAddress("b5_L1_DoubleEG8er2p5_HTT280er", &b5_L1_DoubleEG8er2p5_HTT280er, &b_b5_L1_DoubleEG8er2p5_HTT280er);
   fChain->SetBranchAddress("b5_L1_DoubleEG8er2p5_HTT300er", &b5_L1_DoubleEG8er2p5_HTT300er, &b_b5_L1_DoubleEG8er2p5_HTT300er);
   fChain->SetBranchAddress("b5_L1_DoubleEG8er2p5_HTT320er", &b5_L1_DoubleEG8er2p5_HTT320er, &b_b5_L1_DoubleEG8er2p5_HTT320er);
   fChain->SetBranchAddress("b5_L1_DoubleEG8er2p5_HTT340er", &b5_L1_DoubleEG8er2p5_HTT340er, &b_b5_L1_DoubleEG8er2p5_HTT340er);
   fChain->SetBranchAddress("b5_L1_DoubleEG_15_10_er2p5", &b5_L1_DoubleEG_15_10_er2p5, &b_b5_L1_DoubleEG_15_10_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_20_10_er2p5", &b5_L1_DoubleEG_20_10_er2p5, &b_b5_L1_DoubleEG_20_10_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_22_10_er2p5", &b5_L1_DoubleEG_22_10_er2p5, &b_b5_L1_DoubleEG_22_10_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_25_12_er2p5", &b5_L1_DoubleEG_25_12_er2p5, &b_b5_L1_DoubleEG_25_12_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_25_14_er2p5", &b5_L1_DoubleEG_25_14_er2p5, &b_b5_L1_DoubleEG_25_14_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_27_14_er2p5", &b5_L1_DoubleEG_27_14_er2p5, &b_b5_L1_DoubleEG_27_14_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_LooseIso20_10_er2p5", &b5_L1_DoubleEG_LooseIso20_10_er2p5, &b_b5_L1_DoubleEG_LooseIso20_10_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_LooseIso22_10_er2p5", &b5_L1_DoubleEG_LooseIso22_10_er2p5, &b_b5_L1_DoubleEG_LooseIso22_10_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_LooseIso22_12_er2p5", &b5_L1_DoubleEG_LooseIso22_12_er2p5, &b_b5_L1_DoubleEG_LooseIso22_12_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleEG_LooseIso25_12_er2p5", &b5_L1_DoubleEG_LooseIso25_12_er2p5, &b_b5_L1_DoubleEG_LooseIso25_12_er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleIsoTau32er2p1", &b5_L1_DoubleIsoTau32er2p1, &b_b5_L1_DoubleIsoTau32er2p1);
   fChain->SetBranchAddress("b5_L1_DoubleIsoTau34er2p1", &b5_L1_DoubleIsoTau34er2p1, &b_b5_L1_DoubleIsoTau34er2p1);
   fChain->SetBranchAddress("b5_L1_DoubleIsoTau36er2p1", &b5_L1_DoubleIsoTau36er2p1, &b_b5_L1_DoubleIsoTau36er2p1);
   fChain->SetBranchAddress("b5_L1_DoubleJet100er2p3_dEta_Max1p6", &b5_L1_DoubleJet100er2p3_dEta_Max1p6, &b_b5_L1_DoubleJet100er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("b5_L1_DoubleJet100er2p5", &b5_L1_DoubleJet100er2p5, &b_b5_L1_DoubleJet100er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet112er2p3_dEta_Max1p6", &b5_L1_DoubleJet112er2p3_dEta_Max1p6, &b_b5_L1_DoubleJet112er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("b5_L1_DoubleJet120er2p5", &b5_L1_DoubleJet120er2p5, &b_b5_L1_DoubleJet120er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet150er2p5", &b5_L1_DoubleJet150er2p5, &b_b5_L1_DoubleJet150er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &b5_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5, &b_b5_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &b5_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5, &b_b5_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &b5_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5, &b_b5_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &b5_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5, &b_b5_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &b5_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5, &b_b5_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &b5_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5, &b_b5_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &b5_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp, &b_b5_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp);
   fChain->SetBranchAddress("b5_L1_DoubleJet40er2p5", &b5_L1_DoubleJet40er2p5, &b_b5_L1_DoubleJet40er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &b5_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620, &b_b5_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620);
   fChain->SetBranchAddress("b5_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &b5_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620, &b_b5_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620);
   fChain->SetBranchAddress("b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620, &b_b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620);
   fChain->SetBranchAddress("b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28, &b_b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28);
   fChain->SetBranchAddress("b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620, &b_b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620);
   fChain->SetBranchAddress("b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28, &b_b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28);
   fChain->SetBranchAddress("b5_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &b5_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ, &b_b5_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ);
   fChain->SetBranchAddress("b5_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &b5_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp, &b_b5_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp);
   fChain->SetBranchAddress("b5_L1_DoubleJet_80_30_Mass_Min420_Mu8", &b5_L1_DoubleJet_80_30_Mass_Min420_Mu8, &b_b5_L1_DoubleJet_80_30_Mass_Min420_Mu8);
   fChain->SetBranchAddress("b5_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &b5_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620, &b_b5_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620);
   fChain->SetBranchAddress("b5_L1_DoubleLooseIsoEG22er2p1", &b5_L1_DoubleLooseIsoEG22er2p1, &b_b5_L1_DoubleLooseIsoEG22er2p1);
   fChain->SetBranchAddress("b5_L1_DoubleLooseIsoEG24er2p1", &b5_L1_DoubleLooseIsoEG24er2p1, &b_b5_L1_DoubleLooseIsoEG24er2p1);
   fChain->SetBranchAddress("b5_L1_DoubleMu0", &b5_L1_DoubleMu0, &b_b5_L1_DoubleMu0);
   fChain->SetBranchAddress("b5_L1_DoubleMu0_Mass_Min1", &b5_L1_DoubleMu0_Mass_Min1, &b_b5_L1_DoubleMu0_Mass_Min1);
   fChain->SetBranchAddress("b5_L1_DoubleMu0_OQ", &b5_L1_DoubleMu0_OQ, &b_b5_L1_DoubleMu0_OQ);
   fChain->SetBranchAddress("b5_L1_DoubleMu0_SQ", &b5_L1_DoubleMu0_SQ, &b_b5_L1_DoubleMu0_SQ);
   fChain->SetBranchAddress("b5_L1_DoubleMu0_SQ_OS", &b5_L1_DoubleMu0_SQ_OS, &b_b5_L1_DoubleMu0_SQ_OS);
   fChain->SetBranchAddress("b5_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &b5_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8, &b_b5_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8);
   fChain->SetBranchAddress("b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4, &b_b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("b5_L1_DoubleMu0er1p5_SQ", &b5_L1_DoubleMu0er1p5_SQ, &b_b5_L1_DoubleMu0er1p5_SQ);
   fChain->SetBranchAddress("b5_L1_DoubleMu0er1p5_SQ_OS", &b5_L1_DoubleMu0er1p5_SQ_OS, &b_b5_L1_DoubleMu0er1p5_SQ_OS);
   fChain->SetBranchAddress("b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4, &b_b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4", &b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4, &b_b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4);
   fChain->SetBranchAddress("b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4, &b_b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4", &b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4, &b_b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4);
   fChain->SetBranchAddress("b5_L1_DoubleMu10_SQ", &b5_L1_DoubleMu10_SQ, &b_b5_L1_DoubleMu10_SQ);
   fChain->SetBranchAddress("b5_L1_DoubleMu18er2p1", &b5_L1_DoubleMu18er2p1, &b_b5_L1_DoubleMu18er2p1);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &b5_L1_DoubleMu3_OS_DoubleEG7p5Upsilon, &b_b5_L1_DoubleMu3_OS_DoubleEG7p5Upsilon);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_SQ_ETMHF50_HTT60er", &b5_L1_DoubleMu3_SQ_ETMHF50_HTT60er, &b_b5_L1_DoubleMu3_SQ_ETMHF50_HTT60er);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5, &b_b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5, &b_b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &b5_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5, &b_b5_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_SQ_HTT220er", &b5_L1_DoubleMu3_SQ_HTT220er, &b_b5_L1_DoubleMu3_SQ_HTT220er);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_SQ_HTT240er", &b5_L1_DoubleMu3_SQ_HTT240er, &b_b5_L1_DoubleMu3_SQ_HTT240er);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_SQ_HTT260er", &b5_L1_DoubleMu3_SQ_HTT260er, &b_b5_L1_DoubleMu3_SQ_HTT260er);
   fChain->SetBranchAddress("b5_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &b5_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8, &b_b5_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8);
   fChain->SetBranchAddress("b5_L1_DoubleMu4_SQ_EG9er2p5", &b5_L1_DoubleMu4_SQ_EG9er2p5, &b_b5_L1_DoubleMu4_SQ_EG9er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleMu4_SQ_OS", &b5_L1_DoubleMu4_SQ_OS, &b_b5_L1_DoubleMu4_SQ_OS);
   fChain->SetBranchAddress("b5_L1_DoubleMu4_SQ_OS_dR_Max1p2", &b5_L1_DoubleMu4_SQ_OS_dR_Max1p2, &b_b5_L1_DoubleMu4_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("b5_L1_DoubleMu4p5_SQ_OS", &b5_L1_DoubleMu4p5_SQ_OS, &b_b5_L1_DoubleMu4p5_SQ_OS);
   fChain->SetBranchAddress("b5_L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &b5_L1_DoubleMu4p5_SQ_OS_dR_Max1p2, &b_b5_L1_DoubleMu4p5_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("b5_L1_DoubleMu4p5er2p0_SQ_OS", &b5_L1_DoubleMu4p5er2p0_SQ_OS, &b_b5_L1_DoubleMu4p5er2p0_SQ_OS);
   fChain->SetBranchAddress("b5_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &b5_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18, &b_b5_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18);
   fChain->SetBranchAddress("b5_L1_DoubleMu5Upsilon_OS_DoubleEG3", &b5_L1_DoubleMu5Upsilon_OS_DoubleEG3, &b_b5_L1_DoubleMu5Upsilon_OS_DoubleEG3);
   fChain->SetBranchAddress("b5_L1_DoubleMu5_SQ_EG9er2p5", &b5_L1_DoubleMu5_SQ_EG9er2p5, &b_b5_L1_DoubleMu5_SQ_EG9er2p5);
   fChain->SetBranchAddress("b5_L1_DoubleMu9_SQ", &b5_L1_DoubleMu9_SQ, &b_b5_L1_DoubleMu9_SQ);
   fChain->SetBranchAddress("b5_L1_DoubleMu_12_5", &b5_L1_DoubleMu_12_5, &b_b5_L1_DoubleMu_12_5);
   fChain->SetBranchAddress("b5_L1_DoubleMu_15_5_SQ", &b5_L1_DoubleMu_15_5_SQ, &b_b5_L1_DoubleMu_15_5_SQ);
   fChain->SetBranchAddress("b5_L1_DoubleMu_15_7", &b5_L1_DoubleMu_15_7, &b_b5_L1_DoubleMu_15_7);
   fChain->SetBranchAddress("b5_L1_DoubleMu_15_7_Mass_Min1", &b5_L1_DoubleMu_15_7_Mass_Min1, &b_b5_L1_DoubleMu_15_7_Mass_Min1);
   fChain->SetBranchAddress("b5_L1_DoubleMu_15_7_SQ", &b5_L1_DoubleMu_15_7_SQ, &b_b5_L1_DoubleMu_15_7_SQ);
   fChain->SetBranchAddress("b5_L1_DoubleTau70er2p1", &b5_L1_DoubleTau70er2p1, &b_b5_L1_DoubleTau70er2p1);
   fChain->SetBranchAddress("b5_L1_ETM120", &b5_L1_ETM120, &b_b5_L1_ETM120);
   fChain->SetBranchAddress("b5_L1_ETM150", &b5_L1_ETM150, &b_b5_L1_ETM150);
   fChain->SetBranchAddress("b5_L1_ETMHF100", &b5_L1_ETMHF100, &b_b5_L1_ETMHF100);
   fChain->SetBranchAddress("b5_L1_ETMHF100_HTT60er", &b5_L1_ETMHF100_HTT60er, &b_b5_L1_ETMHF100_HTT60er);
   fChain->SetBranchAddress("b5_L1_ETMHF110", &b5_L1_ETMHF110, &b_b5_L1_ETMHF110);
   fChain->SetBranchAddress("b5_L1_ETMHF110_HTT60er", &b5_L1_ETMHF110_HTT60er, &b_b5_L1_ETMHF110_HTT60er);
   fChain->SetBranchAddress("b5_L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &b5_L1_ETMHF110_HTT60er_NotSecondBunchInTrain, &b_b5_L1_ETMHF110_HTT60er_NotSecondBunchInTrain);
   fChain->SetBranchAddress("b5_L1_ETMHF120", &b5_L1_ETMHF120, &b_b5_L1_ETMHF120);
   fChain->SetBranchAddress("b5_L1_ETMHF120_HTT60er", &b5_L1_ETMHF120_HTT60er, &b_b5_L1_ETMHF120_HTT60er);
   fChain->SetBranchAddress("b5_L1_ETMHF120_NotSecondBunchInTrain", &b5_L1_ETMHF120_NotSecondBunchInTrain, &b_b5_L1_ETMHF120_NotSecondBunchInTrain);
   fChain->SetBranchAddress("b5_L1_ETMHF130", &b5_L1_ETMHF130, &b_b5_L1_ETMHF130);
   fChain->SetBranchAddress("b5_L1_ETMHF130_HTT60er", &b5_L1_ETMHF130_HTT60er, &b_b5_L1_ETMHF130_HTT60er);
   fChain->SetBranchAddress("b5_L1_ETMHF140", &b5_L1_ETMHF140, &b_b5_L1_ETMHF140);
   fChain->SetBranchAddress("b5_L1_ETMHF150", &b5_L1_ETMHF150, &b_b5_L1_ETMHF150);
   fChain->SetBranchAddress("b5_L1_ETMHF90_HTT60er", &b5_L1_ETMHF90_HTT60er, &b_b5_L1_ETMHF90_HTT60er);
   fChain->SetBranchAddress("b5_L1_ETT1200", &b5_L1_ETT1200, &b_b5_L1_ETT1200);
   fChain->SetBranchAddress("b5_L1_ETT1600", &b5_L1_ETT1600, &b_b5_L1_ETT1600);
   fChain->SetBranchAddress("b5_L1_ETT2000", &b5_L1_ETT2000, &b_b5_L1_ETT2000);
   fChain->SetBranchAddress("b5_L1_FirstBunchAfterTrain", &b5_L1_FirstBunchAfterTrain, &b_b5_L1_FirstBunchAfterTrain);
   fChain->SetBranchAddress("b5_L1_FirstBunchBeforeTrain", &b5_L1_FirstBunchBeforeTrain, &b_b5_L1_FirstBunchBeforeTrain);
   fChain->SetBranchAddress("b5_L1_FirstBunchInTrain", &b5_L1_FirstBunchInTrain, &b_b5_L1_FirstBunchInTrain);
   fChain->SetBranchAddress("b5_L1_FirstCollisionInOrbit", &b5_L1_FirstCollisionInOrbit, &b_b5_L1_FirstCollisionInOrbit);
   fChain->SetBranchAddress("b5_L1_FirstCollisionInTrain", &b5_L1_FirstCollisionInTrain, &b_b5_L1_FirstCollisionInTrain);
   fChain->SetBranchAddress("b5_L1_HCAL_LaserMon_Trig", &b5_L1_HCAL_LaserMon_Trig, &b_b5_L1_HCAL_LaserMon_Trig);
   fChain->SetBranchAddress("b5_L1_HCAL_LaserMon_Veto", &b5_L1_HCAL_LaserMon_Veto, &b_b5_L1_HCAL_LaserMon_Veto);
   fChain->SetBranchAddress("b5_L1_HTT120er", &b5_L1_HTT120er, &b_b5_L1_HTT120er);
   fChain->SetBranchAddress("b5_L1_HTT160er", &b5_L1_HTT160er, &b_b5_L1_HTT160er);
   fChain->SetBranchAddress("b5_L1_HTT200er", &b5_L1_HTT200er, &b_b5_L1_HTT200er);
   fChain->SetBranchAddress("b5_L1_HTT255er", &b5_L1_HTT255er, &b_b5_L1_HTT255er);
   fChain->SetBranchAddress("b5_L1_HTT280er", &b5_L1_HTT280er, &b_b5_L1_HTT280er);
   fChain->SetBranchAddress("b5_L1_HTT280er_QuadJet_70_55_40_35_er2p4", &b5_L1_HTT280er_QuadJet_70_55_40_35_er2p4, &b_b5_L1_HTT280er_QuadJet_70_55_40_35_er2p4);
   fChain->SetBranchAddress("b5_L1_HTT320er", &b5_L1_HTT320er, &b_b5_L1_HTT320er);
   fChain->SetBranchAddress("b5_L1_HTT320er_QuadJet_70_55_40_40_er2p4", &b5_L1_HTT320er_QuadJet_70_55_40_40_er2p4, &b_b5_L1_HTT320er_QuadJet_70_55_40_40_er2p4);
   fChain->SetBranchAddress("b5_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &b5_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3, &b_b5_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3);
   fChain->SetBranchAddress("b5_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &b5_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3, &b_b5_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3);
   fChain->SetBranchAddress("b5_L1_HTT360er", &b5_L1_HTT360er, &b_b5_L1_HTT360er);
   fChain->SetBranchAddress("b5_L1_HTT400er", &b5_L1_HTT400er, &b_b5_L1_HTT400er);
   fChain->SetBranchAddress("b5_L1_HTT450er", &b5_L1_HTT450er, &b_b5_L1_HTT450er);
   fChain->SetBranchAddress("b5_L1_IsoEG32er2p5_Mt40", &b5_L1_IsoEG32er2p5_Mt40, &b_b5_L1_IsoEG32er2p5_Mt40);
   fChain->SetBranchAddress("b5_L1_IsoEG32er2p5_Mt44", &b5_L1_IsoEG32er2p5_Mt44, &b_b5_L1_IsoEG32er2p5_Mt44);
   fChain->SetBranchAddress("b5_L1_IsoEG32er2p5_Mt48", &b5_L1_IsoEG32er2p5_Mt48, &b_b5_L1_IsoEG32er2p5_Mt48);
   fChain->SetBranchAddress("b5_L1_IsoTau40er2p1_ETMHF100", &b5_L1_IsoTau40er2p1_ETMHF100, &b_b5_L1_IsoTau40er2p1_ETMHF100);
   fChain->SetBranchAddress("b5_L1_IsoTau40er2p1_ETMHF110", &b5_L1_IsoTau40er2p1_ETMHF110, &b_b5_L1_IsoTau40er2p1_ETMHF110);
   fChain->SetBranchAddress("b5_L1_IsoTau40er2p1_ETMHF120", &b5_L1_IsoTau40er2p1_ETMHF120, &b_b5_L1_IsoTau40er2p1_ETMHF120);
   fChain->SetBranchAddress("b5_L1_IsoTau40er2p1_ETMHF90", &b5_L1_IsoTau40er2p1_ETMHF90, &b_b5_L1_IsoTau40er2p1_ETMHF90);
   fChain->SetBranchAddress("b5_L1_IsolatedBunch", &b5_L1_IsolatedBunch, &b_b5_L1_IsolatedBunch);
   fChain->SetBranchAddress("b5_L1_LastBunchInTrain", &b5_L1_LastBunchInTrain, &b_b5_L1_LastBunchInTrain);
   fChain->SetBranchAddress("b5_L1_LastCollisionInTrain", &b5_L1_LastCollisionInTrain, &b_b5_L1_LastCollisionInTrain);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &b5_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3, &b_b5_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &b5_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3, &b_b5_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG24er2p1_HTT100er", &b5_L1_LooseIsoEG24er2p1_HTT100er, &b_b5_L1_LooseIsoEG24er2p1_HTT100er);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &b5_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3, &b_b5_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG26er2p1_HTT100er", &b5_L1_LooseIsoEG26er2p1_HTT100er, &b_b5_L1_LooseIsoEG26er2p1_HTT100er);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &b5_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3, &b_b5_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG28er2p1_HTT100er", &b5_L1_LooseIsoEG28er2p1_HTT100er, &b_b5_L1_LooseIsoEG28er2p1_HTT100er);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &b5_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3, &b_b5_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG30er2p1_HTT100er", &b5_L1_LooseIsoEG30er2p1_HTT100er, &b_b5_L1_LooseIsoEG30er2p1_HTT100er);
   fChain->SetBranchAddress("b5_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &b5_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3, &b_b5_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("b5_L1_MinimumBiasHF0_AND_BptxAND", &b5_L1_MinimumBiasHF0_AND_BptxAND, &b_b5_L1_MinimumBiasHF0_AND_BptxAND);
   fChain->SetBranchAddress("b5_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &b5_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6, &b_b5_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("b5_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &b5_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6, &b_b5_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6);
   fChain->SetBranchAddress("b5_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &b5_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6, &b_b5_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("b5_L1_Mu18er2p1_Tau24er2p1", &b5_L1_Mu18er2p1_Tau24er2p1, &b_b5_L1_Mu18er2p1_Tau24er2p1);
   fChain->SetBranchAddress("b5_L1_Mu18er2p1_Tau26er2p1", &b5_L1_Mu18er2p1_Tau26er2p1, &b_b5_L1_Mu18er2p1_Tau26er2p1);
   fChain->SetBranchAddress("b5_L1_Mu20_EG10er2p5", &b5_L1_Mu20_EG10er2p5, &b_b5_L1_Mu20_EG10er2p5);
   fChain->SetBranchAddress("b5_L1_Mu22er2p1_IsoTau32er2p1", &b5_L1_Mu22er2p1_IsoTau32er2p1, &b_b5_L1_Mu22er2p1_IsoTau32er2p1);
   fChain->SetBranchAddress("b5_L1_Mu22er2p1_IsoTau34er2p1", &b5_L1_Mu22er2p1_IsoTau34er2p1, &b_b5_L1_Mu22er2p1_IsoTau34er2p1);
   fChain->SetBranchAddress("b5_L1_Mu22er2p1_IsoTau36er2p1", &b5_L1_Mu22er2p1_IsoTau36er2p1, &b_b5_L1_Mu22er2p1_IsoTau36er2p1);
   fChain->SetBranchAddress("b5_L1_Mu22er2p1_IsoTau40er2p1", &b5_L1_Mu22er2p1_IsoTau40er2p1, &b_b5_L1_Mu22er2p1_IsoTau40er2p1);
   fChain->SetBranchAddress("b5_L1_Mu22er2p1_Tau70er2p1", &b5_L1_Mu22er2p1_Tau70er2p1, &b_b5_L1_Mu22er2p1_Tau70er2p1);
   fChain->SetBranchAddress("b5_L1_Mu3_Jet120er2p5_dR_Max0p4", &b5_L1_Mu3_Jet120er2p5_dR_Max0p4, &b_b5_L1_Mu3_Jet120er2p5_dR_Max0p4);
   fChain->SetBranchAddress("b5_L1_Mu3_Jet120er2p5_dR_Max0p8", &b5_L1_Mu3_Jet120er2p5_dR_Max0p8, &b_b5_L1_Mu3_Jet120er2p5_dR_Max0p8);
   fChain->SetBranchAddress("b5_L1_Mu3_Jet16er2p5_dR_Max0p4", &b5_L1_Mu3_Jet16er2p5_dR_Max0p4, &b_b5_L1_Mu3_Jet16er2p5_dR_Max0p4);
   fChain->SetBranchAddress("b5_L1_Mu3_Jet30er2p5", &b5_L1_Mu3_Jet30er2p5, &b_b5_L1_Mu3_Jet30er2p5);
   fChain->SetBranchAddress("b5_L1_Mu3_Jet35er2p5_dR_Max0p4", &b5_L1_Mu3_Jet35er2p5_dR_Max0p4, &b_b5_L1_Mu3_Jet35er2p5_dR_Max0p4);
   fChain->SetBranchAddress("b5_L1_Mu3_Jet60er2p5_dR_Max0p4", &b5_L1_Mu3_Jet60er2p5_dR_Max0p4, &b_b5_L1_Mu3_Jet60er2p5_dR_Max0p4);
   fChain->SetBranchAddress("b5_L1_Mu3_Jet80er2p5_dR_Max0p4", &b5_L1_Mu3_Jet80er2p5_dR_Max0p4, &b_b5_L1_Mu3_Jet80er2p5_dR_Max0p4);
   fChain->SetBranchAddress("b5_L1_Mu3er1p5_Jet100er2p5_ETMHF40", &b5_L1_Mu3er1p5_Jet100er2p5_ETMHF40, &b_b5_L1_Mu3er1p5_Jet100er2p5_ETMHF40);
   fChain->SetBranchAddress("b5_L1_Mu3er1p5_Jet100er2p5_ETMHF50", &b5_L1_Mu3er1p5_Jet100er2p5_ETMHF50, &b_b5_L1_Mu3er1p5_Jet100er2p5_ETMHF50);
   fChain->SetBranchAddress("b5_L1_Mu5_EG23er2p5", &b5_L1_Mu5_EG23er2p5, &b_b5_L1_Mu5_EG23er2p5);
   fChain->SetBranchAddress("b5_L1_Mu5_LooseIsoEG20er2p5", &b5_L1_Mu5_LooseIsoEG20er2p5, &b_b5_L1_Mu5_LooseIsoEG20er2p5);
   fChain->SetBranchAddress("b5_L1_Mu6_DoubleEG10er2p5", &b5_L1_Mu6_DoubleEG10er2p5, &b_b5_L1_Mu6_DoubleEG10er2p5);
   fChain->SetBranchAddress("b5_L1_Mu6_DoubleEG12er2p5", &b5_L1_Mu6_DoubleEG12er2p5, &b_b5_L1_Mu6_DoubleEG12er2p5);
   fChain->SetBranchAddress("b5_L1_Mu6_DoubleEG15er2p5", &b5_L1_Mu6_DoubleEG15er2p5, &b_b5_L1_Mu6_DoubleEG15er2p5);
   fChain->SetBranchAddress("b5_L1_Mu6_DoubleEG17er2p5", &b5_L1_Mu6_DoubleEG17er2p5, &b_b5_L1_Mu6_DoubleEG17er2p5);
   fChain->SetBranchAddress("b5_L1_Mu6_HTT240er", &b5_L1_Mu6_HTT240er, &b_b5_L1_Mu6_HTT240er);
   fChain->SetBranchAddress("b5_L1_Mu6_HTT250er", &b5_L1_Mu6_HTT250er, &b_b5_L1_Mu6_HTT250er);
   fChain->SetBranchAddress("b5_L1_Mu7_EG23er2p5", &b5_L1_Mu7_EG23er2p5, &b_b5_L1_Mu7_EG23er2p5);
   fChain->SetBranchAddress("b5_L1_Mu7_LooseIsoEG20er2p5", &b5_L1_Mu7_LooseIsoEG20er2p5, &b_b5_L1_Mu7_LooseIsoEG20er2p5);
   fChain->SetBranchAddress("b5_L1_Mu7_LooseIsoEG23er2p5", &b5_L1_Mu7_LooseIsoEG23er2p5, &b_b5_L1_Mu7_LooseIsoEG23er2p5);
   fChain->SetBranchAddress("b5_L1_NotBptxOR", &b5_L1_NotBptxOR, &b_b5_L1_NotBptxOR);
   fChain->SetBranchAddress("b5_L1_QuadJet36er2p5_IsoTau52er2p1", &b5_L1_QuadJet36er2p5_IsoTau52er2p1, &b_b5_L1_QuadJet36er2p5_IsoTau52er2p1);
   fChain->SetBranchAddress("b5_L1_QuadJet60er2p5", &b5_L1_QuadJet60er2p5, &b_b5_L1_QuadJet60er2p5);
   fChain->SetBranchAddress("b5_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &b5_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0, &b_b5_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0);
   fChain->SetBranchAddress("b5_L1_QuadMu0", &b5_L1_QuadMu0, &b_b5_L1_QuadMu0);
   fChain->SetBranchAddress("b5_L1_QuadMu0_OQ", &b5_L1_QuadMu0_OQ, &b_b5_L1_QuadMu0_OQ);
   fChain->SetBranchAddress("b5_L1_QuadMu0_SQ", &b5_L1_QuadMu0_SQ, &b_b5_L1_QuadMu0_SQ);
   fChain->SetBranchAddress("b5_L1_SecondBunchInTrain", &b5_L1_SecondBunchInTrain, &b_b5_L1_SecondBunchInTrain);
   fChain->SetBranchAddress("b5_L1_SecondLastBunchInTrain", &b5_L1_SecondLastBunchInTrain, &b_b5_L1_SecondLastBunchInTrain);
   fChain->SetBranchAddress("b5_L1_SingleEG10er2p5", &b5_L1_SingleEG10er2p5, &b_b5_L1_SingleEG10er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG15er2p5", &b5_L1_SingleEG15er2p5, &b_b5_L1_SingleEG15er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG26er2p5", &b5_L1_SingleEG26er2p5, &b_b5_L1_SingleEG26er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG34er2p5", &b5_L1_SingleEG34er2p5, &b_b5_L1_SingleEG34er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG36er2p5", &b5_L1_SingleEG36er2p5, &b_b5_L1_SingleEG36er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG38er2p5", &b5_L1_SingleEG38er2p5, &b_b5_L1_SingleEG38er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG40er2p5", &b5_L1_SingleEG40er2p5, &b_b5_L1_SingleEG40er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG42er2p5", &b5_L1_SingleEG42er2p5, &b_b5_L1_SingleEG42er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG45er2p5", &b5_L1_SingleEG45er2p5, &b_b5_L1_SingleEG45er2p5);
   fChain->SetBranchAddress("b5_L1_SingleEG50", &b5_L1_SingleEG50, &b_b5_L1_SingleEG50);
   fChain->SetBranchAddress("b5_L1_SingleEG60", &b5_L1_SingleEG60, &b_b5_L1_SingleEG60);
   fChain->SetBranchAddress("b5_L1_SingleEG8er2p5", &b5_L1_SingleEG8er2p5, &b_b5_L1_SingleEG8er2p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG24er1p5", &b5_L1_SingleIsoEG24er1p5, &b_b5_L1_SingleIsoEG24er1p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG24er2p1", &b5_L1_SingleIsoEG24er2p1, &b_b5_L1_SingleIsoEG24er2p1);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG26er1p5", &b5_L1_SingleIsoEG26er1p5, &b_b5_L1_SingleIsoEG26er1p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG26er2p1", &b5_L1_SingleIsoEG26er2p1, &b_b5_L1_SingleIsoEG26er2p1);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG26er2p5", &b5_L1_SingleIsoEG26er2p5, &b_b5_L1_SingleIsoEG26er2p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG28er1p5", &b5_L1_SingleIsoEG28er1p5, &b_b5_L1_SingleIsoEG28er1p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG28er2p1", &b5_L1_SingleIsoEG28er2p1, &b_b5_L1_SingleIsoEG28er2p1);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG28er2p5", &b5_L1_SingleIsoEG28er2p5, &b_b5_L1_SingleIsoEG28er2p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG30er2p1", &b5_L1_SingleIsoEG30er2p1, &b_b5_L1_SingleIsoEG30er2p1);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG30er2p5", &b5_L1_SingleIsoEG30er2p5, &b_b5_L1_SingleIsoEG30er2p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG32er2p1", &b5_L1_SingleIsoEG32er2p1, &b_b5_L1_SingleIsoEG32er2p1);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG32er2p5", &b5_L1_SingleIsoEG32er2p5, &b_b5_L1_SingleIsoEG32er2p5);
   fChain->SetBranchAddress("b5_L1_SingleIsoEG34er2p5", &b5_L1_SingleIsoEG34er2p5, &b_b5_L1_SingleIsoEG34er2p5);
   fChain->SetBranchAddress("b5_L1_SingleJet10erHE", &b5_L1_SingleJet10erHE, &b_b5_L1_SingleJet10erHE);
   fChain->SetBranchAddress("b5_L1_SingleJet120", &b5_L1_SingleJet120, &b_b5_L1_SingleJet120);
   fChain->SetBranchAddress("b5_L1_SingleJet120_FWD3p0", &b5_L1_SingleJet120_FWD3p0, &b_b5_L1_SingleJet120_FWD3p0);
   fChain->SetBranchAddress("b5_L1_SingleJet120er2p5", &b5_L1_SingleJet120er2p5, &b_b5_L1_SingleJet120er2p5);
   fChain->SetBranchAddress("b5_L1_SingleJet12erHE", &b5_L1_SingleJet12erHE, &b_b5_L1_SingleJet12erHE);
   fChain->SetBranchAddress("b5_L1_SingleJet140er2p5", &b5_L1_SingleJet140er2p5, &b_b5_L1_SingleJet140er2p5);
   fChain->SetBranchAddress("b5_L1_SingleJet140er2p5_ETMHF80", &b5_L1_SingleJet140er2p5_ETMHF80, &b_b5_L1_SingleJet140er2p5_ETMHF80);
   fChain->SetBranchAddress("b5_L1_SingleJet140er2p5_ETMHF90", &b5_L1_SingleJet140er2p5_ETMHF90, &b_b5_L1_SingleJet140er2p5_ETMHF90);
   fChain->SetBranchAddress("b5_L1_SingleJet160er2p5", &b5_L1_SingleJet160er2p5, &b_b5_L1_SingleJet160er2p5);
   fChain->SetBranchAddress("b5_L1_SingleJet180", &b5_L1_SingleJet180, &b_b5_L1_SingleJet180);
   fChain->SetBranchAddress("b5_L1_SingleJet180er2p5", &b5_L1_SingleJet180er2p5, &b_b5_L1_SingleJet180er2p5);
   fChain->SetBranchAddress("b5_L1_SingleJet200", &b5_L1_SingleJet200, &b_b5_L1_SingleJet200);
   fChain->SetBranchAddress("b5_L1_SingleJet20er2p5_NotBptxOR", &b5_L1_SingleJet20er2p5_NotBptxOR, &b_b5_L1_SingleJet20er2p5_NotBptxOR);
   fChain->SetBranchAddress("b5_L1_SingleJet20er2p5_NotBptxOR_3BX", &b5_L1_SingleJet20er2p5_NotBptxOR_3BX, &b_b5_L1_SingleJet20er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("b5_L1_SingleJet35", &b5_L1_SingleJet35, &b_b5_L1_SingleJet35);
   fChain->SetBranchAddress("b5_L1_SingleJet35_FWD3p0", &b5_L1_SingleJet35_FWD3p0, &b_b5_L1_SingleJet35_FWD3p0);
   fChain->SetBranchAddress("b5_L1_SingleJet35er2p5", &b5_L1_SingleJet35er2p5, &b_b5_L1_SingleJet35er2p5);
   fChain->SetBranchAddress("b5_L1_SingleJet43er2p5_NotBptxOR_3BX", &b5_L1_SingleJet43er2p5_NotBptxOR_3BX, &b_b5_L1_SingleJet43er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("b5_L1_SingleJet46er2p5_NotBptxOR_3BX", &b5_L1_SingleJet46er2p5_NotBptxOR_3BX, &b_b5_L1_SingleJet46er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("b5_L1_SingleJet60", &b5_L1_SingleJet60, &b_b5_L1_SingleJet60);
   fChain->SetBranchAddress("b5_L1_SingleJet60_FWD3p0", &b5_L1_SingleJet60_FWD3p0, &b_b5_L1_SingleJet60_FWD3p0);
   fChain->SetBranchAddress("b5_L1_SingleJet60er2p5", &b5_L1_SingleJet60er2p5, &b_b5_L1_SingleJet60er2p5);
   fChain->SetBranchAddress("b5_L1_SingleJet8erHE", &b5_L1_SingleJet8erHE, &b_b5_L1_SingleJet8erHE);
   fChain->SetBranchAddress("b5_L1_SingleJet90", &b5_L1_SingleJet90, &b_b5_L1_SingleJet90);
   fChain->SetBranchAddress("b5_L1_SingleJet90_FWD3p0", &b5_L1_SingleJet90_FWD3p0, &b_b5_L1_SingleJet90_FWD3p0);
   fChain->SetBranchAddress("b5_L1_SingleJet90er2p5", &b5_L1_SingleJet90er2p5, &b_b5_L1_SingleJet90er2p5);
   fChain->SetBranchAddress("b5_L1_SingleLooseIsoEG28er1p5", &b5_L1_SingleLooseIsoEG28er1p5, &b_b5_L1_SingleLooseIsoEG28er1p5);
   fChain->SetBranchAddress("b5_L1_SingleLooseIsoEG30er1p5", &b5_L1_SingleLooseIsoEG30er1p5, &b_b5_L1_SingleLooseIsoEG30er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu0_BMTF", &b5_L1_SingleMu0_BMTF, &b_b5_L1_SingleMu0_BMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu0_DQ", &b5_L1_SingleMu0_DQ, &b_b5_L1_SingleMu0_DQ);
   fChain->SetBranchAddress("b5_L1_SingleMu0_EMTF", &b5_L1_SingleMu0_EMTF, &b_b5_L1_SingleMu0_EMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu0_OMTF", &b5_L1_SingleMu0_OMTF, &b_b5_L1_SingleMu0_OMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu10er1p5", &b5_L1_SingleMu10er1p5, &b_b5_L1_SingleMu10er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu12_DQ_BMTF", &b5_L1_SingleMu12_DQ_BMTF, &b_b5_L1_SingleMu12_DQ_BMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu12_DQ_EMTF", &b5_L1_SingleMu12_DQ_EMTF, &b_b5_L1_SingleMu12_DQ_EMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu12_DQ_OMTF", &b5_L1_SingleMu12_DQ_OMTF, &b_b5_L1_SingleMu12_DQ_OMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu12er1p5", &b5_L1_SingleMu12er1p5, &b_b5_L1_SingleMu12er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu14er1p5", &b5_L1_SingleMu14er1p5, &b_b5_L1_SingleMu14er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu15_DQ", &b5_L1_SingleMu15_DQ, &b_b5_L1_SingleMu15_DQ);
   fChain->SetBranchAddress("b5_L1_SingleMu16er1p5", &b5_L1_SingleMu16er1p5, &b_b5_L1_SingleMu16er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu18", &b5_L1_SingleMu18, &b_b5_L1_SingleMu18);
   fChain->SetBranchAddress("b5_L1_SingleMu18er1p5", &b5_L1_SingleMu18er1p5, &b_b5_L1_SingleMu18er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu20", &b5_L1_SingleMu20, &b_b5_L1_SingleMu20);
   fChain->SetBranchAddress("b5_L1_SingleMu22", &b5_L1_SingleMu22, &b_b5_L1_SingleMu22);
   fChain->SetBranchAddress("b5_L1_SingleMu22_BMTF", &b5_L1_SingleMu22_BMTF, &b_b5_L1_SingleMu22_BMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu22_EMTF", &b5_L1_SingleMu22_EMTF, &b_b5_L1_SingleMu22_EMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu22_OMTF", &b5_L1_SingleMu22_OMTF, &b_b5_L1_SingleMu22_OMTF);
   fChain->SetBranchAddress("b5_L1_SingleMu25", &b5_L1_SingleMu25, &b_b5_L1_SingleMu25);
   fChain->SetBranchAddress("b5_L1_SingleMu3", &b5_L1_SingleMu3, &b_b5_L1_SingleMu3);
   fChain->SetBranchAddress("b5_L1_SingleMu5", &b5_L1_SingleMu5, &b_b5_L1_SingleMu5);
   fChain->SetBranchAddress("b5_L1_SingleMu6er1p5", &b5_L1_SingleMu6er1p5, &b_b5_L1_SingleMu6er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu7", &b5_L1_SingleMu7, &b_b5_L1_SingleMu7);
   fChain->SetBranchAddress("b5_L1_SingleMu7_DQ", &b5_L1_SingleMu7_DQ, &b_b5_L1_SingleMu7_DQ);
   fChain->SetBranchAddress("b5_L1_SingleMu7er1p5", &b5_L1_SingleMu7er1p5, &b_b5_L1_SingleMu7er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu8er1p5", &b5_L1_SingleMu8er1p5, &b_b5_L1_SingleMu8er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMu9er1p5", &b5_L1_SingleMu9er1p5, &b_b5_L1_SingleMu9er1p5);
   fChain->SetBranchAddress("b5_L1_SingleMuCosmics", &b5_L1_SingleMuCosmics, &b_b5_L1_SingleMuCosmics);
   fChain->SetBranchAddress("b5_L1_SingleMuCosmics_BMTF", &b5_L1_SingleMuCosmics_BMTF, &b_b5_L1_SingleMuCosmics_BMTF);
   fChain->SetBranchAddress("b5_L1_SingleMuCosmics_EMTF", &b5_L1_SingleMuCosmics_EMTF, &b_b5_L1_SingleMuCosmics_EMTF);
   fChain->SetBranchAddress("b5_L1_SingleMuCosmics_OMTF", &b5_L1_SingleMuCosmics_OMTF, &b_b5_L1_SingleMuCosmics_OMTF);
   fChain->SetBranchAddress("b5_L1_SingleMuOpen", &b5_L1_SingleMuOpen, &b_b5_L1_SingleMuOpen);
   fChain->SetBranchAddress("b5_L1_SingleMuOpen_NotBptxOR", &b5_L1_SingleMuOpen_NotBptxOR, &b_b5_L1_SingleMuOpen_NotBptxOR);
   fChain->SetBranchAddress("b5_L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &b5_L1_SingleMuOpen_er1p1_NotBptxOR_3BX, &b_b5_L1_SingleMuOpen_er1p1_NotBptxOR_3BX);
   fChain->SetBranchAddress("b5_L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &b5_L1_SingleMuOpen_er1p4_NotBptxOR_3BX, &b_b5_L1_SingleMuOpen_er1p4_NotBptxOR_3BX);
   fChain->SetBranchAddress("b5_L1_SingleTau120er2p1", &b5_L1_SingleTau120er2p1, &b_b5_L1_SingleTau120er2p1);
   fChain->SetBranchAddress("b5_L1_SingleTau130er2p1", &b5_L1_SingleTau130er2p1, &b_b5_L1_SingleTau130er2p1);
   fChain->SetBranchAddress("b5_L1_TOTEM_1", &b5_L1_TOTEM_1, &b_b5_L1_TOTEM_1);
   fChain->SetBranchAddress("b5_L1_TOTEM_2", &b5_L1_TOTEM_2, &b_b5_L1_TOTEM_2);
   fChain->SetBranchAddress("b5_L1_TOTEM_3", &b5_L1_TOTEM_3, &b_b5_L1_TOTEM_3);
   fChain->SetBranchAddress("b5_L1_TOTEM_4", &b5_L1_TOTEM_4, &b_b5_L1_TOTEM_4);
   fChain->SetBranchAddress("b5_L1_TripleEG16er2p5", &b5_L1_TripleEG16er2p5, &b_b5_L1_TripleEG16er2p5);
   fChain->SetBranchAddress("b5_L1_TripleEG_16_12_8_er2p5", &b5_L1_TripleEG_16_12_8_er2p5, &b_b5_L1_TripleEG_16_12_8_er2p5);
   fChain->SetBranchAddress("b5_L1_TripleEG_16_15_8_er2p5", &b5_L1_TripleEG_16_15_8_er2p5, &b_b5_L1_TripleEG_16_15_8_er2p5);
   fChain->SetBranchAddress("b5_L1_TripleEG_18_17_8_er2p5", &b5_L1_TripleEG_18_17_8_er2p5, &b_b5_L1_TripleEG_18_17_8_er2p5);
   fChain->SetBranchAddress("b5_L1_TripleEG_18_18_12_er2p5", &b5_L1_TripleEG_18_18_12_er2p5, &b_b5_L1_TripleEG_18_18_12_er2p5);
   fChain->SetBranchAddress("b5_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &b5_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5, &b_b5_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5);
   fChain->SetBranchAddress("b5_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &b5_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5, &b_b5_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5);
   fChain->SetBranchAddress("b5_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &b5_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5, &b_b5_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5);
   fChain->SetBranchAddress("b5_L1_TripleMu0", &b5_L1_TripleMu0, &b_b5_L1_TripleMu0);
   fChain->SetBranchAddress("b5_L1_TripleMu0_OQ", &b5_L1_TripleMu0_OQ, &b_b5_L1_TripleMu0_OQ);
   fChain->SetBranchAddress("b5_L1_TripleMu0_SQ", &b5_L1_TripleMu0_SQ, &b_b5_L1_TripleMu0_SQ);
   fChain->SetBranchAddress("b5_L1_TripleMu3", &b5_L1_TripleMu3, &b_b5_L1_TripleMu3);
   fChain->SetBranchAddress("b5_L1_TripleMu3_SQ", &b5_L1_TripleMu3_SQ, &b_b5_L1_TripleMu3_SQ);
   fChain->SetBranchAddress("b5_L1_TripleMu_5SQ_3SQ_0OQ", &b5_L1_TripleMu_5SQ_3SQ_0OQ, &b_b5_L1_TripleMu_5SQ_3SQ_0OQ);
   fChain->SetBranchAddress("b5_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &b5_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_b5_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("b5_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &b5_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_b5_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("b5_L1_TripleMu_5_3_3", &b5_L1_TripleMu_5_3_3, &b_b5_L1_TripleMu_5_3_3);
   fChain->SetBranchAddress("b5_L1_TripleMu_5_3_3_SQ", &b5_L1_TripleMu_5_3_3_SQ, &b_b5_L1_TripleMu_5_3_3_SQ);
   fChain->SetBranchAddress("b5_L1_TripleMu_5_3p5_2p5", &b5_L1_TripleMu_5_3p5_2p5, &b_b5_L1_TripleMu_5_3p5_2p5);
   fChain->SetBranchAddress("b5_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &b5_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_b5_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("b5_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &b5_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17, &b_b5_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17);
   fChain->SetBranchAddress("b5_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &b5_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_b5_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("b5_L1_TripleMu_5_5_3", &b5_L1_TripleMu_5_5_3, &b_b5_L1_TripleMu_5_5_3);
   fChain->SetBranchAddress("b5_L1_UnpairedBunchBptxMinus", &b5_L1_UnpairedBunchBptxMinus, &b_b5_L1_UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("b5_L1_UnpairedBunchBptxPlus", &b5_L1_UnpairedBunchBptxPlus, &b_b5_L1_UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("b5_L1_ZeroBias", &b5_L1_ZeroBias, &b_b5_L1_ZeroBias);
   fChain->SetBranchAddress("b5_L1_ZeroBias_copy", &b5_L1_ZeroBias_copy, &b_b5_L1_ZeroBias_copy);
   fChain->SetBranchAddress("b5_L1_UnprefireableEvent", &b5_L1_UnprefireableEvent, &b_b5_L1_UnprefireableEvent);
   fChain->SetBranchAddress("b5_Flag_HBHENoiseFilter", &b5_Flag_HBHENoiseFilter, &b_b5_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("b5_Flag_HBHENoiseIsoFilter", &b5_Flag_HBHENoiseIsoFilter, &b_b5_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("b5_Flag_CSCTightHaloFilter", &b5_Flag_CSCTightHaloFilter, &b_b5_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("b5_Flag_CSCTightHaloTrkMuUnvetoFilter", &b5_Flag_CSCTightHaloTrkMuUnvetoFilter, &b_b5_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("b5_Flag_CSCTightHalo2015Filter", &b5_Flag_CSCTightHalo2015Filter, &b_b5_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("b5_Flag_globalTightHalo2016Filter", &b5_Flag_globalTightHalo2016Filter, &b_b5_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("b5_Flag_globalSuperTightHalo2016Filter", &b5_Flag_globalSuperTightHalo2016Filter, &b_b5_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("b5_Flag_HcalStripHaloFilter", &b5_Flag_HcalStripHaloFilter, &b_b5_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("b5_Flag_hcalLaserEventFilter", &b5_Flag_hcalLaserEventFilter, &b_b5_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("b5_Flag_EcalDeadCellTriggerPrimitiveFilter", &b5_Flag_EcalDeadCellTriggerPrimitiveFilter, &b_b5_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("b5_Flag_EcalDeadCellBoundaryEnergyFilter", &b5_Flag_EcalDeadCellBoundaryEnergyFilter, &b_b5_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("b5_Flag_ecalBadCalibFilter", &b5_Flag_ecalBadCalibFilter, &b_b5_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("b5_Flag_goodVertices", &b5_Flag_goodVertices, &b_b5_Flag_goodVertices);
   fChain->SetBranchAddress("b5_Flag_eeBadScFilter", &b5_Flag_eeBadScFilter, &b_b5_Flag_eeBadScFilter);
   fChain->SetBranchAddress("b5_Flag_ecalLaserCorrFilter", &b5_Flag_ecalLaserCorrFilter, &b_b5_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("b5_Flag_trkPOGFilters", &b5_Flag_trkPOGFilters, &b_b5_Flag_trkPOGFilters);
   fChain->SetBranchAddress("b5_Flag_chargedHadronTrackResolutionFilter", &b5_Flag_chargedHadronTrackResolutionFilter, &b_b5_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("b5_Flag_muonBadTrackFilter", &b5_Flag_muonBadTrackFilter, &b_b5_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("b5_Flag_BadChargedCandidateFilter", &b5_Flag_BadChargedCandidateFilter, &b_b5_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("b5_Flag_BadPFMuonFilter", &b5_Flag_BadPFMuonFilter, &b_b5_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("b5_Flag_BadPFMuonDzFilter", &b5_Flag_BadPFMuonDzFilter, &b_b5_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("b5_Flag_hfNoisyHitsFilter", &b5_Flag_hfNoisyHitsFilter, &b_b5_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("b5_Flag_BadChargedCandidateSummer16Filter", &b5_Flag_BadChargedCandidateSummer16Filter, &b_b5_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("b5_Flag_BadPFMuonSummer16Filter", &b5_Flag_BadPFMuonSummer16Filter, &b_b5_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("b5_Flag_trkPOG_manystripclus53X", &b5_Flag_trkPOG_manystripclus53X, &b_b5_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("b5_Flag_trkPOG_toomanystripclus53X", &b5_Flag_trkPOG_toomanystripclus53X, &b_b5_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("b5_Flag_trkPOG_logErrorTooManyClusters", &b5_Flag_trkPOG_logErrorTooManyClusters, &b_b5_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("b5_Flag_METFilters", &b5_Flag_METFilters, &b_b5_Flag_METFilters);
   fChain->SetBranchAddress("b5_L1Reco_step", &b5_L1Reco_step, &b_b5_L1Reco_step);
   fChain->SetBranchAddress("b5_HLTriggerFirstPath", &b5_HLTriggerFirstPath, &b_b5_HLTriggerFirstPath);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet360_TrimMass30", &b5_HLT_AK8PFJet360_TrimMass30, &b_b5_HLT_AK8PFJet360_TrimMass30);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet380_TrimMass30", &b5_HLT_AK8PFJet380_TrimMass30, &b_b5_HLT_AK8PFJet380_TrimMass30);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet400_TrimMass30", &b5_HLT_AK8PFJet400_TrimMass30, &b_b5_HLT_AK8PFJet400_TrimMass30);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet420_TrimMass30", &b5_HLT_AK8PFJet420_TrimMass30, &b_b5_HLT_AK8PFJet420_TrimMass30);
   fChain->SetBranchAddress("b5_HLT_AK8PFHT750_TrimMass50", &b5_HLT_AK8PFHT750_TrimMass50, &b_b5_HLT_AK8PFHT750_TrimMass50);
   fChain->SetBranchAddress("b5_HLT_AK8PFHT800_TrimMass50", &b5_HLT_AK8PFHT800_TrimMass50, &b_b5_HLT_AK8PFHT800_TrimMass50);
   fChain->SetBranchAddress("b5_HLT_AK8PFHT850_TrimMass50", &b5_HLT_AK8PFHT850_TrimMass50, &b_b5_HLT_AK8PFHT850_TrimMass50);
   fChain->SetBranchAddress("b5_HLT_AK8PFHT900_TrimMass50", &b5_HLT_AK8PFHT900_TrimMass50, &b_b5_HLT_AK8PFHT900_TrimMass50);
   fChain->SetBranchAddress("b5_HLT_CaloJet500_NoJetID", &b5_HLT_CaloJet500_NoJetID, &b_b5_HLT_CaloJet500_NoJetID);
   fChain->SetBranchAddress("b5_HLT_CaloJet550_NoJetID", &b5_HLT_CaloJet550_NoJetID, &b_b5_HLT_CaloJet550_NoJetID);
   fChain->SetBranchAddress("b5_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL", &b5_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL, &b_b5_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon", &b5_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon, &b_b5_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon);
   fChain->SetBranchAddress("b5_HLT_Trimuon5_3p5_2_Upsilon_Muon", &b5_HLT_Trimuon5_3p5_2_Upsilon_Muon, &b_b5_HLT_Trimuon5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("b5_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon", &b5_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon, &b_b5_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("b5_HLT_DoubleEle25_CaloIdL_MW", &b5_HLT_DoubleEle25_CaloIdL_MW, &b_b5_HLT_DoubleEle25_CaloIdL_MW);
   fChain->SetBranchAddress("b5_HLT_DoubleEle27_CaloIdL_MW", &b5_HLT_DoubleEle27_CaloIdL_MW, &b_b5_HLT_DoubleEle27_CaloIdL_MW);
   fChain->SetBranchAddress("b5_HLT_DoubleEle33_CaloIdL_MW", &b5_HLT_DoubleEle33_CaloIdL_MW, &b_b5_HLT_DoubleEle33_CaloIdL_MW);
   fChain->SetBranchAddress("b5_HLT_DoubleEle24_eta2p1_WPTight_Gsf", &b5_HLT_DoubleEle24_eta2p1_WPTight_Gsf, &b_b5_HLT_DoubleEle24_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350, &b_b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350);
   fChain->SetBranchAddress("b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("b5_HLT_Ele27_Ele37_CaloIdL_MW", &b5_HLT_Ele27_Ele37_CaloIdL_MW, &b_b5_HLT_Ele27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("b5_HLT_Mu27_Ele37_CaloIdL_MW", &b5_HLT_Mu27_Ele37_CaloIdL_MW, &b_b5_HLT_Mu27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("b5_HLT_Mu37_Ele27_CaloIdL_MW", &b5_HLT_Mu37_Ele27_CaloIdL_MW, &b_b5_HLT_Mu37_Ele27_CaloIdL_MW);
   fChain->SetBranchAddress("b5_HLT_Mu37_TkMu27", &b5_HLT_Mu37_TkMu27, &b_b5_HLT_Mu37_TkMu27);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_3_Bs", &b5_HLT_DoubleMu4_3_Bs, &b_b5_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_3_Jpsi", &b5_HLT_DoubleMu4_3_Jpsi, &b_b5_HLT_DoubleMu4_3_Jpsi);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_JpsiTrk_Displaced", &b5_HLT_DoubleMu4_JpsiTrk_Displaced, &b_b5_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced, &b_b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_Trk_Tau3mu", &b5_HLT_DoubleMu3_Trk_Tau3mu, &b_b5_HLT_DoubleMu3_Trk_Tau3mu);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_TkMu_DsTau3Mu", &b5_HLT_DoubleMu3_TkMu_DsTau3Mu, &b_b5_HLT_DoubleMu3_TkMu_DsTau3Mu);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_PsiPrimeTrk_Displaced", &b5_HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_b5_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_Mass3p8_DZ_PFHT350", &b5_HLT_DoubleMu4_Mass3p8_DZ_PFHT350, &b_b5_HLT_DoubleMu4_Mass3p8_DZ_PFHT350);
   fChain->SetBranchAddress("b5_HLT_Mu3_PFJet40", &b5_HLT_Mu3_PFJet40, &b_b5_HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_L2Mu2_Jpsi", &b5_HLT_Mu7p5_L2Mu2_Jpsi, &b_b5_HLT_Mu7p5_L2Mu2_Jpsi);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_L2Mu2_Upsilon", &b5_HLT_Mu7p5_L2Mu2_Upsilon, &b_b5_HLT_Mu7p5_L2Mu2_Upsilon);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_Track2_Jpsi", &b5_HLT_Mu7p5_Track2_Jpsi, &b_b5_HLT_Mu7p5_Track2_Jpsi);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_Track3p5_Jpsi", &b5_HLT_Mu7p5_Track3p5_Jpsi, &b_b5_HLT_Mu7p5_Track3p5_Jpsi);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_Track7_Jpsi", &b5_HLT_Mu7p5_Track7_Jpsi, &b_b5_HLT_Mu7p5_Track7_Jpsi);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_Track2_Upsilon", &b5_HLT_Mu7p5_Track2_Upsilon, &b_b5_HLT_Mu7p5_Track2_Upsilon);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_Track3p5_Upsilon", &b5_HLT_Mu7p5_Track3p5_Upsilon, &b_b5_HLT_Mu7p5_Track3p5_Upsilon);
   fChain->SetBranchAddress("b5_HLT_Mu7p5_Track7_Upsilon", &b5_HLT_Mu7p5_Track7_Upsilon, &b_b5_HLT_Mu7p5_Track7_Upsilon);
   fChain->SetBranchAddress("b5_HLT_DoublePhoton33_CaloIdL", &b5_HLT_DoublePhoton33_CaloIdL, &b_b5_HLT_DoublePhoton33_CaloIdL);
   fChain->SetBranchAddress("b5_HLT_DoublePhoton70", &b5_HLT_DoublePhoton70, &b_b5_HLT_DoublePhoton70);
   fChain->SetBranchAddress("b5_HLT_DoublePhoton85", &b5_HLT_DoublePhoton85, &b_b5_HLT_DoublePhoton85);
   fChain->SetBranchAddress("b5_HLT_Ele20_WPTight_Gsf", &b5_HLT_Ele20_WPTight_Gsf, &b_b5_HLT_Ele20_WPTight_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele15_WPLoose_Gsf", &b5_HLT_Ele15_WPLoose_Gsf, &b_b5_HLT_Ele15_WPLoose_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele17_WPLoose_Gsf", &b5_HLT_Ele17_WPLoose_Gsf, &b_b5_HLT_Ele17_WPLoose_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele20_WPLoose_Gsf", &b5_HLT_Ele20_WPLoose_Gsf, &b_b5_HLT_Ele20_WPLoose_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele20_eta2p1_WPLoose_Gsf", &b5_HLT_Ele20_eta2p1_WPLoose_Gsf, &b_b5_HLT_Ele20_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("b5_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &b5_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG, &b_b5_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
   fChain->SetBranchAddress("b5_HLT_Ele27_WPTight_Gsf", &b5_HLT_Ele27_WPTight_Gsf, &b_b5_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele32_WPTight_Gsf", &b5_HLT_Ele32_WPTight_Gsf, &b_b5_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele35_WPTight_Gsf", &b5_HLT_Ele35_WPTight_Gsf, &b_b5_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele35_WPTight_Gsf_L1EGMT", &b5_HLT_Ele35_WPTight_Gsf_L1EGMT, &b_b5_HLT_Ele35_WPTight_Gsf_L1EGMT);
   fChain->SetBranchAddress("b5_HLT_Ele38_WPTight_Gsf", &b5_HLT_Ele38_WPTight_Gsf, &b_b5_HLT_Ele38_WPTight_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele40_WPTight_Gsf", &b5_HLT_Ele40_WPTight_Gsf, &b_b5_HLT_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("b5_HLT_Ele32_WPTight_Gsf_L1DoubleEG", &b5_HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_b5_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   fChain->SetBranchAddress("b5_HLT_HT450_Beamspot", &b5_HLT_HT450_Beamspot, &b_b5_HLT_HT450_Beamspot);
   fChain->SetBranchAddress("b5_HLT_HT300_Beamspot", &b5_HLT_HT300_Beamspot, &b_b5_HLT_HT300_Beamspot);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_Beamspot", &b5_HLT_ZeroBias_Beamspot, &b_b5_HLT_ZeroBias_Beamspot);
   fChain->SetBranchAddress("b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1", &b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1, &b_b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1", &b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1, &b_b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1", &b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1, &b_b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1", &b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1", &b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1", &b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1, &b_b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu20", &b5_HLT_IsoMu20, &b_b5_HLT_IsoMu20);
   fChain->SetBranchAddress("b5_HLT_IsoMu24", &b5_HLT_IsoMu24, &b_b5_HLT_IsoMu24);
   fChain->SetBranchAddress("b5_HLT_IsoMu24_eta2p1", &b5_HLT_IsoMu24_eta2p1, &b_b5_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("b5_HLT_IsoMu27", &b5_HLT_IsoMu27, &b_b5_HLT_IsoMu27);
   fChain->SetBranchAddress("b5_HLT_IsoMu30", &b5_HLT_IsoMu30, &b_b5_HLT_IsoMu30);
   fChain->SetBranchAddress("b5_HLT_UncorrectedJetE30_NoBPTX", &b5_HLT_UncorrectedJetE30_NoBPTX, &b_b5_HLT_UncorrectedJetE30_NoBPTX);
   fChain->SetBranchAddress("b5_HLT_UncorrectedJetE30_NoBPTX3BX", &b5_HLT_UncorrectedJetE30_NoBPTX3BX, &b_b5_HLT_UncorrectedJetE30_NoBPTX3BX);
   fChain->SetBranchAddress("b5_HLT_UncorrectedJetE60_NoBPTX3BX", &b5_HLT_UncorrectedJetE60_NoBPTX3BX, &b_b5_HLT_UncorrectedJetE60_NoBPTX3BX);
   fChain->SetBranchAddress("b5_HLT_UncorrectedJetE70_NoBPTX3BX", &b5_HLT_UncorrectedJetE70_NoBPTX3BX, &b_b5_HLT_UncorrectedJetE70_NoBPTX3BX);
   fChain->SetBranchAddress("b5_HLT_L1SingleMu18", &b5_HLT_L1SingleMu18, &b_b5_HLT_L1SingleMu18);
   fChain->SetBranchAddress("b5_HLT_L1SingleMu25", &b5_HLT_L1SingleMu25, &b_b5_HLT_L1SingleMu25);
   fChain->SetBranchAddress("b5_HLT_L2Mu10", &b5_HLT_L2Mu10, &b_b5_HLT_L2Mu10);
   fChain->SetBranchAddress("b5_HLT_L2Mu10_NoVertex_NoBPTX3BX", &b5_HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_b5_HLT_L2Mu10_NoVertex_NoBPTX3BX);
   fChain->SetBranchAddress("b5_HLT_L2Mu10_NoVertex_NoBPTX", &b5_HLT_L2Mu10_NoVertex_NoBPTX, &b_b5_HLT_L2Mu10_NoVertex_NoBPTX);
   fChain->SetBranchAddress("b5_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &b5_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX, &b_b5_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("b5_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &b5_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_b5_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("b5_HLT_L2Mu50", &b5_HLT_L2Mu50, &b_b5_HLT_L2Mu50);
   fChain->SetBranchAddress("b5_HLT_L2Mu23NoVtx_2Cha", &b5_HLT_L2Mu23NoVtx_2Cha, &b_b5_HLT_L2Mu23NoVtx_2Cha);
   fChain->SetBranchAddress("b5_HLT_L2Mu23NoVtx_2Cha_CosmicSeed", &b5_HLT_L2Mu23NoVtx_2Cha_CosmicSeed, &b_b5_HLT_L2Mu23NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4", &b5_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4, &b_b5_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4", &b5_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4, &b_b5_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu50", &b5_HLT_DoubleL2Mu50, &b_b5_HLT_DoubleL2Mu50);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed", &b5_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed, &b_b5_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed", &b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed, &b_b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4", &b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4, &b_b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu23NoVtx_2Cha", &b5_HLT_DoubleL2Mu23NoVtx_2Cha, &b_b5_HLT_DoubleL2Mu23NoVtx_2Cha);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu25NoVtx_2Cha", &b5_HLT_DoubleL2Mu25NoVtx_2Cha, &b_b5_HLT_DoubleL2Mu25NoVtx_2Cha);
   fChain->SetBranchAddress("b5_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4", &b5_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4, &b_b5_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4);
   fChain->SetBranchAddress("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("b5_HLT_Mu25_TkMu0_Onia", &b5_HLT_Mu25_TkMu0_Onia, &b_b5_HLT_Mu25_TkMu0_Onia);
   fChain->SetBranchAddress("b5_HLT_Mu30_TkMu0_Onia", &b5_HLT_Mu30_TkMu0_Onia, &b_b5_HLT_Mu30_TkMu0_Onia);
   fChain->SetBranchAddress("b5_HLT_Mu20_TkMu0_Phi", &b5_HLT_Mu20_TkMu0_Phi, &b_b5_HLT_Mu20_TkMu0_Phi);
   fChain->SetBranchAddress("b5_HLT_Mu25_TkMu0_Phi", &b5_HLT_Mu25_TkMu0_Phi, &b_b5_HLT_Mu25_TkMu0_Phi);
   fChain->SetBranchAddress("b5_HLT_Mu12", &b5_HLT_Mu12, &b_b5_HLT_Mu12);
   fChain->SetBranchAddress("b5_HLT_Mu15", &b5_HLT_Mu15, &b_b5_HLT_Mu15);
   fChain->SetBranchAddress("b5_HLT_Mu20", &b5_HLT_Mu20, &b_b5_HLT_Mu20);
   fChain->SetBranchAddress("b5_HLT_Mu27", &b5_HLT_Mu27, &b_b5_HLT_Mu27);
   fChain->SetBranchAddress("b5_HLT_Mu50", &b5_HLT_Mu50, &b_b5_HLT_Mu50);
   fChain->SetBranchAddress("b5_HLT_Mu55", &b5_HLT_Mu55, &b_b5_HLT_Mu55);
   fChain->SetBranchAddress("b5_HLT_OldMu100", &b5_HLT_OldMu100, &b_b5_HLT_OldMu100);
   fChain->SetBranchAddress("b5_HLT_TkMu100", &b5_HLT_TkMu100, &b_b5_HLT_TkMu100);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve40", &b5_HLT_DiPFJetAve40, &b_b5_HLT_DiPFJetAve40);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve60", &b5_HLT_DiPFJetAve60, &b_b5_HLT_DiPFJetAve60);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve80", &b5_HLT_DiPFJetAve80, &b_b5_HLT_DiPFJetAve80);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve140", &b5_HLT_DiPFJetAve140, &b_b5_HLT_DiPFJetAve140);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve200", &b5_HLT_DiPFJetAve200, &b_b5_HLT_DiPFJetAve200);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve260", &b5_HLT_DiPFJetAve260, &b_b5_HLT_DiPFJetAve260);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve320", &b5_HLT_DiPFJetAve320, &b_b5_HLT_DiPFJetAve320);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve400", &b5_HLT_DiPFJetAve400, &b_b5_HLT_DiPFJetAve400);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve500", &b5_HLT_DiPFJetAve500, &b_b5_HLT_DiPFJetAve500);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve15_HFJEC", &b5_HLT_DiPFJetAve15_HFJEC, &b_b5_HLT_DiPFJetAve15_HFJEC);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve25_HFJEC", &b5_HLT_DiPFJetAve25_HFJEC, &b_b5_HLT_DiPFJetAve25_HFJEC);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve60_HFJEC", &b5_HLT_DiPFJetAve60_HFJEC, &b_b5_HLT_DiPFJetAve60_HFJEC);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve80_HFJEC", &b5_HLT_DiPFJetAve80_HFJEC, &b_b5_HLT_DiPFJetAve80_HFJEC);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve100_HFJEC", &b5_HLT_DiPFJetAve100_HFJEC, &b_b5_HLT_DiPFJetAve100_HFJEC);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve160_HFJEC", &b5_HLT_DiPFJetAve160_HFJEC, &b_b5_HLT_DiPFJetAve160_HFJEC);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve220_HFJEC", &b5_HLT_DiPFJetAve220_HFJEC, &b_b5_HLT_DiPFJetAve220_HFJEC);
   fChain->SetBranchAddress("b5_HLT_DiPFJetAve300_HFJEC", &b5_HLT_DiPFJetAve300_HFJEC, &b_b5_HLT_DiPFJetAve300_HFJEC);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet15", &b5_HLT_AK8PFJet15, &b_b5_HLT_AK8PFJet15);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet25", &b5_HLT_AK8PFJet25, &b_b5_HLT_AK8PFJet25);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet40", &b5_HLT_AK8PFJet40, &b_b5_HLT_AK8PFJet40);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet60", &b5_HLT_AK8PFJet60, &b_b5_HLT_AK8PFJet60);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet80", &b5_HLT_AK8PFJet80, &b_b5_HLT_AK8PFJet80);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet140", &b5_HLT_AK8PFJet140, &b_b5_HLT_AK8PFJet140);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet200", &b5_HLT_AK8PFJet200, &b_b5_HLT_AK8PFJet200);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet260", &b5_HLT_AK8PFJet260, &b_b5_HLT_AK8PFJet260);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet320", &b5_HLT_AK8PFJet320, &b_b5_HLT_AK8PFJet320);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet400", &b5_HLT_AK8PFJet400, &b_b5_HLT_AK8PFJet400);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet450", &b5_HLT_AK8PFJet450, &b_b5_HLT_AK8PFJet450);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet500", &b5_HLT_AK8PFJet500, &b_b5_HLT_AK8PFJet500);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet550", &b5_HLT_AK8PFJet550, &b_b5_HLT_AK8PFJet550);
   fChain->SetBranchAddress("b5_HLT_PFJet15", &b5_HLT_PFJet15, &b_b5_HLT_PFJet15);
   fChain->SetBranchAddress("b5_HLT_PFJet25", &b5_HLT_PFJet25, &b_b5_HLT_PFJet25);
   fChain->SetBranchAddress("b5_HLT_PFJet40", &b5_HLT_PFJet40, &b_b5_HLT_PFJet40);
   fChain->SetBranchAddress("b5_HLT_PFJet60", &b5_HLT_PFJet60, &b_b5_HLT_PFJet60);
   fChain->SetBranchAddress("b5_HLT_PFJet80", &b5_HLT_PFJet80, &b_b5_HLT_PFJet80);
   fChain->SetBranchAddress("b5_HLT_PFJet140", &b5_HLT_PFJet140, &b_b5_HLT_PFJet140);
   fChain->SetBranchAddress("b5_HLT_PFJet200", &b5_HLT_PFJet200, &b_b5_HLT_PFJet200);
   fChain->SetBranchAddress("b5_HLT_PFJet260", &b5_HLT_PFJet260, &b_b5_HLT_PFJet260);
   fChain->SetBranchAddress("b5_HLT_PFJet320", &b5_HLT_PFJet320, &b_b5_HLT_PFJet320);
   fChain->SetBranchAddress("b5_HLT_PFJet400", &b5_HLT_PFJet400, &b_b5_HLT_PFJet400);
   fChain->SetBranchAddress("b5_HLT_PFJet450", &b5_HLT_PFJet450, &b_b5_HLT_PFJet450);
   fChain->SetBranchAddress("b5_HLT_PFJet500", &b5_HLT_PFJet500, &b_b5_HLT_PFJet500);
   fChain->SetBranchAddress("b5_HLT_PFJet550", &b5_HLT_PFJet550, &b_b5_HLT_PFJet550);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd15", &b5_HLT_PFJetFwd15, &b_b5_HLT_PFJetFwd15);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd25", &b5_HLT_PFJetFwd25, &b_b5_HLT_PFJetFwd25);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd40", &b5_HLT_PFJetFwd40, &b_b5_HLT_PFJetFwd40);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd60", &b5_HLT_PFJetFwd60, &b_b5_HLT_PFJetFwd60);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd80", &b5_HLT_PFJetFwd80, &b_b5_HLT_PFJetFwd80);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd140", &b5_HLT_PFJetFwd140, &b_b5_HLT_PFJetFwd140);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd200", &b5_HLT_PFJetFwd200, &b_b5_HLT_PFJetFwd200);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd260", &b5_HLT_PFJetFwd260, &b_b5_HLT_PFJetFwd260);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd320", &b5_HLT_PFJetFwd320, &b_b5_HLT_PFJetFwd320);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd400", &b5_HLT_PFJetFwd400, &b_b5_HLT_PFJetFwd400);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd450", &b5_HLT_PFJetFwd450, &b_b5_HLT_PFJetFwd450);
   fChain->SetBranchAddress("b5_HLT_PFJetFwd500", &b5_HLT_PFJetFwd500, &b_b5_HLT_PFJetFwd500);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd15", &b5_HLT_AK8PFJetFwd15, &b_b5_HLT_AK8PFJetFwd15);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd25", &b5_HLT_AK8PFJetFwd25, &b_b5_HLT_AK8PFJetFwd25);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd40", &b5_HLT_AK8PFJetFwd40, &b_b5_HLT_AK8PFJetFwd40);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd60", &b5_HLT_AK8PFJetFwd60, &b_b5_HLT_AK8PFJetFwd60);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd80", &b5_HLT_AK8PFJetFwd80, &b_b5_HLT_AK8PFJetFwd80);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd140", &b5_HLT_AK8PFJetFwd140, &b_b5_HLT_AK8PFJetFwd140);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd200", &b5_HLT_AK8PFJetFwd200, &b_b5_HLT_AK8PFJetFwd200);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd260", &b5_HLT_AK8PFJetFwd260, &b_b5_HLT_AK8PFJetFwd260);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd320", &b5_HLT_AK8PFJetFwd320, &b_b5_HLT_AK8PFJetFwd320);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd400", &b5_HLT_AK8PFJetFwd400, &b_b5_HLT_AK8PFJetFwd400);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd450", &b5_HLT_AK8PFJetFwd450, &b_b5_HLT_AK8PFJetFwd450);
   fChain->SetBranchAddress("b5_HLT_AK8PFJetFwd500", &b5_HLT_AK8PFJetFwd500, &b_b5_HLT_AK8PFJetFwd500);
   fChain->SetBranchAddress("b5_HLT_PFHT180", &b5_HLT_PFHT180, &b_b5_HLT_PFHT180);
   fChain->SetBranchAddress("b5_HLT_PFHT250", &b5_HLT_PFHT250, &b_b5_HLT_PFHT250);
   fChain->SetBranchAddress("b5_HLT_PFHT370", &b5_HLT_PFHT370, &b_b5_HLT_PFHT370);
   fChain->SetBranchAddress("b5_HLT_PFHT430", &b5_HLT_PFHT430, &b_b5_HLT_PFHT430);
   fChain->SetBranchAddress("b5_HLT_PFHT510", &b5_HLT_PFHT510, &b_b5_HLT_PFHT510);
   fChain->SetBranchAddress("b5_HLT_PFHT590", &b5_HLT_PFHT590, &b_b5_HLT_PFHT590);
   fChain->SetBranchAddress("b5_HLT_PFHT680", &b5_HLT_PFHT680, &b_b5_HLT_PFHT680);
   fChain->SetBranchAddress("b5_HLT_PFHT780", &b5_HLT_PFHT780, &b_b5_HLT_PFHT780);
   fChain->SetBranchAddress("b5_HLT_PFHT890", &b5_HLT_PFHT890, &b_b5_HLT_PFHT890);
   fChain->SetBranchAddress("b5_HLT_PFHT1050", &b5_HLT_PFHT1050, &b_b5_HLT_PFHT1050);
   fChain->SetBranchAddress("b5_HLT_PFHT500_PFMET100_PFMHT100_IDTight", &b5_HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_b5_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFHT500_PFMET110_PFMHT110_IDTight", &b5_HLT_PFHT500_PFMET110_PFMHT110_IDTight, &b_b5_HLT_PFHT500_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFHT700_PFMET85_PFMHT85_IDTight", &b5_HLT_PFHT700_PFMET85_PFMHT85_IDTight, &b_b5_HLT_PFHT700_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFHT700_PFMET95_PFMHT95_IDTight", &b5_HLT_PFHT700_PFMET95_PFMHT95_IDTight, &b_b5_HLT_PFHT700_PFMET95_PFMHT95_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFHT800_PFMET75_PFMHT75_IDTight", &b5_HLT_PFHT800_PFMET75_PFMHT75_IDTight, &b_b5_HLT_PFHT800_PFMET75_PFMHT75_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFHT800_PFMET85_PFMHT85_IDTight", &b5_HLT_PFHT800_PFMET85_PFMHT85_IDTight, &b_b5_HLT_PFHT800_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMET110_PFMHT110_IDTight", &b5_HLT_PFMET110_PFMHT110_IDTight, &b_b5_HLT_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMET120_PFMHT120_IDTight", &b5_HLT_PFMET120_PFMHT120_IDTight, &b_b5_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMET130_PFMHT130_IDTight", &b5_HLT_PFMET130_PFMHT130_IDTight, &b_b5_HLT_PFMET130_PFMHT130_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMET140_PFMHT140_IDTight", &b5_HLT_PFMET140_PFMHT140_IDTight, &b_b5_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1", &b5_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1, &b_b5_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1);
   fChain->SetBranchAddress("b5_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1", &b5_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1, &b_b5_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1);
   fChain->SetBranchAddress("b5_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1", &b5_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1, &b_b5_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1);
   fChain->SetBranchAddress("b5_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1", &b5_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1, &b_b5_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1);
   fChain->SetBranchAddress("b5_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1", &b5_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1, &b_b5_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1);
   fChain->SetBranchAddress("b5_HLT_PFMET120_PFMHT120_IDTight_PFHT60", &b5_HLT_PFMET120_PFMHT120_IDTight_PFHT60, &b_b5_HLT_PFMET120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("b5_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60", &b5_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60, &b_b5_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("b5_HLT_PFMETTypeOne110_PFMHT110_IDTight", &b5_HLT_PFMETTypeOne110_PFMHT110_IDTight, &b_b5_HLT_PFMETTypeOne110_PFMHT110_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMETTypeOne120_PFMHT120_IDTight", &b5_HLT_PFMETTypeOne120_PFMHT120_IDTight, &b_b5_HLT_PFMETTypeOne120_PFMHT120_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMETTypeOne130_PFMHT130_IDTight", &b5_HLT_PFMETTypeOne130_PFMHT130_IDTight, &b_b5_HLT_PFMETTypeOne130_PFMHT130_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMETTypeOne140_PFMHT140_IDTight", &b5_HLT_PFMETTypeOne140_PFMHT140_IDTight, &b_b5_HLT_PFMETTypeOne140_PFMHT140_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &b5_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_b5_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &b5_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_b5_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("b5_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &b5_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_b5_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("b5_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &b5_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_b5_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("b5_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &b5_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_b5_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("b5_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight", &b5_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_b5_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("b5_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight", &b5_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_b5_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("b5_HLT_L1ETMHadSeeds", &b5_HLT_L1ETMHadSeeds, &b_b5_HLT_L1ETMHadSeeds);
   fChain->SetBranchAddress("b5_HLT_CaloMHT90", &b5_HLT_CaloMHT90, &b_b5_HLT_CaloMHT90);
   fChain->SetBranchAddress("b5_HLT_CaloMET80_NotCleaned", &b5_HLT_CaloMET80_NotCleaned, &b_b5_HLT_CaloMET80_NotCleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET90_NotCleaned", &b5_HLT_CaloMET90_NotCleaned, &b_b5_HLT_CaloMET90_NotCleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET100_NotCleaned", &b5_HLT_CaloMET100_NotCleaned, &b_b5_HLT_CaloMET100_NotCleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET110_NotCleaned", &b5_HLT_CaloMET110_NotCleaned, &b_b5_HLT_CaloMET110_NotCleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET250_NotCleaned", &b5_HLT_CaloMET250_NotCleaned, &b_b5_HLT_CaloMET250_NotCleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET70_HBHECleaned", &b5_HLT_CaloMET70_HBHECleaned, &b_b5_HLT_CaloMET70_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET80_HBHECleaned", &b5_HLT_CaloMET80_HBHECleaned, &b_b5_HLT_CaloMET80_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET90_HBHECleaned", &b5_HLT_CaloMET90_HBHECleaned, &b_b5_HLT_CaloMET90_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET100_HBHECleaned", &b5_HLT_CaloMET100_HBHECleaned, &b_b5_HLT_CaloMET100_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET250_HBHECleaned", &b5_HLT_CaloMET250_HBHECleaned, &b_b5_HLT_CaloMET250_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET300_HBHECleaned", &b5_HLT_CaloMET300_HBHECleaned, &b_b5_HLT_CaloMET300_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_CaloMET350_HBHECleaned", &b5_HLT_CaloMET350_HBHECleaned, &b_b5_HLT_CaloMET350_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_PFMET200_NotCleaned", &b5_HLT_PFMET200_NotCleaned, &b_b5_HLT_PFMET200_NotCleaned);
   fChain->SetBranchAddress("b5_HLT_PFMET200_HBHECleaned", &b5_HLT_PFMET200_HBHECleaned, &b_b5_HLT_PFMET200_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_PFMET250_HBHECleaned", &b5_HLT_PFMET250_HBHECleaned, &b_b5_HLT_PFMET250_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_PFMET300_HBHECleaned", &b5_HLT_PFMET300_HBHECleaned, &b_b5_HLT_PFMET300_HBHECleaned);
   fChain->SetBranchAddress("b5_HLT_PFMET200_HBHE_BeamHaloCleaned", &b5_HLT_PFMET200_HBHE_BeamHaloCleaned, &b_b5_HLT_PFMET200_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("b5_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned", &b5_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned, &b_b5_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("b5_HLT_MET105_IsoTrk50", &b5_HLT_MET105_IsoTrk50, &b_b5_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("b5_HLT_MET120_IsoTrk50", &b5_HLT_MET120_IsoTrk50, &b_b5_HLT_MET120_IsoTrk50);
   fChain->SetBranchAddress("b5_HLT_SingleJet30_Mu12_SinglePFJet40", &b5_HLT_SingleJet30_Mu12_SinglePFJet40, &b_b5_HLT_SingleJet30_Mu12_SinglePFJet40);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71", &b5_HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71, &b_b5_HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71", &b5_HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71, &b_b5_HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71", &b5_HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71, &b_b5_HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71", &b5_HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71, &b_b5_HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &b5_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71, &b_b5_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &b5_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71, &b_b5_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &b5_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71, &b_b5_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_DoublePFJets40_CaloBTagDeepCSV_p71", &b5_HLT_DoublePFJets40_CaloBTagDeepCSV_p71, &b_b5_HLT_DoublePFJets40_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_DoublePFJets100_CaloBTagDeepCSV_p71", &b5_HLT_DoublePFJets100_CaloBTagDeepCSV_p71, &b_b5_HLT_DoublePFJets100_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_DoublePFJets200_CaloBTagDeepCSV_p71", &b5_HLT_DoublePFJets200_CaloBTagDeepCSV_p71, &b_b5_HLT_DoublePFJets200_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_DoublePFJets350_CaloBTagDeepCSV_p71", &b5_HLT_DoublePFJets350_CaloBTagDeepCSV_p71, &b_b5_HLT_DoublePFJets350_CaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &b5_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71, &b_b5_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71", &b5_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71, &b_b5_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71);
   fChain->SetBranchAddress("b5_HLT_Photon300_NoHE", &b5_HLT_Photon300_NoHE, &b_b5_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("b5_HLT_Mu8_TrkIsoVVL", &b5_HLT_Mu8_TrkIsoVVL, &b_b5_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, &b_b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, &b_b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("b5_HLT_Mu17_TrkIsoVVL", &b5_HLT_Mu17_TrkIsoVVL, &b_b5_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("b5_HLT_Mu19_TrkIsoVVL", &b5_HLT_Mu19_TrkIsoVVL, &b_b5_HLT_Mu19_TrkIsoVVL);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK4DiJet20_Mu5", &b5_HLT_BTagMu_AK4DiJet20_Mu5, &b_b5_HLT_BTagMu_AK4DiJet20_Mu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK4DiJet40_Mu5", &b5_HLT_BTagMu_AK4DiJet40_Mu5, &b_b5_HLT_BTagMu_AK4DiJet40_Mu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK4DiJet70_Mu5", &b5_HLT_BTagMu_AK4DiJet70_Mu5, &b_b5_HLT_BTagMu_AK4DiJet70_Mu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK4DiJet110_Mu5", &b5_HLT_BTagMu_AK4DiJet110_Mu5, &b_b5_HLT_BTagMu_AK4DiJet110_Mu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK4DiJet170_Mu5", &b5_HLT_BTagMu_AK4DiJet170_Mu5, &b_b5_HLT_BTagMu_AK4DiJet170_Mu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK4Jet300_Mu5", &b5_HLT_BTagMu_AK4Jet300_Mu5, &b_b5_HLT_BTagMu_AK4Jet300_Mu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK8DiJet170_Mu5", &b5_HLT_BTagMu_AK8DiJet170_Mu5, &b_b5_HLT_BTagMu_AK8DiJet170_Mu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK8Jet170_DoubleMu5", &b5_HLT_BTagMu_AK8Jet170_DoubleMu5, &b_b5_HLT_BTagMu_AK8Jet170_DoubleMu5);
   fChain->SetBranchAddress("b5_HLT_BTagMu_AK8Jet300_Mu5", &b5_HLT_BTagMu_AK8Jet300_Mu5, &b_b5_HLT_BTagMu_AK8Jet300_Mu5);
   fChain->SetBranchAddress("b5_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL", &b5_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL, &b_b5_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu12_DoublePhoton20", &b5_HLT_Mu12_DoublePhoton20, &b_b5_HLT_Mu12_DoublePhoton20);
   fChain->SetBranchAddress("b5_HLT_TriplePhoton_20_20_20_CaloIdLV2", &b5_HLT_TriplePhoton_20_20_20_CaloIdLV2, &b_b5_HLT_TriplePhoton_20_20_20_CaloIdLV2);
   fChain->SetBranchAddress("b5_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL", &b5_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL, &b_b5_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("b5_HLT_TriplePhoton_30_30_10_CaloIdLV2", &b5_HLT_TriplePhoton_30_30_10_CaloIdLV2, &b_b5_HLT_TriplePhoton_30_30_10_CaloIdLV2);
   fChain->SetBranchAddress("b5_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL", &b5_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL, &b_b5_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("b5_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL", &b5_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL, &b_b5_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL);
   fChain->SetBranchAddress("b5_HLT_Photon20", &b5_HLT_Photon20, &b_b5_HLT_Photon20);
   fChain->SetBranchAddress("b5_HLT_Photon33", &b5_HLT_Photon33, &b_b5_HLT_Photon33);
   fChain->SetBranchAddress("b5_HLT_Photon50", &b5_HLT_Photon50, &b_b5_HLT_Photon50);
   fChain->SetBranchAddress("b5_HLT_Photon75", &b5_HLT_Photon75, &b_b5_HLT_Photon75);
   fChain->SetBranchAddress("b5_HLT_Photon90", &b5_HLT_Photon90, &b_b5_HLT_Photon90);
   fChain->SetBranchAddress("b5_HLT_Photon120", &b5_HLT_Photon120, &b_b5_HLT_Photon120);
   fChain->SetBranchAddress("b5_HLT_Photon150", &b5_HLT_Photon150, &b_b5_HLT_Photon150);
   fChain->SetBranchAddress("b5_HLT_Photon175", &b5_HLT_Photon175, &b_b5_HLT_Photon175);
   fChain->SetBranchAddress("b5_HLT_Photon200", &b5_HLT_Photon200, &b_b5_HLT_Photon200);
   fChain->SetBranchAddress("b5_HLT_Photon100EB_TightID_TightIso", &b5_HLT_Photon100EB_TightID_TightIso, &b_b5_HLT_Photon100EB_TightID_TightIso);
   fChain->SetBranchAddress("b5_HLT_Photon110EB_TightID_TightIso", &b5_HLT_Photon110EB_TightID_TightIso, &b_b5_HLT_Photon110EB_TightID_TightIso);
   fChain->SetBranchAddress("b5_HLT_Photon120EB_TightID_TightIso", &b5_HLT_Photon120EB_TightID_TightIso, &b_b5_HLT_Photon120EB_TightID_TightIso);
   fChain->SetBranchAddress("b5_HLT_Photon100EBHE10", &b5_HLT_Photon100EBHE10, &b_b5_HLT_Photon100EBHE10);
   fChain->SetBranchAddress("b5_HLT_Photon100EEHE10", &b5_HLT_Photon100EEHE10, &b_b5_HLT_Photon100EEHE10);
   fChain->SetBranchAddress("b5_HLT_Photon100EE_TightID_TightIso", &b5_HLT_Photon100EE_TightID_TightIso, &b_b5_HLT_Photon100EE_TightID_TightIso);
   fChain->SetBranchAddress("b5_HLT_Photon50_R9Id90_HE10_IsoM", &b5_HLT_Photon50_R9Id90_HE10_IsoM, &b_b5_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("b5_HLT_Photon75_R9Id90_HE10_IsoM", &b5_HLT_Photon75_R9Id90_HE10_IsoM, &b_b5_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3", &b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3, &b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3);
   fChain->SetBranchAddress("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3", &b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3, &b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3);
   fChain->SetBranchAddress("b5_HLT_Photon90_R9Id90_HE10_IsoM", &b5_HLT_Photon90_R9Id90_HE10_IsoM, &b_b5_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("b5_HLT_Photon120_R9Id90_HE10_IsoM", &b5_HLT_Photon120_R9Id90_HE10_IsoM, &b_b5_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("b5_HLT_Photon165_R9Id90_HE10_IsoM", &b5_HLT_Photon165_R9Id90_HE10_IsoM, &b_b5_HLT_Photon165_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   fChain->SetBranchAddress("b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95, &b_b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
   fChain->SetBranchAddress("b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, &b_b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55);
   fChain->SetBranchAddress("b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55, &b_b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Jpsi_L1_NoOS", &b5_HLT_Dimuon0_Jpsi_L1_NoOS, &b_b5_HLT_Dimuon0_Jpsi_L1_NoOS);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &b5_HLT_Dimuon0_Jpsi_NoVertexing_NoOS, &b_b5_HLT_Dimuon0_Jpsi_NoVertexing_NoOS);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Jpsi", &b5_HLT_Dimuon0_Jpsi, &b_b5_HLT_Dimuon0_Jpsi);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Jpsi_NoVertexing", &b5_HLT_Dimuon0_Jpsi_NoVertexing, &b_b5_HLT_Dimuon0_Jpsi_NoVertexing);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &b5_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R, &b_b5_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R, &b_b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Jpsi3p5_Muon2", &b5_HLT_Dimuon0_Jpsi3p5_Muon2, &b_b5_HLT_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_L1_4p5", &b5_HLT_Dimuon0_Upsilon_L1_4p5, &b_b5_HLT_Dimuon0_Upsilon_L1_4p5);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_L1_5", &b5_HLT_Dimuon0_Upsilon_L1_5, &b_b5_HLT_Dimuon0_Upsilon_L1_5);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_L1_4p5NoOS", &b5_HLT_Dimuon0_Upsilon_L1_4p5NoOS, &b_b5_HLT_Dimuon0_Upsilon_L1_4p5NoOS);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0", &b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0, &b_b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0M, &b_b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0M);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_NoVertexing", &b5_HLT_Dimuon0_Upsilon_NoVertexing, &b_b5_HLT_Dimuon0_Upsilon_NoVertexing);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_L1_5M", &b5_HLT_Dimuon0_Upsilon_L1_5M, &b_b5_HLT_Dimuon0_Upsilon_L1_5M);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_LowMass_L1_0er1p5R", &b5_HLT_Dimuon0_LowMass_L1_0er1p5R, &b_b5_HLT_Dimuon0_LowMass_L1_0er1p5R);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_LowMass_L1_0er1p5", &b5_HLT_Dimuon0_LowMass_L1_0er1p5, &b_b5_HLT_Dimuon0_LowMass_L1_0er1p5);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_LowMass", &b5_HLT_Dimuon0_LowMass, &b_b5_HLT_Dimuon0_LowMass);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_LowMass_L1_4", &b5_HLT_Dimuon0_LowMass_L1_4, &b_b5_HLT_Dimuon0_LowMass_L1_4);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_LowMass_L1_4R", &b5_HLT_Dimuon0_LowMass_L1_4R, &b_b5_HLT_Dimuon0_LowMass_L1_4R);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_LowMass_L1_TM530", &b5_HLT_Dimuon0_LowMass_L1_TM530, &b_b5_HLT_Dimuon0_LowMass_L1_TM530);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_Muon_L1_TM0", &b5_HLT_Dimuon0_Upsilon_Muon_L1_TM0, &b_b5_HLT_Dimuon0_Upsilon_Muon_L1_TM0);
   fChain->SetBranchAddress("b5_HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &b5_HLT_Dimuon0_Upsilon_Muon_NoL1Mass, &b_b5_HLT_Dimuon0_Upsilon_Muon_NoL1Mass);
   fChain->SetBranchAddress("b5_HLT_TripleMu_5_3_3_Mass3p8_DZ", &b5_HLT_TripleMu_5_3_3_Mass3p8_DZ, &b_b5_HLT_TripleMu_5_3_3_Mass3p8_DZ);
   fChain->SetBranchAddress("b5_HLT_TripleMu_10_5_5_DZ", &b5_HLT_TripleMu_10_5_5_DZ, &b_b5_HLT_TripleMu_10_5_5_DZ);
   fChain->SetBranchAddress("b5_HLT_TripleMu_12_10_5", &b5_HLT_TripleMu_12_10_5, &b_b5_HLT_TripleMu_12_10_5);
   fChain->SetBranchAddress("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15, &b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15);
   fChain->SetBranchAddress("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1, &b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1);
   fChain->SetBranchAddress("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15, &b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15);
   fChain->SetBranchAddress("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1, &b_b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &b5_HLT_DoubleMu3_DZ_PFMET50_PFMHT60, &b_b5_HLT_DoubleMu3_DZ_PFMET50_PFMHT60);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &b5_HLT_DoubleMu3_DZ_PFMET70_PFMHT70, &b_b5_HLT_DoubleMu3_DZ_PFMET70_PFMHT70);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &b5_HLT_DoubleMu3_DZ_PFMET90_PFMHT90, &b_b5_HLT_DoubleMu3_DZ_PFMET90_PFMHT90);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &b5_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass, &b_b5_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_Jpsi_Displaced", &b5_HLT_DoubleMu4_Jpsi_Displaced, &b_b5_HLT_DoubleMu4_Jpsi_Displaced);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_Jpsi_NoVertexing", &b5_HLT_DoubleMu4_Jpsi_NoVertexing, &b_b5_HLT_DoubleMu4_Jpsi_NoVertexing);
   fChain->SetBranchAddress("b5_HLT_DoubleMu4_JpsiTrkTrk_Displaced", &b5_HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_b5_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("b5_HLT_DoubleMu43NoFiltersNoVtx", &b5_HLT_DoubleMu43NoFiltersNoVtx, &b_b5_HLT_DoubleMu43NoFiltersNoVtx);
   fChain->SetBranchAddress("b5_HLT_DoubleMu48NoFiltersNoVtx", &b5_HLT_DoubleMu48NoFiltersNoVtx, &b_b5_HLT_DoubleMu48NoFiltersNoVtx);
   fChain->SetBranchAddress("b5_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &b5_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL, &b_b5_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL);
   fChain->SetBranchAddress("b5_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &b5_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL, &b_b5_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL);
   fChain->SetBranchAddress("b5_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", &b5_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL, &b_b5_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL);
   fChain->SetBranchAddress("b5_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL", &b5_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL, &b_b5_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL);
   fChain->SetBranchAddress("b5_HLT_DoubleMu33NoFiltersNoVtxDisplaced", &b5_HLT_DoubleMu33NoFiltersNoVtxDisplaced, &b_b5_HLT_DoubleMu33NoFiltersNoVtxDisplaced);
   fChain->SetBranchAddress("b5_HLT_DoubleMu40NoFiltersNoVtxDisplaced", &b5_HLT_DoubleMu40NoFiltersNoVtxDisplaced, &b_b5_HLT_DoubleMu40NoFiltersNoVtxDisplaced);
   fChain->SetBranchAddress("b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4", &b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4, &b_b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4);
   fChain->SetBranchAddress("b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG", &b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG, &b_b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG);
   fChain->SetBranchAddress("b5_HLT_HT425", &b5_HLT_HT425, &b_b5_HLT_HT425);
   fChain->SetBranchAddress("b5_HLT_HT430_DisplacedDijet40_DisplacedTrack", &b5_HLT_HT430_DisplacedDijet40_DisplacedTrack, &b_b5_HLT_HT430_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("b5_HLT_HT500_DisplacedDijet40_DisplacedTrack", &b5_HLT_HT500_DisplacedDijet40_DisplacedTrack, &b_b5_HLT_HT500_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("b5_HLT_HT430_DisplacedDijet60_DisplacedTrack", &b5_HLT_HT430_DisplacedDijet60_DisplacedTrack, &b_b5_HLT_HT430_DisplacedDijet60_DisplacedTrack);
   fChain->SetBranchAddress("b5_HLT_HT400_DisplacedDijet40_DisplacedTrack", &b5_HLT_HT400_DisplacedDijet40_DisplacedTrack, &b_b5_HLT_HT400_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("b5_HLT_HT650_DisplacedDijet60_Inclusive", &b5_HLT_HT650_DisplacedDijet60_Inclusive, &b_b5_HLT_HT650_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("b5_HLT_HT550_DisplacedDijet60_Inclusive", &b5_HLT_HT550_DisplacedDijet60_Inclusive, &b_b5_HLT_HT550_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("b5_HLT_DiJet110_35_Mjj650_PFMET110", &b5_HLT_DiJet110_35_Mjj650_PFMET110, &b_b5_HLT_DiJet110_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("b5_HLT_DiJet110_35_Mjj650_PFMET120", &b5_HLT_DiJet110_35_Mjj650_PFMET120, &b_b5_HLT_DiJet110_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("b5_HLT_DiJet110_35_Mjj650_PFMET130", &b5_HLT_DiJet110_35_Mjj650_PFMET130, &b_b5_HLT_DiJet110_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("b5_HLT_TripleJet110_35_35_Mjj650_PFMET110", &b5_HLT_TripleJet110_35_35_Mjj650_PFMET110, &b_b5_HLT_TripleJet110_35_35_Mjj650_PFMET110);
   fChain->SetBranchAddress("b5_HLT_TripleJet110_35_35_Mjj650_PFMET120", &b5_HLT_TripleJet110_35_35_Mjj650_PFMET120, &b_b5_HLT_TripleJet110_35_35_Mjj650_PFMET120);
   fChain->SetBranchAddress("b5_HLT_TripleJet110_35_35_Mjj650_PFMET130", &b5_HLT_TripleJet110_35_35_Mjj650_PFMET130, &b_b5_HLT_TripleJet110_35_35_Mjj650_PFMET130);
   fChain->SetBranchAddress("b5_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1", &b5_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1, &b_b5_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1);
   fChain->SetBranchAddress("b5_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1", &b5_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1, &b_b5_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1);
   fChain->SetBranchAddress("b5_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1", &b5_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1, &b_b5_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1);
   fChain->SetBranchAddress("b5_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &b5_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned, &b_b5_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
   fChain->SetBranchAddress("b5_HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &b5_HLT_Ele28_eta2p1_WPTight_Gsf_HT150, &b_b5_HLT_Ele28_eta2p1_WPTight_Gsf_HT150);
   fChain->SetBranchAddress("b5_HLT_Ele28_HighEta_SC20_Mass55", &b5_HLT_Ele28_HighEta_SC20_Mass55, &b_b5_HLT_Ele28_HighEta_SC20_Mass55);
   fChain->SetBranchAddress("b5_HLT_DoubleMu20_7_Mass0to30_Photon23", &b5_HLT_DoubleMu20_7_Mass0to30_Photon23, &b_b5_HLT_DoubleMu20_7_Mass0to30_Photon23);
   fChain->SetBranchAddress("b5_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &b5_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5, &b_b5_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5);
   fChain->SetBranchAddress("b5_HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &b5_HLT_Ele15_IsoVVVL_PFHT450_PFMET50, &b_b5_HLT_Ele15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("b5_HLT_Ele15_IsoVVVL_PFHT450", &b5_HLT_Ele15_IsoVVVL_PFHT450, &b_b5_HLT_Ele15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("b5_HLT_Ele50_IsoVVVL_PFHT450", &b5_HLT_Ele50_IsoVVVL_PFHT450, &b_b5_HLT_Ele50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("b5_HLT_Ele15_IsoVVVL_PFHT600", &b5_HLT_Ele15_IsoVVVL_PFHT600, &b_b5_HLT_Ele15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("b5_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &b5_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_b5_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   fChain->SetBranchAddress("b5_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &b5_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_b5_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   fChain->SetBranchAddress("b5_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &b5_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60, &b_b5_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
   fChain->SetBranchAddress("b5_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5", &b5_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5, &b_b5_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5);
   fChain->SetBranchAddress("b5_HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &b5_HLT_Mu15_IsoVVVL_PFHT450_PFMET50, &b_b5_HLT_Mu15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("b5_HLT_Mu15_IsoVVVL_PFHT450", &b5_HLT_Mu15_IsoVVVL_PFHT450, &b_b5_HLT_Mu15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("b5_HLT_Mu50_IsoVVVL_PFHT450", &b5_HLT_Mu50_IsoVVVL_PFHT450, &b_b5_HLT_Mu50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("b5_HLT_Mu15_IsoVVVL_PFHT600", &b5_HLT_Mu15_IsoVVVL_PFHT600, &b_b5_HLT_Mu15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight);
   fChain->SetBranchAddress("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight", &b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight, &b_b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight);
   fChain->SetBranchAddress("b5_HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &b5_HLT_Dimuon10_PsiPrime_Barrel_Seagulls, &b_b5_HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
   fChain->SetBranchAddress("b5_HLT_Dimuon20_Jpsi_Barrel_Seagulls", &b5_HLT_Dimuon20_Jpsi_Barrel_Seagulls, &b_b5_HLT_Dimuon20_Jpsi_Barrel_Seagulls);
   fChain->SetBranchAddress("b5_HLT_Dimuon12_Upsilon_y1p4", &b5_HLT_Dimuon12_Upsilon_y1p4, &b_b5_HLT_Dimuon12_Upsilon_y1p4);
   fChain->SetBranchAddress("b5_HLT_Dimuon14_Phi_Barrel_Seagulls", &b5_HLT_Dimuon14_Phi_Barrel_Seagulls, &b_b5_HLT_Dimuon14_Phi_Barrel_Seagulls);
   fChain->SetBranchAddress("b5_HLT_Dimuon18_PsiPrime", &b5_HLT_Dimuon18_PsiPrime, &b_b5_HLT_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("b5_HLT_Dimuon25_Jpsi", &b5_HLT_Dimuon25_Jpsi, &b_b5_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("b5_HLT_Dimuon18_PsiPrime_noCorrL1", &b5_HLT_Dimuon18_PsiPrime_noCorrL1, &b_b5_HLT_Dimuon18_PsiPrime_noCorrL1);
   fChain->SetBranchAddress("b5_HLT_Dimuon24_Upsilon_noCorrL1", &b5_HLT_Dimuon24_Upsilon_noCorrL1, &b_b5_HLT_Dimuon24_Upsilon_noCorrL1);
   fChain->SetBranchAddress("b5_HLT_Dimuon24_Phi_noCorrL1", &b5_HLT_Dimuon24_Phi_noCorrL1, &b_b5_HLT_Dimuon24_Phi_noCorrL1);
   fChain->SetBranchAddress("b5_HLT_Dimuon25_Jpsi_noCorrL1", &b5_HLT_Dimuon25_Jpsi_noCorrL1, &b_b5_HLT_Dimuon25_Jpsi_noCorrL1);
   fChain->SetBranchAddress("b5_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8", &b5_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8, &b_b5_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8);
   fChain->SetBranchAddress("b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("b5_HLT_DoubleIsoMu20_eta2p1", &b5_HLT_DoubleIsoMu20_eta2p1, &b_b5_HLT_DoubleIsoMu20_eta2p1);
   fChain->SetBranchAddress("b5_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &b5_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx, &b_b5_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx);
   fChain->SetBranchAddress("b5_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx", &b5_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx, &b_b5_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx);
   fChain->SetBranchAddress("b5_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &b5_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx, &b_b5_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
   fChain->SetBranchAddress("b5_HLT_Mu8", &b5_HLT_Mu8, &b_b5_HLT_Mu8);
   fChain->SetBranchAddress("b5_HLT_Mu17", &b5_HLT_Mu17, &b_b5_HLT_Mu17);
   fChain->SetBranchAddress("b5_HLT_Mu19", &b5_HLT_Mu19, &b_b5_HLT_Mu19);
   fChain->SetBranchAddress("b5_HLT_Mu17_Photon30_IsoCaloId", &b5_HLT_Mu17_Photon30_IsoCaloId, &b_b5_HLT_Mu17_Photon30_IsoCaloId);
   fChain->SetBranchAddress("b5_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &b5_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_b5_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("b5_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &b5_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_b5_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("b5_HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30", &b5_HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_b5_HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("b5_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &b5_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_b5_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("b5_HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &b5_HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_b5_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("b5_HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &b5_HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_b5_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("b5_HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &b5_HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_b5_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("b5_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &b5_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_b5_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   fChain->SetBranchAddress("b5_HLT_Ele115_CaloIdVT_GsfTrkIdT", &b5_HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_b5_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("b5_HLT_Ele135_CaloIdVT_GsfTrkIdT", &b5_HLT_Ele135_CaloIdVT_GsfTrkIdT, &b_b5_HLT_Ele135_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("b5_HLT_Ele145_CaloIdVT_GsfTrkIdT", &b5_HLT_Ele145_CaloIdVT_GsfTrkIdT, &b_b5_HLT_Ele145_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("b5_HLT_Ele200_CaloIdVT_GsfTrkIdT", &b5_HLT_Ele200_CaloIdVT_GsfTrkIdT, &b_b5_HLT_Ele200_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("b5_HLT_Ele250_CaloIdVT_GsfTrkIdT", &b5_HLT_Ele250_CaloIdVT_GsfTrkIdT, &b_b5_HLT_Ele250_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("b5_HLT_Ele300_CaloIdVT_GsfTrkIdT", &b5_HLT_Ele300_CaloIdVT_GsfTrkIdT, &b_b5_HLT_Ele300_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5", &b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5, &b_b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5);
   fChain->SetBranchAddress("b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40", &b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40, &b_b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40);
   fChain->SetBranchAddress("b5_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2", &b5_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2, &b_b5_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2);
   fChain->SetBranchAddress("b5_HLT_PFHT380_SixPFJet32", &b5_HLT_PFHT380_SixPFJet32, &b_b5_HLT_PFHT380_SixPFJet32);
   fChain->SetBranchAddress("b5_HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5", &b5_HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5, &b_b5_HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5);
   fChain->SetBranchAddress("b5_HLT_PFHT430_SixPFJet40", &b5_HLT_PFHT430_SixPFJet40, &b_b5_HLT_PFHT430_SixPFJet40);
   fChain->SetBranchAddress("b5_HLT_PFHT350", &b5_HLT_PFHT350, &b_b5_HLT_PFHT350);
   fChain->SetBranchAddress("b5_HLT_PFHT350MinPFJet15", &b5_HLT_PFHT350MinPFJet15, &b_b5_HLT_PFHT350MinPFJet15);
   fChain->SetBranchAddress("b5_HLT_Photon60_R9Id90_CaloIdL_IsoL", &b5_HLT_Photon60_R9Id90_CaloIdL_IsoL, &b_b5_HLT_Photon60_R9Id90_CaloIdL_IsoL);
   fChain->SetBranchAddress("b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL", &b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL, &b_b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL);
   fChain->SetBranchAddress("b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15", &b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15, &b_b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15);
   fChain->SetBranchAddress("b5_HLT_ECALHT800", &b5_HLT_ECALHT800, &b_b5_HLT_ECALHT800);
   fChain->SetBranchAddress("b5_HLT_DiSC30_18_EIso_AND_HE_Mass70", &b5_HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_b5_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   fChain->SetBranchAddress("b5_HLT_Physics", &b5_HLT_Physics, &b_b5_HLT_Physics);
   fChain->SetBranchAddress("b5_HLT_Physics_part0", &b5_HLT_Physics_part0, &b_b5_HLT_Physics_part0);
   fChain->SetBranchAddress("b5_HLT_Physics_part1", &b5_HLT_Physics_part1, &b_b5_HLT_Physics_part1);
   fChain->SetBranchAddress("b5_HLT_Physics_part2", &b5_HLT_Physics_part2, &b_b5_HLT_Physics_part2);
   fChain->SetBranchAddress("b5_HLT_Physics_part3", &b5_HLT_Physics_part3, &b_b5_HLT_Physics_part3);
   fChain->SetBranchAddress("b5_HLT_Physics_part4", &b5_HLT_Physics_part4, &b_b5_HLT_Physics_part4);
   fChain->SetBranchAddress("b5_HLT_Physics_part5", &b5_HLT_Physics_part5, &b_b5_HLT_Physics_part5);
   fChain->SetBranchAddress("b5_HLT_Physics_part6", &b5_HLT_Physics_part6, &b_b5_HLT_Physics_part6);
   fChain->SetBranchAddress("b5_HLT_Physics_part7", &b5_HLT_Physics_part7, &b_b5_HLT_Physics_part7);
   fChain->SetBranchAddress("b5_HLT_Random", &b5_HLT_Random, &b_b5_HLT_Random);
   fChain->SetBranchAddress("b5_HLT_ZeroBias", &b5_HLT_ZeroBias, &b_b5_HLT_ZeroBias);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part0", &b5_HLT_ZeroBias_part0, &b_b5_HLT_ZeroBias_part0);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part1", &b5_HLT_ZeroBias_part1, &b_b5_HLT_ZeroBias_part1);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part2", &b5_HLT_ZeroBias_part2, &b_b5_HLT_ZeroBias_part2);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part3", &b5_HLT_ZeroBias_part3, &b_b5_HLT_ZeroBias_part3);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part4", &b5_HLT_ZeroBias_part4, &b_b5_HLT_ZeroBias_part4);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part5", &b5_HLT_ZeroBias_part5, &b_b5_HLT_ZeroBias_part5);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part6", &b5_HLT_ZeroBias_part6, &b_b5_HLT_ZeroBias_part6);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_part7", &b5_HLT_ZeroBias_part7, &b_b5_HLT_ZeroBias_part7);
   fChain->SetBranchAddress("b5_HLT_AK4CaloJet30", &b5_HLT_AK4CaloJet30, &b_b5_HLT_AK4CaloJet30);
   fChain->SetBranchAddress("b5_HLT_AK4CaloJet40", &b5_HLT_AK4CaloJet40, &b_b5_HLT_AK4CaloJet40);
   fChain->SetBranchAddress("b5_HLT_AK4CaloJet50", &b5_HLT_AK4CaloJet50, &b_b5_HLT_AK4CaloJet50);
   fChain->SetBranchAddress("b5_HLT_AK4CaloJet80", &b5_HLT_AK4CaloJet80, &b_b5_HLT_AK4CaloJet80);
   fChain->SetBranchAddress("b5_HLT_AK4CaloJet100", &b5_HLT_AK4CaloJet100, &b_b5_HLT_AK4CaloJet100);
   fChain->SetBranchAddress("b5_HLT_AK4CaloJet120", &b5_HLT_AK4CaloJet120, &b_b5_HLT_AK4CaloJet120);
   fChain->SetBranchAddress("b5_HLT_AK4PFJet30", &b5_HLT_AK4PFJet30, &b_b5_HLT_AK4PFJet30);
   fChain->SetBranchAddress("b5_HLT_AK4PFJet50", &b5_HLT_AK4PFJet50, &b_b5_HLT_AK4PFJet50);
   fChain->SetBranchAddress("b5_HLT_AK4PFJet80", &b5_HLT_AK4PFJet80, &b_b5_HLT_AK4PFJet80);
   fChain->SetBranchAddress("b5_HLT_AK4PFJet100", &b5_HLT_AK4PFJet100, &b_b5_HLT_AK4PFJet100);
   fChain->SetBranchAddress("b5_HLT_AK4PFJet120", &b5_HLT_AK4PFJet120, &b_b5_HLT_AK4PFJet120);
   fChain->SetBranchAddress("b5_HLT_SinglePhoton10_Eta3p1ForPPRef", &b5_HLT_SinglePhoton10_Eta3p1ForPPRef, &b_b5_HLT_SinglePhoton10_Eta3p1ForPPRef);
   fChain->SetBranchAddress("b5_HLT_SinglePhoton20_Eta3p1ForPPRef", &b5_HLT_SinglePhoton20_Eta3p1ForPPRef, &b_b5_HLT_SinglePhoton20_Eta3p1ForPPRef);
   fChain->SetBranchAddress("b5_HLT_SinglePhoton30_Eta3p1ForPPRef", &b5_HLT_SinglePhoton30_Eta3p1ForPPRef, &b_b5_HLT_SinglePhoton30_Eta3p1ForPPRef);
   fChain->SetBranchAddress("b5_HLT_Photon20_HoverELoose", &b5_HLT_Photon20_HoverELoose, &b_b5_HLT_Photon20_HoverELoose);
   fChain->SetBranchAddress("b5_HLT_Photon30_HoverELoose", &b5_HLT_Photon30_HoverELoose, &b_b5_HLT_Photon30_HoverELoose);
   fChain->SetBranchAddress("b5_HLT_EcalCalibration", &b5_HLT_EcalCalibration, &b_b5_HLT_EcalCalibration);
   fChain->SetBranchAddress("b5_HLT_HcalCalibration", &b5_HLT_HcalCalibration, &b_b5_HLT_HcalCalibration);
   fChain->SetBranchAddress("b5_HLT_L1UnpairedBunchBptxMinus", &b5_HLT_L1UnpairedBunchBptxMinus, &b_b5_HLT_L1UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("b5_HLT_L1UnpairedBunchBptxPlus", &b5_HLT_L1UnpairedBunchBptxPlus, &b_b5_HLT_L1UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("b5_HLT_L1NotBptxOR", &b5_HLT_L1NotBptxOR, &b_b5_HLT_L1NotBptxOR);
   fChain->SetBranchAddress("b5_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &b5_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_b5_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("b5_HLT_HcalNZS", &b5_HLT_HcalNZS, &b_b5_HLT_HcalNZS);
   fChain->SetBranchAddress("b5_HLT_HcalPhiSym", &b5_HLT_HcalPhiSym, &b_b5_HLT_HcalPhiSym);
   fChain->SetBranchAddress("b5_HLT_HcalIsolatedbunch", &b5_HLT_HcalIsolatedbunch, &b_b5_HLT_HcalIsolatedbunch);
   fChain->SetBranchAddress("b5_HLT_IsoTrackHB", &b5_HLT_IsoTrackHB, &b_b5_HLT_IsoTrackHB);
   fChain->SetBranchAddress("b5_HLT_IsoTrackHE", &b5_HLT_IsoTrackHE, &b_b5_HLT_IsoTrackHE);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_FirstCollisionAfterAbortGap", &b5_HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_b5_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_IsolatedBunches", &b5_HLT_ZeroBias_IsolatedBunches, &b_b5_HLT_ZeroBias_IsolatedBunches);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_FirstCollisionInTrain", &b5_HLT_ZeroBias_FirstCollisionInTrain, &b_b5_HLT_ZeroBias_FirstCollisionInTrain);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_LastCollisionInTrain", &b5_HLT_ZeroBias_LastCollisionInTrain, &b_b5_HLT_ZeroBias_LastCollisionInTrain);
   fChain->SetBranchAddress("b5_HLT_ZeroBias_FirstBXAfterTrain", &b5_HLT_ZeroBias_FirstBXAfterTrain, &b_b5_HLT_ZeroBias_FirstBXAfterTrain);
   fChain->SetBranchAddress("b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1", &b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1, &b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1", &b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1, &b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1", &b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1, &b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1);
   fChain->SetBranchAddress("b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1", &b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1", &b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1", &b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1, &b_b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1);
   fChain->SetBranchAddress("b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg", &b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg", &b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg", &b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg, &b_b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg", &b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg, &b_b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg", &b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, &b_b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg", &b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, &b_b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90", &b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90, &b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", &b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100, &b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", &b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110, &b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", &b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120, &b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", &b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130, &b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140", &b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140, &b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr", &b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr, &b_b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr", &b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr, &b_b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1", &b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1, &b_b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1", &b5_HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1, &b_b5_HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1);
   fChain->SetBranchAddress("b5_HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1", &b5_HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1, &b_b5_HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1);
   fChain->SetBranchAddress("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1", &b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1, &b_b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1", &b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1, &b_b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("b5_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &b5_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_b5_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("b5_HLT_Rsq0p35", &b5_HLT_Rsq0p35, &b_b5_HLT_Rsq0p35);
   fChain->SetBranchAddress("b5_HLT_Rsq0p40", &b5_HLT_Rsq0p40, &b_b5_HLT_Rsq0p40);
   fChain->SetBranchAddress("b5_HLT_RsqMR300_Rsq0p09_MR200", &b5_HLT_RsqMR300_Rsq0p09_MR200, &b_b5_HLT_RsqMR300_Rsq0p09_MR200);
   fChain->SetBranchAddress("b5_HLT_RsqMR320_Rsq0p09_MR200", &b5_HLT_RsqMR320_Rsq0p09_MR200, &b_b5_HLT_RsqMR320_Rsq0p09_MR200);
   fChain->SetBranchAddress("b5_HLT_RsqMR300_Rsq0p09_MR200_4jet", &b5_HLT_RsqMR300_Rsq0p09_MR200_4jet, &b_b5_HLT_RsqMR300_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("b5_HLT_RsqMR320_Rsq0p09_MR200_4jet", &b5_HLT_RsqMR320_Rsq0p09_MR200_4jet, &b_b5_HLT_RsqMR320_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("b5_HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1", &b5_HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1, &b_b5_HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1", &b5_HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1, &b_b5_HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1", &b5_HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1, &b_b5_HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu27_MET90", &b5_HLT_IsoMu27_MET90, &b_b5_HLT_IsoMu27_MET90);
   fChain->SetBranchAddress("b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1", &b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1, &b_b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1);
   fChain->SetBranchAddress("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1", &b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1, &b_b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1);
   fChain->SetBranchAddress("b5_HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", &b5_HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg, &b_b5_HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("b5_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50", &b5_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50, &b_b5_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50);
   fChain->SetBranchAddress("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   fChain->SetBranchAddress("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3", &b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3, &b_b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3);
   fChain->SetBranchAddress("b5_HLT_PFMET100_PFMHT100_IDTight_PFHT60", &b5_HLT_PFMET100_PFMHT100_IDTight_PFHT60, &b_b5_HLT_PFMET100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("b5_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &b5_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60, &b_b5_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
   fChain->SetBranchAddress("b5_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60", &b5_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60, &b_b5_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60);
   fChain->SetBranchAddress("b5_HLT_Mu18_Mu9_SameSign", &b5_HLT_Mu18_Mu9_SameSign, &b_b5_HLT_Mu18_Mu9_SameSign);
   fChain->SetBranchAddress("b5_HLT_Mu18_Mu9_SameSign_DZ", &b5_HLT_Mu18_Mu9_SameSign_DZ, &b_b5_HLT_Mu18_Mu9_SameSign_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu18_Mu9", &b5_HLT_Mu18_Mu9, &b_b5_HLT_Mu18_Mu9);
   fChain->SetBranchAddress("b5_HLT_Mu18_Mu9_DZ", &b5_HLT_Mu18_Mu9_DZ, &b_b5_HLT_Mu18_Mu9_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu20_Mu10_SameSign", &b5_HLT_Mu20_Mu10_SameSign, &b_b5_HLT_Mu20_Mu10_SameSign);
   fChain->SetBranchAddress("b5_HLT_Mu20_Mu10_SameSign_DZ", &b5_HLT_Mu20_Mu10_SameSign_DZ, &b_b5_HLT_Mu20_Mu10_SameSign_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu20_Mu10", &b5_HLT_Mu20_Mu10, &b_b5_HLT_Mu20_Mu10);
   fChain->SetBranchAddress("b5_HLT_Mu20_Mu10_DZ", &b5_HLT_Mu20_Mu10_DZ, &b_b5_HLT_Mu20_Mu10_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu23_Mu12_SameSign", &b5_HLT_Mu23_Mu12_SameSign, &b_b5_HLT_Mu23_Mu12_SameSign);
   fChain->SetBranchAddress("b5_HLT_Mu23_Mu12_SameSign_DZ", &b5_HLT_Mu23_Mu12_SameSign_DZ, &b_b5_HLT_Mu23_Mu12_SameSign_DZ);
   fChain->SetBranchAddress("b5_HLT_Mu23_Mu12", &b5_HLT_Mu23_Mu12, &b_b5_HLT_Mu23_Mu12);
   fChain->SetBranchAddress("b5_HLT_Mu23_Mu12_DZ", &b5_HLT_Mu23_Mu12_DZ, &b_b5_HLT_Mu23_Mu12_DZ);
   fChain->SetBranchAddress("b5_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", &b5_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05, &b_b5_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05);
   fChain->SetBranchAddress("b5_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", &b5_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi, &b_b5_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi);
   fChain->SetBranchAddress("b5_HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &b5_HLT_DoubleMu3_DCA_PFMET50_PFMHT60, &b_b5_HLT_DoubleMu3_DCA_PFMET50_PFMHT60);
   fChain->SetBranchAddress("b5_HLT_TripleMu_5_3_3_Mass3p8_DCA", &b5_HLT_TripleMu_5_3_3_Mass3p8_DCA, &b_b5_HLT_TripleMu_5_3_3_Mass3p8_DCA);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &b5_HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1, &b_b5_HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &b5_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1, &b_b5_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet105_90_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &b5_HLT_QuadPFJet105_90_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1, &b_b5_HLT_QuadPFJet105_90_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1", &b5_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1, &b_b5_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2", &b5_HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2, &b_b5_HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2", &b5_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2, &b_b5_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2", &b5_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2, &b_b5_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2", &b5_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2, &b_b5_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet98_83_71_15", &b5_HLT_QuadPFJet98_83_71_15, &b_b5_HLT_QuadPFJet98_83_71_15);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet103_88_75_15", &b5_HLT_QuadPFJet103_88_75_15, &b_b5_HLT_QuadPFJet103_88_75_15);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet105_88_76_15", &b5_HLT_QuadPFJet105_88_76_15, &b_b5_HLT_QuadPFJet105_88_76_15);
   fChain->SetBranchAddress("b5_HLT_QuadPFJet111_90_80_15", &b5_HLT_QuadPFJet111_90_80_15, &b_b5_HLT_QuadPFJet111_90_80_15);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17", &b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17, &b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1", &b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1, &b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02", &b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02, &b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2", &b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2, &b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2);
   fChain->SetBranchAddress("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4", &b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4, &b_b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4);
   fChain->SetBranchAddress("b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55", &b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55, &b_b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55);
   fChain->SetBranchAddress("b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto", &b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto, &b_b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto);
   fChain->SetBranchAddress("b5_HLT_Mu8p5_IP3p5_part0", &b5_HLT_Mu8p5_IP3p5_part0, &b_b5_HLT_Mu8p5_IP3p5_part0);
   fChain->SetBranchAddress("b5_HLT_Mu8p5_IP3p5_part1", &b5_HLT_Mu8p5_IP3p5_part1, &b_b5_HLT_Mu8p5_IP3p5_part1);
   fChain->SetBranchAddress("b5_HLT_Mu8p5_IP3p5_part2", &b5_HLT_Mu8p5_IP3p5_part2, &b_b5_HLT_Mu8p5_IP3p5_part2);
   fChain->SetBranchAddress("b5_HLT_Mu8p5_IP3p5_part3", &b5_HLT_Mu8p5_IP3p5_part3, &b_b5_HLT_Mu8p5_IP3p5_part3);
   fChain->SetBranchAddress("b5_HLT_Mu8p5_IP3p5_part4", &b5_HLT_Mu8p5_IP3p5_part4, &b_b5_HLT_Mu8p5_IP3p5_part4);
   fChain->SetBranchAddress("b5_HLT_Mu8p5_IP3p5_part5", &b5_HLT_Mu8p5_IP3p5_part5, &b_b5_HLT_Mu8p5_IP3p5_part5);
   fChain->SetBranchAddress("b5_HLT_Mu10p5_IP3p5_part0", &b5_HLT_Mu10p5_IP3p5_part0, &b_b5_HLT_Mu10p5_IP3p5_part0);
   fChain->SetBranchAddress("b5_HLT_Mu10p5_IP3p5_part1", &b5_HLT_Mu10p5_IP3p5_part1, &b_b5_HLT_Mu10p5_IP3p5_part1);
   fChain->SetBranchAddress("b5_HLT_Mu10p5_IP3p5_part2", &b5_HLT_Mu10p5_IP3p5_part2, &b_b5_HLT_Mu10p5_IP3p5_part2);
   fChain->SetBranchAddress("b5_HLT_Mu10p5_IP3p5_part3", &b5_HLT_Mu10p5_IP3p5_part3, &b_b5_HLT_Mu10p5_IP3p5_part3);
   fChain->SetBranchAddress("b5_HLT_Mu10p5_IP3p5_part4", &b5_HLT_Mu10p5_IP3p5_part4, &b_b5_HLT_Mu10p5_IP3p5_part4);
   fChain->SetBranchAddress("b5_HLT_Mu10p5_IP3p5_part5", &b5_HLT_Mu10p5_IP3p5_part5, &b_b5_HLT_Mu10p5_IP3p5_part5);
   fChain->SetBranchAddress("b5_HLT_Mu9_IP6_part0", &b5_HLT_Mu9_IP6_part0, &b_b5_HLT_Mu9_IP6_part0);
   fChain->SetBranchAddress("b5_HLT_Mu9_IP6_part1", &b5_HLT_Mu9_IP6_part1, &b_b5_HLT_Mu9_IP6_part1);
   fChain->SetBranchAddress("b5_HLT_Mu9_IP6_part2", &b5_HLT_Mu9_IP6_part2, &b_b5_HLT_Mu9_IP6_part2);
   fChain->SetBranchAddress("b5_HLT_Mu9_IP6_part3", &b5_HLT_Mu9_IP6_part3, &b_b5_HLT_Mu9_IP6_part3);
   fChain->SetBranchAddress("b5_HLT_Mu9_IP6_part4", &b5_HLT_Mu9_IP6_part4, &b_b5_HLT_Mu9_IP6_part4);
   fChain->SetBranchAddress("b5_HLT_Mu9_IP6_part5", &b5_HLT_Mu9_IP6_part5, &b_b5_HLT_Mu9_IP6_part5);
   fChain->SetBranchAddress("b5_HLT_Mu8_IP3_part0", &b5_HLT_Mu8_IP3_part0, &b_b5_HLT_Mu8_IP3_part0);
   fChain->SetBranchAddress("b5_HLT_Mu8_IP3_part1", &b5_HLT_Mu8_IP3_part1, &b_b5_HLT_Mu8_IP3_part1);
   fChain->SetBranchAddress("b5_HLT_Mu8_IP3_part2", &b5_HLT_Mu8_IP3_part2, &b_b5_HLT_Mu8_IP3_part2);
   fChain->SetBranchAddress("b5_HLT_Mu8_IP3_part3", &b5_HLT_Mu8_IP3_part3, &b_b5_HLT_Mu8_IP3_part3);
   fChain->SetBranchAddress("b5_HLT_Mu8_IP3_part4", &b5_HLT_Mu8_IP3_part4, &b_b5_HLT_Mu8_IP3_part4);
   fChain->SetBranchAddress("b5_HLT_Mu8_IP3_part5", &b5_HLT_Mu8_IP3_part5, &b_b5_HLT_Mu8_IP3_part5);
   fChain->SetBranchAddress("b5_HLTriggerFinalPath", &b5_HLTriggerFinalPath, &b_b5_HLTriggerFinalPath);
   fChain->SetBranchAddress("bG_run", &bG_run, &b_bG_run);
   fChain->SetBranchAddress("bG_event", &bG_event, &b_bG_event);
   fChain->SetBranchAddress("bG_lumis", &bG_lumis, &b_bG_lumis);
   fChain->SetBranchAddress("bG_isData", &bG_isData, &b_bG_isData);
   fChain->SetBranchAddress("bG_nSC", &bG_nSC, &b_bG_nSC);
   fChain->SetBranchAddress("bG_nscE", &bG_nscE, &b_bG_nscE);
   fChain->SetBranchAddress("bG_scE", bG_scE, &b_bG_scE);
   fChain->SetBranchAddress("bG_nscEt", &bG_nscEt, &b_bG_nscEt);
   fChain->SetBranchAddress("bG_scEt", bG_scEt, &b_bG_scEt);
   fChain->SetBranchAddress("bG_nscRawE", &bG_nscRawE, &b_bG_nscRawE);
   fChain->SetBranchAddress("bG_scRawE", bG_scRawE, &b_bG_scRawE);
   fChain->SetBranchAddress("bG_nscEta", &bG_nscEta, &b_bG_nscEta);
   fChain->SetBranchAddress("bG_scEta", bG_scEta, &b_bG_scEta);
   fChain->SetBranchAddress("bG_nscPhi", &bG_nscPhi, &b_bG_nscPhi);
   fChain->SetBranchAddress("bG_scPhi", bG_scPhi, &b_bG_scPhi);
   fChain->SetBranchAddress("bG_nscX", &bG_nscX, &b_bG_nscX);
   fChain->SetBranchAddress("bG_scX", bG_scX, &b_bG_scX);
   fChain->SetBranchAddress("bG_nscY", &bG_nscY, &b_bG_nscY);
   fChain->SetBranchAddress("bG_scY", bG_scY, &b_bG_scY);
   fChain->SetBranchAddress("bG_nscZ", &bG_nscZ, &b_bG_nscZ);
   fChain->SetBranchAddress("bG_scZ", bG_scZ, &b_bG_scZ);
   fChain->SetBranchAddress("bG_nscEtaWidth", &bG_nscEtaWidth, &b_bG_nscEtaWidth);
   fChain->SetBranchAddress("bG_scEtaWidth", bG_scEtaWidth, &b_bG_scEtaWidth);
   fChain->SetBranchAddress("bG_nscPhiWidth", &bG_nscPhiWidth, &b_bG_nscPhiWidth);
   fChain->SetBranchAddress("bG_scPhiWidth", bG_scPhiWidth, &b_bG_scPhiWidth);
   fChain->SetBranchAddress("bG_nscRawEt", &bG_nscRawEt, &b_bG_nscRawEt);
   fChain->SetBranchAddress("bG_scRawEt", bG_scRawEt, &b_bG_scRawEt);
   fChain->SetBranchAddress("bG_nscMinDrWithGsfElectornSC_", &bG_nscMinDrWithGsfElectornSC_, &b_bG_nscMinDrWithGsfElectornSC_);
   fChain->SetBranchAddress("bG_scMinDrWithGsfElectornSC_", bG_scMinDrWithGsfElectornSC_, &b_bG_scMinDrWithGsfElectornSC_);
   fChain->SetBranchAddress("bG_nscFoundGsfMatch_", &bG_nscFoundGsfMatch_, &b_bG_nscFoundGsfMatch_);
   fChain->SetBranchAddress("bG_scFoundGsfMatch_", bG_scFoundGsfMatch_, &b_bG_scFoundGsfMatch_);
   fChain->SetBranchAddress("bG_nscE5x5", &bG_nscE5x5, &b_bG_nscE5x5);
   fChain->SetBranchAddress("bG_scE5x5", bG_scE5x5, &b_bG_scE5x5);
   fChain->SetBranchAddress("bG_nscE2x2Ratio", &bG_nscE2x2Ratio, &b_bG_nscE2x2Ratio);
   fChain->SetBranchAddress("bG_scE2x2Ratio", bG_scE2x2Ratio, &b_bG_scE2x2Ratio);
   fChain->SetBranchAddress("bG_nscE3x3Ratio", &bG_nscE3x3Ratio, &b_bG_nscE3x3Ratio);
   fChain->SetBranchAddress("bG_scE3x3Ratio", bG_scE3x3Ratio, &b_bG_scE3x3Ratio);
   fChain->SetBranchAddress("bG_nscEMaxRatio", &bG_nscEMaxRatio, &b_bG_nscEMaxRatio);
   fChain->SetBranchAddress("bG_scEMaxRatio", bG_scEMaxRatio, &b_bG_scEMaxRatio);
   fChain->SetBranchAddress("bG_nscE2ndRatio", &bG_nscE2ndRatio, &b_bG_nscE2ndRatio);
   fChain->SetBranchAddress("bG_scE2ndRatio", bG_scE2ndRatio, &b_bG_scE2ndRatio);
   fChain->SetBranchAddress("bG_nscETopRatio", &bG_nscETopRatio, &b_bG_nscETopRatio);
   fChain->SetBranchAddress("bG_scETopRatio", bG_scETopRatio, &b_bG_scETopRatio);
   fChain->SetBranchAddress("bG_nscERightRatio", &bG_nscERightRatio, &b_bG_nscERightRatio);
   fChain->SetBranchAddress("bG_scERightRatio", bG_scERightRatio, &b_bG_scERightRatio);
   fChain->SetBranchAddress("bG_nscEBottomRatio", &bG_nscEBottomRatio, &b_bG_nscEBottomRatio);
   fChain->SetBranchAddress("bG_scEBottomRatio", bG_scEBottomRatio, &b_bG_scEBottomRatio);
   fChain->SetBranchAddress("bG_nscELeftRatio", &bG_nscELeftRatio, &b_bG_nscELeftRatio);
   fChain->SetBranchAddress("bG_scELeftRatio", bG_scELeftRatio, &b_bG_scELeftRatio);
   fChain->SetBranchAddress("bG_nscE2x5MaxRatio", &bG_nscE2x5MaxRatio, &b_bG_nscE2x5MaxRatio);
   fChain->SetBranchAddress("bG_scE2x5MaxRatio", bG_scE2x5MaxRatio, &b_bG_scE2x5MaxRatio);
   fChain->SetBranchAddress("bG_nscE2x5TopRatio", &bG_nscE2x5TopRatio, &b_bG_nscE2x5TopRatio);
   fChain->SetBranchAddress("bG_scE2x5TopRatio", bG_scE2x5TopRatio, &b_bG_scE2x5TopRatio);
   fChain->SetBranchAddress("bG_nscE2x5RightRatio", &bG_nscE2x5RightRatio, &b_bG_nscE2x5RightRatio);
   fChain->SetBranchAddress("bG_scE2x5RightRatio", bG_scE2x5RightRatio, &b_bG_scE2x5RightRatio);
   fChain->SetBranchAddress("bG_nscE2x5BottomRatio", &bG_nscE2x5BottomRatio, &b_bG_nscE2x5BottomRatio);
   fChain->SetBranchAddress("bG_scE2x5BottomRatio", bG_scE2x5BottomRatio, &b_bG_scE2x5BottomRatio);
   fChain->SetBranchAddress("bG_nscE2x5LeftRatio", &bG_nscE2x5LeftRatio, &b_bG_nscE2x5LeftRatio);
   fChain->SetBranchAddress("bG_scE2x5LeftRatio", bG_scE2x5LeftRatio, &b_bG_scE2x5LeftRatio);
   fChain->SetBranchAddress("bG_nscSwissCross", &bG_nscSwissCross, &b_bG_nscSwissCross);
   fChain->SetBranchAddress("bG_scSwissCross", bG_scSwissCross, &b_bG_scSwissCross);
   fChain->SetBranchAddress("bG_nscR9", &bG_nscR9, &b_bG_nscR9);
   fChain->SetBranchAddress("bG_scR9", bG_scR9, &b_bG_scR9);
   fChain->SetBranchAddress("bG_nscSigmaIetaIeta", &bG_nscSigmaIetaIeta, &b_bG_nscSigmaIetaIeta);
   fChain->SetBranchAddress("bG_scSigmaIetaIeta", bG_scSigmaIetaIeta, &b_bG_scSigmaIetaIeta);
   fChain->SetBranchAddress("bG_nscSigmaIetaIphi", &bG_nscSigmaIetaIphi, &b_bG_nscSigmaIetaIphi);
   fChain->SetBranchAddress("bG_scSigmaIetaIphi", bG_scSigmaIetaIphi, &b_bG_scSigmaIetaIphi);
   fChain->SetBranchAddress("bG_nscSigmaIphiIphi", &bG_nscSigmaIphiIphi, &b_bG_nscSigmaIphiIphi);
   fChain->SetBranchAddress("bG_scSigmaIphiIphi", bG_scSigmaIphiIphi, &b_bG_scSigmaIphiIphi);
   fChain->SetBranchAddress("bG_nscFull5x5_e5x5", &bG_nscFull5x5_e5x5, &b_bG_nscFull5x5_e5x5);
   fChain->SetBranchAddress("bG_scFull5x5_e5x5", bG_scFull5x5_e5x5, &b_bG_scFull5x5_e5x5);
   fChain->SetBranchAddress("bG_nscFull5x5_e2x2Ratio", &bG_nscFull5x5_e2x2Ratio, &b_bG_nscFull5x5_e2x2Ratio);
   fChain->SetBranchAddress("bG_scFull5x5_e2x2Ratio", bG_scFull5x5_e2x2Ratio, &b_bG_scFull5x5_e2x2Ratio);
   fChain->SetBranchAddress("bG_nscFull5x5_e3x3Ratio", &bG_nscFull5x5_e3x3Ratio, &b_bG_nscFull5x5_e3x3Ratio);
   fChain->SetBranchAddress("bG_scFull5x5_e3x3Ratio", bG_scFull5x5_e3x3Ratio, &b_bG_scFull5x5_e3x3Ratio);
   fChain->SetBranchAddress("bG_nscFull5x5_eMaxRatio", &bG_nscFull5x5_eMaxRatio, &b_bG_nscFull5x5_eMaxRatio);
   fChain->SetBranchAddress("bG_scFull5x5_eMaxRatio", bG_scFull5x5_eMaxRatio, &b_bG_scFull5x5_eMaxRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_e2ndRatio", &bG_nscFull5x5_e2ndRatio, &b_bG_nscFull5x5_e2ndRatio);
   fChain->SetBranchAddress("bG_scFull5x5_e2ndRatio", bG_scFull5x5_e2ndRatio, &b_bG_scFull5x5_e2ndRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_eTopRatio", &bG_nscFull5x5_eTopRatio, &b_bG_nscFull5x5_eTopRatio);
   fChain->SetBranchAddress("bG_scFull5x5_eTopRatio", bG_scFull5x5_eTopRatio, &b_bG_scFull5x5_eTopRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_eRightRatio", &bG_nscFull5x5_eRightRatio, &b_bG_nscFull5x5_eRightRatio);
   fChain->SetBranchAddress("bG_scFull5x5_eRightRatio", bG_scFull5x5_eRightRatio, &b_bG_scFull5x5_eRightRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_eBottomRatio", &bG_nscFull5x5_eBottomRatio, &b_bG_nscFull5x5_eBottomRatio);
   fChain->SetBranchAddress("bG_scFull5x5_eBottomRatio", bG_scFull5x5_eBottomRatio, &b_bG_scFull5x5_eBottomRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_eLeftRatio", &bG_nscFull5x5_eLeftRatio, &b_bG_nscFull5x5_eLeftRatio);
   fChain->SetBranchAddress("bG_scFull5x5_eLeftRatio", bG_scFull5x5_eLeftRatio, &b_bG_scFull5x5_eLeftRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_e2x5MaxRatio", &bG_nscFull5x5_e2x5MaxRatio, &b_bG_nscFull5x5_e2x5MaxRatio);
   fChain->SetBranchAddress("bG_scFull5x5_e2x5MaxRatio", bG_scFull5x5_e2x5MaxRatio, &b_bG_scFull5x5_e2x5MaxRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_e2x5TopRatio", &bG_nscFull5x5_e2x5TopRatio, &b_bG_nscFull5x5_e2x5TopRatio);
   fChain->SetBranchAddress("bG_scFull5x5_e2x5TopRatio", bG_scFull5x5_e2x5TopRatio, &b_bG_scFull5x5_e2x5TopRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_e2x5RightRatio", &bG_nscFull5x5_e2x5RightRatio, &b_bG_nscFull5x5_e2x5RightRatio);
   fChain->SetBranchAddress("bG_scFull5x5_e2x5RightRatio", bG_scFull5x5_e2x5RightRatio, &b_bG_scFull5x5_e2x5RightRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_e2x5BottomRatio", &bG_nscFull5x5_e2x5BottomRatio, &b_bG_nscFull5x5_e2x5BottomRatio);
   fChain->SetBranchAddress("bG_scFull5x5_e2x5BottomRatio", bG_scFull5x5_e2x5BottomRatio, &b_bG_scFull5x5_e2x5BottomRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_e2x5LeftRatio", &bG_nscFull5x5_e2x5LeftRatio, &b_bG_nscFull5x5_e2x5LeftRatio);
   fChain->SetBranchAddress("bG_scFull5x5_e2x5LeftRatio", bG_scFull5x5_e2x5LeftRatio, &b_bG_scFull5x5_e2x5LeftRatio);
   fChain->SetBranchAddress("bG_nscFull5x5_swissCross", &bG_nscFull5x5_swissCross, &b_bG_nscFull5x5_swissCross);
   fChain->SetBranchAddress("bG_scFull5x5_swissCross", bG_scFull5x5_swissCross, &b_bG_scFull5x5_swissCross);
   fChain->SetBranchAddress("bG_nscFull5x5_r9", &bG_nscFull5x5_r9, &b_bG_nscFull5x5_r9);
   fChain->SetBranchAddress("bG_scFull5x5_r9", bG_scFull5x5_r9, &b_bG_scFull5x5_r9);
   fChain->SetBranchAddress("bG_nscFull5x5_sigmaIetaIeta", &bG_nscFull5x5_sigmaIetaIeta, &b_bG_nscFull5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("bG_scFull5x5_sigmaIetaIeta", bG_scFull5x5_sigmaIetaIeta, &b_bG_scFull5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("bG_nscFull5x5_sigmaIetaIphi", &bG_nscFull5x5_sigmaIetaIphi, &b_bG_nscFull5x5_sigmaIetaIphi);
   fChain->SetBranchAddress("bG_scFull5x5_sigmaIetaIphi", bG_scFull5x5_sigmaIetaIphi, &b_bG_scFull5x5_sigmaIetaIphi);
   fChain->SetBranchAddress("bG_nscFull5x5_sigmaIphiIphi", &bG_nscFull5x5_sigmaIphiIphi, &b_bG_nscFull5x5_sigmaIphiIphi);
   fChain->SetBranchAddress("bG_scFull5x5_sigmaIphiIphi", bG_scFull5x5_sigmaIphiIphi, &b_bG_scFull5x5_sigmaIphiIphi);
   fChain->SetBranchAddress("bG_nscNHcalRecHitInDIEta5IPhi5", &bG_nscNHcalRecHitInDIEta5IPhi5, &b_bG_nscNHcalRecHitInDIEta5IPhi5);
   fChain->SetBranchAddress("bG_scNHcalRecHitInDIEta5IPhi5", bG_scNHcalRecHitInDIEta5IPhi5, &b_bG_scNHcalRecHitInDIEta5IPhi5);
   fChain->SetBranchAddress("bG_nscEFromHcalRecHitInDIEta5IPhi5", &bG_nscEFromHcalRecHitInDIEta5IPhi5, &b_bG_nscEFromHcalRecHitInDIEta5IPhi5);
   fChain->SetBranchAddress("bG_scEFromHcalRecHitInDIEta5IPhi5", bG_scEFromHcalRecHitInDIEta5IPhi5, &b_bG_scEFromHcalRecHitInDIEta5IPhi5);
   fChain->SetBranchAddress("bG_nscNHcalRecHitInDIEta2IPhi2", &bG_nscNHcalRecHitInDIEta2IPhi2, &b_bG_nscNHcalRecHitInDIEta2IPhi2);
   fChain->SetBranchAddress("bG_scNHcalRecHitInDIEta2IPhi2", bG_scNHcalRecHitInDIEta2IPhi2, &b_bG_scNHcalRecHitInDIEta2IPhi2);
   fChain->SetBranchAddress("bG_nscEFromHcalRecHitInDIEta2IPhi2", &bG_nscEFromHcalRecHitInDIEta2IPhi2, &b_bG_nscEFromHcalRecHitInDIEta2IPhi2);
   fChain->SetBranchAddress("bG_scEFromHcalRecHitInDIEta2IPhi2", bG_scEFromHcalRecHitInDIEta2IPhi2, &b_bG_scEFromHcalRecHitInDIEta2IPhi2);
   fChain->SetBranchAddress("bG_nscPFChIso1", &bG_nscPFChIso1, &b_bG_nscPFChIso1);
   fChain->SetBranchAddress("bG_scPFChIso1", bG_scPFChIso1, &b_bG_scPFChIso1);
   fChain->SetBranchAddress("bG_nscPFChIso2", &bG_nscPFChIso2, &b_bG_nscPFChIso2);
   fChain->SetBranchAddress("bG_scPFChIso2", bG_scPFChIso2, &b_bG_scPFChIso2);
   fChain->SetBranchAddress("bG_nscPFChIso3", &bG_nscPFChIso3, &b_bG_nscPFChIso3);
   fChain->SetBranchAddress("bG_scPFChIso3", bG_scPFChIso3, &b_bG_scPFChIso3);
   fChain->SetBranchAddress("bG_nscPFChIso4", &bG_nscPFChIso4, &b_bG_nscPFChIso4);
   fChain->SetBranchAddress("bG_scPFChIso4", bG_scPFChIso4, &b_bG_scPFChIso4);
   fChain->SetBranchAddress("bG_nscPFChIso5", &bG_nscPFChIso5, &b_bG_nscPFChIso5);
   fChain->SetBranchAddress("bG_scPFChIso5", bG_scPFChIso5, &b_bG_scPFChIso5);
   fChain->SetBranchAddress("bG_nscPFPhoIso1", &bG_nscPFPhoIso1, &b_bG_nscPFPhoIso1);
   fChain->SetBranchAddress("bG_scPFPhoIso1", bG_scPFPhoIso1, &b_bG_scPFPhoIso1);
   fChain->SetBranchAddress("bG_nscPFPhoIso2", &bG_nscPFPhoIso2, &b_bG_nscPFPhoIso2);
   fChain->SetBranchAddress("bG_scPFPhoIso2", bG_scPFPhoIso2, &b_bG_scPFPhoIso2);
   fChain->SetBranchAddress("bG_nscPFPhoIso3", &bG_nscPFPhoIso3, &b_bG_nscPFPhoIso3);
   fChain->SetBranchAddress("bG_scPFPhoIso3", bG_scPFPhoIso3, &b_bG_scPFPhoIso3);
   fChain->SetBranchAddress("bG_nscPFPhoIso4", &bG_nscPFPhoIso4, &b_bG_nscPFPhoIso4);
   fChain->SetBranchAddress("bG_scPFPhoIso4", bG_scPFPhoIso4, &b_bG_scPFPhoIso4);
   fChain->SetBranchAddress("bG_nscPFPhoIso5", &bG_nscPFPhoIso5, &b_bG_nscPFPhoIso5);
   fChain->SetBranchAddress("bG_scPFPhoIso5", bG_scPFPhoIso5, &b_bG_scPFPhoIso5);
   fChain->SetBranchAddress("bG_nscPFNeuIso1", &bG_nscPFNeuIso1, &b_bG_nscPFNeuIso1);
   fChain->SetBranchAddress("bG_scPFNeuIso1", bG_scPFNeuIso1, &b_bG_scPFNeuIso1);
   fChain->SetBranchAddress("bG_nscPFNeuIso2", &bG_nscPFNeuIso2, &b_bG_nscPFNeuIso2);
   fChain->SetBranchAddress("bG_scPFNeuIso2", bG_scPFNeuIso2, &b_bG_scPFNeuIso2);
   fChain->SetBranchAddress("bG_nscPFNeuIso3", &bG_nscPFNeuIso3, &b_bG_nscPFNeuIso3);
   fChain->SetBranchAddress("bG_scPFNeuIso3", bG_scPFNeuIso3, &b_bG_scPFNeuIso3);
   fChain->SetBranchAddress("bG_nscPFNeuIso4", &bG_nscPFNeuIso4, &b_bG_nscPFNeuIso4);
   fChain->SetBranchAddress("bG_scPFNeuIso4", bG_scPFNeuIso4, &b_bG_scPFNeuIso4);
   fChain->SetBranchAddress("bG_nscPFNeuIso5", &bG_nscPFNeuIso5, &b_bG_nscPFNeuIso5);
   fChain->SetBranchAddress("bG_scPFNeuIso5", bG_scPFNeuIso5, &b_bG_scPFNeuIso5);
   fChain->SetBranchAddress("bG_nPrimaryVertex", &bG_nPrimaryVertex, &b_bG_nPrimaryVertex);
   fChain->SetBranchAddress("bG_primaryVertex_isFake", bG_primaryVertex_isFake, &b_bG_primaryVertex_isFake);
   fChain->SetBranchAddress("bG_primaryVertex_x", bG_primaryVertex_x, &b_bG_primaryVertex_x);
   fChain->SetBranchAddress("bG_primaryVertex_y", bG_primaryVertex_y, &b_bG_primaryVertex_y);
   fChain->SetBranchAddress("bG_primaryVertex_z", bG_primaryVertex_z, &b_bG_primaryVertex_z);
   fChain->SetBranchAddress("bG_primaryVertex_t", bG_primaryVertex_t, &b_bG_primaryVertex_t);
   fChain->SetBranchAddress("bG_primaryVertex_covXX", bG_primaryVertex_covXX, &b_bG_primaryVertex_covXX);
   fChain->SetBranchAddress("bG_primaryVertex_covXY", bG_primaryVertex_covXY, &b_bG_primaryVertex_covXY);
   fChain->SetBranchAddress("bG_primaryVertex_covXZ", bG_primaryVertex_covXZ, &b_bG_primaryVertex_covXZ);
   fChain->SetBranchAddress("bG_primaryVertex_covYY", bG_primaryVertex_covYY, &b_bG_primaryVertex_covYY);
   fChain->SetBranchAddress("bG_primaryVertex_covYZ", bG_primaryVertex_covYZ, &b_bG_primaryVertex_covYZ);
   fChain->SetBranchAddress("bG_primaryVertex_covZZ", bG_primaryVertex_covZZ, &b_bG_primaryVertex_covZZ);
   fChain->SetBranchAddress("bG_primaryVertex_x_error", bG_primaryVertex_x_error, &b_bG_primaryVertex_x_error);
   fChain->SetBranchAddress("bG_primaryVertex_y_error", bG_primaryVertex_y_error, &b_bG_primaryVertex_y_error);
   fChain->SetBranchAddress("bG_primaryVertex_z_error", bG_primaryVertex_z_error, &b_bG_primaryVertex_z_error);
   fChain->SetBranchAddress("bG_primaryVertex_t_error", bG_primaryVertex_t_error, &b_bG_primaryVertex_t_error);
   fChain->SetBranchAddress("bG_primaryVertex_ntracks", bG_primaryVertex_ntracks, &b_bG_primaryVertex_ntracks);
   fChain->SetBranchAddress("bG_primaryVertex_ndof", bG_primaryVertex_ndof, &b_bG_primaryVertex_ndof);
   fChain->SetBranchAddress("bG_primaryVertex_chi2", bG_primaryVertex_chi2, &b_bG_primaryVertex_chi2);
   fChain->SetBranchAddress("bG_primaryVertex_normalizedChi2", bG_primaryVertex_normalizedChi2, &b_bG_primaryVertex_normalizedChi2);
   Notify();
}

Bool_t MergedBMMX2018Data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MergedBMMX2018Data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MergedBMMX2018Data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MergedBMMX2018Data_cxx
