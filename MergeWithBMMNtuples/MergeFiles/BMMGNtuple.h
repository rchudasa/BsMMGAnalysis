//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 12 18:19:40 2021 by ROOT version 6.14/09
// from TTree EventTree/Event data
// found on file: /eos/home-a/athachay/workarea/data/BsToMuMuGamma/Run2Studies/BsToMuMuGammaNtuples/Charmonium/crab_charmonium2018A_ntuplizer/0000/Charmonium_Run2018A12Nov2019_UL2018_rsbv1_hadd1.root
//////////////////////////////////////////////////////////

#ifndef BMMGNtuple_h
#define BMMGNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class BMMGNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxscMinDrWithGsfElectornSC = 1;
   static constexpr Int_t kMaxscFoundGsfMatch = 1;

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       event;
   UInt_t          lumis;
   Bool_t          isData;
   Double_t        beamspot_x;
   Double_t        beamspot_y;
   Double_t        beamspot_z;
   Double_t        beamspot_x_error;
   Double_t        beamspot_y_error;
   Double_t        beamspot_z_error;
   Double_t        beamspot_covXX;
   Double_t        beamspot_covXY;
   Double_t        beamspot_covXZ;
   Double_t        beamspot_covYY;
   Double_t        beamspot_covYZ;
   Double_t        beamspot_covZZ;
   Double_t        beamspot_dxdz;
   Double_t        beamspot_dydz;
   Double_t        beamspot_sigmaZ;
   Double_t        beamspot_dxdz_error;
   Double_t        beamspot_dydz_error;
   Double_t        beamspot_sigmaZError;
   Double_t        beamspot_beamWidthX;
   Double_t        beamspot_beamWidthY;
   Double_t        beamspot_beamWidthX_error;
   Double_t        beamspot_beamWidthY_error;
   Int_t           nPrimaryVertex;
   vector<bool>    *primaryVertex_isFake;
   vector<double>  *primaryVertex_x;
   vector<double>  *primaryVertex_y;
   vector<double>  *primaryVertex_z;
   vector<double>  *primaryVertex_t;
   vector<double>  *primaryVertex_covXX;
   vector<double>  *primaryVertex_covXY;
   vector<double>  *primaryVertex_covXZ;
   vector<double>  *primaryVertex_covYY;
   vector<double>  *primaryVertex_covYZ;
   vector<double>  *primaryVertex_covZZ;
   vector<double>  *primaryVertex_x_error;
   vector<double>  *primaryVertex_y_error;
   vector<double>  *primaryVertex_z_error;
   vector<double>  *primaryVertex_t_error;
   vector<double>  *primaryVertex_ntracks;
   vector<double>  *primaryVertex_ndof;
   vector<double>  *primaryVertex_chi2;
   vector<double>  *primaryVertex_normalizedChi2;
   vector<bool>    *trigResult;
   vector<int>     *trigPrescales;
   vector<string>  *l1Table;
   vector<int>     *l1Prescales;
   vector<bool>    *HLT_DoubleMu4_3_Bs_v14_result;
   vector<int>     *HLT_DoubleMu4_3_Bs_v14_prescale;
   vector<bool>    *HLT_DoubleMu4_3_Jpsi_v2_result;
   vector<int>     *HLT_DoubleMu4_3_Jpsi_v2_prescale;
   vector<bool>    *HLT_DoubleMu4_JpsiTrk_Displaced_v15_result;
   vector<int>     *HLT_DoubleMu4_JpsiTrk_Displaced_v15_prescale;
   vector<bool>    *HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_result;
   vector<int>     *HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_prescale;
   vector<bool>    *HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_result;
   vector<int>     *HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_prescale;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoSigmaE;
   vector<float>   *phoCalibE;
   vector<float>   *phoCalibEt;
   vector<float>   *phoSCE;
   vector<float>   *phoSCEt;
   vector<float>   *phoSCRawE;
   vector<float>   *phoESEnP1;
   vector<float>   *phoESEnP2;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   vector<float>   *phoE2x2Full5x5;
   vector<float>   *phoE5x5Full5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoEcalPFClusterIso;
   vector<float>   *phoHcalPFClusterIso;
   vector<float>   *phoIDMVA;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<float>   *phoMIPTotEnergy;
   vector<float>   *phoMIPChi2;
   vector<float>   *phoMIPSlope;
   vector<float>   *phoMIPIntercept;
   vector<float>   *phoMIPNhitCone;
   vector<float>   *phoMIPIsHalo;
   Int_t           nPFPho;
   vector<float>   *phoPFE;
   vector<float>   *phoPFEt;
   vector<float>   *phoPFEta;
   vector<float>   *phoPFPhi;
   Int_t           nSC;
   vector<float>   *scE;
   vector<float>   *scRawE;
   vector<float>   *scEta;
   vector<float>   *scPhi;
   vector<float>   *scX;
   vector<float>   *scY;
   vector<float>   *scZ;
   vector<float>   *scEtaWidth;
   vector<float>   *scPhiWidth;
   vector<float>   *scRawEt;
   vector<float>   *scMinDrWithGsfElectornSC_;
   vector<bool>    *scFoundGsfMatch_;
   vector<float>   *superCluster_e5x5;
   vector<float>   *superCluster_e2x2Ratio;
   vector<float>   *superCluster_e3x3Ratio;
   vector<float>   *superCluster_eMaxRatio;
   vector<float>   *superCluster_e2ndRatio;
   vector<float>   *superCluster_eTopRatio;
   vector<float>   *superCluster_eRightRatio;
   vector<float>   *superCluster_eBottomRatio;
   vector<float>   *superCluster_eLeftRatio;
   vector<float>   *superCluster_e2x5MaxRatio;
   vector<float>   *superCluster_e2x5TopRatio;
   vector<float>   *superCluster_e2x5RightRatio;
   vector<float>   *superCluster_e2x5BottomRatio;
   vector<float>   *superCluster_e2x5LeftRatio;
   vector<float>   *superCluster_swissCross;
   vector<float>   *superCluster_r9;
   vector<float>   *superCluster_sigmaIetaIeta;
   vector<float>   *superCluster_sigmaIetaIphi;
   vector<float>   *superCluster_sigmaIphiIphi;
   vector<float>   *superCluster_full5x5_e5x5;
   vector<float>   *superCluster_full5x5_e2x2Ratio;
   vector<float>   *superCluster_full5x5_e3x3Ratio;
   vector<float>   *superCluster_full5x5_eMaxRatio;
   vector<float>   *superCluster_full5x5_e2ndRatio;
   vector<float>   *superCluster_full5x5_eTopRatio;
   vector<float>   *superCluster_full5x5_eRightRatio;
   vector<float>   *superCluster_full5x5_eBottomRatio;
   vector<float>   *superCluster_full5x5_eLeftRatio;
   vector<float>   *superCluster_full5x5_e2x5MaxRatio;
   vector<float>   *superCluster_full5x5_e2x5TopRatio;
   vector<float>   *superCluster_full5x5_e2x5RightRatio;
   vector<float>   *superCluster_full5x5_e2x5BottomRatio;
   vector<float>   *superCluster_full5x5_e2x5LeftRatio;
   vector<float>   *superCluster_full5x5_swissCross;
   vector<float>   *superCluster_full5x5_r9;
   vector<float>   *superCluster_full5x5_sigmaIetaIeta;
   vector<float>   *superCluster_full5x5_sigmaIetaIphi;
   vector<float>   *superCluster_full5x5_sigmaIphiIphi;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_beamspot_x;   //!
   TBranch        *b_beamspot_y;   //!
   TBranch        *b_beamspot_z;   //!
   TBranch        *b_beamspot_x_error;   //!
   TBranch        *b_beamspot_y_error;   //!
   TBranch        *b_beamspot_z_error;   //!
   TBranch        *b_beamspot_covXX;   //!
   TBranch        *b_beamspot_covXY;   //!
   TBranch        *b_beamspot_covXZ;   //!
   TBranch        *b_beamspot_covYY;   //!
   TBranch        *b_beamspot_covYZ;   //!
   TBranch        *b_beamspot_covZZ;   //!
   TBranch        *b_beamspot_dxdz;   //!
   TBranch        *b_beamspot_dydz;   //!
   TBranch        *b_beamspot_sigmaZ;   //!
   TBranch        *b_beamspot_dxdz_error;   //!
   TBranch        *b_beamspot_dydz_error;   //!
   TBranch        *b_beamspot_sigmaZError;   //!
   TBranch        *b_beamspot_beamWidthX;   //!
   TBranch        *b_beamspot_beamWidthY;   //!
   TBranch        *b_beamspot_beamWidthX_error;   //!
   TBranch        *b_beamspot_beamWidthY_error;   //!
   TBranch        *b_nPrimaryVertex;   //!
   TBranch        *b_primaryVertex_isFake;   //!
   TBranch        *b_primaryVertex_x;   //!
   TBranch        *b_primaryVertex_y;   //!
   TBranch        *b_primaryVertex_z;   //!
   TBranch        *b_primaryVertex_t;   //!
   TBranch        *b_primaryVertex_covXX;   //!
   TBranch        *b_primaryVertex_covXY;   //!
   TBranch        *b_primaryVertex_covXZ;   //!
   TBranch        *b_primaryVertex_covYY;   //!
   TBranch        *b_primaryVertex_covYZ;   //!
   TBranch        *b_primaryVertex_covZZ;   //!
   TBranch        *b_primaryVertex_x_error;   //!
   TBranch        *b_primaryVertex_y_error;   //!
   TBranch        *b_primaryVertex_z_error;   //!
   TBranch        *b_primaryVertex_t_error;   //!
   TBranch        *b_primaryVertex_ntracks;   //!
   TBranch        *b_primaryVertex_ndof;   //!
   TBranch        *b_primaryVertex_chi2;   //!
   TBranch        *b_primaryVertex_normalizedChi2;   //!
   TBranch        *b_trigResult;   //!
   TBranch        *b_trigPrescales;   //!
   TBranch        *b_l1Table;   //!
   TBranch        *b_l1Prescales;   //!
   TBranch        *b_HLT_DoubleMu4_3_Bs_v14_result;   //!
   TBranch        *b_HLT_DoubleMu4_3_Bs_v14_prescale;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi_v2_result;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi_v2_prescale;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrk_Displaced_v15_result;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrk_Displaced_v15_prescale;   //!
   TBranch        *b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_result;   //!
   TBranch        *b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_prescale;   //!
   TBranch        *b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_result;   //!
   TBranch        *b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_prescale;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSigmaE;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoCalibEt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_phoE2x2Full5x5;   //!
   TBranch        *b_phoE5x5Full5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoEcalPFClusterIso;   //!
   TBranch        *b_phoHcalPFClusterIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoMIPTotEnergy;   //!
   TBranch        *b_phoMIPChi2;   //!
   TBranch        *b_phoMIPSlope;   //!
   TBranch        *b_phoMIPIntercept;   //!
   TBranch        *b_phoMIPNhitCone;   //!
   TBranch        *b_phoMIPIsHalo;   //!
   TBranch        *b_nPFPho;   //!
   TBranch        *b_phoPFE;   //!
   TBranch        *b_phoPFEt;   //!
   TBranch        *b_phoPFEta;   //!
   TBranch        *b_phoPFPhi;   //!
   TBranch        *b_nSC;   //!
   TBranch        *b_scE;   //!
   TBranch        *b_scRawE;   //!
   TBranch        *b_scEta;   //!
   TBranch        *b_scPhi;   //!
   TBranch        *b_scX;   //!
   TBranch        *b_scY;   //!
   TBranch        *b_scZ;   //!
   TBranch        *b_scEtaWidth;   //!
   TBranch        *b_scPhiWidth;   //!
   TBranch        *b_scRawEt;   //!
   TBranch        *b_scMinDrWithGsfElectornSC_;   //!
   TBranch        *b_scFoundGsfMatch_;   //!
   TBranch        *b_superCluster_e5x5;   //!
   TBranch        *b_superCluster_e2x2Ratio;   //!
   TBranch        *b_superCluster_e3x3Ratio;   //!
   TBranch        *b_superCluster_eMaxRatio;   //!
   TBranch        *b_superCluster_e2ndRatio;   //!
   TBranch        *b_superCluster_eTopRatio;   //!
   TBranch        *b_superCluster_eRightRatio;   //!
   TBranch        *b_superCluster_eBottomRatio;   //!
   TBranch        *b_superCluster_eLeftRatio;   //!
   TBranch        *b_superCluster_e2x5MaxRatio;   //!
   TBranch        *b_superCluster_e2x5TopRatio;   //!
   TBranch        *b_superCluster_e2x5RightRatio;   //!
   TBranch        *b_superCluster_e2x5BottomRatio;   //!
   TBranch        *b_superCluster_e2x5LeftRatio;   //!
   TBranch        *b_superCluster_swissCross;   //!
   TBranch        *b_superCluster_r9;   //!
   TBranch        *b_superCluster_sigmaIetaIeta;   //!
   TBranch        *b_superCluster_sigmaIetaIphi;   //!
   TBranch        *b_superCluster_sigmaIphiIphi;   //!
   TBranch        *b_superCluster_full5x5_e5x5;   //!
   TBranch        *b_superCluster_full5x5_e2x2Ratio;   //!
   TBranch        *b_superCluster_full5x5_e3x3Ratio;   //!
   TBranch        *b_superCluster_full5x5_eMaxRatio;   //!
   TBranch        *b_superCluster_full5x5_e2ndRatio;   //!
   TBranch        *b_superCluster_full5x5_eTopRatio;   //!
   TBranch        *b_superCluster_full5x5_eRightRatio;   //!
   TBranch        *b_superCluster_full5x5_eBottomRatio;   //!
   TBranch        *b_superCluster_full5x5_eLeftRatio;   //!
   TBranch        *b_superCluster_full5x5_e2x5MaxRatio;   //!
   TBranch        *b_superCluster_full5x5_e2x5TopRatio;   //!
   TBranch        *b_superCluster_full5x5_e2x5RightRatio;   //!
   TBranch        *b_superCluster_full5x5_e2x5BottomRatio;   //!
   TBranch        *b_superCluster_full5x5_e2x5LeftRatio;   //!
   TBranch        *b_superCluster_full5x5_swissCross;   //!
   TBranch        *b_superCluster_full5x5_r9;   //!
   TBranch        *b_superCluster_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_superCluster_full5x5_sigmaIetaIphi;   //!
   TBranch        *b_superCluster_full5x5_sigmaIphiIphi;   //!

   BMMGNtuple(TTree *tree=0);
   virtual ~BMMGNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
  // virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BMMGNtuple_cxx
BMMGNtuple::BMMGNtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

   }
   Init(tree);
}

BMMGNtuple::~BMMGNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BMMGNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BMMGNtuple::LoadTree(Long64_t entry)
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

void BMMGNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   primaryVertex_isFake = 0;
   primaryVertex_x = 0;
   primaryVertex_y = 0;
   primaryVertex_z = 0;
   primaryVertex_t = 0;
   primaryVertex_covXX = 0;
   primaryVertex_covXY = 0;
   primaryVertex_covXZ = 0;
   primaryVertex_covYY = 0;
   primaryVertex_covYZ = 0;
   primaryVertex_covZZ = 0;
   primaryVertex_x_error = 0;
   primaryVertex_y_error = 0;
   primaryVertex_z_error = 0;
   primaryVertex_t_error = 0;
   primaryVertex_ntracks = 0;
   primaryVertex_ndof = 0;
   primaryVertex_chi2 = 0;
   primaryVertex_normalizedChi2 = 0;
   trigResult = 0;
   trigPrescales = 0;
   l1Table = 0;
   l1Prescales = 0;
   HLT_DoubleMu4_3_Bs_v14_result = 0;
   HLT_DoubleMu4_3_Bs_v14_prescale = 0;
   HLT_DoubleMu4_3_Jpsi_v2_result = 0;
   HLT_DoubleMu4_3_Jpsi_v2_prescale = 0;
   HLT_DoubleMu4_JpsiTrk_Displaced_v15_result = 0;
   HLT_DoubleMu4_JpsiTrk_Displaced_v15_prescale = 0;
   HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_result = 0;
   HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_prescale = 0;
   HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_result = 0;
   HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_prescale = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoSigmaE = 0;
   phoCalibE = 0;
   phoCalibEt = 0;
   phoSCE = 0;
   phoSCEt = 0;
   phoSCRawE = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoESEffSigmaRR = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoE2x2Full5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoEcalPFClusterIso = 0;
   phoHcalPFClusterIso = 0;
   phoIDMVA = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoMIPTotEnergy = 0;
   phoMIPChi2 = 0;
   phoMIPSlope = 0;
   phoMIPIntercept = 0;
   phoMIPNhitCone = 0;
   phoMIPIsHalo = 0;
   phoPFE = 0;
   phoPFEt = 0;
   phoPFEta = 0;
   phoPFPhi = 0;
   scE = 0;
   scRawE = 0;
   scEta = 0;
   scPhi = 0;
   scX = 0;
   scY = 0;
   scZ = 0;
   scEtaWidth = 0;
   scPhiWidth = 0;
   scRawEt = 0;
   scMinDrWithGsfElectornSC_ = 0;
   scFoundGsfMatch_ = 0;
   superCluster_e5x5 = 0;
   superCluster_e2x2Ratio = 0;
   superCluster_e3x3Ratio = 0;
   superCluster_eMaxRatio = 0;
   superCluster_e2ndRatio = 0;
   superCluster_eTopRatio = 0;
   superCluster_eRightRatio = 0;
   superCluster_eBottomRatio = 0;
   superCluster_eLeftRatio = 0;
   superCluster_e2x5MaxRatio = 0;
   superCluster_e2x5TopRatio = 0;
   superCluster_e2x5RightRatio = 0;
   superCluster_e2x5BottomRatio = 0;
   superCluster_e2x5LeftRatio = 0;
   superCluster_swissCross = 0;
   superCluster_r9 = 0;
   superCluster_sigmaIetaIeta = 0;
   superCluster_sigmaIetaIphi = 0;
   superCluster_sigmaIphiIphi = 0;
   superCluster_full5x5_e5x5 = 0;
   superCluster_full5x5_e2x2Ratio = 0;
   superCluster_full5x5_e3x3Ratio = 0;
   superCluster_full5x5_eMaxRatio = 0;
   superCluster_full5x5_e2ndRatio = 0;
   superCluster_full5x5_eTopRatio = 0;
   superCluster_full5x5_eRightRatio = 0;
   superCluster_full5x5_eBottomRatio = 0;
   superCluster_full5x5_eLeftRatio = 0;
   superCluster_full5x5_e2x5MaxRatio = 0;
   superCluster_full5x5_e2x5TopRatio = 0;
   superCluster_full5x5_e2x5RightRatio = 0;
   superCluster_full5x5_e2x5BottomRatio = 0;
   superCluster_full5x5_e2x5LeftRatio = 0;
   superCluster_full5x5_swissCross = 0;
   superCluster_full5x5_r9 = 0;
   superCluster_full5x5_sigmaIetaIeta = 0;
   superCluster_full5x5_sigmaIetaIphi = 0;
   superCluster_full5x5_sigmaIphiIphi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("beamspot_x", &beamspot_x, &b_beamspot_x);
   fChain->SetBranchAddress("beamspot_y", &beamspot_y, &b_beamspot_y);
   fChain->SetBranchAddress("beamspot_z", &beamspot_z, &b_beamspot_z);
   fChain->SetBranchAddress("beamspot_x_error", &beamspot_x_error, &b_beamspot_x_error);
   fChain->SetBranchAddress("beamspot_y_error", &beamspot_y_error, &b_beamspot_y_error);
   fChain->SetBranchAddress("beamspot_z_error", &beamspot_z_error, &b_beamspot_z_error);
   fChain->SetBranchAddress("beamspot_covXX", &beamspot_covXX, &b_beamspot_covXX);
   fChain->SetBranchAddress("beamspot_covXY", &beamspot_covXY, &b_beamspot_covXY);
   fChain->SetBranchAddress("beamspot_covXZ", &beamspot_covXZ, &b_beamspot_covXZ);
   fChain->SetBranchAddress("beamspot_covYY", &beamspot_covYY, &b_beamspot_covYY);
   fChain->SetBranchAddress("beamspot_covYZ", &beamspot_covYZ, &b_beamspot_covYZ);
   fChain->SetBranchAddress("beamspot_covZZ", &beamspot_covZZ, &b_beamspot_covZZ);
   fChain->SetBranchAddress("beamspot_dxdz", &beamspot_dxdz, &b_beamspot_dxdz);
   fChain->SetBranchAddress("beamspot_dydz", &beamspot_dydz, &b_beamspot_dydz);
   fChain->SetBranchAddress("beamspot_sigmaZ", &beamspot_sigmaZ, &b_beamspot_sigmaZ);
   fChain->SetBranchAddress("beamspot_dxdz_error", &beamspot_dxdz_error, &b_beamspot_dxdz_error);
   fChain->SetBranchAddress("beamspot_dydz_error", &beamspot_dydz_error, &b_beamspot_dydz_error);
   fChain->SetBranchAddress("beamspot_sigmaZError", &beamspot_sigmaZError, &b_beamspot_sigmaZError);
   fChain->SetBranchAddress("beamspot_beamWidthX", &beamspot_beamWidthX, &b_beamspot_beamWidthX);
   fChain->SetBranchAddress("beamspot_beamWidthY", &beamspot_beamWidthY, &b_beamspot_beamWidthY);
   fChain->SetBranchAddress("beamspot_beamWidthX_error", &beamspot_beamWidthX_error, &b_beamspot_beamWidthX_error);
   fChain->SetBranchAddress("beamspot_beamWidthY_error", &beamspot_beamWidthY_error, &b_beamspot_beamWidthY_error);
   fChain->SetBranchAddress("nPrimaryVertex", &nPrimaryVertex, &b_nPrimaryVertex);
   fChain->SetBranchAddress("primaryVertex_isFake", &primaryVertex_isFake, &b_primaryVertex_isFake);
   fChain->SetBranchAddress("primaryVertex_x", &primaryVertex_x, &b_primaryVertex_x);
   fChain->SetBranchAddress("primaryVertex_y", &primaryVertex_y, &b_primaryVertex_y);
   fChain->SetBranchAddress("primaryVertex_z", &primaryVertex_z, &b_primaryVertex_z);
   fChain->SetBranchAddress("primaryVertex_t", &primaryVertex_t, &b_primaryVertex_t);
   fChain->SetBranchAddress("primaryVertex_covXX", &primaryVertex_covXX, &b_primaryVertex_covXX);
   fChain->SetBranchAddress("primaryVertex_covXY", &primaryVertex_covXY, &b_primaryVertex_covXY);
   fChain->SetBranchAddress("primaryVertex_covXZ", &primaryVertex_covXZ, &b_primaryVertex_covXZ);
   fChain->SetBranchAddress("primaryVertex_covYY", &primaryVertex_covYY, &b_primaryVertex_covYY);
   fChain->SetBranchAddress("primaryVertex_covYZ", &primaryVertex_covYZ, &b_primaryVertex_covYZ);
   fChain->SetBranchAddress("primaryVertex_covZZ", &primaryVertex_covZZ, &b_primaryVertex_covZZ);
   fChain->SetBranchAddress("primaryVertex_x_error", &primaryVertex_x_error, &b_primaryVertex_x_error);
   fChain->SetBranchAddress("primaryVertex_y_error", &primaryVertex_y_error, &b_primaryVertex_y_error);
   fChain->SetBranchAddress("primaryVertex_z_error", &primaryVertex_z_error, &b_primaryVertex_z_error);
   fChain->SetBranchAddress("primaryVertex_t_error", &primaryVertex_t_error, &b_primaryVertex_t_error);
   fChain->SetBranchAddress("primaryVertex_ntracks", &primaryVertex_ntracks, &b_primaryVertex_ntracks);
   fChain->SetBranchAddress("primaryVertex_ndof", &primaryVertex_ndof, &b_primaryVertex_ndof);
   fChain->SetBranchAddress("primaryVertex_chi2", &primaryVertex_chi2, &b_primaryVertex_chi2);
   fChain->SetBranchAddress("primaryVertex_normalizedChi2", &primaryVertex_normalizedChi2, &b_primaryVertex_normalizedChi2);
   fChain->SetBranchAddress("trigResult", &trigResult, &b_trigResult);
   fChain->SetBranchAddress("trigPrescales", &trigPrescales, &b_trigPrescales);
   fChain->SetBranchAddress("l1Table", &l1Table, &b_l1Table);
   fChain->SetBranchAddress("l1Prescales", &l1Prescales, &b_l1Prescales);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs_v14_result", &HLT_DoubleMu4_3_Bs_v14_result, &b_HLT_DoubleMu4_3_Bs_v14_result);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs_v14_prescale", &HLT_DoubleMu4_3_Bs_v14_prescale, &b_HLT_DoubleMu4_3_Bs_v14_prescale);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_v2_result", &HLT_DoubleMu4_3_Jpsi_v2_result, &b_HLT_DoubleMu4_3_Jpsi_v2_result);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_v2_prescale", &HLT_DoubleMu4_3_Jpsi_v2_prescale, &b_HLT_DoubleMu4_3_Jpsi_v2_prescale);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced_v15_result", &HLT_DoubleMu4_JpsiTrk_Displaced_v15_result, &b_HLT_DoubleMu4_JpsiTrk_Displaced_v15_result);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced_v15_prescale", &HLT_DoubleMu4_JpsiTrk_Displaced_v15_prescale, &b_HLT_DoubleMu4_JpsiTrk_Displaced_v15_prescale);
   fChain->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_result", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_result, &b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_result);
   fChain->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_prescale", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_prescale, &b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15_prescale);
   fChain->SetBranchAddress("HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_result", &HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_result, &b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_result);
   fChain->SetBranchAddress("HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_prescale", &HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_prescale, &b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8_prescale);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSigmaE", &phoSigmaE, &b_phoSigmaE);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCEt", &phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoEcalPFClusterIso", &phoEcalPFClusterIso, &b_phoEcalPFClusterIso);
   fChain->SetBranchAddress("phoHcalPFClusterIso", &phoHcalPFClusterIso, &b_phoHcalPFClusterIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoMIPTotEnergy", &phoMIPTotEnergy, &b_phoMIPTotEnergy);
   fChain->SetBranchAddress("phoMIPChi2", &phoMIPChi2, &b_phoMIPChi2);
   fChain->SetBranchAddress("phoMIPSlope", &phoMIPSlope, &b_phoMIPSlope);
   fChain->SetBranchAddress("phoMIPIntercept", &phoMIPIntercept, &b_phoMIPIntercept);
   fChain->SetBranchAddress("phoMIPNhitCone", &phoMIPNhitCone, &b_phoMIPNhitCone);
   fChain->SetBranchAddress("phoMIPIsHalo", &phoMIPIsHalo, &b_phoMIPIsHalo);
   fChain->SetBranchAddress("nPFPho", &nPFPho, &b_nPFPho);
   fChain->SetBranchAddress("phoPFE", &phoPFE, &b_phoPFE);
   fChain->SetBranchAddress("phoPFEt", &phoPFEt, &b_phoPFEt);
   fChain->SetBranchAddress("phoPFEta", &phoPFEta, &b_phoPFEta);
   fChain->SetBranchAddress("phoPFPhi", &phoPFPhi, &b_phoPFPhi);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("scE", &scE, &b_scE);
   fChain->SetBranchAddress("scRawE", &scRawE, &b_scRawE);
   fChain->SetBranchAddress("scEta", &scEta, &b_scEta);
   fChain->SetBranchAddress("scPhi", &scPhi, &b_scPhi);
   fChain->SetBranchAddress("scX", &scX, &b_scX);
   fChain->SetBranchAddress("scY", &scY, &b_scY);
   fChain->SetBranchAddress("scZ", &scZ, &b_scZ);
   fChain->SetBranchAddress("scEtaWidth", &scEtaWidth, &b_scEtaWidth);
   fChain->SetBranchAddress("scPhiWidth", &scPhiWidth, &b_scPhiWidth);
   fChain->SetBranchAddress("scRawEt", &scRawEt, &b_scRawEt);
   fChain->SetBranchAddress("scMinDrWithGsfElectornSC_", &scMinDrWithGsfElectornSC_, &b_scMinDrWithGsfElectornSC_);
   fChain->SetBranchAddress("scFoundGsfMatch_", &scFoundGsfMatch_, &b_scFoundGsfMatch_);
   fChain->SetBranchAddress("superCluster_e5x5", &superCluster_e5x5, &b_superCluster_e5x5);
   fChain->SetBranchAddress("superCluster_e2x2Ratio", &superCluster_e2x2Ratio, &b_superCluster_e2x2Ratio);
   fChain->SetBranchAddress("superCluster_e3x3Ratio", &superCluster_e3x3Ratio, &b_superCluster_e3x3Ratio);
   fChain->SetBranchAddress("superCluster_eMaxRatio", &superCluster_eMaxRatio, &b_superCluster_eMaxRatio);
   fChain->SetBranchAddress("superCluster_e2ndRatio", &superCluster_e2ndRatio, &b_superCluster_e2ndRatio);
   fChain->SetBranchAddress("superCluster_eTopRatio", &superCluster_eTopRatio, &b_superCluster_eTopRatio);
   fChain->SetBranchAddress("superCluster_eRightRatio", &superCluster_eRightRatio, &b_superCluster_eRightRatio);
   fChain->SetBranchAddress("superCluster_eBottomRatio", &superCluster_eBottomRatio, &b_superCluster_eBottomRatio);
   fChain->SetBranchAddress("superCluster_eLeftRatio", &superCluster_eLeftRatio, &b_superCluster_eLeftRatio);
   fChain->SetBranchAddress("superCluster_e2x5MaxRatio", &superCluster_e2x5MaxRatio, &b_superCluster_e2x5MaxRatio);
   fChain->SetBranchAddress("superCluster_e2x5TopRatio", &superCluster_e2x5TopRatio, &b_superCluster_e2x5TopRatio);
   fChain->SetBranchAddress("superCluster_e2x5RightRatio", &superCluster_e2x5RightRatio, &b_superCluster_e2x5RightRatio);
   fChain->SetBranchAddress("superCluster_e2x5BottomRatio", &superCluster_e2x5BottomRatio, &b_superCluster_e2x5BottomRatio);
   fChain->SetBranchAddress("superCluster_e2x5LeftRatio", &superCluster_e2x5LeftRatio, &b_superCluster_e2x5LeftRatio);
   fChain->SetBranchAddress("superCluster_swissCross", &superCluster_swissCross, &b_superCluster_swissCross);
   fChain->SetBranchAddress("superCluster_r9", &superCluster_r9, &b_superCluster_r9);
   fChain->SetBranchAddress("superCluster_sigmaIetaIeta", &superCluster_sigmaIetaIeta, &b_superCluster_sigmaIetaIeta);
   fChain->SetBranchAddress("superCluster_sigmaIetaIphi", &superCluster_sigmaIetaIphi, &b_superCluster_sigmaIetaIphi);
   fChain->SetBranchAddress("superCluster_sigmaIphiIphi", &superCluster_sigmaIphiIphi, &b_superCluster_sigmaIphiIphi);
   fChain->SetBranchAddress("superCluster_full5x5_e5x5", &superCluster_full5x5_e5x5, &b_superCluster_full5x5_e5x5);
   fChain->SetBranchAddress("superCluster_full5x5_e2x2Ratio", &superCluster_full5x5_e2x2Ratio, &b_superCluster_full5x5_e2x2Ratio);
   fChain->SetBranchAddress("superCluster_full5x5_e3x3Ratio", &superCluster_full5x5_e3x3Ratio, &b_superCluster_full5x5_e3x3Ratio);
   fChain->SetBranchAddress("superCluster_full5x5_eMaxRatio", &superCluster_full5x5_eMaxRatio, &b_superCluster_full5x5_eMaxRatio);
   fChain->SetBranchAddress("superCluster_full5x5_e2ndRatio", &superCluster_full5x5_e2ndRatio, &b_superCluster_full5x5_e2ndRatio);
   fChain->SetBranchAddress("superCluster_full5x5_eTopRatio", &superCluster_full5x5_eTopRatio, &b_superCluster_full5x5_eTopRatio);
   fChain->SetBranchAddress("superCluster_full5x5_eRightRatio", &superCluster_full5x5_eRightRatio, &b_superCluster_full5x5_eRightRatio);
   fChain->SetBranchAddress("superCluster_full5x5_eBottomRatio", &superCluster_full5x5_eBottomRatio, &b_superCluster_full5x5_eBottomRatio);
   fChain->SetBranchAddress("superCluster_full5x5_eLeftRatio", &superCluster_full5x5_eLeftRatio, &b_superCluster_full5x5_eLeftRatio);
   fChain->SetBranchAddress("superCluster_full5x5_e2x5MaxRatio", &superCluster_full5x5_e2x5MaxRatio, &b_superCluster_full5x5_e2x5MaxRatio);
   fChain->SetBranchAddress("superCluster_full5x5_e2x5TopRatio", &superCluster_full5x5_e2x5TopRatio, &b_superCluster_full5x5_e2x5TopRatio);
   fChain->SetBranchAddress("superCluster_full5x5_e2x5RightRatio", &superCluster_full5x5_e2x5RightRatio, &b_superCluster_full5x5_e2x5RightRatio);
   fChain->SetBranchAddress("superCluster_full5x5_e2x5BottomRatio", &superCluster_full5x5_e2x5BottomRatio, &b_superCluster_full5x5_e2x5BottomRatio);
   fChain->SetBranchAddress("superCluster_full5x5_e2x5LeftRatio", &superCluster_full5x5_e2x5LeftRatio, &b_superCluster_full5x5_e2x5LeftRatio);
   fChain->SetBranchAddress("superCluster_full5x5_swissCross", &superCluster_full5x5_swissCross, &b_superCluster_full5x5_swissCross);
   fChain->SetBranchAddress("superCluster_full5x5_r9", &superCluster_full5x5_r9, &b_superCluster_full5x5_r9);
   fChain->SetBranchAddress("superCluster_full5x5_sigmaIetaIeta", &superCluster_full5x5_sigmaIetaIeta, &b_superCluster_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("superCluster_full5x5_sigmaIetaIphi", &superCluster_full5x5_sigmaIetaIphi, &b_superCluster_full5x5_sigmaIetaIphi);
   fChain->SetBranchAddress("superCluster_full5x5_sigmaIphiIphi", &superCluster_full5x5_sigmaIphiIphi, &b_superCluster_full5x5_sigmaIphiIphi);
   Notify();
}

Bool_t BMMGNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BMMGNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BMMGNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BMMGNtuple_cxx
