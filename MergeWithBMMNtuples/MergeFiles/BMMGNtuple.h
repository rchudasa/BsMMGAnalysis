//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 27 17:00:50 2021 by ROOT version 6.14/09
// from TTree EventTree/Event data
// found on file: Charmonium_Run2018A12Nov2019_UL2018_rsbv1_1.root
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
   vector<float>   *scEt;
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
   vector<float>   *scE5x5;
   vector<float>   *scE2x2Ratio;
   vector<float>   *scE3x3Ratio;
   vector<float>   *scEMaxRatio;
   vector<float>   *scE2ndRatio;
   vector<float>   *scETopRatio;
   vector<float>   *scERightRatio;
   vector<float>   *scEBottomRatio;
   vector<float>   *scELeftRatio;
   vector<float>   *scE2x5MaxRatio;
   vector<float>   *scE2x5TopRatio;
   vector<float>   *scE2x5RightRatio;
   vector<float>   *scE2x5BottomRatio;
   vector<float>   *scE2x5LeftRatio;
   vector<float>   *scSwissCross;
   vector<float>   *scR9;
   vector<float>   *scSigmaIetaIeta;
   vector<float>   *scSigmaIetaIphi;
   vector<float>   *scSigmaIphiIphi;
   vector<float>   *scFull5x5_e5x5;
   vector<float>   *scFull5x5_e2x2Ratio;
   vector<float>   *scFull5x5_e3x3Ratio;
   vector<float>   *scFull5x5_eMaxRatio;
   vector<float>   *scFull5x5_e2ndRatio;
   vector<float>   *scFull5x5_eTopRatio;
   vector<float>   *scFull5x5_eRightRatio;
   vector<float>   *scFull5x5_eBottomRatio;
   vector<float>   *scFull5x5_eLeftRatio;
   vector<float>   *scFull5x5_e2x5MaxRatio;
   vector<float>   *scFull5x5_e2x5TopRatio;
   vector<float>   *scFull5x5_e2x5RightRatio;
   vector<float>   *scFull5x5_e2x5BottomRatio;
   vector<float>   *scFull5x5_e2x5LeftRatio;
   vector<float>   *scFull5x5_swissCross;
   vector<float>   *scFull5x5_r9;
   vector<float>   *scFull5x5_sigmaIetaIeta;
   vector<float>   *scFull5x5_sigmaIetaIphi;
   vector<float>   *scFull5x5_sigmaIphiIphi;
   Int_t           nhcalRechit;
   vector<float>   *hcalRechitIEta;
   vector<float>   *hcalRechitIPhi;
   vector<float>   *hcalRechitEnergy;
   vector<float>   *scPFChIso1;
   vector<float>   *scPFChIso2;
   vector<float>   *scPFChIso3;
   vector<float>   *scPFChIso4;
   vector<float>   *scPFChIso5;
   vector<float>   *scPFPhoIso1;
   vector<float>   *scPFPhoIso2;
   vector<float>   *scPFPhoIso3;
   vector<float>   *scPFPhoIso4;
   vector<float>   *scPFPhoIso5;
   vector<float>   *scPFNeuIso1;
   vector<float>   *scPFNeuIso2;
   vector<float>   *scPFNeuIso3;
   vector<float>   *scPFNeuIso4;
   vector<float>   *scPFNeuIso5;

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
   TBranch        *b_scEt;   //!
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
   TBranch        *b_scE5x5;   //!
   TBranch        *b_scE2x2Ratio;   //!
   TBranch        *b_scE3x3Ratio;   //!
   TBranch        *b_scEMaxRatio;   //!
   TBranch        *b_scE2ndRatio;   //!
   TBranch        *b_scETopRatio;   //!
   TBranch        *b_scERightRatio;   //!
   TBranch        *b_scEBottomRatio;   //!
   TBranch        *b_scELeftRatio;   //!
   TBranch        *b_scE2x5MaxRatio;   //!
   TBranch        *b_scE2x5TopRatio;   //!
   TBranch        *b_scE2x5RightRatio;   //!
   TBranch        *b_scE2x5BottomRatio;   //!
   TBranch        *b_scE2x5LeftRatio;   //!
   TBranch        *b_scSwissCross;   //!
   TBranch        *b_scR9;   //!
   TBranch        *b_scSigmaIetaIeta;   //!
   TBranch        *b_scSigmaIetaIphi;   //!
   TBranch        *b_scSigmaIphiIphi;   //!
   TBranch        *b_scFull5x5_e5x5;   //!
   TBranch        *b_scFull5x5_e2x2Ratio;   //!
   TBranch        *b_scFull5x5_e3x3Ratio;   //!
   TBranch        *b_scFull5x5_eMaxRatio;   //!
   TBranch        *b_scFull5x5_e2ndRatio;   //!
   TBranch        *b_scFull5x5_eTopRatio;   //!
   TBranch        *b_scFull5x5_eRightRatio;   //!
   TBranch        *b_scFull5x5_eBottomRatio;   //!
   TBranch        *b_scFull5x5_eLeftRatio;   //!
   TBranch        *b_scFull5x5_e2x5MaxRatio;   //!
   TBranch        *b_scFull5x5_e2x5TopRatio;   //!
   TBranch        *b_scFull5x5_e2x5RightRatio;   //!
   TBranch        *b_scFull5x5_e2x5BottomRatio;   //!
   TBranch        *b_scFull5x5_e2x5LeftRatio;   //!
   TBranch        *b_scFull5x5_swissCross;   //!
   TBranch        *b_scFull5x5_r9;   //!
   TBranch        *b_scFull5x5_sigmaIetaIeta;   //!
   TBranch        *b_scFull5x5_sigmaIetaIphi;   //!
   TBranch        *b_scFull5x5_sigmaIphiIphi;   //!
   TBranch        *b_nhcalRechit;   //!
   TBranch        *b_hcalRechitIEta;   //!
   TBranch        *b_hcalRechitIPhi;   //!
   TBranch        *b_hcalRechitEnergy;   //!
   TBranch        *b_scPFChIso1;   //!
   TBranch        *b_scPFChIso2;   //!
   TBranch        *b_scPFChIso3;   //!
   TBranch        *b_scPFChIso4;   //!
   TBranch        *b_scPFChIso5;   //!
   TBranch        *b_scPFPhoIso1;   //!
   TBranch        *b_scPFPhoIso2;   //!
   TBranch        *b_scPFPhoIso3;   //!
   TBranch        *b_scPFPhoIso4;   //!
   TBranch        *b_scPFPhoIso5;   //!
   TBranch        *b_scPFNeuIso1;   //!
   TBranch        *b_scPFNeuIso2;   //!
   TBranch        *b_scPFNeuIso3;   //!
   TBranch        *b_scPFNeuIso4;   //!
   TBranch        *b_scPFNeuIso5;   //!

   BMMGNtuple(TTree *tree=0);
   virtual ~BMMGNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
 //  virtual void     Loop();
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Charmonium_Run2018A12Nov2019_UL2018_rsbv1_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Charmonium_Run2018A12Nov2019_UL2018_rsbv1_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Charmonium_Run2018A12Nov2019_UL2018_rsbv1_1.root:/Ntuples");
      dir->GetObject("EventTree",tree);

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
   scEt = 0;
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
   scE5x5 = 0;
   scE2x2Ratio = 0;
   scE3x3Ratio = 0;
   scEMaxRatio = 0;
   scE2ndRatio = 0;
   scETopRatio = 0;
   scERightRatio = 0;
   scEBottomRatio = 0;
   scELeftRatio = 0;
   scE2x5MaxRatio = 0;
   scE2x5TopRatio = 0;
   scE2x5RightRatio = 0;
   scE2x5BottomRatio = 0;
   scE2x5LeftRatio = 0;
   scSwissCross = 0;
   scR9 = 0;
   scSigmaIetaIeta = 0;
   scSigmaIetaIphi = 0;
   scSigmaIphiIphi = 0;
   scFull5x5_e5x5 = 0;
   scFull5x5_e2x2Ratio = 0;
   scFull5x5_e3x3Ratio = 0;
   scFull5x5_eMaxRatio = 0;
   scFull5x5_e2ndRatio = 0;
   scFull5x5_eTopRatio = 0;
   scFull5x5_eRightRatio = 0;
   scFull5x5_eBottomRatio = 0;
   scFull5x5_eLeftRatio = 0;
   scFull5x5_e2x5MaxRatio = 0;
   scFull5x5_e2x5TopRatio = 0;
   scFull5x5_e2x5RightRatio = 0;
   scFull5x5_e2x5BottomRatio = 0;
   scFull5x5_e2x5LeftRatio = 0;
   scFull5x5_swissCross = 0;
   scFull5x5_r9 = 0;
   scFull5x5_sigmaIetaIeta = 0;
   scFull5x5_sigmaIetaIphi = 0;
   scFull5x5_sigmaIphiIphi = 0;
   hcalRechitIEta = 0;
   hcalRechitIPhi = 0;
   hcalRechitEnergy = 0;
   scPFChIso1 = 0;
   scPFChIso2 = 0;
   scPFChIso3 = 0;
   scPFChIso4 = 0;
   scPFChIso5 = 0;
   scPFPhoIso1 = 0;
   scPFPhoIso2 = 0;
   scPFPhoIso3 = 0;
   scPFPhoIso4 = 0;
   scPFPhoIso5 = 0;
   scPFNeuIso1 = 0;
   scPFNeuIso2 = 0;
   scPFNeuIso3 = 0;
   scPFNeuIso4 = 0;
   scPFNeuIso5 = 0;
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
   fChain->SetBranchAddress("scEt", &scEt, &b_scEt);
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
   fChain->SetBranchAddress("scE5x5", &scE5x5, &b_scE5x5);
   fChain->SetBranchAddress("scE2x2Ratio", &scE2x2Ratio, &b_scE2x2Ratio);
   fChain->SetBranchAddress("scE3x3Ratio", &scE3x3Ratio, &b_scE3x3Ratio);
   fChain->SetBranchAddress("scEMaxRatio", &scEMaxRatio, &b_scEMaxRatio);
   fChain->SetBranchAddress("scE2ndRatio", &scE2ndRatio, &b_scE2ndRatio);
   fChain->SetBranchAddress("scETopRatio", &scETopRatio, &b_scETopRatio);
   fChain->SetBranchAddress("scERightRatio", &scERightRatio, &b_scERightRatio);
   fChain->SetBranchAddress("scEBottomRatio", &scEBottomRatio, &b_scEBottomRatio);
   fChain->SetBranchAddress("scELeftRatio", &scELeftRatio, &b_scELeftRatio);
   fChain->SetBranchAddress("scE2x5MaxRatio", &scE2x5MaxRatio, &b_scE2x5MaxRatio);
   fChain->SetBranchAddress("scE2x5TopRatio", &scE2x5TopRatio, &b_scE2x5TopRatio);
   fChain->SetBranchAddress("scE2x5RightRatio", &scE2x5RightRatio, &b_scE2x5RightRatio);
   fChain->SetBranchAddress("scE2x5BottomRatio", &scE2x5BottomRatio, &b_scE2x5BottomRatio);
   fChain->SetBranchAddress("scE2x5LeftRatio", &scE2x5LeftRatio, &b_scE2x5LeftRatio);
   fChain->SetBranchAddress("scSwissCross", &scSwissCross, &b_scSwissCross);
   fChain->SetBranchAddress("scR9", &scR9, &b_scR9);
   fChain->SetBranchAddress("scSigmaIetaIeta", &scSigmaIetaIeta, &b_scSigmaIetaIeta);
   fChain->SetBranchAddress("scSigmaIetaIphi", &scSigmaIetaIphi, &b_scSigmaIetaIphi);
   fChain->SetBranchAddress("scSigmaIphiIphi", &scSigmaIphiIphi, &b_scSigmaIphiIphi);
   fChain->SetBranchAddress("scFull5x5_e5x5", &scFull5x5_e5x5, &b_scFull5x5_e5x5);
   fChain->SetBranchAddress("scFull5x5_e2x2Ratio", &scFull5x5_e2x2Ratio, &b_scFull5x5_e2x2Ratio);
   fChain->SetBranchAddress("scFull5x5_e3x3Ratio", &scFull5x5_e3x3Ratio, &b_scFull5x5_e3x3Ratio);
   fChain->SetBranchAddress("scFull5x5_eMaxRatio", &scFull5x5_eMaxRatio, &b_scFull5x5_eMaxRatio);
   fChain->SetBranchAddress("scFull5x5_e2ndRatio", &scFull5x5_e2ndRatio, &b_scFull5x5_e2ndRatio);
   fChain->SetBranchAddress("scFull5x5_eTopRatio", &scFull5x5_eTopRatio, &b_scFull5x5_eTopRatio);
   fChain->SetBranchAddress("scFull5x5_eRightRatio", &scFull5x5_eRightRatio, &b_scFull5x5_eRightRatio);
   fChain->SetBranchAddress("scFull5x5_eBottomRatio", &scFull5x5_eBottomRatio, &b_scFull5x5_eBottomRatio);
   fChain->SetBranchAddress("scFull5x5_eLeftRatio", &scFull5x5_eLeftRatio, &b_scFull5x5_eLeftRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5MaxRatio", &scFull5x5_e2x5MaxRatio, &b_scFull5x5_e2x5MaxRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5TopRatio", &scFull5x5_e2x5TopRatio, &b_scFull5x5_e2x5TopRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5RightRatio", &scFull5x5_e2x5RightRatio, &b_scFull5x5_e2x5RightRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5BottomRatio", &scFull5x5_e2x5BottomRatio, &b_scFull5x5_e2x5BottomRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5LeftRatio", &scFull5x5_e2x5LeftRatio, &b_scFull5x5_e2x5LeftRatio);
   fChain->SetBranchAddress("scFull5x5_swissCross", &scFull5x5_swissCross, &b_scFull5x5_swissCross);
   fChain->SetBranchAddress("scFull5x5_r9", &scFull5x5_r9, &b_scFull5x5_r9);
   fChain->SetBranchAddress("scFull5x5_sigmaIetaIeta", &scFull5x5_sigmaIetaIeta, &b_scFull5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("scFull5x5_sigmaIetaIphi", &scFull5x5_sigmaIetaIphi, &b_scFull5x5_sigmaIetaIphi);
   fChain->SetBranchAddress("scFull5x5_sigmaIphiIphi", &scFull5x5_sigmaIphiIphi, &b_scFull5x5_sigmaIphiIphi);
   fChain->SetBranchAddress("nhcalRechit", &nhcalRechit, &b_nhcalRechit);
   fChain->SetBranchAddress("hcalRechitIEta", &hcalRechitIEta, &b_hcalRechitIEta);
   fChain->SetBranchAddress("hcalRechitIPhi", &hcalRechitIPhi, &b_hcalRechitIPhi);
   fChain->SetBranchAddress("hcalRechitEnergy", &hcalRechitEnergy, &b_hcalRechitEnergy);
   fChain->SetBranchAddress("scPFChIso1", &scPFChIso1, &b_scPFChIso1);
   fChain->SetBranchAddress("scPFChIso2", &scPFChIso2, &b_scPFChIso2);
   fChain->SetBranchAddress("scPFChIso3", &scPFChIso3, &b_scPFChIso3);
   fChain->SetBranchAddress("scPFChIso4", &scPFChIso4, &b_scPFChIso4);
   fChain->SetBranchAddress("scPFChIso5", &scPFChIso5, &b_scPFChIso5);
   fChain->SetBranchAddress("scPFPhoIso1", &scPFPhoIso1, &b_scPFPhoIso1);
   fChain->SetBranchAddress("scPFPhoIso2", &scPFPhoIso2, &b_scPFPhoIso2);
   fChain->SetBranchAddress("scPFPhoIso3", &scPFPhoIso3, &b_scPFPhoIso3);
   fChain->SetBranchAddress("scPFPhoIso4", &scPFPhoIso4, &b_scPFPhoIso4);
   fChain->SetBranchAddress("scPFPhoIso5", &scPFPhoIso5, &b_scPFPhoIso5);
   fChain->SetBranchAddress("scPFNeuIso1", &scPFNeuIso1, &b_scPFNeuIso1);
   fChain->SetBranchAddress("scPFNeuIso2", &scPFNeuIso2, &b_scPFNeuIso2);
   fChain->SetBranchAddress("scPFNeuIso3", &scPFNeuIso3, &b_scPFNeuIso3);
   fChain->SetBranchAddress("scPFNeuIso4", &scPFNeuIso4, &b_scPFNeuIso4);
   fChain->SetBranchAddress("scPFNeuIso5", &scPFNeuIso5, &b_scPFNeuIso5);
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
