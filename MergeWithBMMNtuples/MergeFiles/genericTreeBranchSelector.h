#ifndef __GENERICTREEMAKER__
#define __GENERICTREEMAKER__

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "iostream"
#define MAX_ARRAY_SIZE 100

#include "Bmm5Ntuple.h"
#include "BMMGNtuple.h"

class genericTreeBranchSelector   {

    public:
    
    Bmm5Ntuple   b5;
    BMMGNtuple   bG;

//      DATAMEMBERS   

	Double_t  bG_primaryVertex_x[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_y[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_z[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_t[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_covXX[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_covXY[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_covXZ[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_covYY[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_covYZ[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_covZZ[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_x_error[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_y_error[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_z_error[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_t_error[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_ntracks[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_ndof[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_chi2[MAX_ARRAY_SIZE] ;
	Double_t  bG_primaryVertex_normalizedChi2[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoE[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoEt[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoEta[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPhi[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCE[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCEt[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCRawE[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoESEnP1[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoESEnP2[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCEta[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCPhi[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCEtaWidth[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCPhiWidth[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSCBrem[MAX_ARRAY_SIZE] ;
	Int_t  bG_phohasPixelSeed[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoR9[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoHoverE[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoESEffSigmaRR[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSigmaIEtaIEtaFull5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSigmaIEtaIPhiFull5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSigmaIPhiIPhiFull5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoE2x2Full5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoE5x5Full5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoR9Full5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPFChIso[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPFPhoIso[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPFNeuIso[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoEcalPFClusterIso[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoHcalPFClusterIso[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSeedTime[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoSeedEnergy[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoMIPTotEnergy[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoMIPChi2[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoMIPSlope[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoMIPIntercept[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoMIPNhitCone[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoMIPIsHalo[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPFE[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPFEt[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPFEta[MAX_ARRAY_SIZE] ;
	Float_t  bG_phoPFPhi[MAX_ARRAY_SIZE] ;
	Float_t  bG_scE[MAX_ARRAY_SIZE] ;
	Float_t  bG_scRawE[MAX_ARRAY_SIZE] ;
	Float_t  bG_scEta[MAX_ARRAY_SIZE] ;
	Float_t  bG_scPhi[MAX_ARRAY_SIZE] ;
	Float_t  bG_scX[MAX_ARRAY_SIZE] ;
	Float_t  bG_scY[MAX_ARRAY_SIZE] ;
	Float_t  bG_scZ[MAX_ARRAY_SIZE] ;
	Float_t  bG_scEtaWidth[MAX_ARRAY_SIZE] ;
	Float_t  bG_scPhiWidth[MAX_ARRAY_SIZE] ;
	Float_t  bG_scRawEt[MAX_ARRAY_SIZE] ;
	Float_t  bG_scMinDrWithGsfElectornSC_[MAX_ARRAY_SIZE] ;
	Bool_t  bG_scFoundGsfMatch_[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e2x2Ratio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e3x3Ratio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_eMaxRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e2ndRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_eTopRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_eRightRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_eBottomRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_eLeftRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e2x5MaxRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e2x5TopRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e2x5RightRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e2x5BottomRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_e2x5LeftRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_swissCross[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_r9[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_sigmaIetaIeta[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_sigmaIetaIphi[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_sigmaIphiIphi[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e5x5[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e2x2Ratio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e3x3Ratio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_eMaxRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e2ndRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_eTopRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_eRightRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_eBottomRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_eLeftRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e2x5MaxRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e2x5TopRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e2x5RightRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e2x5BottomRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_e2x5LeftRatio[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_swissCross[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_r9[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_sigmaIetaIeta[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_sigmaIetaIphi[MAX_ARRAY_SIZE] ;
	Float_t  bG_superCluster_full5x5_sigmaIphiIphi[MAX_ARRAY_SIZE] ;
    //      END DATA MEMBERS   
  
  
  //  Int_t bmmgTree_nMuMu;



    TTree * _compiledTree; 
    void FillTheArraysFromVectors();
    void setCompiledTreeBranches();
    void Fill();
    void Write();
    
    genericTreeBranchSelector( ){};
    
    genericTreeBranchSelector( TChain * motherTree, TChain* searchTreeChain):
        bG(motherTree),
        b5(searchTreeChain)
    {
        std::cout<<"Entries in BmmgTree : "<<bG.fChain->GetEntries()<<"\n";
        std::cout<<"Entries in Bmm5Tree : "<<b5.fChain->GetEntries()<<"\n";
        _compiledTree = new TTree("mergedTree","bmm5 and bmmg merged tree") ;
    //      HACK
	//    _compiledTree->Branch("bmmgTree_nMuMu",	&( bmmgTree_nMuMu ));
	//      Hack
        setCompiledTreeBranches();
    }



};

#endif //  __GENERICTREEMAKER__

#ifdef __GENERICTREEMAKERCXX__

void genericTreeBranchSelector::Fill()
{
    // Hack
    //bmmgTree_nMuMu = bG.mumuPt->size();
    // Hack
    
    FillTheArraysFromVectors();
    _compiledTree->Fill();
}

void genericTreeBranchSelector::Write()
{
    _compiledTree->Write();
}

#endif //  __GENERICTREEMAKERCXX__



