#ifndef __GENERICTREEMAKER__
#define __GENERICTREEMAKER__

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "iostream"
#define MAX_ARRAY_SIZE 2000

#define BMMGNtuple_cxx
#define BMM5Ntuple_cxx
#include "BMM5Ntuple.h"
#include "BMMGNtuple.h"

class genericTreeBranchSelector   {

    public:
    
    BMM5Ntuple   b5;
    BMMGNtuple   bG;


   #include  "DataMembers.h"
  
    TTree * _compiledTree; 
    void FillTheArraysFromVectors();
    void setCompiledTreeBranches();
    void Fill();
    void Write();

    
    
    genericTreeBranchSelector( ) : _compiledTree(nullptr){};
    
    genericTreeBranchSelector( TChain * motherTree, TChain* searchTreeChain):
        bG(motherTree),
        b5(searchTreeChain),
        _compiledTree(nullptr)
    {
        std::cout<<"Entries in BmmgTree : "<<bG.fChain->GetEntries()<<"\n";
        std::cout<<"Entries in Bmm5Tree : "<<b5.fChain->GetEntries()<<"\n";
        _compiledTree = new TTree("mergedTree","bmm5 and bmmg merged tree") ;
        if( not _compiledTree) exit(1);
    //    _compiledTree->SetAutoSave(30);
    //    _compiledTree->SetCircular(200000);
    //      HACK
	//    _compiledTree->Branch("bmmgTree_nMuMu",	&( bmmgTree_nMuMu ));
	//      Hack
        setCompiledTreeBranches();
    }



};


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


#include "MemberFunction.h"

#endif //  __GENERICTREEMAKER__


