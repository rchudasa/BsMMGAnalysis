#include "Util.h"

/*

Int_t setupDimuonBranches( TTree* outTree ,std::map<string, Int_t > &candidateMap , double * storageArray , Int_t offset=0 , Int_t size =0)
{
    #include  "setupDimu.cc"
}

void assignDimuonBMMGCandidates( std::map<string,Int_t > candidateMap, MergedBMMX &ntupleRawTree, int candIdx,double * storageArray,Int_t offset)
{
   #include "assignDimu.cc"
}

Int_t setupMuonBranches( TTree* outTree ,std::map<string, Int_t > &candidateMap , double * storageArray, Int_t offset=0 ,Int_t size =0)
{

    #include  "setupMuon.cc"

}

void assignMuons( std::map<string, Int_t > candidateMap, MergedBMMX &ntupleRawTree,double * storageArray)
{
    Int_t offset=0;
    Int_t candIdx=0;
    for(offset=0;offset < ntupleRawTree.b5_nMuon;offset++,candIdx++)
    {
        #include "assignMuon.cc"
    }
}

Int_t setupPhotonSCBranches( TTree* outTree ,std::map<string, Int_t > &candidateMap , double * storageArray, Int_t offset=0 , Int_t size =0 )
{

    #include  "setupSc.cc"

}

void assignSC( std::map<string, Int_t > candidateMap, MergedBMMX &ntupleRawTree, double * storageArray)
{
    Int_t offset=0;
    Int_t candIdx=0;
    for(  ;offset < ntupleRawTree.bG_nSC;offset++,candIdx++)
    {
        #include "assignSc.cc"
    }
}
*/
