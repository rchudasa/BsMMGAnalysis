#include "Util.h"

Int_t doMuonSelection(MergedBMMX &ntupleRawTree,Int_t muIdx, bool isLead)
{
    
   // auto muGblId= ntupleRawTree.b5_mm_mu1_index[muIdx];
    // Soft Muon ID   : Cut Based
    
    /*

    if( not ntupleRawTree.Muon_softId[muGblId] ) return 1;
    
    */

    // Soft Muon MVA  : Cut Based
    
    if( not ntupleRawTree.b5_Muon_softMvaId[muIdx] ) return 2;
    
    return 0;

    // BMM5 MVA

    const Double_t BDTLooseWorkingPoint (0.0); 
    const Double_t BDTTightWorkingPoint (0.0);

    if(ntupleRawTree.b5_MuonId_newSoftMuonMva[muIdx] < BDTTightWorkingPoint ) return 3;

    return 0;

}
