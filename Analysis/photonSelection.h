#include "Util.h"

Int_t doPhotonSelection(MergedBMMX &ntupleRawTree , Int_t scIdx)
{
    
    if(ntupleRawTree.bG_scEta[scIdx] > 1.0 ) return 1;
    return 0;
}
