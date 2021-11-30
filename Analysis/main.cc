/*
 *          USAGE
 *
 * root -b -q 'bmmX_distributionStudy.cc("YOURCONFIG.CFG")'
 *
 * */

#include "main.h"
#include "BMMGAnalysis.h" 

int  main()
{
    string cfgFile="analysis2018.cfg";
    Long64_t maxEvents(-1000);
    
    //   TTree * outTree   = new TTree("AnalysisTree","Reduced branches from the Merged [bmm5+bmmX] trees.");
    
    BMMGAnalysis analyzer2018;
    analyzer2018.Init(cfgFile);
    
    analyzer2018.SetupAnalysis();
    analyzer2018.Analyze();
    analyzer2018.SaveFile();

    
}
