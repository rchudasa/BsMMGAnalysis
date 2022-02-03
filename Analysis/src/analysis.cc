/*
 *          USAGE
 *
 * root -b -q 'bmmX_distributionStudy.cc("YOURCONFIG.CFG")'
 *
 * */

#include "main.h"
#include "BMMGAnalysis.h" 

int  main(int argc,char *argv[])
{

    if(argc<2)
    {
        std::cout<<" Usage : \n"
                 <<"         ./main.exe <configfile.cfg> <option>\n"
                 <<" Option 0 : GenAnalysis \n"
                 <<" Option 1 : DataAnalysis \n"
                 <<endl;
        exit(1);
    }
    int val(1);
    
    if(argc > 2)
    {
        val=atoi(argv[2]);        
    }
    std::cout<<"\n VAL = "<<val<<"\n";
    string cfgFile(argv[1]);
    Long64_t maxEvents(-1000);
    
    //   TTree * outTree   = new TTree("AnalysisTree","Reduced branches from the Merged [bmm5+bmmX] trees.");
    
    BMMGAnalysis analyzer2018;
    std::cout<<"Initializing with cfg : "<<cfgFile<<"\n";
    analyzer2018.Init(cfgFile);
    
    analyzer2018.SetupAnalysis();
    #ifdef __MCANALYSIS__
    if(val==0) {
        std::cout<<" Doing GEN Analysis "<<"\n";
        analyzer2018.GenAnalyze();
    }
    #endif
    if(val==1) {
        std::cout<<" Doing Data Analysis "<<"\n";
        analyzer2018.Analyze();
    }

    analyzer2018.SaveFile();
    
}
