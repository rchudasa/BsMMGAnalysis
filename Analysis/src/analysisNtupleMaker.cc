/*
 *                        USAGE
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
                 <<"         ./main.exe <configfile.cfg>\n\n";
        exit(1);
    }

    string cfgFile(argv[1]);
    Long64_t maxEvents(-1000);
    
    //   TTree * outTree   = new TTree("AnalysisTree","Reduced branches from the Merged [bmm5+bmmX] trees.");
    
    BMMGAnalysis analyzer2018;
    analyzer2018.Init(cfgFile);
    
    analyzer2018.SetupAnalysis(true);
    analyzer2018.SkimData();
    analyzer2018.SaveFile();

}
