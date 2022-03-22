#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "sstream"
#include "fstream"
#include "iostream"
using namespace std;



/*
 *          USAGE
 *
 * root -b -q 'bmmX_distributionStudy.cc("YOURCONFIG.CFG")'
 *
 * */

#include "TreeMaker.h" 

int  main(int argc,char *argv[])
{

    if(argc<2)
    {
        std::cout<<" Usage : \n"
                 <<"         ./main.exe <configfile.cfg> <DO_GEN>\n\n";
        exit(1);
    }
    int val(0);
    
    if(argc > 2)
    {
        val=atoi(argv[2]);        
    }

    std::cout<<"\n VAL = "<<val<<"\n";
    string cfgFile(argv[1]);
    
    TreeMaker aTreeeMaker;
    aTreeeMaker.Init(cfgFile);
    aTreeeMaker.SetupAnalysis(true);
    
    if(val == 0) aTreeeMaker.DataDimuonTreeMaker();
    //if(val == 1) aTreeeMaker.genMatchedDimuons();
    aTreeeMaker.SaveFile();
}
