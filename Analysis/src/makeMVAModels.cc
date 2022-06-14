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


#include "MVATrainer.h" 

int  main(int argc,char *argv[])
{

    if(argc<2)
    {
        std::cout<<" Usage : \n"
                 <<"         ./makeMVAModels.exe <configfile.cfg> \n";
        exit(1);
    }
    
    string cfgFile(argv[1]);
    
    MVATrainer aMVAModelTrainer;
    aMVAModelTrainer.Init(cfgFile);
    aMVAModelTrainer.SetupTMVAOptions();
    aMVAModelTrainer.trainAndTest();
    aMVAModelTrainer.SaveOutputs();
}
