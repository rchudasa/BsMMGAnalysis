#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"

class MVATrainer 
{
    public  :
    
    // Data members
    std::vector<string> InFileList;
    vector<string> treeNames;
    string ofileName;
    string prefix;
    TFile * outputFile;

    string testTrainConfig;
    
    std::map<std::string,int> Use;
    std::map<std::string,string> mvaMethordOptions;
 
    string modelFactoryName;

    std::vector<string> signalFnames;
    std::vector<string> signalTreeNames;
    std::vector<TFile*> signalFiles ;
    std::vector<TTree*> signalTrees ;
    std::vector<Double_t> signalWeight;
    string sigCuts;

    std::vector<string> bkgFnames;
    std::vector<string> bkgTreeNames;
    std::vector<TFile*> bkgFiles ;
    std::vector<TTree*> bkgTrees ;
    std::vector<Double_t> bkgWeight;
    string bkgCuts;

    string factoryOptions;
    string mvaMethords;
    std::vector<string> mvaModelsToTrain;
    std::vector<string> mvaTrainVars;
    std::vector<string> spectatorVars;

    Long64_t nentries, maxEvents ;
    Bool_t initDone ;

    TMVA::Factory *factory;
    TMVA::DataLoader *dataloader;
 // Constructor
    MVATrainer();
    ~MVATrainer();
    
    void Init( string cfgFileName);

    // Helper funtions
    void readParameters(string fname);
    
    void InitializeTMVAModel();
    void SetupTMVAOptions();
    void SaveOutputs();
    void trainAndTest();
    void getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose=false);
};

MVATrainer::MVATrainer()
{
  initDone=false;
  factory=nullptr;
  dataloader=nullptr;
  outputFile=nullptr;

 
}
MVATrainer::~MVATrainer()
{

   delete factory;
   delete dataloader;
   
}


void MVATrainer::Init(string cfgFileName)
{
    prefix      =  "workarea/";
    ofileName   =  "output.root";
    maxEvents   =  1000;
    mvaMethords =  "" ;

    sigCuts  =  "";
    bkgCuts     =  "";
    factoryOptions="auto";
    testTrainConfig="nTrain_Signal=0.5:nTrain_Background=0.5:SplitMode=Random" ;
    
    modelFactoryName="aModel";

    readParameters(cfgFileName);
    
    InitializeTMVAModel();
    initDone=true;
}

void MVATrainer::InitializeTMVAModel()
{
  TMVA::Tools::Instance();
  
  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  outputFile = TFile::Open( (prefix+ofileName).c_str(), "RECREATE" );
  if(factoryOptions=="auto")
    factoryOptions="ROC:V=True:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";   
  std::cout<<factoryOptions<<"\n";
  factory = new TMVA::Factory( modelFactoryName.c_str(), outputFile, factoryOptions.c_str());
  dataloader=new TMVA::DataLoader(modelFactoryName.c_str());
}

void MVATrainer::SaveOutputs()
{

   // Save the output
   outputFile->Close();

}

void MVATrainer::readParameters(string fname)
{
    fstream cfgFile(fname,ios::in);

    Double_t aDouble;
    
    cfgFile.clear();
    cfgFile.seekg(0,ios::beg);
	bool cfgModeFlag=false;
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
	string line;
    Int_t tmpI;
    
	while(std::getline(cfgFile,line))
	{
	   if(line=="#PARAMS_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#PARAMS_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
       strStream.clear();
       strStream.str(line);
       while (getline(strStream, field,'='))
       {
           std::cout<<field<<"\n";
           if(field.compare("OutputFile")==0){
                 getline(strStream, field);
                 ofileName=field;
                 std::cout<<" setting ofileName = "<<ofileName<<"\n";
            }
            if(field.compare("OutputPrefix")==0){
                 getline(strStream, field);
                 prefix=field;
                 cout<<" setting prefix = "<<prefix<<"\n";
            }
           if(field.compare("MVAMethords")==0){
                 getline(strStream, field);
                 mvaMethords=field;
                 cout<<" setting mvaMethords = "<<mvaMethords<<"\n";
            }
            if(field.compare("MaxEvents")==0){
                 getline(strStream, field);
                 maxEvents=std::atoi(field.c_str());
                 cout<<" setting maxEvents  = "<<maxEvents<<"\n";
            }
            if(field.compare("FactoryOptions")==0){
                 getline(strStream, field);
                 factoryOptions=field;
                 cout<<" setting factoryOptions = "<<factoryOptions<<"\n";
            }
            if(field.compare("BkgCuts")==0){
                 getline(strStream, field);
                 bkgCuts=field;
                 cout<<" setting bkgCuts = "<<bkgCuts<<"\n";
            }
            if(field.compare("TestTrainConfig")==0){
                 getline(strStream, field);
                 testTrainConfig=field;
                 cout<<" setting testTrainConfig = "<<testTrainConfig<<"\n";
            }
            if(field.compare("SigCuts")==0){
                 getline(strStream, field);
                 sigCuts=field;
                 cout<<" setting siggCuts = "<<sigCuts<<"\n";
            }
            if(field.compare("ModelFactoryName")==0){
                 getline(strStream, field);
                 modelFactoryName=field;
                 cout<<" setting modelFactoryName = "<<modelFactoryName<<"\n";
            }

       }
    }

	// getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	while(std::getline(cfgFile,line))
	{
	   if(line=="#FILELIST_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#FILELIST_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   InFileList.push_back(line);
	}

    getVetorFillledFromConfigFile(cfgFile, mvaModelsToTrain   , "#MODELS_BEG", "#MODELS_END", true);
    mvaMethords="";
    for(int i=0;i<mvaModelsToTrain.size();i++)
    {
        if(mvaMethords=="")
            {
                mvaMethords=mvaModelsToTrain[i];
                continue;
            }
        mvaMethords+= ","+mvaModelsToTrain[i];
    }
    std::cout<<"setting mvaMethords = "<<mvaMethords<<"\n";
    std::vector<string> vect;
    vect.clear();
    getVetorFillledFromConfigFile(cfgFile, vect   , "#MODELPARAMS_BEG", "#MODELPARAMS_END", true);
    if( mvaModelsToTrain.size() != vect.size()) {std::cout<<__LINE__<<" config error !! \n";exit(1) ;}
    for(int i=0;i<mvaModelsToTrain.size();i++)
    {
            mvaMethordOptions[mvaModelsToTrain[i]]=vect[i];
    }

    
    getVetorFillledFromConfigFile(cfgFile, signalFnames   , "#SIGFILELIST_BEG", "#SIGFILELIST_END", true);
    getVetorFillledFromConfigFile(cfgFile, signalTreeNames, "#SIGTREELIST_BEG", "#SIGTREELIST_END", true);
    getVetorFillledFromConfigFile(cfgFile, bkgFnames      , "#BKGFILELIST_BEG", "#BKGFILELIST_END", true);
    getVetorFillledFromConfigFile(cfgFile, bkgTreeNames   , "#BKGTREELIST_BEG", "#BKGTREELIST_END", true);
    
    vect.clear();
    getVetorFillledFromConfigFile(cfgFile, vect     , "#SIGWEIGHT_BEG", "#SIGWEIGHT_END", true);
    for(int i=0;i<vect.size();i++)
    {
            signalWeight.push_back(atof(vect[i].c_str()));
    }
    std::cout<<" setting  signal weights as  : ";
    for(int i=0;i<signalWeight.size();i++) { std::cout<<" "<<signalWeight[i]<<" " ;}
    std::cout<<"\n";

    vect.clear();
    getVetorFillledFromConfigFile(cfgFile, vect     , "#BKGWEIGHT_BEG", "#BKGWEIGHT_END", true);
    for(int i=0;i<vect.size();i++)
    {
            bkgWeight.push_back(atof(vect[i].c_str()));
    }
    std::cout<<" setting  bkg weights as  : ";
    for(int i=0;i<bkgWeight.size();i++) { std::cout<<" "<<bkgWeight[i]<<" " ;}
    std::cout<<"\n";


    getVetorFillledFromConfigFile(cfgFile, mvaTrainVars   , "#MVAVARLIST_BEG", "#MVAVARLIST_END", true);
    getVetorFillledFromConfigFile(cfgFile, spectatorVars  , "#SPECTATORLIST_BEG", "#SPECTATORLIST_END", true);
    
    if(signalFnames.size() != signalTreeNames.size() )  {  std::cout<<__LINE__<<" config error !! \n"; exit(1) ;}
    if(signalFnames.size() != signalWeight.size() )     {  std::cout<<__LINE__<<" config error !! \n"; exit(1) ;}
    if(bkgFnames.size() != bkgTreeNames.size() )        {  std::cout<<__LINE__<<" config error !! \n"; exit(1) ;}
    if(bkgFnames.size() != bkgWeight.size() )           {  std::cout<<__LINE__<<" config error !! \n"; exit(1) ;}
}

void MVATrainer::getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose)
{
	
    bool cfgModeFlag=false;
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
	string line;
    
    // getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	int nItems(0);
    while(std::getline(cfgFile,line))
	{
	   if(line==beginTag) {cfgModeFlag=true;continue;}
	   if(line==endTag) {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   vecList.push_back(line);
	   nItems++;
    }

    if(verbose)
    {
       std::cout<<" Added "<<nItems<<" between "<<beginTag<<" and "<< endTag<<"\n";
       for( auto &name : vecList)
       {
           std::cout<<"\t"<<name<<"\n";
       }
    }

}


#include "../src/MVATrainer.cc"

