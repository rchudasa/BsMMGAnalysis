#include <iostream>
#include "Util.h"
#include "chrono"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"



#define PHO_MASS 0.0
// PDG 2012
#define  Mu_MASS 0.10565837

#define NSTORAGE_ARRAY_MAX 200000
#define NDIMU_MAX 20
#define NBMMG_MAX 40
#define NSC_MAX   50
#define NMUONS_MAX  20


#define MergedBMMX2018Data_cxx
#include "MergedBMMX2018Data.h"
typedef MergedBMMX2018Data MergedBMMX ;

struct PhotonMVAvariables
{

  float energy ;
  float et ;
  float eta ;
  float rawE ;
  float full5x5_e5x5  ;
  float full5x5_r9  ;
  float sigmaIetaIeta  ;
  float sigmaIetaIphi  ;
  float etaWidth  ;
  float phiWidth  ;
  
  float r9  ;
  float swisross  ;
  float PFPhoIso ;
  float PFNeuIso ;
  float PFChIso ;
  float FoundGsfMatch;
  float e2x5_MaxRatio; 

};

class BMMGAnalysis 
{
    public  :
    

    // Data members

    // file handling vars
    std::vector<string> InFileList;
    string ofileName;
    string treeName;
    string prefix;
    TFile * outputFile ;
    TChain *treeChain ;
    TTree *outTree;
    MergedBMMX ntupleRawTree;
    
    // Storage Vars
    // Defenition of the storage of type Double_t
    Double_t *storageArrayDouble;
    std::map<string, Int_t > candidateMapDouble;
    Int_t storageIdxFilledInt;
    // Defenition of the storage of type Int_t
    Int_t* storageArrayInt ;
    std::map<string, Int_t > candidateMapInt;
    Int_t storageIdxFilledDouble; 
    Int_t* photonSelectionCheck;
    
    // OutPut Tree vars
    Int_t nDiMuCandidates;
    Int_t nBMMGCandidates;
    Int_t nBMMGCandidatesPerDimu;
    bool isTriggerd;
    
    // book keeping vars
    bool initDone;
    std::map<string,string> string_parameters;
    std::map<string,Double_t> double_parameters;
    Long64_t nentries, maxEvents ;
 
    // Photon MVA vars 
    bool doPhotonMVA;
    bool hasWeightFiles;
    bool hasSetupPhotonMVA;
    TMVA::Reader *reader ;
    PhotonMVAvariables photonMVAdata;
    TString photonIdxMVAWeightFile;
    Float_t photonMVAValue;



    // Analysis Cuts
    Double_t maxMuMuDr          ;
    Double_t maxDimuPhotonDr    ;
    Double_t maxDimuMass        ;
    Double_t minDimuMass        ;
    Double_t maxMMGMass         ;
    Double_t minMMGMass         ;
    // Member functions 
    
    // Constructor
    BMMGAnalysis();
    
    void Init( string cfgFileName);
    // Helper funtions
    void readParameters(string fname);
    void setupInputTree();
    void AllocateMemory();
    void AllocateBMMGBranches();
    void setupBranchStatus();
    void SaveFile();
    void setupOutPuts();
    void SetupAnalysis();
    
    // Analysis Functions
    void doPhotonMVAScores();
    void setUpPhotonMVA();
    
    Int_t doMuonSelection(Int_t muIdx, bool isLead);
    Int_t doPhotonSelection(Int_t scIdx);
    void Analyze();
};

BMMGAnalysis::BMMGAnalysis()
{
    initDone=false;
}

void BMMGAnalysis::Init(string cfgFileName)
{
    treeName="mergedTree";
    prefix="workarea/";
    ofileName="output.root";
    InFileList.clear();
    ofileName="output.root";
    maxEvents=1000;
    
    maxMuMuDr          =1.4;
    maxDimuPhotonDr    =1.4;
    maxDimuMass        =6.0;
    minDimuMass        =1.0;
    maxMMGMass         =6.2;
    minMMGMass         =4.4;
    doPhotonMVA        =false;
    photonIdxMVAWeightFile="";
    hasWeightFiles=false;

    readParameters(cfgFileName);
    
    if(doPhotonMVA)
    {
       setUpPhotonMVA();
    }
        
    initDone=true;
}

void BMMGAnalysis::SetupAnalysis()
{

    if( not initDone) 
    {
        std::cout<<"\n Init ur vars before setupAnalysis !! \n";
        return;
    }
    
    setupInputTree();
    setupOutPuts();

}

void BMMGAnalysis::setupInputTree()
{   
    if( not initDone) 
    {
        std::cout<<"\n Init ur vars before setupInputTree !! \n";
        return;
    }

    treeChain = new TChain(treeName.c_str());
    for(auto i=0;i<InFileList.size();i++)
    {
        treeChain->Add(InFileList[i].c_str());
    }
    
    nentries = treeChain->GetEntries();
    if( maxEvents < 0) maxEvents = nentries;
    maxEvents = nentries < maxEvents ? nentries: maxEvents;
    cout<<"Available total number of events "<<nentries<<" \n";
    
    ntupleRawTree.Init(treeChain);

}

void BMMGAnalysis::setupBranchStatus()
{
    //ntupleRawTree.fChain->SetBranchStatus("b5_Electron_",0);
    //ntupleRawTree.fChain->SetBranchStatus("*MET*",0);
    //ntupleRawTree.fChain->SetBranchStatus("*Jet*",0);
    //ntupleRawTree.fChain->SetBranchStatus("*Tau*",0);
    //ntupleRawTree.fChain->SetBranchStatus("*b5_L1_*",0);
    //ntupleRawTree.fChain->SetBranchStatus("*b5_HLT_*",0);

    //ntupleRawTree.fChain->SetBranchStatus("b5_HLT_DoubleMu4_3_Bs",1);
}

void BMMGAnalysis::setupOutPuts()
{
   outputFile = new  TFile((prefix+ofileName).c_str(),"recreate");    
    
   setupBranchStatus();
   outTree   = treeChain->CloneTree(0);
   AllocateMemory();
   AllocateBMMGBranches();
}

void BMMGAnalysis::AllocateMemory()
{

    // Defenition of the storage of type Double_t
    
    storageArrayDouble=new Double_t[NSTORAGE_ARRAY_MAX];
    std::cout<<"Allocated "<<sizeof(Double_t)*NSTORAGE_ARRAY_MAX/1024<<" kB of storage for "<<NSTORAGE_ARRAY_MAX<<" Doubles \n";
    storageIdxFilledDouble=0; 
    
    storageArrayInt =  new Int_t[NSTORAGE_ARRAY_MAX]; 
    std::cout<<"Allocated "<<sizeof(Int_t)*NSTORAGE_ARRAY_MAX/1024<<" kB of storage for "<<NSTORAGE_ARRAY_MAX<<" Ints\n";
    storageIdxFilledInt=0;

    photonSelectionCheck = new Int_t[NSC_MAX];  
    std::cout<<"Allocated "<<sizeof(Int_t)*NSC_MAX/1024<<" kB of storage for "<<NSTORAGE_ARRAY_MAX<<" Ints\n";
}

void BMMGAnalysis::SaveFile()
{

    outTree->Write();
    outputFile->Write();
    outputFile->Purge();
    outputFile->Close();

}

void BMMGAnalysis::readParameters(string fname)
{
    fstream cfgFile(fname,ios::in);
	string line;
	bool cfgModeFlag=false;

    cfgFile.clear();
    cfgFile.seekg(0,ios::beg);
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
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
            if(field.compare("MaxEvents")==0){
                 getline(strStream, field);
                 maxEvents=std::atoi(field.c_str());
                 cout<<" setting maxEvents  = "<<maxEvents<<"\n";
            }
            if(field.compare("MaxMuMuDr")==0){
                 getline(strStream, field);
                 maxMuMuDr=std::atof(field.c_str());
                 cout<<" setting maxMuMuDr  = "<<maxMuMuDr<<"\n";
            }
            if(field.compare("MaxDimuPhotonDr")==0){
                 getline(strStream, field);
                 maxDimuPhotonDr=std::atof(field.c_str());
                 cout<<" setting maxMuMuDr  = "<<maxDimuPhotonDr<<"\n";
            }
            if(field.compare("MinDimuMass")==0){
                 getline(strStream, field);
                 minDimuMass=std::atof(field.c_str());
                 cout<<" setting maxDimuMass  = "<<maxDimuMass<<"\n";
            }
            if(field.compare("MaxDimuMass")==0){
                 getline(strStream, field);
                 maxDimuMass=std::atof(field.c_str());
                 cout<<" setting maxDimuMass  = "<<maxDimuMass<<"\n";
            }
             if(field.compare("MaxMMGMass")==0){
                 getline(strStream, field);
                 maxMMGMass=std::atof(field.c_str());
                 cout<<" setting maxMMGMass  = "<<maxMMGMass<<"\n";
            }
             if(field.compare("MinMMGMass")==0){
                 getline(strStream, field);
                 minMMGMass=std::atof(field.c_str());
                 cout<<" setting minMMGMass  = "<<minMMGMass<<"\n";
            }
             if(field.compare("DoPhotonMVAID")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 doPhotonMVA= tmpI >0 ? 1 : 0;
                 cout<<" setting DoPhotonMVAID  = "<<doPhotonMVA<<"\n";
            }
             if(field.compare("PhotonIDWeightFile")==0){
                 hasWeightFiles=true;
                 getline(strStream, field);
                 photonIdxMVAWeightFile=field;
                 std::cout<<" setting photonIdxMVAWeightFile = "<<photonIdxMVAWeightFile<<"\n";
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
    std::cout<<"File List has the following files : \n";
    for( auto name : InFileList)
    {
        std::cout<<"\t"<<name<<"\n";
    }

    if(doPhotonMVA)
    {
        setUpPhotonMVA();
    }
}

#include "BMMGAnalysis.cc"

