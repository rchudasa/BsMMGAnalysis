
#ifndef __MCANALYSIS__
  #include "MergedBMMX2018Data.h"
typedef MergedBMMX2018Data MergedBMMX ;
#else
  # include "MergedBMMXMC.h"  
typedef MergedBMMXMC MergedBMMX ;
#endif
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/Vector3D.h"

#include "map"
#include "iostream"
#include "Util.h"
#include "chrono"
#define NSTORAGE_ARRAY_MAX 5000

class TreeMaker 
{
    public  :
    
    // Data members
    std::vector<string> InFileList;
    string ofileName;
    string treeName;
    string prefix;
    TFile * outputFile ;
    TChain *treeChain ;
    MergedBMMX ntupleRawTree;
    Double_t drGenMatchMin;
    Float_t maxMuonEta     ;
    Double_t genParticlePtMin;
    Double_t eventGenMultiplicity;
    Bool_t genParticleIsStable;   
    Float_t minMuIsolation;
    
    // Storage Vars
    // Defenition of the storage of type Double_t
    Int_t storageIdxFilledDouble; 
    std::map<string, Float_t*> storageArrayDoubleMap;
    Double_t *storageArrayDouble;
    std::map<string, Int_t > candidateMapDouble;
    
    // Defenition of the storage of type Int_t
    Int_t storageIdxFilledInt;
    Int_t* storageArrayInt ;
    std::map<TString, Int_t > candidateMapInt;
    std::map<TString, Float_t > storageFloat;
    TTree* outTree;
    
    Int_t* photonSelectionCheck;
    
    // Histograms to Store
    std::map<TString,TH1F*> th1fStore;
    
    // Trees to store
    std::map<TString,TTree*> treeStore;

    // isMC
    bool isMC;
    bool doGenMatching;
    
    // book keeping vars
    bool initDone;
    std::vector<TTree*> treesToStore;
    Long64_t nentries, maxEvents ;
    Int_t reportEvery;
                    // Member functions 
    // Constructor
    TreeMaker();
    
    void Init( string cfgFileName);

    // Helper funtions
    void readParameters(string fname);
    void setupInputTree();
    void AllocateMemory();
    
    void SaveFile();
    void setupOutPuts(bool makeTreeCopy=false);
    void SetupAnalysis(bool makeTreeCopy=false);
    
    void AddSCHistos(TString tag="");
    void AddSCTree(TString SCTreeName="SCTreeStorage");
    
    // Histogram Related Functions
    void bookHistograms();
    void fill_scHists(Int_t idx =-1);
    void fillSCVariablesToOutTree(Int_t scIDX,TString SCTreeName);
    
    void fill_genHists(Int_t idx);
    void fill_scHists(Int_t scIDX,TString tag="",Double_t dr=0.0);
    void fill_eventHists();
    
    void DataDimuonTreeMaker();
    void fillDimuon(Int_t,TString);
    
    void genMatchedDimuons();
    
    void FillDimuonTree();
    void AddDimuonTree(TString tag);
    void fillDimuonVertexVars(Int_t mumuIdx);
    void fillDimuonVertexHists(TString tag);
    void AddDimuonVertexHistos(TString tag);
    Int_t doMuonSelection(Int_t muIdx, bool isLead);
    
    Bool_t muonHasToBeHighPurity  ;
    Bool_t muonHasToBeGlobal      ;
    Bool_t muonHasToBeLoose       ;
    Bool_t muonHasToBeTracker     ;
    Float_t BDTWorkingPoint ; 
    Float_t drMaxForGenMatchMu;
    Float_t drMaxForGenMatchSC;
    Float_t minMuonPt      ;

  private :
     int   pTBins;
     float pTmin;
     float pTmax;

     int   etaBins;
     float etamin;
     float etamax;

     int   deltaRNBins;
     float deltaRMin;
     float deltaRMax;



};

TreeMaker::TreeMaker()
{    
}
void TreeMaker::Init(string cfgFileName)
{
    treeName="mergedTree";
    prefix="workarea/";
    ofileName="output.root";
    InFileList.clear();
    ofileName="output.root";
    maxEvents=1000;
    isMC=false;
    doGenMatching=false;
    drGenMatchMin=0.1;
    reportEvery=10000;
    maxMuonEta =  1e9;
    minMuonPt=0.0;
    BDTWorkingPoint=0.58;
    minMuIsolation=0.8;

    readParameters(cfgFileName);
    AllocateMemory();
    
    initDone=true;
}

void TreeMaker::SetupAnalysis(bool makeTreeCopy)
{

    if( not initDone) 
    {
        std::cout<<"\n Init ur vars before setupAnalysis !! \n";
        return;
    }
    
    setupInputTree();
    bookHistograms();
    setupOutPuts(makeTreeCopy);

}

void TreeMaker::setupInputTree()
{   
    if( not initDone) 
    {
        std::cout<<"\n Init ur vars before setupInputTree !! \n";
        return;
    }

    treeChain = new TChain(treeName.c_str());
    auto rslt=0;
    for(auto i=0;i<InFileList.size();i++)
    {
        rslt=treeChain->AddFile(InFileList[i].c_str(),-1);
        if(rslt!=1) exit(112);
    }
    
    nentries = treeChain->GetEntries();
    if( maxEvents < 0) maxEvents = nentries;
    maxEvents = nentries < maxEvents ? nentries: maxEvents;
    cout<<"Available total number of events "<<nentries<<" \n";
    
    //ntupleRawTree.Init(treeChain,isMC);
    ntupleRawTree.Init(treeChain);
    std::cout<<"OUT OF HERE\n";

}

void TreeMaker::setupOutPuts(bool makeTreeCopy)
{
   outputFile = new  TFile((prefix+ofileName).c_str(),"recreate");    
    
   
   if(makeTreeCopy)
   {
   }

}

void TreeMaker::AllocateMemory()
{

}

void TreeMaker::SaveFile()
{
    if(not outputFile)
    {
        std::cout<<"\n OUT PUT FILE NOT SET \n";
    }

    outputFile->cd();

    for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
    {
        
  //      std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &ahist = *(it->second); 
        ahist.Write();
    }
       
    for (std::map<TString,TTree *>::iterator it=treeStore.begin() ; it!=treeStore.end(); ++it)
    {
        
 //       std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &atree = *(it->second); 
        atree.Write();
    }
    outputFile->Write();
    outputFile->Purge();
    outputFile->Close();

}

void TreeMaker::readParameters(string fname)
{
    fstream cfgFile(fname,ios::in);
	string line;
	bool cfgModeFlag=false;

    Double_t aDouble;
    
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
           if(field.compare("InputTreeName")==0){
                 getline(strStream, field);
                 treeName=field;
                 cout<<" setting treeName  = "<<prefix<<"\n";
            }
            if(field.compare("ReportEvery")==0){
                 getline(strStream, field);
                 reportEvery=std::atoi(field.c_str());
                 cout<<" setting reportEvery  = "<<reportEvery<<"\n";
            }
             if(field.compare("MaxEvents")==0){
                 getline(strStream, field);
                 maxEvents=std::atoi(field.c_str());
                 cout<<" setting maxEvents  = "<<maxEvents<<"\n";
            }
            if(field.compare("DrGenMatchMin")==0){
                 getline(strStream, field);
                 drGenMatchMin=std::atof(field.c_str());
                 cout<<" setting drGenMatchMin  = "<<drGenMatchMin<<"\n";
            }
            if(field.compare("IsMC")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 isMC = tmpI >0 ? 1 : 0;
                 cout<<" setting isMC  = "<<isMC<<"\n";
            }
           
           getBoolFromTag ("MuonHasToBeTracker"		  ,strStream , field , muonHasToBeTracker);
           getBoolFromTag ("MuonHasToBeGlobal"		  ,strStream , field , muonHasToBeGlobal);
           getBoolFromTag ("MuonHasToBeLoose"		  ,strStream , field , muonHasToBeLoose);
           getBoolFromTag ("MuonHasToBeHighPurity"	  ,strStream , field , muonHasToBeHighPurity);
           getFloatFromTag("DrMaxForGenMatchMu"       ,strStream, field  ,drMaxForGenMatchMu);
           getFloatFromTag("DrMaxForGenMatchSC"       ,strStream, field  ,drMaxForGenMatchSC);
           getFloatFromTag("BDTWorkingPoint"          ,strStream, field  ,BDTWorkingPoint);
           getFloatFromTag("MaxMuonEta"                ,strStream, field  ,maxMuonEta);
           getFloatFromTag("MinMuonPt"                ,strStream, field  ,minMuonPt);
           getFloatFromTag("MinMuIsolation"           ,strStream, field  ,minMuIsolation);

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

}

#include "../src/TreeMaker.cc"


