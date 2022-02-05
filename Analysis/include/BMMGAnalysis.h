#include <iostream>
#include "Util.h"
#include "chrono"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/Vector3D.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "RunLumiSelector.h"


#define PHO_MASS 0.0
// PDG 2012
#define  Mu_MASS 0.10565837

#define NSTORAGE_ARRAY_MAX 200000
#define NDIMU_MAX 20
#define NBMMG_MAX 40
#define NSC_MAX   50
#define NMUONS_MAX  20



//#define __MCANALYSIS__

#ifndef __MCANALYSIS__
  #include "MergedBMMX2018Data.h"
typedef MergedBMMX2018Data MergedBMMX ;
#else
  # include "MergedBMMXMC.h"  
typedef MergedBMMXMC MergedBMMX ;
#endif

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
    
    std::map<string,Double_t> storageDouble;
    std::map<string,Float_t> storageFloat;
    std::vector<string> mvaTrainVars;
    std::vector<string> spectatorVars;
    
    // Histograms to Store
    
    std::map<TString,TH1F*> th1fStore;

    // OutPut Tree vars
    Int_t nDiMuCandidates;
    Int_t nBMMGCandidates;
    Int_t nBMMGCandidatesPerDimu;
    bool isTriggerd;
     
    // isMC
    bool isMC;
    bool doGenMatching;
    bool doReducedTree;
    bool doTriggerFiltering;
    // RunLumiMask
    bool doRunLumiMask;
    TString runLumiMaskFname;
    RunLumiSelector runLumiMask;

    // book keeping vars
    bool initDone;
    std::map<string,string> string_parameters;
    std::map<string,Double_t> double_parameters;
    std::vector<TTree*> treesToStore;
    Long64_t nentries, maxEvents ;
 
    // Photon MVA vars 
    bool doPhotonMVA;
    bool hasWeightFiles;
    bool hasSetupPhotonMVA;
    TMVA::Reader *reader ;
    PhotonMVAvariables photonMVAdata;
    TString photonIdxMVAWeightFile;
    Float_t photonMVAValue;
    Float_t photonIDcut;

    // Muon MVA CUT
    const Double_t BDTWorkingPoint ; 
    
    // Analysis Cuts
    Double_t maxMuMuDr          ;
    Double_t maxDimuPhotonDr    ;
    Double_t maxMMGMass         ;
    Double_t minMMGMass         ;
    std::vector<Double_t> minDimuMass;
    std::vector<Double_t> maxDimuMass;

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
    void SetupAnalysis(bool makeTreeCopy=false);
    void setupOutPuts(bool makeTreeCopy=false);
    Int_t getMuonMatch(Double_t muEta,Double_t muPhi);
    // Histogram Related Functions
    void bookHistograms();
    Double_t getDCAGammaToDimuVertex(Int_t mumuIdx,Int_t phoId);
    void fill_muonHists(Int_t idx=-1);
    void fill_scHists(Int_t idx =-1);
    void fill_photonHists(Int_t idx  = -1);
    void fill_dimuonPassHists(Int_t idx=-1);
    void fill_dimuonHists(Int_t mumuIdx=-1);
    void fill_dimuonEnvironmentHists(Int_t mumuIdx);
    void fill_bmmgHists(TLorentzVector &bmmgLV,Int_t mumuIdx, Int_t phoSCIdx);
    void fill_globalEventHists();
    void getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose);
    // Analysis Functions
    void doPhotonMVAScores();
    void setUpPhotonMVA();
    
    Int_t doMuonSelection(Int_t muIdx, bool isLead);
    Int_t doDimuonSelection(Int_t mumuIdx);
    Int_t doVertexSelection(Int_t mumuIdx);
    Int_t doBMMGSelection(Int_t mumuIdx , Int_t phoSCIdx);
    Int_t doPhotonSelection(Int_t scIdx);
    Int_t getPVMatch(Int_t mumuIdx);
    void setupReducedAnalysisTreeBranches();
    void Analyze();
    #ifdef __MCANALYSIS__
    void GenAnalyze();
    #endif

    // Temporary vars with class scope
    TLorentzVector bmmgLV,diMuLV,photonLV;
    ROOT::Math::XYZVector svDisplacementVecor, bmmg3Momentum ;
};

BMMGAnalysis::BMMGAnalysis():
BDTWorkingPoint(0.58)
{
    initDone=false;
}

void BMMGAnalysis::Init(string cfgFileName)
{
    treeName="mergedTree";
    prefix="workarea/";
    ofileName="output.root";
    runLumiMaskFname="";
    InFileList.clear();
    ofileName="output.root";
    maxEvents=1000;
    photonIDcut=0.0;

    outTree=nullptr;

    maxMuMuDr          =1.4;
    maxDimuPhotonDr    =1.4;
    maxMMGMass         =6.2;
    minMMGMass         =4.4;
    doPhotonMVA        =false;
    doTriggerFiltering=false;
    photonIdxMVAWeightFile="";
    hasWeightFiles=false;
    
    doRunLumiMask=false;
    doReducedTree=true;
    minDimuMass.clear(); minDimuMass.push_back(0.0);
    maxDimuMass.clear(); maxDimuMass.push_back(6.0);

    readParameters(cfgFileName);
    
    if(doPhotonMVA)
    {
       setUpPhotonMVA();
    }
    if(doRunLumiMask)
    {
        runLumiMask.loadRunLumiMask(runLumiMaskFname);
    }
    initDone=true;
}

void BMMGAnalysis::SetupAnalysis(bool makeTreeCopy)
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

void BMMGAnalysis::setupInputTree()
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

void BMMGAnalysis::setupOutPuts(bool makeTreeCopy)
{
   outputFile = new  TFile((prefix+ofileName).c_str(),"recreate");    
    
   setupBranchStatus();
   AllocateMemory();
   if(makeTreeCopy)
   {
    if(doReducedTree)
     {
     outTree   = new TTree("mergedTree","Reduced Merged tree");
     setupReducedAnalysisTreeBranches();
    }
    else {
    outTree   = treeChain->CloneTree(0);
    }
    
    AllocateBMMGBranches();
   }

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
    outputFile->cd();

    if(outTree)
        outTree->Write();
       
    for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
    {
        
        //std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &ahist = *(it->second); 
        ahist.Write();
    }


    outputFile->Write();
    outputFile->Purge();
    outputFile->Close();

}

void BMMGAnalysis::readParameters(string fname)
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
            if(field.compare("DoTriggerFiltering")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 doTriggerFiltering= tmpI >0 ? 1 : 0;
                 cout<<" setting doTriggerFiltering  = "<<doTriggerFiltering<<"\n";
            }
            if(field.compare("DoRunLumiMask")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 doRunLumiMask= tmpI >0 ? 1 : 0;
                 cout<<" setting doRunLumiMask  = "<<doRunLumiMask<<"\n";
            }
            if(field.compare("RunLumiMask")==0){
                 getline(strStream, field);
                 runLumiMaskFname=field;
                 cout<<" setting runLumiMaskFname = "<<runLumiMaskFname<<"\n";
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
            if(field.compare("PhotonIDcut")==0){
                 getline(strStream, field);
                 photonIDcut=std::atof(field.c_str());
                 cout<<" setting photonIDcut  = "<<photonIDcut<<"\n";
            }
            if(field.compare("MaxDimuPhotonDr")==0){
                 getline(strStream, field);
                 maxDimuPhotonDr=std::atof(field.c_str());
                 cout<<" setting maxDimuPhotonDr  = "<<maxDimuPhotonDr<<"\n";
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
            
            if(field.compare("MinDimuMass")==0){
                 minDimuMass.clear();
                 cout<<" setting minDimuMass  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    minDimuMass.push_back(aDouble);
                 }
                 cout<<" }\n";
            }            
            if(field.compare("MaxDimuMass")==0){
                 maxDimuMass.clear();
                 cout<<" setting maxDimuMass  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    maxDimuMass.push_back(aDouble);
                 }
                 cout<<" }\n";
            } 
       }
    }

    getVetorFillledFromConfigFile(cfgFile, mvaTrainVars   , "#MVAVARLIST_BEG", "#MVAVARLIST_END", true);
    getVetorFillledFromConfigFile(cfgFile, spectatorVars  , "#SPECTATORLIST_BEG", "#SPECTATORLIST_END", true);
    

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

void BMMGAnalysis::getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose)
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


#include "../src/BMMGAnalysis.cc"

