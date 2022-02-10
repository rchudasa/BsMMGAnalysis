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
#include "RunLumiLogger.h"


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
    
    // cutflow to store
    
    std::map<string,Int_t> cutFlowOffsets;
    std::map<Int_t,TString> cutFlowIdxNamesMap;

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

    bool doRunLumiLog;
    RunLumiLogger runLumiLogger;
    string runLumiJsonName;
    
    // book keeping vars
    bool initDone;
    std::map<string,string> string_parameters;
    std::map<string,Double_t> double_parameters;
    std::vector<TTree*> treesToStore;
    Long64_t nentries, maxEvents ;
    Int_t reportEvery;
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
    
    // GEN MATCH CUTs

    Float_t drMaxForGenMatchMu;
    Float_t drMaxForGenMatchSC;

    // Analysis Cuts
    Float_t svMinDocaTrack;
    Int_t  svMaxNTracksClose;
    Int_t  svMaxN1SigmaTracksClose;
    Int_t  svMaxN2SigmaTracksClose;
    Int_t  svMaxN3SigmaTracksClose;

    Float_t minMuonPt      ;
    Float_t maxMuonEta     ;
    Float_t minSCEt        ;
    Float_t maxSCEta       ;
    Float_t minLeadMuIsolation;
    Float_t minSubLeadMuIsolation;
    Double_t maxMuMuDr          ;
    Double_t maxDimuPhotonDr    ;
    Double_t maxMMGMass         ;
    Double_t minMMGMass         ;
    std::vector<Double_t> minDimuMass;
    std::vector<Double_t> maxDimuMass;

    Bool_t muonHasToBeHighPurity;
    Bool_t muonHasToBeGlobal;
    Bool_t muonHasToBeLoose;
    Bool_t muonHasToBeTracker;

    Float_t muMuVtxProbabilityMin;
    Float_t muMuMaxDOCA;
    Float_t muMuMaxNTracksClose1Sigma;
    Float_t svGammaMaxDCA;
    

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
    void SetupCutFlowOffsets();
    Int_t getMuonMatch(Double_t muEta,Double_t muPhi);
    // Histogram Related Functions
    void bookHistograms();
    Double_t getDCAGammaToDimuVertex(Int_t mumuIdx,Int_t phoId);
    void FillCutFlow(Int_t);
    void fill_muonHists(Int_t idx=-1);
    void fill_scHists(Int_t idx =-1);
    void fill_photonHists(Int_t idx  = -1);
    void fill_dimuonPassHists(Int_t idx=-1);
    void fill_dimuonHists(Int_t mumuIdx=-1);
    void fill_dimuonEnvironmentHists(Int_t mumuIdx);
    void fill_bmmgHists(TLorentzVector &bmmgLV,Int_t mumuIdx, Int_t phoSCIdx);
    void fill_globalEventHists();
    void getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose);
    void getFloatFromTag(string tag,std::istringstream &strStream, string inputStr ,Float_t &var);
    void getIntFromTag(string tag,std::istringstream &strStream, string inputStr ,Int_t &var);
    void getBoolFromTag(string tag,std::istringstream &strStream, string inputStr ,Bool_t &var);
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
    runLumiJsonName="runLumi.json";
    runLumiMaskFname="";
    InFileList.clear();
    ofileName="output.root";
    maxEvents=1000;
    reportEvery=1;
    photonIDcut=0.0;
    muMuVtxProbabilityMin=0.0;
    muMuMaxDOCA=1e3;
    muMuMaxNTracksClose1Sigma=100;
    svGammaMaxDCA = 1e9;
    drMaxForGenMatchSC = 1e9;
    drMaxForGenMatchMu = 1e9;

    minMuonPt          = -1.0;
    maxMuonEta         =  1e9;
    minSCEt            = -1.0;
    maxSCEta           =  1e3;
    minLeadMuIsolation =  0.0;
    minSubLeadMuIsolation = 0.0 ;
    svMinDocaTrack = -0.1;
    svMaxNTracksClose        = 1000 ;
    svMaxN1SigmaTracksClose  = 1000 ;
    svMaxN2SigmaTracksClose  = 1000 ;
    svMaxN3SigmaTracksClose  = 1000 ;
    
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
    doRunLumiLog=false;
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
    SetupCutFlowOffsets();
    initDone=true;
}

void BMMGAnalysis::SetupCutFlowOffsets()
{
    Int_t i(0),j(0),offset(50);
    
    // Need to 
    i=0;
    j=0;
    j=1;    cutFlowIdxNamesMap[cutFlowOffsets["evtSelection"] + j ] = "allEventSelection";

    cutFlowOffsets["genSelection"]=i ;
    j=1;    cutFlowIdxNamesMap[cutFlowOffsets["genSelection"] + j ] = "gen_dimuonNoMuP";
    j=2;    cutFlowIdxNamesMap[cutFlowOffsets["genSelection"] + j ] = "gen_dimuonNoMuM";
    j=3;    cutFlowIdxNamesMap[cutFlowOffsets["genSelection"] + j ] = "gen_noPho";

    // Should always come the end
    i+=offset;
    cutFlowOffsets["basicCuts"]=i ;
    j=0;
    j=1;    cutFlowIdxNamesMap[ cutFlowOffsets["basicCuts"] + j ] = "noDimuons";
    j=2;    cutFlowIdxNamesMap[ cutFlowOffsets["basicCuts"] + j ] = "noSCs";
    j=3;    cutFlowIdxNamesMap[ cutFlowOffsets["basicCuts"] + j ] = "pfMuCut";
    j=4;    cutFlowIdxNamesMap[ cutFlowOffsets["basicCuts"] + j ] = "mu1NotInMuonCollection";
    j=5;    cutFlowIdxNamesMap[ cutFlowOffsets["basicCuts"] + j ] = "mu2NotInMuonCollection";
    
    i+=offset;
    cutFlowOffsets["muonSelection"]=i ;
    j=0;
    j++;    cutFlowIdxNamesMap[cutFlowOffsets["muonSelection"] + j ] = "muon_pT";
    j++;    cutFlowIdxNamesMap[cutFlowOffsets["muonSelection"] + j ] = "muon_eta";
    j++;    cutFlowIdxNamesMap[cutFlowOffsets["muonSelection"] + j ] = "muon_isTraker";
    j++;    cutFlowIdxNamesMap[cutFlowOffsets["muonSelection"] + j ] = "muon_isGlobal";
    j++;    cutFlowIdxNamesMap[cutFlowOffsets["muonSelection"] + j ] = "muon_isLooseId";
    j++;    cutFlowIdxNamesMap[cutFlowOffsets["muonSelection"] + j ] = "muon_isHighPurity";
    j++;    cutFlowIdxNamesMap[cutFlowOffsets["muonSelection"] + j ] = "muon_isNewSoftMVA";

    i+=offset;
    cutFlowOffsets["diMuonSelection"]=i ;
    j=0;
    j=1;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonSelection"] + j ] = "diMuon_m1Iso";
    j=2;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonSelection"] + j ] = "diMuon_m2Iso";
    j=3;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonSelection"] + j ] = "diMuon_mumuDr";
    j=4;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonSelection"] + j ] = "diMuon_mass";
    j=5;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonSelection"] + j ] = "diMuon_nMM";
    
    i+=offset;
    cutFlowOffsets["diMuonVertexSelection"]=i ;
    j=0;
    j=1;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonVertexSelection"] + j ] = "diMuonVtx_muMuDoca";
    j=2;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonVertexSelection"] + j ] = "diMuonVtx_vertexProbablility";
    j=3;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonVertexSelection"] + j ] = "diMuonVtx_trackDoca";
    j=4;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonVertexSelection"] + j ] = "diMuonVtx_nTrackClose";
    j=5;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonVertexSelection"] + j ] = "diMuonVtx_n1SigmaTrackClose";
    j=6;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonVertexSelection"] + j ] = "diMuonVtx_n2SigmaTrackClose";
    j=7;    cutFlowIdxNamesMap[cutFlowOffsets["diMuonVertexSelection"] + j ] = "diMuonVtx_n3SigmaTrackClose";
    
    i+=offset;
    cutFlowOffsets["photonSelection"]=i ;
    j=0;
    j=1;    cutFlowIdxNamesMap[cutFlowOffsets["photonSelection"] + j ] = "photon_Et";
    j=2;    cutFlowIdxNamesMap[cutFlowOffsets["photonSelection"] + j ] = "photon_Eta";
    j=3;    cutFlowIdxNamesMap[cutFlowOffsets["photonSelection"] + j ] = "photon_scMVA";
    j=4;    cutFlowIdxNamesMap[cutFlowOffsets["photonSelection"] + j ] = "photon_genMatch";
    
    i+=offset;
    cutFlowOffsets["bmmgSelection"]=i ;
    j=0;
    j=1;    cutFlowIdxNamesMap[cutFlowOffsets["bmmgSelection"] + j ] = "bmmg_PhotonPt";
    j=2;    cutFlowIdxNamesMap[cutFlowOffsets["bmmgSelection"] + j ] = "bmmg_PhotonEta";
    j=3;    cutFlowIdxNamesMap[cutFlowOffsets["bmmgSelection"] + j ] = "bmmg_DimuPhotonDr";
    j=4;    cutFlowIdxNamesMap[cutFlowOffsets["bmmgSelection"] + j ] = "bmmg_SvGammaDCA";
    
    for( auto &x : cutFlowIdxNamesMap) 
    {
            std::cout<<x.first<<" -> "<<x.second<<"\n";
    }
}

void BMMGAnalysis::FillCutFlow(Int_t cutIdx)
{
    //std::cout<<"Cut idx : "<<cutIdx<<"\n";
    th1fStore["CutFlowSummary"]->Fill(cutFlowIdxNamesMap[cutIdx],1);
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

    if(doRunLumiLog) runLumiLogger.saveAsJson((prefix + runLumiJsonName));
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
           getIntFromTag("ReportEvery"        ,strStream , field , reportEvery);
           getFloatFromTag("MuMuVtxProbabilityMin"	  ,strStream , field , muMuVtxProbabilityMin);
           getFloatFromTag("MuMuMaxDOCA"			  ,strStream , field , muMuMaxDOCA);
           getFloatFromTag("SvGammaMaxDCA"		      ,strStream , field , svGammaMaxDCA);
           getFloatFromTag("MinMuonPt"			      ,strStream , field , minMuonPt);
           getFloatFromTag("MinSCEt"			      ,strStream , field , minSCEt);
           getFloatFromTag("MaxSCEta"			      ,strStream , field , maxSCEta);
           getFloatFromTag("SvMinDocaTrack"			  ,strStream , field , svMinDocaTrack);
           getFloatFromTag("DrMaxForGenMatchSC"	      ,strStream , field , drMaxForGenMatchSC);
           getFloatFromTag("DrMaxForGenMatchMu"	      ,strStream , field , drMaxForGenMatchMu);
           getIntFromTag("SvMaxNTracksClose"        ,strStream , field , svMaxNTracksClose);
           getIntFromTag("SvMaxN1SigmaTracksClose"  ,strStream , field , svMaxN1SigmaTracksClose);
           getIntFromTag("SvMaxN2SigmaTracksClose"  ,strStream , field , svMaxN2SigmaTracksClose);
           getIntFromTag("SvMaxN3SigmaTracksClose"  ,strStream , field , svMaxN3SigmaTracksClose);
           getFloatFromTag("MinLeadMuIsolation"		  ,strStream , field , minLeadMuIsolation);
           getFloatFromTag("MinSubLeadMuIsolation"	  ,strStream , field , minSubLeadMuIsolation);
           getBoolFromTag ("MuonHasToBeTracker"		  ,strStream , field , muonHasToBeTracker);
           getBoolFromTag ("MuonHasToBeGlobal"		  ,strStream , field , muonHasToBeGlobal);
           getBoolFromTag ("MuonHasToBeHighPurity"	  ,strStream , field , muonHasToBeHighPurity);

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
            if(field.compare("DoRunLumiLog")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 doRunLumiLog= tmpI >0 ? 1 : 0;
                 cout<<" setting doRunLumiLog  = "<<doRunLumiLog<<"\n";
            }
           if(field.compare("RunLumiJsonName")==0){
                 getline(strStream, field);
                 runLumiJsonName=field;
                 std::cout<<" setting runLumiJsonName = "<<runLumiJsonName<<"\n";
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

void BMMGAnalysis::getBoolFromTag(string tag,std::istringstream &strStream, string inputStr  , Bool_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
           getline(strStream, field);
           var =  std::atoi(field.c_str()) > 0 ? 1 : 0;
           cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}

void BMMGAnalysis::getFloatFromTag(string tag,std::istringstream &strStream, string inputStr  , Float_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
             getline(strStream, field);
             var=std::atof(field.c_str());
             cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
       }
}

void BMMGAnalysis::getIntFromTag(string tag,std::istringstream &strStream, string inputStr  , Int_t &var)
{   
    std::string field;
    if(inputStr.compare(tag)==0){
             getline(strStream, field);
             var=std::atoi(field.c_str());
             cout<<" SETTING  "<<tag<<" = "<<var<<"\n";
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

