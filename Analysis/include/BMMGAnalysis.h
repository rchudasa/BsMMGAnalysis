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
#define NDIMU_MAX 200
#define NBMMG_MAX 100
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
    std::map<TString,Float_t> storageFloat;
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
    std::map<TString,TTree*> treesToStore;
    Long64_t nentries, maxEvents ;
    Int_t reportEvery;
    
    // Dimuon MVA vars 
    bool doDimuonMVA;
    bool doDimuonMVASelection;
    bool hasDimuonWeightFiles;
    bool hasSetupDimuonMVA;
    TString dimuonMVAWeightFile;
    Float_t dimuonMVAThreshold;
    
    // Photon MVA vars 
    bool doPhotonMVA;
    bool hasWeightFiles;
    bool hasSetupPhotonMVA;
    TMVA::Reader *reader ;
    std::map<string , TMVA::Reader*> readerStore ;
    PhotonMVAvariables photonMVAdata;
    TString photonIdxMVAWeightFile;
    TString photonIdxECapMVAWeightFile;
    TString photonIdxBarrelMVAWeightFile;
    Float_t photonMVAValue;
    Float_t photonIDcut;
    Float_t photonIDcutBarrel;
    Float_t photonIDcutECap;

    

    // BMMG MVA Vars
    bool doBMMGMVA;
    bool doBMMGMVASelection;
    bool hasBMMGWeightFiles;
    bool hasSetupBMMGMVA;
    TString bmmgMVAWeightFile;
    Float_t bmmgMVAThreshold;

    Float_t mmgSideBandRBeg ; 
    Float_t mmgSideBandREnd ;
    Float_t mmgSideBandLBeg ;
    Float_t mmgSideBandLEnd ;
    Float_t mmgSigRegionBeg ;
    Float_t mmgSigRegionEnd ;

    std::vector<Float_t> mvaControlRegionSplits;
    
    void setUpBMMGMVA();
    Float_t doBMMGMVAScores(Int_t mumuIdx,Int_t scIdx);
    
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
    Float_t minMuonEta     ;
    Float_t maxMuonEta     ;
    Float_t minSCEt        ;
    Float_t maxSCEta       ;
    Float_t minLeadMuIsolation;
    Float_t minSubLeadMuIsolation;
    Double_t maxMuMuDr          ;
    Double_t maxDimuPhotonDr    ;
    Double_t maxMMGMass         ;
    Double_t minMMGMass         ;
    std::vector<Double_t> minBMMGMass;
    std::vector<Double_t> maxBMMGMass;
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
    void treeFill(TString tag);
    void SaveFile();
    void SetupAnalysis(bool makeTreeCopy=false);
    void setupOutPuts(bool makeTreeCopy=false);
    void SetupCutFlowOffsets();
    Int_t getMuonMatch(Double_t muEta,Double_t muPhi);
    // Histogram Related Functions
    void bookHistograms();
    void bookBMMGHists(TString tag);

    Double_t getDCAGammaToDimuVertex(Int_t mumuIdx,Int_t phoId);
    void FillCutFlow(Int_t);
    void fill_muonHists(Int_t idx=-1);
    void fill_scHists(Int_t idx =-1);
    void fill_photonHists(Int_t idx  = -1);
    void fill_dimuonPassHists(Int_t idx=-1);
    void fill_dimuonHists(Int_t mumuIdx=-1);
    void fill_dimuonEnvironmentHists(Int_t mumuIdx);
    void fill_bmmgHists(TLorentzVector &bmmgLV,Int_t mumuIdx, Int_t phoSCIdx,TString tag="bmmg_");
    void fill_globalEventHists();

    //void getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose);
    //void getFloatFromTag(string tag,std::istringstream &strStream, string inputStr ,Float_t &var);
    //void getIntFromTag(string tag,std::istringstream &strStream, string inputStr ,Int_t &var);
    //void getBoolFromTag(string tag,std::istringstream &strStream, string inputStr ,Bool_t &var);
    // Analysis Functions
    void doPhotonMVAScores();
    void setUpPhotonMVA();
    
    void doMuonMVAScores();
    
    void setUpDimuonMVA();
    void doDimuonMVAScores();

    void FillMVAControlHistograms(Float_t mvaVal,Float_t mmgMass);
    void FillDimuMVAControlHistograms(Float_t mvaVal,Int_t mumuIdx);
    
    bool doTriggerSelection();
    Int_t doMuonSelection(Int_t muIdx, bool isLead);
    Int_t doDimuonSelection(Int_t mumuIdx);
    Int_t doVertexSelection(Int_t mumuIdx);
    Int_t doBMMGSelection(TLorentzVector &bmmgLV,Int_t mumuIdx , Int_t phoSCIdx);
    Int_t doPhotonSelection(Int_t scIdx);
    Int_t getPVMatch(Int_t mumuIdx);
    void setupReducedAnalysisTreeBranches();
    void SkimData();
    void MVAAnalyze();
    void Analyze();
    
    void makeMVATree();
    void fillBMMGTreeVars(Int_t mumuIdx,Int_t phoSCIdx);
    void AddBMMGmvaTree(TString);
    
    void makeDimuonMVATree();
    void makeGenDimuonMVATree();
    void fillDimuonVertexVars( Int_t mumuIdx);
    void AddDimuonTree(TString tag);


    #ifdef __MCANALYSIS__
    void GenAnalyze();
    void makeGenMVATree();
    #endif

    // Temporary vars with class scope
    TLorentzVector mu1LV,mu2LV,bmmgLV,diMuLV,photonLV;
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
    reader=nullptr;
    readerStore["photonIDMVA"]=nullptr;
    readerStore["dimuonMVA"]=nullptr;
    dimuonMVAThreshold=-2.0;


    mmgSideBandRBeg  =  0.0 ;  
    mmgSideBandREnd  =  4.5 ;
    mmgSideBandLBeg  =  4.5 ;
    mmgSideBandLEnd  =  6.22 ;
    mmgSigRegionBeg  =  6.22 ;
    mmgSigRegionEnd  =  10.0 ;

    minMuonPt          = -1.0;
    minMuonEta         = -1e9;
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

    doTriggerFiltering=false;
    
    doPhotonMVA        =false;
    photonIdxMVAWeightFile="";
    photonIdxBarrelMVAWeightFile="";
    photonIdxECapMVAWeightFile="";
    hasWeightFiles=false;
    photonIDcutBarrel=-1.5;
    photonIDcutECap=-1.5;

    doDimuonMVA=false;
    doDimuonMVASelection=false;
    dimuonMVAWeightFile="";
    hasDimuonWeightFiles=false;
    dimuonMVAThreshold=-2.0;
    hasSetupDimuonMVA=false;
    
    doBMMGMVA=false;
    doBMMGMVASelection=false;
    bmmgMVAWeightFile="";
    hasBMMGWeightFiles=false;
    bmmgMVAThreshold=-2.0;
    hasSetupBMMGMVA=false;
    
    doRunLumiMask=false;
    doReducedTree=false;
    doRunLumiLog=false;
    mvaControlRegionSplits.clear() ; 
        mvaControlRegionSplits.push_back(-1.0);
        mvaControlRegionSplits.push_back(1.0);
    minBMMGMass.clear(); minBMMGMass.push_back(-1.0);
    maxBMMGMass.clear(); maxBMMGMass.push_back(1e9);
    minDimuMass.clear(); minDimuMass.push_back(0.0);
    maxDimuMass.clear(); maxDimuMass.push_back(6.0);

    readParameters(cfgFileName);
    
    if(doDimuonMVA)
    {
       setUpDimuonMVA();
    }
    if(doBMMGMVA)
    {
       setUpBMMGMVA();
    }
 
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
   cutFlowIdxNamesMap[0]="sucess"; 
    // Need to 
    i=0;
    cutFlowOffsets["evtSelection"] =i;
    j=1;    cutFlowIdxNamesMap[cutFlowOffsets["evtSelection"] + j ] = "allEventSelection";
    j=2;    cutFlowIdxNamesMap[cutFlowOffsets["evtSelection"] + j ] = "notGoodEvent";

    i+=offset;
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

void BMMGAnalysis::treeFill(TString tag)
{
    treesToStore[tag]->Fill();
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

    for (std::map<TString,TTree *>::iterator it=treesToStore.begin() ; it!=treesToStore.end(); ++it)
    {
        
        //std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &a= *(it->second); 
        a.Write();
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
           getFloatFromTag("MinMuonEta"			      ,strStream , field , minMuonEta);
           getFloatFromTag("MaxMuonEta"			      ,strStream , field , maxMuonEta);
           getFloatFromTag("MinSCEt"			      ,strStream , field , minSCEt);
           getFloatFromTag("MaxSCEta"			      ,strStream , field , maxSCEta);
           getFloatFromTag("SvMinDocaTrack"			  ,strStream , field , svMinDocaTrack);
           getFloatFromTag("DrMaxForGenMatchSC"	      ,strStream , field , drMaxForGenMatchSC);
           getFloatFromTag("DrMaxForGenMatchMu"	      ,strStream , field , drMaxForGenMatchMu);
           getIntFromTag("SvMaxNTracksClose"          ,strStream , field , svMaxNTracksClose);
           getIntFromTag("SvMaxN1SigmaTracksClose"    ,strStream , field , svMaxN1SigmaTracksClose);
           getIntFromTag("SvMaxN2SigmaTracksClose"    ,strStream , field , svMaxN2SigmaTracksClose);
           getIntFromTag("SvMaxN3SigmaTracksClose"    ,strStream , field , svMaxN3SigmaTracksClose);
           getFloatFromTag("MinLeadMuIsolation"		  ,strStream , field , minLeadMuIsolation);
           getFloatFromTag("MinSubLeadMuIsolation"	  ,strStream , field , minSubLeadMuIsolation);
           getBoolFromTag("MuonHasToBeTracker"		  ,strStream , field , muonHasToBeTracker);
           getBoolFromTag("MuonHasToBeGlobal"		  ,strStream , field , muonHasToBeGlobal);
           getBoolFromTag("MuonHasToBeLoose"		  ,strStream , field , muonHasToBeLoose);
           getBoolFromTag("MuonHasToBeHighPurity"	  ,strStream , field , muonHasToBeHighPurity);
           getBoolFromTag("DoReducedTree"	          ,strStream , field , doReducedTree);
           getBoolFromTag("DoDimuonMVA"	              ,strStream , field , doDimuonMVA);
           getBoolFromTag("DoDiMuonMVASelection"      ,strStream , field , doDimuonMVASelection);
           getBoolFromTag("DoBMMGMVA"	              ,strStream , field , doBMMGMVA);
           getBoolFromTag("DoBMMGMVASelction"	      ,strStream , field , doBMMGMVASelection);
           getFloatFromTag("DimuonMVAThreshold"	      ,strStream , field , dimuonMVAThreshold);
           getFloatFromTag("BmmgMVAThreshold"	      ,strStream , field , bmmgMVAThreshold);
           //getTStringFromTag("BMMGMVAWeightFile"	  ,strStream , field , bmmgMVAWeightFile);
           getFloatFromTag("MmgSideBandRBeg"	      ,strStream , field , mmgSideBandRBeg);
           getFloatFromTag("MmgSideBandREnd"	      ,strStream , field , mmgSideBandREnd);
           getFloatFromTag("MmgSideBandLBeg"	      ,strStream , field , mmgSideBandLBeg);
           getFloatFromTag("MmgSideBandLEnd"	      ,strStream , field , mmgSideBandLEnd);
           getFloatFromTag("MmgSigRegionBeg"	      ,strStream , field , mmgSigRegionBeg);
           getFloatFromTag("MmgSigRegionEnd"	      ,strStream , field , mmgSigRegionEnd);
           getFloatFromTag("PhotonIDcutBarrel"	      ,strStream , field , photonIDcutBarrel);
           getFloatFromTag("PhotonIDcutECap"	      ,strStream , field , photonIDcutECap);

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
             if(field.compare("BmmgMVAWeightFile")==0){
                 hasBMMGWeightFiles=true;
                 getline(strStream, field);
                 bmmgMVAWeightFile=field;
                 std::cout<<" setting BmmgMVAWeightFile = "<<bmmgMVAWeightFile<<"\n";
            }
             if(field.compare("DimuonMVAWeightFile")==0){
                 hasDimuonWeightFiles=true;
                 getline(strStream, field);
                 dimuonMVAWeightFile=field;
                 std::cout<<" setting DimuonMVAWeightFile = "<<dimuonMVAWeightFile<<"\n";
            }
             if(field.compare("PhotonIDWeightFile")==0){
                 hasWeightFiles=true;
                 getline(strStream, field);
                 photonIdxMVAWeightFile=field;
                 std::cout<<" setting photonIdxMVAWeightFile = "<<photonIdxMVAWeightFile<<"\n";
            }
             if(field.compare("PhotonIDBarrelWeightFile")==0){
                 hasWeightFiles=true;
                 getline(strStream, field);
                 photonIdxBarrelMVAWeightFile=field;
                 std::cout<<" setting photonIdxBarrelMVAWeightFile = "<<photonIdxBarrelMVAWeightFile<<"\n";
            }
             if(field.compare("PhotonIDECapWeightFile")==0){
                 hasWeightFiles=true;
                 getline(strStream, field);
                 photonIdxECapMVAWeightFile=field;
                 std::cout<<" setting photonIdxECapMVAWeightFile = "<<photonIdxECapMVAWeightFile<<"\n";
            }
            
            if(field.compare("MvaControlRegionSplits")==0){
                 mvaControlRegionSplits.clear();
                 cout<<" setting mvaControlRegionSplits  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    mvaControlRegionSplits.push_back(aDouble);
                 }
                 cout<<" }\n";
            }            

           if(field.compare("MinBMMGMass")==0){
                 minBMMGMass.clear();
                 cout<<" setting minBMMGMass  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    minBMMGMass.push_back(aDouble);
                 }
                 cout<<" }\n";
            }            
            if(field.compare("MaxBMMGMass")==0){
                 maxBMMGMass.clear();
                 cout<<" setting maxBMMGMass  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    maxBMMGMass.push_back(aDouble);
                 }
                 cout<<" }\n";
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


#include "../src/BMMGAnalysis.cc"
#include "../src/BMMGmvaTreeMaker.cc"


