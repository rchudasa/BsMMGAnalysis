/*
 *      author : aravind t s 
 *      date   : 09 Dec 2021
 *      mail   : aravindsugunan@gmail.com
 */


#include "vector"
#include "map"

struct LumiSet
{
    int min,max;  
};

class RunLumiSelector {
    
    std::map<TString,vector<LumiSet> *> runLumiMask;
    bool loadedMask;
    public :
        RunLumiSelector()
        {
            loadedMask=false;
            
        }
        void loadRunLumiMask(TString fname, TString treeName);
        bool checkRunLumi(Int_t run, Int_t lumi);

};

void RunLumiSelector::loadRunLumiMask(TString fname,TString treeName)
{
    auto file=TFile(fname,"READ");
    auto runLumiTree=(TTree * ) file.Get(treeName);
    Int_t run;
    Int_t nLumiSets;

    Int_t minLumi[200];
    Int_t maxLumi[200];
    
    runLumiTree->SetBranchAddress("run"      , &run );
    runLumiTree->SetBranchAddress("nLumiSets", &nLumiSets);
    runLumiTree->SetBranchAddress("minLumi"  , minLumi);
    runLumiTree->SetBranchAddress("maxLumi"  , maxLumi);
    
    Long64_t nentries = runLumiTree->GetEntries();
    
    for(Int_t i=0; i < nentries;i++)
    {
        runLumiTree->GetEntry(i);
        runLumiMask[run]=new std::vector<LumiSet>;
 //       std::cout<<"Adding run : "<<run<<"\n";    
        for(Int_t j=0;j<nLumiSets;j++)
        {
            
            LumiSet aLumiSet;
            aLumiSet.min=minLumi[j];
            aLumiSet.max=maxLumi[j];
//            std::cout<<"\t [ "<<aLumiSet.min<<" , " <<aLumiSet.max<<" ] "<<"\n";
            runLumiMask[run]->push_back(aLumiSet);
        }
        
    }

    file.Close();
    loadedMask=true;
    
}

bool RunLumiSelector::checkRunLumi(Int_t run,Int_t lumi)
{
    
    if ( not loadedMask) {
    return false;
    }
    if ( runLumiMask.find(run) == runLumiMask.end() ){ 
    return false;
    }
    for( Int_t i=0;i<runLumiMask[run]->size();i++)
    {   

        if( lumi >= (runLumiMask[run]->at(i)).min && lumi<= (runLumiMask[run]->at(i)).max )
            return true;
    }
    return false;
}

