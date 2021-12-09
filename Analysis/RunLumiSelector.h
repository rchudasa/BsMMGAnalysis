
#include "vector"
#include "map"

struct LumiSet
{
    int min,max;  
};

class RunLumiSelector {
    
    std::map<TString,vector<LumiSet>> runLumiMask;
    bool loadedMask;
    public :
        RunLumiSelector()
        {
            loadedMask=false;
            
        }
        void loadRunLumiMask(TString fname);
        bool checkRunLumi(Int_t run, Int_t lumi);

}

void RunLumiSelector::loadRunLumiMask(TString fname,TString treeName)
{
    auto file=TFile(fname,"READ");
    auto runLumiTree=(TTree * ) file.Get(treeName);
    Int_t run;
    Int_t nLumiSets;

    Int_t minLumi[200];
    Int_t maxLumi[200];
    
    runLumiTree->Branch("run"      , &run );
    runLumiTree->Branch("nLumiSets", &nLumiSets);
    runLumiTree->Branch("minLumi"  , minLumi);
    runLumiTree->Branch("maxLumi"  , maxLumi);
    
    Long64_t nentries = runLumiTree->GetEntries();
    
    for(Int_t i=0; i < nentries;i++)
    {
        runLumiTree->GetEntry(i);
        runLumiMask[run]=std::vector<LumiSet>;
        
        for(Int_t j=0;j<nLumiSets;j++)
        {
            LumiSet aLumiSet;
            aLumiSet.min=minLumi[j];
            aLumiSet.max=minLumi[j];

            runLumiMask[run].push_back(aLumiSet);
        }
        
    }

    file.Close();
    loadedMask=true;
    
}

bool RunLumiSelector::checkRunLumi(Int_t run,Int_t lumi)
{

    if ( not loadedMask)  return false;
    if ( runLumiMask.find(run) == runLumiMask.end() ) return false;
    
    for( Int_t i=0;i<runLumiMask[run].size();i++)
    {
        if( lumi >= runLumiMask[run][i].min && lumi<= runLumiMask[run][i].max )
            return true;
    }
    return false;
}

