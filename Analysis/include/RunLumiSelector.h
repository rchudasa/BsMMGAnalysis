/*
 *      author : aravind t s 
 *      date   : 09 Dec 2021
 *      mail   : aravindsugunan@gmail.com
 */

#include "vector"
#include "map"
#include "TFile.h"
struct LumiSet
{
    int min,max;  
};

class RunLumiSelector {
    
    std::map<Int_t,vector<LumiSet> *> runLumiMask;
    bool loadedMask;
    public :
        RunLumiSelector()
        {
            loadedMask=false;
            
        }
        void loadRunLumiMask(TString fname, TString treeName="RunLumiMask");
        bool checkRunLumi(Int_t run, Int_t lumi);
        void printAllActiveLumis();

};

void RunLumiSelector::printAllActiveLumis()
{

       
 /*   for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
    {
        
        //std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &ahist = *(it->second); 
        ahist.Write();
    }
*/

    for(auto i : runLumiMask)
    {
        std::cout<<"\nRun : "<<i.first<<" : ";
        for(auto j :  *(i.second))
        {
                std::cout<<"["<<j.min<<","<<j.max<<"] , ";
        }       
    }
}

void RunLumiSelector::loadRunLumiMask(TString fname,TString treeName)
{
    TFile file(fname,"READ");
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

