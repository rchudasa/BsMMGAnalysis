/*
 *      author : aravind t s 
 *      date   : 09 Dec 2021
 *      mail   : aravindsugunan@gmail.com
 */

#include "set"
#include "map"
#include "fstream"

class RunLumiLogger {
    
    std::map<Int_t,set<int, less<int> >> runLumiMap;
    public :
        RunLumiLogger()
        {
            
        }
        void clear()
        {   
            runLumiMap.clear();
        }
        void addRunLumi(Int_t run, Int_t lumi);
        void printAllRunLumis();
        void printJson();
        void saveAsJson(string);

};

void RunLumiLogger::printAllRunLumis()
{
    for(auto i : runLumiMap)
    {
        std::cout<<"\nRun : "<<i.first<<" : ";
        for(auto j :  i.second )
        {   
            std::cout<<j<<",";
        }          
        std::cout<<"\n";
    }

}
void RunLumiLogger::addRunLumi(Int_t run,Int_t lumi)
{
       runLumiMap[run].insert(lumi);
}

void RunLumiLogger::printJson()
{
    std::cout<<"{";
    int lastSetBegLumi=-1;
    int lastSetEndLumi=-1;
    string runPrefix="";
    string lumiPrefix="";
    for(auto i : runLumiMap)
    {
        std::cout<<runPrefix<<"\n  \""<<i.first<<"\": [\n";
        runPrefix=",";

        lastSetBegLumi=-1;
        lastSetEndLumi=-1;
        lumiPrefix="";
        for(auto j :  i.second )
        {   
            if( ( lastSetEndLumi+1 ) != j)
            {
                if(lastSetBegLumi!=-1) {
                    std::cout<<lumiPrefix<<"\n      [\n        "<<lastSetBegLumi<<",\n        "<<lastSetEndLumi<<"\n      ]";
                    lumiPrefix=",";
                  }
                lastSetBegLumi=j;
            }
            lastSetEndLumi=j;
        }          
        if(lastSetBegLumi!=-1) std::cout<<lumiPrefix<<"\n      [\n        "<<lastSetBegLumi<<",\n        "<<lastSetEndLumi<<"\n      ]\n";
        std::cout<<"  ]";
    }
    std::cout<<"\n}\n";
}

void RunLumiLogger::saveAsJson(string fname)
{
    fstream fileOut(fname,ios::out);
    fileOut<<"{";
    int lastSetBegLumi=-1;
    int lastSetEndLumi=-1;
    string runPrefix="";
    string lumiPrefix="";
    for(auto i : runLumiMap)
    {
        fileOut<<runPrefix<<"\n  \""<<i.first<<"\": [\n";
        runPrefix=",";

        lastSetBegLumi=-1;
        lastSetEndLumi=-1;
        lumiPrefix="";
        for(auto j :  i.second )
        {   
            if( ( lastSetEndLumi+1 ) != j)
            {
                if(lastSetBegLumi!=-1) {
                    fileOut<<lumiPrefix<<"\n      [\n        "<<lastSetBegLumi<<",\n        "<<lastSetEndLumi<<"\n      ]";
                    lumiPrefix=",";
                  }
                lastSetBegLumi=j;
            }
            lastSetEndLumi=j;
        }          
        if(lastSetBegLumi!=-1) fileOut<<lumiPrefix<<"\n      [\n        "<<lastSetBegLumi<<",\n        "<<lastSetEndLumi<<"\n      ]\n";
        fileOut<<"  ]";
    }
    fileOut<<"\n}\n";
    fileOut.close();
}

