/*

usage : root -b -q 'getEventListFromBMM5.cc(\"@@CFGFILENAME\")'

////   Config

#FLIST_BEG
/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/Charmonium+Run2018A-12Nov2019_UL2018_rsb-v1+MINIAOD/05906F81-6ECC-B34A-BADC-5C47DD466F33.root
#FLIST_END

#PARAMS_BEG
OutputPrefix=
OutputFile=bmm5EvntSelection_0.root
OutputFileTxt=bmm5EvntSelection_0.txt
#PARAMS_END

#EOF


*/
#include "genericTreeBranchSelector.h"
R__LOAD_LIBRARY(genericTreeBranchSelector_C.so)



void doMerging (std::vector<string> fMotherList,
                std::vector<string> fSonList,
                ULong64_t maxEvent,
                string ofilename,
                string prifix
                );


void mergeBmmXTrees(string fname)
{
	fstream cfgFile(fname,ios::in);
	string line;
	bool cfgModeFlag=false;
	
	std::vector<string> fMotherList;
	std::vector<string> fSonList;
	
//	for(int i=0;i<fList.size();i++)
//	{
//		std::cout<<"  "<<i+1<<" "<<fList.at(i)<<"\n";
//	}

    ULong64_t maxEvent(10);
    string ofilename("output.root");
    string Bmm5FileIndexFname("output.txt");
    string prefix("");
    
    cfgFile.clear();
    cfgFile.seekg(0,ios::beg);
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
    std::cout<<"Reading Params "<<"\n";
    int gg=0;
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
                 ofilename=field;
                 std::cout<<" setting ofilename = "<<ofilename<<"\n";
            }
            if(field.compare("OutputPrefix")==0){
                 getline(strStream, field);
                 prefix=field;
                 std::cout<<" setting prefix = "<<prefix<<"\n";
            }
            if(field.compare("MaxEvents")==0){
                 getline(strStream, field);
                 maxEvent=std::atoi(field.c_str());
                 std::cout<<" setting maxEvent  = "<<maxEvent<<"\n";
            }
       }
	}

	// getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	while(std::getline(cfgFile,line))
	{
	   if(line=="#MOTHER_FILELIST_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#MOTHER_FILELIST_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   //std::cout<<line<<"|\n";
	   fMotherList.push_back(line);
	}

	std::cout<<"Read Mother List with  "<<fMotherList.size()<<" entries \n";
	
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	while(std::getline(cfgFile,line))
	{
	   if(line=="#SON_FILELIST_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#SON_FILELIST_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   //std::cout<<line<<"|\n";
	   fSonList.push_back(line);
	}

	std::cout<<"Read Son List with  "<<fSonList.size()<<" entries \n";

	doMerging(fMotherList,fSonList,maxEvent,ofilename,prefix);
	
}

void doMerging(   std::vector<string> fMotherList,
                  std::vector<string> fSonList,   
                  ULong64_t maxEvent,
                  string ofilename,
                  string prefix
              )
{

 
  // TODO
  //    Write the write functionality to the new made tree
  

  UInt_t          run(0),s_run(0);
  UInt_t          lumi(0),s_lumi(0);
  ULong64_t       event(0),s_event(0);
  TString         filename;
 


  TFile * mergedTreeFile = new  TFile((prefix+ofilename).c_str(),"recreate");

   // Setting up the TChain for the BMMG Ntuples
   TChain *motherTreeChain = new TChain("Ntuples/EventTree");
   for (int i=0;i<fMotherList.size();i++)
   {
	    motherTreeChain->Add(fMotherList.at(i).c_str());
   }

   // Setting up the TChain for the bmm5 Ntuples

   TChain *searchTreeChain = new TChain("Events");
   for (int i=0;i<fSonList.size();i++)
   {
	    searchTreeChain->Add(fSonList.at(i).c_str());
   }
   
  genericTreeBranchSelector treeBranchSelector(motherTreeChain,searchTreeChain);

   treeBranchSelector.b5.fChain->BuildIndex("run","event");
   treeBranchSelector.bG.fChain->BuildIndex("run","event");
  
   if(treeBranchSelector.b5.fChain == 0 ) { std::cout<<" source fchain not found" ; return ; }
   if(treeBranchSelector.bG.fChain == 0 ) { std::cout<<" search fchain not found !! "; return   ; }

  Long64_t nentries = treeBranchSelector.bG.fChain->GetEntries();
  cout<<"Total Available Entries :  "<< nentries << endl;
  nentries = nentries < maxEvent ? nentries : maxEvent;
  cout<<"Processing Total Entries :  "<< nentries << endl;

  Long64_t jentry=0,eventCount=0,evt=0;
  Long64_t nb = 0,nbytes=0 ;
  
  for ( ; jentry<nentries;jentry++) 
    {
      //if (jentry%2500==0) 
        { 
	        cout << "--> " << jentry << "/" << nentries<<endl;
        }
      Long64_t ientry_evt = treeBranchSelector.bG.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = treeBranchSelector.bG.fChain->GetEntry(jentry);   nbytes += nb;
       
      run   = treeBranchSelector.bG.run;
      lumi   = treeBranchSelector.bG.lumis;
      event = treeBranchSelector.bG.event;

      nb = treeBranchSelector.b5.fChain->GetEntryWithIndex(run,event);

      s_run   = treeBranchSelector.b5.run;
      s_lumi = treeBranchSelector.b5.luminosityBlock;
      s_event = treeBranchSelector.b5.event;

      std::cout<<jentry<<"  FLG : "<< bool(nb>0) <<" :  nb : "<<nb<<"  BMMG {  "<<run<<","<<lumi<<","<<event<<" }  ,{ "<<s_run<<"."<<s_lumi<<","<<s_event<<"}\n";
      if ( nb < 0 ) continue;

      treeBranchSelector.Fill();
        
    }
 
  mergedTreeFile->cd();
  treeBranchSelector.Write();
  mergedTreeFile->Close();

  return;


}

