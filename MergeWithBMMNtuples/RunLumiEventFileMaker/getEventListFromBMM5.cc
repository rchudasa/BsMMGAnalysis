#include "Bmm5Ntuple.h"

Double_t getDETA(Double_t eta1, Double_t eta2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
void getRunDetails(std::vector<string> fList,
                string ofilename,
                string ofilenameTxt,
                string prifix,
                Long64_t maxEvents
                );

void getEventListFromBMM5(string fname)
{
	fstream cfgFile(fname,ios::in);
	string line;
	
	std::vector<string> fList;
	
	// getting flist
	cfgFile.seekp(ios::beg);
	bool cfgModeFlag=false;
	while(std::getline(cfgFile,line))
	{
	   if(line=="#FLIST_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#FLIST_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   std::cout<<line<<"|\n";
	   fList.push_back(line);
	}

	std::cout<<"Read File List as : "<<"\n";
	for(int i=0;i<fList.size();i++)
	{
		std::cout<<"  "<<i+1<<" "<<fList.at(i)<<"\n";
	}

    string ofilename("output.root");
    string ofilenameTxt("output.txt");
    string prefix("");
    Long64_t maxEvents(10);
    bool doSelection(0);
    
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
            if(field.compare("OutputFileTxt")==0){
                 getline(strStream, field);
                 ofilenameTxt=field;
                 std::cout<<" setting ofilename = "<<ofilename<<"\n";
            }
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
                 maxEvents=std::atoi(field.c_str());
                 cout<<" setting maxEvents  = "<<maxEvents<<"\n";
            }
            if(field.compare("DoSelection")==0){
                 getline(strStream, field);
                 if(std::atoi(field.c_str())==1)  doSelection=1;
                 else doSelection = 0;
                 cout<<" setting doSelection  = "<<doSelection<<"\n";
            }
       }
	}


//	getRunDetails(fList,ofilename,ofilenameTxt,prefix, maxEvents);
//	
//}
//
//void getRunDetails(   std::vector<string> fList,
//                           string ofilename,
//                           string ofilenameTxt,
//                           string prefix,
//                           Long64_t maxEvents
//                   )
//{

  Long64_t NMAX=1e9;
  
  TChain *treeChain = new TChain("Events");
  TFile * ofile= new TFile((prefix+ofilename).c_str(),"recreate");

   for (int i=0;i<fList.size();i++)
   {
	    treeChain->Add(fList.at(i).c_str());
   }

  Bmm5Ntuple ntupleRawTree(treeChain);

  if (ntupleRawTree.fChain == 0) return;

  UInt_t          run(0),run_o(1);
  UInt_t          lumi(0),lumi_o(1);
  ULong64_t       event;
  TString         filename;
  
  TTree * runLumiTree = new  TTree("RunLumi","Run LumiBlock and Filenames");
  TFile * runLumiFile = new  TFile((prefix+ofilename).c_str(),"recreate");
  fstream runLumiTxt((prefix+ofilenameTxt).c_str(),ios::out);

  runLumiTree->Branch("run",&run);
  runLumiTree->Branch("lumi",&lumi);
  runLumiTree->Branch("event",&event);
  runLumiTree->Branch("file",&filename);

  

  // need to setup proper Branches for SPEEDUP TODO
  ntupleRawTree.fChain->SetBranchStatus("*",1);
  
  Long64_t nentries = ntupleRawTree.fChain->GetEntries();
  std::cout<<"nentries : "<<nentries<<"\n";
  if (maxEvents >0 ) nentries = nentries > maxEvents ? maxEvents : nentries;
  cout<<"Processing Total Entries :  "<< nentries << endl;
  
  Long64_t jentry=0,eventCount=0,evt=0;
  Long64_t nb = 0,nbytes=0 ;
 
  UInt_t mu1_Idx,mu2_Idx;

  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();

  for (; jentry<nentries;jentry++) 
    {

      if (jentry%2500==0) 
      {
         t_end = std::chrono::high_resolution_clock::now();
         cout<<"Processing Entry in event loop : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/jentry* 0.001
             <<"\n";
        
      }
      
      Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
      if (ientry_evt < 0) break;
      evt++;
      nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;

     run=ntupleRawTree.run;
     lumi=ntupleRawTree.luminosityBlock;
     event=ntupleRawTree.event;
     filename=(treeChain->GetFile())->GetName();
     
     eventCount++;
     runLumiTree->Fill();
     if(run_o != run or lumi_o != lumi) {
        runLumiTxt<<"\n@ "<<run<<","<<lumi<<","<<filename<<"\n\t";
        runLumiTxt<<event;
        run_o=run;
        std::cout<<run<<" , "<<lumi<<" : "<< eventCount<<" / "<<evt <<"\n";
        eventCount=0;
        lumi_o=lumi;
     }
     else{
        //runLumiTxt<<","<<event;
    }
  }

  runLumiTree->Write();
  runLumiFile->Purge();
  runLumiFile->Write();
  runLumiFile->Close();
  runLumiTxt.close();

}


Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;
  
  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;
  
  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }
  
  return TMath::Abs(dphi);
  //return dphi;
}



Double_t getDETA(Double_t eta1, Double_t eta2){
  return TMath::Abs(eta1 - eta2);
}

