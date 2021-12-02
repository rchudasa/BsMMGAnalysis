#include "EventTree.h"


#define BS_MASS 5.366
void getIt()
{

  TChain *nTupleTree = new TChain("Ntuples/EventTree");
  nTupleTree->Add("/home/aravind/store/BsToMMG/ntuples/bmmg/muonNtuplizer_bsmmg_pVInfo.root");

  cout << "file " << nTupleTree->GetEntries()  << endl;
  EventTree ntupleRawTree(nTupleTree);
  ntupleRawTree.fChain->SetBranchStatus("*",1);
  
  TFile *output = new TFile("bSRapidityPsudoRapidity.root","recreate");
 
  auto eta   = TH1D("eta","eta",400,-5.0,5.0);
  auto rapi  = TH1D("rapi","rapi",400,-5.0,5.0);
  auto pz  = TH1D("pz","pz",100,0.0,100.0);
  auto e  = TH1D("e","e",40,0.0,100.0);
  auto p  = TH1D("p","p",40,0.0,100.0);
  
  if (ntupleRawTree.fChain == 0) return;
  Long64_t nentries = ntupleRawTree.fChain->GetEntriesFast();
  cout << nentries << endl;
  

 Long64_t entryMax=1e5;
 Long64_t jentry=0;
 Long64_t nbytes(0) ;
 Long64_t nb(0) ;
 TLorentzVector aVec;

 for (; jentry<nentries;jentry++) 
    {
      if (jentry>entryMax) break;
      if (jentry%1000==0) 
      {
        std::cout<<"doing : "<<jentry<<"\n";

      }
      
      Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;

      if(  ntupleRawTree.gen_Bs_pt->size() != 1 ) continue;
      if(  ntupleRawTree.nMC  != 3 ) continue;
     
        aVec.SetPtEtaPhiM(ntupleRawTree.gen_Bs_pt->at(0),ntupleRawTree.gen_Bs_eta->at(0),ntupleRawTree.gen_Bs_phi->at(0), BS_MASS );
        
        eta.Fill(aVec.Eta());
        rapi.Fill(aVec.Rapidity());
        pz.Fill(abs(aVec.Pz()));
        p.Fill(aVec.P());
        e.Fill(aVec.E());
    }
    output->cd();
    eta.Write();
    rapi.Write();
    pz.Write();
    p.Write();
    e.Write();
    output->Purge();
    output->Write();
    output->Close();
}
