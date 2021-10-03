/*
 *
 *          USAGE
 *
 * root -b -q 'bmmX_distributionStudy.cc("YOURCONFIG.CFG")'
 *
 *
 * */

#include "MergedBMMX.h"
R__LOAD_LIBRARY(MergedBMMX_C.so)

const Double_t TWO_PI(3.141592653589*2) ;
const Double_t PI(3.141592653589)       ;

#define PHO_MASS 0.0

Int_t setupDimuonBranches( TTree* outTree ,std::map<string, Int_t > &candidateMap , double * storageArray , Int_t offset=0 );
void assigndimuonBMMGCandidates( std::map<string,Int_t > &candidateMap, MergedBMMX &ntupleRawTree, int candIdx, double * storageArray);

Int_t setupPhotonSCBranches( TTree* outTree ,std::map<string, Int_t > &candidateMap , double * storageArray, Int_t offset=0 );
void assignSCBMMGCandidates( std::map<string,Int_t > &candidateMap, MergedBMMX &ntupleRawTree, int candIdx, double * storageArray);

Double_t getDR( Double_t eta1, Double_t phi1,Double_t eta2 ,Double_t phi2) ;
Double_t getDPHI( Double_t phi1, Double_t phi2) ;
Double_t getDETA(Double_t eta1, Double_t eta2) ;


void  bmmX_distributionStudy(string fname)
{
	fstream cfgFile(fname,ios::in);
	string line;
	bool cfgModeFlag=false;
	
	std::vector<string> InFileList;
	std::vector<string> branchList;
	std::vector<double> branchMin;
	std::vector<double> branchMax;
	std::vector<int> branchNBins;
	
    Double_t aDouble;
    
    const int nq=4;
    const int NBINS=20;
    int nbins,n;
    Double_t xq[nq]={0.01,0.98,0.25,0.75};  // position where to compute the quantiles in [0,1]
    Double_t yq[nq];              // array to contain the quantiles


    Long64_t maxEvents(10);
    string ofileName("output.root");
    string   treeName("mergedTree");
    string prefix("");
    double PhoDimuDrMax(1.4);
    double PhoPtMin(3.99);
        
    /**********************************************/
    
                    #include "readConfig.h"
    
    /**********************************************/
    
    std::map<string, TH1D * > histMap;
    
    TChain *treeChain = new TChain(treeName.c_str());
    for(auto i=0;i<InFileList.size();i++)
    {
        treeChain->Add(InFileList[i].c_str());
    }
    
    TFile * outputFile = new  TFile((prefix+ofileName).c_str(),"recreate");    
    

    MergedBMMX ntupleRawTree(treeChain);

	Long64_t nentries = treeChain->GetEntries();
    cout<<" Available total "<<nentries<<" \n";
    if (maxEvents >0 ) nentries = nentries > maxEvents ? maxEvents : nentries;
    cout<<" Processing total "<<nentries<<" \n";
    
    std::cout<<"Making the hists for the branches : \n";
    string branchName,tmp;

    float min,max,bwidth,iqr;
    TH1F* ahist;

    for (UInt_t i=0;i<branchList.size();i++) 
    {
        break;
        branchName=branchList[i];
        std::cout<<" \t doing :  "<< branchName<<"  "<<i<<" / "<<nentries<<"\n"; 
        tmp=branchName+" >> th1_tmp";
        treeChain->Draw( tmp.c_str() ,"","",nentries);
        auto hist_raw=(TH1*) gROOT->FindObject("th1_tmp");
        if(hist_raw == nullptr)
        {
            std::cout<<"\t"<<branchName<<" not found ! "<<"\n";
            continue;
        }
        hist_raw->GetQuantiles(nq,yq,xq);
        n= hist_raw->Integral() ;
        delete hist_raw;
        

        if(n<1) n=1;
        min = yq[0]> 0 ? yq[0]*0.80 : yq[0]*1.20;
        max = yq[1]> 0 ? yq[1]*1.20 : yq[1]*0.80;
        iqr = yq[3] - yq[2];
        if (iqr==0) iqr=1;
        bwidth = 2 * iqr / pow(n,1.0/3.0);
        nbins=(max-min)/bwidth + 1;
        if (nbins < NBINS) nbins=NBINS;
        // TODO : Need to float NBINS
        ahist = new TH1F(branchName.c_str(),branchName.c_str(),nbins,min,max);
        treeChain->Project(ahist->GetName(),branchName.c_str());
        
        if (ahist==nullptr)
         {
            std::cout<<"\t\t trouble !! \n";
            continue;
         }
         ahist->Write();     
    }

    TTree * outTree   = new TTree("bmmgCandidates","bmmgCandidates");
    

    UInt_t nCands(0);

    UInt_t photonMultiplicity(0);
    UInt_t   bmmgMultiplicity(0);
    
    //ROOT::Math::PtEtaPhiMVector aLV,bLV,cLV;
    TLorentzVector aLV,bLV,cLV;

    outTree->Branch("run",&(ntupleRawTree.b5_run));
    outTree->Branch("lumi",&(ntupleRawTree.b5_luminosityBlock));
    outTree->Branch("event",&(ntupleRawTree.b5_event));
    outTree->Branch("nCands",&(nCands));
    
    std::map<string, Int_t > candidateMap;
    double storageArray[1024];
    Int_t storageIdxFilled=0; 
    storageIdxFilled=setupDimuonBranches( outTree   , candidateMap ,storageArray, storageIdxFilled);
    storageIdxFilled=setupPhotonSCBranches( outTree , candidateMap , storageArray,storageIdxFilled);
    double dr(1e9);    

    candidateMap["mumu_dr"]   = storageIdxFilled ;
    outTree->Branch("mumu_dr",&(storageArray[storageIdxFilled]));

    candidateMap["dimuGamma_dr"]  = storageIdxFilled ;
    outTree->Branch("dimuGamma_dr",&(storageArray[storageIdxFilled]));
    storageIdxFilled++;

    candidateMap["bmmg_pt"]       = storageIdxFilled ;
    outTree->Branch("bmmg_pt",&(storageArray[storageIdxFilled]));
    storageIdxFilled++;

    candidateMap["bmmg_eta"]       = storageIdxFilled ;
    outTree->Branch("bmmg_eta",&(storageArray[storageIdxFilled]));
    storageIdxFilled++;

    candidateMap["bmmg_phi"]       = storageIdxFilled ;
    outTree->Branch("bmmg_phi",&(storageArray[storageIdxFilled]));
    storageIdxFilled++;

    candidateMap["bmmg_mass"]       = storageIdxFilled ;
    outTree->Branch("bmmg_mass",&(storageArray[storageIdxFilled]));
    storageIdxFilled++;

    candidateMap["bmmg_massErr"]       = storageIdxFilled ;
    outTree->Branch("bmmg_massErr",&(storageArray[storageIdxFilled]));
    storageIdxFilled++;

    outTree->Branch("bmmgMultiplicity",&(bmmgMultiplicity));
    outTree->Branch("photonMultiplicity",&(photonMultiplicity));
    
    bool isBMMGCandidate;
    if (maxEvents >0)
        nentries = nentries < maxEvents ? nentries : maxEvents;
   
    Long64_t EventCount=0;
    Long64_t nb = 0,nbytes=0 ;
    for (Long64_t jentry=0; jentry<nentries; jentry++)
    {
       
       if(jentry%1000 ==0 )
       {
            cout<<"Processing jentry : "<<jentry<<"\n";
       }
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(not ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ) continue;
       
       bmmgMultiplicity=0;
       for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    
           
           isBMMGCandidate=false;
          
           aLV.SetPtEtaPhiM(   ntupleRawTree.b5_mm_kin_pt[mumuIdx],\
                               ntupleRawTree.b5_mm_kin_eta[mumuIdx],\
                               ntupleRawTree.b5_mm_kin_phi[mumuIdx],\
                               ntupleRawTree.b5_mm_kin_mass[mumuIdx]) ;
            
          // aLV.SetCoordinates( ntupleRawTree.b5_mm_kin_pt[mumuIdx],\
                               ntupleRawTree.b5_mm_kin_eta[mumuIdx],\
                               ntupleRawTree.b5_mm_kin_phi[mumuIdx],\
                               ntupleRawTree.b5_mm_kin_mass[mumuIdx]) ;
           
           dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
                    ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
           candidateMap["mumu_dr"]=dr;
           photonMultiplicity=0;
            
           assigndimuonBMMGCandidates(candidateMap,ntupleRawTree,mumuIdx,storageArray);

           for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
           {
               
               auto et= ntupleRawTree.bG_scE[phoSCIdx]/cosh(ntupleRawTree.bG_scEta[phoSCIdx]);
               if( et < PhoPtMin) continue;
               dr=getDR(ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx],
                        ntupleRawTree.b5_mm_kin_eta[mumuIdx],ntupleRawTree.b5_mm_kin_phi[mumuIdx]);
               candidateMap["dimuGamma_dr"]=dr;
               if(dr > PhoDimuDrMax ) continue;
    
               bLV.SetPtEtaPhiM( et,
                                   ntupleRawTree.bG_scEta[phoSCIdx],
                                   ntupleRawTree.bG_scPhi[phoSCIdx],
                                   PHO_MASS
                                  );
               //bLV.SetCoordinates( et,\
                                   ntupleRawTree.bG_scEta[phoSCIdx],\
                                   ntupleRawTree.bG_scPhi[phoSCIdx],\
                                   PHO_MASS\
                                  );

               assignSCBMMGCandidates(candidateMap, ntupleRawTree,phoSCIdx,storageArray);

               candidateMap["bmmg_pt"]   = bLV.Pt();
               candidateMap["bmmg_eta"]  = bLV.Eta();
               candidateMap["bmmg_phi"]  = bLV.Phi();
               candidateMap["bmmg_mass"] = bLV.M();
                
               bmmgMultiplicity++;
               photonMultiplicity++;

               outTree->Fill();
                
           }

            
       }

       EventCount++;
        
    }
    
    outTree->Write();
    outputFile->Write();
    outputFile->Purge();
    outputFile->Close();

}

Double_t getDR( Double_t eta1, Double_t phi1,Double_t eta2 ,Double_t phi2) {

    Double_t de = getDETA(eta1,eta2);
    Double_t dp = getDPHI(phi1,phi2); 
    return sqrt(de*de + dp*dp);

}


Double_t getDPHI( Double_t phi1, Double_t phi2) {

  Double_t dphi = phi1 - phi2;
  
  while( dphi > PI)
        dphi-= TWO_PI; 
  while( dphi < -PI)
        dphi += TWO_PI; 
  //std::cout<<"phi1  : "<<phi1<<" , phi2 : "<<phi2<<" dphi : "<<dphi<<"\n";
  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 "<< dphi << endl;
  }
  
  return TMath::Abs(dphi);
  //return dphi;
}



Double_t getDETA(Double_t eta1, Double_t eta2){
  return TMath::Abs(eta1 - eta2);
}
Int_t setupDimuonBranches( TTree* outTree ,std::map<string, Int_t > &candidateMap , double * storageArray , Int_t offset=0 )
{
    #include  "setupDimu.cc"
    
}

void assigndimuonBMMGCandidates( std::map<string,Int_t > &candidateMap, MergedBMMX &ntupleRawTree, int candIdx,double * storageArray)
{
   #include "assignDimu.cc"
}


Int_t setupPhotonSCBranches( TTree* outTree ,std::map<string, Int_t > &candidateMap , double * storageArray, Int_t offset=0 )
{

    #include  "setupSc.cc"

}

void assignSCBMMGCandidates( std::map<string, Int_t > &candidateMap, MergedBMMX &ntupleRawTree, int candIdx,double * storageArray)
{


   #include "assignSc.cc"
}
