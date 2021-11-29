/*
 *          USAGE
 *
 * root -b -q 'bmmX_distributionStudy.cc("YOURCONFIG.CFG")'
 *
 * */

#include "MergedBMMX2018Data.h"

typedef MergedBMMX2018Data MergedBMMX ;

#include "main.h"

#include "selections.h"
#include "muonSelection.h"
#include "photonSelection.h"

#define PHO_MASS 0.0
#define  Mu_MASS 0.0


#define NSTORAGE_ARRAY_MAX 200000
#define NDIMU_MAX 20
#define NBMMG_MAX 40
#define NSC_MAX   50
#define NMUONS_MAX  20

int  main()
{
    Long64_t maxEvents(-1000);
    std::vector<string> InFileList;
    string ofileName("output.root");
    string   treeName("mergedTree");
    string prefix("");
    
    std::map<string, TH1D * > histMap;
    
    TChain *treeChain = new TChain(treeName.c_str());

    InFileList.push_back("bmmXMerged2018data_0_0.root");
    for(auto i=0;i<InFileList.size();i++)
    {
        treeChain->Add(InFileList[i].c_str());
    }
    
    TFile * outputFile = new  TFile((prefix+ofileName).c_str(),"recreate");    
    

    MergedBMMX ntupleRawTree(treeChain);
    
    TTree * outTree   = new TTree("AnalysisTree","Reduced branches from the Merged [bmm5+bmmX] trees.");
    
    bool isTriggerd;
    outTree->Branch("run",&(ntupleRawTree.b5_run));
    outTree->Branch("lumi",&(ntupleRawTree.b5_luminosityBlock));
    outTree->Branch("event",&(ntupleRawTree.b5_event));
    outTree->Branch("isTriggerd",&(isTriggerd));
    
    
    // Defenition of the storage of type Double_t
    std::map<string, Int_t > candidateMapDouble;
    Double_t storageArrayDouble[NSTORAGE_ARRAY_MAX];
    std::cout<<"Allocated "<<sizeof(Int_t)*NSTORAGE_ARRAY_MAX/1024<<" kB of storage for "<<NSTORAGE_ARRAY_MAX<<" Doubles \n";
    Int_t storageIdxFilledDouble=0; 
    
    outTree->Branch("nMuons",&(ntupleRawTree.b5_nMuon));
    outTree->Branch("nSCPhotons",&(ntupleRawTree.bG_nSC));
    Int_t nDiMuCandidates=0;
    outTree->Branch("nDiMuCandidates",&(nDiMuCandidates));
    Int_t nBMMGCandidates=0;
    outTree->Branch("nBMMGCandidates",&(nBMMGCandidates));
    Int_t nBMMGCandidatesPerDimu(0);
    Int_t photonSelectionCheck[NSC_MAX];

    // Defenition of the storage of type Int_t
    std::map<string, Int_t > candidateMapInt;
    Int_t storageArrayInt[NSTORAGE_ARRAY_MAX]; std::cout<<"Allocated "<<sizeof(Int_t)*NSTORAGE_ARRAY_MAX/1024<<" kB of storage for "<<NSTORAGE_ARRAY_MAX<<" Ints\n";
    Int_t storageIdxFilledInt=0;

    candidateMapInt["bmmg_dimuon_idx"]       = storageIdxFilledInt ;
    outTree->Branch("bmmg_dimuon_idx",&(storageArrayInt[storageIdxFilledInt]),"dimuon_idx[nBMMGCandidates]/I" )   ;   storageIdxFilledInt+=NBMMG_MAX;
    
    candidateMapInt["bmmg_photonSC_idx"]       = storageIdxFilledInt ;
    outTree->Branch("bmmg_photonSC_idx",&(storageArrayInt[storageIdxFilledInt]),"photonSC_idx[nBMMGCandidates]/I" )   ;   storageIdxFilledInt+=NBMMG_MAX;

    candidateMapInt["nBMMGCandidatesPerDimu"]   = storageIdxFilledInt ;
    outTree->Branch("nBMMGCandidatesPerDimu",&(storageArrayInt[storageIdxFilledInt]),"nBMMGCandidatesPerDimu[nDiMuCandidates]/I");storageIdxFilledInt+=NDIMU_MAX;
    
    candidateMapDouble["mumu_dr"]   = storageIdxFilledDouble ;
    outTree->Branch("mumu_dr",&(storageArrayDouble[storageIdxFilledDouble]),"mumu_dr[nDiMuCandidates]/D");storageIdxFilledDouble+=NDIMU_MAX;
    
    
    candidateMapDouble["dimuGamma_dr"]  = storageIdxFilledDouble                ;
    outTree->Branch( "dimuGamma_dr",&(storageArrayDouble[storageIdxFilledDouble]),"dimuGamma_dr[nBMMGCandidates]/D" )   ;   storageIdxFilledDouble+=NBMMG_MAX;

    candidateMapDouble["bmmg_pt"]       = storageIdxFilledDouble ;
    outTree->Branch( "bmmg_pt",&(storageArrayDouble[storageIdxFilledDouble]),"bmmg_pt[nBMMGCandidates]/D" )   ;   storageIdxFilledDouble+=NBMMG_MAX;

    candidateMapDouble["bmmg_eta"]       = storageIdxFilledDouble ;
    outTree->Branch( "bmmg_eta",&(storageArrayDouble[storageIdxFilledDouble]),"bmmg_eta[nBMMGCandidates]/D" )   ;   storageIdxFilledDouble+=NBMMG_MAX;

    candidateMapDouble["bmmg_phi"]       = storageIdxFilledDouble ;
    outTree->Branch( "bmmg_phi",&(storageArrayDouble[storageIdxFilledDouble]),"bmmg_phi[nBMMGCandidates]/D" )   ;   storageIdxFilledDouble+=NBMMG_MAX;

    candidateMapDouble["bmmg_mass"]       = storageIdxFilledDouble ;
    outTree->Branch( "bmmg_mass",&(storageArrayDouble[storageIdxFilledDouble]),"bmmg_mass[nBMMGCandidates]/D" )   ;   storageIdxFilledDouble+=NBMMG_MAX;

    candidateMapDouble["bmmg_massErr"]       = storageIdxFilledDouble ;
    outTree->Branch( "bmmg_massErr",&(storageArrayDouble[storageIdxFilledDouble]),"bmmg_massErr[nBMMGCandidates]/D" )   ;   storageIdxFilledDouble+=NBMMG_MAX;
    
    
    storageIdxFilledDouble  =  setupDimuonBranches  ( outTree , candidateMapDouble , storageArrayDouble,storageIdxFilledDouble , NDIMU_MAX );
    storageIdxFilledDouble  =  setupMuonBranches    ( outTree , candidateMapDouble , storageArrayDouble,storageIdxFilledDouble , NMUONS_MAX);
    storageIdxFilledDouble  =  setupPhotonSCBranches( outTree , candidateMapDouble , storageArrayDouble,storageIdxFilledDouble , NSC_MAX   );
    
	
    TLorentzVector diMuLV,photonLV,bmmgLV;
    Double_t dr;

    Long64_t nentries = treeChain->GetEntries();
    cout<<" Available total "<<nentries<<" \n";
    if (maxEvents >0 ) nentries = nentries > maxEvents ? maxEvents : nentries;
    cout<<" Processing total "<<nentries<<" \n";
   
    Long64_t EventCount=0;
    Long64_t nb = 0,nbytes=0 ;
    for (Long64_t jentry=0; jentry<nentries; jentry++)
    {   

       nDiMuCandidates=0;
       isTriggerd=false;
       
       if(jentry%500 == 0 )
       {
            cout<<"Processing jentry : "<<jentry<<"\n";
       }
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);

       if (ientry_evt < 0) break;

       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       // Trigger Selection
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ) isTriggerd=true;

       if(not isTriggerd) continue;

        for(Int_t i=0;i<NSTORAGE_ARRAY_MAX;i++)
            storageArrayDouble[i]=0;
       for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
            photonSelectionCheck[phoSCIdx]=-1;

       // BMMG Selection

       int rslt=0;
       float dimuEta,dimuPhi;
       

       nBMMGCandidates=0;
       for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    
           // Muon Selection
		   rslt=doMuonSelection( ntupleRawTree , ntupleRawTree.b5_mm_mu1_index[mumuIdx], true);
           if(rslt > 0) continue;
		   rslt=doMuonSelection( ntupleRawTree , ntupleRawTree.b5_mm_mu2_index[mumuIdx], false);
           if(rslt > 0) continue;


           // Dimuon Selection
           diMuLV.SetPtEtaPhiM(     ntupleRawTree.b5_mm_kin_pt[mumuIdx],   \
                                    ntupleRawTree.b5_mm_kin_eta[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_phi[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_mass[mumuIdx]  );
           
           dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
                    ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
           if (dr < minMuMuDr ) continue;
           
           storageArrayDouble[nDiMuCandidates + candidateMapDouble["mumu_dr"] ]   = dr ;
            
            /*
                    VERTEX SELECTION STUFF
            */


           dimuEta=ntupleRawTree.b5_mm_kin_eta[mumuIdx]  ;
           dimuPhi=ntupleRawTree.b5_mm_kin_phi[mumuIdx]  ;
           nBMMGCandidatesPerDimu=0;

           for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
           {
               // photon selection
               if(photonSelectionCheck[phoSCIdx] < 0)
               {
                    photonSelectionCheck[phoSCIdx]=doPhotonSelection( ntupleRawTree , phoSCIdx );
               }

               if( photonSelectionCheck[phoSCIdx] > 0) continue;

               dr=getDR(ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx],dimuEta,dimuPhi);

               if(dr > PhoDimuDrMax ) continue;

               auto et= ntupleRawTree.bG_scE[phoSCIdx]/cosh(ntupleRawTree.bG_scEta[phoSCIdx]);

               photonLV.SetPtEtaPhiM(  et       ,     ntupleRawTree.bG_scEta[phoSCIdx],
                                                       ntupleRawTree.bG_scPhi[phoSCIdx], PHO_MASS );
               bmmgLV = diMuLV + photonLV;
               
              // assignBMMGCandidates(candidateMapDouble,ntupleRawTree,dimuIDx,storageArrayDouble);
               storageArrayDouble[  nBMMGCandidates + candidateMapDouble["dimuGamma_dr"]          ]   = dr           ;
               storageArrayDouble[  nBMMGCandidates + candidateMapDouble["bmmg_pt"]               ]   = bmmgLV.Pt()  ;
               storageArrayDouble[  nBMMGCandidates + candidateMapDouble["bmmg_eta"]              ]   = bmmgLV.Eta() ;
               storageArrayDouble[  nBMMGCandidates + candidateMapDouble["bmmg_phi"]              ]   = bmmgLV.Phi() ;
               storageArrayDouble[  nBMMGCandidates + candidateMapDouble["bmmg_mass"]             ]   = bmmgLV.M()   ;
               storageArrayDouble[  nBMMGCandidates + candidateMapDouble["bmmg_massErr"]          ]   = bmmgLV.M()*-1; // TODO
               
               storageArrayInt[nBMMGCandidates + candidateMapInt["bmmg_dimuon_idx"]   ]        = mumuIdx       ;
               storageArrayInt[nBMMGCandidates + candidateMapInt["bmmg_photonSC_idx"] ]        = phoSCIdx      ;
               
               nBMMGCandidatesPerDimu++;
               nBMMGCandidates++;

               if(nBMMGCandidates >= NBMMG_MAX)
               {
                    std::cout<<" BMMG count per event above NBMMG_MAX , Aborting !! \n";
                    exit(9);
               }          
           }
           
           
           storageArrayInt[ nDiMuCandidates + candidateMapInt["nBMMGCandidatesPerDimu"]  ]   = nBMMGCandidatesPerDimu ;
           assignDimuonBMMGCandidates(candidateMapDouble,ntupleRawTree,mumuIdx,storageArrayDouble,nDiMuCandidates);
           
           nDiMuCandidates++;

           if(nDiMuCandidates >= NDIMU_MAX)
           {
                std::cout<<" Dimuon count per event above NDIMU_MAX , Aborting !! \n";
                exit(9);
           }
            
       }
       
       assignSC(candidateMapDouble, ntupleRawTree,storageArrayDouble);
       assignMuons(candidateMapDouble, ntupleRawTree,storageArrayDouble);

       EventCount++;
       outTree->Fill();
    }
    
    outTree->Write();
    outputFile->Write();
    outputFile->Purge();
    outputFile->Close();

}
