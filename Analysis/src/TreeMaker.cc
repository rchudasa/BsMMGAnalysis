//#include "TreeMaker.h"

void TreeMaker::DataDimuonTreeMaker()
{

    AddDimuonVertexHistos("data");
    AddDimuonTree("data");
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nCands=0,rslt;
    
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();

    Double_t dr=-1.0;
    Double_t drMin;

    Bool_t hasACand=false;
    
    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {  
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       EventCount++;
       hasACand=false;
       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
       
       for(UInt_t mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    

                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] < 0 or ntupleRawTree.b5_mm_mu2_index[mumuIdx] <0 ) 
	            {
                    continue;
                }
                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    continue;
	            }
                if(ntupleRawTree.b5_mm_mu2_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    continue;
	            }
           
           // Muon Selection
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu1_index[mumuIdx], true);
           if(rslt > 0) continue;
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu2_index[mumuIdx], false);
           if(rslt > 0) continue;
           if( ntupleRawTree.b5_mm_m1iso[mumuIdx]  < minMuIsolation ) 
             continue;
           if( ntupleRawTree.b5_mm_m2iso[mumuIdx]  < minMuIsolation ) 
             continue;

           fillDimuon(mumuIdx,"data");
           nCands++;
           hasACand=true;
      }
      if(hasACand) EventCountWithCand++;

    }

    std::cout<<" Number of Evnets processed        : "<<EventCount<<"\n";
    std::cout<<" Number of Evnets with candidates  : "<<EventCountWithCand<<"\n";
    std::cout<<" Number of candidates              : "<<nCands<<"\n";

}

void TreeMaker::fillDimuon(Int_t mumuIdx,TString tag)
{
           fillDimuonVertexVars(mumuIdx);
           fillDimuonVertexHists(tag);
           
           treeStore[tag]->Fill();
}

#ifdef __MCANALYSIS__
void TreeMaker::genMatchedDimuons()
{
    AddDimuonVertexHistos("genMatchedDimuons");
    AddDimuonTree("genMatchedDimuons");
    
    Double_t dr;
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nCands=0;
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    bool goodRunLumi = false;

    Double_t drMin;
    Int_t scMatchIdx;
    Bool_t foundmatch=false;
    Int_t mumMatchIdx(-1),mupMatchIdx(-1),phoMatchIdx(-1);
    Int_t mumFoundIdx(-1),mupFoundIdx(-1),phoFoundIdx(-1);
    Int_t motherIdx(-1);
    Long64_t mumMatchCount=0;
    Long64_t mupMatchCount=0;
    Long64_t scMatchCount=0;
    Long64_t phoMatchCount=0;
    Bool_t hasCand=false;

    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {
        hasCand=false;
        
       eventGenMultiplicity=0;
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
       EventCount++;

//==================================================================================
    std::cout<<__LINE__<<"\n";

       for(Int_t i=0; i< ntupleRawTree.b5_nGenPart; i++)
       {
            if( abs(ntupleRawTree.b5_GenPart_pdgId[i]) == 531 )
            {
                   mumMatchIdx=-1; mupMatchIdx=-1; phoMatchIdx=-1;
                   mumFoundIdx=0;  mupFoundIdx=0;  phoFoundIdx=0;
                   motherIdx=i;
                   for(Int_t j=0; j< ntupleRawTree.b5_nGenPart; j++)
                    {
                          if( ntupleRawTree.b5_GenPart_genPartIdxMother[j] == motherIdx)
                          {
                                if(ntupleRawTree.b5_GenPart_pdgId[j]==13) 
                                {
                                    mumFoundIdx++;
                                    mumMatchIdx=j;
                                }
                                
                                if(ntupleRawTree.b5_GenPart_pdgId[j]==-13) 
                                {
                                    mupFoundIdx++;
                                    mupMatchIdx=j;
                                }

                                if(ntupleRawTree.b5_GenPart_pdgId[j]==22) 
                                {
                                    phoFoundIdx++;
                                    phoMatchIdx=j;
                                }
                          }
                    }
                    if(mupFoundIdx==1 and mumFoundIdx==1)
                    {
                        break;
                    }
                    else
                    {
                        mumMatchIdx=-1; mupMatchIdx=-1; phoMatchIdx=-1;
                    }

            }
       }
    std::cout<<__LINE__<<"\n";
       if(mumFoundIdx != 1 or mupFoundIdx !=1 )
       {
            std::cout<<" No Gen Bs Decay found for jentry =  "<<jentry<<" !!\n ";
            continue;
       }
   
       TLorentzVector mumLV,mupLV,bmmgLV,photonLV;
    std::cout<<__LINE__<<"\n";

       if( phoFoundIdx == 0 )
       {
            bmmgLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[motherIdx], ntupleRawTree.b5_GenPart_eta[motherIdx],
                                ntupleRawTree.b5_GenPart_phi[motherIdx], ntupleRawTree.b5_GenPart_mass[motherIdx] )    ;
            mumLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[mumMatchIdx], ntupleRawTree.b5_GenPart_eta[mumMatchIdx],
                                 ntupleRawTree.b5_GenPart_phi[mumMatchIdx], ntupleRawTree.b5_GenPart_mass[mumMatchIdx] ) ;    
            mupLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[mupMatchIdx], ntupleRawTree.b5_GenPart_eta[mupMatchIdx],
                                ntupleRawTree.b5_GenPart_phi[mupMatchIdx], ntupleRawTree.b5_GenPart_mass[mupMatchIdx] )   ; 
            photonLV=bmmgLV - ( mumLV  +  mupLV );
       //     std::cout<<"Recalculating Pho P4 pt  = "<<photonLV.Pt()<<" , eta = "<<photonLV.Eta()<<"\n" ;
       }
       else if( phoFoundIdx==1)
       {
            photonLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[phoMatchIdx], ntupleRawTree.b5_GenPart_eta[phoMatchIdx],
                                ntupleRawTree.b5_GenPart_phi[phoMatchIdx], ntupleRawTree.b5_GenPart_mass[phoMatchIdx] )   ; 
       }

       Int_t mumRecoMatchIdx(-1),mupRecoMatchIdx(-1),phoRecoMatchIdx(-1);
       Float_t mumRecoMatchDr(drMaxForGenMatchMu),mupRecoMatchDr(drMaxForGenMatchMu),phoRecoMatchDr(drMaxForGenMatchSC);
       for(int i=0;i<ntupleRawTree.b5_nMuon ; i++)
       {
            if(ntupleRawTree.b5_Muon_pdgId[i] == 13 )
            {
                dr= getDR(ntupleRawTree.b5_Muon_eta[i], ntupleRawTree.b5_Muon_phi[i],
                ntupleRawTree.b5_GenPart_eta[mumMatchIdx],ntupleRawTree.b5_GenPart_phi[mumMatchIdx] );
                if( dr < mumRecoMatchDr) 
                {
                    mumRecoMatchDr=dr;
                    mumRecoMatchIdx=i;
                }
            }
            if(ntupleRawTree.b5_Muon_pdgId[i] == -13 )
            {
                dr= getDR(ntupleRawTree.b5_Muon_eta[i], ntupleRawTree.b5_Muon_phi[i], 
                          ntupleRawTree.b5_GenPart_eta[mupMatchIdx] ,ntupleRawTree.b5_GenPart_phi[mupMatchIdx] );
                if( dr < mupRecoMatchDr) 
                {
                    mupRecoMatchDr=dr;
                    mupRecoMatchIdx=i;
                }
            }
       }
        
       for(int i=0;i<ntupleRawTree.bG_nSC ; i++)
       {
                dr= getDR(ntupleRawTree.bG_scEta[i], ntupleRawTree.bG_scPhi[i], 
                 photonLV.Eta() ,photonLV.Phi() );
                if( dr < phoRecoMatchDr) 
                {
                    phoRecoMatchDr=dr;
                    phoRecoMatchIdx=i;
                }
       }
       
       if(mumRecoMatchDr < drMaxForGenMatchMu) { 
        mumMatchCount++;
	  }
      else
      {
        mumRecoMatchIdx=-1;
      }
       if(mupRecoMatchDr < drMaxForGenMatchMu) { 
          mupMatchCount++;
	  }
      else
      {
        mupRecoMatchIdx=-1;
      }

      if(phoRecoMatchDr < drMaxForGenMatchSC ) { 
        scMatchCount++;
        //if(doPhotonSelection(phoRecoMatchIdx) ==0)
        //{   
        //    th1fStore["ProcessingSummary"]->Fill("nGenMatchedPho",1);
        //    phoMatchCount++;
        //}
		//fill_photonHists(phoRecoMatchIdx);
     }
      else
      {
         phoRecoMatchIdx=-1;
      }
     

    std::cout<<__LINE__<<"\n";

       if(mumRecoMatchDr < drMaxForGenMatchMu and mupRecoMatchDr < drMaxForGenMatchMu)
       {
       }
       if(mumRecoMatchDr < drMaxForGenMatchMu and mupRecoMatchDr <  drMaxForGenMatchMu and phoRecoMatchDr < drMaxForGenMatchSC)
       {
           // fullEventMatchesFound++;
       }
       

       if(mumRecoMatchIdx < 0 )
       {
                continue;
       }
       else if ( mupRecoMatchIdx <0 )
       {
                continue;
       }
       else 
       {
    std::cout<<__LINE__<<"\n";
            for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
            {   

                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] < 0 or ntupleRawTree.b5_mm_mu2_index[mumuIdx] <0 ) 
	            {
                    continue;
                }
                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    continue;
	            }
                if(ntupleRawTree.b5_mm_mu2_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    continue;
	            }
    std::cout<<__LINE__<<"\n";
               
            if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] != mumRecoMatchIdx and ntupleRawTree.b5_mm_mu2_index[mumuIdx] != mumRecoMatchIdx ) {
                  continue;
                }
           if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] != mupRecoMatchIdx and ntupleRawTree.b5_mm_mu2_index[mumuIdx] != mupRecoMatchIdx ) { 
                  continue;  
                }
    std::cout<<__LINE__<<"\n";
           
           if( ntupleRawTree.b5_mm_m1iso[mumuIdx]  < minMuIsolation ) 
             continue;
           if( ntupleRawTree.b5_mm_m2iso[mumuIdx]  < minMuIsolation ) 
                continue;
    std::cout<<__LINE__<<"\n";
              fillDimuon(mumuIdx,"genMatchedDimuons");
              hasCand=true;
              nCands++;

            }
        }
        if(hasCand)
            EventCountWithCand++;

//==================================================================================
       
    }

    std::cout<<" Number of Evnets processed        : "<<EventCount<<"\n";
    std::cout<<" Number of Evnets with candidates  : "<<EventCountWithCand<<"\n";
    std::cout<<" Number of candidates              : "<<nCands<<"\n";
}
#endif



void TreeMaker::AddDimuonTree(TString tag)
{
    outTree= new TTree(tag+"_tree",tag);
    treeStore[tag]=outTree;
    TString tg;
    tg="dimuon_lXY";storageFloat[tg]=0;  outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_bdt";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_doca";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_docatrk";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_iso";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kal_lxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kal_mass";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kal_slxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kal_vtx_prob";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_alpha";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_alphaBS";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_alphaBSErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_alphaErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_eta";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_l3d";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_lxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_mass";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_massErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_mu1eta";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_mu1phi";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_mu1pt";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_mu2eta";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_mu2phi";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_mu2pt";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_phi";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pt";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pv2ip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pv2ipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pv2lip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pv2lipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pv2lipSig";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pv_z";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pv_zErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pvip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pvipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pvlip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pvlipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_pvlipSig";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_sl3d";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_slxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_spv2ip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_spvip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_tau";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_taue";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_tauxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_tauxye";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_chi2dof";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_prob";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_x";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_xErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_y";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_yErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_z";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_vtx_zErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_alpha";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_alphaBS";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_alphaBSErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_alphaErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_eta";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_l3d";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_lxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_mass";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_massErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_phi";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pt";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pv2ip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pv2ipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pv2lip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pv2lipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pv2lipSig";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pv_z";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pv_zErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pvip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pvipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pvlip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pvlipErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_pvlipSig";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_sl3d";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_slxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_spv2ip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_spvip";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_tau";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_taue";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_tauxy";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_tauxye";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_chi2dof";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_prob";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_x";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_xErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_y";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_yErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_z";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_vtx_zErr";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_m1iso";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_m2iso";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mass";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu1_eta";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu1_phi";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu1_pt";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu2_eta";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu2_phi";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu2_pt";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mva";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_otherVtxMaxProb";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_otherVtxMaxProb1";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_otherVtxMaxProb2";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_closetrk";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_closetrks1";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_closetrks2";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_closetrks3";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kal_valid";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kin_valid";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_kinpc_valid";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu1_index";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu1_pdgId";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu2_index";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_mu2_pdgId";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_nBMTrks";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_nDisTrks";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
    tg="mm_nTrks";  storageFloat[tg]=0; outTree->Branch(tg,&storageFloat[tg]);
}


void TreeMaker::fillDimuonVertexVars( Int_t mumuIdx)
{
    storageFloat["mm_bdt"]=ntupleRawTree.b5_mm_bdt[mumuIdx] ;
    storageFloat["mm_doca"]=ntupleRawTree.b5_mm_doca[mumuIdx] ;
    storageFloat["mm_docatrk"]=ntupleRawTree.b5_mm_docatrk[mumuIdx] ;
    storageFloat["mm_iso"]=ntupleRawTree.b5_mm_iso[mumuIdx] ;
    storageFloat["mm_kal_lxy"]=ntupleRawTree.b5_mm_kal_lxy[mumuIdx] ;
    storageFloat["mm_kal_mass"]=ntupleRawTree.b5_mm_kal_mass[mumuIdx] ;
    storageFloat["mm_kal_slxy"]=ntupleRawTree.b5_mm_kal_slxy[mumuIdx] ;
    storageFloat["mm_kal_vtx_prob"]=ntupleRawTree.b5_mm_kal_vtx_prob[mumuIdx] ;
    storageFloat["mm_kin_alpha"]=ntupleRawTree.b5_mm_kin_alpha[mumuIdx] ;
    storageFloat["mm_kin_alphaBS"]=ntupleRawTree.b5_mm_kin_alphaBS[mumuIdx] ;
    storageFloat["mm_kin_alphaBSErr"]=ntupleRawTree.b5_mm_kin_alphaBSErr[mumuIdx] ;
    storageFloat["mm_kin_alphaErr"]=ntupleRawTree.b5_mm_kin_alphaErr[mumuIdx] ;
    storageFloat["mm_kin_eta"]=ntupleRawTree.b5_mm_kin_eta[mumuIdx] ;
    storageFloat["mm_kin_l3d"]=ntupleRawTree.b5_mm_kin_l3d[mumuIdx] ;
    storageFloat["mm_kin_lxy"]=ntupleRawTree.b5_mm_kin_lxy[mumuIdx] ;
    storageFloat["mm_kin_mass"]=ntupleRawTree.b5_mm_kin_mass[mumuIdx] ;
    storageFloat["mm_kin_massErr"]=ntupleRawTree.b5_mm_kin_massErr[mumuIdx] ;
    storageFloat["mm_kin_mu1eta"]=ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx] ;
    storageFloat["mm_kin_mu1phi"]=ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ;
    storageFloat["mm_kin_mu1pt"]=ntupleRawTree.b5_mm_kin_mu1pt[mumuIdx] ;
    storageFloat["mm_kin_mu2eta"]=ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx] ;
    storageFloat["mm_kin_mu2phi"]=ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] ;
    storageFloat["mm_kin_mu2pt"]=ntupleRawTree.b5_mm_kin_mu2pt[mumuIdx] ;
    storageFloat["mm_kin_phi"]=ntupleRawTree.b5_mm_kin_phi[mumuIdx] ;
    storageFloat["mm_kin_pt"]=ntupleRawTree.b5_mm_kin_pt[mumuIdx] ;
    storageFloat["mm_kin_pv2ip"]=ntupleRawTree.b5_mm_kin_pv2ip[mumuIdx] ;
    storageFloat["mm_kin_pv2ipErr"]=ntupleRawTree.b5_mm_kin_pv2ipErr[mumuIdx] ;
    storageFloat["mm_kin_pv2lip"]=ntupleRawTree.b5_mm_kin_pv2lip[mumuIdx] ;
    storageFloat["mm_kin_pv2lipErr"]=ntupleRawTree.b5_mm_kin_pv2lipErr[mumuIdx] ;
    storageFloat["mm_kin_pv2lipSig"]=ntupleRawTree.b5_mm_kin_pv2lipSig[mumuIdx] ;
    storageFloat["mm_kin_pv_z"]=ntupleRawTree.b5_mm_kin_pv_z[mumuIdx] ;
    storageFloat["mm_kin_pv_zErr"]=ntupleRawTree.b5_mm_kin_pv_zErr[mumuIdx] ;
    storageFloat["mm_kin_pvip"]=ntupleRawTree.b5_mm_kin_pvip[mumuIdx] ;
    storageFloat["mm_kin_pvipErr"]=ntupleRawTree.b5_mm_kin_pvipErr[mumuIdx] ;
    storageFloat["mm_kin_pvlip"]=ntupleRawTree.b5_mm_kin_pvlip[mumuIdx] ;
    storageFloat["mm_kin_pvlipErr"]=ntupleRawTree.b5_mm_kin_pvlipErr[mumuIdx] ;
    storageFloat["mm_kin_pvlipSig"]=ntupleRawTree.b5_mm_kin_pvlipSig[mumuIdx] ;
    storageFloat["mm_kin_sl3d"]=ntupleRawTree.b5_mm_kin_sl3d[mumuIdx] ;
    storageFloat["mm_kin_slxy"]=ntupleRawTree.b5_mm_kin_slxy[mumuIdx] ;
    storageFloat["mm_kin_spv2ip"]=ntupleRawTree.b5_mm_kin_spv2ip[mumuIdx] ;
    storageFloat["mm_kin_spvip"]=ntupleRawTree.b5_mm_kin_spvip[mumuIdx] ;
    storageFloat["mm_kin_tau"]=ntupleRawTree.b5_mm_kin_tau[mumuIdx] ;
    storageFloat["mm_kin_taue"]=ntupleRawTree.b5_mm_kin_taue[mumuIdx] ;
    storageFloat["mm_kin_tauxy"]=ntupleRawTree.b5_mm_kin_tauxy[mumuIdx] ;
    storageFloat["mm_kin_tauxye"]=ntupleRawTree.b5_mm_kin_tauxye[mumuIdx] ;
    storageFloat["mm_kin_vtx_chi2dof"]=ntupleRawTree.b5_mm_kin_vtx_chi2dof[mumuIdx] ;
    storageFloat["mm_kin_vtx_prob"]=ntupleRawTree.b5_mm_kin_vtx_prob[mumuIdx] ;
    storageFloat["mm_kin_vtx_x"]=ntupleRawTree.b5_mm_kin_vtx_x[mumuIdx] ;
    storageFloat["mm_kin_vtx_xErr"]=ntupleRawTree.b5_mm_kin_vtx_xErr[mumuIdx] ;
    storageFloat["mm_kin_vtx_y"]=ntupleRawTree.b5_mm_kin_vtx_y[mumuIdx] ;
    storageFloat["mm_kin_vtx_yErr"]=ntupleRawTree.b5_mm_kin_vtx_yErr[mumuIdx] ;
    storageFloat["mm_kin_vtx_z"]=ntupleRawTree.b5_mm_kin_vtx_z[mumuIdx] ;
    storageFloat["mm_kin_vtx_zErr"]=ntupleRawTree.b5_mm_kin_vtx_zErr[mumuIdx] ;
    storageFloat["mm_kinpc_alpha"]=ntupleRawTree.b5_mm_kinpc_alpha[mumuIdx] ;
    storageFloat["mm_kinpc_alphaBS"]=ntupleRawTree.b5_mm_kinpc_alphaBS[mumuIdx] ;
    storageFloat["mm_kinpc_alphaBSErr"]=ntupleRawTree.b5_mm_kinpc_alphaBSErr[mumuIdx] ;
    storageFloat["mm_kinpc_alphaErr"]=ntupleRawTree.b5_mm_kinpc_alphaErr[mumuIdx] ;
    storageFloat["mm_kinpc_eta"]=ntupleRawTree.b5_mm_kinpc_eta[mumuIdx] ;
    storageFloat["mm_kinpc_l3d"]=ntupleRawTree.b5_mm_kinpc_l3d[mumuIdx] ;
    storageFloat["mm_kinpc_lxy"]=ntupleRawTree.b5_mm_kinpc_lxy[mumuIdx] ;
    storageFloat["mm_kinpc_mass"]=ntupleRawTree.b5_mm_kinpc_mass[mumuIdx] ;
    storageFloat["mm_kinpc_massErr"]=ntupleRawTree.b5_mm_kinpc_massErr[mumuIdx] ;
    storageFloat["mm_kinpc_phi"]=ntupleRawTree.b5_mm_kinpc_phi[mumuIdx] ;
    storageFloat["mm_kinpc_pt"]=ntupleRawTree.b5_mm_kinpc_pt[mumuIdx] ;
    storageFloat["mm_kinpc_pv2ip"]=ntupleRawTree.b5_mm_kinpc_pv2ip[mumuIdx] ;
    storageFloat["mm_kinpc_pv2ipErr"]=ntupleRawTree.b5_mm_kinpc_pv2ipErr[mumuIdx] ;
    storageFloat["mm_kinpc_pv2lip"]=ntupleRawTree.b5_mm_kinpc_pv2lip[mumuIdx] ;
    storageFloat["mm_kinpc_pv2lipErr"]=ntupleRawTree.b5_mm_kinpc_pv2lipErr[mumuIdx] ;
    storageFloat["mm_kinpc_pv2lipSig"]=ntupleRawTree.b5_mm_kinpc_pv2lipSig[mumuIdx] ;
    storageFloat["mm_kinpc_pv_z"]=ntupleRawTree.b5_mm_kinpc_pv_z[mumuIdx] ;
    storageFloat["mm_kinpc_pv_zErr"]=ntupleRawTree.b5_mm_kinpc_pv_zErr[mumuIdx] ;
    storageFloat["mm_kinpc_pvip"]=ntupleRawTree.b5_mm_kinpc_pvip[mumuIdx] ;
    storageFloat["mm_kinpc_pvipErr"]=ntupleRawTree.b5_mm_kinpc_pvipErr[mumuIdx] ;
    storageFloat["mm_kinpc_pvlip"]=ntupleRawTree.b5_mm_kinpc_pvlip[mumuIdx] ;
    storageFloat["mm_kinpc_pvlipErr"]=ntupleRawTree.b5_mm_kinpc_pvlipErr[mumuIdx] ;
    storageFloat["mm_kinpc_pvlipSig"]=ntupleRawTree.b5_mm_kinpc_pvlipSig[mumuIdx] ;
    storageFloat["mm_kinpc_sl3d"]=ntupleRawTree.b5_mm_kinpc_sl3d[mumuIdx] ;
    storageFloat["mm_kinpc_slxy"]=ntupleRawTree.b5_mm_kinpc_slxy[mumuIdx] ;
    storageFloat["mm_kinpc_spv2ip"]=ntupleRawTree.b5_mm_kinpc_spv2ip[mumuIdx] ;
    storageFloat["mm_kinpc_spvip"]=ntupleRawTree.b5_mm_kinpc_spvip[mumuIdx] ;
    storageFloat["mm_kinpc_tau"]=ntupleRawTree.b5_mm_kinpc_tau[mumuIdx] ;
    storageFloat["mm_kinpc_taue"]=ntupleRawTree.b5_mm_kinpc_taue[mumuIdx] ;
    storageFloat["mm_kinpc_tauxy"]=ntupleRawTree.b5_mm_kinpc_tauxy[mumuIdx] ;
    storageFloat["mm_kinpc_tauxye"]=ntupleRawTree.b5_mm_kinpc_tauxye[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_chi2dof"]=ntupleRawTree.b5_mm_kinpc_vtx_chi2dof[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_prob"]=ntupleRawTree.b5_mm_kinpc_vtx_prob[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_x"]=ntupleRawTree.b5_mm_kinpc_vtx_x[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_xErr"]=ntupleRawTree.b5_mm_kinpc_vtx_xErr[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_y"]=ntupleRawTree.b5_mm_kinpc_vtx_y[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_yErr"]=ntupleRawTree.b5_mm_kinpc_vtx_yErr[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_z"]=ntupleRawTree.b5_mm_kinpc_vtx_z[mumuIdx] ;
    storageFloat["mm_kinpc_vtx_zErr"]=ntupleRawTree.b5_mm_kinpc_vtx_zErr[mumuIdx] ;
    storageFloat["mm_m1iso"]=ntupleRawTree.b5_mm_m1iso[mumuIdx] ;
    storageFloat["mm_m2iso"]=ntupleRawTree.b5_mm_m2iso[mumuIdx] ;
    storageFloat["mm_mass"]=ntupleRawTree.b5_mm_mass[mumuIdx] ;
    storageFloat["mm_mu1_eta"]=ntupleRawTree.b5_mm_mu1_eta[mumuIdx] ;
    storageFloat["mm_mu1_phi"]=ntupleRawTree.b5_mm_mu1_phi[mumuIdx] ;
    storageFloat["mm_mu1_pt"]=ntupleRawTree.b5_mm_mu1_pt[mumuIdx] ;
    storageFloat["mm_mu2_eta"]=ntupleRawTree.b5_mm_mu2_eta[mumuIdx] ;
    storageFloat["mm_mu2_phi"]=ntupleRawTree.b5_mm_mu2_phi[mumuIdx] ;
    storageFloat["mm_mu2_pt"]=ntupleRawTree.b5_mm_mu2_pt[mumuIdx] ;
    storageFloat["mm_mva"]=ntupleRawTree.b5_mm_mva[mumuIdx] ;
    storageFloat["mm_otherVtxMaxProb"]=ntupleRawTree.b5_mm_otherVtxMaxProb[mumuIdx] ;
    storageFloat["mm_otherVtxMaxProb1"]=ntupleRawTree.b5_mm_otherVtxMaxProb1[mumuIdx] ;
    storageFloat["mm_otherVtxMaxProb2"]=ntupleRawTree.b5_mm_otherVtxMaxProb2[mumuIdx] ;
    storageFloat["mm_closetrk"]=ntupleRawTree.b5_mm_closetrk[mumuIdx] ;
    storageFloat["mm_closetrks1"]=ntupleRawTree.b5_mm_closetrks1[mumuIdx] ;
    storageFloat["mm_closetrks2"]=ntupleRawTree.b5_mm_closetrks2[mumuIdx] ;
    storageFloat["mm_closetrks3"]=ntupleRawTree.b5_mm_closetrks3[mumuIdx] ;
    storageFloat["mm_kal_valid"]=ntupleRawTree.b5_mm_kal_valid[mumuIdx] ;
    storageFloat["mm_kin_valid"]=ntupleRawTree.b5_mm_kin_valid[mumuIdx] ;
    storageFloat["mm_kinpc_valid"]=ntupleRawTree.b5_mm_kinpc_valid[mumuIdx] ;
    storageFloat["mm_mu1_index"]=ntupleRawTree.b5_mm_mu1_index[mumuIdx] ;
    storageFloat["mm_mu1_pdgId"]=ntupleRawTree.b5_mm_mu1_pdgId[mumuIdx] ;
    storageFloat["mm_mu2_index"]=ntupleRawTree.b5_mm_mu2_index[mumuIdx] ;
    storageFloat["mm_mu2_pdgId"]=ntupleRawTree.b5_mm_mu2_pdgId[mumuIdx] ;
    storageFloat["mm_nBMTrks"]=ntupleRawTree.b5_mm_nBMTrks[mumuIdx] ;
    storageFloat["mm_nDisTrks"]=ntupleRawTree.b5_mm_nDisTrks[mumuIdx] ;
    storageFloat["mm_nTrks"]=ntupleRawTree.b5_mm_nTrks[mumuIdx] ;
}

void TreeMaker::bookHistograms()
{
     // All sc 
     for( TString tag : {"allSC_"} )
     {
     }

     for( TString tag : {"gen_" })
     {
        th1fStore[tag+"Multiplicity"            ]  = new TH1F(tag + "Multiplicity" ,"Multiplicity", 200 , -0.5  , 199.5  );
        th1fStore[tag+"Pt"                      ]  = new TH1F(tag + "Pt" ,"Pt of GEN", 100 , 0.0  , 50.0  );
        th1fStore[tag+"Eta"                     ]  = new TH1F(tag + "Eta","Eta of GEN", 100 , -5.0 , 5.0  );
        th1fStore[tag+"phi"                     ]  = new TH1F(tag + "phi","Phi of GEN", 64 , -3.20  , 3.20  );
        th1fStore[tag+"E"		                ]  = new TH1F(tag + "E"  ,"E", 600 , 0.0 ,150.0 ) ;                  
     }
     
     th1fStore["sc_Multiplicity"            ]  = new TH1F("scMultiplicity" ,"Multiplicity", 200 , -0.5  , 199.5  );
     th1fStore["pv_Multiplicity"            ]  = new TH1F("pvMultiplicity" ,"Multiplicity", 200 , -0.5  , 199.5  );
}

void TreeMaker::AddDimuonVertexHistos(TString tag)
{
    
     TString tg;
     tg="mm_bdt"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 202 , -1.01  , 1.01  )  ;
     tg="mm_doca"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 1.0  )  ;
     tg="mm_docatrk"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 20.0  )  ;
     tg="mm_iso"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 1.0  )  ;
     tg="mm_kal_lxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kal_mass"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 800 , 0.0  , 100.0  )  ;
     tg="mm_kal_slxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 1600 , 0.0  , 400.0  )  ;
     tg="mm_kal_vtx_prob"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 120 , -0.1  , 1.1  )  ;
     tg="mm_kin_alpha"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64 ,-3.2  , 3.2  )  ;
     tg="mm_kin_alphaBS"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64, -3.2 , 3.2  )  ;
     tg="mm_kin_alphaBSErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2000 , 0.0  , 4000.0  )  ;
     tg="mm_kin_alphaErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2000 , 0.0  , 4000.0  )  ;
     tg="mm_kin_eta"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , -5.0  , 5.0  )  ;
     tg="mm_kin_l3d"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 240 , 0.0  , 6.0  )  ;
     tg="mm_kin_lxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 240 , 0.0  , 6.0  )  ;
     tg="mm_kin_mass"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 800 , 0.0  , 100.0  )  ;
     tg="mm_kin_massErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 50 , 0.0  , 2.0  )  ;
     tg="mm_kin_mu1eta"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , -5.0  , 5.0  )  ;
     tg="mm_kin_mu1phi"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64 , -3.2  , 3.2  )  ;
     tg="mm_kin_mu1pt"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 500 , 0.0  , 50.0  )  ;
     tg="mm_kin_mu2eta"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , -5.0  , 5.0  )  ;
     tg="mm_kin_mu2phi"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64, -3.2  , 3.2  )  ;
     tg="mm_kin_mu2pt"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 500 , 0.0  , 50.0  )  ;
     tg="mm_kin_phi"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64 , -3.2  , 3.2  )  ;
     tg="mm_kin_pt"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 500 , 0.0  , 50.0  )  ;
     tg="mm_kin_pv2ip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 4.0  )  ;
     tg="mm_kin_pv2ipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 50 , 0.0  ,0.5  )  ;
     tg="mm_kin_pv2lip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg,  200 , -20.0  , 20.0  )  ;
     tg="mm_kin_pv2lipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -20.0  , 20.0  )  ;
     tg="mm_kin_pv2lipSig"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -1500.0  , 1500.0  )  ;
     tg="mm_kin_pv_z"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 400 , -20  , 20.0  )  ;
     tg="mm_kin_pv_zErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 0.05  )  ;
     tg="mm_kin_pvip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 0.5  )  ;
     tg="mm_kin_pvipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 0.2  )  ;
     tg="mm_kin_pvlip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -5.0  , 5.0  )  ;
     tg="mm_kin_pvlipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -2.0  , 2.0  )  ;
     tg="mm_kin_pvlipSig"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 500 , -500.0  , 500.0  )  ;
     tg="mm_kin_sl3d"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 200.0  )  ;
     tg="mm_kin_slxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 200.0  )  ;
     tg="mm_kin_spv2ip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 5.0  )  ;
     tg="mm_kin_spvip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 100.0  )  ;
     tg="mm_kin_tau"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kin_taue"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kin_tauxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kin_tauxye"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kin_vtx_chi2dof"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 800 , 0.0  , 400.0  )  ;
     tg="mm_kin_vtx_prob"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 1.0  )  ;
     tg="mm_kin_vtx_x"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 1000 , -5.0  , 5.0  )  ;
     tg="mm_kin_vtx_xErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 5.0  )  ;
     tg="mm_kin_vtx_y"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -5.0  , 5.0  )  ;
     tg="mm_kin_vtx_yErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 5.0  )  ;
     tg="mm_kin_vtx_z"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -25.0 , 25.0  )  ;
     tg="mm_kin_vtx_zErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 1000 , 0.0  , 100.0  )  ;

     tg="mm_kinpc_alpha"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64 ,-3.2  , 3.2  )  ;
     tg="mm_kinpc_alphaBS"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64, -3.2 , 3.2  )  ;
     tg="mm_kinpc_alphaBSErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2000 , 0.0  , 4000.0  )  ;
     tg="mm_kinpc_alphaErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2000 , 0.0  , 4000.0  )  ;
     tg="mm_kinpc_eta"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , -5.0  , 5.0  )  ;
     tg="mm_kinpc_l3d"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 240 , 0.0  , 6.0  )  ;
     tg="mm_kinpc_lxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 240 , 0.0  , 6.0  )  ;
     tg="mm_kinpc_mass"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 800 , 0.0  , 100.0  )  ;
     tg="mm_kinpc_massErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 50 , 0.0  , 2.0  )  ;
     tg="mm_kinpc_phi"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64 , -3.2  , 3.2  )  ;
     tg="mm_kinpc_pt"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 500 , 0.0  , 50.0  )  ;
     tg="mm_kinpc_pv2ip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 4.0  )  ;
     tg="mm_kinpc_pv2ipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 50 , 0.0  ,0.5  )  ;
     tg="mm_kinpc_pv2lip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg,  200 , -20.0  , 20.0  )  ;
     tg="mm_kinpc_pv2lipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -20.0  , 20.0  )  ;
     tg="mm_kinpc_pv2lipSig"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -1500.0  , 1500.0  )  ;
     tg="mm_kinpc_pv_z"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 400 , -20  , 20.0  )  ;
     tg="mm_kinpc_pv_zErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 0.05  )  ;
     tg="mm_kinpc_pvip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 0.5  )  ;
     tg="mm_kinpc_pvipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 0.2  )  ;
     tg="mm_kinpc_pvlip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -5.0  , 5.0  )  ;
     tg="mm_kinpc_pvlipErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -2.0  , 2.0  )  ;
     tg="mm_kinpc_pvlipSig"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 500 , -500.0  , 500.0  )  ;
     tg="mm_kinpc_sl3d"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 200.0  )  ;
     tg="mm_kinpc_slxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 200.0  )  ;
     tg="mm_kinpc_spv2ip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 5.0  )  ;
     tg="mm_kinpc_spvip"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 100.0  )  ;
     tg="mm_kinpc_tau"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kinpc_taue"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kinpc_tauxy"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kinpc_tauxye"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 10.0  )  ;
     tg="mm_kinpc_vtx_chi2dof"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 800 , 0.0  , 400.0  )  ;
     tg="mm_kinpc_vtx_prob"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 1.0  )  ;
     tg="mm_kinpc_vtx_x"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 1000 , -5.0  , 5.0  )  ;
     tg="mm_kinpc_vtx_xErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 5.0  )  ;
     tg="mm_kinpc_vtx_y"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -5.0  , 5.0  )  ;
     tg="mm_kinpc_vtx_yErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 5.0  )  ;
     tg="mm_kinpc_vtx_z"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -25.0 , 25.0  )  ;
     tg="mm_kinpc_vtx_zErr"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 1000 , 0.0  , 100.0  )  ;
     tg="mm_m1iso"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 1.0  )  ;
     tg="mm_m2iso"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 1.0  )  ;
     tg="mm_mass"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 1000 , 0.0  , 100.0  )  ;
     tg="mm_mu1_eta"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -4.0  , 4.0  )  ;
     tg="mm_mu1_phi"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64 , -3.2  , 3.2 );     
     tg="mm_mu1_pt"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 50.0  )  ;
     tg="mm_mu2_eta"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -4.0  , 4.0  )  ;
     tg="mm_mu2_phi"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 64 , -3.2  , 3.2  )  ;
     tg="mm_mu2_pt"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 50.0  )  ;
     tg="mm_mva"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , -1.0  , 1.0  )  ;
     tg="mm_otherVtxMaxProb"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 1.0  )  ;
     tg="mm_otherVtxMaxProb1"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 1.0  )  ;
     tg="mm_otherVtxMaxProb2"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , 0.0  , 1.0  )  ;
     tg="mm_closetrk"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 101 , -0.5  ,100.5  )  ;
     tg="mm_closetrks1"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 101 , -0.5  , 100.0  )  ;
     tg="mm_closetrks2"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 101 , -0.5  , 100.0  )  ;
     tg="mm_closetrks3"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 101 , -0.5  , 100.0  )  ;
     tg="mm_kal_valid"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2.0 , -0.50  ,1.50  )  ;
     tg="mm_kin_valid"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2.0 , -0.50 , 1.5  )  ;
     tg="mm_kinpc_valid"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2 , -0.5  , 1.5  )  ;
     tg="mm_mu1_pdgId"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 601 , -300.5  , 300.5  )  ;
     tg="mm_mu2_pdgId"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 200 , 0.0  , 20.0  )  ;
     tg="mm_nBMTrks"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , -0.5  , 99.5  )  ;
     tg="mm_nDisTrks"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 100 , -0.50  ,99.5   )  ;
     tg="mm_nTrks"; th1fStore[tag+"_"+tg]  = new TH1F(tag+"_"+tg ,tg, 2001 , -0.5  , 2000.5  )  ;
}

void TreeMaker::fill_eventHists()
{
        th1fStore["gen_Multiplicity"    ]->Fill(eventGenMultiplicity);   
        th1fStore["sc_Multiplicity"    ]->Fill(ntupleRawTree.bG_nSC);   
        th1fStore["pv_Multiplicity"    ]->Fill(ntupleRawTree.bG_nPrimaryVertex);   
}

void TreeMaker::fillDimuonVertexHists(TString tag)
{

   th1fStore[tag+"_"+"mm_bdt"]->Fill(storageFloat["mm_bdt"]);
   th1fStore[tag+"_"+"mm_doca"]->Fill(storageFloat["mm_doca"]);
   th1fStore[tag+"_"+"mm_docatrk"]->Fill(storageFloat["mm_docatrk"]);
   th1fStore[tag+"_"+"mm_iso"]->Fill(storageFloat["mm_iso"]);
   th1fStore[tag+"_"+"mm_kal_lxy"]->Fill(storageFloat["mm_kal_lxy"]);
   th1fStore[tag+"_"+"mm_kal_mass"]->Fill(storageFloat["mm_kal_mass"]);
   th1fStore[tag+"_"+"mm_kal_slxy"]->Fill(storageFloat["mm_kal_slxy"]);
   th1fStore[tag+"_"+"mm_kal_vtx_prob"]->Fill(storageFloat["mm_kal_vtx_prob"]);
   th1fStore[tag+"_"+"mm_kin_alpha"]->Fill(storageFloat["mm_kin_alpha"]);
   th1fStore[tag+"_"+"mm_kin_alphaBS"]->Fill(storageFloat["mm_kin_alphaBS"]);
   th1fStore[tag+"_"+"mm_kin_alphaBSErr"]->Fill(storageFloat["mm_kin_alphaBSErr"]);
   th1fStore[tag+"_"+"mm_kin_alphaErr"]->Fill(storageFloat["mm_kin_alphaErr"]);
   th1fStore[tag+"_"+"mm_kin_eta"]->Fill(storageFloat["mm_kin_eta"]);
   th1fStore[tag+"_"+"mm_kin_l3d"]->Fill(storageFloat["mm_kin_l3d"]);
   th1fStore[tag+"_"+"mm_kin_lxy"]->Fill(storageFloat["mm_kin_lxy"]);
   th1fStore[tag+"_"+"mm_kin_mass"]->Fill(storageFloat["mm_kin_mass"]);
   th1fStore[tag+"_"+"mm_kin_massErr"]->Fill(storageFloat["mm_kin_massErr"]);
   th1fStore[tag+"_"+"mm_kin_mu1eta"]->Fill(storageFloat["mm_kin_mu1eta"]);
   th1fStore[tag+"_"+"mm_kin_mu1phi"]->Fill(storageFloat["mm_kin_mu1phi"]);
   th1fStore[tag+"_"+"mm_kin_mu1pt"]->Fill(storageFloat["mm_kin_mu1pt"]);
   th1fStore[tag+"_"+"mm_kin_mu2eta"]->Fill(storageFloat["mm_kin_mu2eta"]);
   th1fStore[tag+"_"+"mm_kin_mu2phi"]->Fill(storageFloat["mm_kin_mu2phi"]);
   th1fStore[tag+"_"+"mm_kin_mu2pt"]->Fill(storageFloat["mm_kin_mu2pt"]);
   th1fStore[tag+"_"+"mm_kin_phi"]->Fill(storageFloat["mm_kin_phi"]);
   th1fStore[tag+"_"+"mm_kin_pt"]->Fill(storageFloat["mm_kin_pt"]);
   th1fStore[tag+"_"+"mm_kin_pv2ip"]->Fill(storageFloat["mm_kin_pv2ip"]);
   th1fStore[tag+"_"+"mm_kin_pv2ipErr"]->Fill(storageFloat["mm_kin_pv2ipErr"]);
   th1fStore[tag+"_"+"mm_kin_pv2lip"]->Fill(storageFloat["mm_kin_pv2lip"]);
   th1fStore[tag+"_"+"mm_kin_pv2lipErr"]->Fill(storageFloat["mm_kin_pv2lipErr"]);
   th1fStore[tag+"_"+"mm_kin_pv2lipSig"]->Fill(storageFloat["mm_kin_pv2lipSig"]);
   th1fStore[tag+"_"+"mm_kin_pv_z"]->Fill(storageFloat["mm_kin_pv_z"]);
   th1fStore[tag+"_"+"mm_kin_pv_zErr"]->Fill(storageFloat["mm_kin_pv_zErr"]);
   th1fStore[tag+"_"+"mm_kin_pvip"]->Fill(storageFloat["mm_kin_pvip"]);
   th1fStore[tag+"_"+"mm_kin_pvipErr"]->Fill(storageFloat["mm_kin_pvipErr"]);
   th1fStore[tag+"_"+"mm_kin_pvlip"]->Fill(storageFloat["mm_kin_pvlip"]);
   th1fStore[tag+"_"+"mm_kin_pvlipErr"]->Fill(storageFloat["mm_kin_pvlipErr"]);
   th1fStore[tag+"_"+"mm_kin_pvlipSig"]->Fill(storageFloat["mm_kin_pvlipSig"]);
   th1fStore[tag+"_"+"mm_kin_sl3d"]->Fill(storageFloat["mm_kin_sl3d"]);
   th1fStore[tag+"_"+"mm_kin_slxy"]->Fill(storageFloat["mm_kin_slxy"]);
   th1fStore[tag+"_"+"mm_kin_spv2ip"]->Fill(storageFloat["mm_kin_spv2ip"]);
   th1fStore[tag+"_"+"mm_kin_spvip"]->Fill(storageFloat["mm_kin_spvip"]);
   th1fStore[tag+"_"+"mm_kin_tau"]->Fill(storageFloat["mm_kin_tau"]);
   th1fStore[tag+"_"+"mm_kin_taue"]->Fill(storageFloat["mm_kin_taue"]);
   th1fStore[tag+"_"+"mm_kin_tauxy"]->Fill(storageFloat["mm_kin_tauxy"]);
   th1fStore[tag+"_"+"mm_kin_tauxye"]->Fill(storageFloat["mm_kin_tauxye"]);
   th1fStore[tag+"_"+"mm_kin_vtx_chi2dof"]->Fill(storageFloat["mm_kin_vtx_chi2dof"]);
   th1fStore[tag+"_"+"mm_kin_vtx_prob"]->Fill(storageFloat["mm_kin_vtx_prob"]);
   th1fStore[tag+"_"+"mm_kin_vtx_x"]->Fill(storageFloat["mm_kin_vtx_x"]);
   th1fStore[tag+"_"+"mm_kin_vtx_xErr"]->Fill(storageFloat["mm_kin_vtx_xErr"]);
   th1fStore[tag+"_"+"mm_kin_vtx_y"]->Fill(storageFloat["mm_kin_vtx_y"]);
   th1fStore[tag+"_"+"mm_kin_vtx_yErr"]->Fill(storageFloat["mm_kin_vtx_yErr"]);
   th1fStore[tag+"_"+"mm_kin_vtx_z"]->Fill(storageFloat["mm_kin_vtx_z"]);
   th1fStore[tag+"_"+"mm_kin_vtx_zErr"]->Fill(storageFloat["mm_kin_vtx_zErr"]);
   th1fStore[tag+"_"+"mm_kinpc_alpha"]->Fill(storageFloat["mm_kinpc_alpha"]);
   th1fStore[tag+"_"+"mm_kinpc_alphaBS"]->Fill(storageFloat["mm_kinpc_alphaBS"]);
   th1fStore[tag+"_"+"mm_kinpc_alphaBSErr"]->Fill(storageFloat["mm_kinpc_alphaBSErr"]);
   th1fStore[tag+"_"+"mm_kinpc_alphaErr"]->Fill(storageFloat["mm_kinpc_alphaErr"]);
   th1fStore[tag+"_"+"mm_kinpc_eta"]->Fill(storageFloat["mm_kinpc_eta"]);
   th1fStore[tag+"_"+"mm_kinpc_l3d"]->Fill(storageFloat["mm_kinpc_l3d"]);
   th1fStore[tag+"_"+"mm_kinpc_lxy"]->Fill(storageFloat["mm_kinpc_lxy"]);
   th1fStore[tag+"_"+"mm_kinpc_mass"]->Fill(storageFloat["mm_kinpc_mass"]);
   th1fStore[tag+"_"+"mm_kinpc_massErr"]->Fill(storageFloat["mm_kinpc_massErr"]);
   th1fStore[tag+"_"+"mm_kinpc_phi"]->Fill(storageFloat["mm_kinpc_phi"]);
   th1fStore[tag+"_"+"mm_kinpc_pt"]->Fill(storageFloat["mm_kinpc_pt"]);
   th1fStore[tag+"_"+"mm_kinpc_pv2ip"]->Fill(storageFloat["mm_kinpc_pv2ip"]);
   th1fStore[tag+"_"+"mm_kinpc_pv2ipErr"]->Fill(storageFloat["mm_kinpc_pv2ipErr"]);
   th1fStore[tag+"_"+"mm_kinpc_pv2lip"]->Fill(storageFloat["mm_kinpc_pv2lip"]);
   th1fStore[tag+"_"+"mm_kinpc_pv2lipErr"]->Fill(storageFloat["mm_kinpc_pv2lipErr"]);
   th1fStore[tag+"_"+"mm_kinpc_pv2lipSig"]->Fill(storageFloat["mm_kinpc_pv2lipSig"]);
   th1fStore[tag+"_"+"mm_kinpc_pv_z"]->Fill(storageFloat["mm_kinpc_pv_z"]);
   th1fStore[tag+"_"+"mm_kinpc_pv_zErr"]->Fill(storageFloat["mm_kinpc_pv_zErr"]);
   th1fStore[tag+"_"+"mm_kinpc_pvip"]->Fill(storageFloat["mm_kinpc_pvip"]);
   th1fStore[tag+"_"+"mm_kinpc_pvipErr"]->Fill(storageFloat["mm_kinpc_pvipErr"]);
   th1fStore[tag+"_"+"mm_kinpc_pvlip"]->Fill(storageFloat["mm_kinpc_pvlip"]);
   th1fStore[tag+"_"+"mm_kinpc_pvlipErr"]->Fill(storageFloat["mm_kinpc_pvlipErr"]);
   th1fStore[tag+"_"+"mm_kinpc_pvlipSig"]->Fill(storageFloat["mm_kinpc_pvlipSig"]);
   th1fStore[tag+"_"+"mm_kinpc_sl3d"]->Fill(storageFloat["mm_kinpc_sl3d"]);
   th1fStore[tag+"_"+"mm_kinpc_slxy"]->Fill(storageFloat["mm_kinpc_slxy"]);
   th1fStore[tag+"_"+"mm_kinpc_spv2ip"]->Fill(storageFloat["mm_kinpc_spv2ip"]);
   th1fStore[tag+"_"+"mm_kinpc_spvip"]->Fill(storageFloat["mm_kinpc_spvip"]);
   th1fStore[tag+"_"+"mm_kinpc_tau"]->Fill(storageFloat["mm_kinpc_tau"]);
   th1fStore[tag+"_"+"mm_kinpc_taue"]->Fill(storageFloat["mm_kinpc_taue"]);
   th1fStore[tag+"_"+"mm_kinpc_tauxy"]->Fill(storageFloat["mm_kinpc_tauxy"]);
   th1fStore[tag+"_"+"mm_kinpc_tauxye"]->Fill(storageFloat["mm_kinpc_tauxye"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_chi2dof"]->Fill(storageFloat["mm_kinpc_vtx_chi2dof"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_prob"]->Fill(storageFloat["mm_kinpc_vtx_prob"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_x"]->Fill(storageFloat["mm_kinpc_vtx_x"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_xErr"]->Fill(storageFloat["mm_kinpc_vtx_xErr"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_y"]->Fill(storageFloat["mm_kinpc_vtx_y"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_yErr"]->Fill(storageFloat["mm_kinpc_vtx_yErr"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_z"]->Fill(storageFloat["mm_kinpc_vtx_z"]);
   th1fStore[tag+"_"+"mm_kinpc_vtx_zErr"]->Fill(storageFloat["mm_kinpc_vtx_zErr"]);
   th1fStore[tag+"_"+"mm_m1iso"]->Fill(storageFloat["mm_m1iso"]);
   th1fStore[tag+"_"+"mm_m2iso"]->Fill(storageFloat["mm_m2iso"]);
   th1fStore[tag+"_"+"mm_mass"]->Fill(storageFloat["mm_mass"]);
   th1fStore[tag+"_"+"mm_mu1_eta"]->Fill(storageFloat["mm_mu1_eta"]);
   th1fStore[tag+"_"+"mm_mu1_phi"]->Fill(storageFloat["mm_mu1_phi"]);
   th1fStore[tag+"_"+"mm_mu1_pt"]->Fill(storageFloat["mm_mu1_pt"]);
   th1fStore[tag+"_"+"mm_mu2_eta"]->Fill(storageFloat["mm_mu2_eta"]);
   th1fStore[tag+"_"+"mm_mu2_phi"]->Fill(storageFloat["mm_mu2_phi"]);
   th1fStore[tag+"_"+"mm_mu2_pt"]->Fill(storageFloat["mm_mu2_pt"]);
   th1fStore[tag+"_"+"mm_mva"]->Fill(storageFloat["mm_mva"]);
   th1fStore[tag+"_"+"mm_otherVtxMaxProb"]->Fill(storageFloat["mm_otherVtxMaxProb"]);
   th1fStore[tag+"_"+"mm_otherVtxMaxProb1"]->Fill(storageFloat["mm_otherVtxMaxProb1"]);
   th1fStore[tag+"_"+"mm_otherVtxMaxProb2"]->Fill(storageFloat["mm_otherVtxMaxProb2"]);
   th1fStore[tag+"_"+"mm_closetrk"]->Fill(storageFloat["mm_closetrk"]);
   th1fStore[tag+"_"+"mm_closetrks1"]->Fill(storageFloat["mm_closetrks1"]);
   th1fStore[tag+"_"+"mm_closetrks2"]->Fill(storageFloat["mm_closetrks2"]);
   th1fStore[tag+"_"+"mm_closetrks3"]->Fill(storageFloat["mm_closetrks3"]);
   th1fStore[tag+"_"+"mm_kal_valid"]->Fill(storageFloat["mm_kal_valid"]);
   th1fStore[tag+"_"+"mm_kin_valid"]->Fill(storageFloat["mm_kin_valid"]);
   th1fStore[tag+"_"+"mm_kinpc_valid"]->Fill(storageFloat["mm_kinpc_valid"]);
   th1fStore[tag+"_"+"mm_mu1_pdgId"]->Fill(storageFloat["mm_mu1_pdgId"]);
   th1fStore[tag+"_"+"mm_mu2_pdgId"]->Fill(storageFloat["mm_mu2_pdgId"]);
   th1fStore[tag+"_"+"mm_nBMTrks"]->Fill(storageFloat["mm_nBMTrks"]);
   th1fStore[tag+"_"+"mm_nDisTrks"]->Fill(storageFloat["mm_nDisTrks"]);
   th1fStore[tag+"_"+"mm_nTrks"]->Fill(storageFloat["mm_nTrks"]);

}

Int_t TreeMaker::doMuonSelection( Int_t muIdx, bool isLead)
{

    Int_t rslt(0);
    rslt++;
    /*
    std::cout<<" Mu "<<isLead<<" pt,eta,isTraker , isGlobal,isLoose,isHighP,mva : "
                     <<ntupleRawTree.b5_Muon_pt[muIdx]<<" / "<<minMuonPt<<","
                     <<ntupleRawTree.b5_Muon_eta[muIdx]<<" / "<<maxMuonEta<<","
                     <<ntupleRawTree.b5_Muon_isGlobal[muIdx]<<","
                     <<ntupleRawTree.b5_Muon_looseId[muIdx]<<","
                     <<ntupleRawTree.b5_Muon_highPurity[muIdx]<<","
                     <<ntupleRawTree.b5_MuonId_newSoftMuonMva[muIdx]<<" / "<<BDTWorkingPoint
                     <<"\n";
    */
    if(  ntupleRawTree.b5_Muon_pt[muIdx] < minMuonPt )  
		 return  rslt ;
	rslt++;
    if(  abs(ntupleRawTree.b5_Muon_eta[muIdx]) > maxMuonEta )  
		 return  rslt ;
	rslt++;
    if( ( not ntupleRawTree.b5_Muon_isTracker[muIdx] ) and muonHasToBeTracker) 
		 return  rslt ;
	rslt++;
    if( ( not ntupleRawTree.b5_Muon_isGlobal[muIdx]) and muonHasToBeGlobal )
		 return  rslt ;
	rslt++;
    if( (not ntupleRawTree.b5_Muon_looseId[muIdx] ) and muonHasToBeLoose )
		 return  rslt ;
	rslt++;
    if( (not ntupleRawTree.b5_MuonId_highPurity[muIdx]) and muonHasToBeHighPurity )  
 		 return  rslt ;
	rslt++;
    if(ntupleRawTree.b5_MuonId_newSoftMuonMva[muIdx] < BDTWorkingPoint ) 
		 return  rslt ;
    rslt++;
	return  0 ;
}
