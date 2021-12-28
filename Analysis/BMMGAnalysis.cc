void BMMGAnalysis::Analyze()
{
 
    /************************************************************************************

 Make sure the branches used here are not turned off to 0 by BMMGAnalysis::setupBranchStatus()

    *************************************************************************************/

    Double_t dr;
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t EventCountWithDimuCand=0;
    Long64_t EventCountWithDimuVertexCand=0;
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    auto nDiMuNoVertexCandidates=0;
    Int_t prevRun(-1),prevLumi(-1);
    bool goodRunLumi = false;


    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {   

       nDiMuCandidates=0;
       isTriggerd=false;
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);

       if (ientry_evt < 0) break;

       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%500 == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }

       if(doRunLumiMask)
       {
        if( (prevRun != ntupleRawTree.b5_run) or (prevLumi != ntupleRawTree.b5_luminosityBlock) )
        {
             prevRun  = ntupleRawTree.b5_run;
             prevLumi = ntupleRawTree.b5_luminosityBlock;

             goodRunLumi =  runLumiMask.checkRunLumi(prevRun,prevLumi);

             if(not goodRunLumi)
             {
                 std::cout<<"Masking : ( "<<prevRun<<" , "<<prevLumi<<" ) \n";
             }
        }

        if(not goodRunLumi) continue;
       }
       
    //   std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "<<"\n";
       
       // Checks if atleast 1 PV is there .. by default there will always be  one pV , the beamspot
       if(ntupleRawTree.bG_nPrimaryVertex < 1 ) continue;
       
       // Trigger Selection
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Jpsi ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing ) isTriggerd=true;
       //if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R ) isTriggerd=true;
    
       if(not isTriggerd) continue;
       for(Int_t i=0;i<NSTORAGE_ARRAY_MAX;i++)
            storageArrayDouble[i]=0;
       for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
            photonSelectionCheck[phoSCIdx]=-1;
       if(doPhotonMVA)
       {
            doPhotonMVAScores();
       }
                

       int rslt=0;
       nBMMGCandidates=0;

       fill_muonHists();
       fill_scHists();
       fill_photonHists();
       nDiMuNoVertexCandidates=0;
       //std::cout<<"\nnMM = "<<ntupleRawTree.b5_nmm<<" , nmu : "<< ntupleRawTree.b5_nMuon<<"\n";
      
       for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    
           if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] < 0 or ntupleRawTree.b5_mm_mu2_index[mumuIdx] <0 ) continue;  
           if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) continue;
           if(ntupleRawTree.b5_mm_mu2_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) continue;

           fill_dimuonHists(mumuIdx);
           //std::cout<<"\tmumuIdx : "<<mumuIdx<<" : m1,m2 : "<<ntupleRawTree.b5_mm_mu1_index[mumuIdx]<<", "<<ntupleRawTree.b5_mm_mu2_index[mumuIdx]<<"\n";
           //std::cout<<"\tmumuIdx : "<<mumuIdx<<" : m1pt,m2pt : "<<ntupleRawTree.b5_mm_mu1_pt[mumuIdx]<<", "<<ntupleRawTree.b5_mm_mu2_pt[mumuIdx]<<"\n";
           //std::cout<<"\tmumu    : "<<mumuIdx<<" : MASS : "<<ntupleRawTree.b5_mm_mass[mumuIdx]<<"\n";
           // Muon Selection
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu1_index[mumuIdx], true);
           // std::cout<<"\tmuon1 rslt : "<<rslt<<"\n";
           if(rslt > 0) continue;
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu2_index[mumuIdx], false);
           // std::cout<<"\tmuon2 rslt : "<<rslt<<"\n";
           if(rslt > 0) continue;
            
           // Dimuon Selection
           rslt=doDimuonSelection(mumuIdx);
           //std::cout<<"\tDimuon rslt : "<<rslt<<"\n";
           if(rslt > 0) continue;        
           nDiMuNoVertexCandidates++;
           rslt = doVertexSelection(mumuIdx);
           //std::cout<<"\tDimuon Vertex rslt : "<<rslt<<"\n";
           if(rslt > 0) continue;
           

           auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
                      ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
           storageArrayDouble[nDiMuCandidates + candidateMapDouble["mumu_dr"] ]   = dr ;
           
           // VERTEX SELECTION STUFF



           fill_dimuonPassHists(mumuIdx);
           fill_dimuonEnvironmentHists(mumuIdx);


           diMuLV.SetPtEtaPhiM(     ntupleRawTree.b5_mm_kin_pt[mumuIdx],   \
                                    ntupleRawTree.b5_mm_kin_eta[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_phi[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_mass[mumuIdx]  );
           nBMMGCandidatesPerDimu=0;

           for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
           {
               // photon selection
               if(photonSelectionCheck[phoSCIdx] < 0)
               {
                    photonSelectionCheck[phoSCIdx]=doPhotonSelection(phoSCIdx );
               }

               if( photonSelectionCheck[phoSCIdx] > 0) continue;
               
               rslt=doBMMGSelection(mumuIdx,phoSCIdx);
               if(rslt > 0) continue;
               
               auto et= ntupleRawTree.bG_scE[phoSCIdx]/cosh(ntupleRawTree.bG_scEta[phoSCIdx]);
               photonLV.SetPtEtaPhiM( et,ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx], PHO_MASS );
               bmmgLV = diMuLV + photonLV;
                
               fill_bmmgHists( bmmgLV , mumuIdx, phoSCIdx );

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
           nDiMuCandidates++;

           if(nDiMuCandidates >= NDIMU_MAX)
           {
                std::cout<<" Dimuon count per event above NDIMU_MAX , Aborting !! \n";
                exit(9);
           }
            
       }
       
    
       EventCount++;
       if(nDiMuNoVertexCandidates > 0)
       {
            EventCountWithDimuCand++;
       }
       if(nDiMuCandidates > 0)
       {
            EventCountWithDimuVertexCand++;
       }
       if(nBMMGCandidates >0)
       {
        // std::cout<<"\t "<<ntupleRawTree.b5_run<<" , "<<ntupleRawTree.b5_event<<" : Number of nBMMGCandidates = "<<nBMMGCandidates<<"\n";
        EventCountWithCand++;
       }
       
       fill_globalEventHists();

     if(outTree)  outTree->Fill();
    }
    std::cout<<"\n\n"
            <<"  Number of Events with trigger = "<<EventCount<<"\n"
            <<"  Number of Events with dimu candidates = "<<EventCountWithDimuCand<<"\n"
            <<"  Number of Events with dimu vtx candidates = "<<EventCountWithDimuVertexCand<<"\n"
            <<"  Number of Events with candidates = "<<EventCountWithCand<<"\n"
            <<"  Analysis loop completed"
            <<"  \n\n";

}

void BMMGAnalysis::GenAnalyze()
{

 
    /************************************************************************************

 Make sure the branches used here are not turned off to 0 by BMMGAnalysis::setupBranchStatus()

    *************************************************************************************/

    Double_t dr;
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithTrigger=0;
    Long64_t EventCountWithCand=0;
    Long64_t EventCountWithDimuCand=0;
    Long64_t EventCountWithDimuVertexCand=0;
    Long64_t nb = 0,nbytes=0 ;

    Long64_t mumMatchCount=0;
    Long64_t mupMatchCount=0;
    Long64_t scMatchCount=0;
    Long64_t phoMatchCount=0;
    Long64_t dimuMatchCount=0;
    Long64_t fullEventMatchesFound=0;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    auto nDiMuNoVertexCandidates=0;
    Int_t prevRun(-1),prevLumi(-1);
    bool goodRunLumi = false;


    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {   

       nDiMuCandidates=0;
       isTriggerd=false;
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);

       if (ientry_evt < 0) break;

       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%1000 == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
      
      
         //std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "<<"\n";
       
       // Checks if atleast 1 PV is there .. by default there will always be  one pV , the beamspot
       if(ntupleRawTree.bG_nPrimaryVertex < 1 ) continue;
       
    //   Trigger Selection
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Jpsi ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R ) isTriggerd=true;
       //if(not isTriggerd) continue;
       if( isTriggerd ) EventCountWithTrigger++;
      for(Int_t i=0;i<NSTORAGE_ARRAY_MAX;i++)
            storageArrayDouble[i]=0;
       for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
            photonSelectionCheck[phoSCIdx]=-1;
       if(doPhotonMVA)
       {
            doPhotonMVAScores();
       }
                

       int rslt=0;
       nBMMGCandidates=0;
       
       Int_t mumMatchIdx(-1),mupMatchIdx(-1),phoMatchIdx(-1);
       Int_t mumFoundIdx(-1),mupFoundIdx(-1),phoFoundIdx(-1);
       Int_t motherIdx(-1);
       /*for(Int_t i=0; i< ntupleRawTree.b5_nGenPart; i++)
       {
            std::cout<<"\t"<<i<<"/"<<ntupleRawTree.b5_nGenPart<<" , "<< ntupleRawTree.b5_GenPart_pdgId[i]  <<" from ";
            if (ntupleRawTree.b5_GenPart_genPartIdxMother[i] >0 )
                std::cout<<ntupleRawTree.b5_GenPart_pdgId[ ntupleRawTree.b5_GenPart_genPartIdxMother[i]  ];
            else
            std::cout<<ntupleRawTree.b5_GenPart_genPartIdxMother[i];
            std::cout<<"\n";
        
       }*/

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
                    //if(phoFoundIdx==1 and mupFoundIdx==1 and mumFoundIdx==1)
                    if(mupFoundIdx==1 and mumFoundIdx==1)
                    {
                     //   std::cout<<" Match Idx s : "<<mupMatchIdx<<", "<<mumMatchIdx<<", "<<phoMatchIdx<<"\n";
                        break;
                    }
                    else
                    {
                        mumMatchIdx=-1; mupMatchIdx=-1; phoMatchIdx=-1;
                    }

            }
       }
       if(mumFoundIdx != 1 or mupFoundIdx !=1 )
       {
            std::cout<<" No Gen Bs Decay found for jentry =  "<<jentry<<" !!\n ";
            continue;
       }

        TLorentzVector mumLV,mupLV;

       if( phoFoundIdx == 0 )
       {
            bmmgLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[motherIdx], ntupleRawTree.b5_GenPart_eta[motherIdx],
                                ntupleRawTree.b5_GenPart_phi[motherIdx], ntupleRawTree.b5_GenPart_mass[motherIdx] )    ;
            mumLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[mumMatchIdx], ntupleRawTree.b5_GenPart_eta[mumMatchIdx],
                                 ntupleRawTree.b5_GenPart_phi[mumMatchIdx], ntupleRawTree.b5_GenPart_mass[mumMatchIdx] ) ;    
            mupLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[mupMatchIdx], ntupleRawTree.b5_GenPart_eta[mupMatchIdx],
                                ntupleRawTree.b5_GenPart_phi[mupMatchIdx], ntupleRawTree.b5_GenPart_mass[mupMatchIdx] )   ; 
            photonLV=bmmgLV- ( mumLV  +  mupLV );
            
       //     std::cout<<"Recalculating Pho P4 pt  = "<<photonLV.Pt()<<" , eta = "<<photonLV.Eta()<<"\n" ;
       }
       else if( phoFoundIdx==1)
       {

            photonLV.SetPtEtaPhiM( ntupleRawTree.b5_GenPart_pt[phoMatchIdx], ntupleRawTree.b5_GenPart_eta[phoMatchIdx],
                                ntupleRawTree.b5_GenPart_phi[phoMatchIdx], ntupleRawTree.b5_GenPart_mass[phoMatchIdx] )   ; 
       }
      // std::cout<<"Mu - eta, phi : "<<ntupleRawTree.b5_GenPart_eta[mumMatchIdx]<<" , "<<ntupleRawTree.b5_GenPart_phi[mumMatchIdx]<<"\n";
      // std::cout<<"Mu + eta, phi : "<<ntupleRawTree.b5_GenPart_eta[mupMatchIdx]<<" , "<<ntupleRawTree.b5_GenPart_phi[mupMatchIdx]<<"\n";
      // std::cout<<"\n";

       Int_t mumRecoMatchIdx(-1),mupRecoMatchIdx(-1),phoRecoMatchIdx(-1);
       Float_t mumRecoMatchDr(1e9),mupRecoMatchDr(1e9),phoRecoMatchDr(1e9);
       for(int i=0;i<ntupleRawTree.b5_nMuon ; i++)
       {
        //    std::cout<<"i = "<<i<<"/"<<ntupleRawTree.b5_nMuon<<" : "<<ntupleRawTree.b5_Muon_pdgId[i]<<" [ "<<ntupleRawTree.b5_Muon_eta[i]<<" , "<<ntupleRawTree.b5_Muon_phi[i]<<" ] "<<"\n";
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
       
       th1fStore["gen_mumDeltaR"]->Fill(mumRecoMatchDr);
       th1fStore["gen_mupDeltaR"]->Fill(mupRecoMatchDr);
       th1fStore["gen_phoDeltaR"]->Fill(phoRecoMatchDr);
      // std::cout<<"mup , mum , pho Dr : "<<mumRecoMatchDr<<" , "<<mupRecoMatchDr<<" , "<<phoRecoMatchDr<<"\n";
       if(mumRecoMatchDr < 0.1 ) { 
          mumMatchCount++;
		fill_muonHists(mumRecoMatchIdx);
	 //	std::cout<<"Match Found for : mu-"<<"\n";
	  }
       if(mupRecoMatchDr < 0.1 ) { 
          mupMatchCount++;
		fill_muonHists(mupRecoMatchIdx);
	//	std::cout<<"Match Found for : mu+"<<"\n";
	  }
       if(phoRecoMatchDr < 0.1 ) { 
          scMatchCount++;
		fill_scHists(phoRecoMatchIdx);
	//	std::cout<<"Match Found for : pho"<<"\n";
	  }
       if(phoRecoMatchDr < 0.1 )  
       {
        if(doPhotonSelection(phoRecoMatchIdx) ==0)
        {   
            phoMatchCount++;
        }
		fill_photonHists(phoRecoMatchIdx);
       }
       if(mumRecoMatchDr < 0.1 and mupRecoMatchDr <0.1 and phoRecoMatchDr < 0.1)
       {
            fullEventMatchesFound++;
       }
       nDiMuNoVertexCandidates=0;
      
       for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    
           if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] < 0 or ntupleRawTree.b5_mm_mu2_index[mumuIdx] <0 ) continue;  
           if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) continue;
           if(ntupleRawTree.b5_mm_mu2_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) continue;

           if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] != mumRecoMatchIdx and ntupleRawTree.b5_mm_mu2_index[mumuIdx] != mumRecoMatchIdx ) continue;  
           if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] != mupRecoMatchIdx and ntupleRawTree.b5_mm_mu2_index[mumuIdx] != mupRecoMatchIdx ) continue;  
           
           dimuMatchCount++;
           // Muon Selection
           fill_dimuonHists(mumuIdx);
           fill_dimuonEnvironmentHists(mumuIdx);
		   
                 
           // Dimuon Selection
           rslt=doDimuonSelection(mumuIdx);
           if(rslt > 0) continue;        
           //nDiMuNoVertexCandidates++;
           //rslt = doVertexSelection(mumuIdx);
           //std::cout<<"\tDimuon Vertex rslt : "<<rslt<<"\n";
           //if(rslt > 0) continue;
           

           auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
                      ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
           storageArrayDouble[nDiMuCandidates + candidateMapDouble["mumu_dr"] ]   = dr ;
           
           // VERTEX SELECTION STUFF

           fill_dimuonPassHists(mumuIdx);
           diMuLV.SetPtEtaPhiM(     ntupleRawTree.b5_mm_kin_pt[mumuIdx],   \
                                    ntupleRawTree.b5_mm_kin_eta[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_phi[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_mass[mumuIdx]  );
           nBMMGCandidatesPerDimu=0;

           for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
           {
                if(phoSCIdx != phoRecoMatchIdx) continue;

               // photon selection
               if(photonSelectionCheck[phoSCIdx] < 0)
               {
                    photonSelectionCheck[phoSCIdx]=doPhotonSelection(phoSCIdx );
               }

               if( photonSelectionCheck[phoSCIdx] > 0) continue;
               
               rslt=doBMMGSelection(mumuIdx,phoSCIdx);
               //if(rslt > 0) continue;
               
               auto et= ntupleRawTree.bG_scE[phoSCIdx]/cosh(ntupleRawTree.bG_scEta[phoSCIdx]);
               photonLV.SetPtEtaPhiM( et,ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx], PHO_MASS );
               bmmgLV = diMuLV + photonLV;
                
               fill_bmmgHists( bmmgLV , mumuIdx, phoSCIdx );

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
           nDiMuCandidates++;

           if(nDiMuCandidates >= NDIMU_MAX)
           {
                std::cout<<" Dimuon count per event above NDIMU_MAX , Aborting !! \n";
                exit(9);
           }
            
       }
       
    
       EventCount++;
       if(nDiMuNoVertexCandidates > 0)
       {
            EventCountWithDimuCand++;
       }
       if(nDiMuCandidates > 0)
       {
            EventCountWithDimuVertexCand++;
       }
       if(nBMMGCandidates >0)
       {
        // std::cout<<"\t "<<ntupleRawTree.b5_run<<" , "<<ntupleRawTree.b5_event<<" : Number of nBMMGCandidates = "<<nBMMGCandidates<<"\n";
        EventCountWithCand++;
       }
       
       fill_globalEventHists();

     if(outTree)  outTree->Fill();
    }
    std::cout<<"\n\n"
            <<"  Number of Events                           = "<<EventCount<<"\n"
            <<"  Number of Events with trigger              = "<<EventCountWithTrigger<<"\n"
            <<"  Number of Events with mu+ reconstructed    = "<<mupMatchCount<<"\n"
            <<"  Number of Events with mu- reconstructed    = "<<mumMatchCount<<"\n"
            <<"  Number of Events with sc reconstructed     = "<<scMatchCount<<"\n"
            <<"  Number of Events with pho reconstructed    = "<<phoMatchCount<<"\n"
            <<"  Number of Events with dimu reconstructed   = "<<dimuMatchCount<<"\n"
            <<"  Number of Events with all reconstructed    = "<<fullEventMatchesFound<<"\n"
            <<"  Analysis loop completed"
            <<"  \n\n";




}
void BMMGAnalysis::setUpPhotonMVA()
{
    if(not hasWeightFiles)
    {
        std::cout<<"Weights are not provided exiting ! \n";
        exit(4);
    }

    std::cout<<" ******** Setting up PhotonMVA ID *********\n"
             <<"Weight file : "<<photonIdxMVAWeightFile<<"\n";
    reader =  new TMVA::Reader( "!Color:!Silent" );
    
    reader->AddVariable("et", &(photonMVAdata.et));
    reader->AddVariable("rawE", &(photonMVAdata.rawE));
    reader->AddVariable("sigmaIetaIeta", &(photonMVAdata.sigmaIetaIeta));
    //reader->AddVariable("sigmaIetaIphi", &(photonMVAdata.sigmaIetaIphi));
    reader->AddVariable("FoundGsfMatch", &(photonMVAdata.FoundGsfMatch));
    reader->AddVariable("r9", &(photonMVAdata.r9));
    reader->AddVariable("etaWidth", &(photonMVAdata.etaWidth));
    reader->AddVariable("PFPhoIso", &(photonMVAdata.PFPhoIso));
    reader->AddVariable("PFNeuIso", &(photonMVAdata.PFNeuIso)); 
    reader->AddVariable("PFChIso", &(photonMVAdata.PFChIso)); 
    reader->AddVariable("full5x5_e5x5", &(photonMVAdata.full5x5_e5x5));
    reader->AddVariable("full5x5_r9", &(photonMVAdata.full5x5_r9));
    //reader->AddVariable("e2x5_MaxRatio", &(photonMVAdata.e2x5_MaxRatio));
    reader->BookMVA("LowPtPhotonIdMVA_MLP", photonIdxMVAWeightFile );

    hasSetupPhotonMVA=true;
}

void BMMGAnalysis::doPhotonMVAScores()
{   
    if(not hasSetupPhotonMVA)
    {
        std::cout<<" MVA not setup : "<<hasSetupPhotonMVA<<"\n";
    }
    for(int  i=0  ; i < ntupleRawTree.bG_nSC ; i++)
    {
      photonMVAdata.energy        = ntupleRawTree.bG_scE[i] ;
      photonMVAdata.et            = ntupleRawTree.bG_scEt[i] ;
      photonMVAdata.eta           = ntupleRawTree.bG_scEta[i] ;
      photonMVAdata.rawE          = ntupleRawTree.bG_scRawE[i] ;
      photonMVAdata.full5x5_e5x5  = ntupleRawTree.bG_scFull5x5_e5x5[i] ;
      photonMVAdata.full5x5_r9    = ntupleRawTree.bG_scFull5x5_r9[i] ;
      photonMVAdata.sigmaIetaIeta = ntupleRawTree.bG_scSigmaIetaIeta[i] ;
   //   photonMVAdata.sigmaIetaIphi = ntupleRawTree.bG_scSigmaIetaIphi[i] ;
      photonMVAdata.etaWidth      = ntupleRawTree.bG_scEtaWidth[i] ;
      photonMVAdata.phiWidth      = ntupleRawTree.bG_scPhiWidth[i] ;
      photonMVAdata.r9            = ntupleRawTree.bG_scR9[i] ;
      photonMVAdata.swisross      = ntupleRawTree.bG_scSwissCross[i] ;
      photonMVAdata.PFPhoIso      = ntupleRawTree.bG_scPFPhoIso3[i] ;
      photonMVAdata.PFNeuIso      = ntupleRawTree.bG_scPFNeuIso3[i] ;
      photonMVAdata.PFChIso       = ntupleRawTree.bG_scPFChIso3[i] ;
      photonMVAdata.FoundGsfMatch = ntupleRawTree.bG_scFoundGsfMatch_[i] ;
  //    photonMVAdata.e2x5_MaxRatio = ntupleRawTree.bG_sc ; 
     
      photonMVAValue = reader->EvaluateMVA("LowPtPhotonIdMVA_MLP");
      storageArrayDouble[ i + candidateMapInt["scPhotonMVAScore"]  ]   = photonMVAValue ;
   }
     

}

void BMMGAnalysis::AllocateBMMGBranches()
{
    outTree->Branch("nDiMuCandidates",&(nDiMuCandidates));
    outTree->Branch("nBMMGCandidates",&(nBMMGCandidates));
    
    outTree->Branch("isTriggerd",&(isTriggerd));

    candidateMapInt["nBMMGCandidatesPerDimu"]   = storageIdxFilledInt ;
    outTree->Branch("nBMMGCandidatesPerDimu",&(storageArrayInt[storageIdxFilledInt]),"nBMMGCandidatesPerDimu[nDiMuCandidates]/I");storageIdxFilledInt+=NDIMU_MAX;
    
    if(doPhotonMVA)
    {
        candidateMapDouble["scPhotonMVAScore"]   = storageIdxFilledDouble ;
        outTree->Branch("scPhotonMVAScore",&storageArrayDouble[storageIdxFilledDouble],"scPhotonMVAScore[bG_nSC]/D"); storageIdxFilledDouble+=NSC_MAX;
    }

    candidateMapDouble["mumu_dr"]   = storageIdxFilledDouble ;
    outTree->Branch("mumu_dr",&(storageArrayDouble[storageIdxFilledDouble]),"mumu_dr[nDiMuCandidates]/D");storageIdxFilledDouble+=NDIMU_MAX;

    candidateMapInt["bmmg_dimuon_idx"]       = storageIdxFilledInt ;
    outTree->Branch("bmmg_dimuon_idx",&(storageArrayInt[storageIdxFilledInt]),"dimuon_idx[nBMMGCandidates]/I" )   ;   storageIdxFilledInt+=NBMMG_MAX;
    
    candidateMapInt["bmmg_photonSC_idx"]       = storageIdxFilledInt ;
    outTree->Branch("bmmg_photonSC_idx",&(storageArrayInt[storageIdxFilledInt]),"photonSC_idx[nBMMGCandidates]/I" )   ;   storageIdxFilledInt+=NBMMG_MAX;
    
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
    
    std::cout<<" Storage usage : " <<storageIdxFilledDouble<<" doubles and "<<storageIdxFilledInt<<" ints.\n ";
}

void BMMGAnalysis::setupReducedAnalysisTreeBranches()
{
   if(ntupleRawTree.fChain) 
   {
        std::cout<<"input Tree  Not Setup !!\n";
        exit(14);

   }
   outTree->Branch("bG_run", &(ntupleRawTree.bG_run));
   outTree->Branch("bG_event", &(ntupleRawTree.bG_event));
   outTree->Branch("bG_lumis", &(ntupleRawTree.bG_lumis));
   outTree->Branch("bG_isData", &(ntupleRawTree.bG_isData));
   outTree->Branch("bG_nPrimaryVertex", &(ntupleRawTree.bG_nPrimaryVertex));
   outTree->Branch("bG_primaryVertex_isFake", (ntupleRawTree.bG_primaryVertex_isFake));
   outTree->Branch("bG_primaryVertex_x", (ntupleRawTree.bG_primaryVertex_x));
   outTree->Branch("bG_primaryVertex_y", (ntupleRawTree.bG_primaryVertex_y));
   outTree->Branch("bG_primaryVertex_z", (ntupleRawTree.bG_primaryVertex_z));
   outTree->Branch("bG_primaryVertex_t", (ntupleRawTree.bG_primaryVertex_t));
   outTree->Branch("bG_primaryVertex_covXX", (ntupleRawTree.bG_primaryVertex_covXX));
   outTree->Branch("bG_primaryVertex_covXY", (ntupleRawTree.bG_primaryVertex_covXY));
   outTree->Branch("bG_primaryVertex_covXZ", (ntupleRawTree.bG_primaryVertex_covXZ));
   outTree->Branch("bG_primaryVertex_covYY", (ntupleRawTree.bG_primaryVertex_covYY));
   outTree->Branch("bG_primaryVertex_covYZ", (ntupleRawTree.bG_primaryVertex_covYZ));
   outTree->Branch("bG_primaryVertex_covZZ", (ntupleRawTree.bG_primaryVertex_covZZ));
   outTree->Branch("bG_primaryVertex_x_error", (ntupleRawTree.bG_primaryVertex_x_error));
   outTree->Branch("bG_primaryVertex_y_error", (ntupleRawTree.bG_primaryVertex_y_error));
   outTree->Branch("bG_primaryVertex_z_error", (ntupleRawTree.bG_primaryVertex_z_error));
   outTree->Branch("bG_primaryVertex_t_error", (ntupleRawTree.bG_primaryVertex_t_error));
   outTree->Branch("bG_primaryVertex_ntracks", (ntupleRawTree.bG_primaryVertex_ntracks));
   outTree->Branch("bG_primaryVertex_ndof", (ntupleRawTree.bG_primaryVertex_ndof));
   outTree->Branch("bG_primaryVertex_chi2", (ntupleRawTree.bG_primaryVertex_chi2));
   outTree->Branch("bG_primaryVertex_normalizedChi2", (ntupleRawTree.bG_primaryVertex_normalizedChi2));
   outTree->Branch("bG_nPho", &(ntupleRawTree.bG_nPho));
   outTree->Branch("bG_phoE", (ntupleRawTree.bG_phoE));
   outTree->Branch("bG_phoEt", (ntupleRawTree.bG_phoEt));
   outTree->Branch("bG_phoEta", (ntupleRawTree.bG_phoEta));
   outTree->Branch("bG_phoPhi", (ntupleRawTree.bG_phoPhi));
   outTree->Branch("bG_phoSCE", (ntupleRawTree.bG_phoSCE));
   outTree->Branch("bG_phoSCEt", (ntupleRawTree.bG_phoSCEt));
   outTree->Branch("bG_phoSCRawE", (ntupleRawTree.bG_phoSCRawE));
   outTree->Branch("bG_phoESEnP1", (ntupleRawTree.bG_phoESEnP1));
   outTree->Branch("bG_phoESEnP2", (ntupleRawTree.bG_phoESEnP2));
   outTree->Branch("bG_phoSCEta", (ntupleRawTree.bG_phoSCEta));
   outTree->Branch("bG_phoSCPhi", (ntupleRawTree.bG_phoSCPhi));
   outTree->Branch("bG_phoSCEtaWidth", (ntupleRawTree.bG_phoSCEtaWidth));
   outTree->Branch("bG_phoSCPhiWidth", (ntupleRawTree.bG_phoSCPhiWidth));
   outTree->Branch("bG_phoSCBrem", (ntupleRawTree.bG_phoSCBrem));
   outTree->Branch("bG_phohasPixelSeed", (ntupleRawTree.bG_phohasPixelSeed));
   outTree->Branch("bG_phoR9", (ntupleRawTree.bG_phoR9));
   outTree->Branch("bG_phoHoverE", (ntupleRawTree.bG_phoHoverE));
   outTree->Branch("bG_phoESEffSigmaRR", (ntupleRawTree.bG_phoESEffSigmaRR));
   outTree->Branch("bG_phoSigmaIEtaIEtaFull5x5", (ntupleRawTree.bG_phoSigmaIEtaIEtaFull5x5));
   outTree->Branch("bG_phoSigmaIEtaIPhiFull5x5", (ntupleRawTree.bG_phoSigmaIEtaIPhiFull5x5));
   outTree->Branch("bG_phoSigmaIPhiIPhiFull5x5", (ntupleRawTree.bG_phoSigmaIPhiIPhiFull5x5));
   outTree->Branch("bG_phoE2x2Full5x5", (ntupleRawTree.bG_phoE2x2Full5x5));
   outTree->Branch("bG_phoE5x5Full5x5", (ntupleRawTree.bG_phoE5x5Full5x5));
   outTree->Branch("bG_phoR9Full5x5", (ntupleRawTree.bG_phoR9Full5x5));
   outTree->Branch("bG_phoPFChIso", (ntupleRawTree.bG_phoPFChIso));
   outTree->Branch("bG_phoPFPhoIso", (ntupleRawTree.bG_phoPFPhoIso));
   outTree->Branch("bG_phoPFNeuIso", (ntupleRawTree.bG_phoPFNeuIso));
   outTree->Branch("bG_phoEcalPFClusterIso", (ntupleRawTree.bG_phoEcalPFClusterIso));
   outTree->Branch("bG_phoHcalPFClusterIso", (ntupleRawTree.bG_phoHcalPFClusterIso));
   outTree->Branch("bG_phoSeedTime", (ntupleRawTree.bG_phoSeedTime));
   outTree->Branch("bG_phoSeedEnergy", (ntupleRawTree.bG_phoSeedEnergy));
   outTree->Branch("bG_phoMIPTotEnergy", (ntupleRawTree.bG_phoMIPTotEnergy));
   outTree->Branch("bG_phoMIPChi2", (ntupleRawTree.bG_phoMIPChi2));
   outTree->Branch("bG_phoMIPSlope", (ntupleRawTree.bG_phoMIPSlope));
   outTree->Branch("bG_phoMIPIntercept", (ntupleRawTree.bG_phoMIPIntercept));
   outTree->Branch("bG_phoMIPNhitCone", (ntupleRawTree.bG_phoMIPNhitCone));
   outTree->Branch("bG_phoMIPIsHalo", (ntupleRawTree.bG_phoMIPIsHalo));
   outTree->Branch("bG_nPFPho", &(ntupleRawTree.bG_nPFPho));
   outTree->Branch("bG_phoPFE", (ntupleRawTree.bG_phoPFE));
   outTree->Branch("bG_phoPFEt", (ntupleRawTree.bG_phoPFEt));
   outTree->Branch("bG_phoPFEta", (ntupleRawTree.bG_phoPFEta));
   outTree->Branch("bG_phoPFPhi", (ntupleRawTree.bG_phoPFPhi));
   outTree->Branch("bG_nSC", &(ntupleRawTree.bG_nSC));
   outTree->Branch("bG_scE", (ntupleRawTree.bG_scE));
   outTree->Branch("bG_scEt", (ntupleRawTree.bG_scEt));
   outTree->Branch("bG_scRawE", (ntupleRawTree.bG_scRawE));
   outTree->Branch("bG_scEta", (ntupleRawTree.bG_scEta));
   outTree->Branch("bG_scPhi", (ntupleRawTree.bG_scPhi));
   outTree->Branch("bG_scX", (ntupleRawTree.bG_scX));
   outTree->Branch("bG_scY", (ntupleRawTree.bG_scY));
   outTree->Branch("bG_scZ", (ntupleRawTree.bG_scZ));
   outTree->Branch("bG_scEtaWidth", (ntupleRawTree.bG_scEtaWidth));
   outTree->Branch("bG_scPhiWidth", (ntupleRawTree.bG_scPhiWidth));
   outTree->Branch("bG_scRawEt", (ntupleRawTree.bG_scRawEt));
   outTree->Branch("bG_scMinDrWithGsfElectornSC_", (ntupleRawTree.bG_scMinDrWithGsfElectornSC_));
   outTree->Branch("bG_scFoundGsfMatch_", (ntupleRawTree.bG_scFoundGsfMatch_));
   outTree->Branch("bG_scE5x5", (ntupleRawTree.bG_scE5x5));
   outTree->Branch("bG_scE2x2Ratio", (ntupleRawTree.bG_scE2x2Ratio));
   outTree->Branch("bG_scE3x3Ratio", (ntupleRawTree.bG_scE3x3Ratio));
   outTree->Branch("bG_scEMaxRatio", (ntupleRawTree.bG_scEMaxRatio));
   outTree->Branch("bG_scE2ndRatio", (ntupleRawTree.bG_scE2ndRatio));
   outTree->Branch("bG_scETopRatio", (ntupleRawTree.bG_scETopRatio));
   outTree->Branch("bG_scERightRatio", (ntupleRawTree.bG_scERightRatio));
   outTree->Branch("bG_scEBottomRatio", (ntupleRawTree.bG_scEBottomRatio));
   outTree->Branch("bG_scELeftRatio", (ntupleRawTree.bG_scELeftRatio));
   outTree->Branch("bG_scE2x5MaxRatio", (ntupleRawTree.bG_scE2x5MaxRatio));
   outTree->Branch("bG_scE2x5TopRatio", (ntupleRawTree.bG_scE2x5TopRatio));
   outTree->Branch("bG_scE2x5RightRatio", (ntupleRawTree.bG_scE2x5RightRatio));
   outTree->Branch("bG_scE2x5BottomRatio", (ntupleRawTree.bG_scE2x5BottomRatio));
   outTree->Branch("bG_scE2x5LeftRatio", (ntupleRawTree.bG_scE2x5LeftRatio));
   outTree->Branch("bG_scSwissCross", (ntupleRawTree.bG_scSwissCross));
   outTree->Branch("bG_scR9", (ntupleRawTree.bG_scR9));
   outTree->Branch("bG_scSigmaIetaIeta", (ntupleRawTree.bG_scSigmaIetaIeta));
   outTree->Branch("bG_scSigmaIetaIphi", (ntupleRawTree.bG_scSigmaIetaIphi));
   outTree->Branch("bG_scSigmaIphiIphi", (ntupleRawTree.bG_scSigmaIphiIphi));
   outTree->Branch("bG_scFull5x5_e5x5", (ntupleRawTree.bG_scFull5x5_e5x5));
   outTree->Branch("bG_scFull5x5_e2x2Ratio", (ntupleRawTree.bG_scFull5x5_e2x2Ratio));
   outTree->Branch("bG_scFull5x5_e3x3Ratio", (ntupleRawTree.bG_scFull5x5_e3x3Ratio));
   outTree->Branch("bG_scFull5x5_eMaxRatio", (ntupleRawTree.bG_scFull5x5_eMaxRatio));
   outTree->Branch("bG_scFull5x5_e2ndRatio", (ntupleRawTree.bG_scFull5x5_e2ndRatio));
   outTree->Branch("bG_scFull5x5_eTopRatio", (ntupleRawTree.bG_scFull5x5_eTopRatio));
   outTree->Branch("bG_scFull5x5_eRightRatio", (ntupleRawTree.bG_scFull5x5_eRightRatio));
   outTree->Branch("bG_scFull5x5_eBottomRatio", (ntupleRawTree.bG_scFull5x5_eBottomRatio));
   outTree->Branch("bG_scFull5x5_eLeftRatio", (ntupleRawTree.bG_scFull5x5_eLeftRatio));
   outTree->Branch("bG_scFull5x5_e2x5MaxRatio", (ntupleRawTree.bG_scFull5x5_e2x5MaxRatio));
   outTree->Branch("bG_scFull5x5_e2x5TopRatio", (ntupleRawTree.bG_scFull5x5_e2x5TopRatio));
   outTree->Branch("bG_scFull5x5_e2x5RightRatio", (ntupleRawTree.bG_scFull5x5_e2x5RightRatio));
   outTree->Branch("bG_scFull5x5_e2x5BottomRatio", (ntupleRawTree.bG_scFull5x5_e2x5BottomRatio));
   outTree->Branch("bG_scFull5x5_e2x5LeftRatio", (ntupleRawTree.bG_scFull5x5_e2x5LeftRatio));
   outTree->Branch("bG_scFull5x5_swissCross", (ntupleRawTree.bG_scFull5x5_swissCross));
   outTree->Branch("bG_scFull5x5_r9", (ntupleRawTree.bG_scFull5x5_r9));
   outTree->Branch("bG_scFull5x5_sigmaIetaIeta", (ntupleRawTree.bG_scFull5x5_sigmaIetaIeta));
   outTree->Branch("bG_scFull5x5_sigmaIetaIphi", (ntupleRawTree.bG_scFull5x5_sigmaIetaIphi));
   outTree->Branch("bG_scFull5x5_sigmaIphiIphi", (ntupleRawTree.bG_scFull5x5_sigmaIphiIphi));
   outTree->Branch("bG_scPFChIso1", (ntupleRawTree.bG_scPFChIso1));
   outTree->Branch("bG_scPFChIso2", (ntupleRawTree.bG_scPFChIso2));
   outTree->Branch("bG_scPFChIso3", (ntupleRawTree.bG_scPFChIso3));
   outTree->Branch("bG_scPFChIso4", (ntupleRawTree.bG_scPFChIso4));
   outTree->Branch("bG_scPFChIso5", (ntupleRawTree.bG_scPFChIso5));
   outTree->Branch("bG_scPFPhoIso1", (ntupleRawTree.bG_scPFPhoIso1));
   outTree->Branch("bG_scPFPhoIso2", (ntupleRawTree.bG_scPFPhoIso2));
   outTree->Branch("bG_scPFPhoIso3", (ntupleRawTree.bG_scPFPhoIso3));
   outTree->Branch("bG_scPFPhoIso4", (ntupleRawTree.bG_scPFPhoIso4));
   outTree->Branch("bG_scPFPhoIso5", (ntupleRawTree.bG_scPFPhoIso5));
   outTree->Branch("bG_scPFNeuIso1", (ntupleRawTree.bG_scPFNeuIso1));
   outTree->Branch("bG_scPFNeuIso2", (ntupleRawTree.bG_scPFNeuIso2));
   outTree->Branch("bG_scPFNeuIso3", (ntupleRawTree.bG_scPFNeuIso3));
   outTree->Branch("bG_scPFNeuIso4", (ntupleRawTree.bG_scPFNeuIso4));
   outTree->Branch("bG_scPFNeuIso5", (ntupleRawTree.bG_scPFNeuIso5));
   outTree->Branch("bG_nhcalRechit", &(ntupleRawTree.bG_nhcalRechit));
   outTree->Branch("bG_hcalRechitIEta", (ntupleRawTree.bG_hcalRechitIEta));
   outTree->Branch("bG_hcalRechitIPhi", (ntupleRawTree.bG_hcalRechitIPhi));
   outTree->Branch("bG_hcalRechitEnergy", (ntupleRawTree.bG_hcalRechitEnergy));
   outTree->Branch("b5_run", &(ntupleRawTree.b5_run));
   outTree->Branch("b5_luminosityBlock", &(ntupleRawTree.b5_luminosityBlock));
   outTree->Branch("b5_event", &(ntupleRawTree.b5_event));
   outTree->Branch("b5_nMuonId", &(ntupleRawTree.b5_nMuonId));
   outTree->Branch("b5_MuonId_chi2LocalPosition", (ntupleRawTree.b5_MuonId_chi2LocalPosition));
   outTree->Branch("b5_MuonId_glbNormChi2", (ntupleRawTree.b5_MuonId_glbNormChi2));
   outTree->Branch("b5_MuonId_glbTrackProbability", (ntupleRawTree.b5_MuonId_glbTrackProbability));
   outTree->Branch("b5_MuonId_match1_dX", (ntupleRawTree.b5_MuonId_match1_dX));
   outTree->Branch("b5_MuonId_match1_dY", (ntupleRawTree.b5_MuonId_match1_dY));
   outTree->Branch("b5_MuonId_match1_pullDxDz", (ntupleRawTree.b5_MuonId_match1_pullDxDz));
   outTree->Branch("b5_MuonId_match1_pullDyDz", (ntupleRawTree.b5_MuonId_match1_pullDyDz));
   outTree->Branch("b5_MuonId_match1_pullX", (ntupleRawTree.b5_MuonId_match1_pullX));
   outTree->Branch("b5_MuonId_match1_pullY", (ntupleRawTree.b5_MuonId_match1_pullY));
   outTree->Branch("b5_MuonId_match2_dX", (ntupleRawTree.b5_MuonId_match2_dX));
   outTree->Branch("b5_MuonId_match2_dY", (ntupleRawTree.b5_MuonId_match2_dY));
   outTree->Branch("b5_MuonId_match2_pullDxDz", (ntupleRawTree.b5_MuonId_match2_pullDxDz));
   outTree->Branch("b5_MuonId_match2_pullDyDz", (ntupleRawTree.b5_MuonId_match2_pullDyDz));
   outTree->Branch("b5_MuonId_match2_pullX", (ntupleRawTree.b5_MuonId_match2_pullX));
   outTree->Branch("b5_MuonId_match2_pullY", (ntupleRawTree.b5_MuonId_match2_pullY));
   outTree->Branch("b5_MuonId_newSoftMuonMva", (ntupleRawTree.b5_MuonId_newSoftMuonMva));
   outTree->Branch("b5_MuonId_simProdRho", (ntupleRawTree.b5_MuonId_simProdRho));
   outTree->Branch("b5_MuonId_simProdZ", (ntupleRawTree.b5_MuonId_simProdZ));
   outTree->Branch("b5_MuonId_trkKink", (ntupleRawTree.b5_MuonId_trkKink));
   outTree->Branch("b5_MuonId_trkValidFrac", (ntupleRawTree.b5_MuonId_trkValidFrac));
   outTree->Branch("b5_MuonId_highPurity", (ntupleRawTree.b5_MuonId_highPurity));
   outTree->Branch("b5_MuonId_nLostHitsInner", (ntupleRawTree.b5_MuonId_nLostHitsInner));
   outTree->Branch("b5_MuonId_nLostHitsOn", (ntupleRawTree.b5_MuonId_nLostHitsOn));
   outTree->Branch("b5_MuonId_nLostHitsOuter", (ntupleRawTree.b5_MuonId_nLostHitsOuter));
   outTree->Branch("b5_MuonId_nPixels", (ntupleRawTree.b5_MuonId_nPixels));
   outTree->Branch("b5_MuonId_nValidHits", (ntupleRawTree.b5_MuonId_nValidHits));
   outTree->Branch("b5_MuonId_simExtType", (ntupleRawTree.b5_MuonId_simExtType));
   outTree->Branch("b5_MuonId_simMotherPdgId", (ntupleRawTree.b5_MuonId_simMotherPdgId));
   outTree->Branch("b5_MuonId_simPdgId", (ntupleRawTree.b5_MuonId_simPdgId));
   outTree->Branch("b5_MuonId_simType", (ntupleRawTree.b5_MuonId_simType));
   outTree->Branch("b5_MuonId_trkLayers", (ntupleRawTree.b5_MuonId_trkLayers));
   outTree->Branch("b5_MuonId_trkLostLayersInner", (ntupleRawTree.b5_MuonId_trkLostLayersInner));
   outTree->Branch("b5_MuonId_trkLostLayersOn", (ntupleRawTree.b5_MuonId_trkLostLayersOn));
   outTree->Branch("b5_MuonId_trkLostLayersOuter", (ntupleRawTree.b5_MuonId_trkLostLayersOuter));
   outTree->Branch("b5_nmm", &(ntupleRawTree.b5_nmm));
   outTree->Branch("b5_mm_bdt", (ntupleRawTree.b5_mm_bdt));
   outTree->Branch("b5_mm_doca", (ntupleRawTree.b5_mm_doca));
   outTree->Branch("b5_mm_docatrk", (ntupleRawTree.b5_mm_docatrk));
   outTree->Branch("b5_mm_gen_alpha_ip", (ntupleRawTree.b5_mm_gen_alpha_ip));
   outTree->Branch("b5_mm_gen_alpha_p_phi", (ntupleRawTree.b5_mm_gen_alpha_p_phi));
   outTree->Branch("b5_mm_gen_alpha_p_theta", (ntupleRawTree.b5_mm_gen_alpha_p_theta));
   outTree->Branch("b5_mm_gen_alpha_vtx", (ntupleRawTree.b5_mm_gen_alpha_vtx));
   outTree->Branch("b5_mm_gen_doca", (ntupleRawTree.b5_mm_gen_doca));
   outTree->Branch("b5_mm_gen_l3d", (ntupleRawTree.b5_mm_gen_l3d));
   outTree->Branch("b5_mm_gen_lxy", (ntupleRawTree.b5_mm_gen_lxy));
   outTree->Branch("b5_mm_gen_mass", (ntupleRawTree.b5_mm_gen_mass));
   outTree->Branch("b5_mm_gen_mu1_pt", (ntupleRawTree.b5_mm_gen_mu1_pt));
   outTree->Branch("b5_mm_gen_mu2_pt", (ntupleRawTree.b5_mm_gen_mu2_pt));
   outTree->Branch("b5_mm_gen_prod_z", (ntupleRawTree.b5_mm_gen_prod_z));
   outTree->Branch("b5_mm_gen_pt", (ntupleRawTree.b5_mm_gen_pt));
   outTree->Branch("b5_mm_gen_tau", (ntupleRawTree.b5_mm_gen_tau));
   outTree->Branch("b5_mm_gen_vtx_x", (ntupleRawTree.b5_mm_gen_vtx_x));
   outTree->Branch("b5_mm_gen_vtx_y", (ntupleRawTree.b5_mm_gen_vtx_y));
   outTree->Branch("b5_mm_gen_vtx_z", (ntupleRawTree.b5_mm_gen_vtx_z));
   outTree->Branch("b5_mm_iso", (ntupleRawTree.b5_mm_iso));
   outTree->Branch("b5_mm_kal_lxy", (ntupleRawTree.b5_mm_kal_lxy));
   outTree->Branch("b5_mm_kal_mass", (ntupleRawTree.b5_mm_kal_mass));
   outTree->Branch("b5_mm_kal_slxy", (ntupleRawTree.b5_mm_kal_slxy));
   outTree->Branch("b5_mm_kal_vtx_prob", (ntupleRawTree.b5_mm_kal_vtx_prob));
   outTree->Branch("b5_mm_kin_alpha", (ntupleRawTree.b5_mm_kin_alpha));
   outTree->Branch("b5_mm_kin_alphaBS", (ntupleRawTree.b5_mm_kin_alphaBS));
   outTree->Branch("b5_mm_kin_alphaBSErr", (ntupleRawTree.b5_mm_kin_alphaBSErr));
   outTree->Branch("b5_mm_kin_alphaErr", (ntupleRawTree.b5_mm_kin_alphaErr));
   outTree->Branch("b5_mm_kin_eta", (ntupleRawTree.b5_mm_kin_eta));
   outTree->Branch("b5_mm_kin_l3d", (ntupleRawTree.b5_mm_kin_l3d));
   outTree->Branch("b5_mm_kin_lxy", (ntupleRawTree.b5_mm_kin_lxy));
   outTree->Branch("b5_mm_kin_mass", (ntupleRawTree.b5_mm_kin_mass));
   outTree->Branch("b5_mm_kin_massErr", (ntupleRawTree.b5_mm_kin_massErr));
   outTree->Branch("b5_mm_kin_mu1eta", (ntupleRawTree.b5_mm_kin_mu1eta));
   outTree->Branch("b5_mm_kin_mu1phi", (ntupleRawTree.b5_mm_kin_mu1phi));
   outTree->Branch("b5_mm_kin_mu1pt", (ntupleRawTree.b5_mm_kin_mu1pt));
   outTree->Branch("b5_mm_kin_mu2eta", (ntupleRawTree.b5_mm_kin_mu2eta));
   outTree->Branch("b5_mm_kin_mu2phi", (ntupleRawTree.b5_mm_kin_mu2phi));
   outTree->Branch("b5_mm_kin_mu2pt", (ntupleRawTree.b5_mm_kin_mu2pt));
   outTree->Branch("b5_mm_kin_phi", (ntupleRawTree.b5_mm_kin_phi));
   outTree->Branch("b5_mm_kin_pt", (ntupleRawTree.b5_mm_kin_pt));
   outTree->Branch("b5_mm_kin_pv2ip", (ntupleRawTree.b5_mm_kin_pv2ip));
   outTree->Branch("b5_mm_kin_pv2ipErr", (ntupleRawTree.b5_mm_kin_pv2ipErr));
   outTree->Branch("b5_mm_kin_pv2lip", (ntupleRawTree.b5_mm_kin_pv2lip));
   outTree->Branch("b5_mm_kin_pv2lipErr", (ntupleRawTree.b5_mm_kin_pv2lipErr));
   outTree->Branch("b5_mm_kin_pv2lipSig", (ntupleRawTree.b5_mm_kin_pv2lipSig));
   outTree->Branch("b5_mm_kin_pv_z", (ntupleRawTree.b5_mm_kin_pv_z));
   outTree->Branch("b5_mm_kin_pv_zErr", (ntupleRawTree.b5_mm_kin_pv_zErr));
   outTree->Branch("b5_mm_kin_pvip", (ntupleRawTree.b5_mm_kin_pvip));
   outTree->Branch("b5_mm_kin_pvipErr", (ntupleRawTree.b5_mm_kin_pvipErr));
   outTree->Branch("b5_mm_kin_pvlip", (ntupleRawTree.b5_mm_kin_pvlip));
   outTree->Branch("b5_mm_kin_pvlipErr", (ntupleRawTree.b5_mm_kin_pvlipErr));
   outTree->Branch("b5_mm_kin_pvlipSig", (ntupleRawTree.b5_mm_kin_pvlipSig));
   outTree->Branch("b5_mm_kin_sl3d", (ntupleRawTree.b5_mm_kin_sl3d));
   outTree->Branch("b5_mm_kin_slxy", (ntupleRawTree.b5_mm_kin_slxy));
   outTree->Branch("b5_mm_kin_spv2ip", (ntupleRawTree.b5_mm_kin_spv2ip));
   outTree->Branch("b5_mm_kin_spvip", (ntupleRawTree.b5_mm_kin_spvip));
   outTree->Branch("b5_mm_kin_tau", (ntupleRawTree.b5_mm_kin_tau));
   outTree->Branch("b5_mm_kin_taue", (ntupleRawTree.b5_mm_kin_taue));
   outTree->Branch("b5_mm_kin_tauxy", (ntupleRawTree.b5_mm_kin_tauxy));
   outTree->Branch("b5_mm_kin_tauxye", (ntupleRawTree.b5_mm_kin_tauxye));
   outTree->Branch("b5_mm_kin_vtx_chi2dof", (ntupleRawTree.b5_mm_kin_vtx_chi2dof));
   outTree->Branch("b5_mm_kin_vtx_prob", (ntupleRawTree.b5_mm_kin_vtx_prob));
   outTree->Branch("b5_mm_kin_vtx_x", (ntupleRawTree.b5_mm_kin_vtx_x));
   outTree->Branch("b5_mm_kin_vtx_xErr", (ntupleRawTree.b5_mm_kin_vtx_xErr));
   outTree->Branch("b5_mm_kin_vtx_y", (ntupleRawTree.b5_mm_kin_vtx_y));
   outTree->Branch("b5_mm_kin_vtx_yErr", (ntupleRawTree.b5_mm_kin_vtx_yErr));
   outTree->Branch("b5_mm_kin_vtx_z", (ntupleRawTree.b5_mm_kin_vtx_z));
   outTree->Branch("b5_mm_kin_vtx_zErr", (ntupleRawTree.b5_mm_kin_vtx_zErr));
   outTree->Branch("b5_mm_kinpc_alpha", (ntupleRawTree.b5_mm_kinpc_alpha));
   outTree->Branch("b5_mm_kinpc_alphaBS", (ntupleRawTree.b5_mm_kinpc_alphaBS));
   outTree->Branch("b5_mm_kinpc_alphaBSErr", (ntupleRawTree.b5_mm_kinpc_alphaBSErr));
   outTree->Branch("b5_mm_kinpc_alphaErr", (ntupleRawTree.b5_mm_kinpc_alphaErr));
   outTree->Branch("b5_mm_kinpc_eta", (ntupleRawTree.b5_mm_kinpc_eta));
   outTree->Branch("b5_mm_kinpc_l3d", (ntupleRawTree.b5_mm_kinpc_l3d));
   outTree->Branch("b5_mm_kinpc_lxy", (ntupleRawTree.b5_mm_kinpc_lxy));
   outTree->Branch("b5_mm_kinpc_mass", (ntupleRawTree.b5_mm_kinpc_mass));
   outTree->Branch("b5_mm_kinpc_massErr", (ntupleRawTree.b5_mm_kinpc_massErr));
   outTree->Branch("b5_mm_kinpc_phi", (ntupleRawTree.b5_mm_kinpc_phi));
   outTree->Branch("b5_mm_kinpc_pt", (ntupleRawTree.b5_mm_kinpc_pt));
   outTree->Branch("b5_mm_kinpc_pv2ip", (ntupleRawTree.b5_mm_kinpc_pv2ip));
   outTree->Branch("b5_mm_kinpc_pv2ipErr", (ntupleRawTree.b5_mm_kinpc_pv2ipErr));
   outTree->Branch("b5_mm_kinpc_pv2lip", (ntupleRawTree.b5_mm_kinpc_pv2lip));
   outTree->Branch("b5_mm_kinpc_pv2lipErr", (ntupleRawTree.b5_mm_kinpc_pv2lipErr));
   outTree->Branch("b5_mm_kinpc_pv2lipSig", (ntupleRawTree.b5_mm_kinpc_pv2lipSig));
   outTree->Branch("b5_mm_kinpc_pv_z", (ntupleRawTree.b5_mm_kinpc_pv_z));
   outTree->Branch("b5_mm_kinpc_pv_zErr", (ntupleRawTree.b5_mm_kinpc_pv_zErr));
   outTree->Branch("b5_mm_kinpc_pvip", (ntupleRawTree.b5_mm_kinpc_pvip));
   outTree->Branch("b5_mm_kinpc_pvipErr", (ntupleRawTree.b5_mm_kinpc_pvipErr));
   outTree->Branch("b5_mm_kinpc_pvlip", (ntupleRawTree.b5_mm_kinpc_pvlip));
   outTree->Branch("b5_mm_kinpc_pvlipErr", (ntupleRawTree.b5_mm_kinpc_pvlipErr));
   outTree->Branch("b5_mm_kinpc_pvlipSig", (ntupleRawTree.b5_mm_kinpc_pvlipSig));
   outTree->Branch("b5_mm_kinpc_sl3d", (ntupleRawTree.b5_mm_kinpc_sl3d));
   outTree->Branch("b5_mm_kinpc_slxy", (ntupleRawTree.b5_mm_kinpc_slxy));
   outTree->Branch("b5_mm_kinpc_spv2ip", (ntupleRawTree.b5_mm_kinpc_spv2ip));
   outTree->Branch("b5_mm_kinpc_spvip", (ntupleRawTree.b5_mm_kinpc_spvip));
   outTree->Branch("b5_mm_kinpc_tau", (ntupleRawTree.b5_mm_kinpc_tau));
   outTree->Branch("b5_mm_kinpc_taue", (ntupleRawTree.b5_mm_kinpc_taue));
   outTree->Branch("b5_mm_kinpc_tauxy", (ntupleRawTree.b5_mm_kinpc_tauxy));
   outTree->Branch("b5_mm_kinpc_tauxye", (ntupleRawTree.b5_mm_kinpc_tauxye));
   outTree->Branch("b5_mm_kinpc_vtx_chi2dof", (ntupleRawTree.b5_mm_kinpc_vtx_chi2dof));
   outTree->Branch("b5_mm_kinpc_vtx_prob", (ntupleRawTree.b5_mm_kinpc_vtx_prob));
   outTree->Branch("b5_mm_kinpc_vtx_x", (ntupleRawTree.b5_mm_kinpc_vtx_x));
   outTree->Branch("b5_mm_kinpc_vtx_xErr", (ntupleRawTree.b5_mm_kinpc_vtx_xErr));
   outTree->Branch("b5_mm_kinpc_vtx_y", (ntupleRawTree.b5_mm_kinpc_vtx_y));
   outTree->Branch("b5_mm_kinpc_vtx_yErr", (ntupleRawTree.b5_mm_kinpc_vtx_yErr));
   outTree->Branch("b5_mm_kinpc_vtx_z", (ntupleRawTree.b5_mm_kinpc_vtx_z));
   outTree->Branch("b5_mm_kinpc_vtx_zErr", (ntupleRawTree.b5_mm_kinpc_vtx_zErr));
   outTree->Branch("b5_mm_m1iso", (ntupleRawTree.b5_mm_m1iso));
   outTree->Branch("b5_mm_m2iso", (ntupleRawTree.b5_mm_m2iso));
   outTree->Branch("b5_mm_mass", (ntupleRawTree.b5_mm_mass));
   outTree->Branch("b5_mm_mu1_eta", (ntupleRawTree.b5_mm_mu1_eta));
   outTree->Branch("b5_mm_mu1_phi", (ntupleRawTree.b5_mm_mu1_phi));
   outTree->Branch("b5_mm_mu1_pt", (ntupleRawTree.b5_mm_mu1_pt));
   outTree->Branch("b5_mm_mu2_eta", (ntupleRawTree.b5_mm_mu2_eta));
   outTree->Branch("b5_mm_mu2_phi", (ntupleRawTree.b5_mm_mu2_phi));
   outTree->Branch("b5_mm_mu2_pt", (ntupleRawTree.b5_mm_mu2_pt));
   outTree->Branch("b5_mm_mva", (ntupleRawTree.b5_mm_mva));
   outTree->Branch("b5_mm_otherVtxMaxProb", (ntupleRawTree.b5_mm_otherVtxMaxProb));
   outTree->Branch("b5_mm_otherVtxMaxProb1", (ntupleRawTree.b5_mm_otherVtxMaxProb1));
   outTree->Branch("b5_mm_otherVtxMaxProb2", (ntupleRawTree.b5_mm_otherVtxMaxProb2));
   outTree->Branch("b5_mm_closetrk", (ntupleRawTree.b5_mm_closetrk));
   outTree->Branch("b5_mm_closetrks1", (ntupleRawTree.b5_mm_closetrks1));
   outTree->Branch("b5_mm_closetrks2", (ntupleRawTree.b5_mm_closetrks2));
   outTree->Branch("b5_mm_closetrks3", (ntupleRawTree.b5_mm_closetrks3));
   outTree->Branch("b5_mm_gen_cpdgId", (ntupleRawTree.b5_mm_gen_cpdgId));
   outTree->Branch("b5_mm_gen_mpdgId", (ntupleRawTree.b5_mm_gen_mpdgId));
   outTree->Branch("b5_mm_gen_mu1_mpdgId", (ntupleRawTree.b5_mm_gen_mu1_mpdgId));
   outTree->Branch("b5_mm_gen_mu1_pdgId", (ntupleRawTree.b5_mm_gen_mu1_pdgId));
   outTree->Branch("b5_mm_gen_mu2_mpdgId", (ntupleRawTree.b5_mm_gen_mu2_mpdgId));
   outTree->Branch("b5_mm_gen_mu2_pdgId", (ntupleRawTree.b5_mm_gen_mu2_pdgId));
   outTree->Branch("b5_mm_gen_pdgId", (ntupleRawTree.b5_mm_gen_pdgId));
   outTree->Branch("b5_mm_kal_valid", (ntupleRawTree.b5_mm_kal_valid));
   outTree->Branch("b5_mm_kin_valid", (ntupleRawTree.b5_mm_kin_valid));
   outTree->Branch("b5_mm_kinpc_valid", (ntupleRawTree.b5_mm_kinpc_valid));
   outTree->Branch("b5_mm_mu1_index", (ntupleRawTree.b5_mm_mu1_index));
   outTree->Branch("b5_mm_mu1_pdgId", (ntupleRawTree.b5_mm_mu1_pdgId));
   outTree->Branch("b5_mm_mu2_index", (ntupleRawTree.b5_mm_mu2_index));
   outTree->Branch("b5_mm_mu2_pdgId", (ntupleRawTree.b5_mm_mu2_pdgId));
   outTree->Branch("b5_mm_nBMTrks", (ntupleRawTree.b5_mm_nBMTrks));
   outTree->Branch("b5_mm_nDisTrks", (ntupleRawTree.b5_mm_nDisTrks));
   outTree->Branch("b5_mm_nTrks", (ntupleRawTree.b5_mm_nTrks));
   outTree->Branch("b5_ngensummary", &(ntupleRawTree.b5_ngensummary));
   outTree->Branch("b5_gensummary_n_anti_b", (ntupleRawTree.b5_gensummary_n_anti_b));
   outTree->Branch("b5_gensummary_n_b", (ntupleRawTree.b5_gensummary_n_b));
   outTree->Branch("b5_gensummary_process_type", (ntupleRawTree.b5_gensummary_process_type));
   outTree->Branch("b5_ngenbmm", &(ntupleRawTree.b5_ngenbmm));
   outTree->Branch("b5_genbmm_dau3_eta", (ntupleRawTree.b5_genbmm_dau3_eta));
   outTree->Branch("b5_genbmm_dau3_phi", (ntupleRawTree.b5_genbmm_dau3_phi));
   outTree->Branch("b5_genbmm_dau3_pt", (ntupleRawTree.b5_genbmm_dau3_pt));
   outTree->Branch("b5_genbmm_dau3_reco_eta", (ntupleRawTree.b5_genbmm_dau3_reco_eta));
   outTree->Branch("b5_genbmm_dau3_reco_phi", (ntupleRawTree.b5_genbmm_dau3_reco_phi));
   outTree->Branch("b5_genbmm_dau3_reco_pt", (ntupleRawTree.b5_genbmm_dau3_reco_pt));
   outTree->Branch("b5_genbmm_dau4_eta", (ntupleRawTree.b5_genbmm_dau4_eta));
   outTree->Branch("b5_genbmm_dau4_phi", (ntupleRawTree.b5_genbmm_dau4_phi));
   outTree->Branch("b5_genbmm_dau4_pt", (ntupleRawTree.b5_genbmm_dau4_pt));
   outTree->Branch("b5_genbmm_dau4_reco_eta", (ntupleRawTree.b5_genbmm_dau4_reco_eta));
   outTree->Branch("b5_genbmm_dau4_reco_phi", (ntupleRawTree.b5_genbmm_dau4_reco_phi));
   outTree->Branch("b5_genbmm_dau4_reco_pt", (ntupleRawTree.b5_genbmm_dau4_reco_pt));
   outTree->Branch("b5_genbmm_dimuon_mass", (ntupleRawTree.b5_genbmm_dimuon_mass));
   outTree->Branch("b5_genbmm_eta", (ntupleRawTree.b5_genbmm_eta));
   outTree->Branch("b5_genbmm_mass", (ntupleRawTree.b5_genbmm_mass));
   outTree->Branch("b5_genbmm_mu1_eta", (ntupleRawTree.b5_genbmm_mu1_eta));
   outTree->Branch("b5_genbmm_mu1_phi", (ntupleRawTree.b5_genbmm_mu1_phi));
   outTree->Branch("b5_genbmm_mu1_pt", (ntupleRawTree.b5_genbmm_mu1_pt));
   outTree->Branch("b5_genbmm_mu2_eta", (ntupleRawTree.b5_genbmm_mu2_eta));
   outTree->Branch("b5_genbmm_mu2_phi", (ntupleRawTree.b5_genbmm_mu2_phi));
   outTree->Branch("b5_genbmm_mu2_pt", (ntupleRawTree.b5_genbmm_mu2_pt));
   outTree->Branch("b5_genbmm_phi", (ntupleRawTree.b5_genbmm_phi));
   outTree->Branch("b5_genbmm_pt", (ntupleRawTree.b5_genbmm_pt));
   outTree->Branch("b5_genbmm_rad_eta", (ntupleRawTree.b5_genbmm_rad_eta));
   outTree->Branch("b5_genbmm_rad_p", (ntupleRawTree.b5_genbmm_rad_p));
   outTree->Branch("b5_genbmm_rad_phi", (ntupleRawTree.b5_genbmm_rad_phi));
   outTree->Branch("b5_genbmm_rad_pt", (ntupleRawTree.b5_genbmm_rad_pt));
   outTree->Branch("b5_genbmm_dau3_pdgId", (ntupleRawTree.b5_genbmm_dau3_pdgId));
   outTree->Branch("b5_genbmm_dau4_pdgId", (ntupleRawTree.b5_genbmm_dau4_pdgId));
   outTree->Branch("b5_genbmm_mu1_good", (ntupleRawTree.b5_genbmm_mu1_good));
   outTree->Branch("b5_genbmm_mu1_index", (ntupleRawTree.b5_genbmm_mu1_index));
   outTree->Branch("b5_genbmm_mu1_pdgId", (ntupleRawTree.b5_genbmm_mu1_pdgId));
   outTree->Branch("b5_genbmm_mu2_good", (ntupleRawTree.b5_genbmm_mu2_good));
   outTree->Branch("b5_genbmm_mu2_index", (ntupleRawTree.b5_genbmm_mu2_index));
   outTree->Branch("b5_genbmm_mu2_pdgId", (ntupleRawTree.b5_genbmm_mu2_pdgId));
   outTree->Branch("b5_genbmm_pdgId", (ntupleRawTree.b5_genbmm_pdgId));
   outTree->Branch("b5_genbmm_signature", (ntupleRawTree.b5_genbmm_signature));
   outTree->Branch("b5_nGenPart", &(ntupleRawTree.b5_nGenPart));
   outTree->Branch("b5_GenPart_eta", (ntupleRawTree.b5_GenPart_eta));
   outTree->Branch("b5_GenPart_mass", (ntupleRawTree.b5_GenPart_mass));
   outTree->Branch("b5_GenPart_phi", (ntupleRawTree.b5_GenPart_phi));
   outTree->Branch("b5_GenPart_pt", (ntupleRawTree.b5_GenPart_pt));
   outTree->Branch("b5_GenPart_genPartIdxMother", (ntupleRawTree.b5_GenPart_genPartIdxMother));
   outTree->Branch("b5_GenPart_pdgId", (ntupleRawTree.b5_GenPart_pdgId));
   outTree->Branch("b5_GenPart_status", (ntupleRawTree.b5_GenPart_status));
   outTree->Branch("b5_GenPart_statusFlags", (ntupleRawTree.b5_GenPart_statusFlags));
   outTree->Branch("b5_nMuon", &(ntupleRawTree.b5_nMuon));
   outTree->Branch("b5_Muon_dxy", (ntupleRawTree.b5_Muon_dxy));
   outTree->Branch("b5_Muon_dxyErr", (ntupleRawTree.b5_Muon_dxyErr));
   outTree->Branch("b5_Muon_dxybs", (ntupleRawTree.b5_Muon_dxybs));
   outTree->Branch("b5_Muon_dz", (ntupleRawTree.b5_Muon_dz));
   outTree->Branch("b5_Muon_dzErr", (ntupleRawTree.b5_Muon_dzErr));
   outTree->Branch("b5_Muon_eta", (ntupleRawTree.b5_Muon_eta));
   outTree->Branch("b5_Muon_ip3d", (ntupleRawTree.b5_Muon_ip3d));
   outTree->Branch("b5_Muon_jetPtRelv2", (ntupleRawTree.b5_Muon_jetPtRelv2));
   outTree->Branch("b5_Muon_jetRelIso", (ntupleRawTree.b5_Muon_jetRelIso));
   outTree->Branch("b5_Muon_mass", (ntupleRawTree.b5_Muon_mass));
   outTree->Branch("b5_Muon_miniPFRelIso_all", (ntupleRawTree.b5_Muon_miniPFRelIso_all));
   outTree->Branch("b5_Muon_miniPFRelIso_chg", (ntupleRawTree.b5_Muon_miniPFRelIso_chg));
   outTree->Branch("b5_Muon_pfRelIso03_all", (ntupleRawTree.b5_Muon_pfRelIso03_all));
   outTree->Branch("b5_Muon_pfRelIso03_chg", (ntupleRawTree.b5_Muon_pfRelIso03_chg));
   outTree->Branch("b5_Muon_pfRelIso04_all", (ntupleRawTree.b5_Muon_pfRelIso04_all));
   outTree->Branch("b5_Muon_phi", (ntupleRawTree.b5_Muon_phi));
   outTree->Branch("b5_Muon_pt", (ntupleRawTree.b5_Muon_pt));
   outTree->Branch("b5_Muon_ptErr", (ntupleRawTree.b5_Muon_ptErr));
   outTree->Branch("b5_Muon_segmentComp", (ntupleRawTree.b5_Muon_segmentComp));
   outTree->Branch("b5_Muon_sip3d", (ntupleRawTree.b5_Muon_sip3d));
   outTree->Branch("b5_Muon_softMva", (ntupleRawTree.b5_Muon_softMva));
   outTree->Branch("b5_Muon_tkRelIso", (ntupleRawTree.b5_Muon_tkRelIso));
   outTree->Branch("b5_Muon_tunepRelPt", (ntupleRawTree.b5_Muon_tunepRelPt));
   outTree->Branch("b5_Muon_mvaLowPt", (ntupleRawTree.b5_Muon_mvaLowPt));
   outTree->Branch("b5_Muon_mvaTTH", (ntupleRawTree.b5_Muon_mvaTTH));
   outTree->Branch("b5_Muon_charge", (ntupleRawTree.b5_Muon_charge));
   outTree->Branch("b5_Muon_jetIdx", (ntupleRawTree.b5_Muon_jetIdx));
   outTree->Branch("b5_Muon_nStations", (ntupleRawTree.b5_Muon_nStations));
   outTree->Branch("b5_Muon_nTrackerLayers", (ntupleRawTree.b5_Muon_nTrackerLayers));
   outTree->Branch("b5_Muon_pdgId", (ntupleRawTree.b5_Muon_pdgId));
   outTree->Branch("b5_Muon_tightCharge", (ntupleRawTree.b5_Muon_tightCharge));
   outTree->Branch("b5_Muon_fsrPhotonIdx", (ntupleRawTree.b5_Muon_fsrPhotonIdx));
   outTree->Branch("b5_Muon_highPtId", (ntupleRawTree.b5_Muon_highPtId));
   outTree->Branch("b5_Muon_highPurity", (ntupleRawTree.b5_Muon_highPurity));
   outTree->Branch("b5_Muon_inTimeMuon", (ntupleRawTree.b5_Muon_inTimeMuon));
   outTree->Branch("b5_Muon_isGlobal", (ntupleRawTree.b5_Muon_isGlobal));
   outTree->Branch("b5_Muon_isPFcand", (ntupleRawTree.b5_Muon_isPFcand));
   outTree->Branch("b5_Muon_isTracker", (ntupleRawTree.b5_Muon_isTracker));
   outTree->Branch("b5_Muon_jetNDauCharged", (ntupleRawTree.b5_Muon_jetNDauCharged));
   outTree->Branch("b5_Muon_looseId", (ntupleRawTree.b5_Muon_looseId));
   outTree->Branch("b5_Muon_mediumId", (ntupleRawTree.b5_Muon_mediumId));
   outTree->Branch("b5_Muon_mediumPromptId", (ntupleRawTree.b5_Muon_mediumPromptId));
   outTree->Branch("b5_Muon_miniIsoId", (ntupleRawTree.b5_Muon_miniIsoId));
   outTree->Branch("b5_Muon_multiIsoId", (ntupleRawTree.b5_Muon_multiIsoId));
   outTree->Branch("b5_Muon_mvaId", (ntupleRawTree.b5_Muon_mvaId));
   outTree->Branch("b5_Muon_mvaLowPtId", (ntupleRawTree.b5_Muon_mvaLowPtId));
   outTree->Branch("b5_Muon_pfIsoId", (ntupleRawTree.b5_Muon_pfIsoId));
   outTree->Branch("b5_Muon_puppiIsoId", (ntupleRawTree.b5_Muon_puppiIsoId));
   outTree->Branch("b5_Muon_softId", (ntupleRawTree.b5_Muon_softId));
   outTree->Branch("b5_Muon_softMvaId", (ntupleRawTree.b5_Muon_softMvaId));
   outTree->Branch("b5_Muon_tightId", (ntupleRawTree.b5_Muon_tightId));
   outTree->Branch("b5_Muon_tkIsoId", (ntupleRawTree.b5_Muon_tkIsoId));
   outTree->Branch("b5_SV_chi2", (ntupleRawTree.b5_SV_chi2));
   outTree->Branch("b5_SV_eta", (ntupleRawTree.b5_SV_eta));
   outTree->Branch("b5_SV_mass", (ntupleRawTree.b5_SV_mass));
   outTree->Branch("b5_SV_ndof", (ntupleRawTree.b5_SV_ndof));
   outTree->Branch("b5_SV_phi", (ntupleRawTree.b5_SV_phi));
   outTree->Branch("b5_SV_pt", (ntupleRawTree.b5_SV_pt));
   outTree->Branch("b5_SV_x", (ntupleRawTree.b5_SV_x));
   outTree->Branch("b5_SV_y", (ntupleRawTree.b5_SV_y));
   outTree->Branch("b5_SV_z", (ntupleRawTree.b5_SV_z));
   outTree->Branch("b5_SV_ntracks", (ntupleRawTree.b5_SV_ntracks));
   outTree->Branch("b5_Muon_triggerIdLoose", (ntupleRawTree.b5_Muon_triggerIdLoose));
   outTree->Branch("b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &(ntupleRawTree.b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4));
   outTree->Branch("b5_L1_DoubleMu0er1p5_SQ", &(ntupleRawTree.b5_L1_DoubleMu0er1p5_SQ));
   outTree->Branch("b5_L1_DoubleMu0er1p5_SQ_OS", &(ntupleRawTree.b5_L1_DoubleMu0er1p5_SQ_OS));
   outTree->Branch("b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &(ntupleRawTree.b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4));
   outTree->Branch("b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4", &(ntupleRawTree.b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4));
   outTree->Branch("b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &(ntupleRawTree.b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4));
   outTree->Branch("b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4", &(ntupleRawTree.b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4));
   outTree->Branch("b5_L1_DoubleMu4_SQ_OS" , &(ntupleRawTree.b5_L1_DoubleMu4_SQ_OS));
   outTree->Branch("b5_L1_DoubleMu4_SQ_OS_dR_Max1p2" , &(ntupleRawTree.b5_L1_DoubleMu4_SQ_OS_dR_Max1p2));
   outTree->Branch("b5_HLT_DoubleMu4_3_Bs" , &(ntupleRawTree.b5_HLT_DoubleMu4_3_Bs));
   outTree->Branch("b5_HLT_DoubleMu4_3_Jpsi" , &(ntupleRawTree.b5_HLT_DoubleMu4_3_Jpsi));
}

void BMMGAnalysis::bookHistograms()
{
     const int   pTBins=50;
     const float pTmin(0.0);
     const float pTmax(50.0);

     const int   etaBins=70;
     const float etamin(-3.5);
     const float etamax(3.5);

     const int   massNBins=20*100;
     const float massMin(0.0);
     const float massMax(20.0);
     
     const int   deltaRNBins=30;
     const float deltaRMin(0.0);
     const float deltaRMax(3.0);
    
    //All Mu Hist
     th1fStore["mu_Pt"]     = new TH1F("mu_Pt","Pt of mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["mu_Eta"]     = new TH1F("mu_Eta","Eta of mu Candidate", etaBins , etamin  , etamax  );
     th1fStore["mu_Charge"]     = new TH1F("mu_charge","charge of mu Candidate", 3 , -1.5  , 1.5  );
     th1fStore["mu_dz"]     = new TH1F("mu_dz","dz of mu Candidate", 60*4 , -30.0  , 30.0  );
     th1fStore["mu_mva"]     = new TH1F("mu_mva","bmm5 mva score of mu Candidate", 60*4 , -1.5  , 1.5  );
    
    //All SC Hists
     th1fStore["sc_Pt"]     = new TH1F("sc_Pt","Pt of sc Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["sc_Eta"]     = new TH1F("sc_Eta","Eta of sc Candidate", etaBins , etamin  , etamax  );
     th1fStore["sc_mva"]     = new TH1F("sc_mva","MVA score of sc Candidate", 80 , -2.0  , 2.0  );

    //All Photon Hists
     th1fStore["photon_Pt"]     = new TH1F("photon_Pt","Pt of photon Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["photon_Eta"]     = new TH1F("photon_Eta","Eta of photon Candidate", etaBins , etamin  , etamax  );
     th1fStore["photon_mva"]     = new TH1F("photon_mva","MVA photonore of photon Candidate", 80 , -2.0  , 2.0  );

    //All Dimu Hists
     th1fStore["dimu_nDimuons"]   = new TH1F("dimu_nDimuons"," number of dimuons in events", 20 , -0.5  , 19.5  );
     th1fStore["dimu_mumuDeltaR"]   = new TH1F("dimu_mumuDeltaR","#Delta(#mu,#mu)", deltaRNBins , deltaRMin  , deltaRMax  );
     th1fStore["dimu_muLeadPt"]     = new TH1F("dimu_muLeadPt","Pt of lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["dimu_muSubLeadPt"]  = new TH1F("dimu_muSubLeadPt","Pt of sub-lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["dimu_muLeadIsolation"]  = new TH1F("dimu_muLeadIsolation","I(#mu_{1}) for #mu#mu Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["dimu_muSubLeadIsolation"]  = new TH1F("dimu_muSubLeadIsolation","I(#mu_{2}) for #mu#mu Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["dimu_mass"]         = new TH1F("dimu_mass","InvMass of mumu Candidate", massNBins , massMin , massMax  );
     th1fStore["dimu_dxy"]          = new TH1F("dimu_dxy","dxy of dimu  Candidate", 20 ,0.0,2.50  );
     th1fStore["dimu_vtxProba"]     = new TH1F("dimu_vtxProbay","vtx Proba. of dimu  Candidate", 20 ,0.0,2.50  );
     th1fStore["dimu_alpha"]        = new TH1F("dimu_alpha","#alpha of dimu  Candidate", 20 ,0.0,3.50  );
     th1fStore["dimu_cosalpha"]        = new TH1F("dimu_cosalpha","cos(#alpha) of dimu  Candidate", 22 ,-1.10 ,1.10  );
     th1fStore["dimu_VertexChi2overNDOF"] = new TH1F("dimu_VertexChi2overNDOF","#Chi^{2}/ndof of #mu#mu vertex for #mu#mu Candidate" ,200 , 0.0 , 50.0  );
     th1fStore["dimu_DCA"]  = new TH1F("dimu_Dca","DCA(#mu#mu) for #mu#mu Candidate" ,100 , 0.0 , 0.050  );
     th1fStore["dimu_VertexOtherTrackDOCA"]  = new TH1F("dimu_VertexOtherTrackDOCA","DCA(#mu#mu vertex, other track) for #mu#mu Candidate" ,100 , 0.0 , 5.0 );
     th1fStore["dimu_L3D"]              = new TH1F("dimu_L3D","l_{3D}(#mu#mu) for #mu#mu Candidate" ,120 , 0.0 , 6.0  );
     th1fStore["dimu_L3Dsignificance"]  = new TH1F("dimu_L3Dsignificance","l_{3D}(#mu#mu)/#sigma(l_{3D}) for #mu#mu Candidate" ,400 , 0.0 , 100  );
     th1fStore["dimu_LXY"]              = new TH1F("dimu_LXY","l_{XY}(#mu#mu) for #mu#mu Candidate" ,50 , 0.0 , 2.0  );
     th1fStore["dimu_LXYsignificance"]  = new TH1F("dimu_LXYsignificance","l_{XY}(#mu#mu)/#sigma(l_{XY}) for #mu#mu Candidate" ,240 , 0.0 , 120.0 );
     th1fStore["dimu_IPwrtPV"]              = new TH1F("dimu_IPwrtPV","IP(#mu#mu) wrt. PV for #mu#mu Candidate" , 80 , 0.0 , 0.40  );
     th1fStore["dimu_IPwrtPVsignificance"]  = new TH1F("dimu_IPwrtPVsignificance","IP(#mu#mu) wrt. PV for #mu#mu Candidate" ,160 , 0.0 , 40.0 );
     th1fStore["dimu_LIPwrtPV"]              = new TH1F("dimu_LIPwrtPV","Longitudinal IP(#mu#mu) wrt. PV for #mu#mu Candidate" , 200 , -2.0 , 2.0  );
     th1fStore["dimu_LIPwrtPVsignificance"]  = new TH1F("dimu_LIPwrtPVsignificance","Longitudinal IP(#mu#mu) wrt. PV for #mu#mu Candidate" , 400 , -50.0 , 50.0 );
     th1fStore["dimu_NTrakClose"]       = new TH1F("dimu_NTrakClose","#Tracks doca(trk,sv) < 0.3 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_NTrakCloseSig1"]   = new TH1F("dimu_NTrakCloseSig1","#Tracks d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 1 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_NTrakCloseSig2"]   = new TH1F("dimu_NTrakCloseSig2","#Tracks d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 2 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_NTrakCloseSig3"]   = new TH1F("dimu_NTrakCloseSig3","#Tracks d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 3 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_Isolation"]  = new TH1F("dimu_Isolation","I(#mu#mu) for #mu#mu Candidate" ,120 , -0.10 , 1.1  );

    // Selected Dimu Vertex Hists
     th1fStore["dimuPass_mumuDeltaR"]   = new TH1F("mumuDeltaR","#Delta(#mu,#mu)", deltaRNBins , deltaRMin  , deltaRMax  );
     th1fStore["dimuPass_muLeadPt"]     = new TH1F("dimuPass_muLeadPt","Pt of lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["dimuPass_muSubLeadPt"]  = new TH1F("dimuPass_muSubLeadPt","Pt of sub-lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["dimuPass_muLeadIsolation"]  = new TH1F("dimuPass_muLeadIsolation","I(#mu_{1}) for #mu#mu#gamma Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["dimuPass_muSubLeadIsolation"]  = new TH1F("dimuPass_muSubLeadIsolation","I(#mu_{2}) for #mu#mu#gamma Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["dimuPass_mass"]         = new TH1F("dimuPass_mass","InvMass of mumu Candidate", massNBins , massMin , massMax  );
     th1fStore["dimuPass_dxy"]          = new TH1F("dimuPass_dxy","dxy of dimuPass  Candidate", 20 ,0.0,2.50  );
     th1fStore["dimuPass_vtxProba"]     = new TH1F("dimuPass_vtxProbay","vtx Proba. of dimuPass  Candidate", 20 ,0.0,2.50  );
     th1fStore["dimuPass_alpha"]        = new TH1F("dimuPass_alpha","#alpha of dimu  Candidate", 20 ,0.0,3.50  );
     th1fStore["dimuPass_cosalpha"]        = new TH1F("dimuPass_cosalpha","cos(#alpha) of dimu  Candidate", 22 ,-1.10 ,1.10  );
     th1fStore["dimuPass_VertexChi2overNDOF"] = new TH1F("dimuPass_VertexChi2overNDOF","#Chi^{2}/ndof of #mu#mu vertex for #mu#mu#gamma Candidate" ,200 , 0.0 , 50.0  );
     th1fStore["dimuPass_DCA"]  = new TH1F("dimuPass_Dca","DCA(#mu#mu) for #mu#mu#gamma Candidate" ,100 , 0.0 , 0.050  );
     th1fStore["dimuPass_VertexOtherTrackDOCA"]  = new TH1F("dimuPass_VertexOtherTrackDOCA","DCA(#mu#mu vertex, other track) for #mu#mu#gamma Candidate" ,100 , 0.0 , 5.0 );
     th1fStore["dimuPass_L3D"]              = new TH1F("dimuPass_L3D","l_{3D}(#mu#mu) for #mu#mu#gamma Candidate" ,120 , 0.0 , 6.0  );
     th1fStore["dimuPass_L3Dsignificance"]  = new TH1F("dimuPass_L3Dsignificance","l_{3D}(#mu#mu)/#sigma(l_{3D}) for #mu#mu#gamma Candidate" ,50 , 0.0 , 2.5  );
     th1fStore["dimuPass_LXY"]              = new TH1F("dimuPass_LXY","l_{XY}(#mu#mu) for #mu#mu#gamma Candidate" ,50 , 0.0 , 2.0  );
     th1fStore["dimuPass_LXYsignificance"]  = new TH1F("dimuPass_LXYsignificance","l_{XY}(#mu#mu)/#sigma(l_{XY}) for #mu#mu#gamma Candidate" ,240 , 0.0 , 120.0 );
     th1fStore["dimuPass_IPwrtPV"]              = new TH1F("dimuPass_IPwrtPV","IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" , 80 , 0.0 , 0.40  );
     th1fStore["dimuPass_IPwrtPVsignificance"]  = new TH1F("dimuPass_IPwrtPVsignificance","IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" ,160 , 0.0 , 40.0 );
     th1fStore["dimuPass_LIPwrtPV"]              = new TH1F("dimuPass_LIPwrtPV","Longitudinal IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" , 200 , -2.0 , 2.0  );
     th1fStore["dimuPass_LIPwrtPVsignificance"]  = new TH1F("dimuPass_LIPwrtPVsignificance","Longitudinal IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" , 400 , -50.0 , 50.0 );
     th1fStore["dimuPass_NTrakClose"]       = new TH1F("dimuPass_NTrakClose","#Tracks doca(trk,sv) < 0.3 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimuPass_NTrakCloseSig1"]= new TH1F("dimuPass_NTrakCloseSig1","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 1 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimuPass_NTrakCloseSig2"]= new TH1F("dimuPass_NTrakCloseSig2","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 2 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimuPass_NTrakCloseSig3"]= new TH1F("dimuPass_NTrakCloseSig3","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 3 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimuPass_Isolation"]  = new TH1F("dimuPass_Isolation","I(#mu#mu) for #mu#mu#gamma Candidate" ,120 , -0.10 , 1.1  );
 
    // Photon Environment Hists after Dimuon Selection
    th1fStore["dimu_photonMultiplicityDrMax0p5"] = new TH1F("dimu_photonMultiplicityDrMax0p5","Multiplicity of photon around #mu#mu , #Delta R < 0.5",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax0p7"] = new TH1F("dimu_photonMultiplicityDrMax0p7","Multiplicity of photon around #mu#mu , #Delta R < 0.7",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax0p9"] = new TH1F("dimu_photonMultiplicityDrMax0p9","Multiplicity of photon around #mu#mu , #Delta R < 0.9",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax1p1"] = new TH1F("dimu_photonMultiplicityDrMax1p1","Multiplicity of photon around #mu#mu , #Delta R < 1.1",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax1p3"] = new TH1F("dimu_photonMultiplicityDrMax1p3","Multiplicity of photon around #mu#mu , #Delta R < 1.3",15, -0.5,14.5);
 
    // BMMG HISTS
     th1fStore["bmmg_photonPt"]         = new TH1F("bmmg_photonPt","Pt of photon  Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["bmmg_mmMass"]           = new TH1F("bmmg_mmMass","InvMass of mumu Candidate", massNBins , massMin , massMax  );
     th1fStore["bmmg_mmgMass"]          = new TH1F("bmmg_mmgMass","InvMass of mmgCandidate" , massNBins , massMin , massMax  );
     th1fStore["bmmg_mmgPt"]            = new TH1F("bmmg_mmgPt","p_{T} of #mu#mu#gamma Candidate" , pTBins , pTmin , pTmax  );
     th1fStore["bmmg_beta"]             = new TH1F("bmmg_beta","#beta(#mu#mu) of #mu#mu#gamma  Candidate", 20 ,0.0,3.50  );
     th1fStore["bmmg_cosbeta"]          = new TH1F("bmmg_cosbeta","cos(#beta(#mu#mu)) of #mu#mu#gamma  Candidate", 22 ,-1.10 ,1.10  );
     th1fStore["bmmg_dimuGammaDeltaR"] = new TH1F("bmmg_dimuGammaDeltaR","#Delta(#mu#mu,#gamma) for #mu#mu#gamma", deltaRNBins , deltaRMin  , deltaRMax  );
     
     th1fStore["bmmg_mumuDeltaR"]   = new TH1F("bmmg_mumuDeltaR","#Delta(#mu,#mu) for #mu#mu#gamma", deltaRNBins , deltaRMin  , deltaRMax  );
     th1fStore["bmmg_muLeadPt"]     = new TH1F("bmmg_muLeadPt","Pt of lead mu Candidate for #mu#mu#gamma", pTBins , pTmin  , pTmax  );
     th1fStore["bmmg_muSubLeadPt"]  = new TH1F("bmmg_muSubLeadPt","Pt of sub-lead #mu#mu#gamma Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["bmmg_muLeadIsolation"]  = new TH1F("bmmg_muLeadIsolation","I(#mu_{1}) for #mu#mu#gamma Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["bmmg_muSubLeadIsolation"]  = new TH1F("bmmg_muSubLeadIsolation","I(#mu_{2}) for #mu#mu#gamma Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["bmmg_mass"]         = new TH1F("bmmg_mass","InvMass of #mu#mu#gamma Candidate", massNBins , massMin , massMax  );
     th1fStore["bmmg_dxy"]          = new TH1F("bmmg_dxy","dxy(#mu#mu)of #mu#mu#gamma  Candidate", 20 ,0.0,2.50  );
     th1fStore["bmmg_vtxProba"]     = new TH1F("bmmg_vtxProbay","vtx Proba. of bmmg  Candidate", 20 ,0.0,2.50  );
     th1fStore["bmmg_alpha"]        = new TH1F("bmmg_alpha","#alpha(#mu#mu) of #mu#mu#gamma  Candidate", 20 ,0.0,3.50  );
     th1fStore["bmmg_cosalpha"]        = new TH1F("bmmg_cosalpha","cos(#alpha(#mu#mu)) of #mu#mu#gamma  Candidate", 22 ,-1.10 ,1.10  );
     th1fStore["bmmg_VertexChi2overNDOF"] = new TH1F("bmmg_VertexChi2overNDOF","#Chi^{2}/ndof(#mu#mu) vertex for #mu#mu#gamma Candidate" ,200 , 0.0 , 50.0  );
     th1fStore["bmmg_DCA"]  = new TH1F("bmmg_Dca","DCA(#mu#mu) for #mu#mu#gamma Candidate" ,100 , 0.0 , 0.050  );
     th1fStore["bmmg_VertexOtherTrackDOCA"]  = new TH1F("bmmg_VertexOtherTrackDOCA","DCA(#mu#mu vertex, other track) for #mu#mu#gamma Candidate" ,100 , 0.0 , 5.0 );
     th1fStore["bmmg_L3D"]              = new TH1F("bmmg_L3D","l_{3D}(#mu#mu) for #mu#mu#gamma Candidate" ,120 , 0.0 , 6.0  );
     th1fStore["bmmg_L3Dsignificance"]  = new TH1F("bmmg_L3Dsignificance","l_{3D}(#mu#mu)/#sigma(l_{3D}) for #mu#mu#gamma Candidate" ,50 , 0.0 , 2.5  );
     th1fStore["bmmg_LXY"]              = new TH1F("bmmg_LXY","l_{XY}(#mu#mu) for #mu#mu#gamma Candidate" ,50 , 0.0 , 2.0  );
     th1fStore["bmmg_LXYsignificance"]  = new TH1F("bmmg_LXYsignificance","l_{XY}(#mu#mu)/#sigma(l_{XY}) for #mu#mu#gamma Candidate" ,240 , 0.0 , 120.0 );
     th1fStore["bmmg_IPwrtPV"]              = new TH1F("bmmg_IPwrtPV","IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" , 80 , 0.0 , 0.40  );
     th1fStore["bmmg_IPwrtPVsignificance"]  = new TH1F("bmmg_IPwrtPVsignificance","IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" ,160 , 0.0 , 40.0 );
     th1fStore["bmmg_LIPwrtPV"]              = new TH1F("bmmg_LIPwrtPV","Longitudinal IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" , 200 , -2.0 , 2.0  );
     th1fStore["bmmg_LIPwrtPVsignificance"]  = new TH1F("bmmg_LIPwrtPVsignificance","Longitudinal IP(#mu#mu) wrt. PV for #mu#mu#gamma Candidate" , 400 , -50.0 , 50.0 );
     th1fStore["bmmg_NTrakClose"]       = new TH1F("bmmg_NTrakClose","#Tracks doca(trk,sv) < 0.3 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_NTrakCloseSig1"]= new TH1F("bmmg_NTrakCloseSig1","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 1 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_NTrakCloseSig2"]= new TH1F("bmmg_NTrakCloseSig2","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 2 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_NTrakCloseSig3"]= new TH1F("bmmg_NTrakCloseSig3","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 3 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_Isolation"]  = new TH1F("bmmg_Isolation","I(#mu#mu) for #mu#mu#gamma Candidate" ,120 , -0.10 , 1.1  );
 
    // Global Event Hists
    th1fStore["dimuPass_bmmgCandidateMultiplicity"]   = new TH1F("dimuPass_bmmgCandidateMultiplicity","Multiplicity of BMMG Candidates Per dimuon",15, -0.5,14.5);
    th1fStore["dimuPass_Multiplicity"]       = new TH1F("dimuPass_Multiplicity","Multiplicity of dimu Candidates in an event",15, -0.5,14.5);
    th1fStore["bmmg_Multiplicity"]       = new TH1F("bmmg_Multiplicity","Multiplicity of dimu Candidates in an event",15, -0.5,14.5);

    // MC : Gen Match Hists 
     th1fStore["gen_mumDeltaR"]   = new TH1F("gen_mumDeltaR","#Delta(#mu^-,#mu_{Gen}^-)", 200 , 0.0  , 2.0  );
     th1fStore["gen_mupDeltaR"]   = new TH1F("gen_mupDeltaR","#Delta(#mu^+,#mu_{Gen}^+)", 200 , 0.0  , 2.0  );
     th1fStore["gen_phoDeltaR"]   = new TH1F("gen_phoDeltaR","#Delta(#gamma,#gamma_{Gen})", 200 , 0.0  , 2.0  );
}



void BMMGAnalysis::fill_muonHists(Int_t idx)
{   
    if(idx < 0)
    for(int i=0;i<ntupleRawTree.b5_nMuon ; i++)
    {
        th1fStore["mu_Pt"]     ->Fill(ntupleRawTree.b5_Muon_pt[i]);
        th1fStore["mu_Eta"]    ->Fill(ntupleRawTree.b5_Muon_eta[i]);
        th1fStore["mu_Charge"] ->Fill(ntupleRawTree.b5_Muon_charge[i]);
        th1fStore["mu_dz"]     ->Fill(ntupleRawTree.b5_Muon_dz[i]);
        th1fStore["mu_mva"]    ->Fill(ntupleRawTree.b5_MuonId_newSoftMuonMva[i] );
    }
    else if(idx < ntupleRawTree.b5_nMuon)
    {
        
        th1fStore["mu_Pt"]     ->Fill(ntupleRawTree.b5_Muon_pt[idx]);
        th1fStore["mu_Eta"]    ->Fill(ntupleRawTree.b5_Muon_eta[idx]);
        th1fStore["mu_Charge"] ->Fill(ntupleRawTree.b5_Muon_charge[idx]);
        th1fStore["mu_dz"]     ->Fill(ntupleRawTree.b5_Muon_dz[idx]);
        th1fStore["mu_mva"]    ->Fill(ntupleRawTree.b5_MuonId_newSoftMuonMva[idx] );
    }
    
}

void BMMGAnalysis::fill_scHists(Int_t idx  )
{
    if (idx <0 )
    for(int i=0;i<ntupleRawTree.bG_nSC ; i++)
    {
        th1fStore["sc_Pt"]    ->Fill(ntupleRawTree.bG_scEt[i]);
        th1fStore["sc_Eta"]   ->Fill(ntupleRawTree.bG_scEta[i]); 
        th1fStore["sc_mva"]   ->Fill(storageArrayDouble[ i + candidateMapInt["scPhotonMVAScore"]  ]); 
    }
    else if (idx < ntupleRawTree.bG_nSC)
    {
        th1fStore["sc_Pt"]    ->Fill(ntupleRawTree.bG_scEt[idx]);
        th1fStore["sc_Eta"]   ->Fill(ntupleRawTree.bG_scEta[idx]); 
        th1fStore["sc_mva"]   ->Fill(storageArrayDouble[ idx + candidateMapInt["scPhotonMVAScore"]  ]); 

    }
}

void BMMGAnalysis::fill_photonHists(Int_t idx)
{
    if(idx <0)
     for(int i=0;i<ntupleRawTree.bG_nSC ; i++)
        {
            if(doPhotonSelection(i)) continue;
            
            th1fStore["photon_Pt"]    ->Fill(ntupleRawTree.bG_scEt[i]);
            th1fStore["photon_Eta"]   ->Fill(ntupleRawTree.bG_scEta[i]); 
            th1fStore["photon_mva"]   ->Fill(storageArrayDouble[ i + candidateMapInt["scPhotonMVAScore"]  ]); 
        }
     else if(idx < ntupleRawTree.bG_nSC )   
        {
            if(doPhotonSelection(idx)< 1) {
            
            th1fStore["photon_Pt"]    ->Fill(ntupleRawTree.bG_scEt[idx]);
            th1fStore["photon_Eta"]   ->Fill(ntupleRawTree.bG_scEta[idx]); 
            th1fStore["photon_mva"]   ->Fill(storageArrayDouble[ idx + candidateMapInt["scPhotonMVAScore"]  ]); 
            }
        }
}
void BMMGAnalysis::fill_dimuonHists(Int_t idx)
{
    th1fStore["dimu_nDimuons"]->Fill(ntupleRawTree.b5_nmm);
  for(Int_t mumuIdx=0;mumuIdx<ntupleRawTree.b5_nmm ; mumuIdx++)
   {
        if(idx >= 0  and mumuIdx != idx ) continue;
     auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx],
                     ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
     th1fStore["dimu_mumuDeltaR"]  ->Fill(dr);
     th1fStore["dimu_muLeadPt"]    ->Fill(ntupleRawTree.b5_mm_kin_mu1pt[mumuIdx]);
     th1fStore["dimu_muSubLeadPt"] ->Fill(ntupleRawTree.b5_mm_kin_mu2pt[mumuIdx]);
     th1fStore["dimu_muLeadIsolation"]    ->Fill(ntupleRawTree.b5_mm_m1iso[mumuIdx]);
     th1fStore["dimu_muSubLeadIsolation"] ->Fill(ntupleRawTree.b5_mm_m2iso[mumuIdx]);
     th1fStore["dimu_mass"]        ->Fill(ntupleRawTree.b5_mm_kin_mass[mumuIdx]);
     th1fStore["dimu_dxy"]         ->Fill(ntupleRawTree.b5_mm_kin_lxy[mumuIdx]);
     th1fStore["dimu_vtxProba"]    ->Fill(ntupleRawTree.b5_mm_kin_vtx_prob[mumuIdx]);
     th1fStore["dimu_alpha"]       ->Fill(ntupleRawTree.b5_mm_kin_alpha[mumuIdx]);
     th1fStore["dimu_cosalpha"]        ->Fill(cos(ntupleRawTree.b5_mm_kin_alpha[mumuIdx]));
     th1fStore["dimu_VertexChi2overNDOF"] ->Fill(ntupleRawTree.b5_mm_kin_vtx_chi2dof[mumuIdx]);
     th1fStore["dimu_DCA"]  ->Fill(ntupleRawTree.b5_mm_doca[mumuIdx]);
     th1fStore["dimu_VertexOtherTrackDOCA"]  ->Fill(ntupleRawTree.b5_mm_docatrk[mumuIdx]);
     th1fStore["dimu_L3D"]              ->Fill(ntupleRawTree.b5_mm_kin_l3d[mumuIdx]);
     th1fStore["dimu_L3Dsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_sl3d[mumuIdx]);
     th1fStore["dimu_LXY"]              ->Fill(ntupleRawTree.b5_mm_kin_lxy[mumuIdx]);
     th1fStore["dimu_LXYsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_slxy[mumuIdx]);
     th1fStore["dimu_IPwrtPV"]              ->Fill(ntupleRawTree.b5_mm_kin_pvip[mumuIdx]);
     th1fStore["dimu_IPwrtPVsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_spvip[mumuIdx]);
     th1fStore["dimu_LIPwrtPV"]              ->Fill(ntupleRawTree.b5_mm_kin_pvlip[mumuIdx]);
     th1fStore["dimu_LIPwrtPVsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_pvlipSig[mumuIdx]);
     th1fStore["dimu_NTrakClose"]       ->Fill(ntupleRawTree.b5_mm_closetrk[mumuIdx]);
     th1fStore["dimu_NTrakCloseSig1"]   ->Fill(ntupleRawTree.b5_mm_closetrks1[mumuIdx]);
     th1fStore["dimu_NTrakCloseSig2"]   ->Fill(ntupleRawTree.b5_mm_closetrks2[mumuIdx]);
     th1fStore["dimu_NTrakCloseSig3"]   ->Fill(ntupleRawTree.b5_mm_closetrks3[mumuIdx]);
     th1fStore["dimu_Isolation"]  ->Fill(ntupleRawTree.b5_mm_iso[mumuIdx]);
   }
}
void BMMGAnalysis::fill_dimuonPassHists(Int_t mumuIdx)
{
    
     auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx],
                     ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
     th1fStore["dimuPass_mumuDeltaR"]  ->Fill(dr);
     th1fStore["dimuPass_muLeadPt"]    ->Fill(ntupleRawTree.b5_mm_kin_mu1pt[mumuIdx]);
     th1fStore["dimuPass_muSubLeadPt"] ->Fill(ntupleRawTree.b5_mm_kin_mu2pt[mumuIdx]);
     th1fStore["dimuPass_muLeadIsolation"]    ->Fill(ntupleRawTree.b5_mm_m1iso[mumuIdx]);
     th1fStore["dimuPass_muSubLeadIsolation"] ->Fill(ntupleRawTree.b5_mm_m2iso[mumuIdx]);
     th1fStore["dimuPass_mass"]        ->Fill(ntupleRawTree.b5_mm_kin_mass[mumuIdx]);
     th1fStore["dimuPass_dxy"]         ->Fill(ntupleRawTree.b5_mm_kin_lxy[mumuIdx]);
     th1fStore["dimuPass_vtxProba"]    ->Fill(ntupleRawTree.b5_mm_kin_vtx_prob[mumuIdx]);
     th1fStore["dimuPass_alpha"]       ->Fill(ntupleRawTree.b5_mm_kin_alpha[mumuIdx]);
     th1fStore["dimuPass_cosalpha"]        ->Fill(cos(ntupleRawTree.b5_mm_kin_alpha[mumuIdx]));
     th1fStore["dimuPass_VertexChi2overNDOF"] ->Fill(ntupleRawTree.b5_mm_kin_vtx_chi2dof[mumuIdx]);
     th1fStore["dimuPass_DCA"]  ->Fill(ntupleRawTree.b5_mm_doca[mumuIdx]);
     th1fStore["dimuPass_VertexOtherTrackDOCA"]  ->Fill(ntupleRawTree.b5_mm_docatrk[mumuIdx]);
     th1fStore["dimuPass_L3D"]              ->Fill(ntupleRawTree.b5_mm_kin_l3d[mumuIdx]);
     th1fStore["dimuPass_L3Dsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_sl3d[mumuIdx]);
     th1fStore["dimuPass_LXY"]              ->Fill(ntupleRawTree.b5_mm_kin_lxy[mumuIdx]);
     th1fStore["dimuPass_LXYsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_slxy[mumuIdx]);
     th1fStore["dimuPass_IPwrtPV"]              ->Fill(ntupleRawTree.b5_mm_kin_pvip[mumuIdx]);
     th1fStore["dimuPass_IPwrtPVsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_spvip[mumuIdx]);
     th1fStore["dimuPass_LIPwrtPV"]              ->Fill(ntupleRawTree.b5_mm_kin_pvlip[mumuIdx]);
     th1fStore["dimuPass_LIPwrtPVsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_pvlipSig[mumuIdx]);
     th1fStore["dimuPass_NTrakClose"]       ->Fill(ntupleRawTree.b5_mm_closetrk[mumuIdx]);
     th1fStore["dimuPass_NTrakCloseSig1"]   ->Fill(ntupleRawTree.b5_mm_closetrks1[mumuIdx]);
     th1fStore["dimuPass_NTrakCloseSig2"]   ->Fill(ntupleRawTree.b5_mm_closetrks2[mumuIdx]);
     th1fStore["dimuPass_NTrakCloseSig3"]   ->Fill(ntupleRawTree.b5_mm_closetrks3[mumuIdx]);
     th1fStore["dimuPass_Isolation"]  ->Fill(ntupleRawTree.b5_mm_iso[mumuIdx]);

}

void BMMGAnalysis::fill_dimuonEnvironmentHists(Int_t mumuIdx)
{
     auto eta = ntupleRawTree.b5_mm_kin_eta[mumuIdx];
     auto phi = ntupleRawTree.b5_mm_kin_phi[mumuIdx];
     auto dr =0.0;   
     int array[5]={0,0,0,0,0};
     for(int i=0;i<ntupleRawTree.bG_nSC ; i++)
        {
            if(doPhotonSelection(i)) continue;
            
            dr = getDR(eta,phi , ntupleRawTree.bG_scEta[i], ntupleRawTree.bG_scPhi[i]);
            
            if(dr < 0.5 ) array[0]++;
            if(dr < 0.7 ) array[1]++;
            if(dr < 0.9 ) array[2]++;
            if(dr < 1.1 ) array[3]++;
            if(dr < 1.3 ) array[4]++;
        }
    
    th1fStore["dimu_photonMultiplicityDrMax0p5"]->Fill(array[0]); 
    th1fStore["dimu_photonMultiplicityDrMax0p7"]->Fill(array[1]); 
    th1fStore["dimu_photonMultiplicityDrMax0p9"]->Fill(array[2]); 
    th1fStore["dimu_photonMultiplicityDrMax1p1"]->Fill(array[3]); 
    th1fStore["dimu_photonMultiplicityDrMax1p3"]->Fill(array[4]); 


}


void BMMGAnalysis::fill_bmmgHists(  TLorentzVector &bmmgLV , Int_t mumuIdx, Int_t phoSCIdx)
{
    auto dr=getDR(ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx],
                  ntupleRawTree.b5_mm_kin_eta[mumuIdx],ntupleRawTree.b5_mm_kin_phi[mumuIdx]);
    th1fStore["bmmg_dimuGammaDeltaR"]  ->Fill(dr);
    
    dr=getDR(ntupleRawTree.b5_mm_mu1_eta[mumuIdx],ntupleRawTree.b5_mm_mu1_phi[mumuIdx],
                  ntupleRawTree.b5_mm_mu2_eta[mumuIdx],ntupleRawTree.b5_mm_mu2_phi[mumuIdx]);
    th1fStore["bmmg_mumuDeltaR"]       ->Fill(dr);

    th1fStore["bmmg_photonPt"]         ->Fill(ntupleRawTree.bG_scEt[phoSCIdx]);
    th1fStore["bmmg_mmMass"]           ->Fill(ntupleRawTree.b5_mm_kin_mass[mumuIdx]);
    th1fStore["bmmg_mmgMass"]          ->Fill(bmmgLV.M());
    
    th1fStore["bmmg_mumuDeltaR"]  ->Fill(dr);
    th1fStore["bmmg_muLeadPt"]    ->Fill(ntupleRawTree.b5_mm_kin_mu1pt[mumuIdx]);
    th1fStore["bmmg_muSubLeadPt"] ->Fill(ntupleRawTree.b5_mm_kin_mu2pt[mumuIdx]);
    th1fStore["bmmg_muLeadIsolation"]    ->Fill(ntupleRawTree.b5_mm_m1iso[mumuIdx]);
    th1fStore["bmmg_muSubLeadIsolation"] ->Fill(ntupleRawTree.b5_mm_m2iso[mumuIdx]);
    th1fStore["bmmg_mass"]        ->Fill(ntupleRawTree.b5_mm_kin_mass[mumuIdx]);
    th1fStore["bmmg_dxy"]         ->Fill(ntupleRawTree.b5_mm_kin_lxy[mumuIdx]);
    th1fStore["bmmg_vtxProba"]    ->Fill(ntupleRawTree.b5_mm_kin_vtx_prob[mumuIdx]);
    th1fStore["bmmg_alpha"]       ->Fill(ntupleRawTree.b5_mm_kin_alpha[mumuIdx]);
    th1fStore["bmmg_cosalpha"]        ->Fill(cos(ntupleRawTree.b5_mm_kin_alpha[mumuIdx]));
    th1fStore["bmmg_VertexChi2overNDOF"] ->Fill(ntupleRawTree.b5_mm_kin_vtx_chi2dof[mumuIdx]);
    th1fStore["bmmg_DCA"]  ->Fill(ntupleRawTree.b5_mm_doca[mumuIdx]);
    th1fStore["bmmg_VertexOtherTrackDOCA"]  ->Fill(ntupleRawTree.b5_mm_docatrk[mumuIdx]);
    th1fStore["bmmg_L3D"]              ->Fill(ntupleRawTree.b5_mm_kin_l3d[mumuIdx]);
    th1fStore["bmmg_L3Dsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_sl3d[mumuIdx]);
    th1fStore["bmmg_LXY"]              ->Fill(ntupleRawTree.b5_mm_kin_lxy[mumuIdx]);
    th1fStore["bmmg_LXYsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_slxy[mumuIdx]);
    th1fStore["bmmg_IPwrtPV"]              ->Fill(ntupleRawTree.b5_mm_kin_pvip[mumuIdx]);
    th1fStore["bmmg_IPwrtPVsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_spvip[mumuIdx]);
    th1fStore["bmmg_LIPwrtPV"]              ->Fill(ntupleRawTree.b5_mm_kin_pvlip[mumuIdx]);
    th1fStore["bmmg_LIPwrtPVsignificance"]  ->Fill(ntupleRawTree.b5_mm_kin_pvlipSig[mumuIdx]);
    th1fStore["bmmg_NTrakClose"]       ->Fill(ntupleRawTree.b5_mm_closetrk[mumuIdx]);
    th1fStore["bmmg_NTrakCloseSig1"]   ->Fill(ntupleRawTree.b5_mm_closetrks1[mumuIdx]);
    th1fStore["bmmg_NTrakCloseSig2"]   ->Fill(ntupleRawTree.b5_mm_closetrks2[mumuIdx]);
    th1fStore["bmmg_NTrakCloseSig3"]   ->Fill(ntupleRawTree.b5_mm_closetrks3[mumuIdx]);
    th1fStore["bmmg_Isolation"]  ->Fill(ntupleRawTree.b5_mm_iso[mumuIdx]);

   auto pvMatch=getPVMatch(mumuIdx);

//   TODO  : NEED TO USES THE PROPER ERROR FINDING APPROCH FOR BETA // try out implimnetation in saras Util.cc 
   
   svDisplacementVecor.SetCoordinates( 
                                        ntupleRawTree.b5_mm_kin_vtx_x[mumuIdx] - ntupleRawTree.bG_primaryVertex_x[pvMatch],
                                        ntupleRawTree.b5_mm_kin_vtx_y[mumuIdx] - ntupleRawTree.bG_primaryVertex_y[pvMatch],
                                        ntupleRawTree.b5_mm_kin_vtx_z[mumuIdx] - ntupleRawTree.bG_primaryVertex_z[pvMatch]
                                     );
   bmmg3Momentum.SetCoordinates(
                                    bmmgLV.Px(),
                                    bmmgLV.Py(),
                                    bmmgLV.Pz()
                               );
    
   auto cos_beta       =  bmmg3Momentum.Dot(svDisplacementVecor)/sqrt(bmmg3Momentum.Mag2())/sqrt(svDisplacementVecor.Mag2());
   auto beta           =  acos(cos_beta);  
   auto cos_beta_error =  cos_beta*-1;
   

    th1fStore["bmmg_beta"]           ->Fill(beta);
    th1fStore["bmmg_cosbeta"]        ->Fill(cos_beta);

}

void BMMGAnalysis::fill_globalEventHists()
{
    th1fStore["bmmg_Multiplicity"]->Fill(nBMMGCandidates);
    th1fStore["bmmg_Multiplicity"]->Fill(nDiMuCandidates); 
    
    for(int i=0;i<nDiMuCandidates;i++)
    {
        if(i>=NDIMU_MAX) break;
        th1fStore["dimuPass_bmmgCandidateMultiplicity"]->Fill(i+storageArrayInt[candidateMapInt["nBMMGCandidatesPerDimu"]  ] );
    }
}



Int_t BMMGAnalysis::doMuonSelection(Int_t muIdx, bool isLead)
{
    
 //   std::cout<<"\t\t  muIdx : "<<muIdx<<"\n"; 
    if(  ntupleRawTree.b5_Muon_pt[muIdx] < 4.0 )  return 1;
    if(  abs(ntupleRawTree.b5_Muon_eta[muIdx]) > 1.4 )  return 2;
    if(not ntupleRawTree.b5_Muon_isTracker[muIdx]) return 4;
    if(not ntupleRawTree.b5_Muon_isGlobal[muIdx]) return 5;
    if(not ntupleRawTree.b5_Muon_looseId[muIdx]) return 6;
    if(not ntupleRawTree.b5_MuonId_highPurity[muIdx]) return 7;
    // auto muGblId= ntupleRawTree.b5_mm_mu1_index[muIdx];
    // Soft Muon ID   : Cut Based
    
    /*

    if( not ntupleRawTree.Muon_softId[muGblId] ) return 1;
    
    */

    // Soft Muon MVA 
    
    //if( not ntupleRawTree.b5_Muon_softMvaId[muIdx] ) return 2;

    // BMM5 MVA


    if(ntupleRawTree.b5_MuonId_newSoftMuonMva[muIdx] < BDTWorkingPoint ) return 8;

    return 0;

}

Int_t BMMGAnalysis::doDimuonSelection(Int_t mumuIdx)
{
    if( ntupleRawTree.b5_mm_m1iso[mumuIdx]  < 0.80 ) return 1;
    if( ntupleRawTree.b5_mm_m2iso[mumuIdx]  < 0.80 ) return 2;


    diMuLV.SetPtEtaPhiM(     ntupleRawTree.b5_mm_kin_pt[mumuIdx],   \
                             ntupleRawTree.b5_mm_kin_eta[mumuIdx],  \
                             ntupleRawTree.b5_mm_kin_phi[mumuIdx],  \
                             ntupleRawTree.b5_mm_kin_mass[mumuIdx]  );
    
    auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
               ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );

    if (dr > maxMuMuDr ) return 3;
    
    auto diMuMass = diMuLV.M();
    bool selection=false;

    for(int i=0;i<minDimuMass.size();i++)
    {
         if(diMuMass > minDimuMass[i] and diMuMass < maxDimuMass[i]) selection=true;
         if(selection) break;
    }
    
    if( not selection ) return 4;   

    
    return 0;
           
}

Int_t BMMGAnalysis::doVertexSelection(Int_t mumuIdx)
{
  //  if( cos(ntupleRawTree.b5_mm_kin_alpha[mumuIdx])  > 0.80 ) return 1;
    if( ntupleRawTree.b5_mm_kin_vtx_chi2dof[mumuIdx] > 2.20  ) return 2;
    if( ntupleRawTree.b5_mm_kin_sl3d[mumuIdx]  < 13.0  ) return 3;
    if( ntupleRawTree.b5_mm_kin_slxy[mumuIdx]  < 3.0   ) return 4;
 //   if( ntupleRawTree.b5_mm_kin_pvip[mumuIdx]  > 0.008 ) return 5;
 //   if( ntupleRawTree.b5_mm_kin_spvip[mumuIdx] < 2.0   ) return 6;
    if( ntupleRawTree.b5_mm_doca[mumuIdx]  > 0.100 ) return 7;
    if( ntupleRawTree.b5_mm_iso[mumuIdx]   < 0.80  ) return 8;
    if( ntupleRawTree.b5_mm_closetrk[mumuIdx]  >= 2.0 ) return 9;
    if( ntupleRawTree.b5_mm_docatrk[mumuIdx]  > 0.0150 ) return 10;

    return 0;
}

Int_t BMMGAnalysis::doPhotonSelection(Int_t scIdx)
{
    if( ntupleRawTree.bG_scEt[scIdx] < 4.0 )  return 2;

    if( storageArrayDouble[ scIdx + candidateMapInt["scPhotonMVAScore"]  ] < 0.43 ) return 1;
    
    return 0;
}

Int_t BMMGAnalysis::doBMMGSelection(Int_t mumuIdx, Int_t phoSCIdx)
{
   auto dr=getDR(ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx],
                 ntupleRawTree.b5_mm_kin_eta[mumuIdx] , ntupleRawTree.b5_mm_kin_phi[mumuIdx]  );
   
   if(dr > maxDimuPhotonDr ) return 1;

/*   

   auto pvMatch=getPVMatch(mumuIdx);

//   TODO  : NEED TO USES THE PROPER ERROR FINDING APPROCH FOR BETA // try out implimnetation in saras Util.cc 
   
   svDisplacementVecor.SetCoordinates( 
                                        ntupleRawTree.b5_mm_kin_vtx_x[mumuIdx] - ntupleRawTree.bG_primaryVertex_x[pvMatch],
                                        ntupleRawTree.b5_mm_kin_vtx_y[mumuIdx] - ntupleRawTree.bG_primaryVertex_y[pvMatch],
                                        ntupleRawTree.b5_mm_kin_vtx_z[mumuIdx] - ntupleRawTree.bG_primaryVertex_z[pvMatch]
                                     );
   bmmg3Momentum.SetCoordinates(
                                    bmmgLV.Px(),
                                    bmmgLV.Py(),
                                    bmmgLV.Pz()
                               );
    
   auto cos_beta       =  bmmg3Momentum.Dot(svDisplacementVecor)/sqrt(bmmg3Momentum.Mag2())/sqrt(svDisplacementVecor.Mag2());
   auto cos_beta_error =  cos_beta*-1;
   
   if(cos_beta > 0.8 ) return 2;
*/
   return 0;
}

Int_t BMMGAnalysis::getPVMatch(Int_t mumuIdx)
{

    auto zMuMu= ntupleRawTree.b5_mm_kin_pv_z[mumuIdx];
    Int_t pvMatch=0;
    auto dz = fabs(ntupleRawTree.bG_primaryVertex_z[0] -zMuMu) ; 
    
    for(int pvIdx=0; pvIdx < ntupleRawTree.bG_nPrimaryVertex ; pvIdx++)
    {
        if(dz > fabs(ntupleRawTree.bG_primaryVertex_z[pvIdx] -zMuMu) )
        {
            dz = fabs(ntupleRawTree.bG_primaryVertex_z[pvIdx] -zMuMu) ; 
            pvMatch=pvIdx;
        }
        if( dz < 1e-3) break;
    }
    
    return pvMatch;

}

