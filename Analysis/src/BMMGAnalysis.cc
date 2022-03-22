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
    Long64_t EventCountWith2GoodMuons=0;
    Long64_t EventCountWithCand=0;
    Long64_t EventCountWithDimuCand=0;
    Long64_t EventCountWithDimuVertexCand=0;
    Long64_t nb = 0,nbytes=0 ;
    Int_t eventLostAt(0);
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    auto nDiMuNoVertexCandidates=0;
    Int_t prevRun(-1),prevLumi(-1);
    bool goodRunLumi = false;
    Int_t n2GoodMuonInDimuon;
    th1fStore["ProcessingSummary"]->Fill("DataAnalysis",1);

    int rslt;

    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {   

       eventLostAt=0; 
       rslt=0;
       nDiMuCandidates=0;
       isTriggerd=false;
       n2GoodMuonInDimuon =0;
       nDiMuNoVertexCandidates =0;
       nDiMuCandidates =0;
       nBMMGCandidates =0;
 	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);

       if (ientry_evt < 0) break;

       th1fStore["ProcessingSummary"]->Fill("TotalEvents",1);
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
        
       if(doRunLumiLog) { runLumiLogger.addRunLumi(ntupleRawTree.b5_run,ntupleRawTree.b5_luminosityBlock); }


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

        if(not goodRunLumi) rslt++;
       }

       th1fStore["ProcessingSummary"]->Fill("GoodLumiEvents",1);
       // std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "<<"\n";
       // Checks if atleast 1 PV is there .. by default there will always be  one pV , the beamspot
       if(ntupleRawTree.bG_nPrimaryVertex < 1 ) rslt++;
       
       // Trigger Selection
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Jpsi ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing ) isTriggerd=true;
       //if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R ) isTriggerd=true;
    
       if( (not isTriggerd) and doTriggerFiltering ) rslt++;


       th1fStore["ProcessingSummary"]->Fill("TriggerPassEvents",1);

       for(Int_t i=0;i<NSTORAGE_ARRAY_MAX;i++)
            storageArrayDouble[i]=0;
       for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
            photonSelectionCheck[phoSCIdx]=-1;
       if(doPhotonMVA)
       {
            doPhotonMVAScores();
       }
       if(doDimuonMVA)
       {
            doDimuonMVAScores();
       }
        
                

       nBMMGCandidates=0;

       fill_muonHists();
       fill_scHists();
       fill_photonHists();
       nDiMuNoVertexCandidates=0;
       
       if(rslt > 0 )
       {
            eventLostAt=cutFlowOffsets["evtSelection"]+2;
       }
       else { 
       eventLostAt=cutFlowOffsets["basicCuts"];

       //std::cout<<"Setiing eventLostAt = "<<eventLostAt<<"\n";
       rslt=0;
       if(ntupleRawTree.b5_nmm <1 ) 
       { 
            if(eventLostAt < ( 1 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+1;
       }
       if(ntupleRawTree.bG_nSC<1)
       {
              if(eventLostAt < ( 2 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+2;
       }

       n2GoodMuonInDimuon=0;

       for(UInt_t mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    

                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] < 0 or ntupleRawTree.b5_mm_mu2_index[mumuIdx] <0 ) 
	            {
                    if(eventLostAt < ( 3 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+3;
                    continue;
                }
                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    if(eventLostAt < ( 4 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+4;
                    continue;
	            }
                if(ntupleRawTree.b5_mm_mu2_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    if(eventLostAt < ( 5 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+5;
                    continue;
	            }
           
           fill_dimuonHists(mumuIdx);
           // Muon Selection
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu1_index[mumuIdx], true);
           if(eventLostAt < rslt) eventLostAt=rslt;
           if(rslt > 0) continue;
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu2_index[mumuIdx], false);
           if(eventLostAt < rslt) eventLostAt=rslt;
           if(rslt > 0) continue;
           n2GoodMuonInDimuon++;

           // Dimuon Selection
           rslt=doDimuonSelection(mumuIdx);
           if(eventLostAt < rslt) eventLostAt=rslt;
           if(rslt > 0) { continue ; };
           nDiMuNoVertexCandidates++;

           // Dimuon Vertex selection
           rslt = doVertexSelection(mumuIdx);
           if(eventLostAt < rslt) eventLostAt=rslt;
           if(rslt > 0) { continue ; };

           auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
                      ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
           storageArrayDouble[nDiMuCandidates + candidateMapDouble["mumu_dr"] ]   = dr ;
           
           // Filling the histograms after dimuon selection
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
               
               rslt=photonSelectionCheck[phoSCIdx];
               if(eventLostAt < rslt) eventLostAt=rslt;
               if( rslt > 0) continue;

               rslt=doBMMGSelection(mumuIdx,phoSCIdx);
               if(eventLostAt < rslt) eventLostAt=rslt;
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
               //std::cout<<"Found Bmmg for "<<phoSCIdx<<" , "<<mumuIdx<<"\n";
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

       }
       
    
       EventCount++;
       if(n2GoodMuonInDimuon > 0)
       {
            th1fStore["ProcessingSummary"]->Fill("EventCountWith2GoodMuons",1);
            EventCountWith2GoodMuons++;
       }
       if(nDiMuNoVertexCandidates > 0)
       {
            th1fStore["ProcessingSummary"]->Fill("EventCountWithDimuCand",1);
            EventCountWithDimuCand++;
       }
       if(nDiMuCandidates > 0)
       {
            EventCountWithDimuVertexCand++;
            th1fStore["ProcessingSummary"]->Fill("EventCountWithDimuVertexCand",1);
       }
       if(nBMMGCandidates >0)
       {
   //   std::cout<<"\t "<<ntupleRawTree.b5_run<<" , "<<ntupleRawTree.b5_event<<" : Number of nBMMGCandidates = "<<nBMMGCandidates<<"\n";
        EventCountWithCand++;
        th1fStore["ProcessingSummary"]->Fill("EventCountWithBMMGCand",1);
        eventLostAt=0;
       }
   //         std::cout<<"EventCountWith2GoodMuons : "<<EventCountWith2GoodMuons<<"\n";
   //         std::cout<<"EventCountWithDimuCand : "<<EventCountWithDimuCand<<"\n";
   //         std::cout<<"EventCountWithDimuVertexCand : "<<EventCountWithDimuVertexCand<<"\n";
   //         std::cout<<"EventCountWithBMMGCand : "<<EventCountWithCand<<"\n";
       
       fill_globalEventHists();
   //    std::cout<<"nBMMGCandidates  "<< nBMMGCandidates <<" , "<<eventLostAt<<"\n";
       FillCutFlow(eventLostAt);
       
       if(outTree)  outTree->Fill();
    }
    std::cout<<"\n\n"
            <<"  Number of Events with trigger = "<<EventCount<<"\n"
            <<"  Number of Events with 2 good muon = "<<EventCountWith2GoodMuons<<"\n"
            <<"  Number of Events with dimu candidates = "<<EventCountWithDimuCand<<"\n"
            <<"  Number of Events with dimu vtx candidates = "<<EventCountWithDimuVertexCand<<"\n"
            <<"  Number of Events with candidates = "<<EventCountWithCand<<"\n"
            <<"  Analysis loop completed"
            <<"  \n\n";

}

#ifdef __MCANALYSIS__
void BMMGAnalysis::GenAnalyze()
{
 
    /************************************************************************************

 Make sure the branches used here are not turned off to 0 by BMMGAnalysis::setupBranchStatus()

    *************************************************************************************/
    Double_t dr;
    
    std::cout<<"\nBegining GEN Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithTrigger=0;
    Long64_t EventCountWithCand=0;
    Long64_t EventCountWith2GoodMuonsInDimu=0;
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
    Int_t eventLostAt=0,n2GoodMuonInDimuon;

    th1fStore["ProcessingSummary"]->Fill("GenAnalyze",1);

    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {   
       eventLostAt=0;

       nDiMuCandidates=0;
       isTriggerd=false;
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);

       if (ientry_evt < 0) break;

       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       th1fStore["ProcessingSummary"]->Fill("TotalEvents",1);
       
       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
      
       // Checks if atleast 1 PV is there .. by default there will always be  one pV , the beamspot
       if(ntupleRawTree.bG_nPrimaryVertex < 1 ) continue;
       
    //   Trigger Selection
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Jpsi ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R ) isTriggerd=true;
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
       n2GoodMuonInDimuon=0;
       Int_t mumMatchIdx(-1),mupMatchIdx(-1),phoMatchIdx(-1);
       Int_t mumFoundIdx(-1),mupFoundIdx(-1),phoFoundIdx(-1);
       Int_t motherIdx(-1);

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
   
       th1fStore["ProcessingSummary"]->Fill("nGenDecay",1);
       TLorentzVector mumLV,mupLV;

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
      // std::cout<<"Mu - eta, phi : "<<ntupleRawTree.b5_GenPart_eta[mumMatchIdx]<<" , "<<ntupleRawTree.b5_GenPart_phi[mumMatchIdx]<<"\n";
      // std::cout<<"Mu + eta, phi : "<<ntupleRawTree.b5_GenPart_eta[mupMatchIdx]<<" , "<<ntupleRawTree.b5_GenPart_phi[mupMatchIdx]<<"\n";
      // std::cout<<"\n";

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
       
       th1fStore["gen_mumDeltaR"]->Fill(mumRecoMatchDr);
       th1fStore["gen_mupDeltaR"]->Fill(mupRecoMatchDr);
       th1fStore["gen_phoDeltaR"]->Fill(phoRecoMatchDr);
      // std::cout<<"mup , mum , pho Dr : "<<mumRecoMatchDr<<" , "<<mupRecoMatchDr<<" , "<<phoRecoMatchDr<<"\n";
       if(mumRecoMatchDr < drMaxForGenMatchMu) { 
        th1fStore["ProcessingSummary"]->Fill("nGenMatchedMuM",1);
        mumMatchCount++;
		fill_muonHists(mumRecoMatchIdx);
	 //	std::cout<<"Match Found for : mu-"<<"\n";
	  }
      else
      {
        mumRecoMatchIdx=-1;
      }
       if(mupRecoMatchDr < drMaxForGenMatchMu) { 
          mupMatchCount++;
        th1fStore["ProcessingSummary"]->Fill("nGenMatchedMuP",1);
		fill_muonHists(mupRecoMatchIdx);
	//	std::cout<<"Match Found for : mu+"<<"\n";
	  }
      else
      {
        mupRecoMatchIdx=-1;
      }

      if(phoRecoMatchDr < drMaxForGenMatchSC ) { 

        th1fStore["ProcessingSummary"]->Fill("nGenMatchedSC",1);
        scMatchCount++;
		fill_scHists(phoRecoMatchIdx);
        if(doPhotonSelection(phoRecoMatchIdx) ==0)
        {   
            th1fStore["ProcessingSummary"]->Fill("nGenMatchedPho",1);
            phoMatchCount++;
        }

		fill_photonHists(phoRecoMatchIdx);
     }
      else
      {
         phoRecoMatchIdx=-1;
      }
     


       if(mumRecoMatchDr < drMaxForGenMatchMu and mupRecoMatchDr < drMaxForGenMatchSC )
       {
            th1fStore["ProcessingSummary"]->Fill("nGenMatchedMM",1);
       }

       if(mumRecoMatchDr < drMaxForGenMatchMu and mupRecoMatchDr <  drMaxForGenMatchSC and phoRecoMatchDr < drMaxForGenMatchSC)
       {
            fullEventMatchesFound++;
            th1fStore["ProcessingSummary"]->Fill("nGenMatchedMMG",1);
       }
       

       if(mumRecoMatchIdx < 0 )
       {
            eventLostAt=cutFlowOffsets["genSelection"]+2;    
       }
       else if ( mupRecoMatchIdx <0 )
       {
            eventLostAt=cutFlowOffsets["genSelection"]+3;    
       }
       else 
       {
            nDiMuNoVertexCandidates=0;
            if(ntupleRawTree.b5_nmm <1 ) 
            { 
                 rslt=cutFlowOffsets["basicCuts"]+1;
                 if(rslt > eventLostAt ) eventLostAt=rslt;
            }
            if(ntupleRawTree.bG_nSC < 1 ) 
            { 
                 rslt= cutFlowOffsets["basicCuts"]+2;
                 if(rslt > eventLostAt ) eventLostAt=rslt;
            }
            for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
            {   

                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] < 0 or ntupleRawTree.b5_mm_mu2_index[mumuIdx] <0 ) 
	            {
                    if(eventLostAt < ( 3 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+3;
                    continue;
                }
                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    if(eventLostAt < ( 4 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+4;
                    continue;
	            }
                if(ntupleRawTree.b5_mm_mu2_index[mumuIdx] >= ntupleRawTree.b5_nMuon ) 
	            {
                    if(eventLostAt < ( 5 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+5;
                    continue;
	            }
               
               if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] != mumRecoMatchIdx and ntupleRawTree.b5_mm_mu2_index[mumuIdx] != mumRecoMatchIdx ) {
                  rslt=cutFlowOffsets["genSelection"]+1;
                  if(rslt >eventLostAt ) eventLostAt=rslt;
                  continue;
                }
                if(ntupleRawTree.b5_mm_mu1_index[mumuIdx] != mupRecoMatchIdx and ntupleRawTree.b5_mm_mu2_index[mumuIdx] != mupRecoMatchIdx ) { 
                  rslt=cutFlowOffsets["genSelection"]+2;
                  if(rslt >eventLostAt ) eventLostAt=rslt;
                  continue;  
                }
                n2GoodMuonInDimuon++;
                // Muon Selection
                fill_dimuonHists(mumuIdx);
                fill_dimuonEnvironmentHists(mumuIdx);
	            
                      
                // Dimuon Selection
                rslt=doDimuonSelection(mumuIdx);
                if(eventLostAt < rslt) eventLostAt=rslt;
                if(rslt > 0) continue;        
                nDiMuNoVertexCandidates++;
                
                //Dimuon Vertex Selection
                //rslt = doVertexSelection(mumuIdx);
                // std::cout<<EventCountWithDimuVertexCand<<" doVertexSelection : "<<rslt<<"\n";
                if(eventLostAt < rslt) eventLostAt=rslt;
                if(rslt > 0) continue;
                

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

                if(phoRecoMatchIdx < 0 ){
                    rslt=cutFlowOffsets["gen"]+3;
                    if(rslt > eventLostAt) eventLostAt=rslt;
                }
                else 
                {

                    for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
                    {
                         if(phoSCIdx != phoRecoMatchIdx) { 
                                rslt=cutFlowOffsets["photonSelection"]+4;
                                if(rslt > eventLostAt ) eventLostAt=rslt;
                                continue;
                         }
                        // photon selection
                        if(photonSelectionCheck[phoSCIdx] < 0)
                        {
                             photonSelectionCheck[phoSCIdx]=doPhotonSelection(phoSCIdx );
                        }
                        
                        rslt = photonSelectionCheck[phoSCIdx] ;
                        //std::cout<<"photonSelectionCheck : "<<rslt<<"\n";
                        if(eventLostAt < rslt) eventLostAt=rslt;
                        if( rslt > 0) continue;
                        
                        rslt=doBMMGSelection(mumuIdx,phoSCIdx);
                        //std::cout<<"doBMMGSelection : "<<rslt<<"\n";
                        if(eventLostAt < rslt) eventLostAt=rslt;
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
                
                }
                storageArrayInt[ nDiMuCandidates + candidateMapInt["nBMMGCandidatesPerDimu"]  ]   = nBMMGCandidatesPerDimu ;
                nDiMuCandidates++;

                if(nDiMuCandidates >= NDIMU_MAX)
                {
                     std::cout<<" Dimuon count per event above NDIMU_MAX , Aborting !! \n";
                     exit(9);
                }
                 
            }
       } 
       
       //std::cout<<"nDiMuNoVertexCandidates "<<nDiMuNoVertexCandidates <<"\n";
       EventCount++;
       if(n2GoodMuonInDimuon > 0)
       {
           dimuMatchCount++;
           th1fStore["ProcessingSummary"]->Fill("EventCountWith2GoodMuons",1);
       }
       if(nDiMuNoVertexCandidates > 0)
       {
            EventCountWithDimuCand++;
            th1fStore["ProcessingSummary"]->Fill("EventCountWithGenMatchedDimu",1);
       }
       if(nDiMuCandidates > 0)
       {
            EventCountWithDimuVertexCand++;
             th1fStore["ProcessingSummary"]->Fill("EventCountWithGenMatchedDimuAfterVertexSelection",1);
       }
       if(nBMMGCandidates >0)
       {
        // std::cout<<"\t "<<ntupleRawTree.b5_run<<" , "<<ntupleRawTree.b5_event<<" : Number of nBMMGCandidates = "<<nBMMGCandidates<<"\n";
        th1fStore["ProcessingSummary"]->Fill("EventCountWithGenMatchedMMG",1);
        eventLostAt=0;
        EventCountWithCand++;
       }
       
     fill_globalEventHists();
     FillCutFlow(eventLostAt);
     if(outTree)  outTree->Fill();
    }
    std::cout<<"\n\n"
            <<"  Number of Events                                    = "<<EventCount<<"\n"
            <<"  Number of Events with trigger                       = "<<EventCountWithTrigger<<"\n"
            <<"  Number of Events with mu+ reconstructed             = "<<mupMatchCount<<"\n"
            <<"  Number of Events with mu- reconstructed             = "<<mumMatchCount<<"\n"
            <<"  Number of Events with sc reconstructed              = "<<scMatchCount<<"\n"
            <<"  Number of Events with pho reconstructed             = "<<phoMatchCount<<"\n"
            <<"  Number of Events with mu+/mu- reconstructed as dimu = "<<dimuMatchCount<<"\n"
            <<"  Number of Events with all reconstructed             = "<<fullEventMatchesFound<<"\n"
            <<"  Number of Events with dimu candidates               = "<<EventCountWithDimuCand<<"\n"
            <<"  Number of Events with dimu vtx candidates           = "<<EventCountWithDimuVertexCand<<"\n"
            <<"  Number of Events with candidates                    = "<<EventCountWithCand<<"\n"
            <<"  Analysis loop completed"
            <<"  \n\n";




}

#endif

void BMMGAnalysis::setUpDimuonMVA()
{
    if(not hasDimuonWeightFiles)
    {
        std::cout<<"Weights are not provided exiting ! \n";
        exit(4);
    }

    std::cout<<" ******** Setting up DimuonMVA *********\n"
             <<"Weight file : "<<dimuonMVAWeightFile<<"\n";
    
    readerStore["dimuonMVA"] =  new TMVA::Reader( "!Color:!Silent" );

    storageFloat["mm_bdt"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_bdt" ,  &(storageFloat["mm_bdt"]));
    storageFloat["mm_doca"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_doca" ,  &(storageFloat["mm_doca"]));
    storageFloat["mm_docatrk"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_docatrk" ,  &(storageFloat["mm_docatrk"]));
    storageFloat["mm_iso"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_iso" ,  &(storageFloat["mm_iso"]));
    storageFloat["mm_kal_lxy"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kal_lxy" ,  &(storageFloat["mm_kal_lxy"]));
    storageFloat["mm_kal_mass"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kal_mass" ,  &(storageFloat["mm_kal_mass"]));
    storageFloat["mm_kal_slxy"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kal_slxy" ,  &(storageFloat["mm_kal_slxy"]));
    storageFloat["mm_kal_vtx_prob"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kal_vtx_prob" ,  &(storageFloat["mm_kal_vtx_prob"]));
    storageFloat["mm_kin_pt"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pt" ,  &(storageFloat["mm_kin_pt"]));
    storageFloat["mm_kin_mu1pt"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_mu1pt" ,  &(storageFloat["mm_kin_mu1pt"]));
    storageFloat["mm_kin_mu2pt"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_mu2pt" ,  &(storageFloat["mm_kin_mu2pt"]));
    storageFloat["mm_kin_l3d"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_l3d" ,  &(storageFloat["mm_kin_l3d"]));
    storageFloat["mm_kin_lxy"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_lxy" ,  &(storageFloat["mm_kin_lxy"]));
    storageFloat["mm_kin_mass"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_mass" ,  &(storageFloat["mm_kin_mass"]));
    storageFloat["mm_kin_massErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_massErr" ,  &(storageFloat["mm_kin_massErr"]));
    storageFloat["mm_kin_pv2ip"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pv2ip" ,  &(storageFloat["mm_kin_pv2ip"]));
    storageFloat["mm_kin_pv2lip"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pv2lip" ,  &(storageFloat["mm_kin_pv2lip"]));
    storageFloat["mm_kin_pv2lipErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pv2lipErr" ,  &(storageFloat["mm_kin_pv2lipErr"]));
    storageFloat["mm_kin_pv2lipSig"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pv2lipSig" ,  &(storageFloat["mm_kin_pv2lipSig"]));
    storageFloat["mm_kin_pv_zErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pv_zErr" ,  &(storageFloat["mm_kin_pv_zErr"]));
    storageFloat["mm_kin_pvip"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pvip" ,  &(storageFloat["mm_kin_pvip"]));
    storageFloat["mm_kin_pvipErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pvipErr" ,  &(storageFloat["mm_kin_pvipErr"]));
    storageFloat["mm_kin_pvlip"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pvlip" ,  &(storageFloat["mm_kin_pvlip"]));
    storageFloat["mm_kin_pvlipErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pvlipErr" ,  &(storageFloat["mm_kin_pvlipErr"]));
    storageFloat["mm_kin_pvlipSig"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_pvlipSig" ,  &(storageFloat["mm_kin_pvlipSig"]));
    storageFloat["mm_kin_sl3d"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_sl3d" ,  &(storageFloat["mm_kin_sl3d"]));
    storageFloat["mm_kin_slxy"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_slxy" ,  &(storageFloat["mm_kin_slxy"]));
    storageFloat["mm_kin_spv2ip"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_spv2ip" ,  &(storageFloat["mm_kin_spv2ip"]));
    storageFloat["mm_kin_spvip"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_spvip" ,  &(storageFloat["mm_kin_spvip"]));
    storageFloat["mm_kin_vtx_chi2dof"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_vtx_chi2dof" ,  &(storageFloat["mm_kin_vtx_chi2dof"]));
    storageFloat["mm_kin_vtx_prob"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_vtx_prob" ,  &(storageFloat["mm_kin_vtx_prob"]));
    storageFloat["mm_kin_vtx_xErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_vtx_xErr" ,  &(storageFloat["mm_kin_vtx_xErr"]));
    storageFloat["mm_kin_vtx_yErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_vtx_yErr" ,  &(storageFloat["mm_kin_vtx_yErr"]));
    storageFloat["mm_kin_vtx_zErr"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_kin_vtx_zErr" ,  &(storageFloat["mm_kin_vtx_zErr"]));
    storageFloat["mm_mva"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_mva" ,  &(storageFloat["mm_mva"]));
    storageFloat["mm_otherVtxMaxProb1"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_otherVtxMaxProb1" ,  &(storageFloat["mm_otherVtxMaxProb1"]));
    storageFloat["mm_closetrk"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_closetrk" ,  &(storageFloat["mm_closetrk"]));
    storageFloat["mm_closetrks1"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_closetrks1" ,  &(storageFloat["mm_closetrks1"]));
    storageFloat["mm_nBMTrks"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_nBMTrks" ,  &(storageFloat["mm_nBMTrks"]));
    storageFloat["mm_nDisTrks"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_nDisTrks" ,  &(storageFloat["mm_nDisTrks"]));
    storageFloat["mm_nTrks"] = 0.0; readerStore["dimuonMVA"]->AddVariable("mm_nTrks" ,  &(storageFloat["mm_nTrks"]));

    readerStore["dimuonMVA"]->BookMVA("DimuonMVA", dimuonMVAWeightFile );
    
    hasSetupDimuonMVA=true;
}

void BMMGAnalysis::doDimuonMVAScores()
{   
    if(not hasSetupDimuonMVA)
    {
        std::cout<<" MVA not setup : "<<hasSetupPhotonMVA<<"\n";
        exit(2);
    }

    Double_t dimuonMVAScore;
    for(int  i=0  ; i < ntupleRawTree.b5_nmm ; i++)
    {
    
      storageFloat["mm_bdt"                  ]  = ntupleRawTree.b5_mm_bdt[i] ;
      storageFloat["mm_doca"                  ]  = ntupleRawTree.b5_mm_doca[i] ;
      storageFloat["mm_docatrk"                  ]  = ntupleRawTree.b5_mm_docatrk[i] ;
      storageFloat["mm_iso"                  ]  = ntupleRawTree.b5_mm_iso[i] ;
      storageFloat["mm_kal_lxy"                  ]  = ntupleRawTree.b5_mm_kal_lxy[i] ;
      storageFloat["mm_kal_mass"                  ]  = ntupleRawTree.b5_mm_kal_mass[i] ;
      storageFloat["mm_kal_slxy"                  ]  = ntupleRawTree.b5_mm_kal_slxy[i] ;
      storageFloat["mm_kal_vtx_prob"                  ]  = ntupleRawTree.b5_mm_kal_vtx_prob[i] ;
      storageFloat["mm_kin_pt"                  ]  = ntupleRawTree.b5_mm_kin_pt[i] ;
      storageFloat["mm_kin_mu1pt"                  ]  = ntupleRawTree.b5_mm_kin_mu1pt[i] ;
      storageFloat["mm_kin_mu2pt"                  ]  = ntupleRawTree.b5_mm_kin_mu2pt[i] ;
      storageFloat["mm_kin_l3d"                  ]  = ntupleRawTree.b5_mm_kin_l3d[i] ;
      storageFloat["mm_kin_lxy"                  ]  = ntupleRawTree.b5_mm_kin_lxy[i] ;
      storageFloat["mm_kin_mass"                  ]  = ntupleRawTree.b5_mm_kin_mass[i] ;
      storageFloat["mm_kin_massErr"                  ]  = ntupleRawTree.b5_mm_kin_massErr[i] ;
      storageFloat["mm_kin_pv2ip"                  ]  = ntupleRawTree.b5_mm_kin_pv2ip[i] ;
      storageFloat["mm_kin_pv2lip"                  ]  = ntupleRawTree.b5_mm_kin_pv2lip[i] ;
      storageFloat["mm_kin_pv2lipErr"                  ]  = ntupleRawTree.b5_mm_kin_pv2lipErr[i] ;
      storageFloat["mm_kin_pv2lipSig"                  ]  = ntupleRawTree.b5_mm_kin_pv2lipSig[i] ;
      storageFloat["mm_kin_pv_zErr"                  ]  = ntupleRawTree.b5_mm_kin_pv_zErr[i] ;
      storageFloat["mm_kin_pvip"                  ]  = ntupleRawTree.b5_mm_kin_pvip[i] ;
      storageFloat["mm_kin_pvipErr"                  ]  = ntupleRawTree.b5_mm_kin_pvipErr[i] ;
      storageFloat["mm_kin_pvlip"                  ]  = ntupleRawTree.b5_mm_kin_pvlip[i] ;
      storageFloat["mm_kin_pvlipErr"                  ]  = ntupleRawTree.b5_mm_kin_pvlipErr[i] ;
      storageFloat["mm_kin_pvlipSig"                  ]  = ntupleRawTree.b5_mm_kin_pvlipSig[i] ;
      storageFloat["mm_kin_sl3d"                  ]  = ntupleRawTree.b5_mm_kin_sl3d[i] ;
      storageFloat["mm_kin_slxy"                  ]  = ntupleRawTree.b5_mm_kin_slxy[i] ;
      storageFloat["mm_kin_spv2ip"                  ]  = ntupleRawTree.b5_mm_kin_spv2ip[i] ;
      storageFloat["mm_kin_spvip"                  ]  = ntupleRawTree.b5_mm_kin_spvip[i] ;
      storageFloat["mm_kin_vtx_chi2dof"                  ]  = ntupleRawTree.b5_mm_kin_vtx_chi2dof[i] ;
      storageFloat["mm_kin_vtx_prob"                  ]  = ntupleRawTree.b5_mm_kin_vtx_prob[i] ;
      storageFloat["mm_kin_vtx_xErr"                  ]  = ntupleRawTree.b5_mm_kin_vtx_xErr[i] ;
      storageFloat["mm_kin_vtx_yErr"                  ]  = ntupleRawTree.b5_mm_kin_vtx_yErr[i] ;
      storageFloat["mm_kin_vtx_zErr"                  ]  = ntupleRawTree.b5_mm_kin_vtx_zErr[i] ;
      storageFloat["mm_mva"                  ]  = ntupleRawTree.b5_mm_mva[i] ;
      storageFloat["mm_otherVtxMaxProb1"                  ]  = ntupleRawTree.b5_mm_otherVtxMaxProb1[i] ;
      storageFloat["mm_closetrk"                  ]  = ntupleRawTree.b5_mm_closetrk[i] ;
      storageFloat["mm_closetrks1"                  ]  = ntupleRawTree.b5_mm_closetrks1[i] ;
      storageFloat["mm_nBMTrks"                  ]  = ntupleRawTree.b5_mm_nBMTrks[i] ;
      storageFloat["mm_nDisTrks"                  ]  = ntupleRawTree.b5_mm_nDisTrks[i] ;
      storageFloat["mm_nTrks"                  ]  = ntupleRawTree.b5_mm_nTrks[i] ;
      
      dimuonMVAScore = readerStore["dimuonMVA"]->EvaluateMVA("DimuonMVA");
      storageArrayDouble[ i + candidateMapInt["dimuonMVAScore"]  ]   = dimuonMVAScore ;
      //std::cout<<" i = "<<dimuonMVAScore<<"\n";
    }
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
    if(not reader) delete reader;
    reader =  new TMVA::Reader( "!Color:!Silent" );
    
    for(Int_t i=0;i<mvaTrainVars.size();i++)
    {
        std::cout<<"Adding mva vars "<<mvaTrainVars[i]<<" to "<<" reader  \n";
        storageDouble[mvaTrainVars[i]]=0.0;
        storageFloat[mvaTrainVars[i]] = storageDouble[mvaTrainVars[i]];
        reader->AddVariable(mvaTrainVars[i].c_str(), &(storageFloat[mvaTrainVars[i]]));
    }
    
  for(int i =0;i<spectatorVars.size();i++)
  {
        std::cout<<"Adding spectator "<<spectatorVars[i]<<" to "<<" reader \n";
        storageDouble[spectatorVars[i]]=0.0;
        storageFloat[spectatorVars[i]] = storageDouble[spectatorVars[i]];
        reader->AddSpectator(spectatorVars[i].c_str(),  &(storageFloat[spectatorVars[i]]));
  }
    //reader->BookMVA("LowPtPhotonIdMVA_MLP", photonIdxMVAWeightFile );
    std::cout<<"Barrel File : "<<photonIdxBarrelMVAWeightFile<<"\n";
    reader->BookMVA("LowPtPhotonIdBarrelMVA_MLP", photonIdxBarrelMVAWeightFile );
    std::cout<<"ECap File : "<<photonIdxECapMVAWeightFile<<"\n";
    reader->BookMVA("LowPtPhotonIdECapMVA_MLP", photonIdxECapMVAWeightFile );
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
      storageFloat["scR9"                  ]  = ntupleRawTree.bG_scR9[i] ;
      storageFloat["scSigmaIetaIeta"       ]  = ntupleRawTree.bG_scSigmaIetaIeta[i] ;
      storageFloat["scFull5x5_E2x5MaxRatio"]  = ntupleRawTree.bG_scFull5x5_e2x5MaxRatio[i] ;
      storageFloat["scPFPhoIso3"           ]  = ntupleRawTree.bG_scPFPhoIso3[i] ;
      storageFloat["scSigmaIetaIphi"       ]  = ntupleRawTree.bG_scSigmaIetaIphi[i] ;
      storageFloat["scFull5x5_SwissCross"  ]  = ntupleRawTree.bG_scFull5x5_swissCross[i] ;
      storageFloat["scEtaWidth"            ]  = ntupleRawTree.bG_scEtaWidth[i] ;
      storageFloat["scE2ndRatio"           ]  = ntupleRawTree.bG_scE2ndRatio[i] ;
      storageFloat["scEMaxRatio"           ]  = ntupleRawTree.bG_scEMaxRatio[i] ;
      storageFloat["scPFChIso1"            ]  = ntupleRawTree.bG_scPFChIso1[i] ;
      storageFloat["scPFChIso2"            ]  = ntupleRawTree.bG_scPFChIso2[i] ;
      storageFloat["scPFChIso3"            ]  = ntupleRawTree.bG_scPFChIso3[i] ;
      storageFloat["scPFChIso4"            ]  = ntupleRawTree.bG_scPFChIso4[i] ;
      storageFloat["scPFChIso5"            ]  = ntupleRawTree.bG_scPFChIso5[i] ;
     
      storageFloat["scRawE"            ]  = ntupleRawTree.bG_scRawE[i];
      storageFloat["scE"               ]  = ntupleRawTree.bG_scE[i];
      storageFloat["scRawE"            ]  = ntupleRawTree.bG_scRawE[i];
      storageFloat["scPhi"             ]  = ntupleRawTree.bG_scPhi[i];
      storageFloat["scX"               ]  = ntupleRawTree.bG_scX[i];
      storageFloat["scY"               ]  = ntupleRawTree.bG_scY[i];
      storageFloat["scZ"               ]  = ntupleRawTree.bG_scZ[i];
      
      if(abs(ntupleRawTree.bG_scEta[i]) <1.48) 
      {
        photonMVAValue = reader->EvaluateMVA("LowPtPhotonIdBarrelMVA_MLP");
        //std::cout<<"barrel "<<ntupleRawTree.bG_scEta[i]<<" : "<<photonMVAValue<<"\n";
      }
      else if(abs(ntupleRawTree.bG_scEta[i]) < 2.5)
      {
        photonMVAValue = reader->EvaluateMVA("LowPtPhotonIdECapMVA_MLP");
        //std::cout<<"ecap "<<ntupleRawTree.bG_scEta[i]<<" : "<<photonMVAValue<<"\n";
      }
      else
      {
        photonMVAValue = -1.0;
      }
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
    
    if(doDimuonMVA)
    {
        candidateMapDouble["dimuonMVAScore"]   = storageIdxFilledDouble ;
        outTree->Branch("dimuonMVAScore",&storageArrayDouble[storageIdxFilledDouble],"dimuonMVAScore[b5_nmm]/D"); storageIdxFilledDouble+=NDIMU_MAX;
    }

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
    AllocateBMMGBranches();
   if(not ntupleRawTree.fChain) 
   {
        std::cout<<"input Tree  Not Setup !!\n";
        exit(14);

   }

	outTree->Branch("b5_run",&	( ntupleRawTree.b5_run ));
	outTree->Branch("b5_luminosityBlock",&	( ntupleRawTree.b5_luminosityBlock ));
	outTree->Branch("b5_event",&	( ntupleRawTree.b5_event ));
	outTree->Branch("b5_nMuonId",&	( ntupleRawTree.b5_nMuonId ));
	outTree->Branch("b5_MuonId_chi2LocalPosition",	( ntupleRawTree.b5_MuonId_chi2LocalPosition ),"b5_MuonId_chi2LocalPosition[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_glbNormChi2",	( ntupleRawTree.b5_MuonId_glbNormChi2 ),"b5_MuonId_glbNormChi2[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_glbTrackProbability",	( ntupleRawTree.b5_MuonId_glbTrackProbability ),"b5_MuonId_glbTrackProbability[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match1_dX",	( ntupleRawTree.b5_MuonId_match1_dX ),"b5_MuonId_match1_dX[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match1_dY",	( ntupleRawTree.b5_MuonId_match1_dY ),"b5_MuonId_match1_dY[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match1_pullDxDz",	( ntupleRawTree.b5_MuonId_match1_pullDxDz ),"b5_MuonId_match1_pullDxDz[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match1_pullDyDz",	( ntupleRawTree.b5_MuonId_match1_pullDyDz ),"b5_MuonId_match1_pullDyDz[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match1_pullX",	( ntupleRawTree.b5_MuonId_match1_pullX ),"b5_MuonId_match1_pullX[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match1_pullY",	( ntupleRawTree.b5_MuonId_match1_pullY ),"b5_MuonId_match1_pullY[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match2_dX",	( ntupleRawTree.b5_MuonId_match2_dX ),"b5_MuonId_match2_dX[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match2_dY",	( ntupleRawTree.b5_MuonId_match2_dY ),"b5_MuonId_match2_dY[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match2_pullDxDz",	( ntupleRawTree.b5_MuonId_match2_pullDxDz ),"b5_MuonId_match2_pullDxDz[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match2_pullDyDz",	( ntupleRawTree.b5_MuonId_match2_pullDyDz ),"b5_MuonId_match2_pullDyDz[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match2_pullX",	( ntupleRawTree.b5_MuonId_match2_pullX ),"b5_MuonId_match2_pullX[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_match2_pullY",	( ntupleRawTree.b5_MuonId_match2_pullY ),"b5_MuonId_match2_pullY[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_newSoftMuonMva",	( ntupleRawTree.b5_MuonId_newSoftMuonMva ),"b5_MuonId_newSoftMuonMva[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_trkKink",	( ntupleRawTree.b5_MuonId_trkKink ),"b5_MuonId_trkKink[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_trkValidFrac",	( ntupleRawTree.b5_MuonId_trkValidFrac ),"b5_MuonId_trkValidFrac[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_highPurity",	( ntupleRawTree.b5_MuonId_highPurity ),"b5_MuonId_highPurity[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_nLostHitsInner",	( ntupleRawTree.b5_MuonId_nLostHitsInner ),"b5_MuonId_nLostHitsInner[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_nLostHitsOn",	( ntupleRawTree.b5_MuonId_nLostHitsOn ),"b5_MuonId_nLostHitsOn[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_nLostHitsOuter",	( ntupleRawTree.b5_MuonId_nLostHitsOuter ),"b5_MuonId_nLostHitsOuter[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_nPixels",	( ntupleRawTree.b5_MuonId_nPixels ),"b5_MuonId_nPixels[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_nValidHits",	( ntupleRawTree.b5_MuonId_nValidHits ),"b5_MuonId_nValidHits[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_trkLayers",	( ntupleRawTree.b5_MuonId_trkLayers ),"b5_MuonId_trkLayers[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_trkLostLayersInner",	( ntupleRawTree.b5_MuonId_trkLostLayersInner ),"b5_MuonId_trkLostLayersInner[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_trkLostLayersOn",	( ntupleRawTree.b5_MuonId_trkLostLayersOn ),"b5_MuonId_trkLostLayersOn[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_trkLostLayersOuter",	( ntupleRawTree.b5_MuonId_trkLostLayersOuter ),"b5_MuonId_trkLostLayersOuter[b5_nMuonId]/I");
	outTree->Branch("b5_nmm",&	( ntupleRawTree.b5_nmm ));
	outTree->Branch("b5_mm_bdt",	( ntupleRawTree.b5_mm_bdt ),"b5_mm_bdt[b5_nmm]/F");
	outTree->Branch("b5_mm_doca",	( ntupleRawTree.b5_mm_doca ),"b5_mm_doca[b5_nmm]/F");
	outTree->Branch("b5_mm_docatrk",	( ntupleRawTree.b5_mm_docatrk ),"b5_mm_docatrk[b5_nmm]/F");
	outTree->Branch("b5_mm_iso",	( ntupleRawTree.b5_mm_iso ),"b5_mm_iso[b5_nmm]/F");
	outTree->Branch("b5_mm_kal_lxy",	( ntupleRawTree.b5_mm_kal_lxy ),"b5_mm_kal_lxy[b5_nmm]/F");
	outTree->Branch("b5_mm_kal_mass",	( ntupleRawTree.b5_mm_kal_mass ),"b5_mm_kal_mass[b5_nmm]/F");
	outTree->Branch("b5_mm_kal_slxy",	( ntupleRawTree.b5_mm_kal_slxy ),"b5_mm_kal_slxy[b5_nmm]/F");
	outTree->Branch("b5_mm_kal_vtx_prob",	( ntupleRawTree.b5_mm_kal_vtx_prob ),"b5_mm_kal_vtx_prob[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_alpha",	( ntupleRawTree.b5_mm_kin_alpha ),"b5_mm_kin_alpha[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_alphaBS",	( ntupleRawTree.b5_mm_kin_alphaBS ),"b5_mm_kin_alphaBS[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_alphaBSErr",	( ntupleRawTree.b5_mm_kin_alphaBSErr ),"b5_mm_kin_alphaBSErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_alphaErr",	( ntupleRawTree.b5_mm_kin_alphaErr ),"b5_mm_kin_alphaErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_eta",	( ntupleRawTree.b5_mm_kin_eta ),"b5_mm_kin_eta[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_l3d",	( ntupleRawTree.b5_mm_kin_l3d ),"b5_mm_kin_l3d[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_lxy",	( ntupleRawTree.b5_mm_kin_lxy ),"b5_mm_kin_lxy[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_mass",	( ntupleRawTree.b5_mm_kin_mass ),"b5_mm_kin_mass[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_massErr",	( ntupleRawTree.b5_mm_kin_massErr ),"b5_mm_kin_massErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_mu1eta",	( ntupleRawTree.b5_mm_kin_mu1eta ),"b5_mm_kin_mu1eta[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_mu1phi",	( ntupleRawTree.b5_mm_kin_mu1phi ),"b5_mm_kin_mu1phi[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_mu1pt",	( ntupleRawTree.b5_mm_kin_mu1pt ),"b5_mm_kin_mu1pt[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_mu2eta",	( ntupleRawTree.b5_mm_kin_mu2eta ),"b5_mm_kin_mu2eta[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_mu2phi",	( ntupleRawTree.b5_mm_kin_mu2phi ),"b5_mm_kin_mu2phi[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_mu2pt",	( ntupleRawTree.b5_mm_kin_mu2pt ),"b5_mm_kin_mu2pt[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_phi",	( ntupleRawTree.b5_mm_kin_phi ),"b5_mm_kin_phi[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pt",	( ntupleRawTree.b5_mm_kin_pt ),"b5_mm_kin_pt[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pv2ip",	( ntupleRawTree.b5_mm_kin_pv2ip ),"b5_mm_kin_pv2ip[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pv2ipErr",	( ntupleRawTree.b5_mm_kin_pv2ipErr ),"b5_mm_kin_pv2ipErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pv2lip",	( ntupleRawTree.b5_mm_kin_pv2lip ),"b5_mm_kin_pv2lip[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pv2lipErr",	( ntupleRawTree.b5_mm_kin_pv2lipErr ),"b5_mm_kin_pv2lipErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pv2lipSig",	( ntupleRawTree.b5_mm_kin_pv2lipSig ),"b5_mm_kin_pv2lipSig[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pv_z",	( ntupleRawTree.b5_mm_kin_pv_z ),"b5_mm_kin_pv_z[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pv_zErr",	( ntupleRawTree.b5_mm_kin_pv_zErr ),"b5_mm_kin_pv_zErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pvip",	( ntupleRawTree.b5_mm_kin_pvip ),"b5_mm_kin_pvip[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pvipErr",	( ntupleRawTree.b5_mm_kin_pvipErr ),"b5_mm_kin_pvipErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pvlip",	( ntupleRawTree.b5_mm_kin_pvlip ),"b5_mm_kin_pvlip[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pvlipErr",	( ntupleRawTree.b5_mm_kin_pvlipErr ),"b5_mm_kin_pvlipErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_pvlipSig",	( ntupleRawTree.b5_mm_kin_pvlipSig ),"b5_mm_kin_pvlipSig[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_sl3d",	( ntupleRawTree.b5_mm_kin_sl3d ),"b5_mm_kin_sl3d[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_slxy",	( ntupleRawTree.b5_mm_kin_slxy ),"b5_mm_kin_slxy[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_spv2ip",	( ntupleRawTree.b5_mm_kin_spv2ip ),"b5_mm_kin_spv2ip[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_spvip",	( ntupleRawTree.b5_mm_kin_spvip ),"b5_mm_kin_spvip[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_tau",	( ntupleRawTree.b5_mm_kin_tau ),"b5_mm_kin_tau[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_taue",	( ntupleRawTree.b5_mm_kin_taue ),"b5_mm_kin_taue[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_tauxy",	( ntupleRawTree.b5_mm_kin_tauxy ),"b5_mm_kin_tauxy[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_tauxye",	( ntupleRawTree.b5_mm_kin_tauxye ),"b5_mm_kin_tauxye[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_chi2dof",	( ntupleRawTree.b5_mm_kin_vtx_chi2dof ),"b5_mm_kin_vtx_chi2dof[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_prob",	( ntupleRawTree.b5_mm_kin_vtx_prob ),"b5_mm_kin_vtx_prob[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_x",	( ntupleRawTree.b5_mm_kin_vtx_x ),"b5_mm_kin_vtx_x[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_xErr",	( ntupleRawTree.b5_mm_kin_vtx_xErr ),"b5_mm_kin_vtx_xErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_y",	( ntupleRawTree.b5_mm_kin_vtx_y ),"b5_mm_kin_vtx_y[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_yErr",	( ntupleRawTree.b5_mm_kin_vtx_yErr ),"b5_mm_kin_vtx_yErr[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_z",	( ntupleRawTree.b5_mm_kin_vtx_z ),"b5_mm_kin_vtx_z[b5_nmm]/F");
	outTree->Branch("b5_mm_kin_vtx_zErr",	( ntupleRawTree.b5_mm_kin_vtx_zErr ),"b5_mm_kin_vtx_zErr[b5_nmm]/F");
	outTree->Branch("b5_mm_m1iso",	( ntupleRawTree.b5_mm_m1iso ),"b5_mm_m1iso[b5_nmm]/F");
	outTree->Branch("b5_mm_m2iso",	( ntupleRawTree.b5_mm_m2iso ),"b5_mm_m2iso[b5_nmm]/F");
	outTree->Branch("b5_mm_mass",	( ntupleRawTree.b5_mm_mass ),"b5_mm_mass[b5_nmm]/F");
	outTree->Branch("b5_mm_mu1_eta",	( ntupleRawTree.b5_mm_mu1_eta ),"b5_mm_mu1_eta[b5_nmm]/F");
	outTree->Branch("b5_mm_mu1_phi",	( ntupleRawTree.b5_mm_mu1_phi ),"b5_mm_mu1_phi[b5_nmm]/F");
	outTree->Branch("b5_mm_mu1_pt",	( ntupleRawTree.b5_mm_mu1_pt ),"b5_mm_mu1_pt[b5_nmm]/F");
	outTree->Branch("b5_mm_mu2_eta",	( ntupleRawTree.b5_mm_mu2_eta ),"b5_mm_mu2_eta[b5_nmm]/F");
	outTree->Branch("b5_mm_mu2_phi",	( ntupleRawTree.b5_mm_mu2_phi ),"b5_mm_mu2_phi[b5_nmm]/F");
	outTree->Branch("b5_mm_mu2_pt",	( ntupleRawTree.b5_mm_mu2_pt ),"b5_mm_mu2_pt[b5_nmm]/F");
	outTree->Branch("b5_mm_mva",	( ntupleRawTree.b5_mm_mva ),"b5_mm_mva[b5_nmm]/F");
	outTree->Branch("b5_mm_otherVtxMaxProb",	( ntupleRawTree.b5_mm_otherVtxMaxProb ),"b5_mm_otherVtxMaxProb[b5_nmm]/F");
	outTree->Branch("b5_mm_otherVtxMaxProb1",	( ntupleRawTree.b5_mm_otherVtxMaxProb1 ),"b5_mm_otherVtxMaxProb1[b5_nmm]/F");
	outTree->Branch("b5_mm_otherVtxMaxProb2",	( ntupleRawTree.b5_mm_otherVtxMaxProb2 ),"b5_mm_otherVtxMaxProb2[b5_nmm]/F");
	outTree->Branch("b5_mm_closetrk",	( ntupleRawTree.b5_mm_closetrk ),"b5_mm_closetrk[b5_nmm]/I");
	outTree->Branch("b5_mm_closetrks1",	( ntupleRawTree.b5_mm_closetrks1 ),"b5_mm_closetrks1[b5_nmm]/I");
	outTree->Branch("b5_mm_closetrks2",	( ntupleRawTree.b5_mm_closetrks2 ),"b5_mm_closetrks2[b5_nmm]/I");
	outTree->Branch("b5_mm_closetrks3",	( ntupleRawTree.b5_mm_closetrks3 ),"b5_mm_closetrks3[b5_nmm]/I");
	//outTree->Branch("b5_mm_gen_cpdgId",	( ntupleRawTree.b5_mm_gen_cpdgId ),"b5_mm_gen_cpdgId[b5_nmm]/I");
	//outTree->Branch("b5_mm_gen_mpdgId",	( ntupleRawTree.b5_mm_gen_mpdgId ),"b5_mm_gen_mpdgId[b5_nmm]/I");
	//outTree->Branch("b5_mm_gen_mu1_mpdgId",	( ntupleRawTree.b5_mm_gen_mu1_mpdgId ),"b5_mm_gen_mu1_mpdgId[b5_nmm]/I");
	//outTree->Branch("b5_mm_gen_mu1_pdgId",	( ntupleRawTree.b5_mm_gen_mu1_pdgId ),"b5_mm_gen_mu1_pdgId[b5_nmm]/I");
	//outTree->Branch("b5_mm_gen_mu2_mpdgId",	( ntupleRawTree.b5_mm_gen_mu2_mpdgId ),"b5_mm_gen_mu2_mpdgId[b5_nmm]/I");
	//outTree->Branch("b5_mm_gen_mu2_pdgId",	( ntupleRawTree.b5_mm_gen_mu2_pdgId ),"b5_mm_gen_mu2_pdgId[b5_nmm]/I");
	//outTree->Branch("b5_mm_gen_pdgId",	( ntupleRawTree.b5_mm_gen_pdgId ),"b5_mm_gen_pdgId[b5_nmm]/I");
	outTree->Branch("b5_mm_kal_valid",	( ntupleRawTree.b5_mm_kal_valid ),"b5_mm_kal_valid[b5_nmm]/I");
	outTree->Branch("b5_mm_kin_valid",	( ntupleRawTree.b5_mm_kin_valid ),"b5_mm_kin_valid[b5_nmm]/I");
	outTree->Branch("b5_mm_kinpc_valid",	( ntupleRawTree.b5_mm_kinpc_valid ),"b5_mm_kinpc_valid[b5_nmm]/I");
	outTree->Branch("b5_mm_mu1_index",	( ntupleRawTree.b5_mm_mu1_index ),"b5_mm_mu1_index[b5_nmm]/I");
	outTree->Branch("b5_mm_mu1_pdgId",	( ntupleRawTree.b5_mm_mu1_pdgId ),"b5_mm_mu1_pdgId[b5_nmm]/I");
	outTree->Branch("b5_mm_mu2_index",	( ntupleRawTree.b5_mm_mu2_index ),"b5_mm_mu2_index[b5_nmm]/I");
	outTree->Branch("b5_mm_mu2_pdgId",	( ntupleRawTree.b5_mm_mu2_pdgId ),"b5_mm_mu2_pdgId[b5_nmm]/I");
	outTree->Branch("b5_mm_nBMTrks",	( ntupleRawTree.b5_mm_nBMTrks ),"b5_mm_nBMTrks[b5_nmm]/I");
	outTree->Branch("b5_mm_nDisTrks",	( ntupleRawTree.b5_mm_nDisTrks ),"b5_mm_nDisTrks[b5_nmm]/I");
	outTree->Branch("b5_mm_nTrks",	( ntupleRawTree.b5_mm_nTrks ),"b5_mm_nTrks[b5_nmm]/I");
	//outTree->Branch("b5_ngensummary",&	( ntupleRawTree.b5_ngensummary ));
	//outTree->Branch("b5_gensummary_n_anti_b",	( ntupleRawTree.b5_gensummary_n_anti_b ),"b5_gensummary_n_anti_b[b5_ngensummary]/I");
	//outTree->Branch("b5_gensummary_n_b",	( ntupleRawTree.b5_gensummary_n_b ),"b5_gensummary_n_b[b5_ngensummary]/I");
	//outTree->Branch("b5_gensummary_process_type",	( ntupleRawTree.b5_gensummary_process_type ),"b5_gensummary_process_type[b5_ngensummary]/I");
	//outTree->Branch("b5_ngenbmm",&	( ntupleRawTree.b5_ngenbmm ));
	//outTree->Branch("b5_genbmm_dau3_eta",	( ntupleRawTree.b5_genbmm_dau3_eta ),"b5_genbmm_dau3_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau3_phi",	( ntupleRawTree.b5_genbmm_dau3_phi ),"b5_genbmm_dau3_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau3_pt",	( ntupleRawTree.b5_genbmm_dau3_pt ),"b5_genbmm_dau3_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau3_reco_eta",	( ntupleRawTree.b5_genbmm_dau3_reco_eta ),"b5_genbmm_dau3_reco_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau3_reco_phi",	( ntupleRawTree.b5_genbmm_dau3_reco_phi ),"b5_genbmm_dau3_reco_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau3_reco_pt",	( ntupleRawTree.b5_genbmm_dau3_reco_pt ),"b5_genbmm_dau3_reco_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau4_eta",	( ntupleRawTree.b5_genbmm_dau4_eta ),"b5_genbmm_dau4_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau4_phi",	( ntupleRawTree.b5_genbmm_dau4_phi ),"b5_genbmm_dau4_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau4_pt",	( ntupleRawTree.b5_genbmm_dau4_pt ),"b5_genbmm_dau4_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau4_reco_eta",	( ntupleRawTree.b5_genbmm_dau4_reco_eta ),"b5_genbmm_dau4_reco_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau4_reco_phi",	( ntupleRawTree.b5_genbmm_dau4_reco_phi ),"b5_genbmm_dau4_reco_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau4_reco_pt",	( ntupleRawTree.b5_genbmm_dau4_reco_pt ),"b5_genbmm_dau4_reco_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dimuon_mass",	( ntupleRawTree.b5_genbmm_dimuon_mass ),"b5_genbmm_dimuon_mass[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_eta",	( ntupleRawTree.b5_genbmm_eta ),"b5_genbmm_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_mass",	( ntupleRawTree.b5_genbmm_mass ),"b5_genbmm_mass[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_mu1_eta",	( ntupleRawTree.b5_genbmm_mu1_eta ),"b5_genbmm_mu1_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_mu1_phi",	( ntupleRawTree.b5_genbmm_mu1_phi ),"b5_genbmm_mu1_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_mu1_pt",	( ntupleRawTree.b5_genbmm_mu1_pt ),"b5_genbmm_mu1_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_mu2_eta",	( ntupleRawTree.b5_genbmm_mu2_eta ),"b5_genbmm_mu2_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_mu2_phi",	( ntupleRawTree.b5_genbmm_mu2_phi ),"b5_genbmm_mu2_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_mu2_pt",	( ntupleRawTree.b5_genbmm_mu2_pt ),"b5_genbmm_mu2_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_phi",	( ntupleRawTree.b5_genbmm_phi ),"b5_genbmm_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_pt",	( ntupleRawTree.b5_genbmm_pt ),"b5_genbmm_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_rad_eta",	( ntupleRawTree.b5_genbmm_rad_eta ),"b5_genbmm_rad_eta[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_rad_p",	( ntupleRawTree.b5_genbmm_rad_p ),"b5_genbmm_rad_p[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_rad_phi",	( ntupleRawTree.b5_genbmm_rad_phi ),"b5_genbmm_rad_phi[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_rad_pt",	( ntupleRawTree.b5_genbmm_rad_pt ),"b5_genbmm_rad_pt[b5_ngenbmm]/F");
	//outTree->Branch("b5_genbmm_dau3_pdgId",	( ntupleRawTree.b5_genbmm_dau3_pdgId ),"b5_genbmm_dau3_pdgId[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_dau4_pdgId",	( ntupleRawTree.b5_genbmm_dau4_pdgId ),"b5_genbmm_dau4_pdgId[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_mu1_good",	( ntupleRawTree.b5_genbmm_mu1_good ),"b5_genbmm_mu1_good[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_mu1_index",	( ntupleRawTree.b5_genbmm_mu1_index ),"b5_genbmm_mu1_index[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_mu1_pdgId",	( ntupleRawTree.b5_genbmm_mu1_pdgId ),"b5_genbmm_mu1_pdgId[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_mu2_good",	( ntupleRawTree.b5_genbmm_mu2_good ),"b5_genbmm_mu2_good[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_mu2_index",	( ntupleRawTree.b5_genbmm_mu2_index ),"b5_genbmm_mu2_index[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_mu2_pdgId",	( ntupleRawTree.b5_genbmm_mu2_pdgId ),"b5_genbmm_mu2_pdgId[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_pdgId",	( ntupleRawTree.b5_genbmm_pdgId ),"b5_genbmm_pdgId[b5_ngenbmm]/I");
	//outTree->Branch("b5_genbmm_signature",	( ntupleRawTree.b5_genbmm_signature ),"b5_genbmm_signature[b5_ngenbmm]/I");
	outTree->Branch("b5_nMuon",&	( ntupleRawTree.b5_nMuon ));
	outTree->Branch("b5_Muon_dxy",	( ntupleRawTree.b5_Muon_dxy ),"b5_Muon_dxy[b5_nMuon]/F");
	outTree->Branch("b5_Muon_dxyErr",	( ntupleRawTree.b5_Muon_dxyErr ),"b5_Muon_dxyErr[b5_nMuon]/F");
	outTree->Branch("b5_Muon_dxybs",	( ntupleRawTree.b5_Muon_dxybs ),"b5_Muon_dxybs[b5_nMuon]/F");
	outTree->Branch("b5_Muon_dz",	( ntupleRawTree.b5_Muon_dz ),"b5_Muon_dz[b5_nMuon]/F");
	outTree->Branch("b5_Muon_dzErr",	( ntupleRawTree.b5_Muon_dzErr ),"b5_Muon_dzErr[b5_nMuon]/F");
	outTree->Branch("b5_Muon_eta",	( ntupleRawTree.b5_Muon_eta ),"b5_Muon_eta[b5_nMuon]/F");
	outTree->Branch("b5_Muon_ip3d",	( ntupleRawTree.b5_Muon_ip3d ),"b5_Muon_ip3d[b5_nMuon]/F");
	outTree->Branch("b5_Muon_jetPtRelv2",	( ntupleRawTree.b5_Muon_jetPtRelv2 ),"b5_Muon_jetPtRelv2[b5_nMuon]/F");
	outTree->Branch("b5_Muon_jetRelIso",	( ntupleRawTree.b5_Muon_jetRelIso ),"b5_Muon_jetRelIso[b5_nMuon]/F");
	outTree->Branch("b5_Muon_mass",	( ntupleRawTree.b5_Muon_mass ),"b5_Muon_mass[b5_nMuon]/F");
	outTree->Branch("b5_Muon_miniPFRelIso_all",	( ntupleRawTree.b5_Muon_miniPFRelIso_all ),"b5_Muon_miniPFRelIso_all[b5_nMuon]/F");
	outTree->Branch("b5_Muon_miniPFRelIso_chg",	( ntupleRawTree.b5_Muon_miniPFRelIso_chg ),"b5_Muon_miniPFRelIso_chg[b5_nMuon]/F");
	outTree->Branch("b5_Muon_pfRelIso03_all",	( ntupleRawTree.b5_Muon_pfRelIso03_all ),"b5_Muon_pfRelIso03_all[b5_nMuon]/F");
	outTree->Branch("b5_Muon_pfRelIso03_chg",	( ntupleRawTree.b5_Muon_pfRelIso03_chg ),"b5_Muon_pfRelIso03_chg[b5_nMuon]/F");
	outTree->Branch("b5_Muon_pfRelIso04_all",	( ntupleRawTree.b5_Muon_pfRelIso04_all ),"b5_Muon_pfRelIso04_all[b5_nMuon]/F");
	outTree->Branch("b5_Muon_phi",	( ntupleRawTree.b5_Muon_phi ),"b5_Muon_phi[b5_nMuon]/F");
	outTree->Branch("b5_Muon_pt",	( ntupleRawTree.b5_Muon_pt ),"b5_Muon_pt[b5_nMuon]/F");
	outTree->Branch("b5_Muon_ptErr",	( ntupleRawTree.b5_Muon_ptErr ),"b5_Muon_ptErr[b5_nMuon]/F");
	outTree->Branch("b5_Muon_segmentComp",	( ntupleRawTree.b5_Muon_segmentComp ),"b5_Muon_segmentComp[b5_nMuon]/F");
	outTree->Branch("b5_Muon_sip3d",	( ntupleRawTree.b5_Muon_sip3d ),"b5_Muon_sip3d[b5_nMuon]/F");
	outTree->Branch("b5_Muon_softMva",	( ntupleRawTree.b5_Muon_softMva ),"b5_Muon_softMva[b5_nMuon]/F");
	outTree->Branch("b5_Muon_tkRelIso",	( ntupleRawTree.b5_Muon_tkRelIso ),"b5_Muon_tkRelIso[b5_nMuon]/F");
	outTree->Branch("b5_Muon_tunepRelPt",	( ntupleRawTree.b5_Muon_tunepRelPt ),"b5_Muon_tunepRelPt[b5_nMuon]/F");
	outTree->Branch("b5_Muon_mvaLowPt",	( ntupleRawTree.b5_Muon_mvaLowPt ),"b5_Muon_mvaLowPt[b5_nMuon]/F");
	outTree->Branch("b5_Muon_mvaTTH",	( ntupleRawTree.b5_Muon_mvaTTH ),"b5_Muon_mvaTTH[b5_nMuon]/F");
	outTree->Branch("b5_Muon_charge",	( ntupleRawTree.b5_Muon_charge ),"b5_Muon_charge[b5_nMuon]/I");
	outTree->Branch("b5_Muon_jetIdx",	( ntupleRawTree.b5_Muon_jetIdx ),"b5_Muon_jetIdx[b5_nMuon]/I");
	outTree->Branch("b5_Muon_nStations",	( ntupleRawTree.b5_Muon_nStations ),"b5_Muon_nStations[b5_nMuon]/I");
	outTree->Branch("b5_Muon_nTrackerLayers",	( ntupleRawTree.b5_Muon_nTrackerLayers ),"b5_Muon_nTrackerLayers[b5_nMuon]/I");
	outTree->Branch("b5_Muon_pdgId",	( ntupleRawTree.b5_Muon_pdgId ),"b5_Muon_pdgId[b5_nMuon]/I");
	outTree->Branch("b5_Muon_tightCharge",	( ntupleRawTree.b5_Muon_tightCharge ),"b5_Muon_tightCharge[b5_nMuon]/I");
	outTree->Branch("b5_Muon_fsrPhotonIdx",	( ntupleRawTree.b5_Muon_fsrPhotonIdx ),"b5_Muon_fsrPhotonIdx[b5_nMuon]/I");
	outTree->Branch("b5_Muon_highPtId",	( ntupleRawTree.b5_Muon_highPtId ),"b5_Muon_highPtId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_highPurity",	( ntupleRawTree.b5_Muon_highPurity ),"b5_Muon_highPurity[b5_nMuon]/O");
	outTree->Branch("b5_Muon_inTimeMuon",	( ntupleRawTree.b5_Muon_inTimeMuon ),"b5_Muon_inTimeMuon[b5_nMuon]/O");
	outTree->Branch("b5_Muon_isGlobal",	( ntupleRawTree.b5_Muon_isGlobal ),"b5_Muon_isGlobal[b5_nMuon]/O");
	outTree->Branch("b5_Muon_isPFcand",	( ntupleRawTree.b5_Muon_isPFcand ),"b5_Muon_isPFcand[b5_nMuon]/O");
	outTree->Branch("b5_Muon_isTracker",	( ntupleRawTree.b5_Muon_isTracker ),"b5_Muon_isTracker[b5_nMuon]/O");
	outTree->Branch("b5_Muon_jetNDauCharged",	( ntupleRawTree.b5_Muon_jetNDauCharged ),"b5_Muon_jetNDauCharged[b5_nMuon]/b");
	outTree->Branch("b5_Muon_looseId",	( ntupleRawTree.b5_Muon_looseId ),"b5_Muon_looseId[b5_nMuon]/O");
	outTree->Branch("b5_Muon_mediumId",	( ntupleRawTree.b5_Muon_mediumId ),"b5_Muon_mediumId[b5_nMuon]/O");
	outTree->Branch("b5_Muon_mediumPromptId",	( ntupleRawTree.b5_Muon_mediumPromptId ),"b5_Muon_mediumPromptId[b5_nMuon]/O");
	outTree->Branch("b5_Muon_miniIsoId",	( ntupleRawTree.b5_Muon_miniIsoId ),"b5_Muon_miniIsoId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_multiIsoId",	( ntupleRawTree.b5_Muon_multiIsoId ),"b5_Muon_multiIsoId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_mvaId",	( ntupleRawTree.b5_Muon_mvaId ),"b5_Muon_mvaId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_mvaLowPtId",	( ntupleRawTree.b5_Muon_mvaLowPtId ),"b5_Muon_mvaLowPtId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_pfIsoId",	( ntupleRawTree.b5_Muon_pfIsoId ),"b5_Muon_pfIsoId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_puppiIsoId",	( ntupleRawTree.b5_Muon_puppiIsoId ),"b5_Muon_puppiIsoId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_softId",	( ntupleRawTree.b5_Muon_softId ),"b5_Muon_softId[b5_nMuon]/O");
	outTree->Branch("b5_Muon_softMvaId",	( ntupleRawTree.b5_Muon_softMvaId ),"b5_Muon_softMvaId[b5_nMuon]/O");
	outTree->Branch("b5_Muon_tightId",	( ntupleRawTree.b5_Muon_tightId ),"b5_Muon_tightId[b5_nMuon]/O");
	outTree->Branch("b5_Muon_tkIsoId",	( ntupleRawTree.b5_Muon_tkIsoId ),"b5_Muon_tkIsoId[b5_nMuon]/b");
	outTree->Branch("b5_Muon_triggerIdLoose",	( ntupleRawTree.b5_Muon_triggerIdLoose ),"b5_Muon_triggerIdLoose[b5_nMuon]/O");
	outTree->Branch("b5_nTrigObj",&	( ntupleRawTree.b5_nTrigObj ));
	outTree->Branch("b5_TrigObj_pt",	( ntupleRawTree.b5_TrigObj_pt ),"b5_TrigObj_pt[b5_nTrigObj]/F");
	outTree->Branch("b5_TrigObj_eta",	( ntupleRawTree.b5_TrigObj_eta ),"b5_TrigObj_eta[b5_nTrigObj]/F");
	outTree->Branch("b5_TrigObj_phi",	( ntupleRawTree.b5_TrigObj_phi ),"b5_TrigObj_phi[b5_nTrigObj]/F");
	outTree->Branch("b5_TrigObj_l1pt",	( ntupleRawTree.b5_TrigObj_l1pt ),"b5_TrigObj_l1pt[b5_nTrigObj]/F");
	outTree->Branch("b5_TrigObj_l1pt_2",	( ntupleRawTree.b5_TrigObj_l1pt_2 ),"b5_TrigObj_l1pt_2[b5_nTrigObj]/F");
	outTree->Branch("b5_TrigObj_l2pt",	( ntupleRawTree.b5_TrigObj_l2pt ),"b5_TrigObj_l2pt[b5_nTrigObj]/F");
	outTree->Branch("b5_TrigObj_id",	( ntupleRawTree.b5_TrigObj_id ),"b5_TrigObj_id[b5_nTrigObj]/I");
	outTree->Branch("b5_TrigObj_l1iso",	( ntupleRawTree.b5_TrigObj_l1iso ),"b5_TrigObj_l1iso[b5_nTrigObj]/I");
	outTree->Branch("b5_TrigObj_l1charge",	( ntupleRawTree.b5_TrigObj_l1charge ),"b5_TrigObj_l1charge[b5_nTrigObj]/I");
	outTree->Branch("b5_TrigObj_filterBits",	( ntupleRawTree.b5_TrigObj_filterBits ),"b5_TrigObj_filterBits[b5_nTrigObj]/I");
	outTree->Branch("b5_nOtherPV",&	( ntupleRawTree.b5_nOtherPV ));
	outTree->Branch("b5_OtherPV_z",	( ntupleRawTree.b5_OtherPV_z ),"b5_OtherPV_z[b5_nOtherPV]/F");
	outTree->Branch("b5_PV_ndof",&	( ntupleRawTree.b5_PV_ndof ));
	outTree->Branch("b5_PV_x",&	( ntupleRawTree.b5_PV_x ));
	outTree->Branch("b5_PV_y",&	( ntupleRawTree.b5_PV_y ));
	outTree->Branch("b5_PV_z",&	( ntupleRawTree.b5_PV_z ));
	outTree->Branch("b5_PV_chi2",&	( ntupleRawTree.b5_PV_chi2 ));
	outTree->Branch("b5_PV_score",&	( ntupleRawTree.b5_PV_score ));
	outTree->Branch("b5_PV_npvs",&	( ntupleRawTree.b5_PV_npvs ));
	outTree->Branch("b5_PV_npvsGood",&	( ntupleRawTree.b5_PV_npvsGood ));
	outTree->Branch("b5_nSV",&	( ntupleRawTree.b5_nSV ));
	outTree->Branch("b5_SV_dlen",	( ntupleRawTree.b5_SV_dlen ),"b5_SV_dlen[b5_nSV]/F");
	outTree->Branch("b5_SV_dlenSig",	( ntupleRawTree.b5_SV_dlenSig ),"b5_SV_dlenSig[b5_nSV]/F");
	outTree->Branch("b5_SV_dxy",	( ntupleRawTree.b5_SV_dxy ),"b5_SV_dxy[b5_nSV]/F");
	outTree->Branch("b5_SV_dxySig",	( ntupleRawTree.b5_SV_dxySig ),"b5_SV_dxySig[b5_nSV]/F");
	outTree->Branch("b5_SV_chi2",	( ntupleRawTree.b5_SV_chi2 ),"b5_SV_chi2[b5_nSV]/F");
	outTree->Branch("b5_SV_eta",	( ntupleRawTree.b5_SV_eta ),"b5_SV_eta[b5_nSV]/F");
	outTree->Branch("b5_SV_mass",	( ntupleRawTree.b5_SV_mass ),"b5_SV_mass[b5_nSV]/F");
	outTree->Branch("b5_SV_ndof",	( ntupleRawTree.b5_SV_ndof ),"b5_SV_ndof[b5_nSV]/F");
	outTree->Branch("b5_SV_phi",	( ntupleRawTree.b5_SV_phi ),"b5_SV_phi[b5_nSV]/F");
	outTree->Branch("b5_SV_pt",	( ntupleRawTree.b5_SV_pt ),"b5_SV_pt[b5_nSV]/F");
	outTree->Branch("b5_SV_x",	( ntupleRawTree.b5_SV_x ),"b5_SV_x[b5_nSV]/F");
	outTree->Branch("b5_SV_y",	( ntupleRawTree.b5_SV_y ),"b5_SV_y[b5_nSV]/F");
	outTree->Branch("b5_SV_z",	( ntupleRawTree.b5_SV_z ),"b5_SV_z[b5_nSV]/F");
	outTree->Branch("b5_SV_ntracks",	( ntupleRawTree.b5_SV_ntracks ),"b5_SV_ntracks[b5_nSV]/b");
	outTree->Branch("b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",&	( ntupleRawTree.b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 ));
	outTree->Branch("b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",&	( ntupleRawTree.b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 ));
	outTree->Branch("b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4",&	( ntupleRawTree.b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4 ));
	outTree->Branch("b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4",&	( ntupleRawTree.b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4 ));
	outTree->Branch("b5_L1_DoubleMu4_SQ_OS_dR_Max1p2",&	( ntupleRawTree.b5_L1_DoubleMu4_SQ_OS_dR_Max1p2 ));
	outTree->Branch("b5_L1_DoubleMu4p5_SQ_OS",&	( ntupleRawTree.b5_L1_DoubleMu4p5_SQ_OS ));
	outTree->Branch("b5_HLT_DoubleMu4_3_Bs",&	( ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ));
	outTree->Branch("b5_HLT_DoubleMu4_3_Jpsi",&	( ntupleRawTree.b5_HLT_DoubleMu4_3_Jpsi ));
	outTree->Branch("b5_HLT_DoubleMu4_JpsiTrk_Displaced",&	( ntupleRawTree.b5_HLT_DoubleMu4_JpsiTrk_Displaced ));
	outTree->Branch("b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",&	( ntupleRawTree.b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced ));
	outTree->Branch("bG_run",&	( ntupleRawTree.bG_run ));
	outTree->Branch("bG_event",&	( ntupleRawTree.bG_event ));
	outTree->Branch("bG_lumis",&	( ntupleRawTree.bG_lumis ));
	outTree->Branch("bG_isData",&	( ntupleRawTree.bG_isData ));
	outTree->Branch("bG_nPrimaryVertex",&	( ntupleRawTree.bG_nPrimaryVertex ));
	outTree->Branch("bG_primaryVertex_isFake",	( ntupleRawTree.bG_primaryVertex_isFake ),"bG_primaryVertex_isFake[bG_nPrimaryVertex]/O");
	outTree->Branch("bG_primaryVertex_x",	( ntupleRawTree.bG_primaryVertex_x ),"bG_primaryVertex_x[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_y",	( ntupleRawTree.bG_primaryVertex_y ),"bG_primaryVertex_y[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_z",	( ntupleRawTree.bG_primaryVertex_z ),"bG_primaryVertex_z[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_t",	( ntupleRawTree.bG_primaryVertex_t ),"bG_primaryVertex_t[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_covXX",	( ntupleRawTree.bG_primaryVertex_covXX ),"bG_primaryVertex_covXX[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_covXY",	( ntupleRawTree.bG_primaryVertex_covXY ),"bG_primaryVertex_covXY[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_covXZ",	( ntupleRawTree.bG_primaryVertex_covXZ ),"bG_primaryVertex_covXZ[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_covYY",	( ntupleRawTree.bG_primaryVertex_covYY ),"bG_primaryVertex_covYY[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_covYZ",	( ntupleRawTree.bG_primaryVertex_covYZ ),"bG_primaryVertex_covYZ[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_covZZ",	( ntupleRawTree.bG_primaryVertex_covZZ ),"bG_primaryVertex_covZZ[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_x_error",	( ntupleRawTree.bG_primaryVertex_x_error ),"bG_primaryVertex_x_error[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_y_error",	( ntupleRawTree.bG_primaryVertex_y_error ),"bG_primaryVertex_y_error[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_z_error",	( ntupleRawTree.bG_primaryVertex_z_error ),"bG_primaryVertex_z_error[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_t_error",	( ntupleRawTree.bG_primaryVertex_t_error ),"bG_primaryVertex_t_error[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_ntracks",	( ntupleRawTree.bG_primaryVertex_ntracks ),"bG_primaryVertex_ntracks[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_ndof",	( ntupleRawTree.bG_primaryVertex_ndof ),"bG_primaryVertex_ndof[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_chi2",	( ntupleRawTree.bG_primaryVertex_chi2 ),"bG_primaryVertex_chi2[bG_nPrimaryVertex]/D");
	outTree->Branch("bG_primaryVertex_normalizedChi2",	( ntupleRawTree.bG_primaryVertex_normalizedChi2 ),"bG_primaryVertex_normalizedChi2[bG_nPrimaryVertex]/D");
	//outTree->Branch("bG_nPho",&	( ntupleRawTree.bG_nPho ));
	//outTree->Branch("bG_phoE",	( ntupleRawTree.bG_phoE ),"bG_phoE[bG_nPho]/F");
	//outTree->Branch("bG_phoEt",	( ntupleRawTree.bG_phoEt ),"bG_phoEt[bG_nPho]/F");
	//outTree->Branch("bG_phoEta",	( ntupleRawTree.bG_phoEta ),"bG_phoEta[bG_nPho]/F");
	//outTree->Branch("bG_phoPhi",	( ntupleRawTree.bG_phoPhi ),"bG_phoPhi[bG_nPho]/F");
	//outTree->Branch("bG_phoSCE",	( ntupleRawTree.bG_phoSCE ),"bG_phoSCE[bG_nPho]/F");
	//outTree->Branch("bG_phoSCEt",	( ntupleRawTree.bG_phoSCEt ),"bG_phoSCEt[bG_nPho]/F");
	//outTree->Branch("bG_phoSCRawE",	( ntupleRawTree.bG_phoSCRawE ),"bG_phoSCRawE[bG_nPho]/F");
	//outTree->Branch("bG_phoESEnP1",	( ntupleRawTree.bG_phoESEnP1 ),"bG_phoESEnP1[bG_nPho]/F");
	//outTree->Branch("bG_phoESEnP2",	( ntupleRawTree.bG_phoESEnP2 ),"bG_phoESEnP2[bG_nPho]/F");
	//outTree->Branch("bG_phoSCEta",	( ntupleRawTree.bG_phoSCEta ),"bG_phoSCEta[bG_nPho]/F");
	//outTree->Branch("bG_phoSCPhi",	( ntupleRawTree.bG_phoSCPhi ),"bG_phoSCPhi[bG_nPho]/F");
	//outTree->Branch("bG_phoSCEtaWidth",	( ntupleRawTree.bG_phoSCEtaWidth ),"bG_phoSCEtaWidth[bG_nPho]/F");
	//outTree->Branch("bG_phoSCPhiWidth",	( ntupleRawTree.bG_phoSCPhiWidth ),"bG_phoSCPhiWidth[bG_nPho]/F");
	//outTree->Branch("bG_phoSCBrem",	( ntupleRawTree.bG_phoSCBrem ),"bG_phoSCBrem[bG_nPho]/F");
	//outTree->Branch("bG_phohasPixelSeed",	( ntupleRawTree.bG_phohasPixelSeed ),"bG_phohasPixelSeed[bG_nPho]/I");
	//outTree->Branch("bG_phoR9",	( ntupleRawTree.bG_phoR9 ),"bG_phoR9[bG_nPho]/F");
	//outTree->Branch("bG_phoHoverE",	( ntupleRawTree.bG_phoHoverE ),"bG_phoHoverE[bG_nPho]/F");
	//outTree->Branch("bG_phoESEffSigmaRR",	( ntupleRawTree.bG_phoESEffSigmaRR ),"bG_phoESEffSigmaRR[bG_nPho]/F");
	//outTree->Branch("bG_phoSigmaIEtaIEtaFull5x5",	( ntupleRawTree.bG_phoSigmaIEtaIEtaFull5x5 ),"bG_phoSigmaIEtaIEtaFull5x5[bG_nPho]/F");
	//outTree->Branch("bG_phoSigmaIEtaIPhiFull5x5",	( ntupleRawTree.bG_phoSigmaIEtaIPhiFull5x5 ),"bG_phoSigmaIEtaIPhiFull5x5[bG_nPho]/F");
	//outTree->Branch("bG_phoSigmaIPhiIPhiFull5x5",	( ntupleRawTree.bG_phoSigmaIPhiIPhiFull5x5 ),"bG_phoSigmaIPhiIPhiFull5x5[bG_nPho]/F");
	//outTree->Branch("bG_phoE2x2Full5x5",	( ntupleRawTree.bG_phoE2x2Full5x5 ),"bG_phoE2x2Full5x5[bG_nPho]/F");
	//outTree->Branch("bG_phoE5x5Full5x5",	( ntupleRawTree.bG_phoE5x5Full5x5 ),"bG_phoE5x5Full5x5[bG_nPho]/F");
	//outTree->Branch("bG_phoR9Full5x5",	( ntupleRawTree.bG_phoR9Full5x5 ),"bG_phoR9Full5x5[bG_nPho]/F");
	//outTree->Branch("bG_phoPFChIso",	( ntupleRawTree.bG_phoPFChIso ),"bG_phoPFChIso[bG_nPho]/F");
	//outTree->Branch("bG_phoPFPhoIso",	( ntupleRawTree.bG_phoPFPhoIso ),"bG_phoPFPhoIso[bG_nPho]/F");
	//outTree->Branch("bG_phoPFNeuIso",	( ntupleRawTree.bG_phoPFNeuIso ),"bG_phoPFNeuIso[bG_nPho]/F");
	//outTree->Branch("bG_phoEcalPFClusterIso",	( ntupleRawTree.bG_phoEcalPFClusterIso ),"bG_phoEcalPFClusterIso[bG_nPho]/F");
	//outTree->Branch("bG_phoHcalPFClusterIso",	( ntupleRawTree.bG_phoHcalPFClusterIso ),"bG_phoHcalPFClusterIso[bG_nPho]/F");
	//outTree->Branch("bG_phoSeedTime",	( ntupleRawTree.bG_phoSeedTime ),"bG_phoSeedTime[bG_nPho]/F");
	//outTree->Branch("bG_phoSeedEnergy",	( ntupleRawTree.bG_phoSeedEnergy ),"bG_phoSeedEnergy[bG_nPho]/F");
	//outTree->Branch("bG_phoMIPTotEnergy",	( ntupleRawTree.bG_phoMIPTotEnergy ),"bG_phoMIPTotEnergy[bG_nPho]/F");
	//outTree->Branch("bG_phoMIPChi2",	( ntupleRawTree.bG_phoMIPChi2 ),"bG_phoMIPChi2[bG_nPho]/F");
	//outTree->Branch("bG_phoMIPSlope",	( ntupleRawTree.bG_phoMIPSlope ),"bG_phoMIPSlope[bG_nPho]/F");
	//outTree->Branch("bG_phoMIPIntercept",	( ntupleRawTree.bG_phoMIPIntercept ),"bG_phoMIPIntercept[bG_nPho]/F");
	//outTree->Branch("bG_phoMIPNhitCone",	( ntupleRawTree.bG_phoMIPNhitCone ),"bG_phoMIPNhitCone[bG_nPho]/F");
	//outTree->Branch("bG_phoMIPIsHalo",	( ntupleRawTree.bG_phoMIPIsHalo ),"bG_phoMIPIsHalo[bG_nPho]/F");
	//outTree->Branch("bG_nPFPho",&	( ntupleRawTree.bG_nPFPho ));
	//outTree->Branch("bG_phoPFE",	( ntupleRawTree.bG_phoPFE ),"bG_phoPFE[bG_nPFPho]/F");
	//outTree->Branch("bG_phoPFEt",	( ntupleRawTree.bG_phoPFEt ),"bG_phoPFEt[bG_nPFPho]/F");
	//outTree->Branch("bG_phoPFEta",	( ntupleRawTree.bG_phoPFEta ),"bG_phoPFEta[bG_nPFPho]/F");
	//outTree->Branch("bG_phoPFPhi",	( ntupleRawTree.bG_phoPFPhi ),"bG_phoPFPhi[bG_nPFPho]/F");
	outTree->Branch("bG_nSC",&	( ntupleRawTree.bG_nSC ));
	outTree->Branch("bG_scE",	( ntupleRawTree.bG_scE ),"bG_scE[bG_nSC]/F");
	outTree->Branch("bG_scEt",	( ntupleRawTree.bG_scEt ),"bG_scEt[bG_nSC]/F");
	outTree->Branch("bG_scRawE",	( ntupleRawTree.bG_scRawE ),"bG_scRawE[bG_nSC]/F");
	outTree->Branch("bG_scEta",	( ntupleRawTree.bG_scEta ),"bG_scEta[bG_nSC]/F");
	outTree->Branch("bG_scPhi",	( ntupleRawTree.bG_scPhi ),"bG_scPhi[bG_nSC]/F");
	outTree->Branch("bG_scX",	( ntupleRawTree.bG_scX ),"bG_scX[bG_nSC]/F");
	outTree->Branch("bG_scY",	( ntupleRawTree.bG_scY ),"bG_scY[bG_nSC]/F");
	outTree->Branch("bG_scZ",	( ntupleRawTree.bG_scZ ),"bG_scZ[bG_nSC]/F");
	outTree->Branch("bG_scEtaWidth",	( ntupleRawTree.bG_scEtaWidth ),"bG_scEtaWidth[bG_nSC]/F");
	outTree->Branch("bG_scPhiWidth",	( ntupleRawTree.bG_scPhiWidth ),"bG_scPhiWidth[bG_nSC]/F");
	outTree->Branch("bG_scRawEt",	( ntupleRawTree.bG_scRawEt ),"bG_scRawEt[bG_nSC]/F");
	outTree->Branch("bG_scMinDrWithGsfElectornSC_",	( ntupleRawTree.bG_scMinDrWithGsfElectornSC_ ),"bG_scMinDrWithGsfElectornSC_[bG_nSC]/F");
	outTree->Branch("bG_scFoundGsfMatch_",	( ntupleRawTree.bG_scFoundGsfMatch_ ),"bG_scFoundGsfMatch_[bG_nSC]/O");
	outTree->Branch("bG_scE5x5",	( ntupleRawTree.bG_scE5x5 ),"bG_scE5x5[bG_nSC]/F");
	outTree->Branch("bG_scE2x2Ratio",	( ntupleRawTree.bG_scE2x2Ratio ),"bG_scE2x2Ratio[bG_nSC]/F");
	outTree->Branch("bG_scE3x3Ratio",	( ntupleRawTree.bG_scE3x3Ratio ),"bG_scE3x3Ratio[bG_nSC]/F");
	outTree->Branch("bG_scEMaxRatio",	( ntupleRawTree.bG_scEMaxRatio ),"bG_scEMaxRatio[bG_nSC]/F");
	outTree->Branch("bG_scE2ndRatio",	( ntupleRawTree.bG_scE2ndRatio ),"bG_scE2ndRatio[bG_nSC]/F");
	outTree->Branch("bG_scETopRatio",	( ntupleRawTree.bG_scETopRatio ),"bG_scETopRatio[bG_nSC]/F");
	outTree->Branch("bG_scERightRatio",	( ntupleRawTree.bG_scERightRatio ),"bG_scERightRatio[bG_nSC]/F");
	outTree->Branch("bG_scEBottomRatio",	( ntupleRawTree.bG_scEBottomRatio ),"bG_scEBottomRatio[bG_nSC]/F");
	outTree->Branch("bG_scELeftRatio",	( ntupleRawTree.bG_scELeftRatio ),"bG_scELeftRatio[bG_nSC]/F");
	outTree->Branch("bG_scE2x5MaxRatio",	( ntupleRawTree.bG_scE2x5MaxRatio ),"bG_scE2x5MaxRatio[bG_nSC]/F");
	outTree->Branch("bG_scE2x5TopRatio",	( ntupleRawTree.bG_scE2x5TopRatio ),"bG_scE2x5TopRatio[bG_nSC]/F");
	outTree->Branch("bG_scE2x5RightRatio",	( ntupleRawTree.bG_scE2x5RightRatio ),"bG_scE2x5RightRatio[bG_nSC]/F");
	outTree->Branch("bG_scE2x5BottomRatio",	( ntupleRawTree.bG_scE2x5BottomRatio ),"bG_scE2x5BottomRatio[bG_nSC]/F");
	outTree->Branch("bG_scE2x5LeftRatio",	( ntupleRawTree.bG_scE2x5LeftRatio ),"bG_scE2x5LeftRatio[bG_nSC]/F");
	outTree->Branch("bG_scSwissCross",	( ntupleRawTree.bG_scSwissCross ),"bG_scSwissCross[bG_nSC]/F");
	outTree->Branch("bG_scR9",	( ntupleRawTree.bG_scR9 ),"bG_scR9[bG_nSC]/F");
	outTree->Branch("bG_scSigmaIetaIeta",	( ntupleRawTree.bG_scSigmaIetaIeta ),"bG_scSigmaIetaIeta[bG_nSC]/F");
	outTree->Branch("bG_scSigmaIetaIphi",	( ntupleRawTree.bG_scSigmaIetaIphi ),"bG_scSigmaIetaIphi[bG_nSC]/F");
	outTree->Branch("bG_scSigmaIphiIphi",	( ntupleRawTree.bG_scSigmaIphiIphi ),"bG_scSigmaIphiIphi[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e5x5",	( ntupleRawTree.bG_scFull5x5_e5x5 ),"bG_scFull5x5_e5x5[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e2x2Ratio",	( ntupleRawTree.bG_scFull5x5_e2x2Ratio ),"bG_scFull5x5_e2x2Ratio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e3x3Ratio",	( ntupleRawTree.bG_scFull5x5_e3x3Ratio ),"bG_scFull5x5_e3x3Ratio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_eMaxRatio",	( ntupleRawTree.bG_scFull5x5_eMaxRatio ),"bG_scFull5x5_eMaxRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e2ndRatio",	( ntupleRawTree.bG_scFull5x5_e2ndRatio ),"bG_scFull5x5_e2ndRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_eTopRatio",	( ntupleRawTree.bG_scFull5x5_eTopRatio ),"bG_scFull5x5_eTopRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_eRightRatio",	( ntupleRawTree.bG_scFull5x5_eRightRatio ),"bG_scFull5x5_eRightRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_eBottomRatio",	( ntupleRawTree.bG_scFull5x5_eBottomRatio ),"bG_scFull5x5_eBottomRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_eLeftRatio",	( ntupleRawTree.bG_scFull5x5_eLeftRatio ),"bG_scFull5x5_eLeftRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e2x5MaxRatio",	( ntupleRawTree.bG_scFull5x5_e2x5MaxRatio ),"bG_scFull5x5_e2x5MaxRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e2x5TopRatio",	( ntupleRawTree.bG_scFull5x5_e2x5TopRatio ),"bG_scFull5x5_e2x5TopRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e2x5RightRatio",	( ntupleRawTree.bG_scFull5x5_e2x5RightRatio ),"bG_scFull5x5_e2x5RightRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e2x5BottomRatio",	( ntupleRawTree.bG_scFull5x5_e2x5BottomRatio ),"bG_scFull5x5_e2x5BottomRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_e2x5LeftRatio",	( ntupleRawTree.bG_scFull5x5_e2x5LeftRatio ),"bG_scFull5x5_e2x5LeftRatio[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_swissCross",	( ntupleRawTree.bG_scFull5x5_swissCross ),"bG_scFull5x5_swissCross[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_r9",	( ntupleRawTree.bG_scFull5x5_r9 ),"bG_scFull5x5_r9[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_sigmaIetaIeta",	( ntupleRawTree.bG_scFull5x5_sigmaIetaIeta ),"bG_scFull5x5_sigmaIetaIeta[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_sigmaIetaIphi",	( ntupleRawTree.bG_scFull5x5_sigmaIetaIphi ),"bG_scFull5x5_sigmaIetaIphi[bG_nSC]/F");
	outTree->Branch("bG_scFull5x5_sigmaIphiIphi",	( ntupleRawTree.bG_scFull5x5_sigmaIphiIphi ),"bG_scFull5x5_sigmaIphiIphi[bG_nSC]/F");
	outTree->Branch("bG_scPFChIso1",	( ntupleRawTree.bG_scPFChIso1 ),"bG_scPFChIso1[bG_nSC]/F");
	outTree->Branch("bG_scPFChIso2",	( ntupleRawTree.bG_scPFChIso2 ),"bG_scPFChIso2[bG_nSC]/F");
	outTree->Branch("bG_scPFChIso3",	( ntupleRawTree.bG_scPFChIso3 ),"bG_scPFChIso3[bG_nSC]/F");
	outTree->Branch("bG_scPFChIso4",	( ntupleRawTree.bG_scPFChIso4 ),"bG_scPFChIso4[bG_nSC]/F");
	outTree->Branch("bG_scPFChIso5",	( ntupleRawTree.bG_scPFChIso5 ),"bG_scPFChIso5[bG_nSC]/F");
	outTree->Branch("bG_scPFPhoIso1",	( ntupleRawTree.bG_scPFPhoIso1 ),"bG_scPFPhoIso1[bG_nSC]/F");
	outTree->Branch("bG_scPFPhoIso2",	( ntupleRawTree.bG_scPFPhoIso2 ),"bG_scPFPhoIso2[bG_nSC]/F");
	outTree->Branch("bG_scPFPhoIso3",	( ntupleRawTree.bG_scPFPhoIso3 ),"bG_scPFPhoIso3[bG_nSC]/F");
	outTree->Branch("bG_scPFPhoIso4",	( ntupleRawTree.bG_scPFPhoIso4 ),"bG_scPFPhoIso4[bG_nSC]/F");
	outTree->Branch("bG_scPFPhoIso5",	( ntupleRawTree.bG_scPFPhoIso5 ),"bG_scPFPhoIso5[bG_nSC]/F");
	outTree->Branch("bG_scPFNeuIso1",	( ntupleRawTree.bG_scPFNeuIso1 ),"bG_scPFNeuIso1[bG_nSC]/F");
	outTree->Branch("bG_scPFNeuIso2",	( ntupleRawTree.bG_scPFNeuIso2 ),"bG_scPFNeuIso2[bG_nSC]/F");
	outTree->Branch("bG_scPFNeuIso3",	( ntupleRawTree.bG_scPFNeuIso3 ),"bG_scPFNeuIso3[bG_nSC]/F");
	outTree->Branch("bG_scPFNeuIso4",	( ntupleRawTree.bG_scPFNeuIso4 ),"bG_scPFNeuIso4[bG_nSC]/F");
	outTree->Branch("bG_scPFNeuIso5",	( ntupleRawTree.bG_scPFNeuIso5 ),"bG_scPFNeuIso5[bG_nSC]/F");
	//outTree->Branch("bG_nhcalRechit",&	( ntupleRawTree.bG_nhcalRechit ));
	//outTree->Branch("bG_hcalRechitIEta",	( ntupleRawTree.bG_hcalRechitIEta ),"bG_hcalRechitIEta[bG_nhcalRechit]/F");
	//outTree->Branch("bG_hcalRechitIPhi",	( ntupleRawTree.bG_hcalRechitIPhi ),"bG_hcalRechitIPhi[bG_nhcalRechit]/F");
    //outTree->Branch("bG_hcalRechitEnergy",	( ntupleRawTree.bG_hcalRechitEnergy ),"bG_hcalRechitEnergy[bG_nhcalRechit]/F");
	#ifdef __MCANALYSIS__
	outTree->Branch("b5_mm_gen_alpha_ip",	( ntupleRawTree.b5_mm_gen_alpha_ip ),"b5_mm_gen_alpha_ip[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_alpha_p_phi",	( ntupleRawTree.b5_mm_gen_alpha_p_phi ),"b5_mm_gen_alpha_p_phi[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_alpha_p_theta",	( ntupleRawTree.b5_mm_gen_alpha_p_theta ),"b5_mm_gen_alpha_p_theta[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_alpha_vtx",	( ntupleRawTree.b5_mm_gen_alpha_vtx ),"b5_mm_gen_alpha_vtx[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_doca",	( ntupleRawTree.b5_mm_gen_doca ),"b5_mm_gen_doca[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_l3d",	( ntupleRawTree.b5_mm_gen_l3d ),"b5_mm_gen_l3d[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_lxy",	( ntupleRawTree.b5_mm_gen_lxy ),"b5_mm_gen_lxy[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_mass",	( ntupleRawTree.b5_mm_gen_mass ),"b5_mm_gen_mass[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_mu1_pt",	( ntupleRawTree.b5_mm_gen_mu1_pt ),"b5_mm_gen_mu1_pt[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_mu2_pt",	( ntupleRawTree.b5_mm_gen_mu2_pt ),"b5_mm_gen_mu2_pt[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_prod_z",	( ntupleRawTree.b5_mm_gen_prod_z ),"b5_mm_gen_prod_z[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_pt",	( ntupleRawTree.b5_mm_gen_pt ),"b5_mm_gen_pt[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_tau",	( ntupleRawTree.b5_mm_gen_tau ),"b5_mm_gen_tau[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_vtx_x",	( ntupleRawTree.b5_mm_gen_vtx_x ),"b5_mm_gen_vtx_x[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_vtx_y",	( ntupleRawTree.b5_mm_gen_vtx_y ),"b5_mm_gen_vtx_y[b5_nmm]/F");
	outTree->Branch("b5_mm_gen_vtx_z",	( ntupleRawTree.b5_mm_gen_vtx_z ),"b5_mm_gen_vtx_z[b5_nmm]/F");
    outTree->Branch("b5_MuonId_simProdRho",	( ntupleRawTree.b5_MuonId_simProdRho ),"b5_MuonId_simProdRho[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_simProdZ",	( ntupleRawTree.b5_MuonId_simProdZ ),"b5_MuonId_simProdZ[b5_nMuonId]/F");
	outTree->Branch("b5_MuonId_simExtType",	( ntupleRawTree.b5_MuonId_simExtType ),"b5_MuonId_simExtType[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_simMotherPdgId",	( ntupleRawTree.b5_MuonId_simMotherPdgId ),"b5_MuonId_simMotherPdgId[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_simPdgId",	( ntupleRawTree.b5_MuonId_simPdgId ),"b5_MuonId_simPdgId[b5_nMuonId]/I");
	outTree->Branch("b5_MuonId_simType",	( ntupleRawTree.b5_MuonId_simType ),"b5_MuonId_simType[b5_nMuonId]/I");
	outTree->Branch("b5_nGenPart",&	( ntupleRawTree.b5_nGenPart ));
	outTree->Branch("b5_GenPart_eta",	( ntupleRawTree.b5_GenPart_eta ),"b5_GenPart_eta[b5_nGenPart]/F");
	outTree->Branch("b5_GenPart_mass",	( ntupleRawTree.b5_GenPart_mass ),"b5_GenPart_mass[b5_nGenPart]/F");
	outTree->Branch("b5_GenPart_phi",	( ntupleRawTree.b5_GenPart_phi ),"b5_GenPart_phi[b5_nGenPart]/F");
	outTree->Branch("b5_GenPart_pt",	( ntupleRawTree.b5_GenPart_pt ),"b5_GenPart_pt[b5_nGenPart]/F");
	outTree->Branch("b5_GenPart_genPartIdxMother",	( ntupleRawTree.b5_GenPart_genPartIdxMother ),"b5_GenPart_genPartIdxMother[b5_nGenPart]/I");
	outTree->Branch("b5_GenPart_pdgId",	( ntupleRawTree.b5_GenPart_pdgId ),"b5_GenPart_pdgId[b5_nGenPart]/I");
	outTree->Branch("b5_GenPart_status",	( ntupleRawTree.b5_GenPart_status ),"b5_GenPart_status[b5_nGenPart]/I");
	outTree->Branch("b5_GenPart_statusFlags",	( ntupleRawTree.b5_GenPart_statusFlags ),"b5_GenPart_statusFlags[b5_nGenPart]/I");
	outTree->Branch("b5_genTtbarId",&	( ntupleRawTree.b5_genTtbarId ));
    #endif


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
     th1fStore["dimu_docatrk"]       = new TH1F("dimu_dcaToCloseTrack","min(doca,0.03) to all pf tracks  #mu#mu Candidate" ,100 , 0.0 , 0.04  );
     th1fStore["dimu_NTrakClose"]       = new TH1F("dimu_NTrakClose","#Tracks doca(trk,sv) < 0.3 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_NTrakCloseSig1"]   = new TH1F("dimu_NTrakCloseSig1","#Tracks d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) > 1 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_NTrakCloseSig2"]   = new TH1F("dimu_NTrakCloseSig2","#Tracks d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) > 2 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_NTrakCloseSig3"]   = new TH1F("dimu_NTrakCloseSig3","#Tracks d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) > 3 for #mu#mu Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["dimu_Isolation"]  = new TH1F("dimu_Isolation","I(#mu#mu) for #mu#mu Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["dimu_mva"]     = new TH1F("dimu_mva","MVA score of dimu Candidate", 80 , -2.0  , 2.0  );

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
     th1fStore["dimuPass_docatrk"]       = new TH1F("dimuPass_dcaToCloseTrack","min(doca,0.03) to all pf tracks  #mu#mu Candidate" ,100 , 0.0 , 0.04  );
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
    
    th1fStore["dimu_dcaGammaDimuSVDrMax0p5"] = new TH1F("dimu_dcaGammaDimuSVDrMax0p5","DCA(#mu#mu,#gamma) for #mu#mu#gamma Candidate" ,200 , 0.0 , 2.00 );
    th1fStore["dimu_dcaGammaDimuSVDrMax0p7"] = new TH1F("dimu_dcaGammaDimuSVDrMax0p7","DCA(#mu#mu,#gamma) for #mu#mu#gamma Candidate" ,200 , 0.0 , 2.00 );
    th1fStore["dimu_dcaGammaDimuSVDrMax0p9"] = new TH1F("dimu_dcaGammaDimuSVDrMax0p9","DCA(#mu#mu,#gamma) for #mu#mu#gamma Candidate" ,200 , 0.0 , 2.00 );
    th1fStore["dimu_dcaGammaDimuSVDrMax1p1"] = new TH1F("dimu_dcaGammaDimuSVDrMax1p1","DCA(#mu#mu,#gamma) for #mu#mu#gamma Candidate" ,200 , 0.0 , 2.00 );
    th1fStore["dimu_dcaGammaDimuSVDrMax1p3"] = new TH1F("dimu_dcaGammaDimuSVDrMax1p3","DCA(#mu#mu,#gamma) for #mu#mu#gamma Candidate" ,200 , 0.0 , 2.00 );
 
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
     th1fStore["bmmg_docatrk"]       = new TH1F("bmmg_dcaToCloseTrack","min(doca,0.03) to all pf tracks  #mu#mu Candidate" ,100 , 0.0 , 0.04  );
     th1fStore["bmmg_NTrakClose"]       = new TH1F("bmmg_NTrakClose","#Tracks doca(trk,sv) < 0.3 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_NTrakCloseSig1"]= new TH1F("bmmg_NTrakCloseSig1","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 1 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_NTrakCloseSig2"]= new TH1F("bmmg_NTrakCloseSig2","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 2 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_NTrakCloseSig3"]= new TH1F("bmmg_NTrakCloseSig3","#trk d(trk,sv) < 0.3 & d(trk,sv)/#sigma_{d}(trk,sv) < 3 for #mu#mu#gamma Candidate" ,100 , -0.50 , 99.5  );
     th1fStore["bmmg_Isolation"]  = new TH1F("bmmg_Isolation","I(#mu#mu) for #mu#mu#gamma Candidate" ,120 , -0.10 , 1.1  );
     th1fStore["bmmg_dcaGammaDimuSV"] = new TH1F("bmmg_dcaGammaDimuSV","DCA(#mu#mu,#gamma) for #mu#mu#gamma Candidate" ,200 , 0.0 , 2.00 );

    // Global Event Hists
    th1fStore["dimuPass_bmmgCandidateMultiplicity"]   = new TH1F("dimuPass_bmmgCandidateMultiplicity","Multiplicity of BMMG Candidates Per dimuon",15, -0.5,14.5);
    th1fStore["dimuPass_Multiplicity"]       = new TH1F("dimuPass_Multiplicity","Multiplicity of dimu Candidates in an event",15, -0.5,14.5);
    th1fStore["bmmg_Multiplicity"]       = new TH1F("bmmg_Multiplicity","Multiplicity of bmmg Candidates in an event",15, -0.5,14.5);

    // MC : Gen Match Hists 
     th1fStore["gen_mumDeltaR"]   = new TH1F("gen_mumDeltaR","#Delta(#mu^-,#mu_{Gen}^-)", 200 , 0.0  , 2.0  );
     th1fStore["gen_mupDeltaR"]   = new TH1F("gen_mupDeltaR","#Delta(#mu^+,#mu_{Gen}^+)", 200 , 0.0  , 2.0  );
     th1fStore["gen_phoDeltaR"]   = new TH1F("gen_phoDeltaR","#Delta(#gamma,#gamma_{Gen})", 200 , 0.0  , 2.0  );
   // Processing Summary Hists 
     th1fStore["ProcessingSummary"] = new TH1F("processing_summary","",3,0.0,3.0);
     th1fStore["ProcessingSummary"]->SetCanExtend(TH1::kAllAxes);
     
     th1fStore["CutFlowSummary"] = new TH1F("cutFlow_summary","",3,0.0,3.0);
     th1fStore["ProcessingSummary"]->SetCanExtend(TH1::kAllAxes);
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
     th1fStore["dimu_docatrk"]       ->Fill(min(ntupleRawTree.b5_mm_docatrk[mumuIdx],Float_t(0.04)));
     th1fStore["dimu_NTrakClose"]       ->Fill(ntupleRawTree.b5_mm_closetrk[mumuIdx]);
     th1fStore["dimu_NTrakCloseSig1"]   ->Fill(ntupleRawTree.b5_mm_closetrks1[mumuIdx]);
     th1fStore["dimu_NTrakCloseSig2"]   ->Fill(ntupleRawTree.b5_mm_closetrks2[mumuIdx]);
     th1fStore["dimu_NTrakCloseSig3"]   ->Fill(ntupleRawTree.b5_mm_closetrks3[mumuIdx]);
     th1fStore["dimu_Isolation"]  ->Fill(ntupleRawTree.b5_mm_iso[mumuIdx]);
     th1fStore["dimu_mva"]   ->Fill(storageArrayDouble[ mumuIdx + candidateMapInt["dimuonMVAScore"]  ]); 
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
     th1fStore["dimuPass_docatrk"]       ->Fill(min(ntupleRawTree.b5_mm_docatrk[mumuIdx],Float_t(0.04)));
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
            
            if(dr < 0.5 )  
            { 
                array[0]++;
		        th1fStore["dimu_dcaGammaDimuSVDrMax0p5"]->Fill(getDCAGammaToDimuVertex(mumuIdx,i));
			}
            if(dr < 0.7 )  
            { 
                array[1]++;
		        th1fStore["dimu_dcaGammaDimuSVDrMax0p7"]->Fill(getDCAGammaToDimuVertex(mumuIdx,i));
			}
            if(dr < 0.9 )  
            { 
                array[2]++;
		        th1fStore["dimu_dcaGammaDimuSVDrMax0p9"]->Fill(getDCAGammaToDimuVertex(mumuIdx,i));
			}
            if(dr < 1.1 )  
            { 
                array[3]++;
		        th1fStore["dimu_dcaGammaDimuSVDrMax1p1"]->Fill(getDCAGammaToDimuVertex(mumuIdx,i));
			}
            if(dr < 1.3 )  
            { 
                array[4]++;
		        th1fStore["dimu_dcaGammaDimuSVDrMax1p3"]->Fill(getDCAGammaToDimuVertex(mumuIdx,i));
			}
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
    th1fStore["bmmg_docatrk"]       ->Fill(min(ntupleRawTree.b5_mm_docatrk[mumuIdx],Float_t(0.04)));
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
    th1fStore["bmmg_dcaGammaDimuSV"]        ->Fill(getDCAGammaToDimuVertex(mumuIdx,phoSCIdx));

}

void BMMGAnalysis::fill_globalEventHists()
{
    th1fStore["bmmg_Multiplicity"]->Fill(nBMMGCandidates);
    
    
    for(int i=0;i<nDiMuCandidates;i++)
    {
        if(i>=NDIMU_MAX) break;
        th1fStore["dimuPass_bmmgCandidateMultiplicity"]->Fill(i+storageArrayInt[candidateMapInt["nBMMGCandidatesPerDimu"]  ] );
    }
}



Int_t BMMGAnalysis::doMuonSelection(Int_t muIdx, bool isLead)
{
    Int_t rslt(0);

    rslt++;
    if(  ntupleRawTree.b5_Muon_pt[muIdx] < minMuonPt )  
		 return cutFlowOffsets["muonSelection"] + rslt ;
	rslt++;
    if(  abs(ntupleRawTree.b5_Muon_eta[muIdx]) > maxMuonEta )  
		 return cutFlowOffsets["muonSelection"] + rslt ;
	rslt++;
    if( ( not ntupleRawTree.b5_Muon_isTracker[muIdx] ) and muonHasToBeTracker) 
		 return cutFlowOffsets["muonSelection"] + rslt ;
	rslt++;
    if( ( not ntupleRawTree.b5_Muon_isGlobal[muIdx]) and muonHasToBeGlobal )
		 return cutFlowOffsets["muonSelection"] + rslt ;
	rslt++;
    if( (not ntupleRawTree.b5_Muon_looseId[muIdx] ) and muonHasToBeLoose )
		 return cutFlowOffsets["muonSelection"] + rslt ;
	rslt++;
    if( (not ntupleRawTree.b5_MuonId_highPurity[muIdx]) and muonHasToBeHighPurity )  
		 return cutFlowOffsets["muonSelection"] + rslt ;
	rslt++;
    if(ntupleRawTree.b5_MuonId_newSoftMuonMva[muIdx] < BDTWorkingPoint ) 
		 return cutFlowOffsets["muonSelection"] + rslt ;
		 
         return  0;
}

Int_t BMMGAnalysis::doDimuonSelection(Int_t mumuIdx)
{
    Int_t rslt(0);
	rslt++;   
    //std::cout<<" "<<ntupleRawTree.b5_mm_m1iso[mumuIdx] <<" < "<<minLeadMuIsolation<<"\n";
    if( ntupleRawTree.b5_mm_m1iso[mumuIdx]  < minLeadMuIsolation ) 
		 return cutFlowOffsets["diMuonSelection"] + rslt ;
	rslt++;
    //std::cout<<" "<<ntupleRawTree.b5_mm_m2iso[mumuIdx] <<" < "<<minLeadMuIsolation<<"\n";
    if( ntupleRawTree.b5_mm_m2iso[mumuIdx]  < minSubLeadMuIsolation ) 
		 return cutFlowOffsets["diMuonSelection"] + rslt ;

    diMuLV.SetPtEtaPhiM(     ntupleRawTree.b5_mm_kin_pt[mumuIdx],   \
                             ntupleRawTree.b5_mm_kin_eta[mumuIdx],  \
                             ntupleRawTree.b5_mm_kin_phi[mumuIdx],  \
                             ntupleRawTree.b5_mm_kin_mass[mumuIdx]  );
    
    auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
               ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );

	rslt++;
    if (dr > maxMuMuDr ) 
		 return cutFlowOffsets["diMuonSelection"] + rslt ;
    
    auto diMuMass = diMuLV.M();
    bool selection=false;
    
    //std::cout<<"MuMu mass : "<<diMuMass<<"\n";
    for(int i=0;i<minDimuMass.size();i++)
    {
         if(diMuMass > minDimuMass[i] and diMuMass < maxDimuMass[i]) selection=true;
         if(selection) break;
    }
	rslt++;   
    if( not selection ) 
		 return cutFlowOffsets["diMuonSelection"] + rslt ;
    //std::cout<<"MuMu mass : pass"<<"\n";

    
    return 0;
           
}

Int_t BMMGAnalysis::doVertexSelection(Int_t mumuIdx)
{
    Int_t rslt(0);
    
	rslt++;
    if( ntupleRawTree.b5_mm_doca[mumuIdx]  > muMuMaxDOCA ) 
		 return cutFlowOffsets["diMuonVertexSelection"] + rslt ;
	rslt++;
    if( ntupleRawTree.b5_mm_kal_vtx_prob[mumuIdx] <  muMuVtxProbabilityMin ) 
		 return cutFlowOffsets["diMuonVertexSelection"] + rslt ;
	rslt++;
   // std::cout<<ntupleRawTree.b5_mm_docatrk[mumuIdx]<<"  < "<<svMinDocaTrack<<"\n";
    if( ntupleRawTree.b5_mm_docatrk[mumuIdx]  < svMinDocaTrack ) 
		 return cutFlowOffsets["diMuonVertexSelection"] + rslt ;
	rslt++;
    if( ntupleRawTree.b5_mm_closetrk[mumuIdx] > svMaxNTracksClose ) 
		 return cutFlowOffsets["diMuonVertexSelection"] + rslt ;
	rslt++;
    if( ntupleRawTree.b5_mm_closetrks1[mumuIdx] > svMaxN1SigmaTracksClose ) 
		 return cutFlowOffsets["diMuonVertexSelection"] + rslt ;
	rslt++;
    if( ntupleRawTree.b5_mm_closetrks2[mumuIdx] > svMaxN2SigmaTracksClose ) 
		 return cutFlowOffsets["diMuonVertexSelection"] + rslt ;
	rslt++;
    if( ntupleRawTree.b5_mm_closetrks3[mumuIdx] > svMaxN3SigmaTracksClose ) 
		 return cutFlowOffsets["diMuonVertexSelection"] + rslt ;

 //   if( cos(ntupleRawTree.b5_mm_kin_alpha[mumuIdx])  > 0.80 ) return 1;
 //   if( ntupleRawTree.b5_mm_kin_vtx_chi2dof[mumuIdx] > 2.20  ) return 2;
 //   if( ntupleRawTree.b5_mm_kin_sl3d[mumuIdx]  < 13.0  ) return 3;
 //   if( ntupleRawTree.b5_mm_kin_slxy[mumuIdx]  < 3.0   ) return 4;
 //   if( ntupleRawTree.b5_mm_kin_pvip[mumuIdx]  > 0.008 ) return 5;
 //   if( ntupleRawTree.b5_mm_kin_spvip[mumuIdx] < 2.0   ) return 6;
 //   if( ntupleRawTree.b5_mm_doca[mumuIdx]  > 0.100 ) return 7;
 //   if( ntupleRawTree.b5_mm_iso[mumuIdx]   < 0.80  ) return 8;
 //   if( ntupleRawTree.b5_mm_closetrk[mumuIdx]  >= 2.0 ) return 9;
    return 0;
}

Int_t BMMGAnalysis::doPhotonSelection(Int_t scIdx)
{
    Int_t rslt=0;
    rslt++;
    if( ntupleRawTree.bG_scEt[scIdx] < minSCEt )  return  cutFlowOffsets["photonSelection"]+rslt;
    rslt++;
    if( abs(ntupleRawTree.bG_scEta[scIdx]) > maxSCEta )  return  cutFlowOffsets["photonSelection"]+rslt;
    
    rslt++;
    if( storageArrayDouble[ scIdx + candidateMapInt["scPhotonMVAScore"]  ] < photonIDcut  ) 
    return cutFlowOffsets["photonSelection"]+rslt;
    
    return 0;
}

Int_t BMMGAnalysis::doBMMGSelection(Int_t mumuIdx, Int_t phoSCIdx)
{
   Int_t rslt(0);
   rslt++;
   if( ntupleRawTree.bG_scEt[phoSCIdx]  < minSCEt )  
        return cutFlowOffsets["bmmgSelection"] + rslt;
   rslt++;
   if( ntupleRawTree.bG_scEta[phoSCIdx]  > maxSCEta)  
        return cutFlowOffsets["bmmgSelection"] + rslt;
   
   auto dr=getDR(ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx],
                 ntupleRawTree.b5_mm_kin_eta[mumuIdx] , ntupleRawTree.b5_mm_kin_phi[mumuIdx]  );
   rslt++;
   if(dr > maxDimuPhotonDr ) 
        return cutFlowOffsets["bmmgSelection"] + rslt;
   rslt++;
   if(getDCAGammaToDimuVertex(mumuIdx,phoSCIdx) > svGammaMaxDCA) 
        return cutFlowOffsets["bmmgSelection"] + rslt;

   return 0;
}

Int_t BMMGAnalysis::getMuonMatch(Double_t muEta,Double_t muPhi)
{

    Double_t dr,drMin( drMaxForGenMatchMu );
    Int_t muMatch=-1;
    
    std::cout<<"nMM = "<<ntupleRawTree.b5_nmm<<" , nmu : "<< ntupleRawTree.b5_nMuon<<"\n";
    for(UInt_t idx =0; idx < ntupleRawTree.b5_nMuon ; idx++)
    {
        dr=getDR(muEta,muPhi,ntupleRawTree.b5_Muon_eta[idx],ntupleRawTree.b5_Muon_phi[idx]);
        std::cout<<"dr for "<<ntupleRawTree.b5_Muon_eta[idx]<<" , "<<ntupleRawTree.b5_Muon_phi[idx]<<" = "<<dr<<"\n";
        if(dr < drMin )
        {
            drMin = dr;
            muMatch=idx;
        }
        if( dr < 1e-3) break;
    }
    std::cout<<"\r\rdrMin for "<<muEta<<","<<muPhi<<" : "<<drMin<<"\n";
    
    return muMatch;

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

Double_t BMMGAnalysis::getDCAGammaToDimuVertex(Int_t mumuIdx,Int_t phoId)
{
   
    Double_t x1[3],x2[3],p[3];
    
    x1[0] = ntupleRawTree.b5_mm_kin_vtx_x[mumuIdx];
    x1[1] = ntupleRawTree.b5_mm_kin_vtx_y[mumuIdx];
    x1[2] = ntupleRawTree.b5_mm_kin_vtx_z[mumuIdx];
    
    x2[0] = 0.0; x2[1] = 0.0 ; x2[2] =0.0;

    p[0] = ntupleRawTree.b5_SV_x[mumuIdx];
    p[1] = ntupleRawTree.b5_SV_y[mumuIdx];
    p[2] = ntupleRawTree.b5_SV_z[mumuIdx];
    
    return  getDCALineAndPoint( x1 , x2 , p );
}

void BMMGAnalysis::SkimData()
{
 
    /************************************************************************************

 Make sure the branches used here are not turned off to 0 by BMMGAnalysis::setupBranchStatus()

    *************************************************************************************/

    Double_t dr;
    
    std::cout<<"\nBegining Skimming Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t EventCountWithDimuCand=0;
    Long64_t EventCountWithDimuVertexCand=0;
    Long64_t nb = 0,nbytes=0 ;
    Int_t eventLostAt(0);
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    auto nDiMuNoVertexCandidates=0;
    Int_t prevRun(-1),prevLumi(-1);
    bool goodRunLumi = false;
    bool skimAdd=true;
    th1fStore["ProcessingSummary"]->Fill("Skimmer",1);
    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {   
    
       EventCount++;
       skimAdd=false;
       eventLostAt=0; 
       nDiMuCandidates=0;
       isTriggerd=false;
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);

       if (ientry_evt < 0) break;

       th1fStore["ProcessingSummary"]->Fill("TotalEvents",1);
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
        
       if(doRunLumiLog) { runLumiLogger.addRunLumi(ntupleRawTree.b5_run,ntupleRawTree.b5_luminosityBlock); }


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

       th1fStore["ProcessingSummary"]->Fill("GoodLumiEvents",1);
       // std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "<<"\n";
       // Checks if atleast 1 PV is there .. by default there will always be  one pV , the beamspot
       if(ntupleRawTree.bG_nPrimaryVertex < 1 ) continue;
       
       // Trigger Selection
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Bs ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_DoubleMu4_3_Jpsi ) isTriggerd=true;
       if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing ) isTriggerd=true;
       //if( ntupleRawTree.b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R ) isTriggerd=true;
    
       if( (not isTriggerd) and doTriggerFiltering ) continue;

       th1fStore["ProcessingSummary"]->Fill("TriggerPassEvents",1);

       for(Int_t i=0;i<NSTORAGE_ARRAY_MAX;i++)
            storageArrayDouble[i]=0;
       for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
            photonSelectionCheck[phoSCIdx]=-1;
       if(doPhotonMVA)
       {
            doPhotonMVAScores();
       }
       if(doDimuonMVA)
       {
            doDimuonMVAScores();
       }
                
       int rslt=0;
       nBMMGCandidates=0;

       fill_muonHists();
       fill_scHists();
       fill_photonHists();
       nDiMuNoVertexCandidates=0;
       eventLostAt=0;
       rslt=0;
       
       if(ntupleRawTree.b5_nmm <1 ) 
       { 
            if(eventLostAt < ( 1 + cutFlowOffsets["basicCuts"]) )   eventLostAt=cutFlowOffsets["basicCuts"]+1;
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
           
           fill_dimuonHists(mumuIdx);
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu1_index[mumuIdx], true);
           if(eventLostAt < rslt) eventLostAt=rslt;
           if(rslt > 0) continue;
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu2_index[mumuIdx], false);
           if(eventLostAt < rslt) eventLostAt=rslt;
           if(rslt > 0) continue;
           
           th1fStore["ProcessingSummary"]->Fill("EventCountWith2GoodMuons",1);

           diMuLV.SetPtEtaPhiM(     ntupleRawTree.b5_mm_kin_pt[mumuIdx],   \
                                    ntupleRawTree.b5_mm_kin_eta[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_phi[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_mass[mumuIdx]  );
           nBMMGCandidatesPerDimu=0;
           
           if(ntupleRawTree.bG_nSC<1)
           {
                continue;
           }

           for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
           {
               // photon selection
               if(photonSelectionCheck[phoSCIdx] < 0)
               {
                    photonSelectionCheck[phoSCIdx]=doPhotonSelection(phoSCIdx );
               }
               
               auto et= ntupleRawTree.bG_scE[phoSCIdx]/cosh(ntupleRawTree.bG_scEta[phoSCIdx]);
               photonLV.SetPtEtaPhiM( et,ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx], PHO_MASS );
               dr=photonLV.DeltaR(diMuLV);
               if( dr > 1.3 ) continue;
               skimAdd=true;
           }
           
           
           storageArrayInt[ nDiMuCandidates + candidateMapInt["nBMMGCandidatesPerDimu"]  ]   = nBMMGCandidatesPerDimu ;
           nDiMuCandidates++;

           if(nDiMuCandidates >= NDIMU_MAX)
           {
                std::cout<<" Dimuon count per event above NDIMU_MAX , Aborting !! \n";
                exit(9);
           }
            
       }
       
    
       if(nDiMuNoVertexCandidates > 0)
       {
            th1fStore["ProcessingSummary"]->Fill("EventCountWithDimuCand",1);
            EventCountWithDimuCand++;
       }
       if(nDiMuCandidates > 0)
       {
            EventCountWithDimuVertexCand++;
            th1fStore["ProcessingSummary"]->Fill("EventCountWithDimuVertexCand",1);
       }
       if(nBMMGCandidates >0)
       {
        // std::cout<<"\t "<<ntupleRawTree.b5_run<<" , "<<ntupleRawTree.b5_event<<" : Number of nBMMGCandidates = "<<nBMMGCandidates<<"\n";
        EventCountWithCand++;
        th1fStore["ProcessingSummary"]->Fill("EventCountWithBMMGCand",1);
       }
       
       fill_globalEventHists();
       FillCutFlow(eventLostAt);
       
       if(skimAdd)
       { 
         EventCountWithCand++;
         if(outTree)  outTree->Fill();
       }
    }
    std::cout<<"\n\n"
            <<"  Total Number of Events           = "<<EventCount<<"\n"
            <<"  Number of Events with candidates = "<<EventCountWithCand<<"\n"
            <<"  Analysis loop completed"
            <<"  \n\n";

}
