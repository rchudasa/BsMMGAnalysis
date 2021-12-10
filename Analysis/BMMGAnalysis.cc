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
       fill_dimuonHists();
       nDiMuNoVertexCandidates=0;
       for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    

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



           fill_dimuonHists(mumuIdx);
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
     th1fStore["dimu_L3Dsignificance"]  = new TH1F("dimu_L3Dsignificance","l_{3D}(#mu#mu)/#sigma(l_{3D}) for #mu#mu Candidate" ,50 , 0.0 , 2.5  );
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
     th1fStore["bmmg_dimuGammaDeltaR"] = new TH1F("bmmg_mumuDeltaR","#Delta(#mu,#mu) for #mu#mu#gamma", deltaRNBins , deltaRMin  , deltaRMax  );
     
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
}



void BMMGAnalysis::fill_muonHists()
{   
    for(int i=0;i<ntupleRawTree.b5_nMuon ; i++)
    {
        if( doMuonSelection(i,false) ) continue ;
        
        th1fStore["mu_Pt"]     ->Fill(ntupleRawTree.b5_Muon_pt[i]);
        th1fStore["mu_Eta"]    ->Fill(ntupleRawTree.b5_Muon_eta[i]);
        th1fStore["mu_Charge"] ->Fill(ntupleRawTree.b5_Muon_charge[i]);
        th1fStore["mu_dz"]     ->Fill(ntupleRawTree.b5_Muon_dz[i]);
    }
    
}

void BMMGAnalysis::fill_scHists()
{
    for(int i=0;i<ntupleRawTree.bG_nSC ; i++)
    {
        th1fStore["sc_Pt"]    ->Fill(ntupleRawTree.bG_scEt[i]);
        th1fStore["sc_Eta"]   ->Fill(ntupleRawTree.bG_scEta[i]); 
        th1fStore["sc_mva"]   ->Fill(storageArrayDouble[ i + candidateMapInt["scPhotonMVAScore"]  ]); 
    }
}

void BMMGAnalysis::fill_photonHists()
{
     for(int i=0;i<ntupleRawTree.bG_nSC ; i++)
        {
            if(doPhotonSelection(i)) continue;

            th1fStore["photon_Pt"]    ->Fill(ntupleRawTree.bG_scEt[i]);
            th1fStore["photon_Eta"]   ->Fill(ntupleRawTree.bG_scEta[i]); 
            th1fStore["photon_mva"]   ->Fill(storageArrayDouble[ i + candidateMapInt["scPhotonMVAScore"]  ]); 
        }

}
void BMMGAnalysis::fill_dimuonHists()
{
    th1fStore["dimu_nDimuons"]->Fill(ntupleRawTree.b5_nmm);
  for(Int_t mumuIdx=0;mumuIdx<ntupleRawTree.b5_nmm ; mumuIdx++)
   {
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
void BMMGAnalysis::fill_dimuonHists(Int_t mumuIdx)
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
    
    
    if(  ntupleRawTree.b5_Muon_pt[muIdx] < 4.0 )  return 1;
   
    // auto muGblId= ntupleRawTree.b5_mm_mu1_index[muIdx];
    // Soft Muon ID   : Cut Based
    
    /*

    if( not ntupleRawTree.Muon_softId[muGblId] ) return 1;
    
    */

    // Soft Muon MVA 
    
    //if( not ntupleRawTree.b5_Muon_softMvaId[muIdx] ) return 2;

    // BMM5 MVA


    if(ntupleRawTree.b5_MuonId_newSoftMuonMva[muIdx] < BDTWorkingPoint ) return 3;

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


