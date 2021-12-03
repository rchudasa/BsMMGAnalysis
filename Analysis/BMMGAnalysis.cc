
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

     const int   massNBins=80;
     const float massMin(0.0);
     const float massMax(20.0);
     
     const int   deltaRNBins=30;
     const float deltaRMin(0.0);
     const float deltaRMax(3.0);
    
    // Mu Hist
     th1fStore["mu_Pt"]     = new TH1F("mu_Pt","Pt of mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["mu_Eta"]     = new TH1F("mu_Eta","Eta of mu Candidate", etaBins , etamin  , etamax  );
     th1fStore["mu_Charge"]     = new TH1F("mu_charge","charge of mu Candidate", 3 , -1.5  , 1.5  );
     th1fStore["mu_dz"]     = new TH1F("mu_dz","dz of mu Candidate", 60*4 , -30.0  , 30.0  );
    

    // SC Hists
     th1fStore["sc_Pt"]     = new TH1F("sc_Pt","Pt of sc Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["sc_Eta"]     = new TH1F("sc_Eta","Eta of sc Candidate", etaBins , etamin  , etamax  );
     th1fStore["sc_mva"]     = new TH1F("sc_mva","MVA score of sc Candidate", 80 , -2.0  , 2.0  );

    // Photon Hists
     th1fStore["photon_Pt"]     = new TH1F("photon_Pt","Pt of photon Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["photon_Eta"]     = new TH1F("photon_Eta","Eta of photon Candidate", etaBins , etamin  , etamax  );
     th1fStore["photon_mva"]     = new TH1F("photon_mva","MVA photonore of photon Candidate", 80 , -2.0  , 2.0  );

    // Mu Mu Hist
     th1fStore["dimu_mumuDeltaR"]   = new TH1F("mumuDeltaR","#Delta(#mu,#mu)", deltaRNBins , deltaRMin  , deltaRMax  );
     th1fStore["dimu_muLeadPt"]     = new TH1F("dimu_muLeadPt","Pt of lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["dimu_muSubLeadPt"]  = new TH1F("dimu_muSubLeadPt","Pt of sub-lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["dimu_mass"]         = new TH1F("dimu_mass","InvMass of mumu Candidate", massNBins , massMin , massMax  );
     th1fStore["dimu_dxy"]          = new TH1F("dimu_dxy","dxy of dimu  Candidate", 20 ,0.0,2.50  );
     th1fStore["dimu_vtxProba"]     = new TH1F("dimu_vtxProbay","vtx Proba. of dimu  Candidate", 20 ,0.0,2.50  );
     th1fStore["dimu_alpha"]        = new TH1F("dimu_alpha","alpha of dimu  Candidate", 20 ,0.0,3.50  );

    // Dimuon Environment Hists
    th1fStore["dimu_photonMultiplicityDrMax0p5"] = new TH1F("dimu_photonMultiplicityDrMax0p5","Multiplicity of photon around #mu#mu , #Delta R < 0.5",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax0p7"] = new TH1F("dimu_photonMultiplicityDrMax0p7","Multiplicity of photon around #mu#mu , #Delta R < 0.7",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax0p9"] = new TH1F("dimu_photonMultiplicityDrMax0p9","Multiplicity of photon around #mu#mu , #Delta R < 0.9",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax1p1"] = new TH1F("dimu_photonMultiplicityDrMax1p1","Multiplicity of photon around #mu#mu , #Delta R < 1.1",15, -0.5,14.5);
    th1fStore["dimu_photonMultiplicityDrMax1p3"] = new TH1F("dimu_photonMultiplicityDrMax1p3","Multiplicity of photon around #mu#mu , #Delta R < 1.3",15, -0.5,14.5);
 
    // BMMG HISTS

     th1fStore["bmmg_muLeadPt"]         = new TH1F("bmmg_muLeadPt","Pt of lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["bmmg_muSubLeadPt"]      = new TH1F("bmmg_muSubLeadPt","Pt of sub-lead mu Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["bmmg_photonPt"]         = new TH1F("bmmg_photonPt","Pt of photon  Candidate", pTBins , pTmin  , pTmax  );
     th1fStore["bmmg_mumuDeltaR"]       = new TH1F("bmmg_mumuDeltaR","#Delta(#mu,#mu)", deltaRNBins , deltaRMin  , deltaRMax  );
     th1fStore["bmmg_dimuGammaDeltaR"]  = new TH1F("bmmg_dimuGammaDeltaR","#Delta(#mu#mu,#gamma)", deltaRNBins , deltaRMin  , deltaRMax  );
     th1fStore["bmmg_mmMass"]           = new TH1F("bmmg_mmMass","InvMass of mumu Candidate", massNBins , massMin , massMax  );
     th1fStore["bmmg_mmgMass"]          = new TH1F("bmmg_mmgMass","InvMass of mmgCandidate" , massNBins , massMin , massMax  );

    // Global Event Hists
    th1fStore["dimu_bmmgCandidateMultiplicity"]   = new TH1F("dimu_bmmgCandidateMultiplicity","Multiplicity of BMMG Candidates Per dimuon",15, -0.5,14.5);
    th1fStore["dimu_Multiplicity"]       = new TH1F("dimu_Multiplicity","Multiplicity of dimu Candidates in an event",15, -0.5,14.5);
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

void BMMGAnalysis::fill_dimuonHists(Int_t mumuIdx)
{
    
     auto dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx],
                     ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
     th1fStore["dimu_mumuDeltaR"]  ->Fill(dr);
     th1fStore["dimu_muLeadPt"]    ->Fill(ntupleRawTree.b5_mm_kin_mu1pt[mumuIdx]);
     th1fStore["dimu_muSubLeadPt"] ->Fill(ntupleRawTree.b5_mm_kin_mu2pt[mumuIdx]);
     th1fStore["dimu_mass"]        ->Fill(ntupleRawTree.b5_mm_kin_mass[mumuIdx]);
     th1fStore["dimu_dxy"]         ->Fill(ntupleRawTree.b5_mm_kin_lxy[mumuIdx]);
     th1fStore["dimu_vtxProba"]    ->Fill(ntupleRawTree.b5_mm_kin_vtx_prob[mumuIdx]);
     th1fStore["dimu_alpha"]       ->Fill(ntupleRawTree.b5_mm_kin_alpha[mumuIdx]);

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

    th1fStore["bmmg_muLeadPt"]         ->Fill(ntupleRawTree.b5_mm_mu1_pt[mumuIdx]); 
    th1fStore["bmmg_muSubLeadPt"]      ->Fill(ntupleRawTree.b5_mm_mu1_pt[mumuIdx]);
    th1fStore["bmmg_photonPt"]         ->Fill(ntupleRawTree.bG_scEt[phoSCIdx]);
    th1fStore["bmmg_mmMass"]           ->Fill(ntupleRawTree.b5_mm_kin_mass[mumuIdx]);
    th1fStore["bmmg_mmgMass"]          ->Fill(bmmgLV.M());

}

void BMMGAnalysis::fill_globalEventHists()
{
    th1fStore["bmmg_Multiplicity"]->Fill(nBMMGCandidates);
    th1fStore["dimu_Multiplicity"]->Fill(nDiMuCandidates); 

    for(int i=0;i<nDiMuCandidates;i++)
        th1fStore["dimu_bmmgCandidateMultiplicity"]->Fill(i+storageArrayInt[candidateMapInt["nBMMGCandidatesPerDimu"]  ] );
}

void BMMGAnalysis::Analyze()
{
 
    /************************************************************************************

 Make sure the branches used here are not turned off to 0 by BMMGAnalysis::setupBranchStatus()

    *************************************************************************************/

    TLorentzVector diMuLV,photonLV,bmmgLV;
    Double_t dr;
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();

    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {   

       nDiMuCandidates=0;
       isTriggerd=false;
       
       if(jentry%500 == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       }
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);

       if (ientry_evt < 0) break;

       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
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
        

       // BMMG Selection
       int rslt=0;
       float dimuEta,dimuPhi;
       

       nBMMGCandidates=0;

       if(doPhotonMVA)
       {
            doPhotonMVAScores();
       }
        
       fill_muonHists();
       fill_scHists();
       fill_photonHists();

       for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    
           fill_dimuonHists(mumuIdx);

           // Muon Selection
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu1_index[mumuIdx], true);
           if(rslt > 0) continue;
		   rslt=doMuonSelection( ntupleRawTree.b5_mm_mu2_index[mumuIdx], false);
           if(rslt > 0) continue;
           
           // Dimuon Selection
           diMuLV.SetPtEtaPhiM(     ntupleRawTree.b5_mm_kin_pt[mumuIdx],   \
                                    ntupleRawTree.b5_mm_kin_eta[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_phi[mumuIdx],  \
                                    ntupleRawTree.b5_mm_kin_mass[mumuIdx]  );
           
           dr= getDR( ntupleRawTree.b5_mm_kin_mu1eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu1phi[mumuIdx] ,
                      ntupleRawTree.b5_mm_kin_mu2eta[mumuIdx], ntupleRawTree.b5_mm_kin_mu2phi[mumuIdx] );
           if (dr > maxMuMuDr ) continue;
           if( diMuLV.M() < minDimuMass or diMuLV.M() > maxDimuMass) continue;   
           storageArrayDouble[nDiMuCandidates + candidateMapDouble["mumu_dr"] ]   = dr ;
            

            /*
                    VERTEX SELECTION STUFF
            */

            
            fill_dimuonEnvironmentHists(mumuIdx);

           dimuEta=ntupleRawTree.b5_mm_kin_eta[mumuIdx]  ;
           dimuPhi=ntupleRawTree.b5_mm_kin_phi[mumuIdx]  ;
           nBMMGCandidatesPerDimu=0;

           for(int phoSCIdx=0;phoSCIdx < ntupleRawTree.bG_nSC ; phoSCIdx++)
           {
               // photon selection
               if(photonSelectionCheck[phoSCIdx] < 0)
               {
                    photonSelectionCheck[phoSCIdx]=doPhotonSelection(phoSCIdx );
               }

               if( photonSelectionCheck[phoSCIdx] > 0) continue;

               dr=getDR(ntupleRawTree.bG_scEta[phoSCIdx],ntupleRawTree.bG_scPhi[phoSCIdx],dimuEta,dimuPhi);
               auto et= ntupleRawTree.bG_scE[phoSCIdx]/cosh(ntupleRawTree.bG_scEta[phoSCIdx]);
               photonLV.SetPtEtaPhiM( et       ,     ntupleRawTree.bG_scEta[phoSCIdx],
                                                       ntupleRawTree.bG_scPhi[phoSCIdx], PHO_MASS );
               bmmgLV = diMuLV + photonLV;
               if(dr > maxDimuPhotonDr ) continue;

               
               //if(bmmgLV.M() < minMMGMass or bmmgLV.M() > maxMMGMass ) continue;

               
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
           // assignDimuonBMMGCandidates(candidateMapDouble,ntupleRawTree,mumuIdx,storageArrayDouble,nDiMuCandidates);
           nDiMuCandidates++;

           if(nDiMuCandidates >= NDIMU_MAX)
           {
                std::cout<<" Dimuon count per event above NDIMU_MAX , Aborting !! \n";
                exit(9);
           }
            
       }
       
    
       EventCount++;
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
            <<"  Number of Events with candidates = "<<EventCountWithCand<<"\n"
            <<"  Analysis loop completed"
            <<"  \n\n";

}

Int_t BMMGAnalysis::doMuonSelection(Int_t muIdx, bool isLead)
{
    
   // auto muGblId= ntupleRawTree.b5_mm_mu1_index[muIdx];
    // Soft Muon ID   : Cut Based
    
    /*

    if( not ntupleRawTree.Muon_softId[muGblId] ) return 1;
    
    */

    // Soft Muon MVA  : Cut Based
    
    if( not ntupleRawTree.b5_Muon_softMvaId[muIdx] ) return 2;
    
    return 0;

    // BMM5 MVA

    const Double_t BDTLooseWorkingPoint (0.0); 
    const Double_t BDTTightWorkingPoint (0.0);

    if(ntupleRawTree.b5_MuonId_newSoftMuonMva[muIdx] < BDTTightWorkingPoint ) return 3;

    return 0;

}

Int_t BMMGAnalysis::doPhotonSelection(Int_t scIdx)
{
    
    if( storageArrayDouble[ scIdx + candidateMapInt["scPhotonMVAScore"]  ] < 0.43 ) return 1;
    
    return 0;
}
