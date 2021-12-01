
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
    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
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
       for(int mumuIdx=0; mumuIdx < ntupleRawTree.b5_nmm;mumuIdx++)
       {    
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

               
               if(bmmgLV.M() < minMMGMass or bmmgLV.M() > maxMMGMass ) continue;

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
           // assignDimuonBMMGCandidates(candidateMapDouble,ntupleRawTree,mumuIdx,storageArrayDouble,nDiMuCandidates);
           
           nDiMuCandidates++;

           if(nDiMuCandidates >= NDIMU_MAX)
           {
                std::cout<<" Dimuon count per event above NDIMU_MAX , Aborting !! \n";
                exit(9);
           }
            
       }
       
       //assignSC(candidateMapDouble, ntupleRawTree,storageArrayDouble);
       //assignMuons(candidateMapDouble, ntupleRawTree,storageArrayDouble);
    
       EventCount++;
       if(nBMMGCandidates >0)
       {
        std::cout<<"\t "<<ntupleRawTree.b5_run<<" , "<<ntupleRawTree.b5_event<<" : Number of nBMMGCandidates = "<<nBMMGCandidates<<"\n";
        EventCountWithCand++;
       }
       outTree->Fill();
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
    
    if( storageArrayDouble[ scIdx + candidateMapInt["scPhotonMVAScore"]  ] > 0.43 ) return 1;
    
    return 0;
}
