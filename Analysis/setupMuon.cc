
	auto idx = offset ;
 
	candidateMap["b5_MuonId_chi2LocalPosition"] = idx; 
	outTree->Branch("b5_MuonId_chi2LocalPosition" , &(storageArray[idx] ) , "b5_MuonId_chi2LocalPosition[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_glbNormChi2"] = idx; 
	outTree->Branch("b5_MuonId_glbNormChi2" , &(storageArray[idx] ) , "b5_MuonId_glbNormChi2[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_glbTrackProbability"] = idx; 
	outTree->Branch("b5_MuonId_glbTrackProbability" , &(storageArray[idx] ) , "b5_MuonId_glbTrackProbability[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match1_dX"] = idx; 
	outTree->Branch("b5_MuonId_match1_dX" , &(storageArray[idx] ) , "b5_MuonId_match1_dX[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match1_dY"] = idx; 
	outTree->Branch("b5_MuonId_match1_dY" , &(storageArray[idx] ) , "b5_MuonId_match1_dY[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match1_pullDxDz"] = idx; 
	outTree->Branch("b5_MuonId_match1_pullDxDz" , &(storageArray[idx] ) , "b5_MuonId_match1_pullDxDz[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match1_pullDyDz"] = idx; 
	outTree->Branch("b5_MuonId_match1_pullDyDz" , &(storageArray[idx] ) , "b5_MuonId_match1_pullDyDz[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match1_pullX"] = idx; 
	outTree->Branch("b5_MuonId_match1_pullX" , &(storageArray[idx] ) , "b5_MuonId_match1_pullX[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match1_pullY"] = idx; 
	outTree->Branch("b5_MuonId_match1_pullY" , &(storageArray[idx] ) , "b5_MuonId_match1_pullY[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match2_dX"] = idx; 
	outTree->Branch("b5_MuonId_match2_dX" , &(storageArray[idx] ) , "b5_MuonId_match2_dX[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match2_dY"] = idx; 
	outTree->Branch("b5_MuonId_match2_dY" , &(storageArray[idx] ) , "b5_MuonId_match2_dY[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match2_pullDxDz"] = idx; 
	outTree->Branch("b5_MuonId_match2_pullDxDz" , &(storageArray[idx] ) , "b5_MuonId_match2_pullDxDz[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match2_pullDyDz"] = idx; 
	outTree->Branch("b5_MuonId_match2_pullDyDz" , &(storageArray[idx] ) , "b5_MuonId_match2_pullDyDz[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match2_pullX"] = idx; 
	outTree->Branch("b5_MuonId_match2_pullX" , &(storageArray[idx] ) , "b5_MuonId_match2_pullX[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_match2_pullY"] = idx; 
	outTree->Branch("b5_MuonId_match2_pullY" , &(storageArray[idx] ) , "b5_MuonId_match2_pullY[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_newSoftMuonMva"] = idx; 
	outTree->Branch("b5_MuonId_newSoftMuonMva" , &(storageArray[idx] ) , "b5_MuonId_newSoftMuonMva[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_trkKink"] = idx; 
	outTree->Branch("b5_MuonId_trkKink" , &(storageArray[idx] ) , "b5_MuonId_trkKink[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_trkValidFrac"] = idx; 
	outTree->Branch("b5_MuonId_trkValidFrac" , &(storageArray[idx] ) , "b5_MuonId_trkValidFrac[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_highPurity"] = idx; 
	outTree->Branch("b5_MuonId_highPurity" , &(storageArray[idx] ) , "b5_MuonId_highPurity[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_nLostHitsInner"] = idx; 
	outTree->Branch("b5_MuonId_nLostHitsInner" , &(storageArray[idx] ) , "b5_MuonId_nLostHitsInner[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_nLostHitsOn"] = idx; 
	outTree->Branch("b5_MuonId_nLostHitsOn" , &(storageArray[idx] ) , "b5_MuonId_nLostHitsOn[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_nLostHitsOuter"] = idx; 
	outTree->Branch("b5_MuonId_nLostHitsOuter" , &(storageArray[idx] ) , "b5_MuonId_nLostHitsOuter[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_nPixels"] = idx; 
	outTree->Branch("b5_MuonId_nPixels" , &(storageArray[idx] ) , "b5_MuonId_nPixels[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_nValidHits"] = idx; 
	outTree->Branch("b5_MuonId_nValidHits" , &(storageArray[idx] ) , "b5_MuonId_nValidHits[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_trkLayers"] = idx; 
	outTree->Branch("b5_MuonId_trkLayers" , &(storageArray[idx] ) , "b5_MuonId_trkLayers[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_trkLostLayersInner"] = idx; 
	outTree->Branch("b5_MuonId_trkLostLayersInner" , &(storageArray[idx] ) , "b5_MuonId_trkLostLayersInner[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_trkLostLayersOn"] = idx; 
	outTree->Branch("b5_MuonId_trkLostLayersOn" , &(storageArray[idx] ) , "b5_MuonId_trkLostLayersOn[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_MuonId_trkLostLayersOuter"] = idx; 
	outTree->Branch("b5_MuonId_trkLostLayersOuter" , &(storageArray[idx] ) , "b5_MuonId_trkLostLayersOuter[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_dxy"] = idx; 
	outTree->Branch("b5_Muon_dxy" , &(storageArray[idx] ) , "b5_Muon_dxy[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_dxyErr"] = idx; 
	outTree->Branch("b5_Muon_dxyErr" , &(storageArray[idx] ) , "b5_Muon_dxyErr[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_dxybs"] = idx; 
	outTree->Branch("b5_Muon_dxybs" , &(storageArray[idx] ) , "b5_Muon_dxybs[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_dz"] = idx; 
	outTree->Branch("b5_Muon_dz" , &(storageArray[idx] ) , "b5_Muon_dz[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_dzErr"] = idx; 
	outTree->Branch("b5_Muon_dzErr" , &(storageArray[idx] ) , "b5_Muon_dzErr[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_eta"] = idx; 
	outTree->Branch("b5_Muon_eta" , &(storageArray[idx] ) , "b5_Muon_eta[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_ip3d"] = idx; 
	outTree->Branch("b5_Muon_ip3d" , &(storageArray[idx] ) , "b5_Muon_ip3d[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_jetPtRelv2"] = idx; 
	outTree->Branch("b5_Muon_jetPtRelv2" , &(storageArray[idx] ) , "b5_Muon_jetPtRelv2[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_jetRelIso"] = idx; 
	outTree->Branch("b5_Muon_jetRelIso" , &(storageArray[idx] ) , "b5_Muon_jetRelIso[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_mass"] = idx; 
	outTree->Branch("b5_Muon_mass" , &(storageArray[idx] ) , "b5_Muon_mass[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_miniPFRelIso_all"] = idx; 
	outTree->Branch("b5_Muon_miniPFRelIso_all" , &(storageArray[idx] ) , "b5_Muon_miniPFRelIso_all[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_miniPFRelIso_chg"] = idx; 
	outTree->Branch("b5_Muon_miniPFRelIso_chg" , &(storageArray[idx] ) , "b5_Muon_miniPFRelIso_chg[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_pfRelIso03_all"] = idx; 
	outTree->Branch("b5_Muon_pfRelIso03_all" , &(storageArray[idx] ) , "b5_Muon_pfRelIso03_all[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_pfRelIso03_chg"] = idx; 
	outTree->Branch("b5_Muon_pfRelIso03_chg" , &(storageArray[idx] ) , "b5_Muon_pfRelIso03_chg[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_pfRelIso04_all"] = idx; 
	outTree->Branch("b5_Muon_pfRelIso04_all" , &(storageArray[idx] ) , "b5_Muon_pfRelIso04_all[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_phi"] = idx; 
	outTree->Branch("b5_Muon_phi" , &(storageArray[idx] ) , "b5_Muon_phi[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_pt"] = idx; 
	outTree->Branch("b5_Muon_pt" , &(storageArray[idx] ) , "b5_Muon_pt[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_ptErr"] = idx; 
	outTree->Branch("b5_Muon_ptErr" , &(storageArray[idx] ) , "b5_Muon_ptErr[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_segmentComp"] = idx; 
	outTree->Branch("b5_Muon_segmentComp" , &(storageArray[idx] ) , "b5_Muon_segmentComp[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_sip3d"] = idx; 
	outTree->Branch("b5_Muon_sip3d" , &(storageArray[idx] ) , "b5_Muon_sip3d[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_softMva"] = idx; 
	outTree->Branch("b5_Muon_softMva" , &(storageArray[idx] ) , "b5_Muon_softMva[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_tkRelIso"] = idx; 
	outTree->Branch("b5_Muon_tkRelIso" , &(storageArray[idx] ) , "b5_Muon_tkRelIso[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_tunepRelPt"] = idx; 
	outTree->Branch("b5_Muon_tunepRelPt" , &(storageArray[idx] ) , "b5_Muon_tunepRelPt[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_mvaLowPt"] = idx; 
	outTree->Branch("b5_Muon_mvaLowPt" , &(storageArray[idx] ) , "b5_Muon_mvaLowPt[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_mvaTTH"] = idx; 
	outTree->Branch("b5_Muon_mvaTTH" , &(storageArray[idx] ) , "b5_Muon_mvaTTH[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_charge"] = idx; 
	outTree->Branch("b5_Muon_charge" , &(storageArray[idx] ) , "b5_Muon_charge[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_jetIdx"] = idx; 
	outTree->Branch("b5_Muon_jetIdx" , &(storageArray[idx] ) , "b5_Muon_jetIdx[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_nStations"] = idx; 
	outTree->Branch("b5_Muon_nStations" , &(storageArray[idx] ) , "b5_Muon_nStations[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_nTrackerLayers"] = idx; 
	outTree->Branch("b5_Muon_nTrackerLayers" , &(storageArray[idx] ) , "b5_Muon_nTrackerLayers[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_pdgId"] = idx; 
	outTree->Branch("b5_Muon_pdgId" , &(storageArray[idx] ) , "b5_Muon_pdgId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_tightCharge"] = idx; 
	outTree->Branch("b5_Muon_tightCharge" , &(storageArray[idx] ) , "b5_Muon_tightCharge[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_fsrPhotonIdx"] = idx; 
	outTree->Branch("b5_Muon_fsrPhotonIdx" , &(storageArray[idx] ) , "b5_Muon_fsrPhotonIdx[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_highPtId"] = idx; 
	outTree->Branch("b5_Muon_highPtId" , &(storageArray[idx] ) , "b5_Muon_highPtId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_highPurity"] = idx; 
	outTree->Branch("b5_Muon_highPurity" , &(storageArray[idx] ) , "b5_Muon_highPurity[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_inTimeMuon"] = idx; 
	outTree->Branch("b5_Muon_inTimeMuon" , &(storageArray[idx] ) , "b5_Muon_inTimeMuon[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_isGlobal"] = idx; 
	outTree->Branch("b5_Muon_isGlobal" , &(storageArray[idx] ) , "b5_Muon_isGlobal[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_isPFcand"] = idx; 
	outTree->Branch("b5_Muon_isPFcand" , &(storageArray[idx] ) , "b5_Muon_isPFcand[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_isTracker"] = idx; 
	outTree->Branch("b5_Muon_isTracker" , &(storageArray[idx] ) , "b5_Muon_isTracker[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_jetNDauCharged"] = idx; 
	outTree->Branch("b5_Muon_jetNDauCharged" , &(storageArray[idx] ) , "b5_Muon_jetNDauCharged[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_looseId"] = idx; 
	outTree->Branch("b5_Muon_looseId" , &(storageArray[idx] ) , "b5_Muon_looseId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_mediumId"] = idx; 
	outTree->Branch("b5_Muon_mediumId" , &(storageArray[idx] ) , "b5_Muon_mediumId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_mediumPromptId"] = idx; 
	outTree->Branch("b5_Muon_mediumPromptId" , &(storageArray[idx] ) , "b5_Muon_mediumPromptId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_miniIsoId"] = idx; 
	outTree->Branch("b5_Muon_miniIsoId" , &(storageArray[idx] ) , "b5_Muon_miniIsoId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_multiIsoId"] = idx; 
	outTree->Branch("b5_Muon_multiIsoId" , &(storageArray[idx] ) , "b5_Muon_multiIsoId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_mvaId"] = idx; 
	outTree->Branch("b5_Muon_mvaId" , &(storageArray[idx] ) , "b5_Muon_mvaId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_mvaLowPtId"] = idx; 
	outTree->Branch("b5_Muon_mvaLowPtId" , &(storageArray[idx] ) , "b5_Muon_mvaLowPtId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_pfIsoId"] = idx; 
	outTree->Branch("b5_Muon_pfIsoId" , &(storageArray[idx] ) , "b5_Muon_pfIsoId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_softId"] = idx; 
	outTree->Branch("b5_Muon_softId" , &(storageArray[idx] ) , "b5_Muon_softId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_softMvaId"] = idx; 
	outTree->Branch("b5_Muon_softMvaId" , &(storageArray[idx] ) , "b5_Muon_softMvaId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_tightId"] = idx; 
	outTree->Branch("b5_Muon_tightId" , &(storageArray[idx] ) , "b5_Muon_tightId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_tkIsoId"] = idx; 
	outTree->Branch("b5_Muon_tkIsoId" , &(storageArray[idx] ) , "b5_Muon_tkIsoId[nMuons]/D"); 
	idx+=size; 

	candidateMap["b5_Muon_triggerIdLoose"] = idx; 
	outTree->Branch("b5_Muon_triggerIdLoose" , &(storageArray[idx] ) , "b5_Muon_triggerIdLoose[nMuons]/D"); 
	idx+=size; 

	return idx;