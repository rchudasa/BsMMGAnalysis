
	auto idx = offset ;
 
	candidateMap["bG_scE"] = idx; 
	outTree->Branch("bG_scE" , &(storageArray[idx] ) , "bG_scE[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scEt"] = idx; 
	outTree->Branch("bG_scEt" , &(storageArray[idx] ) , "bG_scEt[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scRawE"] = idx; 
	outTree->Branch("bG_scRawE" , &(storageArray[idx] ) , "bG_scRawE[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scEta"] = idx; 
	outTree->Branch("bG_scEta" , &(storageArray[idx] ) , "bG_scEta[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPhi"] = idx; 
	outTree->Branch("bG_scPhi" , &(storageArray[idx] ) , "bG_scPhi[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scX"] = idx; 
	outTree->Branch("bG_scX" , &(storageArray[idx] ) , "bG_scX[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scY"] = idx; 
	outTree->Branch("bG_scY" , &(storageArray[idx] ) , "bG_scY[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scZ"] = idx; 
	outTree->Branch("bG_scZ" , &(storageArray[idx] ) , "bG_scZ[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scEtaWidth"] = idx; 
	outTree->Branch("bG_scEtaWidth" , &(storageArray[idx] ) , "bG_scEtaWidth[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPhiWidth"] = idx; 
	outTree->Branch("bG_scPhiWidth" , &(storageArray[idx] ) , "bG_scPhiWidth[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scRawEt"] = idx; 
	outTree->Branch("bG_scRawEt" , &(storageArray[idx] ) , "bG_scRawEt[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scMinDrWithGsfElectornSC_"] = idx; 
	outTree->Branch("bG_scMinDrWithGsfElectornSC_" , &(storageArray[idx] ) , "bG_scMinDrWithGsfElectornSC_[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFoundGsfMatch_"] = idx; 
	outTree->Branch("bG_scFoundGsfMatch_" , &(storageArray[idx] ) , "bG_scFoundGsfMatch_[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE5x5"] = idx; 
	outTree->Branch("bG_scE5x5" , &(storageArray[idx] ) , "bG_scE5x5[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE2x2Ratio"] = idx; 
	outTree->Branch("bG_scE2x2Ratio" , &(storageArray[idx] ) , "bG_scE2x2Ratio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE3x3Ratio"] = idx; 
	outTree->Branch("bG_scE3x3Ratio" , &(storageArray[idx] ) , "bG_scE3x3Ratio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scEMaxRatio"] = idx; 
	outTree->Branch("bG_scEMaxRatio" , &(storageArray[idx] ) , "bG_scEMaxRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE2ndRatio"] = idx; 
	outTree->Branch("bG_scE2ndRatio" , &(storageArray[idx] ) , "bG_scE2ndRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scETopRatio"] = idx; 
	outTree->Branch("bG_scETopRatio" , &(storageArray[idx] ) , "bG_scETopRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scERightRatio"] = idx; 
	outTree->Branch("bG_scERightRatio" , &(storageArray[idx] ) , "bG_scERightRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scEBottomRatio"] = idx; 
	outTree->Branch("bG_scEBottomRatio" , &(storageArray[idx] ) , "bG_scEBottomRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scELeftRatio"] = idx; 
	outTree->Branch("bG_scELeftRatio" , &(storageArray[idx] ) , "bG_scELeftRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE2x5MaxRatio"] = idx; 
	outTree->Branch("bG_scE2x5MaxRatio" , &(storageArray[idx] ) , "bG_scE2x5MaxRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE2x5TopRatio"] = idx; 
	outTree->Branch("bG_scE2x5TopRatio" , &(storageArray[idx] ) , "bG_scE2x5TopRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE2x5RightRatio"] = idx; 
	outTree->Branch("bG_scE2x5RightRatio" , &(storageArray[idx] ) , "bG_scE2x5RightRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE2x5BottomRatio"] = idx; 
	outTree->Branch("bG_scE2x5BottomRatio" , &(storageArray[idx] ) , "bG_scE2x5BottomRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scE2x5LeftRatio"] = idx; 
	outTree->Branch("bG_scE2x5LeftRatio" , &(storageArray[idx] ) , "bG_scE2x5LeftRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scSwissCross"] = idx; 
	outTree->Branch("bG_scSwissCross" , &(storageArray[idx] ) , "bG_scSwissCross[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scR9"] = idx; 
	outTree->Branch("bG_scR9" , &(storageArray[idx] ) , "bG_scR9[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scSigmaIetaIeta"] = idx; 
	outTree->Branch("bG_scSigmaIetaIeta" , &(storageArray[idx] ) , "bG_scSigmaIetaIeta[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scSigmaIetaIphi"] = idx; 
	outTree->Branch("bG_scSigmaIetaIphi" , &(storageArray[idx] ) , "bG_scSigmaIetaIphi[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scSigmaIphiIphi"] = idx; 
	outTree->Branch("bG_scSigmaIphiIphi" , &(storageArray[idx] ) , "bG_scSigmaIphiIphi[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e5x5"] = idx; 
	outTree->Branch("bG_scFull5x5_e5x5" , &(storageArray[idx] ) , "bG_scFull5x5_e5x5[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e2x2Ratio"] = idx; 
	outTree->Branch("bG_scFull5x5_e2x2Ratio" , &(storageArray[idx] ) , "bG_scFull5x5_e2x2Ratio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e3x3Ratio"] = idx; 
	outTree->Branch("bG_scFull5x5_e3x3Ratio" , &(storageArray[idx] ) , "bG_scFull5x5_e3x3Ratio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_eMaxRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_eMaxRatio" , &(storageArray[idx] ) , "bG_scFull5x5_eMaxRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e2ndRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_e2ndRatio" , &(storageArray[idx] ) , "bG_scFull5x5_e2ndRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_eTopRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_eTopRatio" , &(storageArray[idx] ) , "bG_scFull5x5_eTopRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_eRightRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_eRightRatio" , &(storageArray[idx] ) , "bG_scFull5x5_eRightRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_eBottomRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_eBottomRatio" , &(storageArray[idx] ) , "bG_scFull5x5_eBottomRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_eLeftRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_eLeftRatio" , &(storageArray[idx] ) , "bG_scFull5x5_eLeftRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e2x5MaxRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_e2x5MaxRatio" , &(storageArray[idx] ) , "bG_scFull5x5_e2x5MaxRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e2x5TopRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_e2x5TopRatio" , &(storageArray[idx] ) , "bG_scFull5x5_e2x5TopRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e2x5RightRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_e2x5RightRatio" , &(storageArray[idx] ) , "bG_scFull5x5_e2x5RightRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e2x5BottomRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_e2x5BottomRatio" , &(storageArray[idx] ) , "bG_scFull5x5_e2x5BottomRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_e2x5LeftRatio"] = idx; 
	outTree->Branch("bG_scFull5x5_e2x5LeftRatio" , &(storageArray[idx] ) , "bG_scFull5x5_e2x5LeftRatio[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_swissCross"] = idx; 
	outTree->Branch("bG_scFull5x5_swissCross" , &(storageArray[idx] ) , "bG_scFull5x5_swissCross[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_r9"] = idx; 
	outTree->Branch("bG_scFull5x5_r9" , &(storageArray[idx] ) , "bG_scFull5x5_r9[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_sigmaIetaIeta"] = idx; 
	outTree->Branch("bG_scFull5x5_sigmaIetaIeta" , &(storageArray[idx] ) , "bG_scFull5x5_sigmaIetaIeta[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_sigmaIetaIphi"] = idx; 
	outTree->Branch("bG_scFull5x5_sigmaIetaIphi" , &(storageArray[idx] ) , "bG_scFull5x5_sigmaIetaIphi[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scFull5x5_sigmaIphiIphi"] = idx; 
	outTree->Branch("bG_scFull5x5_sigmaIphiIphi" , &(storageArray[idx] ) , "bG_scFull5x5_sigmaIphiIphi[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFChIso1"] = idx; 
	outTree->Branch("bG_scPFChIso1" , &(storageArray[idx] ) , "bG_scPFChIso1[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFChIso2"] = idx; 
	outTree->Branch("bG_scPFChIso2" , &(storageArray[idx] ) , "bG_scPFChIso2[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFChIso3"] = idx; 
	outTree->Branch("bG_scPFChIso3" , &(storageArray[idx] ) , "bG_scPFChIso3[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFChIso4"] = idx; 
	outTree->Branch("bG_scPFChIso4" , &(storageArray[idx] ) , "bG_scPFChIso4[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFChIso5"] = idx; 
	outTree->Branch("bG_scPFChIso5" , &(storageArray[idx] ) , "bG_scPFChIso5[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFPhoIso1"] = idx; 
	outTree->Branch("bG_scPFPhoIso1" , &(storageArray[idx] ) , "bG_scPFPhoIso1[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFPhoIso2"] = idx; 
	outTree->Branch("bG_scPFPhoIso2" , &(storageArray[idx] ) , "bG_scPFPhoIso2[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFPhoIso3"] = idx; 
	outTree->Branch("bG_scPFPhoIso3" , &(storageArray[idx] ) , "bG_scPFPhoIso3[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFPhoIso4"] = idx; 
	outTree->Branch("bG_scPFPhoIso4" , &(storageArray[idx] ) , "bG_scPFPhoIso4[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFPhoIso5"] = idx; 
	outTree->Branch("bG_scPFPhoIso5" , &(storageArray[idx] ) , "bG_scPFPhoIso5[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFNeuIso1"] = idx; 
	outTree->Branch("bG_scPFNeuIso1" , &(storageArray[idx] ) , "bG_scPFNeuIso1[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFNeuIso2"] = idx; 
	outTree->Branch("bG_scPFNeuIso2" , &(storageArray[idx] ) , "bG_scPFNeuIso2[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFNeuIso3"] = idx; 
	outTree->Branch("bG_scPFNeuIso3" , &(storageArray[idx] ) , "bG_scPFNeuIso3[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFNeuIso4"] = idx; 
	outTree->Branch("bG_scPFNeuIso4" , &(storageArray[idx] ) , "bG_scPFNeuIso4[nSCPhotons]/D"); 
	idx+=size; 

	candidateMap["bG_scPFNeuIso5"] = idx; 
	outTree->Branch("bG_scPFNeuIso5" , &(storageArray[idx] ) , "bG_scPFNeuIso5[nSCPhotons]/D"); 
	idx+=size; 

	return idx;