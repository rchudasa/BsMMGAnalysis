void genericTreeBranchSelector::FillTheArraysFromVectors()
{
	bG_nscE = bG.scE->size() ;
	bG_nscE = bG_nscE > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE ;
	bG_nscEt = bG.scEt->size() ;
	bG_nscEt = bG_nscEt > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscEt ;
	bG_nscRawE = bG.scRawE->size() ;
	bG_nscRawE = bG_nscRawE > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscRawE ;
	bG_nscEta = bG.scEta->size() ;
	bG_nscEta = bG_nscEta > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscEta ;
	bG_nscPhi = bG.scPhi->size() ;
	bG_nscPhi = bG_nscPhi > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPhi ;
	bG_nscX = bG.scX->size() ;
	bG_nscX = bG_nscX > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscX ;
	bG_nscY = bG.scY->size() ;
	bG_nscY = bG_nscY > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscY ;
	bG_nscZ = bG.scZ->size() ;
	bG_nscZ = bG_nscZ > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscZ ;
	bG_nscEtaWidth = bG.scEtaWidth->size() ;
	bG_nscEtaWidth = bG_nscEtaWidth > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscEtaWidth ;
	bG_nscPhiWidth = bG.scPhiWidth->size() ;
	bG_nscPhiWidth = bG_nscPhiWidth > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPhiWidth ;
	bG_nscRawEt = bG.scRawEt->size() ;
	bG_nscRawEt = bG_nscRawEt > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscRawEt ;
	bG_nscMinDrWithGsfElectornSC_ = bG.scMinDrWithGsfElectornSC_->size() ;
	bG_nscMinDrWithGsfElectornSC_ = bG_nscMinDrWithGsfElectornSC_ > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscMinDrWithGsfElectornSC_ ;
	bG_nscFoundGsfMatch_ = bG.scFoundGsfMatch_->size() ;
	bG_nscFoundGsfMatch_ = bG_nscFoundGsfMatch_ > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFoundGsfMatch_ ;
	bG_nscE5x5 = bG.scE5x5->size() ;
	bG_nscE5x5 = bG_nscE5x5 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE5x5 ;
	bG_nscE2x2Ratio = bG.scE2x2Ratio->size() ;
	bG_nscE2x2Ratio = bG_nscE2x2Ratio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE2x2Ratio ;
	bG_nscE3x3Ratio = bG.scE3x3Ratio->size() ;
	bG_nscE3x3Ratio = bG_nscE3x3Ratio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE3x3Ratio ;
	bG_nscEMaxRatio = bG.scEMaxRatio->size() ;
	bG_nscEMaxRatio = bG_nscEMaxRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscEMaxRatio ;
	bG_nscE2ndRatio = bG.scE2ndRatio->size() ;
	bG_nscE2ndRatio = bG_nscE2ndRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE2ndRatio ;
	bG_nscETopRatio = bG.scETopRatio->size() ;
	bG_nscETopRatio = bG_nscETopRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscETopRatio ;
	bG_nscERightRatio = bG.scERightRatio->size() ;
	bG_nscERightRatio = bG_nscERightRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscERightRatio ;
	bG_nscEBottomRatio = bG.scEBottomRatio->size() ;
	bG_nscEBottomRatio = bG_nscEBottomRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscEBottomRatio ;
	bG_nscELeftRatio = bG.scELeftRatio->size() ;
	bG_nscELeftRatio = bG_nscELeftRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscELeftRatio ;
	bG_nscE2x5MaxRatio = bG.scE2x5MaxRatio->size() ;
	bG_nscE2x5MaxRatio = bG_nscE2x5MaxRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE2x5MaxRatio ;
	bG_nscE2x5TopRatio = bG.scE2x5TopRatio->size() ;
	bG_nscE2x5TopRatio = bG_nscE2x5TopRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE2x5TopRatio ;
	bG_nscE2x5RightRatio = bG.scE2x5RightRatio->size() ;
	bG_nscE2x5RightRatio = bG_nscE2x5RightRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE2x5RightRatio ;
	bG_nscE2x5BottomRatio = bG.scE2x5BottomRatio->size() ;
	bG_nscE2x5BottomRatio = bG_nscE2x5BottomRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE2x5BottomRatio ;
	bG_nscE2x5LeftRatio = bG.scE2x5LeftRatio->size() ;
	bG_nscE2x5LeftRatio = bG_nscE2x5LeftRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscE2x5LeftRatio ;
	bG_nscSwissCross = bG.scSwissCross->size() ;
	bG_nscSwissCross = bG_nscSwissCross > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscSwissCross ;
	bG_nscR9 = bG.scR9->size() ;
	bG_nscR9 = bG_nscR9 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscR9 ;
	bG_nscSigmaIetaIeta = bG.scSigmaIetaIeta->size() ;
	bG_nscSigmaIetaIeta = bG_nscSigmaIetaIeta > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscSigmaIetaIeta ;
	bG_nscSigmaIetaIphi = bG.scSigmaIetaIphi->size() ;
	bG_nscSigmaIetaIphi = bG_nscSigmaIetaIphi > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscSigmaIetaIphi ;
	bG_nscSigmaIphiIphi = bG.scSigmaIphiIphi->size() ;
	bG_nscSigmaIphiIphi = bG_nscSigmaIphiIphi > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscSigmaIphiIphi ;
	bG_nscFull5x5_e5x5 = bG.scFull5x5_e5x5->size() ;
	bG_nscFull5x5_e5x5 = bG_nscFull5x5_e5x5 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e5x5 ;
	bG_nscFull5x5_e2x2Ratio = bG.scFull5x5_e2x2Ratio->size() ;
	bG_nscFull5x5_e2x2Ratio = bG_nscFull5x5_e2x2Ratio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e2x2Ratio ;
	bG_nscFull5x5_e3x3Ratio = bG.scFull5x5_e3x3Ratio->size() ;
	bG_nscFull5x5_e3x3Ratio = bG_nscFull5x5_e3x3Ratio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e3x3Ratio ;
	bG_nscFull5x5_eMaxRatio = bG.scFull5x5_eMaxRatio->size() ;
	bG_nscFull5x5_eMaxRatio = bG_nscFull5x5_eMaxRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_eMaxRatio ;
	bG_nscFull5x5_e2ndRatio = bG.scFull5x5_e2ndRatio->size() ;
	bG_nscFull5x5_e2ndRatio = bG_nscFull5x5_e2ndRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e2ndRatio ;
	bG_nscFull5x5_eTopRatio = bG.scFull5x5_eTopRatio->size() ;
	bG_nscFull5x5_eTopRatio = bG_nscFull5x5_eTopRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_eTopRatio ;
	bG_nscFull5x5_eRightRatio = bG.scFull5x5_eRightRatio->size() ;
	bG_nscFull5x5_eRightRatio = bG_nscFull5x5_eRightRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_eRightRatio ;
	bG_nscFull5x5_eBottomRatio = bG.scFull5x5_eBottomRatio->size() ;
	bG_nscFull5x5_eBottomRatio = bG_nscFull5x5_eBottomRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_eBottomRatio ;
	bG_nscFull5x5_eLeftRatio = bG.scFull5x5_eLeftRatio->size() ;
	bG_nscFull5x5_eLeftRatio = bG_nscFull5x5_eLeftRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_eLeftRatio ;
	bG_nscFull5x5_e2x5MaxRatio = bG.scFull5x5_e2x5MaxRatio->size() ;
	bG_nscFull5x5_e2x5MaxRatio = bG_nscFull5x5_e2x5MaxRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e2x5MaxRatio ;
	bG_nscFull5x5_e2x5TopRatio = bG.scFull5x5_e2x5TopRatio->size() ;
	bG_nscFull5x5_e2x5TopRatio = bG_nscFull5x5_e2x5TopRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e2x5TopRatio ;
	bG_nscFull5x5_e2x5RightRatio = bG.scFull5x5_e2x5RightRatio->size() ;
	bG_nscFull5x5_e2x5RightRatio = bG_nscFull5x5_e2x5RightRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e2x5RightRatio ;
	bG_nscFull5x5_e2x5BottomRatio = bG.scFull5x5_e2x5BottomRatio->size() ;
	bG_nscFull5x5_e2x5BottomRatio = bG_nscFull5x5_e2x5BottomRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e2x5BottomRatio ;
	bG_nscFull5x5_e2x5LeftRatio = bG.scFull5x5_e2x5LeftRatio->size() ;
	bG_nscFull5x5_e2x5LeftRatio = bG_nscFull5x5_e2x5LeftRatio > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_e2x5LeftRatio ;
	bG_nscFull5x5_swissCross = bG.scFull5x5_swissCross->size() ;
	bG_nscFull5x5_swissCross = bG_nscFull5x5_swissCross > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_swissCross ;
	bG_nscFull5x5_r9 = bG.scFull5x5_r9->size() ;
	bG_nscFull5x5_r9 = bG_nscFull5x5_r9 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_r9 ;
	bG_nscFull5x5_sigmaIetaIeta = bG.scFull5x5_sigmaIetaIeta->size() ;
	bG_nscFull5x5_sigmaIetaIeta = bG_nscFull5x5_sigmaIetaIeta > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_sigmaIetaIeta ;
	bG_nscFull5x5_sigmaIetaIphi = bG.scFull5x5_sigmaIetaIphi->size() ;
	bG_nscFull5x5_sigmaIetaIphi = bG_nscFull5x5_sigmaIetaIphi > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_sigmaIetaIphi ;
	bG_nscFull5x5_sigmaIphiIphi = bG.scFull5x5_sigmaIphiIphi->size() ;
	bG_nscFull5x5_sigmaIphiIphi = bG_nscFull5x5_sigmaIphiIphi > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscFull5x5_sigmaIphiIphi ;
	bG_nscNHcalRecHitInDIEta5IPhi5 = bG.scNHcalRecHitInDIEta5IPhi5->size() ;
	bG_nscNHcalRecHitInDIEta5IPhi5 = bG_nscNHcalRecHitInDIEta5IPhi5 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscNHcalRecHitInDIEta5IPhi5 ;
	bG_nscEFromHcalRecHitInDIEta5IPhi5 = bG.scEFromHcalRecHitInDIEta5IPhi5->size() ;
	bG_nscEFromHcalRecHitInDIEta5IPhi5 = bG_nscEFromHcalRecHitInDIEta5IPhi5 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscEFromHcalRecHitInDIEta5IPhi5 ;
	bG_nscNHcalRecHitInDIEta2IPhi2 = bG.scNHcalRecHitInDIEta2IPhi2->size() ;
	bG_nscNHcalRecHitInDIEta2IPhi2 = bG_nscNHcalRecHitInDIEta2IPhi2 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscNHcalRecHitInDIEta2IPhi2 ;
	bG_nscEFromHcalRecHitInDIEta2IPhi2 = bG.scEFromHcalRecHitInDIEta2IPhi2->size() ;
	bG_nscEFromHcalRecHitInDIEta2IPhi2 = bG_nscEFromHcalRecHitInDIEta2IPhi2 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscEFromHcalRecHitInDIEta2IPhi2 ;
	bG_nscPFChIso1 = bG.scPFChIso1->size() ;
	bG_nscPFChIso1 = bG_nscPFChIso1 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFChIso1 ;
	bG_nscPFChIso2 = bG.scPFChIso2->size() ;
	bG_nscPFChIso2 = bG_nscPFChIso2 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFChIso2 ;
	bG_nscPFChIso3 = bG.scPFChIso3->size() ;
	bG_nscPFChIso3 = bG_nscPFChIso3 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFChIso3 ;
	bG_nscPFChIso4 = bG.scPFChIso4->size() ;
	bG_nscPFChIso4 = bG_nscPFChIso4 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFChIso4 ;
	bG_nscPFChIso5 = bG.scPFChIso5->size() ;
	bG_nscPFChIso5 = bG_nscPFChIso5 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFChIso5 ;
	bG_nscPFPhoIso1 = bG.scPFPhoIso1->size() ;
	bG_nscPFPhoIso1 = bG_nscPFPhoIso1 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFPhoIso1 ;
	bG_nscPFPhoIso2 = bG.scPFPhoIso2->size() ;
	bG_nscPFPhoIso2 = bG_nscPFPhoIso2 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFPhoIso2 ;
	bG_nscPFPhoIso3 = bG.scPFPhoIso3->size() ;
	bG_nscPFPhoIso3 = bG_nscPFPhoIso3 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFPhoIso3 ;
	bG_nscPFPhoIso4 = bG.scPFPhoIso4->size() ;
	bG_nscPFPhoIso4 = bG_nscPFPhoIso4 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFPhoIso4 ;
	bG_nscPFPhoIso5 = bG.scPFPhoIso5->size() ;
	bG_nscPFPhoIso5 = bG_nscPFPhoIso5 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFPhoIso5 ;
	bG_nscPFNeuIso1 = bG.scPFNeuIso1->size() ;
	bG_nscPFNeuIso1 = bG_nscPFNeuIso1 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFNeuIso1 ;
	bG_nscPFNeuIso2 = bG.scPFNeuIso2->size() ;
	bG_nscPFNeuIso2 = bG_nscPFNeuIso2 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFNeuIso2 ;
	bG_nscPFNeuIso3 = bG.scPFNeuIso3->size() ;
	bG_nscPFNeuIso3 = bG_nscPFNeuIso3 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFNeuIso3 ;
	bG_nscPFNeuIso4 = bG.scPFNeuIso4->size() ;
	bG_nscPFNeuIso4 = bG_nscPFNeuIso4 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFNeuIso4 ;
	bG_nscPFNeuIso5 = bG.scPFNeuIso5->size() ;
	bG_nscPFNeuIso5 = bG_nscPFNeuIso5 > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : bG_nscPFNeuIso5 ;


	for( UInt_t i=0; i< UInt_t ( bG_nscE) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE[i] = bG.scE->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscEt) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scEt[i] = bG.scEt->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscRawE) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scRawE[i] = bG.scRawE->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscEta) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scEta[i] = bG.scEta->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPhi) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPhi[i] = bG.scPhi->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscX) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scX[i] = bG.scX->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscY) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scY[i] = bG.scY->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscZ) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scZ[i] = bG.scZ->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscEtaWidth) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scEtaWidth[i] = bG.scEtaWidth->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPhiWidth) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPhiWidth[i] = bG.scPhiWidth->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscRawEt) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scRawEt[i] = bG.scRawEt->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscMinDrWithGsfElectornSC_) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scMinDrWithGsfElectornSC_[i] = bG.scMinDrWithGsfElectornSC_->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFoundGsfMatch_) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFoundGsfMatch_[i] = bG.scFoundGsfMatch_->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE5x5) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE5x5[i] = bG.scE5x5->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE2x2Ratio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE2x2Ratio[i] = bG.scE2x2Ratio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE3x3Ratio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE3x3Ratio[i] = bG.scE3x3Ratio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscEMaxRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scEMaxRatio[i] = bG.scEMaxRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE2ndRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE2ndRatio[i] = bG.scE2ndRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscETopRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scETopRatio[i] = bG.scETopRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscERightRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scERightRatio[i] = bG.scERightRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscEBottomRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scEBottomRatio[i] = bG.scEBottomRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscELeftRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scELeftRatio[i] = bG.scELeftRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE2x5MaxRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE2x5MaxRatio[i] = bG.scE2x5MaxRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE2x5TopRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE2x5TopRatio[i] = bG.scE2x5TopRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE2x5RightRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE2x5RightRatio[i] = bG.scE2x5RightRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE2x5BottomRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE2x5BottomRatio[i] = bG.scE2x5BottomRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscE2x5LeftRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scE2x5LeftRatio[i] = bG.scE2x5LeftRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscSwissCross) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scSwissCross[i] = bG.scSwissCross->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscR9) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scR9[i] = bG.scR9->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscSigmaIetaIeta) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scSigmaIetaIeta[i] = bG.scSigmaIetaIeta->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscSigmaIetaIphi) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scSigmaIetaIphi[i] = bG.scSigmaIetaIphi->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscSigmaIphiIphi) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scSigmaIphiIphi[i] = bG.scSigmaIphiIphi->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e5x5) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e5x5[i] = bG.scFull5x5_e5x5->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e2x2Ratio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e2x2Ratio[i] = bG.scFull5x5_e2x2Ratio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e3x3Ratio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e3x3Ratio[i] = bG.scFull5x5_e3x3Ratio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_eMaxRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_eMaxRatio[i] = bG.scFull5x5_eMaxRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e2ndRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e2ndRatio[i] = bG.scFull5x5_e2ndRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_eTopRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_eTopRatio[i] = bG.scFull5x5_eTopRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_eRightRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_eRightRatio[i] = bG.scFull5x5_eRightRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_eBottomRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_eBottomRatio[i] = bG.scFull5x5_eBottomRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_eLeftRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_eLeftRatio[i] = bG.scFull5x5_eLeftRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e2x5MaxRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e2x5MaxRatio[i] = bG.scFull5x5_e2x5MaxRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e2x5TopRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e2x5TopRatio[i] = bG.scFull5x5_e2x5TopRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e2x5RightRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e2x5RightRatio[i] = bG.scFull5x5_e2x5RightRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e2x5BottomRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e2x5BottomRatio[i] = bG.scFull5x5_e2x5BottomRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_e2x5LeftRatio) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_e2x5LeftRatio[i] = bG.scFull5x5_e2x5LeftRatio->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_swissCross) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_swissCross[i] = bG.scFull5x5_swissCross->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_r9) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_r9[i] = bG.scFull5x5_r9->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_sigmaIetaIeta) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_sigmaIetaIeta[i] = bG.scFull5x5_sigmaIetaIeta->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_sigmaIetaIphi) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_sigmaIetaIphi[i] = bG.scFull5x5_sigmaIetaIphi->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscFull5x5_sigmaIphiIphi) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scFull5x5_sigmaIphiIphi[i] = bG.scFull5x5_sigmaIphiIphi->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscNHcalRecHitInDIEta5IPhi5) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scNHcalRecHitInDIEta5IPhi5[i] = bG.scNHcalRecHitInDIEta5IPhi5->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscEFromHcalRecHitInDIEta5IPhi5) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scEFromHcalRecHitInDIEta5IPhi5[i] = bG.scEFromHcalRecHitInDIEta5IPhi5->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscNHcalRecHitInDIEta2IPhi2) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scNHcalRecHitInDIEta2IPhi2[i] = bG.scNHcalRecHitInDIEta2IPhi2->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscEFromHcalRecHitInDIEta2IPhi2) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scEFromHcalRecHitInDIEta2IPhi2[i] = bG.scEFromHcalRecHitInDIEta2IPhi2->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFChIso1) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFChIso1[i] = bG.scPFChIso1->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFChIso2) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFChIso2[i] = bG.scPFChIso2->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFChIso3) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFChIso3[i] = bG.scPFChIso3->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFChIso4) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFChIso4[i] = bG.scPFChIso4->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFChIso5) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFChIso5[i] = bG.scPFChIso5->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFPhoIso1) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFPhoIso1[i] = bG.scPFPhoIso1->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFPhoIso2) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFPhoIso2[i] = bG.scPFPhoIso2->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFPhoIso3) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFPhoIso3[i] = bG.scPFPhoIso3->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFPhoIso4) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFPhoIso4[i] = bG.scPFPhoIso4->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFPhoIso5) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFPhoIso5[i] = bG.scPFPhoIso5->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFNeuIso1) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFNeuIso1[i] = bG.scPFNeuIso1->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFNeuIso2) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFNeuIso2[i] = bG.scPFNeuIso2->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFNeuIso3) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFNeuIso3[i] = bG.scPFNeuIso3->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFNeuIso4) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFNeuIso4[i] = bG.scPFNeuIso4->at(i);

	}

	for( UInt_t i=0; i< UInt_t ( bG_nscPFNeuIso5) ; i++ ) { 
		 if( i>= MAX_ARRAY_SIZE ) {    break;}
 		bG_scPFNeuIso5[i] = bG.scPFNeuIso5->at(i);

	}


}

void genericTreeBranchSelector::setCompiledTreeBranches() 
{

	_compiledTree->Branch("b5_run",&	( b5.run ));
	_compiledTree->Branch("b5_luminosityBlock",&	( b5.luminosityBlock ));
	_compiledTree->Branch("b5_event",&	( b5.event ));
	_compiledTree->Branch("b5_nMuonId",&	( b5.nMuonId ));
	_compiledTree->Branch("b5_MuonId_chi2LocalPosition",	( b5.MuonId_chi2LocalPosition ),"b5_MuonId_chi2LocalPosition[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_glbNormChi2",	( b5.MuonId_glbNormChi2 ),"b5_MuonId_glbNormChi2[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_glbTrackProbability",	( b5.MuonId_glbTrackProbability ),"b5_MuonId_glbTrackProbability[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match1_dX",	( b5.MuonId_match1_dX ),"b5_MuonId_match1_dX[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match1_dY",	( b5.MuonId_match1_dY ),"b5_MuonId_match1_dY[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match1_pullDxDz",	( b5.MuonId_match1_pullDxDz ),"b5_MuonId_match1_pullDxDz[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match1_pullDyDz",	( b5.MuonId_match1_pullDyDz ),"b5_MuonId_match1_pullDyDz[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match1_pullX",	( b5.MuonId_match1_pullX ),"b5_MuonId_match1_pullX[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match1_pullY",	( b5.MuonId_match1_pullY ),"b5_MuonId_match1_pullY[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match2_dX",	( b5.MuonId_match2_dX ),"b5_MuonId_match2_dX[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match2_dY",	( b5.MuonId_match2_dY ),"b5_MuonId_match2_dY[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match2_pullDxDz",	( b5.MuonId_match2_pullDxDz ),"b5_MuonId_match2_pullDxDz[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match2_pullDyDz",	( b5.MuonId_match2_pullDyDz ),"b5_MuonId_match2_pullDyDz[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match2_pullX",	( b5.MuonId_match2_pullX ),"b5_MuonId_match2_pullX[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_match2_pullY",	( b5.MuonId_match2_pullY ),"b5_MuonId_match2_pullY[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_newSoftMuonMva",	( b5.MuonId_newSoftMuonMva ),"b5_MuonId_newSoftMuonMva[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_trkKink",	( b5.MuonId_trkKink ),"b5_MuonId_trkKink[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_trkValidFrac",	( b5.MuonId_trkValidFrac ),"b5_MuonId_trkValidFrac[b5_nMuonId]/F");
	_compiledTree->Branch("b5_MuonId_highPurity",	( b5.MuonId_highPurity ),"b5_MuonId_highPurity[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_nLostHitsInner",	( b5.MuonId_nLostHitsInner ),"b5_MuonId_nLostHitsInner[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_nLostHitsOn",	( b5.MuonId_nLostHitsOn ),"b5_MuonId_nLostHitsOn[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_nLostHitsOuter",	( b5.MuonId_nLostHitsOuter ),"b5_MuonId_nLostHitsOuter[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_nPixels",	( b5.MuonId_nPixels ),"b5_MuonId_nPixels[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_nValidHits",	( b5.MuonId_nValidHits ),"b5_MuonId_nValidHits[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_trkLayers",	( b5.MuonId_trkLayers ),"b5_MuonId_trkLayers[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_trkLostLayersInner",	( b5.MuonId_trkLostLayersInner ),"b5_MuonId_trkLostLayersInner[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_trkLostLayersOn",	( b5.MuonId_trkLostLayersOn ),"b5_MuonId_trkLostLayersOn[b5_nMuonId]/I");
	_compiledTree->Branch("b5_MuonId_trkLostLayersOuter",	( b5.MuonId_trkLostLayersOuter ),"b5_MuonId_trkLostLayersOuter[b5_nMuonId]/I");
	_compiledTree->Branch("b5_nmm",&	( b5.nmm ));
	_compiledTree->Branch("b5_mm_bdt",	( b5.mm_bdt ),"b5_mm_bdt[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_doca",	( b5.mm_doca ),"b5_mm_doca[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_docatrk",	( b5.mm_docatrk ),"b5_mm_docatrk[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_iso",	( b5.mm_iso ),"b5_mm_iso[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kal_lxy",	( b5.mm_kal_lxy ),"b5_mm_kal_lxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kal_mass",	( b5.mm_kal_mass ),"b5_mm_kal_mass[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kal_slxy",	( b5.mm_kal_slxy ),"b5_mm_kal_slxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kal_vtx_prob",	( b5.mm_kal_vtx_prob ),"b5_mm_kal_vtx_prob[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_alpha",	( b5.mm_kin_alpha ),"b5_mm_kin_alpha[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_alphaBS",	( b5.mm_kin_alphaBS ),"b5_mm_kin_alphaBS[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_alphaBSErr",	( b5.mm_kin_alphaBSErr ),"b5_mm_kin_alphaBSErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_alphaErr",	( b5.mm_kin_alphaErr ),"b5_mm_kin_alphaErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_eta",	( b5.mm_kin_eta ),"b5_mm_kin_eta[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_l3d",	( b5.mm_kin_l3d ),"b5_mm_kin_l3d[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_lxy",	( b5.mm_kin_lxy ),"b5_mm_kin_lxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_mass",	( b5.mm_kin_mass ),"b5_mm_kin_mass[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_massErr",	( b5.mm_kin_massErr ),"b5_mm_kin_massErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_mu1eta",	( b5.mm_kin_mu1eta ),"b5_mm_kin_mu1eta[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_mu1phi",	( b5.mm_kin_mu1phi ),"b5_mm_kin_mu1phi[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_mu1pt",	( b5.mm_kin_mu1pt ),"b5_mm_kin_mu1pt[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_mu2eta",	( b5.mm_kin_mu2eta ),"b5_mm_kin_mu2eta[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_mu2phi",	( b5.mm_kin_mu2phi ),"b5_mm_kin_mu2phi[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_mu2pt",	( b5.mm_kin_mu2pt ),"b5_mm_kin_mu2pt[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_phi",	( b5.mm_kin_phi ),"b5_mm_kin_phi[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pt",	( b5.mm_kin_pt ),"b5_mm_kin_pt[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pv2ip",	( b5.mm_kin_pv2ip ),"b5_mm_kin_pv2ip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pv2ipErr",	( b5.mm_kin_pv2ipErr ),"b5_mm_kin_pv2ipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pv2lip",	( b5.mm_kin_pv2lip ),"b5_mm_kin_pv2lip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pv2lipErr",	( b5.mm_kin_pv2lipErr ),"b5_mm_kin_pv2lipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pv2lipSig",	( b5.mm_kin_pv2lipSig ),"b5_mm_kin_pv2lipSig[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pv_z",	( b5.mm_kin_pv_z ),"b5_mm_kin_pv_z[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pv_zErr",	( b5.mm_kin_pv_zErr ),"b5_mm_kin_pv_zErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pvip",	( b5.mm_kin_pvip ),"b5_mm_kin_pvip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pvipErr",	( b5.mm_kin_pvipErr ),"b5_mm_kin_pvipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pvlip",	( b5.mm_kin_pvlip ),"b5_mm_kin_pvlip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pvlipErr",	( b5.mm_kin_pvlipErr ),"b5_mm_kin_pvlipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_pvlipSig",	( b5.mm_kin_pvlipSig ),"b5_mm_kin_pvlipSig[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_sl3d",	( b5.mm_kin_sl3d ),"b5_mm_kin_sl3d[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_slxy",	( b5.mm_kin_slxy ),"b5_mm_kin_slxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_spv2ip",	( b5.mm_kin_spv2ip ),"b5_mm_kin_spv2ip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_spvip",	( b5.mm_kin_spvip ),"b5_mm_kin_spvip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_tau",	( b5.mm_kin_tau ),"b5_mm_kin_tau[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_taue",	( b5.mm_kin_taue ),"b5_mm_kin_taue[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_tauxy",	( b5.mm_kin_tauxy ),"b5_mm_kin_tauxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_tauxye",	( b5.mm_kin_tauxye ),"b5_mm_kin_tauxye[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_chi2dof",	( b5.mm_kin_vtx_chi2dof ),"b5_mm_kin_vtx_chi2dof[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_prob",	( b5.mm_kin_vtx_prob ),"b5_mm_kin_vtx_prob[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_x",	( b5.mm_kin_vtx_x ),"b5_mm_kin_vtx_x[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_xErr",	( b5.mm_kin_vtx_xErr ),"b5_mm_kin_vtx_xErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_y",	( b5.mm_kin_vtx_y ),"b5_mm_kin_vtx_y[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_yErr",	( b5.mm_kin_vtx_yErr ),"b5_mm_kin_vtx_yErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_z",	( b5.mm_kin_vtx_z ),"b5_mm_kin_vtx_z[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kin_vtx_zErr",	( b5.mm_kin_vtx_zErr ),"b5_mm_kin_vtx_zErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_alpha",	( b5.mm_kinpc_alpha ),"b5_mm_kinpc_alpha[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_alphaBS",	( b5.mm_kinpc_alphaBS ),"b5_mm_kinpc_alphaBS[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_alphaBSErr",	( b5.mm_kinpc_alphaBSErr ),"b5_mm_kinpc_alphaBSErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_alphaErr",	( b5.mm_kinpc_alphaErr ),"b5_mm_kinpc_alphaErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_eta",	( b5.mm_kinpc_eta ),"b5_mm_kinpc_eta[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_l3d",	( b5.mm_kinpc_l3d ),"b5_mm_kinpc_l3d[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_lxy",	( b5.mm_kinpc_lxy ),"b5_mm_kinpc_lxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_mass",	( b5.mm_kinpc_mass ),"b5_mm_kinpc_mass[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_massErr",	( b5.mm_kinpc_massErr ),"b5_mm_kinpc_massErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_phi",	( b5.mm_kinpc_phi ),"b5_mm_kinpc_phi[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pt",	( b5.mm_kinpc_pt ),"b5_mm_kinpc_pt[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pv2ip",	( b5.mm_kinpc_pv2ip ),"b5_mm_kinpc_pv2ip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pv2ipErr",	( b5.mm_kinpc_pv2ipErr ),"b5_mm_kinpc_pv2ipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pv2lip",	( b5.mm_kinpc_pv2lip ),"b5_mm_kinpc_pv2lip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pv2lipErr",	( b5.mm_kinpc_pv2lipErr ),"b5_mm_kinpc_pv2lipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pv2lipSig",	( b5.mm_kinpc_pv2lipSig ),"b5_mm_kinpc_pv2lipSig[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pv_z",	( b5.mm_kinpc_pv_z ),"b5_mm_kinpc_pv_z[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pv_zErr",	( b5.mm_kinpc_pv_zErr ),"b5_mm_kinpc_pv_zErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pvip",	( b5.mm_kinpc_pvip ),"b5_mm_kinpc_pvip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pvipErr",	( b5.mm_kinpc_pvipErr ),"b5_mm_kinpc_pvipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pvlip",	( b5.mm_kinpc_pvlip ),"b5_mm_kinpc_pvlip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pvlipErr",	( b5.mm_kinpc_pvlipErr ),"b5_mm_kinpc_pvlipErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_pvlipSig",	( b5.mm_kinpc_pvlipSig ),"b5_mm_kinpc_pvlipSig[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_sl3d",	( b5.mm_kinpc_sl3d ),"b5_mm_kinpc_sl3d[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_slxy",	( b5.mm_kinpc_slxy ),"b5_mm_kinpc_slxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_spv2ip",	( b5.mm_kinpc_spv2ip ),"b5_mm_kinpc_spv2ip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_spvip",	( b5.mm_kinpc_spvip ),"b5_mm_kinpc_spvip[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_tau",	( b5.mm_kinpc_tau ),"b5_mm_kinpc_tau[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_taue",	( b5.mm_kinpc_taue ),"b5_mm_kinpc_taue[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_tauxy",	( b5.mm_kinpc_tauxy ),"b5_mm_kinpc_tauxy[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_tauxye",	( b5.mm_kinpc_tauxye ),"b5_mm_kinpc_tauxye[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_chi2dof",	( b5.mm_kinpc_vtx_chi2dof ),"b5_mm_kinpc_vtx_chi2dof[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_prob",	( b5.mm_kinpc_vtx_prob ),"b5_mm_kinpc_vtx_prob[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_x",	( b5.mm_kinpc_vtx_x ),"b5_mm_kinpc_vtx_x[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_xErr",	( b5.mm_kinpc_vtx_xErr ),"b5_mm_kinpc_vtx_xErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_y",	( b5.mm_kinpc_vtx_y ),"b5_mm_kinpc_vtx_y[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_yErr",	( b5.mm_kinpc_vtx_yErr ),"b5_mm_kinpc_vtx_yErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_z",	( b5.mm_kinpc_vtx_z ),"b5_mm_kinpc_vtx_z[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_kinpc_vtx_zErr",	( b5.mm_kinpc_vtx_zErr ),"b5_mm_kinpc_vtx_zErr[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_m1iso",	( b5.mm_m1iso ),"b5_mm_m1iso[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_m2iso",	( b5.mm_m2iso ),"b5_mm_m2iso[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mass",	( b5.mm_mass ),"b5_mm_mass[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mu1_eta",	( b5.mm_mu1_eta ),"b5_mm_mu1_eta[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mu1_phi",	( b5.mm_mu1_phi ),"b5_mm_mu1_phi[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mu1_pt",	( b5.mm_mu1_pt ),"b5_mm_mu1_pt[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mu2_eta",	( b5.mm_mu2_eta ),"b5_mm_mu2_eta[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mu2_phi",	( b5.mm_mu2_phi ),"b5_mm_mu2_phi[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mu2_pt",	( b5.mm_mu2_pt ),"b5_mm_mu2_pt[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_mva",	( b5.mm_mva ),"b5_mm_mva[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_otherVtxMaxProb",	( b5.mm_otherVtxMaxProb ),"b5_mm_otherVtxMaxProb[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_otherVtxMaxProb1",	( b5.mm_otherVtxMaxProb1 ),"b5_mm_otherVtxMaxProb1[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_otherVtxMaxProb2",	( b5.mm_otherVtxMaxProb2 ),"b5_mm_otherVtxMaxProb2[b5_nmm]/F");
	_compiledTree->Branch("b5_mm_closetrk",	( b5.mm_closetrk ),"b5_mm_closetrk[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_closetrks1",	( b5.mm_closetrks1 ),"b5_mm_closetrks1[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_closetrks2",	( b5.mm_closetrks2 ),"b5_mm_closetrks2[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_closetrks3",	( b5.mm_closetrks3 ),"b5_mm_closetrks3[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_kal_valid",	( b5.mm_kal_valid ),"b5_mm_kal_valid[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_kin_valid",	( b5.mm_kin_valid ),"b5_mm_kin_valid[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_kinpc_valid",	( b5.mm_kinpc_valid ),"b5_mm_kinpc_valid[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_mu1_index",	( b5.mm_mu1_index ),"b5_mm_mu1_index[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_mu1_pdgId",	( b5.mm_mu1_pdgId ),"b5_mm_mu1_pdgId[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_mu2_index",	( b5.mm_mu2_index ),"b5_mm_mu2_index[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_mu2_pdgId",	( b5.mm_mu2_pdgId ),"b5_mm_mu2_pdgId[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_nBMTrks",	( b5.mm_nBMTrks ),"b5_mm_nBMTrks[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_nDisTrks",	( b5.mm_nDisTrks ),"b5_mm_nDisTrks[b5_nmm]/I");
	_compiledTree->Branch("b5_mm_nTrks",	( b5.mm_nTrks ),"b5_mm_nTrks[b5_nmm]/I");
	_compiledTree->Branch("b5_nd0",&	( b5.nd0 ));
	_compiledTree->Branch("b5_d0_doca",	( b5.d0_doca ),"b5_d0_doca[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kaon_eta",	( b5.d0_kaon_eta ),"b5_d0_kaon_eta[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kaon_phi",	( b5.d0_kaon_phi ),"b5_d0_kaon_phi[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kaon_pt",	( b5.d0_kaon_pt ),"b5_d0_kaon_pt[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kaon_sip",	( b5.d0_kaon_sip ),"b5_d0_kaon_sip[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_cosAlphaXY",	( b5.d0_kin_cosAlphaXY ),"b5_d0_kin_cosAlphaXY[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_eta",	( b5.d0_kin_eta ),"b5_d0_kin_eta[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_lxy",	( b5.d0_kin_lxy ),"b5_d0_kin_lxy[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_mass",	( b5.d0_kin_mass ),"b5_d0_kin_mass[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_massErr",	( b5.d0_kin_massErr ),"b5_d0_kin_massErr[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_phi",	( b5.d0_kin_phi ),"b5_d0_kin_phi[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_pt",	( b5.d0_kin_pt ),"b5_d0_kin_pt[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_sipBS",	( b5.d0_kin_sipBS ),"b5_d0_kin_sipBS[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_sipPV",	( b5.d0_kin_sipPV ),"b5_d0_kin_sipPV[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_slxy",	( b5.d0_kin_slxy ),"b5_d0_kin_slxy[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_vtx_chi2dof",	( b5.d0_kin_vtx_chi2dof ),"b5_d0_kin_vtx_chi2dof[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kin_vtx_prob",	( b5.d0_kin_vtx_prob ),"b5_d0_kin_vtx_prob[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_mass",	( b5.d0_mass ),"b5_d0_mass[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_pion_eta",	( b5.d0_pion_eta ),"b5_d0_pion_eta[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_pion_phi",	( b5.d0_pion_phi ),"b5_d0_pion_phi[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_pion_pt",	( b5.d0_pion_pt ),"b5_d0_pion_pt[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_pion_sip",	( b5.d0_pion_sip ),"b5_d0_pion_sip[b5_nd0]/F");
	_compiledTree->Branch("b5_d0_kaon_mu_index",	( b5.d0_kaon_mu_index ),"b5_d0_kaon_mu_index[b5_nd0]/I");
	_compiledTree->Branch("b5_d0_kin_valid",	( b5.d0_kin_valid ),"b5_d0_kin_valid[b5_nd0]/I");
	_compiledTree->Branch("b5_d0_pion_mu_index",	( b5.d0_pion_mu_index ),"b5_d0_pion_mu_index[b5_nd0]/I");
	_compiledTree->Branch("b5_nks",&	( b5.nks ));
	_compiledTree->Branch("b5_ks_doca",	( b5.ks_doca ),"b5_ks_doca[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_cosAlphaXY",	( b5.ks_kin_cosAlphaXY ),"b5_ks_kin_cosAlphaXY[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_eta",	( b5.ks_kin_eta ),"b5_ks_kin_eta[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_lxy",	( b5.ks_kin_lxy ),"b5_ks_kin_lxy[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_mass",	( b5.ks_kin_mass ),"b5_ks_kin_mass[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_massErr",	( b5.ks_kin_massErr ),"b5_ks_kin_massErr[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_phi",	( b5.ks_kin_phi ),"b5_ks_kin_phi[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_pt",	( b5.ks_kin_pt ),"b5_ks_kin_pt[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_sipBS",	( b5.ks_kin_sipBS ),"b5_ks_kin_sipBS[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_sipPV",	( b5.ks_kin_sipPV ),"b5_ks_kin_sipPV[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_slxy",	( b5.ks_kin_slxy ),"b5_ks_kin_slxy[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_vtx_chi2dof",	( b5.ks_kin_vtx_chi2dof ),"b5_ks_kin_vtx_chi2dof[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_vtx_prob",	( b5.ks_kin_vtx_prob ),"b5_ks_kin_vtx_prob[b5_nks]/F");
	_compiledTree->Branch("b5_ks_mass",	( b5.ks_mass ),"b5_ks_mass[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk1_eta",	( b5.ks_trk1_eta ),"b5_ks_trk1_eta[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk1_phi",	( b5.ks_trk1_phi ),"b5_ks_trk1_phi[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk1_pt",	( b5.ks_trk1_pt ),"b5_ks_trk1_pt[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk1_sip",	( b5.ks_trk1_sip ),"b5_ks_trk1_sip[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk2_eta",	( b5.ks_trk2_eta ),"b5_ks_trk2_eta[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk2_phi",	( b5.ks_trk2_phi ),"b5_ks_trk2_phi[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk2_pt",	( b5.ks_trk2_pt ),"b5_ks_trk2_pt[b5_nks]/F");
	_compiledTree->Branch("b5_ks_trk2_sip",	( b5.ks_trk2_sip ),"b5_ks_trk2_sip[b5_nks]/F");
	_compiledTree->Branch("b5_ks_kin_valid",	( b5.ks_kin_valid ),"b5_ks_kin_valid[b5_nks]/I");
	_compiledTree->Branch("b5_ks_trk1_mu_index",	( b5.ks_trk1_mu_index ),"b5_ks_trk1_mu_index[b5_nks]/I");
	_compiledTree->Branch("b5_ks_trk2_mu_index",	( b5.ks_trk2_mu_index ),"b5_ks_trk2_mu_index[b5_nks]/I");
	_compiledTree->Branch("b5_nlambda",&	( b5.nlambda ));
	_compiledTree->Branch("b5_lambda_doca",	( b5.lambda_doca ),"b5_lambda_doca[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_cosAlphaXY",	( b5.lambda_kin_cosAlphaXY ),"b5_lambda_kin_cosAlphaXY[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_eta",	( b5.lambda_kin_eta ),"b5_lambda_kin_eta[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_lxy",	( b5.lambda_kin_lxy ),"b5_lambda_kin_lxy[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_mass",	( b5.lambda_kin_mass ),"b5_lambda_kin_mass[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_massErr",	( b5.lambda_kin_massErr ),"b5_lambda_kin_massErr[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_phi",	( b5.lambda_kin_phi ),"b5_lambda_kin_phi[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_pt",	( b5.lambda_kin_pt ),"b5_lambda_kin_pt[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_sipBS",	( b5.lambda_kin_sipBS ),"b5_lambda_kin_sipBS[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_sipPV",	( b5.lambda_kin_sipPV ),"b5_lambda_kin_sipPV[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_slxy",	( b5.lambda_kin_slxy ),"b5_lambda_kin_slxy[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_vtx_chi2dof",	( b5.lambda_kin_vtx_chi2dof ),"b5_lambda_kin_vtx_chi2dof[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_vtx_prob",	( b5.lambda_kin_vtx_prob ),"b5_lambda_kin_vtx_prob[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_mass",	( b5.lambda_mass ),"b5_lambda_mass[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_pion_eta",	( b5.lambda_pion_eta ),"b5_lambda_pion_eta[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_pion_phi",	( b5.lambda_pion_phi ),"b5_lambda_pion_phi[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_pion_pt",	( b5.lambda_pion_pt ),"b5_lambda_pion_pt[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_pion_sip",	( b5.lambda_pion_sip ),"b5_lambda_pion_sip[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_proton_eta",	( b5.lambda_proton_eta ),"b5_lambda_proton_eta[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_proton_phi",	( b5.lambda_proton_phi ),"b5_lambda_proton_phi[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_proton_pt",	( b5.lambda_proton_pt ),"b5_lambda_proton_pt[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_proton_sip",	( b5.lambda_proton_sip ),"b5_lambda_proton_sip[b5_nlambda]/F");
	_compiledTree->Branch("b5_lambda_kin_valid",	( b5.lambda_kin_valid ),"b5_lambda_kin_valid[b5_nlambda]/I");
	_compiledTree->Branch("b5_lambda_pion_mu_index",	( b5.lambda_pion_mu_index ),"b5_lambda_pion_mu_index[b5_nlambda]/I");
	_compiledTree->Branch("b5_lambda_proton_mu_index",	( b5.lambda_proton_mu_index ),"b5_lambda_proton_mu_index[b5_nlambda]/I");
	_compiledTree->Branch("b5_nphi",&	( b5.nphi ));
	_compiledTree->Branch("b5_phi_doca",	( b5.phi_doca ),"b5_phi_doca[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_cosAlphaXY",	( b5.phi_ds_cosAlphaXY ),"b5_phi_ds_cosAlphaXY[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_eta",	( b5.phi_ds_eta ),"b5_phi_ds_eta[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_lxy",	( b5.phi_ds_lxy ),"b5_phi_ds_lxy[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_mass",	( b5.phi_ds_mass ),"b5_phi_ds_mass[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_massErr",	( b5.phi_ds_massErr ),"b5_phi_ds_massErr[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_phi",	( b5.phi_ds_phi ),"b5_phi_ds_phi[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_pion_eta",	( b5.phi_ds_pion_eta ),"b5_phi_ds_pion_eta[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_pion_mu_index",	( b5.phi_ds_pion_mu_index ),"b5_phi_ds_pion_mu_index[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_pion_phi",	( b5.phi_ds_pion_phi ),"b5_phi_ds_pion_phi[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_pion_pt",	( b5.phi_ds_pion_pt ),"b5_phi_ds_pion_pt[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_pt",	( b5.phi_ds_pt ),"b5_phi_ds_pt[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_sipBS",	( b5.phi_ds_sipBS ),"b5_phi_ds_sipBS[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_sipPV",	( b5.phi_ds_sipPV ),"b5_phi_ds_sipPV[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_slxy",	( b5.phi_ds_slxy ),"b5_phi_ds_slxy[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_vtx_chi2dof",	( b5.phi_ds_vtx_chi2dof ),"b5_phi_ds_vtx_chi2dof[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_ds_vtx_prob",	( b5.phi_ds_vtx_prob ),"b5_phi_ds_vtx_prob[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_cosAlphaXY",	( b5.phi_kin_cosAlphaXY ),"b5_phi_kin_cosAlphaXY[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_eta",	( b5.phi_kin_eta ),"b5_phi_kin_eta[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_lxy",	( b5.phi_kin_lxy ),"b5_phi_kin_lxy[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_mass",	( b5.phi_kin_mass ),"b5_phi_kin_mass[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_massErr",	( b5.phi_kin_massErr ),"b5_phi_kin_massErr[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_phi",	( b5.phi_kin_phi ),"b5_phi_kin_phi[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_pt",	( b5.phi_kin_pt ),"b5_phi_kin_pt[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_sipBS",	( b5.phi_kin_sipBS ),"b5_phi_kin_sipBS[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_sipPV",	( b5.phi_kin_sipPV ),"b5_phi_kin_sipPV[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_slxy",	( b5.phi_kin_slxy ),"b5_phi_kin_slxy[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_vtx_chi2dof",	( b5.phi_kin_vtx_chi2dof ),"b5_phi_kin_vtx_chi2dof[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_vtx_prob",	( b5.phi_kin_vtx_prob ),"b5_phi_kin_vtx_prob[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_mass",	( b5.phi_mass ),"b5_phi_mass[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk1_eta",	( b5.phi_trk1_eta ),"b5_phi_trk1_eta[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk1_phi",	( b5.phi_trk1_phi ),"b5_phi_trk1_phi[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk1_pt",	( b5.phi_trk1_pt ),"b5_phi_trk1_pt[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk1_sip",	( b5.phi_trk1_sip ),"b5_phi_trk1_sip[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk2_eta",	( b5.phi_trk2_eta ),"b5_phi_trk2_eta[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk2_phi",	( b5.phi_trk2_phi ),"b5_phi_trk2_phi[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk2_pt",	( b5.phi_trk2_pt ),"b5_phi_trk2_pt[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_trk2_sip",	( b5.phi_trk2_sip ),"b5_phi_trk2_sip[b5_nphi]/F");
	_compiledTree->Branch("b5_phi_kin_valid",	( b5.phi_kin_valid ),"b5_phi_kin_valid[b5_nphi]/I");
	_compiledTree->Branch("b5_phi_trk1_mu_index",	( b5.phi_trk1_mu_index ),"b5_phi_trk1_mu_index[b5_nphi]/I");
	_compiledTree->Branch("b5_phi_trk2_mu_index",	( b5.phi_trk2_mu_index ),"b5_phi_trk2_mu_index[b5_nphi]/I");
	_compiledTree->Branch("b5_CaloMET_phi",&	( b5.CaloMET_phi ));
	_compiledTree->Branch("b5_CaloMET_pt",&	( b5.CaloMET_pt ));
	_compiledTree->Branch("b5_CaloMET_sumEt",&	( b5.CaloMET_sumEt ));
	_compiledTree->Branch("b5_ChsMET_phi",&	( b5.ChsMET_phi ));
	_compiledTree->Branch("b5_ChsMET_pt",&	( b5.ChsMET_pt ));
	_compiledTree->Branch("b5_ChsMET_sumEt",&	( b5.ChsMET_sumEt ));
	_compiledTree->Branch("b5_nCorrT1METJet",&	( b5.nCorrT1METJet ));
	_compiledTree->Branch("b5_CorrT1METJet_area",	( b5.CorrT1METJet_area ),"b5_CorrT1METJet_area[b5_nCorrT1METJet]/F");
	_compiledTree->Branch("b5_CorrT1METJet_eta",	( b5.CorrT1METJet_eta ),"b5_CorrT1METJet_eta[b5_nCorrT1METJet]/F");
	_compiledTree->Branch("b5_CorrT1METJet_muonSubtrFactor",	( b5.CorrT1METJet_muonSubtrFactor ),"b5_CorrT1METJet_muonSubtrFactor[b5_nCorrT1METJet]/F");
	_compiledTree->Branch("b5_CorrT1METJet_phi",	( b5.CorrT1METJet_phi ),"b5_CorrT1METJet_phi[b5_nCorrT1METJet]/F");
	_compiledTree->Branch("b5_CorrT1METJet_rawPt",	( b5.CorrT1METJet_rawPt ),"b5_CorrT1METJet_rawPt[b5_nCorrT1METJet]/F");
	_compiledTree->Branch("b5_DeepMETResolutionTune_phi",&	( b5.DeepMETResolutionTune_phi ));
	_compiledTree->Branch("b5_DeepMETResolutionTune_pt",&	( b5.DeepMETResolutionTune_pt ));
	_compiledTree->Branch("b5_DeepMETResponseTune_phi",&	( b5.DeepMETResponseTune_phi ));
	_compiledTree->Branch("b5_DeepMETResponseTune_pt",&	( b5.DeepMETResponseTune_pt ));
	_compiledTree->Branch("b5_nElectron",&	( b5.nElectron ));
	_compiledTree->Branch("b5_Electron_deltaEtaSC",	( b5.Electron_deltaEtaSC ),"b5_Electron_deltaEtaSC[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dr03EcalRecHitSumEt",	( b5.Electron_dr03EcalRecHitSumEt ),"b5_Electron_dr03EcalRecHitSumEt[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dr03HcalDepth1TowerSumEt",	( b5.Electron_dr03HcalDepth1TowerSumEt ),"b5_Electron_dr03HcalDepth1TowerSumEt[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dr03TkSumPt",	( b5.Electron_dr03TkSumPt ),"b5_Electron_dr03TkSumPt[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dr03TkSumPtHEEP",	( b5.Electron_dr03TkSumPtHEEP ),"b5_Electron_dr03TkSumPtHEEP[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dxy",	( b5.Electron_dxy ),"b5_Electron_dxy[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dxyErr",	( b5.Electron_dxyErr ),"b5_Electron_dxyErr[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dz",	( b5.Electron_dz ),"b5_Electron_dz[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_dzErr",	( b5.Electron_dzErr ),"b5_Electron_dzErr[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_eCorr",	( b5.Electron_eCorr ),"b5_Electron_eCorr[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_eInvMinusPInv",	( b5.Electron_eInvMinusPInv ),"b5_Electron_eInvMinusPInv[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_energyErr",	( b5.Electron_energyErr ),"b5_Electron_energyErr[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_eta",	( b5.Electron_eta ),"b5_Electron_eta[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_hoe",	( b5.Electron_hoe ),"b5_Electron_hoe[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_ip3d",	( b5.Electron_ip3d ),"b5_Electron_ip3d[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_jetPtRelv2",	( b5.Electron_jetPtRelv2 ),"b5_Electron_jetPtRelv2[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_jetRelIso",	( b5.Electron_jetRelIso ),"b5_Electron_jetRelIso[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_mass",	( b5.Electron_mass ),"b5_Electron_mass[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_miniPFRelIso_all",	( b5.Electron_miniPFRelIso_all ),"b5_Electron_miniPFRelIso_all[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_miniPFRelIso_chg",	( b5.Electron_miniPFRelIso_chg ),"b5_Electron_miniPFRelIso_chg[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_mvaFall17V1Iso",	( b5.Electron_mvaFall17V1Iso ),"b5_Electron_mvaFall17V1Iso[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_mvaFall17V1noIso",	( b5.Electron_mvaFall17V1noIso ),"b5_Electron_mvaFall17V1noIso[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_mvaFall17V2Iso",	( b5.Electron_mvaFall17V2Iso ),"b5_Electron_mvaFall17V2Iso[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_mvaFall17V2noIso",	( b5.Electron_mvaFall17V2noIso ),"b5_Electron_mvaFall17V2noIso[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_pfRelIso03_all",	( b5.Electron_pfRelIso03_all ),"b5_Electron_pfRelIso03_all[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_pfRelIso03_chg",	( b5.Electron_pfRelIso03_chg ),"b5_Electron_pfRelIso03_chg[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_phi",	( b5.Electron_phi ),"b5_Electron_phi[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_pt",	( b5.Electron_pt ),"b5_Electron_pt[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_r9",	( b5.Electron_r9 ),"b5_Electron_r9[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_scEtOverPt",	( b5.Electron_scEtOverPt ),"b5_Electron_scEtOverPt[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_sieie",	( b5.Electron_sieie ),"b5_Electron_sieie[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_sip3d",	( b5.Electron_sip3d ),"b5_Electron_sip3d[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_mvaTTH",	( b5.Electron_mvaTTH ),"b5_Electron_mvaTTH[b5_nElectron]/F");
	_compiledTree->Branch("b5_Electron_charge",	( b5.Electron_charge ),"b5_Electron_charge[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_cutBased",	( b5.Electron_cutBased ),"b5_Electron_cutBased[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_cutBased_Fall17_V1",	( b5.Electron_cutBased_Fall17_V1 ),"b5_Electron_cutBased_Fall17_V1[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_jetIdx",	( b5.Electron_jetIdx ),"b5_Electron_jetIdx[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_pdgId",	( b5.Electron_pdgId ),"b5_Electron_pdgId[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_photonIdx",	( b5.Electron_photonIdx ),"b5_Electron_photonIdx[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_tightCharge",	( b5.Electron_tightCharge ),"b5_Electron_tightCharge[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_vidNestedWPBitmap",	( b5.Electron_vidNestedWPBitmap ),"b5_Electron_vidNestedWPBitmap[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_vidNestedWPBitmapHEEP",	( b5.Electron_vidNestedWPBitmapHEEP ),"b5_Electron_vidNestedWPBitmapHEEP[b5_nElectron]/I");
	_compiledTree->Branch("b5_Electron_convVeto",	( b5.Electron_convVeto ),"b5_Electron_convVeto[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_cutBased_HEEP",	( b5.Electron_cutBased_HEEP ),"b5_Electron_cutBased_HEEP[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_isPFcand",	( b5.Electron_isPFcand ),"b5_Electron_isPFcand[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_jetNDauCharged",	( b5.Electron_jetNDauCharged ),"b5_Electron_jetNDauCharged[b5_nElectron]/b");
	_compiledTree->Branch("b5_Electron_lostHits",	( b5.Electron_lostHits ),"b5_Electron_lostHits[b5_nElectron]/b");
	_compiledTree->Branch("b5_Electron_mvaFall17V1Iso_WP80",	( b5.Electron_mvaFall17V1Iso_WP80 ),"b5_Electron_mvaFall17V1Iso_WP80[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V1Iso_WP90",	( b5.Electron_mvaFall17V1Iso_WP90 ),"b5_Electron_mvaFall17V1Iso_WP90[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V1Iso_WPL",	( b5.Electron_mvaFall17V1Iso_WPL ),"b5_Electron_mvaFall17V1Iso_WPL[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V1noIso_WP80",	( b5.Electron_mvaFall17V1noIso_WP80 ),"b5_Electron_mvaFall17V1noIso_WP80[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V1noIso_WP90",	( b5.Electron_mvaFall17V1noIso_WP90 ),"b5_Electron_mvaFall17V1noIso_WP90[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V1noIso_WPL",	( b5.Electron_mvaFall17V1noIso_WPL ),"b5_Electron_mvaFall17V1noIso_WPL[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V2Iso_WP80",	( b5.Electron_mvaFall17V2Iso_WP80 ),"b5_Electron_mvaFall17V2Iso_WP80[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V2Iso_WP90",	( b5.Electron_mvaFall17V2Iso_WP90 ),"b5_Electron_mvaFall17V2Iso_WP90[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V2Iso_WPL",	( b5.Electron_mvaFall17V2Iso_WPL ),"b5_Electron_mvaFall17V2Iso_WPL[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V2noIso_WP80",	( b5.Electron_mvaFall17V2noIso_WP80 ),"b5_Electron_mvaFall17V2noIso_WP80[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V2noIso_WP90",	( b5.Electron_mvaFall17V2noIso_WP90 ),"b5_Electron_mvaFall17V2noIso_WP90[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_mvaFall17V2noIso_WPL",	( b5.Electron_mvaFall17V2noIso_WPL ),"b5_Electron_mvaFall17V2noIso_WPL[b5_nElectron]/O");
	_compiledTree->Branch("b5_Electron_seedGain",	( b5.Electron_seedGain ),"b5_Electron_seedGain[b5_nElectron]/b");
	_compiledTree->Branch("b5_nFatJet",&	( b5.nFatJet ));
	_compiledTree->Branch("b5_FatJet_area",	( b5.FatJet_area ),"b5_FatJet_area[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagCMVA",	( b5.FatJet_btagCMVA ),"b5_FatJet_btagCMVA[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagCSVV2",	( b5.FatJet_btagCSVV2 ),"b5_FatJet_btagCSVV2[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDBvL",	( b5.FatJet_btagDDBvL ),"b5_FatJet_btagDDBvL[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDBvLV2",	( b5.FatJet_btagDDBvLV2 ),"b5_FatJet_btagDDBvLV2[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDBvL_noMD",	( b5.FatJet_btagDDBvL_noMD ),"b5_FatJet_btagDDBvL_noMD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDCvB",	( b5.FatJet_btagDDCvB ),"b5_FatJet_btagDDCvB[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDCvBV2",	( b5.FatJet_btagDDCvBV2 ),"b5_FatJet_btagDDCvBV2[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDCvB_noMD",	( b5.FatJet_btagDDCvB_noMD ),"b5_FatJet_btagDDCvB_noMD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDCvL",	( b5.FatJet_btagDDCvL ),"b5_FatJet_btagDDCvL[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDCvLV2",	( b5.FatJet_btagDDCvLV2 ),"b5_FatJet_btagDDCvLV2[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDDCvL_noMD",	( b5.FatJet_btagDDCvL_noMD ),"b5_FatJet_btagDDCvL_noMD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagDeepB",	( b5.FatJet_btagDeepB ),"b5_FatJet_btagDeepB[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_btagHbb",	( b5.FatJet_btagHbb ),"b5_FatJet_btagHbb[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_H4qvsQCD",	( b5.FatJet_deepTagMD_H4qvsQCD ),"b5_FatJet_deepTagMD_H4qvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_HbbvsQCD",	( b5.FatJet_deepTagMD_HbbvsQCD ),"b5_FatJet_deepTagMD_HbbvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_TvsQCD",	( b5.FatJet_deepTagMD_TvsQCD ),"b5_FatJet_deepTagMD_TvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_WvsQCD",	( b5.FatJet_deepTagMD_WvsQCD ),"b5_FatJet_deepTagMD_WvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_ZHbbvsQCD",	( b5.FatJet_deepTagMD_ZHbbvsQCD ),"b5_FatJet_deepTagMD_ZHbbvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_ZHccvsQCD",	( b5.FatJet_deepTagMD_ZHccvsQCD ),"b5_FatJet_deepTagMD_ZHccvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_ZbbvsQCD",	( b5.FatJet_deepTagMD_ZbbvsQCD ),"b5_FatJet_deepTagMD_ZbbvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_ZvsQCD",	( b5.FatJet_deepTagMD_ZvsQCD ),"b5_FatJet_deepTagMD_ZvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_bbvsLight",	( b5.FatJet_deepTagMD_bbvsLight ),"b5_FatJet_deepTagMD_bbvsLight[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTagMD_ccvsLight",	( b5.FatJet_deepTagMD_ccvsLight ),"b5_FatJet_deepTagMD_ccvsLight[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTag_H",	( b5.FatJet_deepTag_H ),"b5_FatJet_deepTag_H[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTag_QCD",	( b5.FatJet_deepTag_QCD ),"b5_FatJet_deepTag_QCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTag_QCDothers",	( b5.FatJet_deepTag_QCDothers ),"b5_FatJet_deepTag_QCDothers[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTag_TvsQCD",	( b5.FatJet_deepTag_TvsQCD ),"b5_FatJet_deepTag_TvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTag_WvsQCD",	( b5.FatJet_deepTag_WvsQCD ),"b5_FatJet_deepTag_WvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_deepTag_ZvsQCD",	( b5.FatJet_deepTag_ZvsQCD ),"b5_FatJet_deepTag_ZvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_eta",	( b5.FatJet_eta ),"b5_FatJet_eta[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_mass",	( b5.FatJet_mass ),"b5_FatJet_mass[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_msoftdrop",	( b5.FatJet_msoftdrop ),"b5_FatJet_msoftdrop[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_n2b1",	( b5.FatJet_n2b1 ),"b5_FatJet_n2b1[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_n3b1",	( b5.FatJet_n3b1 ),"b5_FatJet_n3b1[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNetMD_QCD",	( b5.FatJet_particleNetMD_QCD ),"b5_FatJet_particleNetMD_QCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNetMD_Xbb",	( b5.FatJet_particleNetMD_Xbb ),"b5_FatJet_particleNetMD_Xbb[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNetMD_Xcc",	( b5.FatJet_particleNetMD_Xcc ),"b5_FatJet_particleNetMD_Xcc[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNetMD_Xqq",	( b5.FatJet_particleNetMD_Xqq ),"b5_FatJet_particleNetMD_Xqq[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNet_H4qvsQCD",	( b5.FatJet_particleNet_H4qvsQCD ),"b5_FatJet_particleNet_H4qvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNet_HbbvsQCD",	( b5.FatJet_particleNet_HbbvsQCD ),"b5_FatJet_particleNet_HbbvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNet_HccvsQCD",	( b5.FatJet_particleNet_HccvsQCD ),"b5_FatJet_particleNet_HccvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNet_QCD",	( b5.FatJet_particleNet_QCD ),"b5_FatJet_particleNet_QCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNet_TvsQCD",	( b5.FatJet_particleNet_TvsQCD ),"b5_FatJet_particleNet_TvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNet_WvsQCD",	( b5.FatJet_particleNet_WvsQCD ),"b5_FatJet_particleNet_WvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_particleNet_ZvsQCD",	( b5.FatJet_particleNet_ZvsQCD ),"b5_FatJet_particleNet_ZvsQCD[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_phi",	( b5.FatJet_phi ),"b5_FatJet_phi[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_pt",	( b5.FatJet_pt ),"b5_FatJet_pt[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_rawFactor",	( b5.FatJet_rawFactor ),"b5_FatJet_rawFactor[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_tau1",	( b5.FatJet_tau1 ),"b5_FatJet_tau1[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_tau2",	( b5.FatJet_tau2 ),"b5_FatJet_tau2[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_tau3",	( b5.FatJet_tau3 ),"b5_FatJet_tau3[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_tau4",	( b5.FatJet_tau4 ),"b5_FatJet_tau4[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_lsf3",	( b5.FatJet_lsf3 ),"b5_FatJet_lsf3[b5_nFatJet]/F");
	_compiledTree->Branch("b5_FatJet_jetId",	( b5.FatJet_jetId ),"b5_FatJet_jetId[b5_nFatJet]/I");
	_compiledTree->Branch("b5_FatJet_subJetIdx1",	( b5.FatJet_subJetIdx1 ),"b5_FatJet_subJetIdx1[b5_nFatJet]/I");
	_compiledTree->Branch("b5_FatJet_subJetIdx2",	( b5.FatJet_subJetIdx2 ),"b5_FatJet_subJetIdx2[b5_nFatJet]/I");
	_compiledTree->Branch("b5_FatJet_electronIdx3SJ",	( b5.FatJet_electronIdx3SJ ),"b5_FatJet_electronIdx3SJ[b5_nFatJet]/I");
	_compiledTree->Branch("b5_FatJet_muonIdx3SJ",	( b5.FatJet_muonIdx3SJ ),"b5_FatJet_muonIdx3SJ[b5_nFatJet]/I");
	_compiledTree->Branch("b5_nFsrPhoton",&	( b5.nFsrPhoton ));
	_compiledTree->Branch("b5_FsrPhoton_dROverEt2",	( b5.FsrPhoton_dROverEt2 ),"b5_FsrPhoton_dROverEt2[b5_nFsrPhoton]/F");
	_compiledTree->Branch("b5_FsrPhoton_eta",	( b5.FsrPhoton_eta ),"b5_FsrPhoton_eta[b5_nFsrPhoton]/F");
	_compiledTree->Branch("b5_FsrPhoton_phi",	( b5.FsrPhoton_phi ),"b5_FsrPhoton_phi[b5_nFsrPhoton]/F");
	_compiledTree->Branch("b5_FsrPhoton_pt",	( b5.FsrPhoton_pt ),"b5_FsrPhoton_pt[b5_nFsrPhoton]/F");
	_compiledTree->Branch("b5_FsrPhoton_relIso03",	( b5.FsrPhoton_relIso03 ),"b5_FsrPhoton_relIso03[b5_nFsrPhoton]/F");
	_compiledTree->Branch("b5_FsrPhoton_muonIdx",	( b5.FsrPhoton_muonIdx ),"b5_FsrPhoton_muonIdx[b5_nFsrPhoton]/I");
	_compiledTree->Branch("b5_nIsoTrack",&	( b5.nIsoTrack ));
	_compiledTree->Branch("b5_IsoTrack_dxy",	( b5.IsoTrack_dxy ),"b5_IsoTrack_dxy[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_dz",	( b5.IsoTrack_dz ),"b5_IsoTrack_dz[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_eta",	( b5.IsoTrack_eta ),"b5_IsoTrack_eta[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_pfRelIso03_all",	( b5.IsoTrack_pfRelIso03_all ),"b5_IsoTrack_pfRelIso03_all[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_pfRelIso03_chg",	( b5.IsoTrack_pfRelIso03_chg ),"b5_IsoTrack_pfRelIso03_chg[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_phi",	( b5.IsoTrack_phi ),"b5_IsoTrack_phi[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_pt",	( b5.IsoTrack_pt ),"b5_IsoTrack_pt[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_miniPFRelIso_all",	( b5.IsoTrack_miniPFRelIso_all ),"b5_IsoTrack_miniPFRelIso_all[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_miniPFRelIso_chg",	( b5.IsoTrack_miniPFRelIso_chg ),"b5_IsoTrack_miniPFRelIso_chg[b5_nIsoTrack]/F");
	_compiledTree->Branch("b5_IsoTrack_fromPV",	( b5.IsoTrack_fromPV ),"b5_IsoTrack_fromPV[b5_nIsoTrack]/I");
	_compiledTree->Branch("b5_IsoTrack_pdgId",	( b5.IsoTrack_pdgId ),"b5_IsoTrack_pdgId[b5_nIsoTrack]/I");
	_compiledTree->Branch("b5_IsoTrack_isHighPurityTrack",	( b5.IsoTrack_isHighPurityTrack ),"b5_IsoTrack_isHighPurityTrack[b5_nIsoTrack]/O");
	_compiledTree->Branch("b5_IsoTrack_isPFcand",	( b5.IsoTrack_isPFcand ),"b5_IsoTrack_isPFcand[b5_nIsoTrack]/O");
	_compiledTree->Branch("b5_IsoTrack_isFromLostTrack",	( b5.IsoTrack_isFromLostTrack ),"b5_IsoTrack_isFromLostTrack[b5_nIsoTrack]/O");
	_compiledTree->Branch("b5_nJet",&	( b5.nJet ));
	_compiledTree->Branch("b5_Jet_area",	( b5.Jet_area ),"b5_Jet_area[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagCMVA",	( b5.Jet_btagCMVA ),"b5_Jet_btagCMVA[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagCSVV2",	( b5.Jet_btagCSVV2 ),"b5_Jet_btagCSVV2[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepB",	( b5.Jet_btagDeepB ),"b5_Jet_btagDeepB[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepC",	( b5.Jet_btagDeepC ),"b5_Jet_btagDeepC[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepCvB",	( b5.Jet_btagDeepCvB ),"b5_Jet_btagDeepCvB[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepCvL",	( b5.Jet_btagDeepCvL ),"b5_Jet_btagDeepCvL[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepFlavB",	( b5.Jet_btagDeepFlavB ),"b5_Jet_btagDeepFlavB[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepFlavC",	( b5.Jet_btagDeepFlavC ),"b5_Jet_btagDeepFlavC[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepFlavCvB",	( b5.Jet_btagDeepFlavCvB ),"b5_Jet_btagDeepFlavCvB[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepFlavCvL",	( b5.Jet_btagDeepFlavCvL ),"b5_Jet_btagDeepFlavCvL[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_btagDeepFlavQG",	( b5.Jet_btagDeepFlavQG ),"b5_Jet_btagDeepFlavQG[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_chEmEF",	( b5.Jet_chEmEF ),"b5_Jet_chEmEF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_chFPV0EF",	( b5.Jet_chFPV0EF ),"b5_Jet_chFPV0EF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_chFPV1EF",	( b5.Jet_chFPV1EF ),"b5_Jet_chFPV1EF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_chFPV2EF",	( b5.Jet_chFPV2EF ),"b5_Jet_chFPV2EF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_chFPV3EF",	( b5.Jet_chFPV3EF ),"b5_Jet_chFPV3EF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_chHEF",	( b5.Jet_chHEF ),"b5_Jet_chHEF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_eta",	( b5.Jet_eta ),"b5_Jet_eta[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_hfsigmaEtaEta",	( b5.Jet_hfsigmaEtaEta ),"b5_Jet_hfsigmaEtaEta[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_hfsigmaPhiPhi",	( b5.Jet_hfsigmaPhiPhi ),"b5_Jet_hfsigmaPhiPhi[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_mass",	( b5.Jet_mass ),"b5_Jet_mass[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_muEF",	( b5.Jet_muEF ),"b5_Jet_muEF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_muonSubtrFactor",	( b5.Jet_muonSubtrFactor ),"b5_Jet_muonSubtrFactor[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_neEmEF",	( b5.Jet_neEmEF ),"b5_Jet_neEmEF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_neHEF",	( b5.Jet_neHEF ),"b5_Jet_neHEF[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_phi",	( b5.Jet_phi ),"b5_Jet_phi[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_pt",	( b5.Jet_pt ),"b5_Jet_pt[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_puIdDisc",	( b5.Jet_puIdDisc ),"b5_Jet_puIdDisc[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_qgl",	( b5.Jet_qgl ),"b5_Jet_qgl[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_rawFactor",	( b5.Jet_rawFactor ),"b5_Jet_rawFactor[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_bRegCorr",	( b5.Jet_bRegCorr ),"b5_Jet_bRegCorr[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_bRegRes",	( b5.Jet_bRegRes ),"b5_Jet_bRegRes[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_cRegCorr",	( b5.Jet_cRegCorr ),"b5_Jet_cRegCorr[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_cRegRes",	( b5.Jet_cRegRes ),"b5_Jet_cRegRes[b5_nJet]/F");
	_compiledTree->Branch("b5_Jet_electronIdx1",	( b5.Jet_electronIdx1 ),"b5_Jet_electronIdx1[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_electronIdx2",	( b5.Jet_electronIdx2 ),"b5_Jet_electronIdx2[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_hfadjacentEtaStripsSize",	( b5.Jet_hfadjacentEtaStripsSize ),"b5_Jet_hfadjacentEtaStripsSize[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_hfcentralEtaStripSize",	( b5.Jet_hfcentralEtaStripSize ),"b5_Jet_hfcentralEtaStripSize[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_jetId",	( b5.Jet_jetId ),"b5_Jet_jetId[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_muonIdx1",	( b5.Jet_muonIdx1 ),"b5_Jet_muonIdx1[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_muonIdx2",	( b5.Jet_muonIdx2 ),"b5_Jet_muonIdx2[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_nElectrons",	( b5.Jet_nElectrons ),"b5_Jet_nElectrons[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_nMuons",	( b5.Jet_nMuons ),"b5_Jet_nMuons[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_puId",	( b5.Jet_puId ),"b5_Jet_puId[b5_nJet]/I");
	_compiledTree->Branch("b5_Jet_nConstituents",	( b5.Jet_nConstituents ),"b5_Jet_nConstituents[b5_nJet]/b");
	_compiledTree->Branch("b5_L1PreFiringWeight_Dn",&	( b5.L1PreFiringWeight_Dn ));
	_compiledTree->Branch("b5_L1PreFiringWeight_Nom",&	( b5.L1PreFiringWeight_Nom ));
	_compiledTree->Branch("b5_L1PreFiringWeight_Up",&	( b5.L1PreFiringWeight_Up ));
	_compiledTree->Branch("b5_MET_MetUnclustEnUpDeltaX",&	( b5.MET_MetUnclustEnUpDeltaX ));
	_compiledTree->Branch("b5_MET_MetUnclustEnUpDeltaY",&	( b5.MET_MetUnclustEnUpDeltaY ));
	_compiledTree->Branch("b5_MET_covXX",&	( b5.MET_covXX ));
	_compiledTree->Branch("b5_MET_covXY",&	( b5.MET_covXY ));
	_compiledTree->Branch("b5_MET_covYY",&	( b5.MET_covYY ));
	_compiledTree->Branch("b5_MET_phi",&	( b5.MET_phi ));
	_compiledTree->Branch("b5_MET_pt",&	( b5.MET_pt ));
	_compiledTree->Branch("b5_MET_significance",&	( b5.MET_significance ));
	_compiledTree->Branch("b5_MET_sumEt",&	( b5.MET_sumEt ));
	_compiledTree->Branch("b5_MET_sumPtUnclustered",&	( b5.MET_sumPtUnclustered ));
	_compiledTree->Branch("b5_nMuon",&	( b5.nMuon ));
	_compiledTree->Branch("b5_Muon_dxy",	( b5.Muon_dxy ),"b5_Muon_dxy[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_dxyErr",	( b5.Muon_dxyErr ),"b5_Muon_dxyErr[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_dxybs",	( b5.Muon_dxybs ),"b5_Muon_dxybs[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_dz",	( b5.Muon_dz ),"b5_Muon_dz[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_dzErr",	( b5.Muon_dzErr ),"b5_Muon_dzErr[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_eta",	( b5.Muon_eta ),"b5_Muon_eta[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_ip3d",	( b5.Muon_ip3d ),"b5_Muon_ip3d[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_jetPtRelv2",	( b5.Muon_jetPtRelv2 ),"b5_Muon_jetPtRelv2[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_jetRelIso",	( b5.Muon_jetRelIso ),"b5_Muon_jetRelIso[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_mass",	( b5.Muon_mass ),"b5_Muon_mass[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_miniPFRelIso_all",	( b5.Muon_miniPFRelIso_all ),"b5_Muon_miniPFRelIso_all[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_miniPFRelIso_chg",	( b5.Muon_miniPFRelIso_chg ),"b5_Muon_miniPFRelIso_chg[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_pfRelIso03_all",	( b5.Muon_pfRelIso03_all ),"b5_Muon_pfRelIso03_all[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_pfRelIso03_chg",	( b5.Muon_pfRelIso03_chg ),"b5_Muon_pfRelIso03_chg[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_pfRelIso04_all",	( b5.Muon_pfRelIso04_all ),"b5_Muon_pfRelIso04_all[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_phi",	( b5.Muon_phi ),"b5_Muon_phi[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_pt",	( b5.Muon_pt ),"b5_Muon_pt[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_ptErr",	( b5.Muon_ptErr ),"b5_Muon_ptErr[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_segmentComp",	( b5.Muon_segmentComp ),"b5_Muon_segmentComp[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_sip3d",	( b5.Muon_sip3d ),"b5_Muon_sip3d[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_softMva",	( b5.Muon_softMva ),"b5_Muon_softMva[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_tkRelIso",	( b5.Muon_tkRelIso ),"b5_Muon_tkRelIso[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_tunepRelPt",	( b5.Muon_tunepRelPt ),"b5_Muon_tunepRelPt[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_mvaLowPt",	( b5.Muon_mvaLowPt ),"b5_Muon_mvaLowPt[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_mvaTTH",	( b5.Muon_mvaTTH ),"b5_Muon_mvaTTH[b5_nMuon]/F");
	_compiledTree->Branch("b5_Muon_charge",	( b5.Muon_charge ),"b5_Muon_charge[b5_nMuon]/I");
	_compiledTree->Branch("b5_Muon_jetIdx",	( b5.Muon_jetIdx ),"b5_Muon_jetIdx[b5_nMuon]/I");
	_compiledTree->Branch("b5_Muon_nStations",	( b5.Muon_nStations ),"b5_Muon_nStations[b5_nMuon]/I");
	_compiledTree->Branch("b5_Muon_nTrackerLayers",	( b5.Muon_nTrackerLayers ),"b5_Muon_nTrackerLayers[b5_nMuon]/I");
	_compiledTree->Branch("b5_Muon_pdgId",	( b5.Muon_pdgId ),"b5_Muon_pdgId[b5_nMuon]/I");
	_compiledTree->Branch("b5_Muon_tightCharge",	( b5.Muon_tightCharge ),"b5_Muon_tightCharge[b5_nMuon]/I");
	_compiledTree->Branch("b5_Muon_fsrPhotonIdx",	( b5.Muon_fsrPhotonIdx ),"b5_Muon_fsrPhotonIdx[b5_nMuon]/I");
	_compiledTree->Branch("b5_Muon_highPtId",	( b5.Muon_highPtId ),"b5_Muon_highPtId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_highPurity",	( b5.Muon_highPurity ),"b5_Muon_highPurity[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_inTimeMuon",	( b5.Muon_inTimeMuon ),"b5_Muon_inTimeMuon[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_isGlobal",	( b5.Muon_isGlobal ),"b5_Muon_isGlobal[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_isPFcand",	( b5.Muon_isPFcand ),"b5_Muon_isPFcand[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_isTracker",	( b5.Muon_isTracker ),"b5_Muon_isTracker[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_jetNDauCharged",	( b5.Muon_jetNDauCharged ),"b5_Muon_jetNDauCharged[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_looseId",	( b5.Muon_looseId ),"b5_Muon_looseId[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_mediumId",	( b5.Muon_mediumId ),"b5_Muon_mediumId[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_mediumPromptId",	( b5.Muon_mediumPromptId ),"b5_Muon_mediumPromptId[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_miniIsoId",	( b5.Muon_miniIsoId ),"b5_Muon_miniIsoId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_multiIsoId",	( b5.Muon_multiIsoId ),"b5_Muon_multiIsoId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_mvaId",	( b5.Muon_mvaId ),"b5_Muon_mvaId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_mvaLowPtId",	( b5.Muon_mvaLowPtId ),"b5_Muon_mvaLowPtId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_pfIsoId",	( b5.Muon_pfIsoId ),"b5_Muon_pfIsoId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_puppiIsoId",	( b5.Muon_puppiIsoId ),"b5_Muon_puppiIsoId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_softId",	( b5.Muon_softId ),"b5_Muon_softId[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_softMvaId",	( b5.Muon_softMvaId ),"b5_Muon_softMvaId[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_tightId",	( b5.Muon_tightId ),"b5_Muon_tightId[b5_nMuon]/O");
	_compiledTree->Branch("b5_Muon_tkIsoId",	( b5.Muon_tkIsoId ),"b5_Muon_tkIsoId[b5_nMuon]/b");
	_compiledTree->Branch("b5_Muon_triggerIdLoose",	( b5.Muon_triggerIdLoose ),"b5_Muon_triggerIdLoose[b5_nMuon]/O");
	_compiledTree->Branch("b5_nPhoton",&	( b5.nPhoton ));
	_compiledTree->Branch("b5_Photon_eCorr",	( b5.Photon_eCorr ),"b5_Photon_eCorr[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_energyErr",	( b5.Photon_energyErr ),"b5_Photon_energyErr[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_eta",	( b5.Photon_eta ),"b5_Photon_eta[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_hoe",	( b5.Photon_hoe ),"b5_Photon_hoe[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_mass",	( b5.Photon_mass ),"b5_Photon_mass[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_mvaID",	( b5.Photon_mvaID ),"b5_Photon_mvaID[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_mvaID_Fall17V1p1",	( b5.Photon_mvaID_Fall17V1p1 ),"b5_Photon_mvaID_Fall17V1p1[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_pfRelIso03_all",	( b5.Photon_pfRelIso03_all ),"b5_Photon_pfRelIso03_all[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_pfRelIso03_chg",	( b5.Photon_pfRelIso03_chg ),"b5_Photon_pfRelIso03_chg[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_phi",	( b5.Photon_phi ),"b5_Photon_phi[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_pt",	( b5.Photon_pt ),"b5_Photon_pt[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_r9",	( b5.Photon_r9 ),"b5_Photon_r9[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_sieie",	( b5.Photon_sieie ),"b5_Photon_sieie[b5_nPhoton]/F");
	_compiledTree->Branch("b5_Photon_charge",	( b5.Photon_charge ),"b5_Photon_charge[b5_nPhoton]/I");
	_compiledTree->Branch("b5_Photon_cutBased",	( b5.Photon_cutBased ),"b5_Photon_cutBased[b5_nPhoton]/I");
	_compiledTree->Branch("b5_Photon_cutBased_Fall17V1Bitmap",	( b5.Photon_cutBased_Fall17V1Bitmap ),"b5_Photon_cutBased_Fall17V1Bitmap[b5_nPhoton]/I");
	_compiledTree->Branch("b5_Photon_electronIdx",	( b5.Photon_electronIdx ),"b5_Photon_electronIdx[b5_nPhoton]/I");
	_compiledTree->Branch("b5_Photon_jetIdx",	( b5.Photon_jetIdx ),"b5_Photon_jetIdx[b5_nPhoton]/I");
	_compiledTree->Branch("b5_Photon_pdgId",	( b5.Photon_pdgId ),"b5_Photon_pdgId[b5_nPhoton]/I");
	_compiledTree->Branch("b5_Photon_vidNestedWPBitmap",	( b5.Photon_vidNestedWPBitmap ),"b5_Photon_vidNestedWPBitmap[b5_nPhoton]/I");
	_compiledTree->Branch("b5_Photon_electronVeto",	( b5.Photon_electronVeto ),"b5_Photon_electronVeto[b5_nPhoton]/O");
	_compiledTree->Branch("b5_Photon_isScEtaEB",	( b5.Photon_isScEtaEB ),"b5_Photon_isScEtaEB[b5_nPhoton]/O");
	_compiledTree->Branch("b5_Photon_isScEtaEE",	( b5.Photon_isScEtaEE ),"b5_Photon_isScEtaEE[b5_nPhoton]/O");
	_compiledTree->Branch("b5_Photon_mvaID_WP80",	( b5.Photon_mvaID_WP80 ),"b5_Photon_mvaID_WP80[b5_nPhoton]/O");
	_compiledTree->Branch("b5_Photon_mvaID_WP90",	( b5.Photon_mvaID_WP90 ),"b5_Photon_mvaID_WP90[b5_nPhoton]/O");
	_compiledTree->Branch("b5_Photon_pixelSeed",	( b5.Photon_pixelSeed ),"b5_Photon_pixelSeed[b5_nPhoton]/O");
	_compiledTree->Branch("b5_Photon_seedGain",	( b5.Photon_seedGain ),"b5_Photon_seedGain[b5_nPhoton]/b");
	_compiledTree->Branch("b5_PuppiMET_phi",&	( b5.PuppiMET_phi ));
	_compiledTree->Branch("b5_PuppiMET_phiJERDown",&	( b5.PuppiMET_phiJERDown ));
	_compiledTree->Branch("b5_PuppiMET_phiJERUp",&	( b5.PuppiMET_phiJERUp ));
	_compiledTree->Branch("b5_PuppiMET_phiJESDown",&	( b5.PuppiMET_phiJESDown ));
	_compiledTree->Branch("b5_PuppiMET_phiJESUp",&	( b5.PuppiMET_phiJESUp ));
	_compiledTree->Branch("b5_PuppiMET_phiUnclusteredDown",&	( b5.PuppiMET_phiUnclusteredDown ));
	_compiledTree->Branch("b5_PuppiMET_phiUnclusteredUp",&	( b5.PuppiMET_phiUnclusteredUp ));
	_compiledTree->Branch("b5_PuppiMET_pt",&	( b5.PuppiMET_pt ));
	_compiledTree->Branch("b5_PuppiMET_ptJERDown",&	( b5.PuppiMET_ptJERDown ));
	_compiledTree->Branch("b5_PuppiMET_ptJERUp",&	( b5.PuppiMET_ptJERUp ));
	_compiledTree->Branch("b5_PuppiMET_ptJESDown",&	( b5.PuppiMET_ptJESDown ));
	_compiledTree->Branch("b5_PuppiMET_ptJESUp",&	( b5.PuppiMET_ptJESUp ));
	_compiledTree->Branch("b5_PuppiMET_ptUnclusteredDown",&	( b5.PuppiMET_ptUnclusteredDown ));
	_compiledTree->Branch("b5_PuppiMET_ptUnclusteredUp",&	( b5.PuppiMET_ptUnclusteredUp ));
	_compiledTree->Branch("b5_PuppiMET_sumEt",&	( b5.PuppiMET_sumEt ));
	_compiledTree->Branch("b5_RawMET_phi",&	( b5.RawMET_phi ));
	_compiledTree->Branch("b5_RawMET_pt",&	( b5.RawMET_pt ));
	_compiledTree->Branch("b5_RawMET_sumEt",&	( b5.RawMET_sumEt ));
	_compiledTree->Branch("b5_RawPuppiMET_phi",&	( b5.RawPuppiMET_phi ));
	_compiledTree->Branch("b5_RawPuppiMET_pt",&	( b5.RawPuppiMET_pt ));
	_compiledTree->Branch("b5_RawPuppiMET_sumEt",&	( b5.RawPuppiMET_sumEt ));
	_compiledTree->Branch("b5_fixedGridRhoFastjetAll",&	( b5.fixedGridRhoFastjetAll ));
	_compiledTree->Branch("b5_fixedGridRhoFastjetCentral",&	( b5.fixedGridRhoFastjetCentral ));
	_compiledTree->Branch("b5_fixedGridRhoFastjetCentralCalo",&	( b5.fixedGridRhoFastjetCentralCalo ));
	_compiledTree->Branch("b5_fixedGridRhoFastjetCentralChargedPileUp",&	( b5.fixedGridRhoFastjetCentralChargedPileUp ));
	_compiledTree->Branch("b5_fixedGridRhoFastjetCentralNeutral",&	( b5.fixedGridRhoFastjetCentralNeutral ));
	_compiledTree->Branch("b5_nSoftActivityJet",&	( b5.nSoftActivityJet ));
	_compiledTree->Branch("b5_SoftActivityJet_eta",	( b5.SoftActivityJet_eta ),"b5_SoftActivityJet_eta[b5_nSoftActivityJet]/F");
	_compiledTree->Branch("b5_SoftActivityJet_phi",	( b5.SoftActivityJet_phi ),"b5_SoftActivityJet_phi[b5_nSoftActivityJet]/F");
	_compiledTree->Branch("b5_SoftActivityJet_pt",	( b5.SoftActivityJet_pt ),"b5_SoftActivityJet_pt[b5_nSoftActivityJet]/F");
	_compiledTree->Branch("b5_SoftActivityJetHT",&	( b5.SoftActivityJetHT ));
	_compiledTree->Branch("b5_SoftActivityJetHT10",&	( b5.SoftActivityJetHT10 ));
	_compiledTree->Branch("b5_SoftActivityJetHT2",&	( b5.SoftActivityJetHT2 ));
	_compiledTree->Branch("b5_SoftActivityJetHT5",&	( b5.SoftActivityJetHT5 ));
	_compiledTree->Branch("b5_SoftActivityJetNjets10",&	( b5.SoftActivityJetNjets10 ));
	_compiledTree->Branch("b5_SoftActivityJetNjets2",&	( b5.SoftActivityJetNjets2 ));
	_compiledTree->Branch("b5_SoftActivityJetNjets5",&	( b5.SoftActivityJetNjets5 ));
	_compiledTree->Branch("b5_nSubJet",&	( b5.nSubJet ));
	_compiledTree->Branch("b5_SubJet_btagCMVA",	( b5.SubJet_btagCMVA ),"b5_SubJet_btagCMVA[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_btagCSVV2",	( b5.SubJet_btagCSVV2 ),"b5_SubJet_btagCSVV2[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_btagDeepB",	( b5.SubJet_btagDeepB ),"b5_SubJet_btagDeepB[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_eta",	( b5.SubJet_eta ),"b5_SubJet_eta[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_mass",	( b5.SubJet_mass ),"b5_SubJet_mass[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_n2b1",	( b5.SubJet_n2b1 ),"b5_SubJet_n2b1[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_n3b1",	( b5.SubJet_n3b1 ),"b5_SubJet_n3b1[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_phi",	( b5.SubJet_phi ),"b5_SubJet_phi[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_pt",	( b5.SubJet_pt ),"b5_SubJet_pt[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_rawFactor",	( b5.SubJet_rawFactor ),"b5_SubJet_rawFactor[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_tau1",	( b5.SubJet_tau1 ),"b5_SubJet_tau1[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_tau2",	( b5.SubJet_tau2 ),"b5_SubJet_tau2[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_tau3",	( b5.SubJet_tau3 ),"b5_SubJet_tau3[b5_nSubJet]/F");
	_compiledTree->Branch("b5_SubJet_tau4",	( b5.SubJet_tau4 ),"b5_SubJet_tau4[b5_nSubJet]/F");
	_compiledTree->Branch("b5_nTau",&	( b5.nTau ));
	_compiledTree->Branch("b5_Tau_chargedIso",	( b5.Tau_chargedIso ),"b5_Tau_chargedIso[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_dxy",	( b5.Tau_dxy ),"b5_Tau_dxy[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_dz",	( b5.Tau_dz ),"b5_Tau_dz[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_eta",	( b5.Tau_eta ),"b5_Tau_eta[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_leadTkDeltaEta",	( b5.Tau_leadTkDeltaEta ),"b5_Tau_leadTkDeltaEta[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_leadTkDeltaPhi",	( b5.Tau_leadTkDeltaPhi ),"b5_Tau_leadTkDeltaPhi[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_leadTkPtOverTauPt",	( b5.Tau_leadTkPtOverTauPt ),"b5_Tau_leadTkPtOverTauPt[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_mass",	( b5.Tau_mass ),"b5_Tau_mass[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_neutralIso",	( b5.Tau_neutralIso ),"b5_Tau_neutralIso[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_phi",	( b5.Tau_phi ),"b5_Tau_phi[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_photonsOutsideSignalCone",	( b5.Tau_photonsOutsideSignalCone ),"b5_Tau_photonsOutsideSignalCone[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_pt",	( b5.Tau_pt ),"b5_Tau_pt[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_puCorr",	( b5.Tau_puCorr ),"b5_Tau_puCorr[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawAntiEle",	( b5.Tau_rawAntiEle ),"b5_Tau_rawAntiEle[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawAntiEle2018",	( b5.Tau_rawAntiEle2018 ),"b5_Tau_rawAntiEle2018[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawDeepTau2017v2p1VSe",	( b5.Tau_rawDeepTau2017v2p1VSe ),"b5_Tau_rawDeepTau2017v2p1VSe[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawDeepTau2017v2p1VSjet",	( b5.Tau_rawDeepTau2017v2p1VSjet ),"b5_Tau_rawDeepTau2017v2p1VSjet[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawDeepTau2017v2p1VSmu",	( b5.Tau_rawDeepTau2017v2p1VSmu ),"b5_Tau_rawDeepTau2017v2p1VSmu[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawIso",	( b5.Tau_rawIso ),"b5_Tau_rawIso[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawIsodR03",	( b5.Tau_rawIsodR03 ),"b5_Tau_rawIsodR03[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawMVAnewDM2017v2",	( b5.Tau_rawMVAnewDM2017v2 ),"b5_Tau_rawMVAnewDM2017v2[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawMVAoldDM",	( b5.Tau_rawMVAoldDM ),"b5_Tau_rawMVAoldDM[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawMVAoldDM2017v1",	( b5.Tau_rawMVAoldDM2017v1 ),"b5_Tau_rawMVAoldDM2017v1[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawMVAoldDM2017v2",	( b5.Tau_rawMVAoldDM2017v2 ),"b5_Tau_rawMVAoldDM2017v2[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_rawMVAoldDMdR032017v2",	( b5.Tau_rawMVAoldDMdR032017v2 ),"b5_Tau_rawMVAoldDMdR032017v2[b5_nTau]/F");
	_compiledTree->Branch("b5_Tau_charge",	( b5.Tau_charge ),"b5_Tau_charge[b5_nTau]/I");
	_compiledTree->Branch("b5_Tau_decayMode",	( b5.Tau_decayMode ),"b5_Tau_decayMode[b5_nTau]/I");
	_compiledTree->Branch("b5_Tau_jetIdx",	( b5.Tau_jetIdx ),"b5_Tau_jetIdx[b5_nTau]/I");
	_compiledTree->Branch("b5_Tau_rawAntiEleCat",	( b5.Tau_rawAntiEleCat ),"b5_Tau_rawAntiEleCat[b5_nTau]/I");
	_compiledTree->Branch("b5_Tau_rawAntiEleCat2018",	( b5.Tau_rawAntiEleCat2018 ),"b5_Tau_rawAntiEleCat2018[b5_nTau]/I");
	_compiledTree->Branch("b5_Tau_idAntiEle",	( b5.Tau_idAntiEle ),"b5_Tau_idAntiEle[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idAntiEle2018",	( b5.Tau_idAntiEle2018 ),"b5_Tau_idAntiEle2018[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idAntiEleDeadECal",	( b5.Tau_idAntiEleDeadECal ),"b5_Tau_idAntiEleDeadECal[b5_nTau]/O");
	_compiledTree->Branch("b5_Tau_idAntiMu",	( b5.Tau_idAntiMu ),"b5_Tau_idAntiMu[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idDecayMode",	( b5.Tau_idDecayMode ),"b5_Tau_idDecayMode[b5_nTau]/O");
	_compiledTree->Branch("b5_Tau_idDecayModeNewDMs",	( b5.Tau_idDecayModeNewDMs ),"b5_Tau_idDecayModeNewDMs[b5_nTau]/O");
	_compiledTree->Branch("b5_Tau_idDeepTau2017v2p1VSe",	( b5.Tau_idDeepTau2017v2p1VSe ),"b5_Tau_idDeepTau2017v2p1VSe[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idDeepTau2017v2p1VSjet",	( b5.Tau_idDeepTau2017v2p1VSjet ),"b5_Tau_idDeepTau2017v2p1VSjet[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idDeepTau2017v2p1VSmu",	( b5.Tau_idDeepTau2017v2p1VSmu ),"b5_Tau_idDeepTau2017v2p1VSmu[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idMVAnewDM2017v2",	( b5.Tau_idMVAnewDM2017v2 ),"b5_Tau_idMVAnewDM2017v2[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idMVAoldDM",	( b5.Tau_idMVAoldDM ),"b5_Tau_idMVAoldDM[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idMVAoldDM2017v1",	( b5.Tau_idMVAoldDM2017v1 ),"b5_Tau_idMVAoldDM2017v1[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idMVAoldDM2017v2",	( b5.Tau_idMVAoldDM2017v2 ),"b5_Tau_idMVAoldDM2017v2[b5_nTau]/b");
	_compiledTree->Branch("b5_Tau_idMVAoldDMdR032017v2",	( b5.Tau_idMVAoldDMdR032017v2 ),"b5_Tau_idMVAoldDMdR032017v2[b5_nTau]/b");
	_compiledTree->Branch("b5_TkMET_phi",&	( b5.TkMET_phi ));
	_compiledTree->Branch("b5_TkMET_pt",&	( b5.TkMET_pt ));
	_compiledTree->Branch("b5_TkMET_sumEt",&	( b5.TkMET_sumEt ));
	_compiledTree->Branch("b5_nTrigObj",&	( b5.nTrigObj ));
	_compiledTree->Branch("b5_TrigObj_pt",	( b5.TrigObj_pt ),"b5_TrigObj_pt[b5_nTrigObj]/F");
	_compiledTree->Branch("b5_TrigObj_eta",	( b5.TrigObj_eta ),"b5_TrigObj_eta[b5_nTrigObj]/F");
	_compiledTree->Branch("b5_TrigObj_phi",	( b5.TrigObj_phi ),"b5_TrigObj_phi[b5_nTrigObj]/F");
	_compiledTree->Branch("b5_TrigObj_l1pt",	( b5.TrigObj_l1pt ),"b5_TrigObj_l1pt[b5_nTrigObj]/F");
	_compiledTree->Branch("b5_TrigObj_l1pt_2",	( b5.TrigObj_l1pt_2 ),"b5_TrigObj_l1pt_2[b5_nTrigObj]/F");
	_compiledTree->Branch("b5_TrigObj_l2pt",	( b5.TrigObj_l2pt ),"b5_TrigObj_l2pt[b5_nTrigObj]/F");
	_compiledTree->Branch("b5_TrigObj_id",	( b5.TrigObj_id ),"b5_TrigObj_id[b5_nTrigObj]/I");
	_compiledTree->Branch("b5_TrigObj_l1iso",	( b5.TrigObj_l1iso ),"b5_TrigObj_l1iso[b5_nTrigObj]/I");
	_compiledTree->Branch("b5_TrigObj_l1charge",	( b5.TrigObj_l1charge ),"b5_TrigObj_l1charge[b5_nTrigObj]/I");
	_compiledTree->Branch("b5_TrigObj_filterBits",	( b5.TrigObj_filterBits ),"b5_TrigObj_filterBits[b5_nTrigObj]/I");
	_compiledTree->Branch("b5_nOtherPV",&	( b5.nOtherPV ));
	_compiledTree->Branch("b5_OtherPV_z",	( b5.OtherPV_z ),"b5_OtherPV_z[b5_nOtherPV]/F");
	_compiledTree->Branch("b5_PV_ndof",&	( b5.PV_ndof ));
	_compiledTree->Branch("b5_PV_x",&	( b5.PV_x ));
	_compiledTree->Branch("b5_PV_y",&	( b5.PV_y ));
	_compiledTree->Branch("b5_PV_z",&	( b5.PV_z ));
	_compiledTree->Branch("b5_PV_chi2",&	( b5.PV_chi2 ));
	_compiledTree->Branch("b5_PV_score",&	( b5.PV_score ));
	_compiledTree->Branch("b5_PV_npvs",&	( b5.PV_npvs ));
	_compiledTree->Branch("b5_PV_npvsGood",&	( b5.PV_npvsGood ));
	_compiledTree->Branch("b5_nSV",&	( b5.nSV ));
	_compiledTree->Branch("b5_SV_dlen",	( b5.SV_dlen ),"b5_SV_dlen[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_dlenSig",	( b5.SV_dlenSig ),"b5_SV_dlenSig[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_dxy",	( b5.SV_dxy ),"b5_SV_dxy[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_dxySig",	( b5.SV_dxySig ),"b5_SV_dxySig[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_pAngle",	( b5.SV_pAngle ),"b5_SV_pAngle[b5_nSV]/F");
	_compiledTree->Branch("b5_Electron_cleanmask",	( b5.Electron_cleanmask ),"b5_Electron_cleanmask[b5_nElectron]/b");
	_compiledTree->Branch("b5_Jet_cleanmask",	( b5.Jet_cleanmask ),"b5_Jet_cleanmask[b5_nJet]/b");
	_compiledTree->Branch("b5_Muon_cleanmask",	( b5.Muon_cleanmask ),"b5_Muon_cleanmask[b5_nMuon]/b");
	_compiledTree->Branch("b5_Photon_cleanmask",	( b5.Photon_cleanmask ),"b5_Photon_cleanmask[b5_nPhoton]/b");
	_compiledTree->Branch("b5_Tau_cleanmask",	( b5.Tau_cleanmask ),"b5_Tau_cleanmask[b5_nTau]/b");
	_compiledTree->Branch("b5_SV_chi2",	( b5.SV_chi2 ),"b5_SV_chi2[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_eta",	( b5.SV_eta ),"b5_SV_eta[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_mass",	( b5.SV_mass ),"b5_SV_mass[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_ndof",	( b5.SV_ndof ),"b5_SV_ndof[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_phi",	( b5.SV_phi ),"b5_SV_phi[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_pt",	( b5.SV_pt ),"b5_SV_pt[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_x",	( b5.SV_x ),"b5_SV_x[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_y",	( b5.SV_y ),"b5_SV_y[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_z",	( b5.SV_z ),"b5_SV_z[b5_nSV]/F");
	_compiledTree->Branch("b5_SV_ntracks",	( b5.SV_ntracks ),"b5_SV_ntracks[b5_nSV]/b");
	_compiledTree->Branch("b5_L1_AlwaysTrue",&	( b5.L1_AlwaysTrue ));
	_compiledTree->Branch("b5_L1_BPTX_AND_Ref1_VME",&	( b5.L1_BPTX_AND_Ref1_VME ));
	_compiledTree->Branch("b5_L1_BPTX_AND_Ref3_VME",&	( b5.L1_BPTX_AND_Ref3_VME ));
	_compiledTree->Branch("b5_L1_BPTX_AND_Ref4_VME",&	( b5.L1_BPTX_AND_Ref4_VME ));
	_compiledTree->Branch("b5_L1_BPTX_BeamGas_B1_VME",&	( b5.L1_BPTX_BeamGas_B1_VME ));
	_compiledTree->Branch("b5_L1_BPTX_BeamGas_B2_VME",&	( b5.L1_BPTX_BeamGas_B2_VME ));
	_compiledTree->Branch("b5_L1_BPTX_BeamGas_Ref1_VME",&	( b5.L1_BPTX_BeamGas_Ref1_VME ));
	_compiledTree->Branch("b5_L1_BPTX_BeamGas_Ref2_VME",&	( b5.L1_BPTX_BeamGas_Ref2_VME ));
	_compiledTree->Branch("b5_L1_BPTX_NotOR_VME",&	( b5.L1_BPTX_NotOR_VME ));
	_compiledTree->Branch("b5_L1_BPTX_OR_Ref3_VME",&	( b5.L1_BPTX_OR_Ref3_VME ));
	_compiledTree->Branch("b5_L1_BPTX_OR_Ref4_VME",&	( b5.L1_BPTX_OR_Ref4_VME ));
	_compiledTree->Branch("b5_L1_BPTX_RefAND_VME",&	( b5.L1_BPTX_RefAND_VME ));
	_compiledTree->Branch("b5_L1_BptxMinus",&	( b5.L1_BptxMinus ));
	_compiledTree->Branch("b5_L1_BptxOR",&	( b5.L1_BptxOR ));
	_compiledTree->Branch("b5_L1_BptxPlus",&	( b5.L1_BptxPlus ));
	_compiledTree->Branch("b5_L1_BptxXOR",&	( b5.L1_BptxXOR ));
	_compiledTree->Branch("b5_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142",&	( b5.L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 ));
	_compiledTree->Branch("b5_L1_DoubleEG8er2p5_HTT260er",&	( b5.L1_DoubleEG8er2p5_HTT260er ));
	_compiledTree->Branch("b5_L1_DoubleEG8er2p5_HTT280er",&	( b5.L1_DoubleEG8er2p5_HTT280er ));
	_compiledTree->Branch("b5_L1_DoubleEG8er2p5_HTT300er",&	( b5.L1_DoubleEG8er2p5_HTT300er ));
	_compiledTree->Branch("b5_L1_DoubleEG8er2p5_HTT320er",&	( b5.L1_DoubleEG8er2p5_HTT320er ));
	_compiledTree->Branch("b5_L1_DoubleEG8er2p5_HTT340er",&	( b5.L1_DoubleEG8er2p5_HTT340er ));
	_compiledTree->Branch("b5_L1_DoubleEG_15_10_er2p5",&	( b5.L1_DoubleEG_15_10_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_20_10_er2p5",&	( b5.L1_DoubleEG_20_10_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_22_10_er2p5",&	( b5.L1_DoubleEG_22_10_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_25_12_er2p5",&	( b5.L1_DoubleEG_25_12_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_25_14_er2p5",&	( b5.L1_DoubleEG_25_14_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_27_14_er2p5",&	( b5.L1_DoubleEG_27_14_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_LooseIso20_10_er2p5",&	( b5.L1_DoubleEG_LooseIso20_10_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_LooseIso22_10_er2p5",&	( b5.L1_DoubleEG_LooseIso22_10_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_LooseIso22_12_er2p5",&	( b5.L1_DoubleEG_LooseIso22_12_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleEG_LooseIso25_12_er2p5",&	( b5.L1_DoubleEG_LooseIso25_12_er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleIsoTau32er2p1",&	( b5.L1_DoubleIsoTau32er2p1 ));
	_compiledTree->Branch("b5_L1_DoubleIsoTau34er2p1",&	( b5.L1_DoubleIsoTau34er2p1 ));
	_compiledTree->Branch("b5_L1_DoubleIsoTau36er2p1",&	( b5.L1_DoubleIsoTau36er2p1 ));
	_compiledTree->Branch("b5_L1_DoubleJet100er2p3_dEta_Max1p6",&	( b5.L1_DoubleJet100er2p3_dEta_Max1p6 ));
	_compiledTree->Branch("b5_L1_DoubleJet100er2p5",&	( b5.L1_DoubleJet100er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet112er2p3_dEta_Max1p6",&	( b5.L1_DoubleJet112er2p3_dEta_Max1p6 ));
	_compiledTree->Branch("b5_L1_DoubleJet120er2p5",&	( b5.L1_DoubleJet120er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet150er2p5",&	( b5.L1_DoubleJet150er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5",&	( b5.L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5",&	( b5.L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5",&	( b5.L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5",&	( b5.L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5",&	( b5.L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5",&	( b5.L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp",&	( b5.L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp ));
	_compiledTree->Branch("b5_L1_DoubleJet40er2p5",&	( b5.L1_DoubleJet40er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620",&	( b5.L1_DoubleJet_100_30_DoubleJet30_Mass_Min620 ));
	_compiledTree->Branch("b5_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620",&	( b5.L1_DoubleJet_110_35_DoubleJet35_Mass_Min620 ));
	_compiledTree->Branch("b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620",&	( b5.L1_DoubleJet_115_40_DoubleJet40_Mass_Min620 ));
	_compiledTree->Branch("b5_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28",&	( b5.L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28 ));
	_compiledTree->Branch("b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620",&	( b5.L1_DoubleJet_120_45_DoubleJet45_Mass_Min620 ));
	_compiledTree->Branch("b5_L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28",&	( b5.L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28 ));
	_compiledTree->Branch("b5_L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ",&	( b5.L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ ));
	_compiledTree->Branch("b5_L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp",&	( b5.L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp ));
	_compiledTree->Branch("b5_L1_DoubleJet_80_30_Mass_Min420_Mu8",&	( b5.L1_DoubleJet_80_30_Mass_Min420_Mu8 ));
	_compiledTree->Branch("b5_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620",&	( b5.L1_DoubleJet_90_30_DoubleJet30_Mass_Min620 ));
	_compiledTree->Branch("b5_L1_DoubleLooseIsoEG22er2p1",&	( b5.L1_DoubleLooseIsoEG22er2p1 ));
	_compiledTree->Branch("b5_L1_DoubleLooseIsoEG24er2p1",&	( b5.L1_DoubleLooseIsoEG24er2p1 ));
	_compiledTree->Branch("b5_L1_DoubleMu0",&	( b5.L1_DoubleMu0 ));
	_compiledTree->Branch("b5_L1_DoubleMu0_Mass_Min1",&	( b5.L1_DoubleMu0_Mass_Min1 ));
	_compiledTree->Branch("b5_L1_DoubleMu0_OQ",&	( b5.L1_DoubleMu0_OQ ));
	_compiledTree->Branch("b5_L1_DoubleMu0_SQ",&	( b5.L1_DoubleMu0_SQ ));
	_compiledTree->Branch("b5_L1_DoubleMu0_SQ_OS",&	( b5.L1_DoubleMu0_SQ_OS ));
	_compiledTree->Branch("b5_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8",&	( b5.L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8 ));
	_compiledTree->Branch("b5_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",&	( b5.L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 ));
	_compiledTree->Branch("b5_L1_DoubleMu0er1p5_SQ",&	( b5.L1_DoubleMu0er1p5_SQ ));
	_compiledTree->Branch("b5_L1_DoubleMu0er1p5_SQ_OS",&	( b5.L1_DoubleMu0er1p5_SQ_OS ));
	_compiledTree->Branch("b5_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",&	( b5.L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 ));
	_compiledTree->Branch("b5_L1_DoubleMu0er1p5_SQ_dR_Max1p4",&	( b5.L1_DoubleMu0er1p5_SQ_dR_Max1p4 ));
	_compiledTree->Branch("b5_L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4",&	( b5.L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4 ));
	_compiledTree->Branch("b5_L1_DoubleMu0er2p0_SQ_dR_Max1p4",&	( b5.L1_DoubleMu0er2p0_SQ_dR_Max1p4 ));
	_compiledTree->Branch("b5_L1_DoubleMu10_SQ",&	( b5.L1_DoubleMu10_SQ ));
	_compiledTree->Branch("b5_L1_DoubleMu18er2p1",&	( b5.L1_DoubleMu18er2p1 ));
	_compiledTree->Branch("b5_L1_DoubleMu3_OS_DoubleEG7p5Upsilon",&	( b5.L1_DoubleMu3_OS_DoubleEG7p5Upsilon ));
	_compiledTree->Branch("b5_L1_DoubleMu3_SQ_ETMHF50_HTT60er",&	( b5.L1_DoubleMu3_SQ_ETMHF50_HTT60er ));
	_compiledTree->Branch("b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5",&	( b5.L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5",&	( b5.L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5",&	( b5.L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleMu3_SQ_HTT220er",&	( b5.L1_DoubleMu3_SQ_HTT220er ));
	_compiledTree->Branch("b5_L1_DoubleMu3_SQ_HTT240er",&	( b5.L1_DoubleMu3_SQ_HTT240er ));
	_compiledTree->Branch("b5_L1_DoubleMu3_SQ_HTT260er",&	( b5.L1_DoubleMu3_SQ_HTT260er ));
	_compiledTree->Branch("b5_L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8",&	( b5.L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8 ));
	_compiledTree->Branch("b5_L1_DoubleMu4_SQ_EG9er2p5",&	( b5.L1_DoubleMu4_SQ_EG9er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleMu4_SQ_OS",&	( b5.L1_DoubleMu4_SQ_OS ));
	_compiledTree->Branch("b5_L1_DoubleMu4_SQ_OS_dR_Max1p2",&	( b5.L1_DoubleMu4_SQ_OS_dR_Max1p2 ));
	_compiledTree->Branch("b5_L1_DoubleMu4p5_SQ_OS",&	( b5.L1_DoubleMu4p5_SQ_OS ));
	_compiledTree->Branch("b5_L1_DoubleMu4p5_SQ_OS_dR_Max1p2",&	( b5.L1_DoubleMu4p5_SQ_OS_dR_Max1p2 ));
	_compiledTree->Branch("b5_L1_DoubleMu4p5er2p0_SQ_OS",&	( b5.L1_DoubleMu4p5er2p0_SQ_OS ));
	_compiledTree->Branch("b5_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18",&	( b5.L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18 ));
	_compiledTree->Branch("b5_L1_DoubleMu5Upsilon_OS_DoubleEG3",&	( b5.L1_DoubleMu5Upsilon_OS_DoubleEG3 ));
	_compiledTree->Branch("b5_L1_DoubleMu5_SQ_EG9er2p5",&	( b5.L1_DoubleMu5_SQ_EG9er2p5 ));
	_compiledTree->Branch("b5_L1_DoubleMu9_SQ",&	( b5.L1_DoubleMu9_SQ ));
	_compiledTree->Branch("b5_L1_DoubleMu_12_5",&	( b5.L1_DoubleMu_12_5 ));
	_compiledTree->Branch("b5_L1_DoubleMu_15_5_SQ",&	( b5.L1_DoubleMu_15_5_SQ ));
	_compiledTree->Branch("b5_L1_DoubleMu_15_7",&	( b5.L1_DoubleMu_15_7 ));
	_compiledTree->Branch("b5_L1_DoubleMu_15_7_Mass_Min1",&	( b5.L1_DoubleMu_15_7_Mass_Min1 ));
	_compiledTree->Branch("b5_L1_DoubleMu_15_7_SQ",&	( b5.L1_DoubleMu_15_7_SQ ));
	_compiledTree->Branch("b5_L1_DoubleTau70er2p1",&	( b5.L1_DoubleTau70er2p1 ));
	_compiledTree->Branch("b5_L1_ETM120",&	( b5.L1_ETM120 ));
	_compiledTree->Branch("b5_L1_ETM150",&	( b5.L1_ETM150 ));
	_compiledTree->Branch("b5_L1_ETMHF100",&	( b5.L1_ETMHF100 ));
	_compiledTree->Branch("b5_L1_ETMHF100_HTT60er",&	( b5.L1_ETMHF100_HTT60er ));
	_compiledTree->Branch("b5_L1_ETMHF110",&	( b5.L1_ETMHF110 ));
	_compiledTree->Branch("b5_L1_ETMHF110_HTT60er",&	( b5.L1_ETMHF110_HTT60er ));
	_compiledTree->Branch("b5_L1_ETMHF110_HTT60er_NotSecondBunchInTrain",&	( b5.L1_ETMHF110_HTT60er_NotSecondBunchInTrain ));
	_compiledTree->Branch("b5_L1_ETMHF120",&	( b5.L1_ETMHF120 ));
	_compiledTree->Branch("b5_L1_ETMHF120_HTT60er",&	( b5.L1_ETMHF120_HTT60er ));
	_compiledTree->Branch("b5_L1_ETMHF120_NotSecondBunchInTrain",&	( b5.L1_ETMHF120_NotSecondBunchInTrain ));
	_compiledTree->Branch("b5_L1_ETMHF130",&	( b5.L1_ETMHF130 ));
	_compiledTree->Branch("b5_L1_ETMHF130_HTT60er",&	( b5.L1_ETMHF130_HTT60er ));
	_compiledTree->Branch("b5_L1_ETMHF140",&	( b5.L1_ETMHF140 ));
	_compiledTree->Branch("b5_L1_ETMHF150",&	( b5.L1_ETMHF150 ));
	_compiledTree->Branch("b5_L1_ETMHF90_HTT60er",&	( b5.L1_ETMHF90_HTT60er ));
	_compiledTree->Branch("b5_L1_ETT1200",&	( b5.L1_ETT1200 ));
	_compiledTree->Branch("b5_L1_ETT1600",&	( b5.L1_ETT1600 ));
	_compiledTree->Branch("b5_L1_ETT2000",&	( b5.L1_ETT2000 ));
	_compiledTree->Branch("b5_L1_FirstBunchAfterTrain",&	( b5.L1_FirstBunchAfterTrain ));
	_compiledTree->Branch("b5_L1_FirstBunchBeforeTrain",&	( b5.L1_FirstBunchBeforeTrain ));
	_compiledTree->Branch("b5_L1_FirstBunchInTrain",&	( b5.L1_FirstBunchInTrain ));
	_compiledTree->Branch("b5_L1_FirstCollisionInOrbit",&	( b5.L1_FirstCollisionInOrbit ));
	_compiledTree->Branch("b5_L1_FirstCollisionInTrain",&	( b5.L1_FirstCollisionInTrain ));
	_compiledTree->Branch("b5_L1_HCAL_LaserMon_Trig",&	( b5.L1_HCAL_LaserMon_Trig ));
	_compiledTree->Branch("b5_L1_HCAL_LaserMon_Veto",&	( b5.L1_HCAL_LaserMon_Veto ));
	_compiledTree->Branch("b5_L1_HTT120er",&	( b5.L1_HTT120er ));
	_compiledTree->Branch("b5_L1_HTT160er",&	( b5.L1_HTT160er ));
	_compiledTree->Branch("b5_L1_HTT200er",&	( b5.L1_HTT200er ));
	_compiledTree->Branch("b5_L1_HTT255er",&	( b5.L1_HTT255er ));
	_compiledTree->Branch("b5_L1_HTT280er",&	( b5.L1_HTT280er ));
	_compiledTree->Branch("b5_L1_HTT280er_QuadJet_70_55_40_35_er2p4",&	( b5.L1_HTT280er_QuadJet_70_55_40_35_er2p4 ));
	_compiledTree->Branch("b5_L1_HTT320er",&	( b5.L1_HTT320er ));
	_compiledTree->Branch("b5_L1_HTT320er_QuadJet_70_55_40_40_er2p4",&	( b5.L1_HTT320er_QuadJet_70_55_40_40_er2p4 ));
	_compiledTree->Branch("b5_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3",&	( b5.L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 ));
	_compiledTree->Branch("b5_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3",&	( b5.L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 ));
	_compiledTree->Branch("b5_L1_HTT360er",&	( b5.L1_HTT360er ));
	_compiledTree->Branch("b5_L1_HTT400er",&	( b5.L1_HTT400er ));
	_compiledTree->Branch("b5_L1_HTT450er",&	( b5.L1_HTT450er ));
	_compiledTree->Branch("b5_L1_IsoEG32er2p5_Mt40",&	( b5.L1_IsoEG32er2p5_Mt40 ));
	_compiledTree->Branch("b5_L1_IsoEG32er2p5_Mt44",&	( b5.L1_IsoEG32er2p5_Mt44 ));
	_compiledTree->Branch("b5_L1_IsoEG32er2p5_Mt48",&	( b5.L1_IsoEG32er2p5_Mt48 ));
	_compiledTree->Branch("b5_L1_IsoTau40er2p1_ETMHF100",&	( b5.L1_IsoTau40er2p1_ETMHF100 ));
	_compiledTree->Branch("b5_L1_IsoTau40er2p1_ETMHF110",&	( b5.L1_IsoTau40er2p1_ETMHF110 ));
	_compiledTree->Branch("b5_L1_IsoTau40er2p1_ETMHF120",&	( b5.L1_IsoTau40er2p1_ETMHF120 ));
	_compiledTree->Branch("b5_L1_IsoTau40er2p1_ETMHF90",&	( b5.L1_IsoTau40er2p1_ETMHF90 ));
	_compiledTree->Branch("b5_L1_IsolatedBunch",&	( b5.L1_IsolatedBunch ));
	_compiledTree->Branch("b5_L1_LastBunchInTrain",&	( b5.L1_LastBunchInTrain ));
	_compiledTree->Branch("b5_L1_LastCollisionInTrain",&	( b5.L1_LastCollisionInTrain ));
	_compiledTree->Branch("b5_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3",&	( b5.L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3 ));
	_compiledTree->Branch("b5_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3",&	( b5.L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3 ));
	_compiledTree->Branch("b5_L1_LooseIsoEG24er2p1_HTT100er",&	( b5.L1_LooseIsoEG24er2p1_HTT100er ));
	_compiledTree->Branch("b5_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3",&	( b5.L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3 ));
	_compiledTree->Branch("b5_L1_LooseIsoEG26er2p1_HTT100er",&	( b5.L1_LooseIsoEG26er2p1_HTT100er ));
	_compiledTree->Branch("b5_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3",&	( b5.L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 ));
	_compiledTree->Branch("b5_L1_LooseIsoEG28er2p1_HTT100er",&	( b5.L1_LooseIsoEG28er2p1_HTT100er ));
	_compiledTree->Branch("b5_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3",&	( b5.L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 ));
	_compiledTree->Branch("b5_L1_LooseIsoEG30er2p1_HTT100er",&	( b5.L1_LooseIsoEG30er2p1_HTT100er ));
	_compiledTree->Branch("b5_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3",&	( b5.L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 ));
	_compiledTree->Branch("b5_L1_MinimumBiasHF0_AND_BptxAND",&	( b5.L1_MinimumBiasHF0_AND_BptxAND ));
	_compiledTree->Branch("b5_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6",&	( b5.L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6 ));
	_compiledTree->Branch("b5_L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6",&	( b5.L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6 ));
	_compiledTree->Branch("b5_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6",&	( b5.L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6 ));
	_compiledTree->Branch("b5_L1_Mu18er2p1_Tau24er2p1",&	( b5.L1_Mu18er2p1_Tau24er2p1 ));
	_compiledTree->Branch("b5_L1_Mu18er2p1_Tau26er2p1",&	( b5.L1_Mu18er2p1_Tau26er2p1 ));
	_compiledTree->Branch("b5_L1_Mu20_EG10er2p5",&	( b5.L1_Mu20_EG10er2p5 ));
	_compiledTree->Branch("b5_L1_Mu22er2p1_IsoTau32er2p1",&	( b5.L1_Mu22er2p1_IsoTau32er2p1 ));
	_compiledTree->Branch("b5_L1_Mu22er2p1_IsoTau34er2p1",&	( b5.L1_Mu22er2p1_IsoTau34er2p1 ));
	_compiledTree->Branch("b5_L1_Mu22er2p1_IsoTau36er2p1",&	( b5.L1_Mu22er2p1_IsoTau36er2p1 ));
	_compiledTree->Branch("b5_L1_Mu22er2p1_IsoTau40er2p1",&	( b5.L1_Mu22er2p1_IsoTau40er2p1 ));
	_compiledTree->Branch("b5_L1_Mu22er2p1_Tau70er2p1",&	( b5.L1_Mu22er2p1_Tau70er2p1 ));
	_compiledTree->Branch("b5_L1_Mu3_Jet120er2p5_dR_Max0p4",&	( b5.L1_Mu3_Jet120er2p5_dR_Max0p4 ));
	_compiledTree->Branch("b5_L1_Mu3_Jet120er2p5_dR_Max0p8",&	( b5.L1_Mu3_Jet120er2p5_dR_Max0p8 ));
	_compiledTree->Branch("b5_L1_Mu3_Jet16er2p5_dR_Max0p4",&	( b5.L1_Mu3_Jet16er2p5_dR_Max0p4 ));
	_compiledTree->Branch("b5_L1_Mu3_Jet30er2p5",&	( b5.L1_Mu3_Jet30er2p5 ));
	_compiledTree->Branch("b5_L1_Mu3_Jet35er2p5_dR_Max0p4",&	( b5.L1_Mu3_Jet35er2p5_dR_Max0p4 ));
	_compiledTree->Branch("b5_L1_Mu3_Jet60er2p5_dR_Max0p4",&	( b5.L1_Mu3_Jet60er2p5_dR_Max0p4 ));
	_compiledTree->Branch("b5_L1_Mu3_Jet80er2p5_dR_Max0p4",&	( b5.L1_Mu3_Jet80er2p5_dR_Max0p4 ));
	_compiledTree->Branch("b5_L1_Mu3er1p5_Jet100er2p5_ETMHF40",&	( b5.L1_Mu3er1p5_Jet100er2p5_ETMHF40 ));
	_compiledTree->Branch("b5_L1_Mu3er1p5_Jet100er2p5_ETMHF50",&	( b5.L1_Mu3er1p5_Jet100er2p5_ETMHF50 ));
	_compiledTree->Branch("b5_L1_Mu5_EG23er2p5",&	( b5.L1_Mu5_EG23er2p5 ));
	_compiledTree->Branch("b5_L1_Mu5_LooseIsoEG20er2p5",&	( b5.L1_Mu5_LooseIsoEG20er2p5 ));
	_compiledTree->Branch("b5_L1_Mu6_DoubleEG10er2p5",&	( b5.L1_Mu6_DoubleEG10er2p5 ));
	_compiledTree->Branch("b5_L1_Mu6_DoubleEG12er2p5",&	( b5.L1_Mu6_DoubleEG12er2p5 ));
	_compiledTree->Branch("b5_L1_Mu6_DoubleEG15er2p5",&	( b5.L1_Mu6_DoubleEG15er2p5 ));
	_compiledTree->Branch("b5_L1_Mu6_DoubleEG17er2p5",&	( b5.L1_Mu6_DoubleEG17er2p5 ));
	_compiledTree->Branch("b5_L1_Mu6_HTT240er",&	( b5.L1_Mu6_HTT240er ));
	_compiledTree->Branch("b5_L1_Mu6_HTT250er",&	( b5.L1_Mu6_HTT250er ));
	_compiledTree->Branch("b5_L1_Mu7_EG23er2p5",&	( b5.L1_Mu7_EG23er2p5 ));
	_compiledTree->Branch("b5_L1_Mu7_LooseIsoEG20er2p5",&	( b5.L1_Mu7_LooseIsoEG20er2p5 ));
	_compiledTree->Branch("b5_L1_Mu7_LooseIsoEG23er2p5",&	( b5.L1_Mu7_LooseIsoEG23er2p5 ));
	_compiledTree->Branch("b5_L1_NotBptxOR",&	( b5.L1_NotBptxOR ));
	_compiledTree->Branch("b5_L1_QuadJet36er2p5_IsoTau52er2p1",&	( b5.L1_QuadJet36er2p5_IsoTau52er2p1 ));
	_compiledTree->Branch("b5_L1_QuadJet60er2p5",&	( b5.L1_QuadJet60er2p5 ));
	_compiledTree->Branch("b5_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0",&	( b5.L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0 ));
	_compiledTree->Branch("b5_L1_QuadMu0",&	( b5.L1_QuadMu0 ));
	_compiledTree->Branch("b5_L1_QuadMu0_OQ",&	( b5.L1_QuadMu0_OQ ));
	_compiledTree->Branch("b5_L1_QuadMu0_SQ",&	( b5.L1_QuadMu0_SQ ));
	_compiledTree->Branch("b5_L1_SecondBunchInTrain",&	( b5.L1_SecondBunchInTrain ));
	_compiledTree->Branch("b5_L1_SecondLastBunchInTrain",&	( b5.L1_SecondLastBunchInTrain ));
	_compiledTree->Branch("b5_L1_SingleEG10er2p5",&	( b5.L1_SingleEG10er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG15er2p5",&	( b5.L1_SingleEG15er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG26er2p5",&	( b5.L1_SingleEG26er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG34er2p5",&	( b5.L1_SingleEG34er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG36er2p5",&	( b5.L1_SingleEG36er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG38er2p5",&	( b5.L1_SingleEG38er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG40er2p5",&	( b5.L1_SingleEG40er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG42er2p5",&	( b5.L1_SingleEG42er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG45er2p5",&	( b5.L1_SingleEG45er2p5 ));
	_compiledTree->Branch("b5_L1_SingleEG50",&	( b5.L1_SingleEG50 ));
	_compiledTree->Branch("b5_L1_SingleEG60",&	( b5.L1_SingleEG60 ));
	_compiledTree->Branch("b5_L1_SingleEG8er2p5",&	( b5.L1_SingleEG8er2p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG24er1p5",&	( b5.L1_SingleIsoEG24er1p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG24er2p1",&	( b5.L1_SingleIsoEG24er2p1 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG26er1p5",&	( b5.L1_SingleIsoEG26er1p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG26er2p1",&	( b5.L1_SingleIsoEG26er2p1 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG26er2p5",&	( b5.L1_SingleIsoEG26er2p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG28er1p5",&	( b5.L1_SingleIsoEG28er1p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG28er2p1",&	( b5.L1_SingleIsoEG28er2p1 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG28er2p5",&	( b5.L1_SingleIsoEG28er2p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG30er2p1",&	( b5.L1_SingleIsoEG30er2p1 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG30er2p5",&	( b5.L1_SingleIsoEG30er2p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG32er2p1",&	( b5.L1_SingleIsoEG32er2p1 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG32er2p5",&	( b5.L1_SingleIsoEG32er2p5 ));
	_compiledTree->Branch("b5_L1_SingleIsoEG34er2p5",&	( b5.L1_SingleIsoEG34er2p5 ));
	_compiledTree->Branch("b5_L1_SingleJet10erHE",&	( b5.L1_SingleJet10erHE ));
	_compiledTree->Branch("b5_L1_SingleJet120",&	( b5.L1_SingleJet120 ));
	_compiledTree->Branch("b5_L1_SingleJet120_FWD3p0",&	( b5.L1_SingleJet120_FWD3p0 ));
	_compiledTree->Branch("b5_L1_SingleJet120er2p5",&	( b5.L1_SingleJet120er2p5 ));
	_compiledTree->Branch("b5_L1_SingleJet12erHE",&	( b5.L1_SingleJet12erHE ));
	_compiledTree->Branch("b5_L1_SingleJet140er2p5",&	( b5.L1_SingleJet140er2p5 ));
	_compiledTree->Branch("b5_L1_SingleJet140er2p5_ETMHF80",&	( b5.L1_SingleJet140er2p5_ETMHF80 ));
	_compiledTree->Branch("b5_L1_SingleJet140er2p5_ETMHF90",&	( b5.L1_SingleJet140er2p5_ETMHF90 ));
	_compiledTree->Branch("b5_L1_SingleJet160er2p5",&	( b5.L1_SingleJet160er2p5 ));
	_compiledTree->Branch("b5_L1_SingleJet180",&	( b5.L1_SingleJet180 ));
	_compiledTree->Branch("b5_L1_SingleJet180er2p5",&	( b5.L1_SingleJet180er2p5 ));
	_compiledTree->Branch("b5_L1_SingleJet200",&	( b5.L1_SingleJet200 ));
	_compiledTree->Branch("b5_L1_SingleJet20er2p5_NotBptxOR",&	( b5.L1_SingleJet20er2p5_NotBptxOR ));
	_compiledTree->Branch("b5_L1_SingleJet20er2p5_NotBptxOR_3BX",&	( b5.L1_SingleJet20er2p5_NotBptxOR_3BX ));
	_compiledTree->Branch("b5_L1_SingleJet35",&	( b5.L1_SingleJet35 ));
	_compiledTree->Branch("b5_L1_SingleJet35_FWD3p0",&	( b5.L1_SingleJet35_FWD3p0 ));
	_compiledTree->Branch("b5_L1_SingleJet35er2p5",&	( b5.L1_SingleJet35er2p5 ));
	_compiledTree->Branch("b5_L1_SingleJet43er2p5_NotBptxOR_3BX",&	( b5.L1_SingleJet43er2p5_NotBptxOR_3BX ));
	_compiledTree->Branch("b5_L1_SingleJet46er2p5_NotBptxOR_3BX",&	( b5.L1_SingleJet46er2p5_NotBptxOR_3BX ));
	_compiledTree->Branch("b5_L1_SingleJet60",&	( b5.L1_SingleJet60 ));
	_compiledTree->Branch("b5_L1_SingleJet60_FWD3p0",&	( b5.L1_SingleJet60_FWD3p0 ));
	_compiledTree->Branch("b5_L1_SingleJet60er2p5",&	( b5.L1_SingleJet60er2p5 ));
	_compiledTree->Branch("b5_L1_SingleJet8erHE",&	( b5.L1_SingleJet8erHE ));
	_compiledTree->Branch("b5_L1_SingleJet90",&	( b5.L1_SingleJet90 ));
	_compiledTree->Branch("b5_L1_SingleJet90_FWD3p0",&	( b5.L1_SingleJet90_FWD3p0 ));
	_compiledTree->Branch("b5_L1_SingleJet90er2p5",&	( b5.L1_SingleJet90er2p5 ));
	_compiledTree->Branch("b5_L1_SingleLooseIsoEG28er1p5",&	( b5.L1_SingleLooseIsoEG28er1p5 ));
	_compiledTree->Branch("b5_L1_SingleLooseIsoEG30er1p5",&	( b5.L1_SingleLooseIsoEG30er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu0_BMTF",&	( b5.L1_SingleMu0_BMTF ));
	_compiledTree->Branch("b5_L1_SingleMu0_DQ",&	( b5.L1_SingleMu0_DQ ));
	_compiledTree->Branch("b5_L1_SingleMu0_EMTF",&	( b5.L1_SingleMu0_EMTF ));
	_compiledTree->Branch("b5_L1_SingleMu0_OMTF",&	( b5.L1_SingleMu0_OMTF ));
	_compiledTree->Branch("b5_L1_SingleMu10er1p5",&	( b5.L1_SingleMu10er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu12_DQ_BMTF",&	( b5.L1_SingleMu12_DQ_BMTF ));
	_compiledTree->Branch("b5_L1_SingleMu12_DQ_EMTF",&	( b5.L1_SingleMu12_DQ_EMTF ));
	_compiledTree->Branch("b5_L1_SingleMu12_DQ_OMTF",&	( b5.L1_SingleMu12_DQ_OMTF ));
	_compiledTree->Branch("b5_L1_SingleMu12er1p5",&	( b5.L1_SingleMu12er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu14er1p5",&	( b5.L1_SingleMu14er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu15_DQ",&	( b5.L1_SingleMu15_DQ ));
	_compiledTree->Branch("b5_L1_SingleMu16er1p5",&	( b5.L1_SingleMu16er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu18",&	( b5.L1_SingleMu18 ));
	_compiledTree->Branch("b5_L1_SingleMu18er1p5",&	( b5.L1_SingleMu18er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu20",&	( b5.L1_SingleMu20 ));
	_compiledTree->Branch("b5_L1_SingleMu22",&	( b5.L1_SingleMu22 ));
	_compiledTree->Branch("b5_L1_SingleMu22_BMTF",&	( b5.L1_SingleMu22_BMTF ));
	_compiledTree->Branch("b5_L1_SingleMu22_EMTF",&	( b5.L1_SingleMu22_EMTF ));
	_compiledTree->Branch("b5_L1_SingleMu22_OMTF",&	( b5.L1_SingleMu22_OMTF ));
	_compiledTree->Branch("b5_L1_SingleMu25",&	( b5.L1_SingleMu25 ));
	_compiledTree->Branch("b5_L1_SingleMu3",&	( b5.L1_SingleMu3 ));
	_compiledTree->Branch("b5_L1_SingleMu5",&	( b5.L1_SingleMu5 ));
	_compiledTree->Branch("b5_L1_SingleMu6er1p5",&	( b5.L1_SingleMu6er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu7",&	( b5.L1_SingleMu7 ));
	_compiledTree->Branch("b5_L1_SingleMu7_DQ",&	( b5.L1_SingleMu7_DQ ));
	_compiledTree->Branch("b5_L1_SingleMu7er1p5",&	( b5.L1_SingleMu7er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu8er1p5",&	( b5.L1_SingleMu8er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMu9er1p5",&	( b5.L1_SingleMu9er1p5 ));
	_compiledTree->Branch("b5_L1_SingleMuCosmics",&	( b5.L1_SingleMuCosmics ));
	_compiledTree->Branch("b5_L1_SingleMuCosmics_BMTF",&	( b5.L1_SingleMuCosmics_BMTF ));
	_compiledTree->Branch("b5_L1_SingleMuCosmics_EMTF",&	( b5.L1_SingleMuCosmics_EMTF ));
	_compiledTree->Branch("b5_L1_SingleMuCosmics_OMTF",&	( b5.L1_SingleMuCosmics_OMTF ));
	_compiledTree->Branch("b5_L1_SingleMuOpen",&	( b5.L1_SingleMuOpen ));
	_compiledTree->Branch("b5_L1_SingleMuOpen_NotBptxOR",&	( b5.L1_SingleMuOpen_NotBptxOR ));
	_compiledTree->Branch("b5_L1_SingleMuOpen_er1p1_NotBptxOR_3BX",&	( b5.L1_SingleMuOpen_er1p1_NotBptxOR_3BX ));
	_compiledTree->Branch("b5_L1_SingleMuOpen_er1p4_NotBptxOR_3BX",&	( b5.L1_SingleMuOpen_er1p4_NotBptxOR_3BX ));
	_compiledTree->Branch("b5_L1_SingleTau120er2p1",&	( b5.L1_SingleTau120er2p1 ));
	_compiledTree->Branch("b5_L1_SingleTau130er2p1",&	( b5.L1_SingleTau130er2p1 ));
	_compiledTree->Branch("b5_L1_TOTEM_1",&	( b5.L1_TOTEM_1 ));
	_compiledTree->Branch("b5_L1_TOTEM_2",&	( b5.L1_TOTEM_2 ));
	_compiledTree->Branch("b5_L1_TOTEM_3",&	( b5.L1_TOTEM_3 ));
	_compiledTree->Branch("b5_L1_TOTEM_4",&	( b5.L1_TOTEM_4 ));
	_compiledTree->Branch("b5_L1_TripleEG16er2p5",&	( b5.L1_TripleEG16er2p5 ));
	_compiledTree->Branch("b5_L1_TripleEG_16_12_8_er2p5",&	( b5.L1_TripleEG_16_12_8_er2p5 ));
	_compiledTree->Branch("b5_L1_TripleEG_16_15_8_er2p5",&	( b5.L1_TripleEG_16_15_8_er2p5 ));
	_compiledTree->Branch("b5_L1_TripleEG_18_17_8_er2p5",&	( b5.L1_TripleEG_18_17_8_er2p5 ));
	_compiledTree->Branch("b5_L1_TripleEG_18_18_12_er2p5",&	( b5.L1_TripleEG_18_18_12_er2p5 ));
	_compiledTree->Branch("b5_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5",&	( b5.L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5 ));
	_compiledTree->Branch("b5_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5",&	( b5.L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5 ));
	_compiledTree->Branch("b5_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5",&	( b5.L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5 ));
	_compiledTree->Branch("b5_L1_TripleMu0",&	( b5.L1_TripleMu0 ));
	_compiledTree->Branch("b5_L1_TripleMu0_OQ",&	( b5.L1_TripleMu0_OQ ));
	_compiledTree->Branch("b5_L1_TripleMu0_SQ",&	( b5.L1_TripleMu0_SQ ));
	_compiledTree->Branch("b5_L1_TripleMu3",&	( b5.L1_TripleMu3 ));
	_compiledTree->Branch("b5_L1_TripleMu3_SQ",&	( b5.L1_TripleMu3_SQ ));
	_compiledTree->Branch("b5_L1_TripleMu_5SQ_3SQ_0OQ",&	( b5.L1_TripleMu_5SQ_3SQ_0OQ ));
	_compiledTree->Branch("b5_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9",&	( b5.L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9 ));
	_compiledTree->Branch("b5_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9",&	( b5.L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9 ));
	_compiledTree->Branch("b5_L1_TripleMu_5_3_3",&	( b5.L1_TripleMu_5_3_3 ));
	_compiledTree->Branch("b5_L1_TripleMu_5_3_3_SQ",&	( b5.L1_TripleMu_5_3_3_SQ ));
	_compiledTree->Branch("b5_L1_TripleMu_5_3p5_2p5",&	( b5.L1_TripleMu_5_3p5_2p5 ));
	_compiledTree->Branch("b5_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17",&	( b5.L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17 ));
	_compiledTree->Branch("b5_L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17",&	( b5.L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17 ));
	_compiledTree->Branch("b5_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17",&	( b5.L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17 ));
	_compiledTree->Branch("b5_L1_TripleMu_5_5_3",&	( b5.L1_TripleMu_5_5_3 ));
	_compiledTree->Branch("b5_L1_UnpairedBunchBptxMinus",&	( b5.L1_UnpairedBunchBptxMinus ));
	_compiledTree->Branch("b5_L1_UnpairedBunchBptxPlus",&	( b5.L1_UnpairedBunchBptxPlus ));
	_compiledTree->Branch("b5_L1_ZeroBias",&	( b5.L1_ZeroBias ));
	_compiledTree->Branch("b5_L1_ZeroBias_copy",&	( b5.L1_ZeroBias_copy ));
	_compiledTree->Branch("b5_L1_UnprefireableEvent",&	( b5.L1_UnprefireableEvent ));
	_compiledTree->Branch("b5_Flag_HBHENoiseFilter",&	( b5.Flag_HBHENoiseFilter ));
	_compiledTree->Branch("b5_Flag_HBHENoiseIsoFilter",&	( b5.Flag_HBHENoiseIsoFilter ));
	_compiledTree->Branch("b5_Flag_CSCTightHaloFilter",&	( b5.Flag_CSCTightHaloFilter ));
	_compiledTree->Branch("b5_Flag_CSCTightHaloTrkMuUnvetoFilter",&	( b5.Flag_CSCTightHaloTrkMuUnvetoFilter ));
	_compiledTree->Branch("b5_Flag_CSCTightHalo2015Filter",&	( b5.Flag_CSCTightHalo2015Filter ));
	_compiledTree->Branch("b5_Flag_globalTightHalo2016Filter",&	( b5.Flag_globalTightHalo2016Filter ));
	_compiledTree->Branch("b5_Flag_globalSuperTightHalo2016Filter",&	( b5.Flag_globalSuperTightHalo2016Filter ));
	_compiledTree->Branch("b5_Flag_HcalStripHaloFilter",&	( b5.Flag_HcalStripHaloFilter ));
	_compiledTree->Branch("b5_Flag_hcalLaserEventFilter",&	( b5.Flag_hcalLaserEventFilter ));
	_compiledTree->Branch("b5_Flag_EcalDeadCellTriggerPrimitiveFilter",&	( b5.Flag_EcalDeadCellTriggerPrimitiveFilter ));
	_compiledTree->Branch("b5_Flag_EcalDeadCellBoundaryEnergyFilter",&	( b5.Flag_EcalDeadCellBoundaryEnergyFilter ));
	_compiledTree->Branch("b5_Flag_ecalBadCalibFilter",&	( b5.Flag_ecalBadCalibFilter ));
	_compiledTree->Branch("b5_Flag_goodVertices",&	( b5.Flag_goodVertices ));
	_compiledTree->Branch("b5_Flag_eeBadScFilter",&	( b5.Flag_eeBadScFilter ));
	_compiledTree->Branch("b5_Flag_ecalLaserCorrFilter",&	( b5.Flag_ecalLaserCorrFilter ));
	_compiledTree->Branch("b5_Flag_trkPOGFilters",&	( b5.Flag_trkPOGFilters ));
	_compiledTree->Branch("b5_Flag_chargedHadronTrackResolutionFilter",&	( b5.Flag_chargedHadronTrackResolutionFilter ));
	_compiledTree->Branch("b5_Flag_muonBadTrackFilter",&	( b5.Flag_muonBadTrackFilter ));
	_compiledTree->Branch("b5_Flag_BadChargedCandidateFilter",&	( b5.Flag_BadChargedCandidateFilter ));
	_compiledTree->Branch("b5_Flag_BadPFMuonFilter",&	( b5.Flag_BadPFMuonFilter ));
	_compiledTree->Branch("b5_Flag_BadPFMuonDzFilter",&	( b5.Flag_BadPFMuonDzFilter ));
	_compiledTree->Branch("b5_Flag_hfNoisyHitsFilter",&	( b5.Flag_hfNoisyHitsFilter ));
	_compiledTree->Branch("b5_Flag_BadChargedCandidateSummer16Filter",&	( b5.Flag_BadChargedCandidateSummer16Filter ));
	_compiledTree->Branch("b5_Flag_BadPFMuonSummer16Filter",&	( b5.Flag_BadPFMuonSummer16Filter ));
	_compiledTree->Branch("b5_Flag_trkPOG_manystripclus53X",&	( b5.Flag_trkPOG_manystripclus53X ));
	_compiledTree->Branch("b5_Flag_trkPOG_toomanystripclus53X",&	( b5.Flag_trkPOG_toomanystripclus53X ));
	_compiledTree->Branch("b5_Flag_trkPOG_logErrorTooManyClusters",&	( b5.Flag_trkPOG_logErrorTooManyClusters ));
	_compiledTree->Branch("b5_Flag_METFilters",&	( b5.Flag_METFilters ));
	_compiledTree->Branch("b5_L1Reco_step",&	( b5.L1Reco_step ));
	_compiledTree->Branch("b5_HLTriggerFirstPath",&	( b5.HLTriggerFirstPath ));
	_compiledTree->Branch("b5_HLT_AK8PFJet360_TrimMass30",&	( b5.HLT_AK8PFJet360_TrimMass30 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet380_TrimMass30",&	( b5.HLT_AK8PFJet380_TrimMass30 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet400_TrimMass30",&	( b5.HLT_AK8PFJet400_TrimMass30 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet420_TrimMass30",&	( b5.HLT_AK8PFJet420_TrimMass30 ));
	_compiledTree->Branch("b5_HLT_AK8PFHT750_TrimMass50",&	( b5.HLT_AK8PFHT750_TrimMass50 ));
	_compiledTree->Branch("b5_HLT_AK8PFHT800_TrimMass50",&	( b5.HLT_AK8PFHT800_TrimMass50 ));
	_compiledTree->Branch("b5_HLT_AK8PFHT850_TrimMass50",&	( b5.HLT_AK8PFHT850_TrimMass50 ));
	_compiledTree->Branch("b5_HLT_AK8PFHT900_TrimMass50",&	( b5.HLT_AK8PFHT900_TrimMass50 ));
	_compiledTree->Branch("b5_HLT_CaloJet500_NoJetID",&	( b5.HLT_CaloJet500_NoJetID ));
	_compiledTree->Branch("b5_HLT_CaloJet550_NoJetID",&	( b5.HLT_CaloJet550_NoJetID ));
	_compiledTree->Branch("b5_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL",&	( b5.HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon",&	( b5.HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon ));
	_compiledTree->Branch("b5_HLT_Trimuon5_3p5_2_Upsilon_Muon",&	( b5.HLT_Trimuon5_3p5_2_Upsilon_Muon ));
	_compiledTree->Branch("b5_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon",&	( b5.HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon ));
	_compiledTree->Branch("b5_HLT_DoubleEle25_CaloIdL_MW",&	( b5.HLT_DoubleEle25_CaloIdL_MW ));
	_compiledTree->Branch("b5_HLT_DoubleEle27_CaloIdL_MW",&	( b5.HLT_DoubleEle27_CaloIdL_MW ));
	_compiledTree->Branch("b5_HLT_DoubleEle33_CaloIdL_MW",&	( b5.HLT_DoubleEle33_CaloIdL_MW ));
	_compiledTree->Branch("b5_HLT_DoubleEle24_eta2p1_WPTight_Gsf",&	( b5.HLT_DoubleEle24_eta2p1_WPTight_Gsf ));
	_compiledTree->Branch("b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350",&	( b5.HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350 ));
	_compiledTree->Branch("b5_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350",&	( b5.HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350 ));
	_compiledTree->Branch("b5_HLT_Ele27_Ele37_CaloIdL_MW",&	( b5.HLT_Ele27_Ele37_CaloIdL_MW ));
	_compiledTree->Branch("b5_HLT_Mu27_Ele37_CaloIdL_MW",&	( b5.HLT_Mu27_Ele37_CaloIdL_MW ));
	_compiledTree->Branch("b5_HLT_Mu37_Ele27_CaloIdL_MW",&	( b5.HLT_Mu37_Ele27_CaloIdL_MW ));
	_compiledTree->Branch("b5_HLT_Mu37_TkMu27",&	( b5.HLT_Mu37_TkMu27 ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_3_Bs",&	( b5.HLT_DoubleMu4_3_Bs ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_3_Jpsi",&	( b5.HLT_DoubleMu4_3_Jpsi ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_JpsiTrk_Displaced",&	( b5.HLT_DoubleMu4_JpsiTrk_Displaced ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",&	( b5.HLT_DoubleMu4_LowMassNonResonantTrk_Displaced ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_Trk_Tau3mu",&	( b5.HLT_DoubleMu3_Trk_Tau3mu ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_TkMu_DsTau3Mu",&	( b5.HLT_DoubleMu3_TkMu_DsTau3Mu ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_PsiPrimeTrk_Displaced",&	( b5.HLT_DoubleMu4_PsiPrimeTrk_Displaced ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_Mass3p8_DZ_PFHT350",&	( b5.HLT_DoubleMu4_Mass3p8_DZ_PFHT350 ));
	_compiledTree->Branch("b5_HLT_Mu3_PFJet40",&	( b5.HLT_Mu3_PFJet40 ));
	_compiledTree->Branch("b5_HLT_Mu7p5_L2Mu2_Jpsi",&	( b5.HLT_Mu7p5_L2Mu2_Jpsi ));
	_compiledTree->Branch("b5_HLT_Mu7p5_L2Mu2_Upsilon",&	( b5.HLT_Mu7p5_L2Mu2_Upsilon ));
	_compiledTree->Branch("b5_HLT_Mu7p5_Track2_Jpsi",&	( b5.HLT_Mu7p5_Track2_Jpsi ));
	_compiledTree->Branch("b5_HLT_Mu7p5_Track3p5_Jpsi",&	( b5.HLT_Mu7p5_Track3p5_Jpsi ));
	_compiledTree->Branch("b5_HLT_Mu7p5_Track7_Jpsi",&	( b5.HLT_Mu7p5_Track7_Jpsi ));
	_compiledTree->Branch("b5_HLT_Mu7p5_Track2_Upsilon",&	( b5.HLT_Mu7p5_Track2_Upsilon ));
	_compiledTree->Branch("b5_HLT_Mu7p5_Track3p5_Upsilon",&	( b5.HLT_Mu7p5_Track3p5_Upsilon ));
	_compiledTree->Branch("b5_HLT_Mu7p5_Track7_Upsilon",&	( b5.HLT_Mu7p5_Track7_Upsilon ));
	_compiledTree->Branch("b5_HLT_DoublePhoton33_CaloIdL",&	( b5.HLT_DoublePhoton33_CaloIdL ));
	_compiledTree->Branch("b5_HLT_DoublePhoton70",&	( b5.HLT_DoublePhoton70 ));
	_compiledTree->Branch("b5_HLT_DoublePhoton85",&	( b5.HLT_DoublePhoton85 ));
	_compiledTree->Branch("b5_HLT_Ele20_WPTight_Gsf",&	( b5.HLT_Ele20_WPTight_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele15_WPLoose_Gsf",&	( b5.HLT_Ele15_WPLoose_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele17_WPLoose_Gsf",&	( b5.HLT_Ele17_WPLoose_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele20_WPLoose_Gsf",&	( b5.HLT_Ele20_WPLoose_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele20_eta2p1_WPLoose_Gsf",&	( b5.HLT_Ele20_eta2p1_WPLoose_Gsf ));
	_compiledTree->Branch("b5_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG",&	( b5.HLT_DiEle27_WPTightCaloOnly_L1DoubleEG ));
	_compiledTree->Branch("b5_HLT_Ele27_WPTight_Gsf",&	( b5.HLT_Ele27_WPTight_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele32_WPTight_Gsf",&	( b5.HLT_Ele32_WPTight_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele35_WPTight_Gsf",&	( b5.HLT_Ele35_WPTight_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele35_WPTight_Gsf_L1EGMT",&	( b5.HLT_Ele35_WPTight_Gsf_L1EGMT ));
	_compiledTree->Branch("b5_HLT_Ele38_WPTight_Gsf",&	( b5.HLT_Ele38_WPTight_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele40_WPTight_Gsf",&	( b5.HLT_Ele40_WPTight_Gsf ));
	_compiledTree->Branch("b5_HLT_Ele32_WPTight_Gsf_L1DoubleEG",&	( b5.HLT_Ele32_WPTight_Gsf_L1DoubleEG ));
	_compiledTree->Branch("b5_HLT_HT450_Beamspot",&	( b5.HLT_HT450_Beamspot ));
	_compiledTree->Branch("b5_HLT_HT300_Beamspot",&	( b5.HLT_HT300_Beamspot ));
	_compiledTree->Branch("b5_HLT_ZeroBias_Beamspot",&	( b5.HLT_ZeroBias_Beamspot ));
	_compiledTree->Branch("b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",&	( b5.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1",&	( b5.HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1",&	( b5.HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1",&	( b5.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1",&	( b5.HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1",&	( b5.HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu20",&	( b5.HLT_IsoMu20 ));
	_compiledTree->Branch("b5_HLT_IsoMu24",&	( b5.HLT_IsoMu24 ));
	_compiledTree->Branch("b5_HLT_IsoMu24_eta2p1",&	( b5.HLT_IsoMu24_eta2p1 ));
	_compiledTree->Branch("b5_HLT_IsoMu27",&	( b5.HLT_IsoMu27 ));
	_compiledTree->Branch("b5_HLT_IsoMu30",&	( b5.HLT_IsoMu30 ));
	_compiledTree->Branch("b5_HLT_UncorrectedJetE30_NoBPTX",&	( b5.HLT_UncorrectedJetE30_NoBPTX ));
	_compiledTree->Branch("b5_HLT_UncorrectedJetE30_NoBPTX3BX",&	( b5.HLT_UncorrectedJetE30_NoBPTX3BX ));
	_compiledTree->Branch("b5_HLT_UncorrectedJetE60_NoBPTX3BX",&	( b5.HLT_UncorrectedJetE60_NoBPTX3BX ));
	_compiledTree->Branch("b5_HLT_UncorrectedJetE70_NoBPTX3BX",&	( b5.HLT_UncorrectedJetE70_NoBPTX3BX ));
	_compiledTree->Branch("b5_HLT_L1SingleMu18",&	( b5.HLT_L1SingleMu18 ));
	_compiledTree->Branch("b5_HLT_L1SingleMu25",&	( b5.HLT_L1SingleMu25 ));
	_compiledTree->Branch("b5_HLT_L2Mu10",&	( b5.HLT_L2Mu10 ));
	_compiledTree->Branch("b5_HLT_L2Mu10_NoVertex_NoBPTX3BX",&	( b5.HLT_L2Mu10_NoVertex_NoBPTX3BX ));
	_compiledTree->Branch("b5_HLT_L2Mu10_NoVertex_NoBPTX",&	( b5.HLT_L2Mu10_NoVertex_NoBPTX ));
	_compiledTree->Branch("b5_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX",&	( b5.HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX ));
	_compiledTree->Branch("b5_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX",&	( b5.HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX ));
	_compiledTree->Branch("b5_HLT_L2Mu50",&	( b5.HLT_L2Mu50 ));
	_compiledTree->Branch("b5_HLT_L2Mu23NoVtx_2Cha",&	( b5.HLT_L2Mu23NoVtx_2Cha ));
	_compiledTree->Branch("b5_HLT_L2Mu23NoVtx_2Cha_CosmicSeed",&	( b5.HLT_L2Mu23NoVtx_2Cha_CosmicSeed ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4",&	( b5.HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4 ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4",&	( b5.HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4 ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu50",&	( b5.HLT_DoubleL2Mu50 ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed",&	( b5.HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed",&	( b5.HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4",&	( b5.HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4 ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu23NoVtx_2Cha",&	( b5.HLT_DoubleL2Mu23NoVtx_2Cha ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu25NoVtx_2Cha",&	( b5.HLT_DoubleL2Mu25NoVtx_2Cha ));
	_compiledTree->Branch("b5_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4",&	( b5.HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4 ));
	_compiledTree->Branch("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&	( b5.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL ));
	_compiledTree->Branch("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL",&	( b5.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL ));
	_compiledTree->Branch("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&	( b5.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ ));
	_compiledTree->Branch("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ",&	( b5.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ ));
	_compiledTree->Branch("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",&	( b5.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ));
	_compiledTree->Branch("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8",&	( b5.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 ));
	_compiledTree->Branch("b5_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",&	( b5.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 ));
	_compiledTree->Branch("b5_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",&	( b5.HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 ));
	_compiledTree->Branch("b5_HLT_Mu25_TkMu0_Onia",&	( b5.HLT_Mu25_TkMu0_Onia ));
	_compiledTree->Branch("b5_HLT_Mu30_TkMu0_Onia",&	( b5.HLT_Mu30_TkMu0_Onia ));
	_compiledTree->Branch("b5_HLT_Mu20_TkMu0_Phi",&	( b5.HLT_Mu20_TkMu0_Phi ));
	_compiledTree->Branch("b5_HLT_Mu25_TkMu0_Phi",&	( b5.HLT_Mu25_TkMu0_Phi ));
	_compiledTree->Branch("b5_HLT_Mu12",&	( b5.HLT_Mu12 ));
	_compiledTree->Branch("b5_HLT_Mu15",&	( b5.HLT_Mu15 ));
	_compiledTree->Branch("b5_HLT_Mu20",&	( b5.HLT_Mu20 ));
	_compiledTree->Branch("b5_HLT_Mu27",&	( b5.HLT_Mu27 ));
	_compiledTree->Branch("b5_HLT_Mu50",&	( b5.HLT_Mu50 ));
	_compiledTree->Branch("b5_HLT_Mu55",&	( b5.HLT_Mu55 ));
	_compiledTree->Branch("b5_HLT_OldMu100",&	( b5.HLT_OldMu100 ));
	_compiledTree->Branch("b5_HLT_TkMu100",&	( b5.HLT_TkMu100 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve40",&	( b5.HLT_DiPFJetAve40 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve60",&	( b5.HLT_DiPFJetAve60 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve80",&	( b5.HLT_DiPFJetAve80 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve140",&	( b5.HLT_DiPFJetAve140 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve200",&	( b5.HLT_DiPFJetAve200 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve260",&	( b5.HLT_DiPFJetAve260 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve320",&	( b5.HLT_DiPFJetAve320 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve400",&	( b5.HLT_DiPFJetAve400 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve500",&	( b5.HLT_DiPFJetAve500 ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve15_HFJEC",&	( b5.HLT_DiPFJetAve15_HFJEC ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve25_HFJEC",&	( b5.HLT_DiPFJetAve25_HFJEC ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve60_HFJEC",&	( b5.HLT_DiPFJetAve60_HFJEC ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve80_HFJEC",&	( b5.HLT_DiPFJetAve80_HFJEC ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve100_HFJEC",&	( b5.HLT_DiPFJetAve100_HFJEC ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve160_HFJEC",&	( b5.HLT_DiPFJetAve160_HFJEC ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve220_HFJEC",&	( b5.HLT_DiPFJetAve220_HFJEC ));
	_compiledTree->Branch("b5_HLT_DiPFJetAve300_HFJEC",&	( b5.HLT_DiPFJetAve300_HFJEC ));
	_compiledTree->Branch("b5_HLT_AK8PFJet15",&	( b5.HLT_AK8PFJet15 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet25",&	( b5.HLT_AK8PFJet25 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet40",&	( b5.HLT_AK8PFJet40 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet60",&	( b5.HLT_AK8PFJet60 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet80",&	( b5.HLT_AK8PFJet80 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet140",&	( b5.HLT_AK8PFJet140 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet200",&	( b5.HLT_AK8PFJet200 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet260",&	( b5.HLT_AK8PFJet260 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet320",&	( b5.HLT_AK8PFJet320 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet400",&	( b5.HLT_AK8PFJet400 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet450",&	( b5.HLT_AK8PFJet450 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet500",&	( b5.HLT_AK8PFJet500 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet550",&	( b5.HLT_AK8PFJet550 ));
	_compiledTree->Branch("b5_HLT_PFJet15",&	( b5.HLT_PFJet15 ));
	_compiledTree->Branch("b5_HLT_PFJet25",&	( b5.HLT_PFJet25 ));
	_compiledTree->Branch("b5_HLT_PFJet40",&	( b5.HLT_PFJet40 ));
	_compiledTree->Branch("b5_HLT_PFJet60",&	( b5.HLT_PFJet60 ));
	_compiledTree->Branch("b5_HLT_PFJet80",&	( b5.HLT_PFJet80 ));
	_compiledTree->Branch("b5_HLT_PFJet140",&	( b5.HLT_PFJet140 ));
	_compiledTree->Branch("b5_HLT_PFJet200",&	( b5.HLT_PFJet200 ));
	_compiledTree->Branch("b5_HLT_PFJet260",&	( b5.HLT_PFJet260 ));
	_compiledTree->Branch("b5_HLT_PFJet320",&	( b5.HLT_PFJet320 ));
	_compiledTree->Branch("b5_HLT_PFJet400",&	( b5.HLT_PFJet400 ));
	_compiledTree->Branch("b5_HLT_PFJet450",&	( b5.HLT_PFJet450 ));
	_compiledTree->Branch("b5_HLT_PFJet500",&	( b5.HLT_PFJet500 ));
	_compiledTree->Branch("b5_HLT_PFJet550",&	( b5.HLT_PFJet550 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd15",&	( b5.HLT_PFJetFwd15 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd25",&	( b5.HLT_PFJetFwd25 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd40",&	( b5.HLT_PFJetFwd40 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd60",&	( b5.HLT_PFJetFwd60 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd80",&	( b5.HLT_PFJetFwd80 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd140",&	( b5.HLT_PFJetFwd140 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd200",&	( b5.HLT_PFJetFwd200 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd260",&	( b5.HLT_PFJetFwd260 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd320",&	( b5.HLT_PFJetFwd320 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd400",&	( b5.HLT_PFJetFwd400 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd450",&	( b5.HLT_PFJetFwd450 ));
	_compiledTree->Branch("b5_HLT_PFJetFwd500",&	( b5.HLT_PFJetFwd500 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd15",&	( b5.HLT_AK8PFJetFwd15 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd25",&	( b5.HLT_AK8PFJetFwd25 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd40",&	( b5.HLT_AK8PFJetFwd40 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd60",&	( b5.HLT_AK8PFJetFwd60 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd80",&	( b5.HLT_AK8PFJetFwd80 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd140",&	( b5.HLT_AK8PFJetFwd140 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd200",&	( b5.HLT_AK8PFJetFwd200 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd260",&	( b5.HLT_AK8PFJetFwd260 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd320",&	( b5.HLT_AK8PFJetFwd320 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd400",&	( b5.HLT_AK8PFJetFwd400 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd450",&	( b5.HLT_AK8PFJetFwd450 ));
	_compiledTree->Branch("b5_HLT_AK8PFJetFwd500",&	( b5.HLT_AK8PFJetFwd500 ));
	_compiledTree->Branch("b5_HLT_PFHT180",&	( b5.HLT_PFHT180 ));
	_compiledTree->Branch("b5_HLT_PFHT250",&	( b5.HLT_PFHT250 ));
	_compiledTree->Branch("b5_HLT_PFHT370",&	( b5.HLT_PFHT370 ));
	_compiledTree->Branch("b5_HLT_PFHT430",&	( b5.HLT_PFHT430 ));
	_compiledTree->Branch("b5_HLT_PFHT510",&	( b5.HLT_PFHT510 ));
	_compiledTree->Branch("b5_HLT_PFHT590",&	( b5.HLT_PFHT590 ));
	_compiledTree->Branch("b5_HLT_PFHT680",&	( b5.HLT_PFHT680 ));
	_compiledTree->Branch("b5_HLT_PFHT780",&	( b5.HLT_PFHT780 ));
	_compiledTree->Branch("b5_HLT_PFHT890",&	( b5.HLT_PFHT890 ));
	_compiledTree->Branch("b5_HLT_PFHT1050",&	( b5.HLT_PFHT1050 ));
	_compiledTree->Branch("b5_HLT_PFHT500_PFMET100_PFMHT100_IDTight",&	( b5.HLT_PFHT500_PFMET100_PFMHT100_IDTight ));
	_compiledTree->Branch("b5_HLT_PFHT500_PFMET110_PFMHT110_IDTight",&	( b5.HLT_PFHT500_PFMET110_PFMHT110_IDTight ));
	_compiledTree->Branch("b5_HLT_PFHT700_PFMET85_PFMHT85_IDTight",&	( b5.HLT_PFHT700_PFMET85_PFMHT85_IDTight ));
	_compiledTree->Branch("b5_HLT_PFHT700_PFMET95_PFMHT95_IDTight",&	( b5.HLT_PFHT700_PFMET95_PFMHT95_IDTight ));
	_compiledTree->Branch("b5_HLT_PFHT800_PFMET75_PFMHT75_IDTight",&	( b5.HLT_PFHT800_PFMET75_PFMHT75_IDTight ));
	_compiledTree->Branch("b5_HLT_PFHT800_PFMET85_PFMHT85_IDTight",&	( b5.HLT_PFHT800_PFMET85_PFMHT85_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMET110_PFMHT110_IDTight",&	( b5.HLT_PFMET110_PFMHT110_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMET120_PFMHT120_IDTight",&	( b5.HLT_PFMET120_PFMHT120_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMET130_PFMHT130_IDTight",&	( b5.HLT_PFMET130_PFMHT130_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMET140_PFMHT140_IDTight",&	( b5.HLT_PFMET140_PFMHT140_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1",&	( b5.HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1 ));
	_compiledTree->Branch("b5_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1",&	( b5.HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1 ));
	_compiledTree->Branch("b5_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1",&	( b5.HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1 ));
	_compiledTree->Branch("b5_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1",&	( b5.HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1 ));
	_compiledTree->Branch("b5_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1",&	( b5.HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1 ));
	_compiledTree->Branch("b5_HLT_PFMET120_PFMHT120_IDTight_PFHT60",&	( b5.HLT_PFMET120_PFMHT120_IDTight_PFHT60 ));
	_compiledTree->Branch("b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&	( b5.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 ));
	_compiledTree->Branch("b5_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60",&	( b5.HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60 ));
	_compiledTree->Branch("b5_HLT_PFMETTypeOne110_PFMHT110_IDTight",&	( b5.HLT_PFMETTypeOne110_PFMHT110_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMETTypeOne120_PFMHT120_IDTight",&	( b5.HLT_PFMETTypeOne120_PFMHT120_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMETTypeOne130_PFMHT130_IDTight",&	( b5.HLT_PFMETTypeOne130_PFMHT130_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMETTypeOne140_PFMHT140_IDTight",&	( b5.HLT_PFMETTypeOne140_PFMHT140_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",&	( b5.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",&	( b5.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",&	( b5.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight ));
	_compiledTree->Branch("b5_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",&	( b5.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight ));
	_compiledTree->Branch("b5_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight",&	( b5.HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight ));
	_compiledTree->Branch("b5_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight",&	( b5.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight ));
	_compiledTree->Branch("b5_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight",&	( b5.HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight ));
	_compiledTree->Branch("b5_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight",&	( b5.HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight ));
	_compiledTree->Branch("b5_HLT_L1ETMHadSeeds",&	( b5.HLT_L1ETMHadSeeds ));
	_compiledTree->Branch("b5_HLT_CaloMHT90",&	( b5.HLT_CaloMHT90 ));
	_compiledTree->Branch("b5_HLT_CaloMET80_NotCleaned",&	( b5.HLT_CaloMET80_NotCleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET90_NotCleaned",&	( b5.HLT_CaloMET90_NotCleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET100_NotCleaned",&	( b5.HLT_CaloMET100_NotCleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET110_NotCleaned",&	( b5.HLT_CaloMET110_NotCleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET250_NotCleaned",&	( b5.HLT_CaloMET250_NotCleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET70_HBHECleaned",&	( b5.HLT_CaloMET70_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET80_HBHECleaned",&	( b5.HLT_CaloMET80_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET90_HBHECleaned",&	( b5.HLT_CaloMET90_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET100_HBHECleaned",&	( b5.HLT_CaloMET100_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET250_HBHECleaned",&	( b5.HLT_CaloMET250_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET300_HBHECleaned",&	( b5.HLT_CaloMET300_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_CaloMET350_HBHECleaned",&	( b5.HLT_CaloMET350_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_PFMET200_NotCleaned",&	( b5.HLT_PFMET200_NotCleaned ));
	_compiledTree->Branch("b5_HLT_PFMET200_HBHECleaned",&	( b5.HLT_PFMET200_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_PFMET250_HBHECleaned",&	( b5.HLT_PFMET250_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_PFMET300_HBHECleaned",&	( b5.HLT_PFMET300_HBHECleaned ));
	_compiledTree->Branch("b5_HLT_PFMET200_HBHE_BeamHaloCleaned",&	( b5.HLT_PFMET200_HBHE_BeamHaloCleaned ));
	_compiledTree->Branch("b5_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",&	( b5.HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned ));
	_compiledTree->Branch("b5_HLT_MET105_IsoTrk50",&	( b5.HLT_MET105_IsoTrk50 ));
	_compiledTree->Branch("b5_HLT_MET120_IsoTrk50",&	( b5.HLT_MET120_IsoTrk50 ));
	_compiledTree->Branch("b5_HLT_SingleJet30_Mu12_SinglePFJet40",&	( b5.HLT_SingleJet30_Mu12_SinglePFJet40 ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71",&	( b5.HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71",&	( b5.HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71",&	( b5.HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71",&	( b5.HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",&	( b5.HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",&	( b5.HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",&	( b5.HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_DoublePFJets40_CaloBTagDeepCSV_p71",&	( b5.HLT_DoublePFJets40_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_DoublePFJets100_CaloBTagDeepCSV_p71",&	( b5.HLT_DoublePFJets100_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_DoublePFJets200_CaloBTagDeepCSV_p71",&	( b5.HLT_DoublePFJets200_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_DoublePFJets350_CaloBTagDeepCSV_p71",&	( b5.HLT_DoublePFJets350_CaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",&	( b5.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",&	( b5.HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ));
	_compiledTree->Branch("b5_HLT_Photon300_NoHE",&	( b5.HLT_Photon300_NoHE ));
	_compiledTree->Branch("b5_HLT_Mu8_TrkIsoVVL",&	( b5.HLT_Mu8_TrkIsoVVL ));
	_compiledTree->Branch("b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ",&	( b5.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ ));
	_compiledTree->Branch("b5_HLT_Mu8_DiEle12_CaloIdL_TrackIdL",&	( b5.HLT_Mu8_DiEle12_CaloIdL_TrackIdL ));
	_compiledTree->Branch("b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ",&	( b5.HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ ));
	_compiledTree->Branch("b5_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350",&	( b5.HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350 ));
	_compiledTree->Branch("b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",&	( b5.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ));
	_compiledTree->Branch("b5_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",&	( b5.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL ));
	_compiledTree->Branch("b5_HLT_Mu17_TrkIsoVVL",&	( b5.HLT_Mu17_TrkIsoVVL ));
	_compiledTree->Branch("b5_HLT_Mu19_TrkIsoVVL",&	( b5.HLT_Mu19_TrkIsoVVL ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK4DiJet20_Mu5",&	( b5.HLT_BTagMu_AK4DiJet20_Mu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK4DiJet40_Mu5",&	( b5.HLT_BTagMu_AK4DiJet40_Mu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK4DiJet70_Mu5",&	( b5.HLT_BTagMu_AK4DiJet70_Mu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK4DiJet110_Mu5",&	( b5.HLT_BTagMu_AK4DiJet110_Mu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK4DiJet170_Mu5",&	( b5.HLT_BTagMu_AK4DiJet170_Mu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK4Jet300_Mu5",&	( b5.HLT_BTagMu_AK4Jet300_Mu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK8DiJet170_Mu5",&	( b5.HLT_BTagMu_AK8DiJet170_Mu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK8Jet170_DoubleMu5",&	( b5.HLT_BTagMu_AK8Jet170_DoubleMu5 ));
	_compiledTree->Branch("b5_HLT_BTagMu_AK8Jet300_Mu5",&	( b5.HLT_BTagMu_AK8Jet300_Mu5 ));
	_compiledTree->Branch("b5_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL",&	( b5.HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL ));
	_compiledTree->Branch("b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&	( b5.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ));
	_compiledTree->Branch("b5_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&	( b5.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ));
	_compiledTree->Branch("b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&	( b5.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ));
	_compiledTree->Branch("b5_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",&	( b5.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL ));
	_compiledTree->Branch("b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",&	( b5.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL ));
	_compiledTree->Branch("b5_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",&	( b5.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ));
	_compiledTree->Branch("b5_HLT_Mu12_DoublePhoton20",&	( b5.HLT_Mu12_DoublePhoton20 ));
	_compiledTree->Branch("b5_HLT_TriplePhoton_20_20_20_CaloIdLV2",&	( b5.HLT_TriplePhoton_20_20_20_CaloIdLV2 ));
	_compiledTree->Branch("b5_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL",&	( b5.HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL ));
	_compiledTree->Branch("b5_HLT_TriplePhoton_30_30_10_CaloIdLV2",&	( b5.HLT_TriplePhoton_30_30_10_CaloIdLV2 ));
	_compiledTree->Branch("b5_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL",&	( b5.HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL ));
	_compiledTree->Branch("b5_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL",&	( b5.HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL ));
	_compiledTree->Branch("b5_HLT_Photon20",&	( b5.HLT_Photon20 ));
	_compiledTree->Branch("b5_HLT_Photon33",&	( b5.HLT_Photon33 ));
	_compiledTree->Branch("b5_HLT_Photon50",&	( b5.HLT_Photon50 ));
	_compiledTree->Branch("b5_HLT_Photon75",&	( b5.HLT_Photon75 ));
	_compiledTree->Branch("b5_HLT_Photon90",&	( b5.HLT_Photon90 ));
	_compiledTree->Branch("b5_HLT_Photon120",&	( b5.HLT_Photon120 ));
	_compiledTree->Branch("b5_HLT_Photon150",&	( b5.HLT_Photon150 ));
	_compiledTree->Branch("b5_HLT_Photon175",&	( b5.HLT_Photon175 ));
	_compiledTree->Branch("b5_HLT_Photon200",&	( b5.HLT_Photon200 ));
	_compiledTree->Branch("b5_HLT_Photon100EB_TightID_TightIso",&	( b5.HLT_Photon100EB_TightID_TightIso ));
	_compiledTree->Branch("b5_HLT_Photon110EB_TightID_TightIso",&	( b5.HLT_Photon110EB_TightID_TightIso ));
	_compiledTree->Branch("b5_HLT_Photon120EB_TightID_TightIso",&	( b5.HLT_Photon120EB_TightID_TightIso ));
	_compiledTree->Branch("b5_HLT_Photon100EBHE10",&	( b5.HLT_Photon100EBHE10 ));
	_compiledTree->Branch("b5_HLT_Photon100EEHE10",&	( b5.HLT_Photon100EEHE10 ));
	_compiledTree->Branch("b5_HLT_Photon100EE_TightID_TightIso",&	( b5.HLT_Photon100EE_TightID_TightIso ));
	_compiledTree->Branch("b5_HLT_Photon50_R9Id90_HE10_IsoM",&	( b5.HLT_Photon50_R9Id90_HE10_IsoM ));
	_compiledTree->Branch("b5_HLT_Photon75_R9Id90_HE10_IsoM",&	( b5.HLT_Photon75_R9Id90_HE10_IsoM ));
	_compiledTree->Branch("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3",&	( b5.HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3 ));
	_compiledTree->Branch("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3",&	( b5.HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3 ));
	_compiledTree->Branch("b5_HLT_Photon90_R9Id90_HE10_IsoM",&	( b5.HLT_Photon90_R9Id90_HE10_IsoM ));
	_compiledTree->Branch("b5_HLT_Photon120_R9Id90_HE10_IsoM",&	( b5.HLT_Photon120_R9Id90_HE10_IsoM ));
	_compiledTree->Branch("b5_HLT_Photon165_R9Id90_HE10_IsoM",&	( b5.HLT_Photon165_R9Id90_HE10_IsoM ));
	_compiledTree->Branch("b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",&	( b5.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 ));
	_compiledTree->Branch("b5_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95",&	( b5.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95 ));
	_compiledTree->Branch("b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55",&	( b5.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 ));
	_compiledTree->Branch("b5_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55",&	( b5.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Jpsi_L1_NoOS",&	( b5.HLT_Dimuon0_Jpsi_L1_NoOS ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Jpsi_NoVertexing_NoOS",&	( b5.HLT_Dimuon0_Jpsi_NoVertexing_NoOS ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Jpsi",&	( b5.HLT_Dimuon0_Jpsi ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Jpsi_NoVertexing",&	( b5.HLT_Dimuon0_Jpsi_NoVertexing ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",&	( b5.HLT_Dimuon0_Jpsi_L1_4R_0er1p5R ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",&	( b5.HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Jpsi3p5_Muon2",&	( b5.HLT_Dimuon0_Jpsi3p5_Muon2 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_L1_4p5",&	( b5.HLT_Dimuon0_Upsilon_L1_4p5 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_L1_5",&	( b5.HLT_Dimuon0_Upsilon_L1_5 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_L1_4p5NoOS",&	( b5.HLT_Dimuon0_Upsilon_L1_4p5NoOS ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0",&	( b5.HLT_Dimuon0_Upsilon_L1_4p5er2p0 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_L1_4p5er2p0M",&	( b5.HLT_Dimuon0_Upsilon_L1_4p5er2p0M ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_NoVertexing",&	( b5.HLT_Dimuon0_Upsilon_NoVertexing ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_L1_5M",&	( b5.HLT_Dimuon0_Upsilon_L1_5M ));
	_compiledTree->Branch("b5_HLT_Dimuon0_LowMass_L1_0er1p5R",&	( b5.HLT_Dimuon0_LowMass_L1_0er1p5R ));
	_compiledTree->Branch("b5_HLT_Dimuon0_LowMass_L1_0er1p5",&	( b5.HLT_Dimuon0_LowMass_L1_0er1p5 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_LowMass",&	( b5.HLT_Dimuon0_LowMass ));
	_compiledTree->Branch("b5_HLT_Dimuon0_LowMass_L1_4",&	( b5.HLT_Dimuon0_LowMass_L1_4 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_LowMass_L1_4R",&	( b5.HLT_Dimuon0_LowMass_L1_4R ));
	_compiledTree->Branch("b5_HLT_Dimuon0_LowMass_L1_TM530",&	( b5.HLT_Dimuon0_LowMass_L1_TM530 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_Muon_L1_TM0",&	( b5.HLT_Dimuon0_Upsilon_Muon_L1_TM0 ));
	_compiledTree->Branch("b5_HLT_Dimuon0_Upsilon_Muon_NoL1Mass",&	( b5.HLT_Dimuon0_Upsilon_Muon_NoL1Mass ));
	_compiledTree->Branch("b5_HLT_TripleMu_5_3_3_Mass3p8_DZ",&	( b5.HLT_TripleMu_5_3_3_Mass3p8_DZ ));
	_compiledTree->Branch("b5_HLT_TripleMu_10_5_5_DZ",&	( b5.HLT_TripleMu_10_5_5_DZ ));
	_compiledTree->Branch("b5_HLT_TripleMu_12_10_5",&	( b5.HLT_TripleMu_12_10_5 ));
	_compiledTree->Branch("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15",&	( b5.HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15 ));
	_compiledTree->Branch("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1",&	( b5.HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1 ));
	_compiledTree->Branch("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15",&	( b5.HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 ));
	_compiledTree->Branch("b5_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1",&	( b5.HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_DZ_PFMET50_PFMHT60",&	( b5.HLT_DoubleMu3_DZ_PFMET50_PFMHT60 ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_DZ_PFMET70_PFMHT70",&	( b5.HLT_DoubleMu3_DZ_PFMET70_PFMHT70 ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_DZ_PFMET90_PFMHT90",&	( b5.HLT_DoubleMu3_DZ_PFMET90_PFMHT90 ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass",&	( b5.HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_Jpsi_Displaced",&	( b5.HLT_DoubleMu4_Jpsi_Displaced ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_Jpsi_NoVertexing",&	( b5.HLT_DoubleMu4_Jpsi_NoVertexing ));
	_compiledTree->Branch("b5_HLT_DoubleMu4_JpsiTrkTrk_Displaced",&	( b5.HLT_DoubleMu4_JpsiTrkTrk_Displaced ));
	_compiledTree->Branch("b5_HLT_DoubleMu43NoFiltersNoVtx",&	( b5.HLT_DoubleMu43NoFiltersNoVtx ));
	_compiledTree->Branch("b5_HLT_DoubleMu48NoFiltersNoVtx",&	( b5.HLT_DoubleMu48NoFiltersNoVtx ));
	_compiledTree->Branch("b5_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL",&	( b5.HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL ));
	_compiledTree->Branch("b5_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL",&	( b5.HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL ));
	_compiledTree->Branch("b5_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL",&	( b5.HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL ));
	_compiledTree->Branch("b5_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL",&	( b5.HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL ));
	_compiledTree->Branch("b5_HLT_DoubleMu33NoFiltersNoVtxDisplaced",&	( b5.HLT_DoubleMu33NoFiltersNoVtxDisplaced ));
	_compiledTree->Branch("b5_HLT_DoubleMu40NoFiltersNoVtxDisplaced",&	( b5.HLT_DoubleMu40NoFiltersNoVtxDisplaced ));
	_compiledTree->Branch("b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4",&	( b5.HLT_DoubleMu20_7_Mass0to30_L1_DM4 ));
	_compiledTree->Branch("b5_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG",&	( b5.HLT_DoubleMu20_7_Mass0to30_L1_DM4EG ));
	_compiledTree->Branch("b5_HLT_HT425",&	( b5.HLT_HT425 ));
	_compiledTree->Branch("b5_HLT_HT430_DisplacedDijet40_DisplacedTrack",&	( b5.HLT_HT430_DisplacedDijet40_DisplacedTrack ));
	_compiledTree->Branch("b5_HLT_HT500_DisplacedDijet40_DisplacedTrack",&	( b5.HLT_HT500_DisplacedDijet40_DisplacedTrack ));
	_compiledTree->Branch("b5_HLT_HT430_DisplacedDijet60_DisplacedTrack",&	( b5.HLT_HT430_DisplacedDijet60_DisplacedTrack ));
	_compiledTree->Branch("b5_HLT_HT400_DisplacedDijet40_DisplacedTrack",&	( b5.HLT_HT400_DisplacedDijet40_DisplacedTrack ));
	_compiledTree->Branch("b5_HLT_HT650_DisplacedDijet60_Inclusive",&	( b5.HLT_HT650_DisplacedDijet60_Inclusive ));
	_compiledTree->Branch("b5_HLT_HT550_DisplacedDijet60_Inclusive",&	( b5.HLT_HT550_DisplacedDijet60_Inclusive ));
	_compiledTree->Branch("b5_HLT_DiJet110_35_Mjj650_PFMET110",&	( b5.HLT_DiJet110_35_Mjj650_PFMET110 ));
	_compiledTree->Branch("b5_HLT_DiJet110_35_Mjj650_PFMET120",&	( b5.HLT_DiJet110_35_Mjj650_PFMET120 ));
	_compiledTree->Branch("b5_HLT_DiJet110_35_Mjj650_PFMET130",&	( b5.HLT_DiJet110_35_Mjj650_PFMET130 ));
	_compiledTree->Branch("b5_HLT_TripleJet110_35_35_Mjj650_PFMET110",&	( b5.HLT_TripleJet110_35_35_Mjj650_PFMET110 ));
	_compiledTree->Branch("b5_HLT_TripleJet110_35_35_Mjj650_PFMET120",&	( b5.HLT_TripleJet110_35_35_Mjj650_PFMET120 ));
	_compiledTree->Branch("b5_HLT_TripleJet110_35_35_Mjj650_PFMET130",&	( b5.HLT_TripleJet110_35_35_Mjj650_PFMET130 ));
	_compiledTree->Branch("b5_HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1",&	( b5.HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1 ));
	_compiledTree->Branch("b5_HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1",&	( b5.HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1 ));
	_compiledTree->Branch("b5_HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1",&	( b5.HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1 ));
	_compiledTree->Branch("b5_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned",&	( b5.HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned ));
	_compiledTree->Branch("b5_HLT_Ele28_eta2p1_WPTight_Gsf_HT150",&	( b5.HLT_Ele28_eta2p1_WPTight_Gsf_HT150 ));
	_compiledTree->Branch("b5_HLT_Ele28_HighEta_SC20_Mass55",&	( b5.HLT_Ele28_HighEta_SC20_Mass55 ));
	_compiledTree->Branch("b5_HLT_DoubleMu20_7_Mass0to30_Photon23",&	( b5.HLT_DoubleMu20_7_Mass0to30_Photon23 ));
	_compiledTree->Branch("b5_HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5",&	( b5.HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5 ));
	_compiledTree->Branch("b5_HLT_Ele15_IsoVVVL_PFHT450_PFMET50",&	( b5.HLT_Ele15_IsoVVVL_PFHT450_PFMET50 ));
	_compiledTree->Branch("b5_HLT_Ele15_IsoVVVL_PFHT450",&	( b5.HLT_Ele15_IsoVVVL_PFHT450 ));
	_compiledTree->Branch("b5_HLT_Ele50_IsoVVVL_PFHT450",&	( b5.HLT_Ele50_IsoVVVL_PFHT450 ));
	_compiledTree->Branch("b5_HLT_Ele15_IsoVVVL_PFHT600",&	( b5.HLT_Ele15_IsoVVVL_PFHT600 ));
	_compiledTree->Branch("b5_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",&	( b5.HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60 ));
	_compiledTree->Branch("b5_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",&	( b5.HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60 ));
	_compiledTree->Branch("b5_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60",&	( b5.HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60 ));
	_compiledTree->Branch("b5_HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5",&	( b5.HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5 ));
	_compiledTree->Branch("b5_HLT_Mu15_IsoVVVL_PFHT450_PFMET50",&	( b5.HLT_Mu15_IsoVVVL_PFHT450_PFMET50 ));
	_compiledTree->Branch("b5_HLT_Mu15_IsoVVVL_PFHT450",&	( b5.HLT_Mu15_IsoVVVL_PFHT450 ));
	_compiledTree->Branch("b5_HLT_Mu50_IsoVVVL_PFHT450",&	( b5.HLT_Mu50_IsoVVVL_PFHT450 ));
	_compiledTree->Branch("b5_HLT_Mu15_IsoVVVL_PFHT600",&	( b5.HLT_Mu15_IsoVVVL_PFHT600 ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight ));
	_compiledTree->Branch("b5_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight",&	( b5.HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight ));
	_compiledTree->Branch("b5_HLT_Dimuon10_PsiPrime_Barrel_Seagulls",&	( b5.HLT_Dimuon10_PsiPrime_Barrel_Seagulls ));
	_compiledTree->Branch("b5_HLT_Dimuon20_Jpsi_Barrel_Seagulls",&	( b5.HLT_Dimuon20_Jpsi_Barrel_Seagulls ));
	_compiledTree->Branch("b5_HLT_Dimuon12_Upsilon_y1p4",&	( b5.HLT_Dimuon12_Upsilon_y1p4 ));
	_compiledTree->Branch("b5_HLT_Dimuon14_Phi_Barrel_Seagulls",&	( b5.HLT_Dimuon14_Phi_Barrel_Seagulls ));
	_compiledTree->Branch("b5_HLT_Dimuon18_PsiPrime",&	( b5.HLT_Dimuon18_PsiPrime ));
	_compiledTree->Branch("b5_HLT_Dimuon25_Jpsi",&	( b5.HLT_Dimuon25_Jpsi ));
	_compiledTree->Branch("b5_HLT_Dimuon18_PsiPrime_noCorrL1",&	( b5.HLT_Dimuon18_PsiPrime_noCorrL1 ));
	_compiledTree->Branch("b5_HLT_Dimuon24_Upsilon_noCorrL1",&	( b5.HLT_Dimuon24_Upsilon_noCorrL1 ));
	_compiledTree->Branch("b5_HLT_Dimuon24_Phi_noCorrL1",&	( b5.HLT_Dimuon24_Phi_noCorrL1 ));
	_compiledTree->Branch("b5_HLT_Dimuon25_Jpsi_noCorrL1",&	( b5.HLT_Dimuon25_Jpsi_noCorrL1 ));
	_compiledTree->Branch("b5_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8",&	( b5.HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8 ));
	_compiledTree->Branch("b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",&	( b5.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ ));
	_compiledTree->Branch("b5_HLT_DiMu9_Ele9_CaloIdL_TrackIdL",&	( b5.HLT_DiMu9_Ele9_CaloIdL_TrackIdL ));
	_compiledTree->Branch("b5_HLT_DoubleIsoMu20_eta2p1",&	( b5.HLT_DoubleIsoMu20_eta2p1 ));
	_compiledTree->Branch("b5_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx",&	( b5.HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx ));
	_compiledTree->Branch("b5_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx",&	( b5.HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx ));
	_compiledTree->Branch("b5_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx",&	( b5.HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx ));
	_compiledTree->Branch("b5_HLT_Mu8",&	( b5.HLT_Mu8 ));
	_compiledTree->Branch("b5_HLT_Mu17",&	( b5.HLT_Mu17 ));
	_compiledTree->Branch("b5_HLT_Mu19",&	( b5.HLT_Mu19 ));
	_compiledTree->Branch("b5_HLT_Mu17_Photon30_IsoCaloId",&	( b5.HLT_Mu17_Photon30_IsoCaloId ));
	_compiledTree->Branch("b5_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",&	( b5.HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30 ));
	_compiledTree->Branch("b5_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",&	( b5.HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 ));
	_compiledTree->Branch("b5_HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30",&	( b5.HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30 ));
	_compiledTree->Branch("b5_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",&	( b5.HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30 ));
	_compiledTree->Branch("b5_HLT_Ele8_CaloIdM_TrackIdM_PFJet30",&	( b5.HLT_Ele8_CaloIdM_TrackIdM_PFJet30 ));
	_compiledTree->Branch("b5_HLT_Ele17_CaloIdM_TrackIdM_PFJet30",&	( b5.HLT_Ele17_CaloIdM_TrackIdM_PFJet30 ));
	_compiledTree->Branch("b5_HLT_Ele23_CaloIdM_TrackIdM_PFJet30",&	( b5.HLT_Ele23_CaloIdM_TrackIdM_PFJet30 ));
	_compiledTree->Branch("b5_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",&	( b5.HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 ));
	_compiledTree->Branch("b5_HLT_Ele115_CaloIdVT_GsfTrkIdT",&	( b5.HLT_Ele115_CaloIdVT_GsfTrkIdT ));
	_compiledTree->Branch("b5_HLT_Ele135_CaloIdVT_GsfTrkIdT",&	( b5.HLT_Ele135_CaloIdVT_GsfTrkIdT ));
	_compiledTree->Branch("b5_HLT_Ele145_CaloIdVT_GsfTrkIdT",&	( b5.HLT_Ele145_CaloIdVT_GsfTrkIdT ));
	_compiledTree->Branch("b5_HLT_Ele200_CaloIdVT_GsfTrkIdT",&	( b5.HLT_Ele200_CaloIdVT_GsfTrkIdT ));
	_compiledTree->Branch("b5_HLT_Ele250_CaloIdVT_GsfTrkIdT",&	( b5.HLT_Ele250_CaloIdVT_GsfTrkIdT ));
	_compiledTree->Branch("b5_HLT_Ele300_CaloIdVT_GsfTrkIdT",&	( b5.HLT_Ele300_CaloIdVT_GsfTrkIdT ));
	_compiledTree->Branch("b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5",&	( b5.HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 ));
	_compiledTree->Branch("b5_HLT_PFHT330PT30_QuadPFJet_75_60_45_40",&	( b5.HLT_PFHT330PT30_QuadPFJet_75_60_45_40 ));
	_compiledTree->Branch("b5_HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",&	( b5.HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2 ));
	_compiledTree->Branch("b5_HLT_PFHT380_SixPFJet32",&	( b5.HLT_PFHT380_SixPFJet32 ));
	_compiledTree->Branch("b5_HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5",&	( b5.HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5 ));
	_compiledTree->Branch("b5_HLT_PFHT430_SixPFJet40",&	( b5.HLT_PFHT430_SixPFJet40 ));
	_compiledTree->Branch("b5_HLT_PFHT350",&	( b5.HLT_PFHT350 ));
	_compiledTree->Branch("b5_HLT_PFHT350MinPFJet15",&	( b5.HLT_PFHT350MinPFJet15 ));
	_compiledTree->Branch("b5_HLT_Photon60_R9Id90_CaloIdL_IsoL",&	( b5.HLT_Photon60_R9Id90_CaloIdL_IsoL ));
	_compiledTree->Branch("b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL",&	( b5.HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL ));
	_compiledTree->Branch("b5_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15",&	( b5.HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15 ));
	_compiledTree->Branch("b5_HLT_ECALHT800",&	( b5.HLT_ECALHT800 ));
	_compiledTree->Branch("b5_HLT_DiSC30_18_EIso_AND_HE_Mass70",&	( b5.HLT_DiSC30_18_EIso_AND_HE_Mass70 ));
	_compiledTree->Branch("b5_HLT_Physics",&	( b5.HLT_Physics ));
	_compiledTree->Branch("b5_HLT_Physics_part0",&	( b5.HLT_Physics_part0 ));
	_compiledTree->Branch("b5_HLT_Physics_part1",&	( b5.HLT_Physics_part1 ));
	_compiledTree->Branch("b5_HLT_Physics_part2",&	( b5.HLT_Physics_part2 ));
	_compiledTree->Branch("b5_HLT_Physics_part3",&	( b5.HLT_Physics_part3 ));
	_compiledTree->Branch("b5_HLT_Physics_part4",&	( b5.HLT_Physics_part4 ));
	_compiledTree->Branch("b5_HLT_Physics_part5",&	( b5.HLT_Physics_part5 ));
	_compiledTree->Branch("b5_HLT_Physics_part6",&	( b5.HLT_Physics_part6 ));
	_compiledTree->Branch("b5_HLT_Physics_part7",&	( b5.HLT_Physics_part7 ));
	_compiledTree->Branch("b5_HLT_Random",&	( b5.HLT_Random ));
	_compiledTree->Branch("b5_HLT_ZeroBias",&	( b5.HLT_ZeroBias ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part0",&	( b5.HLT_ZeroBias_part0 ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part1",&	( b5.HLT_ZeroBias_part1 ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part2",&	( b5.HLT_ZeroBias_part2 ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part3",&	( b5.HLT_ZeroBias_part3 ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part4",&	( b5.HLT_ZeroBias_part4 ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part5",&	( b5.HLT_ZeroBias_part5 ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part6",&	( b5.HLT_ZeroBias_part6 ));
	_compiledTree->Branch("b5_HLT_ZeroBias_part7",&	( b5.HLT_ZeroBias_part7 ));
	_compiledTree->Branch("b5_HLT_AK4CaloJet30",&	( b5.HLT_AK4CaloJet30 ));
	_compiledTree->Branch("b5_HLT_AK4CaloJet40",&	( b5.HLT_AK4CaloJet40 ));
	_compiledTree->Branch("b5_HLT_AK4CaloJet50",&	( b5.HLT_AK4CaloJet50 ));
	_compiledTree->Branch("b5_HLT_AK4CaloJet80",&	( b5.HLT_AK4CaloJet80 ));
	_compiledTree->Branch("b5_HLT_AK4CaloJet100",&	( b5.HLT_AK4CaloJet100 ));
	_compiledTree->Branch("b5_HLT_AK4CaloJet120",&	( b5.HLT_AK4CaloJet120 ));
	_compiledTree->Branch("b5_HLT_AK4PFJet30",&	( b5.HLT_AK4PFJet30 ));
	_compiledTree->Branch("b5_HLT_AK4PFJet50",&	( b5.HLT_AK4PFJet50 ));
	_compiledTree->Branch("b5_HLT_AK4PFJet80",&	( b5.HLT_AK4PFJet80 ));
	_compiledTree->Branch("b5_HLT_AK4PFJet100",&	( b5.HLT_AK4PFJet100 ));
	_compiledTree->Branch("b5_HLT_AK4PFJet120",&	( b5.HLT_AK4PFJet120 ));
	_compiledTree->Branch("b5_HLT_SinglePhoton10_Eta3p1ForPPRef",&	( b5.HLT_SinglePhoton10_Eta3p1ForPPRef ));
	_compiledTree->Branch("b5_HLT_SinglePhoton20_Eta3p1ForPPRef",&	( b5.HLT_SinglePhoton20_Eta3p1ForPPRef ));
	_compiledTree->Branch("b5_HLT_SinglePhoton30_Eta3p1ForPPRef",&	( b5.HLT_SinglePhoton30_Eta3p1ForPPRef ));
	_compiledTree->Branch("b5_HLT_Photon20_HoverELoose",&	( b5.HLT_Photon20_HoverELoose ));
	_compiledTree->Branch("b5_HLT_Photon30_HoverELoose",&	( b5.HLT_Photon30_HoverELoose ));
	_compiledTree->Branch("b5_HLT_EcalCalibration",&	( b5.HLT_EcalCalibration ));
	_compiledTree->Branch("b5_HLT_HcalCalibration",&	( b5.HLT_HcalCalibration ));
	_compiledTree->Branch("b5_HLT_L1UnpairedBunchBptxMinus",&	( b5.HLT_L1UnpairedBunchBptxMinus ));
	_compiledTree->Branch("b5_HLT_L1UnpairedBunchBptxPlus",&	( b5.HLT_L1UnpairedBunchBptxPlus ));
	_compiledTree->Branch("b5_HLT_L1NotBptxOR",&	( b5.HLT_L1NotBptxOR ));
	_compiledTree->Branch("b5_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142",&	( b5.HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 ));
	_compiledTree->Branch("b5_HLT_HcalNZS",&	( b5.HLT_HcalNZS ));
	_compiledTree->Branch("b5_HLT_HcalPhiSym",&	( b5.HLT_HcalPhiSym ));
	_compiledTree->Branch("b5_HLT_HcalIsolatedbunch",&	( b5.HLT_HcalIsolatedbunch ));
	_compiledTree->Branch("b5_HLT_IsoTrackHB",&	( b5.HLT_IsoTrackHB ));
	_compiledTree->Branch("b5_HLT_IsoTrackHE",&	( b5.HLT_IsoTrackHE ));
	_compiledTree->Branch("b5_HLT_ZeroBias_FirstCollisionAfterAbortGap",&	( b5.HLT_ZeroBias_FirstCollisionAfterAbortGap ));
	_compiledTree->Branch("b5_HLT_ZeroBias_IsolatedBunches",&	( b5.HLT_ZeroBias_IsolatedBunches ));
	_compiledTree->Branch("b5_HLT_ZeroBias_FirstCollisionInTrain",&	( b5.HLT_ZeroBias_FirstCollisionInTrain ));
	_compiledTree->Branch("b5_HLT_ZeroBias_LastCollisionInTrain",&	( b5.HLT_ZeroBias_LastCollisionInTrain ));
	_compiledTree->Branch("b5_HLT_ZeroBias_FirstBXAfterTrain",&	( b5.HLT_ZeroBias_FirstBXAfterTrain ));
	_compiledTree->Branch("b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",&	( b5.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1 ));
	_compiledTree->Branch("b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1",&	( b5.HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1 ));
	_compiledTree->Branch("b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1",&	( b5.HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1 ));
	_compiledTree->Branch("b5_HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1",&	( b5.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1 ));
	_compiledTree->Branch("b5_HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1",&	( b5.HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1 ));
	_compiledTree->Branch("b5_HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1",&	( b5.HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1 ));
	_compiledTree->Branch("b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg",&	( b5.HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg",&	( b5.HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg",&	( b5.HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",&	( b5.HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",&	( b5.HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",&	( b5.HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",&	( b5.HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",&	( b5.HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",&	( b5.HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90",&	( b5.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100",&	( b5.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110",&	( b5.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120",&	( b5.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130",&	( b5.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140",&	( b5.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",&	( b5.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr",&	( b5.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1",&	( b5.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1",&	( b5.HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1 ));
	_compiledTree->Branch("b5_HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1",&	( b5.HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1 ));
	_compiledTree->Branch("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1",&	( b5.HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1",&	( b5.HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1",&	( b5.HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1",&	( b5.HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1 ));
	_compiledTree->Branch("b5_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",&	( b5.HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL ));
	_compiledTree->Branch("b5_HLT_Rsq0p35",&	( b5.HLT_Rsq0p35 ));
	_compiledTree->Branch("b5_HLT_Rsq0p40",&	( b5.HLT_Rsq0p40 ));
	_compiledTree->Branch("b5_HLT_RsqMR300_Rsq0p09_MR200",&	( b5.HLT_RsqMR300_Rsq0p09_MR200 ));
	_compiledTree->Branch("b5_HLT_RsqMR320_Rsq0p09_MR200",&	( b5.HLT_RsqMR320_Rsq0p09_MR200 ));
	_compiledTree->Branch("b5_HLT_RsqMR300_Rsq0p09_MR200_4jet",&	( b5.HLT_RsqMR300_Rsq0p09_MR200_4jet ));
	_compiledTree->Branch("b5_HLT_RsqMR320_Rsq0p09_MR200_4jet",&	( b5.HLT_RsqMR320_Rsq0p09_MR200_4jet ));
	_compiledTree->Branch("b5_HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1",&	( b5.HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1",&	( b5.HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1",&	( b5.HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu27_MET90",&	( b5.HLT_IsoMu27_MET90 ));
	_compiledTree->Branch("b5_HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1",&	( b5.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 ));
	_compiledTree->Branch("b5_HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1",&	( b5.HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1 ));
	_compiledTree->Branch("b5_HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",&	( b5.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg ));
	_compiledTree->Branch("b5_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50",&	( b5.HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50 ));
	_compiledTree->Branch("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3",&	( b5.HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3 ));
	_compiledTree->Branch("b5_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3",&	( b5.HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3 ));
	_compiledTree->Branch("b5_HLT_PFMET100_PFMHT100_IDTight_PFHT60",&	( b5.HLT_PFMET100_PFMHT100_IDTight_PFHT60 ));
	_compiledTree->Branch("b5_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60",&	( b5.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60 ));
	_compiledTree->Branch("b5_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60",&	( b5.HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60 ));
	_compiledTree->Branch("b5_HLT_Mu18_Mu9_SameSign",&	( b5.HLT_Mu18_Mu9_SameSign ));
	_compiledTree->Branch("b5_HLT_Mu18_Mu9_SameSign_DZ",&	( b5.HLT_Mu18_Mu9_SameSign_DZ ));
	_compiledTree->Branch("b5_HLT_Mu18_Mu9",&	( b5.HLT_Mu18_Mu9 ));
	_compiledTree->Branch("b5_HLT_Mu18_Mu9_DZ",&	( b5.HLT_Mu18_Mu9_DZ ));
	_compiledTree->Branch("b5_HLT_Mu20_Mu10_SameSign",&	( b5.HLT_Mu20_Mu10_SameSign ));
	_compiledTree->Branch("b5_HLT_Mu20_Mu10_SameSign_DZ",&	( b5.HLT_Mu20_Mu10_SameSign_DZ ));
	_compiledTree->Branch("b5_HLT_Mu20_Mu10",&	( b5.HLT_Mu20_Mu10 ));
	_compiledTree->Branch("b5_HLT_Mu20_Mu10_DZ",&	( b5.HLT_Mu20_Mu10_DZ ));
	_compiledTree->Branch("b5_HLT_Mu23_Mu12_SameSign",&	( b5.HLT_Mu23_Mu12_SameSign ));
	_compiledTree->Branch("b5_HLT_Mu23_Mu12_SameSign_DZ",&	( b5.HLT_Mu23_Mu12_SameSign_DZ ));
	_compiledTree->Branch("b5_HLT_Mu23_Mu12",&	( b5.HLT_Mu23_Mu12 ));
	_compiledTree->Branch("b5_HLT_Mu23_Mu12_DZ",&	( b5.HLT_Mu23_Mu12_DZ ));
	_compiledTree->Branch("b5_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",&	( b5.HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05 ));
	_compiledTree->Branch("b5_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",&	( b5.HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi ));
	_compiledTree->Branch("b5_HLT_DoubleMu3_DCA_PFMET50_PFMHT60",&	( b5.HLT_DoubleMu3_DCA_PFMET50_PFMHT60 ));
	_compiledTree->Branch("b5_HLT_TripleMu_5_3_3_Mass3p8_DCA",&	( b5.HLT_TripleMu_5_3_3_Mass3p8_DCA ));
	_compiledTree->Branch("b5_HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",&	( b5.HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",&	( b5.HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet105_90_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",&	( b5.HLT_QuadPFJet105_90_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",&	( b5.HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2",&	( b5.HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2",&	( b5.HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2",&	( b5.HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2",&	( b5.HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet98_83_71_15",&	( b5.HLT_QuadPFJet98_83_71_15 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet103_88_75_15",&	( b5.HLT_QuadPFJet103_88_75_15 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet105_88_76_15",&	( b5.HLT_QuadPFJet105_88_76_15 ));
	_compiledTree->Branch("b5_HLT_QuadPFJet111_90_80_15",&	( b5.HLT_QuadPFJet111_90_80_15 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17",&	( b5.HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1",&	( b5.HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02",&	( b5.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2",&	( b5.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2 ));
	_compiledTree->Branch("b5_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4",&	( b5.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 ));
	_compiledTree->Branch("b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55",&	( b5.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55 ));
	_compiledTree->Branch("b5_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto",&	( b5.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto ));
	_compiledTree->Branch("b5_HLT_Mu8p5_IP3p5_part0",&	( b5.HLT_Mu8p5_IP3p5_part0 ));
	_compiledTree->Branch("b5_HLT_Mu8p5_IP3p5_part1",&	( b5.HLT_Mu8p5_IP3p5_part1 ));
	_compiledTree->Branch("b5_HLT_Mu8p5_IP3p5_part2",&	( b5.HLT_Mu8p5_IP3p5_part2 ));
	_compiledTree->Branch("b5_HLT_Mu8p5_IP3p5_part3",&	( b5.HLT_Mu8p5_IP3p5_part3 ));
	_compiledTree->Branch("b5_HLT_Mu8p5_IP3p5_part4",&	( b5.HLT_Mu8p5_IP3p5_part4 ));
	_compiledTree->Branch("b5_HLT_Mu8p5_IP3p5_part5",&	( b5.HLT_Mu8p5_IP3p5_part5 ));
	_compiledTree->Branch("b5_HLT_Mu10p5_IP3p5_part0",&	( b5.HLT_Mu10p5_IP3p5_part0 ));
	_compiledTree->Branch("b5_HLT_Mu10p5_IP3p5_part1",&	( b5.HLT_Mu10p5_IP3p5_part1 ));
	_compiledTree->Branch("b5_HLT_Mu10p5_IP3p5_part2",&	( b5.HLT_Mu10p5_IP3p5_part2 ));
	_compiledTree->Branch("b5_HLT_Mu10p5_IP3p5_part3",&	( b5.HLT_Mu10p5_IP3p5_part3 ));
	_compiledTree->Branch("b5_HLT_Mu10p5_IP3p5_part4",&	( b5.HLT_Mu10p5_IP3p5_part4 ));
	_compiledTree->Branch("b5_HLT_Mu10p5_IP3p5_part5",&	( b5.HLT_Mu10p5_IP3p5_part5 ));
	_compiledTree->Branch("b5_HLT_Mu9_IP6_part0",&	( b5.HLT_Mu9_IP6_part0 ));
	_compiledTree->Branch("b5_HLT_Mu9_IP6_part1",&	( b5.HLT_Mu9_IP6_part1 ));
	_compiledTree->Branch("b5_HLT_Mu9_IP6_part2",&	( b5.HLT_Mu9_IP6_part2 ));
	_compiledTree->Branch("b5_HLT_Mu9_IP6_part3",&	( b5.HLT_Mu9_IP6_part3 ));
	_compiledTree->Branch("b5_HLT_Mu9_IP6_part4",&	( b5.HLT_Mu9_IP6_part4 ));
	_compiledTree->Branch("b5_HLT_Mu9_IP6_part5",&	( b5.HLT_Mu9_IP6_part5 ));
	_compiledTree->Branch("b5_HLT_Mu8_IP3_part0",&	( b5.HLT_Mu8_IP3_part0 ));
	_compiledTree->Branch("b5_HLT_Mu8_IP3_part1",&	( b5.HLT_Mu8_IP3_part1 ));
	_compiledTree->Branch("b5_HLT_Mu8_IP3_part2",&	( b5.HLT_Mu8_IP3_part2 ));
	_compiledTree->Branch("b5_HLT_Mu8_IP3_part3",&	( b5.HLT_Mu8_IP3_part3 ));
	_compiledTree->Branch("b5_HLT_Mu8_IP3_part4",&	( b5.HLT_Mu8_IP3_part4 ));
	_compiledTree->Branch("b5_HLT_Mu8_IP3_part5",&	( b5.HLT_Mu8_IP3_part5 ));
	_compiledTree->Branch("b5_HLTriggerFinalPath",&	( b5.HLTriggerFinalPath ));
	_compiledTree->Branch("bG_run",&	( bG.run ));
	_compiledTree->Branch("bG_event",&	( bG.event ));
	_compiledTree->Branch("bG_lumis",&	( bG.lumis ));
	_compiledTree->Branch("bG_isData",&	( bG.isData ));
	_compiledTree->Branch("bG_nSC",&	( bG.nSC ));
	_compiledTree->Branch("bG_nscE", &bG_nscE);
	_compiledTree->Branch("bG_scE",	( bG_scE ),"bG_scE[bG_nscE]/F");
	_compiledTree->Branch("bG_nscEt", &bG_nscEt);
	_compiledTree->Branch("bG_scEt",	( bG_scEt ),"bG_scEt[bG_nscEt]/F");
	_compiledTree->Branch("bG_nscRawE", &bG_nscRawE);
	_compiledTree->Branch("bG_scRawE",	( bG_scRawE ),"bG_scRawE[bG_nscRawE]/F");
	_compiledTree->Branch("bG_nscEta", &bG_nscEta);
	_compiledTree->Branch("bG_scEta",	( bG_scEta ),"bG_scEta[bG_nscEta]/F");
	_compiledTree->Branch("bG_nscPhi", &bG_nscPhi);
	_compiledTree->Branch("bG_scPhi",	( bG_scPhi ),"bG_scPhi[bG_nscPhi]/F");
	_compiledTree->Branch("bG_nscX", &bG_nscX);
	_compiledTree->Branch("bG_scX",	( bG_scX ),"bG_scX[bG_nscX]/F");
	_compiledTree->Branch("bG_nscY", &bG_nscY);
	_compiledTree->Branch("bG_scY",	( bG_scY ),"bG_scY[bG_nscY]/F");
	_compiledTree->Branch("bG_nscZ", &bG_nscZ);
	_compiledTree->Branch("bG_scZ",	( bG_scZ ),"bG_scZ[bG_nscZ]/F");
	_compiledTree->Branch("bG_nscEtaWidth", &bG_nscEtaWidth);
	_compiledTree->Branch("bG_scEtaWidth",	( bG_scEtaWidth ),"bG_scEtaWidth[bG_nscEtaWidth]/F");
	_compiledTree->Branch("bG_nscPhiWidth", &bG_nscPhiWidth);
	_compiledTree->Branch("bG_scPhiWidth",	( bG_scPhiWidth ),"bG_scPhiWidth[bG_nscPhiWidth]/F");
	_compiledTree->Branch("bG_nscRawEt", &bG_nscRawEt);
	_compiledTree->Branch("bG_scRawEt",	( bG_scRawEt ),"bG_scRawEt[bG_nscRawEt]/F");
	_compiledTree->Branch("bG_nscMinDrWithGsfElectornSC_", &bG_nscMinDrWithGsfElectornSC_);
	_compiledTree->Branch("bG_scMinDrWithGsfElectornSC_",	( bG_scMinDrWithGsfElectornSC_ ),"bG_scMinDrWithGsfElectornSC_[bG_nscMinDrWithGsfElectornSC_]/F");
	_compiledTree->Branch("bG_nscFoundGsfMatch_", &bG_nscFoundGsfMatch_);
	_compiledTree->Branch("bG_scFoundGsfMatch_",	( bG_scFoundGsfMatch_ ),"bG_scFoundGsfMatch_[bG_nscFoundGsfMatch_]/O");
	_compiledTree->Branch("bG_nscE5x5", &bG_nscE5x5);
	_compiledTree->Branch("bG_scE5x5",	( bG_scE5x5 ),"bG_scE5x5[bG_nscE5x5]/F");
	_compiledTree->Branch("bG_nscE2x2Ratio", &bG_nscE2x2Ratio);
	_compiledTree->Branch("bG_scE2x2Ratio",	( bG_scE2x2Ratio ),"bG_scE2x2Ratio[bG_nscE2x2Ratio]/F");
	_compiledTree->Branch("bG_nscE3x3Ratio", &bG_nscE3x3Ratio);
	_compiledTree->Branch("bG_scE3x3Ratio",	( bG_scE3x3Ratio ),"bG_scE3x3Ratio[bG_nscE3x3Ratio]/F");
	_compiledTree->Branch("bG_nscEMaxRatio", &bG_nscEMaxRatio);
	_compiledTree->Branch("bG_scEMaxRatio",	( bG_scEMaxRatio ),"bG_scEMaxRatio[bG_nscEMaxRatio]/F");
	_compiledTree->Branch("bG_nscE2ndRatio", &bG_nscE2ndRatio);
	_compiledTree->Branch("bG_scE2ndRatio",	( bG_scE2ndRatio ),"bG_scE2ndRatio[bG_nscE2ndRatio]/F");
	_compiledTree->Branch("bG_nscETopRatio", &bG_nscETopRatio);
	_compiledTree->Branch("bG_scETopRatio",	( bG_scETopRatio ),"bG_scETopRatio[bG_nscETopRatio]/F");
	_compiledTree->Branch("bG_nscERightRatio", &bG_nscERightRatio);
	_compiledTree->Branch("bG_scERightRatio",	( bG_scERightRatio ),"bG_scERightRatio[bG_nscERightRatio]/F");
	_compiledTree->Branch("bG_nscEBottomRatio", &bG_nscEBottomRatio);
	_compiledTree->Branch("bG_scEBottomRatio",	( bG_scEBottomRatio ),"bG_scEBottomRatio[bG_nscEBottomRatio]/F");
	_compiledTree->Branch("bG_nscELeftRatio", &bG_nscELeftRatio);
	_compiledTree->Branch("bG_scELeftRatio",	( bG_scELeftRatio ),"bG_scELeftRatio[bG_nscELeftRatio]/F");
	_compiledTree->Branch("bG_nscE2x5MaxRatio", &bG_nscE2x5MaxRatio);
	_compiledTree->Branch("bG_scE2x5MaxRatio",	( bG_scE2x5MaxRatio ),"bG_scE2x5MaxRatio[bG_nscE2x5MaxRatio]/F");
	_compiledTree->Branch("bG_nscE2x5TopRatio", &bG_nscE2x5TopRatio);
	_compiledTree->Branch("bG_scE2x5TopRatio",	( bG_scE2x5TopRatio ),"bG_scE2x5TopRatio[bG_nscE2x5TopRatio]/F");
	_compiledTree->Branch("bG_nscE2x5RightRatio", &bG_nscE2x5RightRatio);
	_compiledTree->Branch("bG_scE2x5RightRatio",	( bG_scE2x5RightRatio ),"bG_scE2x5RightRatio[bG_nscE2x5RightRatio]/F");
	_compiledTree->Branch("bG_nscE2x5BottomRatio", &bG_nscE2x5BottomRatio);
	_compiledTree->Branch("bG_scE2x5BottomRatio",	( bG_scE2x5BottomRatio ),"bG_scE2x5BottomRatio[bG_nscE2x5BottomRatio]/F");
	_compiledTree->Branch("bG_nscE2x5LeftRatio", &bG_nscE2x5LeftRatio);
	_compiledTree->Branch("bG_scE2x5LeftRatio",	( bG_scE2x5LeftRatio ),"bG_scE2x5LeftRatio[bG_nscE2x5LeftRatio]/F");
	_compiledTree->Branch("bG_nscSwissCross", &bG_nscSwissCross);
	_compiledTree->Branch("bG_scSwissCross",	( bG_scSwissCross ),"bG_scSwissCross[bG_nscSwissCross]/F");
	_compiledTree->Branch("bG_nscR9", &bG_nscR9);
	_compiledTree->Branch("bG_scR9",	( bG_scR9 ),"bG_scR9[bG_nscR9]/F");
	_compiledTree->Branch("bG_nscSigmaIetaIeta", &bG_nscSigmaIetaIeta);
	_compiledTree->Branch("bG_scSigmaIetaIeta",	( bG_scSigmaIetaIeta ),"bG_scSigmaIetaIeta[bG_nscSigmaIetaIeta]/F");
	_compiledTree->Branch("bG_nscSigmaIetaIphi", &bG_nscSigmaIetaIphi);
	_compiledTree->Branch("bG_scSigmaIetaIphi",	( bG_scSigmaIetaIphi ),"bG_scSigmaIetaIphi[bG_nscSigmaIetaIphi]/F");
	_compiledTree->Branch("bG_nscSigmaIphiIphi", &bG_nscSigmaIphiIphi);
	_compiledTree->Branch("bG_scSigmaIphiIphi",	( bG_scSigmaIphiIphi ),"bG_scSigmaIphiIphi[bG_nscSigmaIphiIphi]/F");
	_compiledTree->Branch("bG_nscFull5x5_e5x5", &bG_nscFull5x5_e5x5);
	_compiledTree->Branch("bG_scFull5x5_e5x5",	( bG_scFull5x5_e5x5 ),"bG_scFull5x5_e5x5[bG_nscFull5x5_e5x5]/F");
	_compiledTree->Branch("bG_nscFull5x5_e2x2Ratio", &bG_nscFull5x5_e2x2Ratio);
	_compiledTree->Branch("bG_scFull5x5_e2x2Ratio",	( bG_scFull5x5_e2x2Ratio ),"bG_scFull5x5_e2x2Ratio[bG_nscFull5x5_e2x2Ratio]/F");
	_compiledTree->Branch("bG_nscFull5x5_e3x3Ratio", &bG_nscFull5x5_e3x3Ratio);
	_compiledTree->Branch("bG_scFull5x5_e3x3Ratio",	( bG_scFull5x5_e3x3Ratio ),"bG_scFull5x5_e3x3Ratio[bG_nscFull5x5_e3x3Ratio]/F");
	_compiledTree->Branch("bG_nscFull5x5_eMaxRatio", &bG_nscFull5x5_eMaxRatio);
	_compiledTree->Branch("bG_scFull5x5_eMaxRatio",	( bG_scFull5x5_eMaxRatio ),"bG_scFull5x5_eMaxRatio[bG_nscFull5x5_eMaxRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_e2ndRatio", &bG_nscFull5x5_e2ndRatio);
	_compiledTree->Branch("bG_scFull5x5_e2ndRatio",	( bG_scFull5x5_e2ndRatio ),"bG_scFull5x5_e2ndRatio[bG_nscFull5x5_e2ndRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_eTopRatio", &bG_nscFull5x5_eTopRatio);
	_compiledTree->Branch("bG_scFull5x5_eTopRatio",	( bG_scFull5x5_eTopRatio ),"bG_scFull5x5_eTopRatio[bG_nscFull5x5_eTopRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_eRightRatio", &bG_nscFull5x5_eRightRatio);
	_compiledTree->Branch("bG_scFull5x5_eRightRatio",	( bG_scFull5x5_eRightRatio ),"bG_scFull5x5_eRightRatio[bG_nscFull5x5_eRightRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_eBottomRatio", &bG_nscFull5x5_eBottomRatio);
	_compiledTree->Branch("bG_scFull5x5_eBottomRatio",	( bG_scFull5x5_eBottomRatio ),"bG_scFull5x5_eBottomRatio[bG_nscFull5x5_eBottomRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_eLeftRatio", &bG_nscFull5x5_eLeftRatio);
	_compiledTree->Branch("bG_scFull5x5_eLeftRatio",	( bG_scFull5x5_eLeftRatio ),"bG_scFull5x5_eLeftRatio[bG_nscFull5x5_eLeftRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_e2x5MaxRatio", &bG_nscFull5x5_e2x5MaxRatio);
	_compiledTree->Branch("bG_scFull5x5_e2x5MaxRatio",	( bG_scFull5x5_e2x5MaxRatio ),"bG_scFull5x5_e2x5MaxRatio[bG_nscFull5x5_e2x5MaxRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_e2x5TopRatio", &bG_nscFull5x5_e2x5TopRatio);
	_compiledTree->Branch("bG_scFull5x5_e2x5TopRatio",	( bG_scFull5x5_e2x5TopRatio ),"bG_scFull5x5_e2x5TopRatio[bG_nscFull5x5_e2x5TopRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_e2x5RightRatio", &bG_nscFull5x5_e2x5RightRatio);
	_compiledTree->Branch("bG_scFull5x5_e2x5RightRatio",	( bG_scFull5x5_e2x5RightRatio ),"bG_scFull5x5_e2x5RightRatio[bG_nscFull5x5_e2x5RightRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_e2x5BottomRatio", &bG_nscFull5x5_e2x5BottomRatio);
	_compiledTree->Branch("bG_scFull5x5_e2x5BottomRatio",	( bG_scFull5x5_e2x5BottomRatio ),"bG_scFull5x5_e2x5BottomRatio[bG_nscFull5x5_e2x5BottomRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_e2x5LeftRatio", &bG_nscFull5x5_e2x5LeftRatio);
	_compiledTree->Branch("bG_scFull5x5_e2x5LeftRatio",	( bG_scFull5x5_e2x5LeftRatio ),"bG_scFull5x5_e2x5LeftRatio[bG_nscFull5x5_e2x5LeftRatio]/F");
	_compiledTree->Branch("bG_nscFull5x5_swissCross", &bG_nscFull5x5_swissCross);
	_compiledTree->Branch("bG_scFull5x5_swissCross",	( bG_scFull5x5_swissCross ),"bG_scFull5x5_swissCross[bG_nscFull5x5_swissCross]/F");
	_compiledTree->Branch("bG_nscFull5x5_r9", &bG_nscFull5x5_r9);
	_compiledTree->Branch("bG_scFull5x5_r9",	( bG_scFull5x5_r9 ),"bG_scFull5x5_r9[bG_nscFull5x5_r9]/F");
	_compiledTree->Branch("bG_nscFull5x5_sigmaIetaIeta", &bG_nscFull5x5_sigmaIetaIeta);
	_compiledTree->Branch("bG_scFull5x5_sigmaIetaIeta",	( bG_scFull5x5_sigmaIetaIeta ),"bG_scFull5x5_sigmaIetaIeta[bG_nscFull5x5_sigmaIetaIeta]/F");
	_compiledTree->Branch("bG_nscFull5x5_sigmaIetaIphi", &bG_nscFull5x5_sigmaIetaIphi);
	_compiledTree->Branch("bG_scFull5x5_sigmaIetaIphi",	( bG_scFull5x5_sigmaIetaIphi ),"bG_scFull5x5_sigmaIetaIphi[bG_nscFull5x5_sigmaIetaIphi]/F");
	_compiledTree->Branch("bG_nscFull5x5_sigmaIphiIphi", &bG_nscFull5x5_sigmaIphiIphi);
	_compiledTree->Branch("bG_scFull5x5_sigmaIphiIphi",	( bG_scFull5x5_sigmaIphiIphi ),"bG_scFull5x5_sigmaIphiIphi[bG_nscFull5x5_sigmaIphiIphi]/F");
	_compiledTree->Branch("bG_nscNHcalRecHitInDIEta5IPhi5", &bG_nscNHcalRecHitInDIEta5IPhi5);
	_compiledTree->Branch("bG_scNHcalRecHitInDIEta5IPhi5",	( bG_scNHcalRecHitInDIEta5IPhi5 ),"bG_scNHcalRecHitInDIEta5IPhi5[bG_nscNHcalRecHitInDIEta5IPhi5]/F");
	_compiledTree->Branch("bG_nscEFromHcalRecHitInDIEta5IPhi5", &bG_nscEFromHcalRecHitInDIEta5IPhi5);
	_compiledTree->Branch("bG_scEFromHcalRecHitInDIEta5IPhi5",	( bG_scEFromHcalRecHitInDIEta5IPhi5 ),"bG_scEFromHcalRecHitInDIEta5IPhi5[bG_nscEFromHcalRecHitInDIEta5IPhi5]/F");
	_compiledTree->Branch("bG_nscNHcalRecHitInDIEta2IPhi2", &bG_nscNHcalRecHitInDIEta2IPhi2);
	_compiledTree->Branch("bG_scNHcalRecHitInDIEta2IPhi2",	( bG_scNHcalRecHitInDIEta2IPhi2 ),"bG_scNHcalRecHitInDIEta2IPhi2[bG_nscNHcalRecHitInDIEta2IPhi2]/F");
	_compiledTree->Branch("bG_nscEFromHcalRecHitInDIEta2IPhi2", &bG_nscEFromHcalRecHitInDIEta2IPhi2);
	_compiledTree->Branch("bG_scEFromHcalRecHitInDIEta2IPhi2",	( bG_scEFromHcalRecHitInDIEta2IPhi2 ),"bG_scEFromHcalRecHitInDIEta2IPhi2[bG_nscEFromHcalRecHitInDIEta2IPhi2]/F");
	_compiledTree->Branch("bG_nscPFChIso1", &bG_nscPFChIso1);
	_compiledTree->Branch("bG_scPFChIso1",	( bG_scPFChIso1 ),"bG_scPFChIso1[bG_nscPFChIso1]/F");
	_compiledTree->Branch("bG_nscPFChIso2", &bG_nscPFChIso2);
	_compiledTree->Branch("bG_scPFChIso2",	( bG_scPFChIso2 ),"bG_scPFChIso2[bG_nscPFChIso2]/F");
	_compiledTree->Branch("bG_nscPFChIso3", &bG_nscPFChIso3);
	_compiledTree->Branch("bG_scPFChIso3",	( bG_scPFChIso3 ),"bG_scPFChIso3[bG_nscPFChIso3]/F");
	_compiledTree->Branch("bG_nscPFChIso4", &bG_nscPFChIso4);
	_compiledTree->Branch("bG_scPFChIso4",	( bG_scPFChIso4 ),"bG_scPFChIso4[bG_nscPFChIso4]/F");
	_compiledTree->Branch("bG_nscPFChIso5", &bG_nscPFChIso5);
	_compiledTree->Branch("bG_scPFChIso5",	( bG_scPFChIso5 ),"bG_scPFChIso5[bG_nscPFChIso5]/F");
	_compiledTree->Branch("bG_nscPFPhoIso1", &bG_nscPFPhoIso1);
	_compiledTree->Branch("bG_scPFPhoIso1",	( bG_scPFPhoIso1 ),"bG_scPFPhoIso1[bG_nscPFPhoIso1]/F");
	_compiledTree->Branch("bG_nscPFPhoIso2", &bG_nscPFPhoIso2);
	_compiledTree->Branch("bG_scPFPhoIso2",	( bG_scPFPhoIso2 ),"bG_scPFPhoIso2[bG_nscPFPhoIso2]/F");
	_compiledTree->Branch("bG_nscPFPhoIso3", &bG_nscPFPhoIso3);
	_compiledTree->Branch("bG_scPFPhoIso3",	( bG_scPFPhoIso3 ),"bG_scPFPhoIso3[bG_nscPFPhoIso3]/F");
	_compiledTree->Branch("bG_nscPFPhoIso4", &bG_nscPFPhoIso4);
	_compiledTree->Branch("bG_scPFPhoIso4",	( bG_scPFPhoIso4 ),"bG_scPFPhoIso4[bG_nscPFPhoIso4]/F");
	_compiledTree->Branch("bG_nscPFPhoIso5", &bG_nscPFPhoIso5);
	_compiledTree->Branch("bG_scPFPhoIso5",	( bG_scPFPhoIso5 ),"bG_scPFPhoIso5[bG_nscPFPhoIso5]/F");
	_compiledTree->Branch("bG_nscPFNeuIso1", &bG_nscPFNeuIso1);
	_compiledTree->Branch("bG_scPFNeuIso1",	( bG_scPFNeuIso1 ),"bG_scPFNeuIso1[bG_nscPFNeuIso1]/F");
	_compiledTree->Branch("bG_nscPFNeuIso2", &bG_nscPFNeuIso2);
	_compiledTree->Branch("bG_scPFNeuIso2",	( bG_scPFNeuIso2 ),"bG_scPFNeuIso2[bG_nscPFNeuIso2]/F");
	_compiledTree->Branch("bG_nscPFNeuIso3", &bG_nscPFNeuIso3);
	_compiledTree->Branch("bG_scPFNeuIso3",	( bG_scPFNeuIso3 ),"bG_scPFNeuIso3[bG_nscPFNeuIso3]/F");
	_compiledTree->Branch("bG_nscPFNeuIso4", &bG_nscPFNeuIso4);
	_compiledTree->Branch("bG_scPFNeuIso4",	( bG_scPFNeuIso4 ),"bG_scPFNeuIso4[bG_nscPFNeuIso4]/F");
	_compiledTree->Branch("bG_nscPFNeuIso5", &bG_nscPFNeuIso5);
	_compiledTree->Branch("bG_scPFNeuIso5",	( bG_scPFNeuIso5 ),"bG_scPFNeuIso5[bG_nscPFNeuIso5]/F");
	_compiledTree->Branch("bG_nPrimaryVertex",&	( bG.nPrimaryVertex ));
	_compiledTree->Branch("bG_primaryVertex_isFake",	( bG.primaryVertex_isFake ),"bG_primaryVertex_isFake[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_x",	( bG.primaryVertex_x ),"bG_primaryVertex_x[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_y",	( bG.primaryVertex_y ),"bG_primaryVertex_y[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_z",	( bG.primaryVertex_z ),"bG_primaryVertex_z[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_t",	( bG.primaryVertex_t ),"bG_primaryVertex_t[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_covXX",	( bG.primaryVertex_covXX ),"bG_primaryVertex_covXX[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_covXY",	( bG.primaryVertex_covXY ),"bG_primaryVertex_covXY[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_covXZ",	( bG.primaryVertex_covXZ ),"bG_primaryVertex_covXZ[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_covYY",	( bG.primaryVertex_covYY ),"bG_primaryVertex_covYY[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_covYZ",	( bG.primaryVertex_covYZ ),"bG_primaryVertex_covYZ[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_covZZ",	( bG.primaryVertex_covZZ ),"bG_primaryVertex_covZZ[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_x_error",	( bG.primaryVertex_x_error ),"bG_primaryVertex_x_error[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_y_error",	( bG.primaryVertex_y_error ),"bG_primaryVertex_y_error[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_z_error",	( bG.primaryVertex_z_error ),"bG_primaryVertex_z_error[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_t_error",	( bG.primaryVertex_t_error ),"bG_primaryVertex_t_error[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_ntracks",	( bG.primaryVertex_ntracks ),"bG_primaryVertex_ntracks[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_ndof",	( bG.primaryVertex_ndof ),"bG_primaryVertex_ndof[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_chi2",	( bG.primaryVertex_chi2 ),"bG_primaryVertex_chi2[bG_nPrimaryVertex]/F");
	_compiledTree->Branch("bG_primaryVertex_normalizedChi2",	( bG.primaryVertex_normalizedChi2 ),"bG_primaryVertex_normalizedChi2[bG_nPrimaryVertex]/F");

}
