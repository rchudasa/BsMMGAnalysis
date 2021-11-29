
	auto idx = offset ;
 
	candidateMap["b5_mm_kal_slxy"] = idx; 
	outTree->Branch("b5_mm_kal_slxy" , &(storageArray[idx] ) , "b5_mm_kal_slxy[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kal_vtx_prob"] = idx; 
	outTree->Branch("b5_mm_kal_vtx_prob" , &(storageArray[idx] ) , "b5_mm_kal_vtx_prob[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_alpha"] = idx; 
	outTree->Branch("b5_mm_kin_alpha" , &(storageArray[idx] ) , "b5_mm_kin_alpha[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_eta"] = idx; 
	outTree->Branch("b5_mm_kin_eta" , &(storageArray[idx] ) , "b5_mm_kin_eta[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_l3d"] = idx; 
	outTree->Branch("b5_mm_kin_l3d" , &(storageArray[idx] ) , "b5_mm_kin_l3d[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_lxy"] = idx; 
	outTree->Branch("b5_mm_kin_lxy" , &(storageArray[idx] ) , "b5_mm_kin_lxy[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_mass"] = idx; 
	outTree->Branch("b5_mm_kin_mass" , &(storageArray[idx] ) , "b5_mm_kin_mass[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_massErr"] = idx; 
	outTree->Branch("b5_mm_kin_massErr" , &(storageArray[idx] ) , "b5_mm_kin_massErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_mu1eta"] = idx; 
	outTree->Branch("b5_mm_kin_mu1eta" , &(storageArray[idx] ) , "b5_mm_kin_mu1eta[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_mu1phi"] = idx; 
	outTree->Branch("b5_mm_kin_mu1phi" , &(storageArray[idx] ) , "b5_mm_kin_mu1phi[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_mu1pt"] = idx; 
	outTree->Branch("b5_mm_kin_mu1pt" , &(storageArray[idx] ) , "b5_mm_kin_mu1pt[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_mu2eta"] = idx; 
	outTree->Branch("b5_mm_kin_mu2eta" , &(storageArray[idx] ) , "b5_mm_kin_mu2eta[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_mu2phi"] = idx; 
	outTree->Branch("b5_mm_kin_mu2phi" , &(storageArray[idx] ) , "b5_mm_kin_mu2phi[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_mu2pt"] = idx; 
	outTree->Branch("b5_mm_kin_mu2pt" , &(storageArray[idx] ) , "b5_mm_kin_mu2pt[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_phi"] = idx; 
	outTree->Branch("b5_mm_kin_phi" , &(storageArray[idx] ) , "b5_mm_kin_phi[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pt"] = idx; 
	outTree->Branch("b5_mm_kin_pt" , &(storageArray[idx] ) , "b5_mm_kin_pt[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pv2ip"] = idx; 
	outTree->Branch("b5_mm_kin_pv2ip" , &(storageArray[idx] ) , "b5_mm_kin_pv2ip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pv2ipErr"] = idx; 
	outTree->Branch("b5_mm_kin_pv2ipErr" , &(storageArray[idx] ) , "b5_mm_kin_pv2ipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pv2lip"] = idx; 
	outTree->Branch("b5_mm_kin_pv2lip" , &(storageArray[idx] ) , "b5_mm_kin_pv2lip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pv2lipErr"] = idx; 
	outTree->Branch("b5_mm_kin_pv2lipErr" , &(storageArray[idx] ) , "b5_mm_kin_pv2lipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pv2lipSig"] = idx; 
	outTree->Branch("b5_mm_kin_pv2lipSig" , &(storageArray[idx] ) , "b5_mm_kin_pv2lipSig[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pv_z"] = idx; 
	outTree->Branch("b5_mm_kin_pv_z" , &(storageArray[idx] ) , "b5_mm_kin_pv_z[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pv_zErr"] = idx; 
	outTree->Branch("b5_mm_kin_pv_zErr" , &(storageArray[idx] ) , "b5_mm_kin_pv_zErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pvip"] = idx; 
	outTree->Branch("b5_mm_kin_pvip" , &(storageArray[idx] ) , "b5_mm_kin_pvip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pvipErr"] = idx; 
	outTree->Branch("b5_mm_kin_pvipErr" , &(storageArray[idx] ) , "b5_mm_kin_pvipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pvlip"] = idx; 
	outTree->Branch("b5_mm_kin_pvlip" , &(storageArray[idx] ) , "b5_mm_kin_pvlip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pvlipErr"] = idx; 
	outTree->Branch("b5_mm_kin_pvlipErr" , &(storageArray[idx] ) , "b5_mm_kin_pvlipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_pvlipSig"] = idx; 
	outTree->Branch("b5_mm_kin_pvlipSig" , &(storageArray[idx] ) , "b5_mm_kin_pvlipSig[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_sl3d"] = idx; 
	outTree->Branch("b5_mm_kin_sl3d" , &(storageArray[idx] ) , "b5_mm_kin_sl3d[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_slxy"] = idx; 
	outTree->Branch("b5_mm_kin_slxy" , &(storageArray[idx] ) , "b5_mm_kin_slxy[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_tau"] = idx; 
	outTree->Branch("b5_mm_kin_tau" , &(storageArray[idx] ) , "b5_mm_kin_tau[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_taue"] = idx; 
	outTree->Branch("b5_mm_kin_taue" , &(storageArray[idx] ) , "b5_mm_kin_taue[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_tauxy"] = idx; 
	outTree->Branch("b5_mm_kin_tauxy" , &(storageArray[idx] ) , "b5_mm_kin_tauxy[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_tauxye"] = idx; 
	outTree->Branch("b5_mm_kin_tauxye" , &(storageArray[idx] ) , "b5_mm_kin_tauxye[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_chi2dof"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_chi2dof" , &(storageArray[idx] ) , "b5_mm_kin_vtx_chi2dof[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_prob"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_prob" , &(storageArray[idx] ) , "b5_mm_kin_vtx_prob[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_x"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_x" , &(storageArray[idx] ) , "b5_mm_kin_vtx_x[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_xErr"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_xErr" , &(storageArray[idx] ) , "b5_mm_kin_vtx_xErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_y"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_y" , &(storageArray[idx] ) , "b5_mm_kin_vtx_y[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_yErr"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_yErr" , &(storageArray[idx] ) , "b5_mm_kin_vtx_yErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_z"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_z" , &(storageArray[idx] ) , "b5_mm_kin_vtx_z[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_vtx_zErr"] = idx; 
	outTree->Branch("b5_mm_kin_vtx_zErr" , &(storageArray[idx] ) , "b5_mm_kin_vtx_zErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_alpha"] = idx; 
	outTree->Branch("b5_mm_kinpc_alpha" , &(storageArray[idx] ) , "b5_mm_kinpc_alpha[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_eta"] = idx; 
	outTree->Branch("b5_mm_kinpc_eta" , &(storageArray[idx] ) , "b5_mm_kinpc_eta[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_l3d"] = idx; 
	outTree->Branch("b5_mm_kinpc_l3d" , &(storageArray[idx] ) , "b5_mm_kinpc_l3d[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_lxy"] = idx; 
	outTree->Branch("b5_mm_kinpc_lxy" , &(storageArray[idx] ) , "b5_mm_kinpc_lxy[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_mass"] = idx; 
	outTree->Branch("b5_mm_kinpc_mass" , &(storageArray[idx] ) , "b5_mm_kinpc_mass[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_massErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_massErr" , &(storageArray[idx] ) , "b5_mm_kinpc_massErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_phi"] = idx; 
	outTree->Branch("b5_mm_kinpc_phi" , &(storageArray[idx] ) , "b5_mm_kinpc_phi[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pt"] = idx; 
	outTree->Branch("b5_mm_kinpc_pt" , &(storageArray[idx] ) , "b5_mm_kinpc_pt[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pv2ip"] = idx; 
	outTree->Branch("b5_mm_kinpc_pv2ip" , &(storageArray[idx] ) , "b5_mm_kinpc_pv2ip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pv2ipErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_pv2ipErr" , &(storageArray[idx] ) , "b5_mm_kinpc_pv2ipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pv2lip"] = idx; 
	outTree->Branch("b5_mm_kinpc_pv2lip" , &(storageArray[idx] ) , "b5_mm_kinpc_pv2lip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pv2lipErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_pv2lipErr" , &(storageArray[idx] ) , "b5_mm_kinpc_pv2lipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pv2lipSig"] = idx; 
	outTree->Branch("b5_mm_kinpc_pv2lipSig" , &(storageArray[idx] ) , "b5_mm_kinpc_pv2lipSig[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pv_z"] = idx; 
	outTree->Branch("b5_mm_kinpc_pv_z" , &(storageArray[idx] ) , "b5_mm_kinpc_pv_z[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pv_zErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_pv_zErr" , &(storageArray[idx] ) , "b5_mm_kinpc_pv_zErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pvip"] = idx; 
	outTree->Branch("b5_mm_kinpc_pvip" , &(storageArray[idx] ) , "b5_mm_kinpc_pvip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pvipErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_pvipErr" , &(storageArray[idx] ) , "b5_mm_kinpc_pvipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pvlip"] = idx; 
	outTree->Branch("b5_mm_kinpc_pvlip" , &(storageArray[idx] ) , "b5_mm_kinpc_pvlip[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pvlipErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_pvlipErr" , &(storageArray[idx] ) , "b5_mm_kinpc_pvlipErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_pvlipSig"] = idx; 
	outTree->Branch("b5_mm_kinpc_pvlipSig" , &(storageArray[idx] ) , "b5_mm_kinpc_pvlipSig[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_sl3d"] = idx; 
	outTree->Branch("b5_mm_kinpc_sl3d" , &(storageArray[idx] ) , "b5_mm_kinpc_sl3d[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_slxy"] = idx; 
	outTree->Branch("b5_mm_kinpc_slxy" , &(storageArray[idx] ) , "b5_mm_kinpc_slxy[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_tau"] = idx; 
	outTree->Branch("b5_mm_kinpc_tau" , &(storageArray[idx] ) , "b5_mm_kinpc_tau[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_taue"] = idx; 
	outTree->Branch("b5_mm_kinpc_taue" , &(storageArray[idx] ) , "b5_mm_kinpc_taue[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_tauxy"] = idx; 
	outTree->Branch("b5_mm_kinpc_tauxy" , &(storageArray[idx] ) , "b5_mm_kinpc_tauxy[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_tauxye"] = idx; 
	outTree->Branch("b5_mm_kinpc_tauxye" , &(storageArray[idx] ) , "b5_mm_kinpc_tauxye[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_chi2dof"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_chi2dof" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_chi2dof[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_prob"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_prob" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_prob[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_x"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_x" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_x[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_xErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_xErr" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_xErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_y"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_y" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_y[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_yErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_yErr" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_yErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_z"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_z" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_z[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_vtx_zErr"] = idx; 
	outTree->Branch("b5_mm_kinpc_vtx_zErr" , &(storageArray[idx] ) , "b5_mm_kinpc_vtx_zErr[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_m1iso"] = idx; 
	outTree->Branch("b5_mm_m1iso" , &(storageArray[idx] ) , "b5_mm_m1iso[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_m2iso"] = idx; 
	outTree->Branch("b5_mm_m2iso" , &(storageArray[idx] ) , "b5_mm_m2iso[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mass"] = idx; 
	outTree->Branch("b5_mm_mass" , &(storageArray[idx] ) , "b5_mm_mass[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu1_eta"] = idx; 
	outTree->Branch("b5_mm_mu1_eta" , &(storageArray[idx] ) , "b5_mm_mu1_eta[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu1_phi"] = idx; 
	outTree->Branch("b5_mm_mu1_phi" , &(storageArray[idx] ) , "b5_mm_mu1_phi[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu1_pt"] = idx; 
	outTree->Branch("b5_mm_mu1_pt" , &(storageArray[idx] ) , "b5_mm_mu1_pt[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu2_eta"] = idx; 
	outTree->Branch("b5_mm_mu2_eta" , &(storageArray[idx] ) , "b5_mm_mu2_eta[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu2_phi"] = idx; 
	outTree->Branch("b5_mm_mu2_phi" , &(storageArray[idx] ) , "b5_mm_mu2_phi[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu2_pt"] = idx; 
	outTree->Branch("b5_mm_mu2_pt" , &(storageArray[idx] ) , "b5_mm_mu2_pt[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mva"] = idx; 
	outTree->Branch("b5_mm_mva" , &(storageArray[idx] ) , "b5_mm_mva[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_otherVtxMaxProb"] = idx; 
	outTree->Branch("b5_mm_otherVtxMaxProb" , &(storageArray[idx] ) , "b5_mm_otherVtxMaxProb[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_otherVtxMaxProb1"] = idx; 
	outTree->Branch("b5_mm_otherVtxMaxProb1" , &(storageArray[idx] ) , "b5_mm_otherVtxMaxProb1[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_otherVtxMaxProb2"] = idx; 
	outTree->Branch("b5_mm_otherVtxMaxProb2" , &(storageArray[idx] ) , "b5_mm_otherVtxMaxProb2[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_closetrk"] = idx; 
	outTree->Branch("b5_mm_closetrk" , &(storageArray[idx] ) , "b5_mm_closetrk[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_closetrks1"] = idx; 
	outTree->Branch("b5_mm_closetrks1" , &(storageArray[idx] ) , "b5_mm_closetrks1[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_closetrks2"] = idx; 
	outTree->Branch("b5_mm_closetrks2" , &(storageArray[idx] ) , "b5_mm_closetrks2[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_closetrks3"] = idx; 
	outTree->Branch("b5_mm_closetrks3" , &(storageArray[idx] ) , "b5_mm_closetrks3[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kal_valid"] = idx; 
	outTree->Branch("b5_mm_kal_valid" , &(storageArray[idx] ) , "b5_mm_kal_valid[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kin_valid"] = idx; 
	outTree->Branch("b5_mm_kin_valid" , &(storageArray[idx] ) , "b5_mm_kin_valid[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_kinpc_valid"] = idx; 
	outTree->Branch("b5_mm_kinpc_valid" , &(storageArray[idx] ) , "b5_mm_kinpc_valid[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu1_index"] = idx; 
	outTree->Branch("b5_mm_mu1_index" , &(storageArray[idx] ) , "b5_mm_mu1_index[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu1_pdgId"] = idx; 
	outTree->Branch("b5_mm_mu1_pdgId" , &(storageArray[idx] ) , "b5_mm_mu1_pdgId[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu2_index"] = idx; 
	outTree->Branch("b5_mm_mu2_index" , &(storageArray[idx] ) , "b5_mm_mu2_index[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_mu2_pdgId"] = idx; 
	outTree->Branch("b5_mm_mu2_pdgId" , &(storageArray[idx] ) , "b5_mm_mu2_pdgId[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_nBMTrks"] = idx; 
	outTree->Branch("b5_mm_nBMTrks" , &(storageArray[idx] ) , "b5_mm_nBMTrks[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_nDisTrks"] = idx; 
	outTree->Branch("b5_mm_nDisTrks" , &(storageArray[idx] ) , "b5_mm_nDisTrks[nDiMuCandidates]/D"); 
	idx+=size; 

	candidateMap["b5_mm_nTrks"] = idx; 
	outTree->Branch("b5_mm_nTrks" , &(storageArray[idx] ) , "b5_mm_nTrks[nDiMuCandidates]/D"); 
	idx+=size; 

	return idx;