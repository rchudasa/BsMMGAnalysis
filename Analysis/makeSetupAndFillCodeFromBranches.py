#!/usr/bin/env python3

def writeSetupAssignFiles(txt,setupFile,fillFile,nCands,candMax):
    setupFile.write("\n")
    setupFile.write("\tauto idx = offset ;\n ")
    fillFile.write( "\tInt_t idx = 0 ;\n ")
    
    for l in txt:
        brname=l[:-1]
        setupFile.write("\n")
        setupFile.write('\tcandidateMap["'+brname+'"] = idx; \n')
        setupFile.write('\toutTree->Branch("'+brname+'" , &(storageArray[idx] ) , "'+brname+'['+nCands+']/D"); \n')
       # setupFile.write('\tstd::cout<< "'+brname+'"<<" idx : "<<idx<<"  candidateMap['+brname+'] = "<<candidateMap["'+brname+'"]<<"\\n";\n')
        setupFile.write('\tidx+='+candMax+'; \n')
    
        fillFile.write('\tstorageArray[ offset + candidateMap["'+brname+'"]] = ntupleRawTree.' + brname+ '[candIdx] ; \n')
        #fillFile.write('\tidx = candidateMap["'+brname+'"] ; \n')
        #fillFile.write('\tstd::cout<< "'+brname+'"<<" idx : "<<idx<<"  val = "<<storageArray[idx]<<"\\n";\n\n')
    setupFile.write("\n\treturn idx;")


f=open("branches/branches.br.dimu")
txt=f.readlines()
f.close()
setupFile=open("setupDimu.cc",'w')
fillFile =open("assignDimu.cc",'w')
nCands="nDiMuCandidates"
candMax="size"
writeSetupAssignFiles(txt,setupFile,fillFile,nCands,candMax)
fillFile.close()
setupFile.close()

f=open("branches/branches.br.sc")
txt=f.readlines()
f.close()
setupFile=open("setupSc.cc",'w')
fillFile =open("assignSc.cc",'w')
nCands="nSCPhotons"
candMax="size"
writeSetupAssignFiles(txt,setupFile,fillFile,nCands,candMax)
fillFile.close()
setupFile.close()


f=open("branches/branches.br.muon")
txt=f.readlines()
f.close()
setupFile=open("setupMuon.cc",'w')
fillFile =open("assignMuon.cc",'w')
nCands="nMuons"
candMax="size"
writeSetupAssignFiles(txt,setupFile,fillFile,nCands,candMax)
fillFile.close()
setupFile.close()


