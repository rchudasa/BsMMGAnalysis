
f=open("branches.br.dimu")
txt=f.readlines()
f.close()

setupFile=open("setupDimu.cc",'w')
fillFile =open("assignDimu.cc",'w')

i=0
setupFile.write("\n")
setupFile.write("\tauto idx = offset ;\n ")
fillFile.write( "\tInt_t idx = 0 ;\n ")

for l in txt:
    brname=l[:-1]
    setupFile.write("\n")
    setupFile.write('\tcandidateMap["'+brname+'"] = idx; \n')
    setupFile.write('\toutTree->Branch("'+brname+'" , &(storageArray[idx] ) ); \n')
   # setupFile.write('\tstd::cout<< "'+brname+'"<<" idx : "<<idx<<"  candidateMap['+brname+'] = "<<candidateMap["'+brname+'"]<<"\\n";\n')
    setupFile.write('\tidx++; \n')

    fillFile.write('\tidx = candidateMap["'+brname+'"] ; \n')
    fillFile.write('\tstorageArray[idx] = ntupleRawTree.' + brname + '[candIdx] ; \n')
  #  fillFile.write('\tstd::cout<< "'+brname+'"<<" idx : "<<idx<<"  val = "<<storageArray[idx]<<"\\n";\n\n')

setupFile.write("\t return idx;")
fillFile.close()
setupFile.close()

f=open("branches.br.sc")
txt=f.readlines()
f.close()

setupFile=open("setupSc.cc",'w')
fillFile =open("assignSc.cc",'w')

i=0
setupFile.write("\n")
setupFile.write("\tauto idx = offset ;\n ")
fillFile.write( "\tInt_t idx = 0 ;\n ")

for l in txt:
    brname=l[:-1]
    setupFile.write("\n")
    setupFile.write('\tcandidateMap["'+brname+'"] = idx; \n')
    setupFile.write('\toutTree->Branch("'+brname+'" , &(storageArray[idx] ) ); \n')
    #setupFile.write('\tstd::cout<< "'+brname+'"<<" idx : "<<idx<<"  candidateMap['+brname+'] = "<<candidateMap["'+brname+'"]<<"\\n";\n')
    setupFile.write('\tidx++; \n')

    fillFile.write('\tidx = candidateMap["'+brname+'"] ; \n')
    fillFile.write('\tstorageArray[idx] = ntupleRawTree.' + brname + '[candIdx] ; \n')
    #fillFile.write('\tstd::cout<< "'+brname+'"<<" idx : "<<idx<<"  val = "<<storageArray[idx]<<"\\n";\n')

setupFile.write("\n\treturn idx;")
fillFile.close()
setupFile.close()

