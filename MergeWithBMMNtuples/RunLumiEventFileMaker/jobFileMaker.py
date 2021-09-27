#!/usr/bin/env python3
import os

NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=25
PU_NEVENTS=-1
destination='/afs/cern.ch/work/a/athachay/private/bs2mumug/run2studies/bmm5/CMSSW_10_6_19_patch2/src/runLumiEvents/'

FileSource ="CharmoniumRun2018A12Nov2019_UL2018_filelist.txt"

pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']
xrdRedirector="root://cms-xrd-global.cern.ch/"

Fnames=open(FileSource,'r')
sourceFileList=Fnames.readlines()
Fnames.close()


configurationTxt="\
#FLIST_BEG\n\
@@FNAMES\n\
#FLIST_END\n\
\n\
#PARAMS_BEG\n\
OutputPrefix=\n\
OutputFile=bmm5EvntSelection_@@IDX.root\n\
OutputFileTxt=bmm5EvntSelection_@@IDX.txt\n\
#PARAMS_END\n\
\n\
#EOF\n\
"

condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
+JobFlavour = \"longlunch\"\n\
"
condorScript=open('jobSubmit.sub','w')
condorScript.write(condorScriptString)



runScriptTxt="\
#!/bin/bash\n\
set -x\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME="+HOME+"\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd "+pwd+"/@@DIRNAME \n\
eval `scramv1 runtime -sh`\n\
cp  /afs/cern.ch/work/a/athachay/private/bs2mumug/run2studies/bmm5/CMSSW_10_6_19_patch2/src/Bmm5Ntuple* .\n\
cp  /afs/cern.ch/work/a/athachay/private/bs2mumug/run2studies/bmm5/CMSSW_10_6_19_patch2/src/getEventListFromBMM5.cc .\n\
root -b -q 'getEventListFromBMM5.cc(\"@@CFGFILENAME\")'\n\
mv *.root "+destination+"\n\
mv *.txt "+destination+"\n\
"

if not os.path.exists('Jobs'):
    os.system('mkdir Jobs')
print("Making ",NJOBS," Jobs ")

for ii in range(NJOBS):
    i=ii+ZERO_OFFSET
    if len(sourceFileList)<FILES_PER_JOB:
       print("fname count less than required .. stoping ")
       FILES_PER_JOB=len(sourceFileList)
    if len(sourceFileList) ==0:
       break 
    dirName= 'Jobs/Job_'+str(i)
    print(i," Job Made")
    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir '+dirName)
    
    cfgFileName='bmm5Selection_'+str(i)+'.cfg'
    cfgFile=open(dirName+'/'+cfgFileName,'w')
    tmp=""
    for j in range(FILES_PER_JOB):
      tmp+= sourceFileList.pop(0)[:-1]+"\n"
    tmp=configurationTxt.replace("@@FNAMES",tmp[:-1])
    tmp=tmp.replace("@@IDX",str(i))
    cfgFile.write(tmp)
    cfgFile.close()   
    
    runScriptName=dirName+'/run'+str(i)+'.sh'
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")

print(" Number of files left : ", len(sourceFileList) )
condorScript.close()
