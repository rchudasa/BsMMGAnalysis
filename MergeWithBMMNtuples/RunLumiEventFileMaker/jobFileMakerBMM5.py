#!/usr/bin/env python3
import os

NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=2
destination='/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/CMSSW_10_6_19_patch2/src/BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/runLumiList/'

FileSource ="bmm5FileList.txt"

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
DoSelection=0\n\
MaxEvents=-1\n\
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
condorScript=open('subCondorBMM5.sub','w')
condorScript.write(condorScriptString)



runScriptTxt="\
#!/bin/bash\n\
set -x\n\
mv @@RUNSCRIPT @@RUNSCRIPT.busy \n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME="+HOME+"\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd @@DIRNAME \n\
eval `scramv1 runtime -sh`\n\
cp  "+pwd+"/Bmm5Ntuple* .\n\
cp  "+pwd+"/getEventListFromBMM5.cc .\n\
root -b -q 'getEventListFromBMM5.cc(\"@@CFGFILENAME\")'\n\
if [ $? -eq 0 ]; then \n\
    mv *.root "+destination+"\n\
    mv *.txt "+destination+"\n\
    mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess \n\
    echo OK\n\
else\n\
    mv @@RUNSCRIPT.busy @@RUNSCRIPT \n\
    echo FAIL\n\
fi\n\
"

if not os.path.exists('JobsBMM5'):
    os.system('mkdir JobsBMM5')
print("Making ",NJOBS," Jobs ")

njobs=0
for ii in range(NJOBS):
    i=ii+ZERO_OFFSET
    
    if len(sourceFileList)<FILES_PER_JOB:
       print("fname count less than required .. stoping ")
       FILES_PER_JOB=len(sourceFileList)
    
    if len(sourceFileList) ==0:
       break 

    dirName= pwd+'/JobsBMM5/Job_'+str(i)
    
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
    tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")
    njobs+=1
print(" Number of jobs made : ", njobs)
print(" Number of files left : ", len(sourceFileList) )
condorScript.close()
