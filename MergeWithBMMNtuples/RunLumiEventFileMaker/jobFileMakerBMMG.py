#!/usr/bin/env python3
import os
import sys

UsageStr = "\
\n./jogFileMakerBMMG.py <FileSource> <destination> <filesPerJob> <njobs> <tag> \
"
NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=20
destination='/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/ntuplizer/data/CMSSW_10_6_4_patch1/src/BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/RunLumiFiles_charmoniumD'

FileSource ="bmmgFileList.txt"
tag=""
if len(sys.argv) > 1:
    FileSource=sys.argv[1]
else:
    print(UsageStr)
    sys.exit(1)
if len(sys.argv)>2:
    destination=os.path.abspath(sys.argv[2])
if len(sys.argv) > 3:
    FILES_PER_JOB=int(sys.argv[3])
if len(sys.argv) > 4:
    NJOBS=int(sys.argv[4])
if len(sys.argv) > 5:
    tag=sys.argv[5]

print("Source  : ",FileSource)
print("Destn : ",destination)
print("Files Per Job : ",FILES_PER_JOB)
print("NJOBS : ",NJOBS)


pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']

Fnames=open(FileSource,'r')
sourceFileList=Fnames.readlines()
Fnames.close()

if not os.path.exists(destination):
    os.system('mkdir -p '+destination)    

print("TAG         : ",tag)
print("N Fnames    : ",len(sourceFileList))
print("File source : ",FileSource)
print("FILES_PER_JOB : ",FILES_PER_JOB)
print("NJobs : ",NJOBS )
print("Destination : ",destination)



configurationTxt="\
#FLIST_BEG\n\
@@FNAMES\n\
#FLIST_END\n\
\n\
#PARAMS_BEG\n\
MaxEvents=-1\n\
OutputPrefix=\n\
OutputFile=bmmGEvntSelection_@@IDX.root\n\
OutputFileTxt=bmmGEvntSelection_@@IDX.txt\n\
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


runScriptTxt="\
#!/bin/bash\n\
set -x\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME="+HOME+"\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd @@DIRNAME \n\
mv @@RUNSCRIPT @@RUNSCRIPT.busy \n\
eval `scramv1 runtime -sh`\n\
TMPDIR=`mktemp -d`\n\
cp @@CFGFILENAME $TMPDIR \n\
cd $TMPDIR\n\
cp  "+pwd+"/BmmGNtuple* .\n\
cp  "+pwd+"/getEventListFromBMMG.cc .\n\
root -b -q 'getEventListFromBMMG.cc(\"@@CFGFILENAME\")'\n\
if [ $? -eq 0 ]; then \n\
    mv *.root "+destination+"\n\
    mv *.txt  "+destination+"\n\
    mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess \n\
    echo OK\n\
else\n\
    mv @@RUNSCRIPT.busy @@RUNSCRIPT \n\
    echo FAIL\n\
fi\n\
rm  BmmGNtuple* \n\
"

head = "Condor/JobsBmmG" + tag
if not os.path.exists(head):
    os.system('mkdir '+ head)
condorScriptName=head+'/subCondorBMMG'+tag+'.sub'
condorScript=open(condorScriptName,'w')
condorScript.write(condorScriptString)



condorScriptName=head+'/subCondorBMMG'+tag+'.sub'
condorScript=open(condorScriptName,'w')
condorScript.write(condorScriptString)



n = int( len(sourceFileList)/FILES_PER_JOB ) + 1
print("Making ",n," Jobs ")

njobs=0
for ii in range(NJOBS):
    i=ii+ZERO_OFFSET
    
    if len(sourceFileList)<FILES_PER_JOB:
       print("\nfname count less than required .. stoping ")
       FILES_PER_JOB=len(sourceFileList)
    
    if len(sourceFileList) ==0:
       break 

    dirName= pwd+'/'+head+'/Job_'+str(i)
    if(ii%25==0):
        print("\n Job Made : ",end="")
    print( ii ,end =" ")

    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir '+dirName)
    
    cfgFileName='bmmGSelection_'+str(i)+'.cfg'
    cfgFile=open(dirName+'/'+cfgFileName,'w')
    tmp=""
    for j in range(FILES_PER_JOB):
      tmp+= sourceFileList.pop(0)[:-1]+"\n"
    tmp=configurationTxt.replace("@@FNAMES",tmp[:-1])
    tmp=tmp.replace("@@IDX",str(i))
    cfgFile.write(tmp)
    cfgFile.close()   
    
    runScriptName=dirName+'/'+tag+'_run'+str(i)+'.sh'
    if os.path.exists(runScriptName+".*"):
        os.system("rm "+runScriptName+'.*')
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
    tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")
    njobs+=1
condorScript.close()
print()
print(" Number of jobs made : ", njobs)
print(" condor submit file  : ",condorScriptName  )
print(" Number of files left : ", len(sourceFileList) )
