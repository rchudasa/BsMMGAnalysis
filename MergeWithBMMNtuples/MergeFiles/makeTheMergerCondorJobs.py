#!/usr/bin/env python3
import os

NJOBS=1000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_MERGE=2
MERGE_PER_JOB=2


destination='/eos/home-a/athachay/workarea/data/BsToMuMuGamma/Run2Studies/MergedNtuples/BsToMuMuGammaNtuples/Charmonium/crab_charmonium2018A_ntuplizer/0000/'
pwd= os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']

xrdRedirector="root://cms-xrd-global.cern.ch/"

FileSource="bmmgFiles.txt"
Fnames=open(FileSource,'r')
sourceFileList=Fnames.readlines()
Fnames.close()

cfgTxt="\
#PARAMS_BEG\n\
OutputPrefix=\n\
OutputFile=mergedCharmonium2018A0000_BMMG_@@IDX.root\n\
MaxEvents=-1\n\
#PARAMS_END\n\
#MOTHER_FILELIST_BEG\n\
@@MOTHERFNAMES\n\
#MOTHER_FILELIST_END\n\
\n\
#SON_FILELIST_BEG\n\
@@SONFNAMES\n\
#SON_FILELIST_BEG\n\
"


condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)cdr.stdout\n\
error = $Fp(filename)cdr.stderr\n\
log = $Fp(filename)cdr.log\n\
+JobFlavour = \"longlunch\"\n\
"
condorScript=open('jobSubmit.sub','w')
condorScript.write(condorScriptString)


exitCheckScript="\
if [ $? -eq 0 ]; then \n\
    echo OK\n\
    mv *.root "+destination+"\n\
else\n\
    SUCCESS=0\n\
    echo FAIL\n\
fi\n\
"

runScriptTxt="\
#!/bin/bash\n\
set -x\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME=/home/athachay\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd "+pwd+"\n\
cp Bmm5Ntuple* @@DIRNAME\n\
cp BMMGNtuple* @@DIRNAME\n\
cp genericTreeBranchSelector* @@DIRNAME\n\
cp mergeBmmXTrees.cc @@DIRNAME\n\
cd "+pwd+"/@@DIRNAME \n\
eval `scramv1 runtime -sh`\n\
SUCCESS=1\n\
@@ROOTCMD\n\
if [ $SUCCESS -eq 1 ]; then \n\
    mv @@RUNSCRIPTNAME @@RUNSCRIPTNAME.success\n\
fi \n\
"
rootCmd="root -b -q 'mergeBmmXTrees.cc(\"@@CFGFILENAME\")' \n"

motherFList=open('bmmgEvntSelection.txt','r')
l=motherFList.readline()
motherFMap={}
while l:
    if l[0]=='@':
        items=l[1:-1].split(',')
        run=int(items[0])
        lumi=int(items[1])
        fname=items[2]
        if fname not in motherFMap:
            motherFMap[fname]={'run':[],'lumi':[]}
        motherFMap[fname]['run'].append(run)
        motherFMap[fname]['lumi'].append(lumi)
    l=motherFList.readline()
motherFList.close()

sonFList=open('bmm5EvntSelection.txt','r')
l=sonFList.readline()
sonFMap={}
while l:
    if l[0]=='@':
        items=l[1:-1].split(',')
        run=int(items[0])
        lumi=int(items[1])
        fname=items[2]
        if run not in sonFMap:
            sonFMap[run]={lumi:[]}
        if lumi not in sonFMap[run]:
            sonFMap[run][lumi]=[]
        sonFMap[run][lumi].append(fname)
    l=sonFList.readline()
sonFList.close()
#print(sonFMap)
#for r in sonFMap:
#    print("\n\n r = ",r," -> ",end =" ")
#    for l in sonFMap[r]:
#        print( l , end=" , ")

if not os.path.exists('Jobs'):
    os.system('mkdir Jobs')
print("Making ",NJOBS," Jobs ")
lostLumis=0
for ii in range(NJOBS):
    if FILES_PER_MERGE==0:
        break
    i=ii+ZERO_OFFSET
    dirName= 'Jobs/Job_'+str(i)
    print(i," Job Made")
    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir '+dirName)
    rootCMDtemp=""
    runScriptBareName='run'+str(i)+'.sh'
    for jj in range(MERGE_PER_JOB):
        if len(sourceFileList)<FILES_PER_MERGE:
            FILES_PER_MERGE=len(sourceFileList)
            print("fname count less than required .. stoping after seting final script to  have ",FILES_PER_MERGE," files ")
        if FILES_PER_MERGE==0:
            print("Filelist empty exiting")
            break

        cfgFileName='merger'+str(jj)+'.cfg'
        cfgFile=open(dirName+'/'+cfgFileName,'w')
        idx=str(i)+'_'+str(jj)
        cfgTxttmp=cfgTxt.replace("@@IDX",idx)
        tmp=""
        srcFiles=[]
        for j in range(FILES_PER_MERGE):
            srcFiles.append(sourceFileList.pop(0)[:-1])
            tmp+=srcFiles[-1]+"\n"
        cfgTxttmp=cfgTxttmp.replace("@@MOTHERFNAMES",tmp)
        
        tmp=""
        searchFiels=[]
        for src in srcFiles:
            for run,lumi in zip(motherFMap[src]['run'],motherFMap[src]['lumi']):
                if run not in sonFMap:
                    print("r  not in search files, run = ",run," , lumi : ",lumi," , parent : ",src)
                    lostLumis+=1
                    continue
                if lumi not in sonFMap[run]:
                    print("l  not in search files, run = ",run," , lumi : ",lumi," , parent : ",src)
                    lostLumis+=1
                    continue
                for fn in sonFMap[run][lumi]:
                    if fn not in searchFiels:
                        searchFiels.append(fn)
        for fn in searchFiels:
            tmp+= fn + "\n"
        cfgTxttmp=cfgTxttmp.replace("@@SONFNAMES",tmp)
        
        cfgTxttmp=cfgTxttmp.replace("@@MAXEVENTS",str(NEVENTS_PER_JOB))
        cfgFile.write(cfgTxttmp)
        cfgFile.close()   
        rootCMDtemp+=rootCmd.replace("@@CFGFILENAME",cfgFileName)
        rootCMDtemp+=exitCheckScript
    runScriptName=dirName+'/'+runScriptBareName
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@ROOTCMD",rootCMDtemp)
    tmp=tmp.replace("@@RUNSCRIPTNAME",runScriptBareName)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")

print("Lost Lumis : " , lostLumis," Number of files left : ", len(sourceFileList) )
condorScript.close()
