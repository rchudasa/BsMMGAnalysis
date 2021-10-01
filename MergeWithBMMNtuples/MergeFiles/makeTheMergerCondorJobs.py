#!/usr/bin/env python3
import os

NJOBS=1000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_MERGE=1
MERGE_PER_JOB=1


destination='/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/CMSSW_10_6_19_patch2/src/BsMMGAnalysis/MergeWithBMMNtuples/MergeFiles/mergedBmmXfiles'
pwd= os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']


## For seeding jobs 
FileSource="bmmgFileList.txt"
Fnames=open(FileSource,'r')
sourceFileList=Fnames.readlines()
Fnames.close()


ferr = open("error.log",'w')


motherFList=open('bmmGEvntSelection.txt','r')
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


cfgTxt="\
#PARAMS_BEG\n\
OutputPrefix=\n\
OutputFile=bmmXMerged_@@IDX.root\n\
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
condorScript=open('subMergerJobs.sub','w')
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
mv @@RUNSCRIPTNAME @@RUNSCRIPTNAME.busy \n\
eval `scramv1 runtime -sh`\n\
SUCCESS=1\n\
@@ROOTCMD\n\
if [ $SUCCESS -eq 1 ]; then \n\
    mv @@RUNSCRIPTNAME.busy @@RUNSCRIPTNAME.success\n\
else\n\
    mv @@RUNSCRIPTNAME.busy @@RUNSCRIPTNAME \n\
    echo FAIL\n\
fi\n\
"
rootCmd="root -b -q 'mergeBmmXTrees.cc(\"@@CFGFILENAME\")' \n"

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
        mF=0
        snF=0
        sF=0
        for src in srcFiles:
            if src not in motherFMap:
                ferr.write( str(i) + '\t'+"job seed file not found in mother filelist !! -> "+src+'\n')
                mF+=1
                continue
            for run,lumi in zip(motherFMap[src]['run'],motherFMap[src]['lumi']):
                if run not in sonFMap:
                    ferr.write(str(i)+'\t'+"run  not in sonMap files, run = "+str(run)+" , lumi : "+str(lumi)+" , parent : "+src+"\n")
                    lostLumis+=1
                    snF+=1
                    continue
                if lumi not in sonFMap[run]:
                    ferr.write(str(i) + '\t' + "lumi  not in sonMap files, run = "+str(run)+" , lumi : "+str(lumi)+" , parent : "+src+"\n")
                    lostLumis+=1
                    snF+=1
                    continue
                for fn in sonFMap[run][lumi]:
                    if fn not in searchFiels:
                        searchFiels.append(fn)
                        sF+=1
        for fn in searchFiels:
            tmp+= fn + "\n"
        cfgTxttmp=cfgTxttmp.replace("@@SONFNAMES",tmp)
        cfgTxttmp=cfgTxttmp.replace("@@MAXEVENTS",str(NEVENTS_PER_JOB))
        cfgFile.write(cfgTxttmp)
        cfgFile.close()   
        rootCMDtemp+=rootCmd.replace("@@CFGFILENAME",cfgFileName)
        rootCMDtemp+=exitCheckScript
        print("\t Job made with nSource = ",len(srcFiles),
              " nMothers failed : ",mF,
              " nSons = ",sF ," , nSons Lost =  ",snF )

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
