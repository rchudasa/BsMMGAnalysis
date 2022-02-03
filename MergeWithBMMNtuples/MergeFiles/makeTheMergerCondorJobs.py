#!/usr/bin/env python3
import os
import sys 
import json

##################


NJOBS=-1
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_MERGE=3
MERGE_PER_JOB=1

destination='/eos/home-a/athachay/workarea/data/BsToMuMuGamma/Run2Studies/MergedNtuples/BsToMuMuGammaNtuples/Charmonium/crab_charmonium2018A_ntuplizer_V2'
bkpDestination=destination
pwd= os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
FileSource="bmmgFileList.txt"
mFileNames="bmmGEvntSelection.txt"
sFileNames="bmm5EvntSelection.txt"
tag=""

if len(sys.argv)<3 :
    print("Usage \n\t\t ./makeTheMergerCondorJobs.py <File source> <mother files list> <son file list> <destination> <tag> <njobs> <bkpDestination> ")
    sys.exit(1)
if len(sys.argv)>1:
    FileSource=sys.argv[1]
if len(sys.argv)>2:
    mFileNames=sys.argv[2]
if len(sys.argv)>3:
    sFileNames=sys.argv[3]
if len(sys.argv)>4:
    destination=sys.argv[4]
if len(sys.argv)>5:
    tag=sys.argv[5]
if len(sys.argv)>6:
    NJOBS=int(sys.argv[6])
if len(sys.argv)>7:
    bkpDestination=sys.argv[7]

#destination=os.path.abspath(destination)
bkpDestination=os.path.abspath(bkpDestination)

print("File Source : ",FileSource)
print("Mother flist : ",mFileNames)
print("Son flist    : ",sFileNames)
print("Destinaltion    : ",destination)
print("bkpDestinaltion    : ",bkpDestination)
print("tag  : ",tag)
print("Njobs : ",NJOBS)
print("Files per Merge : ",FILES_PER_MERGE)
print("Merges per job : ",MERGE_PER_JOB)

if not os.path.exists(destination):
    os.system('mkdir -p '+destination)

if not os.path.exists(destination):
    print("Unable to make path !! [ ",destination," ]")
    exit(1)
## For seeding jobs 
print("Reading the master filelist to sed the jobs !! ")
Fnames=open(FileSource,'r')
sourceFileList=Fnames.readlines()
print("Read ",len(sourceFileList)," seed files")
Fnames.close()
n=int(len(sourceFileList)/MERGE_PER_JOB/FILES_PER_MERGE) + 1
print("Making ",n," Jobs ")

ferr = open("error.log",'w')



print("\nReading the mother FileList :  ",mFileNames)
motherFMap={}
n=0
f=open(mFileNames,'r')
flist=f.readlines()
f.close()
print("\t Files to read  : ",len(flist))
cacheName='.cache/'+mFileNames.split('/')[-1]+'.cache.json'
motherRunLumiMap={}
if  True or not os.path.exists(cacheName):
    for mFileName in flist:
        mFileName=mFileName[:-1]
        motherFList=open(mFileName,'r')
        l=motherFList.readline()
        n+=1
        if(n%20 == 1):
            tstr="\tRead file  : "
            print()
            print(tstr,end=" ")
        tstr+=" " + str(n)
        print("\r"+tstr,end=" ")
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
                if run not in motherRunLumiMap:
                    motherRunLumiMap[run]={lumi:[]}
                if lumi not in motherRunLumiMap[run]:
                    motherRunLumiMap[run][lumi]=[]
                if fname not in motherRunLumiMap[run][lumi]:
                    motherRunLumiMap[run][lumi].append(fname)


            l=motherFList.readline()
        motherFList.close()
    print("cache does not exist ! making one : ",cacheName)
    os.system('mkdir -p .cache')
    f=open(cacheName,'w')
    json.dump(motherFMap,f,indent=4)
    f.close()
else:
    print()
    print("Loading mfile from cache : ",cacheName)
    f=open(cacheName,'r')
    motherFMap=json.load(f)
    #f.close()
print("Mother Map made with ",len(motherFMap)," entries")
#for ky in motherFMap:
#    print("\tfor ",ky," we have #runlumi = ",len(motherFMap[ky]['run']))
print()
print("\nReading the son FileList :  ",sFileNames)
sonFMap={}
n=0
f=open(sFileNames,'r')
flist=f.readlines()
f.close()
print("\t Files to read  : ",len(flist))
cacheName='.cache/'+sFileNames.split('/')[-1]+'.cache.json'

missingFiles=[]
if True or not os.path.exists(cacheName):
    for sFileName in flist:
        sFileName=sFileName[:-1]
        sonFList=open(sFileName,'r')
        l=sonFList.readline()
        n+=1
        if(n%20 == 1):
            print()
            tstr="\r\tRead file  : "
            print(tstr,end=" ")
        tstr+=" " + str(n)
        print(tstr,end="")
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
    #            fname=fname.replace('AOD/516/Cha','AOD/517/Cha')
    #            if( not os.path.exists(fname)):
    #                fname=fname.replace('AOD/517/Cha','AOD/518/Cha')
    #                if fname not in missingFiles:
    #                    print("\n\tSon Fname not found in disk : ",fname)
    #                    missingFiles.append(fname)
    #                    continue
                if fname not in sonFMap[run][lumi]:
                    sonFMap[run][lumi].append(fname)
            l=sonFList.readline()
        sonFList.close()

    print()
    print("cache does not exist ! making one ")
    if not os.path.exists('.cache'):
        os.system('mkdir -p .cache')
    f=open(cacheName,'w')
    json.dump(sonFMap,f,indent=4)
    f.close()
else:
    print("Loading sfile from cache : ",cacheName)
    f=open(cacheName,'r')
    sonFMap=json.load(f)
    f.close()

f=open('testJson.json','w')
json.dump(sonFMap,f,indent=4)
f.close()
print(" ================== Peparatory Summary ================== \n")

print("Mother Map made with ",len(motherRunLumiMap)," runs")
for ky in motherRunLumiMap:
    print("\t r = ",ky," has ",len(motherRunLumiMap[ky]), " lumis ")
print()

print("Son Map made with ",len(sonFMap)," runs")
for ky in sonFMap:
    print("\t r = ",ky," has ",len(sonFMap[ky]), " lumis ")
print("======================================================== \n")

cfgTxt="\
#PARAMS_BEG\n\
OutputPrefix=\n\
OutputFile=bmmXMerged2018data_@@TAG_@@IDX.root\n\
MaxEvents=@@MAXEVENTS\n\
#PARAMS_END\n\
#MOTHER_FILELIST_BEG\n\
@@MOTHERFNAMES\n\
#MOTHER_FILELIST_END\n\
\n\
#SON_FILELIST_BEG\n\
@@SONFNAMES\n\
#SON_FILELIST_END\n\
"


condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)cdr.stdout\n\
error = $Fp(filename)cdr.stderr\n\
log = $Fp(filename)cdr.log\n\
+JobFlavour = \"microcentury\"\n\
"
condorScriptName='subMergerJobs'+tag+'.sub'
condorScript=open(condorScriptName,'w')
condorScript.write(condorScriptString)


exitCheckScript="\
if [ $? -eq 0 ]; then \n\
    xrdcp *.root "+destination+"\n\
    if [ $? -ne 0 ] ; then\n\
        mv *.root "+bkpDestination+"\n\
        if [ $? -ne 0 ] ; then\n\
            echo FAIL\n\
            SUCCESS=0\n\
        fi\n\
    else\n\
        rm *.root\n\
        mv @@CFGFILENAME @@CFGFILENAME.success\n\
        cp @@CFGFILENAME.success "+bkpDestination+"\n\
        echo OK\n\
   fi\n\
else\n\
    SUCCESS=0\n\
    echo FAIL\n\
fi\n\
date\n\
"

runScriptTxt="\
#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME=/home/athachay\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd "+pwd+"\n\
eval `scramv1 runtime -sh`\n\
set -x\n\
TMPDIR=`mktemp -d`\n\
cp Bmm5Ntuple* $TMPDIR\n\
cp BMMGNtuple* $TMPDIR\n\
cp genericTreeBranchSelector* $TMPDIR\n\
cp *.h $TMPDIR\n\
cp mergeBmmXTrees.cc $TMPDIR\n\
cd $TMPDIR \n\
mv @@RUNSCRIPTNAME @@RUNSCRIPTNAME.busy \n\
SUCCESS=1\n\
date\n\
@@ROOTCMD\n\
if [ $SUCCESS -eq 1 ]; then \n\
    mv @@RUNSCRIPTNAME.busy @@RUNSCRIPTNAME.success\n\
else\n\
    mv @@RUNSCRIPTNAME.busy @@RUNSCRIPTNAME \n\
    echo FAIL\n\
fi\n\
"
rootCmd="\
date\n\
root -b -q 'mergeBmmXTrees.cc(\"@@CFGFILENAME\")' \n"

#print(sonFMap)
#for r in sonFMap:
#    print("\n\n r = ",r," -> ",end =" ")
#    for l in sonFMap[r]:
#        print( l , end=" , ")

head = pwd+'/Jobs'+tag
if not os.path.exists(head):
    os.system('mkdir -p '+head)

lostLumis=0
foundLumi={}

njbs=int(len(sourceFileList)/MERGE_PER_JOB/FILES_PER_MERGE) + 1
if len(sourceFileList)%(MERGE_PER_JOB*FILES_PER_MERGE):
    njbs-=1
if njbs==0:
    njbs=1
print("Making ",njbs," Jobs\n")
jobsMade=0
jobsSkipped=0
for ii in range(NJOBS):
    if FILES_PER_MERGE==0:
        break
    i=ii+ZERO_OFFSET
    dirName=head+'/Job_'+str(i)
    if i%10==0:
        print("\nJob Made : ",end =" ")
    print(" ",i,end =" ")

    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir '+dirName)
    rootCMDtemp=""
    runScriptBareName='run'+str(i)+'.sh'
    hasAtleastAMerge=False
    for jj in range(MERGE_PER_JOB):
        if len(sourceFileList)<FILES_PER_MERGE:
            FILES_PER_MERGE=len(sourceFileList)
            print("fname count less than required .. stoping after seting final script to  have ",FILES_PER_MERGE," files ")
        if FILES_PER_MERGE==0:
            print("Filelist empty exiting")
            break

        cfgFileName=dirName+'/merger_'+str(i)+"_"+str(jj)+'.cfg'
        cfgFile=open(cfgFileName,'w')
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
        hasAMotherWithSon=False
        for src in srcFiles:
            if src not in motherFMap:
 #               print( str(i) + '\t'+"job seed file not found in mother filelist !! -> "+src)
                ferr.write( str(i) + '\t'+"job seed file not found in mother filelist !! -> "+src+'\n')
                mF+=1
                continue
            for run,lumi in zip(motherFMap[src]['run'],motherFMap[src]['lumi']):
                if run not in sonFMap:
 #                   print(str(i)+'\t'+"run  not in sonMap files, run = "+str(run)+" , lumi : "+str(lumi)+" , parent : "+src)
                    ferr.write(str(i)+'\t'+"run  not in sonMap files, run = "+str(run)+" , lumi : "+str(lumi)+" , parent : "+src+"\n")
                    lostLumis+=1
                    snF+=1
                    continue
                if lumi not in sonFMap[run]:
 #                   print(str(i) + '\t' + "lumi  not in sonMap files, run = "+str(run)+" , lumi : "+str(lumi)+" , parent : "+src)
                    ferr.write(str(i) + '\t' + "lumi  not in sonMap files, run = "+str(run)+" , lumi : "+str(lumi)+" , parent : "+src+"\n")
                    lostLumis+=1
                    snF+=1
                    continue
                if run not in foundLumi:
                    foundLumi[run]=set()
                foundLumi[run].add(lumi)
                for fn in sonFMap[run][lumi]:
                    if fn not in searchFiels:
                        searchFiels.append(fn)
                        sF+=1
        for fn in searchFiels:
            tmp+= fn + "\n"
            hasAMotherWithSon=True
        
        if not hasAMotherWithSon : continue
        hasAtleastAMerge=True
        cfgTxttmp=cfgTxttmp.replace("@@SONFNAMES",tmp)
        cfgTxttmp=cfgTxttmp.replace("@@MAXEVENTS",str(NEVENTS_PER_JOB))
        cfgTxttmp=cfgTxttmp.replace("@@TAG",tag)
        cfgFile.write(cfgTxttmp)
        cfgFile.close()   
        rootCMDtemp+=rootCmd.replace("@@CFGFILENAME",cfgFileName)
        rootCMDtemp+=exitCheckScript.replace("@@CFGFILENAME",cfgFileName)
        
        jobStr = "\t Job made with nSource = "+str(len(srcFiles))+\
              " nMothers failed : "+str(mF)+\
              " nSons = "+str(sF)+" , nSons Lost =  "+str(snF) +'\n'
        #print(jobStr)
        ferr.write(jobStr)
    if not hasAtleastAMerge :
        jobsSkipped+=1
        continue
    jobsMade+=1
    runScriptName=dirName+'/'+runScriptBareName
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@ROOTCMD",rootCMDtemp)
    tmp=tmp.replace("@@RUNSCRIPTNAME",runScriptName)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")
print()
lumiF=0
print("\n\n================  Merger Preperation Summary ============  \n")
print("Summary of Matches found ")
for r in foundLumi:
    print("\t in run : ",r," found ",len(foundLumi[r])," lumis ")
    lumiF+=len(foundLumi[r])
print("\n")
print("\t Total  found runs : ",len(foundLumi),"  ,  total found lumis : ", lumiF)

print("Lost Lumis : " , lostLumis," Number of files left : ", len(sourceFileList) )
condorScript.close()
 
print("Njbs : ",njbs, "Jobs Made : ",jobsMade," Jobs Skipped : " ,jobsSkipped)
print("Condor script : ",condorScriptName)
print("========================================================\n")
print("\nMissing files ")
for i in missingFiles:
    print(i)
