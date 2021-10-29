#!/usr/bin/env python3

f=open('bmmGEvntSelection.txt','r')
txt=f.readlines()
f.close()

ifnames=[]
for l in txt:
    ifnames.append(l[:-1])

run='0'
lumi='0'
evt='0'
i=0
imax=len(ifnames)

fout    = open('EventListFile.txt','w')

runLumiDict={}
for fname in ifnames:
    print( " doing : ",i," / ",imax,"\t",fname)
    i+=1
    f=open(fname,'r')
    l=f.readline()
    while l:
        if(len(l)<2): 
            l=f.readline()
            continue
        if(l[0]=='@'):
            run,lumi,fname=l[2:-1].split(',')
            rint=int(run)
            lint=int(lumi)
            if rint not in runLumiDict:
                runLumiDict[rint]=[]
            runLumiDict[rint].append(lint)
            l=f.readline()
            continue
        items=l[:-1].strip().split(',')
        for evt in items:
            fout.write(run+':'+lumi+':'+evt+'\n')
        l=f.readline()
    f.close()

sortedRuns=list(runLumiDict.keys())
sortedRuns.sort()

refinedRunLumiDict={}


print(" processing the run-lumi json file")
idx=0
maxLen=len(sortedRuns)
for r in sortedRuns:
    refinedRunLumiDict[r]=[]
    runLumiDict[r].sort()
    print("[",idx+1,"/",maxLen,"]",end="")
    idx+=1

    if len(runLumiDict[r]) <1 : continue
    lumiStart=runLumiDict[r][0]
    prevLumi=lumiStart
    for l in runLumiDict[r]:
#        print(l,prevLumi)
#        print(type(l),type(prevLumi))
        if l==prevLumi or l==(prevLumi+1):
            prevLumi=l
            continue
        refinedRunLumiDict[r].append([lumiStart,prevLumi])
#        print("Adding r = ",r," ls le : ",lumiStart,prevLumi )
        lumiStart=l
        prevLumi=l
    refinedRunLumiDict[r].append([lumiStart,prevLumi])

print()
print("writing out the run-lumi json file")
jsonOut = open('pickEvents.json','w')
jsonOut.write("{ ")
isFirstRun=True
idx=0
for r in refinedRunLumiDict:
    print("[",idx+1,"/",maxLen,"]",end="")
    idx+=1
    if not isFirstRun:
        jsonOut.write(", \n")
    isFirstRun=False
    jsonOut.write('"'+str(r)+'": [')
    isFirst=True
    for l in refinedRunLumiDict[r]:
        if not isFirst:
            jsonOut.write(", ")
        jsonOut.write("[ "+str(l[0])+", "+str(l[1])+"]")
        isFirst=False
    jsonOut.write(" ]")
jsonOut.write("}")
print()
jsonOut.close()
fout.close()
print()
