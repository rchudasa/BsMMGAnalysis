#!/usr/bin/env python3


ifnames=['bmmGEvntSelection.txt']

run='0'
lumi='0'
evt='0'
i=0
imax=len(ifnames)

fout=open('EventListFile.txt','w')
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
            l=f.readline()
            continue
        items=l[:-1].strip().split(',')
        for evt in items:
            fout.write(run+':'+lumi+':'+evt+'\n')
        l=f.readline()
    f.close()
fout.close()
