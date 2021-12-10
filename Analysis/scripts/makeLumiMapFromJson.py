#!/usr/bin/env python
import sys
import os
import ROOT
import json
import array

iFname="goldenJson.txt"
oFname="goldenJson.root"

if len(sys.argv) >1:
    iFname=sys.argv[1]
if len(sys.argv) >2:
    oFname=sys.argv[2]

ofile=ROOT.TFile(oFname,"RECREATE")
oTree=ROOT.TTree("RunLumiMask","Run , nLumiSets , minLumi[nLumiSets] , maxLumi[nLumiSets]" )

run=array.array('i',[0])
nLumiSets=array.array('i',[0])
minLumi=array.array('i',[0]*1000)
maxLumi=array.array('i',[0]*1000)

oTree.Branch('run'      , run , 'run/I')
oTree.Branch('nLumiSets', nLumiSets, 'nLumiSets/I')
oTree.Branch('minLumi'  , minLumi, 'minLumi[nLumiSets]/I')
oTree.Branch('maxLumi'  , maxLumi, 'maxLumi[nLumiSets]/I')

f=open(iFname)
runlumi_json=json.load(f)
f.close()
f.close()

totalRuns=0
totalLumiSets=0
totalLumis=0
for srun in runlumi_json:
    run[0]=int(srun)
    nLumiSets[0]=len(runlumi_json[srun])
    totalRuns+=1
    totalLumiSets+=nLumiSets[0]
    for i in range(len(runlumi_json[srun])):
        minLumi[i]=runlumi_json[srun][i][0]
        maxLumi[i]=runlumi_json[srun][i][1]
        totalLumis+= 1 + (maxLumi[i]-minLumi[i])
    oTree.Fill()

print("Total number of runs     :",totalRuns)
print("Total number of LumiSets :",totalLumiSets)
print("Total number of Lumis    :",totalLumis)


oTree.Write()
ofile.Close()

