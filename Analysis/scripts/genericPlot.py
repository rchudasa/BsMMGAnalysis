#!/usr/bin/env python3
import copy
import ROOT
import PlotUtil as putil
import HistUtil as hutil
import CanvasSetup as cSetup
import sys

UsageString="\
Usage :\n\
    ./genericPlot.py <hName> <Description> <xmin> <xmax> <ymin> <ymax> <rebin> 
"

plots = []
inputFile2018         = ROOT.TFile.Open("../data/bmmSelectionNoAlpha_analysis2018A.root")


cPars=cSetup.getCanvas2018A()

plotDir="plots/2018A"
pname="bmmg_mmgMass"
legend="m(#mu#mu#gamma)"
xtitle="GeV"
xmin=1.0
xmax=8.0
ymin=0.0
ymax=100.0

rebin = -1
logy=False

idx=1

if len(sys.argv)>idx:
    pname=sys.argv[idx]
else:
    print(UsageString)
idx+=1
if len(sys.argv)>idx:
    legend=sys.argv[idx]
idx+=1
if len(sys.argv)>idx:
    xmin=float(sys.argv[idx])
idx+=1
if len(sys.argv)>idx:
    xmax=float(sys.argv[idx])
idx+=1
if len(sys.argv)>idx:
    ymin=float(sys.argv[idx])
idx+=1
if len(sys.argv)>idx:
    ymax=float(sys.argv[idx])
idx+=1
if len(sys.argv)>idx:
    rebin=int(sys.argv[idx])
idx+=1

plots.append(putil.PlotCanvas(prefix=plotDir,canvasPars=cPars))
plots[-1].name = pname
plots[-1].xRange = (xmin,xmax)
plots[-1].yRange = (ymin,ymax)
plots[-1].logy   = logy
plots[-1].yTitle = "Events"
plots[-1].xTitle = xtitle
plots[-1].desc=[legend]
plots[-1].legendPosition = (0.53,0.28,0.93,0.40)
plots[-1].descPosition =(0.6,0.58)

histo0 = inputFile2018.Get(pname)
histo0.__class__ = ROOT.TH1D
if rebin>0:
    histo0.Rebin(rebin)
aplot = putil.Plot(Name="bmmg_mmgMass", Histo=histo0)
plots[-1].addPlot(aplot)


canvas = []
for plot in plots:
    canvas.append(plot.plot())

inputFile2018.Close()
