#!/usr/bin/env python3
import copy
import ROOT
import PlotUtil as putil
import HistUtil as hutil
import CanvasSetup as cSetup

plots = []
inputFile2018         = ROOT.TFile.Open("../data/bmmSelectionNoAlpha_analysis2018A.root")


cPars=cSetup.getCanvas2018A()

plotDir="plots/2018A/"

plots.append(putil.PlotCanvas(prefix=plotDir,canvasPars=cPars))
plots[-1].name = "bmmgMass"
plots[-1].xRange = (1.0,8.0)
plots[-1].yRange = (0.0,100.0)
plots[-1].yTitle = "Candidate/ 100 MeV"
plots[-1].xTitle = "m(#mu#mu#gamma) [GeV]"
plots[-1].desc=["m(#mu#mu#gamma)"]
plots[-1].legendPosition = (0.53,0.28,0.93,0.40)
plots[-1].descPosition =(0.6,0.58)

histo0 = inputFile2018.Get("bmmg_mass")
histo0.__class__ = ROOT.TH1D
histo0.Rebin()
aplot = putil.Plot(Name="bmmg_mmgMass", Histo=histo0)
plots[-1].addPlot(aplot)

plots.append(putil.PlotCanvas(prefix=plotDir,canvasPars=cPars))

plots[-1].name = "bmmg_mumu_Mass"
plots[-1].xRange = (1.0,8.0)
plots[-1].yRange = (0.0,100.0)
plots[-1].yTitle = "Candidate/ 100 MeV"
plots[-1].xTitle = "m(#mu#mu) [GeV]"
plots[-1].desc=["m(#mu#mu)"]
plots[-1].legendPosition = (0.53,0.28,0.93,0.40)
plots[-1].descPosition =(0.6,0.58)

histo0 = inputFile2018.Get("bmmg_mmMass")
histo0.__class__ = ROOT.TH1D
histo0.Rebin(5)
aplot = putil.Plot(Name="bmmg_mass", Histo=histo0)
plots[-1].addPlot(aplot)



canvas = []
for plot in plots:
    canvas.append(plot.plot())

inputFile2018.Close()
