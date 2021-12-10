import ROOT 
import operator 
import array

from CanvasSetup import *

ROOT.gSystem.Load('libRooFit')


class Plot:
    def __init__(self, **args):
        self.name        = args.get("Name", "plot")
        self.legend      = args.get("Legend",None)
        self.histo       = args.get("Histo", None)
        self.fit         = args.get("Fit", None)
        self.markerColor = args.get("MarkerColor", ROOT.kBlue+2)
        self.markerStyle = args.get("MarkerStyle", 22)
        self.lineColor   = args.get("LineColor", ROOT.kBlue+2)
        self.lineStyle   = args.get("LineStyle", 1)
        self.options   = args.get("Options", "")
        self.drawLegend   = args.get("DrawLegend",False)
        self.drawLine     = args.get("DrawLine", False)
        self.histo.SetName(self.name+"_histo")


class PlotCanvas:
    def __init__(self, **args):
        self.name  = ""
        self.plots = []
        self.rates = []
        self.plotDir =  args.get("prefix","plots/")
        self.yRange = (0.0, 1.0)
        self.xRange = (0, 100)
        self.yRange2 = (0.0, 1.0)
        self.xTitle = args.get("xTitle","E_{T}^{e offl} [GeV]")
        self.yTitle = args.get("yTitle","a.u")
        self.yTitle2 = args.get("yTitle","a.u")
        self.legendPosition = (0.4,0.2,0.9,0.6)
        self.descPosition = (0.6,0.58)
        self.desc = args.get("desc",["L1 Efficiency"])
        self.logx = args.get("logx",False)
        self.logy = args.get("logy",False)
        self.logy2 = args.get("logy2",False)
        self.canvasParams = args.get("canvasPars",getDefaultCanvas())    
        self.setPlotStyle()

    def addRate(self, rate):
        self.rates.append(rate)
    def addPlot(self, plot):
        self.plots.append(plot)
    def clearPlots(self):
        self.plots=[]


    def plot(self):
        canvas = ROOT.TCanvas("c_"+self.name, self.name, 900, 800)
        canvas.SetGrid()
        hDummy = ROOT.TH1F("hDummy_"+self.name, self.name, 1, self.xRange[0], self.xRange[1])
        hDummy.SetAxisRange(self.yRange[0], self.yRange[1], "Y")
        hDummy.SetXTitle(self.xTitle)
        hDummy.SetYTitle(self.yTitle)
        hDummy.Draw()
        if self.logx : canvas.SetLogx()
        if self.logy : canvas.SetLogy()
        self.logz= True
        if self.logz : canvas.SetLogz()
        
        
        CMSbox=self.getCMSBox()
        extraTextBox=self.getExtraTextBox()
        lumibox=self.getLumiBox()


        selbox=[]
        for n in range(len(self.desc)):
            selbox.append(ROOT.TLatex  (self.descPosition[0], self.descPosition[1] -n*0.04, self.desc[n]))
            selbox[-1].SetNDC()
            selbox[-1].SetTextSize(0.035)
            selbox[-1].SetTextFont(12)
            selbox[-1].SetTextAlign(13)

        #Line legend
        legend = ROOT.TLegend(self.legendPosition[0],self.legendPosition[1],self.legendPosition[2],self.legendPosition[3])
        legend.SetTextFont(42)
        legend.SetFillColor(0)

        for plot in self.plots:
            histo = plot.histo
            histo.SetMarkerStyle(plot.markerStyle)
            histo.SetMarkerColor(plot.markerColor)
            histo.SetLineColor(plot.markerColor)
            histo.Draw("ple1 same"+ plot.options)
            if(plot.legend):
                legend.AddEntry(histo, plot.legend, "pe")
                legend.SetTextSize(0.025) 
            if(plot.drawLegend) :  legend.Draw()
        
        scale = (self.yRange[1] -self.yRange[0])/(self.yRange2[1] - self.yRange2[0])
        y20=self.yRange2[0]

        if(self.logy2):
            scale = (self.yRange[1] -self.yRange[0])/(ROOT.log(self.yRange2[1]) - ROOT.log(self.yRange2[0]) )
            y20 = ROOT.log(self.yRange2[0])


        #####   begin :            For plotting a graph in the Second Axis         #####

        if len(self.rates) >0 :
            axis=None
            canvas.SetRightMargin(0.18);
            if(self.logy2):
                axis = ROOT.TGaxis(self.xRange[1],self.yRange[0],self.xRange[1], self.yRange[1],self.yRange2[0],self.yRange2[1],510,"+LG")
            else :
                axis = ROOT.TGaxis(self.xRange[1],self.yRange[0],self.xRange[1], self.yRange[1],self.yRange2[0],self.yRange2[1],510,"+L")
            axis.SetLineColor(ROOT.kRed);
            axis.SetTextColor(ROOT.kRed);
            axis.SetTitle(self.yTitle2)
            axis.SetTitleColor(ROOT.kBlack)
            axis.SetTitleOffset(1.2)
            axis.Draw();

        for rate in self.rates:
            histo = rate.histo
            if self.logy2:
                for i in range(histo.GetNbinsX()):
                    content=histo.GetBinContent(i+1)
                    scaledContent = ( ROOT.log(content + 1e-20 )  - y20 )* scale
                    histo.SetBinContent(i+1,scaledContent)
                    #print("\t ",content," -> ",scaledContent,"  = ( ",ROOT.log(content + 1e-20),"  - ",y20," )* ",scale)
            else:
                for i in range(histo.GetNbinsX()):
                    content=histo.GetBinContent(i+1)
                    scaledContent = (content  - y20 )* scale
                    histo.SetBinContent(i+1,scaledContent)
            
            histo.SetMarkerStyle(rate.markerStyle)
            histo.SetMarkerColor(rate.markerColor)
            histo.SetLineColor(rate.markerColor)
            histo.Draw("pl same"+ rate.options)
            legend.AddEntry(histo, rate.legend, "pe")
            legend.SetTextSize(0.025) 
            if(rate.drawLegend) :  legend.Draw()
        
        #####   end :            For plotting a graph in the Second Axis         #####
            
        CMSbox.Draw()
        extraTextBox.Draw()
        lumibox.Draw()
        for selb in selbox:
            selb.Draw()
        
        canvas.Print(self.plotDir+"/"+self.name+".png", "png")

        return canvas

    def getCMSBox(self):
        par=self.canvasParams["CMS"]
        CMSbox       = ROOT.TLatex  (par["xpos"],par["ypos"], par["text"])
        CMSbox.SetNDC()
        CMSbox.SetTextSize(par["textSize"])
        CMSbox.SetTextFont(par["textFont"])
        CMSbox.SetTextColor(ROOT.kBlack)
        CMSbox.SetTextAlign(13) 
        return CMSbox
   
    def getExtraTextBox(self):
        par=self.canvasParams["TAG"]
        extraTextBox = ROOT.TLatex  (par["xpos"], par["ypos"], par["text"])
        extraTextBox.SetNDC()
        extraTextBox.SetTextSize(par["textSize"])
        extraTextBox.SetTextFont(par["textFont"])
        extraTextBox.SetTextColor(ROOT.kBlack)
        extraTextBox.SetTextAlign(13)
        
        return extraTextBox

    def getLumiBox(self):
        par=self.canvasParams["EnLumi"]
        lumibox = ROOT.TLatex  ( par["xpos"], par["ypos"], par["text"])
        lumibox.SetNDC()
        lumibox.SetTextAlign(31)
        lumibox.SetTextSize(par["textSize"])
        lumibox.SetTextFont(par["textFont"])
        lumibox.SetTextColor(ROOT.kBlack)
        
        return lumibox
        
    def setPlotStyle(self):
        ROOT.gROOT.SetStyle("Plain")
        ROOT.gStyle.SetOptStat()
        ROOT.gStyle.SetOptFit(0)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetFrameLineWidth(1)
        ROOT.gStyle.SetPadBottomMargin(0.13)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadTopMargin(0.06)
        ROOT.gStyle.SetPadRightMargin(0.05)

        ROOT.gStyle.SetLabelFont(42,"X")
        ROOT.gStyle.SetLabelFont(42,"Y")
        ROOT.gStyle.SetLabelSize(0.04,"X")
        ROOT.gStyle.SetLabelSize(0.04,"Y")
        ROOT.gStyle.SetLabelOffset(0.01,"Y")
        ROOT.gStyle.SetTickLength(0.02,"X")
        ROOT.gStyle.SetTickLength(0.02,"Y")
        ROOT.gStyle.SetLineWidth(1)
        ROOT.gStyle.SetTickLength(0.02 ,"Z")

        ROOT.gStyle.SetTitleSize(0.1)
        ROOT.gStyle.SetTitleFont(42,"X")
        ROOT.gStyle.SetTitleFont(42,"Y")
        ROOT.gStyle.SetTitleSize(0.05,"X")
        ROOT.gStyle.SetTitleSize(0.05,"Y")
        ROOT.gStyle.SetTitleOffset(1.1,"X")
        ROOT.gStyle.SetTitleOffset(1.4,"Y")
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPaintTextFormat("3.2f")
        ROOT.gROOT.ForceStyle()
