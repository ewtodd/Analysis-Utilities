import ROOT
import numpy as np
import os


class PyPlottingUtils:
    def __init__(self):
        self.SetROOTStyle()
        self.colors = self.GetDefaultColors()

    def SetROOTStyle(self):
        """Apply ROOT style"""
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadRightMargin(0.05)
        ROOT.gStyle.SetPadTopMargin(0.12)
        ROOT.gStyle.SetPadBottomMargin(0.15)
        ROOT.gStyle.SetTitleSize(0.06, "XY")
        ROOT.gStyle.SetLabelSize(0.06, "XY")
        ROOT.gStyle.SetLegendFont(132)
        ROOT.gStyle.SetTitleOffset(1.2, "X")
        ROOT.gStyle.SetTitleOffset(1.2, "Y")
        ROOT.gStyle.SetTextFont(42)
        ROOT.gStyle.SetHistLineWidth(1)
        ROOT.gStyle.SetLineWidth(1)

        # Enable automatic grid with minor divisions
        ROOT.gStyle.SetPadGridX(1)
        ROOT.gStyle.SetPadGridY(1)
        ROOT.gStyle.SetGridStyle(3)  # Dotted
        ROOT.gStyle.SetGridWidth(1)
        ROOT.gStyle.SetGridColor(ROOT.kGray)

        # Key for minor grid lines in log scale
        ROOT.gStyle.SetPadTickX(2)
        ROOT.gStyle.SetPadTickY(2)

        ROOT.gROOT.ForceStyle(ROOT.kTRUE)

    def ConfigureHistogram(self, hist, color, title=""):
        if not hist:
            return

        hist.SetLineColor(color)
        hist.SetTitle(title)
        hist.SetFillColorAlpha(color, 0.2)

        # Key settings for minor divisions in log scale
        hist.GetYaxis().SetMoreLogLabels(ROOT.kFALSE)
        hist.GetYaxis().SetNoExponent(ROOT.kFALSE)
        hist.SetMinimum(10)

        # Division settings
        hist.GetYaxis().SetNdivisions(50109)
        hist.GetXaxis().SetNdivisions(506)
        hist.GetXaxis().SetTitleSize(0.06)
        hist.GetYaxis().SetTitleSize(0.06)
        hist.GetXaxis().SetLabelSize(0.06)
        hist.GetYaxis().SetLabelSize(0.06)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetYaxis().SetTitleOffset(1.2)

    def ConfigureCanvas(self, canvas, logy=False):
        if not canvas:
            return

        canvas.SetGridx(1)
        canvas.SetGridy(1)
        canvas.SetLogy(logy)
        canvas.SetTicks(1, 1)
        ROOT.gPad.SetTicks(1, 1)

    def GetDefaultColors(self):
        return [
            ROOT.kRed + 1,
            ROOT.kBlue + 1,
            ROOT.kGreen + 2,
            ROOT.kMagenta,
            ROOT.kOrange,
            ROOT.kCyan + 1,
            ROOT.kYellow + 2,
        ]

    def GetSourceColor(self, source_id):
        """Get color for source ID"""
        return self.colors[source_id % len(self.colors)]

    def CreateLegend(self, x1=0.6, y1=0.7, x2=0.9, y2=0.85):
        leg = ROOT.TLegend(x1, y1, x2, y2)
        leg.SetBorderSize(1)
        leg.SetFillColor(ROOT.kWhite)
        leg.SetTextSize(0.05)
        leg.SetTextFont(132)
        return leg

    def AddSubplotLabel(self, label, x=0.9, y=0.83):
        """Add subplot label"""
        latex = ROOT.TLatex(x, y, label)
        latex.SetNDC()
        latex.SetTextSize(0.06)
        latex.SetTextAlign(33)  # Right top alignment
        latex.Draw()
        return latex


# Initialize the plotting utilities
plot_utils = PyPlottingUtils()
