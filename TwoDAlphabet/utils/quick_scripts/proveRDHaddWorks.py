#####################################################################################################
# This script is meant to test whether adding RDH_1 to RDH_2 and then rescaling RDH-1 changes RDH-2 #
#####################################################################################################

import ROOT
from ROOT import *

gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)

file1 = TFile.Open('../../2DCombineFitCode/bstar/input_hists/rootfiles/TW2Dalphabetweightedttbar_Trigger_nominal_none_PSET_default.root')
file2 = TFile.Open('../../2DCombineFitCode/bstar/input_hists/rootfiles/TW2DalphabetQCD_Trigger_nominal_none_PSET_default.root')

TH2_1 = file1.Get('MtwvMtPass')
TH2_2 = file2.Get('MtwvMtPass')

xvar = RooRealVar('xvar','xvar',50,350)
yvar = RooRealVar('yvar','yvar',700,4000)

RAL_vars = RooArgList(xvar,yvar)

name1 = TH2_1.GetName()
RDH1 = RooDataHist(name1,name1,RAL_vars,TH2_1)

name2 = TH2_2.GetName()
RDH2 = RooDataHist(name2,name2,RAL_vars,TH2_2)

RDH2.add(RDH1)

RDH2_before_scale = RDH2.createHistogram('RDH2_before_scale',xvar,RooFit.Binning(30,50,350),RooFit.YVar(yvar,RooFit.Binning(33,700,4000)))
RDH1_before_scale = RDH1.createHistogram('RDH1_before_scale',xvar,RooFit.Binning(30,50,350),RooFit.YVar(yvar,RooFit.Binning(33,700,4000)))

print RDH1.isWeighted()
quit()

RDH2_after_scale = RDH2.createHistogram('RDH2_after_scale',xvar,RooFit.Binning(30,50,350),RooFit.YVar(yvar,RooFit.Binning(33,700,4000)))
RDH1_after_scale = RDH1.createHistogram('RDH1_after_scale',xvar,RooFit.Binning(30,50,350),RooFit.YVar(yvar,RooFit.Binning(33,700,4000)))

xCan = TCanvas('x','x',1400,1400)
xCan.Divide(2,2)
xCan.cd(1)
RDH1_before_scale.Draw('lego')
xCan.cd(2)
RDH2_before_scale.Draw('lego')
xCan.cd(3)
RDH1_after_scale.Draw('lego')
xCan.cd(4)
RDH2_after_scale.Draw('lego')

xCan.Print('proveRDHaddWorks.pdf','pdf')