import ROOT
from ROOT import *

file = TFile.Open('../fakeData_QCDMC_plus_TTMC.root')

hist_og = file.Get('MtwvMtPass')

xvar = RooRealVar('xvar','xvar',50,350)
yvar = RooRealVar('yvar','yvar',700,4000)
RAL_myvars = RooArgList(xvar,yvar)
RAS_myvars = RooArgSet(RAL_myvars)

RDH = RooDataHist('RDH','RDH',RAL_myvars,hist_og)

RHP = RooHistPdf('RHP','RHP',RAS_myvars,RDH)

hist_new = RHP.createHistogram('hist_new',xvar,RooFit.Binning(30,50,350),RooFit.YVar(yvar,RooFit.Binning(33,700,4000)))

print 'TH1: ' + str(hist_og.Integral())
print 'RDH (): ' + str(RDH.sum(False))
print 'RDH (true): ' + str(RDH.sum(True))
print 'RHP: ' + str(RHP.analyticalIntegral(0))
print 'hist_new: ' + str(hist_new.Integral())