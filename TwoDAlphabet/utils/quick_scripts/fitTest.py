import ROOT
from ROOT import *

rpf = TF2('rpf','( 0.033098 + 0.000276*x + 4.245e-06*x^2) + y*(-1e-6 - 1e-8*x - 1e-10*x^2)',50,350,700,4000)
print rpf.Eval(50,700)
rpf_hist = rpf.CreateHistogram()

rpf_hist.Draw('lego')
raw_input('waiti')

for xbin in range(1,rpf_hist.GetNbinsX()+1):
	for ybin in range(1,rpf_hist.GetNbinsY()+1):
		if rpf_hist.GetBinContent(xbin,ybin) <= 0:
			print str(xbin) + ' ' + str(ybin)

