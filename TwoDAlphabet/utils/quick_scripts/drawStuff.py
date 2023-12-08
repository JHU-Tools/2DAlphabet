import ROOT
from ROOT import *

def drawStuff(x,y,pdf):
	c = TCanvas('c_'+pdf.GetName(),'c_'+pdf.GetName(),600,500)
	c.cd()
	hist = pdf.createHistogram(pdf.GetName(),x,RooFit.YVar(y)).Draw('lego')
	

myf = TFile.Open('../QCDMCwTTfewerBins/MaxLikelihoodFitResult.root')
myw = myf.Get('MaxLikelihoodFitResult')

x = myw.var('jetmass')
y = myw.var('resmass')

myw.pdf('shapeBkg_fail_ttbar_morph').Print()

drawStuff(x,y,myw.pdf('model_s'))
drawStuff(x,y,myw.pdf('pdf_binfail'))
drawStuff(x,y,myw.pdf('pdf_binfail_nuis'))
drawStuff(x,y,myw.pdf('shapeBkg_fail_ttbar_morph'))
# drawStuff(x,y,myw.pdf('jer_Pdf'))
# drawStuff(x,y,myw.pdf('jmr_Pdf'))

raw_input('waiting')