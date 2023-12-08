import ROOT
from ROOT import *

f = TFile.Open('QCDMCwTT/base_QCDMCwTT.root')
w = f.Get('w_2D')

RAS_vars = w.allVars().selectByName('Fail_bin*',True)
iterator = RAS_vars.createIterator()

var = iterator.Next()

while var:
	var.Print()
	var = iterator.Next()
