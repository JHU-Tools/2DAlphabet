import ROOT
from ROOT import *

myf = TFile.Open('../MaxLikelihoodFitResult.root')
myw = myf.Get('MaxLikelihoodFitResult')

myf1 = TFile.Open('../../2DCombineFitCode/bstar/input_hists/rootfiles/TW2Dalphabetweightedttbar_Trigger_nominal_none_PSET_default.root')
myh1 = myf1.Get('MtwvMtFail')
print myh1.Integral()

mything1 = myw.function('n_exp_final_binfail_proc_ttbar')
mything2 = myw.function('n_exp_binfail_proc_ttbar') 
mything3 = myw.function('shapeBkg_ttbar_fail__norm')
mything4 = myw.function('ttbar_norm')
mything5 = myw.function('ttbar_fail_norm_start')

mything1.Print()
mything2.Print()
mything3.Print()
mything4.Print()
mything5.Print()
