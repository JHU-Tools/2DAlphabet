import ROOT 
from ROOT import *

files = []
files.append(TFile.Open('../dataSidebandFewerbins_tt/MaxLikelihoodFitResult.root'))
files.append(TFile.Open('../dataSidebandFewerbins/MaxLikelihoodFitResult.root'))

workspaces = []
for f in files:
    workspaces.append(f.Get('MaxLikelihoodFitResult'))

out = open('nusiance_comparison.txt','w')

for syst in ['lumi','ttbar_xsec','st_t_xsec','st_tB_xsec','st_tW_xsec','st_tWB_xsec', 'topsf','wtagsf','Extrap','Tpt','Trig','Scale','Pileup','Scale','Pileup','pdf','jer','jes','jmr','jms']:
    line = syst + '\t'
    for w in workspaces:
        try:
            line+=str(w.var(syst).getValV()) + '\t'
        except:
            line+=str('\t-\t')

    out.write(line+'\n')

out.close()