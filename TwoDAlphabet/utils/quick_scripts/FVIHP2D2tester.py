import ROOT
from ROOT import *

allVars = []

# ----- Get files and workspaces ------------------------------
new_f = TFile.Open('../dataSidebandFewerbins/MaxLikelihoodFitResult.root')
new_w = new_f.Get('MaxLikelihoodFitResult')

old_f = TFile.Open('../dataSidebandFewerbins/base_dataSidebandFewerbins.root')
old_w = old_f.Get('w_2D')

# ----- Get vars, pdfs, funcs, hists, etc ---------------------
xvar = new_w.var('jetmass')
yvar = new_w.var('resmass')
allVars.extend([xvar,yvar])

postfit_sidebands_pdf = new_w.pdf('shapeBkg_fail_ttbar_morph')
postfit_sidebands_norm = new_w.function('n_exp_final_binfail_proc_ttbar')
postfit_sidebands = postfit_sidebands_pdf.createHistogram('postfit_sidebands',xvar,RooFit.Binning(7,50,330),RooFit.YVar(yvar,RooFit.Binning(8,800,4000)))
postfit_sidebands.Scale(postfit_sidebands_norm.getValV())
allVars.extend([postfit_sidebands_pdf,postfit_sidebands_norm])

prefit_sidebands_RDH = old_w.data('ttbar_fail')
prefit_sidebands_RHP = RooHistPdf(prefit_sidebands_RDH.GetName()+'Pdf',prefit_sidebands_RDH.GetName()+'Pdf',RooArgSet(xvar,yvar),prefit_sidebands_RDH)
allVars.extend([prefit_sidebands_RDH,prefit_sidebands_RHP])

############################## Begin Morphed Shape Reconstruction ##############################
# ----- Get TList with pdfs and RooArgList with RooRealVars  for FVIHP2D2 construction ----
# RAL is ordered [a,b,c,...]
# TList is ordered [nom, a_hi, a_lo, b_hi, b_lo,...]
coefList = RooArgList()
pdfList = TList()
pdfList.Add(prefit_sidebands_RHP)
for n in ["Scale", "Tpt", "Trig", "pdf","Pileup", "jer", "jes", "jmr", "jms"]:
    coefList.add(new_w.var(n))
    thisRHPup = RooHistPdf('ttbar_fail_'+n+'UpPdf','ttbar_fail_'+n+'UpPdf',RooArgSet(xvar,yvar),old_w.data('ttbar_fail_'+n+'Up'))
    thisRHPdown = RooHistPdf('ttbar_fail_'+n+'DownPdf','ttbar_fail_'+n+'DownPdf',RooArgSet(xvar,yvar),old_w.data('ttbar_fail_'+n+'Down'))
    pdfList.Add(thisRHPup)
    pdfList.Add(thisRHPdown)
    allVars.extend([new_w.var(n),thisRHPup,thisRHPdown])

# Conditional should be FALSE (I checked this)
# I think it normalizes each x slice to 1 independently
# Last two arguments are the values used by Combine
recofit_sidebands_FVIHP = FastVerticalInterpHistPdf2D2('recofit_sidebands_FVIHP', 'postfit_sidebands', xvar, yvar, False, pdfList, coefList, 1, 0)
recofit_sidebands_final = recofit_sidebands_FVIHP.createHistogram('recofit_sidebands_final',xvar,RooFit.Binning(7,50,330),RooFit.YVar(yvar,RooFit.Binning(8,800,4000)))
allVars.append(recofit_sidebands_FVIHP)

# ---------- Recreate the normalizations -------------------------
reco_sidebands_normlist = RooArgList()
# 1 - From shape uncertainties using AsymPow
asympows = {}
for su in ["Scale", "Tpt", "Trig", "pdf", "Pileup", "jer", "jes", "jmr", "jms"]:
    kappaLowVal = old_w.data('ttbar_fail_'+su+'Down').sumEntries()/prefit_sidebands_RDH.sumEntries()
    kappaHighVal = old_w.data('ttbar_fail_'+su+'Up').sumEntries()/prefit_sidebands_RDH.sumEntries()
    if abs(kappaHighVal-1) < 1e-3 and abs(kappaLowVal-1) < 1e-3:    # Combine won't count the systematic unless this condition is satisfied so we shouldn't either
        kappaLow = RooConstVar('ttbar_fail_'+su+'_kappaLow','ttbar_fail_'+su+'_kappaLow',kappaLowVal)
        kappaHigh = RooConstVar('ttbar_fail_'+su+'_kappaHigh','ttbar_fail_'+su+'_kappaHigh',kappaHighVal)
        asympows[su] = AsymPow('systeff_fail_ttbar_'+su, 'systeff_fail_ttbar_'+su, kappaLow, kappaHigh, new_w.var(su))
        allVars.extend([kappaLow,kappaHigh])
        reco_sidebands_normlist.add(asympows[su])


# 2 - THIS STEP UNNECCESSARY - JUST GRAB IT 
# OLD - From lnN (sym and asym) using ProcessNormalization
# Initialize
recofit_sidebands_lnN_norm = ProcessNormalization('recofit_sidebands_lnN_norm','recofit_sidebands_lnN_norm',prefit_sidebands_RDH.sumEntries())
allVars.append(recofit_sidebands_lnN_norm)
# Symmetric
sym = {'lumi':1.026}
for s in sym.keys():
    recofit_sidebands_lnN_norm.addLogNormal(sym[s], new_w.var(s))
# Asymmetric
asym = {'topsf':{'up':1.14,'down': 1.056},'ttbar_xsec':{'up':1.024,'down':1.035}}
for a in asym.keys():
    recofit_sidebands_lnN_norm.addAsymmLogNormal(asym[a]['up'], asym[a]['down'], new_w.var(a))

reco_sidebands_normlist.add(recofit_sidebands_lnN_norm)
reco_sidebands_fullNorm = RooProduct('reco_sidebands_fullNorm','reco_sidebands_fullNorm',reco_sidebands_normlist)
allVars.append(reco_sidebands_fullNorm)

reco_sidebands_fullNorm.Print()
postfit_sidebands_norm.Print()
quit()

############################### End Reconstruction ##############################################
c = TCanvas('c','c',1400,700)
c.Divide(2,1)
c.cd(1)
postfit_sidebands.Draw('lego')
c.cd(2)
recofit_sidebands_final.Draw('lego')

c2 = TCanvas('c2','c2',800,700)
c2.cd()
diff = postfit_sidebands.Clone()
diff.Add(recofit_sidebands_final,-1)
diff.Draw('surf')
raw_input('waiting')
