from TwoDAlphabet import plot
from TwoDAlphabet.twoDalphabet import MakeCard, TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction
import os
import ROOT
import json

'''--------------------------Helper functions---------------------------'''
def _select_signal(row, args):
    signame = args[0]
    poly_order = args[1]
    if row.process_type == 'SIGNAL':
        if signame in row.process:
            return True
        else:
            return False
    elif 'Background_' in row.process:
        if row.process == 'Background_'+poly_order:
            return True
        else:
            return False
    else:
        return True

def _generate_constraints(nparams):
    out = {}
    for i in range(nparams):
        if i == 0:
            out[i] = {"MIN":0,"MAX":1}
        else:
            out[i] = {"MIN":-5,"MAX":5}
    return out

_rpf_options = {
    '0x0': {
        'form': '0.1*(@0)',
        'constraints': _generate_constraints(1)
    },
    '1x0': {
        'form': '0.1*(@0+@1*x)',
        'constraints': _generate_constraints(2)
    },
    '0x1': {
        'form': '0.1*(@0+@1*y)',
        'constraints': _generate_constraints(2)
    },
    '1x1': {
        'form': '0.1*(@0+@1*x)*(1+@2*y)',
        'constraints': _generate_constraints(3)
    }
}

def test_make():
    twoD = TwoDAlphabet('CR_cicd', 'twoDtest_cicd.json', loadPrevious=False)
    qcd_hists = twoD.InitQCDHists()
    f = "CR_fail"
    p = "CR_pass"
    binning_f, _ = twoD.GetBinningFor(f)
    fail_name = 'Background_'+f
    qcd_f = BinnedDistribution(fail_name, qcd_hists[f],binning_f, constant=False)
    twoD.AddAlphaObj('Background',f,qcd_f)
    for opt_name, opt in _rpf_options.items():
        qcd_rpf = ParametricFunction(fail_name.replace('fail','rpf')+'_'+opt_name,binning_f, opt['form'],constraints=opt['constraints'])
        qcd_p = qcd_f.Multiply(fail_name.replace('fail','pass')+'_'+opt_name, qcd_rpf)
        twoD.AddAlphaObj('Background_'+opt_name,p,qcd_p,title='Background')
    twoD.Save()

def test_fit():
    working_area = 'CR_cicd'
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, 'TprimeB-1800_16', '1x1')
    twoD.MakeCard(subset, 'TprimeB-1800_16_area')
    twoD.MLfit('TprimeB-1800_16_area',rMin=-1,rMax=20,verbosity=0)

def test_plot():
    working_area = 'CR_cicd'
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, 'TprimeB-1800_16', '1x1')
    twoD.StdPlots('TprimeB-1800_16_area', subset)

def test_limit():
    poly_order = '1x1'
    working_area = 'CR_cicd'
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)
    signame = "TprimeB-1800_16"
    subset = twoD.ledger.select(_select_signal, signame, poly_order)
    twoD.MakeCard(subset, signame+'_area')
    twoD.Limit(subtag=signame+'_area',blindData=False,verbosity=0,condor=False)

def extract_limits():
    f = ROOT.TFile.Open("CR_cicd/TprimeB-1800_16_area/higgsCombineTest.AsymptoticLimits.mH120.root")
    t = f.Get("limit")

    limits = []
    for entry in t:
        limits.append(entry.limit)

    if len(limits) != 6:
        print(f"WARNING: Expected 6 limit entries, got {len(limits)}")

    out = {
        "exp_m2": limits[0],
        "exp_m1": limits[1],
        "exp":    limits[2],
        "exp_p1": limits[3],
        "exp_p2": limits[4],
        "obs":    limits[5],
    }

    with open("limits_cicd.json", "w") as jf:
        json.dump(out, jf, indent=2)

    return out

def test_compare_limits(tolerance=0.001):
    with open("limits_cicd_reference.json") as f:
        ref = json.load(f)

    with open("limits_cicd.json") as f:
        new = json.load(f)

    print("\n### Comparing limits to reference ###")
    fail = False

    for key in ref:
        r = ref[key]
        p = new[key]

        if r == 0:
            print(f"{key}: reference=0, skipping")
            continue

        rel = abs(p - r) / abs(r)

        print(f"{key}: reference={r:.4f}, new={p:.4f}, rel.diff.={rel*100:.2f}%")

        if rel > tolerance:
            print(f"  -> exceeds tolerance ({tolerance*100:.1f}%)")
            fail = True

    if fail:
        raise AssertionError("Limit comparison failed: differences exceed tolerance")

    print("All limits within tolerance\n")


if __name__ == '__main__':
    test_make()
    test_fit()
    test_plot()
    test_limit()
    extract_limits()
    test_compare_limits()