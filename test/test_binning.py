import itertools
from TwoDAlphabet.binning import *
from ROOT import TH2F
import pytest
from copy import deepcopy

basedict = {
    "HELP": "The binning of the x and y axes should be configured here",
    "X": {
        "NAME": "xaxis",
        "TITLE": "xaxis",
        "MIN": 0,
        "MAX": 24,
        "NBINS": 12,
        "SIGSTART": 14,
        "SIGEND": 16
    },
    "Y": {
        "NAME": "yaxis",
        "TITLE": "yaxis",
        "MIN": 0, 
        "MAX": 20,
        "NBINS": 10
    }
}

def filled_template(h):
    out = h.Clone('filled_test')
    for x,y in itertools.product(range(1,h.GetNbinsX()+1),range(1,h.GetNbinsY()+1)):
        out.SetBinContent(x,y,1)
    return out

template = TH2F('test','test',12,0,24,10,0,20) # NOTE: If you change this, the tests will fail!!!
filled = filled_template(template)

def test_BinningBase():
    b = Binning('test',basedict,template)
    assert(b.xbinList == [0,2,4,6,8,10,12,14,16,18,20,22,24])
    assert(b.xbinByCat == {"LOW":[0,2,4,6,8,10,12,14],"SIG":[14,16],"HIGH":[16,18,20,22,24]})
    assert(b.ybinList == [0,2,4,6,8,10,12,14,16,18,20])

def test_BinningSubsetX():
    test_dict = deepcopy(basedict)
    test_dict['X']['MIN'] = 2
    test_dict['X']['MAX'] = 22
    test_dict['X']['NBINS'] = 10
    b = Binning('test',test_dict,template)
    assert(b.xbinList == [2,4,6,8,10,12,14,16,18,20,22])
    assert(b.xbinByCat == {"LOW":[2,4,6,8,10,12,14],"SIG":[14,16],"HIGH":[16,18,20,22]})

def test_BinningSubsetY():
    test_dict = deepcopy(basedict)
    test_dict['Y']['MIN'] = 2
    test_dict['Y']['MAX'] = 18
    test_dict['Y']['NBINS'] = 8
    b = Binning('test',test_dict,template)
    assert(b.ybinList == [2,4,6,8,10,12,14,16,18])

def test_BinningVariableX():
    test_dict = deepcopy(basedict)
    test_dict['X'].pop('MIN')
    test_dict['X'].pop('MAX')
    test_dict['X'].pop('NBINS')
    test_dict['X']['BINS'] = [0,2,4,6,8,10,12,14,16,18,20,22,24]
    b = Binning('test',test_dict,template)
    assert(b.xbinList == [0,2,4,6,8,10,12,14,16,18,20,22,24])
    assert(b.xbinByCat == {"LOW":[0,2,4,6,8,10,12,14],"SIG":[14,16],"HIGH":[16,18,20,22,24]})

def test_BinningVariableY():
    test_dict = deepcopy(basedict)
    test_dict['Y'].pop('MIN')
    test_dict['Y'].pop('MAX')
    test_dict['Y'].pop('NBINS')
    test_dict['Y']['BINS'] = [0,2,4,6,8,10,12,14,16,18,20]
    b = Binning('test',basedict,template)
    assert(b.ybinList == [0,2,4,6,8,10,12,14,16,18,20])

def test_BinningBreakLowerX():
    test_dict = deepcopy(basedict)
    test_dict['X']['MIN'] = -2

    with pytest.raises(ValueError):
        b = Binning('test',test_dict,template)

def test_BinningBreakLowerY():
    test_dict = deepcopy(basedict)
    test_dict['Y']['MIN'] = -2

    with pytest.raises(ValueError):
        b = Binning('test',test_dict,template)

def test_BinningBreakUpperX():
    test_dict = deepcopy(basedict)
    test_dict['X']['MAX'] = 26

    with pytest.raises(ValueError):
        b = Binning('test',test_dict,template)

def test_BinningBreakUpperY():
    test_dict = deepcopy(basedict)
    test_dict['Y']['MAX'] = 22

    with pytest.raises(ValueError):
        b = Binning('test',test_dict,template)

def test_BinningOrderX():
    test_dict = deepcopy(basedict)
    del test_dict['X']['MIN']
    del test_dict['X']['MAX']
    del test_dict['X']['NBINS']
    test_dict['X']['BINS'] = [0,2,4,6,8,10,16,14,18,20]

    with pytest.raises(ValueError):
        b = Binning('test',test_dict,template)

def test_BinningOrderY():
    test_dict = deepcopy(basedict)
    del test_dict['Y']['MIN']
    del test_dict['Y']['MAX']
    del test_dict['Y']['NBINS']
    test_dict['Y']['BINS'] = [0,2,4,6,8,10,16,14,18,20]

    with pytest.raises(ValueError):
        b = Binning('test',test_dict,template)

def test__concat_bin_lists():
    assert(concat_bin_lists([[0,1,2],[2,3]]) == [0,1,2,3])

def test__concat_bin_lists_VALUE():
    with pytest.raises(ValueError):
        concat_bin_lists([[0,1,2],[3,4]])

def test__histlist_to_binlist():
    assert (histlist_to_binlist('X',[template]) == [0,2,4,6,8,10,12,14,16,18,20,22,24])
    assert (histlist_to_binlist('Y',[template]) == [0,2,4,6,8,10,12,14,16,18,20])

def test_SplittingHists():
    low = copy_hist_with_new_bins('low','X',filled,[0,2,4,6,8])
    sig = copy_hist_with_new_bins('sig','X',filled,[8,10,12])
    high = copy_hist_with_new_bins('high','X',filled,[12,14,16,18,20,22,24])
    test = stitch_hists_in_x(filled,[low,sig,high])
    assert(low.GetNbinsX() == 4)
    assert(sig.GetNbinsX() == 2)
    assert(high.GetNbinsX() == 6)
    assert(test.GetNbinsX() == filled.GetNbinsX())

def test_SplittingHists_VALUE():
    low = copy_hist_with_new_bins('low','X',filled,[2,4,6])
    sig = copy_hist_with_new_bins('sig','X',filled,[8,10,12])
    high = copy_hist_with_new_bins('high','X',filled,[14,16,18,20])
    with pytest.raises(ValueError):
        test = stitch_hists_in_x(filled,[low,sig,high])

def test__make_blinded_hist():
    nbins_start = filled.GetNbinsX()
    h = make_blinded_hist(filled,[8,10])
    assert(h.Integral() == filled.Integral()*(nbins_start-1)/nbins_start)

def test__make_blinded_hist_VALUE():
    with pytest.raises(ValueError):
        make_blinded_hist(filled,[9,10])
    with pytest.raises(ValueError):
        make_blinded_hist(filled,[20,26])

def test__make_blinded_hist_INDEX():
    with pytest.raises(IndexError):
        make_blinded_hist(filled,[8,10,12])

def test__copy_hist_with_new_bins():
    h = copy_hist_with_new_bins('test','Y',filled,[2,4,6,8])
    assert (h.Integral() == 36)
    h = copy_hist_with_new_bins('test','X',filled,[2,4,6,8])
    assert (h.Integral() == 30)

def test__copy_hist_with_new_bins_VALUE1():
    with pytest.raises(ValueError):
        h = copy_hist_with_new_bins('test','Z',template,[2,4,6,8])

def test__copy_hist_with_new_bins_VALUE2():
    with pytest.raises(ValueError):
        h = copy_hist_with_new_bins('test','Y',template,[3,4,6,8])

def test__get_min_bin_width():
    assert(get_min_bin_width(template.ProjectionY()) == 2)

def test__get_min_bin_width_VALUE():
    with pytest.raises(TypeError):
        get_min_bin_width(template) == 2

def test__convert_to_events_per_unit():
    h = convert_to_events_per_unit(filled.ProjectionY())
    assert(h.Integral() == filled.Integral())

def test__zero_negative_bins():
    clone = filled.Clone()
    clone.SetBinContent(1,1,-10)
    h = zero_negative_bins('test',clone)
    assert(h.Integral() == filled.Integral()-1)

def test__remap_binlist_to_unity():
    assert (remap_binlist([0,1,2,3,4],0,4) == [0,1,2,3,4])
    assert (remap_binlist([0,1,2,3,4],0,1) == [0,0.25,0.5,0.75,1])
    assert (remap_binlist([0,1,2,3,4],-1,1) == [-1,-0.5,0,0.5,1])
    assert (remap_binlist([0,1,2,3,4],1,2) == [1,1.25,1.5,1.75,2.0])

def test__remap_to_unity():
    h = remap_hist_axis(template)
    assert(h.GetNbinsX() == template.GetNbinsX())
    assert(h.GetNbinsY() == template.GetNbinsY())
    assert(h.GetXaxis().GetXmin() == 0)
    assert(h.GetXaxis().GetXmax() == 1)
    assert(h.GetYaxis().GetXmin() == 0)
    assert(h.GetYaxis().GetXmax() == 1)