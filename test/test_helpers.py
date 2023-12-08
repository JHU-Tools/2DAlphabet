
from argparse import ArgumentParser
from typing import Type
from ROOT import TH2F, RooArgList, RooRealVar, TH1F
import pytest
from TwoDAlphabet.helpers import is_filled_list, open_json, arg_dict_to_list, parse_arg_dict, make_RDH, dict_copy, nested_dict, roofit_form_to_TF1, set_hist_maximums

test_dict = {
    "NAME": "bare",
    "this": "is",
    "a" : {
        "test": 1,
        "dictionary": False
    },
    "OPTIONS": {},
    "GLOBAL": {}
}

test_dict_zeros = {
    'NAME': 0,
    'this': 0,
    'a' : {
        'test': 0,
        'dictionary': 0
    },
    "OPTIONS": {},
    "GLOBAL": {}
}

def test__open_json():
    assert (open_json("test/testjson.json") == test_dict)
    open_json('test/twoDtest.json')

def test__arg_dict_to_list():
    test = {"key": "value"}
    assert arg_dict_to_list(test) == ["key=value"]

def test__parse_arg_dict():
    parser = ArgumentParser(prog='TEST')
    parser.add_argument("key",type=str,nargs='?',default="notvalue")
    parser.add_argument("key2",type=str,nargs='?',default="notvalue2")
    test = {"key": "value"}
    args = parse_arg_dict(parser, test)
    assert args.key == "value"
    assert args.key2 == "notvalue2"

def test__make_RDH():
    var1 = RooRealVar('axis1','',2)
    var2 = RooRealVar('axis2','',2)
    h = TH2F('name','',10,0,10,10,0,10)
    make_RDH(h, RooArgList(var1,var2))
    assert True

def test__dict_copy():
    assert dict_copy(test_dict, structureOnly=False) == test_dict
    assert dict_copy(test_dict, structureOnly=True) == test_dict_zeros

def test__nested_dict():
    d = nested_dict(3, int)
    d[0] = 0
    d[1][1] = 1
    d[2][2][2] = 2
    assert True

def test__nested_dict_TypeError1():
    d = nested_dict(3, int)
    with pytest.raises(TypeError):
        d[0][0][0][0] = 2

def test__roofit_form_to_TF1():
    rfvform = '(@0*x**2 + @1*x)*(@2*y**3)+@2'
    assert roofit_form_to_TF1(rfvform, shift=0) == '([0]*x**2 + [1]*x)*([2]*y**3)+[2]'
    assert roofit_form_to_TF1(rfvform, shift=1) == '([1]*x**2 + [2]*x)*([3]*y**3)+[3]'

def test__set_hist_maximums():
    h1 = TH1F('h1','',10,0,10)
    h2 = h1.Clone('h2')
    h3 = h1.Clone('h3')
    h1.SetMaximum(100)
    for h in set_hist_maximums([h1,h2,h3], factor=1.1):
        assert int(h.GetMaximum()) == 110

def test__is_filled_list():
    assert is_filled_list({},'test_key') == False
    assert is_filled_list({'no_key':'no_val'},'test_key') == False
    with pytest.raises(TypeError):
        is_filled_list('not a dict','test_key')
    assert is_filled_list({'test_key':'test_val'},'test_key') == False
    assert is_filled_list({'test_key':[]},'test_key') == False
    assert is_filled_list({'test_key':['test_val']},'test_key') == True