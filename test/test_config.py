import pytest
from numpy import nan
from TwoDAlphabet.config import *
from TwoDAlphabet.config import _keyword_replace, _get_syst_attrs, _df_sanity_checks

class TestConfig():
    @classmethod
    def setup_class(cls):
        cls.objBase = Config('test/twoDtest.json','test/unit_test',
                             findreplace={},
                             externalOptions={}
                             )
        cls.objBare = Config('test/twoDtest.json','test/unit_test',
                             findreplace={"ttbar_16":"THIS","lumi":"IS"},
                             externalOptions={"blindedPlots":True}
                             )
        cls.objBase.Construct()

    def test__init(self):
        assert "THIS" in self.objBare.config["PROCESSES"].keys()
        assert "IS" in self.objBare.config["SYSTEMATICS"].keys()
        assert self.objBare.config["OPTIONS"]["blindedPlots"] == True
        assert self.objBare.options.blindedPlots == True
        assert self.objBare.options.blindedFit == False

    def test_Construct(self):
        assert os.path.exists('test/unit_test/CR/organized_hists.root')
        objBase2 = Config('test/twoDtest.json', 'test/unit_test', findreplace={}, externalOptions={}, loadPrevious=True)
        objBase2.Construct()
        assert objBase2.organized_hists.file.GetName() == 'test/unit_test/CR/organized_hists.root'

    def test_Section(self):
        assert set(self.objBase.Section("GLOBAL").keys()) == set(['path','FILE','FILE_UP','FILE_DOWN','HIST','HIST_UP','HIST_DOWN'])

    def test__addFindReplace(self):
        with pytest.raises(ValueError):
            self.objBare._addFindReplace({'path':'anotherPath/'})

    def test__varReplacement(self):
        pass # Done in init

    def test_GetOptions(self):
        assert self.objBase.GetOptions(externalOptions={'freezeFail':True}).freezeFail == True

    def test_SaveOut(self):
        self.objBase.SaveOut()
        assert os.path.exists('test/unit_test/CR/runConfig.json')

    def test_FullTable(self):
        nan = -1000
        df = self.objBase.df.sort_values(by=['process','region']).reset_index().drop(columns='index').fillna(-1000).to_dict()
        assert df == {
            'direction': {
                0: nan, 1: nan, 2: nan, 3: 'Up', 4: 'Down', 5: nan, 6: nan, 7: 'Up',
                8: 'Down', 9: nan, 10: nan, 11: nan, 12: 'Up', 13: 'Down', 14: nan, 
                15: nan, 16: nan, 17: 'Up', 18: 'Down', 19: nan},
            'scale': {
                0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0,
                8: 1.0, 9: 1.0, 10: 1.0, 11: 1.0, 12: 1.0, 13: 1.0, 14: 1.0,
                15: 1.0, 16: 1.0, 17: 1.0, 18: 1.0, 19: 1.0},
            'shape_sigma': {
                0: nan, 1: nan, 2: nan, 3: 1.0, 4: 1.0, 5: nan, 6: nan, 7: 1.0,
                8: 1.0, 9: nan, 10: nan, 11: nan, 12: 1.0, 13: 1.0, 14: nan,
                15: nan, 16: nan, 17: 1.0, 18: 1.0, 19: nan},
            'lnN_asym': {
                0: nan, 1: nan, 2: nan, 3: nan, 4: nan, 5: nan, 6: nan, 7: nan,
                8: nan, 9: nan, 10: nan, 11: [1.2, 0.8], 12: nan, 13: nan, 14: nan,
                15: nan, 16: [1.2, 0.8], 17: nan, 18: nan, 19: nan},
            'process': {
                0: 'Data_Run2', 1: 'Data_Run2', 2: 'TprimeB-1800_16', 3: 'TprimeB-1800_16', 4: 'TprimeB-1800_16',
                5: 'TprimeB-1800_16', 6: 'TprimeB-1800_16', 7: 'TprimeB-1800_16',
                8: 'TprimeB-1800_16', 9: 'TprimeB-1800_16', 10: 'ttbar_16', 11: 'ttbar_16', 12: 'ttbar_16',
                13: 'ttbar_16', 14: 'ttbar_16', 15: 'ttbar_16', 16: 'ttbar_16', 17: 'ttbar_16', 18: 'ttbar_16', 19: 'ttbar_16'},
            'syst_type': {
                0: nan, 1: nan, 2: 'lnN', 3: 'shapes', 4: 'shapes', 5: nan, 6: 'lnN', 7: 'shapes',
                8: 'shapes', 9: nan, 10: 'lnN', 11: 'lnN_asym', 12: 'shapes', 13: 'shapes', 14: nan,
                15: 'lnN', 16: 'lnN_asym', 17: 'shapes', 18: 'shapes', 19: nan},
            'region': {
                0: 'CR_fail', 1: 'CR_pass', 2: 'CR_fail', 3: 'CR_fail', 4: 'CR_fail', 5: 'CR_fail', 6: 'CR_pass', 7: 'CR_pass',
                8: 'CR_pass', 9: 'CR_pass', 10: 'CR_fail', 11: 'CR_fail', 12: 'CR_fail', 13: 'CR_fail', 14: 'CR_fail',
                15: 'CR_pass', 16: 'CR_pass', 17: 'CR_pass', 18: 'CR_pass', 19: 'CR_pass'},
            'variation': {
                0: 'nominal', 1: 'nominal', 2: 'lumi', 3: 'JER', 4: 'JER', 5: 'nominal', 6: 'lumi', 7: 'JER',
                8: 'JER', 9: 'nominal', 10: 'lumi', 11: 'TT_xsec', 12: 'TptReweight', 13: 'TptReweight', 14: 'nominal',
                15: 'lumi', 16: 'TT_xsec', 17: 'TptReweight', 18: 'TptReweight', 19: 'nominal'},
            'color': {
                0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
                8: 0, 9: 0, 10: 2, 11: 2, 12: 2, 13: 2, 14: 2, 15: 2, 16: 2, 17: 2, 18: 2, 19: 2},
            'lnN': {
                0: nan, 1: nan, 2: 1.018, 3: nan, 4: nan, 5: nan, 6: 1.018, 7: nan,
                8: nan, 9: nan, 10: 1.018, 11: nan, 12: nan, 13: nan, 14: nan, 15: 1.018, 16: nan, 17: nan, 18: nan, 19: nan},
            'process_type': {
                0: 'DATA', 1: 'DATA', 2: 'SIGNAL', 3: 'SIGNAL', 4: 'SIGNAL', 5: 'SIGNAL', 6: 'SIGNAL', 7: 'SIGNAL',
                8: 'SIGNAL', 9: 'SIGNAL', 10: 'BKG', 11: 'BKG', 12: 'BKG', 13: 'BKG', 14: 'BKG', 15: 'BKG', 16: 'BKG', 17: 'BKG', 18: 'BKG', 19: 'BKG'},
            'source_filename': {
                0: 'test/data/THselection_Data_Run2.root', 1: 'test/data/THselection_Data_Run2.root',
                2: 'test/data/THselection_TprimeB-1800_16.root', 3: 'test/data/THselection_TprimeB-1800_16_JER_up.root',
                4: 'test/data/THselection_TprimeB-1800_16_JER_down.root', 5: 'test/data/THselection_TprimeB-1800_16.root',
                6: 'test/data/THselection_TprimeB-1800_16.root', 7: 'test/data/THselection_TprimeB-1800_16_JER_up.root',
                8: 'test/data/THselection_TprimeB-1800_16_JER_down.root', 9: 'test/data/THselection_TprimeB-1800_16.root',
                10: 'test/data/THselection_ttbar_16.root', 11: 'test/data/THselection_ttbar_16.root', 
                12: 'test/data/THselection_ttbar_16.root', 13: 'test/data/THselection_ttbar_16.root',
                14: 'test/data/THselection_ttbar_16.root', 15: 'test/data/THselection_ttbar_16.root',
                16: 'test/data/THselection_ttbar_16.root', 17: 'test/data/THselection_ttbar_16.root',
                18: 'test/data/THselection_ttbar_16.root', 19: 'test/data/THselection_ttbar_16.root'},
            'source_histname': {
                0: 'MthvMh_particleNet_CR_fail__nominal', 1: 'MthvMh_particleNet_CR_pass__nominal',
                2: 'MthvMh_particleNet_CR_fail__nominal', 3: 'MthvMh_particleNet_CR_fail__nominal',
                4: 'MthvMh_particleNet_CR_fail__nominal', 5: 'MthvMh_particleNet_CR_fail__nominal',
                6: 'MthvMh_particleNet_CR_pass__nominal', 7: 'MthvMh_particleNet_CR_pass__nominal',
                8: 'MthvMh_particleNet_CR_pass__nominal', 9: 'MthvMh_particleNet_CR_pass__nominal',
                10: 'MthvMh_particleNet_CR_fail__nominal', 11: 'MthvMh_particleNet_CR_fail__nominal',
                12: 'MthvMh_particleNet_CR_fail__TptReweight_up', 13: 'MthvMh_particleNet_CR_fail__TptReweight_down',
                14: 'MthvMh_particleNet_CR_fail__nominal', 15: 'MthvMh_particleNet_CR_pass__nominal',
                16: 'MthvMh_particleNet_CR_pass__nominal', 17: 'MthvMh_particleNet_CR_pass__TptReweight_up',
                18: 'MthvMh_particleNet_CR_pass__TptReweight_down', 19: 'MthvMh_particleNet_CR_pass__nominal'}
            }

    def test__regionTable(self):
        c = Config('test/twoDtest.json','test/unit_test')
        c.config['REGIONS']['CR_pass'].append('FAKEPROC')
        with pytest.raises(RuntimeError):
            c.FullTable()
    
    def test__data_not_included(self):
        c = Config('test/twoDtest.json','test/unit_test')
        del c.config['PROCESSES']['Data_Run2']
        with pytest.raises(RuntimeError):
            c.FullTable()
        c = Config('test/twoDtest.json','test/unit_test')
        c.config['PROCESSES']['DataFAKE'] = {
            "SYSTEMATICS":[],
            "SCALE": 1.0,
            "COLOR": 0,
            "TYPE": "DATA",
            "ALIAS": "data_obs",
            "LOC": "path/FILE:HIST"
        }
        with pytest.raises(RuntimeError):
            c.FullTable()

    def test__processTable(self):
        pass # tested in test__FullTable

    def test__systematicsTable(self):
        pass # tested in test__FullTable

    def test_GetHistMap(self):
        test_dict = {}
        for k,m in self.objBase.GetHistMap().items():
            test_dict[k] = m.to_dict()
            
        assert test_dict == {
            'test/data/THselection_TprimeB-1800_16.root': {
                'source_histname': {8: 'MthvMh_particleNet_CR_fail__nominal', 18: 'MthvMh_particleNet_CR_pass__nominal'},
                'color': {8: 0, 18: 0},
                'out_histname': {8: 'TprimeB-1800_16_CR_fail_nominal_FULL', 18: 'TprimeB-1800_16_CR_pass_nominal_FULL'},
                'scale': {8: 1.0, 18: 1.0}
            }, 'test/data/THselection_TprimeB-1800_16_JER_up.root': {
                'source_histname': {16: 'MthvMh_particleNet_CR_pass__nominal', 6: 'MthvMh_particleNet_CR_fail__nominal'},
                'color': {16: 0, 6: 0},
                'out_histname': {16: 'TprimeB-1800_16_CR_pass_JERUp_FULL', 6: 'TprimeB-1800_16_CR_fail_JERUp_FULL'},
                'scale': {16: 1.0, 6: 1.0}
            }, 'test/data/THselection_Data_Run2.root': {
                'source_histname': {9: 'MthvMh_particleNet_CR_fail__nominal', 19: 'MthvMh_particleNet_CR_pass__nominal'},
                'color': {9: 0, 19: 0},
                'out_histname': {9: 'Data_Run2_CR_fail_nominal_FULL', 19: 'Data_Run2_CR_pass_nominal_FULL'},
                'scale': {9: 1.0, 19: 1.0}
            }, 'test/data/THselection_ttbar_16.root': {
                'source_histname': {2: 'MthvMh_particleNet_CR_fail__TptReweight_up', 3: 'MthvMh_particleNet_CR_fail__TptReweight_down', 4: 'MthvMh_particleNet_CR_fail__nominal', 12: 'MthvMh_particleNet_CR_pass__TptReweight_up', 13: 'MthvMh_particleNet_CR_pass__TptReweight_down', 14: 'MthvMh_particleNet_CR_pass__nominal'},
                'color': {2: 2, 3: 2, 4: 2, 12: 2, 13: 2, 14: 2},
                'out_histname': {2: 'ttbar_16_CR_fail_TptReweightUp_FULL', 3: 'ttbar_16_CR_fail_TptReweightDown_FULL', 4: 'ttbar_16_CR_fail_nominal_FULL', 12: 'ttbar_16_CR_pass_TptReweightUp_FULL', 13: 'ttbar_16_CR_pass_TptReweightDown_FULL', 14: 'ttbar_16_CR_pass_nominal_FULL'},
                'scale': {2: 1.0, 3: 1.0, 4: 1.0, 12: 1.0, 13: 1.0, 14: 1.0}
            }, 'test/data/THselection_TprimeB-1800_16_JER_down.root': {
                'source_histname': {17: 'MthvMh_particleNet_CR_pass__nominal', 7: 'MthvMh_particleNet_CR_fail__nominal'},
                'color': {17: 0, 7: 0},
                'out_histname': {17: 'TprimeB-1800_16_CR_pass_JERDown_FULL', 7: 'TprimeB-1800_16_CR_fail_JERDown_FULL'},
                'scale': {17: 1.0, 7: 1.0}
            }
        }

class TestOrganizedHist():
    @classmethod
    def setup_class(cls):
        cls.c = Config('test/twoDtest.json','test/unit_test',
                             findreplace={},
                             externalOptions={}
                             )
        cls.c.Construct()
        cls.objBase = OrganizedHists(cls.c)

    def test__init(self):
        assert self.objBase.name == self.c.name
        assert self.objBase.filename == self.c.projPath + 'organized_hists.root'
        assert self.objBase.file.GetName() == self.objBase.filename

    def test_Add(self):
        self.objBase.Add()

    def test_Get(self):
        assert self.objBase.Get(process='Data_Run2',region='CR_pass',systematic='nominal') != None
        assert self.objBase.Get(histname='FAKE') == None

    def test_CreateSubRegions(self):
        h = self.objBase.Get(process='Data_Run2',region='CR_pass',systematic='nominal')
        self.objBase.CreateSubRegions(h)
        assert self.objBase.Get(process='Data_Run2',region='CR_pass',systematic='nominal',subspace='LOW') != None
        assert self.objBase.Get(process='Data_Run2',region='CR_pass',systematic='nominal',subspace='SIG') != None
        assert self.objBase.Get(process='Data_Run2',region='CR_pass',systematic='nominal',subspace='HIGH') != None
        with pytest.raises(NameError):
            assert self.objBase.Get(process='Data_Run2',region='CR_pass',systematic='nominal',subspace='TEST') != None

def test__keyword_replace():
    d = {'process':['ttbar'],
         'region':['SR'],
         'variation':['nominal'],
         'TEST':['$process_$region_$syst']
        }
    df = pandas.DataFrame(d)
    df = _keyword_replace(df,['TEST'])
    assert df.to_dict() == {
        'process': {0: 'ttbar'},
        'region': {0: 'SR'},
        'variation': {0: 'nominal'},
        'TEST': {0: 'ttbar_SR_nominal'},
    }

def test__get_syst_attrs():
    d = {"lumi": {
        "VAL": 1.018
    },
    "TT_xsec": {
        "VALUP": 1.20,
        "VALDOWN": 0.8
    },
    "TptReweight":{
        "ALIAS": "Top3",
        "UP": "path/FILE:HIST_UP",
        "DOWN": "path/FILE:HIST_DOWN",
        "SIGMA": 1.0
    },
    "JER": {
        "ALIAS": "JER",
        "UP": "path/FILE_UP:HIST",
        "DOWN": "path/FILE_DOWN:HIST",
        "SIGMA": 1.0
    },
    "FAKE": {}}
    assert _get_syst_attrs('lumi',d['lumi'])[0].at['lnN'] == 1.018
    assert pandas.isna(_get_syst_attrs('lumi',d['lumi'])[0].at['lnN_asym'])
    assert _get_syst_attrs('lumi',d['lumi'])[0].at['syst_type'] == 'lnN'
    assert pandas.isna(_get_syst_attrs('lumi',d['lumi'])[0].at['shape_sigma'])
    assert pandas.isna(_get_syst_attrs('lumi',d['lumi'])[0].at['source_filename'])
    assert pandas.isna(_get_syst_attrs('lumi',d['lumi'])[0].at['source_histname'])
    assert pandas.isna(_get_syst_attrs('lumi',d['lumi'])[0].at['direction'])

    assert pandas.isna(_get_syst_attrs('TT_xsec',d['TT_xsec'])[0].at['lnN'])
    assert _get_syst_attrs('TT_xsec',d['TT_xsec'])[0].at['lnN_asym'] == [1.20, 0.8]
    assert _get_syst_attrs('TT_xsec',d['TT_xsec'])[0].at['syst_type'] == 'lnN_asym'
    assert pandas.isna(_get_syst_attrs('TT_xsec',d['TT_xsec'])[0].at['shape_sigma'])
    assert pandas.isna(_get_syst_attrs('TT_xsec',d['TT_xsec'])[0].at['source_filename'])
    assert pandas.isna(_get_syst_attrs('TT_xsec',d['TT_xsec'])[0].at['source_histname'])
    assert pandas.isna(_get_syst_attrs('TT_xsec',d['TT_xsec'])[0].at['direction'])

    assert pandas.isna(_get_syst_attrs('TptReweight',d['TptReweight'])[0].at['lnN'])
    assert pandas.isna(_get_syst_attrs('TptReweight',d['TptReweight'])[0].at['lnN_asym'])
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[0].at['syst_type'] == 'shapes'
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[0].at['shape_sigma'] == 1.0
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[0].at['source_filename'] == 'path/FILE'
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[0].at['source_histname'] == 'HIST_UP'
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[0].at['direction'] == 'Up'

    assert pandas.isna(_get_syst_attrs('TptReweight',d['TptReweight'])[1].at['lnN'])
    assert pandas.isna(_get_syst_attrs('TptReweight',d['TptReweight'])[1].at['lnN_asym'])
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[1].at['syst_type'] == 'shapes'
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[1].at['shape_sigma'] == 1.0
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[1].at['source_filename'] == 'path/FILE'
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[1].at['source_histname'] == 'HIST_DOWN'
    assert _get_syst_attrs('TptReweight',d['TptReweight'])[1].at['direction'] == 'Down'

    assert pandas.isna(_get_syst_attrs('JER',d['JER'])[0].at['lnN'])
    assert pandas.isna(_get_syst_attrs('JER',d['JER'])[0].at['lnN_asym'])
    assert _get_syst_attrs('JER',d['JER'])[0].at['syst_type'] == 'shapes'
    assert _get_syst_attrs('JER',d['JER'])[0].at['shape_sigma'] == 1.0
    assert _get_syst_attrs('JER',d['JER'])[0].at['source_filename'] == 'path/FILE_UP'
    assert _get_syst_attrs('JER',d['JER'])[0].at['source_histname'] == 'HIST'
    assert _get_syst_attrs('JER',d['JER'])[0].at['direction'] == 'Up'

    assert pandas.isna(_get_syst_attrs('JER',d['JER'])[1].at['lnN'])
    assert pandas.isna(_get_syst_attrs('JER',d['JER'])[1].at['lnN_asym'])
    assert _get_syst_attrs('JER',d['JER'])[1].at['syst_type'] == 'shapes'
    assert _get_syst_attrs('JER',d['JER'])[1].at['shape_sigma'] == 1.0
    assert _get_syst_attrs('JER',d['JER'])[1].at['source_filename'] == 'path/FILE_DOWN'
    assert _get_syst_attrs('JER',d['JER'])[1].at['source_histname'] == 'HIST'
    assert _get_syst_attrs('JER',d['JER'])[1].at['direction'] == 'Down'

    with pytest.raises(RuntimeError):
        _get_syst_attrs('FAKE',d['FAKE'])

    with pytest.raises(RuntimeError):
        _get_syst_attrs('JER',d)

def test__df_condense_nameinfo():
    pass # tested in test__FullTable

def test__df_sanity_checks():
    c = Config('test/twoDtest.json','unit_test',
                findreplace={},
                externalOptions={}
                )
    healthy_df = c.FullTable()
    sick_df = healthy_df.append(healthy_df.iloc[0,:])
    _df_sanity_checks(healthy_df)
    with pytest.raises(RuntimeError):
        _df_sanity_checks(sick_df)

def test_config_loop_replace():
    config = open_json('test/testjson.json')
    assert "THIS" in config_loop_replace(config, "this", "THIS")
    assert config_loop_replace(config, "is", "IS")["THIS"] == "IS"
    assert "DICTIONARY" in config_loop_replace(config, "dictionary", "DICTIONARY")['a']
    with pytest.raises(TypeError):
        config_loop_replace("dummy",1,2)