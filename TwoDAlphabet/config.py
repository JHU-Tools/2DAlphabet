from collections import OrderedDict
import ROOT, json, os, pandas, re, warnings, itertools
from numpy import nan
import pprint
pp = pprint.PrettyPrinter(indent=4)
from TwoDAlphabet.helpers import copy_update_dict, open_json, parse_arg_dict, replace_multi
from TwoDAlphabet.binning import Binning, copy_hist_with_new_bins, get_bins_from_hist

_protected_keys = ["PROCESSES","SYSTEMATICS","REGIONS","BINNING","OPTIONS","GLOBAL","SCALE","COLOR","TYPE","X","Y","TITLE","BINS","NBINS","LOW","HIGH"]
_syst_col_defaults = {
    # 'variation': nan,
    'lnN': nan,
    'shapes': nan, # shape sigma
    'syst_type': nan,
    'source_filename': nan,
    'source_histname': nan,
    'direction': nan,
    'variation_alias': nan
}
class Config:
    '''Class to handle the reading and manipulation of data provided 
    in 2DAlphabet JSON configuration files. Constructor initializes
    a Config object for a given set of JSON files and performs
    all initial checks and manipulations.

    Args:
        jsonPath (str): File name and path.
        findreplace (dict, optional): Find-replace pairs. Defaults to {}.

    Attributes:
        config (dict): JSON config opened as a dict.
        nsignals (int): Number of signal processes. Zero before running Construct().
        nbkgs (int): Number of signal processes. Zero before running Construct().
    '''
    def __init__(self,jsonPath,findreplace={}):
        self.config = open_json(jsonPath)
        self._addFindReplace(findreplace)
        self.iterWorkspaceObjs = {}
        if 'GLOBAL' in self.config.keys(): self._varReplacement()
        self._addedConfigs = []

    def _section(self,key):
        '''Derive the dictionary for a given section of the configuration file
        with HELP keys removed.

        Args:
            key (str): Section name (all capital letters) of the configuration file.

        Returns:
            dict: Section of config.
        '''
        if not isinstance(self.config[key],dict):
            raise TypeError('Section access is only for sub-dictionaries. Try accessing directly reading config[key].')
        return {k:v for k,v in self.config[key].items() if k != 'HELP'}

    def _addFindReplace(self,findreplace):
        '''Add external find-replace pairs to the "GLOBAL" entry of config.

        Args:
            findreplace (dict): Find-replace pairs in non-nested dictionary.

        Raises:
            ValueError: If a "find" already exists in self.config['GLOBAL'].
        '''
        for s in findreplace:
            if s in self._section('GLOBAL'):
                raise ValueError('A command line string replacement (%s) conflicts with one already in the configuration file.' %(s))
            self.config['GLOBAL'][s] = findreplace[s]

    def _varReplacement(self):
        '''Do string substitution for config entries based on the dictionary of find
        and replace values found in self.config['GLOBAL']. Call self._addFindReplace()
        before running this function to add in external find-replace pairs.

        Args:
            findreplace (dict): Non-nested dictionary with key-value pairs to find-replace
                in the internal class configuration dict.

        Returns:
            None.
        '''
        print ("Doing GLOBAL variable replacement in input json...")
        for old_string in self._section('GLOBAL'):
            if old_string == "HELP":
                print ('WARNING: The HELP entry is deprecated and checking for it will be removed in the future. Please remove it from your config.')
                continue
            new_obj = self._section('GLOBAL')[old_string]
            if isinstance(new_obj,list):
                self.iterWorkspaceObjs[old_string] = new_obj
            else:
                self.config = config_loop_replace(self.config, old_string, new_obj)

    def SaveOut(self, projPath): # pragma: no cover
        '''Save the histogram table to the `projPath` directory in csv
        and markdown formats. No copy of the input JSONs is saved since multiple
        could be provided with the only final combination making it to the histogram table.
        '''
        file_out = open(projPath+'runConfig.json', 'w')
        json.dump(self.config,file_out,indent=2,sort_keys=True)
        file_out.close()

    def FullTable(self):
        '''Generate full table of processes, regions, and systematic variations
        to account for, including relevant information for each. The table is
        returned as a pandas DataFrame for convenient manipulation.

        Returns:
            pandas.DataFrame: Table
        '''
        regions = self._regionTable()
        processes = self._processTable()
        systematics = self._systematicsTable()

        for p,group in processes.groupby(processes.index):
            if group.title.nunique() > 1:
                raise RuntimeError('There are multiple titles specified for process %s. Do not use multiple substitution for titles. Use plotting options instead.\n%s'%(p,group))

        for syst in list(processes.variation.unique()):
            if syst not in list(systematics.index.unique())+['nominal']:
                raise NameError('Systematic variation %s does not exist among possible systematics:\n\t%s\nIs there a typo in the JSON?'
                                %(syst,list(systematics.index.unique())))

        proc_syst = processes.merge(systematics,right_index=True,left_on='variation',how='left',suffixes=['','_syst'])
        proc_syst = _df_condense_nameinfo(proc_syst,'source_histname')
        proc_syst = _df_condense_nameinfo(proc_syst,'source_filename')

        final = regions.merge(proc_syst,right_index=True,left_on='process',how='left')
        final = _keyword_replace(final, ['source_filename', 'source_histname']).reset_index(drop=True)
        _df_sanity_checks(final)
        return final

    def _regionTable(self):
        '''Generate the table of region information based on the JSON config.
        Columns are `process`, `region`, and `binning`.

        Returns:
            pandas.DataFrame
        '''
        def _data_not_included(region):
            '''Check if the list of processes associated with a region
            includes the one marked as `type` `DATA` in the `PROCESSES`
            section of the config. If it doesn't included it but it exists,
            return the name so it can be added to the list of processes
            for the region.

            Args:
                region (str): Name of region to check.

            Raises:
                RuntimeError: No process of TYPE == "DATA" provided.
                RuntimeError: Multiple processes of TYPE == "DATA" provided in PROCESSES section.
                RuntimeError: Multiple processes of TYPE == "DATA" provided in REGIONS subsection.

            Returns:
                str: The name of the data key if it's not already included in the list of processes
                for the region.
            '''
            region_df = pandas.DataFrame({'process':self._section('REGIONS')[region]['PROCESSES']})
            process_df = pandas.DataFrame(self._section('PROCESSES')).T[['TYPE']]
            region_df = region_df.merge(process_df,
                                        left_on='process',
                                        right_index=True,
                                        how='left')

            # Check DATA type is even provided in the PROCESSES
            if (process_df['TYPE'] == 'DATA').sum() == 0:
                warnings.warn('No process of TYPE == "DATA" provided. Ignoring...', RuntimeWarning)
                data_key = False
            elif (process_df['TYPE'] == 'DATA').sum() > 1:
                raise RuntimeError('Multiple processes of TYPE == "DATA" provided in PROCESSES section.')
            else:
                data_key = 'data_obs'
            # Check if it's included in the regions
            if (region_df['TYPE'] == 'DATA').sum() == 0:
                out = data_key
            elif (region_df['TYPE'] == 'DATA').sum() == 1:
                out = False
            else:
                raise RuntimeError('Multiple processes of TYPE == "DATA" provided in REGIONS subsection.')

            return out

        out_df = pandas.DataFrame(columns=['process','region','binning'])
        for r in self._section('REGIONS'):
            data_key = _data_not_included(r)
            if data_key:
                out_df = out_df.append(pandas.Series({'process':data_key,'region':r, 'binning':self._section('REGIONS')[r]['BINNING']}),ignore_index=True)

            for p in self._section('REGIONS')[r]['PROCESSES']:
                if p not in self._section('PROCESSES') and len([kglobal for kglobal in self._section('GLOBAL') if kglobal in p]) == 0:
                    raise RuntimeError('Process "%s" listed for region "%s" not defined in PROCESSES section.'%(p,r))
                
                row_format = lambda c: {
                    'process': c['PROCESS'],
                    'region': c['REGION'],
                    'binning':self._section('REGIONS')[c['REGION']]['BINNING']
                }
                rows_to_append = self._iterObjReplaceProducer({'PROCESS':p, 'REGION':r}, row_format)
                for new_row in rows_to_append:
                    out_df = out_df.append(new_row,ignore_index=True)
                
        return out_df

    def _processTable(self):
        '''Generate the table of process information based on the JSON config.
        Columns are `color`, `process_type`, `scale`, `variation`, `source_filename`,
        and `source_histname`.

        Returns:
            pandas.DataFrame
        '''
        out_df = pandas.DataFrame(columns=['color','process_type','scale','variation','source_filename','source_histname','alias','title','combine_idx'])
        for p in self._section('PROCESSES'):
            this_proc_info = self._section('PROCESSES')[p]
            this_proc_info['NAME'] = p
            if this_proc_info['TYPE'] == 'DATA' and p != 'data_obs':
                raise RuntimeError('Any process of type DATA must have section key "data_obs".')
            for s in this_proc_info['SYSTEMATICS']+['nominal']:
                this_proc_info['VARIATION'] = s
                row_format = lambda info: pandas.Series(
                    {'color': nan if 'COLOR' not in info else info['COLOR'],
                    'process_type': info['TYPE'],
                    'scale': 1.0 if 'SCALE' not in info else info['SCALE'],
                    'source_filename': info['LOC'].split(':')[0],
                    'source_histname': info['LOC'].split(':')[1],
                    'alias': info['NAME'] if 'ALIAS' not in info.keys() else info['ALIAS'], #in file name
                    'title': info['NAME'] if 'TITLE' not in info.keys() else info['TITLE'], #in legend entry
                    'variation': info['VARIATION'],
                    }, name=info['NAME']
                )
                rows_to_append = self._iterObjReplaceProducer(this_proc_info, row_format)
                for new_row in rows_to_append:
                    out_df = out_df.append(new_row)

        return out_df

    def _systematicsTable(self):
        '''Generate the table of process information based on the JSON config.
        Columns are  'lnN', 'shapes', 'syst_type', 'source_filename',
        'source_histname', and 'direction' (ie. NaN, 'Up', or 'Down').

        Note that 'shapes' is short for 'shape sigma'.

        Returns:
            pandas.DataFrame
        '''
        out_df = pandas.DataFrame(columns=_syst_col_defaults.keys())
        for s in self._section('SYSTEMATICS'):
            iterations_to_process = self._iterObjReplaceProducer(self._section('SYSTEMATICS')[s], lambda c: c)
            for iteration in iterations_to_process:
                for syst in _get_syst_attrs(s, iteration):
                    out_df = out_df.append(syst)
        
        return out_df

    def _iterObjReplaceProducer(self, obj_package, func):
        '''Pre-processes input to DataFrame in the case that the inputs
        can take multiple values with keyword replacement.

        Args:
            obj_package ([type]): [description]
            func ([type]): [description]
        '''
        def _iterativeReplace(base,find_replace_map):
            out = base
            for f,r in find_replace_map.items():
                out = out.replace(f,r)
            return out

        to_vary = OrderedDict()
        for objKey, obj in obj_package.items(): # do replacement for multiple objects
            if not isinstance(obj, str):
                to_vary[objKey] = [obj]
                continue
            to_vary[objKey] = []

            matches = [] # find all matches to this obj
            for iterKey in self.iterWorkspaceObjs.keys():
                if iterKey in obj:
                    matches.append(iterKey)

            # Build all combinations of N matches
            replacements_for_matches = [self.iterWorkspaceObjs[iterKey] for iterKey in matches]
            for replacement_set in itertools.product(*replacements_for_matches):
                if isinstance(replacement_set, str): replacement_set = set(replacement_set)
                # Loop of replacements, perform replacements, and track the variation
                find_replace_map = {matches[i]:replacement_set[i] for i in range(len(matches))}
                to_vary[objKey].append(_iterativeReplace(obj, find_replace_map))

            if len(to_vary[objKey]) == 0:
                # if no replacements, store original
                to_vary[objKey] = [obj]

        # Use func to plug everything back together
        out = []
        for varied_set in itertools.product(*(to_vary.values())):
            varied_set_map = {key: varied_set[i] for i,key in enumerate(to_vary.keys())}
            out.append(func(varied_set_map))
        
        return out

    def Add(self,cNew,onlyOn=['process','region']):
        raise NotImplementedError('Multiple config support is currently a work in progress. Only the first config will be used.')
        def _drop(row,dupes_list):
            drop = False
            for d in dupes_list:
                to_compare = []
                for i,c in enumerate(onlyOn): # build bool from onlyOn cols
                    to_compare.append(row[c] == d[i])
                drop = pandas.Series(to_compare).all()

                if drop == True: # if we found a match to a
                    break
            return drop

        if self.constructed == True:
            raise RuntimeError('This config has already been constructed so no additions can be made.')
        if isinstance(onlyOn,str):
            if onlyOn not in ['process','region']:
                raise RuntimeError('Can only add configs together on the "process" or "region" information.')
            onlyOn = [onlyOn]
        elif onlyOn != ['process','region']:
            raise RuntimeError('Can only add configs together on the "process" or "region" information.')
        
        df_modified_base         = self.df.append(cNew.df).reset_index(drop=True)
        df_modified_nominal_only = df_modified_base[df_modified_base.variation.eq('nominal')]
        df_modified_dupes        = df_modified_nominal_only[ df_modified_nominal_only.duplicated(subset=onlyOn,keep='first') ]

        dupes_list = set(zip(*(df_modified_dupes[k] for k in onlyOn)))
        if len(dupes_list) > 0: # if duplicates, replace old with new
            print ('Found duplicates in attempting to modify base Config. Replacing...')
            for d in dupes_list:
                print('\t(%s)'%(','.join(d)))
            df_final = self.df.loc[
                            ~self.df.apply(_drop,args=[dupes_list],axis='columns')
                        ].append(cNew.df).reset_index(drop=True)

        else: # if no duplicates, just use the appended df
            df_final = df_modified_base

        self._addedConfigs.append(cNew)
        self.df = df_final

    def GetNregions(self):
        return self.df.region.nunique()

    def GetNsystematics(self):
        return self.df.variation.nunique()-1

class OrganizedHists():
    '''Class to store histograms in a consistent data structure and with accompanying
    methods to access the histograms.

    Attributes:
        name (str): Name, taken from input configObj.
        filename (str): Path to `organized_hists.root`.
        hists (dict): Three-level nested dictionary organized as [process][region][systematic variation].
        binning (Binning): Binning object, taken from configObj.
        rebinned (bool): Flag to denote if a rebinning has already occured.
        file (ROOT.TFile): TFile to store histograms on disk.

    Args:
        configObj (Config): Config object.
    '''
    def __init__(self,projPath,binnings,hist_map,readOnly=False):
        self.filename = projPath + 'organized_hists.root'
        self.hist_map = hist_map

        if os.path.exists(self.filename) and readOnly:
            self.file = ROOT.TFile.Open(self.filename,"OPEN")
        else:
            self.file = ROOT.TFile.Open(self.filename,"RECREATE")
            self.Add(binnings)
            self.file.Close()
            self.file = ROOT.TFile.Open(self.filename,"OPEN")

    def Add(self, binnings):
        '''Manipulate all histograms in self.hist_map and save them to organized_hists.root.

        Returns:
            None
        '''
        for infilename,histdf in self.hist_map.items():
            infile = ROOT.TFile.Open(infilename)
            for row in histdf.itertuples():
                if row.source_histname not in [k.GetName() for k in infile.GetListOfKeys()]:
                    raise NameError('Histogram name %s does not exist in file %s.'%(row.source_histname,infile.GetName()))
                h = infile.Get(row.source_histname)
                h.SetDirectory(0)
                h.Scale(row.scale)
                binning = binnings[row.binning]

                if get_bins_from_hist("Y", h) != binning.ybinList:
                    h = copy_hist_with_new_bins(row.out_histname+'_rebinY','Y',h,binning.ybinList)
                if get_bins_from_hist("X", h) != binning.xbinList:
                    h = copy_hist_with_new_bins(row.out_histname,'X',h,binning.xbinList)
                else:
                    h.SetName(row.out_histname)

                h.SetTitle(row.out_histname)
                h.SetFillColor(row.color)

                self.file.WriteTObject(h, row.out_histname)
                self.CreateSubRegions(h, binning)

            infile.Close()

    def Get(self,histname='',process='',region='',systematic='',subspace='FULL'):
        '''Get histogram from the opened TFile. Specify the histogram
        you want via `histname` or by the combination of `process`, `region`,
        and `systematic` options. The `histname` option will take priority.

        Args:
            histname (str, optional): Name of histogram to get. Overrides other three options if specified. Defaults to ''.
            process (str, optional): Name of process to search for. Must be used in conjunction with `region` and `systematic` options. Overridden by `histname`. Defaults to ''.
            region (str, optional): Name of region to search for. Must be used in conjunction with `process` and `systematic` options. Overridden by `histname`. Defaults to ''.
            systematic (str, optional): Name of systematic to search for. Must be used in conjunction with `process` and `region` options. Overridden by `histname`. Defaults to ''.
            subspace (str, optional): Name of subspace. Default is 'FULL' with other options being 'LOW', 'SIG', and 'HIGH'.

        Raises:
            NameError: If subspace option is not 'FULL','LOW','SIG', or 'HIGH'.

        Returns:
            TH2F: Histogram from file.
        '''
        if subspace not in ['FULL','LOW','SIG','HIGH']:
            raise NameError("Subspace '%s' not accepted. Options are 'FULL','LOW','SIG','HIGH'.")
        if histname == '':
            histname = '_'.join([process,region,subspace])
            if systematic != '':
                histname+='_'+systematic

        all_histnames = [k.GetName() for k in self.file.GetListOfKeys()]
        if histname not in all_histnames:
            raise NameError('Histogram %s does not exist.'%(histname))

        return self.file.Get(histname)

    def GetHistNames(self):
        return [hkey.GetName() for hkey in self.file.GetListOfKeys()]

    def BinningLookup(self,histname):
        all_hists = pandas.concat([v[['out_histname','binning']] for v in self.hist_map.values()])
        return all_hists.loc[all_hists.out_histname.eq(histname)].iloc[0].binning

    def CreateSubRegions(self,h,binning):
        '''Sub-divide input histogram along the X axis into the regions specified in the config
        and write the new histogram to organized_hists.root.

        Returns:
            None
        '''
        for sub in binning.xbinByCat.keys():
            hsub = h.Clone()
            hsub = copy_hist_with_new_bins(h.GetName().replace('_FULL','_'+sub),'X',h,binning.xbinByCat[sub])
            hsub.SetTitle(hsub.GetName())
            if hsub.Integral() <= 0:
                print ('WARNING: %s has zero or negative events - %s'%(hsub.GetName(), hsub.Integral()))
                for b in range(1,hsub.GetNbinsX()*hsub.GetNbinsY()+1):
                    hsub.SetBinContent(b,1e-6)
            self.file.WriteObject(hsub, hsub.GetName())

def _keyword_replace(df,col_strs):
    '''Given a DataFrame and list of column names,
    find and replace the three keywords ("$process", "$region$", "$syst") with their
    respective values in the row for the DataFrame.

    Args:
        df (pandas.DataFrame): DataFrame in which to do the find-replace and to find the per-row keyword matches.
        col_strs (list(str)): List of column names to consider for the find-replace.

    Returns:
        pandas.DataFrame: The manipulated DataFrame copy.
    '''
    def _batch_replace(row,s=None):
        if pandas.isna(row[s]):
            return nan
        else:
            return replace_multi(
                row[col_str],
                {'$process': row.alias,
                 '$region':  row.region,
                 '$syst':    row.variation_alias}
            )

    for col_str in col_strs:
        df[col_str] = df.apply(_batch_replace, axis='columns', s=col_str)
    return df

def _get_syst_attrs(name,syst_dict):
    '''Parse an entry in the `"SYSTEMATICS"` section of the JSON config and
    generate the row(s) to append to the systematics DataFrame based on
    that information.

    Args:
        name (str): Name of the systematic variation.
        syst_dict (dict): Dictionary of the config["SYSTEMATICS"][name] section of the JSON config.

    Raises:
        RuntimeError: Systematic variation type could not be determined.

    Returns:
        list(pands.Series): List of new rows to append to the main systematics DataFrame.
    '''
    if 'VAL' in syst_dict:
        out = [{
            'lnN':str(syst_dict['VAL']),
            'syst_type': 'lnN'
        }]
    elif 'VALUP' in syst_dict and 'VALDOWN' in syst_dict:
        out = [{
            'lnN':'%s/%s'%(syst_dict['VALDOWN'], syst_dict['VALUP']),
            'syst_type': 'lnN'
        }]
    elif 'UP' in syst_dict and 'DOWN' in syst_dict:
        out = [
            {
                'shapes':syst_dict['SIGMA'],
                'syst_type': 'shapes',
                'source_filename': syst_dict['UP'].split(':')[0],
                'source_histname': syst_dict['UP'].split(':')[1],
                'direction': 'Up',
                'variation_alias': name if 'ALIAS' not in syst_dict else syst_dict['ALIAS']
            }, {
                'shapes':syst_dict['SIGMA'],
                'syst_type': 'shapes',
                'source_filename': syst_dict['DOWN'].split(':')[0],
                'source_histname': syst_dict['DOWN'].split(':')[1],
                'direction': 'Down',
                'variation_alias': name if 'ALIAS' not in syst_dict else syst_dict['ALIAS']
            }
        ]
    else:
        raise RuntimeError('Systematic variation type could not be determined for "%s".'%name)

    out = [pandas.Series(copy_update_dict(_syst_col_defaults, d), name=name) for d in out]
    return out

def _df_condense_nameinfo(df,baseColName):
    '''Condense information after the left-join of the process and systematics DataFrames which creates duplicates
    of the `source_*` columns. Manipulate `df` and return it so that `baseColName+"_syst"` replaces `baseColName`
    and is then dropped.

    Args:
        df (pandas.DataFrame): Input DataFrame to condense.
        baseColName (str): Name of column to condense into.

    Returns:
        pandas.DataFrame: Condensed DataFrame.
    '''
    df[baseColName] = df.apply(lambda row: row[baseColName+'_syst'] if pandas.notna(row[baseColName+'_syst']) else row[baseColName], axis='columns')
    df.drop(baseColName+'_syst',axis='columns',inplace=True)
    return df

def _df_sanity_checks(df):
    '''Check for duplicated (process,region,variation,source_filename,source_histname).
    Prints any duplicate rows if they exist (and raises RuntimeError).

    Args:
        df (pandas.DataFrame): DataFrame to check.

    Raises:
        RuntimeError: Duplicates exist (duplicated rows are printed).
    '''
    # check for duplicate process+region+variation
    dupes = df[df.duplicated(subset=['process','region','variation','source_filename','source_histname'],keep=False)]
    if dupes.shape[0] > 0:
        raise RuntimeError('Duplicates exist. Printing them...\n%s'%dupes)
      
def config_loop_replace(config,old,new,inGLOBAL=False):
    '''Self-calling function loop to find-replace one pair (old,new)
    in a nested dictionary or list (config). If old, new, and the config entry
    examined are all strings, all instances of <old> will be replaced by <new> in the
    string in the config. If this condition is not satisified, <old> must match the config
    entry in its entirety (ie. <old> == <config dict value>).

    Args:
        config (dict,list): Nested dictionary or list where keys and values will have the
            find-replace algorithm applied.
        old (non-nested obj): Object to replace (of type string, int, float, etc - no lists or dictionaries).
        new (non-nested obj): Object replacement (of type string, int, float, etc) - no lists or dictionaries).

    Raises:
        TypeError: If input is not a dict or list.

    Returns:
        dict,list: Modified dict/list.
    '''
    next_is_global = False
    if isinstance(config,dict):
        for k,v in config.items():
            if k == 'GLOBAL':
                next_is_global = True
            if old in k and not inGLOBAL:
                config[re.sub(r'\b%s\b'%old,new,k)] = config.pop(k)
                k = re.sub(r'\b%s\b'%old,new,k)
            if isinstance(v,dict) or isinstance(v,list):
                config[k] = config_loop_replace(v,old,new,inGLOBAL=next_is_global)
            elif isinstance(v,str) and isinstance(old,str) and isinstance(new,str):
                if old in v:
                    config[k] = re.sub(r'\b%s\b'%old,new,v)
            else:
                if old == v:
                    config[k] = new
            next_is_global = False # never consider next k to be GLOBAL
    elif isinstance(config,list):
        for i,v in enumerate(config):
            if isinstance(v,dict) or isinstance(v,list):
                config[i] = config_loop_replace(v,old,new)
            elif isinstance(v,str):
                if old in v:
                    config[i] = re.sub(r'\b%s\b'%old,new,v)
            else:
                if old == v:
                    config[i] = new
    else:
        raise TypeError('Type "%s" not accepted in config_loop_replace.')

    return config
