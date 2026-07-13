from collections import OrderedDict
import ROOT, json, os, pandas, re, warnings, itertools
import math
from numpy import nan
import pprint
pp = pprint.PrettyPrinter(indent=4)
from TwoDAlphabet.plotstyle import mpl_to_root_colors, root_to_matplotlib_color
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

    def FullTable(self,verbose=False):
        '''Generate full table of processes, regions, and systematic variations
        to account for, including relevant information for each. The table is
        returned as a pandas DataFrame for convenient manipulation.

        Args:
            verbose (bool): If True, prints the regions, processes, and systematics dataframe tables to .txt files.

        Returns:
            pandas.DataFrame: Table
        '''
        regions = self._regionTable()
        processes = self._processTable()
        systematics = self._systematicsTable()

        if verbose:
            regions.to_string('regions.txt')
            processes.to_string('processes.txt')
            systematics.to_string('systematics.txt')

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
                new_row = pandas.DataFrame([{
                    'process': data_key,
                    'region': r,
                    'binning': self._section('REGIONS')[r]['BINNING']
                }])
                out_df = pandas.concat([out_df, new_row], ignore_index=True)


            for p in self._section('REGIONS')[r]['PROCESSES']:
                if p not in self._section('PROCESSES') and len([kglobal for kglobal in self._section('GLOBAL') if kglobal in p]) == 0:
                    raise RuntimeError('Process "%s" listed for region "%s" not defined in PROCESSES section.'%(p,r))

                row_format = lambda c: pandas.DataFrame([{
                    'process': c['PROCESS'],
                    'region': c['REGION'],
                    'binning': self._section('REGIONS')[c['REGION']]['BINNING']
                }])
                rows_to_append = self._iterObjReplaceProducer({'PROCESS':p, 'REGION':r}, row_format)
                for new_row in rows_to_append:
                    out_df = pandas.concat([out_df,new_row],ignore_index=True)
  
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
                    'source_filename': info['LOC'].split(':')[0] if 'root://' not in info['LOC'] else info['LOC'].split(':')[0]+':'+info['LOC'].split(':')[1],
                    'source_histname': info['LOC'].split(':')[1] if 'root://' not in info['LOC'] else info['LOC'].split(':')[2],
                    'alias': info['NAME'] if 'ALIAS' not in info.keys() else info['ALIAS'], #in file name
                    'title': info['NAME'] if 'TITLE' not in info.keys() else info['TITLE'], #in legend entry
                    'variation': info['VARIATION'],
                    }, name=info['NAME']
                )
                rows_to_append = self._iterObjReplaceProducer(this_proc_info, row_format)
                for new_row in rows_to_append:
                    new_row = new_row.to_frame().T
                    out_df = pandas.concat([out_df, new_row], ignore_index=False)

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
                    syst = syst.to_frame().T
                    out_df = pandas.concat([out_df,syst])
        
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
            infilename = infilename[0] #Gets extracted as tuple for some reason
            infile = ROOT.TFile.Open(infilename)
            for row in histdf.itertuples():
                if row.source_histname not in [k.GetName() for k in infile.GetListOfKeys()]:
                    raise NameError('Histogram name %s does not exist in file %s.'%(row.source_histname,infile.GetName()))
                h = infile.Get(row.source_histname)
                h.SetDirectory(0)
                h.Scale(row.scale)
                binning = binnings[row.binning]

                if binning.is3D and get_bins_from_hist("Z", h) != binning.zbinList:
                    h = copy_hist_with_new_bins(row.out_histname+'_rebinZ','Z',h,binning.zbinList)
                if get_bins_from_hist("Y", h) != binning.ybinList:
                    h = copy_hist_with_new_bins(row.out_histname+'_rebinY','Y',h,binning.ybinList)
                if get_bins_from_hist("X", h) != binning.xbinList:
                    h = copy_hist_with_new_bins(row.out_histname,'X',h,binning.xbinList)
                else:
                    h.SetName(row.out_histname)

                h.SetTitle(row.out_histname)
                if row.color not in mpl_to_root_colors.keys():
                    available_colors = '", "'.join(mpl_to_root_colors.keys())
                    raise ValueError(f'Color "{row.color}" not defined. Please add the ROOT TColor code to the "mpl_to_root_colors" dictionary defined in TwoDAlphabet.plotstyle. Available default colors are: "{available_colors}"')
                else:
                    h.SetFillColor(mpl_to_root_colors[row.color])

                self.file.WriteTObject(h, row.out_histname)
                self.CreateSubRegions(h, binning)

            infile.Close()

    def AddMCStatShapes(self, df, binnings, threshold=10, include_signal=False, excluded_procs=[], alpha_min=0.1, name_prefix='mcstat', verbose=True):
        '''
        Generate autoMCStats-style per-bin shape templates from the rebinned nominal
        background hists, write them to organized_hists.root, register
        them in hist_map for BinningLookup(), and return rows to append to the ledger.

        Call after self.Add(binnings). Returns list[dict] (syst_type='shapes').
        Set verbose=True for a per-bin report in the style of Combine's autoMCStats.
        '''
        self.file.Close()
        self.file = ROOT.TFile.Open(self.filename, 'UPDATE')

        types = ['BKG', 'SIGNAL'] if include_signal else ['BKG'] # normally we don't want to add MC stat unc on signal, but keep it just in case
        bkg = df[df.process_type.isin(types)]
        if bkg.empty: # Warn user that MCstats will not be run at all, return
            print('\n'+'='*100)
            print('WARNING: All MC background processes excluded from calculation of MC statistical uncertainties.')
            print(f'\tExcluded processes: {", ".join(excluded_procs)}')
            print('='*100+'\n')
            return []
        existing = self.GetHistNames()
        new_rows = []
        hist_map_rows = []
        n_bblite = n_perproc = n_skipped = 0

        # Loop over all regions in the workspace and get the nominal histograms for all unique processes.
        # Store them in a {process: TH2} dict called `nominal`.
        # Avoid any processes specified by the user in the JSON, passed here as `excluded_procs` list.
        for region in df.region.unique():
            procs = list(bkg[bkg.region.eq(region) & ~bkg.process.isin(excluded_procs)].process.unique())
            if not procs:
                # This occurs when `include_signal=False` and all MC processes are excluded. Just return. 
                print('\n'+'='*100)
                print('WARNING: All MC processes excluded from calculation of MC statistical uncertainties.')
                print(f'\tExcluded processes: {", ".join(excluded_procs + list(df[df.process_type.isin(["SIGNAL"])].process.unique()))}')
                print('='*100+'\n')
                return []
            binning_name = df[df.region.eq(region)].binning.iloc[0]
            binning = binnings[binning_name]

            nominal = {}
            for p in procs:
                hname = '%s_%s_FULL' % (p, region)
                if hname in existing:
                    h = self.file.Get(hname); h.SetDirectory(0)
                    nominal[p] = h
            if not nominal:
                continue

            # For verbose output in the style of Combine: load the excluded non-data nominal histos, for the "total sum" line
            all_nominal = dict(nominal)
            for p in df[df.region.eq(region) & df.process_type.ne('DATA')].process.unique():
                if p in all_nominal:
                    continue
                hname = '%s_%s_FULL' % (p, region)
                if hname in existing:
                    h = self.file.Get(hname); h.SetDirectory(0)
                    all_nominal[p] = h
            excluded = [p for p in all_nominal if p not in nominal]
            meta_data = {p: df[df.process.eq(p) & df.region.eq(region) & df.variation.eq('nominal')].iloc[0] for p in nominal}
            sample = next(iter(nominal.values()))
            nx, ny = sample.GetNbinsX(), sample.GetNbinsY()

            if verbose:
                print('\n' + '=' * 64)
                print('Analyzing bin errors for region: %s' % region)
                print('Poisson cut-off: %d' % threshold)
                print('Processes summed: %s' % ' '.join(nominal.keys()))
                print('Processes excluded from sums: %s' % (' '.join(excluded) if excluded else '(none)'))
                print('=' * 64)
                print('%-12s  %14s  %14s  %s' % ('Bin', 'Contents', 'Error', 'Notes'))

            def _make_mcstat_template(proc, variation, bx, by, rel):
                '''
                Produce a 2D histogram for a given process, where only a single bin (bx, by) is different from the nominal.
                This bin is shifted up and down by some multiplicative factor `rel`, and the down variation is bounded below at zero. 
                The factor `rel` is determined based on the number of effective events. 
                '''
                c = nominal[proc].GetBinContent(bx, by) # Nominal bin contents
                for direction, factor in (('Up', 1.0 + rel), ('Down', max(0.0, 1.0 - rel))):
                    h = nominal[proc].Clone('%s_%s_FULL_%s%s' % (proc, region, variation, direction))
                    h.SetDirectory(0)
                    h.SetBinContent(bx, by, c * factor)
                    h.SetTitle(h.GetName())
                    self.file.WriteTObject(h, h.GetName())
                    self.CreateSubRegions(h, binning)
                    # register the "FULL" template so BinningLookup() can resolve it later. This way, 2DA does the RooDataHist creation for us
                    hist_map_rows.append({
                        'source_histname': h.GetName(),
                        'out_histname':    h.GetName(),
                        'scale':           1.0,
                        'color':           meta_data[proc].get('color', 0),
                        'binning':         binning_name,
                    })

                    row = meta_data[proc].copy()
                    row['variation'] = variation
                    row['variation_alias'] = variation
                    row['direction'] = direction
                    row['syst_type'] = 'shapes'
                    row['shapes'] = 1.0
                    row['lnN'] = nan
                    row['mcstat'] = True
                    new_rows.append(row.to_dict())

            # Loop over all bins, and all nominal MC processes in this region. Perform the BB/BB-lite algorithm
            for bx in range(1, nx + 1):
                for by in range(1, ny + 1):
                    tag = f'{region} ({bx}, {by})' #'bx%d_by%d' % (bx, by)
                    n_tot, e2 = 0.0, 0.0
                    for p in nominal:
                        n_tot += nominal[p].GetBinContent(bx, by)
                        e2    += nominal[p].GetBinError(bx, by) ** 2
                    e_tot = math.sqrt(e2)
                    if verbose:
                        na, e2a = 0.0, 0.0
                        for p in all_nominal:
                            na  += all_nominal[p].GetBinContent(bx, by)
                            e2a += all_nominal[p].GetBinError(bx, by) ** 2
                        print('-' * 64)
                        print('%-12s  %14.6f  %14.6f  %s (%s)' % (tag, na, math.sqrt(e2a), 'total sum', ', '.join([p for p in all_nominal])))
                        print('%-12s  %14.6f  %14.6f  %s (%s)' % (tag, n_tot, e_tot, 'excluding marked processes', ', '.join(excluded)))

                    if e_tot == 0.0:
                        n_skipped += 1
                        if verbose:
                            print('  => Error is zero, ignore')
                        continue

                    N_eff_tot = int(round(n_tot * n_tot / (e_tot * e_tot)))

                    alpha = n_tot / N_eff_tot if N_eff_tot > 0 else 0.0
                    if alpha < alpha_min:
                        print('%-12s  %14.6f  %14.6f  %s' % (
                              tag, float(N_eff_tot), math.sqrt(N_eff_tot),
                              'Effective events, alpha=%.6f' % alpha))
                        print('  => alpha < alpha_min: %.3f < %.3f'%(alpha, alpha_min))
                        print('  => No MC statistical uncertainty parameters will be calculated')
                        continue
                    if verbose:
                        print('%-12s  %14.6f  %14.6f  %s' % (
                            tag, float(N_eff_tot), math.sqrt(N_eff_tot),
                            'Effective events, alpha=%.6f' % alpha))

                    ##################################################################################################
                    # Case 1: n_tot^eff > threshold
                    # In this case, we make a single Gaussian-constrained NP that scales the total yield in the bin.
                    # Since we're doing this with shape templates instead of whatever Combine does, we have to ensure
                    # that every process in this bin shares the same nuisance parameter name (histogram name), and that 
                    # every process has the same scaling value `rel`. In Case 1, we use `rel = e_tot / n_tot`
                    #
                    # This means that, in a given bin, process `i` contributes n_i * (1 + v*rel) for nuisance value `v`.
                    # For all processes in this bin, we get:
                    #   \Sum_i n_i * (1 + v * rel) = (1 + v * rel) * \Sum_i n_i
                    #                              = (1 + v * e_tot / n_tot) * n_tot
                    #                              = n_tot + v * e_tot
                    # 
                    # Thus, we're left with +/- e_tot at for v = +/-1. 
                    # Since e_tot = sqrt(\Sum_i e_i^2), we get the combined MC stat uncertainty in this given bin. 
                    ##################################################################################################
                    if N_eff_tot > threshold:
                        variation = '%s_%s_bx%d_by%d' % (name_prefix, region, bx, by)
                        rel = e_tot / n_tot
                        if verbose:
                            print('  => N_eff > %d : shared gaussian shape %r (rel=%.6f) shared by [%s]'
                                % (threshold, variation, rel, ', '.join(nominal.keys())))
                        for p in nominal:
                            _make_mcstat_template(p, variation, bx, by, rel)
                        n_bblite += 1

                    ##################################################################################################
                    # Case 2: n_tot^eff < threshold
                    # In this case, the number of effective events is less than the threshold. Now, we simply make a 
                    # new TH2 for each process, where every bin is kept the same except for (bx, by) which is shifted
                    # up and down by +/-1 sigma. 
                    # Here, we use `rel == width = e_i / n_i`, where "i" represents the process. Thus, in a given bin,
                    # a process "i" has a yield (governed by the value `v` of the nuisance parameter):
                    #       yield_i(v) = n_i + v * e_i 
                    #                  = n_i·(1 + v * width)
                    #
                    # Since we are providing Combine the +/-1 sigma variation histograms, it will automatically 
                    # assign this nuisance parameter a unit-Gaussian prior. 
                    # NOTE: this is not strictly what Combine does, as Combine will implement a Poisson-constrained 
                    # parameter. However, there's no good/easy way to do this with the TH2/RooDataHist implementation 
                    # of 2DAlphabet, and in general a Gaussian is going to be good enough anyway. 
                    ##################################################################################################
                    else:
                        if verbose:
                            print('  => N_eff <= %d : per-process (Poisson approximated as gaussian)' % threshold)
                        made_one = False # the below-threshold bin where every process has zero content/error doesn't get counted as a per-process bin
                        for p in nominal:
                            c = nominal[p].GetBinContent(bx, by)
                            e = nominal[p].GetBinError(bx, by)
                            if verbose:
                                print('  ' + '-' * 58)
                                print('    %-18s %14.6f %14.6f' % (p, c, e))
                            if c <= 0.0 or e <= 0.0:
                                if verbose:
                                    print('      => Content or error is zero, ignore')
                                continue
                            N_eff_i = int(round(c * c / (e * e)))
                            width = e / c
                            variation = '%s_%s_%s_bx%d_by%d' % (name_prefix, region, p, bx, by)

                            alpha_i = c / N_eff_i if N_eff_i > 0 else 0.0
                            if alpha < alpha_min:
                                print('%-12s  %14.6f  %14.6f  %s' % (
                                    tag, float(N_eff_tot), math.sqrt(N_eff_tot),
                                    'Effective events, alpha=%.6f' % alpha))
                                print('  => alpha < alpha_min: %.3f < %.3f'%(alpha, alpha_min))
                                print('  => No MC statistical uncertainty parameters will be calculated')
                                continue
                            if verbose:
                                print('    %-18s %14.6f %14.6f  %s' % (
                                    '', float(N_eff_i), math.sqrt(N_eff_i),
                                    'Effective events, alpha=%.6f' % alpha_i))
                                print('      => %r [1.00,%.2f,%.2f] to be gaussian constrained (width=%.6f)'
                                    % (variation, max(0.0, 1.0 - 7.0 * width), 1.0 + 7.0 * width, width))
                            _make_mcstat_template(p, variation, bx, by, width)
                            made_one = True
                        if made_one:
                            n_perproc += 1

        if hist_map_rows:
            self.hist_map['__mcstat__'] = pandas.concat(
                [self.hist_map.get('__mcstat__'), pandas.DataFrame(hist_map_rows)],
                ignore_index=True) if '__mcstat__' in self.hist_map else pandas.DataFrame(hist_map_rows)
            hist_map_rows = []

        if verbose:
            n_nuis = len({r['variation'] for r in new_rows})
            print('\n' + '=' * 64)
            print('MC-stat shape summary')
            print('  BB-lite bins        : %d' % n_bblite)
            print('  per-process bins    : %d' % n_perproc)
            print('  skipped (zero error): %d' % n_skipped)
            print('  distinct nuisances  : %d  (%d shape templates)' % (n_nuis, len(new_rows)))
            print('=' * 64 + '\n')

        self.file.Close()
        self.file = ROOT.TFile.Open(self.filename, 'OPEN')
        return new_rows

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
                nbins_tot = hsub.GetNbinsX()*hsub.GetNbinsY()
                if binning.is3D:
                    nbins_tot *= hsub.GetNbinsZ()
                for b in range(1, nbins_tot+1):
                    hsub.SetBinContent(b, 1e-6)
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
        # Check if user requests files be read remotely via xrootd
        source_filename_up = syst_dict['UP'].split(':')[0] if 'root://' not in syst_dict['UP'] else syst_dict['UP'].split(':')[0]+':'+syst_dict['UP'].split(':')[1]
        source_histname_up = syst_dict['UP'].split(':')[1] if 'root://' not in syst_dict['UP'] else syst_dict['UP'].split(':')[2]
        source_filename_dn = syst_dict['DOWN'].split(':')[0] if 'root://' not in syst_dict['DOWN'] else syst_dict['DOWN'].split(':')[0]+':'+syst_dict['DOWN'].split(':')[1]
        source_histname_dn = syst_dict['DOWN'].split(':')[1] if 'root://' not in syst_dict['DOWN'] else syst_dict['DOWN'].split(':')[2]
        out = [
            {
                'shapes':syst_dict['SIGMA'],
                'syst_type': 'shapes',
                'source_filename': source_filename_up,
                'source_histname': source_histname_up,
                'direction': 'Up',
                'variation_alias': name if 'ALIAS' not in syst_dict else syst_dict['ALIAS']
            }, {
                'shapes':syst_dict['SIGMA'],
                'syst_type': 'shapes',
                'source_filename': source_filename_dn,
                'source_histname': source_histname_dn,
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
    dupes = df[df.duplicated(subset=['process','region','variation','source_filename','source_histname','direction'],keep=False)]
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
