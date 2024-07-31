import glob
import ROOT, os, warnings, pandas, math, time
from PIL import Image
from collections import OrderedDict
from TwoDAlphabet.helpers import set_hist_maximums, execute_cmd, cd, hist2array
from TwoDAlphabet.binning import stitch_hists_in_x, convert_to_events_per_unit, get_min_bin_width
from TwoDAlphabet.ext import tdrstyle, CMS_lumi
from TwoDAlphabet.plotstyle import * # dictionaries containing plotting styles for mplhep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep

class Plotter(object):
    '''Class to manage output distributions, manipulate them, and provide access to plotting
    standard groups of distributions.

    Attributes:
        fit_tag (str): Either 's' or 'b'.
        fit_results (RooFitResult): The RooFitResult corresponding to the fit_tag.
        signal_strength (float): The post-fit signal strength.
        twoD (TwoDAlphabet): TwoDAlphabet object storing various meta information needed for access.
        yaxis1D_title (str): Title for "counts" axis of 1D plots. Defaults to 'Events / bin' but can change if plotting events per unit.
        df (pandas.DataFrame): DataFrame organizing the post-fit and pre-fit plots and their 1D projections.
        dir (str): Directory path to save final images.
        slices (dict): Stores edges to slice "x" and "y" axes. 
        root_out (ROOT.TFile): File storing all histograms that are made.
    '''
    def __init__(self,ledger,twoD,fittag,loadExisting=False):
        '''Constructor.

        Args:
            twoD (TwoDAlphabet): Object with meta information about the run.
            fittag (str): Either 's' or 'b'.
            loadExisting (bool, optional): Flag to load existing projections instead of remaking everything. Defaults to False.
        '''

        self.fittag = fittag
        self.twoD = twoD
        self.ledger = ledger
        self.yaxis1D_title = 'Events / bin'
        self.df = pandas.DataFrame(columns=['process','region','process_type','title'])
        self.dir = 'plots_fit_{f}'.format(f=self.fittag)
        self.slices = {'x': {}, 'y': {}}
        self.root_out = None
        self.df_path = None

        if not loadExisting:
            self._make()
        else:
            self._load()

    def __del__(self):
        '''On deletion, save DataFrame to csv and close ROOT file.'''
        if self.df_path is not None:
            self.df.to_csv(self.df_path)
        if self.root_out is not None and hasattr(self.root_out, 'Close'):
            self.root_out.Close()

    def _load(self):
        '''Open pickled DataFrame and output ROOT file
        and reference with `self.df` and `self.root_out` attributes.'''
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name)
        
        df_csv_path = os.path.join(self.dir, 'df.csv')
        self.df = pandas.read_csv(df_csv_path)
        self.df_path = df_csv_path

    def _format_1Dhist(self, hslice, title, xtitle, ytitle, color, proc_type):
        '''Perform some basic formatting of a 1D histogram so that the ROOT.TH1
        already has some of the meta information set like line and fill colors.
        Also renormalizes bins to the minimum bin width if the y-axis is requested
        to be plotted as events per unit rather than events per bin.

        Args:
            hslice (ROOT.TH1): Histogram.
            title (str): Title for the output histogram.
            xtitle (str): X-axis title for the output histogram.
            ytitle (str): Y-axis title for the output histogram. Will be modified if requesting to plot events/unit.
            color (int): ROOT color code.
            proc_type (str): Process type. Either 'BKG', 'SIGNAL', or 'DATA'.

        Raises:
            NameError: If proc_type is not 'BKG', 'SIGNAL', or 'DATA'.

        Returns:
            TH1: Formatted histogram.
        '''
        ytitle = self.yaxis1D_title
        if self.twoD.options.plotEvtsPerUnit:
            hslice = convert_to_events_per_unit(hslice)
            ytitle = 'Events / %s GeV' % get_min_bin_width(hslice)
        hslice.SetMinimum(0)
        hslice.SetTitle(title)
        hslice.GetXaxis().SetTitle(xtitle)
        hslice.GetYaxis().SetTitle(ytitle)

        # Obtain the ROOT color code from the dictionary defined in TwoDAlphabet.plotstyle
        if color not in mpl_to_root_colors.keys():
            available_colors = '", "'.join(mpl_to_root_colors.keys())
            raise ValueError(f'Color "{color}" not defined. Please add the ROOT TColor code to the "mpl_to_root_colors" dictionary defined in TwoDAlphabet.plotstyle. Available default colors are: "{available_colors}"')
        else:
            color = int(mpl_to_root_colors[color]) #ROOT call in C++ sometimes cannot convert it to int

        if proc_type == 'BKG':
            hslice.SetFillColor(color)
            hslice.SetLineColorAlpha(0,0)
        elif proc_type == 'SIGNAL' or proc_type == 'TOTAL':
            hslice.SetLineColor(color)
        elif proc_type == 'DATA':
            hslice.SetLineColor(color)
            hslice.SetMarkerColor(color)
        else:
            raise NameError('Process type "%s" is not supported.'%proc_type)

        return hslice

    def _make(self):
        '''Make the DataFrame and output ROOT file from scratch
        and reference with `self.df` and `self.root_out` attributes.
        
        Loops over all regions and processes from the pre-fit and post-fit shapes
        and tracks/constructs the 2D histograms and six projections (three each for "x" and "y"). 
        '''
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name,'RECREATE')

        shapes_file = ROOT.TFile.Open('postfitshapes_%s.root'%self.fittag)
        loc_base = '{r}_{c}_{t}/{p}'

        proc_reg_pairs = self.ledger.GetProcRegPairs()+[('TotalBkg', r) for r in self.ledger.GetRegions()]

        for region in self.ledger.GetRegions():
            binning,_ = self.twoD.GetBinningFor(region)
            if len(self.twoD.options.blindedPlots) > 0:
                if region in self.twoD.options.blindedPlots:
                    blinding = [1]
                else: blinding = []
            else: blinding = []

            self.slices['x'][region] = {'vals': binning.xSlices,'idxs':binning.xSliceIdx}
            self.slices['y'][region] = {'vals': binning.ySlices,'idxs':binning.ySliceIdx}
            
            for process in self.ledger.GetProcesses()+['TotalBkg']:
                # Skip processes not in this region
                if process not in [pair[0] for pair in proc_reg_pairs if pair[1] == region]:
                    continue
                
                if process != 'TotalBkg':
                    color = self.ledger.GetProcessColor(process)
                    proc_type = self.ledger.GetProcessType(process)
                    proc_title = self.ledger.GetProcessTitle(process)
                else:
                    #color = ROOT.kBlack
                    color = 'black'
                    proc_type = 'TOTAL'
                    proc_title = 'TotalBkg'

                self.df = pandas.concat([self.df,pandas.DataFrame([{'process':process,
                                          'region':region,
                                          'process_type': proc_type,
                                          'title': proc_title}])], ignore_index=True)

                for time in ['prefit','postfit']:
                    # 2D distributions first
                    out2d_name = '%s_%s_%s_2D'%(process,region,time)

                    low_name = loc_base.format(r=region, c='LOW', t=time, p=process)
                    sig_name = loc_base.format(r=region, c='SIG', t=time, p=process)
                    high_name = loc_base.format(r=region, c='HIGH',t=time, p=process)
                    
                    low  = shapes_file.Get(low_name)
                    sig  = shapes_file.Get(sig_name)
                    high = shapes_file.Get(high_name)

                    if low == None: raise IOError('Could not find histogram %s in postfitshapes_%s.root'%(low_name, self.fittag))
                    if sig == None: raise IOError('Could not find histogram %s in postfitshapes_%s.root'%(sig_name, self.fittag))
                    if high == None: raise IOError('Could not find histogram %s in postfitshapes_%s.root'%(high_name, self.fittag))

                    full = stitch_hists_in_x(out2d_name, binning, [low,sig,high], blinded=blinding if process == 'data_obs' else [])
                    full.SetMinimum(0)
                    full.SetTitle('%s, %s, %s'%(proc_title,region,time))

                    self.root_out.WriteTObject(full,full.GetName())

                    # Now do projections using the 2D
                    out_proj_name = '{p}_{r}_{t}_proj{x}{i}'
                    for proj in ['X','Y']:
                        slices = self.slices['x' if proj == 'Y' else 'y'][region]

                        for islice in range(3):
                            hname = out_proj_name.format(p=process,r=region,t=time,x=proj.lower(),i=islice)
                            start,stop = _get_start_stop(islice,slices['idxs'])
                            
                            hslice = getattr(full,'Projection'+proj)(hname,start,stop,'e')
                            hslice_title = '%s, %s, %s, %s-%s'%(proc_title,region,time,slices['vals'][islice],slices['vals'][islice+1])
                            hslice = self._format_1Dhist(
                                hslice, hslice_title,
                                binning.xtitle if proj == 'X' else binning.ytitle,
                                self.yaxis1D_title,
                                color, proc_type)

                            self.root_out.WriteTObject(hslice,hslice.GetName())

        shapes_file.Close()
        self.root_out.Close()
        self.root_out = ROOT.TFile.Open(root_out_name)

    def Get(self,hname=None,row=None,hist_type=None):
        '''Get a histogram by name from the master ROOT file.
        Does quick check for if the histogram is saved.
        
        Args:
            hname (str): Histogram name.

        Raises:
            LookupError: If histogram cannot be found.
        '''
        if hname != None:
            if hname not in [k.GetName() for k in self.root_out.GetListOfKeys()]:
                raise LookupError('Histogram %s not found in %s'%(hname,self.root_out.GetName()))
            name = hname
        else:
            name = '_'.join([row.process,row.region,hist_type])
        return self.root_out.Get(name)

    def _order_df_on_proc_list(self,df,proc_type,proclist=[],alphaBottom=True):
        '''Re-order input dataframe (`df`) based on the ordered list of process names (`proclist`).
        Useful for pre-ordering process before trying to construct a THStack.

        Args:
            df (pandas.DataFrame): Input DataFrame to manipulate.
            proc_type (str): Process type.  Either 'BKG', 'SIGNAL', or 'DATA'.
            proclist (list, optional): Ordered list of processes to order by. Defaults to [] in which case order is determined based on alphaBottom.
            alphaBottom (bool, optional): Only matters if proclist == []. Defaults to True in which case parametric Alphabet objects included first (thus, on the bottom of the THStack).

        Returns:
            pandas.DataFrame: Ordered DataFrame.
        '''
        if proclist == []:
            if alphaBottom:
                process_order = self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=True, includeConfig=False) + self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=False, includeConfig=True)
            else:
                process_order = self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=False, includeConfig=True) + self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=True, includeConfig=False)
        else:
            process_order = proclist

        process_order_df = pandas.DataFrame({'process':process_order})
        return process_order_df.merge(df,on='process',how='inner')

    def plot_2D_distributions(self):
        '''Take the saved 2D distributions and plot them together on sub-pads
        based on process and region groupings.

        Plots are grouped based on process and then the regions and pre-fit/post-fit
        plots share the same canvas as sub-pads.

        Returns:
            None
        '''
        for pr, _ in self.df.groupby(['process','region']):
            process, region = pr[0], pr[1]
            out_file_name = '{d}/base_figs/{p}_{r}_%s_2D'.format(d=self.dir,p=process,r=region)
            make_pad_2D(outname=out_file_name%('prefit'), hist=self.Get('{p}_{r}_{t}'.format(p=process,r=region,t='prefit_2D')),
                            year=self.twoD.options.year, savePDF=True, savePNG=True)
            make_pad_2D(outname=out_file_name%('postfit'), hist=self.Get('{p}_{r}_{t}'.format(p=process,r=region,t='postfit_2D')),
                            year=self.twoD.options.year, savePDF=True, savePNG=True)

            make_can('{d}/{p}_{r}_2D'.format(d=self.dir,p=process,r=region), [out_file_name%('prefit')+'.png', out_file_name%('postfit')+'.png'])

    def plot_projections(self, lumiText=r'138 $fb^{-1}$ (13 TeV)', extraText='Preliminary', subtitles={}, units='GeV', regionsToGroup=[]):
        '''Plot comparisons of data and the post-fit background model and signal
        using the 1D projections. Canvases are grouped based on projection axis.
        The canvas rows are separate selection regions while the columns
        are the different slices of the un-plotted axis.

        Args:
            lumiText (string): LaTeX-formatted string containing luminosity information. Defaults to Run2 conditions.
            extraText (str): Additional text to place after experiment (CMS) text.
            subtitles ({str:str}, optional): Dict of raw strings corresponding to each region specified in the JSON to be placed underneath axis slice text in top left corner of pad. If multiple titles are desired, separate with semicolon character.
                Example: {"SR_fail": r"$ParticleNet_{TvsQCD}$ Pass\n$ParticleNetMD_{Xbb}$ Fail", "SR_pass": r"$ParticleNet_{TvsQCD}$ Pass;$ParticleNetMD_{Xbb}$ Pass"}
            units (str): Units of measurement for observable. Placed on x-axis and in slice string.

            regionsToGroup ([[str]]): 
                List of list of strings representing the desired regions to group. For example if the fit involved
                four regions: CR_fail, CR_pass, SR_fail, SR_pass then 2DAlphabet will try to plot all
                (4 regions) x (3 slices) = 12 plots on the same page, which is greater than the number that can be
                plotted. Instead, pass regionsToGroup = [['CR'],['SR']] to have the CR and SR plotted on separate
                2x3 canvases. 
                If you wanted to plot, e.g. SR_fail, SR_pass, and CR_pass, pass in the following list of lists:
                    [['CR'], ['SR'], ['SR_fail','SR_pass','CR_pass']]
                The order matters if the sub-list contains multiple strings - the regions are plotted in order
                with the first on top and last on bottom of the canvas.

        Returns:
            None
        '''
        axes = pandas.DataFrame() # Book a dataframe for creating a full group of plots
        for region, group in self.df.groupby('region'):
            binning, _ = self.twoD.GetBinningFor(region)

            # Make both regular and logarithmic y-axis plots
            for logyFlag in [False, True]:
                # Get reduced dataframes for all backgrounds and the signals
                ordered_bkgs = self._order_df_on_proc_list(
                                            group[group.process_type.eq('BKG')], proc_type='BKG',
                                            alphaBottom=(not logyFlag))
                signals = group[group.process_type.eq('SIGNAL')]

                for proj in ['postfit_projx','postfit_projy']:
                    for islice in range(3):
                        projn     = f'{proj}{islice}'
                        sig_projn = projn
                        if self.twoD.options.plotPrefitSigInFitB and self.fittag == 'b':
                            sig_projn = projn.replace('postfit','prefit') # Plot prefit signal in b-only plots

                        # TH1s representing the sum of all data and background histograms in the workspace for this projection and slice
                        this_data =      self.Get(row=group.loc[group.process_type.eq('DATA')].iloc[0], hist_type=projn)
                        this_totalbkg =  self.Get(row=group.loc[group.process_type.eq('TOTAL')].iloc[0], hist_type=projn)
                        # Lists of individual TH1s comprising all bkg and sig histograms for this projection and slice.
                        these_bkgs =    [self.Get(row=ordered_bkgs.iloc[irow], hist_type=projn) for irow in range(ordered_bkgs.shape[0])]
                        these_signals = [self.Get(row=signals.iloc[irow], hist_type=sig_projn) for irow in range(signals.shape[0])]

                        slice_edges = (
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice],
                            binning.xtitle if 'y' in proj else binning.ytitle,
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice+1],
                        )
                        slice_str = '%s < %s < %s'%slice_edges

                        out_pad_name = f'{self.dir}/base_figs/{projn}_{region}{"" if not logyFlag else "_logy"}'

                        # If user requested sub-titles, obtain the ones for this specific region
                        if subtitles:
                            subtitle = subtitles[region]
                        else:
                            subtitle = {}

                        # Produce a matplotlib axis containing the full data + stacked bkg histogram for this projection and slice
                        make_ax_1D(
                            out_pad_name, 
                            binning, 
                            data=this_data, 
                            bkgs=these_bkgs, 
                            signals=these_signals, 
                            totalBkg=this_totalbkg, 
                            logyFlag=logyFlag, 
                            slicetitle=slice_str,   # opposite axis slices (top left, in pad)
                            subtitle=subtitle,       # (top left, under slice string)
                            lumiText=lumiText,      # (top right, above pad)
                            extraText=extraText,    # top left, above pad, next to "CMS"
                            units=units,
                            savePDF=True,
                            savePNG=True
                        )

                        # Append projection plot to master dataframe
                        axes = pandas.concat([axes, pandas.DataFrame([{'ax':out_pad_name+'.png', 'region':region, 'proj':projn, 'logy':logyFlag}])], ignore_index=True)

        # Create a canvas with the full set of projection plots
        for logy in ['', '_logy']:
            for proj in ['postfit_projx','postfit_projy']:
                these_axes = axes.loc[axes.proj.str.contains(proj)]
                if logy == '':
                    these_axes = these_axes.loc[these_axes.logy.eq(False)]
                else:
                    these_axes = these_axes.loc[these_axes.logy.eq(True)]

                if (len(these_axes) > 9) and (len(regionsToGroup) == 0):
                    raise RuntimeError('histlist of size %s not currently supported. Instead, call plot_projections() with regionsToGroup list describing the regions you want to group together.'%len(these_axes))
                elif (len(these_axes) > 9) and (len(regionsToGroup) > 0):
                    validRegions = these_axes['region'].to_list()
                    for i in range(len(regionsToGroup)):
                        region = regionsToGroup[i]
                        if (len(region) > 1):	# e.g. ['SR_fail', 'SR_pass', 'ttbarCR_pass']
                            # Ensure that the regions are plotted in order
                            new_axes = []
                            for r in region:
                                if r not in validRegions:
                                    raise ValueError(f'Region "{r}" specified in regionsToGroup is not available in the 2DAlphabet workspace. Available regions:\n\t{validRegions}')
                                r_axes = these_axes[these_axes['region'].str.match(r)].sort_values(by=['region','proj'])['ax'].to_list()
                                new_axes += r_axes
                            '''
                            new_axes = []
                            # Ensure that the regions are plotted in order provided
                            for r in region:
                                for ax in these_axes:
                                    # As an edge case, assume r would be: CR_pass
                                    # pad would be: plots_fit_b/base_figs/postfit_projx2_CR_pass.png
                                    # but also: plots_fit_b/base_figs/postfit_projx2_ttbarCR_pass.png
                                    # so we have to append an underscore to ensure we get the right one 
                                    rNew = '_'+r
                                    if rNew in ax:
                                        new_axes.append(ax)
                            '''
                            out_can_name = '{d}/{reg}_{proj}{logy}'.format(d=self.dir,proj=proj,logy=logy,reg='_and_'.join(region))
                            make_can(out_can_name, new_axes)
                        else:	# e.g. ['SR']
                            new_axes = these_axes[these_axes['region'].str.match(f'{region[0]}_')].sort_values(by=['region','proj'])['ax'].to_list()
                            out_can_name = '{d}/{reg}_{proj}{logy}'.format(d=self.dir,proj=proj,logy=logy,reg=region[0])
                            make_can(out_can_name, new_axes)
                else:
                    these_axes = these_axes.sort_values(by=['region','proj'])['ax'].to_list()
                    out_can_name = '{d}/{proj}{logy}'.format(d=self.dir, proj=proj, logy=logy)
                    make_can(out_can_name, these_axes)


    def plot_transfer_funcs(self):
        raise NotImplementedError()

def _save_pad_generic(outname, pad, ROOTout, savePDF, savePNG):
    if isinstance(ROOTout, ROOT.TFile):
        ROOTout.WriteTObject(pad,outname.split('/')[-1])
    if savePDF:
        pad.Print(outname+'.pdf','pdf')
    if savePNG:
        pad.Print(outname+'.png','png')

def _make_pad_gen(name):
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleAlign(33)
    ROOT.gStyle.SetTitleX(.77)

    pad = ROOT.TCanvas(name, name, 800, 700)
    pad.cd(); pad.SetRightMargin(0.0); pad.SetTopMargin(0.0); pad.SetBottomMargin(0.0)
    return pad

def make_pad_2D(outname, hist, style='lego', logzFlag=False, ROOTout=None,
                savePDF=False, savePNG=False, year=1, extraText='Preliminary'):
    '''Make a pad holding a 2D plot with standardized formatting conventions.

    Args:
        outname (str): Output file path name.
        hist (TH2): Histogram to draw on the pad.
        style (str, optional): ROOT drawing style. Defaults to 'lego'.
        logzFlag (bool, optional): Make log z-axis. Defaults to False.
        saveROOT (bool, optional): Save to master ROOT file. Defaults to True.
        savePDF (bool, optional): Save to PDF. Defaults to False.
        savePNG (bool, optional): Save to PNG. Defaults to False.
        year (int, optional): Luminosity formatting. Options are 16, 17, 18, 1 (full Run 2), 2 (16+17+18). Defaults to 1.
        extraText (str, optional): Prepended to the CMS subtext. Defaults to 'Preliminary'.

    Returns:
        [type]: [description]
    '''
    pad = _make_pad_gen(outname)
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.2)
    pad.SetBottomMargin(0.12)
    pad.SetTopMargin(0.1)
    if logzFlag: pad.SetLogz()

    hist.GetXaxis().SetTitleOffset(1.15); hist.GetXaxis().SetLabelSize(0.05); hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleOffset(1.5);  hist.GetYaxis().SetLabelSize(0.05); hist.GetYaxis().SetTitleSize(0.05)
    hist.GetZaxis().SetTitleOffset(1.5);  hist.GetZaxis().SetLabelSize(0.05); hist.GetZaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetNdivisions(505)
    
    if 'lego' in style.lower():
        hist.GetZaxis().SetTitleOffset(1.4)

    hist.Draw(style)
    
    CMS_lumi.extraText = extraText
    CMS_lumi.CMS_lumi(pad, year, 11, sim=False if 'data' in hist.GetName().lower() else True)

    ROOT.SetOwnership(pad, False)
    _save_pad_generic(outname, pad, ROOTout, savePDF, savePNG)

    return pad


def make_ax_1D(outname, binning, data, bkgs=[], signals=[], title='', subtitle='', slicetitle='',
            totalBkg=None, logyFlag=False, ROOTout=None, savePDF=False, savePNG=False,
            dataOff=False, datastyle='pe X0', year=1, addSignals=True, 
            lumiText=r'$138 fb^{-1} (13 TeV)$', extraText='Preliminary', units='GeV', hspace=0.0):
    '''Create a matplotlib.axis.Axis object holding a 1D plot with standardized CMS formatting conventions
    Args:
        outname (str): Output file path + name.
        binning (Binning): TwoDAlphabet binning object for the given region from which to obtain binning info.
        data (TH1): Data histogram.
        bkgs ([TH1]): List of background histograms (to be stacked).
        signals ([TH1]): List of signal histograms.
        title (str, optional): Title of plot. Only applicable if bkgs is empty. Defaults to ''.
        subtitle (str, optional): Subtitle text for physics information (like slice ranges). Defaults to ''.
        totalBkg (TH1, optional): Total background estimate from fit. Used to get total background uncertianty. Defaults to None.
        logyFlag (bool, optional): Make log y-axis. Defaults to False.
        savePDF (bool, optional): Save to PDF. Defaults to True.
        savePNG (bool, optional): Save to PNG. Defaults to True.
        dataOff (bool, optional): Turn off the data from plotting. Defaults to False.
        extraText (str, optional): Prepended to the CMS subtext. Defaults to 'Preliminary'.
        units (str, optional): Units of measurement for x- and y-axes. Defaults to 'GeV'.
        hspace (float, optional): Spacing between the main plotting axis and the ratio subplot axis. Defaults to 0.0 (no spacing). Original 2DAlphabet spacing would correspond to hspace=0.7.
    Returns:
        ax (matplotlib.axis.Axis): Output axis object for further manipulation.
    '''
    # Convert all histograms to numpy arrays. First, determine the bin edges
    projn  = outname.split('/')[-1].split('_')[1].split('proj')[-1][0]
    islice = outname.split('/')[-1].split('_')[1].split('proj')[-1][1]; islice = int(islice)
    if projn == 'x':
        xbins = binning.xbinByCat
        edges = np.array(xbins['LOW'][:-1]+xbins['SIG'][:-1]+xbins['HIGH'])
    else:
        edges = np.array(binning.ybinList)
    widths = np.diff(edges)     # obtain bin widths

    # Convert all unique bkg and signal hists to arrays and group by process. They will come already ordered, so keep the ordering
    bkgDict = OrderedDict()
    sigDict = OrderedDict()
    bkgNames = list(dict.fromkeys([hist.GetTitle().split(',')[0] for hist in bkgs]))
    sigNames = list(dict.fromkeys([hist.GetTitle().split(',')[0] for hist in signals]))
    # Replace the ROOT latex "#" with standard latex "\" escape character for python rstring
    bkgNamesLatex = [r'${}$'.format(bkgName.replace("#","\\")) for bkgName in bkgNames]
    sigNamesLatex = [r'${}$'.format(sigName.replace("#","\\")) for sigName in sigNames]
    # Sum the common backgrounds and signals
    for bkg in bkgNames:
        bkg_arrs = [hist2array(i, return_errors=True)[0] for i in bkgs if bkg in i.GetTitle()]
        bkg_errs = [hist2array(i, return_errors=True)[1] for i in bkgs if bkg in i.GetTitle()]
        bkgDict[f'{bkg}_arr'] = sum(bkg_arrs)
        bkgDict[f'{bkg}_err'] = sum(bkg_errs)
    for sig in sigNames:
        sigDict[f'{sig}_arr'] = sum([hist2array(i) for i in signals if sig in i.GetTitle()]) # ignore uncertainty
    # For the colors, we need to translate from ROOT TColor to matplotlib 
    bkgColors = list(dict.fromkeys([hist.GetFillColor() for hist in bkgs]))
    bkgColors = [root_to_matplotlib_color(TColor) for TColor in bkgColors]
    sigColors = list(dict.fromkeys([hist.GetLineColor() for hist in signals]))  # signals have no fill, only line color
    sigColors = [root_to_matplotlib_color(TColor) for TColor in sigColors]

    # Get the data array and total bkg
    data_arr = hist2array(data); data_arr = np.array([int(i) for i in data_arr]) # hist2array converts to floats which may leave very small differences when converting TH1 (int) -> array (float). This fixes it
    if totalBkg:
        totalBkg_arr, totalBkg_err = hist2array(totalBkg, return_errors=True)

    # Begin plotting
    plt.style.use([hep.style.CMS])
    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style) # ax is the stacked plot, rax contains the ratio plot
    fig.subplots_adjust(hspace=hspace)

    bkg_stack = np.vstack([arr for key, arr in bkgDict.items() if '_arr' in key]) # stack all unique background processes
    # depending on step option ('pre' or 'post'), the last bin needs be concatenated on one side, so that the edge bin is drawn (annoying)
    bkg_stack = np.hstack([bkg_stack, bkg_stack[:,-1:]])
    bkg_stack = np.hstack([bkg_stack])

    ax.stackplot(edges, bkg_stack, labels=bkgNamesLatex, colors=bkgColors, step='post', **stack_style)
    unc = totalBkg_err # the hist2array() call returns the sqrt of the sumw2 for the histogram
    unc = np.hstack([unc, unc[-1]]) # concatenation of last bin so that edge bin is drawn
    ax.fill_between(x=edges, y1=bkg_stack.sum(axis=0)-unc, y2=bkg_stack.sum(axis=0)+unc, label='Bkg. Unc.', step='post', **hatch_style)
    bin_centers = (edges[:-1] + edges[1:])/2
    # Determine if we have variable binning, and if so, add the CMS-required horizontal bars to denote bin width
    if len(np.unique(widths)) > 1: # Detected variable bin widths
        xerrs = (edges[1:]-edges[:-1])/2
        ax.set_ylabel('Events / bin')
    else:
        ax.set_ylabel(f'Events / {widths[0]} {units}')
        xerrs = None
    ax.errorbar(x=bin_centers, y=data_arr, yerr=np.sqrt(data_arr), xerr=xerrs, label='Data', **errorbar_style)

    # Plot signals
    for i, sig in enumerate(sigNames):
        sigarr = sigDict[f'{sig}_arr']
        ax.step(x=edges, y=np.hstack([sigarr, sigarr[-1]]), where='post', color=sigColors[i], label=sigNamesLatex[i])

    if logyFlag:
        ax.set_ylim(0.01, totalBkg_arr.max()*1e5)
        ax.set_yscale('log')
    else:
        ax.set_ylim(0, totalBkg_arr.max()*1.38)

    # Make sure data and signal(s) come first 
    handles, labels = ax.get_legend_handles_labels()
    data_idx = labels.index('Data') # should always be last entry in legend
    unc_idx  = labels.index('Bkg. Unc.')
    sig_idxs = [i for i in range(unc_idx+1, data_idx)] # indices of all signals in legend
    # Resulting legend ordering will be Data, signal(s), bkgs, then bkg unc. The bkg ordering is already sent to this function properly.
    leg_order = [data_idx] + sig_idxs + [idx for idx in range(len(labels)) if idx not in [data_idx] + sig_idxs]
    ax.legend([handles[idx] for idx in leg_order],[labels[idx] for idx in leg_order])
    ax.autoscale(axis='x', tight=True)
    ax.margins(x=0) # remove white space at left and right margins of plot 

    hep.cms.label(loc=0, ax=ax, data = not dataOff, label=extraText, rlabel='') # CMS + label, where label is typically “Preliminary” “Supplementary”, “Private Work” or “Work in Progress”
    hep.cms.lumitext(lumiText, ax=ax)                       # Typically luminosity + sqrt(s)
    # Can't use the hep.cms.text() wrapper without "CMS" being added, so add the slice text manually
    if slicetitle:
        slicetitle = r'${}$ {}'.format(slicetitle.replace('#','\\'), units)
        ax.text(0.3, 0.95, slicetitle, ha='center', va='top', fontsize='small', transform=ax.transAxes)
    if subtitle:
        # Check if user requested multi-line
        if (len(subtitle.split(';')) == 1):
            ax.text(0.3, 0.89, r'%s'%subtitle, ha='center', va='top', fontsize='small', transform=ax.transAxes)
        else:
            for i, title in enumerate(subtitle.split(';')):
                ax.text(0.3, 0.95-(0.06*(i+1)), r'%s'%title, ha='center', va='top', fontsize='small', transform=ax.transAxes)

    # pull
    dataMinusBkg = data_arr - totalBkg_arr
    sigmas = np.sqrt(np.sqrt(data_arr)*np.sqrt(data_arr) + totalBkg_err*totalBkg_err)
    sigmas[sigmas==0.0] = 1e-5 # avoid division by zero 
    rax.bar(bin_centers, dataMinusBkg/sigmas, width=widths, color='gray')
    rax.set_ylim(-3,3)
    rax.set_ylabel(r'$\frac{Data-Bkg}{\sigma}$')
    axisTitle = binning.xtitle if projn == 'x' else binning.ytitle
    axisTitle = axisTitle.replace("#","\\")
    rax.set_xlabel(r'${}$ [{}]'.format(axisTitle, units))
    rax.autoscale(axis='x', tight=True)
    rax.margins(x=0)

    if savePDF:
        plt.savefig(f'{outname}.pdf')
    if savePNG:
        plt.savefig(f'{outname}.png')
    if ((not savePDF) and (not savePNG)):
        print(f'WARNING: plot "{outname}" has not been saved.')

    print(f'Plotting {outname}')
    plt.close()
    # return ax, rax      ????????????

def make_can(outname, padnames, padx=0, pady=0):
    '''Combine multiple pads/canvases into one canvas for convenience of viewing.
    Input pad order matters.

    Args:
        outname (str): Output file path name.
        pads ([TCanvas,TPad]): List of canvases/pads to plot together on one canvas.
        saveROOT (bool, optional): Save to master ROOT file. Defaults to False.
        savePDF (bool, optional): Save to PDF. Defaults to True.
        savePNG (bool, optional): Save to PNG. Defaults to True.

    Raises:
        RuntimeError: If 10 or more subdivisions are requested.

    Returns:
        None
    '''
    if padx == 0 or pady == 0:
        if len(padnames) == 1:
            padx = 1; pady = 1
        elif len(padnames) == 2:
            padx = 2; pady = 1
        elif len(padnames) == 3:
            padx = 3; pady = 1
        elif len(padnames) == 4:
            padx = 2; pady = 2
        elif len(padnames) <= 6:
            padx = 3; pady = 2
        elif len(padnames) <= 9:
            padx = 3; pady = 3
        else:
            raise RuntimeError('histlist of size %s not currently supported'%len(padnames))
            #raise RuntimeError('histlist of size %s not currently supported: %s'%(len(padnames),[p.GetName() for p in padnames]))

    pads = [Image.open(os.path.abspath(pname)) for pname in padnames]
    w, h = pads[0].size
    grid = Image.new('RGB', size=(padx*w, pady*h))
    
    for i, pad in enumerate(pads):
        grid.paste(pad, box=(i%padx*w, i//padx*h))
    
    print ('Writing grid of images %s.pdf'%outname)
    grid.save(outname+'.pdf')

def _get_start_stop(i,slice_idxs):
    start = slice_idxs[i]+1
    stop  = slice_idxs[i+1]
    return start, stop

def gen_projections(ledger, twoD, fittag, loadExisting=False, lumiText=r'138 $fb^{-1}$ (13 TeV)', extraText='Preliminary', subtitles={}, units='GeV', regionsToGroup=[]):
    plotter = Plotter(ledger, twoD, fittag, loadExisting)
    plotter.plot_2D_distributions()
    plotter.plot_projections(lumiText, extraText, subtitles, units, regionsToGroup)

def make_systematic_plots(twoD):
    '''Make plots of the systematic shape variations of each process based on those
    processes and systematic shapes specified in the config. Shapes are presented
    as projections onto 1D axis where no selection has been made on the axis not
    being plotted. Plots are saved to UncertPlots/.
    '''
    for (p,r), _ in twoD.df.groupby(['process','region']):
        if p == 'data_obs': continue

        nominal_full = twoD.organizedHists.Get(process=p, region=r, systematic='')
        binning, _ = twoD.GetBinningFor(r)

        for axis in ['X','Y']:

            nominal_hist = getattr(nominal_full,'Projection'+axis)('%s_%s_%s_%s'%(p,r,'nom','proj'+axis))
            nominal = hist2array(nominal_hist)

            # Get the bin edges from the 2DAlphabet binning object. Avoid edge duplication in the case of X-axis stitching
            if axis == 'X':
                xbins = binning.xbinByCat
                edges = xbins['LOW'][:-1]+xbins['SIG'][:-1]+xbins['HIGH']
            else:
                edges = binning.ybinList

            # GetShapeSystematics() lists all shape systs, even if the given process is not subject to them.
            check_df = twoD.df.loc[(twoD.df.process == p) & (twoD.df.region == r)] # Get the entries for this proc and region
            proc_vars = check_df.variation.unique() # Get all unique variations for this process

            for s in twoD.ledger.GetShapeSystematics(drop_norms=True):
                if s not in proc_vars: continue

                up_hist = getattr(twoD.organizedHists.Get(process=p,region=r,systematic=s+'Up'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Up','proj'+axis))
                up = hist2array(up_hist)
                down_hist = getattr(twoD.organizedHists.Get(process=p,region=r,systematic=s+'Down'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Down','proj'+axis))
                down = hist2array(down_hist)

                # Begin plotting
                labels = ['Nominal', r'$+1\sigma$', r'$-1\sigma$']
                colors = ['black', 'red', 'blue']
                styles = ['solid', 'dashed', 'dashed']
                histos = [nominal, up, down]

                plt.style.use([hep.style.CMS])
                fig, ax = plt.subplots(figsize=(12, 8), dpi=200)
                hep.cms.text("WiP",loc=1)
                ax.set_title(f'{p}, {s}', pad=12)
                hep.histplot(histos, edges, stack=False, ax=ax, label=labels, histtype='step', linestyle=styles, color=colors)
                handles, labelsproj = ax.get_legend_handles_labels()
                ax.set_xlabel(getattr(binning,axis.lower()+'title'))
                ax.set_ylabel('Events')
                plt.legend(loc='best')

                outname = f'{twoD.tag}/UncertPlots/Uncertainty_{p}_{r}_{s}_proj{axis}.png'
                print(f'[2DAlphabet.plot] INFO: Plotting histogram\n\t{outname}')

                plt.savefig(outname)
                plt.close() # free up memory

def _make_pull_plot(data, bkg, preVsPost=False):
    pull = data.Clone(data.GetName()+"_pull")
    pull.Add(bkg,-1)
    
    sigma = 0.0
    for ibin in range(1,pull.GetNbinsX()+1):
        d = data.GetBinContent(ibin)
        b = bkg.GetBinContent(ibin)
        if d >= b:
            derr = data.GetBinErrorLow(ibin)
            berr = bkg.GetBinErrorUp(ibin)
        elif d < b:
            derr = data.GetBinErrorUp(ibin)
            berr = bkg.GetBinErrorLow(ibin)
        
        if d == 0:
            derr = 1

        sigma = math.sqrt(derr*derr + berr*berr)
        if sigma != 0:
            pull.SetBinContent(ibin, (pull.GetBinContent(ibin))/sigma)
        else:
            pull.SetBinContent(ibin, 0.0 )

    pull.SetFillColor(ROOT.kBlue)
    pull.SetTitle(";"+data.GetXaxis().GetTitle()+";({})/#sigma".format('Post-Pre' if preVsPost else 'Data-Bkg'))
    pull.SetStats(0)

    pull.GetYaxis().SetRangeUser(-2.9,2.9)
    pull.GetYaxis().SetTitleOffset(0.4)                             
    pull.GetYaxis().SetLabelSize(0.13)
    pull.GetYaxis().SetTitleSize(0.12)
    pull.GetYaxis().SetNdivisions(306)

    pull.GetXaxis().SetLabelSize(0.13)
    pull.GetXaxis().SetTitleSize(0.15)
    pull.GetYaxis().SetTitle("({})/#sigma".format('Post-Pre' if preVsPost else 'Data-Bkg'))
    return pull

def _get_good_fit_results(tfile):
    successful_fits = []
    for fittag in ['b','s']:
        if 'fit_'+fittag not in [k.GetName() for k in tfile.GetListOfKeys()]:
            warnings.warn('Unable to find result fit_%s...'%fittag,RuntimeWarning)
        else:
            successful_fits.append(fittag)
    return successful_fits

def nuis_pulls(vtol=0.3, stol=0.1, vtol2=2.0, stol2=0.5, regex='^(?!.*(_bin_|_par))'):
    diffnuis_cmd = f"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnosticsTest.root --abs -g nuisance_pulls.root --vtol={vtol} --stol={stol} --vtol2={vtol2} --stol2={stol2} --all --regex='{regex}'"
    execute_cmd(diffnuis_cmd)
    # Make a PDF of the nuisance_pulls.root
    if os.path.exists('nuisance_pulls.root'):
        nuis_file = ROOT.TFile.Open('nuisance_pulls.root')
        nuis_can = nuis_file.Get('nuisances')
        nuis_can.Print('nuisance_pulls.pdf','pdf')
        nuis_file.Close()

def save_post_fit_parametric_vals():
    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    goodFitTags = _get_good_fit_results(fit_result_file)
    for fittag in goodFitTags:
        coeffs_final = fit_result_file.Get('fit_'+fittag).floatParsFinal()
        all_par_names = [] # get names of anything matching *_par<N>
        for i in range(coeffs_final.getSize()):
            name = coeffs_final.at(i).GetName()
            if name.split('_')[-1].startswith('par'):
                all_par_names.append(name)
        all_par_names.sort()
        # Get unique prefixes to _par<N>
        all_obj_names = set(['_'.join(name.split('_')[:-1]) for name in all_par_names])

        for obj_name in all_obj_names:
            with open('rpf_params_%s_fit%s.txt'%(obj_name,fittag),'w') as param_out:
                for par_name in [p for p in all_par_names if p.startswith(obj_name)]:
                    var = coeffs_final.find(par_name)
                    param_out.write('%s: %s +/- %s\n'%(par_name, var.getValV(), var.getError()))

    fit_result_file.Close()

def gen_post_fit_shapes():
    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    goodFitTags = _get_good_fit_results(fit_result_file)
    for t in goodFitTags:
        if os.path.exists('higgsCombineTest.FitDiagnostics.mH120.root'):
            workspace_file = 'higgsCombineTest.FitDiagnostics.mH120.root'
        else:
            workspace_file = 'higgsCombineTest.FitDiagnostics.mH120.123456.root'
        shapes_cmd = 'PostFit2DShapesFromWorkspace -w {w} --output postfitshapes_{t}.root -f fitDiagnosticsTest.root:fit_{t} --postfit --samples 100 --print > PostFitShapes2D_stderr_{t}.txt'.format(t=t,w=workspace_file)
        execute_cmd(shapes_cmd)
    fit_result_file.Close()

def _reduced_corr_matrix(fit_result, varsToIgnore=[], varsOfInterest=[], threshold=0):
    if threshold < 0:
        raise ValueError('Threshold for correlation matrix values to plot must be a positive number.')

    ROOT.gStyle.SetOptStat(0)
    # ROOT.gStyle.SetPaintTextFormat('.3f')
    CM = fit_result.correlationMatrix()
    finalPars = fit_result.floatParsFinal()

    nParams = CM.GetNcols()
    finalParamsDict = {}
    for cm_index in range(nParams):
        if varsOfInterest == []:
            if finalPars.at(cm_index).GetName() not in varsToIgnore:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index
        else:
            if finalPars.at(cm_index).GetName() in varsOfInterest:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index

    nFinalParams = len(finalParamsDict.keys())
    out = ROOT.TH2D('correlation_matrix','correlation_matrix',nFinalParams,0,nFinalParams,nFinalParams,0,nFinalParams)
    out_txt = ''

    for out_x_index, paramXName in enumerate(sorted(finalParamsDict.keys())):
        cm_index_x = finalParamsDict[paramXName]

        if not any([abs(CM[cm_index_x][cm_index_y]) > threshold for cm_index_y in range(nParams) if cm_index_y != cm_index_x]):
            continue

        for out_y_index, paramYName in enumerate(sorted(finalParamsDict.keys())):
            cm_index_y = finalParamsDict[paramYName]
            if cm_index_x > cm_index_y:
                out_txt += '%s:%s = %s\n'%(paramXName,paramYName,CM[cm_index_x][cm_index_y])
            out.Fill(out_x_index+0.5,out_y_index+0.5,CM[cm_index_x][cm_index_y])

        out.GetXaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
        out.GetYaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
    out.SetMinimum(-1)
    out.SetMaximum(+1)

    return out, out_txt

def plot_correlation_matrix(varsToIgnore, threshold=0, corrText=False):
    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    for fittag in _get_good_fit_results(fit_result_file):
        fit_result = fit_result_file.Get("fit_"+fittag)
        if hasattr(fit_result,'correlationMatrix'):
            corrMtrx, corrTxt = _reduced_corr_matrix(fit_result, varsToIgnore=varsToIgnore, threshold=threshold)
            corrMtrxCan = ROOT.TCanvas('c','c',1400,1000)
            corrMtrxCan.cd()
            corrMtrxCan.SetBottomMargin(0.22)
            corrMtrxCan.SetLeftMargin(0.17)
            corrMtrxCan.SetTopMargin(0.06)

            corrMtrx.GetXaxis().SetLabelSize(0.01)
            corrMtrx.GetYaxis().SetLabelSize(0.01)
            corrMtrx.Draw('colz text' if corrText else 'colz')
            corrMtrxCan.Print('plots_fit_%s/correlation_matrix.png'%fittag,'png')
            corrMtrxCan.Print('plots_fit_%s/correlation_matrix.pdf'%fittag,'pdf')

            with open('plots_fit_%s/correlation_matrix.txt'%fittag,'w') as corrTxtFile:
                corrTxtFile.write(corrTxt)

        else:
            warnings.warn('Not able to produce correlation matrix.',RuntimeWarning)

    fit_result_file.Close()

def plot_gof(tag, subtag, seed=123456, condor=False):
    with cd(tag+'/'+subtag):
        if condor:
            tmpdir = 'notneeded/tmp/'
            execute_cmd('mkdir '+tmpdir) 
            execute_cmd('cat %s_%s_gof_toys_output_*.tgz | tar zxvf - -i --strip-components 2 -C %s'%(tag,subtag,tmpdir))
            toy_limit_tree = ROOT.TChain('limit')
            if len(glob.glob(tmpdir+'higgsCombine_gof_toys.GoodnessOfFit.mH120.*.root')) == 0:
                raise Exception('No files found')
            toy_limit_tree.Add(tmpdir+'higgsCombine_gof_toys.GoodnessOfFit.mH120.*.root') 
            
        else:
            toyOutput = ROOT.TFile.Open('higgsCombine_gof_toys.GoodnessOfFit.mH120.{seed}.root'.format(seed=seed))
            toy_limit_tree = toyOutput.Get('limit')

        # Now to analyze the output
        # Get observation
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetOptStat(False)
        gof_data_file = ROOT.TFile.Open('higgsCombine_gof_data.GoodnessOfFit.mH120.root')
        gof_limit_tree = gof_data_file.Get('limit')
        gof_limit_tree.GetEntry(0)
        gof_data = gof_limit_tree.limit

        # Get toys
        toy_limit_tree.Draw('limit>>hlimit','limit>1.0 && limit<%s && limit != %s'%(gof_data*2.0,gof_data)) 
        htoy_gof = ROOT.gDirectory.Get('hlimit')
        time.sleep(1) # if you don't sleep the code moves too fast and won't perform the fit
        htoy_gof.Fit("gaus")

        # Fit toys and derive p-value
        gaus = htoy_gof.GetFunction("gaus")
        pvalue = 1-(1/gaus.Integral(-float("inf"),float("inf")))*gaus.Integral(-float("inf"),gof_data)

        # Write out for reference
        with open('gof_results.txt','w') as out:
            out.write('Test statistic in data = '+str(gof_data))
            out.write('Mean from toys = '+str(gaus.GetParameter(1)))
            out.write('Width from toys = '+str(gaus.GetParameter(2)))
            out.write('p-value = '+str(pvalue))

        # Extend the axis if needed
        if htoy_gof.GetXaxis().GetXmax() < gof_data:
            print ('Axis limit greater than GOF p value')
            binwidth = htoy_gof.GetXaxis().GetBinWidth(1)
            xmin = htoy_gof.GetXaxis().GetXmin()
            new_xmax = int(gof_data*1.1)
            new_nbins = int((new_xmax-xmin)/binwidth)
            toy_limit_tree.Draw('limit>>hlimitrebin('+str(new_nbins)+', '+str(xmin)+', '+str(new_xmax)+')','limit>0.001 && limit<1500') 
            htoy_gof = ROOT.gDirectory.Get('hlimitrebin')
            htoy_gof.Fit("gaus")
            gaus = htoy_gof.GetFunction("gaus")

        # Arrow for observed
        arrow = ROOT.TArrow(gof_data,0.25*htoy_gof.GetMaximum(),gof_data,0)
        arrow.SetLineWidth(2)

        # Legend
        leg = ROOT.TLegend(0.1,0.7,0.4,0.9)
        leg.SetLineColor(ROOT.kWhite)
        leg.SetLineWidth(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.AddEntry(htoy_gof,"toy data","lep")
        leg.AddEntry(arrow,"observed = %.1f"%gof_data,"l")
        leg.AddEntry(0,"p-value = %.2f"%(pvalue),"")

        # Draw
        cout = ROOT.TCanvas('cout','cout',800,700)
        htoy_gof.SetTitle('')
        htoy_gof.Draw('pez')
        arrow.Draw()
        leg.Draw()

        cout.Print('gof_plot.pdf','pdf')
        cout.Print('gof_plot.png','png')

    if condor:
            execute_cmd('rm -r '+tmpdir)

def plot_signalInjection(tag, subtag, injectedAmount, seed=123456, stats=True, condor=False):
    # if injectedAmount is not an integer, need to look for different file
    # see: https://github.com/lcorcodilos/2DAlphabet/blob/e089ed1da63172770726b3e6f406c11e611e057d/TwoDAlphabet/twoDalphabet.py#L488
    injectedName = str(injectedAmount).replace('.','p')
    with cd(tag+'/'+subtag):
        if condor:
            tmpdir = 'notneeded/tmp/'
            execute_cmd('mkdir '+tmpdir) 
            execute_cmd('cat %s_%s_sigInj_r%s_output_*.tgz | tar zxvf - -i --strip-components 2 -C %s'%(tag,subtag,injectedName,tmpdir))
            tree_fit_sb = ROOT.TChain('tree_fit_sb')
            tree_fit_sb.Add(tmpdir+'fitDiagnostics_sigInj_r{rinj}*.root'.format(rinj=injectedName)) 
            
        else:
            toyOutput = ROOT.TFile.Open('fitDiagnostics_sigInj_r{rinj}_{seed}.root'.format(rinj=injectedName,seed=seed))
            tree_fit_sb = toyOutput.Get('tree_fit_sb')

        ROOT.gROOT.SetBatch(True)
        if stats:
            ROOT.gStyle.SetOptStat(True)
        else:
            ROOT.gStyle.SetOptStat(False)
        # Final plotting
        result_can = ROOT.TCanvas('sigpull_can','sigpull_can',800,700)

        tree_fit_sb.Draw("(r-{rinj})/(rHiErr*(r<{rinj})+rLoErr*(r>{rinj}))>>sigpull(20,-5,5)".format(rinj=injectedAmount),"fit_status>=0")
        hsigpull = ROOT.gDirectory.Get('sigpull')
        tree_fit_sb.Draw("(r-{rinj})>>sigstrength(20,-1,1)".format(rinj=injectedAmount),"fit_status>=0")
        hsignstrength = ROOT.gDirectory.Get('sigstrength')

        hsigpull.Fit("gaus","L")
        hsigpull.SetTitle('')
        hsigpull.GetXaxis().SetTitle('(r-%s)/rErr'%injectedAmount)
        result_can.cd()
        hsigpull.Draw('pe')
        result_can.Print('signalInjection_r%s_pull.png'%(str(injectedAmount).replace('.','p')),'png')

        hsignstrength.Fit("gaus","L")
        hsignstrength.SetTitle('')
        hsignstrength.GetXaxis().SetTitle('r-%s'%injectedAmount)
        result_can.cd()
        hsignstrength.Draw('pe')
        result_can.Print('signalInjection_r%s.png'%(str(injectedAmount).replace('.','p')),'png')

    if condor:
            execute_cmd('rm -r '+tmpdir)
