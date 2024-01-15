import glob
import ROOT, os, warnings, pandas, math, time
from PIL import Image
from TwoDAlphabet.helpers import set_hist_maximums, execute_cmd, cd
from TwoDAlphabet.binning import stitch_hists_in_x, convert_to_events_per_unit, get_min_bin_width
from TwoDAlphabet.ext import tdrstyle, CMS_lumi

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

        if not loadExisting:
            self._make()
        else:
            self._load()

    def __del__(self):
        '''On deletion, save DataFrame to csv and close ROOT file.'''
        self.df.to_csv('%s/df.csv'%self.dir)
        self.root_out.Close()

    def _load(self):
        '''Open pickled DataFrame and output ROOT file
        and reference with `self.df` and `self.root_out` attributes.'''
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name)
        self.df = pandas.read_csv('%s/df.csv'%self.dir)

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

        color = int(color) #ROOT call in C++ sometimes cannot convert it to int

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
                    color = ROOT.kBlack
                    proc_type = 'TOTAL'
                    proc_title = 'TotalBkg'

                self.df = self.df.append({'process':process,
                                          'region':region,
                                          'process_type': proc_type,
                                          'title': proc_title}, ignore_index=True)
                
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

    def plot_projections(self, prefit=False):
        '''Plot comparisons of data and the post-fit background model and signal
        using the 1D projections. Canvases are grouped based on projection axis.
        The canvas rows are separate selection regions while the columns 
        are the different slices of the un-plotted axis.

        Args:
        prefit (bool): If True, will plot the prefit distributions instead of postfit. Defaults to False.
        Returns:
            None
        '''
        pads = pandas.DataFrame()
        for region, group in self.df.groupby('region'):
            binning,_ = self.twoD.GetBinningFor(region)

            for logyFlag in [False, True]:
                ordered_bkgs = self._order_df_on_proc_list(
                                            group[group.process_type.eq('BKG')], proc_type='BKG',
                                            alphaBottom=(not logyFlag))
                signals = group[group.process_type.eq('SIGNAL')]

                for proj in ['prefit_projx','prefit_projy'] if prefit else ['postfit_projx','postfit_projy']:
                    for islice in range(3):
                        projn = proj+str(islice)
                        sig_projn = projn
                        if self.twoD.options.plotPrefitSigInFitB and self.fittag == 'b':
                            sig_projn = projn.replace('postfit','prefit')

                        this_data =      self.Get(row=group.loc[group.process_type.eq('DATA')].iloc[0], hist_type=projn)
                        this_totalbkg =  self.Get(row=group.loc[group.process_type.eq('TOTAL')].iloc[0], hist_type=projn)
                        these_bkgs =    [self.Get(row=ordered_bkgs.iloc[irow], hist_type=projn) for irow in range(ordered_bkgs.shape[0])]
                        these_signals = [self.Get(row=signals.iloc[irow], hist_type=sig_projn) for irow in range(signals.shape[0])]

                        slice_edges = (
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice],
                            binning.xtitle if 'y' in proj else binning.ytitle,
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice+1],
                            'GeV'
                        )
                        slice_str = '%s < %s < %s %s'%slice_edges

                        out_pad_name = '{d}/base_figs/{projn}_{reg}{logy}'.format(
                                            d=self.dir, projn=projn, reg=region,
                                            logy='' if logyFlag == False else '_logy')
                        
                        make_pad_1D(out_pad_name, data=this_data, bkgs=these_bkgs, signals=these_signals,
                                    subtitle=slice_str, totalBkg=this_totalbkg,
                                    logyFlag=logyFlag, year=self.twoD.options.year, preVsPost=False,
                                    extraText='', savePDF=True, savePNG=True, ROOTout=False)
                        pads = pads.append({'pad':out_pad_name+'.png', 'region':region, 'proj':projn, 'logy':logyFlag}, ignore_index=True)

        for logy in ['','_logy']:
            for proj in ['prefit_projx','prefit_projy'] if prefit else ['postfit_projx','postfit_projy']:
                these_pads = pads.loc[pads.proj.str.contains(proj)]
                if logy == '':
                    these_pads = these_pads.loc[these_pads.logy.eq(False)]
                else:
                    these_pads = these_pads.loc[these_pads.logy.eq(True)]

                these_pads = these_pads.sort_values(by=['region','proj'])['pad'].to_list()
                out_can_name = '{d}/{proj}{logy}'.format(d=self.dir, proj=proj, logy=logy)
                make_can(out_can_name, these_pads)
        
    def plot_pre_vs_post(self):
        '''Make comparisons for each background process of pre and post fit projections.
        '''
        for proj in ['projx','projy']:
            pads = pandas.DataFrame()
            for pr, _ in self.df[~self.df.process.isin(['data_obs','TotalBkg'])].groupby(['process','region']):
                process, region = pr[0], pr[1]
                binning,_ = self.twoD.GetBinningFor(region)
                for islice in range(3):
                    projn = proj+str(islice)
                    post = self.Get('%s_%s_postfit_%s'%(process,region,projn))
                    post.SetLineColor(ROOT.kBlack)
                    post.SetTitle('          Postfit,'+process) # spaces are for legend aesthetics

                    pre = self.Get('%s_%s_prefit_%s'%(process,region,projn))
                    pre.SetLineColor(ROOT.kRed)
                    pre.SetTitle('Prefit, '+process)

                    slice_edges = (
                        self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice],
                        binning.xtitle if 'y' in proj else binning.ytitle,
                        self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice+1],
                        'GeV'
                    )
                    slice_str = '%s < %s < %s %s'%slice_edges

                    out_pad_name = '{d}/base_figs/{p}_{reg}_{projn}'.format(d=self.dir,p=process,projn=projn, reg=region)
                    make_pad_1D(
                        out_pad_name, 
                        post, [pre], totalBkg=pre, subtitle=slice_str, savePDF=True, savePNG=True, 
                        datastyle='histe', year=self.twoD.options.year, extraText='',
            preVsPost=True  # This tells make_pad_1D() that we're not passing in data distributions but rather a non-data postfit dist and to relabel the legend
                    )
                    
                    pads = pads.append({'pad':out_pad_name+'.png','process':process,'region':region,'proj':projn}, ignore_index=True)

        for process, padgroup in pads.groupby('process'):
            these_pads = padgroup.sort_values(by=['region','proj'])['pad'].to_list()
        make_can('{d}/{p}_{proj}'.format(d=self.dir, p=process,proj=proj), these_pads)


    def plot_transfer_funcs(self):
        raise NotImplementedError()
        # # Need to sample the space to get the Rp/f with proper errors (1000 samples)
        # rpf_xnbins = len(self.fullXbins)-1
        # rpf_ynbins = len(self.newYbins)-1
        # if self.rpfRatio == False: rpf_zbins = [i/1000000. for i in range(0,1000001)]
        # else: rpf_zbins = [i/1000. for i in range(0,5001)]
        # rpf_samples = TH3F('rpf_samples','rpf_samples',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins), len(rpf_zbins)-1, array.array('d',rpf_zbins))# TH3 to store samples
        # sample_size = 500

        # # Collect all final parameter values
        # param_final = fit_result.floatParsFinal()
        # coeffs_final = RooArgSet()
        # for v in self.rpf.funcVars.keys():
        #     coeffs_final.add(param_final.find(v))

        # # Now sample to generate the Rpf distribution
        # for i in range(sample_size):
        #     sys.stdout.write('\rSampling '+str(100*float(i)/float(sample_size)) + '%')
        #     sys.stdout.flush()
        #     param_sample = fit_result.randomizePars()

        #     # Set params of the Rpf object
        #     coeffIter_sample = param_sample.createIterator()
        #     coeff_sample = coeffIter_sample.Next()
        #     while coeff_sample:
        #         # Set the rpf parameter to the sample value
        #         if coeff_sample.GetName() in self.rpf.funcVars.keys():
        #             self.rpf.setFuncParam(coeff_sample.GetName(), coeff_sample.getValV())
        #         coeff_sample = coeffIter_sample.Next()

        #     # Loop over bins and fill
        #     for xbin in range(1,rpf_xnbins+1):
        #         for ybin in range(1,rpf_ynbins+1):
        #             bin_val = 0

        #             thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
        #             thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)

        #             if self.recycleAll:
        #                 # Remap to [-1,1]
        #                 x_center_mapped = (thisXCenter - self.newXbins['LOW'][0])/(self.newXbins['HIGH'][-1] - self.newXbins['LOW'][0])
        #                 y_center_mapped = (thisYCenter - self.newYbins[0])/(self.newYbins[-1] - self.newYbins[0])

        #                 # And assign it to a RooConstVar 
        #                 x_const = RooConstVar("ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,x_center_mapped)
        #                 y_const = RooConstVar("ConstVar_y_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,y_center_mapped)
                        
        #                 # Now get the Rpf function value for this bin 
        #                 self.allVars.append(x_const)
        #                 self.allVars.append(y_const)
        #                 self.rpf.evalRpf(x_const, y_const,xbin,ybin)

        #             # Determine the category
        #             if thisXCenter > self.newXbins['LOW'][0] and thisXCenter < self.newXbins['LOW'][-1]: # in the LOW category
        #                 thisxcat = 'LOW'
        #             elif thisXCenter > self.newXbins['SIG'][0] and thisXCenter < self.newXbins['SIG'][-1]: # in the SIG category
        #                 thisxcat = 'SIG'
        #             elif thisXCenter > self.newXbins['HIGH'][0] and thisXCenter < self.newXbins['HIGH'][-1]: # in the HIGH category
        #                 thisxcat = 'HIGH'

        #             bin_val = self.rpf.getFuncBinVal(thisxcat,xbin,ybin)

        #             rpf_samples.Fill(thisXCenter,thisYCenter,bin_val)

        # print ('\n')
        # rpf_final = TH2F('rpf_final','rpf_final',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins))
        # # Now loop over all x,y bin in rpf_samples, project onto Z axis, 
        # # get the mean and RMS and set as the bin content and error in rpf_final
        # for xbin in range(1,rpf_final.GetNbinsX()+1):
        #     for ybin in range(1,rpf_final.GetNbinsY()+1):
        #         temp_projz = rpf_samples.ProjectionZ('temp_projz',xbin,xbin,ybin,ybin)
        #         rpf_final.SetBinContent(xbin,ybin,temp_projz.GetMean())
        #         rpf_final.SetBinError(xbin,ybin,temp_projz.GetRMS())

        # rpf_final.SetTitle('')
        # rpf_final.GetXaxis().SetTitle(self.xVarTitle)
        # rpf_final.GetYaxis().SetTitle(self.yVarTitle)
        # rpf_final.GetZaxis().SetTitle('R_{P/F}' if self.rpfRatio == False else 'R_{Ratio}')
        # rpf_final.GetXaxis().SetTitleSize(0.045)
        # rpf_final.GetYaxis().SetTitleSize(0.045)
        # rpf_final.GetZaxis().SetTitleSize(0.045)
        # rpf_final.GetXaxis().SetTitleOffset(1.2)
        # rpf_final.GetYaxis().SetTitleOffset(1.5)
        # rpf_final.GetZaxis().SetTitleOffset(1.3)

        # rpf_c = TCanvas('rpf_c','Post-fit R_{P/F}',1000,700)
        # CMS_lumi.lumiTextSize = 0.75
        # CMS_lumi.cmsTextSize = 0.85
        # CMS_lumi.extraText = 'Preliminary'
        # CMS_lumi.CMS_lumi(rpf_c, self.year, 11)
        # rpf_c.SetRightMargin(0.2)
        # rpf_final.Draw('colz')
        # rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_colz.pdf','pdf')
        # rpf_final.Draw('surf')
        # rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_surf.pdf','pdf')
        # rpf_final.Draw('pe')
        # rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_errs.pdf','pdf')

        # rpf_file = TFile.Open(self.projPath+'/plots/postfit_rpf_fit'+fittag+'.root','RECREATE')
        # rpf_file.cd()
        # rpf_final.Write()
        # rpf_file.Close()

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

def make_pad_1D(outname, data, bkgs=[], signals=[], title='', subtitle='',
            totalBkg=None, logyFlag=False, ROOTout=None, savePDF=False, savePNG=False,
            dataOff=False, preVsPost=False, datastyle='pe X0', year=1, addSignals=True, extraText='Preliminary'):
    '''Make a pad holding a 1D plot with standardized formatting conventions.

    Args:
        outname (str): Output file path name.
        data (TH1): Data histogram.
        bkgs ([TH1]): List of background histograms (will be stacked).
        signals ([TH1]): List of signal histograms.
        title (str, optional): Title of plot. Only applicable if bkgs is empty. Defaults to ''.
        subtitle (str, optional): Subtitle text for physics information (like slice ranges). Defaults to ''.
        totalBkg (TH1, optional): Total background estimate from fit. Used to get total background uncertianty. Defaults to None.
        logyFlag (bool, optional): Make log y-axis. Defaults to False.
        saveROOT (bool, optional): Save to master ROOT file. Defaults to False.
        savePDF (bool, optional): Save to PDF. Defaults to True.
        savePNG (bool, optional): Save to PNG. Defaults to True.
        dataOff (bool, optional): Turn off the data from plotting. Defaults to False.
    preVsPost (bool, optional): Incoming data histogram is postfit distribution of non-data process. If True, renames legend entry. See plot_pre_vs_post().
        datastyle (str, optional): ROOT drawing style for the data. Defaults to 'pe X0'.
        year (int, optional): Luminosity formatting. Options are 16, 17, 18, 1 (full Run 2), 2 (16+17+18). Defaults to 1.
        addSignals (bool, optional): If True, multiple signals will be added together and plotted as one. If False, signals are plotted individually. Defaults to True.
        extraText (str, optional): Prepended to the CMS subtext. Defaults to 'Preliminary'.

    Returns:
        ROOT.TPad: Output pad.
    '''

    def _draw_extralumi_tex():
        lumiE = ROOT.TLatex()
        lumiE.SetNDC()
        lumiE.SetTextAngle(0)
        lumiE.SetTextColor(ROOT.kBlack)
        lumiE.SetTextFont(42)
        lumiE.SetTextAlign(31) 
        lumiE.SetTextSize(0.7*0.1)
        lumiE.DrawLatex(1-0.05,1-0.1+0.2*0.1,"137 fb^{-1} (13 TeV)")

    pad = _make_pad_gen(outname)

    data.SetBinErrorOption(ROOT.TH1.kPoisson)
    data.SetLineColorAlpha(ROOT.kBlack, 0 if dataOff else 1)
    if 'pe' in datastyle.lower():
        data.SetMarkerColorAlpha(ROOT.kBlack,0 if dataOff else 1)
        data.SetMarkerStyle(8)
    if 'hist' in datastyle.lower():
        data.SetFillColorAlpha(0,0)

    data.SetTitleOffset(1.15,"xy")
    data.GetYaxis().SetTitleOffset(1.04)
    data.GetYaxis().SetLabelSize(0.07)
    data.GetYaxis().SetTitleSize(0.09)
    data.GetXaxis().SetLabelSize(0.07)
    data.GetXaxis().SetTitleSize(0.09)
    data.GetXaxis().SetLabelOffset(0.05)
    data.GetYaxis().SetNdivisions(508)
    data.SetMaximum(1.35*data.GetMaximum())

    if len(bkgs) == 0:
        data.SetTitle(title)
        data.SetTitleOffset(1.1)
        data.Draw(datastyle)
        CMS_lumi.CMS_lumi(pad, year, 11)
    else:
        data.SetTitle('')
        if not dataOff:
            main_pad = ROOT.TPad(data.GetName()+'_main',data.GetName()+'_main',0, 0.35, 1, 1)
            sub_pad  = ROOT.TPad(data.GetName()+'_sub',data.GetName()+'_sub',0, 0, 1, 0.35)
        else:
            main_pad = ROOT.TPad(data.GetName()+'_main',data.GetName()+'_main',0, 0.1, 1, 1)
            sub_pad  = ROOT.TPad(data.GetName()+'_sub',data.GetName()+'_sub',0, 0, 0, 0)

        main_pad.SetBottomMargin(0.04)
        main_pad.SetLeftMargin(0.17)
        main_pad.SetRightMargin(0.05)
        main_pad.SetTopMargin(0.1)

        sub_pad.SetLeftMargin(0.17)
        sub_pad.SetRightMargin(0.05)
        sub_pad.SetTopMargin(0)
        sub_pad.SetBottomMargin(0.35)           
        
        pad.cd()
        main_pad.Draw()
        sub_pad.Draw()

        if len(signals) == 0: nsignals = 0
        elif addSignals:      nsignals = 1
        else:                 nsignals = len(signals)

        legend_bottomY = 0.73-0.03*(min(len(bkgs),6)+nsignals+1)
        legend = ROOT.TLegend(0.65,legend_bottomY,0.90,0.88)
        legend.SetBorderSize(0)
        if not dataOff:
            if preVsPost:
                legend.AddEntry(data,'Postfit',datastyle)
            else:
                legend.AddEntry(data,'Data',datastyle)

        totalBkg.SetMarkerStyle(0)
        totalBkg_err = totalBkg.Clone()
        totalBkg.SetLineColor(ROOT.kBlack)
        totalBkg_err.SetLineColor(ROOT.kBlack)
        totalBkg_err.SetLineWidth(0)
        totalBkg_err.SetFillColor(ROOT.kBlack)
        totalBkg_err.SetFillStyle(3354)
    if preVsPost:
        # Determine whether we're plotting uncertainty on Prefit or Postfit distribution (plot_pre_vs_post() uses prefit unc)
        isPrefit = 1 if 'Prefit' in totalBkg.GetTitle() else 0
        legend.AddEntry(totalBkg_err,'Total {} unc.'.format('prefit' if isPrefit else 'postfit'), 'F')
    else:
        legend.AddEntry(totalBkg_err,'Total bkg unc.','F')

        sigs_to_plot = signals
        # Can add together for total signal
        if addSignals and len(signals) > 0:
            totsig = signals[0].Clone()
            for isig in range(1,len(signals)):
                totsig.Add(signals[isig])
            sigs_to_plot = [totsig]
        # Plot either way
        for isig,sig in enumerate(sigs_to_plot):
            sig.SetLineWidth(2)
            legend.AddEntry(sig,sig.GetTitle().split(',')[0],'L')

        stack = ROOT.THStack(outname.split('/')[-1]+'_stack',outname.split('/')[-1]+'_stack')
        # Build the stack
        legend_info = []
        for bkg in bkgs:     # Won't loop if bkglist is empty
            if logyFlag:
                bkg.SetMinimum(1e-3)
            stack.Add(bkg)
            legend_info.append((bkg.GetTitle().split(',')[0], bkg))

        # Deal with legend which needs ordering reversed from stack build
        legend_duplicates = []
        for bname, background in reversed(legend_info):
            if bname not in legend_duplicates:
                legend.AddEntry(background, bname, 'f')
                legend_duplicates.append(bname)
        # set_hist_maximums([stack, totalBkg, data], 2.5-legend_topY+0.03)
        # Go to main pad and draw
        main_pad.cd()
        if logyFlag:
            main_pad.SetLogy()
            data.SetMaximum(50*data.GetMaximum())
            data.SetMinimum(1e-3)
            totalBkg.SetMinimum(1e-3)
            totalBkg_err.SetMinimum(1e-3)
            stack.SetMinimum(1e-3)
            for sig in signals: sig.SetMinimum(1e-3)
            
            # main_pad.RedrawAxis()
        data.Draw(datastyle)
        stack.Draw('hist same') # need to draw twice because the axis doesn't exist for modification until drawing
        try:    stack.GetYaxis().SetNdivisions(508)
        except: stack.GetYaxis().SetNdivisions(8,5,0)
        stack.Draw('hist same')
        for sig in signals:
            sig.Draw('hist same')
        
        # Draw total hist and error
        totalBkg.Draw('hist same')
        totalBkg_err.Draw('e2 same')
        data.Draw(datastyle+' same')
        legend.Draw()
        
        CMS_lumi.extraText = extraText
        CMS_lumi.cmsTextSize = 0.9
        CMS_lumi.cmsTextOffset = 2
        CMS_lumi.lumiTextSize = 0.9
        CMS_lumi.CMS_lumi(main_pad, year, 11)
        
        subtitle_tex = ROOT.TLatex()
        subtitle_tex.SetNDC()
        subtitle_tex.SetTextAngle(0)
        subtitle_tex.SetTextColor(ROOT.kBlack)
        subtitle_tex.SetTextFont(42)
        subtitle_tex.SetTextAlign(12) 
        subtitle_tex.SetTextSize(0.06)
        subtitle_tex.DrawLatex(0.208,0.74,subtitle)

        sub_pad.cd()
        pull = _make_pull_plot(data, totalBkg, preVsPost) # If plotting pre vs postfit dists, ensure the Y axis of pull plot reflects this
        pull.Draw('hist')
        pad.cd()
    
    ROOT.SetOwnership(pad, False)
    _save_pad_generic(outname, pad, ROOTout, savePDF, savePNG)
    return pad

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
            raise RuntimeError('histlist of size %s not currently supported: %s'%(len(padnames),[p.GetName() for p in padnames]))

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

def gen_projections(ledger, twoD, fittag, loadExisting=False, prefit=False):
    '''
    Optional Args:
    loadExisting (bool): Flag to load existing projections instead of remaking everything. Defaults to False.
    prefit       (bool): Flag to plot prefit distributions instead of postfit. Defaults to False.
    '''
    plotter = Plotter(ledger, twoD, fittag, loadExisting)
    plotter.plot_2D_distributions()
    plotter.plot_projections(prefit)
    plotter.plot_pre_vs_post()
    # plotter.plot_transfer_funcs()

def make_systematic_plots(twoD):
    '''Make plots of the systematic shape variations of each process based on those
    processes and systematic shapes specified in the config. Shapes are presented 
    as projections onto 1D axis where no selection has been made on the axis not
    being plotted. Plots are saved to UncertPlots/.
    '''
    c = ROOT.TCanvas('c','c',800,700)

    for (p,r), _ in twoD.df.groupby(['process','region']):
        if p == 'data_obs': continue

        nominal_full = twoD.organizedHists.Get(process=p,region=r,systematic='')
        binning, _ = twoD.GetBinningFor(r)

        for axis in ['X','Y']:
            nominal = getattr(nominal_full,'Projection'+axis)('%s_%s_%s_%s'%(p,r,'nom','proj'+axis))
            for s in twoD.ledger.GetShapeSystematics(drop_norms=True):
                up = getattr(twoD.organizedHists.Get(process=p,region=r,systematic=s+'Up'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Up','proj'+axis))
                down = getattr(twoD.organizedHists.Get(process=p,region=r,systematic=s+'Down'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Down','proj'+axis))

                c.cd()
                nominal.SetLineColor(ROOT.kBlack)
                nominal.SetFillColor(ROOT.kYellow-9)
                up.SetLineColor(ROOT.kRed)
                up.SetFillColorAlpha(ROOT.kWhite, 0)
                down.SetLineColor(ROOT.kBlue)
                down.SetFillColorAlpha(ROOT.kWhite, 0)

                up.SetLineStyle(9)
                down.SetLineStyle(9)
                up.SetLineWidth(2)
                down.SetLineWidth(2)

                nominal,up,down = set_hist_maximums([nominal,up,down])
                nominal.SetXTitle(getattr(binning,axis.lower()+'title'))

                nominal.SetTitle('')
                nominal.GetXaxis().SetTitleOffset(1.0)
                nominal.GetXaxis().SetTitleSize(0.05)
                c.SetRightMargin(0.16)

                nominal.Draw('hist')
                up.Draw('same hist')
                down.Draw('same hist')

                c.Print(twoD.tag+'/UncertPlots/Uncertainty_%s_%s_%s_%s.png'%(p,r,s,'proj'+axis),'png')

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

def nuis_pulls():
    diffnuis_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnosticsTest.root --abs -g nuisance_pulls.root'
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
