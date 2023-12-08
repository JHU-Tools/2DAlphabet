import ROOT, array
from math import sqrt

class Binning:
    '''Class to handle information on and manipulations of binning schemes.'''
    def __init__(self, name, binning_dict, start_template):
        '''Initialize Binning object.

        Args:
            name (str): Name to create unique objects owned by Binning.
            binning_dict (dict): "BINNING" section of configuration json which specifies the X and Y binning schemes.
            start_template (TH2): Histogram to compare against when doing sanity checks.
        '''
        self.name = name
        self.sigStart = binning_dict['X']['SIGSTART']
        self.sigEnd = binning_dict['X']['SIGEND']
        self.xtitle = binning_dict['X']['TITLE']
        self.ytitle = binning_dict['Y']['TITLE']
        self.xbinByCat, self.ybinList = parse_binning_info(binning_dict)
        self.ySlices,self.ySliceIdx = self._getYslices(binning_dict) # x slices defined as properties
        self._checkBinning('X',start_template)
        self._checkBinning('Y',start_template)
        self.xVars, self.yVar = self.CreateRRVs(binning_dict['X'], binning_dict['Y']) 

    def CreateRRVs(self,xdict,ydict):
        '''Create the RooRealVars representing the X and Y axes.
        For the X axis, three RooRealVars are returned in a dictionary with
        each key mapping the RooRealVar to either LOW, SIG, or HIGH categories.

        Args:
            xdict (dict): "X" subsection of the "BINNING" section of the config.
            ydict (dict): "Y" subsection of the "BINNING" section of the config.

        Returns:
            tuple: (0) dict of X axis RooRealVars and (1) Y axis RooRealVar.
        '''
        yRRV = create_RRV_base(ydict['NAME']+'_'+self.name,
                          ydict['TITLE'],
                          self.ybinList)
        xRRVs = {}
        for c in ['LOW','SIG','HIGH']:
            xRRVs[c] = create_RRV_base(xdict['NAME']+'_'+c+'_'+self.name,
                                  xdict['TITLE'],
                                  self.xbinByCat[c])
        return xRRVs,yRRV

    def _checkBinning(self,axis,start_template):
        '''Perform sanity check that new binning scheme is a subset of the
        starting input space.

        Args:
            axis (str): Either "X" or "Y".

        Raises:
            ValueError: If requested binning is not a subset of the input space.
            ValueError: If binning sheme is not in increasing order.
        '''
        input_min = getattr(start_template,'Get%saxis'%axis)().GetXmin()
        input_max = getattr(start_template,'Get%saxis'%axis)().GetXmax()

        if axis == 'X': new_bins = self.xbinList
        else: new_bins = self.ybinList

        if (new_bins[0] < input_min) or (new_bins[-1] > input_max):
            raise ValueError('%s axis requested is larger than input\n\tInput: [%s,%s]\n\tRequested: [%s,%s]'%(axis,input_min,input_max,new_bins[0],new_bins[-1]))
        prev = -1000000
        for b in new_bins:
            if b > prev:
                prev = b
            else:
                raise ValueError('%s axis bin edges must be in increasing order!'%axis)

    def _getYslices(self,binning_dict):
        if 'SLICES' in binning_dict['Y']:
            if len(binning_dict['Y']['SLICES']) != 4:
                raise RuntimeError('Must define Y SLICES as a list of four values which represent the edges of the continuous slices.')
            elif binning_dict['Y']['SLICES'][0] != self.ybinList[0]:
                raise ValueError('First edge of Y SLICES does not match axis (%s vs %s)'%(binning_dict['Y']['SLICES'][0], self.ybinList[0]))
            elif binning_dict['Y']['SLICES'][-1] != self.ybinList[-1]:
                raise ValueError('Last edges of Y SLICES does not match axis (%s vs %s)'%(binning_dict['Y']['SLICES'][-1], self.ybinList[-1]))
            slices = binning_dict['Y']['SLICES']
            idxs = [0, self.ybinList.index(slices[1]), self.ybinList.index(slices[2]), len(self.ybinList)-1]
        else:
            slices, idxs = self._autoYslices()

        return slices, idxs

    def _autoYslices(self):
        nbins = len(self.ybinList)-1
        idxs = [0, int(nbins/4), int(nbins/4)+int(nbins/3), nbins]
        slices = [int(self.ybinList[i]) for i in idxs]

        return slices, idxs

    @property
    def xSliceIdx(self):
        return [0,self.GlobalXbinIdx(0,'SIG'),self.GlobalXbinIdx(-1,'SIG'),len(self.xbinList)-1]

    @property
    def xSlices(self):
        return [int(self.xbinList[i]) for i in self.xSliceIdx]

    def GlobalXbinIdx(self,xbin,c):
        '''Evaluate for the bin - a bit tricky since it was built with separate categories.
        Returns the index of upper global wall - AKA index of bin we want in the the
        histogram (remember those bin indices start at 1 not 0).

        Args:
            xbin (int): Category local bin index
            c ([type]): Category name - LOW, SIG, or HIGH

        Returns:
            int: Global index
        '''
        return self.xbinList.index(self.xbinByCat[c][xbin])

    def xcatFromGlobal(self,xbin):
        n_low_bins = len(self.xbinByCat['LOW'])-1
        n_sig_bins = len(self.xbinByCat['SIG'])-1
        if xbin < n_low_bins+1:
            return xbin,'LOW'
        elif xbin < n_low_bins+n_sig_bins+1:
            return xbin-n_low_bins,'SIG'
        else:
            return xbin-n_low_bins-n_sig_bins,'HIGH'

    @property
    def xbinList(self):
        '''
        Returns:
            list: X axis binning dict converted from a dictionary of the regions to
            a continuous list of bin edges for the full X axis.
        '''
        return concat_bin_dicts(self.xbinByCat)

    def GetBinCenterBase(self,ibin,binlist):
        if ibin < 1: raise ValueError('Binning is indexed at 1 for compatibility with ROOT.')
        return (binlist[ibin]+binlist[ibin-1])/2

    def GetBinCenterX(self,ibin,cat):
        return self.GetBinCenterBase(ibin,self.xbinByCat[cat])

    def GetBinCenterY(self,ibin):
        return self.GetBinCenterBase(ibin,self.ybinList)

    def CreateHist(self,name,cat=''):
        if cat != '':
            xbins = self.xbinByCat[cat]
        else:
            xbins = self.xbinList

        return ROOT.TH2F(name,name,
                        len(xbins)-1,
                        array.array('d',xbins),
                        len(self.ybinList)-1,
                        array.array('d',self.ybinList)
        )

def create_RRV_base(name,title,bins):
    '''Generically create a RooRealVar with the specified bin edges.

    Args:
        name (str): Object name.
        title (str): Object title.
        bins ([type]): List of bin edges.

    Returns:
        RooRealVar: 
    '''
    RRV = ROOT.RooRealVar(name,title,bins[0],bins[-1])
    bin_array = array.array('d',bins)
    RooBinning = ROOT.RooBinning(len(bins)-1,bin_array)
    RRV.setBinning(RooBinning)
    return RRV

def parse_binning_info(binDict):
    '''If running a blinded fit, then we want to do a combined fit over 
    two categories: below and above the signal region. This requires
    generating histograms in those three regions and it's useful
    to have different binning for all of those. If the signal region
    is not blinded then we can fit the entire region but it's convenient
    to still do three categories for later parts of the code. So here are
    the options.
    1) It may be desired or even required to bin the fit in three categories
    each with its own binning structure (say statistics are good in 
    region below the signal region but bad above it so you'd like to
    use more and fewer bins, respectively). 
    2) Additionally, variable binning can be used for each category. 
    3) Single binning strategy across all three regions and only defined
    once in the configuration.
    The only requirement for any of this is that the bin walls of the new
    binning match a bin wall of the input histograms (you can't make bins smaller or split bins!)

    For config input, this breaks down into
    - Standard bins over one category - one NBINS,MIN,MAX
    - Standard bins over three categories - three NBINS,MIN,MAX (organized by dict)
    - Custom bins over one category - list of bin walls
    - Custom bins over three categories - three lists of bin walls (organized by dict)

    Args:
        binDict (dict): Usually the "BINNING" section of the configuration file.

    Returns:
        tuple: In order - new bins in the X axis, new bins in the Y axis
    '''
    for v in ['X','Y']:
        axis = binDict[v]
        if (v == 'X') and ('LOW' in axis.keys()) and ('SIG' in axis.keys()) and ('HIGH' in axis.keys()):
            new_bins = {c:parse_axis_info(axis[c]) for c in ['LOW','SIG','HIGH']}
        else:
            new_bins = parse_axis_info(axis)
            
        if v == 'X':
            if isinstance(new_bins,list):
                newXbins = binlist_to_bindict(new_bins,axis['SIGSTART'],axis['SIGEND'])
            else:
                newXbins = new_bins
        elif v == 'Y': newYbins = new_bins

    return newXbins,newYbins
    
def parse_axis_info(axisDict):
    '''Return list of bin edges based off of the scheme specified in the "BINNING"
    section of the json config.

    Args:
        axisDict (dict): config["BINNING"]["X" or "Y"]

    Raises:
        RuntimeError: If binning scheme not specified with the correct syntax.

    Returns:
        list: Bin edges.
    '''
    if 'BINS' in axisDict.keys():
        new_bins = axisDict['BINS']
    elif ('MIN' in axisDict.keys()) and ('MAX' in axisDict.keys()) and ('NBINS' in axisDict.keys()):
        new_width = float(axisDict['MAX']-axisDict['MIN'])/float(axisDict['NBINS'])
        new_bins = [axisDict['MIN'] + new_width*i for i in range(axisDict['NBINS'])] + [axisDict['MAX']]
    else:
        raise RuntimeError('Bins not specified correctly in BINNING section of config.')
    return new_bins

def binlist_to_bindict(binList, sigLow, sigHigh):
    '''Convert a list of bins into a dictionary with keys LOW, SIG, and HIGH
    where the three regions are separated by the values sigLow and sigHigh.
    The sigLow and sigHigh values should be bin edges.

    Args:
        binList ([type]): 
        sigLow (int): [description]
        sigHigh (int): [description]

    Raises:
        ValueError: If sigLow or sigHig is not in binList.

    Returns:
        dict: Dictionary of shape {'LOW':[...],'SIG':[...],'HIGH':[...]}
    '''
    return_bins = {'LOW':[],'SIG':[],'HIGH':[]}
    for s in [sigLow,sigHigh]:
        if s not in binList:
            raise ValueError('The signal region edges must be in the list of bin edges. The value %s is not in the provided list of bin edges (%s).'%(s,binList))
    for b in binList:
        if b <= sigLow:
            return_bins['LOW'].append(b)
        if b >= sigLow and b <= sigHigh:
            return_bins['SIG'].append(b)
        if b >= sigHigh:
            return_bins['HIGH'].append(b)

    return return_bins 

def concat_bin_dicts(binDict):
    '''Convert a dictionary of shape {'LOW':[...],'SIG':[...],'HIGH':[...]}
    to a list, avoiding overlapping values and the start and ends of each dict entry.

    Args:
        binDict (dict): Input dictionary

    Returns:
        list: List of bin edges concatenated from the dictionary entries.
    '''
    bins_list = list(binDict['LOW']) # need list() to make a copy - not a reference
    for c in ['SIG','HIGH']:
        bins_list.extend(binDict[c][1:])
    return bins_list

def concat_bin_lists(binLists):
    '''Convert a list of separate bin edges to one continuous list.
    Will check that the start and ends of each sublist in binList are the same (continuous).

    Args:
        binList (list): Input list of bin edges.

    Raises:
        ValueError: If bins in binList are not continuous along axis.

    Returns:
        list: List of bin edges concatenated from the binList entries.
    '''
    bins_list = binLists[0]
    for b in binLists[1:]:
        if b[0] != bins_list[-1]:
            raise ValueError('Bins in binLists are not continuous along axis.')
        bins_list.extend(b[1:])
    return bins_list

def get_bins_from_hist(XYZ,h):
    '''Get list of all bin edges from a histogram.
    Specify which axis in the histogram via XYZ.

    Args:
        XYZ (str): "X", "Y", or "Z".
        h (TH1): Histogram to get binning from.

    Returns:
        list: List of all bin edges.
    '''
    nbins = getattr(h,'GetNbins'+XYZ)()
    axis = getattr(h,'Get%saxis'%XYZ)()
    bins = [axis.GetBinLowEdge(i+1) for i in range(nbins)]+[axis.GetBinUpEdge(nbins)]
    return bins

def histlist_to_binlist(XYZ,histList):
    '''Input a list of separate histograms and return one continuous list of bins along the XYZ axis.
    Will check that the start and ends of each sublist in binList are the same (continuous).

    Args:
        XYZ (str): "X", "Y", or "Z".
        histList (list(TH1)): List of histograms.

    Returns:
        list: List of bin edges. 
    '''
    binList = [get_bins_from_hist(XYZ,h) for h in histList]
    return concat_bin_lists(binList)

def stitch_hists_in_x(name,binning,histList,blinded=[]):
    '''Required that histList be in order of desired stitching
    `blinded` is a list of the index of regions you wish to skip/blind.

    Args:
        name (str): Name of output histogram.
        binning (Binning): Binning storage object.
        histList (list(TH2)): List of histograms to stitch together. Binning must be continious.
        blinded (list(int), optional): List of indexes of histList which should be dropped/blinded. Defaults to [].

    Raises:
        ValueError: If X axis bins stitched together from histList are not the same as the input template.

    Returns:
        TH2: Output stitched histograms.
    '''
    stitched_hist = binning.CreateHist(name)
    stitched_hist.Reset()
    # Sanity checks
    histListBins = histlist_to_binlist("X",histList)
    if histListBins != get_bins_from_hist("X",stitched_hist):
        raise ValueError('X axis bins stitched together from histList are not the same as the input template.\n%s vs %s'%(histListBins,get_bins_from_hist("X",stitched_hist)))
    # Stitch
    bin_jump = 0
    for i,h in enumerate(histList):
        if i in blinded:
            bin_jump += histList[i].GetNbinsX()
            continue
        
        for ybin in range(1,h.GetNbinsY()+1):
            for xbin in range(1,h.GetNbinsX()+1):
                stitched_xindex = xbin + bin_jump
                stitched_hist.SetBinContent(stitched_xindex,ybin,h.GetBinContent(xbin,ybin))
                stitched_hist.SetBinError(stitched_xindex,ybin,h.GetBinError(xbin,ybin))

        bin_jump += histList[i].GetNbinsX()

    return stitched_hist

def make_blinded_hist(h,sigregion):
    '''Clone histogram (h) and set the bins in range
    sigregion[0] to sigregion[1] on the X axis to zero.

    Args:
        h (TH2): Input unblinded histogram
        sigregion (list(float)): Axis range to blind.
        Edges must line up with h's bin edges. Must have length of 2 (lower and upper bound).

    Raises:
        ValueError: If signal region edges do not line up with histogram bin edges.
        IndexError: If sigregion is not of length 2.

    Returns:
        TH2: Modified version of h with specified region blinded.
    '''
    blindedHist = h.Clone()
    blindedHist.Reset(); blindedHist.Sumw2()
    h.SetName(blindedHist.GetName()+'_unblinded') # Need to change nominal hist name or we'll get a memory leak

    if (sigregion[0] not in get_bins_from_hist("X",h)) or (sigregion[1] not in get_bins_from_hist("X",h)):
        raise ValueError('Signal region edges %s do not line up with histogram bin edges.')
    if len(sigregion) != 2:
        raise IndexError('Signal region must be specified by list of length 2.')

    for binY in range(1,h.GetNbinsY()+1):
        for binX in range(1,h.GetNbinsX()+1):
            if h.GetXaxis().GetBinUpEdge(binX) <= sigregion[0] or h.GetXaxis().GetBinLowEdge(binX) >= sigregion[1]:
                if h.GetBinContent(binX,binY) > 0:
                    blindedHist.SetBinContent(binX,binY,h.GetBinContent(binX,binY))
                    blindedHist.SetBinError(binX,binY,h.GetBinError(binX,binY))

    return blindedHist

def copy_hist_with_new_bins(copyName,XorY,inHist,new_bins):
    '''Make a copy of a 2D histogram with new bins specified for a given axis (X or Y).
    New bins must be larger than the old bins and the edges of new bins must line up with 
    existing edges (no finer binning and no splitting bins).

    Args:
        copyName (str): Name of copy.
        XorY (str): "X" or "Y" to change which axis is rebinned.
        inHist (TH2): Input histogram to rebin.
        new_bins (list): New list of bin edges.

    Raises:
        ValueError: If XorY is not "X" or "Y".
        ValueError: If the requested rebinning does not align bin edges with the available input bin edges.

    Returns:
        TH2: Copy of histogram with new binning scheme. Note that the number of entries
            will not be correct but integrated yield will be. 
    '''
    if XorY not in ["X","Y"]:
        raise ValueError('Arg XorY is not "X" or "Y".')
    axis_to_rebin = XorY
    axis_to_hold = "X" if XorY=="Y" else "Y"
    
    static_array = array.array('f',get_bins_from_hist(axis_to_hold,inHist))
    static_nbins = len(static_array)-1
    rebin_array = array.array('f',new_bins)
    rebin_nbins = len(rebin_array)-1 

    # Use copyName with _temp to avoid overwriting if inHist has the same name
    # We can do this at the end but not before we're finished with inHist
    if XorY == "X":
        hist_copy = ROOT.TH2F(copyName+'_temp',copyName+'_temp',rebin_nbins,rebin_array,static_nbins,static_array)
    else:
        hist_copy = ROOT.TH2F(copyName+'_temp',copyName+'_temp',static_nbins,static_array,rebin_nbins,rebin_array)
    hist_copy.Sumw2()
    hist_copy.GetXaxis().SetName(inHist.GetXaxis().GetName())
    hist_copy.GetYaxis().SetName(inHist.GetYaxis().GetName())
    old_axis = getattr(inHist,'Get%saxis'%axis_to_rebin)()
    rebin_axis = getattr(hist_copy,'Get%saxis'%axis_to_rebin)()

    # Loop through the old bins
    for static_bin in range(1,static_nbins+1):
        # print 'Bin y: ' + str(binY)
        for rebin in range(1,rebin_nbins+1):
            new_bin_content = 0
            new_bin_errorsq = 0
            new_bin_min = rebin_axis.GetBinLowEdge(rebin)
            new_bin_max = rebin_axis.GetBinUpEdge(rebin)

            # print '\t New bin x: ' + str(newBinX) + ', ' + str(newBinXlow) + ', ' + str(newBinXhigh)
            for old_bin in range(1,old_axis.GetNbins()+1):
                old_bin_min = old_axis.GetBinLowEdge(old_bin)
                old_bin_max = old_axis.GetBinUpEdge(old_bin)
                if old_bin_min >= new_bin_max:
                    break
                elif old_bin_min >= new_bin_min and old_bin_min < new_bin_max:
                    if old_bin_max <= new_bin_max:
                        if axis_to_rebin == "X":
                            new_bin_content += inHist.GetBinContent(old_bin,static_bin)
                            new_bin_errorsq += inHist.GetBinError(old_bin,static_bin)**2
                        else:
                            new_bin_content += inHist.GetBinContent(static_bin,old_bin)
                            new_bin_errorsq += inHist.GetBinError(static_bin,old_bin)**2
                    elif old_bin_max > new_bin_max:
                        raise ValueError(
                            '''The requested %s rebinning does not align bin edges with the input bin edge.
                            Cannot split input bin [%s,%s] with output bin [%s,%s]'''%(axis_to_rebin,old_bin_min,old_bin_max,new_bin_min,new_bin_max))
                elif old_bin_min <= new_bin_min and old_bin_max > new_bin_min:
                    raise ValueError(
                        '''The requested %s rebinning does not align bin edges with the input bin edge.
                        Cannot split input bin [%s,%s] with output bin [%s,%s]'''%(axis_to_rebin,old_bin_min,old_bin_max,new_bin_min,new_bin_max))

            # print '\t Setting content ' + str(newBinContent) + '+/-' + str(sqrt(newBinErrorSq))
            if new_bin_content > 0:
                if axis_to_rebin == "X":
                    hist_copy.SetBinContent(rebin,static_bin,new_bin_content)
                    hist_copy.SetBinError(rebin,static_bin,sqrt(new_bin_errorsq))
                else:
                    hist_copy.SetBinContent(static_bin,rebin,new_bin_content)
                    hist_copy.SetBinError(static_bin,rebin,sqrt(new_bin_errorsq))

    # Will now set the copyName which will overwrite inHist if it has the same name
    hist_copy.SetName(copyName)
    hist_copy.SetTitle(copyName)
    return hist_copy

def get_min_bin_width(hist):
    '''Get the minimum width among all bins in a 1D histogram.

    Args:
        hist (TH1): 1D histogram to analyze.

    Raises:
        TypeError: If number of dimensions in hist is != 1.

    Returns:
        int: Minimum bin width.
    '''
    if hist.GetDimension() != 1:
        raise TypeError('Only 1D histograms can be analyzed for minimum bin width.')
    use_width = 10**6
    for ibin in range(1,hist.GetNbinsX()+1):
        if hist.GetBinWidth(ibin) < use_width:
            use_width = hist.GetBinWidth(ibin)
    return int(use_width)

def convert_to_events_per_unit(hist,width=None):
    '''Convert the bin contents of a 1D histogram so they are normalized to the
    width of the narrowest bin. Only useful if the bin widths are variable.

    Args:
        hist (TH1): 1D histogram to manipulate.
        width (int, optional): Override the automatic determination of the minimum bin width.
        Defaults to None in which case get_min_bin_width is used.

    Raises:
        TypeError: If number of dimensions in hist is != 1.

    Returns:
        TH1: Histogram clone with the proper bin contents renormalized.
    '''
    if hist.GetDimension() != 1:
        raise TypeError('Only 1D histograms can be analyzed for minimum bin width.')
    if width == None:
        use_width = get_min_bin_width(hist)
    else:
        use_width = width

    converted = hist.Clone()
    for ibin in range(1,hist.GetNbinsX()+1):
        if hist.GetBinWidth(ibin) == use_width:
            continue
        else:
            factor = use_width/hist.GetBinWidth(ibin)
            new_content = factor * converted.GetBinContent(ibin)
            new_error = factor * converted.GetBinError(ibin)
            converted.SetBinContent(ibin,new_content)
            converted.SetBinError(ibin,new_error)
    
    return converted

def zero_negative_bins(name,inhist):
    '''Set all negative bins in 2D inhist to zero with zero error.

    Args:
        name (str): Name of returned histogram.
        inhist (TH2): Input histogram.

    Returns:
        TH2: Clone of input histogram with negative bins set to zero.
    '''
    outhist = inhist.Clone(name)
    for ix in range(1,inhist.GetNbinsX()+1):
        for iy in range(1,inhist.GetNbinsY()+1):
            if inhist.GetBinContent(ix,iy) < 0:
                outhist.SetBinContent(ix,iy,0)
                outhist.SetBinError(ix,iy,0)

    return outhist

def remap_binlist(binList,new_min=0.0,new_max=1.0):
    '''Remap a list of bin edges to [new_min, new_max].

    Args:
        binList (list(float)): List of bin edges (including first and last).
        new_min (float, optional): New minimum. Defaults to 0.0.
        new_max (float, optional): New maximum. Defaults to 1.0.

    Returns:
        list(float): List of bin edges mapped to new range.
    '''
    minimum = binList[0]
    new_min = float(new_min); new_max = float(new_max)
    new_length = new_max - new_min
    scale = (binList[-1]-binList[0]) / new_length
    new_vals = []
    for old_val in binList:
        new_vals.append( (old_val-minimum)/scale + new_min )
    return new_vals

def remap_hist_axis(hist,new_min=0,new_max=1):
    '''Remap axes of a 2D histogram to [new_min, new_max].

    Args:
        hist (TH2): Histogram to remap.
        new_min (int, optional): New minimum. Defaults to 0.
        new_max (int, optional): New maximum. Defaults to 1.

    Raises:
        TypeError: If number of dimensions in hist is != 1.

    Returns:
        TH1: [description]
    '''
    ybins = array.array('d', remap_binlist(
                                get_bins_from_hist('Y',hist),
                                new_min, new_max
                             )
                        )
    xbins = array.array('d', remap_binlist(
                                get_bins_from_hist('X',hist),
                                new_min, new_max
                             )
                        )

    remap = ROOT.TH2F(hist.GetName()+'_unit',hist.GetName()+'_unit',len(xbins)-1,xbins,len(ybins)-1,ybins)
    remap.Sumw2()

    for xbin in range(1,hist.GetNbinsX()+1):
        for ybin in range(1,hist.GetNbinsY()+1):
            remap.SetBinContent(xbin,ybin,hist.GetBinContent(xbin,ybin))
            remap.SetBinError(xbin,ybin,hist.GetBinError(xbin,ybin))

    return remap
