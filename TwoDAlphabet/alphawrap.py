from collections import OrderedDict
from TwoDAlphabet.helpers import roofit_form_to_TF1
from ROOT import RooRealVar, RooFormulaVar, RooArgList, RooParametricHist2D, RooParametricHist3D, RooConstVar, TFormula, RooAddition
from TwoDAlphabet.binning import copy_hist_with_new_bins
import itertools

class Generic2D(object):
    '''Wraps binned distributions in a common type so that
    distributions can easily be manipulated and compared. While a full distribution
    is always input, this class will actually store three sets of bins, one each
    for the 'LOW', 'SIG', and 'HIGH' categories of the X axis. The three sets are considered
    one "object" to the user. Combining an object of class Generic2D with another
    (see `_manipulate` method) creates a new Generic2D object with new
    RooFormulaVars corresponding to the manipulation. To avoid python's
    garbage collection of RooFit objects, assign each instance of this class to a persistent
    variable.

    Attributes:
        name (str): Unique name of object which will be prepended to all associated RooFit objects.
        binning (Binning): Binning object.
        is3D (bool): Whether the model is 3D (TH3*).
        nuisances (list): All tracked nuisance dictionaries.
        binVars (OrderedDict): Ordered dictionary of all RooAbsArgs representing the bins of the space.
        binArgLists (dict): Dict mapping of the subspaces (LOW, SIG, HIGH) to the RooArgList of the RooAbsArgs in the subspace.
        rph (dict): Dict mapping of the subspaces (LOW, SIG, HIGH) to the RooParametricHist2D objects of the subspaces.
        forcePositive (bool): Option to ensure bin values cannot be negative.
    '''

    def __init__(self,name,binning,forcePositive=True):
        self.name = name
        self.binning = binning
        self.is3D = getattr(binning, 'is3D', False)  # [3D]
        self.nuisances = []
        self.binVars = OrderedDict()
        self.subspaces = [str(c) for c in binning.xbinByCat]
        self.binArgLists = {c:None for c in self.subspaces}
        self.rph = {c:None for c in self.subspaces}
        self.forcePositive = forcePositive
        self._varStorage = [] # only used by AddShapeTemplates

    def _binIter(self,cat):
        '''Provides an iterator over (xbin, ybin, zbin) for every bin in subspace `cat`.
        Uses `yield` so that we can move the iteration to this function rather than 
        checking is3D inside the _manipulate function. 
        '''
        nx = len(self.binning.xbinByCat[cat])
        ny = len(self.binning.ybinList)
        if self.is3D:
            nz = len(self.binning.zbinList)
            for zbin in range(1,nz):
                for ybin in range(1,ny):
                    for xbin in range(1,nx):
                        yield xbin, ybin, zbin
        else:
            for ybin in range(1,ny):
                for xbin in range(1,nx):
                    yield xbin, ybin, None

    def _binSuffix(self,xbin,ybin,zbin=None):
        return '%s-%s-%s'%(xbin,ybin,zbin) if self.is3D else '%s-%s'%(xbin,ybin)

    def _binName(self,cat_name,xbin,ybin,zbin=None):
        return '%s_bin_%s'%(cat_name, self._binSuffix(xbin,ybin,zbin))

    def _binContent(self,hist,xbin,ybin,zbin=None):
        return hist.GetBinContent(xbin,ybin,zbin) if self.is3D else hist.GetBinContent(xbin,ybin)

    def _manipulate(self,name,other,operator=''):
        '''Base method to create a new Generic2D object combining `self` and `other`.'''
        out = Generic2D(name,self.binning,self.forcePositive)
        for cat in self.subspaces:
            new_cat_name = name+'_'+cat
            self_cat_name = self.name+'_'+cat
            other_cat_name = other.name+'_'+cat
            for xbin, ybin, zbin in self._binIter(cat):
                new_bin_name   = self._binName(new_cat_name,xbin,ybin,zbin)
                self_bin_name  = self._binName(self_cat_name,xbin,ybin,zbin)
                other_bin_name = self._binName(other_cat_name,xbin,ybin,zbin)
                out.binVars[new_bin_name] = RooFormulaVar(
                                                new_bin_name, new_bin_name, '@0%s@1'%operator,
                                                RooArgList(
                                                    self.binVars[self_bin_name],
                                                    other.binVars[other_bin_name]))

        all_nuisances = self.nuisances+other.nuisances
        for nuisance in all_nuisances:
            if nuisance['name'] in [n['name'] for n in out.nuisances]:
                raise RuntimeError('Already tracking nuisance %s. Printing all nuisances...\n\t'%(nuisance['name'],all_nuisances))

            out.nuisances.append(nuisance)

        return out

    def Add(self,name,other,factor='1'):
        '''Add `self` with `other`. Optionally change the
        factor in front of `other` (defaults to 1). This option is
        primarly for the case of subtracting `other` from `self`.

        Args:
            name (str): Unique name for the new output Generic2D object.
            other (Generic2D): Object to add to `self`.
            factor (str, optional): Factor to include in front of `other` in the combination. Defaults to '1'.
                Primary use case is "-1" which will subtract `other` from `self`.

        Returns:
            Generic2D: Object containing the addition of `self` and `other`.
        '''
        if factor.startswith('-'):
            op = '%s*'%factor
        elif factor == '1':
            op = '+'
        else:
            op = '+%s*'%factor
        return self._manipulate(name,other,op)

    def Multiply(self,name,other):
        '''Multiply `self` with `other`.

        Args:
            name (str): Unique name for the new output Generic2D object.
            other (Generic2D): Object to multiply `self` by.

        Returns:
            Generic2D: Object containing the multiplication of `self` and `other`.
        '''
        return self._manipulate(name,other,'*')

    def Divide(self,name,other):
        '''Divide `self` by `other`.

        Args:
            name (str): Unique name for the new output Generic2D object.
            other (Generic2D): Object to divide `self` by.

        Returns:
            Generic2D: Object containing the division of `self` by `other`.
        '''
        return self._manipulate(name,other,'/')

    def RooParametricHist(self,name=''):
        '''Produce a RooParametricHist2D (2D mode) or RooParametricHist3D (3D mode)
        filled with this object's binVars, plus the matching RooAddition norm.

        Returns:
            RooParametricHist(2,3)D, RooAddition
        '''
        out_rph = {}
        out_add = {}
        for cat in self.subspaces:
            cat_name = self.name+'_'+cat
            cat_hist = self.binning.CreateHist(cat_name+'_temp',cat)  # TH2F or TH3F depending on binning
            obj_name = '%s_%s'%(name if name != '' else self.name, cat)

            self.binArgLists[cat] = RooArgList()
            for xbin, ybin, zbin in self._binIter(cat):
                self.binArgLists[cat].add(self.binVars[self._binName(cat_name,xbin,ybin,zbin)])

            if self.is3D:
                out_rph[cat] = RooParametricHist3D(
                            obj_name, obj_name,
                            self.binning.xVars[cat],
                            self.binning.yVar,
                            self.binning.zVar,
                            self.binArgLists[cat], cat_hist
                )
            else:
                out_rph[cat] = RooParametricHist2D(
                            obj_name, obj_name,
                            self.binning.xVars[cat],
                            self.binning.yVar,
                            self.binArgLists[cat], cat_hist
                )
            out_add[cat] = RooAddition(obj_name+'_norm',obj_name+'_norm',self.binArgLists[cat])
        return out_rph, out_add

    def getBinVal(self,xbin,ybin,zbin=None):
        '''Get bin value (for the current parameter values) in bin (xbin, ybin[, zbin]).
        Args:
            xbin (int): Indexed at 1 for ROOT compatibility.
            ybin (int): Indexed at 1 for ROOT compatibility.
            zbin (int): (optional) Z-bin if using 3D template.

        Returns:
            float: Current bin value.        
        '''
        return self.getBinVar(xbin,ybin,zbin=zbin).getValV()

    def getBinVar(self,xbin,ybin,c='',zbin=None):
        '''Get the bin variable associated with (xbin,ybin).
        The `xbin` and `ybin` args are assumed to be for global bin
        numbers but can be for a given category ("LOW", "SIG", or "HIGH") 
        if specified with the `c` option.

        Args:
            xbin ([type]): [description]
            ybin ([type]): [description]
            c (str, optional): One of "LOW", "SIG", or "HIGH" which will
                cause xbin and ybin to be interpreted as indexes for the given subspace.
                Defaults to '' in which case xbin and ybin are treated as global.

        Returns:
            RooFormulaVar: RooFit object for the requested bin.
        '''
        if c == '': # using a global xbin that needs to be translated
            xbin, c = self.binning.xcatFromGlobal(xbin)
        return self.binVars[self._binName(self.name+'_'+c,xbin,ybin,zbin)]


class ParametricFunction(Generic2D):
    def __init__(self,name,binning,formula,constraints={},forcePositive=True):
        '''Represents parametric functions as a group of RooFormulaVars which
        create a binned distribution and which change
        as the underlying function parameters change. Set parameter specific
        values by specifying the `constraints` argument with a dict formatted as

        \code{.json}
            {0: {
                'constraint': 'flatParam' or 'param <mu> <sigma>',
                'MIN' = -1000,
                'MAX' = 1000,
                'NOM' = 0,
                'ERROR' = 0.1
            } }
        \endcode

        The 'constraint' can only be 'flatParam' or 'param <mu> <sigma>' (options in the Combine card) 
        which represent "no constraint" and "Gaussian constraint centered at <mu> and with width <sigma>", respectively.

        @param name (str): Unique name for the new object.
        @param formula (str): Must reference fit parameters by ordinal with @. Use "x" and "y" to represent
                the "x" and "y" axes of the space. All other terms are indexed starting at 0. Ex. "@0 + x*@1 +y*@2".
        @param constraints (dict, optional): Map of formula parameters to constraint information. Defaults to {} in which
                case the constraint will be flat, the starting value of the parameter will be 0 with a step size of 0.1,
                and the range of the parameter will be [-1000,1000]. 
            
        @param forcePositive (bool, optional). Defaults to True in which case the bin values will be lower bound by 1e-9.
        '''
        super(ParametricFunction,self).__init__(name,binning,forcePositive)
        self.formula = formula
        self.nuisances = self._createFuncVars(constraints)
        self.arglist = RooArgList()
        for n in self.nuisances: self.arglist.add(n['obj'])

        for cat in self.subspaces:
            cat_name = name+'_'+cat
            for xbin, ybin, zbin in self._binIter(cat):
                bin_name = self._binName(cat_name,xbin,ybin,zbin)
                centers = self.mappedBinCenter(xbin,ybin,cat,zbin)  # (x,y) or (x,y,z)
                if forcePositive: final_formula = "max(1e-9,%s)"%(self._replaceXYZ(*centers))
                else:             final_formula = self._replaceXYZ(*centers)

                self.binVars[bin_name] = RooFormulaVar(
                    bin_name, bin_name,
                    final_formula,
                    self.arglist
                )

    def _replaceXYZ(self,x,y,z=None):
        '''Find and replace "x", "y", and optionally "z" in the input formula
        with this method's arguments (floats) which should
        correspond to bin centers.

        Args:
            x (float): Value along x axis to evaluate.
            y (float): Value along y axis to evaluate.
            z (float): Value along z axis to evaluate (only if 3D).

        Returns:
            str: Formula with "x" and "y" (and optionally "z") replaced by provided numerical values.
        '''

        f = self.formula.replace(' ','')
        tokens = [('x',x),('y',y)]
        if z is not None:
            tokens.append(('z',z))
        for tok,val in tokens:
            f = f.replace('+'+tok,'+%s'%val).replace('*'+tok,'*%s'%val)
            f = f.replace('-'+tok,'-%s'%val).replace('/'+tok,'/%s'%val)
            f = f.replace('('+tok,'(%s'%val)
        return f

    def _replaceXY(self,x,y):  # kept in case of old _replaceXY calls that were missedd
        return self._replaceXYZ(x,y)

    def getNparams(self):
        '''Get the number of parameters in the formula (not counting "x" or "y" or "z").
        Converts the formula formating and uses ROOT's TFormula to count.

        Returns:
            int: Number of parameters in the fit (not counting "x" or "y" or "z").
        '''
        return TFormula('tempFormula',roofit_form_to_TF1(self.formula)).GetNpar()

    def _createFuncVars(self,constraints):
        '''Creates the nuisances list of the function variables (RooRealVars)
        and associated meta data (nuisance name, constraint).

        Args:
            constraints (dict): Information of which parameters to constrain differently from the default.
                By default, the constraint will be flat, the starting value of the parameter will be 0 with a step size of 0.1,
                and the range of the parameter will be [-1000,1000].

        Returns:
            list: List of dictionaries with keys "name" (str), "obj" (RooRealVar), "constraint" (str).
        '''
        out = []
        for i in range(self.getNparams()):
            name = '%s_par%s'%(self.name,i)
            constraint = 'flatParam'; MIN = -1000; MAX = 1000; NOM = 0.1; ERROR = 0.1
            if i in constraints:
                if 'constraint' in constraints[i]: constraint = constraints[i]['constraint']
                if 'MIN' in constraints[i]: MIN = constraints[i]['MIN']
                if 'MAX' in constraints[i]: MAX = constraints[i]['MAX']
                if 'NOM' in constraints[i]: NOM = constraints[i]['NOM']
                if 'ERROR' in constraints[i]: ERROR = constraints[i]['ERROR']

            this_out = {'name':name, 'obj': RooRealVar(name,name,NOM,MIN,MAX), 'constraint': constraint}
            this_out['obj'].setError(ERROR)
            out.append(this_out)
        return out

    def mappedBinCenter(self,xbin,ybin,cat,zbin=None):
        '''Convert bin indices to bin-center values with each axis mapped to [0,1].

        Returns:
            tuple of floats: (x, y) in 2D mode, or (x, y, z) in 3D mode.
        '''
        x_center = self.binning.GetBinCenterX(xbin,cat)
        y_center = self.binning.GetBinCenterY(ybin)

        x_min = self.binning.xbinList[0]
        y_min = self.binning.ybinList[0]
        x_range = self.binning.xbinList[-1] - x_min
        y_range = self.binning.ybinList[-1] - y_min

        # Remap to [-1,1]
        x_center_mapped = float(x_center - x_min)/x_range #float cast prevents returning zero if bin edges are ints
        y_center_mapped = float(y_center - y_min)/y_range

        if self.is3D and zbin is not None:
            z_center = self.binning.GetBinCenterZ(zbin)
            z_min = self.binning.zbinList[0]
            z_range = self.binning.zbinList[-1] - z_min
            z_center_mapped = float(z_center - z_min)/z_range
            return x_center_mapped, y_center_mapped, z_center_mapped

        return x_center_mapped,y_center_mapped

    def setFuncParam(self,parIdx,value):
        '''Set the value of a given ROOT.RooRealVar object within a ParametricFunction

        Args:
            parIdx (int,str): Parameter index to access, or parameter name.
            value (float): Value to assign.

        Raises:
            RuntimeError: If the parameter could not be found.

        Returns:
            None
        '''
        parfound = False
        for i,n in enumerate(self.nuisances):
            # user supplies full parameter name, e.g. 'Background_CR_rpfT_par0'
            if n['name'] == parIdx:
                self.nuisances[i]['obj'].setVal(value)
                parfound = True
                break
            # user supplies only parameter index, e.g. '0', 0
            elif n['name'].endswith('_par%s'%parIdx):
                self.nuisances[i]['obj'].setVal(value)
                parfound = True
                break
        if parfound == False:
            raise RuntimeError('Could not find par%s in set of nuisances:\n\t%s'%(parIdx,[n['name'] for n in self.nuisances]))

class SemiParametricFunction(ParametricFunction,Generic2D):
    def __init__(self,name,inhist,binning,formula,constraints={},forcePositive=True,funcCeiling=10.):
        '''A hybrid of ParametricFunction and BinnedDistribution classes. 
        Uses former (RooFormulaVar) if bin count<funcCeiling, later (RooRealVar) if not.
        Args: 
            same as in Parametric Function and BinnedDistribution classes, except funcCeiling
            funcCeiling (float, optional). Bins with content >funcCeiling will use floating bin parametrization
        instead of a functional form. Enables using functional form only in the tails of the distribution. Defaults to 10
        '''
        Generic2D.__init__(self,name,binning,forcePositive)
        self.formula = formula #This is already done in init
        self.nuisances = self._createFuncVars(constraints)
        self.arglist = RooArgList()
        for n in self.nuisances: self.arglist.add(n['obj'])

        for cat in self.subspaces:
            cat_name = name+'_'+cat
            cat_hist = copy_hist_with_new_bins(cat_name,'X',inhist,self.binning.xbinByCat[cat])
            for xbin, ybin, zbin in self._binIter(cat):
                content = self._binContent(cat_hist,xbin,ybin,zbin)
                bin_name = self._binName(cat_name,xbin,ybin,zbin)
                if(content<funcCeiling):
                    centers = self.mappedBinCenter(xbin,ybin,cat,zbin)
                    if forcePositive:
                        final_formula = "max(1e-9,%s)"%(self._replaceXYZ(*centers))
                    else:
                        final_formula = self._replaceXYZ(*centers)

                    self.binVars[bin_name] = RooFormulaVar(
                        bin_name, bin_name,
                        final_formula,
                        self.arglist
                    )
                else:
                    self.binVars[bin_name] = RooRealVar(bin_name, bin_name, content, 1e-6, 1e9)
                    self.nuisances.append({'name':bin_name, 'constraint':'flatParam', 'obj': self.binVars[bin_name]})

class BinnedDistribution(Generic2D):
    def __init__(self,name,inhist,binning,constant=False,forcePositive=True):
        '''Represents a binned distribution as a group of RooRealVar parameters.
        If constant == False, each bin is considered an unconstrained parameter of the model.

        Args:
            name (str): Unique name for the new object.
            inhist (TH2): Input 2D histogram to build set of variables.
            binning (Binning): Binning object used to create LOW, SIG, HIGH regions along X axis.
            constant (bool, optional): If true, use RooConstVars for bins. Defaults to False and RooRealVars are used.
            forcePositive (bool, optional). Defaults to True in which case the bin values will be lower bound by 1e-9
                and any shape templates will asymptotically approach zero as the associated nuisance increases/decreases.
        '''
        super(BinnedDistribution,self).__init__(name,binning,forcePositive=forcePositive)
        zero_thresh = 3**(3 if self.is3D else 2) - 2  # almost all neighbors zero, replaced the old hard-coded value of 7
        for cat in self.subspaces:
            cat_name = name+'_'+cat
            cat_hist = copy_hist_with_new_bins(cat_name,'X',inhist,self.binning.xbinByCat[cat])
            for xbin, ybin, zbin in self._binIter(cat):
                bin_name = self._binName(cat_name,xbin,ybin,zbin)
                content = self._binContent(cat_hist,xbin,ybin,zbin)
                if constant or self._nSurroundingZeros(cat_hist,xbin,ybin,zbin) > zero_thresh:
                    self.binVars[bin_name] = RooConstVar(bin_name, bin_name, content)
                else:
                    self.binVars[bin_name] = RooRealVar(bin_name, bin_name, max(5,content), 1e-6, 1e6)
                    self.nuisances.append({'name':bin_name, 'constraint':'flatParam', 'obj': self.binVars[bin_name]})
                self._varStorage.append(self.binVars[bin_name]) # For safety if we add shape templates

    def AddShapeTemplates(self,nuis_name,up_shape,down_shape,constraint="param 1 0"):
        '''Add variation shape templates that are used to create a map between
        a new nuisance parameter (named `nuis_name`) and the values for a given bin.
        To accomodate the potential for multiple shape templates, the new parameter
        will control the relative yield of the bin (ie. as a percentage). 

        This means for a nuisance value of 0, the multiplicative term on the bin yield will
        be 1. For nuisance value +1(-1), the multiplicative term on the bin yield will be
        the ratio of the bin value in `up_shape`(`down_shape`) to the nominal value.
        
        If `BinnedDistribution.forcePositive` is True, the parameters will extrapolate bin values above(below)
        nuisance values of +1(-1) using
        exponentials so that the values asymptotically approach 0. When `BinnedDistribution.forcePositive`
        is False, the values are exptrapolated linearly.

        For asymmetric uncertainties in a given nuisance `n`, the region defined by `n > -1` and `n < 1`
        is modeled using sigmoid functions which smoothly turn "on" and "off" the extrapolated pieces.
        This modeling provides a consistent description between -1 and 1, satisifies the boundary conditions
        at `n` of 0, 1, and -1, and is continuous in its first and second derivatives.

        Args:
            nuis_name (str): [description]
            up_shape (TH2): Input 2D histogram representing "up" variation.
            down_shape (TH2): Input 2D histogram representing "down" variation.
            constraint (str, optional): Can only be 'flatParam' or 'param <mu> <sigma>' (options in the Combine card) 
                which represent "no constraint" and "Gaussian constraint centered at <mu> and with width <sigma>", respectively.
                Defaults to "param 1 0".
            forcePositive (bool, optional): If True, shape template mapping will use exponentials so that values asymptotically
                approach zero as the associated nuisance increases/decreases. If False, the mapping will be linear.
        '''
        nuisance_par = RooRealVar(nuis_name,nuis_name,0,-5,5)
        self.nuisances.append({'name':nuis_name, 'constraint':constraint, 'obj': nuisance_par})

        for cat in self.subspaces:
            cat_name = self.name+'_'+cat
            cat_hist_up =   copy_hist_with_new_bins(up_shape.GetName()+'_'+cat,  'X', up_shape,   self.binning.xbinByCat[cat])
            cat_hist_down = copy_hist_with_new_bins(down_shape.GetName()+'_'+cat,'X', down_shape, self.binning.xbinByCat[cat])
            for xbin, ybin, zbin in self._binIter(cat):
                bin_name = '%s_%s_bin_%s'%(cat_name,nuis_name,self._binSuffix(xbin,ybin,zbin))
                upC   = self._binContent(cat_hist_up,xbin,ybin,zbin)
                downC = self._binContent(cat_hist_down,xbin,ybin,zbin)
                self.binVars[bin_name] = singleBinInterp( # change to singleBinInterpQuad to change interpolation method
                                            bin_name, self.getBinVar(xbin,ybin,c=cat,zbin=zbin), nuisance_par,
                                            upC,
                                            downC,
                                            self.forcePositive
                )
                self._varStorage.append(self.binVars[bin_name]) # For safety if we add more shape templates

    def KDESmooth(self):
        raise NotImplementedError()

    def _nSurroundingZeros(self,hist,xbin,ybin,zbin=None):
        '''Count how many cells of the surrounding neighborhood (3x3 in 2D, 3x3x3 in 3D),
        including the center, are <= 0. Returns 0 immediately if the center is > 0.
        '''
        if self.is3D and zbin is not None:
            if hist.GetBinContent(xbin,ybin,zbin) > 0:
                return 0
            nzeros = 0
            for x_,y_,z_ in itertools.product([xbin,xbin-1,xbin+1],[ybin-1,ybin,ybin+1],[zbin-1,zbin,zbin+1]):
                if hist.GetBinContent(x_,y_,z_) <= 0:
                    nzeros += 1
            return nzeros
        else:
            nzeros = 0
            if hist.GetBinContent(xbin,ybin) > 0:
                nzeros = 0
            else:
                for pair in itertools.product([xbin,xbin-1,xbin+1],[ybin-1,ybin,ybin+1]):
                    if hist.GetBinContent(pair[0],pair[1]) <= 0:
                        nzeros += 1
            return nzeros

def singleBinInterp(name, nuis, binVar, upVal, downVal, forcePositive):
    '''Create a RooFormulaVar containing the nuisance parameter that can
    morph the initial `binVar` value between the values of `upVal` and `downVal`.
    
    To accomodate the potential for multiple shape templates, the new parameter
    will control the relative yield of the bin (ie. as a percentage). 

    This means for a nuisance value of 0, the multiplicative term on the bin yield will
    be 1. For nuisance value +1(-1), the multiplicative term on the bin yield will be
    the ratio of the bin value in `up_shape`(`down_shape`) to the starting nominal value.
    
    If `forcePositive` is True, the parameters will extrapolate bin values above(below)
    nuisance values of +1(-1) using
    exponentials so that the values asymptotically approach 0. When `forcePositive`
    is False, the values are exptrapolated linearly.

    For asymmetric uncertainties in a given nuisance `n`, the region defined by `n > -1` and `n < 1`
    is modeled using sigmoid functions which smoothly turn "on" and "off" the extrapolated pieces.
    This modeling provides a consistent description between -1 and 1, satisifies the boundary conditions
    at `n` of 0, 1, and -1, and is continuous in its first and second derivatives.

    Args:
        name (str): Name for output RooFormulaVar.
        nuis (RooRealVar): Parameter to control yield changes across multiple bins.
        binVar (RooAbsReal): Current bin value represented as RooRealVar or RooFormulaVar (derives from RooAbsReal).
        upVal (float): Absolute "up" variation of the bin value.
        downVal (float): Absolute "down" variation of the bin value.
        forcePositive (bool): If True, shape template mapping will use exponentials so that values asymptotically
                approach zero as the associated nuisance increases/decreases. If False, the mapping will be linear.
    Returns:
        RooFormulaVar: New bin value which includes interpolation term.
    '''
    activate_pos = '(1/(1 + exp(-5x)))' # Use sigmoid for activation
    activate_neg = '(1/(1 + exp(5x)))'
    if forcePositive:
        pos_term = '({u}^@0)'.format(u=upVal)
        neg_term = '({d}^(-1*@0))'.format(d=downVal)
    else:
        pos_term = '(1+({u}-1)*@0)'.format(u=upVal)
        neg_term = '(1+(1-{d})*@0)'.format(d=downVal)

    full = '@1*({act_pos}*{pos}+{act_neg}*{neg})/{nom}'.format(act_pos=activate_pos, act_neg=activate_neg, pos=pos_term, neg=neg_term, nom=binVar.getValV())
    return RooFormulaVar(name, name, full, RooArgList(nuis,binVar))
