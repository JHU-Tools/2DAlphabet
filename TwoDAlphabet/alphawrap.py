from collections import OrderedDict
from TwoDAlphabet.helpers import roofit_form_to_TF1
from ROOT import RooRealVar, RooFormulaVar, RooArgList, RooParametricHist2D, RooConstVar, TFormula, RooAddition
from TwoDAlphabet.binning import copy_hist_with_new_bins
import itertools
# import numpy as np
# from numpy.lib.function_base import piecewise

_subspace = ['LOW','SIG','HIGH']
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
        nuisances (list): All tracked nuisance dictionaries.
        binVars (OrderedDict): Ordered dictionary of all RooAbsArgs representing the bins of the space.
        binArgLists (dict): Dict mapping of the subspaces (LOW, SIG, HIGH) to the RooArgList of the RooAbsArgs in the subspace.
        rph (dict): Dict mapping of the subspaces (LOW, SIG, HIGH) to the RooParametricHist2D objects of the subspaces.
        forcePositive (bool): Option to ensure bin values cannot be negative.
    '''
    
    def __init__(self,name,binning,forcePositive=True):
        '''Constructor.

        Args:
            name (str): Unique name of object which will be prepended to all associated RooFit objects.
            binning (TwoDAlphabet.Binning): Binning scheme object.
            forcePositive (bool, optional). Defaults to True in which case the bin values will be lower bound by 1e-9.
        '''
        self.name = name
        self.binning = binning
        self.nuisances = []
        self.binVars = OrderedDict()
        self.binArgLists = {c:None for c in _subspace}
        self.rph = {c:None for c in _subspace}
        self.forcePositive = forcePositive
        self._varStorage = [] # only used by AddShapeTemplates

    def _manipulate(self,name,other,operator=''):
        '''Base method to create a new Generic2D object. When combining
        `self` and `other`, a new set of RooFormulaVars will be created for
        the new Generic2D object that connect `self` and `other` with the
        `operator` string. The associated nuisances of `self` and `other` will
        also be passed to the new object as one set (with potential duplicates removed).
        
        If attempting to add, subtract, multiply, or divide,
        use the dedicated methods. More complex use cases could be built here.

        Args:
            name (str): Unique name for the new output Generic2D object.
            other (Generic2D): Object to combine with self.
            operator (str, optional): Connecting mathematical operator string for the combination. Defaults to '*'
                which causes the method to return self*other.

        Returns:
            Generic2D: Object containing the combination of `self` and `other`.
        '''
        out = Generic2D(name,self.binning,self.forcePositive)
        for cat in _subspace:
            new_cat_name = name+'_'+cat
            for ybin in range(1,len(self.binning.ybinList)):
                for xbin in range(1,len(self.binning.xbinByCat[cat])):
                    new_bin_name   = '%s_bin_%s-%s'%(new_cat_name,xbin,ybin)
                    self_bin_name  = new_bin_name.replace(new_cat_name, self.name+'_'+cat)
                    other_bin_name = new_bin_name.replace(new_cat_name, other.name+'_'+cat)
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
        '''Produce a RooParametricHist2D filled with this object's binVars.

        Returns:
            RooParametricHist2D: Output RooFit object to pass to Combine.
        '''
        out_rph = {}
        out_add = {}
        for cat in _subspace:
            cat_name = self.name+'_'+cat
            cat_hist = self.binning.CreateHist(cat_name+'_temp',cat)
            obj_name = '%s_%s'%(name if name != '' else self.name, cat)

            self.binArgLists[cat] = RooArgList()
            for ybin in range(1,len(self.binning.ybinList)):
                for xbin in range(1,len(self.binning.xbinByCat[cat])):
                    bin_name   = '%s_bin_%s-%s'%(cat_name,xbin,ybin)
                    self.binArgLists[cat].add(self.binVars[bin_name])

            out_rph[cat] = RooParametricHist2D(
                        obj_name, obj_name,
                        self.binning.xVars[cat],
                        self.binning.yVar,
                        self.binArgLists[cat], cat_hist
            )
            out_add[cat] = RooAddition(obj_name+'_norm',obj_name+'_norm',self.binArgLists[cat])
        return out_rph, out_add

    def getBinVal(self,xbin,ybin):
        '''Get bin value (for the current parameter values) in
        bin (xbin, ybin).

        Args:
            xbin (int): Indexed at 1 for ROOT compatibility.
            ybin (int): Indexed at 1 for ROOT compatibility.

        Returns:
            float: Current bin value.
        '''
        return self.getBinVar(xbin,ybin).getValV()

    def getBinVar(self,xbin,ybin,c=''):
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
        formula_name = '%s_bin_%s-%s'%(self.name+'_'+c,xbin,ybin)
        return self.binVars[formula_name]
            

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

        for cat in _subspace:
            cat_name = name+'_'+cat
            for ybin in range(1,len(self.binning.ybinList)):
                for xbin in range(1,len(self.binning.xbinByCat[cat])):
                    bin_name = '%s_bin_%s-%s'%(cat_name,xbin,ybin)
                    xConst,yConst = self.mappedBinCenter(xbin,ybin,cat)
                    if forcePositive: final_formula = "max(1e-9,%s)"%(self._replaceXY(xConst,yConst))
                    else:             final_formula = self._replaceXY(xConst,yConst)

                    self.binVars[bin_name] = RooFormulaVar(
                        bin_name, bin_name,
                        final_formula,
                        self.arglist
                    )

    def _replaceXY(self,x,y):
        '''Find and replace "x" and "y" in the input formula
        with this method's arguments (floats) which should
        correspond to bin centers.

        Args:
            x (float): Value along x axis to evaluate.
            y (float): Value along y axis to evaluate.

        Returns:
            str: Formula with "x" and "y" replaced by provided numerical values.
        '''
        f = self.formula.replace(' ','')
        f = f.replace('+x','+%s'%x).replace('+y','+%s'%y)
        f = f.replace('*x','*%s'%x).replace('*y','*%s'%y)
        f = f.replace('-x','-%s'%x).replace('-y','-%s'%y)
        f = f.replace('/x','/%s'%x).replace('/y','/%s'%y)
        f = f.replace('(x','(%s'%x).replace('(y','(%s'%y)
        return f

    def getNparams(self):
        '''Get the number of parameters in the formula (not counting "x" or "y").
        Converts the formula formating and uses ROOT's TFormula to count.

        Returns:
            int: Number of parameters in the fit (not counting "x" or "y").
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
    
    def mappedBinCenter(self,xbin,ybin,cat):
        '''Convert global x and y bin to values in center of bins
        where the axis has been mapped to range [0,1].

        Args:
            xbin (int): X axis bin number.
            ybin (int): Y axis bin number.

        Returns:
            tuple of floats: x and y values, respectively.
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

        for cat in _subspace:
            cat_name = name+'_'+cat
            cat_hist = copy_hist_with_new_bins(cat_name,'X',inhist,self.binning.xbinByCat[cat])
            for ybin in range(1,cat_hist.GetNbinsY()+1):
                for xbin in range(1,cat_hist.GetNbinsX()+1):
                    content = cat_hist.GetBinContent(xbin,ybin)
                    bin_name = '%s_bin_%s-%s'%(cat_name,xbin,ybin)
                    if(content<funcCeiling):
                        xConst,yConst = self.mappedBinCenter(xbin,ybin,cat)
                        if forcePositive: 
                            final_formula = "max(1e-9,%s)"%(self._replaceXY(xConst,yConst))
                        else:             
                            final_formula = self._replaceXY(xConst,yConst)

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
        for cat in _subspace:
            cat_name = name+'_'+cat
            cat_hist = copy_hist_with_new_bins(cat_name,'X',inhist,self.binning.xbinByCat[cat])
            for ybin in range(1,cat_hist.GetNbinsY()+1):
                for xbin in range(1,cat_hist.GetNbinsX()+1):
                    bin_name = '%s_bin_%s-%s'%(cat_name,xbin,ybin)
                    if constant or self._nSurroundingZeros(cat_hist,xbin,ybin) > 7:
                        self.binVars[bin_name] = RooConstVar(bin_name, bin_name, cat_hist.GetBinContent(xbin,ybin))
                    else:
                        self.binVars[bin_name] = RooRealVar(bin_name, bin_name, max(5,cat_hist.GetBinContent(xbin,ybin)), 1e-6, 1e6)
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

        for cat in _subspace:
            cat_name = self.name+'_'+cat
            cat_hist_up =   copy_hist_with_new_bins(up_shape.GetName()+'_'+cat,  'X', up_shape,   self.binning.xbinByCat[cat])
            cat_hist_down = copy_hist_with_new_bins(down_shape.GetName()+'_'+cat,'X', down_shape, self.binning.xbinByCat[cat])
            for ybin in range(1,cat_hist_up.GetNbinsY()+1):
                for xbin in range(1,cat_hist_up.GetNbinsX()+1):
                    bin_name = '%s_%s_bin_%s-%s'%(cat_name,nuis_name,xbin,ybin)
                    self.binVar[bin_name] = singleBinInterp( # change to singleBinInterpQuad to change interpolation method
                                                bin_name, self.getBinVar(xbin,ybin,cat), nuisance_par,
                                                cat_hist_up.GetBinContent(xbin,ybin),
                                                cat_hist_down.GetBinContent(xbin,ybin),
                                                self.forcePositive
                    )
                    self._varStorage.append(self.binVars[bin_name]) # For safety if we add more shape templates  

    def KDESmooth(self):
        raise NotImplementedError()

    def _nSurroundingZeros(self,hist,xbin,ybin):
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

# def singleBinInterpQuad(name, nuis, binVar, upVal, downVal): # NOT USED
#     nomVal = binVar.getValV()
#     a,b,c = solve_quad([(-1,downVal),(0,nomVal),(1,upVal)])
#     m_0, b_0 = linear_extrap((-1,downVal),[a,b,c])
#     m_1, b_1 = linear_extrap((1,upVal),[a,b,c])
#     pieces = [
#         '(@0 < -1)*(%s*@0+%s)'%(m_0,b_0),
#         '(@0 > -1 && @0 < 1)*(%s*@0**2+%s*@0+%s)'%(a,b,c),
#         '(@0 > 1)*(%s*@0+%s)'%(m_1,b_1)]
    
#     return RooFormulaVar(name, name, '@1*(({0})-{1})/{1}'.format('+'.join(pieces),nomVal), RooArgList(nuis,binVar))

# def singleBinInterpCombine(name, nuis, binVar, upVal, downVal): # NOT USED
#     nomVal = binVar.getValV()
#     pieces = '(abs(@0) < 1)*(0.125*@0*(pow(@0,2)*(3*pow(@0,2) - 10) + 15))-(@0 < -1)+(@0 > 1)'
#     alpha = '@1*( (@0/2) * (({up}-{down})+({up}+{down}-2*{nom})*({piecewise})) )/{nom}'.format(up=upVal,down=downVal,nom=nomVal,piecewise=pieces)
#     return RooFormulaVar(name,name,alpha,RooArgList(nuis,binVar))

# class quadratic:
#     def __init__(self, params):
#         self.a = params[0]
#         self.b = params[1]
#         self.c = params[2]
#     def y(self,x):
#         return self.a * x**2 + self.b * x + self.c
#     def slope(self,x):
#         return 2*self.a*x + self.b

# class line:
#     def __init__(self, slope, intercept):
#         self.m = slope
#         self.b = intercept
#     def y(self,x):
#         return self.m*x + self.b

# def solve_quad(threepts):
#     # A = [[x_1^2, x_1, 1],
#     #      [x_2^2, x_2, 1],
#     #      [x_3^2, x_3, 1]]
#     # X = [[a],[b],[c]]
#     # B = [[y_1], [y_2], [y_3]]
#     A, B = [], []
#     for pt in threepts:
#         x, y = pt[0], pt[1]
#         A.append([x**2, x, 1])
#         B.append([y])
#     A = np.array(A)
#     B = np.array(B)
#     out = [x[0] for x in np.linalg.solve(A,B)]
#     return out[0], out[1], out[2]

# def linear_extrap(point,quad_params):
#     q = quadratic(quad_params)
#     slope = q.slope(point[0])
#     const = -1*slope*point[0]+point[1]
#     return slope,const

# print (solve_quad([(0,0),(-1,1),(1,1)]))
