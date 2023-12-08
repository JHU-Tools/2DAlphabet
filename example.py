from time import time
from TwoDAlphabet import plot
from TwoDAlphabet.twoDalphabet import MakeCard, TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction
from TwoDAlphabet.helpers import make_env_tarball
import os

'''--------------------------Helper functions---------------------------'''
def _get_other_region_names(pass_reg_name):
    return pass_reg_name, pass_reg_name.replace('pass','loose'),pass_reg_name.replace('pass','fail')

def _select_signal(row, args):
    '''Used by the Ledger.select() method to create a subset of a Ledger.
    This function provides the logic to determine which entries/rows of the Ledger
    to keep for the subset. The first argument should always be the row to process.
    The arguments that follow will be the other arguments of Ledger.select().
    This function should ALWAYS return a bool that signals whether to keep (True)
    or drop (False) the row.

    To check if entries in the Ledger pass, we can access a given row's 
    column value via attributes which are named after the columns (ex. row.process
    gets the "process" column). One can also access them as keys (ex. row["process"]).

    In this example, we want to select for signals that have a specific string
    in their name ("process"). Thus, the first element of `args` contains the string
    we want to find.

    We also want to pick a TF to use so the second element of `args` contains a
    string to specify the Background_args[1] process we want to use.

    Args:
        row (pandas.Series): The row to evaluate.
        args (list): Arguments to pass in for the evaluation.

    Returns:
        Bool: True if keeping the row, False if dropping.
    '''
    signame = args[0]
    poly_order = args[1]
    if row.process_type == 'SIGNAL':
        if signame in row.process:
            return True
        else:
            return False
    elif 'Background_' in row.process:
        if row.process == 'Background_'+poly_order:
            return True
        else:
            return False
    else:
        return True

def _load_CR_rpf(poly_order):
    twoD_CRonly = TwoDAlphabet('XHYfits_CR','example_config.json', loadPrevious=True)
    params_to_set = twoD_CRonly.GetParamsOnMatch('rpf.*'+poly_order, 'MX_2000_MY_800_area', 'b')
    return {k:v['val'] for k,v in params_to_set.items()}

def _load_CR_rpf_as_SR(poly_order):
    params_to_set = {}
    for k,v in _load_CR_rpf(poly_order).items():
        params_to_set[k.replace('CR','SR')] = v
    return params_to_set

def _generate_constraints(nparams):
    out = {}
    for i in range(nparams):
        if i == 0:
            out[i] = {"MIN":0,"MAX":1}
        else:
            out[i] = {"MIN":-5,"MAX":5}
    return out

_rpf_options = {
    '0x0': {
        'form': '0.1*(@0)',
        'constraints': _generate_constraints(1)
    },
    '1x0': {
        'form': '0.1*(@0+@1*x)',
        'constraints': _generate_constraints(2)
    },
    '0x1': {
        'form': '0.1*(@0+@1*y)',
        'constraints': _generate_constraints(2)
    },
    '1x1': {
        'form': '0.1*(@0+@1*x)*(1+@2*y)',
        'constraints': _generate_constraints(3)
    }
}

'''---------------Primary functions---------------------------'''
def test_make(SRorCR):
    '''Constructs the workspace for either the CR or SR (a different function
    could build them simultanesouly but in this example, we don't care to fit
    the two simultanesouly so separate treatment is fine).

    Args:
        SRorCR (str): 'SR' or 'CR'.
    '''
    assert SRorCR in ['SR', 'CR'] # quick sanity check

    fr = {}
    if SRorCR == 'SR': # example_config.json is setup for the CR so setup a find-replace dictionary to make the SR version if needed.
        fr = {'CR_loose':'SR_loose', 'CR_fail': 'SR_fail', 'CR_pass':'SR_pass'}

    # Create the twoD object which starts by reading the JSON config and input arguments to
    # grab input simulation and data histograms, rebin them if needed, and save them all
    # in one place (organized_hists.root). The modified JSON config (with find-replaces applied, etc)
    # is also saved as runConfig.json. This means, if you want to share your analysis with
    # someone, they can grab everything they need from this one spot - no need to have access to
    # the original files! (Note though that you'd have to change the config to point to organized_hists.root).
    twoD = TwoDAlphabet('XHYfits_'+SRorCR, 'example_config.json', findreplace=fr, loadPrevious=False)
    qcd_hists = twoD.InitQCDHists() # Create the data - BKGs histograms

    # This loop will only run once since the only regions are CR_pass, CR_fail, CR_loose (or the SR versions)
    # but it future proofs if we want to add more later for simultaneous fits.
    for p, l, f in [_get_other_region_names(r) for r in twoD.ledger.GetRegions() if 'pass' in r]:
        # Gets the Binning object and some meta information (stored in `_`) that we don't care about
        # The Binning object is needed for constructing the Alphabet objects.
        # If one wanted to be very robust, they could get the binning for `p` and `l` as well and check the binning
        # is consistent between the three.
        binning_f, _ = twoD.GetBinningFor(f)

        # Next we construct the Alphabet objects which all inherit from the Generic2D class.
        # This class constructs and stores RooAbsArg objects (RooRealVar, RooFormulaVar, etc)
        # which represent each bin in the space.

        # First we make a BinnedDistribution which is a collection of RooRealVars built from a starting
        # histogram (`qcd_hists[f]`). These can be set to be constants but, if not, they become free floating
        # parameters in the fit. 
        fail_name = 'Background_'+f
        qcd_f = BinnedDistribution(
                    fail_name, qcd_hists[f],
                    binning_f, constant=False
                )

        # We'll then book a flat TF which will be used to transfer between loose and pass
        # We keep it out of the loop below though because this will keep the same form
        # while the fail-to-loose TF changes with the different possible options.
        qcd_rpfT = ParametricFunction(
                        fail_name.replace('fail','rpfT'),
                        binning_f, '0.1*(@0)',
                        constraints={0:{"MIN":0, "MAX": 1}}
                    )

        # We add it to `twoD` so its included when making the RooWorkspace and ledger.
        # We specify the name of the process, the region it lives in, and the object itself.
        # The process is assumed to be a background and colored yellow but this can be changed
        # with optional arguments.
        twoD.AddAlphaObj('Background',f,qcd_f)

        # As a global variables, we've defined some different transfer function (TF) options.
        # We only want to include one of these at the time of fitting but we want to construct
        # them all right now so we can pick and choose later.
        for opt_name, opt in _rpf_options.items():
            # We have two regions determined by a TF, "pass" and "loose" with the "pass"
            # being a simple flat scaling of the loose. The functional form and the 
            # dictionary of constraints is defined in _rpf_options so we just plug
            # these in, being careful to name the objects uniquely (this affects
            # the naming of the RooFormulaVars created, which need to be unique).

            # The ParametricFunction class is the same as the BinnedDistribution except
            # the bins are RooFormulaVars constructed from the input formula with the
            # "x" and "y" taken as the centers of each bin.
            # The constraints option takes as input a dictionary with keys that control
            # the minimum, maximum, and error (initial step) of each parameter. It can
            # also be used to specify if the parameter should be unconstrainted (flatParam)
            # or Gaussian constrained (param <mu> <sigma>).
            qcd_rpfL = ParametricFunction(
                        fail_name.replace('fail','rpfL')+'_'+opt_name,
                        binning_f, opt['form'],
                        constraints=opt['constraints']
                    )
            
            # Of course, what we actually need is these TFs multiplied by something else:
            #     qcd_l = qcd_f*rpfL
            #     qcd_p = qcd_l*rpfT
            # The Multiply method will make a new set of RooFormulaVars defined by multiplying the RooAbsArgs
            # of each object together. Other methods exist for adding and dividing (where Add() can take a optional factor
            # so that subtraction is possible).
            qcd_l = qcd_f.Multiply(fail_name.replace('fail','loose')+'_'+opt_name, qcd_rpfL)
            qcd_p = qcd_l.Multiply(fail_name.replace('fail','pass')+'_'+opt_name, qcd_rpfT)

            # Now add the final models to the `twoD` object for tracking
            # Note that we have unique process names so they are identifiable
            # but we give them different titles so that they look pretty in
            # the final plot legends.
            twoD.AddAlphaObj('Background_'+opt_name,l,qcd_l,title='Background')
            twoD.AddAlphaObj('Background_'+opt_name,p,qcd_p,title='Background')

    # Save() will save the RooWorkspace and the ledgers and other associated pieces
    # so the twoD object can be reconstructed later. If this line doesn't run or
    # if something in the above needs to change, everything will need to be re-run to this point.
    twoD.Save()

def test_fit(SRorCR):
    '''Loads a TwoDAlphabet object from an existing project area, selects
    a subset of objects to run over (a specific signal and TF), makes a sub-directory
    to store the information, and runs the fit in that sub-directory. To make clear
    when a directory/area is being specified vs when a signal is being selected,
    I've redundantly prepended the "subtag" argument with "_area".
    '''
    assert SRorCR == 'CR' # Setup for either SR or CR but don't want to unblind accidentally until ready

    # So that the find-replace in the config doesn't need to be done again if I want
    # the SR (since it would have been performed already by test_make()), I grab
    # the runConfig.json that's already been saved in the created directory.
    working_area = 'XHYfits_'+SRorCR
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)
    # Access the Ledger and perform a selection on it to create a subset
    # from which to build the card. One can modify the Ledger DataFrames
    # manually to do more sophisticated manipulations but the select()
    # method will not modify the Ledger in-place. It always generates a new Ledger
    # which, by itself, is only stored in memory.

    # Create a subset of the primary ledger using the select() method.
    # The select() method takes as a function as its first argument
    # and any args to pass to that function as the remiaining arguments
    # to select(). See _select_signal for how to construct the function.
    subset = twoD.ledger.select(_select_signal, 'MX_2000_MY_800', '0x0')

    # Make card reads the ledger and creates a Combine card from it.
    # The second argument specifices the sub-directory to save the card in.
    # MakeCard() will also save the corresponding Ledger DataFrames as csvs
    # in the sub-directory for later reference/debugging. By default, MakeCard()
    # will reference the base.root workspace in the first level of the project directory
    # (../ relative to the card). However, one can specify another path if a different
    # workspace is desired. Additionally, a different dataset can be supplied via
    # toyData but this requires supplying almost the full Combine card line and
    # is reserved for quick hacks by those who are familiar with Combine cards.
    twoD.MakeCard(subset, 'MX_2000_MY_800_area')

    # Run the fit! Will run in the area specified by the `subtag` (ie. sub-directory) argument
    # and use the card in that area. Via the cardOrW argument, a different card or workspace can be
    # supplied (passed to the -d option of Combine). 
    twoD.MLfit('MX_2000_MY_800_area',rMin=-1,rMax=20,verbosity=0)

def test_plot(SRorCR):
    '''Load the twoD object again and run standard plots for a specific subtag.
    Assumes loading the Ledger in this sub-directory but a different one can
    be provided if desired.
    '''
    assert SRorCR == 'CR' # Setup for either SR or CR but don't want to unblind accidentally until ready
    working_area = 'XHYfits_'+SRorCR
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, 'MX_2000_MY_800', '0x0')
    twoD.StdPlots('MX_2000_MY_800_area', subset)

def test_limit(SRorCR):
    '''Perform a blinded limit. To be blinded, the Combine algorithm (via option `--run blind`)
    will create an Asimov toy dataset from the pre-fit model. Since the TF parameters are meaningless
    in our true "pre-fit", we need to load in the parameter values from a different fit so we have
    something reasonable to create the Asimov toy. 
    '''
    poly_order = '0x0'
    # Returns a dictionary of the TF parameters with the names as keys and the post-fit values as dict values.
    params_to_set = _load_CR_rpf_as_SR(poly_order) if SRorCR == 'SR' else _load_CR_rpf()
    working_area = 'XHYfits_'+SRorCR
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)
    
    # The iterWorkspaceObjs attribute stores the key-value pairs in the JSON config 
    # where the value is a list. This allows for later access like here so the user
    # can loop over the list values without worrying if the config has changed over time
    # (necessitating remembering that it changed and having to hard-code the list here).
    print ('Possible signals: %s'%twoD.iterWorkspaceObjs['SIGNAME'])
    for signame in twoD.iterWorkspaceObjs['SIGNAME']:
        # signame is going too look like <what we want>_18 so drop the last three characters
        signame = signame[:-3]
        print ('Performing limit for %s'%signame)

        # Make a subset and card as in test_fit()
        subset = twoD.ledger.select(_select_signal, signame, poly_order)
        twoD.MakeCard(subset, signame+'_area')
        # Run the blinded limit with our dictionary of TF parameters
        twoD.Limit(
            subtag=signame+'_area',
            blindData=True,
            verbosity=0,
            setParams=params_to_set,
            condor=False
        )

def test_GoF(SRorCR):
    '''Perform a Goodness of Fit test using an existing working area.
    Requires using data so SRorCR is enforced to be 'CR' to avoid accidental unblinding.
    '''
    assert SRorCR == 'CR' # Setup for either SR or CR but don't want to unblind accidentally until ready

    poly_order = '0x0'
    working_area = 'XHYfits_'+SRorCR
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)

    # If the card doesn't exist, make it (in the case that test_fit() wasn't run first).
    signame = 'MX_2000_MY_800'
    if not os.path.exists(twoD.tag+'/'+signame+'_area/card.txt'):
        subset = twoD.ledger.select(_select_signal, signame, poly_order)
        twoD.MakeCard(subset, signame+'_area')
    # Run the Goodness of fit test with 500 toys, r frozen to 0, TF parameters set to prefit.
    # This method always runs the evaluation on data interactively but the toy generation and evaluation
    # can be sent to condor with condor=True and split over several jobs with njobs=<int>.
    # Note that running a GoF test without data is relatively meaningless so by using this method,
    # you must unblind data. If you wish to use a toy dataset instead, you should set that
    # up when making the card.
    twoD.GoodnessOfFit(
        signame+'_area', ntoys=500, freezeSignal=0,
        condor=True, njobs=10
    )

    # Note that no plotting is done here since one needs to wait for the condor jobs to finish first.
    # See test_GoF_plot() for plotting (which will also collect the outputs from the jobs).

def test_SigInj(SRorCR):
    '''Perform a signal injection test'''
    assert SRorCR in ['SR','CR']

    poly_order = '0x0'
    signame = 'MX_2000_MY_800'
    working_area = 'XHYfits_'+SRorCR
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)

    # If the card doesn't exist, make it (in the case that test_fit() wasn't run first).
    if not os.path.exists(twoD.tag+'/'+signame+'_area/card.txt'):
        subset = twoD.ledger.select(_select_signal, signame, poly_order)
        twoD.MakeCard(subset, signame+'_area')

    # Perform the signal injection test with r=0 and with 500 toys split over 10 jobs on condor.
    # Because the data is blinded, we feed in the parameters from a previous fit so that we
    # have a model from which to generate toys.
    twoD.SignalInjection(
        signame+'_area', injectAmount=0,
        ntoys=500,
        blindData=True,
        setParams=_load_CR_rpf_as_SR(poly_order),
        condor=True, njobs=10)

def test_GoF_plot(SRorCR):
    '''Plot the GoF in XHYfits_<SRorCR>/MX_2000_MY_800_area (condor=True indicates that condor jobs need to be unpacked)'''
    plot.plot_gof('XHYfits_'+SRorCR,'MX_2000_MY_800_area', condor=True)

def test_SigInj_plot(SRorCR):
    '''Plot the signal injection test for r=0 injected and stored in XHYfits_<SRorCR>/MX_2000_MY_800_area
    (condor=True indicates that condor jobs need to be unpacked)'''
    plot.plot_signalInjection('XHYfits_'+SRorCR,'MX_2000_MY_800_area', injectedAmount=0,condor=True)

def test_Impacts(SRorCR):
    '''Calculate the nuisance parameter impacts. The parameters corresponding to the unconstrained bins
    of the fail region are ignored. Assumes that a fit has already been performed so that the post-fit
    uncertainties can be used for the scans. However, another card or workspace can be specified as well
    as a dictionary of parameters to set before running (setParams). With blindData=True, a pre-fit Asimov
    toy is generated for the sake of performing the scans. Since we're only using the CR here, blindData
    is set to False.
    '''
    assert SRorCR == 'SR' # Setup for either SR or CR but don't want to unblind accidentally until ready
    working_area = 'XHYfits_'+SRorCR
    twoD = TwoDAlphabet(working_area, '%s/runConfig.json'%working_area, loadPrevious=True)

    # We need to run impacts in the SR for them to make sense but we can't use the data in the SR while blinded.
    # So we need a toy to play with instead.
    poly_order = '0x0'
    subset = twoD.ledger.select(_select_signal, 'MX_2000_MY_800',poly_order)
    # Make a new area to play in
    twoD.MakeCard(subset, 'MX_2000_MY_800_impactArea')

    # Generate the toy
    toy_file_path = twoD.GenerateToys(
        'impactToy', 'MX_2000_MY_800_impactArea',
        card='card.txt', 
        workspace=None, 
        ntoys=1, seed=123456, expectSignal=0,
        setParams=_load_CR_rpf_as_SR(poly_order)
    )
    # Run the parameter impacts on the toy with the pre-fit workspace/card
    twoD.Impacts(
        'MX_2000_MY_800_impactArea',
        cardOrW='card.txt',
        extra='-t 1 --toysFile %s'%toy_file_path.split('/')[-1]
    )

def test_generate_for_SR():
    '''NOTE: This is an expert-level manipulation that requires understanding the underlying Combine
    commands. Use and change it only if you understand what each step is doing.
    
    Use the CR fit result to generate and fit a toy in the SR (without looking at SR data).
    There are two ways to do this which will be broken up into toyArea1 and toyArea2.'''
    # Load in the SR TwoDAlphabet object
    twoD = TwoDAlphabet('XHYfits_SR', 'XHYfits_SR/runConfig.json', loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, 'MX_2000_MY_800','0x0')
    params_to_set = _load_CR_rpf_as_SR('0x0')

    ###################################
    #-------- Version 1 --------------#
    ###################################
    # We'll make a card for each version to ensure the directory structure is made - they will be identical though to start.
    twoD.MakeCard(subset, 'MX_2000_MY_800_toyArea1')

    # Perform a fit as normal but via the `extra` arg, provide some commands
    # directly to combine to generate 1 toy, with seed 123456, and with r=0.
    # Note that --expectSignal 0 will generate with r=0 AND fit with r=0.
    twoD.MLfit(
        subtag='MX_2000_MY_800_toyArea1',
        setParams=params_to_set,
        rMin=0,rMax=5,verbosity=0,
        extra='-t 1 -s 123456 --expectSignal 0'
    )
    # Plot!
    twoD.StdPlots('MX_2000_MY_800_toyArea1',ledger=subset)

    ###################################
    #-------- Version 2 --------------#
    ###################################
    # We'll make a card for each version to ensure the directory structure is made - they will be identical though to start.
    twoD.MakeCard(subset, 'MX_2000_MY_800_toyArea2')
    # First generate a toy by itself. This means we can set r for *just* this step
    # as opposed to Version 1 where r was set for generation and for fitting.
    # Note that this method will generate frequentist toys but always skip the frequentist fit.
    # So if you'd like to generate from a post-fit workspace, you should fit first
    # and then provide a workspace snapshot
    toy_file_path = twoD.GenerateToys(
        'toys', 'MX_2000_MY_800_toyArea2',
        card='card.txt', workspace=None, # A card or workspace MUST be defined manually or one of these options should be set to True to use a default.
        ntoys=1, seed=123456, expectSignal=0,
        setParams=params_to_set
    )
    # Perform a fit as normal but via the `extra` arg, tell combine
    # to access our already-generated toy.
    # Note that r is now freely floating in this fit again.
    twoD.MLfit(
        subtag='MX_2000_MY_800_toyArea2',
        setParams=params_to_set,
        rMin=0,rMax=5,verbosity=0,
        extra='-t 1 --toysFile=%s'%toy_file_path.split('/')[-1]
    )
    # Plot!
    twoD.StdPlots('MX_2000_MY_800_toyArea2',ledger=subset)

if __name__ == '__main__':
    # Provided for convenience is this function which will package the current CMSSW and store it on the user's EOS (assumes FNAL).
    # This only needs to be run once unless you fundamentally change your working environment.
    # make_env_tarball()
    
    test_make('CR')
    test_make('SR')
    
    test_fit('CR')
    test_plot('CR')

    test_limit('SR')

    test_GoF('CR')
    test_SigInj('SR')

    test_Impacts('SR')
    test_generate_for_SR()

    # Run after condor jobs finish
    # test_GoF_plot('CR')
    # test_SigInj_plot('SR')