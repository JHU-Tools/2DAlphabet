# Guides  {#guides}
\tableofcontents

# Running an ML fit {#mlfit}

This tool takes as input a list of JSON configuration
files (detailed in [configs](configs) via a command such as

```bash
    python run_MLfit.py input_config1.json input_config2.json ... input_configN.json
```

where `input_config*.json` are the configurations for each
Alphabet pair. Based on the input from the JSON files, 2D
Alphabet grabs all distributions, generates the RooFit objects
from them, creates the pass and fail distributions of the
QCD/non-resonant background estimated from data, passes these
all to Combine, calls Combine (`-M FitDiagnostics`), interprets the fit result, and
plots the post-fit distributions. It also outputs additional
plots and read-out for debugging or reference. The `run_MLfit.py`
script also creates the starting point for all other `run_`
scripts because it:
- Takes the model from pre-fit (where we know essentially nothing
about the total background estimate) to post-fit where background
estimate is meaningful,
- Creates a directory where all of the information of the run can
be stored and accessed for later use.

Thus, `run_MLfit.py` should be run before the other scripts whenever
attempting to fit a new Alphabet pair or a new set of pairs.

The main command line options are provided below with longer descriptions
along with the output of `python run_MLfit.py --help`.

- `-q, --tag` Assigns a tag for the run. The tag determines the 
naming of the folder (and some of the nested objects) the will be created
during processing. This folder is often referred to as the "project directory".
This option takes precedent
over tags defined in configuration files and is a good way to ignore the need
for configuration files to have the same tag and to make
sure you don't overwrite an old result on-the-fly.
- `--rMin/--rMax` These are the minimum and maximum bounds of the
signal strength (r) in the fit. They default to 0 and 5, respectively.
It may be useful to loosen or tighten the bounds if the fit is failing
near one of the boundaries.
- `--recycleAll` Recycles everything generated in the previous run
using the given tag. This means the construction of the workspace, rebinned histograms, and
other steps previous to the fit are skipped and the previous run versions
are loaded instead. Use this if you haven't changed your configuration
file and would like to speed things up.
- `--skipFit` Skips running the fit and goes directly to plotting.
- `--skipPlots` Skips running the plots. Sometimes this is the part
that takes the longest because a sampling method is used to estimate errors.

```
Options:
-h, --help            show this help message and exit
-q TAG, --tag=TAG     Assigns a tag for this run
-s V1=1.0,V2=1.0..., --setParameter=V1=1.0,V2=1.0...
                        String of parameters to set pre-fit. Uses same comma
                        separated format as Combine (V1=1.0,V2=1.0...)
--rMin=rMin           Minimum bound on r (signal strength)
--rMax=rMax           Minimum bound on r (signal strength)
--recycleAll          Recycle everything from the previous run with this tag
--skipFit             Skip fit and go directly to plotting (WARNING: Will
                        use previous fit result if it exists and crash
                        otherwise)
--skipPlots           Skip plotting
--fullRun2            Plot sum of years 16, 17, 18
--CL=CL               Command-line options to set for all configs
```

# Running asymptotic limits {#limit}

This tool is for calculating limits on the signal strength (`r`) of
a given simulated signal. The script takes the background only fit result
`fit_b` output from `run_MLfit.py`, morphs a copy of the pre-fit workspace
to this result (thus bypassing having to perform a likelihood fit again),
generates toys from the background estimate of the morphed pre-fit workspace
to create pseudo-data, and fits these toys to derive expected limits for the
simulated signal (blinding can be turned off so that the real data is also
fit to get an observed limit on `r`). The actual pseudo-data generation
and asymptotic limit calculation is handled by the AsymptoticLimit method
of Combine.

Because the background only fit result is used, the signal simulation used
in the configuration file for the `run_MLfit.py` does not matter. However,
the signal for limit setting obviously matters and in all likelihood
(no pun intended), the user will want to scan over several simulated signals to
generate limits as a function of their samples (perhaps each one was generated
at a different mass). Note that the model is rebuilt each time a signal 
(and the associated uncertainty templates) needs to be changed.

You do NOT have to write a configuration file for each simulated sample. At
any point in the python call to `run_Limit.py`, one can provide a string
swap specified by the syntax `oldString:newString` which will replace all
instances of `oldString` in the configuration files provided with `newString`.
In fact this can be used in `run_MLfit.py` as well.

Running `python run_Limit.py --help` returns the following:

```
Options:
-h, --help            show this help message and exit
-q <tag>, --tag=<tag>
                        Assigns a tag for this run
-d <dir>, --projDir=<dir>
                        Points to the directory where the b-only fit result is
                        located
--unblindData         Unblind the observation and calculate the observed
                        limit
--recycleAll          Recycle everything from the previous run with this
                        tag. Note that this does not allow for string
                        substitution.
```

# Running statistical tests {#stats}

The `run_Stats.py` tool hosts the infrastructure to run the majority
of necessary statistical tests on the model and fit. Each call
to `run_Stats.py` can only be run with a single "mode". The test
must be run on a project directory.
The available modes are represented by the following options.

- `--gof` Runs the goodness of fit test for the model in the 
provided `--projDir`. If the fit in `--projDir` is blinded, the 
goodness of fit test will be blinded as well. The saturated test
statistic is used. The KS and AD test statistics are not available
because they rely on CDFs which are not well defined for two dimensional
distributions.
- `--signalInjection <r>` Injects the designated amount of signal
`<r>` on top of the background-only model from `--projDir`.
- `--ftest <option>` Runs a Fischer test comparing the model in `--projDir`
against a model with more parameters, `--altDir`. The process of running
the test is split into several options so that pieces can be recycled 
between comparisons:
- **generate** Generates the toys `--projDir`. The `--altDir` does not
need to be specified.
- **fitMain** Fits the base/main model to the toys generated from itself
and calculates the saturated test statistic for each.
- **fitAlt** Fits the alternate model to the toys generated from the 
base model and calculates the saturated test statistic for each.
- **post** Collects the outputs and fitAlt and fitMain, plots, and calculates the 
pvalue when comparing the fits to data against the fits to the toys.
- **pvalue** Skips toys entirely and assumes toys follow an F-distribution.
The pvalue is calculated when comparing data against the F-distribution. 
This option is appropriate if you've first confirmed that the toys follow the
F-distribution) which is plotted in the `post` option. 

@note Statistical tests require that toys (pseudo-data) be generated and run. 
The number of toys to use can be specified with the  `-t/--toys` option.
Typically, 500 toys are a good number to generate and analyze.
The intricacies of the toy generation are handled by 2D Alphabet but 
in general, the toys are generated from a post-fit model (`fit_b` or a 
modified version depending on the "mode") and fit with the pre-fit model.
The seed for the toys can be changed with the `--seed` option.

@note For `run_Stats.py` specifically, there is the `--condor` option which
will run the Combine commands on condor batch nodes. When `--condor` 
is used, the jobs should be monitored by the user. Once the jobs are
finished, the same command should be run but with `--condor` replaced
by `--post` which will collect the returned tarballs, untar them 
temporarily, and generate the plots. 

The main advantage of using `--condor` is to split the toys among separate
jobs with the `--toyJobs` option which specifies into how many jobs to split
the toys. For example `-t 500 --toyJobs 50` would put 10 toys into each job.
The commands run in each job will have different seeds so that the toys
are not identical.

Note that even though the job outputs remain untarred (even after running
`--post`), the folders are not small. Depending on the model, they could be
between 1-2 GB. This can eat up space quickly on disk so be cognizant of how
many old tests you are saving!

Finally, the 2D Alphabet environment is already tarred and stored on
EOS for access to the jobs. If you make changes, you'll need to create
a new tarball of the environment and put it on EOS for the condor nodes
to access.

The `CondorHelper.py` can be used to submit jobs for other tools such as 
`run_Limit.py`. See [Condor Helper](condorhelper) for more information.

Running `python run_Stats.py --help` returns the following:

```
Options:
-h, --help            show this help message and exit
-d <dir>, --projDir=<dir>
                        Home of the project - has the cards, fit results, etc
-a <dir>, --altDir=<dir>
                        Home of the alternative model that you'd like to
                        compare against the one in projDir
-t <N>, --toys=<N>    Number of toys to generate - fed to combine
--toyJobs=<N>         Number of jobs to split toy fitting into - fed to
                        condor job building
--seed=<N>            Seed - fed to combine
--rMin=<rMin>         Minimum bound on r (signal strength)
--rMax=<rMax>         Minimum bound on r (signal strength)
-w <name>             Model workspace (from text2workspace)
--plotOnly            Only plot
--dryrun              Dry run the combine commands to console
--gof                 Perform goodness of fit test
--signalInjection=<r>
                        Perform signal injection test
--ftest=<option>      Perform F test. Options are 'generate', 'fitAlt',
                        'fitMain', 'post', 'pvalue'
--post                Run in conjunction with diagnosticsWithToys or
                        signalInjection and condor jobs to process output
                        files
--condor              Submit a condor job (only for signal injection and
                        diagnosticsWithToys currently
--skipSnapshot        Use if you've already made a snapshot that you trust
                        for this project directory
```
    
