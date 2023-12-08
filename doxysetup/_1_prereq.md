Getting Started {#getting-started}
===============

\tableofcontents

# Installation {#install}

2D Alphabet can only be used in a CMSSW environment where the Higgs Analysis
Combine Tool must be installed. Please follow the instructions below to
checkout a 2D Alphabet-friendly version of Combine. These instructions are
based off of the [Combine documentation](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) 
for 102X setup. Please cross-check the instructions here with the official
instructions. The only difference should be in the cloned repository and lack
of branch change. Note also that the CMSSW release is different from the
release recommended by the Combine tool documentation in order to maintain
compatibility with Fermilab's LPC changing to SL7 by September 1st, 2020.

For `csh`:
```sh
    set SCRAM_ARCH=slc7_amd64_gcc700
    cmsrel CMSSW_10_6_14
    cd CMSSW_10_6_14/src
    cmsenv
    git clone https://github.com/lcorcodilos/2DAlphabet.git
    git clone https://github.com/lcorcodilos/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    bash <(curl -s https://raw.githubusercontent.com/lcorcodilos/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
    scram b clean; scram b -j 10
    cmsenv
```

For `bash`:
```bash
    export SCRAM_ARCH=slc7_amd64_gcc700
    cmsrel CMSSW_10_6_14
    cd CMSSW_10_6_14/src
    cmsenv
    git clone https://github.com/lcorcodilos/2DAlphabet.git
    git clone https://github.com/lcorcodilos/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit/
    curl -s https://raw.githubusercontent.com/lcorcodilos/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh | bash
    scram b clean; scram b -j 10
    cmsenv
```

Run `combine --help` and check it returns the help menu to confirm you've successfully setup combine.

Try to call RooParametricHist2D in interactive python if you're feeling
uneasy and to ensure everything is working. 

```python
import ROOT
r = RooParametricHist2D()
```

# JSON Configuration Files {#config}

The goal of the JSON configuration file is to have an easily understandable and configurable input
that allows the
2D Alphabet software to read and organize analysis files and histograms to its
liking while also giving the user the ability to easily configure values like
bin sizes and ranges without lots of command line options. This means that, while
the user is encouraged to add what they need, there are some keys and values
that must stay the same. These static strings are always in capital letters
to make them easy to distinguish. The six static sections are described below.

## Describing existing histograms {#config-prs}

2D Alphabet assumes that any analysis histogram can be described using three
descriptors:

- the physics "process" (ex. data, ttbar, signal, etc),
- the selection region (ex. signal region, control region, validation region, etc),
- the variation in a systematic uncertainty (ex. nominal/no variation, uncertainty varied up, uncertainty varied down, etc)

Of course, some regions have different physics processes than others and some
systematic uncertainties only apply for a subset of processes or regions. Thus,
the set of all histograms for the analysis is not simply the cross product of
all processes, regions, and systematic variations.

To account for this in the config, there needs to be meta information provided that describes
exactly the combinations that should be considered. Thus, as the three sections for these descriptors
are shown, you will see information cross section boundaries to inform the construction of
the final statistical model.

Let's start with `PROCESSES`.

### `PROCESSES`

In this section, the user can define as many processes as they need. At a minimum,
one must define exactly one data source and at least one signal. The data 
process name should be `"data_obs"` which is the name expected by Combine.
The signal can be named what you please. Background can likewise be defined
but are not necessary for 2D Alphabet to work. 

Please note two important things: 
1. If you plan to define a data-driven background (defined by a transfer function or similar parametric form),
   you do not include that background here. The config should only reference processes for which you intend
   to use histogram templates.
2. Combine always requires that there be an observation (data), a background
   estimate, and a signal sample. In addition to the above, that means that not
   including a background in the config necessitates that you define a parameteric
   background model later.

Each key in `PROCESSES` is the name of each process of interest and is paired with a sub-dictionary
of information about the process.
- `"TYPE"`: one of `"DATA"`, `"BKG"`, or `"SIGNAL"`
- `"TITLE"`: the "pretty" title to include in plot legends (support's ROOT's LaTeX formatting),
- `"COLOR"`: the ROOT color code to use for plotting (fill color for backgrounds, line color for signals),
- `"SCALE"`: a factor by which to scale the normalization of all templates of this process,
- `"LOC"`: the path to the histogram, ex. `"path/file.root:histo_name"`,
- `"SYSTEMATICS"`: a list of strings that correspond to a subset of the names of systematic variations
  in the `SYSTEMATICS` section;

The following is an example `PROCESSES` section:

```json
{
  "PROCESSES": {
    "data_obs": {
      "TYPE": "DATA",
      "TITLE": "Data",
      "COLOR": 1,
      "SCALE": 1.0,
      "LOC": "/to/my/rootfiles/selection_data_obs.root:hist2D_$region",
      "SYSTEMATICS":[]
    },
    "signal": {
      "TYPE": "SIGNAL",
      "TITLE": "m_{X} = 2000, m_{Y} = 800 (GeV)",
      "COLOR": 1,
      "SCALE": 1.0,
      "LOC":"/to/my/rootfiles/selection_signal.root:hist2D_$region",
      "SYSTEMATICS":["lumi","TriggerEff18","Pileup18","jer18","jes18","jmr18","jms18"]
    },
    "ttbar_18": {
      "TYPE": "BKG",
      "TITLE": "t#bar{t}",
      "COLOR": 2,
      "SCALE": 1.0,
      "LOC":"/to/my/rootfiles/selection_ttbar_18.root:hist2D_$region",
      "SYSTEMATICS":["lumi","ttbar_xsec","TriggerEff18","Pileup18","jer18","jes18","jmr18","jms18"],
    } 
  }
}
```

Most of this should seem self-explanatory except for the `"LOC"` values which have some strange strings.
Specifically, they use a reserved substitiuion string - `$region` - which is a stand-in to mean that 
the histogram names change per region and the region names (defined as the keys in the `REGIONS` section)
should be subsituted to find the actual histograms.

There are two more of these reserved substitution strings - `$process` and `$syst`. We'll see `$syst` later
but we can use them here as well to avoid the redundancy of the `PROCESSES` keys if we'd like. For example,

```python
"LOC":"/to/my/rootfiles/selection_ttbar_18.root:hist2D_$region",
```
can change to
```python
"LOC":"/to/my/rootfiles/selection_$process.root:hist2D_$region",
```
with no loss of information and a potential reduction in possible mistakes.

Now that you know this, I can explain the last option that exists for the sub-dictionaries - `"ALIAS"`.
By default, 2D Alphabet will use the process key (ie. "data_obs", "signal", "ttbar_18") for the substitution
to `$process`. If you'd prefer something different, you can specify an `"ALIAS"` which is used instead.
By specifying an `"ALIAS"`,
one can specify separate input and output names for the same process - the input (what appears in histogram paths)
will be the `"ALIAS"` while the output will be the outer key for the process name.

For example, I sometimes run on a subset of data (say, just 2018 data) and the file for that subset would
require the `$process` to be `Data_18`.
To avoid having to rename the file to include `data_obs` just to please Combine, the outer key
can be `data_obs` (so Combine is happy) and the `"ALIAS"` can be `Data_18`.

### `REGIONS` {#config-regions}
The `REGIONS` section is far simpler in structure and, now that the reserved substitution strings
have been introduced, perhaps a bit more intuitive. Let's start with a simple example building off
of the `PROCESSES` section example:

```json
{
  "REGIONS": {
    "SignalRegion_fail": {
      "PROCESSES": ["signal"],
      "BINNING": "default"
    },
    "SignalRegion_loose": {
      "PROCESSES": ["signal"],
      "BINNING": "default"
    },
    "SignalRegion_pass": {
      "PROCESSES": ["signal", "ttbar"],
      "BINNING": "default"
    }
  }
}
```
You should think of this section as a way of listing out the regions, what processes appear
in each, and the binning scheme you want. If it weren't for all of the meta information like
color, titles, etc plus the systematic variations, this section would describe the analysis
in its entirety!

As you can see, each region is named with a key and these keys are what will be used for the
`$region` substituions we saw in the `PROCESSES` `"LOC"` areas. However, they will only be used for the
substitution when there is a match with a process listed in the sub-dictionary! Looking at the example above,
there will only be one histogram loaded for `ttbar` - the one for region "SignalRegion_pass"!
Whereas the `signal` will have histograms for all three regions.
In this way, you can control what appears in each region.

The new concept introduced here is the `"BINNING"` key which names a binning scheme described in the
dedicated `"BINNING"` section described below. In this case, each region will be binned according to
what has been defined in the `"default"` scheme. However, multiple schemes can be defined and assigned
to different regions (for example, if you have a validation region which populates a different phase space or has better/worse statistics).

### `SYSTEMATICS` {#config-syst}

Finally, the `SYSTEMATICS` section of the config describes the how to account for variations
in parameters due to systematic uncertainties and connects back to the `PROCESSES` and `REGIONS`
sections with some of the tools we've learned already.

Note that the systematic uncertainties corresponding to each key are considered
entirely uncorrelated with each other in the final model!

The simplest type of variation is a symmetric, log-normal. A good example
is the Run 2 luminosity uncertinaty of 1.6% which can be defined with:
```json
{
  "SYSTEMATICS": {
    "lumi": {
      "VAL": 1.016
    }  
}
```
The `"VAL"` key can map to any floating point number but note that it is on a log scale
and should be a value greater than one.

An asymmetric, log-normal can be added by specifying up and down variations like so (+5%/-7%):
```json
{
  "SYSTEMATICS": {
    "lumi": {
      "VAL": 1.016
    },
    "asym": {
      "VALUP": 1.05,
      "VALDOWN": 0.93
    }
}
```

For both log-normal types, the values will be used directly in the Combine card so any rules
set by Combine must be followed.

The third type of variation is shape-based. In other words, the bins of the histogram are not
fully correlated as with the normalization uncertainties and the user will input up and down
variations of the shape that correspond to changing the underlying parameter by \$f\sigma\$f
standard deviation of uncertainty. Typically, \$f\sigma\$f is 1. Combine will map these
variations to a unit Gaussian penalty term and create functional forms of the bin values
that smoothly interpolate between and extrapolate beyond the templates.
For 2D Alphabet, a shape uncertainty may be implemented like the following:

```json
{
  "SYSTEMATICS": {
    "shape_syst": {
      "UP":   "/to/my/rootfiles/selection_$process.root:hist2D_$region__$syst_up",
      "DOWN": "/to/my/rootfiles/selection_$process.root:hist2D_$region__$syst_down",
      "SIGMA": 1.0
    }
}
```

Note that the `"UP"` and `"DOWN"` keys map to strings with a scheme to find the template histograms.
Specifically, `$process`, `$region`, and `$syst` are reserved substitution keys that will be replaced
by actual sets or process, region, and systematic variation just as before.

Note also that `"SIGMA"` currently maps to 1 since the templates are always expected to be the result of
one standard deviation variations of the underlying parameter.
There are cases where it may be useful to change the underlying
constraint on-the-fly. To accomplish this, the user may change `"SIGMA"` to a different value which will change
*the factor by which to scale the unit Gaussian constraint*. This is a convention chosen by Combine and
can be read about [here](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/settinguptheanalysis/#binned-shape-analysis).
This means, if one wants to treat the templates as 2 standard deviation variations, `"SIGMA"`
should be 1/2 - or 0.5.

Finally, note that in the above scheme, the user is saying that the
histograms can be found inside the same file as the nominal shape. One could likewise specify
```json
{
  "SYSTEMATICS": {
    "shape_syst": {
      "UP": "/to/my/rootfiles/selection_$process_$syst_up.root:hist2D_$region",
      "DOWN": "/to/my/rootfiles/selection_$process_$syst_down.root:hist2D_$region",
      "SIGMA": 1.0
    }
}
```
which means that there are different files for each variation but with the same name for the histogram
as the nominal variation. Again, we just want to describe the systematic variation
generally, not for a specific process or region. If your histograms don't follow such a pattern
or a pattern which is impossible to scheme, you are encouraged to rename them.

In addition to the above, the `SYSTEMATICS` section also makes use of `"ALIAS"`, just like the `PROCESSES` section. Modifying the above example,

```json
{
  "SYSTEMATICS": {
    "shape_syst_16": {
      "ALIAS": "ShapeSyst",
      "UP": "/to/my/rootfiles/selection_$process_$syst_up.root:hist2D_$region",
      "DOWN": "/to/my/rootfiles/selection_$process_$syst_down.root:hist2D_$region",
      "SIGMA": 1.0
    },
    "shape_syst_17": {
      "ALIAS": "ShapeSyst",
      "UP": "/to/my/rootfiles/selection_$process_$syst_up.root:hist2D_$region",
      "DOWN": "/to/my/rootfiles/selection_$process_$syst_down.root:hist2D_$region",
      "SIGMA": 1.0
    },
    "shape_syst_18": {
      "ALIAS": "ShapeSyst",
      "UP": "/to/my/rootfiles/selection_$process_$syst_up.root:hist2D_$region",
      "DOWN": "/to/my/rootfiles/selection_$process_$syst_down.root:hist2D_$region",
      "SIGMA": 1.0
    }
}
```
This will use the `"ALIAS"` as the substitution for `$syst` but name the uncertainty as the key
(`"shape_syst_16"`, `"shape_syst_17"`, `"shape_syst_18"`). Since the systematics are tied to the
processes in the `PROCESSES` section, one can split the systematic uncertainties to be uncorrelated
per process (per year in this case).

## `GLOBAL` {#config-global}
This section is designed to help users with large configuration file
names by allowing them to create JSON-wide variables. For example,
if all of your files are located in `/long/path/to/my/files/`, you 
can store this string in the GLOBAL dictionary with a custom key 
(let's say `path`). Now instead of having to write the full directory
path for every process and systematic, the user can just write `path`.
This simplifies the JSON and also has the standard advantages of using
variables over several instances of the same object.

This functionality works by searching all strings in the JSON for instances
of each key in `GLOBAL` and replacing the key with its corresponding dictionary value.
Note that it does not work for booleans or integers/floats.

The user must be careful they don't accidentally use strings in the JSON
that are identical to keys in `GLOBAL` so accidental substitutions don't happen.
This means keys in `GLOBAL` should be at least partially descriptive 
(single character keys would be a bad idea). 

This also means that many of the histogram path schemes can be standardized
for later simplification. Consider the following `GLOBAL` section:

```json
"GLOBAL": {
    "FILE": "selection_$process.root",
    "FILE_UP": "selection_$process_$syst_up.root",
    "FILE_DOWN": "selection_$process_$syst_down.root",
    "HIST": "hist2D_$region__nominal",
    "HIST_UP": "hist2D_$region__$syst_up",
    "HIST_DOWN": "hist2D_$region__$syst_down",
    "path": "/go/to/my/rootfiles"
}
```

Note that these are all variables defined by the **user**.
With this, the previous `SYSTEMATIC` section example can simply change to,

```json
{
  "SYSTEMATICS": {
    "shape_syst": {
      "UP": "path/FILE_UP:HIST",
      "DOWN": "path/FILE_DOWN:HIST",
      "SIGMA": 1.0
    }
}
```

## `BINNING` {#config-binning}
This dictionary is the opportunity to define the axis binning of the user's space.
The binning values are split into x and y axis definitions where the x-axis describes
the variable whose signal region is blinded. Note that 2D Alphabet can rebin
and reduce the ranges of the input axes for you. The binning in each axis cannot be
beyond the range of the input histogram (hopefully, this is obvious) but it can be a
subset of the input range. The number of bins are restricted to be equal to or less
than the input number of bins (hopefully, this is also obvious). Additionally, any newly
defined bin edges must line up with input bin edges (ie. there's no bin splitting).

The signal bounds only apply to the `X` axis and must exist within the `X` axis range
defined in the configuration file so that there are always three regions defined
with non-zero width. Only the signal region will be blinded during the fitting
and plotting. The options to perform this blinding are described in the [`OPTIONS` section](OPTIONS).

For either axis, there are two options for binning: constant width and variable width.
To use variable bins, provided a list of bin edges to the `BINS` key. To have 
2D Alphabet calculate the bin edges with constant bin width, use the `MIN`, `MAX`, and
`NBINS` keys (as you would in a `TH1` definition).

An example formatting of the `BINNING` section is provided here.

```json
"X":
    "NAME": <name of internal variable on the x-axis>, # string
    "TITLE": <title to appear in plots with this axis>, # string
    "MIN": <lower bound of x-axis>, # float
    "MAX": <upper bound of x-axis>, # float
    "NBINS": <number of x bins from MIN to MAX>, # int
    "SIGSTART": <lower bound of signal region of x-axis>, # int
    "SIGEND": <upper bound of signal region of x-axis> # int
"Y"`
    "NAME": name of your variable on the y-axis, # string
    "TITLE": title that you'd like to appear in plots with this axis, # string
    "BINS": [<edge1>, <edge2>,...<edgeN>] # [float]
```

Because the `X` axis is split into three portions (`LOW`, `SIG`, `HIGH`), one can also
define binning per region of the `X` axis. `BINS` can be used as well as `MIN`, `MAX`, `NBINS` like so:

```json
"X":
    "NAME": <name of internal variable on the x-axis>, # string
    "TITLE": <title to appear in plots with this axis>, # string
    "LOW": {
        "MIN": <lower bound of x-axis>, # float
        "MAX": <upper bound of x-axis>, # float
        "NBINS": <number of x bins from MIN to MAX> # int
    }
    "SIG": {
        "MIN": <lower bound of x-axis>, # float
        "MAX": <upper bound of x-axis>, # float
        "NBINS": <number of x bins from MIN to MAX> # int
    }
    "HIGH": {
        "MIN": <lower bound of x-axis>, # float
        "MAX": <upper bound of x-axis>, # float
        "NBINS": <number of x bins from MIN to MAX> # int
    }
    "SIGSTART": <lower bound of signal region of x-axis>, # int
    "SIGEND": <upper bound of signal region of x-axis> # int
```



## `OPTIONS` {#config-options}

This section is dedicated to providing per-config options to 2D Alphabet. 
Unless noted, these options only affect `run_MLfit.py` since
that tool is responsible for building the model performing
the initial fit that is read by the other `run_*.py` tools.

Required among these options the `name` for the configuration file
so that the pass-fail pair can be uniquely identified among those provided
to 2D Alphabet.
Truly optional options include those for blinding, plotting, and changing 
fit functionality. The most common possible options are described in detail below.

- **name** (`string`, default=`None`): Required. Provides a unique name to the configuration file.
- **tag** (`string`, default=`None`): Provides a tag to group configuration files and 
  create an output directory of the same name ("project directory"). Overriden by
  the `-q/--tag` command-line option to `run_MLfit.py` which will set the provided tag
  to all configuration files.
- **year** (`int`, default=1): Default `1` indicates this is full Run II. Other options
  are `16`, `17`, and `18`. This option determines what luminosity will be included in the plots
  and helps group configuration files when plotting the full Run II together (when years have been fit
  separately).
- **blindedFit** (`bool`, default=`true`): Blinds the fit to the signal region defined by `SIGSTART` and `SIGEND`.
- **blindedPlots** (`bool`, default=`true`): Blind the plots to the signal region defined by `SIGSTART` and `SIGEND`.
- **plotUncerts** (`bool`, default=`true`): Turns on plotting of the `X` and `Y` projections of the
  shape uncertainty templates provided (after rebinning). Useful for comparing nominal templates
  with uncertainty variations.
- **prerun** (`bool`, default=`true`): Runs the ML fit for this configuration file
  by itself and feeds the resulting transfer function parameters into the simultaneous
  fit with the other configuration files. This option can be effective if the fit
  cannot converge to a minimum because of too few steps.
- **plotPrefitSigInFitB** (`bool`, default=`false`): By default, the background-only
  fit result plots will include no signal since `r = 0`. This option set to `true` will
  plot the pre-fit signal (normalized to `r = 1`) instead. The s+b fit results always
  plot with the post-fit value of `r`.
- **plotTitles** (`bool`, default=`true`): By default, plot titles are always included in the plotting
  (this is useful for determining the boundaries of slices that are projected onto the opposite axis).
  Set to `false` this option turns off plot titles so that they are more appropriate for publication.
- **addSignals** (`bool`, default=`true`): By default, this option adds the signals for the sake of plotting so that
  they are considered as "one signal" in the legend and plots. This is useful if you have defined a separate signal for each
  year in a Run II configuration file. When set to `false`, the plotting will consider the signals
  as separate and they will be drawn separately and with different legend entries. This is useful
  if you have multiple unique signals.

The following are "developer/advanced options" that would mainly be useful for anyone
working on the 2D Alphabet code or debugging a fit.

- **draw** (`bool', default=`false`): This will live draw plots as they are created
  (sets `gROOT.SetBatch(kTRUE)`).
  Useful if one has pauses in the code to view the plots. Since 2D Alphabet most likely
  is running remotely, this will be very slow!
- **verbosity** (`int', default=`0`): Sets Combine verbosity. Setting in the first of
  the provided configuration files will determine behavior of the group run.
- **fitGuesses** ('dict`, default=`None`): Not recommended for use since it has gone
  untested for quite some time. This is because it is largely unnecessary and will
  most likely be removed. For the sake of documentation, the option creates slices
  along `Y` and projects them onto `X` to perform 1D fits in the slices. The parameters
  of the 1D fits are then fit as a function of `Y` to determine the `Y` dependence.
  The output of this pseudo-2D fit is then fed into the true 2D fit with Combine as
  a better pre-fit starting point.
- **overwrite** (`bool`, default=`false`): Explicitly deletes the configuration file's
  corresponding sub-folder in the project directory before running. This is usually 
  not necessary unless you've done something such as remove a MC-based background
  from the model without changing the tag. The plots of that process from the previous
  run would remain since no new plots would be made to replace them.
- **recycle** - (`[string]`, default=`[]`): Not recommended for use since the introduction
  of the `--recycleAll` option of `run_MLfit.py`. Recycles named pieces of saveOut.p. 
  Available pieces are "name", "tag", "xVarName", "yVarName", "xVarTitle", "yVarTitle",
  "sigStart", "sigEnd", "freezeFail", "blindedFit", "blindedPlots", "newXbins", "full_x_bins",
  "newYbins", "rpf", "rpfVarNames", "organizedDict", "floatingBins". I cannot think
  of a good use for this option that would not also be dangerous or overcome by just re-running!
     


# Python Interface {#python}

## Constructing TwoDAlphabet {#python-twoD}

## Defining Parametric Shapes {#python-parametric}