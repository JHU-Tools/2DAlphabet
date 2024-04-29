# Installation instructions
Instructions to set up 2DAlphabet on `el8` and `el9` architectures (tested on FNAL LPC). **NOTE:** 
* [Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#combine-v9-recommended-version) uses `CMSSW_11_3_X` which runs only on `slc7` and therefore has to be setup using apptainer on `el8/9`. 
* Once you have performed the setup via apptainer once, you *do not* have to repeat it each time. Simply go to the directory in which you performed this setup, type `cmssw-el7`, and navigate to the `CMSSW_11_3_4/src/` directory within the container to do your work

## One-time setup
```
# Set up Combine
cmssw-el7
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.2.1

# Set up Combine Harvester
cd ../../
git clone git@github.com:JHU-Tools/CombineHarvester.git
cd CombineHarvester/
git fetch origin
git checkout CMSSW_11_3_X
cd ../

# Build
scramv1 b clean
scramv1 b -j 8

# install 2DAlphabet
git clone https://github.com/JHU-Tools/2DAlphabet.git
python3 -m virtualenv twoD-env
source twoD-env/bin/activate
cd 2DAlphabet/
git fetch origin
git checkout el8-el9-py3
python setup.py develop
```

You can test that everything works by opening a python shell and typing:
```
import ROOT
r = ROOT.RooParametricHist()
```

## Setup after initial installation in apptainer

You do not need to re-compile the packages every time you log back in to the node. Instead:

* Navigate to the directory where you started the apptainer container
* Start the container via `cmssw-el7`
* `cd CMSSW_11_3_4/src/`
* `cmsenv`
* `source twoD-env/bin/activate`

And then you are ready to work with 2DAlphabet + Combine
