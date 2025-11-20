# 2DAlphabet 
2DAlphabet - A framework for performing 2D binned-likelihood fits with one background source derived from a parametric transfer function.
# Installation instructions
The `master` branch of 2DAlphabet is built for el9, corresponding to the following clusters:
* FNAL LPC: `cmslpc-el9.fnal.gov`
* CERN LXPLUS: `lxplus9.cern.ch`
```
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.0.1
cd ../../
git clone --branch CMSSW_14_1_0_pre4 git@github.com:JHU-Tools/CombineHarvester.git
scramv1 b clean
scramv1 b -j 16
git clone git@github.com:JHU-Tools/2DAlphabet.git
python3 -m virtualenv twoD-env
source twoD-env/bin/activate
cd 2DAlphabet/
python setup.py develop
```
