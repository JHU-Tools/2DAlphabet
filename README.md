# Installation instructions
Instructions to set up 2DAlphabet on `el8` and `el9` architectures (tested on FNAL LPC). **NOTE:** [Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#combine-v9-recommended-version) uses`CMSSW_11_3_X` runs on `slc7`, which has to be setup using apptainer on `el8/9`. 

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

