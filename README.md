# Installation instructions
```
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.1.0
cd ../../
git clone git@github.com:JHU-Tools/CombineHarvester.git
cd CombineHarvester/
git fetch origin
git checkout CMSSW_11_3_X
cd ../
scramv1 b clean
scramv1 b -j 8

git clone https://github.com/JHU-Tools/2DAlphabet.git
python3 -m virtualenv twoD-env
source twoD-env/bin/activate
cd 2DAlphabet/
git fetch origin
git checkout py3
python setup.py develop

```

