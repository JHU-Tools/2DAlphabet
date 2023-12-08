# Installation instructions
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_14
cd CMSSW_10_6_14/src
cmsenv
git clone https://github.com/ammitra/2DAlphabet.git
git clone --branch 102x https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
curl -s https://raw.githubusercontent.com/lcorcodilos/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh | bash
scram b clean; scram b -j 8
cmsenv
python -m virtualenv twoD-env
source twoD-env/bin/activate
cd 2DAlphabet/
python setup.py develop
```

