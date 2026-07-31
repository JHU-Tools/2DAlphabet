#!/usr/bin/env bash

_twod_previous_dir="${PWD}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "/home/cmsusr/${CMSSW_VERSION}/src"
eval "$(scramv1 runtime -sh)"
source /home/cmsusr/twoD-env/bin/activate
cd "${_twod_previous_dir}"

unset _twod_previous_dir
