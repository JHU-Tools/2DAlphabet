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

## CVMFS container image

A prebuilt image containing CMSSW, Combine, CombineHarvester, and 2DAlphabet is
available through CVMFS on supported clusters such as CERN LXPLUS and FNAL
LPC:

```text
/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/jhu-tools/2dalphabet:v2.0
```

### Interactive shell

Mount the current directory at `/work` and enter the container:

```bash
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/jhu-tools/2dalphabet:v2.0

singularity shell \
    --bind /cvmfs \
    --bind "$(pwd -P):/work" \
    --pwd /work \
    "$IMAGE"
```

Inside the container, activate 2DAlphabet and enter the working directory:

```bash
source /home/cmsusr/setup_2dalphabet.sh
cd /work
```

Type `exit` to leave the container.

### Run commands

The following shell functions start the container, activate the 2DAlphabet
environment, and run commands from the current directory. Add them to `~/.bashrc`:

```bash
twod() {
    local workdir
    workdir="$(pwd -P)"

    singularity exec \
        --bind /cvmfs \
        --bind "${workdir}:/work" \
        --pwd /work \
        /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/jhu-tools/2dalphabet:v2.0 \
        bash -c 'source /home/cmsusr/setup_2dalphabet.sh && exec "$@"' \
        bash "$@"
}

twod_py() {
    twod python3 "$@"
}
```

Load the functions and use `twod_py` for Python scripts or `twod` for other
commands provided by the image:

```bash
source ~/.bashrc
twod_py my_analysis.py
twod root -l fitResults.root
```

Arguments are passed to the selected command. Files written in the
current directory remain on the host after the container exits.

### Example control-region fit
The example in `test/` performs a fit using two control-region categories,
`CR_fail` and `CR_pass`. QCD is estimated from data in the fail category after
subtracting the simulated backgrounds. Then, transfer function is used to
predict the QCD contribution in the pass category. The example is configured
using two files:

- `test/twoDtest_cicd.json` describes the input ROOT histograms, processes,
  regions, systematic uncertainties, two-dimensional binning, and run options.
- `test/test_CR.py` constructs the statistical model and runs the fit, plots,
  limit calculation, and result check.

The main JSON sections are:

| Section | Purpose |
| --- | --- |
| `GLOBAL` | File and histogram naming patterns and the input data directory. |
| `PROCESSES` | Observed data, the top quark background, and the Tprime signal. |
| `REGIONS` | Process content and binning for `CR_fail` and `CR_pass`. |
| `SYSTEMATICS` | Luminosity, cross-section, reweighting, and JES uncertainties. |
| `BINNING` | Ten x-axis bins from 60 to 260 and 22 y-axis bins from 800 to 3000. |
| `OPTIONS` | MC statistical uncertainties, plotting, signal handling, and year. |

`test_CR.py` performs the following steps:

1. Builds the fail-region QCD distribution from data minus simulated
   backgrounds and constructs `0x0`, `1x0`, `0x1`, and `1x1` transfer-function
   models for the pass region.
2. Selects the `1x1` transfer function and the `TprimeB-1800_16` signal, then
   creates a RooFit workspace and Combine datacard.
3. Runs a maximum-likelihood fit and produces background-only and
   signal-plus-background projections.
4. Runs `combine -M AsymptoticLimits` and saves the expected and observed
   95% confidence-level upper limits on the signal-strength parameter `r`.
5. Compares the limits with `test/limits_cicd_reference.json` as a final
   regression check.

Run the example from a clean checkout:

```bash
cd test
twod_py test_CR.py
```

A successful run ends with:

```text
All limits within tolerance
```

The most useful outputs are:

| Output | Contents |
| --- | --- |
| `CR_cicd/base.root` | RooFit workspace containing the full statistical model. |
| `CR_cicd/TprimeB-1800_16_area/card.txt` | Combine datacard for the selected signal and transfer function. |
| `CR_cicd/TprimeB-1800_16_area/fitDiagnosticsTest.root` | Maximum-likelihood fit results. |
| `CR_cicd/TprimeB-1800_16_area/plots_fit_b/` | Background-only post-fit projections and diagnostics. |
| `CR_cicd/TprimeB-1800_16_area/plots_fit_s/` | Signal-plus-background post-fit projections and diagnostics. |
| `CR_cicd/TprimeB-1800_16_area/higgsCombineTest.AsymptoticLimits.mH120.root` | Combine asymptotic-limit output. |
| `limits_cicd.json` | Expected and observed limits in a readable JSON file. |

Inspect the numerical limits with:

```bash
cat limits_cicd.json
```

For the bundled input sample, the result should be close to the
reference values:

| Limit | Signal-strength upper limit |
| --- | ---: |
| Expected, -2 standard deviations | 6.2813 |
| Expected, -1 standard deviation | 8.4705 |
| Expected, median | 12.0000 |
| Expected, +1 standard deviation | 17.1659 |
| Expected, +2 standard deviations | 23.5822 |
| Observed | 12.2201 |

The keys `exp_m2`, `exp_m1`, `exp`, `exp_p1`, and `exp_p2` are the expected
limits at minus two, minus one, median, plus one, and plus two standard
deviations. The `obs` entry is the observed limit. Here, `r = 1` represents the
nominal signal normalization used by the model. The script compares these
values with `test/limits_cicd_reference.json` and reports a failure if their
relative difference is greater than 0.1%.
