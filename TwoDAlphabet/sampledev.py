'''Correlation-aware uncertainties for 1D projections.

`PostFit3DShapesFromWorkspace --sample-deviations` writes, next to each shape, a
`<proc>_sampledev` TH2 holding the deviation of every bin from the nominal for
every toy sample. Projecting the *sample* keeps the bin-to-bin correlations that summing the
per-bin errors in quadrature (what ROOT's TH3::Projection does) throws away.

Everything here degrades gracefully: if the deviations are absent (an older
postfitshapes file, a pre-fit shape, or data_obs) the loader returns None and
callers leave the ROOT-computed errors alone.
'''

import ROOT
import numpy as np

SAMPLEDEV_SUFFIX = '_sampledev'

# The tool flattens 3D bins as b = ((i-1)*n_y + (j-1))*n_z + k, which is C order
# for shape (n_x, n_y, n_z), so the array below reshapes without any reordering.
_DTYPES = {'TH2F': np.float32, 'TH2D': np.float64}


def _th2_to_array(h):
    '''(n_x, n_y) float64 copy of a TH2, read at the histogram's own precision.'''
    cls = h.ClassName()
    if cls not in _DTYPES:
        raise TypeError('%s: expected TH2F or TH2D, got %s' % (h.GetName(), cls))
    h.BufferEmpty()
    # ROOT stores the buffer y-major and includes under/overflow on both axes.
    shape = (h.GetNbinsY() + 2, h.GetNbinsX() + 2)
    view = np.ndarray(shape, dtype=_DTYPES[cls], buffer=h.GetArray(), order='C')
    return np.array(view[1:-1, 1:-1].T, dtype=np.float64)


def load_sample_deviations(shapes_file, loc_base, region, subregions, time,
                           process, binning, blinded=[]):
    '''Per-sample deviations for one stitched 3D shape.

    Returns an (n_x, n_y, n_z, n_samples) array in the same bin order as the
    histogram produced by stitch_hists_in_x, or None if the deviations were not
    written for this shape.
    '''
    if not binning.is3D:
        return None

    blocks, n_samples = [], None
    for i, subregion in enumerate(subregions):
        base = loc_base.format(r=region, c=subregion, t=time, p=process)
        dev = shapes_file.Get(base + SAMPLEDEV_SUFFIX)
        shape = shapes_file.Get(base)
        if dev == None or shape == None:
            return None

        n_x, n_y, n_z = shape.GetNbinsX(), shape.GetNbinsY(), shape.GetNbinsZ()
        arr = _th2_to_array(dev)
        if arr.shape[0] != n_x * n_y * n_z:
            raise ValueError(
                '%s: %d deviation bins but the shape has %dx%dx%d=%d'
                % (base, arr.shape[0], n_x, n_y, n_z, n_x * n_y * n_z))
        if n_samples is None:
            n_samples = arr.shape[1]
        elif arr.shape[1] != n_samples:
            raise ValueError('%s: %d samples, expected %d from an earlier subregion'
                             % (base, arr.shape[1], n_samples))

        block = arr.reshape(n_x, n_y, n_z, n_samples)
        # stitch_hists_in_x leaves blinded subregions empty; do the same here so
        # the errors cannot disagree with the contents they belong to.
        blocks.append(np.zeros_like(block) if i in blinded else block)

    return np.concatenate(blocks, axis=0)


def _axis_range(lo, hi, n):
    '''ROOT (first,last) bin pair -> slice. plot.py passes (0,-1) for "all bins".'''
    if lo is None or hi is None or lo <= 0 or hi < lo:
        return slice(0, n)
    return slice(lo - 1, hi)


def projection_errors(devs, axis, lo=0, hi=-1):
    '''1-sigma errors of a projection onto `axis`, correlations included.

    `lo`/`hi` restrict the sliced axis, matching the ROOT calls plot.py makes:
    for an X projection they are Y bins, for a Y projection they are X bins, and
    for Z they are ignored (Z is always projected in full). Z is always summed
    over in full for X and Y, again matching plot.py.
    '''
    n_x, n_y, n_z, _ = devs.shape
    if axis == 'X':
        proj = devs[:, _axis_range(lo, hi, n_y), :, :].sum(axis=(1, 2))
    elif axis == 'Y':
        proj = devs[_axis_range(lo, hi, n_x), :, :, :].sum(axis=(0, 2))
    elif axis == 'Z':
        proj = devs.sum(axis=(0, 1))
    else:
        raise ValueError('axis must be X, Y or Z, got %r' % axis)
    return np.sqrt((proj ** 2).mean(axis=-1))


def set_projection_errors(hslice, devs, axis, lo=0, hi=-1):
    '''Overwrite the bin errors of a projection made by ROOT with correlated ones.'''
    errors = projection_errors(devs, axis, lo, hi)
    if hslice.GetNbinsX() != len(errors):
        raise ValueError('%s has %d bins but the projection gives %d'
                         % (hslice.GetName(), hslice.GetNbinsX(), len(errors)))
    for i, err in enumerate(errors):
        hslice.SetBinError(i + 1, float(err))
    return hslice


def max_bin_error_mismatch(devs, shape):
    '''Largest relative disagreement between sqrt(mean(d^2)) and the TH3's own
    per-bin errors. A sanity check on the bin ordering and the stitching: the
    deviations must reproduce the errors the tool already stored.'''
    sampled = np.sqrt((devs ** 2).mean(axis=-1))
    worst = 0.0
    for i in range(shape.GetNbinsX()):
        for j in range(shape.GetNbinsY()):
            for k in range(shape.GetNbinsZ()):
                stored = shape.GetBinError(i + 1, j + 1, k + 1)
                if stored > 0:
                    worst = max(worst, abs(sampled[i, j, k] - stored) / stored)
    return worst
