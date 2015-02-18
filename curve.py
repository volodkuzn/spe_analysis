import numpy as np

__author__ = 'volod_kuzn'
__version__ = 0.1


def gauss(x, *args):
    """
    Defines multiple Gaussians with characteristics defined by 'params':
    params = ymax, stddev, x0...
    For multiple Gaussians, string together parameters into one long list (so we can fit using leastsq).
    For a single Gaussian, params can be a dict containing the keys ymax, x0, and halfwidth.
    """
    # args = np.abs(args)
    freeparams = 3  # number free parameters per lorentzian, not counting yoffset
    if np.mod(len(args), freeparams) != 0:
        print args
        raise NameError('Incorrect number of parameters supplied to function gaussians. N = %d' % len(args))
    total = np.zeros(len(x))
    for ii in range(0, len(args), freeparams):
        ymax = args[ii]
        stddev = args[ii + 1]
        x0 = args[ii + 2]
        total += ymax * np.exp(-((x - x0) / stddev) ** 2 / 2.0)
    total = np.minimum(total, 60000 * np.ones(len(x)))
    return total