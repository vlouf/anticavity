import numpy as np
from numba import jit


@jit(nopython=True)
def anticavity(reflectivity, mask, window_len = 3, refl_thld = 10):
    """
    Anticavity plugs the "speckle holes" created by the filtering of noise.

    Parameters:
    ===========
    reflectivity: ndarray
        Reflectivity fields.
    mask: ndarray
        Filtered mask.
    window_len: int
        Size of the moving window.
    refl_thld: float
        Reflectivity threshold.

    Returns:
    ========
    outmask: ndarray
        Filtered mask.
    """
    thrld = window_len // 3
    nazi, ngate = reflectivity.shape
    outmask = mask.copy()
    window = np.zeros((window_len))
    refl_vec = np.zeros((window_len))
    cnt = 0
    for na in range(nazi):
        for ng in range(ngate - window_len):
            window = mask[na, ng: ng + window_len]
            w = window.sum()
            if w == window_len:
                continue
            if w == 0:
                continue
            if w > thrld:
                continue
            refl_vec = reflectivity[na, ng: ng + window_len]
            mean_refl = np.mean(refl_vec[window])
            for k in range(window_len):
                if window[k] == 1:
                    continue
                if np.abs(refl_vec[k] - mean_refl) <= refl_thld:
                    outmask[na, ng + k] = 1
                    cnt += 1

    return outmask