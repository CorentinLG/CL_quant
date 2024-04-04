def extract_spectra(s, mask, plot=True, mode = 'mean'):
    """This function allows to extract a spectra from a masked signal
    ** s : the signal
    ** M : the mask to apply
    ** mode : the spectra can be either 'mean' or 'sum' of all unmasked pixels"""

    import copy
    import numpy as np
    s_M = copy.deepcopy(s)
    s_M.data[mask.data, :] = np.nan
    if mode == 'mean':
        sp = s_M._get_signal_signal(np.nanmean(s_M.data, axis=(0, 1)))
    
    elif mode== 'sum':
        sp = s_M._get_signal_signal(np.nansum(s_M.data, axis=(0, 1)))
    
    sp.metadata=s.metadata
    if plot ==True: sp.plot(True)
    return sp
