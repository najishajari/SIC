"""Tools for spectral analysis"""

from __future__ import division, print_function, absolute_import
from scipy.fftpack import fft,fftfreq
import numpy as np
from math import ceil,floor
import sys
from scipy.signal import get_window,welch,signaltools
from matplotlib import pyplot as plt
import warnings
from scipy.lib.six import string_types

def win_sig(x,nperseg,padded='False'):
    """A function just to cut a multidimensional time series into pieces of specific length (nperseg) """
    
    #checking whether the size of time series are smaller than the window lenght
    if nperseg>=x.shape[-1]:
        N=nperseg
        win_num=1
    elif int(nperseg/2.)==nperseg/2.:
        N=int(nperseg/2.)
    else:
        N=int(nperseg/2.)+1
    if padded:
        win_num=ceil(x.shape[-1]/float(N))
    elif nperseg<x.shape[-1]:
        win_num=int(x.shape[-1]/float(N))-1

    #index set manipulation for generating the splitted version of the signals faster                        
    idx_temp=np.indices((win_num,nperseg))
    idx_temp=idx_temp[0]*N+idx_temp[1]

    #padding zeros for the last window when the last window is longer than the remaining of signal
    pad_len=(win_num-1)*N+nperseg-x.shape[-1]
    if padded:
        pad_mat=np.zeros((x.shape[0],pad_len))
        x=np.concatenate((x,pad_mat),axis=-1)

    return x.reshape(-1,x.shape[-1])[:,[idx_temp]].reshape((x.shape[0],win_num,-1))
    
def ndim_welch(x,nperseg=256,window='hanning',scaling = 'density',detrend='constant',fs=1.0,axis=-1,padded=False):
    """Multidimensional Welch method: Calculating Power Spectral Density for time series of multiple [and one] dimension."""

    ##################################################################
    # checking whether the time series lengths are longer than window 
    # length for welch method. If negative the size has been set to time series 
    # length. (Taken from original scipy.signal.welch)
    
    if (not padded) and x.shape[-1] < nperseg:
        warnings.warn('nperseg = %d, is greater than x.shape[%d] = %d, using '
                    'nperseg = x.shape[%d]'
                    % (nperseg, axis, x.shape[axis], axis))
        nperseg = x.shape[-1]

    #########################################################
    # setting the window as is done in original scipy.signal

    if isinstance(window, string_types) or type(window) is tuple:
            win = get_window(window, nperseg)
    else:
        win = np.asarray(window)
        if len(win.shape) != 1:
            raise ValueError('window must be 1-D')
        if win.shape[0] > x.shape[-1]:
            raise ValueError('window is longer than x.')
        nperseg = win.shape[0]     
        
    ######################################################
    #setting the scale as is done in original scipy.signal

    if scaling == 'density':
        scale = 1.0 / (fs * (win*win).sum())
    elif scaling == 'spectrum':
        scale = 1.0 / win.sum()**2
    else:
        raise ValueError('Unknown scaling: %r' % scaling)    

    #########################################################
    #  windowing the signal    
    # (turning the multidimensional time series into multiple t
    # time series of windowed sections using the function 'win_seg')

    windowed_sig=win_sig(x,nperseg,padded)


    ##################################    
    # detrending step
    
    if not detrend:
        detrend_func = lambda seg: seg
    elif not hasattr(detrend, '__call__'):
        detrend_func = lambda seg: signaltools.detrend(seg, type=detrend)
    elif axis != -1:
        # Wrap this function so that it receives a shape that it could
        # reasonably expect to receive.
        def detrend_func(seg):
            seg = np.rollaxis(seg, -1, axis)
            seg = detrend(seg)
            return np.rollaxis(seg, axis, len(seg.shape))
    else:
        detrend_func = detrend
    
    windowed_sig=detrend_func(windowed_sig)
    
    
    #######################################
    # spectral density estimation
    

    # multiplying by window
    windowed_sig=np.multiply(win,windowed_sig)
    
    # calculating the fourier transform    
    windowed_seg_fft=fft(windowed_sig)
    windowed_fft=windowed_seg_fft.T
    
    # returning the spectral density with calcualting outerproducts to get the 
    # crossspectrum matrix and also returning the frequency set

    spec_density=np.mean(np.einsum('...i,...j->...ij',windowed_fft,windowed_fft.conjugate())*scale,axis=1)
    spec_freq=fftfreq(nperseg)
    return spec_freq,np.squeeze(spec_density).real
    