'''2D discrete sine transform convenience functions'''

# In this module, all frequencies are measured in "per array index"
# e.g. the k=0.5 fourier component corresponds to x[i+2] = x[i] for all i

# DST-I means you extend a,b,c to (a,b,c,0,-c,-b,-a,0) and then take
# the FFT Returning only the imaginary parts of the half of the
# coefficients which determine the others.

# Scipy supports the DST but there are some "gotchas" like:
# * composition of DST and inverse DST is not the identity
# * DST in multiple dimensions is not supported directly
# * There is no analogue of "fftfreq" which returns the vector (or matrix)
#   of frequencies of the components returned by the transform itself.

# This small module fixes these and adds the missing 2D DST support.

from __future__ import absolute_import
import numpy as np
import scipy.fftpack

def dst(x,axis=-1):
    '''DST-I of an array'''
    return scipy.fftpack.dst(x,type=1,axis=axis)

def idst(x,axis=-1):
    '''Inverse DST-I of an array'''
    if hasattr(x,'shape'):
        N = x.shape[axis]
    else:
        N = len(x)
    return scipy.fftpack.idst(x,type=1,axis=axis)/float(2*N+2)

def dstfreq(x):
    '''Frequency vector of DST-I for an array like x (if x is iterable) or
       of length x (if not)'''
    try:
        N = len(x)
    except TypeError:
        N = x
    return scipy.fftpack.fftfreq(2*N+2)[1:N+1]

def dst2(x,axes=(-2,-1)):
    '''2D DST-I of an array'''
    return dst( dst(x,axis=axes[0]), axis=axes[1] )

def idst2(x,axes=(-2,-1)):
    '''2D Inverse DST-I of an array'''
    return idst( idst(x,axis=axes[1]), axis=axes[0] )

def dst2freq(x,axes=(-2,-1)):
    '''Return frequency grid of dst2 output for a 2D array like x (if x
       has "shape" attribute) or of shape x (if not)'''
    if hasattr(x,'shape'):
        S = x.shape
    else:
        S = x
    if len(S) != 2:
        raise TypeError('shape must be 2D array or tuple of length 2')
    
    return dstfreq(S[axes[1]])[np.newaxis], dstfreq(S[axes[0]])[np.newaxis].transpose()
