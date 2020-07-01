'''Test that framedata and integralequations can write, read, and get same data back again'''
from __future__ import absolute_import
from __future__ import print_function

import os
import warnings
import sys
import numpy as np
import datetime
import fileinfo
from six import string_types
import pytest

# TODO: Some tests give warnings under python 2 and certain recent numpy/scipy versions
# Investigate a better solution than ignoring them.

# TODO: The decorators to suppress warnings are due to incompatibility of python 2.7 
# with some recent numpy/scipy combinations.  Should be fixed properly.

def isprimitive(x):
    return x==None or isinstance(x,(str,int,float,complex,bool))

def recursive_weak_eq(a,b):
    '''Recursively compare objects by comparing all arrays and primitive types encountered, but only testing other types for equality of type, not value'''
    if isinstance(a,dict):
        if not recursive_weak_eq(list(sorted(a.keys())),list(sorted(b.keys()))):
            return False
        return all( recursive_weak_eq(a[k],b[k]) for k in a )
    if isinstance(a,list):
        return all( recursive_weak_eq(x,y) for x,y in zip(a,b) )
    if isinstance(a,np.ndarray):
        res = (a==b).all()
        if not res:
            print('Unequal arrays:',a,b)
        return res
    if isinstance(a,string_types):
        res = a.encode('utf-8')==b.encode('utf-8')
        if not res:
            print('Unequal strings:',a,b)
        return a.encode('utf-8')==b.encode('utf-8')
    if isprimitive(a) or isinstance(a,datetime.datetime):
        if a!=b:
            print('Unequal primitive types:',a,b)
        return a==b
    print('Only comparing types, not values "{}" to "{}"'.format(a,b))
    return type(a)==type(b)

def instances_weak_eq(a,b):
    '''Check that non-callcable non-special attributes of a,b are equal in the sense of recursive_weak_eq'''
    if dir(a) != dir(b):
        return False
    for x in dir(a):
        if x[0] == '_':
            continue
        if callable(getattr(a,x)):
            continue
        if not recursive_weak_eq(getattr(a,x),getattr(b,x)):
            return False
    return True

@pytest.mark.filterwarnings("ignore:numpy.dtype size changed","ignore:numpy.ufunc size changed")
def test_framedata_hitchin_io_identity(tmpdir):
    '''Write Hitchin section frame data (using filename generator logic), read, check equality'''
    for k in fileinfo.PATH_VARS:
        setattr(fileinfo,k,str(tmpdir))
    import framedata
    print('Computing and saving frame data (Hitchin)')
    fout, ffn = framedata.computeAndSaveFrames('A1A2',R=0.1,theta=0.0,pde_nmesh=63,pde_thresh=1e-5,ode_thresh=1e-5,oper=False)
    print(ffn)
    if not isinstance(fout,framedata.framedata):
        raise RuntimeError('computeAndSaveFrames return value is not an instance of framedata')
    print('Loading saved frame data (Hitchin)')
    fin = framedata.loadFrameData('A1A2',R=0.1,theta=0.0,pde_nmesh=63,oper=False)
    assert instances_weak_eq(fout,fin)

@pytest.mark.filterwarnings("ignore:numpy.dtype size changed","ignore:numpy.ufunc size changed")
def test_framedata_oper_io_identity(tmpdir):
    '''Write oper frame data (using filename generator logic), read, check equality'''
    for k in fileinfo.PATH_VARS:
        setattr(fileinfo,k,str(tmpdir))
    import framedata
    print('Computing and saving frame data (oper)')
    fout, ffn = framedata.computeAndSaveFrames('A2A1',absh=10.0,theta=0.0,oper=True,ode_thresh=1e-5)
    print(ffn)
    if not isinstance(fout,framedata.framedata):
        raise RuntimeError('computeAndSaveFrames return value is not an instance of framedata')
    print('Loading saved frame data (oper)')
    fin = framedata.loadFrameData('A2A1',1.0,absh=10.0,theta=0.0,oper=True)
    assert instances_weak_eq(fout,fin)

@pytest.mark.filterwarnings("ignore:numpy.dtype size changed","ignore:numpy.ufunc size changed")
def test_xar_io_identity(tmpdir):
    '''Write xarray (using filename generator logic), read, check equality'''
    for k in fileinfo.PATH_VARS:
        setattr(fileinfo,k,str(tmpdir))
    import integralequations
    print('Computing and saving xar data')
    xout, xfn = integralequations.computeAndSaveXar('A1A2',R=0.1,L=2**5,oper=False,savexar=True,saveclusters=False)
    print(xfn)
    if not isinstance(xout,integralequations.xarray):
        raise RuntimeError('computeAndSaveXar return value is not an instance of xarray')
    print('Loading saved xar data')
    xin = integralequations.loadXar('A1A2',R=0.1,oper=False)
    assert instances_weak_eq(xout,xin)
    
@pytest.mark.filterwarnings("ignore:numpy.dtype size changed","ignore:numpy.ufunc size changed")
def test_xarcluster_io_identity(tmpdir):
    '''Write clusters from xarray (using filename generator logic), read, check equality'''
    for k in fileinfo.PATH_VARS:
        setattr(fileinfo,k,str(tmpdir))
    import integralequations
    print('Computing and saving xarcluster data')
    xout, xfn = integralequations.computeAndSaveXar('A1A2',R=0.1,L=2**5,oper=False,savexar=False,saveclusters=True)
    print(xfn)
    if not isinstance(xout,integralequations.xarray):
        raise RuntimeError('computeAndSaveXar return value is not an instance of xarray')
    print('Loading saved xarcluster data')
    xin = integralequations.loadXarCluster('A1A2',R=0.1,oper=False)
    assert instances_weak_eq(xout.getCluster(),xin)
