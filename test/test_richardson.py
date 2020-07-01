import pytest
import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.DEBUG)


from richardson import *

def _richardson_exponent_recovery(p=2.7,hlist=[0.02,0.05,0.08],ytrue=1.0,coef=1.0,tol=1e-3):
    '''When given data of the form y0 + c*h^p, check that extrapolate finds y0 and p correctly'''
    ylist = [ ytrue + coef*(h**p) for h in hlist ]
    d = extrapolate(ylist=ylist,hlist=hlist,order=None)
    assert abs(d['extrapolated'] - ytrue) < tol
    assert abs(d['extrapolation_order'] - p) < tol

def test_richardson_exponent_recovery():
    _richardson_exponent_recovery(p=2.7,hlist=[0.02,0.05,0.08],ytrue=2.0,coef=-1.0)
    _richardson_exponent_recovery(p=1.5,hlist=[0.1,0.2,0.4],ytrue=-5.0,coef=3.0)
    _richardson_exponent_recovery(p=2.2,hlist=[0.01,0.05,0.09],ytrue=10.0,coef=-2.0)

def test_richardson_rejects_non_monotone():
    with pytest.raises(Exception) as e_info:
        extrapolate(ylist=[1.1,1.0,1.2],hlist=[0.1,0.2,0.4],order_seed=2.0)
