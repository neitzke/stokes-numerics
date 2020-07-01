"""Richardson extrapolation"""

import numpy as np
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def extrapolate(ylist,hlist,order=None,order_seed=2.0):
    '''Given three values (f1,f2,f3) at three step sizes (h1,h2,h3) with
    h1<h2<h3.  If `order` is given, it is assumed to be the theoretical order (in
    h) of the method and is used for extrapolation.  Otherwise, an observed order
    is computed and used.  The observed order calculation uses Newton's method,
    and `order_seed` specifies the starting value in that search.

    Return is a dictionary with attributes

      "extrapolated" : The extrapolated f-value
      "delta" : Absolute difference between f1 and extrapolated value
      "extrapolation_order" : The order in h used in the extrapolation
      "forced_order" : True if the order was specified rather than observed.
    '''

    # Following p320 of Oberkampf & Roy, "Verification and Validation in Scientific Computing"
    
    logger.debug('Applying Richardson extrapolation to:\nhlist={}\nylist={}'.format(hlist,ylist))

    if not all( h.real == h for h in hlist ):
        raise ValueError('hlist entries must be real')

    if not all( h < hnext for h,hnext in zip(hlist,hlist[1:]) ):
        raise ValueError('hlist must be strictly increasing')

    if not all( y.real == y for y in ylist ):
        raise ValueError('ylist entries must be real')

    h1,h2,h3=hlist
    f1,f2,f3=ylist
    r12 = h2/h1
    r23 = h3/h2
    if f2==f1:
        raise ValueError('ylist must have distinct elements')
    diffquot = (f3-f2)/(f2-f1)

    if diffquot <= 0:
        raise ValueError('ylist must be monotone')

    if order:
        logger.debug('Using fixed order p={}'.format(order))
        p = order
    else:
        logger.debug('Using order recovery with seed p0={}'.format(order_seed))
        p = order_seed
        dp = 1.0
        n = 0
        while (dp > 0.0001) and (n < 100):
            p_old,p = p, np.log( (r12**p - 1.0)*diffquot + r12**p ) / np.log(r12*r23)
            dp = np.abs(p - p_old)
            if np.isnan(p) or np.isinf(p):
                raise Exception('Failed to compute observed order (invalid float value)')

        if n == 100:
            raise Exception('Failed to compute observed order (exceeded 100 iterations)')

    # Check that the solution p satisfies the required equation
    LHS = (f3-f2)/(r23**p -1.0)
    RHS = r12**p * (f2-f1)/(r12**p -1.0)
    if abs(LHS-RHS) > 0.0001*(abs(LHS)+abs(RHS)):
        raise Exception('Failed to compute observed order (convergence to non-solution)')

    logger.debug('Found observed order p={} ({} iterations)'.format(p,n+1))

    # Extrapolated value
    fbar = f1 + (f1-f2)/(r12**p - 1.0)
    
    return {"extrapolated": fbar, "delta": abs(fbar-f1), "extrapolation_order": p, "forced_order": (order!=None)}
