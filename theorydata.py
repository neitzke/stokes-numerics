from __future__ import absolute_import
from cmath import pi,log,exp

import copy
import logging
logger = logging.getLogger(__name__)

def xarXclusterA1A2(co, theta):
    if -pi/6 < theta < pi/6:
        return [co([1,0]), co([0,1])]
    raise NotImplementedError("Requested an undetermined cluster (A1A2, theta > pi/6 not available)")

def xarXclusterA1A3(co, theta):
    if 0 < theta < pi/4:
        return [co([1,0,0]), co([0,1,0]), co([0,0,1])]
    if theta == 0:
        return [co([1,0,0]), co([0,1,0])*(1 - co([-1,0,0]))**(0.5), co([0,0,1])*(1 - co([-1,0,0]))**(0.5)]
    if -pi/4 < theta < 0:
        return [co([1,0,0]), co([0,1,0])*(1 - co([-1,0,0])), co([0,0,1])*(1 - co([-1,0,0]))]
    raise NotImplementedError("Requested an undetermined cluster (A1A3, theta out of range)")

def xarXclusterA2A1(co, theta):
    if -pi/6 < theta < pi/6:
        return [co([1,0]), co([0,1])]
    raise NotImplementedError("Requested an undetermined cluster (A2A1, theta > pi/6 not available)")

def xarXclusterA2A2(co, theta):
    if 0 < theta < 0.3601:
        return [co([1,0,0,0]), co([0,1,0,0]), co([0,0,1,0]), co([0,0,0,1])]
    if theta == 0:
        return [co([1,0,0,0]), co([0,1,0,0])*(1 - co([-1,0,0,0]))**(0.5), co([0,0,1,0]), co([0,0,0,1])]
    if -0.3601 < theta < 0:
        return [co([1,0,0,0]), co([0,1,0,0])*(1 - co([-1,0,0,0])), co([0,0,0,1])]
    raise NotImplementedError("Requested an undetermined cluster (A2A2, theta out of range)")

import numpy as np

def plucker(vertexlist, vertices):
    """Plucker coordinate with vertices numbered from 1 to nvert, cyclically"""
    submatrix = [vertexlist[((n-1) % len(vertexlist))] for n in vertices]
    coord = np.linalg.det(submatrix)
    return coord

def nonplucker36(vertexlist, vertexpairs):
    """Get non-Plucker coordinate for Gr(3,6), with vertices given as 3 pairs"""
    def cross(a,b):
        return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    def crossvertices(n):
        return cross(vertexlist[(vertexpairs[n][0]-1) % 6], vertexlist[(vertexpairs[n][1]-1) % 6])
    return np.linalg.det([crossvertices(n) for n in range(3)])

def fdclusterA23(vertexlist, offset=0):
    def a(n1,n2):
        return plucker(vertexlist, [n1+offset,n2+offset])
    return [a(1,2),a(2,3),a(3,1)]

def fdclusterA24(vertexlist, offset=0):
    def a(n1,n2):
        return plucker(vertexlist, [n1+offset,n2+offset])
    return [a(1,3)/a(1,2),a(1,2),a(2,3),a(3,4),a(4,1)]

def fdclusterA25(vertexlist, offset=0):
    def a(n1,n2):
        return plucker(vertexlist, [n1+offset,n2+offset])
    return [a(1,3),a(3,5),a(1,2),a(2,3),a(3,4),a(4,5),a(5,1)]

def fdclusterA26(vertexlist, offset=0):
    def a(n1,n2):
        return plucker(vertexlist, [n1+offset,n2+offset])
    return [a(1,3),a(1,4),a(4,6),a(1,2),a(2,3),a(3,4),a(4,5),a(5,6),a(6,1)]

def fdclusterA34(vertexlist, offset=0):
    def a(n1,n2,n3):
        return plucker(vertexlist, [n1+offset,n2+offset,n3+offset])
    return [a(1,2,3),a(2,3,4),a(3,4,1),a(4,1,2)]

def fdclusterA35(vertexlist, offset=0):
    def a(n1,n2,n3):
        return plucker(vertexlist, [n1+offset,n2+offset,n3+offset])
    return [a(1,2,4),a(2,4,5),a(1,2,3),a(2,3,4),a(3,4,5),a(4,5,1),a(5,1,2)]

def fdclusterA36(vertexlist, offset=0):
    def a(n1,n2,n3):
        return plucker(vertexlist, [n1+offset,n2+offset,n3+offset])
    npc = nonplucker36(vertexlist, [[3+offset,4+offset],[2+offset,1+offset],[6+offset,5+offset]])
    return [npc,a(1,2,5),a(2,5,6),a(3,5,6),a(1,2,3),a(2,3,4),a(3,4,5),a(4,5,6),a(5,6,1),a(6,1,2)]


def minusone(charge):
    return -1

def A1A2rayf1(gamma, charges, xlists):
    [x,y] = xlists

    # write the ray function in terms of quantities which are small, to avoid
    # overflow problems in intermediate formulas
    Xi = np.exp(-x)
    XY = np.exp(x+y)
    Y = np.exp(y)

    if gamma == [1,0]:
        f = - np.log(1 + Y + XY)
    if gamma == [0,1]:
        f = - np.log(1 + Xi) + np.log(1 + Y + XY)

    return f

def A1A2rayf2(gamma, charges, xlists):
    [x,y] = xlists

    # write the ray function in terms of quantities which are small, to avoid
    # overflow problems in intermediate formulas
    X = np.exp(x)
    XYi = np.exp(-x-y)
    Yi = np.exp(-y)

    if gamma == [1,0]:
        f = np.log(1 + Yi + XYi)
    if gamma == [0,1]:
        f = np.log(1 + X) - np.log(1 + Yi + XYi)

    return f

# background on the lattices and bases:

# the data we give is adapted to a specific cluster

# nodes of quiver form a basis for A-lattice
# inside here have "unfrozen part" as a sublattice
# "B" is B-matrix, a partially defined bilinear form on A-lattice
# "Xbasis" gives a basis for the X-lattice as sublattice of A-lattice: 
# each basis element is some linear combination of columns of the B-matrix
# (i.e. take linear combination of the basis vectors 
# obtained by duality from the unfrozen nodes of the quiver)
# "bpshypers" gives charges of BPS hypermultiplets relative to the chosen X-lattice basis
# "bpsvectors" gives charges of BPS vectormultiplets

# "fdcluster" function should give the A-cluster variables in order
# X-cluster is then automatically computed from A-cluster as needed

theorydata = {

    "A1A2": {

                "name": "A1A2",

                "K": 2,
                "Pcoefs": [
                           [-1,0,0,1]
                          ],
                "Xrank": 2,
                "Arank": 7,
                "bpshypers": [[1,0], [0,1], [1,1]],

                # central charge computed by:
                # w = Exp[2 Pi I/3]
                # 2*NIntegrate[Sqrt[1 - z^3], {z, 1, w}, PrecisionGoal -> 50, MaxRecursion -> 20, WorkingPrecision -> 100] // FullForm
                # 2*w*NIntegrate[Sqrt[1 - z^3], {z, 1, w}, PrecisionGoal -> 50, MaxRecursion -> 20, WorkingPrecision -> 100] // FullForm
                "centralcharges": [-2.52392778958581767011503434+1.45719038873254896709196501j, -2.9143807774650979341839300j],
                "angles": [0,0],

                "rays": [
                          {"phase": 1j, "charges": [[1,0],[0,1]], "function": A1A2rayf1},
                          {"phase": -1j, "charges": [[1,0],[0,1]], "function": A1A2rayf2}
                        ],

                "cpshift": 7.0/40.0,

                "Rlist": [exp(0.375*t).real for t in range(-54,8)],
                "abshlist": [exp(0.375*t).real for t in range(-8,30)],
#                "abshlist": [exp(0.1*t).real for t in range(-40,20)],
                "thetalist": [0.0],

                "Xclusternames": ["X1", "X2"],

                "B": [[ 0,-1],
                      [ 1, 0],
                      [ 1, 0],
                      [-1, 0],
                      [ 0, 1],
                      [ 0,-1],
                      [-1, 1]
                     ],

                "Xbasis": [[ -1,  0],
                           [  1,  1]
                          ],

                "Xsigns": [-1, 1],

                "fdAclusterfunc": fdclusterA25,
                "xarXclusterfunc": xarXclusterA1A2,

                "sigma": minusone
            },

    "A1A3": {
                "name": "A1A3",

                "K": 2,
                "Pcoefs": [
                           [-1,0,0,0,1]
                          ],

                "Xrank": 3,
                "Arank": 9,
                "bpshypers": [[1,0,0],[0,1,0],[0,0,1],[1,-1,0],[1,0,-1],[1,-1,-1]],

                "centralcharges": [3.4960767390562006, 1.7480383695281214+1.7480383695281214j, 1.7480383695281214+1.7480383695281214j],
                "angles": [0,0,0],

                "cpshift": None,

                "Rlist": [1e-8,2e-8,5e-8,1e-7,2e-7,5e-7,1e-6,2e-6,5e-6,1e-5,2e-5,3e-5,5e-5,7e-5,1e-4,2e-4,3e-4,5e-4,7e-4,1e-3,2e-3,3e-3,5e-3,7e-3,1e-2,2e-2,3e-2,5e-2,7e-2,1e-1,2e-1,3e-1,5e-1,7e-1,1,2,3,5,7],
                "thetalist": [0.0],
                "abshlist": [exp(0.1*t).real for t in range(-30,20)],

                "Xclusternames": ["X1", "X2", "X3"],

                "B": [[ 0, 1, 0],
                      [-1, 0,-1],
                      [ 0, 1, 0],
                      [ 1, 0, 0],
                      [-1, 0, 0],
                      [ 1,-1, 0],
                      [ 0, 0, 1],
                      [ 0, 0,-1],
                      [ 0,-1, 1]
                     ],

                "Xbasis": [[ 0, 1, 0],
                           [-1, 0, 0],
                           [ 0, 0,-1],
                          ],

                "Xsigns": [-1, 1, -1],

                "fdAclusterfunc": fdclusterA26,
                "xarXclusterfunc": xarXclusterA1A3,

                "sigma": minusone
            },

    "A2A1": {
                "name": "A2A1",

                "K": 3,
                "Pcoefs": [
                           [0],
                           [-0.5,0,0.5]
                          ],

                "Xrank": 2,
                "Arank": 7,
                "bpshypers": [[1,0], [0,1], [1,1]],

                "centralcharges": [-2.0032428141401496922 + 1.1565727779959988786j, -2.3131455559919977573j],
                "angles": [0,0],

                "cpshift": 13.0/60.0,

                "Rlist": [1e-8,3e-8,5e-8,1e-7,3e-7,5e-7,1e-6,3e-6,5e-6,1e-5,2e-5,4e-5,6e-5,8e-5,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.0010,0.0011,0.0012,0.0014,0.0016,0.0018,0.0020,0.0030,0.0040,0.0050,0.0060,0.0070,0.0080,0.0090,0.0100,0.0200,0.0300,0.0400,0.0500,0.0600,0.0700,0.0800,0.0900,0.1000,0.1500,0.2000,0.2500,0.3000,0.3500,0.4000,0.4500,0.5000,0.5500,0.6000,0.6500,0.7000,0.7500,0.8000,0.8500,0.9000,0.9500,1.0000,1.1000,1.2000,1.3000,1.4000,1.5000,1.6000,1.7000,1.8000,1.9000,2.0000,2.5000,3.0000,3.5000,4.0000,5.0000,6.0000,7.0000,8.0000,9.5000],
                "abshlist": [ exp(0.3*t).real for t in range(-10,50) ],
                "thetalist": [0.0],

                "Xclusternames": ["X1", "X2"],

                "B": [[ 0,-1],
                      [ 1, 0],
                      [ 1, 0],
                      [-1, 1],
                      [ 0,-1],
                      [ 0, 1],
                      [-1, 0]
                     ],

                "Xbasis": [[-1, 0],
                           [ 1, 1]
                          ],

                "Xsigns": [1, 1],

                "fdAclusterfunc": fdclusterA35,
                "xarXclusterfunc": xarXclusterA2A1,

                "sigma": minusone
            },

    "A2A2": {
                "name": "A2A2",

                "K": 3,
                "Pcoefs": [
                           [0],
                           [-1,0,-1.5,0.5]
                          ],

                "Xrank": 4,
                "Arank": 10,
                "bpshypers": [[1,0,0,0],
                        [0,-1,-1,-1],
                        [-1,1,1,1],
                        [0,1,0,0],
                        [1,-1,-1,0],
                        [-1,0,1,0],
                        [0,-1,-1,0],
                        [1,0,-1,-1],
                        [-1,1,2,1],
                        [1,1,0,0],
                        [1,-2,-2,-1],
                        [-2,1,2,1]
                       ],
                # central charges computed numerically, may not have as much accuracy as number of places given here
                "centralcharges": [2.30298111643, 5.47033108653 + 4.48792389285j, -4.31884052829 + 2.49348374162j, -4.98696748321j],
                "angles": [0,0,0,0],

                "cpshift": None,

                "Rlist": [1e-8,3e-8,5e-8,1e-7,3e-7,5e-7,1e-6,3e-6,5e-6,1e-5,2e-5,4e-5,6e-5,8e-5,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.0010,0.0011,0.0012,0.0014,0.0015,0.0016,0.0018,0.0020,0.0025,0.0030,0.0040,0.0050,0.0060,0.0070,0.0080,0.0090,0.0100,0.0150,0.0200,0.0250,0.0300,0.0350,0.0400,0.0500,0.0600,0.0700,0.0800,0.0900,0.1000,0.1500,0.2000,0.2500,0.3000,0.3500,0.4000,0.4500,0.5000,0.5500,0.6000,0.6500,0.7000,0.7500,0.8000,0.8500,0.9000,0.9500,1.0000,1.1000,1.2000,1.3000,1.4000,1.5000,1.6000,1.7000,1.8000,1.9000,2.0000,2.5000,3.0000,3.5000,4.0000,5.0000,6.0000,7.0000],
                "abshlist": [ exp(0.1*t).real for t in range(-33,50)],
                "thetalist": [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3],
        
                "Xclusternames": ["X1", "X2", "X3", "X4"],

                "B": [[ 0, 1,-1, 1],
                      [-1, 0, 0, 0],
                      [ 1, 0, 0, 0],
                      [-1, 0, 0, 0],
                      [ 0,-1, 0, 0],
                      [ 0, 0, 0,-1],
                      [ 1, 0, 1, 0],
                      [-1,-1, 0, 0],
                      [ 1, 0, 0,-1],
                      [ 0, 0, 1, 0]
                     ],

                "Xbasis": [
                           [ 0, 1, 0, 0],
                           [-1, 0, 0, 0],
                           [ 0, 1, 1, 0],
                           [ 0, 0,-1,-1]],

                "Xsigns": [-1, 1, 1, 1],

                "fdAclusterfunc": fdclusterA36,
                "xarXclusterfunc": xarXclusterA2A2,

                "sigma": minusone                
            },

}


def A1A2theorydeformed(theoryname, cr = 1, ci = 0.0, Lambda = 0, theta1 = 0, theta2 = 0): 
    '''the example of P = z^3 - Lambda z - c'''

    c = cr + 1j*ci
    Pcoefs = [[-c, -Lambda, 0, 1]]

    if Lambda == 0:
        centralcharges = [(-2.52392778958581767011503434+1.45719038873254896709196501j)*(c**(5.0/6.0)), (-2.9143807774650979341839300j)*(c**(5.0/6.0))]
    else:
        from periods import turningpoints,allintegrals
        # compute central charges by direct integration
        # should work for c near 1 and Lambda small, at least

        w = exp(2*pi*1j/3.0)
        tp = turningpoints(2, Pcoefs, [1,w,w**2])
        logger.info("Computed turning points: %s" % tp)

        zbase = 0
        xbase = [c**0.5,-(c**0.5)]
        logger.info("Computed base sheets: %s" % xbase)

        periods01 = allintegrals(2, Pcoefs, tp[0], tp[1], 0, xbase)
        periods12 = allintegrals(2, Pcoefs, tp[1], tp[2], 0, xbase)

        centralcharges = [2*periods01[0], 2*periods12[0]]

    logger.info("Computed central charges: %s" % centralcharges)

    th = copy.deepcopy(theorydata["A1A2"])
    th["name"] = theoryname
    th["Pcoefs"] = Pcoefs
    th["centralcharges"] = centralcharges
    th["angles"] = [theta1, theta2]

    return th

def A2A1theorydeformed(theoryname, c = 1, theta1 = 0, theta2 = 0):
    '''the example P2 = c, P3 = 0.5z^2 - 0.5''' 
    Pcoefs = [
               [c],
               [-0.5,0,0.5]
             ]

    if c == 0:
        centralcharges = [-2.0032428141401496922 + 1.1565727779959988786j, -2.3131455559919977573j]
    else:
        from periods import turningpoints,spectralpoints,allintegrals
        # compute central charges by direct integration
        # should work for c near 1 at least
        tp = turningpoints(3, Pcoefs, [-1-0.4j,-1+0.4j,1+0.4j,1-0.4j])
        logger.info("Computed turning points: %s" % tp)

        zbase = 0
        xbase = spectralpoints(3, Pcoefs, zbase, [0.4238, -0.2119 - 1.0652j, -0.2119 + 1.0652j])
        logger.info("Computed base sheets: %s" % xbase)

        periods02 = allintegrals(3, Pcoefs, tp[0], tp[2], zbase, xbase)
        periods13 = allintegrals(3, Pcoefs, tp[1], tp[3], zbase, xbase)
        Z1 = periods02[1]-periods02[0]
        Z2 = periods13[0]-periods13[2]

        centralcharges = [-Z2, Z1+Z2]

    logger.info("Computed central charges: %s" % centralcharges)

    th = copy.deepcopy(theorydata["A2A1"])
    th["name"] = theoryname
    th["Pcoefs"] = Pcoefs
    th["centralcharges"] = centralcharges
    th["angles"] = [theta1, theta2]

    return th

# TODO: make this mechanism more general
def getdata(theoryname):
    tnlist = theoryname.split("_")
    args = [i.split("=") for i in tnlist[1:]]
    argdict = {}
    # things past the first underscore are interpreted as floating-point params
    # eg A1A2_c=0.75_Lambda=0.3
    for arg in args:
        argdict[arg[0]] = complex(arg[1])
    if len(tnlist) == 1:
        # no extra arguments
        return theorydata[theoryname]
    if tnlist[0] == "A2A1":
        return A2A1theorydeformed(theoryname, **argdict)
    if tnlist[0] == "A1A2":
        return A1A2theorydeformed(theoryname, **argdict)
    raise NotImplementedError("gettheorydata() didn't recognize its arguments")
