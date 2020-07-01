from framedata import *
from theorydata import plucker,nonplucker36

def eq(x,y):
    return abs(x-y) < 1e-10

def test_signs_A1A2():
    def a(i,j):
        return plucker(vl, [i,j])
    fd = computeFrames(theoryname="A1A2", R=1, theta=0, pde_nmesh=127, oper=False)
    vl = fd.vertexlist

    x1 = (a(2,3)*a(1,5))/(a(1,2)*a(3,5))
    assert x1>0, "x1 = %s not positive in A1A2" % x1
    x2 = (a(2,3)*a(4,5)*a(1,3))/(a(1,2)*a(3,4)*a(3,5))
    assert x2>0, "x2 = %s not positive in A1A2" % x2

    [y1,y2] = fd.getCluster()
    assert eq(x1,-y1), "x1 != -y1 in A1A2: x1 = %s, -y1 = %s" % (x1,-y1)
    assert eq(x2,-y2), "x2 != -y2 in A1A2: x2 = %s, -y2 = %s" % (x2,-y2)

def test_signs_A1A3():
    def a(i,j):
        return plucker(vl, [i,j])
    fd = computeFrames(theoryname="A1A3", R=1, theta=0, pde_nmesh=127, oper=False)
    vl = fd.vertexlist

    x1 = (a(1,3)*a(4,6))/(a(1,6)*a(3,4))
    assert x1>0, "x1 = %s not positive in A1A3" % x1
    x2 = (a(2,3)*a(1,4))/(a(1,2)*a(3,4))
    assert x2>0, "x2 = %s not positive in A1A3" % x2
    x3 = (a(1,4)*a(5,6))/(a(4,5)*a(1,6))
    assert x3>0, "x3 = %s not positive in A1A3" % x3
    xf = (a(2,3)*a(4,5)*a(1,6))/(a(1,2)*a(3,4)*a(5,6))
    assert xf>0, "xf = %s not positive in A1A3" % xf

    [y1,y2,y3] = fd.getCluster()
    assert eq(x1,-y1), "x1 != -y1 in A1A3: x1 = %s, -y1 = %s" % (x1,-y1)
    assert eq(x2,-y2), "x2 != -y2 in A1A3: x2 = %s, -y2 = %s" % (x2,-y2)
    assert eq(x3,-y3), "x3 != -y3 in A1A3: x3 = %s, -y3 = %s" % (x3,-y3)

def test_signs_A2A1():
    def a(i,j,k):
        return plucker(vl, [i,j,k])
    fd = computeFrames(theoryname="A2A1", R=1, theta=0, pde_nmesh=127, oper=False)
    vl = fd.vertexlist

    x1 = (a(2,3,4)*a(1,4,5))/(a(1,2,4)*a(3,4,5))
    assert x1>0, "x1 = %s not positive in A2A1" % x1
    x2 = (a(2,4,5)*a(1,2,3)*a(1,4,5))/(a(1,2,5)*a(1,2,4)*a(3,4,5))
    assert x2>0, "x2 = %s not positive in A2A1" % x2

    [y1,y2] = fd.getCluster()
    assert eq(x1,-y1), "x1 != -y1 in A2A1: x1 = %s, -y1 = %s" % (x1,-y1)
    assert eq(x2,-y2), "x2 != -y2 in A2A1: x2 = %s, -y2 = %s" % (x2,-y2)

def test_signs_A2A2():
    def a(i,j,k):
        return plucker(vl, [i,j,k])
    def q(i,j,k,l,m,n):
        return nonplucker36(vl, [[i,j],[k,l],[m,n]])
    fd = computeFrames(theoryname="A2A2", R=1, theta=0, pde_nmesh=127, oper=False)
    vl = fd.vertexlist

    x1 = (q(1,2,3,4,5,6) / (a(1,2,3)*a(4,5,6)))
    assert x1>0, "x1 = %s not positive in A2A2" % x1
    x2 = (a(1,2,5)*a(3,5,6)*a(4,5,6))/(a(1,5,6)*a(2,5,6)*a(3,4,5))
    assert x2>0, "x2 = %s not positive in A2A2" % x2
    x3 = (a(1,2,6)*a(3,4,5))/(a(1,2,3)*a(4,5,6))
    assert x3>0, "x3 = %s not positive in A2A2" % x3
    x4 = (a(1,5,6)*a(2,3,4))/(a(1,2,6)*a(3,4,5))
    assert x4>0, "x4 = %s not positive in A2A2" % x4

    [y1,y2,y3,y4] = fd.getCluster()
    assert eq(x1,-y1), "x1 != -y1 in A2A2: x1 = %s, -y1 = %s" % (x1,-y1)
    assert eq(x2,-y2), "x2 != -y2 in A2A2: x2 = %s, -y2 = %s" % (x2,-y2)
    assert eq(x3,-y3), "x3 != -y3 in A2A2: x3 = %s, -y3 = %s" % (x3,-y3)
    assert eq(x4,-y4), "x4 != -y4 in A2A2: x4 = %s, -y4 = %s" % (x4,-y4)
