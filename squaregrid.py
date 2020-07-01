from __future__ import absolute_import
from __future__ import print_function
import numpy as np

class SquareGrid(object):
    """Square mesh centered at the origin in the complex plane.  This is
       essentially just a specialized/enhanced version of
       numpy.meshgrid which bundles the coordinate vectors and grids
       together with spacing data, and wich supports access to the
       grid points as either (float x, float y) or (complex z).

    """
    def __init__(self,r,n):
        """r : inradius of the square
           n : number of mesh points along each side
        """
        self.n = int(n)
        self.nx = self.ny = n

        self.N = n*n
        self.r = float(r)
        self.xmin = self.ymin = -self.r
        self.xmax = self.ymax = self.r
        
        self.dx = (self.xmax - self.xmin) / float(n-1)
        self.dy = (self.xmax - self.xmin) / float(n-1)

        self.x = np.linspace(self.xmin,self.xmax,num=n,endpoint=True)
        self.y = np.linspace(self.ymin,self.ymax,num=n,endpoint=True)
        
        self.xm, self.ym = np.meshgrid(self.x,self.y)

        self.xv = self.xm.ravel()
        self.yv = self.ym.ravel()

        self.zm = self.xm + 1j*self.ym
        self.zv = self.xv + 1j*self.yv


    def __str__(self):
        return "<SquareGrid(r=%f,n=%d),id=0x%x>" % (self.r,self.n,id(self))

    def toidx(self,i,j):
        return i + j*self.n

    def fromidx(self,k):
        return k % self.n, int(k/self.n)
