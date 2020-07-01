'''Hyperbolic ideal polygon class that can convert itself to matplotlib patches '''
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

def _near(a,b):
    return abs(a-b) < 1e-6

class H2IdealPoint:
    def __init__(self, theta):
        self.theta = theta
        self.p = np.array([np.cos(self.theta), np.sin(self.theta)])

    def patch(self, **kwargs):
        return plt.Circle(self.p, **kwargs)

class H2Line:
    def __init__(self, theta0=-0.25, theta1=0.25):
        self.theta0 = np.mod(theta0,2*np.pi)
        self.theta1 = np.mod(theta1,2*np.pi)

        dt0 = abs(self.theta0 - self.theta1)
        good_dt = min(dt0, 2*np.pi - dt0)
        #print(f"dt0={dt0} good_dt={good_dt} theta0={self.theta0} theta1={self.theta1}")
        if _near(self.theta1-self.theta0, good_dt):
            self.positive = True
        elif _near(self.theta1+2*np.pi-self.theta0, good_dt):
            self.theta1 = self.theta1 + 2*np.pi
            self.positive = True
        elif _near(self.theta0-self.theta1, good_dt):
            self.positive = False
        elif _near(self.theta0-self.theta1+2*np.pi, good_dt):
            self.theta1 = self.theta1 - 2*np.pi
            self.positive = False
        else:
            raise ValueError("Failed to linearize!")
        self.mid = 0.5*(self.theta1 + self.theta0)
        self.dt = self.theta1 - self.theta0
        self.p0 = H2IdealPoint(self.theta0)
        self.p1 = H2IdealPoint(self.theta1)
        if abs(self.dt - np.pi) < 0.0001:
            self.is_line = True
            self.crad = np.inf
            self.ctr = np.array([np.nan, np.nan])
            self.radius = np.inf
        else:
            self.is_line = False
            self.crad = 1.0 / abs(np.cos(0.5*self.dt))
            self.ctr = self.crad * np.array([np.cos(self.mid),np.sin(self.mid)])
            self.radius = abs(np.tan(0.5*self.dt))

    def endpoints(self,**kwargs):
        return [self.p0.patch(**kwargs), self.p1.patch(**kwargs)]

    def points_along(self,n=100):
        v0 = self.p0.p - self.ctr
        v1 = self.p1.p - self.ctr
        t0 = np.arctan2(v0[1],v0[0])
        t1 = np.arctan2(v1[1],v1[0])
        if abs(t1-t0) > np.pi:
            t1 += np.sign(t0-t1)*2.0*np.pi
        return [ self.ctr + self.radius*np.array([np.cos(t),np.sin(t)]) for t in np.linspace(t0,t1,n) ]

    def patch(self, **kwargs):
        if self.is_line:
            return patches.Polygon([self.p0.p,self.p1.p],closed=False,**kwargs)
        else:
            return patches.Polygon(self.points_along(),closed=False,**kwargs)

class H2IdealPolygon:
    def __init__(self,angles):
        self.angles = list(angles)
        self.n = len(self.angles)
        
    def vertices(self):
        return (H2IdealPoint(t) for t in self.angles)

    def edges(self):
        return (H2Line(self.angles[i],self.angles[(i+1) % self.n]) for i in range(self.n))
    
    def edge_patches(self,**kwargs):
        return [L.patch(**kwargs) for L in self.edges()]

    def patch(self,**kwargs):
        verts = []
        for e in self.edges():
            verts += e.points_along()
        return patches.Polygon(verts,closed=True,**kwargs)

    def vertex_patches(self,**kwargs):
        return [v.patch(**kwargs) for v in self.vertices()]
