from __future__ import absolute_import
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import integralequations
import numpy as np

def plotData(datalist, part = "real", progressive = True, color = None, clip = False, tcutoff = None):
    """Plot real or imaginary parts of a given list of functions.
    arguments:
      datalist (list of tuples, each tuple of form (xlist,ylist)): data to plot;
        xlist should be real numbers, ylist can be complex
      part (str): "real" or "imag"; determines which part of ylist we use
      progressive (bool): if True, alphalevel starts at 0 and increases linearly
        as we go through datalist
      color (str): a color supported by axes.plot; if None then gets set to "blue" or "red"
        for real/imag parts respectively; if "rainbow" then color varies from blue to green
        as we go through datalist
      clip (bool): remove first and last items before plotting"""

    fig = plt.figure()
    axes = fig.add_subplot(1,1,1)
    alphalevel = 0.0

    if color is None:
        if part == "real":
            color = "blue"
        if part == "imag":
            color = "red"

    if color == "rainbow":
        red = 0.5
        green = 0.0
        blue = 1.0

    # iterate over list of functions
    n = len(datalist)
    for k,data in enumerate(datalist):
        if tcutoff is not None:
            tlist,zlist = zip(*[[t,z] for t,z in zip(*data) if abs(t) <= tcutoff])
        else:
            tlist,zlist = data

        if clip:
            tlist = tlist[1:-1]
            zlist = zlist[1:-1]

        if part == "real":
            ylist = [z.real for z in zlist]
        if part == "imag":
            ylist = [z.imag for z in zlist]

        # if "progressive" is set then vary alphalevel as we go through datalist
        if progressive:
            alphalevel += 1.0 / n
            alphalevel = min(alphalevel, 1.0)
        else:
            alphalevel = 1.0

        if color == "rainbow":
            hue = float(k) / n
            currentcolor = colors.hsv_to_rgb([hue, 0.8, 0.8])
        else:
            currentcolor = color

        # now plot the points
        axes.plot(tlist, ylist, color=currentcolor, alpha=alphalevel)

    return fig


def xarrayrayplotF(self, rayid, gamma, part="real"):
    datalist = [[self.tlist(rayid = rayid), self.getrayfunction(rayid, gamma)]]
    fig = plotData(datalist, part=part, color=None, progressive = False)
    return fig

def xarrayrayplots(self, rayid, part="real", addsf=False, tcutoff=None):
    # make ray plots for all charges attached to the ray
    ray = self.getRaydatum(rayid)

    figs = []
    for charge,xinstlist,xsflist in zip(ray.charges, ray.xinstlists, ray.xsflists):
        if addsf:
            datalist = [[ray.tlist, xinstlist+xsflist]]
        else:
            datalist = [[ray.tlist, xinstlist]]
        fig = plotData(datalist, part=part, color=None, progressive=False, tcutoff=tcutoff)
        figs.append(fig)

    return figs

# Warning: Monkey patching the xarray class!
integralequations.xarray.rayplotF = xarrayrayplotF
integralequations.xarray.rayplots = xarrayrayplots

