"""Generate systematic data filenames following conventions of existing collection"""

from __future__ import absolute_import
import fileinfo
import os
import numpy as np

def ie_filename(theoryname, R = 1.0, oper=False, tag="", fullpath = False, clusteronly = False, theta = 0.0, absh=1.0):
    if fullpath:
        dirname = fileinfo.XARPATH
    else:
        dirname = ""
    if oper:
        if clusteronly:
            ending = "xarcluster"
            hstring = "absh%0.8fth%0.4f" % (absh,theta)
        else:
            ending = "xar"
            hstring = ""
        basename = theoryname + hstring + tag + "oper" + "." + ending
    else:
        if clusteronly:
            ending = "xarcluster"
            anglestring = "th%0.4f" % theta
        else:
            ending = "xar"
            anglestring = ""
        basename = theoryname + ("R%0.10f" % R) + anglestring + tag + "." + ending
    return os.path.join(dirname,basename)

def fd_filename(theoryname, R=1.0, theta=0.0, oper=False, tag = "", fullpath = False, absh = 1.0, pde_nmesh=None):
    if fullpath:
        dirname= fileinfo.FRAMEPATH
    else:
        dirname = ""
    if oper:
        basename = theoryname + ("absh%0.8fth%0.4f" % (absh,theta)) + tag + "oper.frames"
    else:
        basename = theoryname + ("N%dR%0.10fth%0.4f" % (pde_nmesh,R,theta)) + tag + ".frames"
    return os.path.join(dirname, basename)

def ieq_metric_filename(c, Lambda, fullpath = False):
    if fullpath:
        dirname= fileinfo.HKMETRICPATH
    else:
        dirname = ""
    basename = ("c%0.8fLambda%0.8f" % (c,Lambda)) + ".ieqmetric"
    return os.path.join(dirname, basename)

def fd_metric_filename(c, Lambda, fullpath = False):
    if fullpath:
        dirname= fileinfo.HKMETRICPATH
    else:
        dirname = ""
    basename = ("c%0.8fLambda%0.8f" % (c,Lambda)) + ".fdmetric"
    return os.path.join(dirname, basename)

# Is not currently used; usually we need both the check and the filename
def ie_file_exists(theoryname, R, oper=False, tag="", clusteronly = False, theta = 0.0, absh=1.0):
    return os.path.isfile(ie_filename(theoryname, R, oper=oper, tag=tag, fullpath=True, clusteronly=clusteronly, theta = theta, absh = absh))

# Is not currently used; usually we need both the check and the filename
def fd_file_exists(theoryname, R, theta = 0.0, oper=False, tag="", absh=1.0, pde_nmesh=None):
    return os.path.isfile(fd_filename(theoryname, R, oper=oper, theta=theta, tag=tag, fullpath=True, absh=absh, ode_thresh=ode_thresh, pde_nmesh=pde_nmesh))
