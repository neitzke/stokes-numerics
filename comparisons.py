from __future__ import absolute_import
from __future__ import print_function
import framedata
import integralequations
import theory
import namegen
import richardson
import numpy as np
from cmath import log,pi,phase

# Utility functions for describing differences

def phase_diff(x, y):
    '''Compute difference in phases between two complex numbers x,y,
       normalized to lie between -pi and pi'''
    a = phase(x)-phase(y)
    if -pi <= a <= pi:
        return a
    elif a < -pi:
        return a+2*pi
    elif a > pi:
        return a-2*pi

def abs_diff(x,y):
    return abs(x-y)

def rel_diff(x,y):
    return 2*abs(x-y) / (abs(x) + abs(y))

def log_diff(x,y):
    return (log(abs(x)) - log(abs(y))).real

def abs_to_rel(x,delta):
    '''Convert absolute error delta in quantity x to maximum corresponding relative error'''
    x = abs(x)
    return rel_diff(x-delta,x)


# get table of latest xars in a fixed theory for varying R
def xartable(theoryname, Rlist = None, tagxar = "", tagframes = "", theta = 0.0, oper = False):
    if Rlist is None:
        Rlist = theory.theory(theoryname=theoryname).data["Rlist"]

    outtable = []
    for R in Rlist:
        xar = integralequations.loadXar(theoryname, R, oper = oper)
        outtable.append([R,xar])

    return outtable

# get table of clusters in a fixed theory
def xarclustertable(theoryname, Rlist = None, tagxar = "", tagframes = "", theta = 0.0, oper = False):
    if Rlist is None:
        Rlist = theory.theory(theoryname=theoryname).data["Rlist"]

    outtable = []
    for R in Rlist:
        outtable.append({"R": R, "xarcluster": integralequations.loadXarCluster(theoryname, R, tag = tagxar, theta = theta, oper = oper)})

    return outtable

# get table of framedata in a fixed theory for varying R
def fdtable(theoryname, Rlist = None, tagxar = "", tagframes = "", theta = 0.0, oper = False, pde_nmesh = framedata.PDE_NMESH):
    if oper:
        raise NotImplementedError("fdtable doesn't work for oper")

    if Rlist is None:
        Rlist = theory.theory(theoryname=theoryname).data["Rlist"]

    outtable = []
    for R in Rlist:
        fd = framedata.loadFrameData(theoryname, R, theta = theta, oper = oper, pde_nmesh = pde_nmesh)
        outtable.append([R,fd])
    return outtable

# get table of clusters in a fixed theory
def fdclustertable(theoryname, Rlist = None, tagxar = "", tagframes = "", theta = 0.0, oper = False, pde_nmesh = framedata.PDE_NMESH):
    if oper:
        raise NotImplementedError("fdclustertable doesn't work for oper")

    if Rlist is None:
        Rlist = theory.theory(theoryname=theoryname).data["Rlist"]

    outtable = []
    for R in Rlist:
        fd = framedata.loadFrameData(theoryname, R, theta = theta, oper = oper, pde_nmesh = pde_nmesh)
        outtable.append({"R": R, "fdcluster": fd.getCluster()})

    return outtable

def compareClusters(theoryname, R=1.0, absh=1.0, theta=0, oper = False, scratch = False, make_monic = True, pde_nmesh = framedata.PDE_NMESH, steps = integralequations.IEQ_STEPS, L = integralequations.IEQ_L, damping = integralequations.IEQ_DAMPING, tolerance = integralequations.IEQ_TOLERANCE, rmax = None, fdmethod = "fourier"):
    '''Load framedata and integralequations data for same theory and return a
    dict with the clusters and comparison data.

    For Hitchin section (opers=False), pde_nmesh can be an integer or a
    sequence of three integers.  In the latter case compareClustersRichardson
    is called and a PDE error estimate will be returned.  Otherwise, the
    return object only contains an estimate of ODE error.'''

    if not oper:
        # Decide whether to do Richardson interpolation
        try:
            pde_nmesh = int(pde_nmesh)
            # integer, proceed
        except TypeError:
            # not integer, assume a list, call Richardson version
            return compareClustersRichardson(theoryname,R=R,theta=theta,oper=oper,scratch=scratch,make_monic=make_monic,pde_nmesh_list=pde_nmesh,steps=steps,L=L,rmax=rmax,fdmethod=fdmethod,damping=damping,tolerance=tolerance)

    if scratch:
        xar = integralequations.computeXar(theoryname, R = R, oper = oper, steps = steps, L = L, damping= damping, tolerance=tolerance)
        frames = framedata.computeFrames(theoryname, R = R, theta = theta, oper = oper, make_monic = make_monic, pde_nmesh = pde_nmesh, rmax = rmax, method = fdmethod, absh = absh)
        xarcluster = xar.getCluster(theta = theta, absh = absh)
        fdcluster = frames.getCluster()
    else:
        xarcluster = integralequations.loadXarCluster(theoryname = theoryname, R = R, oper = oper, absh = absh, theta = theta)
        frames = framedata.loadFrameData(theoryname = theoryname, R = R, oper = oper, absh = absh, theta = theta, pde_nmesh = pde_nmesh)
        fdcluster = frames.getCluster()

    sfcluster = integralequations.getApproxCluster(theoryname, R = R, absh = absh, theta = theta, oper = oper)

    absdiff = [ abs_diff(x,f) for x,f in zip(xarcluster, fdcluster) ]
    logdiff   = [ log_diff(x,f)   for x,f in zip(xarcluster, fdcluster) ]
    phasediff = [ phase_diff(x,f) for x,f in zip(xarcluster, fdcluster) ]
    reldiff   = [ rel_diff(x,f)   for x,f in zip(xarcluster, fdcluster) ]

    errest_absode = list(frames.estimateODEclusterXerror())
    errest_relode = [ abs_to_rel(x,delta) for x,delta in zip(fdcluster, errest_absode) ]

    return {
        "xarcluster": xarcluster,
        "fdcluster": fdcluster,
        "sfcluster": sfcluster,
        "absdiff" : absdiff,
        "logdiff": logdiff,
        "phasediff": phasediff,
        "reldiff": reldiff,
        "frames": frames,
        "errest": {
            "absode": errest_absode,
            "relode": errest_relode,
        },
    }

def compareClustersRichardson(theoryname, R=1.0, absh=1.0, theta=0, oper = False, scratch = False, make_monic = True, pde_nmesh_list = framedata.PDE_NMESH_LIST, steps = integralequations.IEQ_STEPS, L = integralequations.IEQ_L, damping = integralequations.IEQ_DAMPING, tolerance = integralequations.IEQ_TOLERANCE, rmax = None, fdmethod = "fourier"):
    '''Load framedata (for several pde_nmesh values) and integralequations
    data for same theory and return a dict with the clusters and comparison data.

    For clusters from the framedata, Richardson interpolation is used to get
    an empirical PDE error estimate.'''

    if oper:
        raise NotImplementedError("compareClustersRichardson is only for the Hitchin section")
    if len(pde_nmesh_list) != 3:
        raise NotImplementedError("pde_nmesh_list must be a list of length 3")

    pde_nmesh_list.sort(reverse=True)

    if scratch:
        xar = integralequations.computeXar(theoryname, R = R, oper = oper, steps = steps, L = L, damping = damping, tolerance = tolerance)
        frames_list = [ framedata.computeFrames(theoryname, R = R, theta = theta, oper = oper, make_monic = make_monic, pde_nmesh = N, rmax = rmax, method = fdmethod, absh = absh) for N in pde_nmesh_list ]
        xarcluster = xar.getCluster(theta = theta, absh = absh)
        fdcluster_list = [ f.getCluster() for f in frames_list ]
    else:
        xarcluster = integralequations.loadXarCluster(theoryname = theoryname, R = R, oper = oper, absh = absh, theta = theta)
        frames_list = [ framedata.loadFrameData(theoryname = theoryname, R = R, oper = oper, absh = absh, theta = theta, pde_nmesh = N) for N in pde_nmesh_list ]
        fdcluster_list = [ f.getCluster() for f in frames_list ]

    sfcluster = integralequations.getApproxCluster(theoryname, R = R, absh = absh, theta = theta, oper = oper)
    frames = frames_list[0]
    fdcluster = fdcluster_list[0] # Highest pde_nmesh value is primary return value

    absdiff = [ abs_diff(x,f) for x,f in zip(xarcluster, fdcluster) ]
    logdiff   = [ log_diff(x,f)   for x,f in zip(xarcluster, fdcluster) ]
    phasediff = [ phase_diff(x,f) for x,f in zip(xarcluster, fdcluster) ]
    reldiff   = [ rel_diff(x,f)   for x,f in zip(xarcluster, fdcluster) ]

    dxlist = [ f.params['rmax']*2.0/f.params['pde_nmesh'] for f in frames_list ]

    fdcluster_extrap = []
    exponents = []
    errest_abspde = []
    for row in np.array(fdcluster_list).transpose():
        try:
            extrap = richardson.extrapolate(row,dxlist,order=None,order_seed=2.0)
            fdcluster_extrap.append(extrap['extrapolated'])
            exponents.append(extrap['extrapolation_order'])
            errest_abspde.append(2.0*extrap['delta'])
        except Exception:
            fdcluster_extrap.append(float('nan'))
            exponents.append(float('nan'))
            errest_abspde.append(float('nan'))

    errest_relpde = [ abs_to_rel(x,e) for x,e in zip(fdcluster,errest_abspde) ]

    errest_absode = list(frames.estimateODEclusterXerror())
    errest_relode = [ abs_to_rel(x,delta) for x,delta in zip(fdcluster, errest_absode) ]

    return {
        "xarcluster": xarcluster,
        "fdcluster": fdcluster_list[0],
        "sfcluster": sfcluster,
        "absdiff" : absdiff,
        "logdiff": logdiff,
        "phasediff": phasediff,
        "reldiff": reldiff,
        "frames": frames_list[0],
        "fdrichardson": {
            "pde_nmesh": pde_nmesh_list,
            "dx": dxlist,
            "computed_cluster": fdcluster_list,
            "extrapolated_cluster": fdcluster_extrap,
            "observed_exponent": exponents,
            "model_exponent": [ 2.0 for _ in exponents ],
            "frames": frames_list,
            },
        "errest": {
            "absode": errest_absode,
            "relode": errest_relode,
            "abspde": errest_abspde,
            "relpde": errest_relpde,
        },
    }


# Return data comparing integral equation and cross-ratio results
# for a given theory and list of R values.
def rawtables(theoryname, Rlist=None, theta=0, tagxar = "" ,tagframes = "", method = "simps", oper = False, pde_nmesh = framedata.PDE_NMESH):
    if oper:
        raise NotImplementedError("rawtables doesn't work for oper")

    thd = theory.theory(theoryname=theoryname).data

    if Rlist is None:
        Rlist = thd["Rlist"]

    rank = thd["Xrank"]
    names = thd["Xclusternames"]
    # include phase errors in data table in oper case, else don't
    if oper:
        dataheaders = ["R"] + sum([[name+"f", name+"i", name+"logdiff", name+"phasediff"] for name in names], [])
    else:
        dataheaders = ["R"] + sum([[name+"f", name+"i", name+"logdiff"] for name in names], [])
    absdiffheaders = ["R"] + [name+"absdiff" for name in names]
    logdiffheaders = ["R"] + [name+"logdiff" for name in names]
    reldiffheaders = ["R"] + [name+"reldiff" for name in names]
    phasediffheaders = ["R"] + [name+"phasediff" for name in names]

    if oper:
        Xtype = "complex"
    else:
        Xtype = "real"
    datatypes = ["real"] + sum([[Xtype, Xtype, "real", "real"] for name in names], [])
    absdifftypes = ["real"] + ["real" for name in names]
    logdifftypes = ["real"] + ["real" for name in names]
    reldifftypes = ["real"] + ["real" for name in names]
    phasedifftypes = ["real"] + ["real" for name in names]

    datarows = []
    absdiffrows = []
    logdiffrows = []
    phasediffrows = []
    reldiffrows = []

    for R in Rlist:
        iecluster = integralequations.loadXarCluster(theoryname = theoryname, R = R, theta = theta, oper = oper)
        fdcluster = framedata.loadFrameDataCluster(theoryname=theoryname,R=R,theta=theta,tag=tagframes,pde_nmesh=pde_nmesh,oper=oper)
        absdiff   = [ abs_diff(x,f)   for x,f in zip(iecluster, fdcluster) ]
        logdiff   = [ log_diff(x,f)   for x,f in zip(iecluster, fdcluster) ]
        phasediff = [ phase_diff(x,f) for x,f in zip(iecluster, fdcluster) ]
        reldiff   = [ rel_diff(x,f)   for x,f in zip(iecluster, fdcluster) ]
        
        data_row = [R]
        absdiff_row = [R]
        logdiff_row = [R]
        phasediff_row = [R]
        reldiff_row = [R]
        for item in range(rank):
            if oper:
                data_row += [fdcluster[item], iecluster[item], logdiff[item], phasediff[item]]
            else:
                data_row += [fdcluster[item], iecluster[item], logdiff[item]]

            absdiff_row.append(absdiff[item])
            logdiff_row.append(logdiff[item])
            phasediff_row.append(phasediff[item])
            reldiff_row.append(reldiff[item])
        datarows.append(data_row)
        absdiffrows.append(absdiff_row)
        logdiffrows.append(logdiff_row)
        phasediffrows.append(phasediff_row)
        reldiffrows.append(reldiff_row)

    dataout = {"rows": datarows, "headers": dataheaders, "types": datatypes}
    absdiffout = {"rows": absdiffrows, "headers": absdiffheaders, "types": absdifftypes}
    logdiffout = {"rows": logdiffrows, "headers": logdiffheaders, "types": logdifftypes}
    phasediffout = {"rows": phasediffrows, "headers": phasediffheaders, "types": phasedifftypes}
    reldiffout = {"rows": reldiffrows, "headers": reldiffheaders, "types": reldifftypes}

    return {'data': dataout, 'absdiff': absdiffout, 'logdiff': logdiffout, 'phasediff': phasediffout, 'reldiff': reldiffout}

# Comparison tables
def plaintables(theoryname, Rlist=None, theta=0, precision=4, pde_nmesh=framedata.PDE_NMESH, tagxar = "", tagframes = "", method = "simps", oper=False):
    T = rawtables(theoryname,Rlist,theta,tagxar,tagframes,method, oper=oper, pde_nmesh = pde_nmesh)

    header_format = {"real": "{:^%ds}" % (precision+10), "complex": "{:^%ds}" % (2*precision+10)}
    data_format = {"real": "{:^%d.%dE}" % (precision+8,precision), "complex": "{:^%d.%dE}" % (precision+8,precision)}
    item_padlength = {"real": 11+precision, "complex": 19+2*precision}

    tables = {}
    for key, v in T.items():
        table = ""
        table += "".join([header_format[datatype].format(header).ljust(item_padlength[datatype]) for datatype, header in zip(v["types"], v["headers"]) ]) + "\n"
        for row in v["rows"]:
            table += "".join([data_format[datatype].format(item).ljust(item_padlength[datatype]) for datatype, item in zip(v["types"], row) ]) + "\n"
        tables[key] = table

    return tables

def texify(varname):
    if varname[0] == "X" and varname[1].isdigit() and len(varname) == 2:
        return "$X_%s$" % varname[1]
    elif varname[0] == "X" and varname[1].isdigit() and len(varname) == 3:
        return "$X_%s^{\\mathrm{%s}}$" % (varname[1],varname[2])
    elif varname[0:9] == "logdiff(X":
        return "$\\mathrm{logdiff}(X_%s)$" % varname[8]
    elif varname[2:9] == "logdiff":
        return "$\\mathrm{logdiff}(X_%s)$" % varname[1]
    elif varname[0:9] == "reldiff(X":
        return "$\\mathrm{reldiff}(X_%s)$" % varname[8]
    elif varname[2:9] == "reldiff":
        return "$\\mathrm{reldiff}(X_%s)$" % varname[1]
    elif varname[2:11] == "phasediff":
        return "$\\mathrm{phasediff}(X_%s)$" % varname[1]
    elif varname == "R":
        return "$R$"
    elif varname == "absh":
        return "$|\\hbar|$"
    elif varname == "abshi":
        return "$|\\hbar|^{-1}$"
    elif varname == "rmax":
        return "rmax"
    elif varname == "A1A2":
        return "$(A_1,A_2)$"
    elif varname == "A1A3":
        return "$(A_1,A_3)$"
    elif varname == "A2A1":
        return "$(A_2,A_1)$"
    elif varname == "A2A2":
        return "$(A_2,A_2)$"
    else:
        raise ValueError("Don't know how to texify %s" % varname)

def latextables(theoryname, Rlist=None, theta=0, precision=4, tagxar = "", tagframes = "", method = "simps", oper=False):
    T = rawtables(theoryname,Rlist,theta,tagxar,tagframes,method, oper=oper)

    header_format = {"real": "{:^%ds}" % (precision+10), "complex": "{:^%ds}" % (2*precision+10)}
    data_format = {"real": "${:^%d.%dE}$" % (precision+8,precision), "complex": "${:^%d.%dE}$" % (precision+8,precision)}
    item_padlength = {"real": 11+precision, "complex": 19+2*precision}

    tables = {}
    for key, v in T.items():
        table = ""
        table += "\\begin{tabular}{%s}\n" % ("|"+"".join("c|"*len(v["headers"])))
        table += "\\hline\n"
        table += " & ".join([header_format[datatype].format(texlabels.texify(header)).ljust(item_padlength[datatype]) for datatype, header in zip(v["types"], v["headers"]) ]) + "\\\\ \n"
        table += "\\hline\n"
        for row in v["rows"]:
            table += " & ".join([data_format[datatype].format(item).replace("E-", "\\mathrm{e-}").replace("E+", "\\mathrm{e+}") for datatype, item in zip(v["types"], row) ]) + "\\\\ \n"
        table += "\\hline\n"
        table += "\\end{tabular}\n"
        tables[key] = table

    return tables
