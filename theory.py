from __future__ import absolute_import
import theorydata
import integralequations
import numpy as np

from cmath import pi,exp,log

import logging
logger = logging.getLogger(__name__)

class theory:

    def __init__(self, theoryname, R = 1.0, nosymmetric = False, oper = False, splitrays = True, usebc = False, zerorayfunctions = False):
        """
        Constructs a new theory object from the given data.
        """

        self.theoryname = theoryname

        data = theorydata.getdata(theoryname)

        bpscounts = []

        # convenience function for building up the "bpscounts" list
        def addbps(charge, count):
            bpscounts.append({'charge': charge, 'count': count})
            anticharge = [-n for n in charge]
            bpscounts.append({'charge': anticharge, 'count': count})

        if "bpshypers" in data:
            for charge in data["bpshypers"]:
                addbps(charge, 1)
        if "bpsvectors" in data:
            for charge in data["bpsvectors"]:
                addbps(charge, -2)
        if "bpsgeneral" in data:
            for [charge,count] in data["bpsgeneral"]:
                addbps(charge, count)

        rank = data["Xrank"]

        try:
            ip = data["intersectionpairing"]
        except KeyError:
            ip = None

        if ip is None:
            # compute pairing on X-charges from B-matrix
            B = data["B"]
            Bunf = B[:rank]
            Xb = data["Xbasis"]
            pairing = np.dot(Xb,np.dot(Bunf,np.transpose(Xb)))
        else:
            # just get pairing directly
            pairing = ip

        sigmafunc = data["sigma"]

        self.rank = rank
        self.bpscounts = bpscounts
        self.pairing = pairing
        self.R = R
        self.oper = oper
        self.data = data
        self.sigmafunc = sigmafunc
        self.angles = data["angles"]
        try:
            self.Aangles = data["Aangles"]
        except KeyError:
            self.Aangles = None

        def zerofunction(gamma, hlist):
            return 0.0*hlist
        if usebc:
            try:
                self.boundaryfout = data["boundaryfout"]
            except KeyError:
                self.boundaryfout = zerofunction
        else:
            self.boundaryfout = zerofunction
        self.usebc = usebc

        self.zerorayfunctions = zerorayfunctions

        self.centralcharges = data["centralcharges"]

        self.symmetric = True
        if nosymmetric:
            self.symmetric = False
        for angle in self.angles:
            if abs(angle) > 1e-8:
                self.symmetric = False

        self.buildchargeequivclasses()

        self.splitrays = splitrays

        self.buildraylist()
        self.buildrayequivclasses()


    def buildraylist(self):
        if self.splitrays:
            self.buildsplitraylist()
        else:
            self.buildnonsplitraylist()

    def buildchargeequivclasses(self):
        """
        Build and store list of equivalence classes of charges.
        Charges in the same equivalence class should have the same function xinst.
        If self.symmetric is True, then each equivalence class contains 2 charges
        {gamma,-gamma}.
        If self.symmetric is False, then each equivalence class contains 1 charge.
        This could be extended for theories with extra symmetry, but so far that
        is not implemented.
        """
        def negate(charge):
            return [-n for n in charge]
        def firstpositive(charge):
            for n in charge:
                if n > 0:
                    return True
                if n < 0:
                    return False
            return None
        chargeequivclasses = []
        for item in self.bpscounts:
            charge = item['charge']
            if not self.symmetric:
                chargeequivclasses.append([charge])
            if self.symmetric:
                if firstpositive(charge):
                    chargeequivclasses.append([charge, negate(charge)])
        self.chargeequivclasses = chargeequivclasses

    def buildrayequivclasses(self):
        """
        Build and store list of equivalence classes of rays.
        """
        def negate(charge):
            return [-n for n in charge]
        rayidsused = []
        rayequivclasses = []
        for item in self.activerays:
            if item.id not in rayidsused:
                if not self.symmetric:
                    rayequivclasses.append([item.id])
                    rayidsused.append(item.id)
                else:
                    # find the opposite ray
                    foundray = False
                    for item2 in self.activerays:
                        if item.charges[0] == negate(item2.charges[0]):
                            rayequivclasses.append([item.id, item2.id])
                            rayidsused.append(item.id)
                            rayidsused.append(item2.id)
                            foundray = True
                            break
                    if not foundray:
                        raise ArithmeticError("Failed to find opposite for a ray")

        self.rayequivclasses = rayequivclasses

    def basischarges(self):
        """
        Convenience function for producing a list of basis charges;
        e.g. if rank = 3 the output list will be [[1,0,0],[0,1,0],[0,0,1]]
        """
        return [[(1 if i==j else 0) for i in range(self.rank)] for j in range(self.rank)]

    def buildsplitraylist(self):
        # build the ray-list in "split-rays" case, one single-charge ray for each BPS particle charge
        self.activerays = []

        # first aggregate BPS charges which are multiples of a given primitive charge

        def ismultiple(charge1, charge2):
            commonratio = None
            eps = 1e-9
            for k1,k2 in zip(charge1,charge2):
                if (k1 == 0 and k2 == 0):
                    continue
                if (k1 == 0 and k2 != 0) or (k1 != 0 and k2 == 0):
                    return False
                if float(k1)/k2 < 0:
                    return False
                if commonratio is None:
                    commonratio = float(k1)/k2
                if abs(float(k1)/k2 - commonratio) > eps:
                    return False
            return True

        rayspecs = []
        curid = 0
        for item in self.bpscounts:
            charge = item["charge"]
            count = item["count"]
            Z = self.Z(charge)
            phase = -Z/abs(Z)
            for rayspec in rayspecs:
                ratio = None
                chargeonray = rayspec["counts"][0][0]
                if ismultiple(charge, chargeonray):
                    rayspec["counts"] += [[charge,Z,count]]
                    break
            else:
                newrayspec = {"phase": phase, "counts": [[charge,Z,count]]}
                rayspecs.append(newrayspec)

        # now for each ray, detect the primitive charge and express others as multiples of that one
        for rayspec in rayspecs:
            minM = None
            counts = rayspec["counts"]
            for charge,Z,count in counts:
                M = abs(Z)
                if (minM is None) or (M < minM):
                    minM = M
                    primcharge = charge
            rayspec["primcharge"] = primcharge
            rayspec["redcounts"] = []
            for charge,Z,count in counts:
                multiple = int(abs(Z) / minM)
                rayspec["redcounts"].append([multiple,count])

        for curid,rayspec in enumerate(rayspecs, 0):
            phase = rayspec["phase"]
            counts = rayspec["redcounts"]
            primcharge = rayspec["primcharge"]

            # little gadget to construct the needed ray-functions
            def rayfunction(primcharge, counts):
                # the ray-function F_gamma in general depends on a list of charges and their corresponding xlists
                def rayf(gamma, charges, xlists):
                    # we just have one xlist for the ray, the one for X_\gamma with \gamma the primitive charge generating the rest
                    # and the ray-function depends only on that xlist
                    # TODO: check sign
                    runningtotal = 0.0
                    primxlist = xlists[0]
                    for multiple,count in counts:
                        charge = [multiple*k for k in primcharge]
                        sigma = self.sigma(charge)
                        contribution = count * self.DSZ(charge,gamma) * np.log(1 - sigma*np.exp(multiple*primxlist))
                        runningtotal += contribution
                    return runningtotal
                # miraculously this seems to work: the returned function
                # will remember the values of "primcharge" and "counts" from 
                # enclosing scope; buzzword is "Python closure"
                return rayf

            # similar gadget to construct the active-detector
            # TODO: deal with case where all counts are zero, in which case this should return false
            def rayactivedetector(primcharge):
                def rad(gamma):
                    return (self.DSZ(gamma,primcharge) != 0)
                return rad

            def zerofunction(gamma, tlist):
                return 0.0*tlist

            def zerorayf(gamma, charges, xlists):
                return 0.0*xlists[0]

            if not self.zerorayfunctions:
                currayfunction = rayfunction(primcharge, counts)
            else:
                currayfunction = zerorayf
            currayactivedetector = rayactivedetector(primcharge)

            curboundaryfin = zerofunction

            charges = [primcharge]
            currentray = integralequations.activeray(id = curid, phase = phase, function = currayfunction, activedetector = currayactivedetector, mutuallylocal = True, boundaryfin = curboundaryfin, charges = charges)
            self.activerays.append(currentray)

    def buildnonsplitraylist(self):
        # build ray-list in non-split-rays case
        self.activerays = []
        try: 
            raylist = self.data["rays"]
        except KeyError:
            raise NotImplementedError("Asked for non-split rays but no ray data is known for this theory")
        for curid,item in enumerate(raylist, 0):
            def truefunc(charge):
                return True
            currentray = integralequations.activeray(id = curid, phase = item["phase"], function = item["function"], activedetector = truefunc, mutuallylocal = False, charges = item["charges"])
            self.activerays.append(currentray)

    def getraybyid(self, rayid):
        for ray in self.activerays:
            if ray.id == rayid:
                return ray
        return None

    def chargeequivclassof(self, charge):
        """Find the equivalence class containing a given charge, or None if there is none.
        arguments: charge (list of ints)
        returns: equivclass (list of lists of ints) or None
        """
        for equivclass in self.chargeequivclasses:
            if charge in equivclass:
                return equivclass
        return None

    def rayequivclassof(self, rayid):
        for equivclass in self.rayequivclasses:
            if rayid in equivclass:
                return equivclass
        return None

    def bpscount(self, charge):
        """Get the BPS count for a given charge.
        arguments: charge (list of ints)
        returns: int
        """
        for item in self.bpscounts:
            if item["charge"] == charge:
                return item["count"]
        return 0

    def Z(self, charge):
        """
        Compute central charge Z for a given charge.
        argument: charge (list of ints)
        returns: float
        """
        return np.dot(charge, self.centralcharges)

    def angle(self, charge):
        """
        Compute angle theta for a given charge.
        argument: charge (list of ints)
        returns: float
        """
        return np.dot(charge, self.angles)

    def sigma(self, charge):
        """
        Compute sigma for a given charge. (pass-through to function given in
        theory data)
        """
        return self.sigmafunc(charge)

    def DSZ(self, charge1, charge2):
        """Compute the DSZ pairing between charges.
        arguments: charge1 (list of int), charge2 (list of int)
        returns: int
        """
        return np.dot([charge1], np.dot(self.pairing, np.transpose([charge2])))[0,0]

    def xsfonray(self, charge, t):
        """Compute xsf along BPS ray, as a function of t where zeta = -Z/|Z| exp(t)
        arguments: charge (list of int), t (float)
        returns: complex
        """
        if not self.oper:
            return -self.R*abs(self.Z(charge))*(exp(-t)+exp(t)) + 1j*self.angle(charge)
        else:
            zeta = -self.Z(charge)/abs(self.Z(charge))*exp(t)
            return -self.Z(charge)/zeta + 1j*self.angle(charge) + self.boundaryfout(charge, zeta)

    def xsf(self, charge, zeta):
        """Compute xsf, as a function of zeta
        arguments: charge (list of int), zeta (complex)
        returns: complex
        """
        if not self.oper:
            return self.R*self.Z(charge) / zeta + self.R*self.Z(charge).conjugate()*zeta + 1j*self.angle(charge)
        else:
            return self.Z(charge) / zeta + 1j*self.angle(charge) + self.boundaryfout(charge, zeta)

    def Xsf(self, charge, zeta):
        """Compute Xsf as a function of zeta (just pass-through to xsf)
        arguments: charge (list of int), zeta (complex)
        returns: complex
        """
        return self.sigma(charge)*exp(self.xsf(charge, zeta))

# TODO: make this function use self.oper?
    def xcorr(self, charge, nterms, zeta, oper = False):
        """
        Compute some part of the large R or small hbar expansion of x(zeta).

        Uses precomputed formulas for this expansion, which are not very well tested
        beyond the first term. When the saddle point is at a pole, we set the
        contribution to zero; this is probably right at leading order but maybe not
        beyond that.

        Arguments:
          charge (list of ints)
          nterms (int): number of terms in the expansion to take
          zeta (complex)
        Returns:
          complex
        """
        if not oper:
            if nterms > 3:
                raise NotImplementedError("Asked for > 3 terms in large R expansion")
            total = 0
            for item in self.bpscounts:
                bpscharge = item['charge']
                count = item['count']
                bpsphase = -self.Z(bpscharge) / abs(self.Z(bpscharge))
                weight = count * self.DSZ(charge, bpscharge)
                if weight != 0 and zeta != bpsphase:
                    M = self.R*abs(self.Z(bpscharge))
                    contributions = []
                    contributions.append(  ((bpsphase + zeta)/(bpsphase - zeta)) * (pi)**(0.5) * exp(-2*M) / (M)**(0.5) )
                    contributions.append( -((bpsphase + zeta)*(bpsphase**2 - 10*bpsphase*zeta + zeta**2)/(16*(bpsphase - zeta)**3)) * (pi)**(0.5) * exp(-2*M) / (M)**(1.5) )
                    contributions.append(  ((bpsphase + zeta)*3*(3*bpsphase**4 - 28*bpsphase**3*zeta + 178*bpsphase**2*zeta**2 - 28*bpsphase*zeta**3 + 3*zeta**4)/(512*(bpsphase - zeta)**5)) * (pi)**(0.5) * exp(-2*M) / (M)**(2.5) )
                    totalcontrib = sum(contributions[:nterms])
                    total += (weight / (4*pi*1j)) * totalcontrib
            return total+self.xsf(charge, zeta)
        if oper:
            if nterms > 1:
                raise NotImplementedError("Asked for > 1 terms in small hbar expansion for oper")
            total = 0
            if nterms != 0:
                if zeta != 1:
                    raise NotImplementedError("Asked for zeta != 1 in small hbar expansion for oper")
                for item in self.bpscounts:
                    bpscharge = item['charge']
                    count = item['count']
                    Z = self.Z(bpscharge)
                    weight = count * self.DSZ(charge, bpscharge)
                    if weight != 0:
                        contributions = []
                        if self.theoryname == "A1A2":
                            contributions.append( -2*pi**2/(15*Z) )
                        elif self.theoryname == "A2A1":
                            contributions.append( -2*pi**2/(15*Z) )
                        else:
                            raise NotImplementedError("Integrals giving the hbar expansion for oper only computed for A1A2 or A2A1")
                        totalcontrib = sum(contributions[:nterms])
                        total += (weight / (4*pi*1j)) * totalcontrib
            return total+self.xsf(charge, zeta)


    def Xcorr(self, charge, zeta, nterms, oper = False):
        """
        Compute part of the large R expansion of X(zeta). Just pass-through to xcorr;
        see xcorr for details of arguments.
        """
        return self.sigma(charge)*exp(self.xcorr(charge = charge, nterms = nterms, zeta = zeta, oper = oper))

    def xcorrestimate(self, charge):
        """
        A very rough estimate of the size of the first correction in the large R expansion of x,
        evaluated at zeta = -Z/|Z|.
        argument: charge (list of int)
        returns: float
        """
        total = 0
        zeta = -self.Z(charge) / abs(self.Z(charge))
        totalestimate = 0.0
        for item in self.bpscounts:
            bpscharge = item['charge']
            count = item['count']
            bpsphase = self.Z(bpscharge) / abs(self.Z(bpscharge))
            weight = count * self.DSZ(charge, bpscharge)
            if weight != 0 and zeta != bpsphase:
                M = self.R*abs(self.Z(bpscharge))
                contribution = ((bpsphase + zeta)/(bpsphase - zeta)) * (pi)**(0.5) * exp(-2*M) / (M)**(0.5)
                totalestimate += (abs(weight) / (4*pi)) * abs(contribution)
        return totalestimate

    def getallphases(self):
        """
        Get an ordered list of the phases of all BPS particles.
        """
        phases = []
        for item in self.bpscounts:
            charge = item['charge']
            phases.append(np.angle(self.Z(charge)))
        phases.sort()
        return phases

    def getallZ(self):
        """
        Get a list of the central charges of all BPS particles.
        """
        Zlist = []
        for item in self.bpscounts:
            charge = item['charge']
            Zlist.append(self.Z(charge))
        return Zlist

    def initialfunctions(self, charge, L, steps, phase):
        """
        Get initial lists tlist, xinstlist, xsflist, to go into an xarray at zeroth iteration.

        xintlist is identically 0.
        tlist is uniformly spaced list running from -L to L (this is the same for all iterations.)
        xsflist is pass-through to the function xsflist() (also the same for all iterations.)

        arguments: charge (list of ints), L (float), steps (int)
        returns: (tlist, xinstlist, xsflist)
          tlist: list of floats
          xinstlist: list of floats (all 0)
          xsflist: list of complex
        """
        tlist = np.linspace(-L, L, steps)
        xinstlist = np.zeros(steps, dtype = np.complex_)
        xsflist = self.xsflist(charge, L, steps, phase)
        return (tlist, xinstlist, xsflist)

    def xsflist(self, charge, L, steps, phase):
        """
        Get a list of values of xsf on BPS ray, sampled between t=L and t=-L.
        arguments: charge (list of ints), L (float), steps (int)
        returns: list of complex
        """
        tlist = np.linspace(-L, L, steps)
        if not self.oper:
            if phase is None:
                return -self.R*abs(self.Z(charge))*(np.exp(-tlist)+np.exp(tlist)) + 1j*self.angle(charge)
            else:
                return self.R*self.Z(charge)*np.exp(-tlist)/phase + self.R*((self.Z(charge)).conjugate())*np.exp(tlist)*phase + 1j*self.angle(charge)
        if self.oper:
            if phase is None:
                phase = -(abs(self.Z(charge))/self.Z(charge))
            hlist = phase*np.exp(tlist)
            if self.boundaryfout is None:
                return self.Z(charge)/hlist + 1j*self.angle(charge)
            else:
                return self.Z(charge)/hlist + 1j*self.angle(charge) + self.boundaryfout(charge, hlist)

    def bpscountsinrange(self, theta1, theta2):
        if theta2 < theta1:
            raise NotImplementedError("theta2 < theta1 not supported")
        def phasesordered(t1, p, t2):
            if t1 < p < t2:
                return True
            if t1 < p + (2*pi) < t2:
                return True
            if t1 < p - (2*pi) < t2:
                return True
            return False
        counts = [[item,cmath.phase(-self.Z(item["charge"]))] for item in self.bpscounts if phasesordered(theta1, cmath.phase(-self.Z(item["charge"])), theta2)]
        sortedcounts = sorted(counts, key = lambda x: x[1])
        return sortedcounts

    def KStransform(self, charge, Omega, xlist):
        # compute the new x's in order, one at a time
        xoutlist = []
        # compute the x for the BPS state
        xcharge = np.dot(charge, xlist)
        for n in range(self.rank):
            curcharge = [1 if i == n else 0 for i in range(self.rank)]
            newx = xlist[n] + self.DSZ(curcharge,charge)*Omega*log(1 - self.sigma(charge)*exp(xcharge))
            xoutlist.append(newx)
        return xoutlist

    # numerical spectrum generator, computed from BPS counts
    def spectrumgenerator(self, theta1, theta2, xlist):
        counts = self.bpscountsinrange(theta1, theta2)
        curxlist = xlist
        for item,phase in counts:
            curxlist = self.KStransform(item["charge"], item["count"], curxlist)
        return curxlist
