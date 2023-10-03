from numpy import *
import json, gzip
from collections import Sequence
from GF2 import setGF2, multGF2
from weights import get_weight, get_weight_params

class Lattice(Sequence):
    pass

def quadrature(lattice, f, transf = None):
    """
    Evaluation of equal-weight quadrature rule for the given lattice and function f.
    lattice: Lattice point set.
    f: integrand function
    transf: point transformation to be applied to x in [0,1]
    """
    S = 0.
    if transf == None:
        for x in lattice: S += f(x)
    else:
        for x in lattice: S += f(transf(x))
    return S / len(lattice)

def vm_expansion(x,P,m):
    """
    Truncated formal Laurent series of x mod P.
    returns first m terms of formal Laurent series (coeffs of x^{-i}, i>0).
    """
    # convert integers x,P to coefficient vectors of polynomials in Z_2
    xp = format(x,'0%db'%m)
    Pp = format(P,'0%db'%(m+1))
    xp = [int(x) for x in xp[::-1]]
    Pp = [int(x) for x in Pp[::-1]]
    # start loop
    w = xp[::-1]
    for l in range(m):
        for j in range(l):
            w[l] ^= w[j]*Pp[m-l+j]
    # result is w[0]*1/2 + w[1]*1/4 + ...
    return w

def interlace_expansion(XI, alpha, m, b=2):
    """Interlaces the 'alpha' points stored in 'X' to one point."""
    assert(len(XI) == alpha)
    # XI: vector of alpha point expansions (coefficients of x^{-i}, i\ge1
    s = 0.
    for a in range(m):
        for j in range(1,alpha+1):
            s += XI[j-1][a] * b**(-j-a*alpha)
    return s

class PolynomialLattice(Lattice):
    def __init__(self,filename=None,smax=1):
        self.is_loaded = False
        self.smax = smax
        if filename != None: self.load(filename)

    def set_smax(self,smax):
        if self.is_loaded:
            assert(smax <= self.s)
        self.smax = smax

    def load(self,f):
        """Load generating vector from a file"""
        if f[-3:] == ".gz":
            with gzip.GzipFile(f,'r') as ff:
                data = json.loads(ff.read().decode('utf8'))
        else:
            with open(f) as ff:
                data = json.load(ff)
        self._load_helper(data)

    @staticmethod
    def loads(genvec_str):
        """Load generating vector directly from a string"""
        data = json.loads(genvec_str)
        lat = PolynomialLattice(smax=data["s"])
        lat._load_helper(data)
        return lat

    def _load_helper(self,data):
        self.m = int(data["m"])
        self.s = int(data["s"])
        if self.s < self.smax:
            raise Exception("Generating vector not long enough! (max. dim: s=%d)"%self.s)
        self.C = data["C"]
        self.P = data["P(x)"]
        self.P2 = int(data["P(2)"])
        self.WCE = data["WCE"]
        self.alpha = int(data["alpha"])
        self.entropy = data["entropy"]
        self.genvec = array(data["genvec"],dtype=int)
        self.weight = get_weight(data["weight"],get_weight_params(data))
        setGF2(self.m,self.P2)
        self.is_loaded = True

    def __getitem__(self,n):
        """need this for Sequence"""
        if not self.is_loaded:
            raise Exception("Must load a lattice before accessing its points!")
        if n >= 2**self.m:
            raise IndexError("Point set only contains 2**%d points."%self.m)
        return array([interlace_expansion([vm_expansion(multGF2(n,self.genvec[j*self.alpha+k]),self.P2,self.m) for k in range(self.alpha)],self.alpha,self.m) for j in range(self.smax)])

    def __len__(self):
        """Number of points in the lattice"""
        # need this for 'Sequence'
        return 2**self.m

class StandardLattice(Lattice):
    def __init__(self,filename=None,smax=1):
        self.is_loaded = False
        self.smax = smax
        if filename != None: self.load(filename)

    def set_smax(self,smax):
        if self.is_loaded:
            assert(smax <= self.s)
        self.smax = smax

    def load(self,f):
        """Load generating vector from a file"""
        if f[-3:] == ".gz":
            with gzip.GzipFile(f,'r') as ff:
                data = json.loads(ff.read().decode('utf8'))
        else:
            with open(f) as ff:
                data = json.load(ff)
        self._load_helper(data)

    @staticmethod
    def loads(genvec_str):
        """Load generating vector directly from a string"""
        data = json.loads(genvec_str)
        lat = StandardLattice(smax=data["s"])
        lat._load_helper(data)
        return lat

    def _load_helper(self,data):
        self.m = int(data["m"])
        self.s = int(data["s"])
        if self.s < self.smax:
            raise Exception("Generating vector not long enough! (max. dim: s=%d)"%self.s)
        self.alpha = int(data["alpha"])
        self.a1 = data["a1"]
        self.a2 = data["a2"]
        self.a3 = data["a3"]
        self.d1 = data["d1"]
        self.p = data["p"]
        self.WCE = data["WCE"]
        #self.entropy = data["entropy"]
        self.genvec = array(data["genvec"],dtype=int)
        #self.weight = get_weight(data["weight"],get_weight_params(data))
        self.is_loaded = True

    def __getitem__(self,n):
        """need this for Sequence"""
        if not self.is_loaded:
            raise Exception("Must load a lattice before accessing its points!")
        if n >= 2**self.m:
            raise IndexError("Point set only contains p**%d points."%(self.p,self.m))
        N = 2.**self.m
        y = array([mod(n*self.genvec[j]/N,1) for j in range(self.smax)])
        assert(all(y >= 0))
        assert(all(y <= 1))
        return y

    def __len__(self):
        """Number of points in the lattice"""
        # need this for 'Sequence'
        return 2**self.m
