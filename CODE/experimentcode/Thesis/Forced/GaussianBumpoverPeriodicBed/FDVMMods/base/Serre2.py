# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Serre2', [dirname(__file__)])
        except ImportError:
            import _Serre2
            return _Serre2
        if fp is not None:
            try:
                _mod = imp.load_module('_Serre2', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Serre2 = swig_import_helper()
    del swig_import_helper
else:
    import _Serre2
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0



def mallocPy(n):
    return _Serre2.mallocPy(n)
mallocPy = _Serre2.mallocPy

def writetomem(x, i, f):
    return _Serre2.writetomem(x, i, f)
writetomem = _Serre2.writetomem

def readfrommem(x, i):
    return _Serre2.readfrommem(x, i)
readfrommem = _Serre2.readfrommem

def deallocPy(x):
    return _Serre2.deallocPy(x)
deallocPy = _Serre2.deallocPy

def minmod(a, b, c):
    return _Serre2.minmod(a, b, c)
minmod = _Serre2.minmod

def getufromG(h, G, bed, u0, u1, h0, h1, b0, b1, dx, n, u):
    return _Serre2.getufromG(h, G, bed, u0, u1, h0, h1, b0, b1, dx, n, u)
getufromG = _Serre2.getufromG

def edgevaluesSplit(h, G, bed, hMbeg, GMbeg, wMbeg, bMbeg, duMbeg, uEbeg, ddbCbeg, hMend, GMend, wMend, bMend, duMend, uEend, ddbCend, nMBC, nEBC, nCBC, n, nMbc, nEbc, nCbc, hMbc, GMbc, wMbc, bMbc, duMbc, uEbc, ddbCbc, dx, theta):
    return _Serre2.edgevaluesSplit(h, G, bed, hMbeg, GMbeg, wMbeg, bMbeg, duMbeg, uEbeg, ddbCbeg, hMend, GMend, wMend, bMend, duMend, uEend, ddbCend, nMBC, nEBC, nCBC, n, nMbc, nEbc, nCbc, hMbc, GMbc, wMbc, bMbc, duMbc, uEbc, ddbCbc, dx, theta)
edgevaluesSplit = _Serre2.edgevaluesSplit

def evolvewrapBC(h, G, bed, hMbeg, GMbeg, wMbeg, bMbeg, duEbeg, uEbeg, ddbCbeg, hMend, GMend, wMend, bMend, duEend, uEend, ddbCend, hMbeg1, GMbeg1, wMbeg1, duEbeg1, uEbeg1, hMend1, GMend1, wMend1, duEend1, uEend1, nMBC, nEBC, nCBC, n, nMbc, nEbc, nCbc, dx, dt, g, theta, x, t, a0, a1, a2, a3, a4, a5, a6, a7):
    return _Serre2.evolvewrapBC(h, G, bed, hMbeg, GMbeg, wMbeg, bMbeg, duEbeg, uEbeg, ddbCbeg, hMend, GMend, wMend, bMend, duEend, uEend, ddbCend, hMbeg1, GMbeg1, wMbeg1, duEbeg1, uEbeg1, hMend1, GMend1, wMend1, duEend1, uEend1, nMBC, nEBC, nCBC, n, nMbc, nEbc, nCbc, dx, dt, g, theta, x, t, a0, a1, a2, a3, a4, a5, a6, a7)
evolvewrapBC = _Serre2.evolvewrapBC

def HankEnergyall(x, h, u, b, g, n, nBC, dx):
    return _Serre2.HankEnergyall(x, h, u, b, g, n, nBC, dx)
HankEnergyall = _Serre2.HankEnergyall

def uhall(x, h, u, n, nBC, dx):
    return _Serre2.uhall(x, h, u, n, nBC, dx)
uhall = _Serre2.uhall

def hall(x, h, n, nBC, dx):
    return _Serre2.hall(x, h, n, nBC, dx)
hall = _Serre2.hall

def Gall(x, G, n, nBC, dx):
    return _Serre2.Gall(x, G, n, nBC, dx)
Gall = _Serre2.Gall
# This file is compatible with both classic and new-style classes.


