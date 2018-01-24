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
            fp, pathname, description = imp.find_module('_matrixNR', [dirname(__file__)])
        except ImportError:
            import _matrixNR
            return _matrixNR
        if fp is not None:
            try:
                _mod = imp.load_module('_matrixNR', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _matrixNR = swig_import_helper()
    del swig_import_helper
else:
    import _matrixNR
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
    return _matrixNR.mallocPy(n)
mallocPy = _matrixNR.mallocPy

def malloc22Py(r, c):
    return _matrixNR.malloc22Py(r, c)
malloc22Py = _matrixNR.malloc22Py

def conc(a, b, c, n, m, k, d):
    return _matrixNR.conc(a, b, c, n, m, k, d)
conc = _matrixNR.conc

def writetomem(x, i, f):
    return _matrixNR.writetomem(x, i, f)
writetomem = _matrixNR.writetomem

def readfrommem(x, i):
    return _matrixNR.readfrommem(x, i)
readfrommem = _matrixNR.readfrommem

def deallocPy(x):
    return _matrixNR.deallocPy(x)
deallocPy = _matrixNR.deallocPy

def writeto2mem(x, i, j, f):
    return _matrixNR.writeto2mem(x, i, j, f)
writeto2mem = _matrixNR.writeto2mem

def readfrom2Dmem(x, i, j):
    return _matrixNR.readfrom2Dmem(x, i, j)
readfrom2Dmem = _matrixNR.readfrom2Dmem

def mallocLongPy(n):
    return _matrixNR.mallocLongPy(n)
mallocLongPy = _matrixNR.mallocLongPy

def banmul(a, n, m1, m2, x0, b0):
    return _matrixNR.banmul(a, n, m1, m2, x0, b0)
banmul = _matrixNR.banmul

def bandec(a, n, m1, m2, al, indx0):
    return _matrixNR.bandec(a, n, m1, m2, al, indx0)
bandec = _matrixNR.bandec

def banbks(a, n, m1, m2, al, indx0, b0):
    return _matrixNR.banbks(a, n, m1, m2, al, indx0, b0)
banbks = _matrixNR.banbks

def getufromG(h, G, bed, hMbeg, hMend, wMbeg, wMend, GMbeg, GMend, uMbeg, uMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbhc, u, hhbc, whbc, Ghbc, bedhbc):
    return _matrixNR.getufromG(h, G, bed, hMbeg, hMend, wMbeg, wMend, GMbeg, GMend, uMbeg, uMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbhc, u, hhbc, whbc, Ghbc, bedhbc)
getufromG = _matrixNR.getufromG

def getufromGBAND(h, G, bed, hMbeg, hMend, wMbeg, wMend, GMbeg, GMend, uMbeg, uMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbhc, u, hhbc, whbc, Ghbc, bedhbc):
    return _matrixNR.getufromGBAND(h, G, bed, hMbeg, hMend, wMbeg, wMend, GMbeg, GMend, uMbeg, uMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbhc, u, hhbc, whbc, Ghbc, bedhbc)
getufromGBAND = _matrixNR.getufromGBAND

def dmatrix(nrl, nrh, ncl, nch):
    return _matrixNR.dmatrix(nrl, nrh, ncl, nch)
dmatrix = _matrixNR.dmatrix
# This file is compatible with both classic and new-style classes.

