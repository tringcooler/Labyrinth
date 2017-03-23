
def roprop(prop):
    def _mod(cls):
        setattr(cls, prop,
                property(lambda self: getattr(self, '_' + prop)))
        return cls
    return _mod

def neq(cls):
    def _ne(self, dst):
        return not self == dst
    cls.__ne__ = _ne
    return cls

def lazyprop(hndl):
    nm = '_' + hndl.__name__
    def _getter(self):
        if not hasattr(self, nm):
            setattr(self, nm, hndl(self))
        return getattr(self, nm)
    return property(_getter)

def iseq(hndl):
    def _eq(self, dst):
        if self is dst:
            return True
        else:
            return hndl(self, dst)
    return _eq

def lazyeq(hndl):
    def _eq(self, dst):
        if self is dst:
            return True
        if not hasattr(self, '_eq_vchain'):
            self._eq_vchain = vchain()
        if not hasattr(dst, '_eq_vchain'):
            dst._eq_vchain = vchain()
        if self._eq_vchain == dst._eq_vchain:
            return True
        else:
            rslt = hndl(self, dst)
            if rslt:
                self._eq_vchain.vlink(dst._eq_vchain)
            return rslt
    return _eq

@neq
class vchain(object):

    __succ = None

    def __getattr__(self, a):
        if self.__succ is None:
            raise AttributeError(
                "vchain has no attribute '{0}'".format(a))
        else:
            return getattr(self.__succ, a)

    def __setattr__(self, a, v):
        if self.__succ is None:
            super(vchain, self).__setattr__(a, v)
        else:
            self.__succ.__setattr__(a, v)

    def __eq__(self, dest):
        if not type(self) == type(dest):
            return False
        if self.__succ is None and dest.__succ is None:
            return self is dest
        if not self.__succ is None:
            self = self.__succ
        if not dest.__succ is None:
            dest = dest.__succ
        return self == dest

    def __nonzero__(self):
        return True

    def vlink(self, dst):
        if self.__succ is None:
            self.__dict__ = {}
            self.__succ = dst
        else:
            if not self.__succ == dst:
                self.__succ.vlink(dst)
