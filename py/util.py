
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

def lazypropds(hndl):
    nm = '_' + hndl.__name__
    dnm = nm + '_dirty'
    def _getter(self):
        if not hasattr(self, nm) or getattr(self, dnm):
            setattr(self, nm, hndl(self))
            setattr(self, dnm, False)
        return getattr(self, nm)
    def _setter(self, v):
        setattr(self, dnm, True)
    return property(_getter, _setter)

def lazypropdd(hndl):
    nm = '_' + hndl.__name__
    def _getter(self):
        if not hasattr(self, nm):
            setattr(self, nm, hndl(self))
        return getattr(self, nm)
    def _setter(self, v):
        if hasattr(self, nm):
            delattr(self, nm)
    return property(_getter, _setter)

def lazypropdh(hndl):
    nm = '_' + hndl.__name__
    dnm = hndl.__name__ + '_dirty'
    def _getter(self):
        if not hasattr(self, nm) or (
            hasattr(self, dnm) and getattr(self, dnm)()):
            setattr(self, nm, hndl(self))
        return getattr(self, nm)
    return property(_getter)

def iseq(hndl):
    def _iseq(self, dst):
        if self is dst:
            return True
        else:
            return hndl(self, dst)
    return _iseq

def lazyeq(hndl):
    def _lazyeq(self, dst):
        if self is dst:
            return True
        if not type(dst).__eq__.__name__ == '_lazyeq':
            return hndl(self, dst)
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
    return _lazyeq

@neq
class vchain(object):

    __succ = None

    def __getattr__(self, a):
        if self.__succ is None:
            raise AttributeError(
                "vchain has no attribute '{0}'".format(a))
        else:
            #return getattr(self.__succ, a)
            node = self
            while not node.__succ is None:
                node = node.__succ
            return getattr(node, a)

    def __setattr__(self, a, v):
        if self.__succ is None:
            super(vchain, self).__setattr__(a, v)
        else:
            #self.__succ.__setattr__(a, v)
            node = self
            while not node.__succ is None:
                node = node.__succ
            return setattr(node, a, v)

    def __eq__(self, dest):
        if not type(self) == type(dest):
            return False
        if self is dest:
            return True
        if self.__succ is None and dest.__succ is None:
            return self is dest
        snode = self
        while not snode.__succ is None:
            snode = snode.__succ
        dnode = dest
        while not dnode.__succ is None:
            dnode = dnode.__succ
        return snode == dnode

    def __nonzero__(self):
        return True

    def vlink(self, dst):
        if self.__succ is None:
            node = dst
            while not node.__succ is None:
                node = node.__succ
            if not node is self:
                self.__dict__ = {}
                self.__succ = dst
        else:
            node = self
            while not node.__succ is None:
                node = node.__succ
            if not node == dst:
                node.vlink(dst)
