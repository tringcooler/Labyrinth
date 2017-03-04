
def roprop(prop):
    def _mod(cls):
        setattr(cls, prop,
                property(lambda self: getattr(self, '_' + prop)))
        return cls
    return _mod

def neq(cls):
    def _ne(self, dst):
        return not self.__eq__(dst)
    cls.__ne__ = _ne
    return cls

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
