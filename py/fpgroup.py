
import re

class grp_word_basis(object):

    def __init__(self, names):
        elms = []
        for n in names:
            if type(n) in [int, long]:
                raise TypeError("name shouldn't be int")
            elms.append(n)
            elms.append(self._invers_pres(n))
        self._elms = elms

    @property
    def elms(self):
        return self._elms

    def _ident_pres(self):
        return '<id>'

    def _invers_pres(self, n):
        return n + '**-1'

    def _invers_invers_pres(self, n):
        return n[:-4]

    def is_invers(self, n):
        return len(n) > 4 and n[-4:] == '**-1'

    def position(self, n):
        return self.elms.index(n)

    def invers(self, n):
        if self.is_invers(n):
            r = self.position(self._invers_invers_pres(n))
        else:
            r = self.position(self._invers_pres(n))
        return self.elms[r]

    def __eq__(self, dst):
        return self.elms == dst.elms

    def __ne__(self, dst):
        return not self.__eq__(dst)

    def __add__(self, dst):
        if self == dst:
            return self
        else:
            elms = list(self.elms)
            for i in dst.elms:
                if not i in elms:
                    elms.append(i)
            r = type(self)([])
            r._elms = elms
            return r

    def gen(self, n):
        if n == self._ident_pres():
            return grp_word([], self)
        else:
            try:
                self.position(n)
            except:
                raise
            return grp_word([n], self)

    def gens(self, inv = False):
        if not hasattr(self, '_gens'):
            self._gens = []
            self._inv_gens = []
            for elm in self.elms:
                if self.is_invers(elm):
                    self._inv_gens.append(self.gen(elm))
                else:
                    self._gens.append(self.gen(elm))
        if inv:
            return self._gens + self._inv_gens
        else:
            return list(self._gens)

class grp_word(object):

    def __init__(self, wl, basis):
        self.basis = basis
        self.seq = wl

    def exp_resp(self):
        wr = []
        cr = []
        cur = None
        cnt = 0
        def _upd():
            if self.basis.is_invers(cur):
                _cur = self.basis.invers(cur)
                _cnt = -cnt
            else:
                _cur = cur
                _cnt = cnt
            wr.append(_cur)
            cr.append(_cnt)
        for w in self.seq:
            if not cur == w:
                if not cur == None:
                    _upd()
                cur = w
                cnt = 0
            cnt += 1
        _upd()
        return wr, cr

    def mapped(self, gens):
        pass

    def __mul__(self, dst):
        return type(self)(self.seq + dst.seq, self.basis + dst.basis)

    def __pow__(self, dst):
        if type(dst) in [int, long]:
            if dst > 0:
                return type(self)(self.seq * dst, self.basis)
            elif dst == 0:
                return type(self)([], self.basis)
            elif dst < 0:
                seq = list(self.seq)
                seq.reverse()
                return type(self)(
                    [self.basis.invers(i) for i in seq] * (-dst), self.basis)
        else:
            return dst ** -1 * self * dst

    def __str__(self):
        er = self.exp_resp()
        r = []
        for w, c in zip(*er):
            if c == 1:
                wc = w
            else:
                wc = w + '**' + str(c)
            r.append(wc)
        return '*'.join(r)

    def __repr__(self):
        return "'" + self.__str__() + "'"

class grp_element(object):
    pass

class fpgrp_element(grp_element):

    def __init__(self, grp, word):
        self.group = grp
        

class group(object):
    pass

class fp_group(group):

    def __init__(self, gens, rels):
        pass

class valid_fp_group(fp_group):

    def __init__(self):
        pass


