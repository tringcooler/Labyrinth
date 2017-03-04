
import re

from util import *

@neq
@roprop('elms')
class grp_word_basis(object):

    def __init__(self, names):
        elms = []
        for n in set(names):
            if type(n) in [int, long]:
                raise TypeError("name shouldn't be int")
            elms.append(n)
            elms.append(self._invers_pres(n))
        self._elms = elms

    @property
    def id_elm(self):
        return self._ident_pres()

    @property
    def gen_elms(self):
        return [self._elms[i] for i in xrange(0, len(self._elms), 2)]

    @property
    def inv_elms(self):
        assert len(self._elms) % 2 == 0
        return [self._elms[i + 1] for i in xrange(0, len(self._elms), 2)]

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

    def __add__(self, dst):
        if self == dst:
            return self
        else:
            r = type(self)(self.gen_elms + dst.gen_elms)
            return r

    def __len__(self):
        assert len(self._elms) % 2 == 0
        return len(self._elms) / 2

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

@neq
@roprop('basis')
@roprop('seq')
class grp_word(object):

    def __init__(self, wl, basis):
        self._basis = basis
        self._seq = wl

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
        if len(gens) < len(self.basis):
            raise ValueError(
                "gens shoud have at least {0}".format(len(self.basis)))
        trans_dict = {
            syl: gen for gen, syl in zip(gens, self.basis.gen_elms)}
        rslt = None
        for syl, cnt in zip(*self.exp_resp()):
            word = trans_dict[syl] ** cnt
            if rslt == None:
                rslt = word
            else:
                rslt *= word
        return rslt

    def __eq__(self, dst):
        return self.basis == dst.basis and self.seq == dst.seq

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
        self.word = word

    def __mul__(self, dst):
        pass

    def __pow__(self, dst):
        pass

class base_group(object):
    pass

@neq
@roprop('frgroup')
@roprop('gens')
@roprop('rels')
class fp_group(base_group):

    def __new__(cls, frgrp, rels):
        if not rels:
            return frgrp
        else:
            return super(fp_group, cls).__new__(cls)

    def __init__(self, frgrp, rels):
        self._frgroup = frgrp
        self._gens = [fpgrp_element(self, g) for g in frgrp.gens]
        for rel in rels:
            if not rel in frgrp:
                raise ValueError("relator {0} is invalid".format(rel))
        self._rels = list(set(rels))

    def __eq__(self, dst):
        return (
            self.frgroup == dst.frgroup ) and (
            self.gens == dst.gens ) and (
            self.rels == dst.rels )

@neq
@roprop('gens')
class free_group(base_group):

    def __init__(self, basis, gens = None):
        if not isinstance(basis, grp_word_basis):
            basis = grp_word_basis(basis)
        if not gens:
            self._gens = basis.gens()
        else:
            self._gens = []
            for gen in gens:
                if not isinstance(gen, grp_word):
                    gen = basis.gen(gen)
                self._gens.append(gen)
        self._gens = list(set(self._gens))

    def __eq__(self, dst):
        return self.gens == dst.gens

    def __contains__(self, dst):
        pass

class valid_fp_group(fp_group):

    def __init__(self):
        pass


