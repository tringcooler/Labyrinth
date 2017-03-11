
import re

from util import *
from todd_coxeter import coset_table

@neq
@roprop('elms')
class grp_word_basis(object):

    def __init__(self, names):
        elms = []
        for n in sorted(names):
            if type(n) in [int, long]:
                raise TypeError("name shouldn't be int")
            if not n in elms:
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
        idx = self.elms.index(n)
        pos = (idx // 2) + 1
        if idx % 2:
            pos = - pos
        return pos

    def invers(self, n):
        if not n in self.elms:
            raise ValueError('element {0} invalid'.format(n))
        if self.is_invers(n):
            return self._invers_invers_pres(n)
        else:
            return self._invers_pres(n)

    @iseq
    def __eq__(self, dst):
        return isinstance(dst, grp_word_basis) and self.elms == dst.elms

    def __contains__(self, dst):
        return isinstance(dst, grp_word) and self == dst.basis

    def __add__(self, dst):
        if self == dst:
            return self
        else:
            r = type(self)(self.gen_elms + dst.gen_elms)
            return r

    def __len__(self):
        assert len(self._elms) % 2 == 0
        return len(self._elms) / 2

    def gen_num(self, n):
        if n == 0:
            return grp_word([], self)
        else:
            pos = (abs(n) - 1) * 2
            if n < 0:
                pos += 1
            if pos >= len(self.elms):
                raise ValueError('element num {0} invalid'.format(n))
            return grp_word([self.elms[pos]], self)

    def gen(self, n):
        if n == self._ident_pres():
            return grp_word([], self)
        else:
            if not n in self.elms:
                raise ValueError('element {0} invalid'.format(n))
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

    def __len__(self):
        return len(self.seq)

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

    def num_resp(self):
        return [self.basis.position(w) for w in self.seq]

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

    @iseq
    def __eq__(self, dst):
        return isinstance(dst, grp_word) and (
            self.basis == dst.basis and self.seq == dst.seq)

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
        if not self.seq:
            return self.basis.id_elm
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

@neq
@roprop('group')
@roprop('word')
class fpgrp_element(object):

    def __init__(self, grp, word):
        if not grp.basis == word.basis:
            raise ValueError("invalid basis of '{0}'".format(word))
        self._group = grp
        self._word = word

    @property
    def seq(self):
        return self.word.seq
        
    @lazyprop
    def state(self):
        try:
            return self.group.trans.state(self)
        except AttributeError:
            raise TypeError('only for fp-group and words')

    @iseq
    def __eq__(self, dst):
        return isinstance(dst, fpgrp_element) and (
            self.group == dst.group and self.state == dst.state)

    def __mul__(self, dst):
        if not self.group == dst.group:
            raise ValueError('from different groups')
        return type(self)(
            self.group, self.word * dst.word)

    def __pow__(self, dst):
        if type(dst) in [int, long]:
            return type(self)(
                self.group, self.word ** dst)
        elif self.group == dst.group:
            return type(self)(
                self.group, self.word ** dst.word)
        else:
            raise ValueError('from different groups')

    def __str__(self):
        return str(self.word)
    
    def __repr__(self):
        return repr(self.word)

@neq
@roprop('basis')
@roprop('rels')
@roprop('subs')
@roprop('tbl')
@roprop('finished')
class grp_coset(object):

    def __init__(self, basis, rels, subs = None):
        if subs == None:
            subs = []
        self._basis = basis
        self._rels = rels
        self._subs = subs
        self._calc(rels, subs)

    def __len__(self):
        if self.finished:
            return len(self.tbl)
        else:
            #return float('inf')
            raise OverflowError('infinity length')

    def _calc(self, rels, subs):
        #print 'calc coset'
        coset_tbl = coset_table(
            len(self.basis),
            [rel.num_resp() for rel in rels],
            [sub.num_resp() for sub in subs])
        try:
            coset_tbl.run()
        except coset_tbl.noresult as ex:
            self._finished = False
        else:
            self._finished = True
        self._tbl = []
        for i in xrange(coset_tbl.gentbl.statid):
            tbl, tbli = coset_tbl.gentbl.states[i].tbl
            tra = {}
            for g in xrange(len(tbl)):
                #gen = str(self.basis.gen_num(g + 1))
                gen = self.basis.gen_num(g + 1).seq[0]
                if tbl[g] == None:
                    tra[gen] = None
                else:
                    tra[gen] = tbl[g].id
            for g in xrange(len(tbli)):
                #gen = str(self.basis.gen_num(- g - 1))
                gen = self.basis.gen_num(- g - 1).seq[0]
                if tbli[g] == None:
                    tra[gen] = None
                else:
                    tra[gen] = tbli[g].id
            self.tbl.append(tra)

    def state(self, word, st = 0):
        #print 'parse word'
        sta = st
        for w in word.seq:
            sta = self.tbl[sta][w]
            if sta == None:
                return None
        return sta

    def trav_word(self, word, sta = None, result = None):
        if result == None:
            result = [None] * len(self)
        if sta == None:
            sta = self.state(word)
        trans = self.tbl[sta]
        for w in trans:
            nsta = trans[w]
            nwd = word * word.basis.gen(w)
            if result[nsta] == None:
                result[nsta] = nwd
                self.trav_word(nwd, nsta, result)
        return result

    def __contains__(self, dst):
        if not (isinstance(dst, grp_coset) and (
            self.basis == dst.basis and dst.finished)):
            return False
        if len(dst) < 2:
            return True
        for rel in self.rels:
            if not dst.state(rel, 1) == 1:
                return False
        for sub in self.subs:
            if not dst.state(sub, 0) == 0:
                return False
        return True

    @iseq
    def __eq__(self, dst):
        return isinstance(dst, grp_coset) and self in dst and dst in self

@neq
class base_fp_group(object):

    @property
    def rels(self):
        return []

    @property
    def trans(self):
        return None

    @property
    def filt(self):
        return None

    @property
    def elems(self):
        raise TypeError('infinity group')

    def subgroup(self, gens):
        return subgroup_of_fpgrp(self, gens)

    def quogroup(self, ker):
        raise TypeError('quotient invalid')

    @iseq
    def __eq__(self, dst):
        if not isinstance(dst, base_fp_group):
            return False
        elif not (self.trans or self.filt or dst.trans or dst.filt):
            return self.basis == dst.basis
        else:
            return self.trans == dst.trans and self.filt == dst.filt

    def __contains__(self, dst):
        return self.has_element(dst)

    # subgroup
    def __getitem__(self, gens):
        if type(gens) == slice:
            raise TypeError('slice unsupported')
        elif type(gens) in [tuple, list]:
            return self.subgroup(gens)
        else:
            return self.subgroup([gens])

    # quotient group
    def __div__(self, ker):
        return self.quogroup(ker)

@roprop('basis')
class free_group(base_fp_group):

    def __init__(self, basis):
        if not isinstance(basis, grp_word_basis):
            basis = grp_word_basis(basis)
        self._basis = basis

    def __len__(self):
        #return float('inf')
        raise OverflowError('infinity length')

    @property
    def gens(self):
        return self.basis.gens()

    @lazyprop
    def one(self):
        return self.basis.gen_num(0)

    def has_element(self, dst):
        return dst in self.basis

    def has_subgroup(self, dst):
        return isinstance(dst, subgroup_of_fpgrp) and self == dst.fpgroup

    def quogroup(self, ker):
        if type(ker) in [tuple, list]:
            return fp_group(self, ker)
        elif self.has_subgroup(ker):
            rels = ker.genwds
            return fp_group(self, rels)
        else:
            raise TypeError('quotient invalid')

@roprop('frgroup')
@roprop('rels')
class fp_group(base_fp_group):

    def __new__(cls, frgrp, rels):
        if not rels:
            return frgrp
        else:
            return super(fp_group, cls).__new__(cls)

    def __init__(self, frgrp, rels):
        self._frgroup = frgrp
        self._rels = []
        for rel in rels:
            if not rel in frgrp:
                raise ValueError("relator {0} is invalid".format(rel))
            if not rel in self.rels:
                self.rels.append(rel)
        self.rels.sort(key = lambda w: w.seq)

    def __len__(self):
        return len(self.trans)

    @lazyprop
    def elems(self):
        if not self.trans.finished:
            raise TypeError('unfinished fp-group')
        return [fpgrp_element(
            self, w) for w in self.trans.trav_word(self.one.word)]

    def has_element(self, dst):
        return isinstance(dst, fpgrp_element) and self == dst.group

    @property
    def basis(self):
        return self.frgroup.basis

    @lazyprop
    def gens(self):
        return [fpgrp_element(self, g) for g in self.frgroup.gens]

    @lazyprop
    def one(self):
        return fpgrp_element(self, self.frgroup.one)

    @lazyprop
    def trans(self):
        return grp_coset(self.basis, self.rels)

@roprop('fpgroup')
@roprop('gens')
@roprop('genwds')
class subgroup_of_fpgrp(base_fp_group):

    def __new__(cls, fpgrp, gens):
        if not gens:
            return fpgrp
        else:
            return super(subgroup_of_fpgrp, cls).__new__(cls)

    def __init__(self, fpgrp, gens):
        self._fpgroup = fpgrp
        self._gens = []
        for gen in gens:
            if not gen in fpgrp:
                raise ValueError("generator {0} is invalid".format(gen))
            if not gen in self.gens:
                self.gens.append(gen)
        self.gens.sort(key = lambda w: w.seq)
        if isinstance(fpgrp, fp_group):
            self._genwds = [g.word for g in self.gens]
        else:
            self._genwds = self.gens

    def __len__(self):
        wlen = len(self.fpgroup)
        index = self.index
        assert wlen % index == 0
        return wlen / index

    @property
    def index(self):
        return len(self.filt)

    @lazyprop
    def elems(self):
        return [e for e in self.fpgroup.elems if e in self]

    @property
    def basis(self):
        return self.fpgroup.basis

    @property
    def rels(self):
        return self.fpgroup.rels

    @property
    def one(self):
        return self.fpgroup.one

    @lazyprop
    def filt(self):
        return grp_coset(self.basis, self.rels, self.genwds)

    def has_element(self, dst):
        return dst in self.fpgroup and (
            self.filt.state(dst) == 0)

    def subgroup(self, gens):
        return type(self)(
            self.fpgroup, self.gens + list(gens))

class valid_fp_group(fp_group):

    def __init__(self):
        pass


