
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
        self._reduce()

    def _reduce(self):
        i = 0
        while i < len(self.seq):
            w = self.seq[i]
            if not w in self.basis.elms:
                raise ValueError("invalid word seq {0}".format(w))
            if i + 1 < len(self.seq):
                wn = self.seq[i + 1]
                if self.basis.invers(w) == wn:
                    self.seq.pop(i)
                    self.seq.pop(i)
                    if i > 0:
                        i -= 2
            i += 1

    def __len__(self):
        return len(self.seq)

    @property
    def underlying(self):
        return self

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

    def __getitem__(self, slc):
        return type(self)(self.seq[slc], self.basis)

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

    @property
    def underlying(self):
        return self.word
        
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
@roprop('stopped')
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
            [rel.num_resp() for rel in rels if len(rel) > 0],
            [sub.num_resp() for sub in subs if len(sub) > 0])
        try:
            coset_tbl.run()
        except coset_tbl.noresult as ex:
            self._stopped = False
        else:
            self._stopped = True
        self._tbl = []
        self._finished = True
        for i in xrange(coset_tbl.gentbl.statid):
            tbl, tbli = coset_tbl.gentbl.states[i].tbl
            tra = {}
            for g in xrange(len(tbl)):
                #gen = str(self.basis.gen_num(g + 1))
                gen = self.basis.gen_num(g + 1).seq[0]
                if tbl[g] == None:
                    if self.finished:
                        self._finished = False
                    tra[gen] = None
                else:
                    tra[gen] = tbl[g].id
            for g in xrange(len(tbli)):
                #gen = str(self.basis.gen_num(- g - 1))
                gen = self.basis.gen_num(- g - 1).seq[0]
                if tbli[g] == None:
                    if self.finished:
                        self._finished = False
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

    #TODO breadth-first order
    def trav_word(self, word = None, sta = None, result = None,
                  loop_hndl = None, trac_inv = False):
        if word == None:
            word = self.basis.gen_num(0)
        if sta == None:
            sta = self.state(word)
        if result == None:
            result = [None] * len(self.tbl)
            result[sta] = word
        gens = self.basis.gens(trac_inv)
        trans = self.tbl[sta]
        for w in gens:
            w = w.seq[0]
            nsta = trans[w]
            if nsta == None:
                continue
            nwd = word * word.basis.gen(w)
            #print 'touch1', nwd,
            if result[nsta] == None:
                #print 'node'
                result[nsta] = nwd
                self.trav_word(nwd, nsta, result, loop_hndl, trac_inv)
            elif not loop_hndl == None and not result[nsta] == nwd:
                #print 'loop'
                loop_hndl(result[nsta], nwd, nsta, sta, w)
            else:
                #print 'back'
                pass
        return result

    def trav_loop(self, word = None, sta = None):
        #print "trav loop"
        result = []
        loop_tbl = [{} for _ in xrange(len(self.tbl))]
        def _add_loop(wds, wdd, nsta, sta, w):
            wdr = wdd * wds ** -1
            iw = self.basis.invers(w)
            assert sta == self.tbl[nsta][iw]
            assert not w in loop_tbl[sta]
            loop_tbl[sta][w] = (len(result), len(wdd), 1)
            assert not iw in loop_tbl[nsta]
            loop_tbl[nsta][iw] = (len(result), len(wdr) - len(wdd) + 1, -1)
            result.append(wdr)
        self.trav_word(word, sta, loop_hndl = _add_loop, trac_inv = False)
        return result, loop_tbl

    def replace_loop_words(self, words):
        if not len(words) > 0:
            return
        key_sets = [set() for _ in words]
        def _add_key(word, key_set):
            sta = 0
            for w in word.seq:
                if w in self.lp_tbl[1][sta]:
                    key_set.add(self.lp_tbl[1][sta][w])
                sta = self.tbl[sta][w]
                if sta == None:
                    raise ValueError("invalid word")
                elif sta == 0:
                    raise ValueError("word more than 1 loop")
        for wd, ks in zip(words, key_sets):
            _add_key(wd, ks)
        all_kidx = set(xrange(len(key_sets)))
        edge_valid_kidx = set()
        uni_sets = []
        for i in xrange(len(key_sets)):
            uset = set()
            for j in xrange(len(key_sets)):
                if i == j:
                    continue
                uset = uset.union(key_sets[j])
            uni_set = key_sets[i].difference(uset)
            if uni_set:
                edge_valid_kidx.add(i)
            uni_sets.append(uni_set)
        if not edge_valid_kidx:
            raise RuntimeError("bone shouldn't be closed")
        key_choosen = [None] * len(key_sets)
        while not edge_valid_kidx == all_kidx:
            kidx = all_kidx.difference(edge_valid_kidx).pop()
            for vkidx in edge_valid_kidx:
                iset = key_sets[kidx].intersection(key_sets[vkidx])
                if iset:
                    key_choosen[kidx] = iset.pop()
                    edge_valid_kidx.add(kidx)
                    break
            else:
                raise RuntimeError("unsolved loop relation")
        for i in xrange(len(key_sets)):
            if key_choosen[i] == None:
                key_choosen[i] = uni_sets[i].pop()
        def _replace_widx(word, widx):
            sta = 0
            for i in xrange(len(word)):
                w = word.seq[i]
                if (w in self.lp_tbl[1][sta] and
                    self.lp_tbl[1][sta][w] == widx):
                    self.lp_tbl[1][sta][w] = (widx[0], i+1, 1)
                    self.lp_tbl[0][widx[0]] = word
                    iw = self.basis.invers(w)
                    nsta = self.tbl[sta][w]
                    self.lp_tbl[1][nsta][iw] = (widx[0], len(word) - i, -1)
                    break
                sta = self.tbl[sta][w]
                if sta == None:
                    raise ValueError("invalid word")
        for i in xrange(len(key_sets)):
            _replace_widx(words[i], key_choosen[i])

    @property
    def lp_tbl(self):
        if not hasattr(self, '_lp_tbl'):
            self._lp_tbl = self.trav_loop(sta = 0)
            self.replace_loop_words(self.subs)
        return self._lp_tbl

    def recalc_subs(self, subs, rels = None):
        self.replace_loop_words(subs)
        self._subs = self.lp_tbl[0]
        if not rels == None:
            self._rels = rels

    def loop_resolve(self, word, hndl = None):
        sta = 0
        for i in xrange(len(word)):
            w = word.seq[i]
            if w in self.lp_tbl[1][sta]:
                loop_widx = self.lp_tbl[1][sta][w]
                loop_word = self.lp_tbl[0][loop_widx[0]] ** loop_widx[2]
                assert loop_widx[1] > 0
                comb_word1 = word[:i] * loop_word[:loop_widx[1]-1] ** -1
                comb_word2 = loop_word[loop_widx[1]:] ** -1 * word[i+1:]
                if hndl == None:
                    rslt = loop_word
                else:
                    rslt = hndl(loop_word, loop_widx)
                return (self.loop_resolve(comb_word1, hndl) +
                        [rslt] +
                        self.loop_resolve(comb_word2, hndl))
            sta = self.tbl[sta][w]
            if sta == None:
                raise ValueError("invalid word")
        if len(word) > 0:
            if hndl == None:
                return [word]
            else:
                return [hndl(word, None)]
        else:
            return []

    def genwds_maptbl(self, genwds):
        if not len(genwds) >= len(self.lp_tbl[0]):
            raise ValueError("too few genwds")
        maptbl = [None] * len(self.lp_tbl[0])
        for i in xrange(len(self.lp_tbl[0])):
            w = self.lp_tbl[0][i]
            try:
                gidx = genwds.index(w)
            except ValueError:
                raise ValueError("invalid genwds")
            maptbl[i] = gidx
        return maptbl

    def show(self):
        gens = self.basis.gens(True)
        def _s(s):
            if s == None:
                r = 'X'
            else:
                r = 'S' + str(s)
            return r + ' ' * (4 - len(r))
        for i in xrange(len(self.tbl)):
            print _s(i) + '->',
            for w in gens:
                w = w.seq[0]
                print w + ': ' +  _s(self.tbl[i][w]),
            print

    def is_normal(self, genwds):
        if not genwds:
            genwds = self.basis.gens()
        for h in self.subs:
            for gen in genwds:
                sta = self.state(gen)
                if sta == None or not sta == self.state(h, sta):
                    return False
        return True

    def is_trivial(self):
        return self.finished and len(self.tbl) == 1

    def is_open(self):
        for gen in self.basis.gens():
            if not self.state(gen) == None:
                return False
        return True

    def __contains__(self, dst):
        if not (isinstance(dst, grp_coset) and self.basis == dst.basis):
            return False
        if len(dst.tbl) < 2:
            if dst.finished:
                return True
        for sub in self.subs:
            if not dst.state(sub, 0) == 0:
                return False
        for rcs in xrange(len(dst.tbl)):
            for rel in self.rels:
                if not dst.state(rel, rcs) == rcs:
                    return False
        return True

    @lazyeq
    def __eq__(self, dst):
        return isinstance(dst, grp_coset) and self in dst and dst in self

    def __nonzero__(self):
        return True

@neq
class base_fp_group(object):

    @property
    def genwds(self):
        return []

    @property
    def rels(self):
        return []

    def has_subgroup(self, dst):
        return self == dst or (
            isinstance(dst, subgroup_of_fpgrp) and self == dst.fpgroup)

    def subgroup(self, gens):
        return subgroup_of_fpgrp(self, gens)

    def normalclosure(self, gens):
        return normalclosure_of_subgrp(self, gens)

    def quogroup(self, ker):
        raise TypeError('quotient invalid')

    def resolve_elem(self, elem):
        return elem.underlying

    def rebuild_elem(self, re_word):
        return re_word.mapped(self.gens)

    @iseq
    def __eq__(self, dst):
        if not isinstance(dst, base_fp_group):
            return False
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
        if self.trivial:
            return 1
        else:
            raise OverflowError('infinity length')

    @property
    def finite(self):
        return self.trivial

    @lazyprop
    def trivial(self):
        return not len(self.basis) > 0

    @property
    def gens(self):
        return self.basis.gens()

    @lazyprop
    def one(self):
        return self.basis.gen_num(0)

    @lazyprop
    def elems(self):
        if self.trivial:
            return [self.one]
        else:
            raise OverflowError('infinity group')

    @property
    def trans(self):
        return None

    @property
    def filt(self):
        return grp_coset(self.basis, self.gens)

    def has_element(self, dst):
        return dst in self.basis

    def quogroup(self, ker):
        if type(ker) in [tuple, list]:
            return fp_group(self, ker)
        elif self.has_subgroup(ker) and ker.is_normal(self):
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

    @property
    def finite(self):
        return self.trans.finished

    @lazyprop
    def trivial(self):
        return self.trans.is_trivial()

    @lazyprop
    def elems(self):
        if not self.trans.finished:
            raise TypeError('unfinished fp-group')
        return [fpgrp_element(
            self, w) for w in self.trans.trav_word(self.one.word)]

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

    @lazyprop
    def filt(self):
        return grp_coset(self.basis, self.frgroup.gens)

    def has_element(self, dst):
        return isinstance(dst, fpgrp_element) and self == dst.group

    def quogroup(self, ker):
        if self.has_subgroup(ker) and ker.is_normal(self):
            return fp_group(self.frgroup, self.rels + ker.genwds)
        else:
            raise TypeError('quotient invalid')

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
        if not fpgrp.filt.is_trivial():
            raise TypeError("fp-group shouldn't be subgroup resp")
        self._fpgroup = fpgrp
        self._gens = []
        for gen in gens:
            if not gen in fpgrp:
                raise ValueError("generator {0} is invalid".format(gen))
            if not gen in self.gens:
                self.gens.append(gen)
        self.gens.sort(key = lambda w: w.seq)
        self._genwds = [g.underlying for g in self.gens]

    def __len__(self):
        if self.trivial:
            return 1
        else:
            wlen = len(self.fpgroup)
            index = self.index
            assert wlen % index == 0
            return wlen / index

    @property
    def finite(self):
        return self.trivial or self.fpgroup.finite

    @lazyprop
    def trivial(self):
        if self.filt.is_open():
            return True
        elif self.fpgroup.finite and self.filt.finished:
            return len(self.fpgroup) == self.index
        else:
            return False

    @property
    def index(self):
        return len(self.filt)

    @lazyprop
    def elems(self):
        if self.trivial:
            return [self.one]
        else:
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

    @property
    def trans(self):
        return self.fpgroup.trans

    @lazyprop
    def filt(self):
        return grp_coset(self.basis, self.rels, self.genwds)

    def has_element(self, dst):
        return dst in self.fpgroup and (
            self.filt.state(dst) == 0)

    def has_subgroup(self, dst):
        if not (isinstance(
            dst, subgroup_of_fpgrp) and self.fpgroup == dst.fpgroup):
            return False
        return self.filt in dst.filt

    def is_normal(self, dst):
        if not self.basis == dst.basis:
            return False
        return self.filt.is_normal(dst.genwds)

    def subgroup(self, gens):
        for gen in gens:
            if not gen in self:
                raise ValueError("generator {0} is invalid".format(gen))
        return self.fpgroup.subgroup(gens)

    def normalclosure(self, gens):
        for gen in gens:
            if not gen in self:
                raise ValueError("generator {0} is invalid".format(gen))
        return self.fpgroup.normalclosure(gens)

    def quogroup(self, ker):
        if self.has_subgroup(ker) and ker.is_normal(self):
            fpgrp = self.fpgroup.quogroup(ker)
            gens = [fpgrp_element(fpgrp, w) for w in self.genwds]
            return fpgrp.subgroup(gens)
        else:
            raise TypeError('quotient invalid')

    @lazyprop
    def _genwds_maptbl(self):
        return self.filt.genwds_maptbl(self.genwds)

    @lazyprop
    def _genwds_basis(self):
        return grp_word_basis(
            ['gen' + str(i) for i in xrange(len(self.genwds))])

    def resolve_elem(self, elem):
        if not elem in self:
            raise ValueError("{0} is not elem of group".format(elem))
        def _hndl(word, widx):
            n = (self._genwds_maptbl[widx[0]] + 1) * widx[2]
            return self._genwds_basis.gen_num(n)
        r = self.filt.loop_resolve(elem.underlying, _hndl)
        return reduce(lambda a, b: a * b, r)

class normalclosure_of_subgrp(subgroup_of_fpgrp):

    def __init__(self, fpgrp, gens):
        self._need_calc = False
        super(normalclosure_of_subgrp, self).__init__(fpgrp, gens)
        self._need_calc = True

    def _recalc_gens(self):
        self.filt.recalc_subs(self._genwds, rels = self.rels)
        genwds = list(self.filt.subs)
        genwds.sort(key = lambda w: w.seq)
        self._genwds = genwds
        self._gens = [w.mapped(self.fpgroup.gens) for w in genwds]
        self._need_calc = False

    @property
    def gens(self):
        if self._need_calc:
            self._recalc_gens()
        return self._gens

    @property
    def genwds(self):
        if self._need_calc:
            self._recalc_gens()
        return self._genwds

    @lazyprop
    def filt(self):
        #print 'prev filt'
        return grp_coset(self.basis, self.rels + self._genwds)

