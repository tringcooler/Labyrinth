
from util import *

@roprop('gens')
@roprop('states')
@roprop('statid')
class gen_table(object):

    def __init__(self, gens_num):
        self._gens = gens_num
        self._statid = 0
        self._states = {}

    def __len__(self):
        return len(self.states)

    def alloc(self):
        state = vchain()
        state.tbl = ([None] * self.gens, [None] * self.gens)
        state.id = self._statid
        self._statid += 1
        self.states[state.id] = state
        return state

    def trans(self, sta, gen, chk = True):
        if gen == 0:
            raise ValueError('gen is based on 1')
        elif gen > 0:
            gen = gen - 1
            tra = 0
        elif gen < 0:
            gen = - gen - 1
            tra = 1
        rsta = sta.tbl[tra][gen]
        assert chk == False or rsta == None or rsta.tbl[tra^1][gen] == sta
        return rsta

    def record(self, sta, gen, rsta):
        if not (sta and rsta):
            return
        #print 'rec', sta.id, gen, rsta.id
        if gen == 0:
            raise ValueError('gen is based on 1')
        elif gen > 0:
            gen = gen - 1
            ssta = sta
            dsta = rsta
        elif gen < 0:
            gen = - gen - 1
            ssta = rsta
            dsta = sta
        odsta = ssta.tbl[0][gen]
        ossta = dsta.tbl[1][gen]
        if not (odsta == None or odsta == dsta):
            self.confluent(odsta, dsta)
        if not (ossta == None or ossta == ssta):
            self.confluent(ossta, ssta)
        ssta.tbl[0][gen] = dsta
        dsta.tbl[1][gen] = ssta

    def confluent(self, sta1, sta2):
        #print 'con', sta1.id, sta2.id
        record_req = []
        for gen_pair in xrange(self.gens):
            for gen in [gen_pair + 1, - gen_pair - 1]:
                nsta = self.trans(sta2, gen, False)
                #print 'ow', sta1.id, sta2.id, gen, nsta
                record_req.append((sta1, gen, nsta))
        del self.states[sta2.id]
        sta2.vlink(sta1)
        for r in record_req:
            self.record(*r)

    def show(self):
        def _s(s):
            if s:
                return 'S' + str(s.id)
            else:
                return 'X'
        for i in sorted(self.states.keys()):
            sta = self.states[i]
            print i, _s(sta), '->',
            for gen_pair in xrange(self.gens):
                for gen in [gen_pair + 1, - gen_pair - 1]:
                    print str(gen) + ':', _s(self.trans(sta, gen)),
            print

    def resort(self):
        stid = 0
        nstates = {}
        for i in sorted(self.states.keys()):
            sta = self.states[i]
            sta.id = stid
            nstates[stid] = sta
            stid += 1
        self._states = nstates
        self._statid = stid

@roprop('rel')
@roprop('gtbl')
class rel_chain(object):

    def __init__(self, gen_tbl, rel):
        assert len(rel) > 0
        self._rel = rel
        self._gtbl = gen_tbl

    def deduce(self, sta):
        cur = sta
        for gen in self.rel[:-1]:
            nxt = self.gtbl.trans(cur, gen)
            if not nxt:
                nxt = self.gtbl.alloc()
                self.gtbl.record(cur, gen, nxt)
            cur = nxt
        self.gtbl.record(cur, self.rel[-1], sta)

@roprop('gentbl')
@roprop('reltbl')
@roprop('subtbl')
class coset_table(object):

    noresult = type('noresult', (Exception,), {})

    def __init__(self, gens_num, rels, subs):
        self._gentbl = gen_table(gens_num)
        self._reltbl = [
            rel_chain(self.gentbl, rel) for rel in rels]
        self._subtbl = [
            rel_chain(self.gentbl, sub) for sub in subs]

    def run(self, max_idx = 50000):
        if len(self.gentbl):
            raise RuntimeError('already run')
        self.gentbl.alloc()
        sta_idx = 0
        while sta_idx < self.gentbl.statid:
            if sta_idx > max_idx:
                raise self.noresult(sta_idx)
            try:
                sta = self.gentbl.states[sta_idx]
            except KeyError:
                continue
            finally:
                sta_idx += 1
            if sta_idx == 1:
                for chain in self.subtbl:
                    chain.deduce(sta)
            for chain in self.reltbl:
                chain.deduce(sta)
        self.gentbl.resort()

