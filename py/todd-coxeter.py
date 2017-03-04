
from util import *

@roprop('gens')
@roprop('states')
class gen_table(object):

    def __init__(self, gens_num):
        self._gens = gens_num
        self._statid = 0
        self._states = {}

    def __len__(self):
        return len(self.tbl)

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

class rel_table(object):

    def __init__(self):
        pass

class coset_table(object):

    def __init__(self, gens_num, rels, subs):
        self.gentbl = gen_table(gens_num)
        #self.reltbl =
        #self.subtbl = 
        
    
