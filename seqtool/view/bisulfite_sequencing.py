from __future__ import absolute_import

from collections import defaultdict

__all__ = ['BisulfiteSequencingResult']

class BsaResult(object):
    def __init__(self, pcr, result):
        self.pcr = pcr
        self.result = result

        p = pcr.bs_met_products
        if not len(p)==1:
            raise ValueError('number of pcr products of %s must be 1 but %s'%(pcr_name,len(products)))
        self.product = p[0]
        if not all(i in 'MUP?' for i in result):
            raise ValueError('bsa result must contain only M,U,P or ?')

        self.num_cpg = self.product.detectable_cpg()
        if len(result)!=self.num_cpg:
            raise ValueError('%s has %s detectable CpG sites, but result gives %s'%(pcr_name,self.num_cpg,len(result)))

        self.cpg_sites = self.product.cpg_sites()

        self.bsa_map = [(n,result[i]) for i,n in enumerate(self.cpg_sites)]

class BisulfiteSequencingResult(object):
    def __init__(self, name):
        self.results = []
        self.name = name

    def add_bsa_result(self, pcr, result):
        self.results.append(BsaResult(pcr, result))

    def get_map(self):
        start_i = min(e.product.start_i for e in self.results)
        end_i = max(e.product.end_i for e in self.results)

        results = defaultdict(str)
        for e in self.results:
            for i,n in enumerate(e.cpg_sites):
                if e.result[i]!='?':
                    results[n] += e.result[i]

        bsa_map = []
        for index, result in results.items():
            if all(r in 'M' for r in result):
                c = 'M'
            elif all(r in 'U' for r in result):
                c = 'U'
            else:
                c = 'P'
            bsa_map.append( (index, c) )
            
        return bsa_map, start_i, end_i
