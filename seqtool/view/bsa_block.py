

from ..util.namedlist import DefaultNamedList
from ..util.parser import TreekvParser
from collections import defaultdict
from . import block

__all__ = ['BsaBlock']

class BsaResult(object):
    def __init__(self, pcr, result):
        self.pcr = pcr
        self.result = result

        p = pcr.bs_products(True, True)
        if not len(p)==1:
            raise ValueError('number of pcr products of %s must be 1 but %s'%(pcr.name,len(p)))
        self.product = p[0]
        if not all(i in 'MUP?' for i in result):
            raise ValueError('bsa result must contain only M,U,P or ?')

        self.cpg_sites = self.product.cpg_sites()

        self.num_cpg = self.product.num_cpg()
        if len(result)!=self.num_cpg:
            raise ValueError('%s has %s detectable CpG sites, but result gives %s'%(pcr.name,self.num_cpg,len(result)))

        self.bsa_map = [(n,result[i]) for i,n in enumerate(self.cpg_sites)]

class BsaCombined(object):
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
        for index, result in list(results.items()):
            if all(r in 'M' for r in result):
                c = 'M'
            elif all(r in 'U' for r in result):
                c = 'U'
            else:
                c = 'P'
            bsa_map.append( (index, c) )
            
        return bsa_map, start_i, end_i



class BsaBlock(block.BaseBlock):
    def __init__(self, bs_pcrs):
        self.bsas = DefaultNamedList(BsaCombined)
        self.celllines = None
        self.bs_pcrs = bs_pcrs

    def set_celllines(self, celllines):
        self.celllines = celllines

    def add(self, cellline, pcr_name, result):
        assert isinstance(cellline,str)
        assert isinstance(pcr_name,str)
        assert isinstance(result,str)

        pcr = self.bs_pcrs.pcrs.get(pcr_name)
        if not pcr:
            raise ValueError('no such pcr: %s'%pcr_name)

        self.bsas[cellline].add_bsa_result(pcr, result)


    def read(self, filename):
        self.readfp(open(filename,'rU'),filename)

    def readfp(self, fileobj, filename=None):
        parser = TreekvParser()
        tree = parser.readfp(fileobj, filename)

        for kv in list(tree.items()):
            n = [n.strip() for n in kv.key.split(',')]
            if not len(n)>=2:
                print('each bsa must have at least 2 key; cell line name and pcr name')
                continue
            pcrname = n[0].strip()
            cellline = n[1].strip().upper()
            #annotations = n[2:]
            if not pcrname or not cellline:
                print('empty pcr or cellline name: %s, %s'%(pcrname,cellline))
                continue

            self.add(cellline, pcrname, kv.value.strip().upper())

    def svg_genome(self, t):
        if self.celllines:
            for name in self.celllines:
                bsa_map, start, end = self.bsas[name].get_map()
                t.add_bsa_track(name, bsa_map, start, end)
        else:
            for name, bsa in list(self.bsas.items()):
                bsa_map, start, end = bsa.get_map()
                t.add_bsa_track(name, bsa_map, start, end)
