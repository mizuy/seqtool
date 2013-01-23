from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from collections import defaultdict
from math import ceil, log, log10

from .nucleotide import bisulfite, base_color, cpg_sites, seq_cpg_analysis, seq_cpg_obs_per_exp
from .pcr import PCR, Primer
from . import xmlwriter
from .svg import *

__all__ = ['SeqviewTrack']

def get_feature_name(feature):
    t = feature.type
    try:
        if t in ['CDS','mRNA']:
            name = "%s: %s"%(t,feature.qualifiers['product'][0])
        elif t=='gene':
            name = 'gene: %s'%feature.qualifiers['gene'][0]
        elif t=='source':
            organ = feature.qualifiers['organism'][0]
            ch = feature.qualifiers['chromosome'][0]
            name = '%s chromosome %s'%(organ,ch)
        else:
            name = t
    except KeyError:
        name = feature.type
    return name

class SeqviewTrack(SvgMatrix):
    def __init__(self):
        super(SeqviewTrack,self).__init__()

    def add_named(self, name, track):
        self.add_row([SvgText(name, 0,0), track], ['right',None])

    def add(self, track):
        self.add_row([SvgBase(), track])

    def add_padding(self, height):
        self.add(SvgItemsFixedHeight(height))

    def add_hline(self, length, height=5):
        t = SvgItemsFixedHeight(height)
        t.add(SvgHline(0, length, height/2))
        self.add(t)

    def add_measure_track(self, length, start=0, step=500, substep=100, subsubstep=10):
        length = length
        start = start
        step = step
        substep = substep
        subsubstep = subsubstep

        numbers = SvgItemsFixedHeight(20)
        ss = step * int(ceil(1.0*start/step))
        for i in range(ss, start + length, step):
            x = i - start
            numbers.add(SvgText(str(i), x=x, y=10, color='black', anchor='middle'))

        bars = SvgItemsFixedHeight(20)
        bars.add(SvgHline(0, length, 20, color='black'))
        for i in range(ss, start + length, step):
            x = i - start
            bars.add(SvgVline(x, 10, 20, color='black'))
        # substep
        ss = substep * int(ceil(1.0*start/substep))
        for i in range(ss, start + length, substep):
            x = i - start
            bars.add(SvgVline(x, 13, 20, color='black'))

        ss = subsubstep * int(ceil(1.0*start/subsubstep))
        for i in range(ss, start + length, subsubstep):
            x = i - start
            bars.add(SvgVline(x, 18, 20, color='black'))

        self.add(numbers)
        self.add(bars)

    def add_sequence_track(self, seq, features, start=0):
        length = len(seq)

        msize = 500*(max(1,int(length/10)/500))

        self.add_measure_track(length, start, msize)
        self.add_padding(10)
        for f in features:
            self.add_feature_track(f)
        self.add_cpgisland_track(seq)

    def add_transcript_track(self, name, seq, feature):
        length = len(seq)

        msize = 500*(max(1,int(length/10)/500))

        self.add_measure_track(length, 0, msize)

        et = SvgItemsFixedHeight(10)

        def loc_len(loc):
            return loc.nofuzzy_end - loc.nofuzzy_start
        exons = [loc_len(sf.location) for sf in feature.sub_features]

        count = 0
        for e in exons:
            et.add(SvgVline(count, 0, 10, color='black'))
            count += e
        et.add(SvgVline(count, 0, 10, color='black'))
        et.add(SvgHline(0, count, 5, color='black'))

        self.add_named(name, et)

    def add_feature_track(self, feature):
        name = get_feature_name(feature)

        height = 14
        t = SvgItemsFixedHeight(height)
        loc = feature.location

        if feature.sub_features:
            t.add(SvgHline( loc.nofuzzy_start, loc.nofuzzy_end, t.height/2, color='black'))
            for sf in feature.sub_features:
                t.add(SvgHbar( sf.location.nofuzzy_start, sf.location.nofuzzy_end, t.height/2, 3., color='blue'))
        else:
            t.add(SvgHbar( loc.nofuzzy_start, loc.nofuzzy_end, t.height/2, 3., color='green'))


        self.add_named(name, t)

    def add_atgc_color_track(self, seq):
        t = SvgItemsFixedHeight(10)
        for i,c in enumerate(seq):
            color = base_color(c)
            t.add(SvgVline(i, 0, t.height, color=color))

        self.add_named('ATGC', t)


    def add_dbtss_track(self, rt, maxtag, width):
        t = SvgItemsFixedHeight(50)

        maxtag = maxtag
        name = rt.name

        w = width
        h = t.height
        vh = 1.*min(500,maxtag)
        t.add(SvgHline( 0, w, h, color='gray'))
        for x,v in rt.items():
            v = h*v/vh
            if v <= h:
                t.add(SvgVline(x, h-v, h, color='red'))
            else:
                st = 1.*v/h
                t.add(SvgVline(x, 0, h, color='blue', stroke=st))

        self.add_named(rt.name, t)
            
    def add_cpgbar_track(self, seq):
        height = 20
        t = SvgItemsFixedHeight(height)

        name = 'CpG site'
        length = len(seq)
        cpg = cpg_sites(seq)

        t.add(SvgHline(0, length, height/2, color='black'))
        for i in cpg:
            t.add(SvgVline(i, 0, height, color='black'))

        self.add_named('CpG', t)

    def add_cpgisland_track(self, seq, window=200):
        height = 20
        length = len(seq)

        gc_per, obs, cpgi = seq_cpg_analysis(seq, window)

        # GC Percent
        # self.add_named('GC %', SvgGraphline(height, gc_per, bars=[.55], stroke='black', fill='none'))

        # Obs/Exp
        # obs2 = [v/2. for v in obs]
        # self.add_named('Obs/Exp',SvgGraphline(height, obs2, bars=[.65/2], stroke='black', fill='none'))

        self.add_cpgislandbar_track(seq, cpgi)
        self.add_cpgbar_track(seq)
        self.add_padding(5)

    def add_cpgislandbar_track(self, seq, cpg_islands):
        height = 20
        t = SvgItemsFixedHeight(height)

        length = len(seq)
        cpg_islands = cpg_islands

        for start,stop in cpg_islands:
            p,q = max(0,start), min(length,stop)
            t.add(SvgHbar(p, q, height/2, 4, color='black'))
            
            if start < 0:
                t.add(SvgHbar(0, (q-p)/20., height/2, 4, color='red'))
            if length < stop:
                t.add(SvgHbar(length-(q-p)/20., length, height/2, 4, color='red'))
        self.add_named('CpG Island', t)

    def add_pcrs_track(self, name, pcrs):
        def d():
            for pcr in pcrs:
                for p in pcr.products:
                    yield pcr.name, p
        self.add_named(name, PcrsTrack(d()))

    def add_bs_pcrs_track(self, name, pcrs):
        def d():
            for pcr in pcrs:
                for met, unmet, genome in pcr.bisulfite_products:
                    assert(met or unmet or genome)
                    p = None
                    r = ''
                    if met:
                        r += "M"
                        p = met
                    if unmet:
                        r += "U"
                        p = unmet
                    if genome:
                        r += "G"
                        p = genome
                    yield pcr.name+' [%s]'%r, p
        self.add_named(name, PcrsTrack(d()))

    def add_primers_track(self, primers, template):
        self.add_named("Primers", PrimersTrack(primers, template))

    def add_bsa_track(self, name, bsp_map, start, end):
        self.add_named(name, BsaTrack(bsp_map, start, end))

### PCRs

class PcrsTrack(SvgItemsVFree):
    def __init__(self, products):
        super(PcrsTrack, self).__init__()

        self.products = []
        for name, p in products:
            m = (p.start+p.end)/2.
            strlen = len(name)*6
            s = min(p.start, m-strlen)
            e = max(p.end, m+strlen)

            named = SvgItemsVStack()

            bar = SvgItemsFixedHeight(8)
            bar.add(SvgHbar(p.start, p.start_i, 4, 4, color='red'))
            bar.add(SvgHbar(p.start_i, p.end_i, 4,4, color='gray'))
            bar.add(SvgHbar(p.end_i, p.end, 4, 4, color='blue'))
            named.add(bar)

            m = (p.start+p.end)/2.
            named.add(SvgText(name, x=m, y=0, color='black', anchor='middle'))

            self.add(named)

### Primers

class PrimersTrack(SvgItemsVFree):
    def __init__(self, primers, template):
        super(PrimersTrack, self).__init__()
        self.primers = []
        for p in primers:
            f,r = p.search(template)
            for ff in f:
                start,end = ff, ff+len(p)-1
                self.primers.append( (p.name, start, end, True) )
            for rr in r:
                start,end = rr, rr+len(p)-1
                self.primers.append( (p.name, start, end, False) )

        for name, start, end, forward in self.primers:
            named = SvgItemsVStack()
            named.add(SvgHbar(start, end, 1.5, 3, color = '#f00' if forward else '#00f'))
            named.add(SvgText(name, x=(start+end)/2., y=0, color='black', anchor='middle'))
            self.add(named)

### Bisulfite Sequencing Analysis

class BsaTrack(SvgItemsFixedHeight):
    def __init__(self, bsp_map, start, end):
        super(BsaTrack, self).__init__(20)
        self.bsp_map = bsp_map
        self.start = start
        self.end = end

        color = {'M':'#F00','U':'#00F','P':'#AA0'}
        
        self.add(SvgHline(self.start, self.end, 10, color='black'))
        for index, result in self.bsp_map:
            if result != '?':
                self.add(SvgVline(index, 2, 18, color=color[result]))
    
