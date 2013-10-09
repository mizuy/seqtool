

from math import ceil

from ..nucleotide import base_color
from ..nucleotide.cpg import cpg_sites
from ..nucleotide.cpgisland import seq_cpg_analysis
from ..util import svg

__all__ = ['OutlineRenderer']

CSS = '''
rect.subfeature {
    fill: blue; stroke: blue;
}
rect.feature{
    fill: green; stroke: green;
}
.pcr_body{
    fill: gray; stroke: gray;
}
.pcr_fw{
    fill: red; stroke: red;
}
.pcr_rv{
    fill: blue; stroke: blue;
}
.bsa_methyl{
    fill: red; stroke: red;
}
.bsa_unmethyl{
    fill: blue; stroke: blue;
}
.bsa_partial{
    fill: #AA0; stroke: #AA0;
}
.dbtss{
    stroke: red;
}
.dbtss_excess{
    stroke: blue;
}
'''

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
        elif 'note' in feature.qualifiers.keys():
            name = '{}:{}'.format(t,feature.qualifiers['note'][0])
        else:
            name = t
    except KeyError:
        name = feature.type
    return name

class NamedTracks(svg.SvgMatrix):
    def __init__(self, scale=1):
        self.scale = scale
        self.gen = svg.SvgItemGenerator(scale, 1)
        super(NamedTracks,self).__init__()

    def add_named(self, name, track):
        self.add_row([self.gen.text(name, 0,0), self.gen.text(' ',0,0), track], ['right',None,None])

    def add(self, track):
        self.add_row([None, None, track])

    def add_padding(self, height):
        self.add(svg.SvgItemsFixedHeight(height))

    def add_hline(self, length, height=5):
        t = svg.SvgItemsFixedHeight(height)
        t.add(self.gen.hline(0, length, height/2))
        self.add(t)


class OutlineRenderer(NamedTracks):
    def __init__(self, scale):
        super(OutlineRenderer,self).__init__(scale)

    def svg_css(self):
        return super().svg_css() + CSS

    def add_measure_track(self, length, start):
        # 10, 100, 1000

        if 10 / self.scale > 30:
            self.add_measure_index_track(length, start, 100)
            self.add_measure_bar_track(length, start, 100, 10, 1)
        elif 100 / self.scale > 30:
            self.add_measure_index_track(length, start, 100)
            self.add_measure_bar_track(length, start, 100, 10)
        elif 500 / self.scale > 30:
            self.add_measure_index_track(length, start, 500)
            self.add_measure_bar_track(length, start, 500, 100)
        elif 1000 / self.scale > 30:
            self.add_measure_index_track(length, start, 1000)
            self.add_measure_bar_track(length, start, 1000, 100)
        elif 5000 / self.scale > 30:
            self.add_measure_index_track(length, start, 5000)
            self.add_measure_bar_track(length, start, 5000, 1000)


    def add_measure_index_track(self, length, start=0, step=100):
        ss = step * int(ceil(1.0*start/step))

        t = svg.SvgItemsFixedHeight(20)
        for i in range(ss, start + length, step):
            x = i - start
            t.add(self.gen.text(str(i), x=x, y=10, anchor='middle'))
        self.add(t)

    def add_measure_bar_track(self, length, start=0, step=100, substep=None, subsubstep=None):
        ss = step * int(ceil(1.0*start/step))

        t = svg.SvgItemsFixedHeight(20)
        
        t.add(self.gen.hline(0, length, 20))

        lines = []
        for i in range(ss, start + length, step):
            x = i - start
            lines.append((x,10,x,20))

        if substep:
            ss = substep * int(ceil(1.0*start/substep))
            for i in range(ss, start + length, substep):
                x = i - start
                lines.append((x,13,x,20))

        if subsubstep:
            ss = subsubstep * int(ceil(1.0*start/subsubstep))
            for i in range(ss, start + length, subsubstep):
                x = i - start
                lines.append((x,18,x,20))
        if lines:
            t.add(self.gen.lines(lines))

        self.add(t)

    def add_sequence_track(self, seq, features, start=0):
        length = len(seq)

        self.add_measure_track(length, start)
        self.add_padding(10)
        for f in features:
            self.add_feature_track(f)
        self.add_cpgisland_track(seq)
        #self.add_eachbase_color_track(seq)
        #self.add_atgc_color_track(seq)

    def add_transcript_track(self, name, seq, feature):
        length = len(seq)

        self.add_measure_track(length, 0)

        et = svg.SvgItemsFixedHeight(10)

        exons = [len(part) for part in feature.location.parts]

        count = 0
        for e in exons:
            et.add(self.gen.vline(count, 0, 10))
            count += e
        et.add(self.gen.vline(count, 0, 10))
        et.add(self.gen.hline(0, count, 5))

        self.add_named(name, et)

    def add_feature_track(self, feature):
        name = get_feature_name(feature)

        height = 14
        t = svg.SvgItemsFixedHeight(height)
        loc = feature.location

        if len(loc.parts) > 1:
            t.add(self.gen.hline( loc.nofuzzy_start, loc.nofuzzy_end, t.height/2))
            for part in loc.parts:
                t.add(self.gen.hbar( part.nofuzzy_start, part.nofuzzy_end, t.height/2, 3., klass='subfeature'))
        else:
            t.add(self.gen.hbar( loc.nofuzzy_start, loc.nofuzzy_end, t.height/2, 3., klass='feature'))


        self.add_named(name, t)

    def add_atgc_color_track(self, seq):
        t = svg.SvgItemsFixedHeight(10)
        a_lines = []
        t_lines = []
        g_lines = []
        c_lines = []
        atgc_lines = {'A':a_lines, 'T':t_lines, 'G':g_lines, 'C':c_lines}
        for i,c in enumerate(seq):
            atgc_lines[c].append((i,0, i,t.height))
        for base, lines in atgc_lines.items():
            color = base_color(base)
            if lines:
                t.add(self.gen.lines(lines, color=color))

        self.add_named('ATGC', t)

    def add_eachbase_color_track(self, seq):
        for base in 'ATGC':
            color = base_color(base)
            t = svg.SvgItemsFixedHeight(10)
            lines = []
            for i,c in enumerate(seq):
                if c==base:
                    lines.append((i,0, i,t.height))
            if lines:
                t.add(self.gen.lines(lines, color=color))
            self.add_named(base, t)


    def add_dbtss_track(self, rt, maxtag, seq):
        t = svg.SvgItemsFixedHeight(50)

        w = len(seq)
        h = t.height
        vh = 1.*min(500,maxtag)
        t.add(self.gen.hline( 0, w, h, color='gray'))

        lines_red = []
        lines_blue = []
        for x,v in list(rt.items()):
            v = h*v/vh
            if v <= h:
                #t.add(self.gen.vline(x, h-v, h, color='red'))
                lines_red.append((x,h-v,x,h))
            else:
                #st = 1.*v/h
                #t.add(self.gen.vline(x, 0, h, color='blue', stroke=st))
                lines_blue.append((x,0,x,h))
        if lines_red:
            t.add(self.gen.lines(lines_red, klass='dbtss'))
        if lines_blue:
            t.add(self.gen.lines(lines_blue, klass='dbtss_excess'))

        self.add_named(rt.name, t)
            
    def add_cpgbar_track(self, seq):
        height = 20
        t = svg.SvgItemsFixedHeight(height)

        length = len(seq)
        cpg = cpg_sites(seq)

        t.add(self.gen.hline(0, length, height/2))
        lines = []
        for i in cpg:
            lines.append((i,0,i,height))
        if lines:
            t.add(self.gen.lines(lines))

        self.add_named('CpG', t)

    def add_cpgisland_track(self, seq, window=200):
        gc_per, obs, cpgi = seq_cpg_analysis(seq, window)

        """
        # GC Percent
        self.add_named('GC %', svg.SvgGraphline(30, gc_per, bars=[.55], stroke='black', fill='none'))

        # Obs/Exp
        obs2 = [v/2. for v in obs]
        self.add_named('Obs/Exp',svg.SvgGraphline(30, obs2, bars=[.65/2], stroke='black', fill='none'))
        """

        self.add_cpgislandbar_track(seq, cpgi)
        self.add_cpgbar_track(seq)
        self.add_padding(5)

    def add_cpgislandbar_track(self, seq, cpg_islands):
        height = 20
        t = svg.SvgItemsFixedHeight(height)

        length = len(seq)
        cpg_islands = cpg_islands

        for start,stop in cpg_islands:
            p,q = max(0,start), min(length,stop)
            t.add(self.gen.hbar(p, q, height/2, 4, color='black'))
            
            if start < 0:
                t.add(self.gen.hbar(0, (q-p)/20., height/2, 4, color='red'))
            if length < stop:
                t.add(self.gen.hbar(length-(q-p)/20., length, height/2, 4, color='red'))
        self.add_named('CpG Island', t)

    def get_pcrs_track(self, products):
        t = svg.SvgItemsVFree()
        self.products = []
        for name, p in products:
            m = (p.start+p.end)/2.

            named = svg.SvgItemsVStack()

            bar = svg.SvgItemsFixedHeight(8)
            bar.add(self.gen.hbar(p.start, p.start_i, 4, 4, klass='pcr_fw'))
            bar.add(self.gen.hbar(p.start_i, p.end_i, 4,4, klass='pcr_body'))
            bar.add(self.gen.hbar(p.end_i, p.end, 4, 4, klass='pcr_rv'))
            named.add(bar)

            m = (p.start+p.end)/2.
            named.add(self.gen.text(name, x=m, y=0, anchor='middle'))

            t.add(named)
        return t

    def add_pcrs_track(self, name, pcrs):
        def d():
            for pcr in pcrs:
                for p in pcr.products:
                    yield pcr.name, p
        self.add_named(name, self.get_pcrs_track(d()))

    def add_bs_pcrs_track(self, name, pcrs):
        def d():
            for pcr in pcrs:
                for band in pcr.bands():
                    yield '{}[{}]'.format(pcr.name,band.match_str), band.get_product()
        self.add_named(name, self.get_pcrs_track(d()))

    def add_primers_track(self, primers, template):
        t = svg.SvgItemsVFree()
        primers = []
        for p in primers:
            f,r = p.search(template)
            for a in f:
                primers.append( (p.name, a.left, a.right, True) )
            for a in r:
                primers.append( (p.name, a.left, a.right, False) )

        for name, start, end, forward in primers:
            named = svg.SvgItemsVStack()
            named.add(self.gen.hbar(start, end, 1.5, 3, klass = 'pcr_fw' if forward else 'pcr_rv'))
            named.add(self.gen.text(name, x=(start+end)/2., y=0, anchor='middle'))
            t.add(named)
        self.add_named("Primers", t)

    def add_bsa_track(self, name, bsp_map, start, end):
        t = svg.SvgItemsFixedHeight(20)
        bsp_map = bsp_map
        start = start
        end = end

        klass = {'M':'bsa_methyl','U':'bsa_unmethyl','P':'bsa_partial'}
        
        t.add(self.gen.hline(start, end, 10))
        methyl = []
        unmethyl = []
        partial = []
        lines = {'M': methyl, 'U':unmethyl, 'P':partial}
        for index, result in bsp_map:
            if result != '?':
                lines[result].append((index, 2, index, 18))
        for char, ll in lines.items():
            if ll:
                t.add(self.gen.lines(ll, klass=klass[char]))
        self.add_named(name, t)

