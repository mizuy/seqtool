from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from collections import defaultdict
from math import ceil, log, log10

from StringIO import StringIO
from .nucleotide import bisulfite, base_color, cpg_sites, cpg_obs_per_exp
from .pcr import PCR, Primer
from . import xmlwriter


# calc free space
class FreeSpace(object):
    def __init__(self):
        self.lines = []

    def set(self, start, end):
        if start > end:
            raise ValueError('start=%s, end=%s, not start>end'%(start,end))

        for lineno, l in enumerate(self.lines):
            # if there is no common region, return this line
            #if all(max(start, s) > min(end, e) for s,e in l):
            if all(end<s or e<start for s,e in l):
                l.append((start,end))
                return lineno
        # new line
        self.lines.append([(start,end)])
        return len(self.lines)-1

    def num_lines(self):
        return len(self.lines)+1

class TrackBase(object):
    @property
    def name(self):
        return ''

    @property
    def name_length(self):
        return len(self.name)

class TrackGroup(TrackBase):
    def __init__(self):
        self.tracks = []
        self.ymargin = 4

    def __iter__(self):
        return iter(self.tracks)

    def add(self, item):
        self.tracks.append(item)
    
    @property
    def width(self):
        if len(self.tracks) == 0:
            return 0
        return max(i.width for i in self)
    @property
    def height(self):
        return sum(i.height + self.ymargin for i in self)
    @property
    def name_length(self):
        if len(self.tracks) == 0:
            return 0
        return max(i.name_length for i in self)

    def __str__(self):
        "TrackGroup(%s, %s)"%(self.name, ','.join(i for i in self))
        
    def draw(self, b, scale):
        y = 0
        for i in self:
            with b.g(transform='translate(0,%d)'%y):
                i.draw(b, scale)
            y += (i.height + self.ymargin)

    def svg(self, dst_width=None):
        padding_right = 20

        width = self.width
        height = self.height

        nl = (self.name_length+1)*6.0

        if not dst_width:
            dst_width = nl + width + padding_right
        dst_height = height
            
        scalex = (dst_width-nl-padding_right)/width
        view_w = (nl+padding_right)/scalex + width
        view_h = height

        trans_sx = dst_width/view_w
        trans_sy = dst_height/view_h
        
        buff = StringIO()
        b = xmlwriter.builder(xmlwriter.XmlWriter(buff))
        with b.svg(xmlns="http://www.w3.org/2000/svg", 
                   width=dst_width,
                   height=dst_height,
                   viewBox="%d 0 %d %d"%(-nl/scalex, view_w, view_h),
                   preserveAspectRatio='none'):
            self.draw(b, (1./trans_sx,1./trans_sy))

        return \
'''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
         "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
''' \
            + buff.getvalue()
    

class Track(TrackBase):
    def __init__(self, width=0, height=10):
        self.width = width
        self.height = height
        self.scale = (1., 1.)

    def draw(self, b, scale):
        print "drawing %s: %s"%(self.__class__.__name__,self.name)
        self.text(b, self.name, x=-5*scale[0], y=6, anchor='end', scale=scale)
        
    def _sty(self, color):
        if color:
            return {'style':'fill:%s;'%color}
        else:
            return {}
        
    def draw_vline(self, b, x, start, end, color, stroke=1, scale=(1.,1.)):
        with b.g(transform='translate(%s,0) scale(%s,%s)'%(x,scale[0],scale[1])):
            b.line(x1=0, y1=start, x2=0, y2=end, stroke=color)(**{'stroke-width':stroke})

    def draw_vgraph(self, b, x, start, end, value, color):
        f = (end-start)*(1.-value)
        b.line(x1=x, y1=f, x2=x, y2=end, stroke=color)

    def draw_hline(self, b, start, end, y, color):
        b.line(x1=start, y1=y, x2=end, y2=y, stroke=color)
        
    def draw_hbar(self, b, start, end, y, thick, color):
        w = end-start
        b.rect(x=start, y=y-thick/2, width=w, height=thick, **self._sty(color))

    def text(self, b, text, x, y, color='black', anchor='start', scale=(1.,1.)):
        with b.g(transform='translate(%s,%s)'%(x,y)):
            with b.g(transform='scale(%s,%s)'%scale):
                with b['text'](x=0, y=0)(**{'font-size':'10', 'text-anchor':anchor, 'style':'fill:%s;'%color}):
                    b.text(text)

class NamedTrack(Track):
    def __init__(self, name, width=0, height=10):
        super(NamedTrack, self).__init__(width, height)
        self.name_ = name
        
    @property
    def name(self):
        return self.name_


"""
misc
"""
class HbarTrack(NamedTrack):
    def __init__(self, name, length):
        super(HbarTrack, self).__init__(name, length, 5)

    def draw(self, b, scale):
        super(HbarTrack, self).draw(b, scale)
        self.draw_hbar(b, 0, self.width, 2, 2, color='black')

class MeasureTrack(Track):
    def __init__(self, length, start=0, step=500, substep=100, subsubstep=10):
        super(MeasureTrack, self).__init__(length, 20)
        self.length = length
        self.start = start
        self.step = step
        self.substep = substep
        self.subsubstep = subsubstep
        
    def draw(self, b, scale):
        self.draw_hline(b, 0, self.length, 20, color='black')
        # step
        ss = self.step * int(ceil(1.0*self.start/self.step))
        for i in range(ss, self.start + self.length, self.step):
            x = i - self.start
            self.draw_vline(b, x, 10, 20, color='black', scale=scale)
            self.text(b, str(i), x=x-5, y=10, color='black', anchor='middle', scale=scale)
        # substep
        ss = self.substep * int(ceil(1.0*self.start/self.substep))
        for i in range(ss, self.start + self.length, self.substep):
            x = i - self.start
            self.draw_vline(b, x, 13, 20, color='black', scale=scale)

        ss = self.subsubstep * int(ceil(1.0*self.start/self.subsubstep))
        for i in range(ss, self.start + self.length, self.subsubstep):
            x = i - self.start
            self.draw_vline(b, x, 18, 20, color='black', scale=scale)


"""
sequences
"""
class SequenceTrack(TrackGroup):
    def __init__(self, seq, features, start=0):
        super(SequenceTrack, self).__init__()
        self.seq = seq
        self.features = features
        self.length = len(self.seq)

        msize = 500*(max(1,int(self.length/10)/500))

        self.add(MeasureTrack(self.length, start, msize))
        self.add(CpgObsPerExpTrack(self.seq))
        self.add(GcPercentTrack(self.seq))
        #self.add(SequenceTrack(self.seq))
        for f in self.features:
            self.add(FeatureTrack(f))
        self.add(CpgBarTrack(self.seq))

    def draw(self, b, scale):
        super(SequenceTrack, self).draw(b, scale)

class TranscriptTrack(TrackGroup):
    def __init__(self, name, seq, feature):
        super(TranscriptTrack, self).__init__()
        self.seq = seq
        self.length = len(self.seq)

        msize = 500*(max(1,int(self.length/10)/500))

        self.add(MeasureTrack(self.length, 0, msize))
        self.add(ExonTrack(name, feature))

    def draw(self, b, scale):
        super(TranscriptTrack, self).draw(b, scale)

class ExonTrack(NamedTrack):
    def __init__(self, name, feature):
        super(ExonTrack,self).__init__(name, 0, 10)
        self.feature = feature
    def draw(self, b, scale):
        super(ExonTrack,self).draw(b,scale)
        def loc_len(loc):
            return loc.nofuzzy_end - loc.nofuzzy_start
        exons = [loc_len(sf.location) for sf in self.feature.sub_features]

        count = 0
        for e in exons:
            self.draw_vline(b, count, 0, 10, color='black')
            count += e
        self.draw_vline(b, count, 0, 10, color='black')
        self.draw_hline(b, 0, count, 5, color='black')

class FeatureTrack(NamedTrack):
    def __init__(self, feature):
        self.feature = feature
        t = self.feature.type
        try:
            if t in ['CDS','mRNA']:
                name = "%s: %s"%(t,self.feature.qualifiers['product'][0])
            elif t=='gene':
                name = 'gene: %s'%self.feature.qualifiers['gene'][0]
            elif t=='source':
                organ = self.feature.qualifiers['organism'][0]
                ch = self.feature.qualifiers['chromosome'][0]
                name = '%s chromosome %s'%(organ,ch)
            else:
                name = t
        except KeyError:
            name = self.feature.type
        super(FeatureTrack, self).__init__(name)

    def draw(self, b, scale):
        super(FeatureTrack, self).draw(b, scale)
        h = 4
        loc = self.feature.location
        self.draw_hline(b, loc.nofuzzy_start, loc.nofuzzy_end, self.height/2, color='black')
        for sf in self.feature.sub_features:
            self.draw_hbar(b, sf.location.nofuzzy_start, sf.location.nofuzzy_end, self.height/2, 4, color='blue')



class AtgcColorTrack(NamedTrack):
    def __init__(self, seq):
        super(AtgcColorTrack, self).__init__('ATGC',len(seq))
        self.seq = seq

    def draw(self, b, scale):
        super(AtgcColorTrack, self).draw(b, scale)
        s = str(self.seq)
        for i,c in enumerate(self.seq):
            color = base_color(c)
            self.draw_vline(b, i, 0, self.height, color=color)
            

class CpgBarTrack(NamedTrack):
    def __init__(self, seq):
        super(CpgBarTrack, self).__init__('CpG site',len(seq),10)
        self.length = len(seq)
        self.cpg = cpg_sites(seq)

    def draw(self, b, scale):
        super(CpgBarTrack, self).draw(b, scale)
        self.draw_hline(b, 0, self.length, self.height/2, color='black')
        for i in self.cpg:
            self.draw_vline(b, i, 0, self.height, color='black', scale=scale)

class DbtssTrack(NamedTrack):
    def __init__(self, tssfile, seq):
        self.tssfile = tssfile
        super(DbtssTrack, self).__init__(tssfile.name, len(seq), 100)

    def draw(self, b, scale):
        super(DbtssTrack, self).draw(b,scale)
        h = self.height
        self.draw_hline(b, 0, self.width, h, color='gray')
        for x,v in self.tssfile.items():
            v = v/4.
            if v < h:
                self.draw_vline(b, x, h-v, h, color='red', scale=scale)
            else:
                st = 1.*v/h
                self.draw_vline(b, x, 0, h, color='blue', stroke=st, scale=scale)
            
def window_search(seq, window, step=1):
    h = int(window/2)
    l = len(seq)
    for i in xrange(0,l,step):
        yield i, seq[max(0,i-h):min(i+h,l)]


class SeqWindowGraphTrack(NamedTrack):
    def __init__(self, seq, name, calc, maxvalue, window, threshold):
        super(SeqWindowGraphTrack, self).__init__(name,len(seq),20)
        self.seq = seq
        self.window = window
        self.threshold=threshold
        self.calc = calc
        self.maxvalue = maxvalue
        
    def draw(self, b, scale):
        super(SeqWindowGraphTrack, self).draw(b, scale)

        l = len(self.seq)
        h = self.height
        
        def trans(v):
            return h*(1.-v/self.maxvalue)

        step = max(1, int(scale[0]))

        values = ((i, self.calc(subseq)) for i,subseq in window_search(self.seq, self.window, step))
        points = ', '.join("%.2f %.2f"%(i,trans(c)) for i,c in values)
        b.polyline(points=points, stroke='black', fill='none')

        t = trans(self.threshold)
        b.line(x1=0, y1=t, x2=l, y2=t, stroke='red')(**{'stroke-width':0.5,'stroke-dasharray':'30,10'})

class GcPercentTrack(SeqWindowGraphTrack):
    def __init__(self, seq, window=200):
        f = lambda subseq: GC(subseq)/100.
        super(GcPercentTrack, self).__init__(seq, 'GC %', f, 1., window, 0.5)

class CpgObsPerExpTrack(SeqWindowGraphTrack):
    def __init__(self, seq, window=200):
        f = lambda subseq: cpg_obs_per_exp(subseq)
        super(CpgObsPerExpTrack, self).__init__(seq, 'Obs/Exp CpG', f, 1.5, window, 0.65)


"""
pcrs
"""

class BasePcrsTrack(NamedTrack):
    def __init__(self, name, products):
        self.products = []
        fs = FreeSpace()
        for n, p in products:
            y = fs.set(p.start, p.end)
            self.products.append( (n, p, y*20) )
        super(BasePcrsTrack, self).__init__(name, 0, 20*fs.num_lines())

    def draw(self, b, scale):
        super(BasePcrsTrack, self).draw(b, scale)
        for name, p, y in self.products:
            self.draw_hbar(b, p.start, p.start_i, y, 4, color='red')
            self.draw_hbar(b, p.start_i, p.end_i, y,4, color='gray')
            self.draw_hbar(b, p.end_i, p.end, y, 4, color='blue')
            m = (p.start+p.end)/2.
            self.text(b, name, x=m, y=y+11, color='black', anchor='middle', scale=scale)
    
class PcrsTrack(BasePcrsTrack):
    def __init__(self, name, pcrs):
        def d():
            for pcr in pcrs:
                for p in pcr.products:
                    yield pcr.name, p
        super(PcrsTrack, self).__init__(name, d())

class BsPcrsTrack(BasePcrsTrack):
    def __init__(self, name, pcrs):
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

        super(BsPcrsTrack, self).__init__(name, d())


class PcrTrack(NamedTrack):
    def __init__(self, pcr):
        super(PcrTrack, self).__init__(pcr.name, 0, 20)
        self.pcr = pcr

    def draw(self, b, scale):
        super(PcrTrack, self).draw(b, scale)
        for i,p in enumerate(self.pcr.products):
            self.draw_hbar(b, p.start, p.start_i, 4, 4, color='red')
            self.draw_hbar(b, p.start_i, p.end_i, 4,4, color='gray')
            self.draw_hbar(b. p.end_i, p.end, 4, 4, color='blue')
            
            self.text(b, self.pcr.name, x=(p.start+p.end)/2, y=18, color='black', anchor='middle', scale=scale)


"""
primers
"""
class PrimersTrack(Track):
    def __init__(self, primers, template):
        fs = FreeSpace()
        self.primers = []
        for p in primers:
            f,r = p.search(template)
            for ff in f:
                start,end = ff, ff+len(p)-1
                y = fs.set(start,end)
                self.primers.append( (p.name, start, end, y, True) )
            for rr in r:
                start,end = rr, rr+len(p)-1
                y = fs.set(start,end)
                self.primers.append( (p.name, start, end, y, False) )
            
        super(PrimersTrack, self).__init__(0, fs.num_lines() * 20)

    def draw(self, b, scale):
        super(PrimersTrack, self).draw(b, scale)

        for name, start, end, y, forward in self.primers:
            self.draw_hbar(b, start, end, y*20, 3, color = '#f00' if forward else '#00f')
            self.text(b, name, x=(start+end)/2., y=y*20+10, color='black', anchor='middle',scale=scale)

class PrimerTrack(NamedTrack):
    def __init__(self, name, primer, template):
        super(PrimerTrack, self).__init__(name, 0, 10)
        self.primer = primer
        self.template = template
    def draw(self, b, scale):
        super(PrimerTrack, self).draw(b, scale)

        def d(start, end, forward=True):
            self.draw_hbar(b, start, end, 0, 10, color = '#f00' if forward else '#00f')
            
        f,r = self.primer.search(self.template)
        for ff in f:
            d(ff, ff+len(self.primer)-1, True)
        for rr in r:
            d(rr, rr+len(self.primer)-1, False)

class BisulfiteSequenceTrack(NamedTrack):
    def __init__(self, name, bsp_map, start, end):
        super(BisulfiteSequenceTrack, self).__init__(name, 0, 10)
        self.bsp_map = bsp_map
        self.start = start
        self.end = end

    def draw(self, b, scale):
        super(BisulfiteSequenceTrack, self).draw(b, scale)

        color = {'M':'#F00','U':'#00F','P':'#AA0'}
        
        self.draw_hline(b, self.start, self.end, 5, color='black')
        for index, result in self.bsp_map:
            if result != '?':
                self.draw_vline(b, index, 0, 10, color=color[result])
    
class SeqView(object):
    def __init__(self, genbank):
        self.genbank = genbank
        self.seq = genbank.seq
        self.length = len(self.seq)
        #self.bm_template = bisulfite(self.seq)
        #self.bu_template = bisulfite(self.seq)
        self.pcrs = []
        #self.bs_pcrs = []

    def add_pcr(self, name, fw, rv):
        self.pcrs.append(PCR(name, self.seq, fw, rv))
        
    def svg(self, width):
        t = TrackGroup()
        t.add(SequenceTrack(self.genbank.seq, self.genbank.features))
        t.add(PcrsTrack('genome PCR', self.length, self.pcrs))
        return t.svg(width)

if __name__=='__main__':
    import sys
    print sys.argv[1]
    genbank = SeqIO.read(open(sys.argv[1],'r'),'genbank')
    sv = SeqView(genbank)
    with open('test.svg','w') as f:
        f.write(sv.svg(1200))
