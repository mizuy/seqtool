from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from collections import defaultdict

import htmlwriter

from StringIO import StringIO
from seqtool.seqcanvas import SeqCanvas, SeqCanvasOverview, Point
from seqtool.nucleotide import Primer, PCR, bisulfite, is_cpg, search_primer, tm_gc, count_cpg, is_repeat, base_color, cpg_sites, cpg_obs_per_exp

class TrackBase(object):
    @property
    def name(self):
        return 'no name'

class TrackGroup(TrackBase):
    def __init__(self, name):
        self.name_ = name
        self.tracks = []
        self.ymargin = 4

    @property
    def name(self):
        return self.name_

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
        return sum(len(i.name) for i in self)

    def __str__(self):
        "TrackGroup(%s, %s"%(self.name, ','.join(i for i in self))
        
    def draw(self, b, scale):
        y = 0
        for i in self:
            with b.g(transform='translate(0,%d)'%y):
                i.draw(b, scale)
            y += (i.height + self.ymargin)

    def svg(self, dst_width=None):
        width = self.width
        height = self.height

        nl = (self.name_length+1)*7.0

        if not dst_width:
            dst_width = nl + width
        dst_height = height

            
        scalex = (dst_width-nl)/width
        view_w = nl/scalex + width
        view_h = height

        trans_sx = dst_width/view_w
        trans_sy = dst_height/view_h
        
        buff = StringIO()
        b = htmlwriter.builder(htmlwriter.HtmlWriter(buff))
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
        print "drawing %s"%self.name
        self.text(b, self.name, x=-5*scale[0], y=8, anchor='end', scale=scale)
        
    def _sty(self, color):
        if color:
            return {'style':'fill:%s;'%color}
        else:
            return {}
        
    def draw_vline(self, b, x, start, end, color, scale=(1.,1.)):
        with b.g(transform='translate(%s,0) scale(%s,%s)'%(x,scale[0],scale[1])):
            b.line(x1=0, y1=start, x2=0, y2=end, stroke=color)
    
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

class FeatureTrack(NamedTrack):
    def __init__(self, feature):
        self.feature = feature
        super(FeatureTrack, self).__init__(self.feature.type)

    def draw(self, b, scale):
        super(FeatureTrack, self).draw(b, scale)
        h = 4
        loc = self.feature.location
        self.draw_hline(b, loc.nofuzzy_start, loc.nofuzzy_end, self.height/2, color='black')
        for sf in self.feature.sub_features:
            self.draw_hbar(b, sf.location.nofuzzy_start, sf.location.nofuzzy_end, self.height/2, 4, color='blue')


class SequenceTrack(NamedTrack):
    def __init__(self, seq):
        super(SequenceTrack, self).__init__('sequence',len(seq))
        self.seq = seq

    def draw(self, b, scale):
        super(SequenceTrack, self).draw(b, scale)
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
            

def window_search(seq, window, step=1):
    h = int(window/2)
    l = len(seq)
    for i in xrange(0,l,step):
        yield i, seq[i-h:i+h]

class SeqWindowGraphTrack(NamedTrack):
    def __init__(self, seq, name, calc, maxvalue, window, threshold):
        super(SeqWindowGraphTrack, self).__init__(name,len(seq),10)
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
        b.line(x1=0, y1=t, x2=l, y2=t, stroke='red')(**{'stroke-width':0.1,'stroke-dasharray':'30,10'})
    

class GcPercentTrack(SeqWindowGraphTrack):
    def __init__(self, seq, window=100):
        f = lambda subseq: GC(subseq)/100.
        super(GcPercentTrack, self).__init__(seq, 'GC %', f, 1., window, 0.5)

class CpgObsPerExpTrack(SeqWindowGraphTrack):
    def __init__(self, seq, window=200):
        f = lambda subseq: cpg_obs_per_exp(subseq)
        super(CpgObsPerExpTrack, self).__init__(seq, 'CpG o/e', f, 1.5, window, 0.65)
            
class PcrTrack(NamedTrack):
    def __init__(self, pcr):
        super(PcrTrack, self).__init__(pcr.name, 0, 20)
        self.pcr = pcr

    def draw(self, b, scale):
        super(PcrTrack, self).draw(b, scale)
        for i,p in enumerate(self.pcr.get_products()):
            self.draw_hbar(b, p.start, p.start_i, 4, 4, color='red')
            self.draw_hbar(b, p.start_i, p.end_i, 4,4, color='gray')
            self.draw_hbar(b. p.end_i, p.end, 4, 4, color='blue')
            
            self.text(b, self.pcr.name, x=(p.start+p.end)/2, y=18, color='black', anchor='middle', scale=scale)

class MeasureTrack(Track):
    def __init__(self, length, zero, step):
        super(MeasureTrack, self).__init__(length, 18)
        self.length = length
        self.zero = zero
        self.step = step
        
    def draw(self, b, scale):
        for i in range(0, self.length, self.step):
            self.draw_vline(b, i, 10, 15, color='black', scale=scale)
            self.text(b, str(i), x=i-5, y=10, color='black', anchor='middle', scale=scale)


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
        t = TrackGroup(self.genbank.description)
        tt = TrackGroup('template')
        t.add(tt)

        msize = 500*(max(1,int(self.length/10)/500))

        tt.add(MeasureTrack(self.length, 0, msize))
        
        tt.add(CpgBarTrack(self.seq))
        tt.add(CpgObsPerExpTrack(self.seq))
        tt.add(GcPercentTrack(self.seq))
        #tt.add(SequenceTrack(self.seq))
        for f in self.genbank.features:
            tt.add(FeatureTrack(f))
        tt = TrackGroup('pcr')
        t.add(tt)
        for p in self.pcrs:
            tt.add(PcrTrack(p))

        print t.width, t.height
        return t.svg(width)

if __name__=='__main__':
    import sys
    print sys.argv[1]
    genbank = SeqIO.read(open(sys.argv[1],'r'),'genbank')
    sv = SeqView(genbank)
    with open('test.svg','w') as f:
        f.write(sv.svg(1200))
