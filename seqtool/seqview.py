from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from collections import defaultdict

from seqtool.seqcanvas import SeqCanvas, SeqCanvasOverview, Point
from seqtool.nucleotide import Primer, PCR, bisulfite, is_cpg, search_primer, tm_gc, count_cpg, is_repeat, base_color

def load_primer_list_file(fileobj):
    primer = []
    for name, value, em in parse_file(fileobj):
        primer.append(Primer(name, Seq.Seq(value.upper(),IUPAC.ambiguous_dna)))
    return primer

def parse_file(fileobj):
    lineno = 0
    for l in fileobj:
        lineno += 1
        l = l.strip()
        if not l or l.startswith('#') or l.startswith('//'):
            continue
        def error_message(msg):
            print ':%s: %s: "%s"'%(lineno,msg,l)
        ls = l.split(':')
        if len(ls)!=2:
            error_message('unknown line')
            continue
        name = ls[0].strip()
        value = ls[1].strip()
        yield name,value,error_message

def parse_structured_file(fileobj):
    lineno = 0
    category = []
    parameter = {}
    rank=0
    skipping = False
    for l in fileobj:
        lineno += 1
        l = l.strip()
        if not l or l.startswith('#') or l.startswith('//'):
            continue
        def error_message(msg):
            print ':%s: %s: "%s"'%(lineno,msg,l)
        if l.startswith('/+'):
            skipping = True
        elif l.startswith('+/'):
            skipping = False
        else:
            if skipping:
                continue
            if l.startswith('>'):
                rank = len(l)-len(l.lstrip('>'))
                category = category[0:rank-1]+[l.strip('>').strip()]
            elif l.startswith('@'):
                ls = l.lstrip('@').split('=')
                if len(ls)!=2:
                    error_message('unkown parameter line')
                    continue
                parameter[ls[0]]=ls[1]
            else:
                ls = l.split(':')
                if len(ls)!=2:
                    error_message('unknown line')
                    continue
                name = ls[0].strip()
                value = ls[1].strip()
                yield category,parameter,name,value,error_message

# mask for pcr bars and primers
class MaskMap(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.length = end-start
        self.mask = []

    def get(self, start, end):
        p = max(0,start-self.start)
        q = min(self.length,end-self.start)
        if p>q:
            raise 'error'
        else:
            rows = len(self.mask)
            for y in range(rows):
                if not any(self.mask[y][p:q]):
                    self._set(p,q,y)
                    return y
            self.mask.append([0]*self.length)
            self._set(p,q,rows)
            return rows

    def _set(self, start, end, y):
        for x in range(start,end):
            self.mask[y][x] = 1

    def getall(self):
        return len(self.mask)

class ListDict(object):
    def __init__(self):
        self.list = []
        self.dict = {}
    def append(self, value):
        if not self.dict.has_key(value.name):
            self.list.append(value)
        self.dict[value.name] = value
    def __len__(self):
        return len(self.list)

class GenBankAnnotated(object):
    def __init__(self,name=None):
        self.record = None

        self.template_seq = {}
        self.templates = {
            'origin': ('original sequence',None),
            'bs_met': ('Sodium-Bisulfite(Methyl)',lambda s:bisulfite(s,True)),
            'bs_unmet': ('Sodium-Bisulfite(Unmethyl)',lambda s:bisulfite(s, False)),
            }

        self.primer = ListDict()
        self.pcr = defaultdict(ListDict)
        self.pcr[None] = ListDict()

        self.bisulfite_sequence = defaultdict(list)
        self.celllines = []
        self.pattern = ListDict()
        
        self.view_start = 0
        self.view_end = 0
        
        self.cpg_location = None
        self.name = name
        self.general = {}
        
    def load_genbank(self, genbank_record):
        self.record = genbank_record
        self.seq = self.record.seq
        self.length = len(self.seq)
        self.template_seq['origin'] = self.seq
        self.view_end = len(self.record)
        
        self.name = self.record.description

    
        
    def get_templates(self):
        for name,(desc,f) in self.templates.items():
            yield name, desc
    def get_pcrs(self, template):
        for pcr in self.pcr[template].list:
            yield pcr

    def add_primer(self, primer):
        self.primer.append(primer)

    def add_pcr(self, name, fw_primer, rv_primer, tmp_name='origin'):
        if not self.record:
            raise ValueError('load genbank first')
        assert isinstance(name,str)
        assert isinstance(fw_primer,str)
        assert isinstance(rv_primer,str)
        assert (not tmp_name) or isinstance(tmp_name,str)

        if not self.primer.dict.has_key(fw_primer):
            ValueError('no such forward primer: %s'%fw_primer)
        if not self.primer.dict.has_key(rv_primer):
            raise ValueError('no such reverse primer: %s'%rv_primer)

        if tmp_name and (not self.templates.has_key(tmp_name)):
            raise ValueError('no such pcr template: %s'%tmp_name)
        if not self.templates.has_key(tmp_name):
            raise ValueError('no such template: %s'%tmp_name)

        tmp = self._get_template_seq(tmp_name)
        
        self.pcr[tmp_name].append(PCR(name,tmp,self.primer.dict[fw_primer],self.primer.dict[rv_primer]))

    def _get_template_seq(self,tmp_name):
        if not self.template_seq.has_key(tmp_name):
            desc, func = self.templates[tmp_name]
            self.template_seq[tmp_name] = func(self.seq)
        return self.template_seq[tmp_name]

    def add_bisulfite_sequence_result(self, cellline, pcr_name, result):
        if not self.record:
            raise ValueError('load genbank first')
        assert isinstance(cellline,str)
        assert isinstance(pcr_name,str)
        assert isinstance(result,str)
        class BisulfiteSequenceResult(object):
            pass
        pcr_bs_met = self.pcr['bs_met']
        bs = BisulfiteSequenceResult()
        bs.cellline = cellline
        if not pcr_bs_met.dict.has_key(pcr_name):
            raise ValueError('no such pcr: %s'%pcr_name)
        bs.pcr = pcr_bs_met.dict[pcr_name]
        bs.result = result
        p = bs.pcr.get_products()
        if not len(p)==1:
            raise ValueError('%s must have 1 pcr product but %s'%(pcrname,len(products)))
        bs.product = p[0]
        if not all(i in 'MUP?' for i in result):
            raise ValueError('bsp result must contain only M,U,P or ?')
        cpg = bs.product.detectable_cpg()
        if len(result)!=bs.product.detectable_cpg():
            raise ValueError('%s has %s detectable CpG sites, but result gives %s'%(pcr_name,cpg,len(result)))
        self.bisulfite_sequence[cellline].append(bs)
        if not cellline in self.celllines:
            self.celllines.append(cellline)

    def add_pattern(self, name, pattern_seq):
        class PatternSeq(object):
            pass
        ps = PatternSeq()
        ps.name = name
        ps.seq = pattern_seq
        self.pattern.append(ps)

    def has_bisulfite_sequence_result(self):
        return len(self.bisulfite_sequence)>0
    def has_bisulfite_pcr(self):
        return len(self.pcr['bs_met'])>0 or len(self.pcr['bs_unmet'])>0

    def render_overview(self,filename):
        if not self.record:
            raise ValueError('load genbank first')

        if not self.cpg_location:
            self.cpg_location = self._calc_cpg_location()

        show_bs = self.has_bisulfite_pcr()
        w = SeqCanvasOverview(self.view_start, self.view_end)

        y = 1
        y = self._draw_genbank_features(w, y)

        # sequence
        w.put_line(self.view_start, self.view_end, y)
        for i in range(self.view_start,self.view_end):
            w.put_yline(i,y+1,y+1.5,base_color(str(self.seq[i])))
            if show_bs:
                seq_met = self._get_template_seq('bs_met')
                w.put_yline(i,y+1.5,y+2,base_color(str(seq_met[i])))
                if is_repeat(seq_met,i,6):
                    w.put_yline(i,y+2,y+2.3,base_color(str(seq_met[i])))
        if show_bs:
            for i in self.cpg_location:
                w.put_yline(i,y-0.5,y+0.5,'#000')
        y += 4

        # draw bisulfite sequences
        color = {'M':'#F00','U':'#00F','P':'#AA0'}
        for cl in self.celllines:
            if cl not in self.bisulfite_sequence.keys():
                continue
            w.put_line(self.view_start, self.view_end, y)
            w.put_char(Point(self.view_start, y-0.5), cl)
            for bs in self.bisulfite_sequence[cl]:
                for x,i in enumerate(bs.product.cpg_sites()):
                    result = bs.result[x]
                    if result!='?':
                        w.put_yline(i,y-.5,y+.5,color[result])
                    continue
                continue
            y += 2
            continue
        y += 1

        for name,(desc,f) in self.templates.items():
            seq = self._get_template_seq(name)
            pcr = self.pcr[name]
            w.put_line(self.view_start, self.view_end, y, color='#999')
            w.put_char(Point(self.view_start, y), name)
            y += 1
            y = self._draw_pattern(w, y, seq)
            y = self._draw_pcr(w, y, pcr.list)
            y = self._draw_primer(w, y, seq, pcr.list)
        w.put_line(self.view_start, self.view_end, y, color='#999')

        w.save(filename, format='PNG', height=y)
        return True
    
    def _calc_cpg_location(self):
        cpgl = {}
        count = 0
        for i in range(self.view_start,self.view_end):
            if is_cpg(self.seq,i):
                cpgl[i] = count
                count += 1
        return cpgl

    def render_bsp_overview(self, filename, scale=2, init=50):
        '''
        TODO: fix init
        '''
        if not self.record:
            raise ValueError('load genbank first')
        if not self.has_bisulfite_sequence_result():
            raise ValueError('no bisulfite sequence result')

        if not self.cpg_location:
            self.cpg_location = self._calc_cpg_location()

        count = len(self.cpg_location)
        def cx(x):
            return init+x*scale
        w = SeqCanvasOverview(0, cx(count))

        color = {'M':'#F00','U':'#00F','P':'#AA0'}
        
        y = 1
        for cl in self.celllines:
            w.put_line(init, cx(count), y)
            w.put_char(Point(0, y-0.5), cl)
            for bs in self.bisulfite_sequence[cl]:
                for x,i in enumerate(bs.product.cpg_sites()):
                    loc = cx(self.cpg_location[i])
                    result = bs.result[x]
                    if not result=='?':
                        w.put_yline(loc,y-.5,y,color[result])
            y += 1

        w.save(filename, format='PNG',height=y)
        return True

    def render_sequence(self,filename):
        if not self.record:
            raise ValueError('load genbank first')

        show_bs = self.has_bisulfite_pcr()
        w = SeqCanvas(self.view_start, self.view_end, 25, self.view_end)
        y = 1

        # sequence
        for i in range(self.view_start,self.view_end):
            cpg = is_cpg(self.seq,i)
            col = '#000' if not cpg else '#FF0'
            bgcol = '#fff' if not cpg else '#000'
            w.put_char(Point(i,y), self.seq[i], col, bgcol)
            if show_bs:
                seq_met = self._get_template_seq('bs_met')
                seq_unmet = self._get_template_seq('bs_unmet')
                w.put_char(Point(i,y+1), seq_met[i], col, bgcol)
                w.put_char(Point(i,y+2), seq_unmet[i], col, bgcol)
        y += 4

        for t in self.get_templates():
            y = self._draw_pcr(w, y, self.pcr[t].list)

        w.save(filename, format='PNG')
        return True

    def _draw_pcr(self, w, y, pcrs):
        # draw pcrs
        mask = MaskMap(self.view_start, self.view_end)
        for pcr in pcrs:
            for q in pcr.get_products():
                m = mask.get(q.start, q.end)
                w.show_pcr(q.start, q.end, y+m, q.primer_fw, q.primer_rv, '#050', '%s'%(pcr.name))
                #mask.set(q.start, q.end, m+1)
        return y+mask.getall()

    def _draw_pattern(self, w, y, seq):
        mask = MaskMap(self.view_start, self.view_end)
        for p in self.pattern.list:
            def draw(start, end, forward=True):
                m = mask.get(start, end)
                w.put_bar(start, end, y+m, '#f00' if forward else '#00f', p.name, rarrow=forward, larrow=not forward)
            f,r = search_primer(p.seq, seq)
            for ff in f:
                draw(ff, ff+len(p.seq)-1, True)
            for rr in r:
                draw(rr, rr+len(p.seq)-1, False)
        return y + mask.getall()

    def _draw_primer(self, w, y, template, exclude_list=[]):
        # draw primers
        mask = MaskMap(self.view_start, self.view_end)
        pmatch = {}
        for p in self.primer.list:
            def draw(start, end, forward=True):
                m = mask.get(start, end)
                w.put_bar(start, end, y+m, '#f00' if forward else '#00f', p.name, rarrow=forward, larrow=not forward)
            
            f,r = search_primer(p.seq, template)
            for ff in f:
                draw(ff, ff+len(p)-1, True)
            for rr in r:
                draw(rr, rr+len(p)-1, False)
        return y + mask.getall()

    def _draw_genbank_features(self, w, y):
        mask = MaskMap(self.view_start, self.view_end)
        for f in self.record.features:
            q = f.qualifiers
            name = f.type
                name += ': ' + q['product'][0] or ''
            if q.has_key('product'):
            l_p = f.location.nofuzzy_start
            l_q = f.location.nofuzzy_end
            yy = mask.get(l_p,l_q-1)
            w.put_bar(l_p,l_q-1, y+yy, name=name)

            for idx, sf in enumerate(f.sub_features):
                w.put_bar(sf.location.nofuzzy_start, sf.location.nofuzzy_end-1, y+yy, color='#ff0000')
        return y+mask.getall()

    def load_seqview(self, fileobj, relative_path):
        for category, parameter, name, value, em in parse_structured_file(fileobj):
            if category == ['general']:
                if self.general.has_key(name):
                    print 'overwriting general option %s=%s'%(name,value)
                self.general[name]=value
                if name=='genbank':
                    if self.record:
                        em('genbank file(%s) already loaded. skipped: %s'%(self.genbank,value))
                        continue
                    self.genbank = value
                    with open(relative_path(self.genbank)) as f:
                        print 'loading genbank: '+self.genbank,
                        self.load_genbank(SeqIO.read(f, "genbank"))
                        print '...done.'
                elif name=='primers':
                    self.primers_file = value
                    with open(relative_path(self.primers_file)) as f:
                        print 'loading primers: '+self.primers_file,
                        for p in load_primer_list_file(f):
                            self.add_primer(p)
                            print '.',
                        print 'done.'
            
            elif category == ['features']:
                self.add_pattern(name,Seq.Seq(value,IUPAC.ambiguous_dna))
            elif category == ['pcr']:
                name = name.split(',')[0].strip()
                ls = value.split(',')
                if len(ls)!=2:
                    em('you must specify 2 primer names separated by "," for each pcr: %s'%name)
                    continue

                if not parameter.has_key('template'):
                    em('you must set parameter template')
                    continue
                template = parameter['template']
                if not self.templates.has_key(template):
                    em('unkown template: "%s"'%(template))
                    continue
            
                self.add_pcr(name, ls[0].strip(), ls[1].strip(), template)
            elif category == ['bsp']:
                n = [n.strip() for n in name.split(',')]
                if not len(n)>=2:
                    em('each bsp must have at least 2 key; cell line name and pcr name')
                    continue
                pcrname = n[0].strip()
                cellline = n[1].strip().upper()
                annotations = n[2:]
                if not pcrname or not cellline:
                    em('empty pcr or cellline name: %s, %s'%(pcrname,cellline))
                    continue
                self.add_bisulfite_sequence_result(cellline, pcrname, value.strip().upper())
            else:
                em('unkown category: %s'%category)


def main():
    import sys, os
    from optparse import OptionParser
    from htmlwriter import HtmlWriter

    parser = OptionParser('usage: %prog [options] seqviewfile1.seqv ... -o outputfile.html')
    parser.add_option("-o", "--output", dest="output", help="output filename")
    
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error("no input file")
    
    if not options.output:
        parser.error("no output file")

    inputfiles = args
    outputfile = options.output
    output_basename = os.path.basename(outputfile).rpartition('.')[0]
    output_dir = os.path.abspath(os.path.dirname(outputfile))

    gbas = []
    for i in inputfiles:
        gba = GenBankAnnotated()
        print 'loading: ',i
        with open(i,'r') as f:
            def relative_path(name):
                return os.path.join(os.path.dirname(i),name)
            gba.load_seqview(f, relative_path)
        gbas.append(gba)

    with open(outputfile,'w') as f:
        print 'writing html...'
        html = HtmlWriter(f)
        html.push('html')
        html.push('head')

        html.push('style',type='text/css')
        html.text(\
'''
.image{border: solid 1px;}
.template{margin-left: 1em;}
.pcr{margin-left: 4em;}
''')
        html.pop()

        html.pop()
        html.push('body')

        for num,gba in enumerate(gbas):
            bn = output_basename + '%02d'%num
            methyl_r = bn + '__methyl.png'
            seqovw_r = bn + '__seqoverview.png'
            sequen_r = bn + '__sequence.png'
            
            html.insert('h1', gba.name or 'No Name')

            html.push('div',cls='image')
            if gba.has_bisulfite_sequence_result():
                html.insertc('img',src=methyl_r,width='600px')
                html.insertc('br')
            html.push('a',href=seqovw_r)
            html.insertc('img',src=seqovw_r,width='1000px')
            html.pop()

            html.pop()

            for template,desc in gba.get_templates():
                pcrs = list(gba.get_pcrs(template))
                if not len(pcrs):
                    continue
                html.push('div',cls='template')
                html.insert('h2','PCR products, template=%s'%desc)
                for pcr in pcrs:
                    html.push('div',cls='pcr')
                    pcr.debugprint_html(html)
                    html.pop()
                html.pop()
        
            def render(func, name):
                print 'rendering %s..'%name,
                func(os.path.join(output_dir,name))
                print '...done.'
            if gba.has_bisulfite_sequence_result():
                render(gba.render_bsp_overview,methyl_r)
            render(gba.render_overview,seqovw_r)
            render(gba.render_sequence,sequen_r)

        html.pop()
        html.pop()
        print 'done.'

if __name__=='__main__':
    main()
