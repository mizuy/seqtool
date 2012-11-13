from __future__ import absolute_import

import sys, os
from collections import defaultdict
from glob import glob
from . import xmlwriter
from .db import gene

directory = os.getcwd()
input_dir = os.path.join(directory,'input')
output_dir = os.path.join(directory,'_output')

gene_list = os.path.join(input_dir, 'gene_list.txt')
index_file = os.path.join(output_dir, 'index.html')

def input_file(name):
    return os.path.join(input_dir, name)
def output_file(name):
    return os.path.join(output_dir, name)

def seqv_template(gene_name):
    return """
>general
genbank: %s.gb
#primers: primers.txt

>motif
>tss
>tss_count
>pcr
>rt_pcr
>bs_pcr
>bsp
""" % gene_name

def gene_list_of_seqv_file():
    ret = []
    for f in glob(os.path.join(input_dir,'*.seqv')):
        basename = os.path.basename(f)
        gene_text = os.path.splitext(basename)[0]
        #gene_id, gene_symbol = gene.get_gene_from_text(gene_text)
        ret.append(gene_text)
    return ret

class Block(object):
    def __init__(self, name):
        self.name = name
        self.list = []
    def add(self, name):
        self.list.append(name)

class IndexFile(object):
    def __init__(self, fileobj):
        self.blocks = self._parse(fileobj)

        indexgenes = set()
        for b in self.blocks:
            indexgenes |= set(b.list)

        seqvgenes = set(gene_list_of_seqv_file())

        othergenes = seqvgenes - indexgenes
        self.seqvs = list(seqvgenes | indexgenes)
        
        o = Block('Other seqv files')
        o.list = list(othergenes)
        self.blocks.append(o)
                
    def write(self, fileobj):
        html = xmlwriter.XmlWriter(fileobj)
        b = xmlwriter.builder(html)
        with b.html:
            with b.head:
                pass
            with b.body:
                b.h1('Index of Gene List')
                
                for block in self.blocks:
                    b.h2(block.name)
                    with b.ul:
                        for g in block.list:
                            with b.li:
                                b.a(g, href=g+'.html')

    def _parse(self, fileobj):
        blocks = []
        lineno = 0
        skipping = False
        for l in fileobj:
            lineno += 1
            l = l.strip()

            if not l or l.startswith('#') or l.startswith('//'):
                continue

            if l.startswith('/+'):
                skipping = True
            elif l.startswith('+/'):
                skipping = False
            else:
                if skipping:
                    continue
                if l.startswith('>'):
                    blocks.append(Block(l[1:].strip()))
                else:
                    if not blocks:
                        blocks.append(Block('No Name'))
                    blocks[-1].add(l.strip())
        return blocks

def cmd_index(isopen):
    f = IndexFile(open(input_file('index.txt'),'r'))

    html = output_file('index.html')
    f.write(open(html,'w'))

    for seqv in f.seqvs:
        # make html
        pass
    

    if isopen:
        os.system("open %s"%html)

def cmd_add(gene_list):
    for g in gene_list:
        print "adding %s..." % g
        gene_id, gene_symbol = gene.get_gene_from_text(g)

        print "Gene ID is %s, Gene Symbol is %s"%(gene_id, gene_symbol)

        seqv_file = input_file(gene_symbol+'.seqv')
        if os.path.exists(seqv_file):
            print "Already exists. skipping...: ",seqv_file
        else:
            print "writing...: ",seqv_file
            open(seqv_file,'w').write(seqv_template(gene_symbol))

        genbank_file = input_file(gene_symbol+'.gb')
        if os.path.exists(genbank_file):
            print "Already exists. skipping...: ",genbank_file
        else:
            print "retrieving genbank..."
            t = gene.get_genomic_context_genbank(gene_id)
            print "writing...: ",genbank_file
            open(genbank_file,'w').write(t)

def cmd_list():
    for l in gene_list_of_seqv_file():
        print l

def main():
    import sys, os
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='seqvcmd', description='commands for maintain seqv file folder')
    parser.add_argument("command", help="index, add, or list", default='index')
    parser.add_argument("gene", nargs='*', help="list of name or id of gene")
    #parser.add_argument("-o", "--open", action='store_true', help="open file", default=True)

    args = parser.parse_args()

    if args.command=='index':
        cmd_index(True)
        return
    elif args.command=='add':
        cmd_add(args.gene)
        return
    elif args.command=='list':
        cmd_list()
        return
    else:
        print "No such command %s. exit." % args.command
        return

if __name__=='__main__':
    main()
