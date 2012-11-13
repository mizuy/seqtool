from __future__ import absolute_import

from Bio import Entrez
from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

import xml.etree.ElementTree as ET
from collections import defaultdict
import sys, os
from cStringIO import StringIO

directory = os.path.dirname(os.path.abspath(__file__))
default_cache = os.path.join(directory,'../../_cache/')

from . import hgnc

gene_table = hgnc.GeneTable()

Entrez.email = 'mizugy@gmail.com'

class CachedEntrez(object):
    def __init__(self, cache_dir=default_cache):
        self.cache_dir = default_cache

    def cache_file(self, **args):
        return os.path.join(self.cache_dir,'&'.join(["%s=%s"%(k,v) for k,v in args.items()])+'.cache')
        
    def efetch(self, **args):
        fname = self.cache_file(**args)
        if os.path.exists(fname):
            with open(fname,'r') as f:
                return f.read()
        try:
            print "Entrez.efetch %s"%args
            handle = Entrez.efetch(**args)
            ret = handle.read()
        except:
            print "Entrez.efetch error:", sys.exc_info()[0]
            return None
        
        with open(fname, 'w') as f:
            f.write(ret)

        return ret


class GeneXml(object):
    def __init__(self, xml):
        self.xml = ET.fromstring(xml)

        class Locus(object):
            pass
        self.locus = Locus()
        try:
            # Entrezgene-Set / Entrezgene / Entrezgene_locus / Gene-commentary[0] / Gene-commentary_seqs / Seq-loc / Seq-loc-int / Seq-interval /
            #     Seq-interval_from
            #     Seq-interval_to
            #     Seq-interval_strand
            #     Seq-interval_id / Seq-id / Seq-id_gi


            a = self.xml.findall(".//Entrezgene/Entrezgene_locus/Gene-commentary[1]/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval")[0]
            self.locus.start = a.findall('Seq-interval_from')[0].text
            self.locus.end = a.findall('Seq-interval_to')[0].text
            self.locus.strand = True if (a.findall('Seq-interval_strand/Na-strand')[0].attrib['value'] == "plus") else False
            self.locus.id_gi = a.findall('Seq-interval_id/Seq-id/Seq-id_gi')[0].text
        except:
            print "XML parse error: ", sys.exec_info()[0]
    

def genomic_context_genbank(gene_id):
    entrez = CachedEntrez()

    g = GeneXml(entrez.efetch(db='gene', id=gene_id, retmode='xml'))

    return entrez.efetch(db='nuccore', id=g.locus.id_gi, seq_start=g.locus.start, seq_stop=g.locus.end, strand=g.locus.strand, rettype='gb', retmode='text')


def main():
    import sys, os
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='get_genbank', description='Retrieve Genbank from Gene ID or Gene Symbol')
    parser.add_argument("gene", help="Gene ID or Gene Symbol")
    parser.add_argument("-o", "--output", dest="output", help="output filename")

    args = parser.parse_args()

    if args.output:
        output = open(outputfile,'w')
    else:
        output = sys.stdout

    try:
        gene_id = int(args.gene)
    except ValueError:
        gene_id = gene_table.get_gene_id(args.gene)
        if not gene_id:
            print "No such gene: ",args.gene
            return

    genbank = genomic_context_genbank(gene_id)
    if not genbank:
        print "GenBank retrieve error: ", gene_id
    print >>output, genbank

if __name__=='__main__':
    main()

# irx5
#    print genome_context(10265)
