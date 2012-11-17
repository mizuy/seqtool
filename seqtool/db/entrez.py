from __future__ import absolute_import

from Bio import Entrez

import xml.etree.ElementTree as ET
import os

from .locus import Locus
from . import hgnc

directory = os.path.dirname(os.path.abspath(__file__))
default_cache = os.path.join(directory,'../../_cache/')

gene_table = hgnc.GeneTable()

Entrez.email = 'mizugy@gmail.com'

class EntrezEfetchError(Exception):
    def __init__(self, m):
        self.m = str(m)
    def __str__(self):
        return 'Entrez.efetch Error: %s"'%self.m

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
            print "Entrez.efetch..: %s"%(args)
            handle = Entrez.efetch(**args)
            ret = handle.read()
        except:
            raise EntrezEfetchError(args)
        finally:
            print "..efetch finished."
        
        with open(fname, 'w') as f:
            f.write(ret)

        return ret

entrez = CachedEntrez()

def settings(text, sp=','):
    return [[x.strip() for x in c.split(sp)] for c in text.strip().split('\n')]

# "GRCh37.p10 Primary Assembly"
# name, acession, GI
t_chromosome = """
chr1,  NC_000001, 224589800
chr2,  NC_000002, 224589811
chr3,  NC_000003, 224589815
chr4,  NC_000004, 224589816
chr5,  NC_000005, 224589817
chr6,  NC_000006, 224589818
chr7,  NC_000007, 224589819
chr8,  NC_000008, 224589820
chr9,  NC_000009, 224589821
chr10, NC_000010, 224589801
chr11, NC_000011, 224589802
chr12, NC_000012, 224589803
chr13, NC_000013, 224589804
chr14, NC_000014, 224589805
chr15, NC_000015, 224589806
chr16, NC_000016, 224589807
chr17, NC_000017, 224589808
chr18, NC_000018, 224589809
chr19, NC_000019, 224589810
chr20, NC_000020, 224589812
chr21, NC_000021, 224589813
chr22, NC_000022, 224589814
chrX,  NC_000023, 224589822
chrY,  NC_000024, 224589823
"""
l_chromosome = settings(t_chromosome)

class Chromosome(object):
    def __init__(self):
        self.name = {}
        for ch,ac,gid in l_chromosome:
            self.name[ac] = ch

    def get_name(self, accession):
        return self.name[accession]

t_chromosome = Chromosome()

class GeneXml(object):
    def __init__(self, xml):
        self.xml = ET.fromstring(xml)

        try:
            # Entrezgene-Set / Entrezgene / Entrezgene_locus / Gene-commentary[0] / Gene-commentary_seqs / Seq-loc / Seq-loc-int / Seq-interval /
            #     Seq-interval_from
            #     Seq-interval_to
            #     Seq-interval_strand
            #     Seq-interval_id / Seq-id / Seq-id_gi


            accession = self.xml.findall(".//Entrezgene/Entrezgene_locus/Gene-commentary[1]/Gene-commentary_accession")[0].text
            chromosome = t_chromosome.get_name(accession)
            a = self.xml.findall(".//Entrezgene/Entrezgene_locus/Gene-commentary[1]/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval")[0]
            start = int(a.findall('Seq-interval_from')[0].text)
            stop = int(a.findall('Seq-interval_to')[0].text)
            sense = True if (a.findall('Seq-interval_strand/Na-strand')[0].attrib['value'] == "plus") else False

            id_gi = a.findall('Seq-interval_id/Seq-id/Seq-id_gi')[0].text

            lower, higher = start, stop

            self.locus = Locus(chromosome, sense, lower, higher)
            self.locus.id_gi = id_gi
        finally:
            pass

def get_genomic_context_genbank(gene_text):
    gene_id, gene_symbol = get_gene_from_text(gene_text)
    return get_genomic_context_genbank(gene_id)

def get_gene_from_text(text):
    """
    text must be gene-id integer OR gene symbol
    """
    try:
        gene_id = int(text)
        gene_symbol = gene_table.get_gene_symbol(gene_id)
    except ValueError:
        gene_symbol = text
        gene_id = gene_table.get_gene_id(gene_symbol)
    return gene_id, gene_symbol

def get_gene_locus(gene_id, upstream=1000, downstream=1000):
    xml = entrez.efetch(db='gene', id=gene_id, retmode='xml')
    g = GeneXml(xml)
    return g.locus.expand(upstream, downstream)

def get_genomic_context_genbank(gene_id, upstream=1000, downstream=1000):
    xml = entrez.efetch(db='gene', id=gene_id, retmode='xml')

    g = GeneXml(xml)
    locus = g.locus.expand(upstream, downstream)

    return entrez.efetch(db='nuccore',
                         id=g.locus.id_gi,
                         seq_start=locus.pos.start, seq_stop=locus.pos.stop,
                         strand=1 if locus.pos.sense else 2,
                         rettype='gb', retmode='text')
    

