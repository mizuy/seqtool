import os

from .locus import Locus
from . import sql
from . import entrez
from . import dbtss

class NoSuchGene(Exception):
    def __init__(self, m):
        self.gene_text = str(m)
    def __str__(self):
        return 'No Such Gene "%s"'%self.gene_text

def get_gene_symbol(gene_id):
    session = sql.Session()
    try:
        return session.query(sql.GeneTable).filter_by(gene_id=gene_id).one().symbol
    except:
        return None

def get_gene_id(gene_symbol):
    session = sql.Session()
    try:
        return session.query(sql.GeneTable).filter_by(symbol=gene_symbol).one().gene_id
    except:
        return None

def get_gene_from_text(text):
    try:
        gene_id = int(text)
        gene_symbol = get_gene_symbol(gene_id)
    except ValueError:
        gene_symbol = text
        gene_id = get_gene_id(gene_symbol)
    return gene_id, gene_symbol

def get_gene_locus(gene_id):
    symbol = get_gene_symbol(gene_id)
    session = sql.Session()
    # TODO: gene has multiple mRNA.
    # TODO: gene might have even multiple loci !!
    try:
        f = session.query(sql.UcscTable).filter_by(symbol=symbol).first()
        return Locus(f.chrom, f.strand=='+', f.txStart, f.txEnd)
    except:
        return None
    
def get_locus_genbank(locus):
    session = sql.Session()

    try:
        gid_chrom = session.query(sql.Chromosome).filter_by(name=locus.chrom).one().gid

        return entrez.entrez.efetch(db='nuccore',
                             id=gid_chrom,
                             seq_start=locus.pos.start, seq_stop=locus.pos.stop,
                             strand=1 if locus.pos.sense else 2,
                             rettype='gb', retmode='text')
    except:
        return None
    
def get_dbtss_tissues():
    return dbtss.tissues
