import os

directory = os.path.dirname(os.path.abspath(__file__))
hgnc_tab_file = os.path.join(directory,'../../_db/hgnc.tab')
ucsc_tab_file = os.path.join(directory,'../../_db/ucsc_gene_table.tab')
                

from sqlalchemy import ForeignKey
from sqlalchemy import Column, Date, Integer, String, Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.sql import select, or_, and_
from . import Session, engine
from ..prompt import prompt
from collections import defaultdict

Base = declarative_base()

t_chromosome = [[x.strip() for x in c.split(',')] for c in """
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
""".strip().split('\n')]

class Chromosome(Base):
    __tablename__ = 'chromosome'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    accession = Column(String, unique=True)
    gid = Column(Integer, unique=True)

    def __init__(self, ch, ac, gid):
        self.name = ch
        self.accession = ac
        self.gid = gid

    @classmethod
    def load(cls):
        session = Session()
        for ch, ac, gid in t_chromosome:
            session.add(Chromosome(ch, ac, gid))
        session.commit()

class GeneTable(Base):
    __tablename__ = 'gene_table'
    
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, unique=True) # Gene ID
    symbol = Column(String, unique=True)

    def __init__(self, gene_id, symbol):
        self.gene_id = gene_id
        self.symbol = symbol

    @classmethod
    def load(cls, hgnc_tab_file=hgnc_tab_file):
        session = Session()
        with open(hgnc_tab_file, 'r') as f:
            lines = f.readlines()
            title = lines[0].strip().split('\t')
            assert(title[1]=="Approved Symbol")
            assert(title[2]=="Entrez Gene ID (mapped data supplied by NCBI)")
            for line in lines[1:]:
                l = line.strip().split('\t')
                if len(l) < 4:
                    continue
                try:
                    symbol = l[1]
                    gene_id = int(l[2])
                except:
                    continue
                session.add(GeneTable(gene_id, symbol))
        session.commit()

def get_gene_symbol(gene_id):
    session = Session()
    return session.query(GeneTable).filter_by(gene_id=gene_id).one().symbol

def get_gene_id(gene_symbol):
    session = Session()
    return session.query(GeneTable).filter_by(symbol=gene_symbol).one().gene_id

class UcscTable(Base):
    __tablename__ = 'ucsc_table'

    id = Column(Integer, primary_key=True)
    n_bin = Column(Integer)
    accession = Column(String)
    chrom = Column(String)
    strand = Column(String)
    txStart = Column(Integer)
    txEnd = Column(Integer)
    cdsStart = Column(Integer)
    cdsEnd = Column(Integer)
    exonCount = Column(Integer)
    exonStarts  = Column(String)
    exonEnds = Column(String)
    score = Column(Integer)
    symbol = Column(String)
    cdsStartStat = Column(String)
    cdsEndStat = Column(String)
    exonFrames  = Column(String)

    def __init__(self):
        pass

    @classmethod
    def load(cls, ucsc_tab_file=ucsc_tab_file):
        session = Session()
        with open(ucsc_tab_file, 'r') as f:
            lines = f.readlines()
            title = lines[0].strip().split('\t')
            for line in lines[1:]:
                l = line.strip().split('\t')
                assert len(l) == 16
                
                u = UcscTable()
                u.n_bin = l[0]
                u.accession = l[1]
                u.chrom = l[2]
                u.strand = l[3]
                u.txStart = int(l[4])
                u.txEnd = int(l[5])
                u.cdsStart = int(l[6])
                u.cdsEnd = int(l[7])
                u.exonCount = int(l[8])
                u.exonStarts = l[9]
                u.exonEnds = l[10]
                u.score = int(l[11])
                u.symbol = l[12]
                u.cdsStartStat = l[13]
                u.cdsEndStat = l[14]
                u.exonFrames = l[15]
                
                session.add(u)
        session.commit()

t_dbtss = [[x.strip() for x in c.split(':')] for c in """
Brain1: adult-tissue/ambion_brain.bed
Brain2: adult-tissue/ambion_brain2.bed
C-Brain: adult-tissue/clontech_brain.bed
Heart: adult-tissue/ambion_heart.bed
C-Heart: adult-tissue/clontech_heart.bed
Lung: adult-tissue/ambion_lung.bed
Breast: adult-tissue/ambion_breast.bed
Kidney: adult-tissue/ambion_kidney.bed
C-Kidney: adult-tissue/clontech_kidney.bed
Liver: adult-tissue/ambion_liver.bed
Colon: adult-tissue/ambion_colon.bed
Lymph: adult-tissue/ambion_lymph.bed
Adipose: adult-tissue/ambion_adipose.bed
Muscle: adult-tissue/ambion_muscle.bed
Thyroid: adult-tissue/ambion_thyroid.bed
Adrenal: adult-tissue/ambion_adrenal.bed
Ovary: adult-tissue/ambion_ovary.bed
Prostate: adult-tissue/ambion_prostate.bed
Testis: adult-tissue/ambion_testis.bed
""".strip().split('\n')]

database_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'../../_db/dbtss/'))

class DbtssStrand(Base):
    __tablename__ = 'dbtss_strand'
    id = Column(Integer, primary_key=True)

    tissue = Column(String)
    chrom = Column(String)
    sense = Column(Boolean)
    #__table_args__ = (UniqueConstraint('tissue','chrom','sense'),{})

    tags = relationship("DbtssTag")

    def __init__(self, tissue, ch, sense):
        self.tissue = tissue
        self.chrom = ch
        self.sense = sense

    @classmethod
    def get(cls, tissue, ch, sense):
        session = Session()
        return session.query(DbtssStrand).filter_by(tissue=tissue, chrom=ch, sense=sense).one()


class DbtssTag(Base):
    __tablename__ = 'dbtss_tag'

    id = Column(Integer, primary_key=True)
    tag_start = Column(Integer, index=True)
    tag_num = Column(Integer)
    strand_id = Column(Integer, ForeignKey('dbtss_strand.id'), index=True)

    def __init__(self, strand, start, num):
        self.strand_id = strand.id
        self.tag_start = start
        self.tag_num = num

    @classmethod
    def search(cls, tissue, ch, sense, start, stop):
        sid = DbtssStrand.get(tissue,ch,sense).id
        t = cls.__table__
        q = t.select(and_(t.c.strand_id==sid, t.c.tag_start >= start, t.c.tag_start <= stop))
        with engine.begin() as con:
            for v in con.execute(q):
                yield v.tag_start, v.tag_num

    @classmethod
    def load(self):
        for tissue, filename in t_dbtss:
            with prompt('loading %s..'%tissue) as p:
                with open(os.path.join(database_dir,filename), 'r') as fileobj:
                    DbtssTag.load_bed(tissue, fileobj)

    @classmethod
    def load_bed(self, tissue, fileobj):
        storage = defaultdict(list)
        for l in fileobj:
            ll = l.split()
            chrom = ll[0]
            start = int(ll[1])
            sense = ll[3].strip()=='+'
            num = int(ll[4])

            storage[(chrom,sense)].append((start,num))

        session = Session()
        for ch,se in storage.keys():
            strand = DbtssStrand(tissue, ch, se)
            session.add(strand)
        session.commit()

        st = DbtssStrand.__table__
        
        with engine.begin() as con:
            for (ch,se),v in storage.items():
                
                q = st.select().where(and_(st.c.tissue==tissue, st.c.chrom==ch, st.c.sense==se))
                strand, = con.execute(q)
                con.execute(DbtssTag.__table__.insert(),[{
                        'tag_start':a,
                        'tag_num':b,
                        'strand_id':strand.id} for a,b in v])

        
Base.metadata.create_all(engine)

def clear_all():
    con = engine.connect() 

    trans = con.begin() 

    for name, table in Base.metadata.tables.items(): 
        print table.delete() 
        con.execute(table.delete()) 

    trans.commit() 

def load_all():
    clear_all()
    session = Session()
    print "loading Chromosome..."
    Chromosome.load()
    print "loading GeneTable..."
    GeneTable.load()
    print "loading UcscTable..."
    UcscTable.load()
    session.commit()
    print "loading DbtssBed..."
    DbtssTag.load()
    print '...done.'
    session.commit()

