

from Bio import Entrez
import os

def splits(line):
    return [x.strip() for x in line.split(',')]

chromid = { name: (ref,gid) for name,ref,gid in map(splits, '''
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
chrY,  NC_000024, 224589823'''.strip().splitlines()) }


class EntrezEfetchError(Exception):
    def __init__(self, m):
        self.m = str(m)
    def __str__(self):
        return 'Entrez.efetch Error: %s"'%self.m

class CachedEntrez(object):
    def __init__(self, cache_dir, email):
        print("Using email {} for entrez".format(email))
        self.cache_dir = cache_dir
        Entrez.email = email

    def cache_file(self, **args):
        return os.path.join(self.cache_dir,'&'.join(["%s=%s"%(k,v) for k,v in sorted(args.items())])+'.cache')
        
    def efetch(self, **args):
        fname = self.cache_file(**args)
        if os.path.exists(fname):
            with open(fname,'r') as f:
                return f.read()
        try:
            args = {str(k): str(v) for k,v in args.items()}
            print("Entrez.efetch..: %s"%(args))
            handle = Entrez.efetch(**args)
            ret = handle.read()
        except:
            raise EntrezEfetchError(args)
        finally:
            print("..efetch finished.")
        
        with open(fname, 'w') as f:
            f.write(ret)

        return ret

    def get_genbank(self, locus):
        accession, gid = chromid[locus.chrom]
        try:
            return self.efetch(db='nuccore',
                                 id=gid,
                                 seq_start=locus.pos.start, seq_stop=locus.pos.stop,
                                 strand=1 if locus.pos.sense else 2,
                                 rettype='gb', retmode='text')
        except:
            return None
