

from Bio import Entrez
import os

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

    def get_genbank(self, chrom_gid, locus):
        try:
            return self.efetch(db='nuccore',
                                 id=chrom_gid,
                                 seq_start=locus.pos.start, seq_stop=locus.pos.stop,
                                 strand=1 if locus.pos.sense else 2,
                                 rettype='gb', retmode='text')
        except:
            return None
