from __future__ import absolute_import

from Bio import Entrez
import os

directory = os.path.dirname(os.path.abspath(__file__))
default_cache = os.path.join(directory,'../../_cache/')

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


