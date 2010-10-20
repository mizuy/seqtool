
from Bio import SeqIO

import sys

def main():
    g = SeqIO.read(open('example/ndrg2.gb','r'),'genbank')
    for f in g.features:
        if f.type=='mRNA':
            ret = f

    return ret
