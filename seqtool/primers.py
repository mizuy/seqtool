from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

from .xmlwriter import XmlWriter, builder
from .pcr import primers_write_html, primers_write_csv, Primer
from .parser import parse_file
from .prompt import prompt

import sys

__all__ = ['load_primer_list_file']

def load_primer_list_file(fileobj):
    primer = []
    for name, value, em in parse_file(fileobj):
        primer.append(Primer(name, Seq.Seq(value.upper(),IUPAC.ambiguous_dna)))
    return primer

