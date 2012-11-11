
from itertools import chain, imap, islice

"""
% dbtssconv -c 'chr14' -t 21483922 -f 21494935  _private/dbtss/adult-tissue/ambion_brain.bed > _private/ndrg2_brain.tss
% dbtssconv -c 'chr14' -t 21483922 -f 21494935  _private/dbtss/adult-tissue/ambion_colon.bed > _private/ndrg2_colon.tss
"""

def main():
    import sys, os
    from optparse import OptionParser

    parser = OptionParser('usage: %prog [options] dbtss_bed_filename')
    parser.add_option("-o", "--output", dest="output", help="output filename")
    parser.add_option("-c", "--chromosome", dest="chromosome", help="chromosome name")
    parser.add_option("-f", "--from", dest="start", type='int')
    parser.add_option("-t", "--to", dest="to", type="int")

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error("no input file")
        return

    bedfile = args[0]
    t_chromosome = options.chromosome
    t_start = options.start
    t_end = options.to

    if t_start == t_end:
        parser.error("bad 'from' and 'to'")
        return

    if t_start < t_end:
        t_strand = '+'
    else:
        t_strand = '-'
        t_start, t_end = t_end, t_start

    """
    print 'chromosome = ', t_chromosome
    print 'from = ', t_start
    print 'to = ', t_end
    print 'strand = ', t_strand
    """

    if options.output:
        output = open(options.output,'w')
    else:
        output = sys.stdout

    for line in open(bedfile):
        l = line.split()
        chromosome = l[0]

        if t_chromosome != chromosome:
            continue

        strand = l[3]

        if t_strand != strand:
            continue

        start = int(l[1])
        end = int(l[2])
        number = int(l[4])

        if t_start <= start and start <= t_end:
            if t_strand == '+':
                print >>output, start - t_start, start, number
            else:
                print >>output, t_end-start, start, number


"""
bed/
# Definition of each column (tab limited):
-chromosome
-TSStag_start
-tag_end
-strand
-Number of tags
"""
