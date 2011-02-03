from __future__ import absolute_import

from collections import defaultdict
import subprocess
import sys, os, time

bowtie_path = os.path.expanduser('~/work/bowtie-0.12.7/bowtie')

def virtualpcr(template, seqstrs, threshold=None):
    '''
    bowtie manual: http://bowtie-bio.sourceforge.net/manual.shtml
    parser for default bowtie output
    '''

    command = '%s -a -n 2 %s -c %s' %(os.path.abspath(bowtie_path), template, ','.join(seqstrs))
    try:
        print command
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        sys.stdout.write('waiting for bowtie.')
        while p.poll()==None:
            sys.stdout.write('.')
            sys.stdout.flush()
            time.sleep(1)
        print 'finished.'
        retcode = p.wait()
        if retcode < 0:
            print >>sys.stderr, 'bowtie exist with return code ', retcode, '  command:', command
            return None
    except OSError, e:
        print >>sys.stderr, "failed to invoke bowtie. command: ", command, ' error=', e
        return None


    # matched-align list for each template(chromosome). positive means forward strand match.
    matches_positive = defaultdict(list)
    matches_negative = defaultdict(list)

    primer_full_count = defaultdict(int)
    primer_miss_count = defaultdict(int)

    for l in p.stdout:
        '''see http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output for output format. '''
        r = [c.strip() for c in l.split()]
        assert( int(r[0]) in [0,1] )

        primer_i = int(r[0])
        primer = seqstrs[primer_i]
        direction = r[1]
        template = r[2]
        offset = int(r[3])
        try:
            difference = r[7]
        except:
            difference = None

        if difference:
            primer_miss_count[primer_i] += 1
        else:
            primer_full_count[primer_i] += 1

        match = (primer, offset, difference)
        if direction == '+':
            matches_positive[template].append(match)
        else:
            matches_negative[template].append(match)

    print 'primer statistics'
    for i,p in enumerate(seqstrs):
        print 'primer %s:' % p
        print ' full match = %4d' % primer_full_count[i]
        print ' miss match = %4d' % primer_miss_count[i]

    for v in matches_positive.values():
        v.sort(key=lambda x:x[1])
    for v in matches_negative.values():
        v.sort(key=lambda x:x[1])

    templates = set(matches_positive.keys()) & set(matches_negative.keys())

    def pcr_pair(pl, nl):
        if len(pl) < 1 or len(nl) < 1:
            return
        p = 0
        n = 0
        while p < len(pl):
            while not pl[p][1] <= nl[n][1]:
                n += 1
                if not (n < len(nl)):
                    return
            while p+1 < len(pl) and pl[p+1][1] <= nl[n][1]:
                p += 1
            if not threshold or nl[n][1]-pl[p][1] < threshold:
                yield p, n
            p += 1
        

    for t in templates:
        pl = matches_positive[t]
        nl = matches_negative[t]
        for p,n in pcr_pair(pl, nl):
            pp = pl[p]
            nn = nl[n]
            l = nn[1] - pp[1]
            print 'match at %s length=%d(%d->%d) %s(%s) -> %s(%s)' %(t, l, pp[1],nn[1], pp[0],pp[2], nn[0],nn[2])

def main():
    import sys, os
    from optparse import OptionParser

    parser = OptionParser('usage: %prog [options] primer_sequence ...')
    parser.add_option("-t",'--template', dest='template', help="name of indexes of template", default=None)
    parser.add_option('-b', '--bisulfite', dest='bisulfite', action='store_true', default=False)
    parser.add_option('-l', '--length', dest='length', type=int, default=-1)

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error("no sequence input")

    if options.bisulfite and options.template:
        parser.error('bisulfite option and template option are exclusive.')
    
    length = options.length if options.length > 0 else None
    if options.template:
        virtualpcr(options.template, args, length)
    elif options.bisulfite:
        print 'METHYL'
        virtualpcr('bs_sense_met', args, length)
        virtualpcr('bs_antisense_met', args, length)
        print 'UNMETHYL'
        virtualpcr('bs_sense_unmet', args, length)
        virtualpcr('bs_antisense_unmet', args, length)
    else:
        parser.error('no template specified.')

if __name__=='__main__':
    main()
