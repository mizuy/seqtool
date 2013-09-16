

from collections import defaultdict
import subprocess
import sys, os, time
import tempfile
from ..util.prompt import prompt


class PrimerAlign(object):
    def __init__(self, template, seq, mismatch, offset):
        self.template = template
        self.origin = seq
        self.mismatch = mismatch
        self.offset = offset


class PCRFragment(object):
    def __init__(self, template, fw, rv):
        self.template = template
        self.fw = fw
        self.rv = rv

    def __len__(self):
        return self.rv.offset - self.fw.offset + len(self.rv.origin)


def virtualpcr(template, seqstrs, bowtie_path, threshold=None):
    '''
    bowtie manual: http://bowtie-bio.sourceforge.net/manual.shtml
    parser for default bowtie output
    '''

    stdout = tempfile.TemporaryFile()
    stderr = tempfile.TemporaryFile()
    command = '%s %s -c %s' %(os.path.abspath(bowtie_path), template, ','.join(seqstrs))
    try:
        print(command)
        # note. using subprocess.PIPE with p.wait() dead locks if output is large.
        p = subprocess.Popen(command, stdout=stdout, stderr=stderr, shell=True)
        with prompt('waiting for bowtie.') as pr:
            while p.poll() is None:
                pr.progress()
                time.sleep(1)
        retcode = p.wait()
        if retcode < 0:
            print('bowtie exist with return code ', retcode, '  command:', command, file=sys.stderr)
            return None
    except OSError as e:
        print("failed to invoke bowtie. command: ", command, ' error=', e, file=sys.stderr)
        return None

    stdout.seek(0)
    stderr.seek(0)

    # matched-align list for each template(chromosome). positive means forward strand match.
    matches_positive = defaultdict(list)
    matches_negative = defaultdict(list)

    primer_full_count = defaultdict(int)
    primer_miss_count = defaultdict(int)

    for l in stdout:
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

        align = PrimerAlign(template, primer, difference, offset)
        if direction == b'+':
            matches_positive[template].append(align)
        else:
            matches_negative[template].append(align)
        
    print('primer statistics')
    for i,p in enumerate(seqstrs):
        print('primer %s:' % p)
        print(' full match = %4d' % primer_full_count[i])
        print(' miss match = %4d' % primer_miss_count[i])

    for v in matches_positive.values():
        v.sort(key=lambda x:x.offset)
    for v in matches_negative.values():
        v.sort(key=lambda x:x.offset)

    templates = set(matches_positive.keys()) & set(matches_negative.keys())

    def pcr_pair(pl, nl):
        if len(pl) < 1 or len(nl) < 1:
            return
        p = 0
        n = 0
        while p < len(pl):
            while not pl[p].offset <= nl[n].offset:
                n += 1
                if not (n < len(nl)):
                    return
            while p+1 < len(pl) and pl[p+1].offset <= nl[n].offset:
                p += 1
            if not threshold or nl[n].offset-pl[p].offset < threshold:
                yield pl[p], nl[n]
            p += 1

    result = []

    for t in templates:
        pl = matches_positive[t]
        nl = matches_negative[t]
        for p,n in pcr_pair(pl, nl):
            result.append(PCRFragment(t, p, n))

    return result

def main():
    import sys, os
    from optparse import OptionParser

    parser = OptionParser('usage: %prog [options] primer_sequence ...')
    parser.add_option("-t",'--template', dest='template', help="name of indexes of template", default=None)
    parser.add_option('-b', '--bisulfite', dest='bisulfite', action='store_true', default=False)
    parser.add_option('-l', '--length', dest='length', type=int, default=-1)

    (options, primers) = parser.parse_args()

    if len(primers) == 0:
        parser.error("no sequence input")

    if options.bisulfite and options.template:
        parser.error('bisulfite option and template option are exclusive.')

    length = options.length if options.length > 0 else None

    def s(template):
        '''
        -p 2 : number of thread
        --best : matches are ordered on number of mismatches.
        -k 1000 : maximum number of matches.
        -n 2 : allowed number of mismatches per primer.
        '''
        bowtie_path = os.path.expanduser('~/work/bowtie-1.0.0/bowtie -p 2 --best -k 1000 -n 3')
        #bowtie_path = os.path.expanduser('~/work/bowtie-0.12.7/bowtie -p 2 --best -k 1000 -n 2')

        result = virtualpcr(template, primers, bowtie_path, length)
        if result:
            for r in result:
                print('match at %s length=%10d %s(%s) -> %s(%s)' %(str(r.template).ljust(10), len(r), r.fw.origin,r.fw.mismatch, r.rv.origin,r.rv.mismatch))
        else:
            print('no result')
                          
    result = None
    if options.template:
        s(options.template)
    elif options.bisulfite:
        print('METHYL')
        s('bs_pos_met')
        s('bs_neg_met')
        print('UNMETHYL')
        s('bs_pos_unm')
        s('bs_neg_unm')
    else:
        parser.error('no template specified.')

if __name__=='__main__':
    main()
