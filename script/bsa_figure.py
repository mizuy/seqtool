
import seqtool.util.svg as svg
from seqtool.util.parser import TreekvParser
import sys
import math

CSIZE = 10
PADDING = 4
HEIGHT = CSIZE+3

MAX_KEY_LENGTH = 5

def is_unmethyl(base):
    return base=='T' or base=='A'
def is_methyl(base):
    return base=='G' or base=='C'

def get_color(base):
    if is_methyl(base):
        return 'black'
    elif is_unmethyl(base):
        return 'white'
    else:
        return 'gray'

def n_mean_sd(i):
    i = list(i)
    n = len(i)
    mean = sum(i)/n
    sd = math.sqrt(sum(1.*(ii-mean)**2 for ii in i)/n)
    return n,mean,sd


def methyl_ratio(seq):
    m = sum(1 if is_methyl(s) else 0 for s in seq)
    return 1.*m/len(seq)

def main():
    filename = sys.argv[1]

    parser = TreekvParser()

    figure = svg.SvgItemsVStack()

    with open(filename,'r') as fileobj:
        tree = parser.readfp(fileobj, filename)

        for kv in tree.items():
            ff = svg.SvgItemsVStack()

            ll = [(methyl_ratio(c.value), c.key, c.value) for c in kv.items()]
            ll.sort(reverse=True)

            for ratio, key, value in ll:
                bs = value
                gg = svg.SvgItemsFixedHeight(HEIGHT)

                x = 20
                gg.add(svg.SvgText(key, x+(MAX_KEY_LENGTH-len(key))*svg.font_width(), 0))
                x += svg.font_width()*MAX_KEY_LENGTH + 10
                gg.add(svg.SvgHline(x, x+(len(bs)+1)*(CSIZE+PADDING), HEIGHT/2.))
                for site in bs:
                    x += CSIZE+PADDING
                    gg.add(svg.SvgCircle(x, HEIGHT/2., CSIZE/2., fill=get_color(site), stroke='black'))

                ff.add(gg)
            figure.add(svg.SvgText(kv.key, 0, 0))
            figure.add(ff)
            n,mean,sd = n_mean_sd(c[0] for c in ll)
            figure.add(svg.SvgText("n={}, mean={:.2f}, sd={:.2f}".format(n,mean,sd), 20, 0))

    print(svg.SvgPadding(10,10,figure).svg())

if __name__=='__main__':
    main()
