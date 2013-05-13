
import seqtool.view.svg as svg
from seqtool.util.parser import TreekvParser
import sys

CSIZE = 10
PADDING = 4
HEIGHT = CSIZE+3

def main():
    filename = sys.argv[1]

    parser = TreekvParser()

    figure = svg.SvgItemsVStack()

    with open(filename,'r') as fileobj:
        tree = parser.readfp(fileobj, filename)

        for kv in tree.items():
            ff = svg.SvgItemsVStack()
            for clone in kv.items():
                gg = svg.SvgItemsFixedHeight(HEIGHT)
                bs = clone.value
                gg.add(svg.SvgHline(0, (len(bs)+1)*(CSIZE+PADDING), HEIGHT/2.))
                x = 0
                for site in bs:
                    x += CSIZE+PADDING
                    if site=='T' or site=='A':
                        gg.add(svg.SvgCircle(x, HEIGHT/2., CSIZE/2., fill='white'))
                    elif site=='G' or site=='C':
                        gg.add(svg.SvgCircle(x, HEIGHT/2., CSIZE/2., fill='black'))
                    else:
                        gg.add(svg.SvgCircle(x, HEIGHT/2., CSIZE/2., fill='gray'))
                ff.add(gg)
            figure.add(svg.SvgText(kv.key, 0, 0))
            figure.add(ff)

    print svg.SvgPadding(10,10,figure).svg()

if __name__=='__main__':
    main()
