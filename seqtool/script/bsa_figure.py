
from ..util import svg
from ..util.parser import TreekvParser
import sys, os
import math
from collections import OrderedDict
CSIZE = 10
PADDING = 4
HEIGHT = CSIZE+3

def get_color(base):
    if base in 'GC':
        return 'black'
    elif base in 'TA':
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
    m = sum(1 if (s in 'GC') else 0 for s in seq)
    return 1.*m/len(seq)


class BsaFigure:
    def __init__(self, filename=None):
        if filename:
            self.read(filename)
            
    def read(self, filename):
        with open(filename,'r') as fileobj:
            self.readfp(fileobj, filename)

    def readfp(self, fileobj, filename=None):
        parser = TreekvParser()
        tree = parser.readfp(fileobj, filename)
    
        kv = tree.get_one('settings/template')
        if kv:
            self.template = kv.value
            self.sites = [i for i,c in enumerate(self.template) if c in 'YR' ]
            r = tree.get_one('results')
        else:
            self.template = None
            r = tree
            
        self.results = OrderedDict()
        for kv in r.items():
            self.results[kv.key] = [(c.key, c.value) for c in kv.items()]

    def sort(self):
        ret = OrderedDict()

        keys = self.results.keys()
        for key in keys:
            values = self.results[key]
            ret[key] = sorted(values, key = lambda t: methyl_ratio(self.get_condensed(t[1])), reverse=True)

        self.results = ret

    def flatten(self):
        ret = []
        for clone,values in self.results.items():
            ret.extend(values)
        return ret
        
    def max_name_length(self):
        return max(len(name) for name,r in self.flatten())

    def get_condensed(self, result):
        if self.template:
            return ''.join([result[i] for i in self.sites])
        else:
            return result

    def svg_bar_circles(self, result, x, fixed):
        dx = CSIZE+PADDING
        if self.template:
            if not fixed:
                length_original = len(self.template)
                length = len(self.sites) * 2
                ratio = length / length_original
                i_result = [(i*ratio, result[i]) for i in self.sites]
            else:
                length = len(self.sites)
                i_result = enumerate(self.get_condensed(result))
        else:
            length = len(result)
            i_result = enumerate(self.get_condensed(result))

        gg = svg.SvgItems()
        gg.add(svg.SvgHline(x, x+(length+1)*dx, HEIGHT/2.))
        for i, site in i_result:
            if get_color(site)!='gray':
                gg.add(svg.SvgCircle(x+i*dx, HEIGHT/2., CSIZE/2., fill=get_color(site), stroke='black'))
        return gg

    def svg(self, fixed, skip_trash=True, show_names=False):
        max_name_length = self.max_name_length()

        figure = svg.SvgItemsVStack()
        for clone, values in self.results.items():
            if skip_trash and os.path.basename(clone)=='trash':
                continue
                
            figure.add(svg.SvgText(clone, 0, 0))
            figure.add(svg.SvgItemsFixedHeight(HEIGHT/2.))

            ff = svg.SvgItemsVStack()
            for name, result in values:
                gg = svg.SvgItemsFixedHeight(HEIGHT)

                x = 20
                
                if show_names:
                    #gg.add(svg.SvgText(name.rjust(max_name_length), x, 0))
                    gg.add(svg.SvgText(name, x+(max_name_length-len(name))*svg.font_width(), 0))
                    x += svg.font_width()*max_name_length + 10
                
                gg.add(self.svg_bar_circles(result, x, fixed))

                ff.add(gg)
            ff.add(svg.SvgItemsFixedHeight(HEIGHT))
            figure.add(ff)
        return svg.SvgPadding(10,10,figure).svg()
            
def main():
    bf = BsaFigure(sys.argv[1])
    bf.sort()
    print(bf.svg(False))

if __name__=='__main__':
    main()
