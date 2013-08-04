from seqtool.util import svg
from seqtool.util.rectangle import Rectangle
from seqtool.nucleotide import base_color
from seqtool import format

fm = svg.fm



class SvgBaseText(svg.SvgBase):
    def __init__(self, ploc, pbas, scalex = 1., fontsize = svg.DEFAULT_FONTSIZE):
        fw = svg.font_width(fontsize)
        self.x = ' '.join(fm(x * scalex - fw/2.) for x in ploc)
        self.pbas = pbas
        self.y = 0
        self.w = max(ploc) * scalex
        self.h = svg.font_height(fontsize)

        self.style = {}
        if fontsize != svg.DEFAULT_FONTSIZE:
            self.style['fontsize'] = fontsize

        self._rect = Rectangle(0, self.w, 0, self.h)

    def draw(self,b):
        with b['text'](x=self.x, y=fm(self.h), **self.style):
            b.text(self.pbas)

class SvgBasePeaks(svg.SvgItemsVStack):
    def __init__(self, height, peaks, ploc, pbas, scalex = 1.):
        super().__init__()
        self.add(SvgBaseText(ploc, pbas, scalex))
        self.add(SvgPeaks(height, peaks, scalex))


class SvgPeaks(svg.SvgItemsFixedHeight):
    def __init__(self, height, peaks, scalex = 1.):
        super().__init__(height)

        self.scalex = scalex
        self.peaks = peaks

        self.max_peak = 0
        self.length = 0
        for base, values in peaks.items():
            self.max_peak = max(self.max_peak, max(values))
            self.length = max(self.length, len(values))

        self.height = height
        self.scaley = 1. * self.height/self.max_peak

        self.width = self.scalex * self.length
        self._rect = Rectangle(0, self.width, 0, self.height)

    def draw(self, b):
        for base, values in self.peaks.items():
            color = base_color(base)
            style = 'stroke:{}'.format(color)

            v = ['{:.2f} {:.2f}'.format(x * self.scalex, self.height - y * self.scaley) for x, y in enumerate(values)]
            
            b.polyline(points = ','.join(v), style = style)

    @property
    def rect(self):
        return self._rect

