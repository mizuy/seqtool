import struct
import datetime

from seqtool.util import svg, debug
from .render import SvgPeaks, SvgBasePeaks
from ..nucleotide import to_seq

base_index = {0:'G', 1:'A', 2:'T', 3:'C'}
rc = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}

"""
see ABIF_File_Format.pdf
Applied Biosystems 3500/3500xl Genetic Analyzer tags
"""

class InvalidFormat(Exception):
    pass

# elementtype to unpack str
ET_US =  {
    3 : 'H', # unsigned short
    4 : 'h', # short
    5 : 'L', # unsigned int
    6 : 'l', # int
    7 : 'f', # float
    8 : 'd', # double
}

class Entry:
    def __init__(self, buffer, filep):
        self.filep = filep
        assert(len(buffer) == 28)

        self.name =  buffer[:4].decode('utf-8')
        self.number, self.elementtype, self.elementsize, self.numelements, self.datasize = struct.unpack('>lhhll', buffer[4:20])
        
        assert(0 < self.elementtype)
        assert(0 < self.elementsize)
        assert(0 < self.numelements)
        assert(0 < self.datasize)
        assert(self.numelements * self.elementsize <= self.datasize)

        if self.datasize > 4:
            self.dataoffset = struct.unpack('>l', buffer[20:24])[0]
            assert(0 < self.dataoffset)
            self.filep.seek(self.dataoffset)
            content = self.filep.read(self.datasize)
        else:
            self.dataoffset = None
            content =  buffer[20:24]

        self.entries = []
        self.string = None

        if self.elementtype == 18:
            #pascal string
            assert(self.elementsize == 1)
            assert(self.numelements == content[0] + 1)
            s = content[1:].decode('utf-8')
            self.entries.append(s)
            self.string = s
        elif self.elementtype == 19:
            #zero-end string
            assert(self.elementsize == 1)
            l =  self.numelements - 1
            content =  content[:l]
            s =  struct.unpack('{}s'.format(l), content)[0]
            self.entries.append(s)
            self.string = s
        else:
            i = 0
            for t in range(self.numelements):
                self.entries.append(self.read_entry(content[i:i + self.elementsize]))
                i += self.elementsize

    def __repr__(self):
        return 'Entry({},number:{},num:{})'.format(self.name, self.number, self.numelements)

    def read_entry(self, contents):
        if self.elementtype == 1023:
            return Entry(contents, self.filep)
        elif self.elementtype in ET_US.keys():
            up = '>' + ET_US[self.elementtype]
            assert(struct.calcsize(up) == self.elementsize)
            return struct.unpack(up, contents)[0]
        elif self.elementtype == 1:
            return bytes(contents)
        elif self.elementtype == 2:
            return bytes(contents)
        elif self.elementtype == 10:
            year, month, day = struct.unpack('>hBB', contents)
            return datetime.date(year, month, day)
        elif self.elementtype == 11:
            hour, minute, second, hsecond = struct.unpack('>BBBB', contents)
            return datetime.time(hour, minute, second, hsecond)
        elif self.elementtype == 18:
            # pascal string
            pass
        elif self.elementtype == 19:
            # zero-end string
            pass
        else:
            return contents

class AbiPeakViewBase:
    def get_sequence(self):
        return self.pbas
    def get_location(self):
        return self.ploc
    def get_peaks(self):
        return self.peaks
    def get_raw_peaks(self):
        return self.raw

    def get_svg_base_peaks(self, height, scalex=1.):
        return SvgBasePeaks(height, self.peaks, self.ploc, self.pbas, scalex)

    def get_svg_peaks(self, height, scalex=1.):
        return SvgPeaks(height, self.peaks, scalex, self.ploc)

    def get_svg_raw_peaks(self, height, scalex=1.):
        return SvgPeaks(height, self.raw, scalex)
    
class AbiPeakView(AbiPeakViewBase):
    def __init__(self, abifile):
        self.abifile = abifile
        self.pbas = abifile.get_pbas()
        self.ploc = abifile.get_ploc()

        self.peaks = {}
        self.peaks_len = 0
        self.raw = {}
        self.raw_len = 0
        for i in range(4):
            base = base_index[i]
            self.peaks[base] = abifile.get_peaks(i)
            self.peaks_len = max(self.peaks_len, len(self.peaks[base]))
            self.raw[base] = abifile.get_raw_peaks(i)
            self.raw_len = max(self.raw_len, len(self.raw[base]))

    def get_reverse_complement(self):
        return AbiPeakRcView(self)

class AbiPeakRcView(AbiPeakViewBase):
    def __init__(self, original):
        self.original = original
        self.pbas = str(to_seq(self.original.pbas).reverse_complement())
        self.ploc = [(self.original.peaks_len - 1) - i for i in self.original.ploc][::-1]
        assert(len(self.pbas)==len(self.ploc))
        self.peaks = {}
        self.peaks_len = self.original.peaks_len
        self.raw = {}
        self.raw_len = self.original.raw_len

        for p,q in self.original.peaks.items():
            self.peaks[rc[p]] = q[::-1]
        for p,q in self.original.raw.items():
            self.raw[rc[p]] = q[::-1]

    def get_reverse_complement(self):
        return self.original

        
class AbiFormat:
    def __init__(self, filename):
        with open(filename, 'rb') as fp:
            self.readfp(fp)
        self.view = AbiPeakView(self)
        self.rc_view = self.view.get_reverse_complement()

    def get_view(self, sense):
        return self.view if sense else self.rc_view

    def readfp(self, fp):
        header =  fp.read(4)
        if header != b'ABIF':
            raise InvalidFormat()

        self.version =  struct.unpack('h', fp.read(2))[0]
        
        self.head = Entry(fp.read(28), fp)
        self._items = {}
        for i in self.head.entries:
            self._items[(i.name, i.number)] = i

    def get(self, name, number = 1):
        return self._items.get((name, number))

    def get_name(self):
        return self.get('SMPL').string

    def get_pbas(self):
        return ''.join(p.decode('utf-8') for p in self.get('PBAS').entries)
    def get_ploc(self):
        return self.get('PLOC').entries
    def get_peaks(self, i):
        return self.get('DATA', 9 + i).entries
    def get_raw_peaks(self, i):
        return self.get('DATA', 1 + i).entries

    def svg(self):
        t = svg.SvgItemsVStack()
        t.add(svg.SvgText(self.get_name()))
        t.add(svg.SvgText('Original Sequencing Result:'))
        t.add(self.view.get_svg_base_peaks(200))
        t.add(svg.SvgText('Reversed-Complemental Result:'))
        t.add(self.rc_view.get_svg_base_peaks(200))

        t.add(svg.SvgItemsFixedHeight(8))
        t.add(svg.SvgText('Raw Data:'))
        t.add(self.view.get_svg_raw_peaks(500))
        return svg.SvgPadding(20, 20, t).svg()

def write_svg(filename, output):
    with debug.report_exceptions():
        open(output, 'w').write(AbiFormat(filename).svg())

if __name__ == '__main__':
    with debug.report_exceptions():
        a = AbiFormat('files/test.ab1')
        open('test.svg', 'w').write(a.svg())
