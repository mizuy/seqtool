import struct
import datetime

from seqtool.util import svg, debug
from .render import SvgPeaks, SvgBasePeaks

base_index = {0:'G', 1:'A', 2:'T', 3:'C'}

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

        if self.elementtype == 18:
            #pascal string
            assert(self.elementsize == 1)
            assert(self.numelements == content[0] + 1)
            self.entries.append(str(content[1:]))
        elif self.elementtype == 19:
            #zero-end string
            assert(self.elementsize == 1)
            l =  self.numelements - 1
            content =  content[:l]
            s =  struct.unpack('{}s'.format(l), content)[0]
            self.entries.append(s)
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

class AbiFormat:
    def __init__(self, filename):
        with open(filename, 'rb') as fp:
            self.readfp(fp)

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

    def get_sequence(self):
        pbas = self.get('PBAS').entries
        pbas = ''.join(p.decode('utf-8') for p in pbas)
        return pbas

    def get_sequence_loc(self):
        return self.get('PLOC').entries
    
    def get_peaks(self):
        peaks = {}
        for i in range(4):
            base = base_index[i]
            values = self.get('DATA', 9 + i).entries
            peaks[base] = values
        return peaks

    def get_raw(self):
        peaks = {}
        for i in range(4):
            base = base_index[i]
            values = self.get('DATA', 1 + i).entries
            peaks[base] = values
        return peaks

    def svg(self):
        ploc = self.get_sequence_loc()
        pbas = self.get_sequence()
        t = svg.SvgItemsVStack()
        t.add(SvgBasePeaks(200, self.get_peaks(), ploc, pbas))
        t.add(SvgPeaks(200, self.get_raw()))
        return svg.SvgPadding(20, 20, t).svg()

    def svg_raw(self):
        return svg.SvgPadding(20, 20, SvgPeaks(200, self.get_raw())).svg()

def svg_raw(filename, output):
    with debug.report_exceptions():
        a = AbiFormat(filename)
        open(output, 'w').write(a.svg())
        #open(output, 'w').write(a.svg_raw())

def svg_peaks(filename, output):
    with debug.report_exceptions():
        a = AbiFormat(filename)
        open(output, 'w').write(a.svg())
        
if __name__ == '__main__':
    with debug.report_exceptions():
        a = AbiFormat()
        a.read('test.ab1')
        open('test.svg', 'w').write(a.svg())
