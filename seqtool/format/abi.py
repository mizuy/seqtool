import struct
import datetime
from contextlib import contextmanager
from ..nucleotide import base_color
from ..view import svg

base_index = {0:'G', 1:'A', 2:'T', 3:'C'}

@contextmanager
def report_exceptions():
    import traceback
    import pdb
    import sys
    try:
        yield
    except Exception:
        e, m, tb = sys.exc_info()
        print('exception traceback:'.ljust( 80, '=' ))
        for tbi in traceback.format_tb( tb ):
            print(tbi)
        print('  %s' % str( m ))
        print(''.rjust( 80, '=' ))
        pdb.post_mortem(tb)

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
    def __init__(self):
        pass

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
    
    def read(self, filename):
        with open(filename, 'rb') as fp:
            self.readfp(fp)

    def svg(self):
        ploc = self.get('PLOC').entries
        pbas = self.get('PBAS').entries
        pbas = ''.join(p.decode('utf-8') for p in pbas)
        print(pbas)
        t = svg.SvgItemsVStack()
        peak = svg.SvgBasePeaks(200, self.get_peaks(), ploc, pbas)
        raw = svg.SvgPeaks(200, self.get_raw())
        t.add(peak)
        t.add(raw)
        return svg.SvgPadding(20, 20, t).svg()

def main():
    a = AbiFormat()
    a.read('test.ab1')
    open('test.svg', 'w').write(a.svg())

if __name__ == '__main__':
    with report_exceptions():
        main()
