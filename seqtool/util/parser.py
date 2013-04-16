__all__ = ['TreekvParser']

class Logger(object):
    def __init__(self):
        pass

    def write(self, message):
        print message

class KeyValue(object):
    def __init__(self, key, value=None, lineinfo=None, logger=Logger()):
        self.key = key
        self.value = value
        self._list = []
        self.lineinfo = lineinfo
        self.logger = logger
        self._flag = False
        self._parent = None
        self.path = '<unkown>/{}'.format(key)

    def value_list(self):
        return list(v.strip() for v in self.value.split(','))

    def __iter__(self):
        yield self.key
        yield self.value

    def set_parent(self, parent):
        self._parent = parent

        paths = []
        p = self
        while p:
            paths.append(p.key)
            p = p._parent
        self.path = '/'.join(paths[::-1])

    def items(self):
        return iter(self._list)

    def descendants(self):
        yield self
        for l in self._list:
            for ll in l.descendants():
                yield ll

    def get_path(self, path):
        if path=='.':
            return [self]
        elif '/' not in path:
            ret = self.get(path)
            return ret
        else:
            p,sep,q = path.partition('/')
            ret = []
            for rr in self.get(p):
                ret += rr.get_path(q)
            return ret

    def add(self, kv):
        kv.set_parent(self)
        self._list.append(kv)

    def has_item(self):
        return len(self._list)>0

    def last(self):
        try:
            return self._list[-1]
        except:
            raise IndexError()

    def get(self, key):
        ret = []
        for item in self._list:
            if item.key==key:
                ret.append(item)
        return ret

    def __repr__(self):
        return "KeyValue({},{},{})".format(self.key, self.value, self._list)

    def text(self, indent=0):
        ret = ''

        if not self._list:
            ret = "{}: {} {}\n".format(self.key, self.value, self.path)
        else:
            ret = "{}: {} {}\n".format(self.key, self.value, self.path)
            for l in self._list:
                ret += l.text(indent+4)

        return ''.join(' '*indent+l+'\n' for l in ret.splitlines())

    def set_used(self):
        self._flag = True

    def get_one(self, path):
        r = self.get_path(path)
        if not r:
            return None
        elif len(r)==1:
            r[0].set_used()
            return r[0]
        else:
            self.logger.write('Too many items {}'.format(path))
            return None

    def get_many(self,path):
        r = self.get_path(path)
        for rr in r:
            rr.set_used()
        return r

    def get_tree(self, path):
        r = self.get_path(path)
        for rr in r.descendants():
            rr.set_used()
        return r

    def get_unused(self):
        for rr in self.descendants():
            if not rr._flag:
                yield rr



class Lineinfo(object):
    def __init__(self, filename, line, lineno):
        self.filename = filename
        self.line = line
        self.lineno = lineno

    def error_msg(self, msg):
        return '{}\nin "{}", line {}: \n{}'.format(msg, self.filename, self.lineno, self.line)

def count_until(items, value):
    count = 0
    for i in items:
        if i==value:
            count += 1
        else:
            break
    return count

class TreekvParser(object):
    def __init__(self, tab_size=4, logger=Logger()):
        self.tab_size = tab_size
        self.logger = logger

    def read(self, filename):
        return self.readfp(open(filename,'rU'), filename)

    def readfp(self, fileobj, filename=None):
        if not filename:
            try:
                filename = fileobj.name
            except:
                filename = '<unkown>'

        lineno = 0
        
        root = KeyValue('root')
        tab_stops = [root]

        for line in fileobj:
            lineno += 1

            li = Lineinfo(filename, line, lineno)

            tab = count_until(line, ' ')

            if tab % self.tab_size != 0:
                self.logger.write(li.error_msg('Ignoring the line due to unkown tab stop {}. tab stops must be {}*n'.format(tab, self.tab_size)) )
                continue

            l = line.strip()

            if not l or l.startswith('#') or l.startswith('//') or l.startswith(';'):
                continue

            if not l:
                continue
            if ':' not in l:
                self.logger.write(li.error_msg('Unkown line. line format must be "key:value"'))
                continue
            key,sep,value = l.partition(':')
            item = KeyValue(key.strip(), value.strip(), li, self.logger)

            level = tab / self.tab_size
            current_level = len(tab_stops) - 1
            current_parent = tab_stops[-1]

            if level==current_level:
                current_parent.add(item)
            elif level == current_level+1:
                assert(current_parent.has_item())
                new_parent = current_parent.last()
                new_parent.add(item)
                tab_stops.append(new_parent)
            elif level > current_level:
                self.logger.write(li.error_msg('Too many indent spaces. This indent must be less than {}, but {}'.format(self.tab_size*(level+1), self.tab_size*level) ))
                continue
            elif level < current_level:
                tab_stops = tab_stops[:level+1]
                parent = tab_stops[-1]
                parent.add(item)

        return root

sample = """
general:
    gene: NDRG2
    primers: primers.txt
    bsa_data: bsa_data.txt

tss: TSS
    tissues: Brian, Liver, Colon

motifs: motifs
    p53BS: GTGCAAGGTCCGGGGCGCTTGGCA
    TATAbox: TATAWAW
    mir650: TGCCTCC
    BamHI: GGATCC
    XhoI: CTCGAG
    ecorv: GATATC
    ecori: GAATTC
    WT1: GTGTGTGTGTGTG
    HRE3: GCGTG
    HRE2: GCGTG
    HRE1: GCGTCC
    probe: CGGGCGGCTGGACGCTTCCAGGCTCTGCTCGGCTCACCAAAACATTCCAC

pcr: Genomic-PCR
    ChIP1: ChIP1-FW, ChIP1-RV
    ChIP2: ChIP2-FW, ChIP2-RV
    ChIP2-dash: BSP4-FW, ChIP2-RV
    ChIP3: ChIP3-FW, ChIP3-RV
    upstream: genome-up-stream-FW, NDRG2 cDNA 1ab-3 RV

"""

if __name__=='__main__':
    import StringIO
    parser = TreekvParser()
    kv = parser.readfp(StringIO.StringIO(sample))
    print kv.text()
