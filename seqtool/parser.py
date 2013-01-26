
def parse_file(fileobj):
    lineno = 0
    for l in fileobj:
        lineno += 1
        l = l.strip()
        if not l or l.startswith('#') or l.startswith('//'):
            continue
        def error_message(msg):
            print ':%s: %s: "%s"'%(lineno,msg,l)
        ls = l.split(':')
        if len(ls)!=2:
            error_message('unknown line')
            continue
        name = ls[0].strip()
        value = ls[1].strip()
        yield name,value,error_message

class Block(object):
    def __init__(self, name, line, lineno=None):
        self.name = name
        self.line = line
        self.lineno = lineno
        self.list = []
    def add(self, line, lineno=None):
        self.list.append((line,lineno))
    def __iter__(self):
        return iter(self.list)

class SettingFile(object):
    def __init__(self):
        self.blocks = []

    def add_block(self,name,line,lineno):
        self.blocks.append(Block(name,line,lineno))

    def add_line(self,line,lineno):
        if not self.blocks:
            self.blocks.append(Block('No Name',line,lineno))
        self.blocks[-1].add(line,lineno)

    def __iter__(self):
        return iter(self.blocks)

    def parse(self, fileobj):
        lineno = 0
        skipping = False
        for l in fileobj:
            lineno += 1
            l = l.strip()

            if not l or l.startswith('#') or l.startswith('//'):
                continue

            if l.startswith('/+'):
                skipping = True
            elif l.startswith('+/'):
                skipping = False
            else:
                if skipping:
                    continue
                if l.startswith('>'):
                    self.add_block(l[1:].strip(), l, lineno)
                else:
                    self.add_line(l, lineno)
