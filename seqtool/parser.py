

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
