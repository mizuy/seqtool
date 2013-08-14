from contextlib import contextmanager

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
