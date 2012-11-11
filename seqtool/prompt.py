__all__ = ['prompt']

import sys

class prompt(object):
    """How to use
    with prompt('start something!', 'complete!!') as p:
        for something:
            do something
            p.progress()
    """
    def __init__(self, start_msg='start', end_msg='finished.', progress_msg='.'):
        self.s = start_msg
        self.p = progress_msg
        self.e = end_msg

    def __enter__(self):
        print self.s,
        return self

    def __exit__(self, type, value, traceback):
        print self.e
        sys.stdout.flush()

    def progress(self, msg=None):
        if msg:
            sys.stdout.write(msg)
        else:
            sys.stdout.write(self.p)
        sys.stdout.flush()
