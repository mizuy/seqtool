# -*- coding:utf-8 mode:Python -*-
import cStringIO
import string, os, re, difflib, itertools

# http://wiki.python.org/moin/PythonDecoratorLibrary#Memoize
import functools
class memoize(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            value = self.func(*args)
            self.cache[args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def method_call(self, self_, *args):
        try:
            cache = getattr(self_,self.name)
        except:
            self.name = '_cache' + self.func.__name__
            cache = {}
            setattr(self_,self.name,cache)

        try:
            return cache[args]
        except KeyError:
            value = self.func(self_,*args)
            cache[args] = value
            return value
        except TypeError:
            return self.func(self_,*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.method_call, obj)

def curry(func):
    def inner(*args,**kw):
        def inner2(w):
            func(w,*args,**kw)
        return inner2
    return inner

def curry2(func):
    def inner(*args,**kw):
        def inner2(w,ri):
            func(w,ri,*args,**kw)
        return inner2
    return inner
