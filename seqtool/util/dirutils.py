
import os
import inspect

def module_path(local_function):
   ''' returns the module path without the use of __file__.  Requires a function defined 
   locally in the module.
   from http://stackoverflow.com/questions/729583/getting-file-path-of-imported-module'''
   return os.path.abspath(inspect.getsourcefile(local_function))

class Filepath(object):
    def __init__(self, filename):
        self.path = os.path.abspath(filename)
        self.dir = os.path.dirname(self.path)
        self.basename = os.path.basename(self.path)
        self.prefix = self.basename.rpartition('.')[0]

    def change_ext(self, ext):
        return Filepath(os.path.join(self.dir, self.prefix + ext))

    def relative(self, name):
        return os.path.join(self.dir,name)

    def __str__(self):
        return self.path

class LocalFs(object):
    def __init__(self, filename):
        self.p = Filepath(filename)

    def read(self, filename):
        return open(self.p.relative(filename), 'r')

    def write(self, filename, context):
        return open(self.p.relative(filename), 'w')
