from __future__ import absolute_import

from cStringIO import StringIO
import os

class DefaultSubFileSystem(object):
    """just store"""
    def __init__(self):
        self.storage = []

    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return filename

class SubFileSystem(object):
    def __init__(self, output_dir, suffix):
        self.output_dir = output_dir
        self.suffix = suffix
        self.storage = []

    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return self._dirname() + filename

    def _dirname(self):
        return self.suffix + '_files/'

    def finish(self):
        d = os.path.join(self.output_dir,self._dirname())
        if not os.path.exists(d):
            os.makedirs(d)
        for name, content in self.storage:
            fn = os.path.join(self.output_dir, name)
            with open(fn, 'w') as f:
                f.write(content)

    def get_subfs(self, name):
        return SubsubFileSystem(self, name)
        
class SubsubFileSystem(object):
    def __init__(self, parent, suffix):
        self.parent = parent
        self.suffix = suffix
        self.storage = self.parent.storage

    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return self.parent.get_link_path(self.suffix + '_' + filename)

    def get_subfs(self, name):
        return SubsubFileSystem(self, name)
