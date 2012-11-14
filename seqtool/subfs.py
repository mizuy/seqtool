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
    def __init__(self, output_dir, basename):
        self.output_dir = output_dir
        self.basename = basename
        self.storage = []

    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return self.basename + '_' + filename

    def finish(self):
        for name, content in self.storage:
            with open(os.path.join(self.output_dir, name), 'w') as f:
                f.write(content)

    def get_subfs(self, name):
        return SubsubFileSystem(self, name)
        
class SubsubFileSystem(object):
    def __init__(self, parent, basename):
        self.parent = parent
        self.basename = self.parent.basename + '_' + basename
        self.storage = self.parent.storage

    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return self.basename + '_' + filename

    def get_subfs(self, name):
        return SubsubFileSystem(self, name)
