import os

"""
subfs must have 2 methods
def write(self, filename, content_text)
def get_link_path(self, filename)
"""
class DefaultSubFileSystem(object):
    """just store"""
    def __init__(self):
        self.storage = []

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return filename

class SubFileSystem(object):
    def __init__(self, output_dir, prefix):
        self.output_dir = output_dir
        self.prefix = prefix
        self.storage = []

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finish()
        
    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return self._dirname() + filename

    def _dirname(self):
        return self.prefix + '_files/'

    def finish(self):
        d = os.path.join(self.output_dir,self._dirname())
        if not os.path.exists(d):
            os.makedirs(d)
        for name, content in self.storage:
            fn = os.path.join(self.output_dir, name)
            if isinstance(content, bytes):
                with open(fn, 'wb') as f:
                    f.write(content)
            elif isinstance(content, str):
                with open(fn, 'wb') as f:
                    f.write(content.encode('utf-8'))

    def get_subfs(self, name):
        return SubsubFileSystem(self, name)
        
class SubsubFileSystem(object):
    def __init__(self, parent, prefix):
        self.parent = parent
        self.prefix = prefix
        self.storage = self.parent.storage

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def write(self, filename, content_text):
        self.storage.append((self.get_link_path(filename), content_text))

    def get_link_path(self, filename):
        return self.parent.get_link_path(self.prefix + '_' + filename)

    def get_subfs(self, name):
        return SubsubFileSystem(self, name)
