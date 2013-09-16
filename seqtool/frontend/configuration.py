from appdirs import AppDirs
import configparser, os

import readline
import pkg_resources

DEFAULT = """
[general]
"""

class GeneralConfiguration(object):
    def __init__(self):
        version = pkg_resources.require("seqtool")[0].version
        dirs = AppDirs("seqtool", "mizugy", version)
        self.configfile = os.path.join(dirs.user_data_dir,'config.cfg')
        self.cache_dir = dirs.user_cache_dir
        print('configfile: ',self.configfile)
        print('cache_dir: ',self.cache_dir)

        if not os.path.exists(self.configfile):
            os.makedirs(dirs.user_data_dir)
            with open(self.configfile, 'w') as f:
                print('creating...',self.configfile)
                f.write(DEFAULT)

        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

    def get_config(self):
        config = configparser.SafeConfigParser()
        config.read(self.configfile)
        return config

    def get_cache_dir(self):
        return self.cache_dir

    def get_email(self):
        config = self.get_config()
        try:
            email = config.get('general','email')
        except:
            print("Email address is required to access NCBI Entrez.")
            email = input("Email: ")
            if not config.has_section('general'):
                config.add_section('general')
            config.set('general','email',email)
            config.write(open(self.configfile,'w'))
            print("Your email is stored at {}".format(self.configfile))
        return email

