import os,sys
#import ez_setup
#ez_setup.use_setuptools()

from setuptools import setup, find_packages

from Cython.Build import cythonize

# not compatible with distribute
#from distutils.extension import Extension
#from Cython.Distutils import build_ext

import numpy
version = '0.2.1'

README = os.path.join(os.path.dirname(__file__), 'README')
long_description = open(README).read() + '\n\n'

classifiers = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Operating System :: Unix
"""

mainscript = 'seqtool/Application.py'
entry_points = """
[console_scripts]
seqview = seqtool.command:seqview
get_genbank = seqtool.command:get_genbank
geneview = seqtool.command:geneview
tssview = seqtool.command:tssview
primers = seqtool.command:primers
seqvcmd = seqtool.seqvcmd:main
convert_bs = seqtool.bowtie.convert_bs:main
virtualpcr = seqtool.bowtie.virtualpcr:main
bisulfite = seqtool.script.bisulfite:main
rpm = seqtool.script.cf:main
align = seqtool.script.align:main
translate = seqtool.script.translate:main
sequencing = seqtool.command:sequencing
server = seqtool.server.server:main
database_load = seqtool.db.sql:database_load
"""

if sys.platform == 'darwin':
    extra_options = dict(
        setup_requires=['py2app'],
        app=[mainscript],
        # Cross-platform applications generally expect sys.argv to
        # be used for opening files.
        options={'py2app': {'argv_emulation': True}},
        entry_points=entry_points,
    )
elif sys.platform == 'win32':
    extra_options = dict(
        setup_requires=['py2exe'],
        app=[mainscript],
    )
else:
    extra_options = dict(
        # Normally unix-like platforms will use "setup.py install"
        # and install the main script as such
        entry_points=entry_points,
    )

setup(
    data_files=['_db'],

    name='seqtool',
    version=version,
    description=("small scripts visualizing PCR products for molecular biology experimetnts."),
    classifiers = filter(None, classifiers.split("\n")),
    keywords='pcr biology bisulfite',
    author='mizuy',
    author_email='mizugy@gmail.com',
    url='http://github.com/mizuy/seqtool',
    license='MIT',
    packages=find_packages(),
    install_requires=['biopython','numpy', 'sqlalchemy', 'cython', 'bottle'],
    test_suite='nose.collector',
#    test_requires=['Nose'],
    ext_modules = cythonize('seqtool/nucleotide/sw_c.pyx'),
    include_dirs = [numpy.get_include()],
    **extra_options)
