import os,sys
#import ez_setup
#ez_setup.use_setuptools()

from setuptools import setup, find_packages

import numpy
from Cython.Build import cythonize
# not compatible with distribute
#from distutils.extension import Extension
#from Cython.Distutils import build_ext

version = '0.6.0'

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

entry_points = """
[console_scripts]
seqview = seqtool.frontend.command:seqview
get_genbank = seqtool.frontend.command:get_genbank
geneview = seqtool.frontend.command:geneview
tssview = seqtool.frontend.command:tssview
primers = seqtool.frontend.command:primers
sequencing = seqtool.frontend.command:sequencing
seqtooldb = seqtool.db:seqdb_command
abiview = seqtool.frontend.command:abiview

convert_bs = seqtool.bowtie.convert_bs:main
virtualpcr = seqtool.bowtie.virtualpcr:main

primer = seqtool.nucleotide.primer:main
probe = seqtool.nucleotide.primer:probe

pdesign = seqtool.script.pdesign:main
bisulfite = seqtool.script.bisulfite:main
rpm = seqtool.script.cf:main
translate = seqtool.script.translate:main
server = seqtool.server.server:main
"""
setup(
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
    install_requires=['biopython','numpy', 'sqlalchemy', 'cython', 'appdirs', 'mysql-connector-python'],
    test_suite='nose.collector',
    ext_modules = cythonize('seqtool/seqtool/nucleotide/sw_c.pyx'),
    include_dirs = [numpy.get_include()],
    entry_points=entry_points
)
