import os
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

setup(name='seqtool',
    version=version,
    description=("small scripts visualizing PCR products for molecular biology experimetnts."),
    classifiers = filter(None, classifiers.split("\n")),
    keywords='pcr biology bisulfite',
    author='mizuy',
    author_email='mizugy@gmail.com',
    url='http://github.com/mizuy/seqtool',
    license='MIT',
    packages=find_packages(),
    install_requires=['biopython','numpy', 'sqlalchemy', 'cython'],
    test_suite='nose.collector',
#    test_requires=['Nose'],
    ext_modules = cythonize('seqtool/nucleotide/sw_c.pyx'),
    include_dirs = [numpy.get_include()],
    entry_points=\
"""
[console_scripts]
seqview = seqtool.command:seqview
get_genbank = seqtool.command:get_genbank
geneview = seqtool.command:geneview
tssview = seqtool.command:tssview
primers = seqtool.command:primers
seqvcmd = seqtool.seqvcmd:main
convert_bs = seqtool.bowtie.convert_bs:main
virtualpcr = seqtool.bowtie.virtualpcr:main
seqtool_loadsql = seqtool.db.sql:load_all
bisulfite = seqtool.script.bisulfite:main
rpm = seqtool.script.cf:main
align = seqtool.script.align:main
sequencing = seqtool.command:sequencing
""")
