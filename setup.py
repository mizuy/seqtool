import os
from setuptools import setup, find_packages

version = '0.1.0'

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
      install_requires=['biopython','numpy'],
      test_suite='nose.collector',
      test_requires=['Nose'],
      entry_points=\
"""
[console_scripts]
seqview = seqtool.seqview:main
primers = seqtool.primers:main
bisearch = seqtool.bisearch:main
""")