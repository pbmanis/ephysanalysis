from setuptools import setup, find_packages
import os

# Use Semantic Versioning, http://semver.org/
version_info = (0, 1, 9, '')
__version__ = '%d.%d.%d%s' % version_info


setup(name='ephysanalysis',
      version=__version__,
      description='Methods for analysis of elecrophysiology data: IV and FI curves',
      url='http://github.com/pbmanis/cnmodel',
      author='Paul B. Manis',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['ephysanalysis*']),
      zip_safe=False,
      entry_points={
          'console_scripts': [
               'dataSummary=ephysanalysis.dataSummary:main',
               'ma2tiff=ephysanalysis.ma2tiff:convertfiles',
               'bridge=ephysanalysis.bridge:main',
          ]
      },
      classifiers = [
             "Programming Language :: Python :: 3.6+",
             "Development Status ::  Beta",
             "Environment :: Console",
             "Intended Audience :: Neuroscientists",
             "License :: MIT",
             "Operating System :: OS Independent",
             "Topic :: Scientific Software :: Tools :: Python Modules",
             "Topic :: Data Processing :: Neuroscience",
             ],
    )
      