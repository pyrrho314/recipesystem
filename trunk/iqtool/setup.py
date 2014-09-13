#!/usr/bin/env python

"""
Setup script for iqtool for SOS-DA.

In this module:
    A mess.

Usage:
  python setup.py install --prefix=/astro/iraf/i686/gempylocal
  python setup.py sdist
"""

import os.path
from distutils.core import setup

MODULENAME = 'iqtool'

# PACKAGES and PACKAGE_DIR
SUBMODULES = ['iq', 
              'gemplotlib']
PACKAGES = [MODULENAME]
for m in SUBMODULES:
    PACKAGES.append('.'.join([MODULENAME,m]))
PACKAGE_DIRS = {}
PACKAGE_DIRS[MODULENAME] = '.'

#PACKAGE_DATA
PACKAGE_DATA = {}
PACKAGE_DATA[MODULENAME] = []
for s in ['.']+SUBMODULES:
    PACKAGE_DATA[MODULENAME].extend([os.path.join(s,'Copyright'),
                                     os.path.join(s,'ReleaseNote'),
                                     os.path.join(s,'README'),
                                     os.path.join(s,'INSTALL')
                                     ])

# DATA_DIRS and DATA_FILES
DATA_FILES = None

# SCRIPTS
IQTOOL_SCRIPTS = ['iqtool.py'
                  #os.path.join('lib','gemiq'),
#                  os.path.join('lib','pygem')
                  ]
# SOMEOTHER_SCRIPTS.extend(['another', 'script'])
SCRIPTS = []
SCRIPTS.extend(IQTOOL_SCRIPTS)

EXTENSIONS = None

setup ( name='iqtool',
        version='1.1',
        description='Image Quality tools',
        author='Gemini Observatory',
        author_email='jholt@gemini.edu',
        url='http://www.gemini.edu',
        maintainer='Kathleen Labrie',
        maintainer_email='klabrie@gemini.edu',
        packages=PACKAGES,
        package_dir=PACKAGE_DIRS,
        package_data=PACKAGE_DATA,
        data_files=DATA_FILES,
        scripts=SCRIPTS,
        ext_modules=EXTENSIONS,
        classifiers=[
            'Development Status :: Beta',
            'Intended Audience :: Gemini Instrument Scientists',
            'Operating System :: Linux :: RHEL',
            'Programming Language :: Python',
            'Topic :: Instrument Support',
            'Topic :: Instrument Checkouts',
            'Topic :: Engineering',
        ],
      )

