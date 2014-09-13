#!/usr/bin/env python

import os
import sys

if os.path.exists('MANIFEST'): os.remove('MANIFEST')


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)

    # add gspline module
    #config.add_extension('gspline', ['gspline/gspline.pyf','gspline/gspline.c'])

    config.add_data_files('*.cl')
    config.add_data_files('*.par')
    config.add_data_files('bad_r.fits','bad_b.fits')

    #config.add_subpackage('astrodata')
    #config.add_data_dir('astrodata')
    config.add_data_dir('doc')
    config.add_data_dir('doc/build/html')
    config.add_data_dir('test')
    config.add_data_files('README')

    return config

def setup_package():

    from numpy.distutils.core import setup

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0,local_path)

    try:
        setup(
            name = 'nici',
            maintainer = "Nelson Zarate",
            maintainer_email = "nzarate@gemini.edu",
            url = "",
            download_url = "",
            platforms = ["Linux", "Mac OS-X", "Unix"],
            configuration=configuration )
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup_package()
