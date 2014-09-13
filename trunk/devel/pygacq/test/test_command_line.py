import os
import datetime
import shutil

from nose import with_setup

# gacq specific stuff
import testutil
from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_number_as_filename():
    add_cursor_position(0.0, 0.0, "q")

    fullpath = testutil.get_data_file_name("N20130616S0089.fits")
    mdfdir = os.path.dirname(fullpath)
    newfname = "N" + datetime.datetime.utcnow().strftime("%Y%m%d") + "S0123.fits"
    shutil.copyfile(fullpath, newfname)

    ret = main(["gacq_test", "--mdfdir", mdfdir, "123"])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.219, tolerance=0.03)
    testutil.assert_tolerance(offsets[1],  0.102, tolerance=0.01)
    testutil.assert_tolerance(offsets[2],  0.010, tolerance=0.01)

