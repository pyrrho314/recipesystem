import sys
import testutil
from nose import with_setup

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

# each fiber is about .2 arcseconds, can't hope for any more accuracy then that

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_north_ifu_2():
    add_cursor_position(12.0, 97.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(39.0, 14.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20120122S0031.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0],  0.242, tolerance=0.1)
    testutil.assert_tolerance(offsets[1], -0.344, tolerance=0.1)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_north_ifu_R():
    add_cursor_position(10.0, 101.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(20.0, 12.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20111126S0354.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0],  0.142, tolerance=0.1)
    testutil.assert_tolerance(offsets[1],  0.169, tolerance=0.1)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_north_ifu_R_long_object():
    add_cursor_position(12.0, 98.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(19.0, 12.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130605S0153.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0],  0.281, tolerance=0.1)
    testutil.assert_tolerance(offsets[1], -0.161, tolerance=0.1)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_north_extended_fuzzy_object():
    add_cursor_position(8.6, 98.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(8.0, 14.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130716S0160.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.328, tolerance=0.1)
    testutil.assert_tolerance(offsets[1], -0.442, tolerance=0.1)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_ifu_R():
    add_cursor_position(9.0, 101.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(19.0, 16.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0002.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.235, tolerance=0.1)
    testutil.assert_tolerance(offsets[1], -0.920, tolerance=0.1)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_ifu_2_diffuse_object():
    add_cursor_position(10.5, 106.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(20.0, 10.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130421S0035.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], 0.615, tolerance=0.1)
    testutil.assert_tolerance(offsets[1], 0.094, tolerance=0.1)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_ifu_2_bright_star():
    add_cursor_position(10.0, 108.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(19.0, 13.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130227S0313.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.315, tolerance=0.1)
    testutil.assert_tolerance(offsets[1], -0.150, tolerance=0.1)
