import sys
import testutil
from nose import with_setup

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io("4"), clear_cursor_positions)
def test_move_single_object_to_center_of_slit():
    add_cursor_position(498.0, 1146.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20060131S0014.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0], 1.994888, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], 0.929106, tolerance=0.08)

@with_setup(testutil.setup_io("4"), clear_cursor_positions)
def test_gmos_south_longslit():
    add_cursor_position(511.0, 1145.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0003.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0], -0.117, tolerance=0.06)
    testutil.assert_tolerance(offsets[1],  0.632, tolerance=0.06)

@with_setup(testutil.setup_io("4\nyes"), clear_cursor_positions)
def test_target_is_very_far_from_slit():
    add_cursor_position(526.0, 1186.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S00000000S0001.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0], -2.095, tolerance=0.02)
    testutil.assert_tolerance(offsets[1],  4.925, tolerance=0.005)
