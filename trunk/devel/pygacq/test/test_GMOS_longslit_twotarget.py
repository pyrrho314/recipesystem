import sys
import testutil
from nose import with_setup

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io("4"), clear_cursor_positions)
def test_custom_center_with_objests_really_far_apart():
    add_cursor_position(513.0, 1150.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(514.0, 1944.0, "a")
    add_cursor_position(513.0, 1547.0, "c")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N00000000S0002.fits")])

    offsets = testutil.get_offsets()

    # high tolerance due to 1 pixel difference in x dimension in fit of 2nd target
    testutil.assert_tolerance(offsets[0],  -0.079, tolerance=0.1)
    testutil.assert_tolerance(offsets[1],   0.007, tolerance=0.1)
    testutil.assert_tolerance(offsets[2],   0.105, tolerance=0.94)

@with_setup(testutil.setup_io("4\nyes"), clear_cursor_positions)
def test_custom_center():
    add_cursor_position(508.0, 1148.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(522.0, 1099.0, "a")
    add_cursor_position(512.0, 1152.0, "c")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N20060131S0020.fits")])

    offsets = testutil.get_offsets()

    # high tolerance due to 1 pixel difference in x dimension in fit of 2nd target
    testutil.assert_tolerance(offsets[0],  -0.584, tolerance=0.1)
    testutil.assert_tolerance(offsets[1],   5.745, tolerance=0.1)
    testutil.assert_tolerance(offsets[2], -19.375, tolerance=0.94)

@with_setup(testutil.setup_io("4\nyes\nyes"), clear_cursor_positions)
def test_move_far_to_put_midpoint_at_center():
    add_cursor_position(508.0, 1148.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(522.0, 1099.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N20060131S0020.fits")])

    offsets = testutil.get_offsets()

    # high tolerance due to 1 pixel difference in x dimension in fit of 2nd target
    testutil.assert_tolerance(offsets[0],  -0.591, tolerance=0.1)
    testutil.assert_tolerance(offsets[1],   5.741, tolerance=0.1)
    testutil.assert_tolerance(offsets[2], -19.347, tolerance=0.92)

@with_setup(testutil.setup_io("4\nno"), clear_cursor_positions)
def test_first_target_to_center():
    add_cursor_position(508.0, 1148.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(522.0, 1099.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N20060131S0020.fits")])

    offsets = testutil.get_offsets()

    # high tolerance due to 1 pixel difference in x dimension in fit of 2nd target
    testutil.assert_tolerance(offsets[0],   0.799, tolerance=0.2)
    testutil.assert_tolerance(offsets[1],   1.781, tolerance=0.2)
    testutil.assert_tolerance(offsets[2], -19.375, tolerance=0.94)

@with_setup(testutil.setup_io("4\nyes\nyes"), clear_cursor_positions)
def test_move_exact_far_to_put_midpoint_at_center():
    add_cursor_position(508.047, 1148.124, "x")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(522.001, 1099.091, "x")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N20060131S0020.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0],  -0.617, tolerance=0.001)
    testutil.assert_tolerance(offsets[1],   5.402, tolerance=0.001)
    testutil.assert_tolerance(offsets[2], -15.882, tolerance=0.004)

@with_setup(testutil.setup_io("4\nno"), clear_cursor_positions)
def test_move_exact_first_target_to_center():
    add_cursor_position(508.047, 1148.124, "x")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(522.001, 1099.091, "x")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N20060131S0020.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0],   0.512, tolerance=0.002)
    testutil.assert_tolerance(offsets[1],   1.431, tolerance=0.001)
    testutil.assert_tolerance(offsets[2], -15.882, tolerance=0.004)
