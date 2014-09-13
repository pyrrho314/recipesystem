import sys
import testutil
from nose import with_setup

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io("2"), clear_cursor_positions)
def test_gmos_north_ifu_2_imaging():
    add_cursor_position(516.0, 1154.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20120122S0030.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  29.013, tolerance=0.004)
    testutil.assert_tolerance(offsets[1],  -1.457, tolerance=0.01)

@with_setup(testutil.setup_io("2"), clear_cursor_positions)
def test_gmos_north_ifu_2_imaging_exact_center():
    add_cursor_position(516.0, 1154.0, "a")
    add_cursor_position(716.0, 1144.0, "c")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20120122S0030.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  29.017, tolerance=0.002)
    testutil.assert_tolerance(offsets[1],  -1.396, tolerance=0.01)


@with_setup(testutil.setup_io("1"), clear_cursor_positions)
def test_gmos_north_ifu_R_imaging():
    add_cursor_position(516.0, 1154.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20120122S0030.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  30.763, tolerance=0.004)
    testutil.assert_tolerance(offsets[1],  -1.457, tolerance=0.01)

@with_setup(testutil.setup_io("1"), clear_cursor_positions)
def test_gmos_north_ifu_R_imaging_on_a_stamp():
    add_cursor_position(154.0, 152.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20111126S0353.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  31.014, tolerance=0.040)
    testutil.assert_tolerance(offsets[1],  -1.457, tolerance=0.006)

@with_setup(testutil.setup_io("yes\n2"), clear_cursor_positions)
def test_gmos_south_ifu_2_imaging():
    add_cursor_position(302.0, 1160.0, "x")
    add_cursor_position(307.0, 1168.0, "c")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-v", testutil.get_data_file_name("S20090326S0024.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], 0.727, tolerance=0.004)
    testutil.assert_tolerance(offsets[1], 1.165, tolerance=0.006)

@with_setup(testutil.setup_io("yes\n1"), clear_cursor_positions)
def test_gmos_south_ifu_R_imaging():
    add_cursor_position(302.0, 1160.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-v", testutil.get_data_file_name("S20090326S0024.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], 2.471, tolerance=0.004)
    testutil.assert_tolerance(offsets[1], 1.035, tolerance=0.007)
