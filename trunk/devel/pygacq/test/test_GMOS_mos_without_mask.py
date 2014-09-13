from nose import with_setup
from nose.tools import raises
import testutil

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_4_boxes():
    add_cursor_position(98.0, 312.0, "a")
    add_cursor_position(298.0, 110.0, "a")
    add_cursor_position(96.0, 108.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "--mosmasknum", "6", testutil.get_data_file_name("N20121108S0357.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.443, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], -1.236, tolerance=0.01)
    testutil.assert_tolerance(offsets[2],  0.116, tolerance=0.004)

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_should_not_be_able_to_manually_specify_a_box():
    add_cursor_position(98.0, 312.0, "a")
    add_cursor_position(96.0, 108.0, "a")
    add_cursor_position(100.0, 100.0, "b")
    add_cursor_position(298.0, 110.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "--mosmasknum", "6", testutil.get_data_file_name("N20121108S0357.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.443, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], -1.236, tolerance=0.01)
    testutil.assert_tolerance(offsets[2],  0.116, tolerance=0.004)

## too difficult to pick up any of these on different machines
## @with_setup(testutil.setup_io(""), clear_cursor_positions)
## def test_gmos_south_crowded_field_some_specified():
##     add_cursor_position(309.0,  70.0, "a")

##     add_cursor_position(0.0, 0.0, "q")

##     ret = main(["gacq_test", "--mosmasknum", "36", testutil.get_data_file_name("S20121103S0144.fits")])

##     offsets = testutil.get_offsets()
##     testutil.assert_tolerance(offsets[0], -1.457, tolerance=0.030)
##     testutil.assert_tolerance(offsets[1],  4.487, tolerance=0.027)
##     testutil.assert_tolerance(offsets[2],  0.009, tolerance=0.010)

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_south_crowded_field_all_specified():
    # almost the same test as above, but with tighter tolerances since
    # the user is more exactly specifying each star location
    add_cursor_position(110.0,  70.0, "a")
    add_cursor_position(111.0,  270.0, "a")
    add_cursor_position(309.0,  70.0, "a")
    add_cursor_position(311.0, 270.0, "a")

    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "--mosmasknum", "36", testutil.get_data_file_name("S20121103S0144.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -1.457, tolerance=0.01)
    testutil.assert_tolerance(offsets[1],  4.487, tolerance=0.008)
    testutil.assert_tolerance(offsets[2],  0.009, tolerance=0.001)

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_large_objects():
    add_cursor_position(157.0, 140.0, "a")
    add_cursor_position(157.0, 345.0, "a")
    add_cursor_position(358.0, 340.0, "a")

    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "--mosmasknum", "2", testutil.get_data_file_name("N20120219S0004.fits")])

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0], -4.175, tolerance=0.039)
    testutil.assert_tolerance(offsets[1], -3.014, tolerance=0.03)
    testutil.assert_tolerance(offsets[2],  0.075, tolerance=0.008)

@with_setup(testutil.setup_io("yes\nyes"), clear_cursor_positions)
def test_gmos_north_with_poor_mask_alignment():
    add_cursor_position(70.0,   77.0, "a")
    add_cursor_position(57.0,  276.0, "a")
    add_cursor_position(262.0, 266.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "--mosmasknum", "3", testutil.get_data_file_name("N20060131S0017.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  2.714, tolerance=0.001)
    testutil.assert_tolerance(offsets[1],  1.927, tolerance=0.02)
    testutil.assert_tolerance(offsets[2], -0.120, tolerance=0.001)

# Auto detection without the mask isn't good enough to always work
## @with_setup(testutil.setup_io("yes\nyes"), clear_cursor_positions)
## def test_gmos_north_with_poor_mask_alignment_use_one_autodetected():
##     add_cursor_position(70.0,   77.0, "a")
##     add_cursor_position(262.0, 266.0, "a")
##     add_cursor_position(0.0, 0.0, "q")

##     ret = main(["gacq_test", "--mosmasknum", "3", testutil.get_data_file_name("N20060131S0017.fits")])

##     offsets = testutil.get_offsets()
##     testutil.assert_tolerance(offsets[0],  2.714, tolerance=0.006)
##     testutil.assert_tolerance(offsets[1],  1.927, tolerance=0.009)
##     testutil.assert_tolerance(offsets[2], -0.120, tolerance=0.002)

@with_setup(testutil.setup_io("yes\nyes\nyes"), clear_cursor_positions)
def test_gmos_north_with_poor_mask_alignment_discard_first_box():
    add_cursor_position(70.0,   77.0, "a")
    add_cursor_position(57.0,  276.0, "a")
    add_cursor_position(70.0,   77.0, "d")
    add_cursor_position(262.0, 266.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "--mosmasknum", "3", testutil.get_data_file_name("N20060131S0017.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  2.714, tolerance=0.06)
    testutil.assert_tolerance(offsets[1],  1.927, tolerance=0.3)
    testutil.assert_tolerance(offsets[2], -0.120, tolerance=0.021)
