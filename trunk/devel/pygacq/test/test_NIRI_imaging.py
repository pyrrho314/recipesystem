import sys
import testutil
from nose import with_setup

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io(), clear_cursor_positions)
def test_bright_object_that_does_not_converge():
    add_cursor_position(391.0, 539.0, "a")
    add_cursor_position(391.0, 539.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20060131S0011.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_almost_equals(offsets,
                                  (2.779, 0.672),
                                  places=2)

@with_setup(testutil.setup_io(), clear_cursor_positions)
def test_cursor_coordinates_should_round_to_nearest_integer():
    add_cursor_position(489.07, 478.99, "a")
    add_cursor_position(489.07, 478.99, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20060131S0012.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_almost_equals(offsets,
                                  (2.661, -3.874),
                                  places=2)

@with_setup(testutil.setup_io(), clear_cursor_positions)
def test_extra_cursor_coordinates_should_be_ignored():
    add_cursor_position(100.07, 100.99, "a")
    add_cursor_position(489.07, 478.99, "a")
    add_cursor_position(489.07, 478.99, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20060131S0012.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_almost_equals(offsets,
                                  (2.661, -3.874),
                                  places=2)

@with_setup(testutil.setup_io(), clear_cursor_positions)
def test_custom_center():
    add_cursor_position(489.07, 478.99, "a")
    add_cursor_position(512.00, 512.00, "c")
    add_cursor_position(489.07, 478.99, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20060131S0012.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_almost_equals(offsets,
                                  (2.661, -3.874),
                                  places=2)
