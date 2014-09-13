from nose import with_setup
from nose.tools import raises
import testutil

from integratedfieldunit import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_ifu_R():
    add_cursor_position(10.5, 101.5, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(9.6, 11.2, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20111126S0354.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_ifu_2():
    add_cursor_position(12.0, 97.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(19.0, 14.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20120122S0031.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_long_object():
    add_cursor_position(12.0, 99.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(9.0, 13.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130605S0153.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_extended_fuzzy_object():
    add_cursor_position(8.6, 98.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(8.0, 14.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130716S0160.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_ifu_2_flat():
    add_cursor_position(10.3, 102.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(20.0, 10.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130606S0009.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_north_ifu_R_flat():
    add_cursor_position(9.0, 101.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(10.0, 10.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130717S0014.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_south_ifu_R():
    add_cursor_position(9.0, 101.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(10.0, 10.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0002.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_south_ifu_2_flat():
    add_cursor_position(7.0, 108.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(20.0, 10.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130110S0028.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_south_ifu_2_newer_flat():
    add_cursor_position(11.0, 102.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(20.0, 10.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130725S0014.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_south_ifu_2_diffuse_object():
    add_cursor_position(10.5, 106.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(20.0, 10.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130421S0035.fits")])

@with_setup(testutil.setup_io(""), clear_cursor_positions)
def test_gmos_south_ifu_2_bright_star():
    add_cursor_position(10.0, 108.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(19.0, 13.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130227S0313.fits")])
