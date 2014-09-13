from nose import with_setup
from nose.tools import raises
import testutil

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_with_a_weird_artifact():
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130114S0062.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.028, tolerance=0.004)
    testutil.assert_tolerance(offsets[1],  0.171, tolerance=0.01)
    testutil.assert_tolerance(offsets[2],  0.015, tolerance=0.008)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_with_a_weird_artifact_discard_top_left():
    add_cursor_position(20.0, 150.0, "d")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130114S0062.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.028, tolerance=0.008)
    testutil.assert_tolerance(offsets[1],  0.171, tolerance=0.053)
    testutil.assert_tolerance(offsets[2],  0.015, tolerance=0.026)

@raises(SystemExit)
@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_with_a_weird_artifact_discard_too_much():
    add_cursor_position(20.0, 150.0, "d")
    add_cursor_position(150.0, 150.0, "d")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20130114S0062.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.028, tolerance=0.008)
    testutil.assert_tolerance(offsets[1],  0.171, tolerance=0.053)
    testutil.assert_tolerance(offsets[2],  0.015, tolerance=0.026)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_from_entire_field_of_view():
    add_cursor_position(44.0, 64.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20090422S0074.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.271, tolerance=0.018)
    testutil.assert_tolerance(offsets[1], -0.473, tolerance=0.005)
    testutil.assert_tolerance(offsets[2], -0.080, tolerance=0.003)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_from_entire_field_of_view_with_exact_star_location():
    add_cursor_position(45.0, 61.0, "x")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20090422S0074.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.271, tolerance=0.019)
    testutil.assert_tolerance(offsets[1], -0.473, tolerance=0.006)
    testutil.assert_tolerance(offsets[2], -0.080, tolerance=0.002)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_exact_star_location_should_escape_manual_box_specification():
    add_cursor_position(45.0, 61.0, "b")
    add_cursor_position(45.0, 61.0, "x")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20090422S0074.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.271, tolerance=0.019)
    testutil.assert_tolerance(offsets[1], -0.473, tolerance=0.006)
    testutil.assert_tolerance(offsets[2], -0.080, tolerance=0.002)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_exact_star_location_should_escape_discarding_box():
    add_cursor_position(45.0, 61.0, "d")
    add_cursor_position(45.0, 61.0, "x")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20090422S0074.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.271, tolerance=0.019)
    testutil.assert_tolerance(offsets[1], -0.473, tolerance=0.006)
    testutil.assert_tolerance(offsets[2], -0.080, tolerance=0.002)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_north_ignore_offsets_and_start_science():
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130419S0270.fits")])

    # ignore offsets and start science immediately
    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], 0.019, tolerance=0.004)
    testutil.assert_tolerance(offsets[1], 0.073, tolerance=0.010)
    testutil.assert_tolerance(offsets[2], 0.005, tolerance=0.010)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_north_7_boxes():
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20111123S0033.fits")])

    # ignore offsets and start science immediately
    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.030, tolerance=0.03)
    testutil.assert_tolerance(offsets[1], -0.010, tolerance=0.006)
    testutil.assert_tolerance(offsets[2],  0.005, tolerance=0.002)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_old_gmos_north_with_3_boxes():
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20060131S0015.fits")])

    # ignore offsets and start science immediately, looks like the
    # errors in star/box positions offset each other
    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.029, tolerance=0.014)
    testutil.assert_tolerance(offsets[1],  0.024, tolerance=0.008)
    testutil.assert_tolerance(offsets[2], -0.012, tolerance=0.003)

## @with_setup(testutil.setup_io("yes"), clear_cursor_positions)
## def test_old_gmos_north_with_coordinate_file():
##     add_cursor_position(0.0, 0.0, "q")

##     ret = main(["gacq_test", testutil.get_data_file_name("N20060131S0016.fits")])

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_off_center_boxes():
    add_cursor_position(55.0, 58.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0001.fits")])

    # higher tolerance for tough case, new python version gets it "more correct" than automaskdetect

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.512, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], -0.381, tolerance=0.05)
    testutil.assert_tolerance(offsets[2],  0.135, tolerance=0.012)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_off_center_boxes_escape_out_of_manual_box_specification():
    # this test should ignore the attempt to manually specify the box
    add_cursor_position(55.0, 58.0, "b")
    add_cursor_position(55.0, 58.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0001.fits")])

    # higher tolerance for tough case, new python version gets it "more correct" than automaskdetect

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.512, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], -0.381, tolerance=0.05)
    testutil.assert_tolerance(offsets[2],  0.135, tolerance=0.012)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_off_center_boxes_escape_out_of_box_toggling():
    # this test should ignore the attempt to discard the box
    add_cursor_position(55.0, 58.0, "d")
    add_cursor_position(55.0, 58.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0001.fits")])

    # higher tolerance for tough case, new python version gets it "more correct" than automaskdetect

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.512, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], -0.381, tolerance=0.05)
    testutil.assert_tolerance(offsets[2],  0.135, tolerance=0.012)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_port_one():
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20110804S0079.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.321, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], -0.004, tolerance=0.01)
    testutil.assert_tolerance(offsets[2],  0.027, tolerance=0.002)


@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_manual_box_specification():
    add_cursor_position(26.0, 32.0, "b")
    add_cursor_position(62.0, 67.0, "b")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0001.fits")])

    # higher tolerance for tough case, new python version gets it "more correct" than automaskdetect

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.512, tolerance=0.06)
    testutil.assert_tolerance(offsets[1], -0.381, tolerance=0.040)
    testutil.assert_tolerance(offsets[2],  0.135, tolerance=0.020)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_south_manual_box_specification_other_corners():
    add_cursor_position(62.0, 32.0, "b")
    add_cursor_position(26.0, 67.0, "b")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("S20071017S0001.fits")])

    # higher tolerance for tough case, new python version gets it "more correct" than automaskdetect

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.512, tolerance=0.06)
    testutil.assert_tolerance(offsets[1], -0.381, tolerance=0.040)
    testutil.assert_tolerance(offsets[2],  0.135, tolerance=0.020)

from random import Random

@raises(SystemExit)
@with_setup(testutil.setup_io("yes\nno"), clear_cursor_positions)
def test_7_box_monkey_test():
    rand = Random(8799)

    for i in range(100):
        x = rand.choice(range(350))
        y = rand.choice(range(350))
        key = rand.choice(["a", "r", "e", "x", "b", "d"])
        add_cursor_position(x, y, key)

    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20111123S0033.fits")])

# input:
# offsets ok = yes
# continue even though mask is found = yes
# MOS mode = 3
# mask num = 1
@with_setup(testutil.setup_io("yes\nyes\n3\n1"), clear_cursor_positions)
def test_gmos_N20130616S0088():
    add_cursor_position(300.0, 303.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130616S0088.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0],  0.282, tolerance=0.01)
    testutil.assert_tolerance(offsets[1], -0.074, tolerance=0.01)
    testutil.assert_tolerance(offsets[2], -0.019, tolerance=0.02)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_N20130616S0089():
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130616S0089.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.219, tolerance=0.03)
    testutil.assert_tolerance(offsets[1],  0.102, tolerance=0.01)
    testutil.assert_tolerance(offsets[2],  0.010, tolerance=0.01)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions)
def test_gmos_N20130616S0090():
    add_cursor_position(148.0, 158.0, "r")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20130616S0090.fits")])

    offsets = testutil.get_offsets()
    testutil.assert_tolerance(offsets[0], -0.032, tolerance=0.02)
    testutil.assert_tolerance(offsets[1], -0.015, tolerance=0.02)
    testutil.assert_tolerance(offsets[2],  0.001, tolerance=0.01)
