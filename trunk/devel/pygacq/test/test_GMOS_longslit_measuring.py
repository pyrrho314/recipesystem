import sys
import testutil
from nose import with_setup

from gacq import main
from userinterface import add_cursor_position, clear_cursor_positions

##################################################################
# all these tests automatically detect long-slit mode because they
# have a companion slit image
##################################################################

@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_measuring_longslit_from_image():
    add_cursor_position(152.0, 144.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok

    ret = main(["gacq_test", testutil.get_data_file_name("N20111127S0001.fits")]) # N20111127S0002.fits is the image of the slit

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0], -0.210, tolerance=0.08)
    testutil.assert_tolerance(offsets[1],  0.342, tolerance=0.003)

@with_setup(testutil.setup_io("yes"), clear_cursor_positions) 
def test_move_object_far_to_slit():
    add_cursor_position(526.0, 1186.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok

    ret = main(["gacq_test", testutil.get_data_file_name("S20110804S0046.fits")]) # S20110804S0047.fits is the image of the slit

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0], -1.749, tolerance=0.1)
    testutil.assert_tolerance(offsets[1],  4.925, tolerance=0.005)


@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_really_long_longslit_measuring():
    add_cursor_position(513.0, 1150.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok

    ret = main(["gacq_test", testutil.get_data_file_name("N20120216S0231.fits")]) # N20120216S0232 is the image of the slit

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0], -0.121, tolerance=0.06)
    testutil.assert_tolerance(offsets[1],  0.285, tolerance=0.02)

@with_setup(testutil.setup_io("yes\nyes"), clear_cursor_positions)
def test_twotarget_longslit_measuring():
    add_cursor_position(513.0, 1150.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(514.0, 1944.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N20120216S0231.fits")]) # N20120216S0232 is the image of the slit

    offsets = testutil.get_offsets()

    # high tolerance due to 1 pixel difference in x dimension in fit of 2nd target
    # plus the fact that the slit tilt is now taken into account
    testutil.assert_tolerance(offsets[0],  -0.217, tolerance=0.05)
    testutil.assert_tolerance(offsets[1], -57.435, tolerance=0.03)
    testutil.assert_tolerance(offsets[2],   0.105, tolerance=0.06)

@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_twotarget_longslit_measuring_with_custom_center():
    add_cursor_position(513.0,  1150.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(514.0,  1944.0, "a")
    add_cursor_position(513.01, 1544.03, "c")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", "-2", testutil.get_data_file_name("N20120216S0231.fits")]) # N20120216S0232 is the image of the slit

    offsets = testutil.get_offsets()

    testutil.assert_tolerance(offsets[0],  -0.087, tolerance=0.05)
    testutil.assert_tolerance(offsets[1],  -0.434, tolerance=0.03)
    testutil.assert_tolerance(offsets[2],   0.105, tolerance=0.004)

@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_measuring_longslit_from_bright_object():
    add_cursor_position(149.0, 147.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok

    ret = main(["gacq_test", testutil.get_data_file_name("N20110501S0008.fits")]) # N20110501S0009.fits is the image of the slit

    offsets = testutil.get_offsets()

    # high tolerance is because gacq.cl wasn't getting the correct center from automatic detection
    # this calculation is better, but still not perfect, the next test will test the user interactivity
    testutil.assert_tolerance(offsets[0], -0.113, tolerance=0.14)
    testutil.assert_tolerance(offsets[1],  0.153, tolerance=0.001)

@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_manaully_measuring_longslit_from_bright_object():
    add_cursor_position(149.0, 147.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(147.0, 129.0, "a")
    add_cursor_position(0.0, 0.0, "q")

    ret = main(["gacq_test", testutil.get_data_file_name("N20110501S0008.fits")]) # N20110501S0009.fits is the image of the slit

    offsets = testutil.get_offsets()

    # high tolerance is because gacq.cl wasn't getting the correct center from automatic detection
    # this calculation is better, but still not perfect, the next test will test the user interactivity
    testutil.assert_tolerance(offsets[0], -0.211, tolerance=0.004)
    testutil.assert_tolerance(offsets[1],  0.153, tolerance=0.001)

@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_i_scoff_at_the_bright_object_in_the_slit():
    add_cursor_position(149.0, 138.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok

    ret = main(["gacq_test", testutil.get_data_file_name("N20110529S0087.fits")]) # N20110529S0088.fits is the image of the slit

    offsets = testutil.get_offsets()

    # low tolerance is because gacq.cl wasn't actually doing a good
    # job of detecting the slit, gacq.py does much better
    testutil.assert_tolerance(offsets[0], -0.193, tolerance=0.001)
    testutil.assert_tolerance(offsets[1],  0.767, tolerance=0.005)

@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_gmos_south_bright_object_in_the_slit():
    add_cursor_position(505.0, 1145.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(0.0, 0.0, "q") # slit measurement is ok

    ret = main(["gacq_test", testutil.get_data_file_name("S20110428S0015.fits")]) # S20110428S0016.fits is the image of the slit

    offsets = testutil.get_offsets()

    # low tolerance is because gacq.cl wasn't actually doing a good
    # job of detecting the slit, gacq.py does much better
    testutil.assert_tolerance(offsets[0],  0.916, tolerance=0.001)
    testutil.assert_tolerance(offsets[1],  0.979, tolerance=0.01)

@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_gmos_manual_slit_measurement():
    add_cursor_position(505.0, 1145.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(144.073, 109.062, "x")
    add_cursor_position(155.008, 108.969, "x")
    add_cursor_position(0.0, 0.0, "q") 

    ret = main(["gacq_test", testutil.get_data_file_name("S20110428S0015.fits")]) # S20110428S0016.fits is the image of the slit

    offsets = testutil.get_offsets()

    # low tolerance is because gacq.cl wasn't actually doing a good
    # job of detecting the slit, gacq.py does much better
    testutil.assert_tolerance(offsets[0],  0.944, tolerance=0.001)
    testutil.assert_tolerance(offsets[1],  0.979, tolerance=0.01)


@with_setup(testutil.setup_io(""), clear_cursor_positions) 
def test_gmos_south_slit_center_not_correctly_detected_the_first_time():
    add_cursor_position(514.0, 1153.0, "a")
    add_cursor_position(0.0, 0.0, "q")
    add_cursor_position(147.0, 114.0, "a") # slit measurement needs to be adjusted
    add_cursor_position(0.0, 0.0, "q") 

    ret = main(["gacq_test", testutil.get_data_file_name("S20090418S0109.fits")]) # S20090418S0110.fits is the image of the slit

    offsets = testutil.get_offsets()

    # low tolerance is because gacq.cl wasn't actually doing a good
    # job of detecting the slit, gacq.py does much better
    testutil.assert_tolerance(offsets[0], -0.461, tolerance=0.01)
    testutil.assert_tolerance(offsets[1], -0.160, tolerance=0.01)
