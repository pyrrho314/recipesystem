from testutil import assert_tolerance, get_data_file_name
from nose.tools import assert_equals

import userinterface as ui

from acquisitionimage import AcquisitionImage 
from multiobjectspectroscopy import find_optimal_box

#import gacqlogging
#gacqlogging.setup_logging(None, True)

# all centers in this file are from automaskdetect run through gacq.cl
# automaskdetect could only report centers up to half-pixel accuracy,
# so that's all we can hope for here.

def assert_mos_acquisition(acqimage, centers, tolerance):
    print acqimage.unbinned_pixel_scale()
    assert acqimage.has_mos_mask()
    assert_equals(acqimage.get_num_mos_boxes(), len(centers))

    data = acqimage.get_science_data()
    assert_equals(*data.shape)

    boxes = []
    z1 = Ellipsis
    z2 = None
    for box in acqimage.get_mos_boxes():
        acqbox = find_optimal_box(box, display=True)
        boxes.append(acqbox)
        box_z1, box_z2 = acqbox.get_zscale()

        z1 = min(box_z1, z1)
        z2 = max(box_z2, z2)

    ui.display(data, z1=z1, z2=z2, zscale=False)

    for border in acqimage.get_mos_box_borders():
        ui.polyline(points=border, color=ui.BLUE)

    for box in boxes:
        ui.polyline(points=box.get_mosaic_predicted_borders(), color=ui.MAGENTA)

        x_center, y_center = box.get_mosaic_center()
        assert 0 <= x_center < data.shape[1] 
        assert 0 <= y_center < data.shape[0]

    detector_centers = [b.get_detector_center() for b in boxes]
    assert_tolerance(detector_centers, centers, tolerance=tolerance)

def test_gmos_south_with_a_weird_artifact():
    acqimage = AcquisitionImage(get_data_file_name("S20130114S0062.fits"))

    centers = [(3923.0, 1384.5),
               (1708.6, 3424.4),
               (4955.4, 3834.4)]

    assert_mos_acquisition(acqimage, centers, 0.5)

def test_gmos_south_from_entire_field_of_view():
    centers = [(3956.0, 404.0),
               (3466.0, 596.5),
               (2397.0, 1626.5),
               (4949.5, 2795.0)]

    acqimage = AcquisitionImage(get_data_file_name("S20090422S0074.fits"))
    assert_mos_acquisition(acqimage, centers, 0.3)

def test_gmos_north():
    centers = [(4408.0, 1601.0),
               (1387.5, 2168.5),
               (3548.0, 3052.0)]

    acqimage = AcquisitionImage(get_data_file_name("N20130419S0270.fits"))
    assert_mos_acquisition(acqimage, centers, 0.5)

def test_gmos_north_7_boxes():
    centers = [(5008.5, 354.5),
               (1020.0, 611.9),
               (2357.0, 943.1),
               (1553.5, 2046.5),
               (4680.0, 3024.9),
               (1371.5, 3685.0),
               (4512.0, 4390.0)]

    acqimage = AcquisitionImage(get_data_file_name("N20111123S0033.fits"))
    assert_mos_acquisition(acqimage, centers, 0.5)

def test_old_gmos_north_with_3_boxes():
    centers = [(1187.0, 572.0),
               (5138.5, 3694.5),
               (1024.5, 4034.5)]

    acqimage = AcquisitionImage(get_data_file_name("N20060131S0015.fits"))
    assert_mos_acquisition(acqimage, centers, 0.4)

## def test_old_gmos_north_with_coordinate_file():
##     centers = [(),
##                (),
##                ()]

##     acqimage = AcquisitionImage(get_data_file_name("N20060131S0016.fits"))
##     assert_mos_acquisition(acqimage, centers, 0.3)

def test_gmos_south_off_center_boxes():
    centers = [(4780.5, 393.0),
               (4983.0, 3826.0),
               (1744.5, 4198.5)]

    acqimage = AcquisitionImage(get_data_file_name("S20071017S0001.fits"))

    # higher tolerance for tough case, new python version gets it "more correct" than automaskdetect
    assert_mos_acquisition(acqimage, centers, 0.7)

def test_gmos_south_port_one():
    centers = [(1330.0, 4220.0),
               (2708.5, 3518.5),
               (2990.0, 642.0)]

    acqimage = AcquisitionImage(get_data_file_name("S20110804S0079.fits"))

    assert_mos_acquisition(acqimage, centers, 0.2)

