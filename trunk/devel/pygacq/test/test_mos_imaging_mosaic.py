from testutil import assert_tolerance, get_data_file_name
from nose.tools import assert_equals

import userinterface as ui

from acquisitionimage import AcquisitionImage, Box
from multiobjectspectroscopy import InteractiveMosaic, MOSBoxCursorListener

def test_mdf_box_detector_offset():
    acqimage = AcquisitionImage(get_data_file_name("N20121108S0357.fits"), mosmask="6")

    box_mosaic = acqimage.box_mosaic

    pnt = (10, 10)
    detsec = box_mosaic.detsec_finder.find_detector_section(pnt)
    box_size = box_mosaic.detsec_finder.get_box_size()

    box = Box(detsec, pnt, box_size, (0, 0))
    assert_tolerance(box.get_detector_offset(), (2, 2), tolerance=0.00001)
    

def test_mdf_box_locations_are_unchanged_by_acquisition_image():
    acqimage = AcquisitionImage(get_data_file_name("N20121108S0357.fits"), mosmask="6")

    box_centers = [(4370.9898354,   339.18367075),
                   (1948.16935825,  712.01011175),
                   (1384.48997275, 1181.42130025),
                   (1451.827332,   2218.3127036 )]

    assert len(box_centers) == len(acqimage.get_mos_boxes())

    for refcenter, box in zip(box_centers, acqimage.get_mos_boxes()):
        box_detector_center = box.get_mdf_detector_center() 
        assert_tolerance(refcenter, box_detector_center, tolerance=0.0001)


def test_mdf_box_locations_are_unchanged_by_interactive_mosaic():
    acqimage = AcquisitionImage(get_data_file_name("N20121108S0357.fits"), mosmask="6")

    verbose = True
    interactive_mosaic = InteractiveMosaic(acqimage, verbose)

    box_centers = [(4370.9898354,   339.18367075),
                   (1948.16935825,  712.01011175),
                   (1384.48997275, 1181.42130025),
                   (1451.827332,   2218.3127036 )]

    assert len(box_centers) == len(acqimage.get_mos_boxes())

    for refcenter, tile in zip(box_centers, interactive_mosaic.get_tiles()):
        tile_detector_center = tile.get_detector_center() 
        assert_tolerance(refcenter, tile_detector_center, tolerance=0.0001)


# Not reliable across platforms
## def test_automatically_detected_star_locations():
##     acqimage = AcquisitionImage(get_data_file_name("N20121108S0357.fits"), mosmask="6")

##     verbose = True
##     interactive_mosaic = InteractiveMosaic(acqimage, verbose)

##     ui.display(acqimage.get_science_data(), zscale=True)
    
##     for border in acqimage.get_mos_box_borders():
##         ui.polyline(points=border, color=ui.BLUE)

##     listener = MOSBoxCursorListener(interactive_mosaic)

##     ui.add_cursor_position(0.0, 0.0, "q")
##     listener.start()

##     star_centers = [(4361.82,  353.20),
##                     (1939.66,  732.32),
##                     (1374.08, 1201.80),
##                     (1445.88, 2238.26)]

##     assert_tolerance(star_centers, list(listener.get_star_centers()), tolerance=0.25)
