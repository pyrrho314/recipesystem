from nose.tools import assert_equals

from multiobjectspectroscopy import AcquisitionBox
import numpy as np

def floatize(c):
    return [map(float, crds) for crds in c]

def arrayize(points):
    return np.array([np.array(p) for p in points])

contour_level = 1.0 # meaningless for these tests
def test_AcquisitionBox_get_area_easy():
    contour = floatize([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    box = AcquisitionBox(contour, contour_level)
    assert_equals(box.get_area(), 1.0)

    contour.reverse()
    box = AcquisitionBox(contour, contour_level)
    assert_equals(box.get_area(), 1.0)


def test_AcquisitionBox_get_area():
    contour = floatize([(1, 1), (1, 2), (1.5, 2), (2, 2), (2, 1), (1, 1)])
    box = AcquisitionBox(contour, contour_level)
    assert_equals(box.get_area(), 1.0)

    contour.reverse()
    box = AcquisitionBox(contour, contour_level)
    assert_equals(box.get_area(), 1.0)

def test_AcquisitionBox_get_area_numpy_array():
    contour = floatize([(1, 1), (1, 2), (1.5, 2), (2, 2), (2, 1), (1, 1)])
    contour = arrayize(contour)
    print contour
    box = AcquisitionBox(contour, contour_level)
    assert_equals(box.get_area(), 1.0)

    contour = contour[::-1] # reverse direction
    box = AcquisitionBox(contour, contour_level)
    assert_equals(box.get_area(), 1.0)
