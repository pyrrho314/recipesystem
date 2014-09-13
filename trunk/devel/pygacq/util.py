import math
import numpy as np

def fits_filename(fname):
    if not fname.endswith(".fits"):
        fname += ".fits"
    return fname

def is_number(val):
    try:
        int(val)
        return True
    except ValueError:
        return False

def clamp(data, coord, pixel_buffer, dimension):
    assert coord < data.shape[dimension], "%r >= %r" % (coord, data.shape[dimension])
    
    crd = int(round(coord)) - 1

    crdmin = crd - pixel_buffer
    crdmin = max(crdmin, 0)

    crdmax = crd + pixel_buffer + 1
    crdmax = min(crdmax, data.shape[dimension])

    return crdmin, crdmax

def clamp_to_size(data, coord, size, dimension):
    crd = coord - 1.0 # convert from 1-based indicing to 0-based indicing

    crdmin = math.ceil(crd - float(size) / 2.0)
    crdmin = max(crdmin, 0)

    crdmax = crdmin + size
    crdmax = min(crdmax, data.shape[dimension])

    return crdmin, crdmax

def get_window(data, coord, pixel_buffer):
    xcrd, ycrd = coord
    xmin, xmax = clamp(data, xcrd, pixel_buffer, 1)
    ymin, ymax = clamp(data, ycrd, pixel_buffer, 0)
    return xmin, xmax, ymin, ymax

def angle(vec1, vec2):
    """ Calculate the angle to rotate vec2 onto vec1, including the direction. """
    dotprod = vec1.dot(vec2)
    mag = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    rad = np.arccos(dotprod / mag)

    # determine the rotation direction
    if np.cross(vec1, vec2) > 0.0:
        rad = -rad
    
    deg = np.degrees(rad)
    # Give the phase which has the smallest magnitude
    if deg < -180.0:
        deg = deg + 360.0
    if deg >  180.0:
        deg = deg - 360.0

    return deg

def rotate(vec, theta):
    r  = np.radians(theta)
    Rx = np.array([math.cos(r), -math.sin(r)])
    Ry = np.array([math.sin(r),  math.cos(r)])
    return np.array([np.dot(Rx, vec), np.dot(Ry, vec)])

class MeasurementFailed(Exception):
    pass

def make_guess(dmin, dmax):
    next_guess = dmin + ((dmax - dmin) / 2)
    return dmin, next_guess, dmax

def find_optimal(scidata, measure_func, error_func, *args):
    next_dmin, next_guess, next_dmax = make_guess(scidata.min(), scidata.max())

    slits = []
    for i in range(10):
        guess, dmin, dmax = next_guess, next_dmin, next_dmax

        try:
            slit = measure_func(scidata, guess, *args)
            slits.append(slit)

            diff = error_func(slit)
            if diff >= 0.0:
                # Going smaller
                next_dmin, next_guess, next_dmax = make_guess(dmin, guess)
            else:
                # Going larger
                next_dmin, next_guess, next_dmax = make_guess(guess, dmax)
                
        except MeasurementFailed: # less than required number of contours throws an exception
            # Going smaller
            next_dmin, next_guess, next_dmax = make_guess(dmin, guess)

    def compare(s1, s2):
        diff1 = abs(error_func(s1))
        diff2 = abs(error_func(s2))
        return cmp(diff1, diff2)

    slits.sort(cmp=compare)
    return slits[0]
