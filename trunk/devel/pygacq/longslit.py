import os
import math
import scipy
import scipy.stats
import numpy as np
from skimage import measure

import plotting
import userinterface as ui

from util import angle, clamp, find_optimal, MeasurementFailed
from gacqlogging import *

STRIP_SIZE = np.array([50, 40])

class SlitMeasurement:
    def __init__(self, slitdata, contours, contour_level, xmin, ymin, slitxpos, slitypos, left_border, rght_border):
        self.slitdata = slitdata
        self.contours = contours
        self.contour_level = contour_level
        self.origin = np.array([xmin, ymin]) + 1.0 # accounts for slitdata only being a subset, and converting to 1-based indicing
        self.left_border = left_border
        self.rght_border = rght_border
        self.slitxpos = slitxpos # not used for vertical slits
        self.slitypos = slitypos # should not be used for horizontal slits

    def get_data(self):
        return self.slitdata

    def get_contours(self):
        return self.contours

    def get_contour_level(self):
        return self.contour_level

    def get_point_on_line(self, border, ycrd):
        ycrd -= self.origin[1]
        return border.get_point_on_line(ycrd) + self.origin

    def get_point_on_left_line(self, ycrd):
        return self.get_point_on_line(self.left_border, ycrd)

    def get_point_on_rght_line(self, ycrd):
        return self.get_point_on_line(self.rght_border, ycrd)

    def get_width(self):
        left = self.get_point_on_left_line(self.slitypos)
        rght = self.get_point_on_rght_line(self.slitypos)

        return np.linalg.norm(rght - left)

    def get_average_position(self, ycrd):
        left = self.get_point_on_left_line(ycrd)
        rght = self.get_point_on_rght_line(ycrd)

        return (left + rght) / 2.0

    def get_midpoint(self):
        return self.get_average_position(self.slitypos)

    def get_borders(self):
        yield self.left_border
        yield self.rght_border

class SlitBorder:
    def __init__(self, y_length, slope, intercept, r_value, p_value, std_err):
        self.y_length = y_length
        self.slope = slope
        self.intercept = intercept
        self.r_value = r_value
        self.p_value = p_value
        self.std_err = std_err

    def get_point_on_line(self, ycrd):
        xcrd = self.get_slope() * ycrd + self.get_intercept()
        return np.array([xcrd, ycrd])

    def get_intercept(self):
        return self.intercept

    def get_slope(self):
        return self.slope

def compare_y_coordinate(point1, point2):
    return cmp(point1[1], point2[1])

def get_slit_tilt(mid1, mid2):
    points = [mid1, mid2]
    points.sort(cmp=compare_y_coordinate)
    bottom, top = points

    slitvec = top - bottom
    up = np.array([0, 1])

    return angle(up, slitvec)

def get_predicted_slit(midpoint, slitwidth):
    dx = slitwidth / 2
    dy = STRIP_SIZE[1] / 2
    
    top_left    = midpoint + np.array([-dx,  dy])
    bottom_left = midpoint + np.array([-dx, -dy])

    top_rght    = midpoint + np.array([ dx,  dy])
    bottom_rght = midpoint + np.array([ dx, -dy])

    return bottom_left, top_left, top_rght, bottom_rght, bottom_left

def get_contour_length(dim):
    return dim.max() - dim.min()

def compare_slit_borders(lhs, rhs):
    return cmp(lhs.get_intercept(), rhs.get_intercept())

def compare_contours(lhs, rhs):
    return cmp(len(lhs), len(rhs))

def plot_slit_measurement(acqimage, slit):
    import matplotlib.pyplot as plt
    plt.ion()

    # when in debug mode, render into separate windows
    if is_debug_mode():
        fig = plt.figure(figsize=(8, 3))
    else:
        fig = plt.gcf()
        plt.clf()
        
    ax = fig.gca()
        
    plt.subplot(121)
    plt.imshow(slit.get_data(), cmap=plt.cm.gray)

    msg = "mask : %s\ncontour : %.1f" % (str(acqimage.focal_plane_mask()),
                                         slit.get_contour_level())
    
    plt.text(0.0, 0.0,
             msg,
             fontsize=16,
             horizontalalignment='left',
             verticalalignment='top',
             color='white',
             transform=ax.transAxes)

    lines = []
    for contour in slit.get_contours():
        x = contour[:, 1]
        y = contour[:, 0]

        plt.plot(x, y, linewidth=2)

        x, y = y, x
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
        lines.append((slope, intercept))

    plt.axis('image')
    plt.xticks([])
    plt.yticks([])

    plt.subplot(122)
    slitdata = slit.get_data()
    plt.imshow(slitdata, cmap=plt.cm.gray)

    ydim, xdim = slitdata.shape

    for border in slit.get_borders():
        x0, y0 = border.get_point_on_line(0)
        x1, y1 = border.get_point_on_line(ydim)
        plt.plot([x0, x1], [y0, y1], linewidth=2)

    plt.axis('image')
    plt.xticks([])
    plt.yticks([])

    msg = "measured width : %.2f" % (slit.get_width() * acqimage.unbinned_pixel_scale())
    plt.text(0.95, 0.05,
             msg,
             fontsize=16,
             horizontalalignment='right',
             verticalalignment='bottom',
             color='white',
             transform=ax.transAxes)
    plt.draw()
    return fig

def measure_slit(scidata, contour_level, slitxpos, slitypos):
    # cut out the slit data
    pixel_buffer = 50
    xmin, xmax = clamp(scidata, slitxpos, pixel_buffer, 1)

    # cut a strip along the y-axis, small strips are better at avoiding bright objects
    pixel_buffer = STRIP_SIZE[1] / 2
    ymin, ymax = clamp(scidata, slitypos, pixel_buffer, 0)

    slitdata = scidata[ymin:ymax,xmin:xmax]

    contours = measure.find_contours(slitdata, contour_level)
    contours.sort(cmp=compare_contours, reverse=True)

    if len(contours) < 2:
        raise MeasurementFailed("Only one contour found at the level of %f, try somewhere else without a bright object" % contour_level)

    # the first two largest contours should be the borders of the slit
    borders = []
    for contour in contours[:2]:
        x = contour[:, 1]
        y = contour[:, 0]

        y_length = get_contour_length(y)
        limit = (STRIP_SIZE[1] / 2.0)

        if y_length < limit:
            raise MeasurementFailed("One of the edge contours is less than half the height of the strip being searched")

        # for vertical slits, swap the inputs, linregress works better
        x, y = y, x
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

        # adjust for the slicing done
        borders.append(SlitBorder(y_length, slope, intercept, r_value, p_value, std_err))

    borders.sort(cmp=compare_slit_borders)
    slit = SlitMeasurement(slitdata,
                           contours,
                           contour_level,
                           xmin, ymin,
                           slitxpos, slitypos,
                           *borders)

    return slit

def get_error_func(actual_width, pixscale):
    def error_func(slit):
        return actual_width - slit.get_width() * pixscale
    return error_func

def find_optimal_slit(acqimage, slitxpos, slitypos, display):
    actual_width = acqimage.get_mask_width()
    
    scidata = acqimage.get_science_data()

    best_slit = find_optimal(scidata,
                             measure_slit,
                             get_error_func(actual_width, acqimage.unbinned_pixel_scale()),
                             slitxpos,
                             slitypos)

    if display:
        plotting.put_plot(plot_slit_measurement, (acqimage, best_slit))
    
    return best_slit

def find_slit(slitimage_ad, slitxpos, slitypos, display):
    debug("...slitxpos = ", slitxpos)
    debug("...slitypos = ", slitypos)

    # find a long slit
    measurement = find_optimal_slit(slitimage_ad, slitxpos, slitypos, display)
    
    # display the center
    midpoint = measurement.get_midpoint()
    slitxpos, slitypos = midpoint
    info("Measured slit center: %8.2f,%8.2f" % (slitxpos, slitypos))
    
    return measurement

def draw_slit(midpoint, width):
    # mark the center
    mark = ui.marker(midpoint, color=ui.MAGENTA)
    
    # display a box 
    box = get_predicted_slit(midpoint, width)
    line = ui.polyline(points=box, color=ui.MAGENTA)
    return [mark, line]

class ManualSlitMeasurement(object):
    def __init__(self, midpoint):
        self.midpoint = midpoint

    def get_midpoint(self):
        return self.midpoint

class SlitMeasurementCursorListener(ui.CursorListener):
    def __init__(self, slitimage_ad, args):
        self.slitimage_ad = slitimage_ad
        self.args = args
        self.slit_measurement = None
        self.overlays = []
        self.left_side = None

    def maybe_undo(self):
        ui.undo(self.overlays)

    def handle_slit_guess(self, pos):
        if pos.get_keystroke() == "a":
            self.maybe_undo()

        if self.left_side is not None:
            self.left_side = None
        
        xpos, ypos = pos.get_position()

        display = self.args.verbose or pos.get_keystroke() == "e"
        measurement = find_slit(self.slitimage_ad, xpos, ypos, display)

        if pos.get_keystroke() == "a":
            self.slit_measurement = measurement
            self.overlays = draw_slit(measurement.get_midpoint(), self.slitimage_ad.get_mask_width_in_pixels())

    def handle_exact_location(self, pos):
        self.maybe_undo()

        coords = pos.get_position()
        if self.left_side is None:
            self.left_side = coords
            mark = ui.marker(coords, color=ui.MAGENTA)
            self.overlays = [mark]
        else:
            midpoint = (self.left_side + coords) / 2.0
            slit = draw_slit(midpoint, self.slitimage_ad.get_mask_width_in_pixels())
            self.slit_measurement = ManualSlitMeasurement(midpoint)
            
            self.overlays = slit            

    def send_cursor_position(self, pos):
        key = pos.get_keystroke()
        if   key in ["a", "e"]:
            self.handle_slit_guess(pos)
        elif key == "x":
            self.handle_exact_location(pos)
        else:
            warning("Unrecognized keystroke '%s', doing nothing..." % key)
            return

    def get_slit_measurement(self):
        return self.slit_measurement
