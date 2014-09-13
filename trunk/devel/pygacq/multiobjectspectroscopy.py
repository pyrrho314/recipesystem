import os
import math
import scipy
import scipy.stats
import numpy as np
from skimage import measure

import plotting
import userinterface as ui
from gacqlogging import *
from util import find_optimal
from acquisitionimage import ACQUISITION_BOX_SIZE
from selection import get_selection_peak, mark_selection
from longslit import compare_contours


def is_closed_contour(contour):
    return np.all(contour[0] == contour[-1])

def get_borders_from_extents(extents):
    xmin, ymin, xmax, ymax = extents
    points = [(xmin, ymin),
              (xmin, ymax),
              (xmax, ymax),
              (xmax, ymin),
              (xmin, ymin)]
    return np.array(points)

def get_extents_from_center(center, pixscale):
    box_width = float(ACQUISITION_BOX_SIZE) / pixscale
    
    half_width = box_width / 2.0
    offset = np.array([half_width, half_width])
    
    xmin, ymin = center - offset
    xmax, ymax = center + offset
    return xmin, ymin, xmax, ymax


class AcquisitionBox(object):
    def __init__(self, contour, contour_level):
        if not is_closed_contour(contour):
            ValueError("%s expects a 'closed' contour, i.e., must have the same point at the beginning and end")

        if not isinstance(contour, np.ndarray):
            contour = np.array([np.array(p) for p in contour])
        
        self.contour = contour
        self.contour_level = contour_level

    def get_xy(self):
        x = self.contour[:, 1]
        y = self.contour[:, 0]
        return x, y
        
    def get_area(self):
        x, y = self.get_xy()
        
        x1 = x[:-1]
        y1 = y[:-1]

        x2 = x[1:]
        y2 = y[1:]

        signed_area = np.add.reduce(x1*y2 - x2*y1)

        return abs(signed_area) * 0.5

    def get_centroid(self):
        x, y = self.get_xy()
        x = x[:-1]
        y = y[:-1]
        x_center = x.sum() / x.size
        y_center = y.sum() / y.size
        return np.array([x_center, y_center])

    def get_edge_centroid(self):
        p1 = self.contour[:-1]
        p2 = self.contour[1:]

        edges = p1 - p2
        norms = np.sqrt(np.sum(edges**2, axis=-1))

        edge_centers = (p1 + p2) / 2.0
        totals = norms.dot(edge_centers)
        center = totals / norms.sum()

        # flip x and y
        return center[::-1]

    def get_center(self):
        return self.get_edge_centroid()

    def get_predicted_extents(self, pixscale):
        return get_extents_from_center(self.get_center(), pixscale)
        
    def get_predicted_borders(self, pixscale):
        return get_borders_from_extents(self.get_predicted_extents(pixscale))

    def get_bounding_box(self):
        x, y = self.get_xy()
        return min(x), min(y), max(x), max(y)

    def get_bounding_box_corners(self):
        xmin, ymin, xmax, ymax = self.get_bounding_box()
        bottom_left = np.array([xmin, ymin])
        top_right   = np.array([xmax, ymax])
        return bottom_left, top_right

    def get_bounding_box_center(self):
        bottom_left, top_right = self.get_bounding_box_corners()
        total = bottom_left + top_right
        return total / 2.0

    def get_bounding_box_area(self):
        bottom_left, top_right = self.get_bounding_box_corners()
        diagonal = top_right - bottom_left
        return diagonal.dot(diagonal) / 2.0

    def get_bounding_box_borders(self):
        xmin, ymin, xmax, ymax = self.get_bounding_box()
        points = [(xmin, ymin),
                  (xmin, ymax),
                  (xmax, ymax),
                  (xmax, ymin),
                  (xmin, ymin)]
        return points

    def get_contour(self):
        return self.contour

    def get_contour_level(self):
        return self.contour_level

    def bounding_box_contains(self, point):
        xmin, ymin, xmax, ymax = self.get_bounding_box()
        xcrd, ycrd = point

        if not (xmin <= xcrd < xmax):
            return False

        if not (ymin <= ycrd < ymax):
            return False

        return True

def plot_box_measurement(mosaic_box, acq_box):
    import matplotlib.pyplot as plt
    plt.ion()

    # when in debug mode, render into separate windows
    if is_debug_mode():
        fig = plt.figure(figsize=(8, 4))
    else:
        fig = plt.gcf()
        plt.clf()

    data = mosaic_box.get_data()
    z1, z2 = get_zscale(mosaic_box, acq_box)
    clippeddata = np.clip(data, z1, z2 / 2.0)

    plt.subplot(131)
    plt.imshow(clippeddata, cmap=plt.cm.gray)

    # flip the y axis
    ax = fig.gca()
    ax.invert_yaxis()

    x_center, y_center = acq_box.get_centroid()
    plt.plot([x_center], [y_center], 'ro')

    x_edge_center, y_edge_center = acq_box.get_edge_centroid()
    plt.plot([x_edge_center], [y_edge_center], 'bo')

    msg = """vertex centroid : %.1f, %.1f
    edge centroid : %.1f, %.1f
    contour : %.1f
    measured area : %.2f arc^2""" % (x_center, y_center,
                                     x_edge_center, y_edge_center,
                                     acq_box.get_contour_level(),
                                     acq_box.get_area() * mosaic_box.unbinned_pixel_scale()**2)
    
    plt.text(0.0, 0.0,
             msg,
             fontsize=16,
             horizontalalignment='left',
             verticalalignment='bottom',
             color='white',
             transform=ax.transAxes)
    
    contour = acq_box.get_contour()
    x = contour[:, 1]
    y = contour[:, 0]
   
    plt.plot(x, y, linewidth=1)

    plt.axis('image')
    plt.xticks([])
    plt.yticks([])


    plt.subplot(132)
    plt.imshow(clippeddata, cmap=plt.cm.gray)

    # flip the y axis
    ax = fig.gca()
    ax.invert_yaxis()

    x_center, y_center = acq_box.get_edge_centroid()
    plt.plot([x_center], [y_center], 'bo')

    msg = "predicted box" 
   
    plt.text(0.0, 0.0,
             msg,
             fontsize=16,
             horizontalalignment='left',
             verticalalignment='bottom',
             color='white',
             transform=ax.transAxes)
    
    contour = acq_box.get_predicted_borders(mosaic_box.unbinned_pixel_scale())
    x = contour[:, 0]
    y = contour[:, 1]
   
    plt.plot(x, y, linewidth=1)

    plt.axis('image')
    plt.xticks([])
    plt.yticks([])



    plt.subplot(133)
    plt.imshow(clippeddata, cmap=plt.cm.gray)

    # flip the y axis
    ax = fig.gca()
    ax.invert_yaxis()

    x_center, y_center = acq_box.get_bounding_box_center()

    msg = "bounding box center : %.1f, %.1f\nbounding box area : %.2f arc^2" % (x_center,
                                                                                y_center,
                                                                                acq_box.get_bounding_box_area() * mosaic_box.unbinned_pixel_scale()**2)
    plt.text(0.95, 0.05,
             msg,
             fontsize=16,
             horizontalalignment='right',
             verticalalignment='bottom',
             color='white',
             transform=ax.transAxes)

    plt.plot([x_center], [y_center], 'ro')

    xmin, ymin, xmax, ymax = acq_box.get_bounding_box()
    lines = [zip([xmin, ymin], [xmin, ymax]),
             zip([xmin, ymax], [xmax, ymax]),
             zip([xmax, ymax], [xmax, ymin]),
             zip([xmax, ymin], [xmin, ymin])]

    for l in lines:
        plt.plot(*l, linewidth=1)

    plt.axis('image')
    plt.xticks([])
    plt.yticks([])

    plt.draw()
    return fig

def compare_boxes(b1, b2):
    return cmp(b1.get_area(), b2.get_area())

def measure_box(scidata, contour_level):
    contours = measure.find_contours(scidata, contour_level)
    contours.sort(cmp=compare_contours, reverse=True)

    if len(contours) == 0:
        raise MeasurementFailed("No contours found at the level of %f, try somewhere else without a bright object" % contour_level)

    boxes = []
    for contour in contours:
        if not is_closed_contour(contour):
            continue

        boxes.append(AcquisitionBox(contour, contour_level))

    # find the box with the largest area at this contouring
    boxes.sort(cmp=compare_boxes, reverse=True)
    return boxes[0]

def get_box_error_func(actual_area, pixscale):
    def error_func(box):
        boxarea = box.get_area() * pixscale**2
        diff = actual_area - boxarea
        return diff
    return error_func

def is_point(obj):
    try:
        x, y = obj
        return True
    except TypeError:
        return False

def is_set_of_points(obj):
    try:
        for x, y in obj:
            return True
    except TypeError:
        return False

def get_zscale(mosaic_box, acq_box):
    data = mosaic_box.get_data()
    bbox = acq_box.get_bounding_box()
    xmin, ymin, xmax, ymax = bbox

    indices = np.ndarray(data.shape, dtype=bool)
    indices[ymin:ymax,xmin:xmax] = True
    outside = data[~indices]
    
    mode = scipy.stats.mstats.mode(outside, axis=None)[0] # the background makes a good z1
    
    # look at the inside of the acquisition box to determine z2
    insidedata = data[ymin:ymax,xmin:xmax]

    z1 = mode - outside.std()
    z2 = insidedata.mean()
    return float(z1), float(z2)

def assert_integral_offsets(offset):
    for crd in offset:
        icrd = float(int(crd))
        diff = crd - icrd
        assert 0.0 == diff, "offset not integral %r" % offset

class MosaicTileSelection(object):
    def __init__(self, selection, origin, mosaic_box):
        self.selection = selection
        self.origin = origin
        self.mosaic_box = mosaic_box

        assert_integral_offsets(origin)
        assert_integral_offsets(self.mosaic_box.get_mosaic_offset())
        assert_integral_offsets(self.mosaic_box.get_detector_offset())

    def get_fit_center(self):
        # subtract one since get_selection_peak already converts to 1-based indicing a little too early
        return self.selection.get_center() + self.origin - 1.0

    def get_mosaic_center(self):
        return self.get_fit_center() + self.mosaic_box.get_mosaic_offset()

    def get_detector_center(self):
        center = self.get_fit_center() * self.mosaic_box.detector_binning()
        return center + self.mosaic_box.get_detector_offset()

    def get_fwhm(self):
        return self.selection.get_fwhm()

class MosaicTileExactSelection(object):
    def __init__(self, point, mosaic_box):
        self.mosaic_box = mosaic_box
        self.mosaic_location = point
        self.tile_location = point - mosaic_box.get_mosaic_offset()

    def get_mosaic_center(self):
        return self.tile_location + self.mosaic_box.get_mosaic_offset()

    def get_detector_center(self):
        return self.tile_location + self.mosaic_box.get_detector_offset()

    def get_mosaic_location(self):
        return self.mosaic_location

class MosaicTileExactBox(object):
    def __init__(self, point1, point2):
        first_x, first_y = point1
        next_x, next_y   = point2
        
        self.xmin = min(first_x, next_x)
        self.ymin = min(first_y, next_y)
        self.xmax = max(first_x, next_x)
        self.ymax = max(first_y, next_y)

    def get_center(self):
        return np.array([(self.xmin + self.xmax) / 2.0,
                         (self.ymin + self.ymax) / 2.0])
    
    def get_predicted_extents(self, pixscale):
        return self.xmin, self.ymin, self.xmax, self.ymax

    def get_predicted_borders(self, pixscale):
        return get_borders_from_extents(self.get_predicted_extents(pixscale))

class MosaicTileMDFPredictedBox(object):
    def __init__(self, center):
        self.center = center
        
    def get_center(self):
        return self.center
    
    def get_predicted_extents(self, pixscale):
        return get_extents_from_center(self.get_center(), pixscale)

    def get_predicted_borders(self, pixscale):
        return get_borders_from_extents(self.get_predicted_extents(pixscale))

class MosaicTile(object):
    def __init__(self, mosaic_box, acq_box):
        self.mosaic_box = mosaic_box
        self.acq_box = acq_box
        self.discarded = False

    def adjust_points(self, points, offset):
        adjusted = []
        for pnt in points:
            adjusted.append(pnt + offset)

        return np.array(adjusted)

    def get_mosaic_center(self):
        return self.acq_box.get_center() + self.mosaic_box.get_mosaic_offset()

    def get_mosaic_offset(self):
        return self.mosaic_box.get_mosaic_offset()

    def get_detector_center(self):
        center = self.acq_box.get_center() * self.mosaic_box.detector_binning()
        return center + self.mosaic_box.get_detector_offset()

    def contains_cursor_position(self, point):
        return self.mosaic_box.contains_cursor_position(point)

    def get_zscale(self):
        return get_zscale(self.mosaic_box, self.acq_box)

    def get_mosaic_predicted_borders(self):
        points = self.acq_box.get_predicted_borders(self.mosaic_box.binned_pixel_scale())
        return self.adjust_points(points, self.mosaic_box.get_mosaic_offset()) # account for mosaic offset

    def get_fitted_star_selection(self, point, display_plot, cursor_listener=None, radial=False):
        # uses the entire tile of data instead of just the data inside the detected box
        point = point - self.mosaic_box.get_mosaic_offset()
        selection = get_selection_peak(self.mosaic_box.get_data(),
                                       point,
                                       self.mosaic_box.unbinned_pixel_scale(),
                                       verbose=display_plot,
                                       cursor_listener=cursor_listener,
                                       radial=radial)

        # don't need to specify an origin since we're not doing any subsetting in this function
        origin = np.array([0, 0])
        return MosaicTileSelection(selection, origin, self.mosaic_box)

    def automatically_detect_star_in_tile(self):
        center = self.get_mosaic_center()
        debug("...searching for a star in the tile starting at", center)
        display_plot = is_debug_mode()
        selection = self.get_fitted_star_selection(center, display_plot)
        self.handle_select_star(selection)

    def automatically_detect_star_in_box(self):
        pixscale = self.mosaic_box.unbinned_pixel_scale()
        xmin, ymin, xmax, ymax = self.acq_box.get_predicted_extents(pixscale)

        debug("...box extents = %r, %r, %r, %r" % (xmin, ymin, xmax, ymax))

        # make sure the limits are integral and we only get pixels
        # inside the box to not skew the function fitting
        xmin = int(round(xmin + 1.0))
        ymin = int(round(ymin + 1.0))
        xmax = int(round(xmax - 1.0))
        ymax = int(round(ymax - 1.0))
        
        data = self.mosaic_box.get_data()
        boxdata = data[ymin:ymax,xmin:xmax]

        origin = np.array([xmin, ymin])
        debug("...origin of data to be fit = %r" % origin)
        boxcenter = self.acq_box.get_center() - origin

        selection = get_selection_peak(boxdata, 
                                       boxcenter,
                                       pixscale,       
                                       verbose=is_debug_mode(),
                                       pixel_buffer=max(*boxdata.shape)) # force the entire box of data to be used

        self.handle_select_star(MosaicTileSelection(selection, origin, self.mosaic_box))

    ####################
    # specifying the box
    ####################
    def draw_box(self):
        points = self.get_mosaic_predicted_borders()
        border = ui.polyline(points=points, color=ui.MAGENTA)

        point  = self.get_mosaic_center()
        debug("...drawing box with center at %r" % point)

        center = ui.marker(point, color=ui.MAGENTA)
        
        self.box_overlay = [border, center]

    def undo_previous_box(self):
        for overlay in getattr(self, "box_overlay", []):
            ui.undo(overlay)
        self.box_overlay = []

    def is_midway_through_exact_acquisition_box(self):
        return getattr(self, "first_corner", None) is not None

    def undo_midway_through_exact_acquisition_box(self):
        if self.is_midway_through_exact_acquisition_box():
            self.undo_previous_box()
            self.draw_box()
            self.first_corner = None

    def handle_exact_box(self, point):
        self.undo_box_discard()
        self.undo_previous_star_selection()
        self.undo_previous_box()

        if getattr(self, "first_corner", None) is None:
            overlay = ui.marker(point, color=ui.MAGENTA)
            self.first_corner = point - self.mosaic_box.get_mosaic_offset()
            self.box_overlay = [overlay]
        else:
            self.acq_box = MosaicTileExactBox(self.first_corner, point - self.mosaic_box.get_mosaic_offset())
            self.first_corner = None
            self.draw_box()
            self.automatically_detect_star_in_box()

    ##############################
    # specifying the star location
    ##############################
    def undo_previous_star_selection(self):
        for overlay in getattr(self, "star_selection_overlay", []):
            ui.undo(overlay)
        self.star_selection_overlay = []

    def handle_select_star(self, selection):
        self.undo_midway_through_exact_acquisition_box()
        self.undo_box_discard()
        self.undo_previous_star_selection()

        self.star_selection = selection
        radius = (ACQUISITION_BOX_SIZE / self.mosaic_box.unbinned_pixel_scale()) / 4.0
        self.star_selection_overlay = mark_selection(self.star_selection.get_mosaic_center(), radius)

    def handle_exact_star(self, point):
        self.undo_midway_through_exact_acquisition_box()
        self.undo_box_discard()
        self.undo_previous_star_selection()

        self.star_selection = MosaicTileExactSelection(point, self.mosaic_box)
        overlay = ui.marker(self.star_selection.get_mosaic_center(), color=ui.GREEN)
        self.star_selection_overlay = [overlay]

    def get_star_selection(self):
        return self.star_selection

    ###################
    # rejecting the box
    ###################
    def undo_box_discard(self):
        if self.is_discarded():
            for overlay in getattr(self, "box_discard_overlay"):
                ui.undo(overlay)
            self.box_discard_overlay = []
            
            self.draw_box()

            # draw either an exact or fitted star location
            if isinstance(self.star_selection, MosaicTileExactSelection):
                self.handle_exact_star(self.star_selection.get_mosaic_location())
            else:
                assert isinstance(self.star_selection, MosaicTileSelection)
                self.handle_select_star(self.star_selection)            

    def handle_box_discard(self, discard):
        """
        discard == True  => get rid of box
        discard == False => reset box to previous
        """
        self.undo_midway_through_exact_acquisition_box()
        self.undo_previous_star_selection()
        self.undo_previous_box()
        self.undo_box_discard()

        if discard:
            tile_size    = self.mosaic_box.get_size()
            bottom_left  = self.mosaic_box.get_mosaic_offset() - 1.0
            bottom_right = bottom_left + np.array([tile_size, 0])
            top_left     = bottom_left + np.array([0, tile_size])
            top_right    = bottom_left + np.array([tile_size, tile_size])

            backslash = ui.polyline(points=(bottom_left, top_right), color=ui.RED)
            fwdslash  = ui.polyline(points=(top_left, bottom_right), color=ui.RED)
            self.box_discard_overlay = [fwdslash, backslash]

    def is_discarded(self):
        return bool(getattr(self, "box_discard_overlay", []))


def find_optimal_box(box, display=False):
    actual_area = ACQUISITION_BOX_SIZE * ACQUISITION_BOX_SIZE # arcseconds^2
    
    scidata = box.get_data()

    best_box = find_optimal(scidata,
                            measure_box,
                            get_box_error_func(actual_area, box.unbinned_pixel_scale()))

    if display:
        plotting.put_plot(plot_box_measurement, (box, best_box))
    
    return MosaicTile(box, best_box)

class InteractiveMosaic(object):
    def __init__(self, ad, verbose):
        self.ad = ad
        self.verbose = verbose
        if ad.has_mos_mask():
            debug("...detecting box boundaries...")
        else:
            debug("...using box boundaries from MDF file...")

        self.tiles = []
        for box in ad.get_mos_boxes():
            if self.has_mos_mask():
                tile = find_optimal_box(box, verbose)
            else:
                tile_center = box.get_mdf_tile_center()
                acq_box = MosaicTileMDFPredictedBox(tile_center)
                tile = MosaicTile(box, acq_box)
            self.tiles.append(tile)

    def has_mos_mask(self):
        return self.ad.has_mos_mask()

    def get_zscale(self):
        if not self.has_mos_mask():
            return None, None
        
        z1 = Ellipsis
        z2 = None
        for tile in self.tiles:
            t_z1, t_z2 = tile.get_zscale()
            z1 = min(t_z1, z1)
            z2 = max(t_z2, z2)
        return z1, z2

    def get_tiles(self):
        return self.tiles

    def find_tile(self, point):
        for tile in self.get_tiles():
            if tile.contains_cursor_position(point):
                return tile

        raise ValueError("Unable to find a tile that contains the point %r" % point)

    def get_verbose(self):
        return self.verbose
    
class MOSBoxCursorListener(ui.CursorListener):
    def __init__(self, interactive_mosaic):
        self.interactive_mosaic = interactive_mosaic
        for tile in self.interactive_mosaic.get_tiles():
            tile.draw_box()

            if self.interactive_mosaic.has_mos_mask():
                tile.automatically_detect_star_in_box()
            else:
                tile.automatically_detect_star_in_tile()

    def get_offset(self):
        return self.current_tile.get_mosaic_offset()

    def get_tiles(self):
        for tile in self.interactive_mosaic.get_tiles():
            if not tile.is_discarded():
                yield tile

    def get_box_centers(self):
        for tile in self.get_tiles():
            yield tile.get_detector_center()

    def get_star_centers(self):
        for tile in self.get_tiles():
            selection = tile.get_star_selection()
            yield selection.get_detector_center()
        
    def handle_plot_star(self, pos):
        point = pos.get_position()
        tile = self.interactive_mosaic.find_tile(point)
        self.current_tile = tile

        display_plot = self.interactive_mosaic.get_verbose()
        if pos.get_keystroke() in ['e', 'r', 'a']:
            display_plot = True

        radial = False
        if pos.get_keystroke() in ['r']:
            radial = True
            
        return tile, tile.get_fitted_star_selection(point, display_plot, cursor_listener=self, radial=radial)

    def handle_exact_star(self, pos):
        point = pos.get_position()
        tile = self.interactive_mosaic.find_tile(point)
        tile.handle_exact_star(point)

    def handle_select_star(self, pos):
        tile, selection = self.handle_plot_star(pos)
        tile.handle_select_star(selection)

    def handle_box_discard(self, pos):
        point = pos.get_position()
        tile = self.interactive_mosaic.find_tile(point)
        tile.handle_box_discard(not tile.is_discarded())

    def handle_exact_box(self, pos):
        point = pos.get_position()
        tile = self.interactive_mosaic.find_tile(point)
        tile.handle_exact_box(point)

    def send_cursor_position(self, pos):
        try:
            key = pos.get_keystroke()
            if   key == "x":
                self.handle_exact_star(pos)
            elif key in ["a", "r"]:
                self.handle_select_star(pos)
            elif key == "e":
                self.handle_plot_star(pos)
            elif key == "b" and self.interactive_mosaic.has_mos_mask():
                self.handle_exact_box(pos)
            elif key == "d":
                self.handle_box_discard(pos)
            else:
                warning("Unrecognized keystroke '%s', doing nothing..." % key)
        except ValueError, err:
            warning("Unable to process keystroke: %s" % err.message)

    def send_close(self):
        super(self.__class__, self).send_close()

        for tile in self.interactive_mosaic.get_tiles():
            if tile.is_midway_through_exact_acquisition_box():
                msg = "Tile was midway through an exact box specification, that is, only one box corner was given"
                error(msg)
                raise RuntimeError(msg)
