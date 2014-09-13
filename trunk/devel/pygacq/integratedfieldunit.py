from __future__ import print_function

import sys
import math
import copy
import numpy as np
import pyfits as pf

from collections import defaultdict
from astrodata import AstroData, ConfigSpace
from skimage import measure

import plotting
import userinterface as ui

from acquisitionimage import AcquisitionImage
from multiobjectspectroscopy import AcquisitionBox, is_closed_contour
from selection import SelectionCursorListener
from util import get_window, clamp_to_size

from gacqlogging import *

FIBER_WIDTH = 5.0 # pixels
FIBER_RADIUS = FIBER_WIDTH / 2.0
# divide area by 2 to avoid bleeding into neighbor fibers
FIBER_AREA = (math.pi * (FIBER_RADIUS**2)) / 2.0

class IFUImaging(object):
    """
    IFUImaging modes:
    None => not an IFU imaging mode
    IFU-2 => both red and blue
    IFU-R => red IFU
    IFU-B => blue IFU
    """

    MODES = ["IFU-2", "IFU-R", "IFU-B"]
    INSTRUMENT_OFFSETS = {
        ########### IFU offset,            Sky offset (both in arcseconds in both dimensions)
        "GMOS-N" : [( 29.65811, -1.22136), (-30.07177, -1.92131)],
        "GMOS-S" : [(-29.95351,  2.26694), ( 29.47652,  2.01699)]
        }

    SIZES = {
        ########## IFU size,       Sky size (both in arcseconds in both dimensions)
        "IFU-2" : [(7.0,     5.0), (3.5,     5.0)],
        "IFU-R" : [(7.0/2.0, 5.0), (3.5/2.0, 5.0)],
        "IFU-B" : [(7.0/2.0, 5.0), (3.5/2.0, 5.0)]
        }

    MODE_OFFSETS = {
        ########## IFU offset, Sky offset (both in arcseconds in both dimensions)
        "IFU-2" : [(0.0, 0.0), (0.0,   0.0)]
        }
    for mode in ["IFU-R", "IFU-B"]:
        (xifu, yifu), (xsky, ysky) = SIZES[mode]
        ##################### IFU offset, Sky offset (both in arcseconds in both dimensions)
        MODE_OFFSETS[mode] = [(xifu / 2.0, 0.0), (xsky / 2.0, 0.0)]

    def __init__(self, mode=None, acqimage=None):
        self.mode = mode
        if mode is None:
            return

        if mode not in self.MODES:
            raise ValueError("IFU imaging mode '%s' not supported, choose from %r" % (mode, self.MODES))

        self.acqimage = acqimage
        if not self.acqimage.is_gmos():
            raise ValueError("Only GMOS IFUs do imaging for acquisition")

    def get_mode(self):
        return self.mode

    def is_ifu_imaging(self):
        return self.mode is not None

    def get_instrument_offset(self):
        for offset in self.INSTRUMENT_OFFSETS[self.acqimage.instrument()]:
            yield np.array(offset) / self.acqimage.binned_pixel_scale()

    def get_mode_offset(self):
        for offset in self.MODE_OFFSETS[self.mode]:
            yield np.array(offset) / self.acqimage.binned_pixel_scale()

    def get_sky_center(self):
        """ Return the sky center in pixel coordinates. """
        field_center = self.acqimage.get_binned_data_center()
        debug("...field_center =", field_center)

        ifu_offset, sky_offset = self.get_instrument_offset()
        mode_ifu_offset, mode_sky_offset = self.get_mode_offset()

        sky_center = field_center + sky_offset + mode_sky_offset
        return sky_center

    def get_mode_size(self):
        for offset in self.SIZES[self.mode]:
            yield np.array(offset) / self.acqimage.binned_pixel_scale()

    def get_borders(self, center, sizes):
        """ Return a set of pixel coordinates for where to draw the box borders. """
        xdiff, ydiff = (sizes / 2.0)

        bottom_left = center + np.array((-xdiff, -ydiff))
        top_left    = center + np.array((-xdiff,  ydiff))
        bottom_rght = center + np.array(( xdiff, -ydiff))
        top_rght    = center + np.array(( xdiff,  ydiff))
        return [bottom_left,
                top_left,
                top_rght,
                bottom_rght,
                bottom_left]

    def get_sky_borders(self):
        center = self.get_sky_center()
        ifu_size, sky_size = self.get_mode_size()
        return self.get_borders(center, sky_size)

    def get_text_location(self, center, sizes):
        xdiff, ydiff = (sizes / 2.0)
        top_left    = center + np.array((-xdiff + 2,  ydiff + 5))
        return top_left

    def get_sky_text_location(self):
        center = self.get_sky_center()
        ifu_size, sky_size = self.get_mode_size()
        return self.get_text_location(center, sky_size)

    def get_ifu_center(self):
        """ Return the IFU center in pixel coordinates. """
        field_center = self.acqimage.get_binned_data_center()

        ifu_offset, sky_offset = self.get_instrument_offset()
        mode_ifu_offset, mode_sky_offset = self.get_mode_offset()

        ifu_center = field_center + ifu_offset + mode_ifu_offset
        return ifu_center

    def get_ifu_arcsec_size(self):
        return self.SIZES[self.mode][0]
      
    def get_ifu_borders(self):
        center = self.get_ifu_center()
        ifu_size, sky_size = self.get_mode_size()
        return self.get_borders(center, ifu_size)

    def get_ifu_text_location(self):
        center = self.get_ifu_center()
        ifu_size, sky_size = self.get_mode_size()
        return self.get_text_location(center, ifu_size)


FRAME_NUMBER = 1
def get_frame_number():
    global FRAME_NUMBER
    FRAME_NUMBER += 1
    return FRAME_NUMBER

class Fiber(AcquisitionBox):
    def __init__(self, detector_offset, origin, contour, contour_level):
        AcquisitionBox.__init__(self, contour, contour_level)
        self.detector_offset = detector_offset
        self.origin = origin

    def get_detector_center(self):
        return self.get_data_center() + self.detector_offset - 1.0

    def get_data_center(self):
        return self.get_center() + self.origin + 1.0

def plot_fiber_measurement(window, fiber, expected):
    import matplotlib.pyplot as plt
    plt.ion()

    # when in debug mode, render into separate windows
    if is_debug_mode():
        fig = plt.figure()
    else:
        fig = plt.gcf()
        plt.clf()

    plt.imshow(window, cmap=plt.cm.gray)

    # flip the y axis
    ax = fig.gca()
    ax.invert_yaxis()

    x_center, y_center = fiber.get_center() 
    plotted_point = plt.plot([x_center], [y_center], "ro")

    contour_level = fiber.get_contour_level()
    area = fiber.get_area()
    msg = """area = %.2f pixels**2
    expected = %.2f pixels**2
    contour_level = %f""" % (area, expected, contour_level)

    plt.text(0.95, 0.05,
             msg,
             fontsize=16,
             horizontalalignment='right',
             verticalalignment='bottom',
             color='white',
             transform=ax.transAxes)

    contour = fiber.get_contour()
    x = contour[:, 1]
    y = contour[:, 0]

    plt.plot(x, y, linewidth=1)

    plt.axis('image')
    plt.xticks([])
    plt.yticks([])

    plt.draw()

    return fig

def find_fiber(fiberstamp, position, display_plot=False):
    fiber_data = fiberstamp.get_data()
    xmin, xmax, ymin, ymax = get_window(fiber_data, position, 10)
    debug("...searching for the fiber in the window = [%i:%i,%i:%i]" % (ymin, ymax, xmin, xmax))

    window = fiber_data[ymin:ymax,xmin:xmax]
    origin = np.array((xmin, ymin))

    cursor_point = position - origin

    debug("...contour must contain the point =", cursor_point)

    median = np.median(window)
    maximum = window.max()

    fibers = []
    detector_offset = fiberstamp.get_detector_offset()
    for contour_level in np.linspace(maximum, median):
        contours = measure.find_contours(window, contour_level)

        for contour in contours:
            if not is_closed_contour(contour):
                continue
            
            fiber = Fiber(detector_offset, origin, contour, contour_level)
            if not fiber.bounding_box_contains(cursor_point):
                continue

            fibers.append(fiber)

    def compare_fibers(fiber1, fiber2):
        diff1 = abs(fiber1.get_area() - FIBER_AREA)
        diff2 = abs(fiber2.get_area() - FIBER_AREA)
        return cmp(diff1, diff2)

    fibers.sort(cmp=compare_fibers)

    if display_plot:
        for fiber in fibers:
            plotting.put_plot(plot_fiber_measurement, (window, fiber, FIBER_AREA))
            break

    if fibers:
        return fibers[0]
    
    return None

class FiberDetectionListener(ui.CursorListener):
    def __init__(self, fiberstamp, verbose=False):
        self.fiberstamp = fiberstamp
        self.verbose = verbose
        self.fiber = None
        self.undo = []

    def get_fiber(self):
        return self.fiber

    def get_fiber_center(self):
        if isinstance(self.fiber, Fiber):
            return self.fiber.get_detector_center()

        if self.fiber is None:
            return self.fiberstamp.get_data_center() + self.fiberstamp.get_detector_offset()
            
        return self.fiber + self.fiberstamp.get_detector_offset()

    def send_cursor_position(self, pos):
        key = pos.get_keystroke()
        if   key == "a":
            self.handle_fiber_selection(pos)
        elif key == "x":
            self.handle_exact_location(pos)
        else:
            warning("Unrecognized keystroke '%s', doing nothing..." % key)

    def undo_previous_selection(self):
        for overlay in self.undo:
            ui.undo(overlay)
        self.undo = []
        
    def handle_exact_location(self, pos):
        self.undo_previous_selection()
        overlay = ui.marker(pos.get_position(), color=ui.GREEN)
        self.undo = [overlay]
        self.fiber = pos.get_position()

    def handle_fiber_selection(self, pos):
        self.undo_previous_selection()

        fiber = find_fiber(self.fiberstamp, pos.get_position() - 1.0, display_plot=self.verbose)

        # failed to detect a fiber using contours, falling back on exact location
        if fiber is None:
            self.handle_exact_location(pos)
            return

        self.fiber = fiber
        center = fiber.get_data_center() 
        
        overlay = ui.circle(center, radius=FIBER_RADIUS + 1.0, color=ui.GREEN)
        self.undo = [overlay]

class FiberStamp(object):
    def __init__(self, detsec, data, origin):
        self.detsec = detsec
        self.data = data
        self.origin = np.array(origin)

    def get_data(self):
        return self.data

    def get_detector_offset(self):
        return self.origin + self.detsec.get_detector_offset()

    def get_zscale(self):
        z1 = self.data.min()
        z2 = self.data.mean()
        return z1, z2

    def get_data_center(self):
        return np.array(list(reversed(self.data.shape))) / 2.0

class ReconstructedImage(object):
    def __init__(self, ifudata, skydata, pixscale):
        self.ifudata = ifudata
        self.skydata = skydata
        self.pixscale = pixscale

    def get_science_data(self):
        return self.ifudata

    def get_data_center(self):
        return np.array(list(reversed(self.ifudata.shape))) / 2.0

    def get_ifu_and_sky_data(self):
        return np.hstack([self.skydata, self.ifudata])

    def get_sky_offset(self):
        return np.array((self.skydata.shape[1], 0))

    def get_sky_line(self):
        skyoffset = self.get_sky_offset() + np.array((0, 1))
        line = [skyoffset, skyoffset + np.array((0, self.skydata.shape[0] - 1))]
        return line

    def binned_pixel_scale(self):
        x_pixscale, y_pixscale = self.pixscale
        return x_pixscale

    def get_pixel_scale(self):
        return self.pixscale
  
class FiberBundle(object):
    def __init__(self, block, fibers):
        self.block = block
        self.fibers = list(fibers)

        def compare_fibers(f1, f2):
            return cmp(f1.field("NO"), f2.field("NO"))

        self.fibers.sort(cmp=compare_fibers)

    def _get_focal_plane_centroid(self):
        accum = np.array((0.0, 0.0))
        for f in self.fibers:
            focal_plane_pos = get_fiber_focal_plane_position(f)
            accum += focal_plane_pos

        return accum / len(self.fibers)

    def get_focal_plane_centroid(self):
        if not hasattr(self, "focal_plane_centroid"):
            self.focal_plane_centroid = self._get_focal_plane_centroid()
        return self.focal_plane_centroid

    def has_cached_fiber_positions(self):
        try:
            return self.fibers[0].field("XDET") != 0.0
        except KeyError:
            return False

    def get_cached_fiber_positions(self):
        for f in self.fibers:
            yield f.field("XDET"), f.field("YDET")

    def get_block(self):
        return self.block

    def _get_first_fiber_idx(self):
        return min(f.field("NO") for f in self.fibers)

    def get_first_fiber_idx(self):
        if not hasattr(self, "first_fiber_idx"):
            self.first_fiber_idx = self._get_first_fiber_idx()
        return self.first_fiber_idx

    def __cmp__(self, rhs):
        return cmp(self.get_first_fiber_idx(), rhs.get_first_fiber_idx())

    def get_num_fibers(self):
        return len(self.fibers)

    def get_first_fiber(self):
        return self.fibers[0]

    def get_top_fiber(self):
        fibers = ((f.field("YDET"), f) for f in self.fibers)
        fiber_y, fiber = max(fibers)
        return fiber

    def get_bottom_fiber(self):
        fibers = ((f.field("YDET"), f) for f in self.fibers)
        fiber_y, fiber = min(fibers)
        return fiber

    def set_first_fiber_center(self, pos):
        self.first_fiber_center = pos

    def get_last_fiber(self):
        return self.fibers[-1]

    def set_last_fiber_center(self, pos):
        self.last_fiber_center = pos

    def _get_fiber_positions(self):
        n_fibers = self.get_num_fibers()

        x_first, y_first = self.first_fiber_center
        x_last, y_last = self.last_fiber_center
        
        xpositions = np.linspace(x_first, x_last, n_fibers)
        ypositions = np.linspace(y_first, y_last, n_fibers)

        for xpos, ypos in zip(xpositions, ypositions):
            yield xpos, ypos

    def get_fiber_positions(self):
        if (getattr(self, "first_fiber_center", None) is None or
            getattr(self, "last_fiber_center", None) is None):
            if self.has_cached_fiber_positions():
                return self.get_cached_fiber_positions()
            
            col = [0.0] * self.get_num_fibers()
            return zip(col, col)

        if not hasattr(self, "fiber_positions"):
            self.fiber_positions = list(self._get_fiber_positions())

        return self.fiber_positions

    def get_fibers(self):
        return self.fibers

def is_fiber_sky(fiber):
    return not fiber.field("BLOCK").split("_")[0].isdigit()

def get_fiber_position(fiber):
    return np.array((fiber.field("XDET"), fiber.field("YDET")))

def get_fiber_focal_plane_position(fiber):
    return np.array((fiber.field("XINST"), fiber.field("YINST")))

def get_fiber_coordinates(fiber):
    # these appear to be flipped in the MDF file
    return fiber.field("YLDIS") - 1, fiber.field("XLDIS") - 1

class FiberArray(object):
    def __init__(self, fibers, stamp_size):
        self.stamp_size = stamp_size
        fibers = list(fibers)
        debug("...number of fibers =", len(fibers))
        bounding_box = self.get_fiber_coordinates_bounding_box(fibers)
        self.xmin, self.ymin, self.xmax, self.ymax = bounding_box
        debug("...fiber box coordinates bounding box =", bounding_box)
        self.xmax += 1
        self.ymax += 1
        self.xdim = self.xmax - self.xmin
        self.ydim = self.ymax - self.ymin
        debug("...fiber box dimensions =", self.xdim, self.ydim)

        row = [None] * self.xdim
        columns = [copy.deepcopy(row) for _ in range(self.ydim)]

        if is_debug_mode():
            self.stamps = columns
            self.samples = copy.deepcopy(columns)
            self.fluxes = np.zeros((self.ydim, self.xdim))

        self.ydim = self.ydim / 2
        if is_debug_mode():
            self.blocks = [[None] * self.xdim for _ in range(self.ydim)]

        self.ydim += 1

        self.shifted_fluxes = np.zeros((self.ydim, self.xdim))

        self.focal_plane_size = self.get_focal_plane_size(fibers)
        
    @staticmethod
    def get_fiber_coordinates_bounding_box(fibers):
        xmin = Ellipsis
        ymin = Ellipsis
        xmax = None
        ymax = None
        
        for fiber in fibers:
            xpos, ypos = get_fiber_coordinates(fiber)
            
            xmin = min(xpos, xmin)
            ymin = min(ypos, ymin)
            
            xmax = max(xpos, xmax)
            ymax = max(ypos, ymax)

        return xmin, ymin, xmax, ymax

    @staticmethod
    def get_focal_plane_size(fibers):
        positions = []
        for fiber in fibers:
            positions.append(get_fiber_focal_plane_position(fiber))

        xpositions, ypositions = zip(*positions)
        xsize = max(xpositions) - min(xpositions)
        ysize = max(ypositions) - min(ypositions)
        return np.array((xsize, ysize))

    def get_pixel_scale(self):
        return self.focal_plane_size / np.array((self.xdim - 1, self.ydim - 1))

    def get_fiber_indices(self, fiber):
        xidx, yidx = get_fiber_coordinates(fiber)
        xidx = xidx - self.xmin
        yidx = yidx - self.ymin
        return xidx, yidx

    def add_fiber_stamp(self, fiber, stamp):
        xidx, yidx = self.get_fiber_indices(fiber)
        self.stamps[yidx][xidx] = stamp
        
    def add_fiber_sample(self, fiber, sample):
        xidx, yidx = self.get_fiber_indices(fiber)
        self.samples[yidx][xidx] = sample

    def add_fiber_flux(self, fiber, flux):
        xidx, yidx = self.get_fiber_indices(fiber)

        if is_debug_mode():
            self.blocks[yidx / 2][xidx] = fiber.field("BLOCK")
            self.fluxes[yidx,xidx] = flux

        if yidx % 2 == 0:
            yidx = yidx / 2
            self.shifted_fluxes[yidx, xidx] += flux
        else:
            yidx = yidx / 2
            half_flux = flux / 2.0
            self.shifted_fluxes[yidx,     xidx] += half_flux
            self.shifted_fluxes[yidx + 1, xidx] += half_flux

    def _build_image(self, blocks):
        zeros = np.zeros((self.stamp_size, self.stamp_size))
        columns = []
        for col in blocks:
            row = []
            for square in col:
                if square is None:
                    row.append(zeros)
                else:
                    row.append(square)
            
            columns.append(np.hstack(row))

        return np.vstack(columns)

    def get_stamp_image(self):
        return self._build_image(self.stamps)

    def get_samples_image(self):
        return self._build_image(self.samples)

    def get_flux_image(self):
        return self.fluxes

    def get_shifted_image(self):
        return self.shifted_fluxes[1:-1,:]

    def write_block_file(self, fname):
        bfile = open(fname, 'w')
        for row in reversed(self.blocks):
            print('\t'.join(row), file=bfile)

class FiberBundleCollection(object):
    def __init__(self, fname, acqimage):
        if fname is None:
            if acqimage.is_type("GMOS_N"):
                fname = ConfigSpace.lookup_path("Gemini/GMOS/MDF/gnifu_slits_mdf.fits")
            elif acqimage.is_type("GMOS_S"):
                fname = ConfigSpace.lookup_path("Gemini/GMOS/MDF/gsifu_slits_mdf.fits")
            else:
                raise ValueError("Only GMOS North and South supported")

        self.ad = AstroData(fname)

        bundle_map = defaultdict(list)
        self.num_fibers = 0
        for record in self.ad.data:
            block = record.field("BLOCK")
            bundle_idx = block.split("_")[0]
            bundle_map[bundle_idx].append(record)
            self.num_fibers += 1

        self.bundles = []
        for block, fibers in bundle_map.items():
            bundle = FiberBundle(block, fibers)
            self.bundles.append(bundle)

        self.bundles.sort()

    def get_total_num_fibers(self):
        return self.num_fibers

    def get_bundles(self, mask=None):
        if mask is None:
            mask = "IFU-2"

        total_fibers = self.get_total_num_fibers()
        if mask in ["IFU-R", "IFU-2"]:
            minidx = 0
        else:
            minidx = total_fibers / 2

        if mask in ["IFU-B", "IFU-2"]:
            maxidx = total_fibers
        else:
            maxidx = total_fibers / 2

        for bundle in self.bundles:
            if minidx < bundle.get_first_fiber_idx() < maxidx:
                yield bundle

    def write_detector_mapping(self, fname):
        detmap = open(fname, 'w')
        
        for side in ["IFU-R", "IFU-B"]:
            print(side, file=detmap)

            sortedbundles = list(self.get_bundles(side))
            def compare_top_fibers(b1, b2):
                ydet1 = b1.get_top_fiber().field("YDET")
                ydet2 = b2.get_top_fiber().field("YDET")
                return cmp(ydet1, ydet2)

            sortedbundles.sort(cmp=compare_top_fibers, reverse=True)

            xpositions = []
            for bundle in sortedbundles:
                top = bundle.get_top_fiber()
                print(top.field("BLOCK"), file=detmap)
                xdet = top.field("XDET")
                xpositions.append(xdet)

                print("...", file=detmap)

                bottom = bundle.get_bottom_fiber()
                print(bottom.field("BLOCK"), file=detmap)
                xdet = bottom.field("XDET")
                xpositions.append(xdet)
                
            print(sum(xpositions) / len(xpositions), file=detmap)
            print("", file=detmap)
            

    def get_fibers(self, mask=None):
        for bundle in self.get_bundles(mask):
            for fiber in bundle.get_fibers():
                yield fiber

    def get_columns(self):
        for col in self.ad.hdulist[-1].columns:
            if col.name == "XDET":
                continue
            elif col.name == "YDET":
                continue
            yield col

    def get_column_names(self):
        for col in self.get_columns():
            yield col.name

    def get_column_header(self):
        colnames = tuple(self.get_column_names())
        return self.get_format_string() % colnames

    def get_base_number_of_columns(self):
        return len(list(self.get_column_names()))

    def get_format_string(self):
        return "%10s" * self.get_base_number_of_columns()

    def format_record(self, rec):
        return self.get_format_string() % tuple(map(str, rec[:self.get_base_number_of_columns()]))

    def get_new_record_format(self):
        return "%10s" * (len(list(self.get_column_names())) + 2)
    
    def get_new_column_header(self):
        colnames = tuple(self.get_column_names()) + ("XDET", "YDET")
        return self.get_new_record_format() % colnames
    
    def format_fiber(self, fiber, pos):
        args = map(str, fiber[:self.get_base_number_of_columns()])
        for p in pos:
            args.append("%.3f" % p)
        return self.get_new_record_format() % tuple(args)

    def get_fiber_positions_columns(self):
        xpositions = []
        ypositions = []

        for bundle in self.get_bundles():
            xpos, ypos = zip(*bundle.get_fiber_positions())
            xpositions.extend(xpos)
            ypositions.extend(ypos)
            
        xcol = pf.Column(name="XDET", format="E", array=np.array(xpositions))
        ycol = pf.Column(name="YDET", format="E", array=np.array(ypositions))
        return xcol, ycol

    def write_new_table(self, fname):
        cols = list(self.get_columns())
        cols.extend(self.get_fiber_positions_columns())
        
        # Create the table HDU
        tablehdu = pf.new_table(cols)
        
        # Create an AstroData object to contain the table
        # and write to disk.
        new_ad = AstroData(tablehdu)
        new_ad.rename_ext('SCI', 1)
        new_ad.write(fname, clobber=True)

class FiberCoordinates(object):
    def __init__(self, acqimage):
        if acqimage.focal_plane_mask() not in ["IFU-R", "IFU-B", "IFU-2"]:
            raise ValueError("Unable to figure out fiber coordinates for the mask '%s'" % acqimage.focal_plane_mask())
        if acqimage.grating().upper() != "MIRROR":
            raise ValueError("Fiber coordinates are only relevant to non-dispersed images")
        self.acqimage = acqimage

        dname = os.path.dirname(__file__)
        if acqimage.is_type("GMOS_N"):
            fname = os.path.join(dname, "GMOS_NORTH_IFU_detector_positions.fits")
        elif acqimage.is_type("GMOS_S"):
            fname = os.path.join(dname, "GMOS_SOUTH_IFU_detector_positions.fits")
        else:
            raise ValueError("Only GMOS north and south supported")

        self.fiber_bundles = FiberBundleCollection(fname, acqimage)

    def get_fiber_stamp(self):
        # get the fiber bundle closest to the center of the focal plane
        accum = np.array((0.0, 0.0))
        ifubundles = []
        for bundle in self.fiber_bundles.get_bundles(self.acqimage.focal_plane_mask()):
            if not bundle.get_block().isdigit(): # skip sky bundles
                continue

            accum += bundle.get_focal_plane_centroid()
            ifubundles.append(bundle)

        center = accum / len(ifubundles)

        debug("...centroid of all the fiber bundles =", center)
            
        def compare_ifu_bundles(b1, b2):
            dist1 = np.linalg.norm(center - b1.get_focal_plane_centroid())
            dist2 = np.linalg.norm(center - b2.get_focal_plane_centroid())
            return cmp(dist1, dist2)
        ifubundles.sort(cmp=compare_ifu_bundles)

        center_bundle = ifubundles[0]
        fiber = center_bundle.get_top_fiber()
        debug("...displaying fiber =", fiber.field("BLOCK"))
        fiber_center = get_fiber_position(fiber)

        self.reference_fiber_center = fiber_center
        
        return self._get_fiber_stamp(fiber_center, 18, 200)

    def _get_fiber_stamp(self, point, xsize, ysize):
        detector_point = np.array(point)
        detsec = self.acqimage.get_full_field_of_view()
        data = detsec.get_data()

        xifu, yifu = detector_point - detsec.get_detector_offset()

        xmin, xmax = clamp_to_size(data, xifu, xsize, 1)
        ymin, ymax = clamp_to_size(data, yifu, ysize, 0)

        stamp = np.zeros((ysize, xsize))
        dstxmin = max(xsize - xmax, 0)
        dstxmax = dstxmin + (xmax - xmin)

        dstymin = max(ysize - ymax, 0)
        dstymax = dstymin + (ymax - ymin)

        stamp[dstymin:dstymax, dstxmin:dstxmax] = data[ymin:ymax,xmin:xmax]

        return FiberStamp(detsec, stamp, (xmin, ymin))

    def get_fibers(self, is_sky=None):
        for fiber in self.fiber_bundles.get_fibers(self.acqimage.focal_plane_mask()):
            if is_sky is None:
                yield fiber
            elif not (is_sky ^ is_fiber_sky(fiber)):
                yield fiber

    def get_reconstructed_image(self, measured_fiber_center):
        debug("...measured fiber center =", measured_fiber_center)
        debug("...reference fiber center =", self.reference_fiber_center)

        fiber_correction = measured_fiber_center - self.reference_fiber_center
        debug("...fiber correction =", fiber_correction)

        fiber_width = 5.68
        stamp_size = 10

        skyarray = FiberArray(self.get_fibers(is_sky=True), stamp_size)
        ifuarray = FiberArray(self.get_fibers(is_sky=False), stamp_size)

        yvalues, xvalues = np.indices((stamp_size, stamp_size))
        def gaussian(fiber_center):
            xcenter, ycenter = fiber_center
            dist = (xcenter - xvalues)**2 + (ycenter - yvalues)**2
            return np.exp(-dist / fiber_width)

        for fiber in self.get_fibers():
            arr = ifuarray
            if is_fiber_sky(fiber):
                arr = skyarray
            
            detector_position = get_fiber_position(fiber) + fiber_correction

            stamp = self._get_fiber_stamp(detector_position, stamp_size, stamp_size)

            stamp_fiber_center = detector_position - stamp.get_detector_offset()
            sampling = gaussian(stamp_fiber_center - 1.0)
            sample = sampling * stamp.get_data()

            flux = np.sum(sample)
            arr.add_fiber_flux(fiber, flux)

            if is_debug_mode():
                arr.add_fiber_stamp(fiber, stamp.get_data())
                arr.add_fiber_sample(fiber, sample)

        if is_debug_mode():
            ui.display(ifuarray.get_stamp_image(), frame=get_frame_number(), zscale=True)
            ui.display(ifuarray.get_samples_image(), frame=get_frame_number(), zscale=True)
            ui.display(ifuarray.get_flux_image(), frame=get_frame_number(), zscale=True)
            ifuarray.write_block_file(self.acqimage.instrument() + "_ifu_mapping.txt")

            ui.display(skyarray.get_stamp_image(), frame=get_frame_number(), zscale=True)
            ui.display(skyarray.get_samples_image(), frame=get_frame_number(), zscale=True)
            ui.display(skyarray.get_flux_image(), frame=get_frame_number(), zscale=True)
            skyarray.write_block_file(self.acqimage.instrument() + "_sky_mapping.txt")

            self.fiber_bundles.write_detector_mapping(self.acqimage.instrument() + "_detector_mapping.txt")

        ifudata = ifuarray.get_shifted_image()
        skydata = skyarray.get_shifted_image()
        return ReconstructedImage(ifudata, skydata, ifuarray.get_pixel_scale())

def get_integrated_field_unit_offsets(acqimage, verbose=False):
    ifucoords = FiberCoordinates(acqimage)

    # display a stamp around the fiber we would like measured
    fiberstamp = ifucoords.get_fiber_stamp()
    z1, z2 = fiberstamp.get_zscale()
    ui.display(fiberstamp.get_data(), frame=1, z1=z1, z2=z2) #zscale=True)
    ui.circle(fiberstamp.get_data_center(), radius=10, color=ui.MAGENTA)

    print("")
    print("   Put cursor on top fiber of lower set (in circle) - press 'a'.")
    print("   When finished, press 'q'.")
    print("")
    print("   If fibers are too faint, just press 'q' for default values.")
    print("   For short exposures the top set of fibers may not be visible.")
    print("")

    
    listener = FiberDetectionListener(fiberstamp, verbose=verbose)
    listener.start()
    
    fiber_center = listener.get_fiber_center()

    # reconstruct the image based upon that measurement
    reconstructed = ifucoords.get_reconstructed_image(fiber_center)

    data = reconstructed.get_ifu_and_sky_data()
    ui.display(data, frame=1, zscale=True)
    ui.polyline(reconstructed.get_sky_line(), color=ui.BLUE)

    print("")
    print("   Point to the object and press 'a' or 'r', the target will then be fitted and marked in DS9.")
    print("   - press 'x' to use an exact location on the acquisition image")
    print("   - press 'e' to see a plot of the profile (clicking in the plot selects an exact location)")
    print("")
    print("   Press 'q' when you're happy with the selection")


    listener = SelectionCursorListener(reconstructed,
                                       verbose=verbose,
                                       pixel_buffer=max(data.shape),
                                       circle_radius=min(data.shape) / 4,
                                       offset=reconstructed.get_sky_offset())
    
    listener.start()

    objcrds = listener.get_object_coords()
    debug("...object coords =", objcrds)
    
    center = reconstructed.get_data_center()
    debug("...center =", center)

    pixscale = reconstructed.get_pixel_scale()
    debug("...pixscale =", pixscale)
    offsets = (center - objcrds) * pixscale

    if acqimage.is_type("GMOS_S"):
        # Check for N&S mode and adjust offset from rec image cen to field cen
        maskname = acqimage.focal_plane_mask()

        nscorr = np.array((0.52, 0.0))
        if maskname == "IFU-NS-B":
            offsets += nscorr
        elif maskname == "IFU-NS-R":
            offsets -= nscorr
        
    return offsets

def main(argv=[__name__]):
    idx = 1
    verbose = False
    if argv[1] == "-v":
        setup_logging(None, True)
        idx = 2
        verbose = True

    acqimage = AcquisitionImage(argv[idx])

    get_integrated_field_unit_offsets(acqimage, verbose)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
