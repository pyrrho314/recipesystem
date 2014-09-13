import os, sys
import math

import numpy as np

from astrodata import AstroData
from astrodata import Lookups

from gacqlogging import *
from util import fits_filename, get_window, is_number
import userinterface as ui

def cache(method):
    def newmethod(self, *args, **kwargs):
        mname = "_cached_" + method.__name__
        if hasattr(self, mname):
            return getattr(self, mname)
        val = method(self, *args, **kwargs)
        setattr(self, mname, val)
        return val
    return newmethod

def _obtain_unbinned_arraygap(adinput):
    """
    This function was copied directly from primitives_GMOS.py, it
    should be refactored at some point in the future! It was tweaked
    to only return unbinned array gaps. 
    
    This function obtains the raw array gap size for the different GMOS
    detectors and returns it after correcting for binning. There are two
    values in the GMOSArrayGaps.py file in the GMOS
    lookup directory, one for unbinned data and one to be used to calculate
    the chip gap when the data are binned.
    """
    
    # Get the dictionary containing the CCD gaps
    all_arraygaps_dict = Lookups.get_lookup_table("Gemini/GMOS/GMOSArrayGaps.py","gmosArrayGaps")
    
    # Obtain the X binning and detector type for the ad input
    detector_type = adinput.phu_get_key_value("DETTYPE")
    
    # Check the read values
    if detector_type is None:
        if hasattr(ad, "exception_info"):
            raise adinput.exception_info
        
    # We're only interested in the unbinned values for gacq
    binning = "unbinned"
        
    # Form the key
    key = (detector_type, binning)
    
    # Obtain the array gap value 
    if key in all_arraygaps_dict:
        arraygap = all_arraygaps_dict[key] 
    else:
        raise Errors.ScienceError("Array gap value not " +
                                  "found for %s" % (detector_type))
    return arraygap
                                                                                                                        

def subtract_overscan(ad):
    ox1, ox2, oy1, oy2 = ad.overscan_section().as_list()
    overscan_data = ad.data[oy1:oy2,ox1:ox2]
    med = np.median(overscan_data)

    dx1, dx2, dy1, dy2 = ad.data_section().as_list()
    science_data = ad.data[dy1:dy2,dx1:dx2]
    science_data -= med

    return science_data
    
def gmultiamp(ad):
    left_science_data = subtract_overscan(ad[3])
    rght_science_data = subtract_overscan(ad[2])
    
    return np.hstack((left_science_data, rght_science_data))

def extract_values(dim):
    return map(int, dim.split(":"))

def extract_dimensions(header_value):
    xdims, ydims = header_value.strip("[]").split(",")
    xmin, xmax = extract_values(xdims)
    ymin, ymax = extract_values(ydims)
    return xmin, xmax, ymin, ymax

def mosaic(ad):
    detsize = ad.phu.header["DETSIZE"]
    xmin, xdim, ymin, ydim = extract_dimensions(detsize)
    
    data = np.zeros((ydim, xdim))

    for ext in ad:
        xmin, xmax, ymin, ymax = ext.detector_section().as_list()
        data[ymin:ymax,xmin:xmax] = subtract_overscan(ext)

    return data

def parse_box_coords(acqimage, mdffile):
    mdffile_ad = AstroData(mdffile)

    mdfxbin = mdffile_ad.phu.header["ODFXBIN"]
    mdfybin = mdffile_ad.phu.header["ODFYBIN"]
       
    debug("...MDF binning = ", mdfxbin, " x ", mdfybin)
    debug("...reading MDF...")

    box_coords = []

    nboxes = 0
    for row in mdffile_ad["MDF"].data:
        # select the alignment boxes, designated by priority 0
        if row["priority"] != "0":
            continue

        nboxes += 1

        if acqimage.is_new_gmosn_ccd():
            # parse out the coordinates of the acquisition boxes in mm and convert to pixels:
            xb = row["slitpos_mx"]
            yb = row["slitpos_my"]
            debug("...box ", nboxes, ": ", xb, yb, "(mm)")

            # GMOS-N E2V DD detector tranformation from Kathy Roth (2012 Jan 19):
            xb = xb * 22.075 + 3113.1
            yb = yb * 22.075 + 2284.3
        else:
            # parse out the (x,y) coordinates of the acquisition boxes in pixels:
            xb = row["x_ccd"]
            yb = row["y_ccd"]
            
            xb = xb * mdfxbin
            yb = yb * mdfybin

        debug("...box ", nboxes, ": ", xb, yb, "(pix)")

        box_coords.append((xb, yb))

    debug("...MDF has", nboxes, "acquisition boxes defined")

    return box_coords

class DetectorSection(object):
    def __init__(self, ext):
        self.ext = ext
        xmin, xmax, ymin, ymax = ext.detector_section().as_list()
        debug("...detector section = %i, %i, %i, %i" % (xmin, xmax, ymin, ymax))

        chip_num = self.get_chip_number()
        debug("...chip number =", chip_num)

        gap_width = self.get_chip_gap()
        debug("...gap width =", gap_width)

        self.xmin = xmin + (gap_width * chip_num)
        self.ymin = ymin
        
        xsize = xmax - xmin
        self.xmax = self.xmin + xsize

        ysize = ymax - ymin
        self.ymax = self.ymin + ysize

    def get_chip_gap(self):
        return _obtain_unbinned_arraygap(self.ext)

    def get_chip_number(self):
        detsize = self.ext.phu_get_key_value("DETSIZE")
        debug("...DETSIZE =", detsize)
        dxmin, dxmax, dymin, dymax = extract_dimensions(detsize)

        ccdsize = self.ext.get_key_value("CCDSIZE")
        debug("...CCDSIZE =", detsize)
        cxmin, cxmax, cymin, cymax = extract_dimensions(ccdsize)

        # check that all the CCDS are the same size
        if dxmax % cxmax != 0 or dymax % cymax != 0:
            raise ValueError("DETSIZE (%s) is not a multiple of the CCDSIZE (%s) for this extension")

        xmin, xmax, ymin, ymax = self.ext.detector_section().as_list()
        chip_num = xmin / cxmax

        # xmin and xmax for this section should be on the same chip
        if chip_num != (xmax - 1) / cxmax:
            raise ValueError("The detector_section %r doesn't all lie on the same chip %r" %
                             ((xmin,   xmax,  ymin,  ymax),
                              (cxmin, cxmax, cymin, cymax)))
        return chip_num       
        
    def get_detector_offset(self):
        return np.array([self.xmin, self.ymin]) + float(self.detector_binning())

    def contains(self, point):
        x, y = np.array(point) - 1.0
        if not (self.xmin <= x < self.xmax):
            return False

        if not (self.ymin <= y < self.ymax):
            return False

        return True

    @cache
    def get_data(self):
        return subtract_overscan(self.ext)

    def get_shape(self):
        return (self.xmax - self.xmin), (self.ymax - self.ymin)

    def unbinned_pixel_scale(self):
        return self.binned_pixel_scale() / self.detector_binning()

    @cache
    def binned_pixel_scale(self):
        return self.ext.pixel_scale()

    @cache
    def detector_binning(self):
        assert int(self.ext.detector_y_bin()) == int(self.ext.detector_x_bin())
        return int(self.ext.detector_y_bin())

class CombineDetectorSections(object):
    def __init__(self, detsecs):
        def sort_by_xmin(dsec1, dsec2):
            return cmp(dsec1.xmin, dsec2.xmin)

        self.sections = list(detsecs)
        self.sections.sort(cmp=sort_by_xmin)

        self.data_sections = []

        prev_chipnum = None
        for detsec in self.sections:
            data = detsec.get_data()

            chipnum = detsec.get_chip_number()
            if prev_chipnum is not None and prev_chipnum != chipnum:
                ydim = data.shape[0]
                xdim = detsec.get_chip_gap() / detsec.detector_binning()
                gap = np.zeros((ydim, xdim))
                self.data_sections.append(gap)
            prev_chipnum = chipnum
                
            self.data_sections.append(data)

    def _get_left(self):
        return self.sections[0]

    def get_detector_offset(self):
        return self._get_left().get_detector_offset()

    def unbinned_pixel_scale(self):
        return self._get_left().unbinned_pixel_scale()

    def binned_pixel_scale(self):
        return self._get_left().binned_pixel_scale()

    def detector_binning(self):
        return self._get_left().detector_binning()

    @cache
    def get_data(self):
        return np.hstack(self.data_sections)

def clamp_to_fullsize(data, coord, fullsize, dimension):
    dimsize = data.shape[dimension]
    assert coord < dimsize, "%r >= %r" % (coord, dimsize)
    assert fullsize <= dimsize

    pixel_buffer = fullsize / 2
    crd = int(round(coord))

    crdmin = crd - pixel_buffer
    if crdmin < 0:
        crdmin = 0
        crdmax = fullsize
        return crdmin, crdmax

    crdmax = crd + pixel_buffer
    if crdmax > dimsize:
        crdmin = dimsize - fullsize
        crdmax = dimsize
        return crdmin, crdmax

    return crdmin, crdmax

ACQUISITION_BOX_SIZE = 2 # arcseconds
class DetectorSectionFinder(object):
    def __init__(self, acqimage):
        self.acqimage = acqimage
        self.sections = []
        for ext in acqimage.get_extensions():
            detsec = DetectorSection(ext)
            self.sections.append(detsec)

    def get_box_size(self):
        sizes = []
        for detsec in self.sections:
            sizes.extend(detsec.get_shape())
        
        minsize = min(sizes)
        unbinned_pixscale = self.unbinned_pixel_scale()
        if minsize < ACQUISITION_BOX_SIZE / unbinned_pixscale:
            raise ValueError("Detector section smaller than an acquisition box!")

        maximum_box_size = 100 # pixels
        if not self.acqimage.has_mos_mask(): # make the tile size larger whenever there isn't a maskx
            maximum_box_size = 200

        return min(maximum_box_size, minsize)

    def find_detector_section(self, point):
        diff = (self.get_box_size() / 2.0) * self.acqimage.detector_y_bin()
        xmin, ymin = point - np.array([diff, diff])
        xmax, ymax = point + np.array([diff, diff])

        corners = [(xmin, ymin),
                   (xmax, ymax)]
        
        sections = []
        for detsec in self.sections:
            for pnt in corners:
                if detsec.contains(pnt):
                    sections.append(detsec)
                    break
                
        if not sections:
            raise ValueError("No detector section found that contains the point %r" % point)

        if len(sections) == 1:
            return sections[0]

        debug("...multiple detector sections found, combining them...")

        return CombineDetectorSections(sections)

    def unbinned_pixel_scale(self):
        return self.acqimage.unbinned_pixel_scale()

    def binned_pixel_scale(self):
        return self.unbinned_pixel_scale() * self.acqimage.detector_y_bin()

    @cache
    def get_full_field_of_view(self):
        return CombineDetectorSections(self.sections)

class Box(object):
    def __init__(self, detsec, detector_point, box_size, mosaic_position):
        self.detsec = detsec
        self.box_size = box_size
        self.mosaic_position = np.array(mosaic_position) # column_num, row_num

        debug("...box center on detector = %f, %f" % tuple(detector_point))

        point = np.array(detector_point)

        # adjust for where the detsec is
        debug("...origin of this extension  = %f, %f" % tuple(detsec.get_detector_offset()))
        point = point - detsec.get_detector_offset()

        debug("...box center in this extension = %f, %f" % tuple(point))

        point = point / detsec.detector_binning()
        debug("...binned box center in this extension = %f, %f" % tuple(point))

        # adjust to zero based offsets
        point -= 1.0

        # grab as much of the box as possible so the output is always box_size by box_size
        data = detsec.get_data()

        debug("...dimensions of this detector section = %i, %i" % data.shape)

        ydim, xdim = data.shape
        if ydim < box_size:
            raise ValueError("Box of data too small in the y-dimension")
        if xdim < box_size:
            raise ValueError("Box of data too small in the x-dimension")

        xcrd, ycrd = point
        xmin, xmax = clamp_to_fullsize(data, xcrd, box_size, 1)
        assert xmax - xmin == box_size

        ymin, ymax = clamp_to_fullsize(data, ycrd, box_size, 0)
        assert ymax - ymin == box_size

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        window = (ymin, ymax, xmin, xmax)
        debug("...cutting out section around acquisition hole =", window)
        
        self.data = data[ymin:ymax,xmin:xmax]
        self.mdf_tile_center = point - np.array([xmin, ymin])

        debug("...binned tile center =", self.mdf_tile_center)

    def get_mdf_tile_center(self):
        """ Returns where the MDF file predicted this center would be
        in this box. The point is relative to the origin of this box
        in one-based indicing.
        """
        return self.mdf_tile_center + 1.0

    def get_mdf_detector_center(self):
        center = self.get_mdf_tile_center() * self.detector_binning()
        return center + self.get_detector_offset()
        
    def get_data(self):
        return self.data

    def get_detector_offset(self):
        origin = np.array((self.xmin, self.ymin)) * self.detector_binning()
        return self.detsec.get_detector_offset() + origin

    def get_mosaic_offset(self):
        return (self.mosaic_position * self.box_size) + 1.0

    def contains_cursor_position(self, point):
        xmin, ymin = self.get_mosaic_offset()
        xmax = xmin + self.box_size
        ymax = ymin + self.box_size

        xcrd, ycrd = point

        if not (xmin <= xcrd < xmax):
            return False

        if not (ymin <= ycrd < ymax):
            return False

        return True

    def get_size(self):
        assert (self.xmax - self.xmin) == (self.ymax - self.ymin)
        return self.xmax - self.xmin

    def unbinned_pixel_scale(self):
        return self.detsec.unbinned_pixel_scale()

    def binned_pixel_scale(self):
        return self.detsec.binned_pixel_scale()

    def detector_binning(self):
        return self.detsec.detector_binning()
    
class BoxMosaic(object):
    def __init__(self, detsec_finder, points):
        self.detsec_finder = detsec_finder

        # every square will be box_size by box_size
        self.box_size = self.detsec_finder.get_box_size()
        debug("...dimension of each box =", self.box_size)

        # tile all the boxes into a square
        numboxes = len(points)
        debug("...numboxes =", numboxes)
        
        self.numrows = int(math.sqrt(numboxes - 1)) + 1
        debug("...numrows =", self.numrows)

        self.boxes = []

        columns = []
        for imin in range(0, numboxes, self.numrows):
            column_num = imin / self.numrows
            
            rows = []
            imax = imin + self.numrows
            if imax > numboxes: # last column, pad out with zeros
                diff = imax - numboxes
                rows = [np.zeros((self.box_size, self.box_size))] * diff
                imax = numboxes

            for pnt in points[imin:imax]:
                row_num = len(rows)
                detsec = self.detsec_finder.find_detector_section(pnt)
                box = Box(detsec,                # detector section to cull from
                          pnt,                   # point on the detector where the box should be (from MDF)
                          self.box_size,         # the size of the box to cull to
                          (column_num, row_num)) # the location in the final mosaic
            
                assert box.get_size() == self.box_size, "all boxes should be the exact same size"

                rows.append(box.get_data())
                self.boxes.append(box)
            
            debug("...row dimensions =", repr([r.shape for r in rows]))

            columns.append(np.vstack(rows))

        self.scidata = np.hstack(columns)

    def get_num_mos_boxes(self):
        return len(self.boxes)

    def get_boxes(self):
        return self.boxes
        
    def get_science_data(self):
        return self.scidata

    def get_box_borders(self):
        ydim, xdim = self.get_science_data().shape
        
        # the horizontal lines
        for rowid in range(1, self.numrows):
            ylevel = rowid * self.box_size
            yield (0, ylevel), (xdim, ylevel)
        
        # the vertical lines, note, we're assuming the mosaic is a square
        for colid in range(1, self.numrows):
            xlevel = colid * self.box_size
            yield (xlevel, 0), (xlevel, ydim)
        

class AcquisitionImage(object):
    def __init__(self, filename, mosmask=None, mdfdir=None):
        self.ad = AstroData(filename)
        self.mosmask = mosmask
        self.mdfdir = mdfdir

        # Determine extension
        nsci = len(self.ad)
        debug("...nsci = ", nsci)
            
        if nsci > 1:
            l_sci_ext = 1 
        else:
            l_sci_ext = 0

        debug("...using extension [" + str(l_sci_ext) + "]")

        overscan_dv = self.ad[l_sci_ext].overscan_section()

        if self.is_mos_mode():
            self.box_coords = parse_box_coords(self, self.get_mdf_filename())
            self.box_mosaic = BoxMosaic(self, self.box_coords)
            self.scidata = self.box_mosaic.get_science_data()
        elif self.is_new_gmosn_ccd():
            # tile the 2 center parts of the new GMOS image
            self.scidata = gmultiamp(self.ad)
        elif not overscan_dv.is_none():
            # remove the overscan so we don't have to take it into account when guessing the slit location
            self.scidata = subtract_overscan(self.ad[l_sci_ext])

            # it still affects the center of rotation however
            ox1, ox2, oy1, oy2 = overscan_dv.as_list()
            correction = np.array([ox2 - ox1, 0])
            center = self.get_binned_data_center() - correction
            self.fieldcenter = center * self.detector_y_bin()
        else:
            self.scidata = self.ad[l_sci_ext].data

    @cache
    def instrument(self):
        return str(self.ad.instrument())

    def is_new_gmosn_ccd(self):
        header = self.ad.phu.header
        if "DETECTOR" not in header:
            return False
        
        if header["DETECTOR"] == "GMOS + e2v DD CCD42-90":
            return True
        return False

    def get_science_data(self):
        assert self.scidata is not None
        return self.scidata

    @cache
    def unbinned_pixel_scale(self):
        return float(self.ad.pixel_scale()) / self.detector_y_bin()

    @cache
    def binned_pixel_scale(self):
        return float(self.ad.pixel_scale())

    def _check_binning(self):
        if int(self.ad.detector_x_bin()) != int(self.ad.detector_y_bin()):
            error("ERROR: incorrect binning!")
            error("Sorry about that, better luck next time.")
            sys.exit(1)

    @cache
    def detector_x_bin(self):
        self._check_binning()
        return int(self.ad.detector_x_bin())

    @cache
    def detector_y_bin(self):
        self._check_binning()
        return int(self.ad.detector_y_bin())

    @cache
    def program_id(self):
        return str(self.ad.program_id())

    @cache
    def observation_id(self):
        return str(self.ad.observation_id())

    @cache
    def saturation_level(self):
        dv = self.ad.saturation_level()
        return min(dv.as_list())

    @cache
    def focal_plane_mask(self):
        return str(self.ad.focal_plane_mask())

    @cache
    def grating(self):
        return str(self.ad.grating())

    def get_detector_size(self):
        # mos mode acquisitions don't necessarily have the entire
        # field of view in their data sections, so we have to rely on
        # other tricks to figure out the center of rotation.

        detsize = self.ad.phu_get_key_value("DETSIZE")
        xmin, xdim, ymin, ydim = extract_dimensions(detsize)
        
        # adjust for chip gaps
        nccds = int(self.ad.phu_get_key_value("NCCDS"))
        xdim += ((nccds - 1) * _obtain_unbinned_arraygap(self.ad))

        # adjust for un-illuminated pixels
        if self.is_gmos():
            ydim -= 36 # magic number that should be replaced with a lookup table later

        return xdim, ydim
        

    def get_field_center(self):
        """ The center of rotation in pixels. """
        if hasattr(self, "fieldcenter"):
            return self.fieldcenter

        if self.is_mos_mode():
            xdim, ydim = self.get_detector_size()
            return np.array([float(xdim) / 2.0, float(ydim) / 2.0])

        return self.get_data_center()

    def get_data_center(self):
        ydim, xdim = self.get_science_data().shape
        return np.array([float(xdim) / 2.0, float(ydim) / 2.0]) * self.detector_y_bin()

    def get_binned_data_center(self):
        return self.get_data_center() / self.detector_y_bin()

    def set_goal_center(self, center):
        self.goal_center = np.array(center)

    def get_goal_center(self):
        default = self.get_data_center()
        return getattr(self, "goal_center", default)

    def set_binned_custom_center(self, center):
        self.set_binned_goal_center(center)
        self.custom_center = True

    def has_custom_center(self):
        return getattr(self, "custom_center", False)

    def get_binned_goal_center(self):
        return self.get_goal_center() / self.detector_y_bin()

    def set_binned_goal_center(self, center):
        center = np.array(center) * self.detector_y_bin()
        self.set_goal_center(center)

    def get_mask_width(self):
        debug("...finding slit dimensions...")
        
        slitxbin = self.detector_x_bin()
        slitybin = self.detector_y_bin()
        debug("...slit image binning = ", slitxbin, " x ", slitybin)
        if slitxbin > 1 or slitybin > 1:
            warning("! WARNING: Slit image is binned " + slitxbin + " x " + slitybin)

        slitmask = self.focal_plane_mask()
        return float(slitmask.replace("arcsec", ""))

    def get_mask_width_in_pixels(self):
        return self.get_mask_width() / self.unbinned_pixel_scale()

    def get_slit_size_in_pixels(self):
        xsize = self.get_mask_width_in_pixels()
        ysize = self.get_science_data().shape[0]
        return xsize, ysize

    def get_expected_slit_tilt(self):
        if self.is_gmos():
            return 0.0

        error("Instrument is not supported, need to know an expected slit tilt")
        sys.exit(1)

    @property
    def phu(self):
        return self.ad.phu

    @property
    def filename(self):
        return self.ad.filename

    def get_program_id_parts(self):
        gemprgid = str(self.ad.program_id())
        parts = gemprgid.split("-")
        if len(parts) != 4:
            msg = "Cannot parse program id '%s'" % gemprgid
            error(msg)
            raise ValueError(msg)
        observatory, semester, prgtype, queuenum = parts
        return observatory, semester, prgtype, int(queuenum)

    def get_semester(self):
        """ Return something in the form of '2006B' """
        observatory, semester, prgtype, queuenum = self.get_program_id_parts()
        return semester

    def get_observatory_prefix(self):
        """ Return something in the form of 'GN' """
        observatory, semester, prgtype, queuenum = self.get_program_id_parts()
        return observatory

    def is_mos_mode(self):
        return self.mosmask is not None or self.has_mos_mask()

    @cache
    def has_mos_mask(self):
        if not self.is_gmos() and not self.is_f2():
            return False

        maskname = self.focal_plane_mask()
        if ("Q" in maskname or   # Queue program
            "C" in maskname or   # Classical program
            "D" in maskname or   # DD program
            "V" in maskname):    # SV program

            xbin = self.detector_x_bin()
            ybin = self.detector_y_bin()
            if xbin != 1 or ybin != 1:
                error ("MOS acquisition image binning must be 1x1, found %ix%i binning." % (xbin, ybin))
                clean()

            return True
        return False

    def has_mask_in_beam(self):
        maskname = self.focal_plane_mask().lower()
        if "imag" in maskname:
            return False

        slitmask = self.focal_plane_mask()    
        if self.is_gmos() and "arcsec" in slitmask:
            return True
    
        if self.is_gnirs() and "arcsec" in slitmask:
            acqmir = slitimage_ad.phu.header["ACQMIR"]
            debug("...acqmir = ", acqmir)
            if acqmir == "In":
                return True
    
        if self.is_niri() and "cam" not in slitmask:
            return True
    
        if self.is_f2() and "slit" in slitmask:
            return True

        if self.is_f2() and "mos" in slitmask:
            return True

        return self.has_mos_mask()

    def get_mdf_filename(self):
        if hasattr(self, "mdffile"):
            return self.mdffile

        self.mdffile = self._get_mdf_filename()
        return self.mdffile

    def _get_mdf_filename(self):
        if self.mosmask is not None:
            mdffile = self.mosmask

            # Expand the MOS mask number
            if is_number(self.mosmask):
                observatory, semester, prgtype, queuenum = self.get_program_id_parts()
                mdffile = "%s%s%s%03i-%02i" % (observatory, semester, prgtype, queuenum, int(self.mosmask))
                debug("...mosmask =", mdffile)
        else:
            mdffile = self.focal_plane_mask()

        mdffile = fits_filename(mdffile)

        #-----------------------------------------------------------------------
        # Start searching around willy nilly for the MDF file
        if os.path.exists(mdffile):
            return mdffile

        # note, the order in which directories are added to this list gives priority
        dirs = []
        if self.mdfdir is not None:
            dirs.append(self.mdfdir)

        dname = os.path.dirname(self.filename)
        if dname and dname != ".":
            dirs.append(dname)

        dirs.append(os.getcwd())

        # search through semester directories as well
        semester_dir = self.get_observatory_prefix() + self.get_semester()
        directories_to_search = []
        for dname in dirs:
            directories_to_search.append(dname)
            dname = os.path.join(dname, semester_dir)
            directories_to_search.append(dname)

        # now search through the directories
        for dname in directories_to_search:
            fname = os.path.join(dname, mdffile)
            debug("...trying", fname)
            if os.path.exists(fname):
                return fname

        raise ValueError("Unable to find MDF file named '%s'" % mdffile)

    def get_num_mos_boxes(self):
        return self.box_mosaic.get_num_mos_boxes()

    def get_mos_boxes(self):
        return self.box_mosaic.get_boxes()

    def get_mos_box_borders(self):
        for border in self.box_mosaic.get_box_borders():
            yield border

    @cache
    def get_min_slitsize(self):
        mdffile_ad = AstroData(self.get_mdf_filename())

        xsize = Ellipsis
        ysize = Ellipsis
        for row in mdffile_ad["MDF"].data:
            # select the alignment boxes, designated by priority 0
            if row["priority"] not in ["1", "2", "3"]:
                continue

            xsize = min(xsize, row["slitsize_x"])
            ysize = min(ysize, row["slitsize_y"])

        return xsize, ysize

    def get_extensions(self):
        for ext in self.ad:
            yield ext

    def _get_lazy_detector_section_finder(self):
        if not hasattr(self, "detsec_finder"):
            self.detsec_finder = DetectorSectionFinder(self)
        return self.detsec_finder

    def get_box_size(self):
        return self._get_lazy_detector_section_finder().get_box_size()

    def find_detector_section(self, point):
        return self._get_lazy_detector_section_finder().find_detector_section(point)

    def get_full_field_of_view(self):
        return self._get_lazy_detector_section_finder().get_full_field_of_view()

    def is_altair(self):
        aofold = self.ad.phu.header["AOFOLD"]
        if aofold == "IN":
            return True
        return False

    def is_south_port(self):
        inportnum = int(self.ad.phu.header["INPORT"])
        if inportnum == 1:
            return True
        return False

    def is_type(self, typename):
        return self.ad.is_type(typename)

    def is_gmos(self):
        return self.is_type("GMOS")
    
    def is_gmosn(self):
        return self.is_type("GMOS_N")

    def is_gmoss(self):
        return self.is_type("GMOS_S")

    def is_gnirs(self):
        return self.is_type("GNIRS")

    def is_f2(self):
        return self.is_type("F2")

    def is_nifs(self):
        return self.is_type("NIFS")

    def is_niri(self):
        return self.is_type("NIRI")

def main(argv=[__name__]):
    idx = 1
    if argv[1] == "-v":
        setup_logging(None, True)
        idx = 2
    
    ad = AcquisitionImage(argv[idx])

    ui.display(ad.get_science_data(), zscale=True)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
