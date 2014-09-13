import os
import math
import copy

if not os.environ.has_key ("NUMERIX"):
    print "Warning:  environment variable NUMERIX should be set to numpy"
elif os.environ["NUMERIX"] != "numpy":
    print "Warning:  NUMERIX is set to '%s', should be numpy" % \
          os.environ["NUMERIX"]

import numpy as N
from pyraf import iraf
try:
    import stsci.ndimage as ndimage
except ImportError:
    import ndimage
import pyfits

import glog
import gemutil
import irafutil

__version__ = "1.0 (2007 January 25)"

class Imcoadd:
    """Combine science images, includes cleaning for cosmic ray events

       The public methods are:
           find (threshold=None, fwhm=None)     [optional, called by init]
           wcs_map (box=20., rotate=False, scale=False,
                    order=3, sigfit=2.5, niter=5,
                    coolimit=0.3, fl_inter=False, fl_mark=False)
           user_map (box=20., rotate=False, scale=False,
                     order=3, sigfit=2.5, niter=5,
                     fl_refmark=False, dispmag=1.,
                     coolimit=0.3, fl_inter=False, fl_mark=False)
           twodx_map (box=20., rotate=False, scale=False,
                      order=3, sigfit=2.5, niter=5,
                      asection="default", xwindow=181,
                      coolimit=0.3, fl_inter=False, fl_mark=False)
           header_map (box=20., rotate=False, scale=False,
                       order=3, sigfit=2.5, niter=5,
                       key_xoff="default", key_yoff="default",
                       instscale=1., xsign="default", ysign="default",
                       key_pixscale="PIXSCALE", pixscale=1.,
                       coolimit=0.3, fl_inter=False, fl_mark=False)
           transform()                  [optional, called by median()]
           median (output=None, aperture=30.)
           add (limit=15., key_limit=True, lowsigma=7., lowlimit=500.,
                scalenoise=0., growthrad=1)
           average()
           close()

       Output files:
           <reference_image>_pos is output from daofind for reference image
               [written by find()]
           <rootname>_cen is an output from apphot.center
               [written by one of the *_map() functions]
           <rootname>_trn contains x0, y0, x, y;
               x0 & y0 are coordinates in the reference image, while
               x & y are coordinates in the current image (<rootname>.fits)
               [written by one of the *_map() functions]
           <rootname>_trn.fits is the geometrically transformed image
               [written by transform(), which is called by median()]
           <rootname>_trn_mag is output from apphot.phot
               [written by median()]
           <rootname>badpix.fits
               [written by median(), updated by add()]
           <reference_image>_med.fits
               [written by median()]
           <reference_image>_add.fits
               [written by add()]
           <reference_image>_avg.fits
               [written by average()]
           the database file (default name "imcoadd.dat") contains the
               coefficients of the mapping determined by geomap or wcsmap
               in one of the *_map() functions
           the log file (default name "imcoadd.log") contains a summary of
               the processing history

       The following IRAF tasks are used:
           mktemp
           imaccess
           proto.fixpix
           images.immatch.wcsmap
           images.immatch.geomap
           images.immatch.geotran
           images.immatch.geoxytran
           images.immatch.xregister
           images.immatch.imcombine
           images.imutil.imexpr
           images.imutil.imarith
           images.imutil.imstatistics
           images.imutil.imcopy
           images.imutil.imdelete
           images.tv.display
           images.tv.tvmark
           noao.digiphot.apphot.center
           noao.digiphot.apphot.daofind
           noao.digiphot.apphot.phot
           noao.digiphot.ptools.pdump
           noao.imred.ccdred.badpiximage
    """

    def __init__ (self, images, outimage=None,
                  sci_ext="SCI", var_ext="VAR", dq_ext="DQ",
                  immasks="DQ", database="imcoadd.dat",
                  threshold=20., fwhm=7.,               # used by find()
                  key_inst="INSTRUME", key_camera="CAMERA",
                  key_inport="INPORT",
                  geofitgeom="rscale",
                  geointer="linear",
                  geonxblock=2048, geonyblock=2048,
                  key_ron="RDNOISE", key_gain="GAIN", ron=1., gain=1.,
                  datamin=-1000., key_sat="SATURATI", datamax=50000.,
                  statsec="default",
                  badpixfile="default",
                  fl_fixpix=True, fl_scale=True,
                  fl_overwrite=False,
                  logfile="imcoadd.log", verbose=True):

        # Convert to host file name, and open the log file.
        logfile = irafutil.expandFileName (logfile)
        self._gl = glog.GLog (logfile, curtask="imcoadd", verbose=verbose)

        # Copy parameters to attributes.
        self._images = images
        self._n_files = None            # assigned later
        self._fl_mef = None             # assigned later
        self._outimage = outimage
        self._immasks = immasks
        self._tmpmli = None             # assigned later
        self._zero_ext = 0
        # these six may be overridden in _setup()
        self._sci_ext = sci_ext
        self._var_ext = var_ext
        self._dq_ext = dq_ext
        self._l_sci_ext = "[" + str (self._sci_ext) + "]"
        self._l_var_ext = "[" + str (self._var_ext) + "]"
        self._l_dq_ext = "[" + str (self._dq_ext) + "]"
        self._database = database
        self._threshold = threshold     # used only by find()
        self._fwhm = fwhm               # used only by find()
        self._key_inst = key_inst
        self._key_camera = key_camera
        self._key_inport = key_inport
        self._geofitgeom = geofitgeom
        self._geointer = geointer
        self._geonxblock = geonxblock
        self._geonyblock = geonyblock
        if key_ron.strip():
            self._key_ron = key_ron
        else:
            self._key_ron = "RDNOISE"
        if key_gain.strip():
            self._key_gain = key_gain
        else:
            self._key_gain = "GAIN"
        self._default_ron = ron
        self._default_gain = gain
        self._datamin = datamin
        self._key_sat = key_sat
        self._datamax = datamax
        self._statsec = statsec
        self._badpixfile = badpixfile
        self._fl_fixpix = fl_fixpix
        self._fl_scale = fl_scale
        self._fl_overwrite = fl_overwrite
        self._fl_obj = True             # the default is that there are objects
        self._fl_msk = True             # default
        self._tmpres = iraf.mktemp ("tmpres")
        self._verbose = verbose
        self._transform_done = False    # True if transform() has been run
        self._median_done = False       # True if median() has been run

        # These are initial values.
        self._ron = []
        self._gain = []
        self._roneff = 0.0
        self._gaineff = 0.0
        self._ron_ref = 0.0
        self._gain_ref = 0.0

        # assigned later in _setup():
        self._rootnames = []
        self._n_coords = []
        self._image_coords = []
        self._trn_images = []   # names of geometrically transformed images
        self._Xmax = None
        self._Ymax = None
        self._n_xsamp = None
        self._n_ysamp = None

        # average of MED_SKY in _trn headers (assigned by getMedianSky())
        self._median_sky = None

        # Expand comma-separated (or @file) string to a list of file names.
        self._checkImages()

        # Create tmpmli list from immasks or from DQ extensions; run fixpix.
        self._initMasks()

        # The reference image is the first image in the list.  Use the
        # rootname of that file to construct the default output name,
        # <reference>_add.fits
        if outimage is None:
            self._outimage = gemutil.appendSuffix (self._images[0], "_add")
        else:
            # convert to host file name
            self._outimage = irafutil.expandFileName (outimage)

        self._setup()

        # Get the instrument name, camera name and port number from the header.
        self._instrument = gemutil.getkey (self._key_inst,
                               self._images[0], self._zero_ext)
        self._camera = ""
        if self._instrument == "NIRI":
            self._camera = gemutil.getkey (self._key_camera,
                               self._images[0], self._zero_ext, default="f6")
        elif self._instrument not in ["NIRI", "Hokupaa+QUIRC",
                                      "GMOS-N", "GMOS-S"]:
            self._gl.glogprint (
            "Instrument definition in header not recognized", type="warning")
            self._instrument = ""
        self._inport = gemutil.getkey (self._key_inport,
                                       self._images[0], self._zero_ext)

        self._gl.glogprint ("Instrument + camera      : " +
                            self._instrument + "  " + self._camera)
        self._gl.glogprint ("Telescope ISS port number: " + str (self._inport))

        refim = self._images[0]

        # Get datamax if from header
        if self._key_sat:
            self._datamax = gemutil.getkey (self._key_sat, refim,
                                self._zero_ext, default=self._datamax)

        if badpixfile == "default" and self._instrument == "":
            badpixfile = ""

        if badpixfile:
            if badpixfile == "default":
                if self._instrument == "Hokupaa+QUIRC":
                    badpixfile = "quirc$data/quirc_bpm.pl"
                elif self._instrument == "NIRI":
                    badpixfile = "niri$data/niri_bpm.pl"
                elif self._instrument in ["GMOS", "GMOS-N", "GMOS-S"]:
                    ccdsum = gemutil.getkey ("CCDSUM", refim, self._sci_ext,
                                             default=1)
                    words = ccdsum.split()
                    # repeat the first value; see imcoadd.cl line 567 or 569
                    ccdsum2 = words[0] + words[0]
                    if self._instrument in ["GMOS", "GMOS-N"]:
                        badpixfile = "gmos$data/mgmosn_bpm" + ccdsum2 + ".pl"
                    else:
                        badpixfile = "gmos$data/mgmoss_bpm" + ccdsum2 + ".pl"
        else:
            badpixfile = "gemtools$badpix.none"

        self._badpixfile = irafutil.expandFileName (badpixfile)
        if not os.access (self._badpixfile, os.R_OK):
            self._gl.glogprint ("Bad pixel file " + self._badpixfile +
                                " not found", type="warning")
            self._badpixfile = irafutil.expandFileName ("gemtools$badpix.none")
        self._badim = (iraf.imaccess (self._badpixfile) == 1)
        self._gl.glogprint ("Bad pixel file " + self._badpixfile)

        self._gl.glogprint ("Read noise: %6.2f e-   Gain: %6.2f e-/ADU" %
                            (self._ron_ref, self._gain_ref))

        # Find objects in the reference image, unless already done.
        self.find (rerun=False)

    def _setup (self):

        # Cut off filename extensions of .fits or .imh (or .hhh).
        self._rootnames = gemutil.removeExtension (self._images)

        # output from daofind for reference image
        self._n_coords = self._rootnames[0] + "_pos"

        # This is a list of strings, each containing x & y coordinates
        # from the last image processed; this is used in the loop over
        # input images in the _map routines.  The initial value is
        # ref_coords, coordinates in the reference image.
        self._image_coords = []

        # names of geometrically transformed images
        self._trn_images = []
        for image in self._images:
            self._trn_images.append (gemutil.appendSuffix (image, "_trn"))

        # check some parameter values
        if self._geofitgeom not in ["shift", "xyscale", "rotate",
                                    "rscale", "rxyscale", "general"]:
            self._crash (ValueError, 'geofitgeom must be one of ' \
                         '"shift", "xyscale", "rotate", ' \
                         '"rscale", "rxyscale", "general"')

        if self._geointer not in ["nearest", "linear", "poly3",
                                  "poly5", "spline3"]:
            self._crash (ValueError, "geointer is not valid")

        self._gl.glogprint ("Images (and masks) in list")
        self._gl.glogprint (repr (self._images))
        if self._fl_msk:
            self._gl.glogprint (repr (self._tmpmli))

        # fl_mef = True means that the input images are multi-extension FITS
        # (MEF) files; False means simple FITS.
        if not self._fl_mef:
            self._sci_ext = 0
            self._l_sci_ext = "[0]"
            self._var_ext = None
            self._dq_ext = None
            self._l_var_ext = None
            self._l_dq_ext = None

        (self._Xmax, self._Ymax) = gemutil.getkeys (["NAXIS1", "NAXIS2"],
                     self._images[0], self._sci_ext)

        if self._statsec == "default":
            self._statsec = "[100:%d,100:%d]" % (self._Xmax-100, self._Ymax-100)
        self._gl.glogprint ("Statistics section: " + self._statsec)

        # Only some values work for sampling (used by transform and median).
        if self._Xmax > 1000:
            self._n_xsamp = 20
        elif self._Xmax > 500:
            self._n_xsamp = 10
        else:
            self._n_xsamp = 1
        if self._Ymax > 1000:
            self._n_ysamp = 20
        elif self._Ymax > 500:
            self._n_ysamp = 10
        else:
            self._n_ysamp = 1

        # Get and check readout noise (ron) and gain in images.
        for image in self._images:
            l_ron = gemutil.getkey (self._key_ron, image, 0,
                                    default=self._default_ron)
            self._ron.append (l_ron)
            self._roneff += l_ron**2
            l_gain = gemutil.getkey (self._key_gain, image, 0,
                                     default=self._default_gain)
            self._gain.append (l_gain)
            self._gaineff += l_gain

        self._ron_ref = self._ron[0]
        self._gain_ref = self._gain[0]

    def _checkImages (self):
        """Expand comma-separated names, and check for existence.

        This function expands the list of file names in self._images,
        checks that all the specified files exist, and checks whether
        they are all multi-extension FITS or simple FITS.
        The following attributes are updated or assigned:
            _images     list of input images
            _n_files    number of input images
            _fl_mef     True if multi-extension FITS
        """

        # Expand comma-separated (or @file) string to a list of file names.
        self._images = gemutil.expandlist (self._images)
        # Append ".fits" to each file name that lacks an extension.
        self._images = gemutil.appendFits (self._images)
        self._n_files = len (self._images)
        if self._n_files == 0:
            self._crash (ValueError, "No input image specified")
        if self._n_files == 1:
            self._crash (ValueError, "Only one input image in list")

        # Check that the specified images actually exist, and check
        # whether they are multi-extension FITS or something else.
        missing = []
        first = True
        for image in self._images:
            (type, has_dq) = gemutil.gimverify (image, self._sci_ext)
            if type == gemutil.gimverify_does_not_exist:
                missing.append (image)
            elif first:
                l_type = type
                l_has_dq = has_dq
                first = False
            elif type != l_type:
                self._crash (RuntimeError,
                    "Mix of MEF and simple FITS, or other image formats")
        if missing:
            self._gl.glogprint ("The following input images are missing:  " +
                                repr(missing), type="error")
            self._crash (ValueError, "Missing input image(s)")

        if l_type == gemutil.gimverify_simple_FITS:
            self._fl_mef = False
            self._has_dq = False
        elif l_type == gemutil.gimverify_MEF:
            self._fl_mef = True
            self._has_dq = l_has_dq
        else:
            self._crash (ValueError, "Image type is not supported.")

    def _initMasks (self):
        """Create and populate tmpmli, the list of mask names.

        If _immasks is equal to "DQ" and the input images are in
        multi-extension FITS files, the DQ extension for each image
        will be used to create a temporary pixel list file (OK if DQ
        is even and bad if DQ is odd), and the names of these files
        will be saved in the attribute _tmpmli (tmpmli is the name used
        in imcoadd.cl).  If fl_fixpix is True, fixpix will be run to
        fix the input images, using the pixel lists as masks.

        If immasks is a list of file names or an @file list, that
        list will be expanded and assigned to tmpmli.

        If immasks is "DQ" and the input images are simple FITS, or if
        immasks is "none" (a string, not None), fl_msk will be set to
        False and tmpmli will be set to None.

        The following attributes are updated or assigned:
            _fl_msk     True if there are input masks
            _tmpmli     list of input masks
        """

        if self._immasks == "none":
            self._fl_msk = False
            self._tmpmli = None
            return

        if self._fl_mef and not self._has_dq:
            # MEF format, but no DQ extension
            self._fl_msk = False
            self._tmpmli = None

        elif self._immasks == "DQ" and self._fl_mef:
            self._fl_msk = True

            # ---- MEF: If DQ available: fixpix, use DQ as immasks
            # Create a list of file names for pixel lists.
            self._tmpmli = []
            tmpdq = iraf.mktemp ("tmpdq")
            for n in range (self._n_files):
                maskname = tmpdq + str(n+1) + ".pl"
                self._tmpmli.append (maskname)

            # (imcoadd.cl, line 405, also included a test on n_fl_med)
            self._gl.glogprint ("Making individual masks from DQ")
            for n in range (self._n_files):
                image = self._images[n]
                iraf.imexpr ("(int(a/2.) * 2 != a) ? 1 : 0",
                             self._tmpmli[n],
                             a=image+self._l_dq_ext, dims="auto",
                             intype="auto", outtype="auto",
                             bwidth=0, btype="nearest", bpixval=0.,
                             rangecheck=True, verbose=False, exprdb="none")

            if self._fl_fixpix:
                self._gl.glogprint ("Fixpix input images using DQ")
                for n in range (self._n_files):
                    image = self._images[n]
                    iraf.proto.fixpix (image + self._l_sci_ext, self._tmpmli[n],
                                       linterp="INDEF", cinterp="INDEF",
                                       verbose=False, pixels=False)
                    gemutil.gemhedit (image, self._zero_ext, "IMCOFIX", "done",
                                      "Imcoadd: fixpix using DQ")

        elif self._immasks == "DQ" and not self._fl_mef:
            self._fl_msk = False
            self._tmpmli = None

        else:
            self._fl_msk = True

            self._tmpmli = gemutil.expandlist (self._immasks)
            self._tmpmli = gemutil.appendFits (self._tmpmli)
            self._fl_msk = (len (self._tmpmli) > 0)
            missing = []
            for maskname in self._tmpmli:
                if maskname != "none" and not os.access (maskname, os.R_OK):
                    missing.append (maskname)
            if missing:
                self._gl.glogprint ("The following mask images are missing:  " +
                                    repr(missing), type="error")
                self._crash (ValueError, "Missing input image mask(s)")
            if self._fl_msk and (len (self._tmpmli) != self._n_files):
                self._crash (ValueError,
                    "The number of images and image masks must be the same")

    def _checkMapRanges (self, box, order, sigfit, niter, coolimit):
        """Check that *_map parameters are within valid ranges."""

        if box < 5.:
            self._crash (ValueError, "box is less than 5.")

        if order < 2:
            self._crash (ValueError, "order is less than 2")

        if sigfit < 0.:
            self._crash (ValueError, "sigfit is less than 0.")

        if niter < 0:
            self._crash (ValueError, "niter is less than 0")

        if coolimit < 0.:
            self._crash (ValueError, "coolimit is less than 0.")

    def _getScaleInfo (self, key_xoff, key_yoff,
                      instscale, xsign, ysign, key_pixscale, pixscale):
        """Compute scale info if default was given.

        This is called only by header_map().
        """

        refim = self._images[0]

        # Assign default values that may be replaced below.
        if xsign == "default" or xsign == "positive":
            xsign = 1
        else:
            xsign = -1
        if ysign == "default" or ysign == "positive":
            ysign = 1
        else:
            ysign = -1

        if key_pixscale:
            pixscale = gemutil.getkey (key_pixscale,
                                refim, self._zero_ext, default=pixscale)

        self._gl.glogprint ("Pixel scale: " + repr(pixscale))

        if key_xoff == "default" or key_yoff == "default":
            key_xoff = "XOFFSET"        # will be replaced if NIRI
            key_yoff = "YOFFSET"
            if self._instrument == "GMOS" or self._instrument == "GMOS-N":
                instscale = 1.
                if self._inport == 1:
                    xsign = -1
                    ysign = 1
                else:
                    xsign = 1
                    ysign = 1
            elif self._instrument == "GMOS-S":
                instscale = 1.
                if self._inport == 1:
                    xsign = -1
                    ysign = 1
                else:
                    xsign = -1
                    ysign = -1
            elif self._instrument == "Hokupaa+QUIRC":
                instscale = 1.611444
                xsign = -1
                ysign = -1
            elif self._instrument == "NIRI":
                key_xoff = "YOFFSET"
                key_yoff = "XOFFSET"
                instscale = 1.
                if self._camera == "f6":
                    if self._inport == 1:
                        xsign = -1
                        ysign = 1
                    else:
                        xsign = 1
                        ysign = 1
                else:                   # camera != "f6"
                    if self._inport == 1:
                        xsign = 1
                        ysign = -1
                    else:
                        xsign = -1
                        ysign = -1
            else:
                key_xoff = "YOFFSET"
                key_yoff = "XOFFSET"

        return (key_xoff, key_yoff, instscale, xsign, ysign, pixscale)

    def find (self, threshold=None, fwhm=None, rerun=True):
        """Find objects using daofind, or use existing file

        Note that a user would not ordinarily run this with rerun=False.

        @param threshold: if not None, this overrides the original value
        @type threshold: float

        @param fwhm: if not None, this overrides the original value
        @type fwhm: float

        @param rerun: if True, delete output file and rerun;
            if False, only rerun if output doesn't already exist
        @type rerun: boolean
        """

        if rerun:
            # if output already exists, delete it so we can recreate it
            self._deleteOutputFiles (self._n_coords)
        elif os.access (self._n_coords, os.R_OK):
            # output already exists, so we don't need to do anything
            return

        self._gl.glogprint ("Running find()")

        if threshold is not None:
            self._threshold = threshold
        if fwhm is not None:
            self._fwhm = fwhm

        if self._fwhm < 1.:
            self._crash (ValueError, "fwhm is less than 1.")

        # Get sky in order to estimate sigma of the sky -
        #  accurate value not required
        ref_sci_ext = self._images[0] + self._l_sci_ext
        medsky = iraf.imstat (ref_sci_ext + self._statsec,
                              fields="midpt",
                              lower=self._datamin, upper=self._datamax,
                              nclip=0, lsigma=iraf.INDEF, usigma=iraf.INDEF,
                              binwidth=0.1, format=False, cache=False,
                              Stdout=1)
        medsky = float (medsky[0])

        self._gl.glogprint ("Finding objects in " + self._images[0])
        sigma = math.sqrt ((self._ron_ref / self._gain_ref)**2 +
                           medsky / self._gain_ref)
        iraf.daofind (ref_sci_ext, output=self._n_coords,
                      starmap="", skymap="",
                      boundary="nearest", constant=0.,
                      interactive=False, icommands="",
                      gcommands="", wcsout="logical",
                      cache=False, verify=False,
                      update=")_.update", verbose=")_.verbose",
                      graphics=")_.graphics", display=")_.display",
                      scale=1., fwhmpsf=self._fwhm, emission=True,
                      sigma=sigma,
                      datamin=self._datamin, datamax=self._datamax,
                      noise="poisson", ccdread="", gain="",
                      readnoise=0., epadu=1., exposure="",
                      airmass="", filter="", obstime="",
                      itime=1., xairmass=iraf.INDEF,
                      ifilter="iraf.INDEF", otime="INDEF",
                      threshold=self._threshold, nsigma=2.,
                      ratio=1., theta=0.,
                      sharplo=0.2, sharphi=1.,
                      roundlo=-1., roundhi=1.,
                      mkdetections=False)

    def wcs_map (self, box=20., rotate=False, scale=False,
                 order=3, sigfit=2.5, niter=5,
                 coolimit=0.3, fl_inter=False, fl_mark=False):
        """Derive the transformation using GEOMAP and WCS keywords."""

        self._gl.glogprint ("Running wcs_map()")

        self._checkMapRanges (box, order, sigfit, niter, coolimit)

        (rotate, sigfit, niter, ref_coords) = \
                self._map_setup ("wcs", box, rotate, sigfit, niter)

        (cd1_1, cd1_2) = gemutil.getkeys (["CD1_1", "CD1_2"], self._images[0],
                                 0, must_exist=True)
        pixscale = math.sqrt (cd1_1**2 + cd1_2**2) * 3600.
        self._gl.glogprint ("Pixel scale: %8.4f" % pixscale)

        # ---- start of wcs alignment
        # SCAN through the other images
        # This image is the reference for the first image, after which
        # last_image will be replaced with the last image processed.
        last_image = self._images[0]
        # replaced by map_geomap within the loop
        self._image_coords = copy.copy (ref_coords)
        for n in range (1, self._n_files):
            image = self._images[n]
            record_name = self._rootnames[n] + "_trn"   # database record name

            iraf.wcsmap (image+self._l_sci_ext, last_image+self._l_sci_ext,
                         self._tmpres,
                         transforms=record_name, results="",
                         xmin=1, xmax=self._Xmax, ymin=1, ymax=self._Ymax,
                         nx=10, ny=10, wcs="world", transpose=False,
                         xformat="%10.3f", yformat="%10.3f",
                         wxformat="", wyformat="",
                         fitgeometry="rscale", function="polynomial",
                         xxorder=2, xyorder=2, xxterms="half",
                         yxorder=2, yyorder=2, yxterms="half",
                         reject=iraf.INDEF, calctype="real",
                         verbose=self._verbose, interactive=False)

            if niter > 0:
                last_image = image      # this image is reference for the next

            self._map_geomap (n, ref_coords, box, order, sigfit,
                              niter, coolimit, fl_inter, fl_mark)
        # ----- end of wcs alignment

    def user_map (self, box=20., rotate=False, scale=False,
                  order=3, sigfit=2.5, niter=5,
                  fl_refmark=False, dispmag=1.,
                  coolimit=0.3, fl_inter=False, fl_mark=False):
        """Derive the transformation using GEOMAP, with user input."""

        self._gl.glogprint ("Running user_map()")

        self._checkMapRanges (box, order, sigfit, niter, coolimit)

        (rotate, sigfit, niter, ref_coords) = \
                self._map_setup ("user", box, rotate, sigfit, niter)

        # ---- start of user alignment for reference image
        print " ==> Reference image: " + self._images[0]
        ref_sci_ext = self._images[0] + self._l_sci_ext
        iraf.display (ref_sci_ext, 1, xmag=dispmag, ymag=dispmag)
        if fl_refmark:
            print "Marking found objects in reference image"
            # extract the x & y coordinates and the ID
            ref_coords_tvmark = gemutil.fieldsOfTable (ref_coords, "0,1,3")
            iraf.tvmark (1, "STDIN", logfile="", autolog=False, outimage="",
                         deletions="", commands="", mark="point",
                         color=202, pointsize=1, txsize=1,
                         nxoffset=5, nyoffset=5, interactive=False,
                         tolerance=1.5, label=True, number=False,
                         font="raster", Stdin=ref_coords_tvmark)
            del ref_coords_tvmark
            print "Objects found and centered are marked"
        if not rotate and not scale:
            print "Point to one common object in reference image"
            print "    strike any key"
            (x11, y11) = irafutil.imcurXY()
        else:
            print "Point to first common object in reference image"
            print "    strike any key"
            (x11, y11) = irafutil.imcurXY()
            print "Point to second common object in reference image"
            print "    strike any key"
            (x12, y12) = irafutil.imcurXY()
        # ---- end of user alignment for reference image

        l_geofitgeom = self._geofitgeom         # may be modified locally
        if not rotate and not scale:
            l_geofitgeom = "shift"
        else:
            if rotate and not scale:
                l_geofitgeom = "rotate"
            if scale:
                l_geofitgeom = "rscale"

        # ---- start of user alignment
        # SCAN through the other images

        # replaced by map_geomap() within the loop
        self._image_coords = copy.copy (ref_coords)

        for n in range (1, self._n_files):
            image = self._images[n]
            record_name = self._rootnames[n] + "_trn"   # database record name

            print " ==> Image to be transformed:", image
            iraf.display (image+self._l_sci_ext, 1,
                          xmag=dispmag, ymag=dispmag)

            if not rotate and not scale:
                print "Point to one common object in " \
                          "image to be transformed"
                print "    coordinates for last image: %.1f, %.1f" % \
                                  (x11, y11)
                print "    strike any key"
                (x21, y21) = irafutil.imcurXY()
                cursor_positions = ["%7.1f %7.1f %7.1f %7.1f" % \
                                       (x11, y11, x21, y21)]
                # shift the variables
                if niter > 0:
                    x11 = x21
                    y11 = y21
            else:
                print "Point to first common object in " \
                      "image to be transformed"
                print "    coordinates for last image: %.1f, %.1f" % \
                                  (x11, y11)
                print "    strike any key"
                (x21, y21) = irafutil.imcurXY()
                print "Point to second common object in " \
                      "image to be transformed"
                print "    coordinates for last image: %.1f, %.1f" % \
                                  (x12, y12)
                print "    strike any key"
                (x22, y22) = irafutil.imcurXY()
                cursor_positions = \
                        ["%7.1f %7.1f %7.1f %7.1f" % (x11, y11, x21, y21),
                         "%7.1f %7.1f %7.1f %7.1f" % (x12, y12, x22, y22)]
                if niter > 0:
                    x11 = x21
                    y11 = y21
                    x12 = x22
                    y12 = y22

            # transform current image's coordinates to next image's
            # coordinates; allow shift, rotation and/or scale change
            iraf.geomap ("STDIN", self._tmpres,
                         xmin=1, xmax=self._Xmax, ymin=1, ymax=self._Ymax,
                         transforms=record_name, results="",
                         fitgeometry=l_geofitgeom, function="polynomial",
                         xxorder=2, xyorder=2, xxterms="half",
                         yxorder=2, yyorder=2, yxterms="half",
                         maxiter=0, reject=3., calctype="real",
                         verbose=self._verbose, interactive=False,
                         Stdin=cursor_positions)

            self._map_geomap (n, ref_coords, box, order, sigfit,
                              niter, coolimit, fl_inter, fl_mark)
        # ----- end of user alignment

    def twodx_map (self, box=20., rotate=False, scale=False,
                   order=3, sigfit=2.5, niter=5,
                   asection="default", xwindow=181,
                   coolimit=0.3, fl_inter=False, fl_mark=False):
        """Derive the transformation using GEOMAP and cross correlation."""

        self._gl.glogprint ("Running twodx_map()")

        self._checkMapRanges (box, order, sigfit, niter, coolimit)
        if xwindow < 11:
            self._crash (ValueError, "xwindow is less than 11")
        self._checkAsection (asection)          # check the format of asection

        (rotate, sigfit, niter, ref_coords) = \
                self._map_setup ("twodx", box, rotate, sigfit, niter)

        # ---- start of twodx alignment for the reference image
        if asection == "default":
            # [x11:x12,y11:y12] [x21:x22,y21:y22]
            x11 = int (0.25  * self._Xmax) + 1
            x12 = int (0.375 * self._Xmax)
            y11 = int (0.25  * self._Ymax) + 1
            y12 = int (0.375 * self._Ymax)
            x21 = int (0.625 * self._Xmax) + 1
            x22 = int (0.75  * self._Xmax)
            y21 = int (0.625 * self._Ymax) + 1
            y22 = int (0.75  * self._Ymax)
            asection = "[%d:%d,%d:%d] [%d:%d,%d:%d]" % \
                       (x11, x12, y11, y12, x21, x22, y21, y22)
        else:
            (x11, x12, y11, y12, x21, x22, y21, y22) = \
                    irafutil.splitSection (asection)

        xc1 = (x11 + x12) / 2.
        yc1 = (y11 + y12) / 2.
        xc2 = (x21 + x22) / 2.
        yc2 = (y21 + y22) / 2.
        xsize1 = x12 - x11
        ysize1 = y12 - y11
        xsize2 = x22 - x21
        ysize2 = y22 - y21

        if xwindow > min (xsize1, ysize1, xsize2, ysize2):
            xwindow = min (xsize1, ysize1, xsize2, ysize2)
            self._gl.glogprint ("Sections for X-correlation too small",
                                type="warning")
            self._gl.glogprint ("                   Resetting xwindow = " +
                                repr (xwindow))
        # ----- end twodx alignment for the reference image

        # ---- start of twodx alignment for the other images

        l_geofitgeom = self._geofitgeom         # may be modified locally
        if not rotate and not scale:
            l_geofitgeom = "shift"
        if rotate and not scale:
            l_geofitgeom = "rotate"
        if scale:
            l_geofitgeom = "rscale"

        # SCAN through the other images
        # This image is the reference for the first image, after which
        # last_image will be replaced with the last image processed.
        last_image = self._images[0]
        # replaced by map_geomap within the loop
        self._image_coords = copy.copy (ref_coords)
        for n in range (1, self._n_files):

            image = self._images[n]
            # this is a database record name
            record_name = self._rootnames[n] + "_trn"   # database record name

            self._gl.glogprint ("xregister on the sections " + repr (asection))
            self._gl.glogprint ("Window size for search %d pixels" % xwindow)

            if xwindow > 21:
                correlation = "fourier"
            else:
                correlation = "discrete"

            shift_info = iraf.immatch.xregister (image+self._l_sci_ext,
                      last_image+self._l_sci_ext, asection, shifts="STDOUT",
                      output="", databasefmt=True, append=True, records="",
                      coords="", xlag=0, ylag=0, dxlag=0, dylag=0,
                      background="none", border=iraf.INDEF,
                      loreject=iraf.INDEF, hireject=iraf.INDEF, apodize=0.,
                      filter="none", correlation=correlation,
                      xwindow=xwindow, ywindow=xwindow,
                      function="centroid", xcbox=box, ycbox=box,
                      interp_type="linear", boundary_type="nearest",
                      constant=0., interactive=False, verbose=self._verbose,
                      Stdout=1)
            self._gl.glogprint ("Output from xregister:")
            self._gl.glogprint (shift_info, type="list")

            words = shift_info[6].split()
            vx1 = float (words[2]) + xc1
            vy1 = float (words[3]) + yc1
            words = shift_info[7].split()
            vx2 = float (words[2]) + xc2
            vy2 = float (words[3]) + yc2

            cursor_positions = \
                    ["%.5f %.5f %.5f %.5f" % (vx1, vy1, xc1, yc1),
                     "%.5f %.5f %.5f %.5f" % (vx2, vy2, xc2, yc2)]

            iraf.geomap ("STDIN", self._tmpres,
                         xmin=1, xmax=self._Xmax, ymin=1, ymax=self._Ymax,
                         transforms=record_name, results="",
                         fitgeometry=l_geofitgeom, function="polynomial",
                         xxorder=2, xyorder=2, xxterms="half",
                         yxorder=2, yyorder=2, yxterms="half",
                         maxiter=0, reject=3., calctype="real",
                         verbose=self._verbose, interactive=False,
                         Stdin=cursor_positions)

            if niter > 0:
                last_image = image      # this image is reference for the next

            self._map_geomap (n, ref_coords, box, order, sigfit,
                              niter, coolimit, fl_inter, fl_mark)

    def _checkAsection (self, asection):
        """Check the format of the image section for alignment.

        This is used only by twodx_map().
        """

        bad = False             # initial value
        if asection != "default":
            sections = asection.split()
            len_sections = len (sections)
            if len_sections < 2:
                bad = True
                sections.append ("garbage")
            if len_sections > 4:
                bad = True
            if len_sections == 3:
                if sections[1][-1] == "]":
                    sections[0] = sections[0] + sections[1]
                    sections[1] = sections[2]
                else:
                    sections[1] = sections[1] + sections[2]
            elif len_sections == 4:
                sections[0] = sections[0] + sections[1]
                sections[1] = sections[2] + sections[3]
            if sections[0][0] != "[" or sections[0][-1] != "]" or \
               sections[1][0] != "[" or sections[1][-1] != "]":
                bad = True
            else:
                words1 = sections[0].split(",")
                words2 = sections[1].split(",")
                if len (words1) != 2 or len (words2) != 2:
                    bad = True
                else:
                    if words1[0].find(":") < 0 or words1[1].find(":") < 0 or \
                       words2[0].find(":") < 0 or words2[1].find(":") < 0:
                        bad = True

        if bad:
            self._gl.glogprint ("asection must be given as ...")
            self._gl.glogprint (
                    "[x11:x12,y11:y12] [x21:x22,y21:y22] or set to default")
            self._crash (errmess="")
        # ---- end of twodx alignment

    def header_map (self, box=20., rotate=False, scale=False,
                    order=3, sigfit=2.5, niter=5,
                    key_xoff="default", key_yoff="default",
                    instscale=1., xsign="default", ysign="default",
                    key_pixscale="PIXSCALE", pixscale=1.,
                    coolimit=0.3, fl_inter=False, fl_mark=False):
        """Derive the transformation using GEOMAP and header keywords."""

        self._gl.glogprint ("Running header_map()")

        self._checkMapRanges (box, order, sigfit, niter, coolimit)
        if instscale < 0.:
            self._crash (ValueError, "instscale is less than 0.")
        if xsign not in ["default", "negative", "positive"]:
            self._crash (ValueError,
                'xsign must be one of "default", "negative", "positive"')
        if ysign not in ["default", "negative", "positive"]:
            self._crash (ValueError,
                'ysign must be one of "default", "negative", "positive"')
        if pixscale < 0.000001:
            self._crash (ValueError, "pixscale is less than 0.000001")

        # Replace some or all default values, depending on instrument
        # configuration.
        (key_xoff, key_yoff, instscale, xsign, ysign, pixscale) = \
                self._getScaleInfo (key_xoff, key_yoff,
                      instscale, xsign, ysign, key_pixscale, pixscale)

        (rotate, sigfit, niter, ref_coords) = \
                self._map_setup ("header", box, rotate, sigfit, niter)

        # ---- start of header alignment for the reference image
        (x11_ref, y11_ref) = gemutil.getkeys ([key_xoff, key_yoff],
                                     self._images[0], 0, must_exist=True)
        x11 = self._Xmax / 2
        y11 = self._Ymax / 2
        # ----- end of header alignment for the reference image

        # ---- start of header alignment
        # SCAN through the other images
        # This image is the reference for the first image, after which
        # last_image will be replaced with the last image processed.
        last_image = self._images[0]
        # replaced by map_geomap within the loop
        self._image_coords = copy.copy (ref_coords)
        for n in range (1, self._n_files):
            image = self._images[n]
            record_name = self._rootnames[n] + "_trn"   # database record name

            (xoffset, yoffset) = gemutil.getkeys ([key_xoff, key_yoff],
                    image, 0, must_exist=True)
            x21 = xsign * (xoffset - x11_ref) * instscale / pixscale \
                      + self._Xmax / 2.
            y21 = ysign * (yoffset - y11_ref) * instscale / pixscale \
                      + self._Ymax / 2.
            cursor_positions = ["%7.1f %7.1f %7.1f %7.1f" % \
                                   (x11, y11, x21, y21)]

            iraf.geomap ("STDIN", self._tmpres,
                         xmin=1, xmax=self._Xmax, ymin=1, ymax=self._Ymax,
                         transforms=record_name, results="",
                         fitgeometry="shift", function="polynomial",
                         xxorder=2, xyorder=2, xxterms="half",
                         yxorder=2, yyorder=2, yxterms="half",
                         maxiter=0, reject=3., calctype="real",
                         verbose=self._verbose, interactive=False,
                         Stdin=cursor_positions)
            # shift the variables
            if niter > 0:
                x11 = x21
                y11 = y21

            if niter > 0:
                last_image = image      # this image is reference for the next

            self._map_geomap (n, ref_coords, box, order, sigfit,
                              niter, coolimit, fl_inter, fl_mark)
        # ----- end of header alignment

    def _map_setup (self, alignmethod, box, rotate, sigfit, niter):
        """Run apphot.center on the reference image.

        alignmethod is a string indicating which of the four _wcs
        methods called map_setup.  The arguments rotate, sigfit, niter
        can be modified by this method.  The function value is
        (rotate, sigfit, niter, ref_coords), where ref_coords is a
        list of strings, each one containing the X and Y coordinates
        in the reference image, "NoError", and the star ID.
        """

        # (imcoadd.cl, line 831)
        self._gl.glogprint ("Centering objects in " + self._images[0])

        ref_sci_ext = self._images[0] + self._l_sci_ext
        refim_cen = self._rootnames[0] + "_cen"
        self._deleteOutputFiles (refim_cen)
        iraf.apphot.center (ref_sci_ext, coords=self._n_coords,
                output=refim_cen,
                plotfile="", interactive=False, radplots=False,
                icommands="", gcommands="",
                wcsin="logical", wcsout="logical",
                cache=False, verify=False,
                update=")_.update", verbose=")_.verbose",
                graphics=")_.graphics", display=")_.display",
                scale=1., fwhmpsf=2.5,
                emission=True, sigma=iraf.INDEF,
                datamin=self._datamin, datamax=self._datamax,
                noise="poisson", ccdread="",
                gain="", readnoise=0., epadu=1.,
                exposure="", airmass="",
                filter="", obstime="", itime=1.,
                xairmass=iraf.INDEF, ifilter="INDEF",
                otime="INDEF", calgorithm="centroid",
                cbox=box, cthreshold=0.,
                minsnratio=1., cmaxiter=10,
                maxshift=box, clean=False,
                rclean=1., rclip=2.,
                kclean=3., mkcenter=False)

        if not os.access (refim_cen, os.R_OK):
            n_test = 0
        else:
            # changing settings to catch few objects
            ref_coords = iraf.pdump (refim_cen, "xcenter,ycenter,cerror,id",
                                "cerror=='NoError'",
                                headers=False, parameters=False, Stdout=1)
            n_test = len (ref_coords)

        # If no objects, the only working setting is
        # niter=0, fl_scale=no, alignmethod="header"
        # (imcoadd.cl, line 865)
        if n_test == 0:
            self._gl.glogprint ("No objects successfully centered",
                                type="warning")
            if alignmethod == "header" or alignmethod == "wcs":
                self._gl.glogprint ("Setting niter=0 and fl_scale=no",
                                    type="warning")
                niter = 0
                self._fl_scale = False
                self._fl_obj = False
            else:
                self._gl.glogprint ("Only alignmethod=header or wcs may work",
                                    type="error")
                if self._immasks == "DQ" and self._fl_msk:      # clean up
                    for maskname in self._tmpmli:
                        iraf.imdelete (maskname, verify=False)
                raise RuntimeError

        if n_test < 3 and niter > 1:
            self._gl.glogprint ("Too few objects.  Setting niter=1, sigfit=100",
                                type="warning")
            niter = 1
            sigfit = 100

        if n_test < 6 and niter > 2:
            self._gl.glogprint ("Too few objects.  Setting niter=2",
                                type="warning")
            niter = 2

        # Check the fit geometry depending on the number of objects; be safe.
        if n_test == 1:
            if self._geofitgeom != "shift":
                self._gl.glogprint (
                "Too few objects.  Setting geofitgeom=shift", type="warning")
                self._geofitgeom = "shift"
            if rotate:
                self._gl.glogprint ("Cannot determine rotation; " +
                                    "setting rotate=False", type="warning")
                rotate = False

        if n_test == 2 and self._geofitgeom not in ["shift", "rotate"]:
            self._gl.glogprint ("Too few objects.  Setting geofitgeom=rotate",
                                type="warning")
            self._geofitgeom = "rotate"

        if n_test < 6 and \
           self._geofitgeom not in ["shift", "rotate", "xyscale"]:
            if self._geofitgeom != "rscale":
                self._gl.glogprint (
                        "Too few objects.  Setting geofitgeom=rscale",
                        type="warning")
                self._geofitgeom = "rscale"

        self._gl.glogprint ("Number of objects successfully centered " +
                            "in reference image " + repr(n_test))

        if niter == 0:
            self._tmpres = self._database

        return (rotate, sigfit, niter, ref_coords)

    def _map_geomap (self, n, ref_coords, box, order, sigfit,
                     niter, coolimit, fl_inter, fl_mark):

        # this is a text file
        image_trn_file = self._rootnames[n] + "_trn"
        # this is a database record name (but same name as image_trn_file)
        record_name = self._rootnames[n] + "_trn"
        # this is a text file (called tmpcen in imcoadd.cl)
        image_cen_file = self._rootnames[n] + "_cen"
        self._deleteOutputFiles ([image_trn_file, image_cen_file])

        # (imcoadd.cl, line 1196)
        if self._fl_obj:
            image = self._images[n]

            new_image_coords = iraf.geoxytran ("STDIN", "STDOUT", self._tmpres,
                    record_name,
                    geometry="geometric", direction="forward",
                    xref=iraf.INDEF, yref=iraf.INDEF,
                    xmag=iraf.INDEF, ymag=iraf.INDEF,
                    xrotation=iraf.INDEF, yrotation=iraf.INDEF,
                    xout=iraf.INDEF, yout=iraf.INDEF,
                    xshift=iraf.INDEF, yshift=iraf.INDEF,
                    xcolumn=1, ycolumn=2, calctype="real",
                    xformat="", yformat="", min_sigdigit=7,
                    Stdin=self._image_coords, Stdout=1)
            # tmpres was set to database if niter is 0
            if self._tmpres != self._database:
                gemutil.deleteFile (self._tmpres)
            # Use this list of coordinates as the next reference.
            self._image_coords = copy.copy (new_image_coords)

            self._gl.glogprint ("Centering objects in " + image)
            self._gl.glogprint ("", "visual", type="visual", vistype="empty")

            # junk is a list of prompt messages (because coords="STDIN")
            junk = iraf.apphot.center (image+self._l_sci_ext, coords="STDIN",
                    output=image_cen_file, plotfile="",
                    interactive=False, radplots=False,
                    icommands="", gcommands="", wcsin="logical",
                    wcsout="logical", cache=False, verify=False,
                    update=")_.update", verbose=")_.verbose",
                    graphics=")_.graphics", display=")_.display",
                    scale=1., fwhmpsf=2.5, emission=True, sigma=iraf.INDEF,
                    datamin=self._datamin, datamax=self._datamax,
                    noise="poisson", ccdread="", gain="",
                    readnoise=0., epadu=1.,
                    exposure="", airmass="", filter="", obstime="", itime=1.,
                    xairmass=iraf.INDEF, ifilter="INDEF", otime="INDEF",
                    calgorithm="centroid", cbox=2.*box,
                    cthreshold=0., minsnratio=1.,
                    cmaxiter=10, maxshift=box,
                    clean=False, rclean=1., rclip=2., kclean=3.,
                    mkcenter=False, Stdin=self._image_coords, Stdout=1)

            coords = iraf.pdump (image_cen_file, "xcenter,ycenter,id", "yes",
                                 headers=False, parameters=False, Stdout=1)
            gemutil.deleteFile (image_cen_file)

            junk = iraf.apphot.center (image+self._l_sci_ext, coords="STDIN",
                    output=image_cen_file, plotfile="",
                    interactive=False, radplots=False,
                    icommands="", gcommands="", wcsin="logical",
                    wcsout="logical", cache=False, verify=False,
                    update=")_.update", verbose=")_.verbose",
                    graphics=")_.graphics", display=")_.display",
                    scale=1., fwhmpsf=2.5, emission=True, sigma=iraf.INDEF,
                    datamin=self._datamin, datamax=self._datamax,
                    noise="poisson", ccdread="", gain="",
                    readnoise=0., epadu=1.,
                    exposure="", airmass="", filter="", obstime="", itime=1.,
                    xairmass=iraf.INDEF, ifilter="INDEF", otime="INDEF",
                    calgorithm="centroid", cbox=box, cthreshold=0.,
                    minsnratio=1., cmaxiter=10, maxshift=box/1.5,
                    clean=False, rclean=1., rclip=2., kclean=3.,
                    mkcenter=False, Stdin=coords, Stdout=1)
            del junk

            new_coords = iraf.pdump (image_cen_file, "xcenter,ycenter,cerror",
                                     "yes", headers=False, parameters=False,
                                     Stdout=1)

            # ref_plus_new contains x0, y0, cerror, ID, x, y, cerror
            ref_plus_new = gemutil.joinlists (ref_coords, new_coords)
            good_coords = [x for x in ref_plus_new if x.count ("NoError") == 2]
            # coords_trn is a list of strings, corresponding to the text
            # file n_inim // "_trn" in imcoadd.cl
            coords_trn = gemutil.fieldsOfTable (good_coords, fields="0,1,4,5")

            # Save coords_trn to a text file, referred to as
            # n_inim // "_trn" in imcoadd.cl.
            # Note:  the cl script called 'unique' here.  This is handled
            # in _join_trn_files(), but the trn files themselves may still
            # contain duplicate stars.
            fd = open (image_trn_file, "w")
            last_line = ""
            for line in coords_trn:
                if line != last_line:   # skip adjacent, exact duplicates
                    fd.write ("%s\n" % line)
                last_line = copy.copy (line)
            fd.close()

            # Have to check that centering went ok
            # (imcoadd.cl, line 1260)
            if len (coords_trn) == 0 and niter > 0:
                self._gl.glogprint ("Failed to center any objects",
                                    type="error")
                self._gl.glogprint ("Adjust datamax and box as needed. " +
                                    "May have to use niter=0")
                self._gl.glogprint ("Check %s" % image_cen_file)
                self._crash (errmess="")

            if fl_mark:
                print "Marking objects used for transformation"
                coords_tvmark = \
                        gemutil.fieldsOfTable (ref_plus_new, fields="4,5,3")
                iraf.tvmark (1, "STDIN", logfile="", autolog=False,
                        outimage="", deletions="", commands="",
                        mark="point", color=202, pointsize=1, txsize=1,
                        nxoffset=5, nyoffset=5, interactive=False,
                        tolerance=1.5, label=True, number=False,
                        font="raster", Stdin=coords_tvmark)
                del coords_tvmark

        self._gl.glogprint ("Xmax = %5d   Ymax = %d" % (self._Xmax, self._Ymax))

        if niter > 0:
            self._gl.glogprint ("Entering geomap")
            self._gl.glogprint ("Fitting geometry " + self._geofitgeom)
            if self._geofitgeom == "general":
                self._gl.glogprint ("Fitting order %d" % order)
            self._gl.glogprint ("Iterating a maximum of %d times" % niter)
        else:
            self._gl.glogprint ("No iterations done with geomap")

        # Iterate for geomap - geomap cannot do this on its own...
        interactive = False
        break_loop = False
        # tmpdatabase is a scratch file for the database to be written by
        # geomap, and tmpresults is for the "results" file.  For the last
        # iteration, we'll set tmpdatabase to the name of the actual
        # database file and keep it.
        if niter > 0:
            tmpdatabase = iraf.mktemp ("tmpdatabase")
            tmpresults = iraf.mktemp ("tmpresults")
            if order > 2 and niter >= 3:
                ordoff = order - 2  # first 2 iterations with 2nd order fits
            else:
                ordoff = 0

        for count in range (1, niter+1):        # count is one indexed
            if count >= 3 or break_loop:
                ordoff = 0      # return to correct order
            if count == niter or break_loop:
                # optionally set to interactive for last iteration
                interactive = fl_inter
                tmpdatabase = self._database
                break_loop = True
                if fl_inter:
                    print ""
                    print "Last interation, number %d" % count
            if (break_loop or count == niter) and interactive:
                print "---------------------------------------------------"
                print " x   graph x residuals"
                print " y   graph y residuals"
                print " d   delete nearest point"
                print " u   undelete nearest point"
                print " f   make new fit"
                print " g   map data points"
                print " q   exit curve fitting"
                print "---------------------------------------------------"

            iraf.geomap (image_trn_file, tmpdatabase,
                         xmin=1, xmax=self._Xmax, ymin=1, ymax=self._Ymax,
                         transforms=record_name, results=tmpresults,
                         fitgeometry=self._geofitgeom, function="legendre",
                         xxorder=(order-ordoff), xyorder=(order-ordoff),
                         xxterms="half",
                         yxorder=(order-ordoff), yyorder=(order-ordoff),
                         yxterms="half",
                         maxiter=1, reject=sigfit, calctype="real",
                         verbose=self._verbose, interactive=interactive,
                         graphics="stdgraph", cursor="")

            if break_loop:
                gemutil.deleteFile (tmpresults)
                break
            else:
                # scratch file, until the last iteration
                gemutil.deleteFile (tmpdatabase)
                fd = open (tmpresults)
                lines = fd.readlines()
                fd.close()
                # format:  "#     Xin and Yin fit rms: <sigx>  <sigy>"
                words = lines[5].split()
                sigx = float (words[-2])
                sigy = float (words[-1])
                # stop when required accuracy is reached
                if sigx <= coolimit and sigy <= coolimit:
                    break_loop = True
                # check columns 6 and 7 from the file tmpresults, and
                # extract lines with low enough sigmas
                sigmas = gemutil.fieldsOfTable (tmpresults, "0,1,2,3,4,5,6,7")
                n_lines = 0
                low_sigmas = []
                for line in sigmas:
                    words = line.split()
                    if words[6] == "INDEF" or words[7] == "INDEF":
                        continue
                    v0 = float (words[6])
                    v1 = float (words[7])
                    if abs (v0) <= sigx*sigfit and abs (v1) <= sigy*sigfit:
                        low_sigmas.append (line)
                        n_lines += 1
                # (imcoadd.cl, line 1358)
                if n_lines >= 25:
                    final_coords = gemutil.fieldsOfTable (low_sigmas,
                                       "0,1,2,3")
                    gemutil.deleteFile (image_trn_file)
                    fd = open (image_trn_file, "w")
                    for line in final_coords:
                        fd.write ("%s\n" % line)
                    fd.close()
                    del final_coords
                else:
                    break_loop = True
                del sigmas, low_sigmas
                gemutil.deleteFile (tmpresults)

    def transform (self, rerun=True):
        """Transform the images using GEOTRAN.

        This function applies the geometric transformation to each input
        image (except the first, the reference), writing output images
        <rootname>_trn.fits that are accurately aligned with the reference
        image.  These _trn.fits files contain the input images that will be
        combined by median(), add() or average().  Either those files must
        already exist or transform() must be run before the images can be
        combined.  transform() may be called multiple times with rerun set
        to False (as it will be, by median(), add() and average()), and the
        bulk of the code in this function will only be run if the output
        files don't already exist.

        Note that a user would not ordinarily run this with rerun=False.

        @param rerun: if True, execute this function even if it's already
            been done; if False, only rerun if output doesn't already exist
        @type rerun: boolean
        """

        if self._transform_done and not rerun:
            return

        if not rerun:
            # All the output images must already exist; if not, we have to
            # rerun transform().
            for n in range (1, self._n_files):
                if not os.access (self._trn_images[n], os.R_OK):
                    rerun = True
                    break
        if rerun:
            # If some output images already exist, delete them.
            for n in range (1, self._n_files):
                gemutil.deleteFile (self._trn_images[n])
        else:
            # We don't have to rerun transform() because all the output
            # images already exist.
            self._transform_done = True
            return

        self._gl.glogprint ("Running transform()")
        self._gl.glogprint ("Sampling for geotran  : %d  %d" %
                            (self._n_xsamp, self._n_ysamp))
        self._gl.glogprint ("Block size for geotran: %d  %d" %
                            (self._geonxblock, self._geonyblock))

        # Reference image, not actually geometrically transformed,
        # just sky subtracted.
        ref_trn = self._trn_images[0]

        if not os.access (self._database, os.R_OK):
            self._crash (errmess="Database file missing: " + self._database)

        # Get sky for the reference image.
        ref_sci_ext = self._images[0] + self._l_sci_ext
        medsky = iraf.imstat (ref_sci_ext + self._statsec,
                              fields="midpt",
                              lower=self._datamin, upper=self._datamax,
                              nclip=0, lsigma=iraf.INDEF, usigma=iraf.INDEF,
                              binwidth=0.1, format=False, cache=False,
                              Stdout=1)
        medsky = float (medsky[0])
        # changed to get better sky estimate  27.11.95
        sigma = math.sqrt ((self._ron_ref / self._gain_ref)**2 +
                           medsky / self._gain_ref)
        lower = max (self._datamin, medsky - 5. * sigma)
        upper = max (50., medsky + 5. * sigma)
        medsky = iraf.imstat (ref_sci_ext + self._statsec,
                              fields="midpt",
                              lower=lower, upper=upper,
                              nclip=0, lsigma=iraf.INDEF, usigma=iraf.INDEF,
                              binwidth=0.1, format=False, cache=False,
                              Stdout=1)
        medsky = float (medsky[0])
        self._gl.glogprint ("")
        self._gl.glogprint ("Median sky level for %s: %8.1f" %
                            (self._images[0], medsky))

        # Subtract sky level from reference image, putting the result in
        # <ref>_trn.fits
        iraf.imarith (ref_sci_ext, "-", str (medsky), ref_trn,
                      title="", divzero=0., hparams="", pixtype="real",
                      calctype="real", verbose=False, noact=False)

        datetime = gemutil.gemdate()
        gemutil.gemhedit (ref_trn, 0, "GEM-TLM", datetime,
                          "UT Last modification with GEMINI")
        gemutil.gemhedit (ref_trn, 0, "IMCOADD", datetime,
                          "UT Time stamp for imcoadd")
        gemutil.gemhedit (ref_trn, 0, "MED_SKY", medsky, "Median sky level")

        # SCAN through the other images  - transform and subtract off sky
        for n in range (1, self._n_files):
            image = self._images[n]
            image_sci_ext = image + self._l_sci_ext
            image_trn = self._trn_images[n]
            # this is a database record name
            record_name = self._rootnames[n] + "_trn"

            self._gl.glogprint ("Transforming %s to %s" % (image, image_trn))
            # Check if input image is type short !!
            bitpix = gemutil.getkey ("bitpix", image, self._sci_ext)
            if bitpix == 16:
                tmpim = iraf.mktemp ("tmpim") + ".fits"
                (data, phdr, hdr) = gemutil.getdata (image)
                gemutil.putdata (tmpim, data.astype (N.float32), mef=False,
                                 phdr=phdr, clobber=True)
                realimage = tmpim + "[0]"
            else:
                tmpim = None
                realimage = image_sci_ext
            iraf.geotran (realimage, image_trn, self._database, record_name,
                          geometry="geometric",
                          xin=iraf.INDEF, yin=iraf.INDEF,
                          xshift=iraf.INDEF, yshift=iraf.INDEF,
                          xout=iraf.INDEF, yout=iraf.INDEF,
                          xmag=iraf.INDEF, ymag=iraf.INDEF,
                          xrotation=iraf.INDEF, yrotation=iraf.INDEF,
                          xmin=1, xmax=self._Xmax, ymin=1, ymax=self._Ymax,
                          xscale=1., yscale=1.,
                          ncols=iraf.INDEF, nlines=iraf.INDEF,
                          xsample=self._n_xsamp, ysample=self._n_ysamp,
                          interpolant=self._geointer, boundary="nearest",
                          constant=0., fluxconserve=True,
                          nxblock=self._geonxblock, nyblock=self._geonyblock,
                          verbose=self._verbose )
            if tmpim is not None:
                os.remove (tmpim)

            # Get sky for the current image
            medsky = iraf.imstat (image_sci_ext + self._statsec, fields="midpt",
                                  lower=self._datamin, upper=self._datamax,
                                  nclip=0, lsigma=iraf.INDEF, usigma=iraf.INDEF,
                                  binwidth=0.1, format=False, cache=False,
                                  Stdout=1)
            medsky = float (medsky[0])
            # changed to get better sky estimate  27.11.95
            sigma = math.sqrt ((self._ron[n] / self._gain[n])**2 +
                               medsky / self._gain[n])
            lower = max (self._datamin, medsky - 5. * sigma)
            upper = max (50., medsky + 5. * sigma)
            medsky = iraf.imstat (image_sci_ext + self._statsec, fields="midpt",
                                  lower=lower, upper=upper,
                                  nclip=0, lsigma=iraf.INDEF, usigma=iraf.INDEF,
                                  binwidth=0.1, format=False, cache=False,
                                  Stdout=1)
            medsky = float (medsky[0])
            self._gl.glogprint ("Median sky level for %s: %8.1f" %
                                (image, medsky))

            iraf.imarith (image_trn, "-", str (medsky), image_trn,
                          title="", divzero=0., hparams="", pixtype="real",
                          calctype="real", verbose=False, noact=False)
            datetime = gemutil.gemdate()
            gemutil.gemhedit (image_trn, 0, "GEM-TLM", datetime,
                              "UT Last modification with GEMINI")
            gemutil.gemhedit (image_trn, 0, "IMCOADD", datetime,
                              "UT Time stamp for imcoadd")
            gemutil.gemhedit (image_trn, 0, "MED_SKY", medsky,
                              "Median sky level")

        self._transform_done = True

        # end of SCAN through the images
        # end of transformation

    def median (self, output=None, aperture=30., rerun=True):
        """Use imcombine to get the median image.

        Note that a user would not ordinarily run this with rerun=False.

        @param output: if not None, this is the output file name
        @type output: string

        @param aperture: for aperture photometry
        @type aperture: float

        @param rerun: if True, execute this function even if it's already
            been done; if False, only rerun if output doesn't already exist
        @type rerun: boolean
        """

        if self._median_done and not rerun:
            return

        if not self._transform_done:
            self.transform (rerun=False)

        if output is None:
            # default output name, root of reference image + _med.fits
            ref_med = gemutil.appendSuffix (self._images[0], "_med")
        else:
            ref_med = output

        self._gl.glogprint ("Running median()")
        self._gl.glogprint ("Sampling for geotran  : %d  %d" %
                            (self._n_xsamp, self._n_ysamp))
        self._gl.glogprint ("Block size for geotran: %d  %d" %
                            (self._geonxblock, self._geonyblock))

        # Reference image, not actually geometrically transformed,
        # just sky subtracted.
        ref_trn = self._trn_images[0]

        if not os.access (self._database, os.R_OK):
            self._crash (errmess="Database file missing: " + self._database)

        # first subtract backgrounds as midpt for each image
        # the mean of the backgrounds is added to the median filtered image

        ref_badpix = self._rootnames[0] + "badpix.fits"
        ref_sci_ext = self._images[0] + self._l_sci_ext
        if self._badim:
            iraf.imcopy (self._badpixfile, ref_badpix, verbose=False)
        else:
            trimmed = getTrimmed (self._badpixfile) # 'trimmed' or 'untrimmed'
            self._gl.glogprint ("Bad pixel file: %s is for %s images" %
                                (self._badpixfile, trimmed))
            if trimmed == "untrimmed":
                tmpbad = iraf.mktemp ("tmpbad")         # a text file
                # Copy badpixfile to tmpbad, subtracting TRIMSEC offset.
                subtractTrimsec (self._badpixfile, tmpbad,
                                 self._images[0], self._sci_ext)
                self._badpixfile = tmpbad
            else:
                tmpbad = None
            iraf.ccdred.badpiximage (self._badpixfile, ref_sci_ext,
                      ref_badpix, goodvalue=0, badvalue=1)
            if tmpbad is not None:
                os.remove (tmpbad)

        # Get sky for first (reference) image from the header.  We will
        # add all the sky levels together using this variable.
        # (imcoadd.cl, line 1493)
        sum_medsky = gemutil.getkey ("MED_SKY", ref_trn)
        self._gl.glogprint ("Median sky level for %s: %8.1f" %
                            (self._images[0], sum_medsky))

        # put badpixel mask name in first image
        gemutil.gemhedit (ref_trn, 0, "BPM",
                          ref_badpix, "Bad pixel mask")

        # change pixel type to be able to transform this image
        (data, phdr, hdr) = gemutil.getdata (ref_badpix)
        gemutil.putdata (ref_badpix, data.astype (N.float32),
                         phdr=phdr, hdr=hdr, clobber=True)
        del (data, phdr, hdr)

        # input to imarith, if _fl_msk
        refimage_mask = gemutil.appendSuffix (self._images[0], "badpix")
        # sky level for other images, transform their badpixel mask
        # (imcoadd.cl, line 1506)
        for n in range (1, self._n_files):       # skip reference image

            image = self._images[n]
            image_trn = self._trn_images[n]
            image_mask = gemutil.appendSuffix (image, "badpix")
            record_name = self._rootnames[n] + "_trn"    # database record name

            medsky = gemutil.getkey ("MED_SKY", image_trn)
            self._gl.glogprint ("Median sky level for %s: %8.1f" %
                                (image, medsky))
            sum_medsky += medsky

            # make individual mask if required
            if self._fl_msk and self._tmpmli[n] != "none":
                mask = self._tmpmli[n]
                n_mask = iraf.mktemp ("tmpmsk") + ".fits"
                if iraf.imaccess (mask):
                    iraf.imarith (mask, "+", refimage_mask,
                                  n_mask,
                                  title="", divzero=0., hparams="",
                                  pixtype="real", calctype="real",
                                  verbose=False, noact=False)
                elif os.access (mask, os.R_OK):
                    tmpbad = iraf.mktemp ("tmpbad")
                    self._gl.glogprint ("Making individual mask for " + image)
                    trimmed = getTrimmed (mask)
                    self._gl.glogprint ("Bad pixel file: %s is for %s images" %
                                        (mask, trimmed))
                    if trimmed == "untrimmed":
                        # Copy mask to tmpbad, subtracting TRIMSEC offset
                        subtractTrimsec (mask, tmpbad,
                                         self._images[0], self._sci_ext)
                    else:
                        gemutil.copyFile (mask, tmpbad)
                    iraf.ccdred.badpiximage (tmpbad, image, n_mask,
                            goodvalue=0, badvalue=1)
                    os.remove (tmpbad)
                    iraf.imarith (n_mask, "+", ref_badpix,
                                  n_mask, title="", divzero=0., hparams="",
                                  pixtype="real", calctype="real",
                                  verbose=False, noact=False)
                else:
                    n_mask = ref_badpix
            else:
                n_mask = ref_badpix

            # transform the bad pixel mask
            if n_mask == ref_badpix:
                self._gl.glogprint ("Transforming %s to %s" %
                                    (n_mask, image_mask))
            else:
                self._gl.glogprint ("Transforming %s + mask to %s" %
                                    (n_mask, image_mask))
            # (imcoadd.cl, line 1561)
            iraf.geotran (n_mask, image_mask, self._database,
                          record_name, geometry="geometric",
                          xin=iraf.INDEF, yin=iraf.INDEF,
                          xshift=iraf.INDEF, yshift=iraf.INDEF,
                          xout=iraf.INDEF, yout=iraf.INDEF,
                          xmag=iraf.INDEF, ymag=iraf.INDEF,
                          xrotation=iraf.INDEF, yrotation=iraf.INDEF,
                          xmin=1, xmax=self._Xmax, ymin=1, ymax=self._Ymax,
                          xscale=1., yscale=1.,
                          ncols=iraf.INDEF, nlines=iraf.INDEF,
                          xsample=self._n_xsamp, ysample=self._n_ysamp,
                          interpolant="linear", boundary="constant",
                          constant=1., fluxconserve=True,
                          nxblock=self._geonxblock, nyblock=self._geonyblock,
                          verbose=self._verbose)

            # convert the mask to integer
            (data, phdr, hdr) = gemutil.getdata (image_mask)
            data += 0.75
            gemutil.putdata (image_mask, data.astype (N.int16), mef=False,
                             phdr=phdr, hdr=None, clobber=True)
            del (data, phdr, hdr)

            gemutil.gemhedit (image_trn, 0, "BPM",
                              image_mask, "Bad pixel mask")

            if n_mask != ref_badpix:  # delete mask if not ref_badpix
                os.remove (n_mask)

        # average of the median sky values
        medsky = sum_medsky / float (self._n_files)

        # make individual mask file for reference image
        # (imcoadd.cl, line 1591)
        if self._fl_msk and self._tmpmli[0] != "none":
            mask = self._tmpmli[0]
            if iraf.imaccess (mask):
                n_mask = mask
                iraf.imarith (ref_badpix, "+", mask, ref_badpix,
                              title="", divzero=0., hparams="", pixtype="real",
                              calctype="real", verbose=False, noact=False)
            elif os.access (mask, os.R_OK):
                tmpbad = iraf.mktemp ("tmpbad")         # text file
                n_mask = iraf.mktemp ("tmpmsk") + ".fits"
                self._gl.glogprint ("Making individual mask for " +
                                    self._images[0])
                trimmed = getTrimmed (mask)
                self._gl.glogprint ("Bad pixel file: %s is for %s images" %
                                    (mask, trimmed))
                if trimmed == "untrimmed":
                    # Copy mask to tmpbad, subtracting TRIMSEC offset
                    subtractTrimsec (mask, tmpbad, self._images[0],
                                     self._sci_ext)
                else:
                    gemutil.copyFile (mask, tmpbad)
                iraf.ccdred.badpiximage (tmpbad, self._images[0], n_mask,
                        goodvalue=0, badvalue=1)
                os.remove (tmpbad)
                iraf.imarith (ref_badpix, "+", n_mask, ref_badpix,
                              title="", divzero=0., hparams="", pixtype="real",
                              calctype="real", verbose=False, noact=False)
                iraf.imdelete (n_mask, verify=False)

        # ---- start of getting intensities of images from objects
        if self._fl_obj:
            self._fl_scale = True        # may be reset to False below
            # The <image>_trn text files contain two pairs of coordinates,
            # x,y in the reference image and x,y in another image, for
            # objects found in each input image.  Extract the coordinates
            # in the reference image for objects that were found in all
            # input images.  In imcoadd.cl, common_stars was the name of a
            # temporary text file (and it was called tmpcen); here it is a
            # list of strings.
            common_stars = self._join_trn_files()

            maglists = []
            for n in range (self._n_files):     # includes reference image
                image_trn = self._trn_images[n]
                image_trn_mag_file = self._rootnames[n] + "_trn_mag"
                # junk is a list of prompt messages (because coords="STDIN")
                junk = iraf.apphot.phot (image_trn, skyfile="",
                        coords="STDIN", output=image_trn_mag_file,
                        plotfile="", interactive=False,
                        radplots=False, icommands="", gcommands="",
                        wcsin="logical", wcsout="logical",
                        cache=False, verify=False, update=")_.update",
                        verbose=")_.verbose", graphics=")_.graphics",
                        display=")_.display", scale=1., fwhmpsf=2.5,
                        emission=True, sigma=iraf.INDEF,
                        datamin=self._datamin, datamax=self._datamax,
                        noise="poisson", ccdread="", gain="",
                        readnoise=self._ron[n], epadu=self._gain[n],
                        exposure="", airmass="", filter="", obstime="",
                        itime=1., xairmass=iraf.INDEF, ifilter="INDEF",
                        otime="INDEF", calgorithm="none", cbox=5.,
                        cthreshold=0., minsnratio=1., cmaxiter=10,
                        maxshift=1., clean=False, rclean=1., rclip=2.,
                        kclean=3., mkcenter=False, salgorithm="constant",
                        annulus=10., dannulus=10., skyvalue=0.,
                        smaxiter=10, sloclip=0., shiclip=0., snreject=50,
                        sloreject=3., shireject=3., khist=3., binsize=0.1,
                        smooth=False, rgrow=0., mksky=False,
                        weighting="constant", apertures=str(aperture),
                        zmag=22., mkapert=False, Stdin=common_stars, Stdout=1)
                del junk

                if os.access (image_trn_mag_file, os.R_OK):
                    mag = iraf.pdump (image_trn_mag_file,
                                      "xcenter,ycenter,mag,perror", "yes",
                                      headers=False, parameters=False,
                                      Stdout=1)
                    maglists.append (mag)
                else:
                    maglists = None
                    self._gl.glogprint ("No good photometry for %s" % \
                                        image_trn_mag_file, type="warning")
                    self._gl.glogprint ("Setting fl_scale=no", type="warning")
                    self._fl_scale = False
                    break
            if maglists:
                try:
                    (nstars, avg_flux) = computeFlux (maglists)
                except RuntimeError, errmsg:
                    self._gl.glogprint (errmsg, type="warning")
                    self._gl.glogprint ("Setting fl_scale=no", type="warning")
                    self._fl_scale = False

            if self._fl_scale:
                for n in range (self._n_files):  # includes reference image
                    image_trn = self._trn_images[n]
                    self._gl.glogprint (
"Mean intensity for   %s:  (N, absolute, relative)= %4d %8.0f %6.3f" %
(image_trn, nstars, avg_flux[n], avg_flux[n]/avg_flux[0]))

                    # Put the relative intensity in the image header.
                    # The images are *multiplied* by this factor, so it has
                    # to be the inverse relative intensity.
                    gemutil.gemhedit (image_trn, 0, "RELINT",
                             avg_flux[0]/avg_flux[n],
                            "Inverse of the relative intensity")

        # ---- end of getting intensities of images from objects
        else:
            self._gl.glogprint ("No objects. Setting fl_scale=no",
                                type="warning")
            self._fl_scale = False

        if self._fl_scale:
            scale_value = "!RELINT"
            scale_message = "imcoadd combine scaled with RELINT"
        else:
            scale_value = "none"
            scale_message = "imcoadd combine not scaled with RELINT"

        # Threshold not used for combining!
        # The median image is on the scale of the reference image.
        log_info = iraf.imcombine ("@STDIN", ref_med,
                        headers="", bpmasks="",
                        rejmasks="", nrejmasks="", expmasks="", sigmas="",
                        logfile="STDOUT", combine="median", reject="none",
                        project=False, outtype="real", outlimits="",
                        offsets="none", masktype="goodvalue", maskvalue=0.,
                        blank=0., scale=scale_value,
                        zero="none", weight="none", statsec="", expname="",
                        lthreshold=iraf.INDEF, hthreshold=iraf.INDEF,
                        nlow=1, nhigh=1, nkeep=1, mclip=True,
                        lsigma=3., hsigma=3.,
                        rdnoise="0.", gain="1.", snoise="0.",
                        sigscale=0.1, pclip=-0.5, grow=0.,
                        Stdin=self._trn_images, Stdout=1)
        self._gl.glogprint ("Output from imcombine:")
        self._gl.glogprint (log_info, type="list")
        datetime = gemutil.gemdate()

        # add sky to median image
        fd = pyfits.open (ref_med, mode="update")
        fd[0].data += medsky
        fd.close()

        gemutil.gemhedit (ref_med, 0, "GEM-TLM", datetime,
                          "UT Last modification with GEMINI")
        gemutil.gemhedit (ref_med, 0, "GEM-TLM", datetime,
                          "UT Time stamp for imcoadd")
        gemutil.gemhedit (ref_med, 0, "COMBSC", scale_message)
        gemutil.gemhedit (ref_med, 0, "MED_SKY", medsky,
                          comment="Median sky level")
        # remove reference to BPM for reference image
        gemutil.gemhedit (ref_med, 0, "BPM", delete=True)

        # Update gain and readnoise keywords
        if gemutil.getkey ("GAINORIG", ref_med, 0, default=0) == 0:
            gemutil.gemhedit (ref_med, 0, "GAINORIG", self._gain[0],
                              "Input gain")
        if gemutil.getkey ("RONORIG", ref_med, 0, default=0) == 0:
            gemutil.gemhedit (ref_med, 0, "RONORIG", self._ron[0],
                              "Input read-noise")

        # Effective readout noise and gain for combined image
        # readout noise in e-/pix  as expected by other iraf tasks
        gemutil.gemhedit (ref_med, 0, self._key_gain, 2.*self._gaineff/3.,
                          "Effective gain [e-/ADU]")
        gemutil.gemhedit (ref_med, 0, self._key_ron,
                          math.sqrt (2.*self._roneff/3.),
                          "Effective read-noise [e-]")

        self.add_image_key (ref_med)

        self._median_done = True

    def _join_trn_files (self):
        """Join _trn files.

        The <rootname[n]>_trn text files contain contain pixel coordinates
        (x_ref, y_ref, x_n, y_n) for targets found in image n.  Extract
        (x_ref, y_ref) from these sets of coordinates for targets that
        were found in all images.  x_ref and y_ref are the coordinates in
        the reference image, so for a given target these coordinates will
        be the same in all images.  The "join" is done on y coordinate only,
        because that's what imcoadd.cl does (line 1645).  The function value
        is a list of strings containing these coordinates.
        """

        SMALL_DIFF = 0.001              # tiny fraction of a pixel

        coords_all_files = []
        # Exclude reference image, because there is no <rootname>_trn file
        # for it, since there wouldn't be any point in "transforming" the
        # reference image to reference image coordinates.
        num_trn_files = self._n_files - 1
        for n in range (1, self._n_files):
            fd = open (self._rootnames[n] + "_trn")     # image_trn_file
            coords = fd.readlines()
            fd.close()
            newcoords = []
            for line in coords:
                words = line.split()
                newcoords.append ([float (words[0]), float (words[1]), n-1])
            coords_all_files.extend (newcoords)

        sorted_coords = N.array (sorted (coords_all_files, cmp=join_cmp3))

        common_stars = []
        # flags will be all true if the current star is found in all trn files.
        flags = N.zeros (num_trn_files, dtype=N.bool)
        prev_x = sorted_coords[0,0]
        prev_y = sorted_coords[0,1]
        n = sorted_coords[0,2]
        flags[n] = True
        for i in range (len (sorted_coords)):
            x = sorted_coords[i,0]
            y = sorted_coords[i,1]
            n = sorted_coords[i,2]
            if abs (y - prev_y) < SMALL_DIFF:
                flags[n] = True
            else:
                if flags.all():
                    common_stars.append ("%8.3f %8.3f" % (prev_x, prev_y))
                prev_x = x
                prev_y = y
                flags[:] = False
                flags[n] = True
        if flags.all():
            common_stars.append ("%8.3f %8.3f" % (x, y))

        return common_stars

    # end of fl_med IF

    def add (self, limit=15., key_limit=True, lowsigma=7., lowlimit=500.,
             scalenoise=0., growthrad=1):
        """Derive the average image, cleaned for cosmic-ray-events."""

        self._gl.glogprint ("Running add()")

        # Get the average of the median sky from header keywords.
        self._getMedianSky()            # assigns self._median_sky

        # Make masks for individual images by comparison to median image,
        # then calculate mean of all images.

        # check range
        if growthrad < 0 or growthrad > 1:
            self._crash (ValueError, "growthrad is out of range")

        ref_add = gemutil.appendSuffix (self._images[0], "_add")
        ref_med = gemutil.appendSuffix (self._images[0], "_med")
        # Read ref_med into memory (ignore the headers).
        (ref_med_data, junk, junk) = gemutil.getdata (ref_med)
        # Subtract median sky from the in-memory data for ref_med.
        ref_med_data -= self._median_sky

        scalenoise2 = scalenoise**2             # Noise scaling ^2 (fraction)
        n_lim = limit                           # may be changed in the loop
        n_files = len (self._trn_images)
        for n in range (n_files):
            # image_trn is called n_inim in imcoadd.cl
            image_trn = self._trn_images[n]
            ron2 = (self._ron[n] / self._gain[n])**2    # (RON in ADU)^2
            # make limit absolute if key_limit=yes
            if key_limit:
                n_lim = limit * \
                        math.sqrt (ron2 + self._median_sky/self._gain[n])
            self._gl.glogprint ("Masking cosmic ray events in " + image_trn)
            if self._fl_scale:
                n_inint = 1. / gemutil.getkey ("RELINT", image_trn, 0,
                                               default=1.)
            else:
                n_inint = 1.

            # sigma image, get sky value
            n_insky = gemutil.getkey ("MED_SKY", image_trn, 0, default=0.)

            # Read image_trn into memory (ignore the headers).
            (image_trn_data, junk, junk) = gemutil.getdata (image_trn)
            # (imcoadd.cl, line 1910)
            image_trn_diff_data = image_trn_data / n_inint - ref_med_data

            n_sigim_data = N.sqrt (( 1. / n_inint**2 + 1.5625 / n_files) *
                 ron2 + ((image_trn_data + n_insky) / n_inint**2 +
                  (ref_med_data + self._median_sky) * 1.5625 / n_files) /
                 self._gain[n] +
                ref_med_data**2 * scalenoise2)

            # mask image  one: deviation, positive only
            #            zero: no deviation
            image_trn_temp = N.where (
                    N.logical_or (image_trn_diff_data > lowsigma * n_sigim_data,
                                  image_trn_diff_data > lowlimit),
                    1, 0)
            n_mask_data = N.where (
                    N.logical_and (ref_med_data <= n_lim,
                                   image_trn_temp == 1),
                    1, 0)
            del (n_sigim_data, image_trn_diff_data, image_trn_temp)

            # if growthrad=1 shift n_mask around and add
            # (imcoadd.cl, line 1940)
            if growthrad > 0:
                self._gl.glogprint ("Applying growth radius %.1f" % growthrad)

                # Make mask 1 and 0 only
                # Make a mask with huge values in badpix
                n_m0_data = N.where (n_mask_data, 1000., 0.)

                # n_mask and n_m0 are input to imcombine, so they need to be
                # images
                # n_m1 is output from imcombine and is ignored
                # n_m2_pl is output from imcombine and is what is used
                n_mask = iraf.mktemp ("tmpmsk") + ".fits"
                n_m0 = iraf.mktemp ("tmpm0") + ".fits"
                n_m1 = iraf.mktemp ("tmpm1") + ".fits"
                n_m2 = iraf.mktemp ("tmpm2")
                n_m2_fits = n_m2 + ".fits"
                n_m2_pl = n_m2 + ".pl"
                # Write data to temporary images.
                gemutil.putdata (n_mask, n_mask_data.astype(N.int16), mef=False)
                gemutil.putdata (n_m0, n_m0_data, mef=False)
                del (n_m0_data)

                # Combine and use output image as the growth
                iraf.imcombine (n_mask+"," + n_m0,
                        n_m1,
                        headers="", bpmasks="", rejmasks="",
                        nrejmasks=n_m2_pl,
                        expmasks="", sigmas="", logfile="",
                        combine="average", reject="minmax", project=False,
                        outtype="real", outlimits="", offsets="none",
                        masktype="none", maskvalue=0., blank=1.,
                        scale="none", zero="none", weight="none", statsec="",
                        expname="",
                        lthreshold=iraf.INDEF, hthreshold=iraf.INDEF,
                        nlow=0, nhigh=1, nkeep=1, mclip=True,
                        lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
                        snoise="0.", sigscale=0.1, pclip=-0.5, grow=1.)

                # Read temporary image back into memory.
                iraf.imcopy (n_m2_pl, n_m2_fits, verbose=False)
                (n_m2_data, junk, junk) = gemutil.getdata (n_m2_fits)

                # Combine the original mask and the growth
                n_mask_data = N.where (
                        N.logical_or (n_mask_data > 0, n_m2_data > 1),
                        1, 0)
                # Delete temporary images
                iraf.imdelete (n_mask+"," + n_m0+"," + n_m1+"," +
                               n_m2_fits+"," + n_m2_pl)
            # (end growthrad > 0)

            # add the mask image to the badpixmask for the image
            # (imcoadd.cl line 1968)
            image_mask = gemutil.appendSuffix (self._images[n], "badpix")
            (image_mask_data, phdr, hdr) = gemutil.getdata (image_mask)
            image_mask_data += n_mask_data
            gemutil.putdata (image_mask, image_mask_data, phdr=phdr, hdr=hdr,
                             clobber=True)
            del (n_mask_data, image_mask_data)
        # end of masking

        # combine the images with masking
        # threshold not used for combining, scale if fl_scale=yes

        if self._fl_scale:
            scale_value = "!RELINT"
            scale_message = "imcoadd combine scaled with RELINT"
        else:
            scale_value = "none"
            scale_message = "imcoadd combine not scaled with RELINT"

        log_info = iraf.imcombine ("@STDIN", ref_add, headers="", bpmasks="",
                        rejmasks="", nrejmasks="", expmasks="",
                        sigmas="", logfile="STDOUT", combine="average",
                        reject="none", project=False, outtype="real",
                        outlimits="", offsets="none",
                        masktype="goodvalue", maskvalue=0., blank=0.,
                        scale=scale_value, zero="none", weight="none",
                        statsec=self._statsec, expname="",
                        lthreshold=iraf.INDEF, hthreshold=iraf.INDEF,
                        nlow=1, nhigh=1, nkeep=1,
                        mclip=True, lsigma=3., hsigma=3.,
                        rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1,
                        pclip=-0.5, grow=0., Stdin=self._trn_images, Stdout=1)
        self._gl.glogprint ("Output from imcombine:")
        self._gl.glogprint (log_info, type="list")

        if self._fl_mef:
            # Pack up as MEF if input is MEF, using the headers from the
            # reference image.
            gemutil.convertMEF (ref_add, output=ref_add,
                                extname="SCI", template=self._images[0])

        # Update the header, and add the median (sky) back to ref_add.
        fd = pyfits.open (ref_add, mode="update")
        phdr = fd[0].header                     # primary header
        header = fd[self._sci_ext].header       # extension or primary header
        data = fd[self._sci_ext].data

        datetime = gemutil.gemdate()
        phdr.update ("GEM-TLM", datetime, "UT Last modification with GEMINI")
        phdr.update ("IMCOADD", datetime, "UT Time stamp for imcoadd")
        phdr.update ("IMCOBPM", self._badpixfile, "Detector bad pixel file")
        phdr.update ("COMBSC", scale_message)

        # remove reference to BPM for reference image
        del (header["BPM"])

        # Add the average of the median sky levels back to the image.
        data += self._median_sky
        phdr.update ("MED_SKY", self._median_sky)

        # various header info in n_ref // "_add"

        # Input readout noise and gain
        # Update gain and readnoise keywords
        if not header.has_key ("GAINORIG") or header["GAINORIG"] == 0:
            header.update ("GAINORIG", self._gain[0], "Input gain")
        if not header.has_key ("RONORIG") or header["RONORIG"] == 0:
            header.update ("RONORIG", self._ron[0], "Input read-noise")

        # Effective readout noise and gain for combined image
        # readout noise in e-/pix  as expected by other iraf tasks
        # self._ron is currently RON^2 in ADU/pixel
        header.update (self._key_gain, self._gaineff,
                       "Effective gain [e-/ADU]")
        header.update (self._key_ron, (math.sqrt (self._roneff)),
                       "Effective read-noise [e-]")

        if self._sci_ext != self._zero_ext:
            # Also update keywords in primary header.
            if not phdr.has_key ("GAINORIG") or phdr["GAINORIG"] == 0:
                phdr.update ("GAINORIG", self._gain[0], "Input gain")
            if not phdr.has_key ("RONORIG") or phdr["RONORIG"] == 0:
                phdr.update ("RONORIG", self._ron[0], "Input read-noise")
            phdr.update (self._key_gain, self._gaineff,
                         "Effective gain [e-/ADU]")
            phdr.update (self._key_ron, (math.sqrt (self._roneff)),
                         "Effective read-noise [e-]")

        phdr.update ("SCNOISE", scalenoise, "imcoadd scnoise")
        phdr.update ("LIMIT", n_lim, "imcoadd limit for cleaning")
        phdr.update ("LOWSIG", lowsigma, "imcoadd sigma rejection")
        phdr.update ("LOWLIM", lowlimit, "imcoadd absolute rejection limit")
        phdr.update ("GROWTHR", growthrad, "imcoadd growth radius")
        phdr.update ("GEOGEOM", self._geofitgeom, "imcoadd geomap geometry")
        phdr.update ("GEOINTER", self._geointer,
                     "imcoadd geotrans interpolation")

        fd.close()
        self.add_image_key (ref_add)

        # end of fl_add

    def average (self):
        """Derive the average uncleaned image using IMCOMBINE."""

        if not self._transform_done:
            self.transform (rerun=False)

        self._gl.glogprint ("Running average()")

        # Get the average of the median sky from header keywords.
        self._getMedianSky()            # assigns self._median_sky

        ref_avg = gemutil.appendSuffix (self._images[0], "_avg")

        if self._fl_scale:
            scale_value = "!RELINT"
            scale_message = "imcoadd combine scaled with RELINT"
        else:
            scale_value = "none"
            scale_message = "imcoadd combine not scaled with RELINT"
        # Combine the _trn images.
        log_info = iraf.imcombine ("@STDIN", ref_avg, headers="",
                        bpmasks="", rejmasks="", nrejmasks="", expmasks="",
                        sigmas="", logfile="STDOUT",
                        combine="average", reject="none", project=False,
                        outtype="real", outlimits="", offsets="none",
                        masktype="none", maskvalue=0., blank=0.,
                        scale=scale_value, zero="none",
                        weight="none", statsec=self._statsec, expname="",
                        lthreshold=iraf.INDEF, hthreshold=iraf.INDEF,
                        nlow=1, nhigh=1, nkeep=1,
                        mclip=True, lsigma=3., hsigma=3.,
                        rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1,
                        pclip=-0.5, grow=0., Stdin=self._trn_images, Stdout=1)
        self._gl.glogprint ("Output from imcombine:")
        self._gl.glogprint (log_info, type="list")

        if self._fl_mef:
            # Pack up as MEF if input is MEF, using the headers from the
            # reference image.
            gemutil.convertMEF (ref_avg, output=ref_avg,
                                extname="SCI", template=self._images[0])

        # Update the header, and add the median (sky) back to ref_avg.
        fd = pyfits.open (ref_avg, mode="update")
        phdr = fd[0].header                     # primary header
        header = fd[self._sci_ext].header
        # Add the average of the median sky levels back to the image.
        fd[self._sci_ext].data += self._median_sky
        data = fd[self._sci_ext].data

        datetime = gemutil.gemdate()
        phdr.update ("GEM-TLM", datetime, "UT Last modification with GEMINI")
        phdr.update ("IMCOADD", datetime, "UT Time stamp for imcoadd")
        phdr.update ("COMBSC", scale_message)
        # remove reference to BPM for reference image
        del (header["BPM"])

        phdr.update ("MED_SKY", self._median_sky)

        # Input readout noise and gain
        # Update gain and readnoise keywords
        if not header.has_key ("GAINORIG") or header["GAINORIG"] == 0:
            header.update ("GAINORIG", self._gain[0], "Input gain")
        if not header.has_key ("RONORIG") or header["RONORIG"] == 0:
            header.update ("RONORIG", self._ron[0], "Input read-noise")

        # Effective readout noise and gain for combined image
        # readout noise in e-/pix  as expected by other iraf tasks
        header.update (self._key_gain, self._gaineff,
                       "Effective gain [e-/ADU]")
        header.update (self._key_ron, (math.sqrt (self._roneff)),
                       "Effective read-noise [e-]")

        fd.close()
        self.add_image_key (ref_avg)

        # end of fl_avg

    def add_image_key (self, filename):
        """Add IMAGEi = input file keywords to primary header."""

        fd = pyfits.open (filename, mode="update")
        phdr = fd[0].header                     # primary header
        for n in range (self._n_files):         # includes reference image
            image = self._images[n]
            # leading zero for n+1 = 1-9.  can n+1 be > 99 (should be OK)?
            keyword = "IMAGE%02d" % (n+1,)
            phdr.update (keyword, image, "Input image for imcoadd")
        fd.close()

    def close (self):
        """Delete temporary files and images; close log file."""

        if self._immasks == "DQ" and self._tmpmli is not None:
            for maskname in self._tmpmli:
                gemutil.deleteFile (maskname)
        # xxx etc.

        self._gl.glogclose (fl_success=True)


    def _getMedianSky (self):
        """Read MED_SKY from _trn headers, and compute average."""

        if self._median_sky is not None:
            return

        self._median_sky = 0.
        for image_trn in self._trn_images:
            self._median_sky += gemutil.getkey ("MED_SKY", image_trn, 0,
                                                must_exist=True)
        self._median_sky /= float (self._n_files)

        if self._verbose:
            print "Median sky: ", repr (self._median_sky)

    def _crash (self, errtype=RuntimeError, errmess=""):

        if self._gl is not None:
            self._gl.glogprint (errmess, loglevel="status", type="error")
            self._gl.glogclose (fl_success=False)

        self.close()

        if errmess:
            raise errtype, errmess
        else:
            raise errtype

    def _deleteOutputFiles (self, filelist):
        """Delete output files if they already exist.

        filelist is a list of files that the calling routine is about to
        create.  If _fl_overwrite is True, all files in filelist that
        already exist will be deleted.  If _fl_overwrite is False, an
        exception will be raised if any file in filelist already exists.

        @param filelist: the names of the files to be deleted
        @type filelist: a single string or a list of strings
        """

        if isinstance (filelist, str):
            filelist = [filelist]

        if self._fl_overwrite:
            for filename in filelist:
                gemutil.deleteFile (filename)
        else:
            existing_files = []
            for filename in filelist:
                if os.access (filename, os.F_OK):
                    existing_files.append (filename)
            if len (existing_files) > 0:
                if len (existing_files) == 1:
                    errmess = "Output file %s already exists" % \
                          existing_files[0]
                else:
                    errmess = "The following output files already exist:  " + \
                              repr (existing_files)
                self._crash (errmess=errmess)

def join_cmp3 (a, b):
    """Compare a and b for sorting.

    a and b are three-element lists.  The elements are the x and y pixel
    coordinates of a star and the number of the image in which the star
    was found.  This function results in a sort on x, y, image number,
    in that order.

    @param a: first item to compare
    @type a: three-element list
    @param b: second item to compare
    @type b: three-element list
    """

    if a[0] < b[0]:             # x coordinate
        return -1
    elif a[0] > b[0]:
        return 1
    elif a[1] < b[1]:           # y coordinate
        return -1
    elif a[1] > b[1]:
        return 1
    elif a[2] < b[2]:           # image number
        return -1
    elif a[2] > b[2]:
        return 1
    else:
        return 0

def getTrimmed (badpixfile):
    """Read 'trimmed' or 'untrimmed' from the comment line in badpixfile.

    If the bad pixel file is a text file, the first line should be
    of the form "# trimmed" or "# untrimmed" followed by a comment.
    This function returns the word after the "#" sign.

    @param badpixfile: bad pixel file name
    @type badpixfile: string

    @return: should be either 'trimmed' or 'untrimmed'
    @rtype: string
    """

    fd = open (badpixfile)
    line = fd.readline()
    fd.close()
    words = line.split()
    if words[0] == "#":
        trimmed = words[1]
    else:
        trimmed = "???"

    return trimmed


def subtractTrimsec (badpixfile, newbadpixfile, image, extension=0):
    """Subtract the TRIMSEC offset from the elements of a bad pixel file.

    @param badpixfile: bad pixel file name
    @type badpixfile: string

    @param newbadpixfile: copy of badpixfile, with dx or dy offset
        subtraced from each element
    @type newbadpixfile: string

    @param image: name of file containing image with keyword TRIMSEC
    @type image: string

    @param extension: number or name of FITS extension containing image
    @type extension: number, string, or tuple, e.g. ("SCI", 1)
    """

    section = gemutil.getkey ("TRIMSEC", image, extension=extension)
    (x0, x1, y0, y1) = irafutil.splitSection (section)

    bp_list = gemutil.fieldsOfTable (badpixfile, "0,1,2,3")

    new_bp_list = []
    for line in bp_list:
        words = line.split()
        first_x = float (word[0]) + 1. - x0
        last_x  = float (word[1]) + 1. - x0
        first_y = float (word[2]) + 1. - y0
        last_y  = float (word[3]) + 1. - y0
        new_bp_list.append ("%g %g %g %g\n" %
                            (first_x, last_x, first_y, last_y))

    ofd = open (newbadpixfile, mode="w")
    ofd.writelines (new_bp_list)
    ofd.close()


def computeFlux (maglists):
    """Compute the average intensity from the magnitudes.

    @param maglists: One element for each imcoadd input image (including the
        reference image).  Each element is a list, the output of pdump for
        the apphot.phot results.  Each list should be the same length (the
        length is the number of stars), and each should have the same four
        columns:  xcenter, ycenter, mag, perror.
    @type maglists: list of lists

    @return: each element is the average of the flux for the stars that
        had valid photometry; the flux is computed from mag assuming a
        reference magnitude of 22.
    @rtype: array of float64, or None in case of error
    """

    n_files = len (maglists)

    # Make sure every list in maglists has the same length (the same
    # number of stars).
    nstars = len (maglists[0])
    for n in range (n_files):
        if nstars != len (maglists[n]):
            raise RuntimeError, \
"There should be the same number of stars (for photometry) in each image."

    # avg_flux will be returned
    avg_flux = N.zeros (n_files, dtype=N.float64)

    # Some stars in some images may have poor photometry, in which case
    # we have to exclude that star when computing the average flux for
    # every image.  labels is an array of flags, 1 is OK.
    labels = N.ones (nstars, dtype=N.int32)

    # Create a list of 1-D arrays for the magnitudes.  We need to check
    # the 'perror' for all data before we can compute averages.
    fluxarrays = []
    for n in range (n_files):
        photometry = maglists[n]
        flux = N.zeros (nstars, dtype=N.float64)
        for i in range (nstars):
            # x, y, mag, error flag
            xyme = photometry[i].split()
            # mag is OK and perror does not indicate an error?
            if xyme[2] == "INDEF":
                mag = 9999.+1.          # flag as bad
            else:
                mag = float (xyme[2])
            if mag < 9999. and xyme[3] == "NoError":
                flux[i] = 10.**(0.4 * (22.-mag))
            else:
                labels[i] = 0
        fluxarrays.append (flux)

    n_good_stars = labels.sum()
    if n_good_stars < 1:
        raise RuntimeError, "No good photometry for any star."

    # now compute the average of the fluxes for each image
    for n in range (n_files):
        avg_flux[n] = ndimage.mean (fluxarrays[n], labels=labels)

    return (n_good_stars, avg_flux)
