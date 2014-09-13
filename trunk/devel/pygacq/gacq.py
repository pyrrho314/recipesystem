#!/usr/bin/env python
# Copyright(c) 2000-2013 Association of Universities for Research in Astronomy, Inc.

from __future__ import print_function

import os, sys
import datetime
import logging
import argparse

# science stuff
import math
import numpy as np
import scipy
import scipy.stats

# other gemini modules
from astrodata import AstroData

# general gacq modules
import userinterface as ui
from gacqlogging import *
from acquisitionimage import AcquisitionImage
from util import angle, rotate, clamp, is_number

# single object selection modules
from selection import mark_selection, SelectionCursorListener

# long slit modules
from longslit import SlitMeasurementCursorListener, find_slit, get_slit_tilt, draw_slit

# MOS modules
from multiobjectspectroscopy import InteractiveMosaic, MOSBoxCursorListener
from superimpose import superimpose

# IFU modules
from integratedfieldunit import IFUImaging, get_integrated_field_unit_offsets

__version__ = "2013 July 24" # Software version

#-----------------------------------------------------------------------
# Gemini acquisition script
#-----------------------------------------------------------------------
# 
# 2005 Nov 05  IJ  alpha version of generic setup script using gnsetup.cl as base
# 2005 Nov 07  KL  fix indentation, removed unused variables, standardize params
#                  use gemisnumber instead of CL trick, removed unused l_struct
#                  replaced series of 'if' with more efficient 'if-else if'
#                  make l_z1 and l_z2 real. change defaults to INDEF
#                  fixed typo, missing var, and added geomap error handling structure
# 2005 Dec 01  AS  adjust IFU position for Kathy
# 2005 Dec 18  AS  merge with nssetup.cl
# 2005 Dec 20  AS  more fix-ups, add version parameter for user's benefit
# 2005 Dec 21  AS  do not draw slit/ifu if a mask is in
# 2005 Dec 22  AS  set display zscale of slit images properly
# 2005 Dec 23  AS  add human error handler
# 2005 Dec 26  AS  consolidate showslit + showifu into "showoverlay" aka "show"
# 2005 Dec 27  AS  mark center of NIRI frame when show- and no mask in beam
# 2005 Dec 28  AS  draw GMOS slits vertically to match new slit masks
# 2005 Dec 29  AS  flip NIRI f/14 offsets (x & y)
# 2006 Jan 03  AS  set slit measurement rplot = slitwidth
# 2006 Jan 04  AS  implement "observatory" gemlocal parameter (compliments of PG)
# 2006 Jan 05  AS  draw compass on image for Percy
# 2006 Jan 06  AS  fix bug so as to not try to measure IFU; add gflush to compass
# 2006 Jan 07  AS  smarter fitting of slit profile (check that fit sigma is reasonable)
# 2006 Jan 09  AS  fix bug finding slit image specified on command line
# 2006 Jan 10  AS  fix bug in slit image display
# 2006 Jan 22  AS  merge with gsmaskalign.cl
# 2006 Jan 24  AS  merge with gfquick.cl
# 2006 Jan 28  AS  merge with mssetup.cl
# 2006 Jan 29  AS  output RMS of mask alignment fit
# 2006 Jan 30  AS  add catch for MDFs without ODFXBIN/ODFYBIN keywords
# 2006 Jan 31  AS  hardcode michelle to use science extension 1
# 2006 Feb 01  AS  add color option for idiosyncratic observers
# 2006 Feb 17  AS  fix up twotarget mode for NIRI
# 2006 Mar 01  KV  Michelle slit position adjustments
# 2006 Mar 21  AS  changes for GMOS-S IFU
# 2006 Apr 04  AS  fix MOS mask coord file bug, remove vestigial "paste" flag
# 2006 Apr 05  AS  if only one "click" assume object goes to the center of slit (gmos/niri)
# 2006 Apr 11  AS  fix bug with auto centering NIRI imaging field
# 2006 Apr 23  AS  fix MOS mask bug
# 2006 Apr 25  AS  fix bug with getting slit position with only one click on unbinned image
# 2006 May 09  AS  extend smart parsing of positions to twotarget mode (including auto-slit-centering)
# 2006 May 18  AS  fix NIRI f/14 + Altair p/q sign problem
# 2006 May 29  RD  return MOS IAA for GMOS-S and IPA for GMOS-N
# 2006 May 29  AS  check that MOS acquisition image is binned 1x1
# 2006 Jun 19  AS  fix MOS mask name formatting and improve error messages
# 2006 Jun 19  AS  add default slit choice for NIRI + Michelle w/o slit image
# 2006 Jun 19  AS  add Michelle IAA error check & M2 baffle check
# 2006 Jun 19  AS  fill in some missing ybin factors
# 2006 Jul 10  RD  made compass drawing optional and UT filename correction for GS
# 2006 Jul 24  RD  begin merging GNIRS
# 2006 Jul 27  AS  continue merging GNIRS, add GNIRS imaging auto z-scaling
# 2006 Jul 28  AS  add TEXES mode
# 2006 Aug 04  AS  fix mdf path bug, update UT filename correction for GS
# 2006 Aug 06  AS  assimilate NIFS nfimage
# 2006 Aug 08  AS  add forgotten variable definitions (l_nifs, l_flip)
# 2006 Aug 25  AS  add capability for NIRI subarray acquisitions
# 2006 Aug 25  AS  fix bug with yslitpos when no slit image is measured
# 2006 Aug 26  AS  update default GMOS-N slit/IFU positions
# 2006 Sep 01  BM  fix GNIRS support and GS date at month/yr transitions
# 2006 Sep 03  AS  Update Michelle IAA, remove unused code
# 2006 Sep 05  BM  fix neg. slit image bug, add default GNIRS slit position,
#                  fix auto slit image name generation, add prefix parameter 
#                  for running on old data
# 2006 Sep 18  AS  fix bug with NIRI subarray imaging
# 2006 Sep 18  AS  prepare code for sending offsets to TCC
# 2006 Sep 18  AS  fix Michelle imaging Q-offset sign error
# 2006 Sep 19  AS  send offsets to TCC
# 2006 Oct 17  PG  change mode of GNIRS slit acq so that XD is more like a GMOS
#                  mask and the longslit acq only needs one measurement. Also 
#                  cleaned up unneeded gnirs parameters. However, gacq DOES call nfquick.
# 2006 Oct 30  PG  fixed a bug when fl_slit- was called. Now the script writes a
#                  file with the previous GNIRS slit parameters that can be reused.
# 2006 Nov 17  AS  add ability to selectively send x/y/all offsets to TCC
# 2006 Nov 17  AS  assimilate NIFS nfacquire
# 2006 Nov 24  PG  updated GMOS-N and GMOS-S longslit offsets based on filter
# 2006 Nov 25  AS  update GMOS-N slit position & offset values
# 2006 Nov 25  AS  add capability to use GMOS stamp for slit measurement
# 2006 Nov 27  AS  decrease histogram bin size to 0.001 when finding slit center
# 2006 Nov 27  AS  remove ability to selectively send components of offsets
# 2006 Nov 27  AS  fix GMOS stamp image binning bug
# 2006 Nov 30  AS  add ability to click on NIFS target and desired position
# 2006 Dec 03  AS  small bug-fixes/improvements for GMOS stamp acquisitions
# 2006 Dec 05  AS  enable GMOS stamp IFU acquisitions
# 2006 Dec 06  AS  correct for weird GMOS stamp binning feature
# 2006 Dec 13  AS  additional NIFS grating offsets
# 2006 Dec 19  PG  fix the criterion for selecting acq holes for gmos
# 2006 Dec 19  AS  fix MOS mask name parsing for DD and SV programs (DD is untested)
# 2007 Feb 17  KL  added missing closing parenthesis, & replaced a hselect with a keypar
# 2007 Feb 21  PG  fix GMOS-S overscan regions for MOS acquisition
# 2007 Feb 23  KL  additional fixes for PyRAF compatibility
# 2007 Feb 24  AS  fix MOS mask name parsing for classical programs
# 2007 Mar 02  AS  make Michelle acquisitions the same as GMOS & NIRI 
# 2007 Mar 04  KV  update Michelle slit positions
# 2007 Apr 17  AS  fix GMOS classical MOS name, add warning for poorly aligned MOS masks
# 2007 Apr 18  AS  name MOS coo file with mask name, and if exists ask to use as input
# 2007 May 14  AS  draw outline of MOS acquisition box if showoverlay = true
# 2007 May 16  AS  fix minor bug drawing MOS acq box when near edge of chip
# 2007 Jun 07  AS  add separate variable for answer of whether to send coords to tcc
# 2007 Jun 07  AS  add ticks which show region used to measure slit center
# 2007 Jun 07  AS  manual slit measurement now uses the last input value instead of first
# 2007 Jun 17  AS  make showoverlay- measure slit position
# 2007 Jun 17  AS  fix minor bug drawing slit measurement region
# 2007 Jul 05  RM  update Michelle spectroscopic hotspots
# 2007 Oct 17  AS  fix IAA and IPA so that they are always in the range 0-360 deg
# 2007 Nov 20  AS  send formatted string to gacq2tcc to prevent scientific notation
# 2007 Dec 31  AS  fix NIFS reconstructed image size bug
# 2008 Feb  1  JH  Set parameters in imstat and display calls. Add instrument switch check
# 2008 Feb 15  JT  Fix small but unpleasant typo
# 2008 Mar 10  AS  Add check for >1 amplifier per GMOS chip
# 2008 Mar 17  AS  Add check for twotarget combined with GMOS stamp or NIRI subarrays
# 2008 Mar 24  AS  Always ask which NIRI slit to use (even if not showing it)
# 2008 Apr 03  AS  Improve error message when GMOS MDF not found.
# 2008 Apr 10  AS  Check age of file & warn user if older than 5 min
# 2008 May 16  KV  Read the Michelle "IMPA" header keyword for imaging acquisitions
# 2008 May 30  AS  Add infrastructure for TReCS
# 2008 May 31  AS  Improve NIFS flip-mirror instructions
# 2008 Jun 13  KV  Update Michelle hot spots, adding 7.6 pix to all x values
# 2008 Jun 17  AS  Add support for NICI
# 2008 Jun 19  AS  Fix file age check, add missing imexam frame definitions
# 2008 Jun 20  AS  Define uparm to tmp$ to not save param changes, fix 2nd NIFS image size bug
# 2008 Jul 15  AS  NICI: set extension to 1, auto name flat, update default mask position
# 2008 Sep 09  AS  Update default NICI mask positions
# 2008 Oct 05  AS  NICI UI upgrade: search for flat with same mask & offset
# 2008 Nov 25  JT,KL  NICI fixed Cass mode + PyRAF compatibility.
# 2009 Feb 19  AS  Fix to allow combining subarray sizes during NIRI acquisitions
# 2009 Apr 10  AS  Undo uparm trick (doesn't work in pyraf); save local copies instead
# 2009 May 06  AS  Complete removal of tmp uparm and add compass to NIFS image
# 2009 Jul 28  AS  Set default boxsize to 300 pixels to help find stars with large offsets
# 2009 Aug 02  AS  Search for MDF file in ./ then mdfdir then mdfdir/semester (taken from the program ID)
# 2009 Aug 26  AS  Check for correct binning
# 2009 Sep 17  AS  Redirect NIFS coronagraphic observations to use Richard's IDL procedure
# 2009 Sep 19  AS  Minor NICI changes for Fredrik
# 2009 Sep 21  KL  Michelle hotspot changes
# 2009 Nov 09  AS  Allow user to reject MOS alignment boxes
# 2009 Nov 17  AS  Include large offset warning
# 2009 Nov 20  AS  Clarify one- and two-target acquisition instructions
# 2009 Nov 20  AS  Add option to reject MOS alignment star after measuring box
# 2009 Dec 19  AS  Bugfix for MOS alignment star rejection
# 2009 Dec 22  AS  Bugfix for expected position warning when using NIRI subarray & GMOS MOS
# 2010 Jan 23  AS  Input TReCS default slit positions
# 2010 Jan 30  AS  Clarify slit measurement instructions
# 2010 Feb 27  AS  Implement GMOS menu when acquisition type cannot be deduced
# 2010 Mar 10  AS  Bugfix for TReCS call to mireduce with rawdata specified
# 2010 Mar 11  AS  Minor changes to the GMOS menu
# 2010 Apr 01  KL  Update NIFS shifty by +4 pixels as measured by R.McDermid
# 2010 May 15  AS  Fix: MOS slit width parsing, MOS w/o priority-1 slits, MOS bstar-
# 2010 Jun 01  AS  GNIRS integration & remember previously measured slit position
# 2010 Jun 14  AS  Fix naming of newbox for MOS & fix bug passing coo files
# 2010 Jun 28  AS  Fix slit position memory for GMOS stamps
# 2010 Jul 23  AS  GNIRS changes for GN
# 2010 Aug 19  AS  GNIRS + Altair mods
# 2010 Sep 10  AS  GNIRS: extract slit x-position from slit measurement keystroke
# 2010 Sep 19  AS  Slit memory binning fix and give warning if slit image is binned
# 2010 Nov 25  AS  GNIRS procedure change
# 2010 Dec 03  AS  Fix GNIRS 2-target mode
# 2010 Dec 10  AS  Measure GNIRS slit tilt in 2-target mode
# 2010 Dec 11  AS  Remove NIRI spec, catch GMOS central-spectrum slit image
# 2010 Dec 16  AS  Propagate the GNIRS WCS CD matrix from [0] to [1] if compass+
# 2011 Jan 16  AS  Allow GMOS MOS acquisitions to be binned if the mask is not in the beam
# 2011 Feb 18  AS  Fix GNIRS 2-target + XD slit, reverse sign of P for GNIRS on side port + Altair
# 2011 Feb 19  AS  Add GNIRS imaging mode
# 2011 Feb 19  AS  Add warning for large MOS rotations with only 2 alignment stars
# 2011 Feb 25  AS  Consolidate (x,y) -> (p,q) conversion code for one- and two-target modes
# 2011 Feb 25  AS  Include GNIRS slit tilt in (x,y) -> (p,q) conversion
# 2011 Mar 23  AS  Reverse signs of NIRI port-1 f/14 offsets
# 2011 Mar 25  AS  Reverse sign of P for NIRI + Altair port-1 offsets
# 2011 Mar 25  AS  Fix phase error of two-target rotation
# 2011 Mar 26  AS  Tweak GNIRS slit tracing parameters
# 2011 Apr 09  AS  Display GMOS image first
# 2011 Apr 11  AS  Undo hack for shifted GMOS stamp
# 2011 Apr 26  AS  Write log and check for undesired telescope offsets
# 2011 Apr 30  AS  Rename parameters slit{x,y}pos to {x,y}slitpos (so you can say slit=123)
# 2011 May 04  AS  Allow GNIRS slit identification to occur on the left side of the slit.
# 2011 May 04  AS  Measure the x-center of the GNIRS XD slit.
# 2011 Jun 05  AS  Add option for gacquisition advice
# 2011 Jun 11  AS  Fix bug searching for slits when rawpath is given explicitly
# 2011 Jun 16  AS  Add a floor of 0.02 arcsec for advising to take an additional image
# 2011 Jun 17  AS  Tweak GNIRS slit tracing parameters (to compensate for correct slit width parsing)
# 2011 Jul 01  AS  Add more information to advice
# 2011 Jul 01  AS  Increase slit search from +/- 10 pix to +/- 15 pix, add slit measurement filter
# 2011 Jul 02  AS  Replace apall with apfind + aptrace for finer control
# 2011 Jul 24  AS  Include slit-memory for GNIRS two-target mode
# 2011 Aug 03  JT  Merge local changes for GMOS on the bottom port
# 2011 Aug 11  AS  Add warning when GNIRS slit is close to edge of detector
# 2011 Aug 13  AS  MOS advice now takes into account slitlet size
# 2011 Aug 13  AS  Only give poor MOS alignment warning when mask is in the beam
# 2011 Aug 26  AS  Use detector center as slit reference
# 2011 Sep 03  AS  Fix bug when object is off the bottom of the GMOS slit image
# 2011 Sep 05  AS  Implement Mischa's auto MOS box finding algorithm
# 2011 Sep 12  AS  Reinstate poor MOS alignment warning when mask is out of beam
# 2011 Sep 18  AS  Change GNIRS imaging acq ref point so that it works for the long camera
# 2011 Sep 18  AS  Measure GNIRS keyhole position from slitimage
# 2011 Sep 20  AS  Add blind-offset option which modifieds advice, add missing usecoo
# 2011 Sep 22  AS  Add GNIRS advice on when to take the through-slit image
# 2011 Oct 30  AS  Add flag for marking saturated pixels
# 2011 Nov 06  AS  Add onaxis flag to force sending GMOS Q-offset
# 2011 Nov 10  AS  Disable advice for Michelle since slit not configured until science starts
# 2011 Nov 15  AS  Tweak GNIRS slit tracing parameters
# 2011 Nov 21  AS  Changes for new GMOS e2v detectors in 6-amp mode
# 2011 Nov 26  IJ  GMOS-N 6-amp improvements
# 2011 Nov 29  AS  Modify the default GN e2v IFU-R position to Y=2287.6 as reported by Kristin
# 2011 Dec 02  EH  Starting to add support for F2
# 2011 Dec 08  AS  Additional modifications for F2
# 2011 Dec 16  KR  GMOS-N 6-amp mode adjustment for IFU-2 (Kathy & Kristin)
# 2011 Dec 21  AS  Fix on-axis bug
# 2011 Dec 30  AS  Modify GNIRS 2-pixel through-slit advice
# 2012 Jan 04  AS  F2 MOS
# 2012 Jan 15  AS  Fix bug with F2 off-axis long-slit acquisitions
# 2012 Jan 19  AS  Add blind-offset message to IFU
# 2012 Jan 22  AS  Implement gmultiamp / goversub for real GMOS overscan subtraction
# 2012 Jan 27  AS  Remove old GMOS-N IFU hard-coded overscan subtraction
# 2012 Jan 27  AS  Tweak GMOS-N IFU reconstruction parameters
# 2012 Feb 05  AS  Read GMOS MDF mm and convert to pix
# 2012 Feb 11  AS  Add saturation overlay to reconstructed GMOS IFU images
# 2012 Feb 16  AS  Fix bug searching for tiled slit image
# 2012 Feb 17  AS  Only allow off-axis for long-slit
# 2012 Feb 18  AS  Fix bug with GMOS-N 6-amp stamp acquisitions
# 2012 Feb 20  AS  Fix F2 mask name and saturation overlay
# 2012 Mar 15  AS  Check object name for blind offsets
# 2012 Mar 16  AS  Remove extra "$" in slit image raw path
# 2012 Mar 19  AS  Process slit image discovered after displaying image with gmultiamp
# 2012 Mar 21  AS  Fix the P,Q signs for NIRI f/14 w/o AO
# 2012 Mar 25  AS  Merge GMOS-N & GMOS-S IFU acquisition code
# 2012 Mar 26  AS  Get the correct GMOS saturation value from gsat
# 2012 Mar 30  AS  Change two-target advice to use |P+R| like MOS
# 2012 Apr 15  AS  Improve imaging advice
# 2012 Apr 22  AS  Check MOS mask name length.  Minor advice tweaks.
# 2012 May 18  AS  Check that M2 baffles are correct for GMOS
# 2012 May 23  AS  Support GMOS custom ROI MOS acquisitions
# 2012 May 24  AS  Subtract GMOS 6-amp bias offset from central 2 columns
# 2012 Jun 04  AS  Apply rough GMOS-N gain correction
# 2012 Jun 06  AS  Fix bug with x/yslitpos and handle extra characters in TReCS slit name
# 2012 Jun 08  AS  Automatically find sky frames
# 2012 Jun 15  AS  Update off-axis instructions
# 2012 Jun 23  AS  Fix bugs in MOS custom ROI
# 2012 Jun 24  AS  Update GNIRS two-target advice
# 2012 Jul 02  AS  Fix bug in GMOS-S MOS box position calculation
# 2012 Jul 06  AS  Return the EXTVER and INHERIT keywords to bias & gain corrected GMOS images
# 2012 Jul 11  AS  Fix crash when GMOS off-axis acq has slit=""
# 2012 Jul 15  AS  Warn when both offsets are zero
# 2012 Jul 15  AS  Modify MOS advice to have a maximum offset of 0.2 arcsec
# 2012 Jul 29  AS  Fix bug when full-frame binned MOS acq box is too close to edge of detector
# 2012 Jul 29  AS  Determine MDF semester from mask name not program ID (for when a program uses an old mask)
# 2012 Jul 30  AS  Give better instructions for manually specifying a mask name
# 2012 Sep 12  AS  Fix two-target bug.  Make GMOS two-target mode work for all ROIs.  Reduce off-axis limit to 2 arcsec.
# 2012 Sep 18  AS  Allow skipping MOS acq boxes when there are more than GMOS allows (4)
# 2012 Sep 22  AS  Handle merged custom ROIs
# 2012 Sep 24  AS  Fix bug with custom ROIs which cross amplifiers (and have wrong ROI header values)
# 2012 Oct 15  AS  Support a single custom ROI
# 2012 Oct 22  AS  Fix minor bug/typo in GNIRS two-target and improve F2 off-axis instructions
# 2012 Nov 08  AS  Catch boxes too close to edge of binned MOS acquisition & fix bug with binned MOS acquisitions
# 2012 Nov 14  AS  Fix GMOS-S binned MOS acq bug
# 2012 Dec 16  AS  Add logging
# 2012 Dec 24  AS  Hard-code rplot and naver for NIFS flip-mirror in acquisition.  Improvements to logging.
# 2013 Jan 22  AS  New gmultiamp
# 2013 Mar 11  AS  Fix bugs when giving ifuxpos parameter
# 2013 Apr 16  MS  IRAF 2.16 updates; specify table extension names use ".fits"
# 2013 May 01  BC  Using version control software

#-----------------------------------------------------------------------
# To-Do:
# - if binned 1x1 and full-frame, probably MOS acquisition, so ask for mos mask number
#-----------------------------------------------------------------------

class GACQMode(object):
    def __init__(self, args):
        self.args = args
        self.slitimage = "auto" # Slit image to use for measuring slitpos
        self.longslit = False
        self.tryagain = False
        self.measureslit = False
        self.display = True
        self.showslit = True

    def get_site_prefix(self):
        if self.args.gemininorth:
            prefix = "N"
        elif self.args.geminisouth:
            prefix = "S"
        else:
            prefix = "N"

        utdate = datetime.datetime.utcnow().strftime("%Y%m%d")
        
        return prefix + utdate + "S"

    def set_display_image(self, display):
        self.display = display

    def get_display_image(self):
        return self.display

    def set_show_slit(self, showslit):
        self.showslit = showslit

    def get_show_slit(self):
        return self.showslit

    def get_slitimage(self):
        return self.slitimage

    def set_slitimage(self, name):
        self.slitimage = name

    def is_longslit(self):
        return self.longslit

    def set_longslit(self, lslit):
        self.longslit = lslit
    
    def should_try_again(self):
        return self.tryagain

    def set_try_again(self, tryagain):
        self.tryagain = tryagain

    def should_measure_slit(self):
        return self.measureslit

    def set_measure_slit(self, measureslit):
        self.measureslit = measureslit

class ObjectCoords(object):
    def __init__(self, acqimage):
        self.acqimage = acqimage
    
    def set_object_coords(self, coords):
        self.object_coords = coords

        xobj, yobj = coords
        info("Object at %.2f %.2f" % (xobj, yobj))

    def set_first_object_coords(self, coords):
        self.first_object_coords = coords

        xobj, yobj = coords
        info("Object #1 at %.2f %.2f" % (xobj, yobj))
        
    def set_second_object_coords(self, coords):
        self.second_object_coords = coords

        xobj, yobj = coords
        info("Object #2 at %.2f %.2f" % (xobj, yobj))

    def _check_ordering(self):
        # Make sure that the first object is on the left for horizontal slits and on the bottom for vertical slits:
        xobj1, yobj1 = self.first_object_coords
        xobj2, yobj2 = self.second_object_coords
        
        if (((self.acqimage.is_gnirs() or self.acqimage.is_f2()) and xobj1 > xobj2) or
            (not (self.acqimage.is_gnirs() or self.acqimage.is_f2()) and yobj1 > yobj2)):
            debug("...switching input coordinates...")
            self.first_object_coords, self.second_object_coords = self.second_object_coords, self.first_object_coords

    def get_first_object_coords(self):
        self._check_ordering()
        return self.first_object_coords

    def get_second_object_coords(self):
        self._check_ordering()
        return self.second_object_coords
        
    def get_object_coords(self):
        if hasattr(self, "object_coords"):
            return self.object_coords

        self.object_coords = (self.first_object_coords + self.second_object_coords) / 2.0
        debug("...mean obj position =", self.object_coords)

        return self.object_coords
    
def get_timestamps():
    now = datetime.datetime.now()
    datestring = now.strftime("%Y-%m-%d")
    timestring = now.strftime("%H:%M:%S")
    return datestring, timestring

ISO_TIMESTAMP_FORMAT = "%Y-%m-%dT%H:%M:%S"
def parse_iso_timestamp(timestamp):
    ts, partial_seconds = timestamp[:-1].split('.')
    partial_seconds = float("." + partial_seconds)
    time = datetime.datetime.strptime(ts, ISO_TIMESTAMP_FORMAT)
    precisedatetime = time + datetime.timedelta(seconds=partial_seconds)
    return precisedatetime

def total_minutes(timediff):
    return timediff.days * 24 * 60 + timediff.seconds / 60

def ask_question(prompt, default):
    msg = "%s (%s) " % (prompt, default)
    data = raw_input(msg)
    if data.strip() == "":
        data = default
    
    if data.lower()[0] == "y":
        return True
    return False    

def split_image_name(fname):
    rawpath = os.path.dirname(fname) 
    image_root, ext = os.path.splitext(os.path.basename(fname))
    return rawpath, image_root

def nearby_images(acqimage, backwards=False, skipfirst=False):
    rawpath, image_root = split_image_name(acqimage.filename)
    
    step = 1
    if backwards:  # look backwards if requested
        step = -1

    index = image_root.find("S", 1) + 1 # skip the first character, then skip the 'S' as well
    base = image_root[:index] # image root (e.g. S20081003S)

    start = int(image_root[index:]) # image number (e.g. 0065)
    if skipfirst:
        start += step

    how_many_to_check = 5
    stop = start + how_many_to_check * step
    stop = max(stop, 1)

    debug("...start = ", start, "   stop = ", stop, "   step = ", step)

    for fileno in range(start, stop, step):
        imagename = "%s%04d.fits" % (base, fileno)
        imagename = os.path.join(rawpath, imagename)

        yield imagename
 
def search_nearby_for_slit_image(mode, acqimage):
    debug("...looking for a slit image...")

    inst = acqimage.instrument()

    backwards = False
    if acqimage.is_gnirs():  # look backwards for slit image
        backwards = True

    for slitimage in nearby_images(acqimage, backwards):
        
        debug("...checking ", slitimage, "...")

        if os.path.exists(slitimage):
            try:
                slitimage_ad = AcquisitionImage(slitimage)
            except ValueError, e:
                continue
            
            slitinst = slitimage_ad.instrument()
            if slitinst != inst: # we've switched to another instrument
                raise RuntimeError

            if slitimage_ad.observation_id() != acqimage.observation_id():
                raise RuntimeError

            slitmask = slitimage_ad.focal_plane_mask()
            debug("...mask = ", slitmask)

            if slitimage_ad.has_mask_in_beam():
                if slitimage_ad.has_mos_mask():
                    warning("")
                    warning("! WARNING: Found '%s' with a MOS mask in the beam !" % slitimage)

                    if ask_question("Continue running acquisition on '%s'?" % acqimage.filename, "no"):
                        raise RuntimeError

                    error("Try running GACQ on the image with the MOS mask in the beam: %s" % slitimage)
                    sys.exit(1)
                
                debug("...found an image with mask = ", slitmask)
                debug("...turning on longslit mode...")

                mode.set_slitimage(slitimage)
                mode.set_longslit(True)

                return mode
    raise RuntimeError

def find_slit_image(mode, acqimage):
    if mode.get_slitimage() == "auto":
        try:
            slit = search_nearby_for_slit_image(mode, acqimage)
        except RuntimeError:
            warning("! WARNING: No slit image found !")
            mode.set_measure_slit(False)
            mode.set_try_again(True)     # look for the slit image later
            mode.set_slitimage("auto")   # reset and try again later
            debug("...measureslit = ", mode.should_measure_slit())
            debug("...l_slitimage = ", mode.get_slitimage())

    else:    # Slit image explicitly given as parameter
        debug("...expanding supplied slit image...")

        slitimagename = mode.get_slitimage()
        if is_number(slitimagename):
            imstring = "%04d" % int(slitimagename)
            slitimagename = l_prefix + imstring + ".fits"
            debug("...slit image = ", slitimagename)

        rawpath, image_root = split_image_name(acqimage.filename)
        if not os.path.exists(slitimagename):
            slitimagename = os.path.join(rawpath, slitimagename)

        if not os.path.exists(slitimagename):
            error("ERROR: Cannot access slit image " + slitimagename)
            mode.set_measure_slit(False)
            mode.set_try_again(True) # look for the slit image later
        else:
            mode.set_longslit(True)
            mode.set_slitimage(slitimagename)
            
    return mode

def extract_slit_info(mode):
    slitimage = mode.get_slitimage()
    info("Using slit image", slitimage)
    slit_root = os.path.basename(slitimage)

    debug("...slit_root = ", slit_root)

    slitimage_ad = AcquisitionImage(slitimage)

    slitmask = str(slitimage_ad.focal_plane_mask())
    info("Slit image has mask = " + slitmask)

    nyslitimg, nxslitimg = slitimage_ad.get_science_data().shape
    debug("...slit image dimensions = ", nxslitimg, " x ", nyslitimg)

    # Human-error checking
    if slitimage_ad.is_gmos():
        slitgrating = slitimage_ad.phu.header["GRATING"]
        debug("...slit grating = ",slitgrating)
        if "arcsec" not in slitmask:
            info("This is not a slit image!")
            mode.set_measure_slit(False)
        elif slitgrating != "MIRROR":
            error("ERROR: This is a spectrum! (Grating = " + slitgrating + ")")
            mode.set_measure_slit(False)

    elif slitimage_ad.is_niri():
        slitfilter = slitimage_ad.phu.header["FILTER3"]
        debug("...slit filter = ",slitfilter)
        if "cam" in slitmask:
            info("This is not a slit image!")
            mode.set_measure_slit(False)
        elif "grism" in slitfilter:
            error("ERROR: This is a spectrum! (Filter3 = " + slitfilter + ")")
            mode.set_measure_slit(False)

    return (mode, slitimage_ad, nxslitimg, nyslitimg)

def image_identifier(ad):
    field_names = ["OBSID",
                   "INSTRUME",
                   ad.focal_plane_mask(),  
                   "CAMERA",  
                   "FILTER",  
                   "FILTER1", 
                   "FILTER2", 
                   "FILTER3", 
                   "FLIP",    
                   "GRATING", 
                   "GRATWAV", 
                   "DECKER",  
                   "COADDS",  
                   "EXPTIME"
                   ]
    
    fields = [str(ad.get_science_data().shape[0])]
    for name in field_names:
        val = None
        if name in ad.phu.header:
            val = ad.phu.header[name]
        fields.append(val)

    return tuple(fields)

def get_axis_distance(dim, obj, bin):
    return (obj - dim / 2.0) * bin

def get_axis_distance_arcsec(dim, obj, bin, pixscale):
    return abs(get_axis_distance(dim, obj, bin)) * pixscale

def check_off_axis(args, dist):
    info("")
    if args.twotarget:
        info("The midpoint between the two targets is %.1f\" from the middle of the slit." % dist)
        msg = "Should the pair be centered along the length of the slit?"
    else:
        info("The target is %.1f\" from the middle of the slit." % dist)
        msg = "Should the target be centered along the length of the slit?"

    return ask_question(msg, "no")

def calculate_slit_display_scale(slitimage_ad, slitxpos, slitypos):
    scidata = slitimage_ad.get_science_data()
    
    slitxsize, slitysize = slitimage_ad.get_slit_size_in_pixels()
    xbuffer = slitxsize / 2.0
    xmin, xmax = clamp(scidata, slitxpos, xbuffer, 1)

    ybuffer = slitysize / 2.0
    ymin, ymax = clamp(scidata, slitypos, ybuffer, 0)
    
    debug("...calculating statistics from on the slit data [%f:%f,%f:%f] " % (ymin, ymax, xmin, xmax))
    
    slitdata = scidata[ymin:ymax, xmin:xmax]
    
    mode = scipy.stats.mstats.mode(slitdata, axis=None)[0]
    stddev = slitdata.std()
    
    debug("...in slit:  mode=", mode,"  stddev=", stddev)
    
    z1 = mode - 1.0 * stddev
    z2 = mode + 6.0 * stddev

    return z1, z2

def interactive_slit_measurement(slitimage_ad, slitxpos, slitypos, args):
    display_plot = args.verbose
    measurement = find_slit(slitimage_ad, slitxpos, slitypos, display_plot)
    draw_slit(measurement.get_midpoint(), slitimage_ad.get_mask_width_in_pixels())
    
    print("   Press 'q' if this slit center is ok, if not, do one of the following:")
    print("")
    print("   - point to the approximate center of slit and press 'a'")
    print("     (try to avoid any bright objects in the slit)")
    print("   - point to the edge of the slit and press 'x', then point to the other edge and press 'x'")
    print("")
    print("   When done press 'q'.")
    print("")

    listener = SlitMeasurementCursorListener(slitimage_ad, args)
    listener.start()

    user_measurement = listener.get_slit_measurement()

    if user_measurement is not None:
        return user_measurement

    return measurement

def main(argv=[__name__]):
    parser = argparse.ArgumentParser(description="Calculate offsets based upon an acquisition image.")
    parser.add_argument("image")

    parser.add_argument("--verbose", "-v",
                        action="store_true",
                        help="Print lots of debuging information to the log and the console")

    modegroup = parser.add_mutually_exclusive_group()
    modegroup.add_argument("--imaging",
                           action="store_true",
                           help="Force into imaging mode")
    modegroup.add_argument("--mos", "-m",
                           action="store_true",
                           help="Force into multi-object spectroscopy mode")
    modegroup.add_argument("--ifu", "-i",
                           action="store_true",
                           help="Force into IFU mode")

    observatorygroup = parser.add_mutually_exclusive_group()
    observatorygroup.add_argument("--gemininorth", "-N",
                                  action="store_true",
                                  help="Pretend to run at Gemini-North")
    
    observatorygroup.add_argument("--geminisouth", "-S",
                                  action="store_true",
                                  help="Pretend to run at Gemini-South")


    targetgroup = parser.add_mutually_exclusive_group()
    targetgroup.add_argument("--onetarget", "-1",
                             action="store_true",
                             default=True,
                             help="Find the P and Q offset one target")

    targetgroup.add_argument("--twotarget", "-2",
                             action="store_true",
                             help="Find the P and Q offset, as well as the rotation for two targets")

    parser.add_argument("--blindoffset", "-b",
                        action="store_true",
                        help="Whether this is a blind offset acquisition")

    parser.add_argument("--log", "-l",
                        help="File to write logging information to, defaults to gacq.log.YEAR-MO-DAY in UTC")

    parser.add_argument("--noadvice", 
                        action="store_true",
                        help="Whether to print advice")

    parser.add_argument("--onaxis", 
                        action="store_true",
                        help="Force acquisitions to be on-axis, sending all Q offsets")

    parser.add_argument("--z1", 
                        type=int,
                        help="Lower bound on the displayed image")

    parser.add_argument("--z2", 
                        type=int,
                        help="Upper bound on the displayed image")

    parser.add_argument("--mdfdir",
                        help="Directory for where to look for MDF files")

    mdfgroup = parser.add_mutually_exclusive_group()
    mdfgroup.add_argument("--mosmasknum",
                          help="MOS mask number (0=default)")
    mdfgroup.add_argument("--mdffile",
                          help="MOS mask MDF file name")
    
    args = parser.parse_args(argv[1:])

    if args.twotarget:
        args.onetarget = False
    
    mode = GACQMode(args)

    #-----------------------------------------------------------------------

    l_image        = args.image # Image to display, can use number if current UT

    ifuimaging = IFUImaging()

    #-----------------------------------------------------------------------

    logfile, console = setup_logging(args.log, args.verbose)

    info("--------------------------------------------------------")
    info(" ".join(argv))
    info("--------------------------------------------------------")

    debug("...gacq version = ", __version__)
        
    if mode.get_slitimage() == "":
        mode.set_measure_slit(False)
    else:
        mode.set_measure_slit(True)

    debug("...measureslit = ", mode.should_measure_slit())

    #-----------------------------------------------------------------------

    debug("...display = ", mode.get_display_image())
    debug("...mdfdir = ", args.mdfdir)
    debug("...slitimage = ", mode.get_slitimage())

    #-----------------------------------------------------------------------
    # Construct image name

    debug("...l_image = ", l_image)
    if is_number(l_image):
        info("Constructing image name based on today's UT date...")
        prefix = mode.get_site_prefix()
        debug("...prefix = ", prefix)

        imstring = "%04d" % int(l_image)
        l_image = prefix + imstring + ".fits"
        debug("...l_image = ", l_image)


    # make sure it ends with .fits, different from the IRAF version
    if not l_image.endswith(".fits"):
        l_image += ".fits"

    l_image = os.path.abspath(l_image)

    if not os.access(l_image, os.R_OK):
        error("ERROR: Cannot access image ", l_image)
        sys.exit(1)

    mosmask = args.mosmasknum or args.mdffile
    ad = AcquisitionImage(l_image, mosmask, args.mdfdir)

    debug("...l_image = ", l_image)

    info("Using image", l_image)
    full_image_name = l_image  # keep this since l_image may be modified
    debug("...full_image_name = ", full_image_name)

    #-----------------------------------------------------------------------
    # Check file age

    if is_number(l_image):  # input was image number
        debug("...checking age of file...")

        file_date = ad.phu.header["DATE"]
        file_time = ad.phu.header["UT"]
        timestamp = file_date + "T" + file_time

        debug("...file = ", timestamp, "UT")

        file_datetime = parse_iso_timestamp(timestamp)

        now = datetime.datetime.utcnow()
        debug("...now =  ", now.strftime(ISO_TIMESTAMP_FORMAT), "UT")
            
        age = now - file_datetime
        
        debug("...age = ", str(age))

        minutes = total_minutes(age)
            
        if minutes > 5:
            warning("! WARNING: This file is over", minutes, " minutes old !")

            if not ask_question("Is this ok?", "no"):
                info("Sorry about that, better luck next time.")
                sys.exit(1)
    else:
        debug("... did not get a number as the image name so not checking age")



    #-----------------------------------------------------------------------

    info("Observer(s):", ad.phu.header["OBSERVER"])

    objname = ad.phu.header["OBJECT"]
    info("Object(s):", objname)

    if not args.blindoffset:
        debug("...checking object name...")

        if "blind" in objname.lower():
            args.blindoffset = ask_question("Is this a blind offset acquisition?", "no")

    #-----------------------------------------------------------------------
    # Determine instrument

    objectcoords = ObjectCoords(ad)

    l_inst = ad.instrument()
    info("Instrument:", l_inst)    

    # reset the format on the log to include the instrument
    formatter = logging.Formatter(fmt="%(asctime)s %(name)s - " + l_inst + ": %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    logfile.setFormatter(formatter)

    formatter = logging.Formatter(fmt="%(name)s - " + l_inst + ": %(message)s")
    console.setFormatter(formatter)

    if ad.is_gmos():
        detector = ad.phu.header["DETECTOR"]
        debug("...DETECTOR = ", detector)
        if detector == "GMOS + e2v DD CCD42-90":
            namps = ad.phu.header["NAMPS"]
            debug("...NAMPS = ", namps)

            if namps == 1:
                error("ERROR: The GMOS e2v detectors should not use 3-amp mode!")
                error("Please set the GMOS CCD Readout Characteristics to Use 6 Amplifiers.")
                sys.exit(1)

    #-----------------------------------------------------------------------

    info("GEMPRGID = ", ad.program_id())

    obsid = ad.observation_id()
    debug("...obsid = ", obsid)

    siteprefix = obsid[:2]
    debug("...site prefix = ", siteprefix)

    #-----------------------------------------------------------------------
    # Determine configuration & set parameters

    ny, nx = ad.get_science_data().shape
    info("Image dimensions = " + str(nx) + " x " + str(ny))

    l_iaa = ad.phu.header["IAA"]
    info("IAA = ", l_iaa)

    m2baffle = ad.phu.header["M2BAFFLE"]
    debug("...M2BAFFLE = ", m2baffle)

    m2cenbaf = ad.phu.header["M2CENBAF"]
    debug("...M2CENBAF = ", m2cenbaf)

    # Assume that the mask will have the same dimensions as the image:
    nxslitimg = nx
    nyslitimg = ny

    info("Mask = ", ad.focal_plane_mask())

    # Do not draw any overlays if a mask is in the beam:
    if ad.has_mask_in_beam() or ad.is_mos_mode():
        info("Mask is in, turning off slit/mask overlay")
        mode.set_show_slit(False)

    info("Plate scale = " + str(ad.unbinned_pixel_scale()) + " arcsec / unbinned pixel")

    #-----------------------------------------------------------------------
    # Human error handling:

    if ad.is_gmos():
        grating = ad.phu.header["GRATING"]

        debug("...GRATING = ", grating)
            
        if grating != "MIRROR":
            info("Grating = " + grating)
            error("ERROR: This is a spectrum!")
            sys.exit(1)

        l_obstype = ad.phu.header["OBSTYPE"]

        debug("...OBSTYPE = ", l_obstype)
            
        if l_obstype == "BIAS":
            info("OBSTYPE = " + l_obstype)
            error("ERROR: This is a bias!")
            sys.exit(1)

        if m2baffle != "VISIBLE":
            info("M2 BAFFLE = ", m2baffle)
            warning("! WARNING: M2 baffle is not in the visible position !")
            if not ask_question("Is this okay?", "no"):            
                print ("Sorry about that, better luck next time.")
                sys.exit(1)

        if m2cenbaf != "CLOSED":
            info("M2 CENTRAL BAFFLE = ", m2cenbaf)
            warning("! WARNING: M2 central baffle is not closed !")
            if not ask_question("Is this okay?", "no"):
                print ("Sorry about that, better luck next time.")
                sys.exit(1)


    elif ad.is_niri():

        filter = ad.phu.header["FILTER3"]
        debug("...FILTER3 = ", filter)
            
        if "grism" in filter:
            info("Filter3 = " + filter)
            error("ERROR: This is a spectrum!")
            sys.exit(1)

        l_window = ad.phu.header["WINDCOVR"]
        debug("...WINDCOVR = ", l_window)
        if "open" not in l_window.lower():
            info("Window cover = " + l_window)
            error("ERROR: Window cover is %s!" % l_window)
            sys.exit(1)

        if args.twotarget and ny != 1024:
            error("ERROR: two-target mode is not supported with NIRI subarrays.")
            error("Please either set the NIRI array size to full frame or disable two-target mode.")
            sys.exit(1)

    elif ad.is_nifs():

        l_window = ad.phu.header["WINDCOVR"]
        debug("...WINDCOVR = ", l_window)
        if "open" not in l_window.lower():
            info("Window cover = " + l_window)
            error ("ERROR: Window cover is %s!" % l_window)
            sys.exit(1)

    elif ad.is_gnirs():

        if ad.focal_plane_mask() == "pupil viewer":
            error ("ERROR: This is a pupil image!")
            sys.exit(1)

        acqmir = ad.phu.header["ACQMIR"]
        debug("...ACQMIR = ", acqmir)
        if acqmir != "In":
            info("Acquisition Mirror = " + acqmir)
            error ("ERROR: This is a spectrum!")
            sys.exit(1)

    elif ad.is_f2():
        grism = ad.phu.header["GRISM"]
        debug("...GRISM = ", grism)
        if "open" not in grism.lower():
            info("Grism cover = " + grism)
            error ("ERROR: This is a spectrum!")
            sys.exit(1)

        l_window = ad.phu.header["WINDOW"]
        debug("...WINDOW = ", l_window)
        if "open" not in l_window.lower():
            info("Window cover = " + l_window)
            error ("ERROR: Window cover is %s!" % l_window)
            sys.exit(1)


    ratrgoff = ad.phu.header["RATRGOFF"] # Target offset in RA in arcsec
    dectrgof = ad.phu.header["DECTRGOF"] # arget offset in DEC in arcsec
    if abs(ratrgoff) > 10 or abs(dectrgof) > 10:
        warning("! WARNING: Telescope is at an offset position !")
        warning("      RA Target Offset  = %.3f arcsec" % ratrgoff)
        warning("      Dec Target Offset = %.3f arcsec" % dectrgof)
        if not ask_question("Continue anyway?", "yes"):
            info ("Sorry about that, better luck next time.")
            sys.exit(1)
            
    obslog = obsid + ".log"
    debug("...obslog = ", obslog)
    if not os.access(obslog, os.R_OK):
        hdr_poff = ad.phu.header["POFFSET"]
        hdr_qoff = ad.phu.header["QOFFSET"]
        debug("...poffset = ", hdr_poff, "   qoffset = ", hdr_qoff)
        if hdr_poff != 0.0 or hdr_qoff != 0.0:
            warning("! WARNING: Telescope is at an offset position !")
            warning("      P = %.3f arcsec" % hdr_poff)
            warning("      Q = %.3f arcsec" % hdr_qoff)
            if not ask_question("Continue anyway?", "yes"):
                info ("Sorry about that, better luck next time.")
                sys.exit(1)


    #-----------------------------------------------------------------------
    # Check if the IFU is in the beam

    if ad.is_gmos() and "IFU" in ad.focal_plane_mask():
        info("IFU in beam!")

        xbin = ad.detector_x_bin()
        ybin = ad.detector_y_bin()
        if xbin != 1 or ybin != 1:
            error("ERROR: acquisitions with the GMOS IFU must be binned 1x1")
            error("The binning of this image is ", xbin, "x", ybin)
            error("Sorry about that, better luck next time.")
            sys.exit(1)

        info("Turning on IFU mode")
        args.ifu = True
        args.onetarget = False
        args.twotarget = False
        mode.set_measure_slit(False)
        mode.set_display_image(False)
        mode.set_show_slit(False)

    #-----------------------------------------------------------------------
    # Automatic checks for whether this is MOS acquisition
    if ad.is_mos_mode():
        if ad.has_mos_mask():
            info("MOS mask in beam!")
        
        debug("...turning off slit measurement")
        args.mos = True
        mode.set_measure_slit(False)

    #-----------------------------------------------------------------------
    # Look for slit image (first time)

    if mode.should_measure_slit():
        mode = find_slit_image(mode, ad)

    if mode.should_measure_slit():
        (mode, l_slitimage_ad, nxslitimg, nyslitimg) = extract_slit_info(mode)
     
    #-----------------------------------------------------------------------
    # Display GMOS menu when we can't figure out what mode this is

    if ad.is_gmos() and not mode.is_longslit() and not args.mos and not args.ifu:
        info ("")
        info ("Configurations:")
    #       info ("   0 = Imaging")        # hidden option
        info ("   1 = IFU-R")
        info ("   2 = IFU-2")
        info ("   3 = MOS")
        info ("   4 = Long-slit")
    #       info ("   5 = IFU-B")          # hidden option
        info ("   9 = Exit")
        config = raw_input("Which configuration number matches the science observation? ")
        l_config = int(config)  # get configuration from user
        info("")
        if l_config == 0:
            info("Imaging")
            args.imaging = True
            mode.set_show_slit(False)
            
        elif l_config == 1:
            info("IFU Red-Slit")
            ifuimaging = IFUImaging("IFU-R", ad)
            mode.set_show_slit(False)
            mode.set_try_again(False)

        elif l_config == 2:
            info("IFU 2-Slit")
            ifuimaging = IFUImaging("IFU-2", ad)
            mode.set_show_slit(False)
            mode.set_try_again(False)

        elif l_config == 3:
            info("Multi-Object Spectroscopy")
            args.mos = True
            args.onetarget = False
            args.twotarget = False
            mode.set_measure_slit(False)
            mode.set_display_image(False)
            mode.set_show_slit(False)

            info ("")
            config = raw_input("MOS mask number? ")
            ad = AcquisitionImage(l_image, int(config), args.mdfdir)

        elif l_config == 4:
            info("Long-slit spectroscopy")
            mode.set_longslit(True)
            if mode.get_slitimage() != "":
                mode.set_measure_slit(False)
                debug("...measureslit = ", mode.should_measure_slit())

        elif l_config == 5:
            info("IFU Blue-Slit")
            ifuimaging = IFUImaging("IFU-B", ad)
            mode.set_show_slit(False)
            mode.set_try_again(False)

        elif l_config == 9:
            info ("Exiting...")
            sys.exit(1)
        else: 
            error("ERROR: I did not understand your response.")
            error("Sorry about that, better luck next time.")
            sys.exit(1)

    #-----------------------------------------------------------------------
    # Check PA if this is a spectrophotometric standard

    if ad.is_gmos() and mode.is_longslit():
        debug("...checking if this is a spectrophotometric standard...")
        # List of GS + GN specphot standards as of 2012 March
        spectrophotometric_standards = [
            "BD+28 4211",
            "CD-32 9927",
            "CD-34 241", 
            "EG131",     
            "EG21",      
            "EG274",     
            "Feige110",  
            "Feige34",   
            "Feige56",   
            "Feige66",   
            "G158-100",  
            "G191B2B",   
            "GD108",     
            "H600",      
            "Hiltner600",
            "HZ44",      
            "LTT1020",   
            "LTT1788",   
            "LTT2415",   
            "LTT2915",   
            "LTT3218",   
            "LTT 377",   
            "LTT3864",   
            "LTT4364",   
            "LTT4816",   
            "LTT6248",   
            "LTT7379",   
            "LTT7987",   
            "LTT9239",   
            "LTT9491",   
            "VMa2",      
            "Wolf1346",
            ]

        for name in spectrophotometric_standards:
            if name in objname:
                warning ("! WARNING:  If this is a spec-phot standard be sure to set the PA to the parallactic angle !")
                # To-Do: Calculate the Parallactic angle and compare with the IPA; Warn if difference is > ~5 deg


    if ad.is_mos_mode():
        info("Using MDF", ad.get_mdf_filename())
        info("Turning on MOS mode")
        args.mos = True
        args.onetarget = False
        args.twotarget = False
        mode.set_measure_slit(False)
        mode.set_display_image(False)
        mode.set_show_slit(False)

    #-----------------------------------------------------------------------
    # Set the default slit positions

    if ( not args.mos and
         not args.ifu and
         not args.imaging ):
        debug("...finding default slit position...")

        if ad.is_gmos():
            if ad.is_gmosn(): # GMOS-N values updated 2006 Nov 25
                # Filter  dslitpos (pix)  (r-filter is the reference)
                dx_g  =  -0.4
                dx_i  =  -4.0
                dx_z  =   0.8
                # dx_g  =  -0.49  # 2005 Nov
                # dx_i  =  -4.29  # 2005 Nov
                # dx_z  =   0.75  # 2005 Nov

            elif ad.is_gmoss(): # GMOS-S values updated 2006 Nov 22
                # Filter  dslitpos (pix)  (r-filter is the reference)
                dx_g  =  -0.24
                dx_i  =  -2.08
                dx_z  =   0.73

            filt = ad.phu.header["FILTER2"]
            info("Adjusting default position for the " + filt + " filter...")

            firstchar = filt[0]
            if   firstchar == "r":
                dx_filt = 0.0
            elif firstchar == "g":
                dx_filt = dx_g
            elif firstchar == "i":
                dx_filt = dx_i
            elif firstchar == "z":
                dx_filt = dx_z
            else:
                warning("Unknown offset for this filter; assuming 0")
                dx_filt = 0.0

            goal_center = ad.get_goal_center()
            debug("...default goal center =", goal_center)
            goal_center += np.array((dx_filt, 0.0))
            
            if ad.is_gmoss():
                goal_center += np.array([-2.0, 0.0])

            debug("...new goal center =", goal_center)
            ad.set_goal_center(goal_center)

        elif ad.is_niri():
            info("Field center selected")
            mode.set_show_slit(False)
            
            args.imaging = True

        else:
            error("Instrument not supported")
            sys.exit(1)

    #-----------------------------------------------------------------------
    # Display image

    if mode.get_display_image():
        slitxpos, slitypos = ad.get_binned_goal_center()
        debug("...slitxpos = ", slitxpos)
        debug("...slitypos = ", slitypos)

        # If this is a slit image, set the display levels based on the pixels in slit:
        z1 = args.z1
        z2 = args.z2
        zscale = False
        if z1 is None or z2 is None:
            if ad.has_mask_in_beam():
                z1, z2 = calculate_slit_display_scale(l_slitimage_ad, slitxpos, slitypos)
                zscale = False
            else:
                zscale = True

        scidata = ad.get_science_data()
        ui.display(scidata,
                   z1=z1,
                   z2=z2,
                   zscale=zscale)

        # mark saturated pixels in red
        satval = ad.saturation_level()

        indices = list((scidata >= satval).nonzero())
        indices.reverse()
        points = zip(*indices)

        info("Pixels with > %.1f ADU marked in red" % satval)
        ui.mark_points(points)
        
    #-----------------------------------------------------------------------
    # Draw the slit

    if mode.get_show_slit():
        debug("...drawing slit...")

        slitxpos, slitypos = ad.get_binned_goal_center()
        ydim, xdim = ad.get_science_data().shape

        if ad.is_gmos():
            l_x1 = slitxpos
            l_x2 = slitxpos
            l_y1 = 1
            l_y2 = ydim
        else:
            error("Only GMOS supported for showing the slit")
            sys.exit(1)

        debug("...l_x1 = ", l_x1)
        debug("...l_x2 = ", l_x2)
        debug("...l_y1 = ", l_y1)
        debug("...l_y2 = ", l_y2)

        info("Drawing slit at %8.2f, %8.2f" % (slitxpos, slitypos))
        
        points = [(l_x1, l_y1), (l_x2, l_y2)]
        ui.polyline(points=points, color=ui.MAGENTA)

    #-----------------------------------------------------------------------
    # Calculate IFU positions and draw the IFU overlays if requested

    if ifuimaging.is_ifu_imaging():
        debug("...calculating IFU positions...")

        ifucenter = ifuimaging.get_ifu_center()
        debug("...ifu center =", ifucenter)
        ad.set_binned_goal_center(ifucenter)

        ifuborders = ifuimaging.get_ifu_borders()
        ui.polyline(points=ifuborders, color=ui.MAGENTA)

        ui.marker(ifuimaging.get_ifu_text_location(), color=ui.BLUE, mark=ifuimaging.get_mode())

        skycenter = ifuimaging.get_sky_center()
        debug("...sky center =", skycenter)
        ui.marker(skycenter, color=ui.MAGENTA)

        skyborders = ifuimaging.get_sky_borders()
        ui.polyline(points=skyborders, color=ui.MAGENTA)
        ui.marker(ifuimaging.get_sky_text_location(), color=ui.BLUE, mark="SKY")



    #-----------------------------------------------------------------------
    # Single target

    if args.onetarget:
        print ("")
        print ("   Put the cursor on the target and press 'a' or 'r', the target will then be fitted and marked in DS9.")
        print ("   - press 'x' to use an exact location on the acquisition image")
        print ("   - press 'e' to see a plot of the profile (clicking in the plot selects an exact location)")
        print ("   - press 'g', then click on the contour plot to use an exact location (use 'a' or 'e' to bring up the contour plot first)")
        print ("")
        if ifuimaging.is_ifu_imaging():
            print ("   If you would like the target to be centered in the IFU press 'q'.")
        elif mode.is_longslit():
            print ("   If you would like the target to be centered in the slit press 'q'.")
        elif args.imaging:
            print ("   If you would like the target to be centered in the image press 'q'.")
        else:  # this should never happen
            assert False
            print ("   If you would like the target to be centered press 'q'.")
        
        print ("   Otherwise point to where you would like to put the target and press 'c' then 'q'.")
        print ("")

        center = ad.get_binned_goal_center()
        debug("...default center = ", center)

        listener = SelectionCursorListener(ad, args.verbose, mark_center=center)
        listener.start()

        objectcoords.set_object_coords(listener.get_object_coords())

        if listener.has_custom_center(): # were we given custom slit coordinates?
            mode.set_measure_slit(False)
            ad.set_binned_custom_center(listener.get_custom_center())


    #-----------------------------------------------------------------------
    # Two targets

    if args.twotarget:
        print ("")
        print ("   TWO-TARGET MODE: 1st target")
        print ("")
        print ("   Point to the object and press 'a' or 'r', the target will then be fitted and marked in DS9.")
        print ("   - press 'x' to use an exact location on the acquisition image")
        print ("   - press 'e' to see a plot of the profile (clicking in the plot selects an exact location)")
        print ("   - press 'g', then click on the contour plot to use an exact location (use 'a' or 'e' to bring up the contour plot first)")
        print ("")
        print ("   Press 'q' when you're happy with the selection")

        center = ad.get_binned_goal_center()
        listener = SelectionCursorListener(ad, args.verbose, mark_center=center)
        listener.start()
        
        if listener.has_custom_center():
            warning("Ignoring the custom center you just specified!")
            warning("It must be specified with the 2nd target")
        listener.undo_marked_center()

        objectcoords.set_first_object_coords(listener.get_object_coords())

        print ("")
        print ("   TWO-TARGET MODE: 2nd target")
        print ("")
        print ("   Point to the object and press 'a' or 'r', the target will then be fitted and marked in DS9.")
        print ("   - press 'x' to use an exact location on the acquisition image")
        print ("   - press 'e' to see a plot of the profile (clicking in the plot selects an exact location)")
        print ("   - press 'g', then click on the contour plot to use an exact location (use 'a' or 'e' to bring up the contour plot first)")
        print ("")
        print ("   If you would like the pair of objects to be centered in the slit press 'q'.")
        print ("   Otherwise point to where you would like to center the pair and press 'c', then 'q'.")
        print ("")

        listener = SelectionCursorListener(ad, args.verbose, mark_center=center)
        listener.start()

        objectcoords.set_second_object_coords(listener.get_object_coords())

        if listener.has_custom_center(): # were we given custom slit coordinates?
            mode.set_measure_slit(False)
            ad.set_binned_custom_center(listener.get_custom_center())

    #-----------------------------------------------------------------------
    # Look for slit image again if it was not found before

    if mode.should_try_again() and mode.is_longslit():
        mode = find_slit_image(mode, ad)

    if mode.should_try_again() and mode.is_longslit() and mode.should_measure_slit():
        (mode, l_slitimage_ad, nxslitimg, nyslitimg) = extract_slit_info(l_slitimage)

    #-----------------------------------------------------------------------
    # Check for off-axis GMOS acquisitions:

    if ad.is_gmos() and mode.is_longslit() and not ad.has_custom_center():
        xpos, ypos = objectcoords.get_object_coords()
        dist = get_axis_distance_arcsec(ny, ypos, ad.detector_y_bin(), ad.unbinned_pixel_scale())
        if not args.onaxis and dist > 2.0:
            if check_off_axis(args, dist):
                args.onaxis = True
            else:
                ydist = get_axis_distance(ny, ypos, ad.detector_y_bin())
                debug("...ydist = ", ydist)
                
                center = ad.get_goal_center()
                center += np.array((0.0, ydist))
                ad.set_goal_center(center)

    #-----------------------------------------------------------------------
    # Measure slit position on slit image (if we have one)

    if mode.should_measure_slit():
        info("Measuring slit position...")

        diff = ad.get_data_center() - ad.get_goal_center()
        slitcenter = l_slitimage_ad.get_data_center() - diff
        slitxpos, slitypos = slitcenter / l_slitimage_ad.detector_y_bin()

        debug("...slitxpos = ", slitxpos)
        debug("...slitypos = ", slitypos)

        info("Calculating display scale for slit image...")
        z1, z2 = calculate_slit_display_scale(l_slitimage_ad, slitxpos, slitypos)

        debug("...z1 = ", z1, "   z2 = ", z2)
        ui.display(l_slitimage_ad.get_science_data(),
                   z1=z1,
                   z2=z2,
                   zscale=False)
                           

    if mode.should_measure_slit(): # try to automagically find the center of the slit
        if args.twotarget: 
            # measure the slit position at the y-coordinate of each object
            xpos, ypos = objectcoords.get_first_object_coords()
            measurement1 = interactive_slit_measurement(l_slitimage_ad, slitxpos, ypos, args)
            xpos, ypos = objectcoords.get_second_object_coords()
            measurement2 = interactive_slit_measurement(l_slitimage_ad, slitxpos, ypos, args)

            mid1 = measurement1.get_midpoint()
            mid2 = measurement2.get_midpoint()

            # then calculate the slit tilt
            tilt = get_slit_tilt(mid1, mid2)
            info("Measured slit tilt: %.1f degrees (expected %.1f) " % (tilt, l_slitimage_ad.get_expected_slit_tilt()))
            slittilt = tilt

            slitline = [mid1, mid2]
            ui.polyline(slitline, color=ui.MAGENTA)

            # find slitxpos at the given slitypos
            diff = mid2 - mid1
            slope = diff[1] / diff[0]
            slitxpos = ((slitypos - mid1[1]) / slope) + mid1[0]
        else:
            measurement = interactive_slit_measurement(l_slitimage_ad, slitxpos, slitypos, args)
            midpoint = measurement.get_midpoint()
            slitxpos = midpoint[0]

        # tweak the offsets based upon the slit measurement
        diff = l_slitimage_ad.get_binned_data_center() - np.array([slitxpos, slitypos])
        diff = diff * l_slitimage_ad.detector_y_bin()

        center = ad.get_data_center()
        newcenter = center - diff
        ad.set_goal_center(newcenter)
        
    #-----------------------------------------------------------------------
    # Calculate one- and two-target offsets
    l_xslit, l_yslit = ad.get_binned_goal_center()
    if (args.onetarget or args.twotarget) and not ad.has_custom_center():
        if mode.is_longslit() and not args.onaxis:
            xobj, yobj = objectcoords.get_object_coords()

            if ad.is_gmos() and get_axis_distance_arcsec(ny, yobj, ad.detector_y_bin(), ad.unbinned_pixel_scale()) > 5.0:
                warning("! WARNING:  Off-axis acquisition !")
                l_yslit = yobj
            elif ad.is_f2() and get_axis_distance_arcsec(nx, xobj, ad.detector_x_bin(), ad.unbinned_pixel_scale()) > 5.0:
                warning("! WARNING:  Off-axis acquisition !")
                l_xslit = xobj

    debug("...l_xslit = ", l_xslit)
    debug("...l_yslit = ", l_yslit)

    l_poff = 0.0
    l_qoff = 0.0
    l_rot  = 0.0

    if args.onetarget:
        info("Mask =   %7.2f %7.2f" % (l_xslit, l_yslit))

        xobj, yobj = objectcoords.get_object_coords()
        if xobj > 0.0 and yobj > 0.0 and l_xslit > 0.0 and l_yslit > 0.0:
            l_xoff = (l_xslit - xobj) * ad.binned_pixel_scale()
            l_yoff = (l_yslit - yobj) * ad.binned_pixel_scale()
            l_rot  = 0.0
        else:
            error("ERROR: Problem reading coordinates")
            sys.exit(1)

    elif args.twotarget:
        # calculate the unit vector representing the slit
        slitvec = np.array([0, 1]) # vertical slits
        if ad.is_gnirs() or ad.is_f2():  # horizontal slits
            slitvec = np.array([1, 0])

        slittilt = ad.get_expected_slit_tilt()
        debug("...slit tilt = ", slittilt)

        slitvec = rotate(slitvec, slittilt)

        # calculate the rotation
        obj1 = objectcoords.get_first_object_coords()
        obj2 = objectcoords.get_second_object_coords()
        objvec = obj2 - obj1
        l_rot = angle(slitvec, objvec)

        # calculate the translation
        slitcenter = np.array([l_xslit, l_yslit]) 
        debug("...slitcenter =", slitcenter)

        fieldcenter = ad.get_field_center() / ad.detector_y_bin()
        debug("...fieldcenter =", fieldcenter)
        if ad.is_gnirs(): # this is the best we can do until the acq mirror is fixed
            fieldcenter = slitcenter

        obj = objectcoords.get_object_coords()
        
        l_xoff, l_yoff = obj - fieldcenter + rotate(fieldcenter - slitcenter, -l_rot)

        l_xoff = -l_xoff * ad.binned_pixel_scale()
        l_yoff = -l_yoff * ad.binned_pixel_scale()


    #-----------------------------------------------------------------------
    # Transform (X,Y) to (P,Q)

    if args.onetarget or args.twotarget:
        debug("...converting from (X,Y) to (P,Q)")
        debug("...X = ", l_xoff, "   Y = ", l_yoff, "   R = ", l_rot)

        if ad.is_niri():
            l_rot  = -l_rot
            l_poff =  l_xoff
            l_qoff = -l_yoff
            if ad.is_altair() and ad.is_south_port():
                l_poff = -l_poff

        else:
            l_poff = l_xoff
            l_qoff = l_yoff

        if ad.is_south_port():
            if ad.is_gmoss():
                l_qoff = -l_qoff
            else:
                l_poff = -l_poff
                l_rot  = -l_rot

        debug("...P = ", l_poff, "   Q = ", l_qoff, "   R = ", l_rot)

    #-----------------------------------------------------------------------
    # MOS mask alignment
    rexpect = 0.0        # expected rotation
    if args.mos:
        debug("...MOS mode...")

        if ad.is_gmosn():
            rexpect = 0.01

        if ad.get_num_mos_boxes() < 2:
            error("MOS mask alignment requires at least 2 acquisition boxes!")
            error("Sorry about that, better luck next time.")
            sys.exit(1)

        interactive_mosaic = InteractiveMosaic(ad, args.verbose)
        z1, z2 = interactive_mosaic.get_zscale()
        zscale = False
        if z1 is None:
            zscale = True
        ui.display(ad.get_science_data(), z1=z1, z2=z2, zscale=zscale)

        for border in ad.get_mos_box_borders():
            ui.polyline(points=border, color=ui.BLUE)

        print ("")
        print ("   Trying to align boxes with the stars marked in green, you can choose to do any of the following:")
        print ("")
        print ("   - press 'q' to accept what is currently shown")
        print ("   - press 'a' or 'r' to select a different star")
        print ("   - press 'e' to see a contour plot")
        print ("   - press 'x' to use an exact location for a star")
        print ("   - press 'g', then click on the contour plot to use an exact location (use 'a' or 'e' to bring up the contour plot first)")

        if interactive_mosaic.has_mos_mask():
            print ("   - press 'b' at two opposite corners to manually specify an acquisition box")
        print ("   - press 'd' to discard a box entirely")
        print ("")

        listener = MOSBoxCursorListener(interactive_mosaic)
        listener.start()

        box_centers  = list(listener.get_box_centers())
        if len(box_centers) < 2:
            error("")
            error("MOS acquisitions requires at least 2 boxes!")
            error("")
            sys.exit(1)
        
        star_centers = list(listener.get_star_centers())

        num_boxes = len(box_centers)
        assert num_boxes == len(star_centers)
        debug("...box centers =", box_centers)
        debug("...star centers =", star_centers)

        fieldcenter = ad.get_field_center()
        debug("...fieldcenter =", fieldcenter)
        box_centers = np.array(box_centers) - fieldcenter
        star_centers = np.array(star_centers) - fieldcenter

        debug("...box centers relative to the field center =", list(box_centers))
        debug("...star centers relative to the field center =", list(star_centers))

        superposition = superimpose(star_centers, box_centers)

        l_xrms, l_yrms = superposition.get_piecewise_rms()
        
        info("")
        info("RMS of fit:  X: %-6.2f  Y: %-6.2f pixels" % (l_xrms, l_yrms))
        l_xrms = l_xrms * ad.unbinned_pixel_scale()
        l_yrms = l_yrms * ad.unbinned_pixel_scale()
        info("RMS of fit:  X: %-6.3f  Y: %-6.3f arcsec" % (l_xrms, l_yrms))

        slitxsize, slitysize = ad.get_min_slitsize()
        if l_xrms > 0.25 * slitxsize or l_yrms > 0.25 * slitysize:
            warning ("")
            warning ("! Poor mask alignment !")
            warning ("")
            warning ("  RMS greater than 25%% of the slit height (%.2f arcsec) or width (%.2f arcsec)" % (slitysize, slitxsize))
            warning ("")

            if not ask_question("Would you like to continue anyway?", "no"):
                info ("Sorry about that, better luck next time.")
                sys.exit(1)

        xshift, yshift = superposition.get_translation()
        debug("...superposition translation = %f, %f" % (xshift, yshift))

        xshift = -xshift * ad.unbinned_pixel_scale()
        yshift = -yshift * ad.unbinned_pixel_scale()

        rot = -superposition.get_euler_angle()
        debug("...superposition rotation =", rot)

        # Signs are opposite of the uplooking port, IJ
        if ad.is_south_port():
            if ad.is_gmoss():
                yshift = -yshift
            else:
                xshift = -xshift
                rot = -rot

        # For GMOS-S the signs are backwards with GMOS on Port 3:
        if ad.is_gmoss() and not ad.is_south_port():
            rot = -rot

        if ad.is_gmos():
            l_poff = xshift
            l_qoff = yshift
            l_rot  = rot
        elif ad.is_f2():
            l_poff = -yshift
            l_qoff = -xshift
            l_rot  = -rot

        debug("...rexpect = ", rexpect)
        if abs(rot - rexpect) > 0.04 and len(box_centers) == 2:
            warning("")
            warning("! The rotation (%-6.3f deg) is larger than expected !" % rot)
            warning("")
            if not ask_question("Would you like to continue anyway?", "no"):
                info ("Sorry about that, better luck next time.")
                sys.exit(1)

    #-----------------------------------------------------------------------
    # GMOS IFU reconstruction

    if args.ifu:
        debug("...GMOS IFU mode...")

        l_xoff, l_yoff = get_integrated_field_unit_offsets(ad, args.verbose)
        if ad.is_south_port():
            l_xoff = -l_xoff

        l_poff = l_xoff
        l_qoff = l_yoff
        l_rot  = 0.0

    #-----------------------------------------------------------------------
    # Print output & send offsets to TCC

    info ("")

    l_iadd = ad.phu.header["IAA"]
    debug("...IAA = ", l_iaa)

    l_ipa = ad.phu.header["PA"]
    debug("...PA = ", l_ipa)

    if args.twotarget or args.mos:
        info ("OFFSETS (arcsec):  P = %-8.3f  Q = %-8.3f   ROTATION (deg) = %-8.3f" % (l_poff, l_qoff, l_rot))

        if args.twotarget:           # Two-target mode: change the IPA
            l_ipa = l_ipa + l_rot

        elif args.mos and ad.is_gmoss(): # GMOS-S MOS mask alignment: change IAA
            l_iaa = l_iaa + l_rot
            if l_iaa <   0.0:
                l_iaa = l_iaa + 360.
            if l_iaa > 360.0:
                l_iaa = l_iaa - 360.
                
            info ("NEW IAA = %-8.3f" % l_iaa)

        elif args.mos and (ad.is_gmosn() or ad.is_f2()): # GMOS-N MOS mask alignment: change IPA
            l_ipa = l_ipa + l_rot
            if l_ipa <   0.0:
                l_ipa = l_ipa + 360.0
            if l_ipa > 360.0:
                l_ipa = l_ipa - 360.0

            info ("NEW IPA = %-8.3f" % l_ipa)
    else:
        info("OFFSETS (arcsec):  P = %-8.3f  Q = %-8.3f" % (l_poff, l_qoff))

    info ("")

    expected = (ad.get_goal_center() - ad.get_data_center()) * ad.unbinned_pixel_scale()
    pexpect, qexpect = expected
    qexpect = 0.0 # never expect much motion in Q
    if ifuimaging.is_ifu_imaging() and ad.is_gmoss():
        # GMOS south starts with the object in the IFU box
        pexpect = 0.0

    debug("...pexpect = ", pexpect, "    qexpect = ", qexpect, "    rexpect = ", rexpect)

    if abs(l_poff - pexpect) > 5.0 or abs(l_qoff - qexpect) > 5.0:
        warning("! WARNING: These offsets are larger than expected !")
        if not ask_question("Would you like to continue anyway?", "no"):
            info ("Sorry about that, better luck next time.")
            sys.exit(1)
        print ("")

    if mode.is_longslit() and l_poff == 0.0 and l_qoff == 0.0:
        warning("! WARNING: Offsets are zero !")
        warning("You may have clicked the same spot twice.")
        if not ask_question("Would you like to continue anyway?", "no"):
            info("Sorry about that, better luck next time.")
            sys.exit(1)
        info("")

    debug("...P=", l_poff, "  Q=", l_qoff, "  IAA=", l_iaa, "  IPA=", l_ipa)



    #-----------------------------------------------------------------------
    send2tcc = True # prompt="Send offsets to TCC?
    if not args.noadvice:
        debug("...Altair = ", ad.is_altair())
        debug("...GNIRS = ", ad.is_gnirs())
        debug("...NIFS = ", ad.is_nifs())
        
        if args.mos:
            debug("...slitxsize = ", slitxsize, " and slitysize = ", slitysize, "arcsec")

        debug("...MOS mode = ", args.mos)
        debug("...IFU mode = ", args.ifu)
        debug("...One target mode = ", args.onetarget)
        debug("...Two target mode = ", args.twotarget)
        debug("...Imaging mode = ", args.imaging)
        debug("...Long slit mode = ", mode.is_longslit())
        debug("...Mask in beam = ", ad.has_mask_in_beam())

        if args.imaging:
            debug("...giving imaging advice")
            if args.blindoffset:
                if abs(l_poff) > 0.01 * nx * ad.unbinned_pixel_scale() or abs(l_qoff) > 0.01 * ny * ad.unbinned_pixel_scale():
                    info("ADVICE: Apply offsets and take another acquisition image.")
                    info("        ( |P| or |Q| > 1% of imaging field size )")
                    send2tcc = True
                else:
                    info("ADVICE: Ignore offsets and move back to the base position.")
                    info("        Start science when guiding has been restored.")
                    send2tcc = False

            elif abs(l_poff) > 0.1 * nx * ad.unbinned_pixel_scale() or abs(l_qoff) > 0.1 * ny * ad.unbinned_pixel_scale():
                info("ADVICE: Apply offsets and take another acquisition image.")
                info("        ( |P| or |Q| > 10% of imaging field size )")
                send2tcc = True
            elif abs(l_poff) > 0.01 * nx * ad.unbinned_pixel_scale() or abs(l_qoff) > 0.01 * ny * ad.unbinned_pixel_scale():
                info("ADVICE: Apply offsets and start science immediately.")
                info("        ( |P| and |Q| < 10% of imaging field size )")
                send2tcc = True
            else:
                info("ADVICE: Ignore offsets and start science immediately.")
                info("        ( |P| and |Q| < 1% of imaging field size )")
                send2tcc = False

        elif args.ifu or ad.is_nifs():
            debug("...giving IFU advice")
                
            if abs(l_poff) > 0.1 or abs(l_qoff) > 0.1:
                info ("ADVICE: Apply offsets and take another acquisition image.")

                msg = "        ("            
                if abs(l_poff) > 0.1:
                    msg += " |P|>0.1\""
                if abs(l_qoff) > 0.1:
                    msg += " |Q|>0.1\""
                msg += " )"
                send2tcc = True
                
            else:
                if args.blindoffset:
                    info ("ADVICE: Ignore offsets and move back to the base position.")
                    info ("        Start science when guiding has been restored.")
                else:
                    info ("ADVICE: Ignore offsets and start science immediately.")

                send2tcc = False
                info ("        ( |P|<0.1\" and |Q|<0.1\" )")

        elif not ad.has_mask_in_beam():
            debug("...giving no-mask-in-beam advice")

            if ifuimaging.is_ifu_imaging(): # IFU acquisition image
                info ("ADVICE: Apply offsets and take the through-IFU image.")
                info ("Note that ALL acquisitions MUST have a through-IFU image.")
                send2tcc = True

            elif ad.is_gnirs() and args.twotarget:
                debug("...separation between targets = ", (xobj2 - xobj1), "pix")
                    
                l_roff = (xobj2 - xobj1) * ad.unbinned_pixel_scale() * tan(l_rot / 180.0 * math.pi) # linear offset at the targets due to rotation
                
                debug("...l_roff = ", l_roff, "arcsec")
                debug("...l_poff + l_roff = ", abs(l_poff) + abs(l_roff))

                change = abs(l_poff) + abs(l_roff)
                halfpx = ad.unbinned_pixel_scale() / 2.0
                tenpercentwidth = 0.1 * l_slitimage_ad.get_mask_width()
                limit = max(tenpercentwidth, halfpx)

                if change > limit or abs(l_qoff) > 0.5:
                    info ("ADVICE: Apply offsets and take another keyhole acquisition image.")
                    send2tcc = True
                    msg = "        ("
                    if change > limit:
                        msg += " |P|+|R| = %.0f%% of slit width" % (change / l_slitimage_ad.get_mask_width() * 100.0)
                        
                    if abs(l_qoff) > 0.5:
                        msg += " |Q|>0.5\""
                    msg += " )"
                    info (msg)
                    
                else:
                    info ("ADVICE: Ignore offsets and take the through-slit image.")
                    info ("Note that ALL acquisitions MUST have a through-slit image.")
                    send2tcc = False

                    if tenpercentwidth < halfpx:
                        info ("        ( |P|+|R| < pixscale/2.0\" and |Q|<0.5\" )")
                    else:
                        info ("        ( |P|+|R| = %.0f%% of slit width and |Q|<0.5\" )", change / l_slitimage_ad.get_mask_width() * 100.0)

            elif ad.is_gnirs(): # one-target GNIRS
                halfpx = ad.unbinned_pixel_scale() / 2.0
                tenpercentwidth = 0.1 * l_slitimage_ad.get_mask_width()
                limit = max(tenpercentwidth, halfpx)
                
                if abs(l_poff) > limit or abs(l_qoff) > 0.5:
                    info ("ADVICE: Apply offsets and take another keyhole acquisition image.")
                    send2tcc = True
                    msg = "        ("
                    if tenpercentwidth < halfpx and abs(l_poff) > halfpx:
                        msg += " |P|>" + str(halfpx) + "\""

                    if tenpercentwidth > halfpx and abs(l_poff) > tenpercentwidth:
                        msg += " |P|=%.0f%% of slit width" % (abs(l_poff) / l_slitimage_ad.get_mask_width() * 100.0)

                    if abs(l_qoff) > 0.5:
                        msg += " |Q|>0.5\""
                    msg += " )"  
                    info(msg)
                else:
                    info ("ADVICE: Ignore offsets and take the through-slit image.")
                    info ("Note that ALL acquisitions MUST have a through-slit image.")
                    send2tcc = False
                    if tenpercentwidth < halfpx:
                        msg = "        ( |P|< 0.5 pix and |Q|<0.5\""
                    else:
                        msg = "        ( |P|=%.0f%% of slit width and |Q|<0.5\"" % (abs(l_poff) / l_slitimage_ad.get_mask_width() * 100.0)
                        
                    msg += " )"
                    info (msg)
            else:
                info ("ADVICE: Apply offsets and take the through-slit image.")
                info ("Note that ALL acquisitions MUST have a through-slit image.")
                send2tcc = True
        elif args.mos:
            debug("...giving MOS advice")
            l_roff = 165.0 * math.tan(l_rot / 180.0 * math.pi) # linear offset at the field edges (2.75 arcmin) due to rotation
            debug("...l_roff = ", l_roff)

            slitxsize, slitysize = ad.get_min_slitsize()
            if (abs(l_poff) + abs(l_roff) > min(0.1 * slitxsize, 0.2) or
                abs(l_qoff) + abs(l_roff) > min(0.1 * slitysize, 0.2)):
                info ("ADVICE: Apply offsets and take another acquisition image.")
                send2tcc = True
            else:
                info ("ADVICE: Ignore offsets and start science immediately.")
                send2tcc = False

            msg = "        ("
            if abs(l_poff) + abs(l_roff) > 0.2:
                msg += " |P|+|R| > 0.2\""
            else:
                msg += " |P|+|R| = %.0f%% of slit width," % ((abs(l_poff) + abs(l_roff)) / slitxsize * 100.0)
                
            if abs(l_qoff) + abs(l_roff) > 0.2:
                msg += " |Q|+|R| > 0.2\""
            else:
                msg += " |Q|+|R| = %.0f%% of slit length" % ((abs(l_qoff) + abs(l_roff)) / slitysize * 100.0)
            msg += " )"
            info(msg)

        elif args.onetarget:
            debug("...giving one target advice")

            slitwidth = l_slitimage_ad.get_mask_width()
            if slitwidth is not None:
                tenpercentwidth = 0.1 * slitwidth
            
            if slitwidth is None:
                info ("ADVICE: Slit width is unknown so no advice can be given.")
                send2tcc = False
            elif abs(l_poff) > max(tenpercentwidth, 0.02) or abs(l_qoff) > 0.5:
                if ad.is_gnirs() and slitysize < 3:
                    info ("ADVICE: Slit is too narrow to determine reliable offsets.")
                    send2tcc = False
                else:
                    info ("ADVICE: Apply offsets and take another acquisition image.")
                    send2tcc = True

                msg = "        ("
                if tenpercentwidth < 0.02 and abs(l_poff) > 0.02:
                    msg += " |P|>0.02\""
                if tenpercentwidth > 0.02 and abs(l_poff) > tenpercentwidth:
                    msg += " |P|=%.0f%% of slit width" % (abs(l_poff) / slitwidth * 100.0)
                    
                if abs(l_qoff) > 0.5:
                    msg += " |Q|>0.5\""
                msg += " )"
                info(msg)
                
            else:
                if args.blindoffset:
                    info ("ADVICE: Ignore offsets and move back to the base position.")
                    info ("        Start science when guiding has been restored.")
                else:
                    info ("ADVICE: Ignore offsets and start science immediately.")

                if tenpercentwidth < 0.02:
                    info("        ( |P|<0.02\" and |Q|<0.5\" )")
                else:
                    info("        ( |P|=%.0f%% of slit width and |Q|<0.5\" )" % (abs(l_poff) / slitwidth * 100.0))

                send2tcc = False

        elif args.twotarget:
            obj1 = objectcoords.get_first_object_coords()
            obj2 = objectcoords.get_second_object_coords()
            sep = max(obj2 - obj1)
            debug("...giving two target advice")
            # this should be valid for either vertical or horizontal slits:
            debug("...separation between targets = ", sep, "pix")
               
            l_roff = sep * ad.unbinned_pixel_scale() * tan(l_rot / 180.0 * math.pi) # linear offset at the targets due to rotation
            
            debug("...l_roff = ", l_roff, "arcsec")
            debug("...l_poff + l_roff = ", abs(l_poff) + abs(l_roff))

            slitwidth = l_slitimage_ad.get_mask_width()
            if slitwidth is not None:
                tenpercentwidth = 0.1 * slitwidth 
                
            if slitwidth is None:
                info ("ADVICE: Slit width is unknown so no advice can be given.")
                send2tcc = False
            elif (abs(l_poff) + abs(l_roff)) > max(tenpercentwidth, 0.02) or abs(l_qoff) > 0.5:
                if ad.is_gnirs() and slitysize < 3:
                    info ("ADVICE: Slit is too narrow to determine reliable offsets.")
                    send2tcc = False
                else:
                    info ("ADVICE: Apply offsets and take another acquisition image.")
                    send2tcc = True

                msg = "        ("
                if (abs(l_poff) + abs(l_roff)) > max(tenpercentwidth, 0.02):
                    msg += " |P|+|R| = %.0f%% of slit width" % ((abs(l_poff) + abs(l_roff)) / slitwidth * 100.0)
                    
                if abs(l_qoff) > 0.5:
                    msg += " |Q|>0.5\""
                msg += " )"
                info(msg)

            else:
                info ("ADVICE: Ignore offsets and start science immediately.")
                if tenpercentwidth < 0.02:
                    info ("        ( |P|+|R|<0.02\" and |Q|<0.5\" )")
                else:
                    info ("        ( |P|+|R| = %.0f%% of slit width and |Q|<0.5\" )" % ((abs(l_poff)+abs(l_roff)) / slitwidth * 100.0))

                send2tcc = False
        info ("")

    #-----------------------------------------------------------------------

    default = "no"
    if send2tcc:
        default = "yes"
    
    datestring, timestring = get_timestamps()
    if ask_question("Send to TCC?", default):
        raise RuntimeError("gacq2tcc not implemented yet")

        debug("...writing obslog ", obslog)
            
        info("%s %s %s %.4f %.4f %.4f %.4f" % (datestring, timestring, full_image_name, l_poff, l_qoff, l_iaa, l_ipa), file=obslog)

        # Truncate to 4 decimal places to prevent overflows:
        gacq2tccstring = "%.4f %.4f %.4f %.4f" % (l_poff, l_qoff, l_iaa, l_ipa) 

        debug("...gacq2tccstring = ", gacq2tccstring)

        gacq2tcc(gacq2tccstring)

        info("Offsets sent: P=%.4f  Q=%.4f  IAA=%.4f  IPA=%.4f" % (l_poff, l_qoff, l_iaa, l_ipa))

    else:
        info("Offsets NOT sent to TCC")

    info ("")



if __name__ == "__main__":
    sys.exit(main(sys.argv))
