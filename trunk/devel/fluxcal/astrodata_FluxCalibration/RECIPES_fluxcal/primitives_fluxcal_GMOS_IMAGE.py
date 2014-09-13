import sys, StringIO, os

from astrodata.adutils import gemLog
from astrodata import Descriptors
from astrodata.data import AstroData
from gempy import gemini_tools as gemt
from primitives_GMOS import GMOSPrimitives
from astrodata.adutils.gemutil import pyrafLoader
import primitives_GEMINI


import numpy as np
import pyfits as pf
import shutil

# from devel.fluxcal import fluxcal
# import flux_cal as fc
from detectSources import DetectSources
from addReferenceCatalogs import AddReferenceCatalogs
from correlateWithReferenceCatalogs  import CorrelateWithReferenceCatalogs
from calculateZeropoint import CalculateZeropoint


log=gemLog.getGeminiLog()

# NOTE, the sys.stdout stuff is to shut up gemini and gmos startup... some primitives
# don't execute pyraf code and so do not need to print this interactive
# package init display (it shows something akin to the dir(gmos)
import sys, StringIO, os
SAVEOUT = sys.stdout
capture = StringIO.StringIO()
sys.stdout = capture
sys.stdout = SAVEOUT


class FluxCalException:
    """ This is the general exception the classes and functions in the
    Structures.py module raise.
    """
    def __init__(self, msg="Exception Raised in Recipe System"):
        """This constructor takes a message to print to the user."""
        self.message = msg
    def __str__(self):
        """This str conversion member returns the message given by the user (or the default message)
        when the exception is not caught."""
        return self.message

class GMOS_IMAGE_fluxcal_Primitives(GMOSPrimitives):
    """ 
    This is the class of all primitives for the GMOS level of the type 
    hierarchy tree.  It inherits all the primitives to the level above
    , 'GEMINIPrimitives'.

    The list of primitives in this class is: 
    """
    astrotype = 'GMOS_IMAGE'
    
    def init(self, rc):
        GMOSPrimitives.init(self, rc)
        return rc
     

    def detectSources(self, rc):
        """ detectsources primitive
        """
        adOuts = []
        try:
            log.info( "Starting primitive detectSources")
            for ad in rc.get_inputs(style='AD'):
                 ds   = DetectSources   (ad, logfile='', sigma=rc['sigma'], 
                                         threshold=rc['threshold'], fwhm=rc['fwhm'])
                 adout = ds.runDS()
                 adOuts.append(adout)
            rc.report_output(adOuts)
        except:
            raise FluxCalException("Problems with detectSources")
      
        yield rc


    def addReferenceCatalogs(self, rc):
        """ addReferenceCatalogs primitive
        """
        adOuts = []
        try:
            log.info(" primitive addReferenceCatalogs")
            for ad in rc.get_inputs(style='AD'):
                 cc   = AddReferenceCatalogs(ad, logfile='', 
                                             catalogName=rc['catalogName'])
                 adout = cc.getRefs()
                 adOuts.append(adout)
            rc.report_output(adOuts)
        except:
            raise FluxCalException("Problems with addReferenceCatalogs")
                 
        yield rc
        
    def correlateWithReferenceCatalogs(self, rc):
        """ correlateWithReferenceCatalogs primitive
        """
        adOuts = []
        try:
            log.info(" primitive correlateWithReferenceCatalogs")
            for ad in rc.get_inputs(style='AD'):
                 corr = CorrelateWithReferenceCatalogs(ad, logfile='',
                                                       delta=rc['delta'],
                                                       firstPass=rc['firstPass'])
                 adout = corr.runCorr()
                 adOuts.append(adout)
            rc.report_output(adOuts)
        except:
            raise FluxCalException("Problems with correlateWithReferenceCatalogs")
                 
        yield rc
        
    def calculateZeropoint(self, rc):
        """ calculateZeropoint primitive
        """
        adOuts = []
        try:
            log.info(" primitive calculateZeropoint")
            for ad in rc.get_inputs(style='AD'):
                 czp  = CalculateZeropoint(ad, logfile='', 
                                           zplogfile=rc['zplogfile'],
                                           extinction=rc['extinction'])
                 adout = czp.runZP()
                 adOuts.append(adout)
            rc.report_output(adOuts)
        except:
            raise FluxCalException("Problems with calculateZeropoint")

        yield rc

    def fluxcal(self, rc):
        """
          Fluxcal drives:
            - DetectSources
            - AddReferenceCatalogs
            - CorrelateWithReferenceCatalogs
            - CalculateZeropoint
        """

        logfile = 'fluxcal.log'
        adOuts = []

        try:
            log.status('*STARTING*  Running fluxcal functions')
            for ad in rc.get_inputs(style='AD'):

                ds   = DetectSources   (ad, logfile=logfile, sigma=rc['sigma'], 
                                        threshold=rc['threshold'], fwhm=rc['fwhm'])
                adout = ds.runDS()

                cc   = AddReferenceCatalogs(adout, logfile=logfile, 
                                            catalogName=rc['catalogName'])
                adout = cc.getRefs()

                corr = CorrelateWithReferenceCatalogs (adout, logfile=logfile,
                                                       delta=rc['delta'])
                adout = corr.runCorr()

                czp  = CalculateZeropoint(adout, logfile=logfile, 
                                           zplogfile=rc['zplogfile'],
                                           extinction=rc['extinction'])
                adout = czp.runZP()

                adOuts.append(adout)

            rc.report_output(adOuts)

            log.status('Fluxcal completed successfully')
        except:
            raise FluxCalException("Problems with fluxcal")

        yield rc

