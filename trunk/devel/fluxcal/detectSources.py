import os
import re
import imp
from copy import deepcopy

import pyfits as pf
import pywcs
import numpy as np
from astrodata import AstroData
from astrodata.adutils import gemLog
from detSources import detSources


class DetectSources(object):
    """
        Find x,y positions of all the objects in the input
        image. Append a FITS table extension with position
        information plus columns for standard objects to be
        updated when position from addReferenceCatalogs -if any are 
        found for this field.

        **Syntax:**
           ds = detectSources.DetectSources( ad, logLevel=6, logfile='')

        :param ad:
              AstroData object containing a FITS hdulist with at least one
              image extension with EXTNAME='SCI' in the extension header.
        :type ad: AstroData object
        :param sigma: The mean of the background value. If nothing is passed,
              detSources will run background() to determine it.
        :type sigma: Number [0.0]

        :param threshold: Threshold intensity for a point source - should generally
                  be 3 or 4 sigma above background RMS"[1]. It was found that 2.5
                  works best for IQ source detection.
        :type threshold: Number [2.5]

        :param fwhm: FWHM to be used in the convolve filter". This ends up playing
                     a factor in determining the size of the kernel put through 
                     the gaussian convolve.
        :type fwhm: Number [5.5]

        :param logLevel:
              Verbose level
        :type logLevel: integer,[6]
        :param logfile: 
               log file name to replace the default value
        :type logfile: String. ['detectSources.log']

        **INPUTS:**
           An AstroData object containing an IMAGE extension with
           EXTNAME value of 'SCI' and an integer EXTVER value.

        **OBJCAT table:**
          A FITS table ('OBJCAT',extver) is appended to the input 
          Astrodata Object. This table contains the columns:
     
          - 'id':    Unique ID. Simple running number.
          - 'x':     x coordinate of the detected object. 1-based
          - 'y':     x coordinate of the detected object, both x and y are
                     with respect to the lower left corner. They are 1-based
                     values, i.e. as in ds9.
          - 'ra':    ra values, in degrees. Calculated from the fits header WCS
          - 'dec':   dec values in degrees. ditto
          - 'flux':  Flux given by the gauss fit of the object
          - 'refid': Reference ID for the reference star found in the field
          - 'refmag': Reference magnitude. 'refid' and 'refmag' will be fill
                      by the 'correlateWithReferenceCatalogs' function.       

        **Methods:**
            runDS()    - Run the object detection returning the
                         astrodata object with the table.

        **EXAMPLE**
        
        >>> ad = Astrodata('mrgN20100402S0047.fits')
        >>> ds = detectSources.DetectSources(ad)
        >>> adout = ds.runDS()
        >>> ad.info()      # will have the BINTABLE extension

    """
        
    def __init__(self, ad, sigma=0.0, threshold=2.5, fwhm=5.5, logLevel=6, logfile=''):

        if 'CAL' in ad.types:
            raise RuntimeError, " **** AD type 'CAL' not supported  by DetectSources****"

        basename = os.path.basename(ad.filename)

        self.sigma = sigma
        self.threshold = threshold
        self.fwhm = fwhm

        if not logfile:
            logfile = 'detectSources.log'

        self.log = gemLog.getGeminiLog()
        log = self.log
        log.defaultCategory(level='ALL',category='DetecSources')
 
        outad = deepcopy(ad)
        outad.filename = 'ds_'+basename
        self.outad = outad
        log.info( "\n ******  DETECTING  SOURCES  for: %s *************"%ad.filename)

    def runDS(self):        
        """ Do the actual object detection.
            - Create the table OBJCAT
            - return the output Astrodata object

        """
        log = self.log
        extname = 'OBJCAT'

        outad = self.outad

        for scix in outad['SCI']:

            xtver = scix.extver()
  
            # Mask the non illuminated regions.
            if outad['BPM',xtver]:
                sdata = scix.data
                bpmdata = outad['BPM',xtver].data
                if bpmdata.shape != sdata.shape:       # bpmdata is already trimmed.
                    try :                            # See if DATASEC is in the header
                        dsec = scix.data_section()
                    except:
                        log.error("*** ERROR: DATASEC not found in SCI header.")
                        log.error("*** Cannot masked SCI: size(BPM) != size(SCI)."+\
                                   " bpmsize: "+str(bpmdata.shape))
                        log.error(" *** SCI data not masked.")
                    else:
                        # bpmdata is already trimmed.
                        s,e = map(int, dsec.split(',')[0][1:].split(':'))
                        #Trim number of columns to match the bpm data.
                        sdata = sdata[:,s-1:e]
                else:
                    bpmdata = np.where(bpmdata==0,1,0)
                    scix.data = sdata*bpmdata
              
            self.findObjects(scix)

            sciHeader = scix.header
            if len(self.x) == 0:
                log.warning( " **** WARNING: No objects were detected: Table OBJCAT, not created")
                continue
       
            wcs = pywcs.WCS(sciHeader)

            # Convert pixel coordinates to world coordinates
            # The second argument is "origin" -- in this case we're declaring we
            # have 1-based (Fortran-like) coordinates.

            xy = np.array(zip(self.x, self.y),np.float32)
            radec = wcs.wcs_pix2sky(xy, 1)

            ra,dec = radec[:,0],radec[:,1]

            nobjs = len(ra)
            log.info("Found %d sources for field in  %s['SCI',%d]"%\
                (nobjs ,outad.filename,xtver))
            
            if outad[extname,xtver]:
                log.info('Table already exists,updating values.')
                tdata = outad[extname, xtver].data
                theader = outad[extname, xtver].header
                tdata.field('id')[:]    = range(len(ra))
                tdata.field('x')[:]     = self.x
                tdata.field('y')[:]     = self.y
                tdata.field('ra')[:]    = ra
                tdata.field('dec')[:]   = dec
                tdata.field('flux')[:]  = self.flux
            else:
                #colsdef = self.define_Table_cols(ra, dec, flux, ellip, fwhm)
                colsdef = self.define_Table_cols(ra, dec, self.flux)
                tbhdu = pf.new_table(colsdef)         # Creates a BINTABLE

                th = tbhdu.header
                
                tabad = AstroData(tbhdu)
                tabad.rename_ext("OBJCAT", xtver)
            
                outad.append(tabad)

        return outad

    def getBPM(self,userbpm):
        """
          Get the BPM file from the flixcal directory and 
          return its the AstroData object.

          :param userbpm:  BPM file suplied by the user. Its size
                           should match the SCI extension image size.

          :rparam ad: The bpm file astrodata object.
          If 'userbpm' is not empty just return its AstroData object

          NOTE: We still need to create: bpmgmoss11_123.fits.gz   # unbinned mef files
                                         bpmgmosn11_123.fits.gz
        """

        log = self.log
 
        if userbpm:
            return AstroData(userbpm)

        ad = self.ad

        ccdsum = ad['SCI',1].header.get('CCDSUM')
        ccdsum = ccdsum.replace(' ','')     # makes the string '11' or '22'
        nextn = ad.count_exts('SCI')
        indx = ccdsum+str(nextn)         # form {11,22}1 or {11,22}3

        if ad.is_type('GMOS_S'):
            dqfiles = {'223':'bpmgmoss22_123.fits.gz','113':'bpmgmoss11_123.fits.gz',
                       '221':'bpmgmoss22_mosaic.fits.gz','111':'bpmgmoss11_mosaic.fits.gz'}
            try:
                dq = dqfiles[indx]
            except:
                log.warning('getBPM: CCDSUM value not supported with BPM: '+indx)
                return None
        elif ad.is_type('GMOS_N'):
            dqfiles = {'223':'bpmgmosn22_123.fits.gz','113':'bpmgmosn11_123.fits.gz',
                       '221':'bpmgmosn22_mosaic.fits.gz','111':'bpmgmosn11_mosaic.fits.gz'}
            try:
                dq = dqfiles[indx]
            except:
                log.warning('getBPM: CCDSUM value not supported with BPM: '+indx)
                return None
        # Find where are these BPM files

        fp, pathname, description = imp.find_module('detectSources')
        fname = os.path.join(os.path.dirname(pathname),dq)

        return AstroData(fname)

    def findObjects(self, hdu):
        """ Small layer to call detSources.
            returns x,y,flux
        """
        log = self.log
        objs = detSources(hdu=hdu, sigma=self.sigma, threshold=self.threshold,
                          fwhm=self.fwhm)
        if len(objs) == 0:
            self.x, self.y, self.flux = [],[],[]
            return 
         
        self.x, self.y ,self.flux = [np.asarray(objs)[:,k] for k in [0,1,2]]
        log.info("  ... found %d objects."%len(self.x))

    def define_Table_cols(self, ra, dec, flux):
        """
          Define table columns for OBJCAT table

          **INPUT:**
             - ra,dec:    ra,dec array from x,y position in the image field.

             - self.x, self.y and self.flux
        """
        from pyfits import Column

        # Define columns


        nlines = len(ra)

        c1 = Column (name='id',format='J',array=range(nlines))
        c2 = Column (name='x',format='E',array=self.x)
        c3 = Column (name='y',format='E',array=self.y)
        c4 = Column (name='ra',format='E',array=ra)
        c5 = Column (name='dec',format='E',array=dec)
        c6 = Column (name='flux',format='E',array=flux)
        #c7 = Column (name='ellip',format='E',array=ellip)
        #c8 = Column (name='fwhm',format='E',array=fwhm)
        c7 = Column (name='refid',format='22A',array=['']*nlines) # Placeholders
        c8 = Column (name='refmag',format='E',array=[-999]*nlines)  # for corr

        #return pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
        return pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8])

def delete_file(file):
        if os.path.isfile(file):
            os.remove(file)

