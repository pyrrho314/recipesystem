import os
import datetime
import tempfile
from copy import deepcopy

import pyfits as pf
import numpy as np
from pyraf import iraf
from pyraf.iraf import images,noao,digiphot,apphot,daofind,tables
from astrodata import AstroData
from astrodata.adutils import gemLog

class CalculateZeropoint(object):
    """
      **Syntax**
        zp = calculateZeropoint.CalculateZeropoint(
             ad, logLevel=6, logfile='', zplogfile='ZPcorr.log', extinction=None)

      :param ad:  AstroData object containing a FITS hdulist with at least one
                    image extension with EXTNAME='SCI' in the extension header.
      :type ad: AstroData object
      :param logLevel: Verbose level
      :type logLevel: integer,default is 6.
      :param logfile: log file name to replace the default value
      :type logfile: String. Default value is 'calculateZeropoint.log'
      :param zplogfile: log file containing Zero Point correction information
      :type zplogfile: Filename with default value 'ZPcorr.log'
      :param extinction: extinction correction coefficient if the image is not Gemini.
      :type extinction: Float [None].    
  
      **Description.**
      Zero Point correction Calculation.

      - Input Image FITS file with one or more SCI extensions
        having table extensions created with detectSources and 
        addReferenceCatalogs script.
      - Read FITS table extension(S) with extname 'OBJCAT'
        containing columns 'refid' and 'refmag' with at least
        one data point. These have the reference standard object name
        and magnitudes for the filter used in the exposure.
      - Each x,y position  from OBJCAT containing a non-empty standard
        in 'refid' are used to create a text file 'ZPcoo' in the working
        directory.
      - Call iraf.imexam to get FWHM
      - Call iraf.apphot.phot to get Flux, and iMag (instrumental magnitude)
      - Zero point correction: ( 'k' is  the extinction coefficient)

            ZPCorr = Reference_Magnitude - [iMag - k*(airmass-1)]
      - Image Quality (IQ) is calculated as:

            Seeing = FWHM*pixel_scale

            IQ = [mean(Seeing) / (airmass**0.6)]
      - Append keywords ZP_CORR, ZP_ERROR, IMAGE_IQ

      ERROR Propagation:

      - ZP_ERROR: At this time only the GS 'smith.cat' has error values
        for each magnitude. In the meantime we only calculate the standard
        deviation of the ZPcorr array and assign it to ZP_ERROR.

      Output file:

      - ZPcorr.log is the filename to append the new ZP values:
        "ZPCorr  ZPerror IQ     Sample Filter  Date_time". 

      **Example**
    
      >>> ad = Astrodata('corr_mrgN20100402S0047.fits')
      >>> ad.info()    # It should show tables OBJCAT and REFCAT
      >>> czp = calculateZeropoint.CalculateZeropoint(ad)

    """

    def __init__(self, ad, logLevel=6, logfile='', zplogfile='ZPcorr.log', extinction=None):

        if not logfile:
            logfile = 'calculateZeropoint.log'

        self.log = gemLog.getGeminiLog(logName=logfile, logLevel=logLevel)
        log = self.log
        log.defaultCategory(level='ALL',category='calczp')

        basename = os.path.basename(ad.filename)
        outad = deepcopy(ad)
        outad.filename = 'zpcalc_'+basename
        self.outad = outad

        self.zplogfile = zplogfile
        self.extinction = extinction


        log.info( "\n ******  CALCULATE ZERO POINT CORRECTION for: %s.  *********"%ad.filename)

    def runZP(self):

        log = self.log
        extname = 'OBJCAT'

        outad = self.outad
  
        # Check if OBJCAT extension is in the file:

        ztime = datetime.datetime

        # Open text file to append results: ZPCorr,ZPerror,IQ,Sample,Filter,Date_time
        zplogfile = self.zplogfile
        if not zplogfile:
            zplogfile = 'ZPcorr.log'
        zplog = open(zplogfile, mode='a')
       
        if zplog.tell() == 0:
            zplog.write("# ZPCorr  ZPerror IQ     Sample Filter  Date_time\n")
            

        for obj in outad['OBJCAT']:

            xtver = obj.extver()
      
            tb = obj.data
            objx = tb.field('x')
            objy = tb.field('y')
            stdn = tb.field('refid') 
            refmag = tb.field('refmag')  
            flux = tb.field('flux')  

            stdn = np.asarray(stdn)

            g = np.where(refmag != -999)         # The empty value is -999
            if len(g[0]) == 0:
                log.error( "ERROR: Calculate Zero point: No reference data found in" +\
                        " table columns 'refid','refmag' for extension: %d." % xtver)
                continue

            # look for stdn != ' '*stlen and write text file 
            # ZPcoo with x,y.

            ss = zip(objx[g],objy[g])   # For ZPcoo
            refmag = refmag[g]
            flux = flux[g]
            np.savetxt('ZPcoo',ss,fmt='%.3f',delimiter=' ')
            
            log.info("Found %d objects with reference magnitudes. File: %s[%d]"%\
                       (len(g[0]),outad.filename,xtver))
            
            # DO ZP
            
            try:
                ZPCorr, ZPerr, IQ, sample = self.do_phot(outad['SCI',xtver],
                                                         refmag, flux)
            except:
                log.error("***************** Error found in do_phot()") 
                return None

            log.info("SCI,%d: ZPCorr: %.2f error: %.2f IQ: %.2f sample:%d"%
                (xtver, ZPCorr, ZPerr, IQ, sample))
 
            # Write keyword in SCI extension
            outad.ext_set_key_value(('SCI',xtver),'ZP_CORR',ZPCorr, 'Zero Point correction')
            outad.ext_set_key_value(('SCI',xtver),'ZP_ERROR',ZPerr, 'Zero Point standard deviation')
            outad.ext_set_key_value(('SCI',xtver),'IMAGE_IQ',IQ, 'Image quality measurement')
            outad.ext_set_key_value(('SCI',xtver),'ZPSAMPLE',sample, 'ZP sample')
            

            line = '%8.2f %5.2f %7.2f %4d     %s (%s) %s[SCI,%d]\n'%\
                 (ZPCorr, ZPerr, IQ, sample, outad.filter_name(),
                  ztime.now().ctime(), outad.filename,xtver)
            zplog.write(line)

            #pp = gp.GmosPhot()
            #pp.gmosPhotometry(file)

        return outad

    def do_phot(self, hdu, refmag, flux):
        """
          Run iraf.imexamine and iraf.photpars.phot.

          INPUT:
           - hdu:      ad['SCI',extver]. Astrodata Object.
           - 
          OUTPUT:
           - ZPCorr:   Median of ZP correction array
           - ZPerror:  Standard deviation of ZParray values
           - IQ:       Image quality. 
        """  

        log = self.log

        filtername = hdu.filter_name()
        if hdu.is_type('GMOS_N'):
            ExtinctionDict = { 'g_G0301':.14, 'i_G0302':.10, 'r_G0303':.11,
                               'z_G0304':.05, 'u_G0308':.42}
            k = ExtinctionDict[hdu.filter_name()]
        elif hdu.is_type('GMOS_S'): 
            ExtinctionDict = { 'g_G0325':.18, 'i_G0327':.08, 'r_G0326':.10,
                               'z_G0328':.05, 'u_G0332':.38}
            k = ExtinctionDict[hdu.filter_name()]
        else:
            k = self.extinction
            if not k:
               log.warning("***** IMAGE is not GEMINI image and an extinction value was not supplied.")
               log.warning("***** k value set to 0.0")

        exposure_time = hdu.exposure_time()
        airmass = hdu.airmass()

        if not os.path.isfile('ZPcoo'):
            log.error("Cannot find file ZPcoo in "+os.getcwd())
            return

        # ZPcoo is input to imexam. Output is ZPimexam 
        #  
        # Create a FITS file from the current SCI extension
        
        tfile = tempfile.mktemp('.fits')   
        hdu.write(tfile)

        try:
            delete_file('ZPimexam')
            iraf.unlearn(images.imexam)
            images.imexam(input=tfile+'[1]', frame=1, logfile='ZPimexam', 
                    keeplog='yes', defkey='a', imagecur='ZPcoo', 
                    graphics ='stdgraph', use_display='no',Stdout=1)

            delete_file('ZPinfo')
            iraf.unlearn(tables.tdump)
            tables.tdump(table='ZPimexam', datafile='ZPinfo', cdfile='',
                    pfile='', columns='c15,c9,c10')
        except:
            log.error('do_phot:: In iraf.imexam or iraf.tdump')
            return

        # if object is outside GMOS FOV, the imexam will not generate a file
        if not os.path.isfile('ZPinfo'):
            log.error("File ZPinfo is not in "+os.getcwd())
            return
           
        # Read file and turn to floats. INDEF are turned to Nans 
        iFWHM, iMax, iEllip = readcols2fl('ZPinfo',unpack=True)
       
        # After centering objects using phot, using imexam to  
        # determine fwhm and peak, performs additional calculations 
        # to determine seeing and peak/coadd.

        # No ellipiticty correction needed since the images are 
        # taken unguided and will be elliptical, aperture should be big enough.
        
        aperture = 4*np.median(iFWHM) 

        # setting parameters to do photometry on the object 
        # (parameters from Inger's script)

        iraf.unlearn(iraf.apphot.datapars)
        iraf.apphot.datapars.itime     = exposure_time
        iraf.apphot.datapars.xairmass  = airmass
        iraf.apphot.datapars.otime     = hdu.ut_time()
        iraf.apphot.datapars.ifilter   = filtername

        iraf.unlearn(iraf.apphot.photpars)
        iraf.apphot.photpars.weighting = 'constant'

        # using 0 instead of 25 for final zero point calculation
        iraf.apphot.photpars.zmag = 0

        iraf.unlearn(iraf.apphot.fitskypars)
        iraf.unlearn(iraf.apphot.centerpars)
        iraf.apphot.fitskypars.salgorithm = 'gauss'
        iraf.apphot.fitskypars.khist      = 7
        iraf.apphot.centerpars.calgorithm = 'centroid'
        iraf.apphot.centerpars.cbox       = aperture
        iraf.apphot.centerpars.maxshift   = aperture*.5
        iraf.apphot.photpars.apertures    = aperture
        iraf.apphot.fitskypars.annulus    = aperture*3
        iraf.apphot.fitskypars.dannulus   = aperture

        # measuring photometry
        log.info( '*** iraf.apphot.phot: file:'+tfile+'[1]')

        try:
            delete_file('photZP')
            iraf.apphot.phot(image=tfile+'[1]',coords='ZPcoo',output='photZP', 
                 interactive='no', verify='no')
        except:
            log.error( "ERROR in iraf.apphot.phot")
            return

        delete_file('outPhotZP.tab')
        iraf.apphot.pconvert(textfile='photZP', table='outPhotZP.tab',
              fields= 'FLUX,MSKY,MAG,MERR')

        delete_file('fluxZP')
        tables.tdump(table='outPhotZP.tab', datafile='fluxZP',\
            cdfile='',pfile='',columns='')

        # Delete the temp file
        delete_file(tfile)

        iFlux, iBG, mag, merr = readcols2fl('fluxZP',unpack=True)
 
        #iMag = -2.5*np.log10(iFlux/exposure_time) - k*(airmass-1)

        #print 'mag::',mag

        g = ~np.isnan(mag)
        Mag = mag[g] -  k*(airmass-1)

        #print 'Mag::',Mag
        #print 'self.mag[g]::',self.refmag[g]

        nrefmag = len(refmag)
        if nrefmag > 0:
            refmag = refmag[g]

        if len(refmag) > 0:
            ZP = refmag - Mag

            #print 'dMag::',ZPcorr
            #print 'merr::',merr

            ZPCorr = np.median(ZP)
            ZPerror = np.std(ZP)

            log.info('do_phot. ZP,error,filter:%.2f %.2f %s'%(ZPCorr,ZPerror,filtername))
        else:
            ZPCorr = -999
            ZPerror = 0
            log.warning('do_phot. No good magnitudes to get ZP corr')

        iSeeing = iFWHM[~np.isnan(iFWHM)]*hdu.pixel_scale()

        avgS = np.mean(iSeeing)
        IQ = avgS/(airmass**0.6)
        sample = len(iSeeing)

        log.info("Seeing Ave= %s IQ= %s sample=%s"%(avgS,IQ,sample))

        delete_file('ZPcoo')
        delete_file('ZPimexam')
        delete_file('ZPinfo')
        delete_file('photZP')
        delete_file('fluxZP')
        delete_file('outPhotZP.tab')
        
        return ZPCorr, ZPerror ,IQ, sample

def delete_file(dfile):
    """ deletes a list of files"""
    try:
       os.remove(dfile)
    except: 
       pass

def unlearn_czp():
    iraf.unlearn(iraf.digiphot.apphot)
    iraf.unlearn(iraf.digiphot.daophot)
    iraf.unlearn(iraf.digiphot.apphot.phot)
    iraf.unlearn(iraf.digiphot.apphot.photpars)
    iraf.unlearn(iraf.digiphot.apphot.fitskypars)
    iraf.unlearn(iraf.noao.imred.irred.centerpars)
    iraf.unlearn(iraf.noao.imred.irred.datapars)
    iraf.unlearn(iraf.apphot.datapars)
    iraf.unlearn(iraf.apphot.centerpars)

def readcols2fl(file,unpack=False,dval=np.nan):
    """ Read columns to a float.
        if INDEF replace for dval.
    """
    ss=[]
    for ll in open(file):
        ss.append(ll.split())
    if unpack:
        ncols = len(ss[0])
        aa = []
        for s in ss:
            bb = []
            for k in range(ncols):
                if s[k] != 'INDEF':
                    bb.append(float(s[k]))
                else:
                    bb.append(dval)
            aa.append(bb)
        aa = np.asarray(aa)
        return [aa[:,k] for k in range(ncols)]
    else:
        return np.asarray(ss)

