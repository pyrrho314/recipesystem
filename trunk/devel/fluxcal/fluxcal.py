#! /usr/local/bin/python
import os
import sys
import numpy as np
import imp
from copy import deepcopy

from pyraf import iraf
#from pyraf.iraf import display,tvmark
from pyraf.iraf import tvmark
from astrodata import AstroData
from astrodata.adutils import gemLog
import gempy.science.geminiScience as gs
#import geminiScience as gs

from detectSources import DetectSources
from addReferenceCatalogs import AddReferenceCatalogs
from correlateWithReferenceCatalogs  import CorrelateWithReferenceCatalogs
from calculateZeropoint import CalculateZeropoint


class Fluxcal(object):
    """ 
      TODO:
        -A function to look for standards within the field
        - Add a parameter to not write an output file.

      Test Script to drive
        - class DetectSources. Creates [OBJCAT,ext] table extension.
        - class SelectReferences. Create [REFCAT,ext] table extension.
        - class CorrelateObjRef. Look for common position in OBJCAT and REFCAT tables
        - class CalculateZP. Calculate ZP correction from common position array.

        INPUTS:
          - ad             : Input AD object  
          - addBPM=True    : add BPM
          - sigma=0.0      : sigma value for detectSource
          - threshold=2.5, : value for DS
          - fwhm=5.5,      : value for DS
          - logfile='fluxcal.log'

      ::

        Example:
          import flux_cal as fc
          from astrodata import AstroData

          ad = AstroData('zzz.fits')
          ff = Fluxcal(ad)    # Create object

    """


    def __init__(self,ad, addBPM=True, sigma=0.0, threshold=2.5, fwhm=5.5,
               logfile='fluxcal.log'):
       
        """
           Run all the functions to calculate ZP
           with the default parameters. 
           If you want to change some of then, you could
           run the individual functions.
        
           -  DetectSources(ad)
           -  AddReferenceCatalogs(ad)
           -  CorrelateWithReferenceCatalogs (ad)
           -  CalculateZeroPoint(ad)
        """
        
        basename = os.path.basename(ad.filename)
        adout = deepcopy(ad)
        adout.filename = 'fc_'+basename
        self.adout = adout
        self.addBPM=addBPM
        self.sigma=sigma
        self.threshold=threshold
        self.fwhm=fwhm

        self.logfile = logfile
        self.log=gemLog.getGeminiLog(logName=logfile,logLevel=6)
        log = self.log

        log.info( "\n  ==================== FLUXCAL for "+ad.filename+" ====================")

    def runFC(self):
        log = self.log
        adout = self.adout
        logfile = self.logfile

        adout = self._detsources()

        cc   = AddReferenceCatalogs(adout, logfile=logfile)
        adout = cc.getRefs()

        corr = CorrelateWithReferenceCatalogs (adout, logfile=logfile)
        adout = corr.runCorr()

        czp  = CalculateZeropoint(adout, logfile=logfile)
        adout = czp.runZP()
       
        self.adout = adout

        ofile = 'zp_'+adout.filename
        log.info("writing output fits file: "+ofile)
        if adout: adout.write(ofile, clobber=True)        


    def _detsources(self):
        """ Local function
        """ 
        if self.addBPM:
           self.addbpm()
           #print 'AFTER addbpm:',self.adout.info()
        self.log.info("\n DS DS DS DS using: SIGMA THRESHOLD FWHM:: "+\
                 str(self.sigma)+','+\
                 str(self.threshold)+','+str(self.fwhm)+'\n')
        dd = DetectSources   (self.adout, sigma=self.sigma,
                    threshold=self.threshold,fwhm=self.fwhm,
                    logfile=self.logfile)
        adout = dd.runDS()
        self.adout = adout
        return adout

    ds = _detsources        # an alias

    def addbpm(self):
        """Add BPM extension(s) to the AD object.
        """
        bpm = self.getBPM(userbpm='')
        #bpm = AstroData('bpmgmosn22_123.fits')
        self.log.info("Adding BPM extension ("+bpm.filename+") to "+\
                      str(self.adout.filename))

        self.adout=gs.add_bpm(adInputs=self.adout, BPMs=bpm, suffix='bbb', logLevel=6)
        self.addBPM = False    # In case we are running detsources now.

    def getoxy(self,extn=1):
        """Get x,y,mag,err from the current file
           in the table extension ['OBJCAT',extn]; otherwise
           give a filename

           INPUT:
              file: FITS filename, if not default

              extn: Default 1. Extension number, needs
                    to be less than the number of
                    SCI extensions in the file.
           OUTPUT:
              x,ymag,err:   Arrays from the ['OBJCAT',extn] data.
                             

        """
        xtnad = self.adout['OBJCAT',extn]
        if xtnad:
            x = xtnad.data.field('x')
            y = xtnad.data.field('y')
            flux = xtnad.data.field('flux')
        else:
            raise RuntimeError,\
            " **ERROR: FITS extension does not exist: %s['OBJCAT',%d]."%(file,extn)
        
        return x,y,flux
          
    def getrxy(self,extn=1):
        """Get x,y,mag, from the current file
           in the table extension ['REFCAT',extn]; otherwise
           give a filename

           INPUT:
              file: FITS filename, if not default.

              extn: Default 1. Extension number, needs
                        to be less than the number of
                        SCI extensions in the file.
           OUTPUT:
              x,y,mag:   Arrays from the ['REFCAT',extn]
                             

        """
        xtnad = self.adout['REFCAT',extn]
        if xtnad:
            x = xtnad.data.field('x')
            y = xtnad.data.field('y')
            refmag = xtnad.data.field('refmag')
        else:
            raise RuntimeError,\
            " **ERROR: FITS extension does not exist: %s['REFCAT',%d]."%(file,extn)
        
        return x,y,refmag

    def display(self,extn=None):
        """Display SCI extensions from the Astrodata object
           and mark objects and reference stars found in the
           input file tables OBJCAT and REFCAT.

           DS9 must be running.
        """
        try:
            from stsci.numdisplay import display as ndisplay
        except ImportError:
            from numdisplay import display as ndisplay

        #file = self.file
        #ad = AstroData(file)

        ad = self.adout
        nscext = ad.count_exts('SCI')

        if extn != None:
            si,ei = extn,extn+1
        else:
            si,ei = 1,nscext+1
            
        fr = 1
        for xt in range(si,ei):
            #ff = '%s[%d]' % (file,xt)
            ff = ad['SCI',xt].data
            
            iraf.set (stdimage='imt47')
            ndisplay(ff, frame=fr)
    
            xtnad = ad['OBJCAT',xt]
            if xtnad:
                x = xtnad.data.field('x')
                y = xtnad.data.field('y')
                print 'Marking %d objects with red.' % len(x)
                np.savetxt('/tmp/xy.txt', zip(x,y), fmt='%.2f', delimiter=' ')
                tvmark(frame=fr,coords='/tmp/xy.txt',color=204)
            else:
                print "No stars found in extension %d."%xt
        
            xtnad = ad['REFCAT',xt]
            if xtnad:
                x = xtnad.data.field('x')
                y = xtnad.data.field('y')
                print 'Marking %d refstars with green.' % len(x)
                np.savetxt('/tmp/xy.txt', zip(x,y), fmt='%.2f', delimiter=' ')
                tvmark(frame=fr,coords='/tmp/xy.txt',color=205)
            else:
                print "No reference stars found in extension %d."%xt

            fr += 1

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

        ad = self.adout

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


if __name__ == '__main__':

    import sys,glob

    #ad = AstroData('/home/nzarate/zp/zzz.fits')
    #ff = Fluxcal(ad)    # Create object
    #ff.runFC()       # run the scripts
    

    for ifile in glob.glob('[g,m]*.fits'):
        if 'rgS20110125S' in ifile: continue
        if 'zp_' in ifile: continue
        ad = AstroData(ifile)
        if ad.is_type('CAL'):
            print "\nWARNING: 'CAL' type for file:",ifile," not supported."
            continue
        print "\n >>>>>>>>>>>>>>>>>>FluxCAL for:",ifile 
        ff = Fluxcal(ad,addBPM=False)    # Create object
        ff.runFC()
