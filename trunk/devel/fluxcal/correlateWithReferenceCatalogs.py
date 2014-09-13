import os
from math import atan2

import numpy as np
from copy import deepcopy

import pyfits as pf
from pyraf import iraf
from astrodata import AstroData
from astrodata.adutils import gemLog

class CorrelateWithReferenceCatalogs(object):
    """
      **Description**
          Read Table extensions OBJCAT and REFCAT from FITS file
          and writes into the columns 'refid' and 'refmag'.

          Correlate the x,y positions from the field objects (OBJCATs)
          with the Reference positions (REFCATs). If any common positions
          are found, write their reference id and magnitudes in the 
          ('OBJCAT',extv) table extension.

      **Syntax**
          corr = correlateWithReferenceCatalogs.CorrelateWithReferenceCatalogs 
               (ad, logLevel=6, firstPass=50, delta=7, logfile='')

      :param ad:  AstroData object containing a FITS hdulist with at least one
                  image extension with EXTNAME='SCI' in the extension header.
      :type ad: AstroData object
      :param firstPass: Maximum difference allow between object and reference coordinates
                    positions for the first matching pass.
      :type firstPass: Pixels, float. Default value is 50 pixels.
      :param delta: Maximum difference allow between object and reference coordinates
                    positions
      :type delta: Pixels, float. Default value is 6 pixels.
      :param logLevel: Verbose level
      :type logLevel: integer. [6]
      :param logfile: log file name to replace the default value
      :type logfile: String. ['detectSources.log']

      **Example**

      >>> ad = Astrodata('ds_rc_mrgN20100402S0047.fits')
      >>> ad.info()        # Should show OBJCAT and REFCAT extension
      >>> corr = correlateWithReferenceCatalogs.CorrelateWithReferenceCatalogs(ad)
      >>> adout = corr.runCorr()    # Will update the REFCAT extension
      >>> adout.write('corr_mrgN20100402S0047.fits')    # Save the modified file.

    """

    def __init__(self, ad, logLevel=6, firstPass=50, delta=7, logfile=''):

        if not logfile:
            logfile = 'selectreferences.log'

        basename = os.path.basename(ad.filename)
        outad = deepcopy(ad)
        outad.filename = 'corr_'+basename
        self.outad = outad

        self.log = gemLog.getGeminiLog(logName=logfile, logLevel=logLevel)
        log = self.log
        log.defaultCategory(level='ALL',category='corrObjRef')

        log.info( "\n  ******  CORRELATE OBJECT AND REFERENCE POSITIONS for: %s.   *********"%ad.filename)
        self.firstpass = firstPass
        self.delta = delta

    def runCorr(self):        
        """
          Match the reference position in table REFCAT with the
          position in table OBJCAT. 

          This is a 2 pass algorithm.
            - First it look for object positions within a radious of 50 pixels
              around a reference one.

            - Takes the median of the difference between the 2 positions and
              shift the reference positions by that amount.

            - Select the object positions within a radious of 6 pixels
              around a reference one.
        """
        log = self.log

        outad = self.outad

        for scix in outad['SCI']:
            
            xtver = scix.extver()

            # At this time there should be one ('OBJCAT',extver)
            # and one ('REFCAT',exter) table extensions to continue.

            objhdu = outad['OBJCAT',xtver]
            refhdu = outad['REFCAT',xtver]

            if not objhdu:
                log.warning( 'No OBJCAT positions found for SCI extension: %d'%xtver)
                continue

            tbobj = objhdu.data
            objx = tbobj.field('x')
            objy = tbobj.field('y')
            stdn = tbobj.field('refid')

            if  not refhdu:
                log.warning('No REFCAT positions found for SCI extension %d.'%xtver)
                continue

            tbref = refhdu.data
            refid = tbref.field('refid')
            refx = tbref.field('x')
            refy = tbref.field('y')
            refmag = tbref.field('refmag')

            if (len(objx)>0) & (len(refx)>0):

                # Clear these fields before writing
                tbobj.field('refid')[:] = ''      # clear the entries
                tbobj.field('refmag')[:] = [-999]             

                oindx,rindx = match_cxy(objx,refx,objy,refy, 
                              firstPass=self.firstpass, delta=self.delta,log=log)
                
                if len(rindx) == 0:
                    log.warning( "INFO: no common positions found within \
                          %d pixels in extension: %d"%(self.delta,xtver))
                    continue         # No good reference positions found

                log.info( "Found %d correlated positions. "%len(oindx))
                for i,j in zip(oindx,rindx):
                     log.info('obj:: %7.2f %7.2f ref:: %7.2f %7.2f'%\
                                (objx[i],objy[i],refx[j],refy[j]))

                tbobj.field('refid')[oindx] = refid[rindx]
                tbobj.field('refmag')[oindx] = refmag[rindx]
                log.info( "Fields 'refid' and 'refmag' in table"+\
                 " %s['OBJCAT',%d] updated." % (outad.filename,xtver))
    
        return outad     # Return the modified ad.


def match_cxy (xx, sx, yy, sy, firstPass=50, delta=None, log=None):
        """
            Match reference positions (sx,sy) with those of the 
            object catalog (xx,yy). 
            Select those that are within delta pixels from
            the object positions.
         
            firstPass:  (50) First pass delta radius.

            This matching is a 2 pass algorithm. First pass takes the
            larger delta and adjust the reference positions by the median 
            of the x,y offset. The second pass takes 'delta' value and 
            look for those x,y now closer to the reference positions.
            The units are pixels for the deltas. If you are passing
            degrees, the deltas need to be consistent.


            OUTPUT:
                 - obj_index: Index array of the objects matched.

                 - ref_index: Index array of the referecences matched.
        """

        # turn to numpy arrays
        xx, sx, yy, sy = map(np.asarray,(xx,sx,yy,sy))

        def getg(xx, sx, yy, sy, deltax=2.5, deltay=2.5):
            """ Return object(xx) and reference(sx) indices of
                common positions.
                OUTPUT
                    g:    Indices of the object position common to
                    r:    indices of the reference position
            """
            dax=[]; day=[]; g=[]; r=[]
            for k in range(len(sx)):
                gindx,= np.where((abs(xx-sx[k])<deltax) & 
                                 (abs(yy-sy[k])<deltay))
                for i in gindx:
                    dx = xx[i] - sx[k] 
                    dy = yy[i] - sy[k] 

                    # if there are multiple matches, keep only the
                    # closest one
                    if i in g or k in r:
                        if i in g:
                            first_ind = g.index(i)
                        else:
                            first_ind = r.index(k)
                        first_dist = dax[first_ind]**2 + day[first_ind]**2
                        this_dist = dx**2 + dy**2
                        if (first_dist > this_dist):
                            del dax[first_ind]
                            del day[first_ind]
                            del g[first_ind]
                            del r[first_ind]
                            dax.append(dx)
                            day.append(dy)
                            g.append(i)
                            r.append(k)
                    else:
                        dax.append(dx)
                        day.append(dy)
                        g.append(i)
                        r.append(k)
                            
                    #print i,'::','%.2f %.2f  (%.2f %.2f) %d'%(xx[i],yy[i],sx[k],sy[k],k)
            dax,day = map(np.asarray, (dax,day))
            mx = np.median(dax); stdx = np.std(dax)
            my = np.median(day); stdy = np.std(day)
            
            return np.asarray(g),np.asarray(r),mx,my,stdx,stdy 


        # Select only those standards with less than 10 pixels from objects.
        # Get the median values (mx,my) of the differences and add these
        # to the standard positions.

        #NOTE: We are setting a large delta here, we would need to see
        #      median 1st...

        ig,r,mx,my,stdx,stdy = getg(xx, sx, yy, sy, deltax=firstPass,deltay=firstPass)
        log.info( 'Median differences (x,y):%.2f %.2f, %.2f %.2f'%(mx,my,stdx,stdy)+"[First iteration]")

        if len(r) == 0:
            return ig,r
        
        # Now shift reference position by adding the median of the
        # differences. The standards are now closer to the object positions.
        sx = sx + mx   
        sy = sy + my

        # Select only those that are closer than delta or default(6.5) pixels.
        xxx = xx[ig]; yyy = yy[ig] 
        deltax=delta; deltay=delta
        if delta == None:
            deltax=2*stdx; deltay=2*stdy
        g,r,mx,my,stdx,stdy = getg (xxx, sx, yyy, sy, deltax=deltax, deltay=deltay)
        log.info( 'Median differences (x,y):%.2f %.2f %.2f %.2f'%(mx,my,stdx,stdy)+"[Second iteration]")

        if g.size == 0:
            indxy,indr=[],[]
        else:
            indxy = ig[g]
            indr = r
        
        return indxy, indr
