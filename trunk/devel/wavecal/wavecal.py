import os
import sys

import numpy as np
from matplotlib import pyplot as pl

from astrodata import Lookups
from astrodata import AstroData
import spec_utils as spu
import gfit
#import inspect

# Put this in a wcal_config.py
INSTRUMENT = {'ALL': {'clip': 4,'cradius': 12, 'fitfunction': 'chebyshev', 'fitorder': 4,
                      'fwidth': 10, 'match': -3, 'minsep':5, 'nsum': 10, 'ntmax': 50,
                      'bigger_linelist': None,
                     },
              ('GMOS','ls'):  {'linelist':'cuar.dat', 'nbins':6,'bins_step':12, 'best':1.5,
                        'bigger_linelist': 'cuar_big.dat',
                       }, 
              ('GMOS','ifu'):  {'linelist':'cuar.dat', 'nbins':6,'bins_step':12, 'best':1.5,
                        'bigger_linelist': 'cuar_big.dat',
                       }, 
              ('GMOS','mos'):  {'linelist':'cuar.dat', 'nbins':4,'bins_step':12, 'best':1.0,
                         'ntmax':30, 'nsum':5,
                        'bigger_linelist': 'cuar_big.dat',
                       }, 
              'GNIRS': {'linelist':'lowresargon.dat', 'nbins':4,'bins_step':12, 'best':1.0,
                        'bigger_linelist': 'argon.dat','ntmax': 15,
                        },
              ('NIFS','JH'):  {'linelist':'argon.dat', 'nbins':2,'bins_step':12, 'best':1.0, 
                              'match': -10,
                              'ntmax': 25,
                              },
              ('NIFS','HK'):  {'linelist':'ArXe_K.dat', 'nbins':2,'bins_step':12, 'best':1.5, 
                              'clip': 3,
                              'match': -6,
                              'ntmax': 25,
                              'bigger_linelist': 'argon.dat',
                              },
              ('NIFS','ZJ'):  {'linelist':'argon.dat', 'nbins':4,'bins_step':18, 'best':1.0, 
                              'match': -10,
                              'minsep':15,
                              'ntmax': 20,
                              'bigger_linelist': 'argon.dat',
                              },
              ('F2','JH'):  {'linelist':'lowresargon.dat', 'nbins':6,'bins_step':18, 'best':1.0, 
                              'clip':3.5,
                                # This value work  OK for exposure time of
                                # 10 secs or less and MOSPOS value of 2pix-slit
                              'bigger_linelist': 'argon.dat',
                            },
              ('F2','J_'):  {'linelist':'lowresargon.dat', 'nbins':6,'bins_step':12, 'best':1.0, 
                              'clip':3.5,
                              'ntmax':20,
                              'bigger_linelist': 'argon.dat',
                            },
              ('F2','HK'):  {'linelist':'lowresargon.dat', 'nbins':4,'bins_step':12, 'best':1.0, 
                              'match': -6,
                              'clip':3.5, 
                              'bigger_linelist': 'argon.dat',
                            },
              ('F2','H_'):  {'linelist':'lowresargon.dat', 'nbins':4,'bins_step':12, 'best':1.0, 
                              'match': -6,
                              'clip':3.5, 
                              'bigger_linelist': 'argon.dat',
                            },
              ('F2','Ks'):  {'linelist':'argon.dat', 'nbins':4,'bins_step':12, 'best':1.5, 
                              'match': -6,
                              'ntmax': 20,
                              'bigger_linelist': 'argon.dat',
                            },
              ('NIRI'):  {'linelist':'lowresargon.dat', 'nbins':6,'bins_step':12, 'best':1.5, 
                              'bigger_linelist': 'argon.dat',
                         },
             }

class Wavecal(object):
    """The wavelength calibration scripts allows the user to obtain a mapping from 
       pixel to wavelength coordinates given an ARC image from a reduced Gemini FITS
       file. The data requirements to run the Wavecal are:

       - A prepared or reduced GEMINI ARC FITS file. 
       - A set of CRPIX, CRVAL and CDELT (CDn_n, n:1,2) values in the
         dispersion direction.
       - A linelist data file with reference wavelength covering the range
         for the arc spectral lines in the image.


    """

    def __init__(self,ad, linelist=None,  extv=1,
                  fitfunction=None, fitorder=None, ntmax=None, fwidth=None,
                  cradius=None, match=None, minsep=None, clip=None, nsum=None,
                  nbins=None, best=None, bins_step=None, debug=False):
        """
      
        ad: AstroData object containing a Gemini Instrument Arc image.

        extv = 1. 
            The FITS extension number in the SCI extension name list.

        linelist = (None).
             An instrument dependent list of reference wavelengths.
             It should be one entry per line.

        fwidth = 10.
            Full-width at the base (in pixels) of features to be identified.
        
        cradius = 12.
            The maximum distance, in pixels, allowed between a line position
            and  the  initial  estimate  when  defining a new line. The user
            will need to increase cradius if wide slits (1 arcsec or  wider)
            are  used  with  the  detector  unbinned  or  binned by 2 in the
            spectral direction.
        
        minsep = 2.
            The    minimum   separation,   in   pixels,   allowed    between 
            line positions when defining a new line.  The user will need  to
            increase  minsep if wide slits (1 arcsec or wider) are used with
            the detector unbinned or binned by 2 in the spectral direction.

        ntmax = (None). 
            An instrument dependent value. The maximum number of peaks to
            returns from the finding peaks algorithm. The selection is 
            based on peak intensity ranking.
        
        match = -6.  
            The  maximum  difference for a match between the line coordinate
            derived from the dispersion function and  a  coordinate  in  the
            coordinate  list.   Positive values are in user coordinate units
            and  negative values are in units of pixels.  The user will need
            to increase the absolute value of match if wide slits (1  arcsec
            or  wider) are used with the detector unbinned or binned by 2 in
            the spectral direction.
        
        fitfunction = "chebyshev"  (legendre|chebyshev|cubic)
            The function to be fit to user coordinates as a function of  the
            pixel  coordinates.   The  choices  are "chebyshev", "legendre",
            and "cubic" to fot a cubic Spline.
        
        fitorder = 4
            Order of the fitting function.   The  order  is  the  number  of
            polynomial terms (coefficients) or the number of spline pieces.
     
        clip = (none)
            An instrument dependent value. The number of sigma units use to
            reject peak coordinate positions when fitting (peaks,wavelengths)
            with the 'fitfunction'.

        nbins = None  (Value between 4 and 6 depending on the instrument)
            The number of bins to divide the reference range to be used by 
            the matching algorithm.

        best = None (Value between 1.0 and 1.5 depending on the instrument)
            A value related to the matching fraction between the reference lines
            matched and the number of peaks.

        bins_step = None  (Value between 12 and 15 depending on the instrument)
            The step value to subdivide the array of triples composing one bin.
            If the number of triples is greater than 800 we subdivide by 12
            for example.

        *Example*

        >>> from wavecal import Wavecal
        >>> ad = AstroData('wnN20130705S0176.fits')
        >>> wc = Wavecal(ad)   # Create a Wavecal object
        >>> wc.wavecal()  # Perform the wavelength calibration
        >>> wc.info()     # A summary of parameters used an results
        >>> wc.plot_ref() # A subplot of the image middle line plus a reference arc
                          # subplot showing the associated matched lines.
        >>> wc.features   # A list of peaks coordinates, the wavelength fitted value
                          # and the reference wavelength from the reference list.
    
        """

        self.ad = ad
        self.filename = ad.filename

        # Set up a dictionary with parameters values at object creation 
        # time which are replacing the default values if they are 
        # different from None.
        param={}
        param['extv']     = extv
        param['ntmax']    = ntmax
        param['fitfunction'] = fitfunction
        param['fitorder'] = fitorder
        param['match']    = match
        param['minsep']   = minsep
        param['nsum']     = nsum
        param['cradius']  = cradius
        param['fwidth']   = fwidth
        param['clip']     = clip
        param['linelist'] = linelist
        param['best']     = best
        param['nbins']    = nbins
        param['bins_step']= bins_step
        param['debug']= debug

        # Instantiate a subclass for a specific instrument.
        if ad.is_type('F2_LS_ARC'):
            instr_obj = F2(ad,param)
        elif ad.is_type('GMOS_SPECT'):
            instr_obj = GMOS(ad,param)
        elif ad.is_type('GNIRS_SPECT'):
            instr_obj = GNIRS(ad,param)
        elif ad.is_type('NIRI_SPECT'):
            instr_obj = NIRI(ad,param)
        elif ad.is_type('NIFS_SPECT'):
            instr_obj = NIFS(ad,param)
        else:
            raise ValueError('Input file is not supported: '+ad.filename)


        # Copy subclass members to the parent's since they will
        # be used by parent methods.
        self.reference = instr_obj.reference
        self.ntmax =     instr_obj.ntmax
        self.nsum =      instr_obj.nsum
        self.minsep =    instr_obj.minsep
        self.cradius =   instr_obj.cradius
        self.match =     instr_obj.match
        self.dispaxis =  instr_obj.dispaxis
        self.debug =  debug
        wcs = instr_obj.wcs
        self.wcs = wcs
    

        fpath = os.path.dirname(spu.__file__)
        linelist_file = instr_obj.param['linelist']
        linelist_file = os.path.join(fpath,linelist_file)
        #linelist_file = os.path.join('/home/nzarate/wcal/',linelist_file)
        linelist = np.loadtxt(linelist_file,usecols=(0,))

        crpix = wcs['crpix']
        crval = wcs['crval']
        cdelt = wcs['cdelt']

        # Setup a linear pixel--> wavelength transformation
        # Equivalent to   wave = (x-crpix)*cdelt + crval
        z = gfit.Gfit([],[],fitname='polynomial',order=1)
        z.coeff[0] = cdelt
        z.coeff[1] = crval - crpix*cdelt
        instr_obj.z = z
         
        if cdelt < 0:
            linelist = linelist[::-1]    # Change to decreasing order

        self.cuar = linelist 
        instr_obj.cuar = linelist
        self.ad = ad

        self.extv = extv
        self.imdata = instr_obj.imdata
        self.instr_obj = instr_obj


    def info(self):
        """
          Display information about input ARC data, input parameters
          seeting and wavelength calibration information.
        """
        fn = self.ad.filename
        print 'Input File: ',os.path.basename(fn),', Filter: ',self.ad.filter_name()
        mode = self.instr_obj.instrument_mode
        print 'Instrument: ',self.ad.instrument(), ', Mode:',mode
        par = self.instr_obj.param
        s = 'match: '+str(par['match'])
        print s+' '*(32-len(s))+\
            '... Maximum difference (fit - linelist) for a given line.'
        s = 'minsep: '+str(par['minsep'])
        print s+' '*(32-len(s))+\
            '... Minimum separation between lines (pixel).'         
        s = 'nsum:   '+str(par['nsum'])
        print s+' '*(32-len(s))+\
            '... Number of rows to sum about the middle line.'          
        s = 'ntmax:  '+str(par['ntmax'])
        print s+' '*(32-len(s))+\
            '... Maximum number of peaks to find.'
        s = 'clip:   '+str(par['clip'])
        print s+' '*(32-len(s))+\
            '... The number of sigma units fitting (peaks,wavelengths).'
        s =  'linelist: '+par['linelist']
        print s+' '*(max(0,32-len(s)))+\
            '... File with list of reference wavelengths.'
        print 'WCS: ',self.wcs
        print '*** Matching algorithm parameters ***'
        print 'nbins:',par['nbins']
        print 'bins_step:',par['bins_step']
        print 'best:',par['best']
        if hasattr(self, 'z'):
           print '*** Wavecal mapping function ***'
           print 'Function and order:',self.z.fitname,self.z.order
           print 'Coefficients: ',self.z.coeff
           print 'Number of features found:', len(self.xpeaks)
           print 'Best:',self.best

    def mos_info(self):
    
        fn = self.ad.filename
        adm = AstroData(fn)
        for extn in range(1,adm['SCI'].count_exts()+1):
            ad = adm['SCI',extn]
            hdr = ad.header
            sz = ad.data.shape
            crpix = hdr['crpix1']
            crval = hdr['crval1']
            cdelt = hdr['cd1_1']
            p2w = lambda x: (x-crpix)*cdelt + crval
            print extn,sz,'(%.2f %.2f)'%(p2w(1),p2w(sz[1]))
        

    def wavecal(self):
        """ 
          *Algorithm*

             - Given an input image, collapse and take the mean of 
               self.nsum rows from the middle of the image.
             - In this line, find upto self.ntmax peak coordinates 
               corresponfing to the arc spectrum peak positions.
             - Using these peaks and an input linelist with reference
               wavelengths start the matching ratio algorithm to find
               the correct peak position with its linelist wavelength.
             - Fit a function to set of (pix_array, user_array) where
               'user' is an alias for 'wavelength'. The member self.z
               contains the object with a Gfit class members with 
               all the fit information.
               For the Matching ratio algorithm please see the docstring
               for the class MatchLines.
        """

        # Execute the Instrument specific wavecal()
        self.instr_obj.wavecal()

        # Copy the new member values to the parent's. 
        self.xpeaks = self.instr_obj.xpeaks
        self.lpix   = self.instr_obj.lpix
        self.pix    = self.instr_obj.pix
        self.user   = self.instr_obj.user
        self.z      = self.instr_obj.z
        self.fitdata= self.instr_obj.fitdata
        self.fit    = self.instr_obj.fit
        self.rms    = self.instr_obj.rms
        self.best    = self.instr_obj.best

    #wavecal.__doc__ = Instrument.wavecal.__doc__

    def dofit(self):
        self.instr_obj.dofit()

    def features(self):
        """
          Print the tuple (pix,fit,user), The final fitting was done
          using the pixel coordinates and the user (wavelength) values
          lists.

          The pixel and wavelength values are available with the
          class members 'pix' and 'user'.
        """
        self.instr_obj.features()
    #features.__doc__ = Instrument.features.__doc__


    def plotfu(self):
         self.instr_obj.plotfu()
         pl.show()

    def plot_ref(self):
        """ Utility function to plot resulting peaks with
            the associated wavelength making a good
            visual aid to quickly check the mapping.
        """

        import pyfits as pf

        if not pl.isinteractive(): pl.ioff()


        cdelt = self.wcs['cdelt']
        crval = self.wcs['crval']
        crpix = self.wcs['crpix']
        pix_1 = np.min(self.xpeaks)
        pix_2 = np.max(self.xpeaks)
        half = max(crpix - pix_1, pix_2 - crpix)

        pix = self.pix
        user = self.user

        reference = self.reference

        fpath = os.path.dirname(spu.__file__)
        reffile = os.path.join(fpath,reference['reffile'])
        cuarf = pf.getdata(reffile)
        wmin = reference['wmin']
        wmax = reference['wmax']

        pl.clf()


        # Resample cuarf
        xsize = self.lpix.size
        #newf =  spu.resample(cuarf,(xsize,))

        newf = cuarf      # DO NOT RESAMPLE
        nxcuar = np.linspace(wmin,wmax,len(newf))

        # U  Plot reference lines
        subu = pl.subplot(211) 
        subu.plot(nxcuar, newf)       # Plot the reference array.

        # Align the Reference plot to the Arc plot
        xmin =  half*cdelt + crval
        xmax = -half*cdelt + crval
        subu.set_xlim(xmin, xmax)
        subu.set_xlabel('Angstroms')
       
        indx = np.searchsorted(nxcuar,(xmin,xmax))

        #p1=pix_1*cdelt + crval
        #p2=pix_2*cdelt + crval
        #indx=np.searchsorted(nxcuar,[p1,p2])

        p1 = min(indx)
        p2 = max(indx)
        ymax = newf[p1:p2].max()
        #ymax = np.max(newf[np.min(mm):np.max(mm)])
        subu.set_ylim(-ymax*.1, ymax+ymax*.2)
        basename = os.path.basename(reffile)
        subu.set_title('Reference: '+basename+' ('+self.ad.instrument()+')')

        self.subu=subu    # TEMPORARY


        lpix = self.lpix          # The arcs array

        # M Plot the arc array
        # 
        subm = pl.subplot(212)
        subm.plot(lpix)
        subm.set_xlim(crpix+half,crpix-half)
        fname = os.path.basename(self.filename)
        if self.extv > 1:
            fname = '%s[SCI,%d]'%(fname,self.extv)
        label = '(Pixels) '+fname+ ' '+self.ad.filter_name()
        subm.set_xlabel(label)

        # Plot marks
        szf = len(newf)
        cd = (wmax-wmin)/szf
        ydel =np.max(lpix)*.03
        for k,(p,u) in enumerate(zip(pix,user),1):
           ix = (u-wmin)/cd
           if ix >= szf: continue
           height = newf[ix]
           subu.annotate('%.d'%k, xy=(u,height), xycoords='data',xytext=(0, 20),
                 textcoords='offset points',rotation='vertical',
                 horizontalalignment='center',fontsize=10,
                 arrowprops=dict(arrowstyle="->"))
           subm.annotate('%.d'%k, xy=(p,lpix[p]),xycoords='data',xytext=(0, 25),
                 textcoords='offset points',rotation='vertical',
                 horizontalalignment='center',fontsize=10,
                 arrowprops=dict(arrowstyle="->"))
        pl.draw()
        return

    def __title(self):
        ad = self.ad 
        bname = os.path.basename(ad.filename)
        try:
           ccdsum = str(ad['SCI',1].get_key_value('ccdsum'))
        except:
           ccdsum = ''
        filtername = str(ad.filter_name())
        if ad.is_type('F2'):
            line = 'F2:  '+bname+' '+filtername
        else:
            line = bname+' '+filtername+'\n'+str(ad.grating()) +\
                   ' ('+ccdsum+') '+ str(ad.phu_get_key_value('maskname'))
        pl.title(line)


    def pix2ref(self, pixels, delta=None, linelist=None):
        """
          Given a list of pixel coordinates, find the closest reference
          entry in the linelist that is within 'delta' Angstroms 
          from the fitted value.

          *Input*
            pixels: 
                A list of pixel coordinates.
            delta: 
                The  maximum  difference for a match between the line coordinate
                derived from the dispersion function and  a  coordinate  in  the
                coordinate  list.   Positive values are in user coordinate units
                and  negative values are in units of pixels.  The user will need
                to increase the absolute value of match if wide slits (1  arcsec
                or  wider) are used with the detector unbinned or binned by 2 in
                the spectral direction.
            linelist = (None). 
                An instrument dependent list of reference wavelengths.
                It should be one entry per line.

          *Output*

            pixs_waves: 
                A list of tuples. The first element of each tuples
                if the input pixel coordinate and the second is the
                wavelength found that is closest to the fitted value
                within delta Angstroms.
                If no wavelengths are encountered it returs an empty list.

        """
        if not hasattr(self,'z'):
            self.wavecal()

        if delta == None:
            if self.match < 0:
                if not hasattr(self, 'lpix'):
                   npts = self.imdata.shape[1]    # The number of x points
                else:
                   npts = self.lpix.size
                fitdata = self.z(np.arange(npts))
                delta = abs(self.match*(fitdata[0]-fitdata[-1]))/(npts-1)
            else:
                delta = self.match
        else:
            delta = delta
        
        pixels = np.asarray(pixels)
        if pixels.size == 1:
            pixels = [pixels]

        cuar = self.cuar
        if linelist != None:
            cuar = np.loadtxt(linelist,usecols=(0,))

        pix_wave = []
        fit = self.z(pixels)
        count = 1
        for p,f in zip(pixels,fit):
            indx = np.argmin(abs(cuar-f))
            diff = abs(cuar[indx]-f)
            if diff > abs(delta): 
                # Was rejected. See if the neighboors are 
                # farther away from delta. 
                self.print_debug(p,cuar[indx],diff, '>delta')
                out = None
            else: 
                self.print_debug(p,cuar[indx],diff, 'good  ',count)
                count +=1
                pix_wave.append((p,cuar[indx]))
      
        return pix_wave
        
    def Dplot_spatial_effect(self):
        """
          UTIL:
          After linearizing, there is a small effect
          remaining that depends on the spatial 
          coordinate. The dispersion coordinate
          changes slightly as the spatial coordinate
          increases.
        """
        ny,nx = self.imdata.shape
        yy=range(0,ny,20)
        for k in range(0,nx,60):
            # fzz function gives the wavelength at pixel k
            # from row yy[0] to yy[ny/20]
            pl.plot(self.fzz(k,yy)-self.fzz(k,yy[0]))
        pl.show()

    def Dplot_surface(self):
        """      
          UTIL:
          Same as Dplot_spatial_effect but now
          in 3 dimension
        """
        from mpl_toolkits.mplot3d import axes3d

        ny,nx = self.imdata.shape
        Y,X = np.mgrid[0:ny:20,0:nx:20]
        Z = self.fzz(X,Y)- self.fzz(X,Y[0,:])
        fig = pl.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
        pl.show()
                
    def plot_arcs(self):
        """
          After fitting the arcs from the ARC spectrum plots
          the arcs.
        """
 
        pl.clf()
        #from numdisplay import display

        if not hasattr(self,'zpeaks'):
            self.fit_arcs()

        sz = self.imdata.shape
        yy = np.arange(0,sz[0],10)
        #for k,xoff in enumerate(self._zoffset):
        #yy = self.yy
        #yarr = np.asarray(xx,dtype=int)
        for zp in self.zpeaks:
            xvals = zp(yy)
            xarr = np.asarray(zp(yy),dtype=int)
            #self.imdata[yarr,xarr] = 25000
            pl.plot(xvals,yy)

        pl.show()
        return
        #display(self.imdata,frame=2,quiet=True)

    def plot_residual(self):
        """ 
          Plot the fit residuals: (w, w-z(pix))
          where 'w' is the wavelength and 'z'
          if the fit function evaluator.
        """

        pl.clf()
        pix = self.pix
        user = self.user
        y = user-self.z(pix)
        ymin = np.amin(y)
        ymax = np.amax(y)
        pl.plot(user, y,'x')
        umin = np.amin(user)
        umax = np.amax(user)
        pl.xlabel('wavelength')
        pl.ylabel('Residuals (Angstroms)')
        pl.axis([umin,umax,ymin,ymax])
        self.__title()
        pl.show()


    def plot_features(self,peaks=False, reverse=False, save=False):
        """
          Plot the features found along with their wavelengths.
          Notice that the horizontal axis is in pixel units.

          peaks: 
             If True plots dashes on top of the found peaks.
          reverse:
             If True reverses the peaks order.
          save:
             If True saves the plot on a PDF in the '/tmp' directory 
             using the filename as root name.
           
        """

        if not pl.isinteractive(): pl.ioff()

        pl.clf()
        
        lpix = self.lpix
        if reverse: lpix=lpix[::-1]
        pix = self.pix

        if peaks:
            # Plot the peaks, regardless if they have wlen.
            pl.plot(lpix)
            ys = lpix[np.int32(pix)]
            pl.vlines(pix,ys+20,ys+ys/10,color='r')
            pl.axis([0,len(lpix),0,np.max(lpix)])
        else:
            # Scale X axis to wave units.
            npix = lpix.size
            xw = self.z(np.arange(npix))
            #pl.plot(xw,self.lpix)
            pl.plot(lpix)
            # switch the axis
            #pl.xlim(xw[-1],xw[0])
            xy = pl.axis()
            pl.axis([npix,0,xy[2],xy[3]])

            user = self.user
     
            
            za = lpix[np.asarray(np.around(pix),dtype=int)]
            zoff = max(za)/10.   # 6 digits figure plus dot; move upwards
            for k,(w,z) in enumerate(zip(user,za)):
                if w < 1.0: continue
                # Todo. If xa[i] and xa[i+1] are too close then
                #       increase y in yxtest by ~60 points (size of str(y))
                #pl.annotate('%.2f'%w, xy=(w,z), xycoords='data',xytext=(0, 20),
                x = round(pix[k])
                pl.annotate('%.2f'%w, xy=(x,z+z/10), xycoords='data',xytext=(0, 20),
                     textcoords='offset points',rotation='vertical',
                     horizontalalignment='center',fontsize=10)

        pl.xlabel('Pixel')
        pl.ylabel('ADU')
        self.__title()
        pl.draw()
        if save:
            root = bname.split('.')[0]
            pl.savefig('/tmp/'+root+'.pdf')

        return

    def fit_arcs(self):
        """
          For each arc in the ARC spectrum fit an order 2 polynomial.
          The algorithm: Find the peaks at the middle line that 
          is the sum of 'nsum' rows. Traverse to the image top every 
          100 rows or so recording the (x,y) pairs for each arc at this
          position. The same from the the middle line to the bottom.
          Fit each arc's (x,y)s to a 2nd order polynomial and add to a
          list of 'zpeaks'. Finally do a 3-sigma rejection on the fits 
          that have bad coefficients.
        """

        image = self.imdata
        im_ny,im_nx = np.shape(image)
        ym = im_ny/2           # Starting row (middlerow)
        self.ym = ym

        nt = self.ntmax         # Get the ranked nt peaks
        nsum = self.nsum        # Number of rows to average

        
        # Start by looking at the nt most intense peaks in the nsum rows from
        # the middle row of the image. This array 'pcen' is the location 
        # to use for all stripes up and down from the middle row.

 
        # Trace each arc returning 'xx':(x_array,peak_number) and 
        # yy: list of row positions.
        # xx is along the dispersion axis, 0: vertical, 1: horizontal 
        xx,yy = self.trace_arc()
     
        npeaks = len(xx[2])
        # Now fit each arc (npeaks of them)
        fitorder = 2
        fitname = 'polynomial'
        zpeaks=[]
        sigmas=[]
        coeffs=[[] for _ in range(fitorder+1)]
        for k in range(npeaks):
            # Fit spectral line (peak) k
            zp = gfit.Gfit(yy, xx[:,k],fitname=fitname,order=fitorder)
            sig = np.std(xx[:,k] - zp(yy))
            sigmas.append(sig)
            zpeaks.append(zp)

            for i,coe in enumerate(zp.coeff):
                coeffs[i].append(coe)

        zpeaks = np.asarray(zpeaks)
        sigmas = np.asarray(sigmas)
        coeffs = np.asarray(coeffs)

        
        NSIG = 3
        # We have curves that are bad based on the coefficients' stats.
        # Reject those curves with coefficients outside a 3 sigma
        # rejection window.
        if fitname == 'polynomial':    # Check first 2 coeffs
            c0 = abs(coeffs[0]-np.median(coeffs[0]))<NSIG*np.std(coeffs[0])
            c1 = abs(coeffs[1]-np.median(coeffs[1]))<NSIG*np.std(coeffs[1])

        # Eliminate any fit function (zp) that is outside a 3 sigma
        # window from the other zp's sigma.

        mval = np.median(sigmas)
        norm = abs(sigmas-mval)
        stdev = np.std(norm)
        good = norm < NSIG*stdev

        good = ~(~good+~c0+~c1)

        # Select only the good arcs' fit.
        self.zpeaks = zpeaks[good]

        return

    def fit_image(self,fit3d_function="chebyshev",order=4):
        """
          Given that any point of in an arc has the same wavelength,
          we can assemble sets of (x,y,lambda) for each arc in the
          image and fit a 3D function f(x, y, wavelength) to it.

          *Input*
        
            fit3d_function = "chebyshev"
               The surface fitting interpolation function  to
               fit (x_array, y_array, z_array) via the LSQ 
               (LeastSquate) method. Other functions are
               'legendre' and 'cubic'

            order = 4
               The order of the fitting function to use. This
               value applies to both axis.

            self.zpeaks: 
               list of evaluator functions for  most of the arcs
               in the spectrum such that for a given arc number 'i'
               and a row zpeaks[i](row) gives the arc position in the
               dispersion direction.

            self.z: 
               For a given arc position in the dispersion direction
               z(x) gives the wavelegnth associated to x.

          The assumption is that for any row the x coordinates along the 
          dispersion direction corresponding to the arcs will have the 
          same wavelength.

          *Output*
    
            self.zz:
               Function evaluator. Returns a wavelength array
               for and (x_array,y_array) input.

          *Example*

          >>> wc.fit_image()
          >>> wc.zz(2000, 500)  # x: 2000, row: 500
          >>> # For a GMOS LongSlit spectra using the zpeaks
          >>> yy=range(800,4000,200)   # An array of rows
          >>> print wc.zz(zpeaks[10](yy), yy)
          >>> # Will return a list of wavelength (about the same value)
          >>  # for positions along the image arc number 10. 
          
        """
        if not hasattr(self, 'zpeaks'):
            self.fit_arcs()

        zpeaks = self.zpeaks
        npeaks = len(zpeaks)

        # -- FIT a surface (peaks, rows, waves)
        # Each point along a given arc has the same wavelength.
        # Use the ARC functions to create a sparce grid of (x,y) for the
        # rows. Use 10 row intervals.


        x=[];y=[];zlam=[]
        # Wavelength for each peak position:
        #    Get the arcs pixel values at the middle row using the arc
        #    fitting functions, then get the wavelength for each of these
        #    pixel using the 'z' function.
        if not hasattr(self,'z'):
            self.wavecal()

        # The array of wavelengths at the middle line for the arcs positions
        # given by the zpeaks evaluators.
        zvals = self.z([float(zp(self.ym)) for zp in zpeaks])

        ny,nx = self.imdata.shape
        ym = self.ym
        nsum = self.nsum
        # For each row record the triple (x,y,zlam)
        for row in range(ym-nsum,nsum,-nsum):
            # Get pixel position for each arch at row
            xpixs = [float(zp(row)) for zp in zpeaks]
            # Make a list of the same current row position
            y.append([float(row)]*npeaks)
            x.append(xpixs)
            # For each position along an arc, we have the same wavelength
            zlam.append(zvals)

        # Reverse the lists
        y = y[::-1]
        x = x[::-1]
        # From middle row to top
        for row in range(ym+nsum,ny-nsum,nsum):
            xpixs = [float(zp(row)) for zp in zpeaks]
            y.append([float(row)]*npeaks)
            x.append(xpixs)
            zlam.append(zvals)

        x,y,zlam = map(np.asarray,(x,y,zlam))
            
        # Fit a surface (peak,row, wavelength(row,peak))
        x = x.flatten()
        y = y.flatten()
        zlam= zlam.flatten()
        xlim = (x.min(),x.max()); ylim = (y.min(),y.max())
        zz = gfit.Fit3d (x, y, zlam, function=fit3d_function,order=order)
        self.zz = zz    
        zlim = (zlam.min(),zlam.max())
        self.zinv = gfit.Fit3d (zlam, y, x, xlim=zlim,function=fit3d_function,order=order)

        return

    def save_wavecal_fit(self):
        """
          Write the wavecal 3Dfit paremeters in output FITS file
          as a BINTABLE.
          Will have one FITS row per image band
        """
        import pyfits as pf
        from astrodata import AstroData

        # Check that we have a features list already [pix,fit,user]

        if not hasattr(self, 'zz'):
             self.fit_image()
        #    raise RuntimeError('No transformation available (pix-->wlen)')

        # Define columns of the output FITS BINTABLE

        zz = self.zz
        nelem = len(zz.coeff)
        nformat = '%dE'%nelem*nelem
        c1 = pf.Column (name='order', format='J1')
        c2 = pf.Column (name='coeff', format=nformat)
        c3 = pf.Column (name='xlim', format='2E')     # (xmin,xmax)
        c4 = pf.Column (name='ylim', format='2E')     # (ymin,ymax)
                        # Fit function name, 20 char max
        c5 = pf.Column (name='fitname', format='20A')
        c6 = pf.Column (name='fitmode', format='20A')
    
        rrr = 2   
        tabout = pf.new_table([c1,c2,c3,c4,c5,c6], nrows=rrr)

        rowout = tabout.data[0]
        rowout.setfield('order', zz.order)
        rowout.setfield('coeff', zz.coeff)
        rowout.setfield('xlim', zz.xlim)
        rowout.setfield('ylim', zz.ylim)
        rowout.setfield('fitname', zz.fitname)
        rowout.setfield('fitmode', 'zzfit')

        zinv = self.zinv
        rowout = tabout.data[1]
        rowout.setfield('order', zinv.order)
        rowout.setfield('coeff', zinv.coeff)
        rowout.setfield('xlim', zinv.xlim)
        rowout.setfield('ylim', zinv.ylim)
        rowout.setfield('fitname', zinv.fitname)
        rowout.setfield('fitmode', 'zzinvfit')

        tabad = AstroData(tabout)
        tabad.rename_ext("WAVECAL", 1)
        #self.ad.append(tabad)

        return  tabout


    def read_wavecal_table(self):
        """
          read WAVECAL extension from
          input AD.
        """

        data = self.ad['WAVECAL'].data

        ev_function={}
        for row in [0,1]:
            tb = data[row]
            order = tb.field('order')
            coeff = tb.field('coeff')
            xlim = tb.field('xlim')
            ylim = tb.field('ylim')
            fitname = tb.field('fitname')
            mode = tb.field('fitmode')
            if mode == 'zzfit':
               self.zz = gfit.Eval2(fitname,order,coeff,xlim,ylim)
            elif mode == 'zzinvfit':
               self.zinv = gfit.Eval2(fitname,order,coeff,xlim,ylim)
            else:  # Cannot happen
               raise ValueError('..... read_wavecal. Bad fitmode value in table '
                                ' WAVECAL')
         
        return 

    def fit2D_IFU(self):
        """
           TESTING:

           Trace a FLAT IFU image.
           /data2/gmos/ifu/mrgS20120827S0069.fits

           Extracting a fiber from the ARC is done with:

           pix = median(arcdata[int(zp[k](x))-2:int(zp[k](x))+2,x]) for x in xx]
           for k in range(len(zpeaks))
              and xx= range(len(nx))

           Note: Maybe we can interpolate int(zp[k](x))-2:int(zp[k](x))+2
                 to avoid the int().
        """

        image = self.imdata
        ny,nx = image.shape
        xm = nx/2             # Starting column for tracing

        nsum = 10
        cpix = np.mean(image[:,xm-nsum:xm+nsum],axis=1)

        upeaks,uflux,fw = spu.find_upeaks(
                cpix, separation=3, nmax=1000, cradius=3)

        print 'IFU DATA0:',len(upeaks)
        dd = np.diff(upeaks)    # WARNING. diff leaves the last index out
        # The slit peaks are separated by ~5 pixel in spatial direction
        g, = np.where((dd>4.5)&(dd<6.5))
        print 'IFU DATA1:',len(g)
        # (g[i],g[i+1]). Add the last index g[i+1] to each chain.
        # eg. [2,3,6,7] we would need to have [2,3,4,6,7,8]
        ss = [g[0]]
        for i in range(1,len(g)):
            if g[i-1]+1==g[i]:
                ss.append(g[i])
            else:
                ss.extend((g[i-1]+1, g[i]))
        ss.append(g[i]+1)
        g = np.asarray(ss)
        slit_peaks = upeaks[g]
        islit_peaks = np.asarray(slit_peaks,dtype=int)

        yy,xx = self.trace_ifu(image,slit_peaks,fw)

        npeaks = len(slit_peaks)
        # Now fit each slit  (spectrum)
        zpeaks=[]
        sigmas=[]
        for k in range(npeaks):
            # Fit spectral line (peak) k
            zp = gfit.Gfit(xx, yy[:,k],fitname='chebyshev',order=4,xlim=(xx[0],xx[-1]))
            sig = np.std(yy[:,k] - zp(xx))
            sigmas.append(sig)
            zpeaks.append(zp)

            #for i,coe in enumerate(zp.coeff):
            #    coeffs[i].append(coe)
            #if 'f' in self.debug:
            #    for coe in zp.coeff:
            #        print '%8.3f'%coe,      # Coefficients
            #    print 'sig:',k,std

        zpeaks = np.asarray(zpeaks)
        sigmas = np.asarray(sigmas)
        #coeffs = np.asarray(coeffs)

        """  COMMENT
        # We have curves that are bad based on the coefficients' stats.
        # Reject those curves with coefficients outside a 4 sigma rejection
        # window.
        if fitname == 'polynomial':    # Check first 2 coeffs
            c0 = abs(coeffs[0]-np.mean(coeffs[0]))<3*np.std(coeffs[0])
            c1 = abs(coeffs[1]-np.mean(coeffs[1]))<3*np.std(coeffs[1])
            print 'Bad COEFFs:',np.where(~c0),np.where(~c1)

        # Eliminate any fit function (zp) that is outside a 3 sigma
        # window from the other zp's sigma.

        NSIG = 3
        good = abs(np.std(sigmas)) < NSIG*sigmas
        print 'BAD curves:',np.where(~good),np.where(~good+~c0+~c1)
        good = ~(~good+~c0+~c1)
        print '# of goods:',sum(good),sum(c0)

        # Select only the good arcs' fit.
        zpeaks = zpeaks[good]

        # For each row select only the good peaks
        nyy=[]
        for k,ya in enumerate(xx):
            nyy.append(list(np.asarray(ya)[good]))
        xx = nyy
        """
        npeaks = len(zpeaks)
        self.zpeaks = zpeaks
        
        print 'ZPEK:',sig,len(zpeaks)


    def trace_ifu(self, image, peaks, widths):
        """
          TESTING:
          Trace sets of one pixel wide lines along
          the x-axis collapsing in x by an nsum amount
          of columns. The peak positions found are appended to
          the yy.
        """

        # The dispersion axis is horizontal (x-direction). 
        ny,nx = image.shape
        ycen = peaks
        xm = nx/2             # Starting column for tracing
        nsum = 50
        #nsum = self.nsum
        
        # FORM AN ARRAY of arc peaks for each selected row starting from the
        # middle row (ym) and going down every nsum rows. The same from ym
        # up to the last top row in the image every nsum rows.
        yy=[];xx=[]
        xx.append(xm)    # column number
        yy.append(ycen) # append ycen list of peaks position at y row to list yy

        # From middle to lowest row
        #print 'fit2D:',len(ycen),ycen,'width:',widths
        for x in range(xm-nsum,nsum,-nsum):
            cpix = np.mean(image[:,x-nsum:x],axis=1)
            # Notice that we use the peak array from the previous stripe as a
            # rough location for current peak locator.
            
            cen = spu.wrecenter(ycen, cpix, widths)
            ycen = cen.copy()
            xx.append(x)    # row number
            yy.append(ycen) # append ycen list of peaks position to array xx

        

        # Reverse the lists
        xx = xx[::-1]
        yy = yy[::-1]
        # From middle row to top
        ycen = peaks
        for x in range(xm+nsum,nx-nsum,nsum):
            cpix = np.mean(image[:,x:x+nsum],axis=1)

            cen = spu.wrecenter(ycen, cpix, widths)
            ycen = cen.copy()
            xx.append(x)    # row number
            yy.append(ycen) # append xcen list of peaks position to array xx


        yy,xx = map(np.asfarray,(yy,xx))
        return yy,xx

    def trace_arc(self):
        """
          Trace_arc produces an array of (x,y) pairs for each
          arc in the image to a maximum of self.ntmax arcs.
          It starts at the middle row of the image looking for the
          peaks in a collapsed (summed) strip of self.nsum rows in the
          direction of decreasing rows. 
          

          output
          ------
            xx: array[nrows, npeaks]. The x coordinate array for the
                npeaks arcs.
            yy: array[nrows]. The y-coordinate array for the nrows
                position in each arc.
        """
        image = self.imdata
        ym = self.ym
        nsum = self.nsum

        # Look for nt peaks at middle row ym
        pcen,cf,widths = spu.find_upeaks(np.sum(image[ym:ym+nsum,:],axis=0),
                           separation=self.minsep, nmax=self.ntmax)

        # The dispersion axis is horizontal (x-direction). 
        ny,nx = image.shape
        xcen = pcen
        
        nsum = ny/20   
        # FORM AN ARRAY of arc peaks for each selected row starting from the
        # middle row (ym) and going down every nsum rows. The same from ym
        # up to the last top row in the image every nsum rows.
        yy=[];xx=[]
        yy.append(ym)    # row number
        xx.append(xcen) # append xcen list of peaks position at y row to list xx

        # From middle to lowest row
        for y in range(ym-nsum,nsum,-nsum):
            mpix = np.mean(image[y-nsum:y,:],axis=0)
            # Notice that we use the peak array from the previous stripe as a
            # rough location for current peak locator.
            
            cen = spu.wrecenter(xcen, mpix,widths)
            xcen = cen.copy()
            yy.append(y)    # row number
            xx.append(xcen) # append xcen list of peaks position to array yy

        # Reverse the lists
        yy = yy[::-1]
        xx = xx[::-1]
        # From middle row to top
        xcen = pcen
        for y in range(ym+nsum,ny-nsum,nsum):
            mpix = np.mean(image[y:y+nsum,:],axis=0)

            cen = spu.wrecenter(xcen, mpix,widths)
            xcen = cen.copy()
            yy.append(y)    # row number
            xx.append(xcen) # append xcen list of peaks position to array yy

        xx,yy = map(np.asfarray,(xx,yy))
        return xx,yy

    def resample_image_to_linearCoords(self,image=None, fit3d_function="chebyshev",order=4):
        """
         Resample the image to linear wavelength co-ordinates by
         applying the surface function evaluator from the method 
         'fit_image'.         

         Using the inverse function evaluator for each row of the 
         input image:
         pixel = self.zinv(lambda, row)

         Generating a set of lambdas with a dispertion value
         (cdelt = (self.z(nx) - self.z(1))/nx) as 
         (lambdas = (ixx-crpix)*cdelt + crval), where 'ixx' is the
         array of indices (1..nx) along the dispersion axis.
         With the inverse function we obtain the pixel coordinates
         corresponding to each lambda value. Interpolating the
         input image values at each of these new pixel coordinates
         using spline interpolation we linearize the input image.
         
         *Input*
   
            image: (None) 
               If an ndarray is supplied use this one as the input
               image to linearize. It should have the same dispersion
               characteristics since it uses the ARC image dispersion
               and fitting functions.
               
            fit3d_function = "chebyshev"
               The function to fit the input image.  This is use 
               only if the method 'fit_image' has not been ran.
               Other values are: 'legendre' and 'cubic' for a 
               bi-cubic spline function.

            order = 4
               Order of the fitting function

         *Output*
         
            new_image: 
               Ndarray of the same shape as the input image.


        """  
        import ginterpolate as gi

        if not hasattr(self, 'zz'):
            # Fit a function sf(peaks, shift) for each stripe
            self.fit_image(fit3d_function=fit3d_function,order=order)
       
        data = self.imdata 
        if image != None:
            data = image 

        # The output arrays
        ny,nx = np.shape(data)
        outdtype = data.dtype
        linear_im = np.zeros((ny,nx),dtype=outdtype)
        lzz= np.zeros((ny,nx),dtype=np.float)

        # ----- Linearize ------
        
        # Dispersion. We use the already calculated mapping function z
        # such that lambda=z(pixels).  The dispersion direction is 
        # along the x-axis.

        # The z function is equivalent to zz(*,self.ym)
        ym = ny/2
        z = lambda x: self.zz(x,ym)

        cdelt = (z(nx) - z(1))/nx
        crpix = 1
        if cdelt < 0:
            crval = z(nx)
            cdelt = -cdelt
        else:
            crval = z(1)

        # The linear array in wavelength units for the dispersion axis 
        lams = (np.arange(nx)-crpix)*cdelt + crval
        lamf = lambda x: (x-crpix)*cdelt + crval

        xa = []; ya=[]; za=[]
        for y in range(ny):
            # Calculate the pixel position for each lambda
            # using the inverse function.
            pixels = self.zinv(lams,y)
            line = gi.ginterpolate(data[y,:],pixels,order=4)
            linear_im[y,:] = np.asarray(line,dtype=outdtype)    # Keep the input datatype
            #lzz[y,:] = gi.ginterpolate(lamf(pixels),pixels,order=4,mode='nearest')

        # For testing. Generate a 2D function to get wavelength for any
        # point (x,y) in the image.
        #y,x = np.mgrid[0:ny:20,0:nx:20]
        #z = lzz[y,x]
        #x,y,z=(x.flatten(),y.flatten(),z.flatten())
        #fzz = gfit.Fit3d (x, y, z, function=fit3d_function,order=order)
        #self.fzz = fzz

        # Calculate the linear WCS 
        self.linear_wcs((ny,nx))

        del lzz
        return linear_im

        
    def linear_wcs(self,image_shape=None):
        """
           Calculate values of crpix, crval and cdelt
           given a linearized middle line from the resample
           image.
           
           *Input*
           
             image_shape: (None) Tuple (ny,nx)
                If None, it uses the input image.shape tuple.

           *Output*
     
             lwcs: (dictionary) The lineariazed wcs values:
             {'crpix':1,'crval':z(1),'cdelt':(z(nx)-z(1))/nx}
               
        """
        import ginterpolate as gi

        if image_shape == None:
            ny,nx = self.imdata.shape
        elif type(image_shape) is tuple:
            ny,nx = image_shape
        else:
            raise ValueError('..... linear_wcs: image_shape is not a tuple.')


        # Check that we have the zinv function (wave to pix).
        if not hasattr(self, 'zinv'):
            self.fit_image()

        ym = ny/2
        z = lambda x: self.zz(x,ym)

        # FITS standard: 1-base
        cdelt = (z(nx) - z(1))/nx
        crpix = 1
        if cdelt < 0:
           crval = z(nx)
           cdelt = -cdelt
        else:
           crval = z(1)

        # Form a linear array in wavelength units 
        wavef = lambda x: (x-crpix)*cdelt + crval
        ww = np.linspace(wavef(0),wavef(nx-1),nx)

        # Turn the ww array into a pixel unit array using the zinv function;
        # where zinv(wave,pix) is calculated from the original arc image.

        pixels = self.zinv(ww,ym)

        # Generate an interpolating function to calculate wavelength from pixel
        # positions from the original middle line.
        lz = gi.ginterpolate(ww, pixels, order=4, mode='nearest')

        # Given this linearized array calculate the linear wcs values.
        lwcs={'crpix':1,'crval':lz[0],'cdelt':(lz[nx-1]-lz[0])/nx}

        # REDO WCS =======================

        # We want to recalculate the wcs by finding the peaks of the linearized
        # image's middle line.
        # 
        ny = self.nsum
        # Unlinearized middle line
        ulpix = np.median(self.imdata[ym-ny:ym+ny,:],axis=0)

        # Linearize
        linear_lpix = gi.ginterpolate(ulpix,pixels,order=4)

        # Find the peaks
        peaks,uflux,fw = spu.find_upeaks(linear_lpix,
                self.minsep, nmax=self.ntmax, cradius=self.cradius)

        crpix = lwcs['crpix']
        crval = lwcs['crval']
        cdelt = lwcs['cdelt']
        wavef = lambda x: (x-crpix)*cdelt + crval

        # The linelist
        cuar = self.cuar

        npix = []
        nuser = []
        nfit = []

        # Calculate the wavelength for each peak using the linear wcs
        fit = wavef(peaks)

        delta = self.match
        count = 1

        # Find the closest reference line in the cuar (linelist array)
        for p,f in zip(peaks,fit):
            indx = np.argmin(abs(cuar-f))
            diff = abs(cuar[indx]-f)
            if diff > abs(delta):
                continue
            else: 
                count +=1
                npix.append(p)
                nfit.append(f)
                nuser.append(cuar[indx])

        # From this array, calculate the best WCS
        lz=gfit.Gfit(npix,nuser,'polynomial',1)

        wmin = lz(0)[0]
        wmax = lz(nx-1)[0]
        self.lwcs={'crpix':1,'crval':wmin,'cdelt':abs((wmax-wmin)/nx)} 
        self.print_debug( 'matching lines percentage:',1.0*count/len(peaks))

        return 
      

    def im_astrodata(self, linearized_image, ad_in=None):
        """ 
          Generate an AstroData object from the linearized image array
          and the input image phu and hdu['SCI']

          *Input*
     
          linearized_image: ndarray
                 Ndarray with a linearized image data.
          ad_in: None.
                 if not None it uses this AD object's phu and science
                 header, otherwise it will be the minimum phu. 

          *Output*
       
          adout: Astrodata object.
                 The science header will have the linearized WCS values.
        
            
        """
        data = linearized_image    
        if not hasattr(self, 'lwcs'):
            self.linear_wcs(data.shape)
     
        ad = self.ad
        if ad_in != None:
            ad = ad_in

        phu = ad.phu
        adout = AstroData(phu=phu)    # New AD object with input phu
        header = ad['SCI',self.extv].header
        #header = ad['SCI',1].header
        axis = 1
        if self.dispaxis == 2:
           axis = 2
        header.update('crpix%d'%axis, self.lwcs['crpix'])
        header.update('crval%d'%axis, self.lwcs['crval']) 
        header.update('cd%d_%d'%(axis,axis), self.lwcs['cdelt']) 
        adout.append(AstroData(header=ad['SCI',1].header,data=data))
        adout.filename = ad.filename
        
        return adout

    def resample_image_asAstrodata(self):
        # Uses ad['SCI',self.extv]
        linear_image = self.resample_image_to_linearCoords()
        return self.im_astrodata(linear_image)

    def print_debug(self,*args):
        if self.debug:
            for arg in args:
                print arg,
            print

class Instrument(object):
    """
        Base class to incorporate the common methods and
        attributes belonging to all supported instruments.

        The instruments and modes are:
    
        GMOS-S and GMOS-N:
            LongSlit
            IFU
            MOS
        GNIRS:
            LongSlit
            XD
        NIFS:
            IFU
        F2:
            LongSLit
            MOS
        NIRI:
            LongSlit       

    """
    def __init__(self,ad,param,pval):


        """
          Init method for the Instrument base class

          *Input*
          
            ad:    Input AstroData object

            param: (dictionary). Input parameters which are
                   mostly None. The user can change the values,
                   taking precedence over the default ones.

            pval:  (dictionary). The paramerts actual values. The
                   names are the same as param's but the values
                   is a merge between the default and the user supplied
                   ones, where the latter takes precedence.
        """

        filtername = str(ad.filter_name())

        scrpix = 'crpix1'
        scrval = 'crval1'
        scdelt = 'cd1_1'
        self.dispaxis = ad.phu_get_key_value('DISPAXIS')
        if self.dispaxis == 2:
            scrpix = 'crpix2'
            scrval = 'crval2'
            scdelt = 'cd2_2'
        self.ad = ad
        self.filename = ad.filename

        for k,v in param.items():
            if v != None:
                pval[k] = v
 
        ext = pval['extv']
        #if ad.count_exts('SCI') < ext:
        #    raise ValueError("'extv' value:"+str(ext)+" is greater than the number"
        #                     " of 'SCI' extensions: "+str(ad.count_exts('SCI')))

        # Create members out of the input parameter names.
        self.imdata = ad['SCI',ext].data
        self.nsum   = pval['nsum']
        self.minsep = pval['minsep']
        self.ntmax  = pval['ntmax']
        self.match  = pval['match']
        self.cradius= pval['cradius']
        self.bins_step= pval['bins_step']
        self.nbins  = pval['nbins']
        self.best   = pval['best']
        self.clip   = pval['clip']
        self.fitname= pval['fitfunction']
        self.fitorder= pval['fitorder']
        self.fwidth = pval['fwidth']
        self.debug = pval['debug']
        self.bigger_linelist = pval['bigger_linelist']
        self.nfound = 6
        self.param  =  pval
        self.rms = None
      
        hdr = ad['SCI',ext].header
        self.wcs = {'crpix':hdr[scrpix], 'crval':hdr[scrval], 'cdelt':hdr[scdelt]}


    def wavecal(self):
        """ 


        """
        
        self.set_wavecal_data()

        # The current algorithm is to match triples of peak's coordinates
        # with a triple of reference lines.
     
        # Instantiate a MatchLines object
        ml = MatchLines(self.xpeaks, self.cuar, self.wcs['cdelt'], self.bins_step)

        self.ml = ml
        # self.wc
        ml.wc = self      # We want to use this instance with MatchLines object.
                          # Use Inheritance better??

        self.find_wavesolution()
        if not self.solution_found:
            raise ValueError("Solution not found for "+self.ad.filename)

    def set_wavecal_data(self):
        """ 
          Initialize data for wavecal()

          - setup the middle line self.lpix from the image
            as the mean of self.nsum rows.
          - Find upto self.ntmax peak coordinates of the arcs in self.pix
        """
        data = self.imdata 

        nsum = self.nsum

        # Use at least 3 rows above and 3 rows on top of the middle
        # and take the median of these.
        sz = data.shape
        if len(sz) == 2:
            ny = max(3,nsum/2)
            ym = sz[0]/2
            lpix = np.mean(data[ym-ny:ym+ny,:],axis=0)
        elif len(sz) == 1:
            lpix = data
        else:
            raise ValueError("Image dimensions are greater than 2.")

        # Now find the arc peaks pixel coordinates
        upeaks,uflux,fw = spu.find_upeaks( lpix, self.minsep, nmax=self.ntmax,
                     cradius=self.cradius)

        self.xpeaks = upeaks   # A synonim
        self.lpix = lpix

    def find_wavesolution(self):
        """
          Matching ratios algorithm implementation

          - Separate the input wavelength range into
            at most self.nbins subranges.

        """
        
        # Create an array of wavelength  from the initial
        # 'z' function (pix --> wavelength).
        self.fitdata = self.z(np.arange(self.lpix.size))

        cuar = self.cuar
        ml = self.ml

        # Create a new MatchLines member
        ml.fitdata = self.fitdata    # NEW
        pixelUser = 0
        best = 99       # A goodness of fit variable (initial value)

        # Get subranges [wa,wb]
        for ncount,(wa,wb) in enumerate(self.breakWW()):
            if ncount >50: 
                self.print_debug(' <<<!!!! Maximum number of ranges reached:',
                       ' solution might not be optimal.')
                raise ValueError('.....No solution found.')
                break

            # Gest pixel and reference indices to the ratio arrays 
            # corresponding to the subrange [wa,wb].

            for ip1,ip2,ia,ib in ml.subrange_coords(wa,wb):
                # If indices to the pixel ratio are the same, there
                # is no data to match, get the next range.
                if ip1==ip2:
                   continue

                self.print_debug(   '.'*45, ' now yield:(%d, %d)'%(cuar[ia],cuar[ib]))

                # Generate the pixel and reference ratios
                ml.peak_ref_ratios(ip1,ip2,ia,ib)
                pixelUser = ml.voting(ip1,ip2,ia,ib)

                if pixelUser == 99: continue
                ii,jj = np.where(ml.votes > self.ml.MIN_VOTES)

                self.print_debug( 'voting::',len(ii),len(jj))

                if len(ii)>5:
                    best = self.fit_pixel_user(pixelUser)
                    if abs(best) <= self.best:
                        self.print_debug( ' >>>>>>> GOOD SOLUTION....< ',self.best)
                        break
                    #if best == 55:
                    #   print '                   best... Not enough votes, go get somemore'

            if best > self.best:   # No solution was found.
                self.solution_found = False
            else:   # Fit again with a bigget line list if available.

                # Save data before final_fit
                pix_save=self.pix[:]
                fit_save=self.fit[:]
                user_save=self.user[:]
                fx_save = self.z.fx
                coeff_save = self.z.coeff

                big_linelist = self.bigger_linelist
                fpath = os.path.dirname(spu.__file__)
                if big_linelist != None:
                    self.print_debug( 'Final Fit.. >>>>Using reference list....',big_linelist)
                    big_linelist = os.path.join(fpath,big_linelist)
                    self.cuar = np.loadtxt(big_linelist,usecols=(0,))
                rms = np.sum((np.asarray(self.fit) - self.user)**2)
                #self.print_debug( 'Final rms before:',rms)
                #self.rms = rms
                self.clip = 4
                last_best = self.dofit()
                rms_final = np.sum((np.asarray(self.fit) - self.user)**2)

                self.solution_found = True

                # NOTE: maybe it need to be rms_final > rms 
                if rms_final > 1.5:
                    #revert to old
                   self.print_debug(">>>> Reverting to previous solution. RMS > 1.5")
                   self.pix = pix_save
                   self.fit = fit_save
                   self.user = user_save
                   self.z.fx = fx_save
                   self.z.coeff = coeff_save
                #else:
                #   self.rms = rms_final
                return
        self.print_debug('\n  *****NO GOOD SOLUTION FOUND. No more ranges to look for.',
                        ncount)

    def fit_peaks(self,xx,yy,fitname=None,order=None):
         if fitname == None:
             fitname = self.fitname
         if order == None:
             order = self.fitorder
         z = gfit.Gfit(xx,yy,fitname=fitname, order=order)
         self.z = z
         return z

    def dofit(self):
        """
          From a set of candidate identifications features 
          fit and evaluate a dispersion solution.
        """

       
        # Perform an iterative improvement of the fit in order
        # to find more reference lines.
        for k in range(3):

            if (k > 0):
                self.add_linelist()

            z = self.fit_peaks(self.pix,self.user)
            self.fit = list(z(self.pix))

            # Add the peaks
            self.add_peaks(self.xpeaks)

            z = self.fit_peaks(self.pix,self.user)
            self.fit = list(z(self.pix))
            self.sort_pu()
            self.clipData(self.clip)


        nfeatures = len(self.pix)
        nmin = 2
        nfound = self.nfound
       
        cdelt = self.wcs['cdelt']

        fit = z(self.pix)
        rms = np.sum((fit - self.user)**2)
        rms = np.sqrt (rms/ max (1, nfeatures-nmin)) / abs(cdelt)
        rms = rms / self.fwidth
  
        ncandidate = self.ncandidate
        nmatch1 = self.nmatch1 

        # Compute line list matching fraction.
        ncandidate = max (nfound, ncandidate)
        fmatch = 1. - 1.0*nmatch1/ncandidate

        nt = len(self.xpeaks)
        ftmatch = 1 - 1.0*self.ntmatch / min(nt, ncandidate)
        
        # NOTE:: Put this at the top 

        RMSG     = 0.1            # RMS goal (fwidths)
        WRMS     = 0.34           
        FMATCHG  = 0.2            # Matching goal (fraction unmatched)
        WFMATCH  = 0.33           # Line list non-matching fraction
        WFTMATCH = 0.33           # Target line non-matching fraction

        best = WRMS * rms / RMSG
        best += WFMATCH * fmatch / FMATCHG
        best += WFTMATCH * ftmatch / FMATCHG

        return best

    def dofit_small_sample(self):
        """
          From a set of candidate identifications features 
          evaluate a dispersion solution.
          This function in use for small number of peaks (~10)
          and relies on a wcs set (crpix,crval,cdelt) and a linelist
          of reference wavelength for those peaks.

          - If the number of lines fitted is too small you could
            increased 'match' to a large number (~15) for a coarse match
            assuming that the linelist is not very populated as to not
            pick a bad wavelength.

          - A new fit is calculated from the (peaks,users) and 'dofit()'
            function is called to performed a better fit.

        """
        # Start from the beginning. 
        wcs = self.wcs
        p2w = lambda x: (x-wcs['crpix'])*wcs['cdelt'] + wcs['crval']

        if self.match < 0:
                if not hasattr(self, 'lpix'):
                   npts = self.imdata.shape[1]    # The number of x points
                else:
                   npts = self.lpix.size
                fitdata = p2w(np.arange(npts))
                delta = abs(self.match*(fitdata[0]-fitdata[-1]))/(npts-1)
        else:
                delta = self.match

        peaks = self.xpeaks
        fit = p2w(np.asarray(peaks))

        delta = self.match*4    # Double
        cuar = self.cuar
        count = 1

        npix = []
        nuser = []
        # Find the closest reference line in the cuar (linelist array)
        self.print_debug('pix_coor   Lambda   Lambda-List_value')
        for p,f in zip(peaks,fit):
            indx = np.argmin(abs(cuar-f))
            diff = abs(cuar[indx]-f)
            if diff > abs(delta):
                continue
            else: 
                self.print_debug( p,cuar[indx],diff, 'good  ',count)
                count +=1
                npix.append(p)
                nuser.append(cuar[indx])

        # From this array, calculate the best WCS
        self.z =gfit.Gfit(npix,nuser,'polynomial',2)
        self.pix = npix
        self.user = nuser 
        if len(npix) < 5:
           self.fitname = 'polynomial'
           self.fitorder = 2 
        self.dofit() 

        #print '>>>>> F2 H filter. There are very few peaks:',len(self.xpeaks)
        #for p,u in zip(self.pix,self.user):
        #   print p,u
    

    def fit_pixel_user(self,pixWave):
        """
          From the votes array (peaks_array, reference_array)
          fit a parabola if the number of elements in this list is less
          MIN_VOTES else fit with the input fitfunction name and order.   

          pixWave: 
             List of tuples (peak,reference) candidates for the final list.
        
        """
        # Initialize lists
        self.pix = []
        self.fit = []
        self.user = []

        # This is data from votes
        pix = [x for x,y in pixWave]     # Peak positions that have an associated wavel.
        user = [y for x,y in pixWave]     # Corresponding wavelengths (from cuar)
        if len(pix) == 0: return 5.0     # Large value for best

        # Fit low order polynomial becuase we do not have many pairs in votes.
        if len(pix) < self.ml.MIN_VOTES:
            z = self.fit_peaks(pix,user,'polynomial',3)
        else:
            z = self.fit_peaks(pix,user)
        self.z = z 

        fit = list(z(pix))
        self.newfeatures(pix,fit,user)

        self.clipData(self.clip)
        
        best = self.dofit()    

        self.print_debug( '%s best.....%.2f'%(' '*55,best))

        return best

    def add_linelist(self):
        """
          Add lines from a linelist that are inside the range 
          of the current features;e.i: fit[0] < llist < fit[n]

        """
        llist = self.cuar

        ncandidate = 0
        nmatch1= 0

        npts = self.lpix.size    # The input line length containing peaks.

        # Wavelength space for the input arc.
        fitdata = self.z(range(npts))
        self.fitdata = fitdata
        
        if self.match < 0.:
            cd = (fitdata[0] - fitdata[-1])/(npts-1)
        else:
            cd = 1

        fitmin = min (fitdata[0],fitdata[-1])
        fitmax = max (fitdata[0],fitdata[-1])

        # The reference list in increasing order
        slist = np.sort(llist)    # Sort in increasing order

        # Find indexes of reference lines within fitmax - fitmin
        g, = np.where( (slist < fitmax) & (slist > fitmin))
        ng = g.size
        lpix = []
        lfit = []
        luser = []
       
        for nn,i in enumerate(g):

            refline = slist[i] 

            # Get the pixel position of the line 'refline' 

            ncandidate = ncandidate + 1
            # upix is a good aproximation of the peak position
            upix = self.fit_to_pix(refline)

            # Get a precise peak location now
            pix = spu.cen1d(upix, self.lpix, 0.0, 10.)
          
            # If no peak was encountered in self.lpix at upix, continue.
            if (not pix) or pix<0 : 
                 continue
                    
            # turn this position to wavelength unit 
            fit = self.z(pix)
            diff = abs((fit - refline) / cd)
 
            # If the value is far from the reference line, continue.
            if diff > abs(self.match):
                continue

            nmatch1 = nmatch1 + 1
  
            lpix.append(pix)
            lfit.append(fit)
            luser.append(refline)

        # Add these features to self.pix,self.user and self.fit
        self.newfeatures(lpix,lfit,luser)
        self.ncandidate = ncandidate
        self.nmatch1 = nmatch1

        return 


    def features(self):
        """
          Print the tuple (pix,fit,user). The final fitting was done
          using the pixel coordinates and the user (wavelength) values
          lists.

          The pixel and wavelength values are available with the
          class members 'pix' and 'user'.
        """
             #' 16 1067.053 6032.086 6032.127'
       
        print '\nNo  pix      fit      reference'
        for k,(p,f,u) in enumerate(zip(self.pix,self.fit,self.user)):
            print '%3d %8.3f %.3f %.3f'%(k+1,p,f,u)

    def newfeatures(self,lpix,lfit,luser):
        """
           Append the input lists.
           
           *Input*
   
             lpix:
                List of peak positions to be added to self.pix
             lfit:
                List of fit function evaluations of the peak positions
                to be added to self.fit.
             luser:
                List of reference lines to be added to self.user.
        """

        # If output lists are empty, just add the input lists.
        if self.pix == []:
            self.pix.extend(lpix)
            self.fit.extend(lfit)
            self.user.extend(luser)
            return

        for pix,fit,user in zip(lpix,lfit,luser):
            # See if this pix is already in the list (self.pix)
            indx = np.argmin(abs(self.pix-pix))
            if abs(self.pix[indx]-pix) < self.minsep:
               if abs(fit-user) < abs(self.fit[indx]-self.user[indx]):
                   # Replace
                   self.pix[indx] =  pix
                   self.user[indx] = user
                   self.fit[indx] =  fit
            else:
               # Is a new feature
               self.pix.append(pix)
               self.user.append(user)
               self.fit.append(fit)

    def sort_pu(self):
        order = np.argsort(np.asarray(self.pix))
        self.pix = list(np.asarray(self.pix)[order])
        self.fit = list(np.asarray(self.fit)[order])
        self.user = list(np.asarray(self.user)[order])


    def fit_to_pix(self,fitcoord):
        """
          Calculates the pixel position of the input array given
          its wavelength (fitcoord).

          It uses the FITDATA array; i.e.
          FITDATA = (pix - crpix)*cdelt + crval
          in the 1st approximation or when a fitting solution is
          available:  FITDATA = eval_function(pix)
        """
        
        fitdata = self.fitdata
        
        DXMIN = 0.01

        if fitdata[0] < fitdata[-1]:
            if (fitcoord <fitdata[0]) | (fitcoord > fitdata[-1]):
                return (None)

            i = np.argmin(abs(fitcoord-fitdata))
            if fitdata[i] == fitcoord:
                pixcoord = i
                return pixcoord

            pixcoord = i-.5
            dx = 1.
            while dx > DXMIN:
                dx = dx / 2
                if self.z(pixcoord) < fitcoord:
                    pixcoord = pixcoord + dx
                else:
                    pixcoord = pixcoord - dx

        else:
            fmin = min(fitdata[-1],fitdata[0])
            fmax = min(fitdata[-1],fitdata[0])
            #if (fitcoord < fitdata[-1]) | (fitcoord > fitdata[0]):
            if ( (fitcoord > fmax) & (fitcoord < fmin)):
                return (None)


            # get the last fitdata element with condition
            g, = np.where(fitcoord < fitdata)
            if len(g) == 0: return(None)
            i = g[-1]

            if fitdata[i] == fitcoord:
                pixcoord = i
                return pixcoord

            # find the transformation for .5 pixel +/-
            pixcoord = i-.5
            dx = 1.
            while dx > DXMIN:
                dx = dx / 2
                if self.z(pixcoord) < fitcoord:
                    pixcoord = pixcoord - dx
                else:
                    pixcoord = pixcoord + dx

        return pixcoord



    def clipData(self,clip):
        """ Rejects data according to 
            clip*stddev, where stddev is std(y-z(x)).
            If points are rejected, it will refit again.
            Return the new arrays x,y and the new function 'z'.

        """
        x = self.pix
        y = self.user
        z = self.z

        norm = np.abs(y - z(x))
        stdev = np.std(norm)
        condition =  norm < clip * stdev
        
        y = np.asarray(y)[condition]
        x = np.asarray(x)[condition]

        while norm.max() > clip * stdev:
            if len(y) < z.order + 1:
                self.print_debug('Too few points left to fit!:'+str(len(y)))
                return x,y,z
            z = self.fit_peaks(x,y)     # it update self.z also
            norm = np.abs(y - z(x))
            stdev = norm.std()
            condition =  norm < clip * stdev
            
            y = y[condition]
            x = x[condition]

        self.pix =  list(x)
        self.user = list(y)
        self.fit =  list(z(x))
         

    def breakWW(self):

        pix = self.xpeaks
        cuar = np.asarray(self.cuar)

        nbins = self.nbins
        nt = pix[-1]    # last peak pixel location

        crpix=self.wcs['crpix']
        crval=self.wcs['crval']
        cdelt=self.wcs['cdelt']

        dw = nt * cdelt             # A  (Angstrom)
        crsearch = abs (0.1 * dw)         # A
        cdsearch = abs (0.1 * cdelt)      # A/pix
        dwmax = (nt-1)*(abs (cdelt) + cdsearch) + 2*crsearch    # A/pix
        dwmin = (nt-1)*(abs (cdelt) - cdsearch)                 # A/pix
        dwmin = max (0.1, dwmin / dwmax)                        # Scalar

        self.print_debug( '...WCS::',crpix,crval,cdelt)

        crmin = crval - crsearch         # A
        crmax = crval + crsearch
        cdmin = abs(cdelt) - cdsearch
        cdmax = abs(cdelt) + cdsearch
        

        if cdelt > 0.0:
            wa = crmin + (cdelt + cdsearch)*(pix[0]-crpix)
            wb = crmax + (cdelt + cdsearch)*(pix[-1]-crpix)
            wa = crval + cdelt*(pix[0]-crpix)
            wb = crval + cdelt*(pix[-1] -crpix)
            cuar.sort()   # Make sure
        else:
            wa = crmin + (cdelt - cdsearch)*(pix[-1]-crpix)
            wb = crmax + (cdelt - cdsearch)*(pix[0]-crpix)
            if cuar[0] < cuar[1]:
               cuar=cuar[::-1]   # Reverse so that cuar[0] > cuar[1]
            #wa = crval + cdelt*(self.lpix.size-crpix)
            #wb = crval + cdelt*(1 -crpix)
            wa = crval + cdelt*(pix[-1]-crpix)
            wb = crval + cdelt*(pix[0] -crpix)

        # select the valid range

        # Do this for now

        iwa = np.argmin(abs(cuar-wa))
        iwb = np.argmin(abs(cuar-wb))
       
        if nbins == 1:
            a = min(wa,wb)
            b = max(wa,wb)
            if cdelt < 0:
               yield (b,a) 
            else:
               yield (a,b)
        
        wA = cuar[iwa] 
        wB = cuar[iwb] - wA

        self.print_debug( "wa_02 wb_02=",(wA,wB))
        self.print_debug( '.'*60,'nbins:',nbins)

        for i in range(2):
            for j in range(1,nbins+1):  
                if j ==1:
                   nbins = (nbins+2)/2
                elif j%2 == 0:
                   nbins = (nbins+2-j)/2
                else:
                   nbins = (nbins + 1 +j)/2 
                nbins = 2 * nbins - 1
               
                for k in range(1,nbins+1):
                    if k == 1:
                       bin = (nbins+2)/2
                    elif k%2 == 0:
                       bin = (nbins + 2 - k) / 2
                    else: 
                       bin = (nbins + 1 + k) / 2

                    if ((nbins-1)/2)%2 == 0:
                        if bin%2 == i:
                           continue
                    else:
                        if bin%2 != i:
                           continue
      
                    wb = wB / ((nbins + 1) / 2)
                    wa = wA + wb / 2 * (bin - 1)
                    wb = wa + wb

                    a = min(wa,wb)
                    b = max(wa,wb)
                    if cdelt < 0:
                       yield (b,a) 
                    else:
                       yield (a,b)

    def add_peaks(self, pix):
        """
          Match a feature in coordinate 'pix' againt a line list and
          return the closest match within a tolerance given by the
          attribute self.match.

        """

        pix = list(pix)

        if self.match < 0:
            if not hasattr(self, 'lpix'):
                npts = self.imdata.shape[1]    # The length of the spatial axis.
            else:
                npts = self.lpix.size
            fitdata = self.z(np.arange(npts))
            delta = abs(self.match*(fitdata[0]-fitdata[-1]))/(npts-1)
        else:
            delta = self.match

        cuar = self.cuar

        npix = []
        nuser = []
        nfit = []
        fit = self.z(pix)
        for p,f in zip(pix,fit):
            indx = np.argmin(abs(cuar-f))
            diff = abs(cuar[indx]-f)
            if diff > abs(delta): 
                # Was rejected. See if the neighboors are 
                # farther away from delta. 
                out = None
            else: 
                npix.append(p)
                nfit.append(f)
                nuser.append(cuar[indx])
      
        self.newfeatures(npix,nfit,nuser)
        self.ntmatch = len(nuser)
        
    def plotfu(self):

        from matplotlib import pyplot as pl
        import pyfits as pf

        # NOTE: plotfu is not complete when calling before using 'wavecal'
        # method
        if not hasattr(self, 'xpeaks'):
           xpeaks,uflux,fw = spu.find_upeaks(
            self.lpix, self.minsep, nmax=self.ntmax, cradius=self.cradius)

        pix = self.xpeaks

        fpath = os.path.dirname(spu.__file__)
        if self.ad.is_type('GMOS_SPECT'):
            wmin = 3053.    # x*0.25 + 3053. cuar.fits wavelength at cuar[0]
            wmax = 10423.   # cuar.fits wavelength at cuar[-1]
            cuarfile = 'cuar.fits'
        elif self.ad.is_type('F2_LS_ARC') or \
             self.ad.is_type('GNIRS_SPECT') or \
             self.ad.is_type('NIRI_SPECT'):
            wmin = 7035.14    # x*0.2.4321 + 7035.1381. cua.fits wavelength at cuar[0]
            wmax = 25344.54   # cuar.fits wavelength at cuar[-1]
            cuarfile = 'calibrated_arc.fits'
        elif self.ad.is_type('NIFS_SPECT'): 
            #cuarfile = 'calibrated_arc.fits'
            filtername = str(self.ad.filter_name())
            if 'JH' in filtername:
                # 1 14866.2099 1.600
                wmin = 14866.210
                wmax = 18129.449
                cuarfile = 'NIFS_H_arc_1D.fits'
            elif 'HK' in filtername:
                # 1 20140.26  2.13391
                wmin = 20138.12
                wmax = 24491.302
                cuarfile = 'NIFS_K_arc_1D.fits'
            elif 'ZJ' in filtername:
                # 1 9486.033 1.067928
                wmin = 7035.14
                wmax = 25344.54
                cuarfile = 'calibrated_arc.fits'

        cuarf = pf.getdata(os.path.join(fpath,cuarfile))

        # Set up the reference spectrum axis in wave units
        xcuar = np.linspace(wmin,wmax,len(cuarf))         # Get lambda scale

        # Find the indices of the 1st and last arc (xpeaks)
        indx1 = np.searchsorted(xcuar,self.z(pix[0]))
        indx2 = np.searchsorted(xcuar,self.z(pix[-1]))
        imin = min(indx1,indx2)
        imax = max(indx1,indx2)
        spec_max= np.max(self.lpix)
        ymax = cuarf[imin:imax].max()

        # Scale cuar to the same as lpix
        ncuar = cuarf*(spec_max/ymax)

        # Calculate new max,min reference value
        ymax = ncuar[imin:imax].max()
        ymin = ncuar[imin:imax].min()

        pl.clf()

        # 1st subplot.                  # PLOT reference linearized spectra  
        pref = pl.subplot(211)
        pref.plot(xcuar[imin:imax],ncuar[imin:imax])      

        nnx = np.sort(xcuar[imin:imax])    
        nny = np.sort(ncuar[imin:imax])
        ref = np.sort(self.cuar)
        
        if hasattr(self,'z'):
            # Find the indices within ref (self.cuar) array that 
            # encloses the peaks found.
            rindx = np.searchsorted(ref,[self.z(pix[0]),self.z(pix[-1])])

            # Get the indices of peaks in wave units
            refindx = np.searchsorted(nnx,ref[rindx[1]:rindx[0]])
            g = np.where(refindx<len(nny))
            height = nny[refindx[g]]
            pref.vlines(nnx[refindx[g]],ymin=height+20,ymax=height+height/10,color='r')
            xw = self.z(np.arange(self.lpix.size))
            pref.set_xlim((xw.min(),xw.max()))

        pref.set_title(cuarfile)
        pref.set_ylim((-2000,ymax))                    # SET YMIN to -2000


        # 2nd subplot
        sz = self.imdata.shape
        if len(sz) == 2:
            nrows = self.nsum
            ym = self.imdata.shape[0]/2
            lpix = np.median(self.imdata[ym-nrows:ym+nrows],axis=0)
        else:
            lpix = self.imdata
        pspec = pl.subplot(212)
        pspec.plot(lpix)                          # PLOT ARC
        pspec.set_xlim(len(lpix),1)

        dash = max(self.lpix)/10
        height = self.lpix[np.asarray(self.xpeaks,dtype=np.int)]
        pspec.vlines(self.xpeaks,ymin=height+20,ymax=height+dash,color='r')
        pl.show()

    def print_debug(self,*args):
        if self.debug:
            for arg in args:
                print arg,
            print

class NIRI(Instrument):
    def __init__(self,ad,param):
        pval = {}
        pval.update(INSTRUMENT['ALL'])
        pval.update(INSTRUMENT['NIRI'])
        self.instrument_mode = 'Long Slit'

        self.reference = {'reffile':'calibrated_arc.fits','wmin':7035.14,'wmax': 25344.54}
        Instrument.__init__(self,ad,param,pval)
        self.print_debug( "\n>>>>>>>>>>>>>>>>>>>>>>>>>> NIRI: %s %s"%\
                (os.path.basename(ad.filename),ad.filter_name()))

        self.param = pval

class F2(Instrument):
    def __init__(self,ad,param):

        pval = {}
        pval.update(INSTRUMENT['ALL'])
        filtername = ad.filter_name().as_str()
        if 'JH' in filtername:
            pval.update(INSTRUMENT[('F2','JH')])
        if 'J_' in filtername:
            pval.update(INSTRUMENT[('F2','J_')])
        elif 'HK' in filtername:
            pval.update(INSTRUMENT[('F2','HK')])
        elif 'H_' in filtername:
            pval.update(INSTRUMENT[('F2','H_')])
        elif 'Ks' in filtername:
            pval.update(INSTRUMENT[('F2','Ks')])
        self.reference = {'reffile':'calibrated_arc.fits',
                          'wmin':7035.14,'wmax': 25344.54}

        if 'pix-slit' in ad.phu.header['mospos']:
            mode = 'Long Slit'
        elif 'mos' in ad.phu.header['mospos']:
            mode = 'MOS'
        else:
            mode = None 
        self.instrument_mode = mode

        Instrument.__init__(self,ad,param,pval)
        self.print_debug ("\n>>>>>>>>>>>>>>>>>>>>>>>>>> F2: %s %s"
               %(os.path.basename(ad.filename),ad.filter_name()))
        imdata = self.imdata
        sz = np.shape(imdata)
        if len(sz) == 3:
           imdata = imdata.reshape(sz[1],sz[2])
        self.imdata = imdata.transpose()
        self.param = pval


    def wavecal(self):
        """ 
            F2 LongSlit wavelength solution.
        """
        self.print_debug ('=================================== F2 wavecal',
                          '*'*20, '  F2')
        self.set_wavecal_data()

        #============================================================

        # With the above information compute the wavelength solution.

        if 'H_' in self.ad.filter_name().as_str():
            z= self.z
            self.dofit_small_sample()
        else:

            # The current algorithm is to match triples of peak's coordinates
            # with a triple of reference lines.

            # Instantiate a MatchLines object
            ml = MatchLines(self.xpeaks, self.cuar, self.wcs['cdelt'], self.bins_step)

            self.ml = ml
            # self.wc
            #ml.wc = self      # We want to use this instance with MatchLines object.
                              # Use Inheritance better??

            self.find_wavesolution()


class NIFS(Instrument):
    def __init__(self,ad,param):
        pval = {}
        pval.update(INSTRUMENT['ALL'])
        filtername = ad.filter_name().as_str()
        if 'JH' in filtername:
            pval.update(INSTRUMENT[('NIFS','JH')])
            self.reference = {'reffile':'NIFS_H_arc_1D.fits',
                              'wmin':14866.210,'wmax': 18129.449}
        elif 'HK' in filtername:
            pval.update(INSTRUMENT[('NIFS','HK')])
            self.reference = {'reffile':'NIFS_K_arc_1D.fits',
                              'wmin':20138.12,'wmax': 24491.302}
        elif 'ZJ' in filtername:
            pval.update(INSTRUMENT[('NIFS','ZJ')])
            self.reference = {'reffile':'calibrated_arc.fits',
                              'wmin':7035.14,'wmax': 25344.54}

        self.instrument_mode = ad.phu.header['obsmode']
            
        Instrument.__init__(self,ad,param,pval)
        self.print_debug ( "\n>>>>>>>>>>>>>>>>>>>>>>>>>> NIFS: %s %s"%\
                (os.path.basename(ad.filename),ad.filter_name()))

        self.param = pval

    def nifs_wcal(self):
        """ 
          Wavecal a NIFS arc spectrum that has been already reduced
          and have more than one extension with slit data.

          1) Wavecal
          2) fit2D_spectrum
          3) linearize
             - calculate dispersion cd1_1
             - calculate crval for crpix=1 
          3) update WCS in header
        """
        pass

    def nifs_offsets(self):
        """
          From an ARC file find the offset between the
          different slices by looking into a bright line's peak
          position about the middle of the spectrum.
        """

        # Get the MDF data

        for ad in self.ad['SCI']:
            ny,nx = ad.data.shape
            ym = ny/2
            upeaks,uflux,fw = spu.find_upeaks(np.mean(ad.data[ym-5:ym+5,:],axis=0),
            self.minsep, nmax=5, cradius=self.cradius)

class GMOS(Instrument):
    def __init__(self,ad,param):

        long_slit = ad.is_type('GMOS_LS_ARC')
        ifu =       ad.is_type('GMOS_IFU')
        mos =       ad.is_type('GMOS_MOS')

        # Check that we have only one SCI extension
        if ad.count_exts('SCI') > 1:
            if not (long_slit or ifu or mos):
                raise ValueError("........ERROR. Only one SCI extension for GMOS LongSlit allowed.")
        pval = {}
        pval.update(INSTRUMENT['ALL'])
        if long_slit:
            pval.update(INSTRUMENT[('GMOS','ls')])
            mode = 'Long Slit'
        elif mos:
            pval.update(INSTRUMENT[('GMOS','mos')])
            mode = 'MOS'
        elif ifu:
            pval.update(INSTRUMENT[('GMOS','ifu')])
            mode = 'IFU'
        else:
            raise ValueError("ERROR. GMOS spectrum not LongSlit, MOS nor IFU.")
        self.instrument_mode = mode

        # x*0.25 + 3053. cuar.fits wavelength at cuar[0]
        self.reference = {'reffile':'cuar.fits',
                              'wmin':3053.,'wmax': 10423.}

        Instrument.__init__(self,ad,param,pval)
        self.print_debug( "\n>>>>>>>>>>>>>>>>>>>>>>>>>> Class  GMOS: %s %s"%\
            (os.path.basename(ad.filename),ad.filter_name()))
        self.param = pval

    def breakWW(self):
        self.print_debug( 'GMOS breakWW:')
        for a,b in  super(GMOS,self).breakWW():
            yield (a,b)

    def ifu_ident(self):
        """
          Reidentify other rows in the IFU image given a default
          solution for an initial row.

          Given a z(pix) function and an xpeaks array find the solution 
          for another with an xpeaks  
        """
        ny,nx = self.imdata.shape

        # Find a set of zpeaks(arcn,row) such that we get a peak_n for a given
        # row. Fit a function (new_peaks, user), where 'user' are the wavelengths
        # from the initial fit.
        

class GNIRS(Instrument):
    def __init__(self,ad,param):
        pval = {}
        pval.update(INSTRUMENT['ALL'])
        pval.update(INSTRUMENT['GNIRS'])

        mode = None
        if 'Long' in ad.decker().as_pytype():
            mode = 'Long Slit'
        elif 'XD' in ad.decker().as_pytype():
            mode = 'Cross Dispersed'
        self.instrument_mode = mode
        self.reference = {'reffile':'calibrated_arc.fits',
                              'wmin':7035.14,'wmax': 25344.54}
        Instrument.__init__(self,ad,param,pval)
        self.print_debug( "\n>>>>>>>>>>>>>>>>>>>>>>>>>> Class  GNIRS: %s %s"%\
            (os.path.basename(ad.filename),ad.filter_name()))
        self.imdata = self.imdata.transpose()
        self.param = pval
        

    def wavecal(self):
        """ GNIRS LongSlit wavelength solution.
        """
        self.set_wavecal_data()

        #============================================================


        # The current algorithm is to match triples of peak's coordinates
        # with a triple of reference lines.

        # Instantiate a MatchLines object
        ml = MatchLines(self.xpeaks, self.cuar, self.wcs['cdelt'], self.bins_step)

        self.ml = ml
        # self.wc
        ml.wc = self      # We want to use this instance with MatchLines object.
                          # Use Inheritance better??

        self.find_wavesolution()
        if not self.solution_found:
            # Try a straight association
            self.dofit_small_sample()

        
class MatchLines(object):
    """
    Pattern matching schema where the input data is the list of arc pixel
    positions (peaks) and a corresponding list of reference wavelengths (user).

    1) As the algorithm is of the order 6; (number_of_lines**6) we
    divide the list in pieces containing no more than 12 elements. Otherwise
    the number of possible combinations makes the algorithm too slow (from
    tenths of seconds to minutes).

    2) Form a list of triples as a combination of the indices with
    no more than 12 elements. These are indices pointing to the
    different elements in the peaks and user's list.

    3) Form a list a tuples where the first element is the difference of the
    last and first element in the triple and the second element in the tuple is
    the ratio of the difference between second and first element and the last
    and first element of the triple.

    4) Loop:  Subtract each ratio from the peak triple to all the
    ratios in the user's triple. If the difference is less than a given
    tolerance proceed to add all the peaks between the triple indices
    that meet the criteria where a line in the triple neighborhood
    is in the user's list. We put a weight on this position
    and add it to a voting matrix[peak_positions,user's_location].

    5) Extract from the voting matrix the locations containing more than -for
    example, 10 votes giving a list of (peak_positions, reference_wavelength)
    pairs.

    6) If we have at least 5 pairs we attempt to fit a polynomial and then 
    proceed to add more lines from the reference line list. To verify that a
    good solution has been found we calculate a figure of merit based on the 
    IRAF task autoidentify (please see 'help aidpar' on the IRAF cl)

    7) The the value of the figure of merit is larger than a minimum value
    we go back to 2) using a different section of the peaks list and we compute
    values to add to the voting matrix, which give a larger sample and
    a better probability to improve the figure of merit.

    Ref: This work is based on: FOCAS Automatic Catalog Matching Algorithms,
         Valdes, F. G., Campusano, L. E., Velasquez, J. D., & Stetson, P. B.
         Publications of the Astronomical Society of the Pacific, v.107, p.1119
         and on private communications with Mischa Schirmer.
        
    """

    def __init__(self, xpeaks, ref_list, cdelt, bins_step=12, min_votes=10):

        """     
           Instantiates a MatchLines object.

        *Input*
       
           xpeaks:
              Arc peak positions from an Arc image. Generally calculated
              in Wavecal.wavecal() method.

           ref_list: 
              A list of reference wavelengths covering the range of
              arcs in xpeaks.

           cdelt:
              A dispersion value in units of Angstrom/pixel
            
           bins_step: (12)
              The step in the peak indices list when the number of triples
              is greater than 800.

           min_votes: (10)
              The minimum number of votes to add a pair (peak,reference)
              as a candidate matching list. 

        *Attributes*
      
           xpeaks: Same as Input
           cuar:   Same as ref_list.
           cdelt:  Initial dispersion value. Same as Input
           bins_step:
           votes: Matrix to accumulate votes representing the matching
                   triples counting.         
           xpeaks_ratios: 
                  The pixel ratio list.
           ref_ratios:
                  The reference wavelengths ratio list.
           xpeaks_triples: 
                  The list of tuples containing peak coordinates
                  length and ratio.  Example [c-a, (b-a)/(c-a)]
           ref_triples: 
                  The list of tuples containing reference wavelength
                  length and ratio. 
                  
        *Methods*
       
        subrange_coords
        peak_ref_ratios

        """

        self.xpeaks = xpeaks           
        self.cuar = ref_list
        self.cdelt = cdelt
        self.bins_step = bins_step
        self.votes = np.zeros((len(xpeaks), len(ref_list)),dtype=np.int32)
        self.MIN_VOTES = min_votes

    def subrange_coords(self, wa, wb):
        """
          Divide the input wavelength range (wb-wa) in difference pieces
          such that each one of them has no more that 800 different triples.

          wa,wb: 
             Starting and ending wavelength of the current subrange.
        """

        pix = self.xpeaks

        cuar = self.cuar
        # In increasing order for both cases (cdelt > and < 0)
        #cuar = np.sort(cuar)
        cuarsz = len(cuar)-1
        
        # A function to determine the number of combination given n.
        ncomb = lambda n: n*(n-1)*(n-2)/6

        # Get the index of cuar that is closest to wa,wb
        iwa = np.argmin(abs(cuar-wa))
        iwb = np.argmin(abs(cuar-wb)) 
        iwlist = [(iwa,iwb)]    # The initial cuar indices enclosing a wavelength
                                # range to look for matches.
        nw = abs(iwb-iwa+1)     # The number of reference lines in this range
        nc = ncomb(nw)      # The number of ref lines triples (<250)
        self.print_debug( '....... in Subrange_coords:', 
           'nrefs:%d (%d,%d)[nc:%d] NK:%d'%(nw,wa,wb,nc,self.bins_step))
        if nc > 800:   # was 400
           bins_step = self.bins_step 
           iwlist=[]
           #print '+++++++++++++++++++SUR01:',nc,iwa,iwb,'NK:',bins_step
           if self.cdelt < 0:
               for k in range(iwb,iwa+bins_step,-bins_step):
                   # Maybe we break the wavelength range rather than using
                   # cuar's
                   #print '.......SUR02:',(k,k-bins_step),(cuar[k],cuar[k-bins_step])
                   iwlist.append((k,k-bins_step))
           else: 
               for k in range(iwa,iwb, bins_step):
                   # Maybe we break the wavelength range rather than using cuar's
                   #print '.......SUR03:',(k,k+bins_step)
                   iwlist.append((k,min(k+bins_step,cuarsz)))

        for ia,ib in iwlist:
            swa = cuar[ia]    # Reference line at index ia
            swb = cuar[ib]

            # What is the pixel position in the input line of swa and swb.
            #print 'SUBR00',wa,wb,(self.fitdata[0],self.fitdata[-1])
            p1 = self.wc.fit_to_pix(swa)  
            p2 = self.wc.fit_to_pix(swb)
            # If the value is outside de input line coordinates it returns None. 
            if p1==None or p2==None:
                continue
            ip1 = np.argmin(abs(pix-p1))   # Get the closest pix array index to p1
            ip2 = np.argmin(abs(pix-p2))   
            npc = ncomb(abs(ip1-ip2+1))      # The number of peaks between ip1,ip2. 
            nwc = ncomb(abs(ia-ib+1))
            
            self.print_debug( '.....ip,p,w: (%d[%d] %d[%d])%d (%d[%d] %d[%d])%d'%\
                    (p1,ip1,p2,ip2,npc,swa,ia,swb,ib,nwc))
            
            # Make sure they are in increasing order
            par = np.sort([ip1,ip2])
            ip1 = par[0]
            ip2 = par[1]
            nx=ip2-ip1+1

            # return the pixel and reference indices for this iteration.
            yield (ip1,ip2,ia,ib)

    def peak_ref_ratios(self,p1,p2,r1,r2):
        """
          Generate pixel and reference lines (wavelength) ratios.
          From any 3-points (a,b,c) in increasing order form tuples:

          [c-a, (b-a)/(c-a)]

          If n=(p2-p1+1) there are a n: n*(n-1)*(n-2)/6 posible tuples.

          *Input*
         
            p1,p2:
               Indices to the peak coordinates list.
            r1,r2:
               Indices to the reference lines list
     
        """

        decreasing = self.cuar[0] > self.cuar[1]
           
        # Form reference triples 
        if r1>r2:
           ltr=[(k,j,i) for i in range(r1,r2,-1) for j in range(i-1,r2-1,-1) 
                                              for k in range(j-1,r2-1,-1)]
        else:
           ltr=[(i,j,k) for i in range(r1,r2+1) for j in range(i+1,r2+1) 
                                              for k in range(j+1,r2+1)]

        # Form a list of tuples(normalized_length,ratio)
        # The ratio is  (k[1]-k[0])/(k[2]-k[0])

        r=self.cuar

        # Select triples with ratios between 0.1 and 10.
        minr = 0.1      # CONSTANT
        maxr = 1/minr
        
        # Form reference tuples list with the triples
        rl=[[abs(r[l]-r[j]),+(r[k]-r[j])/(r[l]-r[j]),] for (j,k,l) in ltr]
        for rr in rl:
            if rr[1] < minr or rr[1] > maxr:
                rr[1]=1000.      # CONSTANT


        # Form pixel triples with the list of 'pix'
        ltp=[(i,j,k) for i in range(p1,p2+1) for j in range(i+1,p2) 
                                                   for k in range(j+1,p2)]

        p = self.xpeaks
        # Form peak tuples with the triples.
        rp=[[abs(p[l]-p[j]),+(p[k]-p[j])/(p[l]-p[j]),] for (j,k,l) in ltp]
        for xx in rp:
            if xx[1] < minr or xx[1] > maxr:
                xx[1]=1000.     # CONSTANT

        self.xpeaks_ratios = rp    # List of peak ratios
        self.ref_ratios    = rl    # List of reference ratios
        self.xpeaks_triples= ltp   # List of peak tuples
        self.ref_triples   = ltr   # List of reference tuples

    def voting(self,ip1,ip2,iwa,iwb):
        """
          Subtract each ratio from the peak triple to all the
          ratios in the reference's triple. If the difference is less than a given
          tolerance proceed to add all the peaks between the triple indices
          that meet the criteria where a line in the triple neighborhood
          is in the reference's list. We put a weight on this position
          and add it to a voting matrix[peak_positions,ref's_location].

          *Input*
          
            ip1,ip2: 
              Indices to peaks list
            iwa,iwb: 
              Indices to the reference list     

          *Output*
    
            List: 
               List of tuples (peak,reference) with the
               largest votes greater than MIN_VOTES.
        """
        
        pix = self.xpeaks
        cuar = self.cuar

        rp = self.xpeaks_ratios
        # Do not continue if there are no ratios.
        if len(rp) == 0: return 99

        rl = self.ref_ratios
        ltp = self.xpeaks_triples
        ltr = self.ref_triples

        tol_match = 0.01    # matching distance in triples space

        # Get the reference subset list
        if iwa>iwb:
           ref = cuar[iwb:iwa+1]
        else:
           ref = cuar[iwa:iwb+1]

        cdelt = self.wc.wcs['cdelt']

        # Arc pix subset
        rpix = pix[ip1:ip2+1]
        wmax = max(ref)
        wmin = min(ref)

        # Pixel Coordinates subrange
        pix_window = max(rpix)-min(rpix)
        resolution = (wmax-wmin)/pix_window

        # Make a list of lines to check for matching
        # criteria.
        nlines,self.lines_to_check= self.get_nchecks()
  
        # For every pixel ratio tuple see if we have a reference
        # tuple match. 
        for ip,p in enumerate(rp):
            self.ip = ip
            for il,r in enumerate(rl):
                # If the difference of both ratios is less than a 
                # given tolerance, proceed to find matches inside
                # the two reference values
                d = p[1]-r[1]
                if abs(d) < tol_match:
                    # triples match, fit a parabola
                    px=[pix[k]  for k in ltp[ip]]
                    lw=[cuar[k] for k in ltr[il]]
                    z=gfit.Gfit(px,lw,'polynomial',order=2)
                      
                    # WHY 0.4
                    # If the different between the given value of the
                    # dispersion and this fit is off, take the next 
                    # reference tuple.
                    if abs(cdelt - z.coeff[1]) > 0.4:
                        continue
                    residual = z(px)-lw
     
                    if (abs(residual/z.coeff[1])>0.5).any() or (abs(residual)>3.).any():
                        rms = -1
                    else:
                        rms = np.sqrt(np.mean(residual**2))
                  
                    # If the rms is within this range, look for the lines
                    # around the current triple.
                    if rms > 0. and rms < 0.4:
                        nsuccess = self.check_neighbours(z.coeff,resolution,nlines[ip])
                        # If more than two 
                        if nsuccess >= 2:
                            self.votes[ltp[ip],ltr[il]] += 1

        pixelUser = self.pairsByVote()
        return pixelUser

    def pairsByVote(self):
        """
          Find the largest votes for each pair of (peak,reference)
          The votes matrix migth have duplicates for a peak or a
          reference index. From these find the largest.
        
          *Output*
         
            self.pairsByVote: 
               List of tuples (peak,reference) candidates for the
               final list.

          

        """
        pix = self.xpeaks
        cuar = self.cuar
        votes = self.votes

        ii,jj = np.where(votes>self.MIN_VOTES)   

        #print 'votes;;;;;;;;;;'
        #for i,j in zip(ii,jj):
        #    print (pix[i],cuar[j]),votes[i,j]

        px = [pix[i] for i in ii]
        #cu = [cuar[j] for j in jj]
        
        #get index of repeated elements in px  (peaks)
        listmaxs = []

        for peak,pind in sorted(list_duplicates(px)): 
            if len(pind)==1:    # No repeated elements
                listmaxs.append(pind[0])
            else:
                # find the largest vote for pind list of indices
                maxind = np.argmax(votes[ii[pind],jj[pind]])
                # Append the peak index containing the largest vote
                listmaxs.append(pind[maxind])

        # Now delete duplicates cuars with the low votes
        pix= [pix[ii[k]] for k in listmaxs]
        cuar= [cuar[jj[k]] for k in listmaxs]
        votes= [votes[ii[k],jj[k]] for k in listmaxs]
        
        #get index of repeated elements in cu (reference)
        listmaxs = []
        for peak,pind in sorted(list_duplicates(cuar)): 
            if len(pind)==1:    # No repeated elements
               listmaxs.append(pind[0])
            else:
                # find the largest vote for pind list of indices
                maxind = np.argmax(np.asarray(votes)[pind])
                # Append the peak index containing the largest vote
                listmaxs.append(pind[maxind])

        pixelUser = [(pix[k],cuar[k]) for k in listmaxs]
        #print 'votes::',type(pixelUser),pixelUser

        return pixelUser
      

    def get_nchecks(self):
        """
          Make a list of lines to check for matching
          criteria.
     
          *Output*
         
            nlines: 
               Number of lines 
            linesCheck:
               List of peak indices that are within the
               peaks triples peak's position.
        """

        ltp = self.xpeaks_triples
        nlines = []
        linesCheck = []
        pix = self.xpeaks
        nx = len(pix)
        ndat = 3       # We have a triple

        # Now look at each pixel ratio in the input tuple list.
        for indx in ltp:

            l = 0
            imin = min(indx)    # peaks indices
            imax = max(indx)

            # the two lines immediately outside the triples
            low  = imin - 1
            high = imax + 1

            if low  >= 0:  l+=1
            if high <= nx: l+=1

            # find lines inside the triples
            for i in range(low,high):
                u = 0
                for j in range(ndat):
                  if i == indx[j]:
                    u+=1
                    break
                if u==0: l+=1
            ncheck = l

            # Setup up an array to hold 
            lines_to_check = np.zeros(ncheck,dtype=np.int32)
            
            # fill the array containing the line (indices) we have to check
            # (i.e. map from pixels to lambda and check the residuals)
            l = 0
            if low >= 0:
                lines_to_check[l] = pix[low]
                l+=1
            for i in range(low,high):
                u = 0
                for j in range(ndat):
                  if i == indx[j]:
                    u+=1
                    break
                if u==0:
                    lines_to_check[l] = pix[i]
                    l+=1
            
            if high <= nx: 
                lines_to_check[l] = pix[high]

            nlines.append(min(nx,ncheck))
            linesCheck.append(lines_to_check)
            
        return nlines,linesCheck


    def check_neighbours(self,coeff,resolution,nlines):

        nsuccess = 0

        #  evaluate the residuals (projecting the observed arcs to lambda,
        #  looking for the closest match)

        lines_to_check = self.lines_to_check[self.ip]

        for i in range(nlines):
            px = lines_to_check[i]
            lambda_predict = coeff[2] + coeff[1]*px + coeff[0]*px**2
            residuals = min(abs(lambda_predict - self.cuar))
            if residuals < resolution: nsuccess+=1

        return nsuccess
   
    def print_debug(self,*args):
        if self.wc.debug:
            for arg in args:
                print arg,
            print

from collections import defaultdict

def list_duplicates(seq):
        tally = defaultdict(list)
        for i,item in enumerate(seq):
            tally[item].append(i)
        return ((key,locs) for key,locs in tally.items())


