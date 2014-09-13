import os
import numpy as np
from matplotlib import pyplot as pl

from astrodata import AstroData

import wcal

import pyfits as pf
import spec_utils as spu
import gfit
from resample import resample

class Winter(object):

    """
      (In develop version)
         TODO if event.key == 'd':
            # Delete point in fig2 (Residual plot)
            # Note not working at the moment. Need to get events from
            # figure 2.... don't know how to do it yet.

      Interactive Wavelength Calibration program. 

      The program displays a reference ARC spectrum in wavelength units
      with the same scale as the input ARC spectrum's. This latter spectrum
      should have at least a rough aproximation of a WCS (crpix,crval and cdelt)
      and is displayed as a sesond subplot.

      Looking at the 2 plots the user can visually matched the arc lines, and
      easily make a selection of the correspoding lines in the reference and 
      input spectrum. To associate, click a 'u' at the reference line peak 
      and 'm' at the line peak in the input spectrum.

      
      The Winter available keystrokes are: (mouse pointer should be
                                            on the plot area)
    
      'h' : A list of the keystrokes available.
      'u' : Mark a reference line. Mouse pointer on top of the  
            peak in the top subplot. If an "ERROR: Reference line
            not found. Try again." appears you should move the 
            pointer closer to the peak. If the error accurs again
            then there is no reference line and you could choose
            another line.
      'm' : Mark the 'u' associated arc peak in the input spectrum
            on the bottom subplot. 
      'l' : List the (pix, wavelength) accumulated so far.
      'f' : Fit the above list using a second order polynomial if 
            the number of points is 3 or a 4th order Chebyshev polynomial
            if larger.
      'r' : Plot residuals in a new Figure. (Don't know how to delete yet)
      'd' : Delete poins in the residual plot (NOT yet working)
      'q' : To quit the loop.

       The fitting function with the mapping is available in the 'z'
       fit class object member. Submembers are:

       z(pix):    Returns the wavelength for pixel coordinates 'pix'.
       z.coeff:   List of fit function coefficients.
       z.fitname: Name of the fitting function used.
       z.order:   Fitting function order. 
        
       Linelist available for these instruments:
        
         F2:   argon.dat
         GMOS: cuar.dat
         NIRI: lowresargon.dat
         NIFS: argon.dat
          

      Example

      >>> import winter
      >>> reload(winter);wi=winter.Winter(fname)
      >>> wi.interactive()
      >>> # A plot appears with 2 subplots and the user is in a loop
      >>> # until a 'q' is hit.

    """
    def __init__(self,fname,linelist=None,extv=1,nsum=10,match=6,radius=10):

        self.extv = extv
        self.nsum = nsum
        self.match = match
        self.radius = radius

        ad = AstroData(fname)

        self.imdata = ad['SCI',extv].data

        if ad.dispersion_axis() == 2:
            # Transpose to work on the X-axis.
            self.imdata = self.imdata.transpose()

        self.filename = ad.filename
        self.ad = ad

        linelist_file = linelist

        scrpix = 'crpix1'
        scrval = 'crval1'
        scdelt = 'cd1_1'
        fpath=os.path.dirname(spu.__file__)
        if ad.is_type('F2_SPECT') or self.ad.is_type('GNIRS'):
            reffile = 'calibrated_arc.fits'
            wmin = 7035.14    # x*0.24321 + 7035.1381. wavelength at [0]
            wmax = 25344.54   # wavelength at [-1]
            scrpix = 'crpix2'
            scrval = 'crval2'
            scdelt = 'cd2_2'
            sz = np.shape(self.imdata)
            if len(sz) == 3:
               self.imdata = self.imdata.reshape(sz[1],sz[2])

            if linelist_file==None:
                linelist_file = 'argon.dat'
                #linelist_file = 'lowresargon.dat'
        elif ad.is_type('GMOS_SPECT'):
            if linelist_file==None:
                 linelist_file = 'cuar.dat'
            reffile = 'cuar.fits'
            wmin = 3053.    # x*0.25 + 3053. cuar.fits wavelength at cuarf[0]
            wmax = 10423.   # cuar.fits wavelength at cuarf[-1]
        elif ad.is_type('NIRI_SPECT'):
            reffile = 'calibrated_arc.fits'
            wmin = 7032.706
            wmax = 25342.113
            if linelist_file==None:
                 linelist_file = 'lowresargon.dat'
        elif ad.is_type('NIFS_SPECT'):
            filtername = str(self.ad.filter_name())
            if 'JH' in filtername:
                 if linelist_file==None:
                    linelist_file = 'argon.dat'
                 wmin = 7032.706
                 wmax = 25342.113
                 reffile = 'calibrated_arc.fits'
            elif 'HK' in filtername:
                 if linelist_file==None:
                    linelist_file = 'ArXe_K.dat'
                 wmin = 20140.259
                 wmax = 24491.320
                 reffile = 'NIFS_K_arc_1D.fits'
            elif 'ZJ' in filtername:
                 if linelist_file==None:
                    linelist_file = 'argon.dat'
                 wmin = 7035.14
                 wmax = 25344.54
                 reffile = 'calibrated_arc.fits'
        else:
            raise ValueError('Input file is not supported: '+ad.filename)

        self.reffile = os.path.join(fpath,reffile)
        self.wmin = wmin
        self.wmax = wmax
        self.pixs=[]
        self.users=[]

        linelist_file = os.path.join(fpath,linelist_file)
        linelist = np.loadtxt(linelist_file,usecols=(0,))

        hdr = ad['SCI',extv].header
        self.wcs={'crpix':hdr[scrpix], 'crval':hdr[scrval], 'cdelt':hdr[scdelt]}
        cdelt = self.wcs['cdelt']

        if cdelt < 0:
            linelist = linelist[::-1]    # Change to decreasing order

        self.cuar = linelist     # cuar is the IRAF name
        self.ad = ad

        sz = np.shape(self.imdata)

        nsum = 6
        if len(sz) == 2:
            ny = max(3,nsum/2)
            ym = sz[0]/2
            lpix = np.mean(self.imdata[ym-ny:ym+ny],axis=0)
        elif len(sz) == 1:
            lpix = self.imdata
        else:
            raise ValueError("Image dimesions are greater than 2.")

        self.lpix = lpix

        # Now find the arc peaks pixel coordinates
        upeaks,uflux,fw = spu.find_upeaks(self.lpix, 2, nmax=50, cradius=12)

        self.xfluxs = uflux
        self.xpeaks = upeaks   # A synonim

    def interactive(self):
        """
           Start interactive mode.
        
           - See help in the class for key mapping

        """
           

        global ncuarf, nxcuar

        extv = self.extv
        peaks = self.pixs

        pl.clf()

        # Plot the reference lists.  (in a Fits file).

        # Arc file WCS
        cdelt = self.wcs['cdelt']
        crval = self.wcs['crval']
        crpix = self.wcs['crpix']

        # Reffile min, max wavelength
        wmin = self.wmin
        wmax = self.wmax


        # 1st subplot ---------------------      REFFILE
        # Resample cuarf to the ARC size
        xsize = self.lpix.size

        resmp = False
        cuarf = pf.getdata(self.reffile)

        if resmp:
            newf =  resample(cuarf,(xsize,))
        else:
            newf = cuarf

        # Scale cuarf to the same as lpix
        spec_max = np.max(self.lpix)
        ymax     = newf.max()
        ncuarf = newf*(spec_max/ymax)
        self.uu=ncuarf

        print 'WCS:',self.wcs, 'size of reff: ',ncuarf.size

        halfw = max(crval-wmin,wmax-crval)
        #print 'Reff:',crval,halfw,(wmin,wmax),crval-halfw,crval+halfw

        # x-axis values in wavelength units
        nxcuar = np.linspace(wmin+1,wmax+1,len(newf))   

        pref = pl.subplot(211)
        pref.plot(nxcuar,ncuarf)                       
        basename = os.path.basename(self.reffile)
        pref.set_title('Reference: '+basename+' ('+self.ad.instrument()+')')
        self.subu = pref


        # 2nd subplot ----------------------     ARC
        bname = os.path.basename(self.filename)
        if extv>1:
            bname+='[SCI,%d]'%extv
        xlabel = bname+' middle row'

        halfp = max(crpix-0,len(self.lpix)-crpix)

        #print 'Half crpix:',halfp,(crpix+halfp,max(0,crpix-halfp))

        pspec = pl.subplot(212)
        self.subm = pspec
        pspec.plot(self.lpix)                   
        pspec.set_xlim(len(self.lpix)-1,0)
        #pspec.set_xlim(crpix+halfp,max(0,crpix-halfp))
        pspec.set_xlabel(xlabel)

        # Align the Reference plot to the Arc plot
        nwmin =  (0-crpix)*cdelt + crval
        nwmax = (len(self.lpix)-crpix)*cdelt + crval
        #nwmin =  halfp*cdelt + crval
        #nwmax = -halfp*cdelt + crval

        pref.set_xlim(nwmax,nwmin)
        indx = np.searchsorted(nxcuar,(nwmin,nwmax))
        if cdelt > 0:
            indx = indx -1
        
        p1 = min(indx)
        p2 = min(max(indx),len(nxcuar)-1)
        
        #print '.....Align 1st plot:',halfp,(nxcuar[p1],nxcuar[p2])
        #pref.set_xlim(nxcuar[p1],nxcuar[p2])
        #if cdelt < 0:
        #    pref.set_xlim(nwmax,nwmin)
        #else:
        #pref.set_xlim(nwmin,nwmax)

        ymax = ncuarf[p1:p2].max()
        pref.set_ylim(-ymax*.1,ymax+ymax*.2)
        
        # Redefine default keymapping 'l'
        pl.rcParams['keymap.yscale'] = 'i' 

        self.cid = pl.connect('key_press_event',self.key_events)
        pl.draw()



    def key_events(self, event):
        """
          Key Event loop 

          - hit 'q' to quit the loop
        """
        global ncuarf, nxcuar,uu


        self.event = event
        rd = self.radius

        if event.key == 'u':
            # Get a reference point 
            while True:
                try:
                    xm = event.xdata

                    #print 'UUU:',xm,min(nxcuar),max(nxcuar),len(nxcuar)
                    ix = np.searchsorted(nxcuar,xm)
                    y = ncuarf[ix-rd:ix+rd]
                    x = np.arange(ix-rd,ix+rd)
                    try:
                       peakidx = np.argmax(y)
                    except:
                       print 'key needs to be "m" in this subplot.'
                       break
                    center = x[peakidx]
                    width = abs(x[-1]-x[1])/2.
                    height = y[peakidx]
                    #print 'WW',width,len(nxcuar),xm,ix,len(ncuarf)
                    cen = spu.recenter([center],nxcuar,width=5)
                    #print 'xcc:',xm,peakidx,cen,height,nxcuar[cen[0]]
                    # find the reference value from the reference list
                    uu = self.aid_match (nxcuar[cen[0]], self.match)
                    print 'Coord:',xm,', "Line list entry":',uu
                    if not uu: 
                        print '\a \n*** ERROR: Reference line not found. Try again.'
                        break

                    ax = self.subu
                    ax.hold(True)
                    h1 = height+20
                    h2 = h1 + height/10
                    pl.vlines(uu,ymin=h1,ymax=h2,color='r')
                    #ax.vlines(uu,ymin=h1,ymax=h2,color='r')
                    pl.show() 

                    break
                except IndexError:
                    print 'Indx'
                    break
                except TypeError:
                    print 'TYerr'
                    break
                else:
                    break

        if event.key == 'm':
            while True:
                try:
                    
                    xm = event.xdata
                    nx = self.nsum
                    y = self.lpix[xm-20:xm+20]
                    x = np.arange(xm-20,xm+20)
                    try:
                       peakidx = np.argmax(y)
                    except:
                       print 'Error:: key needs to be "u" in this subplot.'
                       break
                    #print 'peakidx',peakidx
                    center = x[peakidx]
                    height = y[peakidx]
                    width = abs(x[-1]-x[1])/2.
                    cen = spu.recenter([center],self.lpix,width=width)
                    if uu<0:
                        print 'Error: Need to mark a ref line first in the reference subplot.'
                        break

                    print 'Center pix::',cen,' Reference:',uu

                    ax = self.subm
                    ax.hold(True)
                    h1 = height+20
                    h2 = h1 + height/10
                    ax.vlines(cen,ymin=h1 ,ymax=h2 ,color='r')
                    #pl.vlines(cen,ymin=height+20,ymax=height+20+height/10,color='r')
                    self.pixs.append(cen[0])
                    self.users.append(uu)
                    #print self.users
                    uu = -99
                    pl.show()
                   
                    break
                except IndexError:
                    print 'Indx'
                    break
                except TypeError:
                    print 'TYerr'
                    break
                else:
                    break

        if event.key == 'r':
            # Residual plot
            pixs = self.pixs
            diff = self.z(pixs)-self.users
            fig2 = pl.figure()
            ax_single = fig2.add_subplot(111)
            line2 = ax_single.plot(pixs,diff,'+')
            fig2.show()
            #pl.plot(pixs,diff)
            # Delete pair in (self.pixs,self.users)

        if event.key == 'd':
            # Delete point in fig2 (Residual plot)
            # Note not working at the moment. Need to get events from
            # figure 2.... don;t know yet how to.
            x = event.xdata
            indx = np.searchsorted(self.pixs,x)
            print indx,self.pixs[indx]

        if event.key == 'f':
            z= self.fit_peaks(self.pixs,self.users)
            if z == None:
                return

            for k,(p,w) in enumerate(zip(self.pixs,self.users)):
                print k,p,z(p),w,'%.2f'%(z(p)-w)

            # Fit the peaks find matches with cuar, refit
            xpeaks = self.xpeaks
            fit = z(xpeaks)
           
            # Find users
            pixs,users = self.fit2user(xpeaks,fit, self.match)

            # Fit again. Use a Chebyshev polynomial since we should have 
            # many (>5) lines.
            z= self.fit_peaks(pixs,users,fitname='chebyshev',order=4)
            print '\n ....New fit.....'
            for k,(p,w) in enumerate(zip(pixs,users)):
                print k,p,z(p),w,'%.2f'%(z(p)-w)

            pixs,users,z = self.clipData(pixs,users,z,clip=4)
            #indx = np.searchsorted(xpeaks,pixs)
            print ' >>>>>>NUMBER OF points after clipping:',len(pixs)
            #hh = self.xfluxs[indx]
            ipix = np.asarray(np.around(pixs),dtype=int)
            hh = self.lpix[ipix+1]

            ax = self.subm
            ax.hold(True)
            ax.vlines(pixs,ymin=hh+20,ymax=hh+20+hh/7.,color='r')

            self.pixs = list(pixs)
            self.users = list(users)
            self.z = z
            pl.draw()

        if event.key == 'h':
            print '\n Winter. Keystrokes available:'
            print '   u: Mark a line in the upper subplot.'
            print '   m: Mark the "u" associated line in the pixel (lower) subplot.' 
            print '   l: List tuples (pix, wavelength)  of lines already marked.'
            print '   r: Plot residuals in a new figure.' 
            print '   f: Fit the (pix,wavelength) tuples.'
            print '   q: Quit'
           
        if event.key == 'a':
            """
               Not working yet
            """
            # Add reference line list
            self.fitdata = self.z(range(self.lpix.size))
            #pixs,users = wcal.add_linelist(self.ll)
            pixs = list(self.pixs)
            user = list(self.users)
            pix.extend(pixs)
            user.extend(users)
            order = np.argsort(pix)
            pix = list(np.asarray(pix)[order])
            user = list(np.asarray(user)[order])
            fit = list(self.z(pix))
            #pix,fit,user = wc.weedoutDups(pix,fit,user)
            self.pixs = pix
            self.users = user
            fit = np.asarray(fit)
            user = np.asarray(user)

            # The new scale is.
            szf = len(ncuarf)
            cd = (10423.-3053.)/szf
            yy = self.lpix
            for p,u in zip(pix,user):
               ix = (u-3053)/cd
               if ix >= szf: continue
               height = ncuarf[ix]
               print p,u
               self.subu.vlines(u,ymin=height,ymax=height+2600,color='r')
               self.subm.vlines(p,ymin=yy[p],ymax=yy[p]+4000,color='r')
            #pl.draw()

            

        if event.key == 'q':
            pl.disconnect(self.cid)
            pl.close('all')
 
        if event.key == 'l':
            # List the current pairs
            k = 1
            for p,w in zip(self.pixs,self.users):
                print k,p,w
                k +=1
               

    def init_linearize(self):
        #if hasattr(self,'nfeatures'):
        #    self.wavecal()

        npix = len(self.pix)
        pix =  self.xpeaks
        user = self.users

        # We need to recalculate z to use the whole x-axis range.
        # linearize 1st
        z = gfit.Gfit(pix, user)
        self.z = z

        lambda1 = z(1)
        lambda2 = z(npix)
        lmin = min(lambda1,lambda2)
        lmax = max(lambda1, lambda2)
        cd11 = (lmax-lmin+1)/npix 

        xmin = min(user[0],user[-1])
        xmax = max(user[0],user[-1])
        # We calculate the inverse z. The xlims are in increasing order
        iz = gfit.Gfit(user, pix, xlim=(xmin,xmax))

        # For testing
        self.iz = iz

    def linearize(self,xx):


        ll = self.z(xx)
        nxx = self.iz(ll)

        return nxx

    def fit2user(self, pixels,features, delta):
        pixs=[]
        users=[]
        count = 1
        cuar = self.cuar
        for p,f in zip(pixels,features):
            indx = np.argmin(abs(cuar-f))
            diff = abs(cuar[indx]-f)
            if diff <= abs(delta):
                #print p,cuar[indx],diff, 'good  ',count
                count +=1
                pixs.append(p)
                users.append(cuar[indx])

        return (pixs,users)


    def aid_match(self, feature, delta):
        """
          Match a feature 'fea' againts a line list and
          return the closest match within 'delta', 
          otherwise None.

          feature: Wavelength units.

        """

        ll = self.cuar

        dars = abs(ll-feature)
        indx = np.argmin(dars)
        diff = dars[indx]

        #print 'AMO0:',feature,indx,ll[indx]
        if diff > delta: 
            out=None
            #print 'AM01:dif>delta:',diff,ll[indx-4:indx+4]
        else:
            out = ll[indx]
            #print 'AM02:',diff

        return out

    def clipData(self,pix,user,z,clip=4):
        """ Rejects data according to
            clip*stddev, where stddev is std(y-z(x)).
            If points are rejected, it will refit again.
            Return the new arrays x,y and the new function 'z'.
            Use the same fitname and order used in 'z'.

        """
        x = pix
        y = user

        norm = np.abs(y - z(x))
        stdev = np.std(norm)
        condition =  norm < clip * stdev

        y = np.asarray(y)[condition]
        x = np.asarray(x)[condition]

        fitname = z.fitname
        fitorder = z.order
        while norm.max() > clip * stdev:
            if len(y) < z.order + 1:
                print('Too few points left to fit!:'+str(len(y)))
                #self.log.warning('Too few points left to fit!:'+str(len(y)))
                return x,y,z
            z = self.fit_peaks(x,y,fitname=fitname, order=fitorder)
            norm = np.abs(y - z(x))
            stdev = norm.std()
            condition =  norm < clip * stdev

            y = y[condition]
            x = x[condition]

        return x,y,z

    def fit_peaks(self,xx,yy,fitname=None,order=None):
         if fitname == None:
             fitname = 'polynomial'
         if len(xx) == 2:
             order = 1
         elif len(xx) == 3:
             order = 2
         elif len(xx)> 3:
             order = 3
         else:
             print('Too few points left to fit!:'+str(len(xx)))
             return None

         z = gfit.Gfit(xx,yy,fitname=fitname, order=order)
         self.z = z
         return z

