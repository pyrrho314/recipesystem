import os
import imp
import time

import numpy as np
import scipy.interpolate as intp
from matplotlib import pyplot as pl

from astrodata import Lookups
import gfit

def find_upeaks(lpix,separation=2,nmax=60,cradius=12):
    """    
      Find spectra line maxima ranked according to peak intensity.

      lpix:     One row of pixel data. It depends on the caller
                if this has been median of several rows or just
                pick a row from an image.
      separation: minimum separation between 2 distinct peaks.
      nmax:     Maximum number of peaks to find.
      cradius:  centering radius in pixels 
    """
     
    #Get the line centers, fluxes (from the continuum) and widths
    gcen, copix, fwth = get_peak_centers(lpix,separation=separation)

    pcen,fmax, fw = findNtmax (gcen, fwth, copix, minsep=separation,
                               ntmax=nmax, radius=cradius)

    # Eliminate lines with zero width
    g, = np.where(fw!=0)
    if g.size != nmax:
        pcen = pcen[g]
        fmax = fmax[g]
        fw = fw[g]
    
    return pcen, fmax, fw

def get_peak_centers(lpix,separation=2):
    """ Get the line centers, fluxes 
        (from the continuum) and widths.
    """

    lpix = np.asarray(lpix)
    npix = lpix.size

    # cut the region in chunks of 60 pixels
    edges = [i for i in range(0,npix,60)]
    xx = np.arange(npix)
    er = np.zeros(npix)

    er[:] = np.std(lpix)
    co,cp,deriv = spline_continuum(xx, lpix, er, edges,debug=False)

    # See gaps:
    
    indl,indr= find_edges(lpix > abs(co))

    # pcen is a pretty good rough estimation of feature peak

    pcen,fmax,fwth = find_peaks(lpix, indl, indr)

    # Eliminate peaks that are too close
    g, = np.where(abs(pcen[1:]-pcen[:-1]) > separation)
    good = np.concatenate((g,[len(pcen)-1]))     # append last value
    pcen = pcen[good]
    fmax = fmax[good]
    fwth = fwth[good]

    return pcen, lpix-co, fwth 

def find_ww_cen1d(cen, lpix,width):
        """
          Find peaks from lpix with a
          rough estimation given by array cen.
          It uses center1d algorithm.

          lpix:   pixel data for the line 
          cen:    Array with rough estimation of the peak center
          width:  array with feature width in pixels at the base
        """
        pcen=[]
        s = ''
        for c,w in zip(cen,width):
            #newc = cen1d(c,lpix,0.0,width=10.0)
            newc = cen1d(c,lpix,0.0,width=w)
            if not newc:
                # Failed. Use a simpler centroid algorithm
                cen,cf,cw=find_peaks(lpix, [c-w/2], [c+w/2])
                if not cen:
                    print 'WW00:','none-->',c,w
                    newc=c
                    s += '%d '%c
                else:
                    newc = cen[0]
                    print 'WW11:',c,'-->',newc
            pcen.append(newc)
        return np.asarray(pcen),s

def findNtmax (pcen, fwth, copix, minsep=2, ntmax=60, radius=10):
    """
       Find upto ntmax peaks from a local_maxima algorithm.
       If one of these values is found in pcen, look for
       its width, otherwise assign a default value of 10.

       INPUTS:
         pcen:   array of peak positions.
         fwth:   width in pixels for each peak.
         copix:  Spectrum with continuum subtracted.
         minsep: Minimum separation between 2 peaks.
         ntmax:  Maximum number of ranked peaks to find.
         radius: Centering radius in pixels.
    """

    # Using local_maxima -since it found many peaks and the
    # widths from pcen, get a good set of peaks centres.

    # Substract the continuum
    
    lmx, = local_maxima(copix,min_distance=minsep)

    fmax = copix[lmx.tolist()]
    fmax=np.where(fmax>0.0,fmax,0.0)

    # NOTE: Needs to be a list, so use range instead of arange
    findx = range(len(fmax))
    findx.sort(lambda x,y: int(fmax[x] - fmax[y]))

    nmax = min(ntmax,len(lmx))    # At most is ntmax 

    cc = lmx[findx[-nmax:]]       # Get up to ntmax centres

    cf = fmax[findx[-nmax:]]      # Same for maximum peaks
    g = np.argsort(cc)            # Center of the maxfeatures ranked

    nn =0
    cen=[]; pw=[]; ffx=[]
    for k in g:
       # Lets look for those peaks in pcen that are in cc
       gg, = np.where(abs(pcen-cc[k])<2)
       if not gg: w = 10
       else:      w = fwth[gg]
       c1d = cen1d(cc[k],copix,0.0,width=w,radius=radius)
       if not c1d:
           for i in range(6):  # let's try 5 times
               if not c1d:
                  w = w-(i+1)
                  c1d = cen1d(cc[k],copix,0.0,width=w,radius=radius)
               else: break
           if i ==5: c1d = cc[k]   # just pick the initial value
       nn +=1
       c1d = float(c1d)
       cen.append(c1d); pw.append(int(w)),ffx.append(copix[c1d])

    return np.asarray(cen),np.asarray(ffx), np.asarray(pw) 

def local_maxima(array, min_distance = 1, periodic=False, edges_allowed=True):
    """Find all local maxima of the array, 
       separated by at least  min_distance.
       If uses ndimage.maximum_filter1d()
    """
    import scipy.ndimage as ndimage
    array = np.asarray(array)
    cval = 0
    if periodic:
       mode = 'wrap'
    elif edges_allowed:
       mode = 'nearest'
    else:
       mode = 'constant'
       cval = array.max()+1
    max_points = array == ndimage.maximum_filter1d(array,
                  1+2*min_distance, mode=mode, cval=cval)

    return [indices[max_points] for indices in  np.indices(array.shape)]


def find_edges(condition):
    """ Finds the indices for the edges of contiguous regions 
        where condition is True.

    Examples
    --------
    >>> a = np.array([3,0,1,4,6,7,8,6,3,2,0,3,4,5,6,4,2,0,2,5,0,3])
    >>> ileft, iright = find_edges_true_regions(a > 2)
    >>> zip(ileft, iright)
    [(0, 0), (3, 8), (11, 15), (19, 19), (21, 21)]

    """
    indices, = condition.nonzero()
    if not len(indices):
        return None, None
    iright, = (indices[1:] - indices[:-1] > 1).nonzero()
    ileft = iright + 1
    iright = np.concatenate( (iright, [len(indices)-1]) )
    ileft = np.concatenate( ([0], ileft) )
    return indices[ileft], indices[iright]

def find_peak (lpix, i0, i1):
    """
      Given the left and right edges of a line list,
      calculates the peak center by simply centroid algorithm
    """

    fl = lpix[i0:i1]
    wa = np.arange(i0,i1)
    try:
        ew = len(wa)*fl
        ewtot = np.sum(ew)
        wa_ewtot = np.sum(wa * ew)
        center = wa_ewtot / ewtot
    except:
        center = (i1-i0)/2.

    return center,max(fl),abs(i1-i0)


def find_peaks(lpix, indl, indr):
    """
      Given the left and right edges of a line list,
      calculates the peak center by simply centroid algorithm
    """

    centres = []
    max_flux = []
    wth = []
    for i0,i1 in zip(indl,indr):
        fl = lpix[i0:i1]
        wa = np.arange(i0,i1)
        if not len(wa): continue
        try:
            ew = len(wa)*fl
            ewtot = np.sum(ew)
            wa_ewtot = np.sum(wa * ew)
            center = wa_ewtot / ewtot
        except:
            center = (i1-i0)/2.

        centres.append(center)
        try:
           if i0==i1:
              print 'FNDPK00:',i0
           max_flux.append(max(fl))
        except:
           print 'FNDPK:',i0,i1
        wth.append(abs(i1-i0))

    return np.asfarray(centres),np.asfarray(max_flux),np.asfarray(wth)

def wrecenter(cen, lpix,width, cradius=10, maxw=15):
        """
          Find peaks from lpix with a
          rough estimation given by array cen.
          It uses center1d algorithm.

          lpix:   pixel data for the line
          cen:    Array with rough estimation of the peak center
          width:  constant in the meantime. Can be le,re
          cradius: 12.
                The maximum distance, in pixels, allowed between a line position
                and  the  initial  estimate  when  defining a new line
        """
        pcen=[]
        for c,w in zip(cen,width):
            if np.isnan(c): 
                pcen.append(c)
                continue
            newc = cen1d(c,lpix,0.0,w,radius=cradius)
            if not newc:
                for i in range(6):  # let's try 5 times
                   if not newc:
                      w = w-(i+1)
                      newc = cen1d(c,lpix,0.0,width=w,radius=10)
                   else: break
                if i ==5:
                    if (w<=1): newc=c
                    else:
                        w = min(maxw,w)
                        newc,f,tt = find_peak(lpix, max(0,c-w/2), c+w/2) 
                        if np.isnan(newc): newc=c
            pcen.append(float(newc))
        return np.asarray(pcen)

def recenter(cen, lpix,width=8, cradius=12):
        """
          Find peaks from lpix with a
          rough estimation given by array cen.
          It uses center1d algorithm.

          lpix:   pixel data for the line
          cen:    Array with rough estimation of the peak center
          width:  constant in the meantime. Can be le,re
          cradius: 12.
                The maximum distance, in pixels, allowed between a line position
                and  the  initial  estimate  when  defining a new line
        """
        pcen=[]
        w = width
        for c in cen:
            newc = cen1d(c,lpix,0.0,w,radius=cradius)

            if not newc:
                # Failed. Use a simpler centroid algorithm
                # Maybe find argmax(lpix[max(0,c-w/2):c+w/2]
                newc,cf,cw = find_peak(lpix, max(0,c-w/2), c+w/2)
            pcen.append(newc)

        return np.asarray(pcen)

def fitSpecLines(xx, yy, function='polynomial', order=3, debug=''):
        """
          Fit each of the spectral lines in a 2D image.

          INPUT:
          xx:    array of row numbers where yy is defined
          yy:    Pixel position for peaks found at rows xx.
                 This is a 2D array containing len(xx) arrays
                 of pixel positions for the spectral lines.
                 The spectral line 'n' pixel position for all
                 its measurements is: yy[:,n]
          function: 'polynomial'.  Fitting functions available:
                     ['polynomial','legendre', 'chebyshev', 'cubic']

          order:  3. Fittting order.

          debug:  If value is 'f', prints the coefficients for each
                  fitted spectral line.

          OUTPUT:
          zpeaks: Array of evaluator functions: zpeaks[row](pix)

        """

        # Fit each spectral line (npeaks of them)
        zpks=[]
        coef1=[]            # Keep the 1st non-constant coefficient
        if 'f' in debug:
            print "\n      off       coef[1]   coef[2]   coef[3]"
        npeaks = np.shape(yy)[1]
      
        x = np.asarray(xx)*1.0
        xmin = min(x[0],x[-1])
        xmax = max(x[0],x[-1])
        for k in range(npeaks):
            # Fit spectral line (peak) k
            y = np.asarray(yy[:,k])
            zp = gfit.Gfit(x, y,fitname=function,
                           order=order, xlim=(xmin,xmax))
            # We migth have some nan's
            if np.isnan(zp.coeff[0]): continue
            zpks.append(zp)
            coef1.append(zp.coeff[1])

            if 'f' in debug:
                for coe in zp.coeff:
                    print '%.3g'%coe,      # Coefficients
                print

        zpks = np.asarray(zpks)

        # Find any misbehaved fit 
        for i in range(1,len(zp.coeff)):
            c = [zpks[k].coeff[i] for k in range(len(zpks))]
            c = np.asarray(c)
            m = np.mean(c)     # Mean
            sig = np.std(c-m)  # Standard deviation
            g = np.where(abs(c-m)< 2*sig)    # Mark those curves that are off
            zpks = zpks[g]
        if 'f' in debug:
            print "number of zpeaks reduced from ",npeaks," to ",len(zpks)
            print "fitSpecLines:",len(coef1)-len(g[0])," fit functions eliminated."

        return zpks

def mapDistortion (zarcs,szim):
    """

       Distortion mapping function. For a given pixel it returns the
       shift value with reference to its position in the middle row.
       This 3-dim function is of the form:
          F(x,y) = A + B*x + C*y + D*x*x + E*y*y + F*x*y
       It is fit using the leastsq method from scipy.optimize. It
       minimizes the differences between F and the shift value given
       by the chebyshev fitting function zarcs[n](row), where
      'n' is the line index and 'row' in the image line number.

      INPUT
      zarcs:   Pixel position for spectral lines.
                This is a function evaluator of the form:
                zarcs[n](row). For any given spectral line index 'n' give its
                pixel position at image line 'row'.    
      szim:     Image size (nrows,npixels)
    """
    nrows = szim[0]
    nsamp = 200                   # LETS PICK 200 rows as max
                                 # (pixel, row, shift)
    x=[];y=[];z=[]
    npeaks = len(zarcs)
    
    nsum = 10

    ym = nrows/2

    # For each x get nrows/nsum y's and z's (shift)
    for k in  range(npeaks):
        zref = float(zarcs[k](ym))
        #for row in range(0,nrows,self.nsum):
        for row in range(0,nrows,10):
             xarc = float(zarcs[k](row))  # x-position of the arc at row
             y.append(row)
             z.append(xarc - zref)   # Evaluate at cen for each row
             x.append(xarc)

    x,y,z = map(np.asarray,(x,y,z))
    #print x.max(),y.max(),len(x),x.size,y.size,z.size,type(z)


    # Fit a surface (peak,row, shift(row,peak))
    xlim = (x.min(),x.max()); ylim = (y.min(),y.max())
    zz = gfit.Fit3d (x, y, z, xlim=xlim, ylim=ylim)

    return zz          # mapdist


def spline_continuum(wa, fl, er, edges, minfrac=0.01, nsig=3.0,
                     resid_std=1.3, debug=False):
    """ Given a section of spectrum, fit a continuum to it very
    loosely based on the method in Aguirre et al. 2002.

    From: Neil Crighton (neilcrighton@gmail.com)

    Parameters
    ----------
    wa               : Wavelengths.
    fl               : Fluxes.
    er               : One sigma errors.
    edges            : Wavelengths giving the chunk edges.
    minfrac = 0.01   : At least this fraction of pixels in a single chunk
                       contributes to the fit.
    nsig = 3.0       : No. of sigma for rejection for clipping.
    resid_std = 1.3  : Maximum residual st. dev. in a given chunk.
    debug = False    : If True, make helpful plots.

    Returns
    -------
    Continuum array, spline points, first derivative at first and last
    spline points

    Examples
    --------
    """

    # Overview:

    # (1) Calculate the median flux value for each wavelength chunk.

    # (2) fit a 1st order spline (i.e. series of straight line
    # segments) through the set of points given by the central
    # wavelength for each chunk and the median flux value for each
    # chunk.

    # (3) Remove any flux values that fall more than nsig*er below
    # the spline.

    # Repeat 1-3 until the continuum converges on a solution (if it
    # doesn't throw hands up in despair! Essential to choose a
    # suitable first guess with small enough chunks).

    if len(edges) < 2:
        raise ValueError('must be at least two bin edges!')

    wa,fl,er = (np.asarray(a) for a in (wa,fl,er))

    if debug:
        ax = pl.gca()
        ax.cla()
        ax.plot(wa,fl)
        ax.plot(wa,er)
        ax.axhline(0, color='0.7')
        good = ~np.isnan(fl) & ~np.isnan(er)
        ymax = 2*sorted(fl[good])[int(len(fl[good])*0.95)]
        ax.set_ylim(-0.1*ymax, ymax)
        ax.set_xlim(min(edges), max(edges))
        ax.set_autoscale_on(0)
        pl.draw()

    npts = len(wa)
    mask = np.ones(npts, bool)
    oldco = np.zeros(npts, float)
    co = np.zeros(npts, float)

    # find indices of chunk edges and central wavelengths of chunks
    indices = wa.searchsorted(edges)
    indices = [(i0,i1) for i0,i1 in zip(indices[:-1],indices[1:])]
    if debug:  print ' indices',indices
    wavc = [0.5*(w1 + w2) for w1,w2 in zip(edges[:-1],edges[1:])]

    # information per chunks
    npts = len(indices)
    mfl = np.zeros(npts, float)     # median fluxes at chunk centres
    goodfit = np.zeros(npts, bool)  # is fit acceptable?
    res_std = np.zeros(npts, float) # residuals standard dev
    res_med = np.zeros(npts, float) # residuals median
    if debug:
        print 'chunk centres',wavc
        cont, = ax.plot(wa,co,'k')
        midpoints, = ax.plot(wavc,mfl,'rx',mew=1.5,ms=8)

    # loop that iterative fits continuum
    while True:
        for i,(j1,j2) in enumerate(indices):
            if goodfit[i]:  continue
            # calculate median flux
            w,f,e,m = (item[j1:j2] for item in (wa,fl,er,mask))
            ercond = e > 0
            cond = m & ercond
            chfl = f[cond]
            chflgood = f[ercond]
            if len(chflgood) == 0: continue
            #print len(chfl), len(chflgood)
            if float(len(chfl)) / len(chflgood) < minfrac:
                f_cutoff = scoreatpercentile(chflgood[ercond], minfrac)
                cond = ercond & (f >= f_cutoff)
            if len(f[cond]) == 0:  continue
            mfl[i] = np.median(f[cond])
        # calculate the spline. add extra points on either end to give
        # a nice slope at the end points.
        extwavc = ([wavc[0]-(wavc[1]-wavc[0])] + wavc +
                   [wavc[-1]+(wavc[-1]-wavc[-2])])
        extmfl = ([mfl[0]-(mfl[1]-mfl[0])] + list(mfl) +
                  [mfl[-1]+(mfl[-1]-mfl[-2])])
        co = np.interp(wa,extwavc,extmfl)
        if debug:
            cont.set_ydata(co)
            midpoints.set_xdata(wavc)
            midpoints.set_ydata(mfl)
            pl.draw()

        # calculate residuals for each chunk
        for i,(j1,j2) in enumerate(indices):
            if goodfit[i]:  continue
            ercond = er[j1:j2] > 0
            cond = ercond & mask[j1:j2]
            chfl = fl[j1:j2][cond]
            chflgood = fl[j1:j2][ercond]
            if len(chflgood) == 0:  continue
            if float(len(chfl)) / len(chflgood) < minfrac:
                f_cutoff = scoreatpercentile(chflgood[ercond], minfrac)
                cond = ercond & (fl[j1:j2] > f_cutoff)
            #print len(co), len(fl), i1, j1, j2
            residuals = (fl[j1:j2][cond] - co[j1:j2][cond]
                         ) / er[j1:j2][cond]
            res_std[i] = residuals.std()
            if len(residuals) == 0:
                continue
            res_med[i] = np.median(residuals)
            # If residuals have std < 1.0 and mean ~1.0, we might have
            # a reasonable fit.
            if res_std[i] <= resid_std:
                goodfit[i] = True

        if debug:
            print 'median and st. dev. of residuals by region - aiming for 0,1'
            for i,(f0,f1) in  enumerate(zip(res_med, res_std)):
                print '%s %.2f %.2f' % (i,f0,f1)
            raw_input('Enter...')

        # (3) Remove flux values that fall more than N*sigma below the
        # spline fit.
        cond = (co - fl) > nsig * er
        if debug:
            print np.nanmax(np.abs(co - oldco)/co)
        # Finish when the biggest change between the new and old
        # medians is smaller than the number below.
        if np.nanmax(np.abs(co - oldco)/co) < 4e-3:
            break
        oldco = co.copy()
        mask[cond] = False

    # finally fit a cubic spline through the median values to
    # get a smooth continuum.
    d1 = (mfl[1] - mfl[0]) / (wavc[1]-wavc[0])
    d2 = (mfl[-1] - mfl[-2]) / (wavc[-1]-wavc[-2])
    final = gfit.InterpCubicSpline(wavc, mfl, firstderiv=d1, lastderiv=d2)

    return final(wa), zip(wavc,mfl), (d1,d2)


# algorithm from scipy
def scoreatpercentile(a, perc):
    """Calculate the score at the given 'perc' percentile of the
    sequence a.  For example, the score at per=50 is the median.

    'perc' can be a sequence of percentile values.

    If the desired quantile lies between two data points, we linearly
    interpolate between them.
    """
    # TODO: this should be a simple wrapper around a well-written quantile
    # function.  GNU R provides 9 quantile algorithms (!), with differing
    # behaviour at, for example, discontinuities.
    vals = np.sort(a, axis=0)

    if not hasattr(perc, '__iter__'):
        perc = [perc]

    out = []
    for p in perc:
        i = 0.01 * p * (vals.shape[0] - 1)
        j = int(i)
        if (i % 1 == 0):
            out.append(vals[j])
        else:
            # linearly interpolate
            out.append(vals[j] + (vals[j+1] - vals[j]) * (i % 1))


    return np.array(out).squeeze()
    
def cen1d(x, data, threshold, width=3., type='emission', radius=12.):
    """
    This is a translation of the iraf/center1d.x procedure.

    Locate the center of a one dimensional feature.
    A value of INDEF is returned in the centering fails for any reason.
    This procedure just sets up the data and adjusts for emission or
    absorption features.  The actual centering is done by C1D_CENTER.
    If twidth <= 1 return the nearest minima or maxima.

    real    x               # Initial guess
    int     npts            # Number of data points
    real    data[npts]      # Data points
    real    width           # Feature width
    type                    # Feature type (emission, absorption)
    real    radius          # Centering radius
    real    threshold       # Minimum range in feature

    """

    # We need to define a class and set this at __init__ time
    MIN_WIDTH = 3.              # Minimum centering width

    if type not in ['emission','absorption']:
       raise RuntimeError,"center1d 'type' value error" 

    wid = max (3., width)
    rad = max (2., radius)

    # Determine the pixel value range around the initial center, including
    # the width and error radius buffer.  Check for a minimum range.
    #npts = len(data)
    npts = data.size - 1
    if x >= npts: return None
    x1 = max (0, x - wid / 2 - rad - wid)
    x2 = min (npts, x + wid / 2 + rad + wid)
    x1,x2 = map(np.int,(x1,x2))
    a = data[x1:x2+1].min()
    b = data[x1:x2+1].max()
    if (b - a) < threshold:
       return (None)

    # calculate the continuum subtracted data vector.  The X
    # range is just large enough to include the error radius and the
    # half width.

    x1 = max (0, x - wid / 2 - rad)
    x2 = min (npts, x + wid / 2 + rad + 1)
    #x2 = min (npts-1, x + wid / 2 + rad + 1)
    x1,x2 = map(np.int,(x1,x2))
    nx = x2 - x1 + 1 


    sdata = data[x1:x2+1].copy()

    # Make the centering data positive, subtract the continuum, and
    # apply a threshold to eliminate noise spikes.

    if type == 'emission':
        a = min (0., a)
        sdata = sdata - (a + threshold)
    elif type == 'absorption':
        sdata = -sdata
        sdata = sdata - (threshold-b)
    smax = sdata.max()
    sdata = np.clip(sdata,0.,smax)
        
    #print x1,x2,x,nx,len(sdata)
    xc = c1d_center (x - x1, sdata, nx, width)

    # Check user centering error radius.
    if xc:       # Value is not None
        xc = xc + x1
        if np.abs(x - xc) > radius:
            xc = None
    return xc    

def c1d_center (x, data, npts, width):
    """
     One dimensional centering algorithm.
     If the width is <= 1. return the nearest local maximum.

     real    x            # Starting guess
     int     npts         # Number of points in data vector
     real    data[npts]   # Data vector
     real    width        # Centering width

    """

    # SOME OF THESE ARE ALREADY DEFINED IN THE CALLING ROUTINE
    MIN_WIDTH  = 3.              # Minimum centering width
    EPSILON    = 0.001           # Accuracy of centering
    EPSILON1   = 0.005           # Tolerance for convergence check 
    MAX_DXCHECK = 3              # Look back for failed convergence
    ITERATIONS = 100             # Maximum number of iterations

    epsilon = EPSILON
    
    # Find the nearest local maxima as the starting point.
    # This is required because the threshold limit may have set
    # large regions of the data to zero and without a gradient
    # the centering will fail.
    
    #ii = int(x+1)    # equivalent to int(ceil(x))
    ii = int(x)
    npm = npts-1
    for i in range(ii,npm):
        if data[i] > data[i+1]:
           break
    while (i > 0) and (data[i]<=data[i-1]): i = i - 1

    #j = round(x+.5)
    j = int(x)
    while (j > 0) and (data[j]<=data[j-1]): j = j - 1 

    while (j < npm) and (data[j]<=data[j+1]): j = j + 1

    if np.abs(i-x) < np.abs(x-j):
        xc = j
    else:
        xc = i

    if width <= 1: 
        return xc
    
    wid = max (width, MIN_WIDTH)

    # Check data range.
    hwidth = wid / 2

    if (xc - hwidth < 1) or (xc + hwidth > npts):
        return (None)

    # Set interpolation functions.

    xx = np.arange(len(data))
    try:
       z1 = intp.InterpolatedUnivariateSpline(xx,data)
    except:
       print "Z1 Exception",xc

    # interpolate the x*y values

    try:
       z2 = intp.InterpolatedUnivariateSpline(xx,data*xx)
    except:
       print "Z2 Exception",xc
 
    # Iterate to find center.  This loop exits when 1) the maximum
    # number of iterations is reached, 2) the delta is less than
    # the required accuracy (criterion for finding a center), 3)
    # there is a problem in the computation, 4) successive steps
    # continue to exceed the minimum delta.

    dxlast = npts
    dxcheck = 0
    iter = 0
    while iter < ITERATIONS:
        # Triangle centering function.
        a = xc - hwidth
        b = xc - hwidth/2
        intgrl1 = z1.integral(a, b)
        intgrl2 = z2.integral(a, b)
        sum1 = (xc - hwidth) * intgrl1 - intgrl2
        sum2 = -intgrl1
        a = b
        b = xc + hwidth/2
        intgrl1 = z1.integral(a, b)
        intgrl2 = z2.integral(a, b)
        sum1 = sum1 - xc * intgrl1 + intgrl2
        sum2 = sum2 + intgrl1
        a = b
        b = xc + hwidth
        intgrl1 = z1.integral(a, b)
        intgrl2 = z2.integral(a, b)
        sum1 = sum1 + (xc + hwidth) * intgrl1 - intgrl2
        sum2 = sum2 - intgrl1

        # Return no center if sum2 is zero.
        if sum2 == 0.:
            break

        # Limit dx change in one iteration to 1 pixel.
        dx = sum1 / np.abs (sum2)
        dxabs = np.abs (dx)
        xc = xc + max (-1., min (1., dx))
        if (xc - hwidth < 0) or (xc + hwidth > npts-1):
            break

        # Convergence tests.
        if dxabs < epsilon:
            return xc
        if dxabs > (dxlast + EPSILON1):
            dxcheck = dxcheck + 1
            if dxcheck > MAX_DXCHECK:
                return None
        elif dxabs > (dxlast - EPSILON1):
            xc = xc - max (-1., min (1., dx)) / 2
            dxcheck = 0
        else:
            dxcheck = 0
            dxlast = dxabs
        iter = iter + 1
    # If we get here, no center was found.
    return None

def _aargs(*pars):
    """ Utility function to turn all the arguments into
        strings and concatenate them into one with a blank
        in between.
    """
    s=''
    for k in pars: s += str(k)
    return s

def findMOSEdges(hdulist):
        """
          findMOSEdges is a take from gmos/gscut.cl to find the
          REFPIX value.
          It also finds the bottom and top edges of each slice by taking
          a vertical bar on CCD2 at 1200:1500
          (vbar = np.mean(bigpix[:,1200:1500],axis=1))
        """

        from math import sin,cos 
        asecmm = 1.611444
        
        pixscale ={'GMOS_N': 0.0727 ,'GMOS-S': 0.073}

        bigpix = hdulist['SCI',1].data
        
        pi = np.pi

        phu = hdulist[0]._header
        header = hdulist['SCI',1]._header
        tb = hdulist['MDF',1].data
     
        inst = phu['INSTRUME']
        
        grating = phu['grating']

        # Check for  data file location
        #fp, pathname, description = imp.find_module('gwavecal')
        #dirn = os.path.dirname(pathname)
        #reffile = os.path.join(dirn, 'StandardGMOSGratings.py')

        print 'fm00: dirni::',dirn,pathname

        GMOSfilters = Lookups.get_lookup_table('Gemini/GMOS/GMOSfilters.py','GMOSfilters')

        StandardGMOSGratings = Lookups.get_lookup_table('Gemini/GMOS/StandardGMOSGratings.py',
                                    'StandardGMOSGratings')

        grule, gblaze, gR, gcoverage, gwave1, gwave2,\
        wavoffset, l_yoff = StandardGMOSGratings[grating]

        filter1 = phu['filter1']
        filter2 = phu['filter2']

        # get filter information
        fwave1 = 0.0     ; wmn1 = 0.0     ; wmn2 = 0.0
        fwave2 = 99999.0 ; wmx1 = 99999.0 ; wmx2 = 99999.0
        if 'open' not in filter1:
            wmn1, wmx1, ffile = GMOSfilters[filter1]
        if 'open' not in filter2:
            wmn2, wmx2, ffile = GMOSfilters[filter2]

        cwave = phu['grwlen']
        tilt = phu['grtilt']
        
        tilt = tilt*pi/180.
        ss = header['ccdsum']     # Read a string 'x y'
        xbin,ybin = float(ss[0]),float(ss[2])

        xscale = pixscale[inst]*xbin
        yscale = pixscale[inst]*ybin

        greq = cwave*grule/1.e6
        gratinfile = os.path.join(dirn, 'gratingeq.dat')
        x,y = np.loadtxt(gratinfile,unpack=True)
        z = gfit.Gfit(x,y,'cubic')
        gtilt = z(greq)

        gtilt = gtilt * pi/180.
        a = sin(gtilt+0.872665) / sin(gtilt)
        gR = 206265. * greq/(0.5*81.0*sin(gtilt))
        nmppx = a*xscale*cwave*81.0*sin(gtilt)/(206265.*greq)
        wave1 = gwave1
        wave2 = gwave2

        fwave1 = max(wmn1,wmn2)
        fwave2 = min(wmx1,wmx2)

        # determine whether filter or grating limits wavelength coverage
        wave1 = max(wave1,fwave1)
        wave2 = min(wave2,fwave2)

        # in pixels
        speclen = round((wave2-wave1)/nmppx)

        # pixel value of central wavelength from left (red) end of spectrum
        #crpix1?
        pixcwave = speclen - (cwave-wave1)/nmppx

        nypix, nxpix = np.shape(bigpix)
        xcen = nxpix/2.
        ycen = nypix/2.



        sx = tb.field('slitpos_mx')
        sy = tb.field('slitpos_my')
        zx = tb.field('slitsize_mx')
        zy = tb.field('slitsize_my')
        pr = tb.field('priority')

        #loop over the slits
        g = np.argsort(sy)
        sx = sx[g]
        sy = sy[g]
        zx = zx[g]
        zy = zy[g]
        pr = pr[g]

        # Determines the edges of the stripes
        sz = np.shape(bigpix)
        cm = sz[1]/2
        vbar = np.mean(bigpix[:,cm-150:cm+150],axis=1)
        lowa, topa = mos_edges(vbar)
        cdif = len(sx)-len(lowa)
        # make these array the same length
        if cdif > 0:
            print "WARNING::::::: The number os slits is less the number of entries in tb:"
            print "               by:",cdif,"\n"
            lowa = np.concatenate(([0],lowa , [0]))
            topa = np.concatenate(([0],topa , [0]))
        print "NUMBER of slits by mos_edges():",len(lowa)


        #TODO 
        #     put tilt
        crpix = [] 
        zzpeaks = []
        #for spos_mx,spos_my,ssize_mx,ssize_my,priority,le,te in zip(sx,sy,zx,zy,pr,lowa,topa):
        k = 1
        for spos_mx,spos_my,ssize_mx,ssize_my,priority in zip(sx,sy,zx,zy,pr):

            k += 1
        #    if k > 3: break

            # Convert from mask to pixel coordinates and correct for
            # binning in both directions
            #if priority == 0: continue

            xccd = spos_mx * asecmm/xscale
            if (inst=='GMOS-S'):
                # yccd=spos_my*asecmm/yscale
                # Not only there is a y-offset (85) but there is also 
                # a distortion.  The solution below is not perfect but
                # it is already much better than the first order
                # solution [Kathleen Labrie]
                yccd = 0.99911*spos_my - 1.7465E-5*spos_my**2 + \
                    3.0494E-7*spos_my**3
                yccd = yccd * asecmm/yscale
            else:
                yccd = 0.99591859227*spos_my + \
                    5.3042211333437E-8*spos_my**2 + \
                    1.7447902551997E-7*spos_my**3
                yccd = yccd * asecmm/yscale

            slitwid = ssize_mx*asecmm
            slitlen = ssize_my*asecmm

            # set slit length if the aperture is a circle,
            # radius=slitwid
            #if (slittype=="circle")
            #    slitlen=2.*slitwid

            xccd = xcen+xccd
            yccd = ycen+yccd

            # simple correction for distortion in x
            y = (yccd/nypix - 0.5)
            dx = nxpix * (0.0014*y - 0.0167*y**2)
            #print(yccd," ",y," ",dx)

            # slit height
            specwid = round(1.05*slitlen/yscale)
            center = specwid/2

            refpix = pixcwave

            # Position of object, take into account that lambda decreases
            # with x
            x1 = round(xcen-(xcen-xccd)/a-pixcwave) + wavoffset/nmppx + dx
            x2 = x1 + speclen-1
            y1 = round(yccd-center+l_yoff)
            y2 = y1 + specwid-1
            #print k,': (%.2f %.2f) (%.2f %.2f)'%(x1,x2,y1,y2),
            # check spectrum isn't off chip
            if x1 < 1:
                refpix = refpix+x1-1.
                x1 = 1

            if x2 > nxpix:
                x2 = nxpix
            if y1 < 1:
                y1 = 1
            if y2 > nypix:
                y2 = nypix

            crpix.append(refpix)
            #ss= '%.2f %.2f %.4f'% (refpix,cwave*10,-10.*nmppx),le,te,te-le
            #self.log.info(ss)

            """
              having the slits lines all vertical, i.e. applied rotation already.
              we now fit each of the slits
              See if we have shift of the lines in each of the slits
            """
        crval = cwave*10
        cdelt = -10.*nmppx

        return np.asarray(crpix), lowa, topa, crval, cdelt

def mos_edges(cpix):
          """
              cpix: one column tru the array, we pick a
                    band of columns and take the mean
              Find edges in gmos mos slit spectra to limit
              the band for processing.
          """
          # NOTE: The bottom slit would not have a sharp edge
          # and it will escape this algorithm. Need to count the
          # number of indices and add the 1st one (bottom) if necessary.  
           
          spix = smooth(cpix)
          dd = spix[1:] - spix[:-1]
          dev = np.std(dd)          # Positive values are the start of edge
                                    # Negative is end of edge
          cond = dd > dev/2
          indices, = cond.nonzero()
          ipos = (indices[1:] - indices[:-1]) > 1
          ipos = np.concatenate( (ipos, [len(indices)-1]) )
          g = np.where(ipos)
          startI = indices[g]

          cond = dd < -dev/2
          indices, = cond.nonzero()
          ipos = (indices[1:] - indices[:-1]) > 1
          ipos = np.concatenate( (ipos, [len(indices)-1]) )
          g = np.where(ipos)
          endI = indices[g] - 2       # We need to be inside the cut
          return startI,endI

def plotpeaks(peaks,lpix,color='g'):
    """Utility function: Plot peaks:

       peaks: peaks pixel position in lpix array
       lpix:  1D spectrum
    """

    pl.clf()
    pl.plot(lpix)
    ys = lpix[np.int32(peaks)]
    pl.vlines(peaks,ys+20,ys+ys/10,color=color)
    pl.axis([0,len(lpix),0,np.max(lpix)])

def impeak(imdata):
    """Utility function: Plot middle row peaks.
       imdata: 2D image containing spectra
    """

    sz = np.shape(imdata)
    ym = sz[0]/2           # Starting row (middlerow)
    nsum = 10
    separation = 2
    lpix = np.mean(imdata[ym-nsum/2:ym+nsum/2,:],axis=0)
    pcen,fmax,fw = find_upeaks(lpix,separation=2,nmax=60,cradius=12)
    
    xcen = pcen.copy()
    low = 80            
    for y in range(ym-nsum,low,-nsum):

        mpix = np.mean(imdata[y-nsum:y,:],axis=0)
        cen = wrecenter(xcen, mpix, fw, cradius=10)
        a,b,c,d=[cen[i]-xcen[i] for i in [0,1,30,50]]
        #e,f,g,h=[cen[i]-fcen[i] for i in [0,1,30,50]]
        print y,'%.2f %.2f %.2f %.2f'%(a,b,c,d),xcen[1]
        #print '   %.2f %.2f %.2f %.2f'%(e,f,g,h)
        xcen = cen.copy()

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming',
                'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead 
          of a string   
    REF: From scipy/cookbook/signalsmooth
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning',"
         "'hamming', 'bartlett', 'blackman'")


    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]


def resample(a,newdims,centre=False):
     '''Arbitrary resampling of source array to new dimension sizes.
        Currently only supports maintaining the same number of dimensions.
        To use 1-D arrays, first promote them to shape (x,1).

        centre:
           True - interpolation points are at the centres of the bins
           False - points are at the front edge of the bin

     '''

     if not a.dtype in [np.float64, np.float32]:
         a = np.cast[float](a)

     minusone = False
     m1 = np.cast[int](minusone)
     ofs = np.cast[int](centre) * 0.5
     old = np.array( a.shape )
     ndims = len( a.shape )
     if len( newdims ) != ndims:
         print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
         return None
     newdims = np.asarray( newdims, dtype=float )
     base = np.arange(newdims)
     dim = (old - m1) / (newdims - m1) * (base + ofs) - ofs 

     # specify old dims
     olddims = np.arange(a.shape[0], dtype = np.float)
     
     # first interpolation - for ndims = any
     mint = intp.interp1d( olddims, a, kind='linear' )
     newa = mint( dim )
     
     return newa
