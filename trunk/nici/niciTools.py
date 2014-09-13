import pyfits as pf
import numpy as np
import re
import glob
import time
import os

def reset_nans(im):
    """
     Reset pixels with Nan value to
     99 (arbitrarily). All other pixels to zero.
    """
    bp = im.copy()
    #bp[np.isfinite(im)]=0
    bp[:] = 0
    gnan = np.where(np.isnan(im))
    bp[gnan]=99.
    im=np.nan_to_num(im)
    
    return bp,im


def restore_nans(im,bp):
    """
      Restore those pixels with values > 50 to nans.
      This procedure works in conjunctio with reset_nans.
    """
    im[bp>5.] = np.nan
    return im
   
def check_dir(dirname,flag):
    """
      Check that directory name have ending '/',
      exist and can be accessed with 'flag' access.
    """
    dirname = os.path.join(dirname.strip(),'')

    if dirname == '': return dirname

    # Verify that directory exist
    if (not os.access(dirname,os.F_OK)):
       errs = '"'+dirname+'"'+' does not exist.'
       raise IOError(errs)
    if (not os.access(dirname,flag)):
       acc = {os.R_OK:'read',os.W_OK:'write'}[flag]
       errs = '"'+dirname+'"'+' does not have '+acc+' access.'
       raise IOError(errs)
    return dirname

def gstat(files,path):
    """
      List information about each file in the input
      list. FILE GCALSHUT CRMODE TYPE MEDR MEDB
      where medr and medb are the median value of the
      red and blue frame respectively.
      
    """

    if (path != '' and path[-1] != '/'):
       path += '/'

    print 'FILE GCALSHUT CRMODE TYPE MEDR MEDB'
    for f in files:
        un = pf.open(path+f)
	if (un[0].header).get('instrume') != 'NICI': 
            print f,'is not NICI file.'
            un.close()
            continue
	if (un[1].header).get('FILTER_R'): 
            exr=1; exb=2
        else:
            exr=2; exb=1
        gcal=(un[0].header).get('GCALSHUT')
        otype=(un[0].header).get('OBSTYPE')
        crmode=(un[0].header).get('crmode')
        medb = np.median(un[exb].data,axis=None)
        medr = np.median(un[exr].data,axis=None)
        print f,gcal,otype,crmode ,'(%.2f)' % medr, '(%.2f)' % medb
        un.close()
        #h=iraf.images.hselect(f+'[0]','$I,GCALSHUT,OBSTYPE','yes',Stdout=1)
        #iraf.images.imstat(f+'[1]',format='no',Stdout=1)
        #iraf.images.imstat(f+'[2]',format='no',Stdout=1)

    return 
def create_list(pattern, file=None):
    """
      Create a gemini file name list. The pattern needs to be of
      the form S20090312S+N-M+(.fits), and/or S20090312S+N+(.fits) and/or
      S20090312S0013(.fits). The extension '.fits' is optional.
      file: Creates a filelist is a file is specified.

      Returns a list and a file if one is giving.

    """
    list = []
    files = []
    list = pattern.split(',')
    # split pattern  root+n-m+ into a list of 
    # [root+str(seq) for seq in range(n,m)]
    rex = r'(?P<root>\w*)\+(?P<r1>\d+)-?(?P<r2>\d*)(\D*)'
    saveroot = ''
    for fpat in list:
        m = re.search(rex, fpat)
        root = m.group('root').strip()
        if len(root) > 1:
            if root[-1] == '+':
                root = root[:-1]
        a = int(m.group('r1'))
        b = m.group('r2')
        if root != '':
            saveroot = root
        if b != '':              # we have a range
            b = int(b)
            for n in range(a,b+1):
                files.append(root + str('%.4d'%n) + '.fits')
        else:
            if saveroot == '':
                print "ERROR: The first element in the input list should \
                       have a root string"
                break
            files.append(saveroot + str('%.4d'%a) + '.fits')
               
    return files

def getFileList(input):
    """
    input: (string) filenames separated by commas or an @file.
           If a filename contains a Unix wildcard (?,*,[-])
           then it will attemp to generate the corresponding
           list of files that matches the template.

    """
    filenames=[]
    if 'str' in str(type(input)):
        if input[0] == '@' and input.find (',') < 0:
            fd = open(input[1:])
            flist = fd.readlines()
            fd.close()
        else:
            flist = input.split(',')
    else:
        if len(input) == 1:
           list = getFileList(input[0])
           filenames += list
           return filenames
        else:   
           flist = input

    for line in flist:
        line = line.strip()
        if len(re.findall(r'\w?\+\d*-\d*\+',line)):     # match 'root+N-N+'
            filenames += create_list(line)
        elif len(re.findall (r'[?,*,-]',line)):
           # We got a wilcard as a filename
           line = glob.glob(line)
           filenames += line
        elif '@' in line:
            list = getFileList(line)
            filenames += list
        else:
            if (len(line.strip()) == 0):
                continue
            filenames.append (line)

    return filenames


def gen_list(wildname, keyname, value):
    """
      list = gen_list(wildname, keyname, value):
      arguments:
      wildname: wild card string: eg: *.fits
      keyname: keyword names
      value: string contained in keyname value
      returns: filename list.      

      Generate a list of filenames that matches the 'value' 
      in a keyword. The keyword is taken from the global extension,
      i.e. the PHU
    """
    import glob

    file_list = []
    flis = glob.glob(wildname)
    for f in flis:
        fits=pf.open(f)
        objval=(fits[0].header).get(keyname)
        if (objval == None):
           continue
        if (value.upper() in objval.upper()):
           file_list.append(f)
        fits.close()
    return np.sort(file_list)


def medbin(im, dimx, dimy):

   naxis = np.size(np.shape(im))
   if (naxis == 2):
      dim = [dimx, dimy]
      sz = np.shape(im)
      bx = sz[1] / dim[1]
      by = sz[0] / dim[0]

      index = np.arange(sz[1]*sz[0]).reshape(sz[1],sz[0])
      iy = index % sz[0]
      ix = index / sz[0]
      ipix = ((index % by) + ((ix) % bx) * by)

      iy = iy / by
      ix = ix / bx
      ibox = iy + ix * dim[1]

      v = np.zeros([bx * by, dim[0] * dim[1]], np.float)
      v[ipix,ibox] = im

      #IDL:v = median(v, dim=2, even=True)
      v = np.median(v,axis=0)   # The 1st python dimension is idl's 2nd.

      out = np.zeros([dim[1], dim[0]], np.float)
      sz = np.shape(out)
      out = v.reshape(sz[1],sz[0])

      return out


   if (naxis == 1):
      out = np.zeros([dimx], np.float)
      for i in np.arange(0, (dimx - 1)+(1)):
         out[i] = np.median(im[i*(sz[1]/dimx) : ((i+1)*(sz[1]/dimx)-1)+1])
      return outflat_file

def rebin_simple( a, newshape ):
    '''
      Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)

    slices = [slice(0,old, np.float(old)/new) for old,new in zip(a.shape,newshape)]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')  #choose the biggest smaller integer index
    return a[tuple(indices)]

def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    lenShape = np.size(np.shape(a))
    factor = np.asarray(np.shape(a))/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    #print ''.join(evList)
    return eval(''.join(evList))

def parangle(ha,dec,lat):
    """
    Return the parallactic angle of a source in degrees.

    HA - the hour angle of the source in decimal hours; a scalar or vector.
    DEC - the declination of the source in decimal degrees; a scalar or
            vector.
    LAT - The latitude of the telescope; a scalar.
    """
    from numpy import pi,sin,cos,tan 

    har = np.radians(15*ha)
    decr = np.radians(dec)
    latr = np.radians(lat)

    ac2 = np.arctan2( -sin(har), \
                 cos(decr)*tan(latr) - sin(decr)*cos(har))

    return -np.degrees(ac2)

def dmstod(dms):
    """ 
    return decimal degrees from +/-dd:mm:ss.dd
    """


    dms=dms.strip()
    if dms[0] == '-':
       sign = -1.
    else:
       sign = 1.
    dar = np.asfarray(dms.split(':'))
    dar[0] = sign*dar[0]

    return sign*(dar[0] + dar[1]/60. + dar[2]/3600.)

def dtodms(deg):
    """
    Turns degrees (a floating) into +/-dd:mm:ss.ddd
    """

   # Separate out the sign.
    if deg < 0.0:
        sign = '-'
        deg = -deg
    else:
        sign = '+'

    # Convert to an integer to avoid inconsistent rounding later.
    # The least significant digit of the result corresponds to
    # centiseconds so multiply the number accordingly.  Note that
    # int rounds towards zero.
    csec = int(deg * 60.0 * 60.0 * 100.0 + 0.5)

    # Split up into the four components.  divmod calculates
    # quotient and remainder simultaneously.
    sec, csec = divmod(csec, 100)
    min, sec = divmod(sec, 60)
    deg, min = divmod(min, 60)

    # Convert the four components and sign into a string.  The
    # % operator on strings works like C's sprintf function.
    return '%s%d:%02d:%2d.%02d' % (sign, deg, min, sec, csec)

def printlog_open(logfile):
    """
    Open logfile for appending
    BETTER TO USE THE Python Logging facility
    """
    lg = open(logfile,'a')
    return lg

def print_error(lg,error):
    pass
    # need to write a class

def nici_noise(imc):
    import copy
    """
     A little function that subtracts the nici readout
     pattern noise every 16 pixels
    """
    im = copy.deepcopy(imc)
    line= np.zeros(16)
    for q in range(2):
        for i in range(128):
            tmp=im[i*8:i*8+8, q*512:q*512+512].ravel()
            for j in range(16): 
                line[j] = np.median(tmp[np.arange(256)*16+j])
            line -= np.median(line)
            for j in range(16):  
               tmp[np.arange(256)*16+j] -= line[j]
            im[i*8:i*8+8, q*512:q*512+512] = tmp.reshape(8,512)

    return im

def fits_utc():
   """
   Return a UTC string in FITS format:
   YYYY-MM-DDThh:mm:ss
   """

   gmt = time.gmtime()
   time.asctime(gmt)
   fitsT = '%d-%02d-%02dT%02d:%02d:%02d' % gmt[:6]

   return fitsT

def order_wcs(header):
   """
   Order the WSC information in the FITS header according to
   RADECSYS, CTYPEn,CRPIXn,CRVALn and CDn_m for n,m [1:3].
   """
   nwcs=[]
   hcc = header.ascard
   for pat in ['RADECSYS','CTYPE.','CRPIX.','CRVAL.', 'CD._.']:
       for h in hcc.filterList(pat):
           nwcs.append(h)
           kw = str(h).split('=')[0].strip()
           header.__delitem__(kw)
   for cc in nwcs:
        hcc.append(cc)


def parallacticAngle (filelist, fdir='',extn=1):
    """
      Calculates the parallactic angle from a FITS header WCS cd matrix.
      Reads the first FITS extension to get the WCS values  ([1])
      
      @filelist:  A list of FITS files
      @type filelist: Python list
      @fdir:      Directory pathname where FITS files are located.
      @type fdir: String (Default value is '')  
      @extn:      FITS extension number to open
      @type extn: int (Default value is 1)

      Returns: A list of corresponding parallactic angles.
    """
    # Get information from the EHU to calculate the Parallactic angle
    #
    fdir = os.path.join(fdir,'')
    radeg = 180/np.pi
    pa = []
    for fi in filelist:
        fi = fdir+fi
        hdr = pf.getheader(fi,extn)
        cd11=hdr['cd1_1']
        cd12=hdr['cd1_2']
        cd21=hdr['cd2_1']
        cd22=hdr['cd2_2']
        dd = 180 - np.arctan2(sum([cd11,cd22]),sum([cd21,-cd12]))*radeg
        pa.append(dd)

    return pa

#import numpy as n
import scipy.interpolate
import scipy.ndimage

def congrid(a, newdims, method='linear', centre=False, minusone=True):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)
    
    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )

    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    dimarg=newdims
    newdimsr = np.asarray( dimarg, dtype=np.float )    

    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdimsr)[i]
            dimlist.append( (old[i] - m1) / (newdimsr[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa
    
    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdimsr[i] )
            dimlist.append( (old[i] - m1) / (newdimsr[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = np.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdimsr) ]
        newcoords = np.mgrid[nslices]

        newcoords_dims = range(np.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs        

        deltas = (np.asarray(old) - m1) / (newdimsr - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None 

#import numpy as np

def robust_sigma(y,zero=None):
    """
    ;FUNCTION  ROBUST_SIGMA,Y, ZERO=REF
    ;
    ;
    ;+
    ; NAME:
    ;	ROBUST_SIGMA  
    ;
    ; PURPOSE:
    ;	Calculate a resistant estimate of the dispersion of a distribution.
    ; EXPLANATION:
    ;	For an uncontaminated distribution, this is identical to the standard
    ;	deviation.
    ;
    ; CALLING SEQUENCE:
    ;	result = ROBUST_SIGMA( Y, [ /ZERO ] )
    ;
    ; INPUT: 
    ;	Y = Vector of quantity for which the dispersion is to be calculated
    ;
    ; OPTIONAL INPUT KEYWORD:
    ;	/ZERO - if set, the dispersion is calculated w.r.t. 0.0 rather than the
    ;		central value of the vector. If Y is a vector of residuals, this
    ;		should be set.
    ;
    ; OUTPUT:
    ;	ROBUST_SIGMA returns the dispersion. In case of failure, returns 
    ;	value of -1.0
    ;
    ; PROCEDURE:
    ;	Use the median absolute deviation as the initial estimate, then weight 
    ;	points using Tukey's Biweight. See, for example, "Understanding Robust
    ;	and Exploratory Data Analysis," by Hoaglin, Mosteller and Tukey, John
    ;	Wiley & Sons, 1983.
    ;
    ; REVSION HISTORY: 
    ;	H. Freudenreich, STX, 8/90
    ;       Replace MED() call with MEDIAN(/EVEN)  W. Landsman   December 2001
      To Python  NZ: 9,2008
    ;
    ;-
    """

    eps = 1.0E-20
    if  (zero != None):
        y0= 0.0
    else:
        y0 = np.median(y,axis=None)

    # First, the median absolute deviation MAD about the median:
    
    mad = np.median( abs(y-y0),axis=None)/0.6745

    # If the MAD=0, try the MEAN absolute deviation:
    if mad < eps: mad = np.average( abs(y-y0) )/.80
    if mad < eps: 
       return 0.0
 

    # Now the biweighted value:
    u   = (y-y0)/(6.*mad)
    uu  = u*u
    q   = np.where(uu <= 1.0)

    if np.size(q) < 3:
       print 'robust_sigma: tHIS DISTRIBUTION IS too weird! rETURNING -1'
       siggma = -1.
       return siggma

    uq=(1-uu[q])
    yq=(y[q]-y0)
    
    arg = yq*yq * uq*uq*uq*uq
    numerator = np.sum(arg)
    n     = np.size(y)

    uuq=uu[q]
    arg=uq * (1.0 - 5*uuq)
    den1  = np.sum( arg )
    siggma = n*numerator/(den1*(den1-1.))
 
    if siggma > 0.0: return np.sqrt(siggma) 
    else: return 0.0



 # $Id: c_correlate.pro,v 1.20 2004/01/21 15:54:48 scottm Exp $
#
# Copyright (c) 1995-2004, Research Systems, Inc.  All rights reserved.
#       Unauthorized reproduction prohibited.
#+
# NAME:
#       C_CORRELATE
#
# PURPOSE:
#       This function computes the cross correlation Pxy(L) or cross
#       covariance Rxy(L) of two sample populations X and Y as a function
#       of the lag (L).
#
# CATEGORY:
#       Statistics.
#
# CALLING SEQUENCE:
#       Result = C_correlate(X, Y, Lag)
#
# INPUTS:
#       X:    An n-element vector of type integer, float or double.
#
#       Y:    An n-element vector of type integer, float or double.
#
#     LAG:    A scalar or n-element vector, in the interval [-(n-2), (n-2)],
#             of type integer that specifies the absolute distance(s) between
#             indexed elements of X.
#
# KEYWORD PARAMETERS:
#       COVARIANCE:    If set to a non-zero value, the sample cross
#                      covariance is computed.
#
#       DOUBLE:        If set to a non-zero value, computations are done in
#                      double precision arithmetic.
#
# EXAMPLE
#       Define two n-element sample populations.
#         x = [3.73, 3.67, 3.77, 3.83, 4.67, 5.87, 6.70, 6.97, 6.40, 5.57]
#         y = [2.31, 2.76, 3.02, 3.13, 3.72, 3.88, 3.97, 4.39, 4.34, 3.95]
#
#       Compute the cross correlation of X and Y for LAG = -5, 0, 1, 5, 6, 7
#         lag = [-5, 0, 1, 5, 6, 7]
#         result = c_correlate(x, y, lag)
#
#       The result should be:
#         [-0.428246, 0.914755, 0.674547, -0.405140, -0.403100, -0.339685]
#
# PROCEDURE:
#       See computational formula published in IDL manual.
#
# REFERENCE:
#       INTRODUCTION TO STATISTICAL TIME SERIES
#       Wayne A. Fuller
#       ISBN 0-471-28715-6
#
# MODIFICATION HISTORY:
#       Written by:  GGS, RSI, October 1994
#       Modified:    GGS, RSI, August 1995
#                    Corrected a condition which excluded the last term of the
#                    time-series.
#                - GGS, RSI, April 1996
#                    Simplified CROSS_COV function. Added DOUBLE keyword.
#                    Modified keyword checking and use of double precision.
#                - W. Biagiotti,  Advanced Testing Technologies
#                Inc., Hauppauge, NY, July 1997, Moved all
#                constant calculations out of main loop for
#                greatly reduced processing time.
#   CT, RSI, September 2002. Further speed improvements, per W. Biagiotti.
#                Now handles large vectors and complex inputs.
#-
def c_correlate(x, y, lag, covariance=None, double=None):

   n_params = 3
   doublein = double
   
   # COMPILE_OPT IDL2
   
   # Compute the sample cross correlation or cross covariance of
   # (Xt, Xt+l) and (Yt, Yt+l) as a function of the lag (l).
   
   # ON_ERROR, 2
   
   x = np.asfarray(x)
   y = np.asfarray(y)
   nx = np.size(x)
   
   if (nx != np.size(y)):
      print("X and Y arrays must have the same number of elements.")
   
   #Check length.
   if (nx < 2):   
       print("X and Y arrays must contain 2 or more elements.")
   
   
       #If the DOUBLE keyword is not set then the internal precision and
       #result are identical to the type of input.
   usedouble = doublein is not None
   tylag = type(lag)
   if usedouble:
     tylag = np.double 
   
   # This will now be in double precision if Double is set.
   xd = x - np.sum(x, dtype=tylag) / nx #Deviations
   yd = y - np.sum(y, dtype=tylag) / nx
   
   nlag = np.size(lag)

   cross = np.zeros(nlag, dtype=tylag)
   
   m = np.absolute(lag)
   for k in range(nlag):
       # Note the reversal of the variables for negative lags.
       cross[k] = (((lag[k] >= 0)) and [np.sum(xd[0:(nx - lag[k] - 1)+1] * yd[lag[k]:])] or [np.sum(yd[0:(nx + lag[k] - 1)+1] * xd[-lag[k]:])])[0]
   
   # Divide by N for covariance, or divide by variance for correlation.
   temp = np.asfarray(cross.copy())
   if covariance is not None:
       cross = temp / nx
   else:
       cross = temp / np.sqrt(np.sum(xd ** 2) * np.sum(yd ** 2))
   del(temp)
   
   return ((usedouble) and [cross] or [np.asfarray(cross)])[0]
   


def gcentroid(img,x,y,fwhm=5,maxgood=None,keepcenter=None):
    """
    Function from IDL routine gcntrd (astrolib):
    http://idlastro.gsfc.nasa.gov/contents.html
    http://idlastro.gsfc.nasa.gov/ftp/pro/idlphot/gcntrd.pro

    Need to run unit test to check if IDL output is the same for a given input data...!!!!!!!!!! 

;+
;  NAME: 
;       GCNTRD
;  PURPOSE:
;       Compute the stellar centroid by Gaussian fits to marginal X,Y, sums 
; EXPLANATION:
;       GCNTRD uses the DAOPHOT "FIND" centroid algorithm by fitting Gaussians
;       to the marginal X,Y distributions.     User can specify bad pixels 
;       (either by using the MAXGOOD keyword or setting them to NaN) to be
;       ignored in the fit.    Pixel values are weighted toward the center to
;       avoid contamination by neighboring stars. 
;
;  CALLING SEQUENCE: 
;       GCNTRD, img, x, y, xcen, ycen, [ fwhm , /SILENT, /DEBUG, MAXGOOD = ,
;                            /KEEPCENTER ]
;
;  INPUTS:     
;       IMG - Two dimensional image array
;       X,Y - Scalar or vector integers giving approximate stellar center
;
;  OPTIONAL INPUT:
;       FWHM - floating scalar; Centroid is computed using a box of half
;               width equal to 1.5 sigma = 0.637* FWHM.  GCNTRD will prompt
;               for FWHM if not supplied
;
;  OUTPUTS:   
;       XCEN - the computed X centroid position, same number of points as X
;       YCEN - computed Y centroid position, same number of points as Y
;
;       Values for XCEN and YCEN will not be computed if the computed
;       centroid falls outside of the box, or if there are too many bad pixels,
;       or if the best-fit Gaussian has a negative height.   If the centroid 
;       cannot be computed, then a  message is displayed (unless /SILENT is 
;       set) and XCEN and YCEN are set to -1.
;
;  OPTIONAL OUTPUT KEYWORDS:
;       MAXGOOD=  Only pixels with values less than MAXGOOD are used to in
;               Gaussian fits to determine the centroid.    For non-integer
;               data, one can also flag bad pixels using NaN values.
;       /SILENT - Normally GCNTRD prints an error message if it is unable
;               to compute the centroid.   Set /SILENT to suppress this.
;       /DEBUG - If this keyword is set, then GCNTRD will display the subarray
;               it is using to compute the centroid.
;       /KeepCenter  By default, GCNTRD finds the maximum pixel in a box 
;              centered on the input X,Y coordinates, and then extracts a new
;              box about this maximum pixel.   Set the /KeepCenter keyword  
;              to skip the step of finding the maximum pixel, and instead use
;              a box centered on the input X,Y coordinates.                          
;  PROCEDURE: 
;       Maximum pixel within distance from input pixel X, Y  determined 
;       from FHWM is found and used as the center of a square, within 
;       which the centroid is computed as the Gaussian least-squares fit
;       to the  marginal sums in the X and Y directions. 
;
;  EXAMPLE:
;       Find the centroid of a star in an image im, with approximate center
;       631, 48.    Assume that bad (saturated) pixels have a value of 4096 or
;       or higher, and that the approximate FWHM is 3 pixels.
;
;       IDL> GCNTRD, IM, 631, 48, XCEN, YCEN, 3, MAXGOOD = 4096       
;  MODIFICATION HISTORY:
;       Written June 2004, W. Landsman  following algorithm used by P. Stetson 
;             in DAOPHOT2.
;-      
       Translated to Python by Sergio Fernandez.
       Modified by NZ Gemini July 2008
    
    """
    
    #Need to run unit test to check if IDL output is the same for a given input data...!!!!!!!!!!
    
    from numpy import fix,round,zeros,arange,exp,transpose,where,mod,float
    from numpy import float64,array,isnan,size,around, isfinite, nansum, nanmax
    from numpy import ones, ravel
    import numpy as np
     
    #print 'Star..ting Gcentroid..'
    sz_image=img.shape
    if len(sz_image) != 2:
        error=True
        print 'Image array (first parameter) must be 2 dimensional!'
        return -1,-1
    xsize=sz_image[0]  # x, y means dimension 0 and dimension 1 for python,
		       # not standard axes!!
    ysize=sz_image[1]
    npts=size(x)
    maxbox=13
    radius=max(0.637*fwhm, 2.001)
    radsq=radius**2
    sigsq=(fwhm/2.35482)**2
    nhalf=int(min(fix(radius),(maxbox - 1)/2))
    nbox=2*nhalf+1     # number of pixels in side of convolution box

    # Turn into a float array of min size 1
    xcen=array(x,dtype=float,ndmin=1)
    ycen=array(y,dtype=float,ndmin=1)

    # Create arrays of size npts
    ix=zeros(npts); iy=zeros(npts)
    if (npts == 1):
        ix[0]=int(round(x))      # central x pixel
        iy[0]=int(round(y))      # central y pixel
    else:
        ix=[int(i) for i in around(x)]
        iy=[int(i) for i in around(y)]

    #Create the Gaussian convolution kernel in variable "g"
    g=zeros((nbox,nbox),dtype=float)
    row2=(arange(float(nbox)) - nhalf)**2
    g[nhalf] = row2
    for i in range(1,nhalf+1):
        temp=row2+i**2
        g[nhalf-i]=temp
        g[nhalf+i]=temp

    g = exp(-0.5*g/sigsq)	#Make c into a Gaussian kernel


# In fitting Gaussians to the marginal sums, pixels will arbitrarily be 
# assigned weights ranging from unity at the corners of the box to 
# NHALF^2 at the center (e.g. if NBOX = 5 or 7, the weights will be
#
#                                 1   2   3   4   3   2   1
#      1   2   3   2   1          2   4   6   8   6   4   2
#      2   4   6   4   2          3   6   9  12   9   6   3
#      3   6   9   6   3          4   8  12  16  12   8   4
#      2   4   6   4   2          3   6   9  12   9   6   3
#      1   2   3   2   1          2   4   6   8   6   4   2
#                                 1   2   3   4   3   2   1
#
# respectively).  This is done to desensitize the derived parameters to 
# possible neighboring, brighter stars.

    #print nbox,arange(10-nhalf,10+nhalf).shape
    x_wt = zeros((nbox,nbox),dtype=float)
    wt = nhalf - abs(arange(float(nbox)) - nhalf) + 1
    for i in range(nbox): 
        x_wt[i] = wt
    y_wt = transpose(x_wt)
    pos= str(x) + ' ' + str(y)

    for i in range(npts):
        if (keepcenter is None):
            if ((ix[i] < nhalf) or ((ix[i]+nhalf) > xsize-1)) or \
               ((iy[i] < nhalf) or ((ix[i]+nhalf) > xsize-1)):
                    print 'WARNING:: position'+pos+' too near edge of image'
                    xcen[i] = -1 ; ycen[i] = -1
		    break

            #OJO!!!!  [y,x]
            d=img[iy[i]-nhalf:iy[i]+nhalf +1,ix[i]-nhalf:ix[i]+nhalf + 1] 
            if  (maxgood is not None):
                if array(maxgood, copy=0).size > 0:
                   ig = where(ravel(d < maxgood))[0]
                   mx = nanmax(d[ig])
            mx = nanmax(ravel(d))       #Maximum pixel value in BIGBOX            
            #How many pixels have maximum value?
            mx_pos = where(ravel(d == mx))[0] # [0] because python returns
                                              # (array([...]),)
            idx= mx_pos % nbox      # X coordinate of Max pixel
            idy= mx_pos / nbox      # Y coordinate of Max pixel

            Nmax=size(mx_pos)
            if Nmax > 1:             # More than 1 pixel at maximum?
                idx=round(idx.sum()/Nmax)
                idy=round(idy.sum()/Nmax)
            else:
                idx=idx[0]
                idy=idy[0]

            xmax = ix[i] - (nhalf) + idx  # X coordinate in original image array
            ymax = iy[i] - (nhalf) + idy  # Y coordinate in original image array
        else:
            xmax = ix[i]
            ymax = iy[i]

        #---------------------------------------------------------------------
        # check *new* center location for range
        # added by Hogg

        if (xmax < nhalf) or ((xmax + nhalf) > (xsize-1)) \
                or (ymax < nhalf) or ((ymax + nhalf) > ysize-1):
                print xmax,nhalf,xsize,ymax,ysize
                print 'position moved too near edge of image'
                #xcen[i] = -1 ; ycen[i] = -1
                return  -1 , -1              

        # Extract  subimage centered on maximum pixel
        ym = round(ymax)
        xm = round(xmax)

        d = img[ym-nhalf : ym+nhalf+1 , xm-nhalf : xm+nhalf+1]

        if maxgood:
            mask=(d<maxgood)*1.     
        else:
            stype= str(img.dtype)
            if ('int' in stype) or ('float' in stype):
                mask = isfinite(d)*1
                mask[-isnan(mask)]=1.
                mask[isnan(mask)]=0.
            else:
                mask = ones([nbox, nbox],dtype=int)

        maskx = (mask.sum(1) > 0)*1
        masky = (mask.sum(0) > 0)*1

        # At least 3 points are needed in the partial sum 
        # to compute the Gaussian

        if maskx.sum() <3 or masky.sum() < 3:
            print 'position has insufficient good points'
            xcen[i]=-1; ycen[i]=-1
            return xcen, ycen

        ywt = y_wt*mask
        xwt = x_wt*mask
        wt1 = wt*maskx
        wt2 = wt*masky

        sd = nansum(d*ywt,0)
        sg = (g*ywt).sum(0)
        sumg = (wt1*sg).sum()
        sumgsq = (wt1*sg*sg).sum()

        sumgd = (wt1*sg*sd).sum()
        sumgx = (wt1*sg).sum()
        sumd = (wt1*sd).sum()
        p = wt1.sum()
        xvec = nhalf - arange(float(nbox)) 
        dgdx = sg*xvec
        sdgdxs = (wt1*dgdx**2).sum()
        sdgdx = (wt1*dgdx).sum()
        sddgdx = (wt1*sd*dgdx).sum()
        sgdgdx = (wt1*sg*dgdx).sum()

        hx = (sumgd - sumg*sumd/p) / (sumgsq - sumg**2/p)

        # HX is the height of the best-fitting marginal Gaussian.
        # If this is not positive then the centroid does not 
        # make sense

        if hx <= 0:
            print '*** Warning: (hx <=0)position cannot be fit by a gaussian',hx
            xcen[i]=-1; ycen[i]=-1
            return xcen,ycen

        skylvl = (sumd - hx*sumg)/p
        dx = (sgdgdx - (sddgdx-sdgdx*(hx*sumg + skylvl*p))) \
                             /(hx*sdgdxs/sigsq)

        #X centroid in original array
        xcen[i] = xmax + dx/(1+abs(dx)) 

        # Now repeat computation for Y centroid

        #sd = (d*xwt).sum(1)
        sd = nansum(d*xwt,1)
        sg = (g*xwt).sum(1)
        sumg = (wt2*sg).sum()
        sumgsq = (wt2*sg*sg).sum()
                               
        sumgd = (wt2*sg*sd).sum()
        sumd = (wt2*sd).sum()
        p = (wt2).sum()
                                
        yvec = nhalf - arange(float(nbox))
        dgdy = sg*yvec
        sdgdys = (wt2*dgdy**2).sum()
        sdgdy = (wt2*dgdy).sum()
        sddgdy = (wt2*sd*dgdy).sum()
        sgdgdy = (wt2*sg*dgdy).sum()
         
        hy = (sumgd - sumg*sumd/p) / (sumgsq - sumg**2/p)

        if (hy <= 0):
            print '*** Warning (hy <=0) position cannot be fit by a gaussian',hy
            xcen[i]=-1; ycen[i]=-1
            return xcen,ycen

        skylvl = (sumd - hy*sumg)/p
        dy = (sgdgdy-(sddgdy-sdgdy*(hy*sumg + skylvl*p)))/(hy*sdgdys/sigsq)
        ycen[i] = ymax + dy/(1+abs(dy))    #X centroid in original array

         #DONE
         #endfor
    return xcen,ycen
