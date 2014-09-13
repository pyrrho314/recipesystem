import numpy as np

from gempy.adlibrary.segmentation import gmos_fplen
from astrodata import Lookups

# Load the standard comments for header keywords that will be updated
# in these functions
keyword_comments = Lookups.get_lookup_table("Gemini/keyword_comments",
                                            "keyword_comments")
def applyApproxWaveCal(ad):
    """
      Calculate an approximate wavelength solution for an ARC spectrum.
      There static tables containing CRVAL and CDELT in the Lookup 
      directory giving the center of the spectrum based on:
      Camera, grating, filter, prism, mask, order

      
    """

    # Camera, grating, filter, prism, mask, order
    # These are keywords in the PHU
    hdrlist = {
       'NIFS':('GRATING', 'FILTER'),
       # For gnirs: the NSCUTSPC is the extension [SCI,(1,2,..)] header
       # keyword corresponding to 'order' number. Get NSCUTSPC
       # in the function 
       'GNIRS':('CAMERA', 'GRATING','FILTER', 'PRISM'),
       # F2: prism is unused
       'F2':('LYOT', 'GRISM', 'FILTER', 'MOSPOS'),
       # NIRI: grating and prism  are unused
       'NIRI':('CAMERA', 'FILTER', 'FPMASK'),
       'GMOS-S':('OBJECT'),   # OBJECT is a dummy name, we use it here
       'GMOS-N':('OBJECT'),   # to not break the design
       }
    # Get header values. The input AD should have been
    # verified to contain a MOS image.

    instr = ad.instrument().as_pytype()
    phu = ad.phu_get_key_value

    kwvalues = [phu(kw) for kw in hdrlist[instr]]
    
    if '_ARC' not in str(ad.types):
        raise ValueError("Input file is not of _ARC type.")

    if ad.is_type('F2'):
        f2_appwave(ad,kwvalues)
    elif ad.is_type('GNIRS'):
        gnirs_appwave(ad,kwvalues)
    elif ad.is_type('GMOS'):
        gmos_appwave(ad)
    elif ad.is_type('NIFS'):
        nifs_appwave(ad,kwvalues)
    elif ad.is_type('NIRI'):
        niri_appwave(ad,kwvalues)
    else:
        raise ValueError("Astrodata type not supported: "+str(ad.types))

def gnirs_appwave(ad,kwvalues):
    

    camera2delta = Lookups.get_lookup_table('Gemini/GNIRS/camera_gnirsspec.py',
                                         'camera2delta')
    camera2wcs = Lookups.get_lookup_table('Gemini/GNIRS/camera_gnirsspec.py',
                                         'camera2wcs')

    kws="%s,%s,%s,%s"%tuple(kwvalues)
    for extn in ['SCI','VAR','DQ']:
        if ad[extn] == None:
            continue
        for xad in ad[extn]:
            # NSCUTSPC is the spectrum order number 
            if xad.header.has_key('NSCUTSPC'):
               # Xdispersed data
               ss = kws+','+str(xad.header['NSCUTSPC'])
               crval,cdelt,t1,t2,t3 = camera2wcs[ss] 
            else:
               # Long Slit data
               ss = kws
               cdelt,resol,hi,lo = camera2delta[ss] 
               crval = xad.phu.header['waveleng']
            update_wcs(xad,crval,cdelt)
 
def f2_appwave(ad,kwvalues):
    
    camera2wcs = Lookups.get_lookup_table('Gemini/F2/camera_f2spec.py',
                                         'camera2wcs')
    ss="%s,%s,%s,%s"%tuple(kwvalues)
    for extn in ['SCI','VAR','DQ']:
        if ad[extn] == None:
            continue
        for xad in ad[extn]:
            crval,cdelt = camera2wcs[ss] 
            update_wcs(xad,crval,cdelt)

def niri_appwave(ad,kwvalues):

    camera2wcs = Lookups.get_lookup_table('Gemini/NIRI/camera_nirispec.py',
                                         'camera2wcs')

    ss="%s,%s,%s"%tuple(kwvalues)
    for extn in ['SCI','VAR','DQ']:
        if ad[extn] == None:
            continue
        for xad in ad[extn]:
            crval,cdelt = camera2wcs[ss] 
            update_wcs(xad,crval,cdelt)
    
def nifs_appwave(ad,kwvalues):

    camera2wcs = Lookups.get_lookup_table('Gemini/NIFS/camera_nifsspec.py',
                                         'camera2wcs')

    ss="%s,%s"%tuple(kwvalues)
    for extn in ['SCI','VAR','DQ']:
        if ad[extn] == None:
            continue
        for xad in ad[extn]:
            crval,cdelt = camera2wcs[ss] 
            update_wcs(xad,crval,cdelt)

def gmos_appwave(ad):
    """
      Similar to gsappwave. Works for GMOS LONG SLIT,
      and MOS. IFU still to be implemented.
      MOS mode: relies on the phu header keyword GSCUT
                and hdu header keyword REFPIX1
      nz Nov 2013
    """

    # For non-IFU modes, calculate reference pixel and length
    # as previously, based on whether the spectra are cut out:

    nypix, nxpix = ad['SCI',1].data.shape

    gscut = ad.phu.header.has_key('GSCUT')
    for extn in ['SCI','VAR','DQ']:
        if ad[extn] == None:
            continue
        for xad in ad[extn]:
            wave1,wave2,wavoffset,nmppx,a,cwave,l_yoff = gmos_fplen(xad)
            if gscut:
                # MOS data
                refpix = xad.header['refpix1']
                #speclen = int(round((wave2-wave1)/nmppx))
                #refpix = float(speclen) - (cwave-wavoffset-wave1)/nmppx
            else:
                speclen = nxpix
                refpix = float(nxpix)/2. + wavoffset/nmppx
    
            crval = cwave*10.
            crpix = refpix
            cdelt = -10.*nmppx
            update_wcs(xad,crval,cdelt,crpix)

def update_wcs(ext,crval,cdelt, crpix=None):

    ny,nx = ext.data.shape
    axis = ext.dispersion_axis()-1
    if crpix != None:
        crpix1,crpix2 = [(crpix, 1.),(1., crpix)][axis]
    else:
        crpix1,crpix2 = [(nx/2, 1.),(1., ny/2)][axis]
         
    crval1,crval2 = [(crval, 1.),(1., crval)][axis]
    cd11,cd22 =     [(cdelt, 1.),(1., cdelt)][axis]
    ext.set_key_value("CTYPE1","LINEAR",
            comment=keyword_comments["CTYPE1"])
    ext.set_key_value("CTYPE2","LINEAR",
            comment=keyword_comments["CTYPE2"])
    ext.set_key_value( "CRVAL1", crval1,
            comment=keyword_comments["CRVAL1"])
    ext.set_key_value( "CRVAL2", crval2,
            comment=keyword_comments["CRVAL2"])
    ext.set_key_value( "CRPIX1", crpix1,
            comment=keyword_comments["CRPIX1"])
    ext.set_key_value( "CRPIX2", crpix2,
            comment=keyword_comments["CRPIX2"])
    ext.set_key_value( "CD1_1", cd11,
            comment=keyword_comments["CD1_1"])
    ext.set_key_value( "CD2_2", cd22,
            comment=keyword_comments["CD2_2"])
    ext.set_key_value( "CD2_1",0.,
            comment=keyword_comments["CD2_1"])
    ext.set_key_value( "CD1_2",0.,
            comment=keyword_comments["CD1_2"])

