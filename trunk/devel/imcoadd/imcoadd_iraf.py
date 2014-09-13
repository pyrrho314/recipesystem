from pyraf import iraf
from pyraf.iraf import gemini, gemtools
import imcoadd

version = "1.0 (2007 January 25)"

def _imcoadd_iraf (images="", outimage="",
                   sci_ext="SCI", var_ext="VAR", dq_ext="DQ",
                   immasks="DQ", database="imcoadd.dat",
                   threshold=20., fwhm=7.,
                   box=20., alignmethod="wcs", asection="default", xwindow=181,
                   key_inst="INSTRUME", key_camera="CAMERA",
                   key_inport="INPORT", key_xoff="default", key_yoff="default",
                   instscale=1., xsign="default", ysign="default",
                   key_pixscale="PIXSCALE", pixscale=1., dispmag=1.0,
                   rotate=iraf.no, scale=iraf.no, geofitgeom="rscale",
                   order=3, sigfit=2.5, niter=5, coolimit=0.3,
                   geointer="linear", geonxblock=2048, geonyblock=2048,
                   key_ron="RDNOISE", key_gain="GAIN", ron=1., gain=1.,
                   datamin=-1000., key_sat="SATURATI", datamax=50000.,
                   aperture=30., limit=15., key_limit=iraf.yes,
                   lowsigma=7., lowlimit=500., scalenoise=0., growthrad=1,
                   statsec="default", badpixfile="default",
                   fl_inter=iraf.no,
                   fl_refmark=iraf.no, fl_mark=iraf.no, fl_fixpix=iraf.yes,
                   fl_find=iraf.yes,
                   fl_map=iraf.yes,
                   fl_trn=iraf.yes,
                   fl_med=iraf.yes,
                   fl_add=iraf.yes,
                   fl_avg=iraf.no,
                   fl_scale=iraf.yes, fl_overwrite=iraf.no,
                   logfile="imcoadd.log", verbose=iraf.yes):
    """Run the Python version of imcoadd with the cl script parameters."""

    if outimage.strip() == "":
        outimage = None

    if rotate == iraf.yes:
        rotate = True
    else:
        rotate = False
    if scale == iraf.yes:
        scale = True
    else:
        scale = False
    if key_limit == iraf.yes:
        key_limit = True
    else:
        key_limit = False
    if fl_inter == iraf.yes:
        fl_inter = True
    else:
        fl_inter = False
    if fl_refmark == iraf.yes:
        fl_refmark = True
    else:
        fl_refmark = False
    if fl_mark == iraf.yes:
        fl_mark = True
    else:
        fl_mark = False
    if fl_fixpix == iraf.yes:
        fl_fixpix = True
    else:
        fl_fixpix = False
    if fl_find == iraf.yes:
        fl_find = True
    else:
        fl_find = False
    if fl_map == iraf.yes:
        fl_map = True
    else:
        fl_map = False
    if fl_trn == iraf.yes:
        fl_trn = True
    else:
        fl_trn = False
    if fl_med == iraf.yes:
        fl_med = True
    else:
        fl_med = False
    if fl_add == iraf.yes:
        fl_add = True
    else:
        fl_add = False
    if fl_avg == iraf.yes:
        fl_avg = True
    else:
        fl_avg = False
    if fl_scale == iraf.yes:
        fl_scale = True
    else:
        fl_scale = False
    if fl_overwrite == iraf.yes:
        fl_overwrite = True
    else:
        fl_overwrite = False
    if verbose == iraf.yes:
        verbose = True
    else:
        verbose = False

    x = imcoadd.Imcoadd (images, outimage=outimage,
                sci_ext=sci_ext, var_ext=var_ext, dq_ext=dq_ext,
                immasks=immasks, database=database,
                threshold=threshold, fwhm=fwhm,
                key_inst=key_inst, key_camera=key_camera,
                key_inport=key_inport,
                geofitgeom=geofitgeom,
                geointer=geointer,
                geonxblock=geonxblock, geonyblock=geonyblock,
                key_ron=key_ron, key_gain=key_gain, ron=ron, gain=gain,
                datamin=datamin, key_sat=key_sat, datamax=datamax,
                statsec=statsec,
                badpixfile=badpixfile,
                fl_fixpix=fl_fixpix, fl_scale=fl_scale,
                fl_overwrite=fl_overwrite,
                logfile=logfile, verbose=verbose)

    # Derive the transformation using GEOMAP (a different version for
    # each alignmethod).
    if fl_map:
        if alignmethod == "wcs":
            x.wcs_map (box=box, rotate=rotate, scale=scale,
                       order=order, sigfit=sigfit, niter=niter,
                       coolimit=coolimit, fl_inter=fl_inter, fl_mark=fl_mark)

        elif alignmethod == "user":
            x.user_map (box=box, rotate=rotate, scale=scale,
                       order=order, sigfit=sigfit, niter=niter,
                       fl_refmark=fl_refmark, dispmag=dispmag,
                       coolimit=coolimit, fl_inter=fl_inter, fl_mark=fl_mark)

        elif alignmethod == "twodx":
            x.twodx_map (box=box, rotate=rotate, scale=scale,
                       order=order, sigfit=sigfit, niter=niter,
                       asection=asection, xwindow=xwindow,
                       coolimit=coolimit, fl_inter=fl_inter, fl_mark=fl_mark)

        elif alignmethod == "header":
            x.header_map (box=box, rotate=rotate, scale=scale,
                       order=order, sigfit=sigfit, niter=niter,
                       key_xoff=key_xoff, key_yoff=key_yoff,
                       instscale=instscale, xsign=xsign, ysign=ysign,
                       key_pixscale=key_pixscale, pixscale=pixscale,
                       coolimit=coolimit, fl_inter=fl_inter, fl_mark=fl_mark)

    # Derive the median image using IMCOMBINE.
    if fl_med:
        x.median (output=outimage, aperture=aperture)

    # Derive the average image, cleaned for cosmic-ray-events, using IMCOMBINE.
    if fl_add:
        x.add (limit=limit, key_limit=key_limit,
               lowsigma=lowsigma, lowlimit=lowlimit,
               scalenoise=scalenoise, growthrad=growthrad)

    # Derive the average uncleaned image using IMCOMBINE.
    if fl_avg:
        x.average()

    x.close()

_parfile = "/data/bowline1/hodge/b5/gemini/test/python/imcoadd.par"
# _parfile = "gemtools$imcoadd.par"

t_imcoadd = iraf.IrafTaskFactory (taskname="imcoadd",
                                  value=iraf.osfn(_parfile),
                                  pkgname=PkgName, pkgbinary=PkgBinary,
                                  function=_imcoadd_iraf)
