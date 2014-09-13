

# Import core Python modules
import time
import os

# Import scientific modules
import pyfits
try:
    import stsci.numdisplay as numdisplay
except ImportError:
    import numdisplay

import numpy as np

# Import IQ modules
import iqUtil
import fit

from astrodata.adutils import gemutil
from astrodata.adutils import mefutil
from astrodata.adutils import filesystem
from astrodata.adutils import gemLog

import iqtool
# This is older imports to have it working with util stuff
# in pygem.
#from iq import fit
#from iq import util
#from pygem import gemutil, mefutil, irafutil, filesystem

from math import sqrt

# Global logger object instantiated in gemiq if verbose=True
log = None
#---------------------------------------------------------------------------
def gemiq(image, outFile='default', function='both', verbose=True,\
          residuals=False, display=True, \
          interactive=False, rawpath='.', prefix='auto', \
          observatory='gemini-north', clip=True, \
          sigma=2.3, pymark=True, niters=4, boxSize=2., mosaic=False,
          debug=False, garbageStat=False, qa=False):
    """get psf measurements of stars

    @param image: Input filename or number if using today's date
    @type image: String or int
    
    @param outFile='default': output file with fit parameters, if default uses image
                               filename+.log
    @type outFile: string
    
    @param function='both': function to fit psf. Currently supported values:
                     moffat/gauss/both
                     where both fits moffat and gaussian functions to data
    @type function: string
    
    @param verbose=True: print information to screen, includes printing psf values
    @type verbose: Boolean
    
    @param residuals=False: create and display residual subtracted images
    @type residuals: Boolean
    
    @param display=True: display images
    @type display: Boolean
    
    @param pymark=True: mark images with X's for all iraf.daofind outputs and
        circles with radius FWHM of fits
    @type pymark: Boolean
    
    @param interactive=False: if interactive = False, allow iraf.daofind to select
        star centers. otherwise, if interactive = True, allow user to select
        object centers using 'space' or 'a'
    @type interactive: Boolean
    
    @param rawpath='.': path for images if not current directory
    @type rawpath: string
    
    @param prefix='auto': if image is a number, then prefix will be of form 'N20080915'
    @type prefix: string
    
    @param observatory='gemini-north': decides whether first letter of filename is
        'N' or 'S'
    @type observatory: string
    
    @param clip=True: sigma clip measured FWHMarcsec and elliticities to remove
        outliers that snuck 
        through. Produces a mean seeing and ellipticity for enire image with standard
        deviation
    @type clip: Boolean
    
    @param sigma=2.5: if clip=True, sigma is the minimum number of standard deviations
        away from the mean that an object is clipped from the final catalogue
    @type sigma: float or int
    
    @param niters=2: integer number of times sigma clipping is iterated through
        parameter list
    @type niter: int
    
    @param boxSize=2.: aperture size around object x,y to make postage stamp for
        the fitting function
    @type boxSize: float or int
    
    @param debug=False: very verbose, print all objects removed from catalogue
        because they were close to detector edge or too close to neighbors.
    @type debug: Boolean
    
    @param qa: A flag to use a grid of sub-windows for detecting the sources in 
               the image frames, rather than the entire frame all at once.
    @type qa: Boolean
    """

    imagelist = []
    
    if not verbose:
        debug = False
    else:
        # instantiate the logger object and put into the global variable 
        global log
        if log==None: 
            log = gemLog.getGeminiLog()
        if debug:   
            # update logger object's file handler level to debug
            log.changeLevels(logLevel=log.loglevel(), debug=True)
        
    
    try:
        image=int(image)
        imagelist.append(image)
    except: 
        if image[0] == '@': # image is list add comma delineated --> glob list
            imagelist = open(image[1:len(image)],'r')
        else:
            imagelist.append(image)

    for image in imagelist:
        #remove the \n from end of each image name if it is there
        try:
            image=int(image)
        except:
            if image[-1]=='\n': image = image[0:len(image)-1]    

        if verbose:
        	log.status("Measuring IQ for image "+str(image))
        
        # imagenorawpath=os.path.basename(image)
        filename, imagenorawpath = gemutil.imageName(image, rawpath, prefix=prefix, 
                                                  observatory=observatory, 
                                                  verbose=verbose)

        if outFile == 'default': 
            outFile = imagenorawpath[0:len(imagenorawpath)-5]+'.dat'
            
        if verbose:
            log.status('Important values looked up and calculated will be stored'+
                       ' in '+outFile)
            
        # Open the output parameter file
        paroutfile=open(outFile, 'w')

        paroutfile.write('\n# filename: '+filename+'\n')
        #open image array
        hdulist = pyfits.open(filename)
        instrument = mefutil.getkey('instrume', filename)
        pixelscale = mefutil.getkey('pixscale', filename)

        ############################## GMOS specific #############################
        if instrument == 'GMOS-N' or instrument == 'GMOS-S':
            
            if instrument == 'GMOS-N': dpixelscale = 0.0727
            if instrument == 'GMOS-S': dpixelscale = 0.073
        
            gmoskeys = ('ELEVATIO','AZIMUTH','UT','DATE','TAMBIENT','WINDSPEE',
                        'WINDDIRE','EXPTIME',
                        'FILTER1','FILTER2','FILTER3', 'CRPA','OBSTYPE','OBSCLASS',
                        'OIWFS_ST',
                        'HUMIDITY','PRESSURE','DTAZEN','PA','PIXSCALE','NCCDS',
                        'RA', 'DEC',
                        'RAOFFSET', 'DECOFFSE', 'OIARA', 'OIADEC')
            
            keylist= mefutil.getkeys(gmoskeys,filename)
            
            #$$$$$$$$$$$ This next line is checking the PHU for EXTNAME which DOES NOT prove it is RAW!!!!!!!!!!!!!$$$$
            #$$$$$$$$$$$ EXTNAME is only found in the image extensions, not the default zero (PHU) extension that getkey uses!!!
            raw = (mefutil.getkey('EXTNAME',filename, extension=1) == 'not found')

            gmoskeysScihdr = ('CRVAL1','CRVAL2')
            if mosaic:
                #$$$$$$$$$ changed extension to 1 from zero, as that key is only in the image extensions
                #$$$$$$$$$ and not the PHU
                CRVAL1, CRVAL2 = mefutil.getkeys(gmoskeysScihdr,filename, extension=1)
            else:
                # Use extension 2 (ie middle chip) if no mosaic as it provides best seeing values
                CRVAL1, CRVAL2 = mefutil.getkeys(gmoskeysScihdr,filename, extension=2)
            #automatically generate this gmoskeydict with a loop
            gmoskeydict = {'ELEVATIO':keylist[0],'AZIMUTH':keylist[1],'UT':keylist[2],
                           'DATE':keylist[3], 'TAMBIENT':keylist[4],
                           'WINDSPEE':keylist[5],'WINDDIRE':keylist[6],
                           'EXPTIME':keylist[7],'FILTER1':keylist[8],'FILTER2':keylist[9],
                           'FILTER3':keylist[10], 'CRPA':keylist[11],'OBSTYPE':keylist[12],
                           'OBSCLASS':keylist[13],'OIWFS_ST':keylist[14],
                           'HUMIDITY':keylist[15],
                           'PRESSURE':keylist[16],'DTAZEN':keylist[17],'PA':keylist[18],
                           'PIXSCALE':keylist[19],'NCCDS':keylist[20],
                           'RA':keylist[21], 'DEC':keylist[22],
                           'RAOFFSET':keylist[23], 'DECOFFSE':keylist[24],
                           'OIARA':keylist[25],
                           'OIADEC':keylist[26], 'CRVAL1':CRVAL1, 'CRVAL2':CRVAL2}

            
            
            for key in gmoskeydict:
                paroutfile.write('# '+key+' = '+str(gmoskeydict[key])+'\n')

            if gmoskeydict['NCCDS'] == 'not found': 
                NCCDS=3
            else: 
                NCCDS = gmoskeydict['NCCDS']
            
            n=1
            allgfit = []
            allmfit = []
            # Loop through GMOS CCDS
            while n <= NCCDS:
                if n==2 and mosaic:
                    break
                scidata = hdulist[n].data # don't open twice
                scihdr = hdulist[n].header
                ccdsum = scihdr['CCDSUM']
                ccdint = int(ccdsum[0])
                if pixelscale == 'not found':
                    if ccdsum != 'not found':
                        #ccdint = int(ccdsum[0])
                        pixelscale = dpixelscale * ccdint
                    else:
                        # Force CCD binning to 2x2
                        ccdint = 2
                        pixelscale = dpixelscale * ccdint
                        log.warning('Could not find PIXELSCALE \
                            in PHU or CCDSUM in science header,\
                            Using pixelscale '+str(pixelscale)) 
                # Setting the min and max X axis values for checking the IQ in
                if n==1:
                    # check for stamp
                    detro1xs = hdulist[0].header['DETRO1XS']
                    # If the region is very small, don't perform IQ check
                    if detro1xs == 300 / ccdint:
                        xmin = None; xmax = None
                    # Else, load up max and min X axis vals for region for 
                    # if 1x1 or 2x2 CCD binning
                    else:
                        if ccdint == 1: xmin =890; xmax=2090
                        else: xmin = 445; xmax = 1045
                # If currently looking at CCD 2 (middle CCD), don't check IQ
                # for this CCD.  $$$$$ WHY ????????
                if n==2:
                    imageSigma = 'default'
                    xmin = None
                    xmax = None
                # If CCD 3, load up max and min X axis vals for region for 
                # if 1x1 or 2x2 CCD binning
                if n==3:
                    if ccdint == 1: xmin = 31; xmax = 1200
                    else: xmin = 31; xmax=600

                # Set the stars (those pixels with values higher than 
                # sigma + mean) and lows (those pixels with values lower than
                # sigma - mean) to the image mean
                imageSigma = iqUtil.starMask(scidata[:,xmin:xmax]).std()
                
                if debug: 
                    log.debug('GETIQ is currently working on GMOS CCD: '+str(n)
                              +' if image '+filename)
                
                # Goes in interactive using pyexam if requested
                if interactive: 
                    pyexam(scidata, function, pixelscale, frame=n,
                                     outFile=paroutfile, verbose=verbose,
                                     pymark=pymark, residuals=residuals,
                                     clip=clip, debug=debug, niters=niters,
                                     boxSize=boxSize)
                else:
                    if raw:
                        gAllstars, mAllstars = pyiq(filename+'['+str(n)+']',
                            scidata, function, paroutfile, pixelscale,
                            frame=n, pverbose=verbose, pymark=pymark,
                            residuals=residuals, clip=False, sigma=sigma,
                            niters=niters, display=display,
                            imageSigma=imageSigma, boxSize=boxSize,
                            debug=debug, xmin=xmin,xmax=xmax, qa=qa)
                    else:
                        gAllstars, mAllstars = pyiq(filename+'[sci,'+str(n)+']',
                            scidata, function, paroutfile, pixelscale,
                            frame=n, pverbose=verbose, pymark=pymark,
                            residuals=residuals, clip=False, sigma=sigma,
                            niters=niters, display=display,
                            imageSigma=imageSigma, boxSize=boxSize,
                            debug=debug, xmin=xmin, xmax=xmax, qa=qa)
                   
                    if gAllstars: 
                        allgfit.append(gAllstars)
                    if mAllstars: 
                        allmfit.append(mAllstars)
                n+=1
            iqdata = []
            if clip: # Clip values from all three CCDS 
                if verbose: 
                    log.fullinfo('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            
                if allgfit:
                    allgfit2 = []
                    # allgfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for ccd in allgfit:
                        for star in ccd:
                            allgfit2.append(star)
                        
    	            if verbose: 
                        log.status('performing Gaussian clipping of all '+
                                   'GMOS ccds')
                            
                    allgfit2, gEllMean, gEllSigma, gFWHMMean, gFWHMSigma = \
    	        	   iqUtil.sigmaClip(allgfit2, paroutfile,
                                                    sigma, verbose, niters, garbageStat=garbageStat)
                    for star in allgfit2: #only write the clipped objects 
                        iqUtil.writePars(star, paroutfile, 'gaussian')
                    
                    g1 = 'All CCDS Total Gauss Ellipticity:\
                    '+str(gEllMean)+' +/- '+str(gEllSigma)
                    g2 = 'ALL CCDS Total Gauss FWHM (arcsec):\
                    '+str(gFWHMMean)+' +/- '+str(gFWHMSigma)+'\n'
                    
                    if verbose:
                        log.stdinfo(g1+'\n'+g2)
                    paroutfile.write(g1+'\n'+g2)
                    iqdata.append((gEllMean,gEllSigma,gFWHMMean,gFWHMSigma))
                if allmfit:
                    allmfit2 = []
                    # allmfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for ccd in allmfit:
                        for star in ccd:
                            allmfit2.append(star)
    	            if verbose: 
                        log.status('Performing Moffat clipping of all GMOS ccds'
                                   )
                    
                    allmfit2, mEllMean, mEllSigma, mFWHMMean, mFWHMSigma = \
                    iqUtil.sigmaClip(allmfit2, paroutfile, sigma,
                                             debug, niters, garbageStat=garbageStat)
                    for star in allmfit2: #only write the clipped objects
                         # to the catalogue
                         iqUtil.writePars(star, paroutfile, 'moffat')
                                     
                    m1 = 'All CCDS Total Moffat Ellipticity:\
                     '+str(mEllMean)+' +/- '+str(mEllSigma)
                    m2 = 'ALL CCDS Total Moffat FWHM (arcsec):\
                     '+str(mFWHMMean)+' +/- '+str(mFWHMSigma)+'\n'
                    if verbose:
                        log.stdinfo(m1+'\n'+m2)
                    paroutfile.write(m1+'\n'+m2)
                    iqdata.append((mEllMean,mEllSigma,mFWHMMean,mFWHMSigma))
            else:
                if allgfit:
                    for ccd in allgfit:
                        for star in ccd:
                            iqUtil.writePars(star, paroutfile, 'gaussian')
                if allmfit:
                    for ccd in allmfit:
                        for star in ccd:
                            iqUtil.writePars(star, paroutfile, 'moffat')

    ############################# NIRI SPECIFIC #############################
                            
        elif instrument == 'NIRI':
            #this is where niri stuff will go
            if verbose: 
                log.fullinfo('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
            if verbose: 
                log.status('Filename= '+filename)
            scidata = hdulist[1].data
            # the following has been re-added from Jen's previous version of IQTool:
            nirikeys = ('CD1_1','CD1_2','CD2_1','CD2_2', 'AIRMASS')
            keylist=mefutil.getkeys(nirikeys,filename)
            nirikeydict = {'CD1_1':keylist[0],'CD1_2':keylist[1],'CD2_1':keylist[2],'CD2_2':keylist[3], 'AIRMASS':keylist[4]}
            for key in nirikeydict:
                paroutfile.write('# '+key+' = '+str(nirikeydict[key])+'\n')
            if pixelscale == 'not found':
                pixelscale = 3600*(sqrt(keylist[0]**2+keylist[1]**2)+sqrt(keylist[2]**2+keylist[3]**2))/2
            imageSigma = 1.8 * iqUtil.starMask(scidata[:,:]).std()
            allgfit = []
            allmfit = []
            gAllstars, mAllstars = pyiq(filename+'[1]', scidata,
                function, paroutfile, pixelscale, frame=1, pverbose=verbose,
                pymark=pymark, residuals=residuals, clip=clip, sigma=sigma,
                niters=niters, display=display, imageSigma=imageSigma,
                boxSize=boxSize,debug=debug, xmin=None, xmax=None, ymin=None, ymax=None, qa=qa)
            if gAllstars: allgfit.append(gAllstars)
            if mAllstars: allmfit.append(mAllstars)

            iqdata = []
            if clip:
                if verbose: 
                    log.fullinfo('\n')
                if allgfit:
                    allgfit2 = []
                    # allgfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for star in allgfit[0]:
                        allgfit2.append(star)
                    if verbose: 
                        log.status('Performing Gaussian clipping')
                    allgfit2, gEllMean, gEllSigma, gFWHMMean, gFWHMSigma = \
                     iqUtil.sigmaClip(allfgit2, paroutfile, sigma, debug, niters, garbageStat=garbageStat)
                    iqdata.append( (gEllMean,gEllSigma,gFWHMMean,gFWHMSigma) )
                if allmfit:
                    allmfit2 = []
                    # allmfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for star in allmfit[0]:
                        allmfit2.append(star)
                    if verbose: 
                        log.status('Performing Moffat clipping')
                    allmfit2, mEllMean, mEllSigma, mFWHMMean, mFWHMSigma = \
                     iqUtil.sigmaClip(allmfit2, paroutfile, sigma, debug, niters, garbageStat=garbageStat)
                    iqdata.append( (mEllMean,mEllSigma,mFWHMMean,mFWHMSigma) )

            else:
                if allgfit:
                    for star in allgfit[0]:
                        iqUtil.writePars(star, paroutfile, 'gaussian')
                if allmfit:
                    for star in allmfit[0]:
                        iqUtil.writePars(star, paroutfile, 'moffat')


        elif instrument == 'GNIRS':
            if verbose: 
                log.fullinfo('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
            if verbose: 
                log.status('Filename= '+filename)
            scidata = hdulist[1].data
            gnirskeys = ('CD1_1','CD1_2','CD2_1','CD2_2', 'AIRMASS')
            keylist=mefutil.getkeys(gnirskeys,filename)
            gnirskeydict = {'CD1_1':keylist[0],'CD1_2':keylist[1],'CD2_1':keylist[2],'CD2_2':keylist[3], 'AIRMASS':keylist[4]}
            for key in gnirskeydict:
                paroutfile.write('# '+key+' = '+str(gnirskeydict[key])+'\n')
            if pixelscale == 'not found':
                pixelscale = 3600*(sqrt(keylist[0]**2+keylist[1]**2)+sqrt(keylist[2]**2+keylist[3]**2))/2
            # only search in (approximate) position of keyhole
            camera = mefutil.getkeys(['CAMERA'],filename)[0]
            if 'Short' in camera:
				xmin = 300
				xmax = 700
				ymin = 400
				ymax = 500
            elif 'Long' in camera:
				xmin = 200
				xmax = 800
				ymin = 300
				ymax = 600
            # this helps detect the right number of objects in low S/N images
            imageSigma = scidata[ymin:ymax,xmin:xmax].std() * 0.3
            allgfit = []
            allmfit = []
            gAllstars, mAllstars = pyiq(filename+'[1]', scidata,
                function, paroutfile, pixelscale, frame=1, pverbose=verbose,
                pymark=pymark, residuals=residuals, clip=clip, sigma=sigma,
                niters=niters, display=display, imageSigma=imageSigma,
                boxSize=boxSize,debug=debug, xmin=None, xmax=None, ymin=ymin, ymax=ymax, qa=qa)
            if gAllstars: allgfit.append(gAllstars)
            if mAllstars: allmfit.append(mAllstars)

            iqdata = []
            if clip:
                if verbose: 
                    log.fullinfo('\n')
                if allgfit:
                    allgfit2 = []
                    # allgfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for star in allgfit[0]:
                        allgfit2.append(star)
                    if verbose: 
                        log.status('Performing Gaussian clipping')
                    allgfit2, gEllMean, gEllSigma, gFWHMMean, gFWHMSigma = \
                     iqUtil.sigmaClip(allfgit2, paroutfile, sigma, debug, niters, garbageStat=garbageStat)
                    iqdata.append( (gEllMean,gEllSigma,gFWHMMean,gFWHMSigma) )
                if allmfit:
                    allmfit2 = []
                    # allmfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for star in allmfit[0]:
                        allmfit2.append(star)
                    if verbose: 
                        log.status('Performing Moffat clipping')
                    allmfit2, mEllMean, mEllSigma, mFWHMMean, mFWHMSigma = \
                     iqUtil.sigmaClip(allmfit2, paroutfile, sigma, debug, niters, garbageStat=garbageStat)
                    iqdata.append( (mEllMean,mEllSigma,mFWHMMean,mFWHMSigma) )

            else:
                if allgfit:
                    for star in allgfit[0]:
                        iqUtil.writePars(star, paroutfile, 'gaussian')
                if allmfit:
                    for star in allmfit[0]:
                        iqUtil.writePars(star, paroutfile, 'moffat')


        elif instrument == 'NIFS':
            if verbose: 
                log.fullinfo('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
            if verbose: 
                log.status('Filename= '+filename)
            scidata = hdulist[1].data
            nifskeys = ('CD1_1','CD1_2','CD2_1','CD2_2', 'AIRMASS')
            keylist=mefutil.getkeys(nifskeys,filename)
            nifskeydict = {'CD1_1':keylist[0],'CD1_2':keylist[1],'CD2_1':keylist[2],'CD2_2':keylist[3], 'AIRMASS':keylist[4]}
            for key in nifskeydict:
                paroutfile.write('# '+key+' = '+str(nifskeydict[key])+'\n')
            # pixel scale for reconstructed NIFS frames:
            pixelscale = 0.02
            imageSigma = 1.8 * iqUtil.starMask(scidata[:,:]).std()
            allgfit = []
            allmfit = []
            gAllstars, mAllstars = pyiq(filename+'[1]', scidata,
                function, paroutfile, pixelscale, frame=1, pverbose=verbose,
                pymark=pymark, residuals=residuals, clip=clip, sigma=sigma,
                niters=niters, display=display, imageSigma=imageSigma,
                boxSize=boxSize,debug=debug, xmin=None, xmax=None, ymin=None, ymax=None, qa=qa)
            if gAllstars: allgfit.append(gAllstars)
            if mAllstars: allmfit.append(mAllstars)

            iqdata = []
            if clip:
                if verbose: 
                    log.fullinfo('\n')
                if allgfit:
                    allgfit2 = []
                    # allgfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for star in allgfit[0]:
                        allgfit2.append(star)
                    if verbose: 
                        log.status('Performing Gaussian clipping')
                    allgfit2, gEllMean, gEllSigma, gFWHMMean, gFWHMSigma = \
                     iqUtil.sigmaClip(allfgit2, paroutfile, sigma, debug, niters, garbageStat=garbageStat)
                    iqdata.append( (gEllMean,gEllSigma,gFWHMMean,gFWHMSigma) )
                if allmfit:
                    allmfit2 = []
                    # allmfit is a tuple of three tuples each containing
                    # a tuple of dictionaries       
                    for star in allmfit[0]:
                        allmfit2.append(star)
                    if verbose: 
                        log.status('Performing Moffat clipping')
                    allmfit2, mEllMean, mEllSigma, mFWHMMean, mFWHMSigma = \
                     iqUtil.sigmaClip(allmfit2, paroutfile, sigma, debug, niters, garbageStat=garbageStat)
                    iqdata.append( (mEllMean,mEllSigma,mFWHMMean,mFWHMSigma) )

            else:
                if allgfit:
                    for star in allgfit[0]:
                        iqUtil.writePars(star, paroutfile, 'gaussian')
                if allmfit:
                    for star in allmfit[0]:
                        iqUtil.writePars(star, paroutfile, 'moffat')


        else:
            #this is where other instrument stuff will go
            log.critical('The instrument in image '+filename+' is not GMOS,'+
                         ' NIRI, or GNIRS')
            pixelscale=1.
            imageSigma='default'
            scidata = hdulist[1].data
            gAllstars, mAllstars = pyiq(filename+'[1]', scidata,
                function, paroutfile, pixelscale, frame=1, pverbose=verbose,
                pymark=pymark, residuals=residuals, clip=clip, sigma=sigma,
                niters=niters, display=display, imageSigma=imageSigma,
                boxSize=boxSize,debug=debug, qa=qa)


        paroutfile.close()
        return iqdata
#---------------------------------------------------------------------------


def pyexam(scidata, function='both', pixelscale=1, frame=1, \
           outFile='testout.txt', verbose=True, pymark=True, \
           residuals=False, clip=True, sigma=3, boxSize=9., \
           debug=False, niters=4.):
        import pylab

        from iqtool.gemplotlib import overlay
        
        # Get box size around each object
        apertureSize = fit.getApSize(pixelscale, boxSize)
 
        outstring = '%10s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n'%('function', 'Cx', 'Cy', 'Bg', 'Peak', 'Wx', 'Wy', 'CooX', 'CooY', 'Ellipticity', 'FWHM_pix','FWHM_arcsec', 'Theta', 'PA_deg', ' FWHMx', 'FWHMy', 'Beta')

        outFile.write(outstring)

        print 'PYEXAM - Please select objects using spacebar, and press \'q\' when finished'
        print '         pressing k, j, r, or  e will make plots along with fits'
        keystroke = ''
        numdisplay.display(scidata, z1=scidata.mean()-scidata.std(),
                           z2=scidata.mean()+scidata.std(), frame=frame)     
        if verbose:
            print '%10s%12s%12s%12s%12s%12s%12s'%('function', 'Coox', 'CooY', 'FWHMpix',
                                                  'FWHMarcsec', 'Ellipticity', 'PAdeg')
      
        gAllstars = []
        mAllstars = []
        done = 'no'
        while done == 'no':
            # Get x,y
            
            cursorInput = numdisplay.readcursor(sample=0)
            components = cursorInput.split()
            xPos = float(components[0])
            yPos = float(components[1]) 
            keystroke = components[3]
            option = ''
            if keystroke == '\\040' or keystroke == 'a' or keystroke == 'k' or keystroke == 'j' or keystroke == 'e' or keystroke == 'r':
                print 'PYEXAM - X:',xPos,'Y:',yPos
                gfitArray=None
                mfitArray=None
                gfxn=None
                mfxn=None
                
                positionCoords = (xPos,yPos)

                stampArray, stampPars = fit.makeStamp(scidata, positionCoords,
                                          apertureSize, outFile, debug=debug)
                
                if stampArray == None: return None, None

                gfitPars, gReturnModel, mfitPars, mReturnModel = \
                    fit.fitprofile(stampArray, function, positionCoords, \
                    outFile, pixelscale, stampPars, debug=debug, frame=frame)

                imageDim = list(scidata.shape) # Gives [y,x] dimensions
                #imageDim.reverse() #give [x,y] dimensions
                
                if gfitPars != None:
                    clipobsPars = ('gauss', gfitPars['CooX'], gfitPars['CooY'],
                                   gfitPars['FWHMpix'],
                                   gfitPars['FWHMarcsec'], gfitPars['Ellip'],
                                   gfitPars['PAdeg'])
                    iqUtil.printPars(clipobsPars, verbose)
                    gAllstars.append(gfitPars)

                    if pymark:
                        overlay.circle(x=gfitPars['CooX']-1,y=gfitPars['CooY']-1, \
                                       radius=gfitPars['FWHMpix'],frame=frame,\
                                       color=204)

                    gfxn = gReturnModel((gfitPars['Bg'], gfitPars['Peak'], \
                                           gfitPars['Cy'], \
                                           gfitPars['Cx'], gfitPars['Wy'], \
                                           gfitPars['Wx'], gfitPars['Theta']+90))

                    gfitArray = np.zeros(imageDim)
                    gfitArray[0:stampPars[3],0:stampPars[1]] = gfxn(*np.indices((stampPars[3], stampPars[1])))
                    
                 
                if mfitPars != None:
                    clipobsPars = ('moffat', mfitPars['CooX'],mfitPars['CooY'],
                                   mfitPars['FWHMpix'],
                                   mfitPars['FWHMarcsec'], mfitPars['Ellip'],
                                   mfitPars['PAdeg'])
                      
                    iqUtil.printPars(clipobsPars, verbose)
                    print mfitPars['CooX']
                    print mfitPars['CooY']
                    print mfitPars['Cx']
                    print mfitPars['Cy']

                    
                    mAllstars.append(mfitPars)
                    if pymark:
                        overlay.circle(x=mfitPars['CooX'],
                                       y=mfitPars['CooY'], frame=frame,
                                       radius=mfitPars['FWHMpix'], color=212)

                    mfxn = mReturnModel((mfitPars['Bg'], mfitPars['Peak'], \
                                        mfitPars['Cy'],\
                                        mfitPars['Cx'], mfitPars['Wy'], \
                                        mfitPars['Wx'], mfitPars['Theta']+90, \
                                        mfitPars['Beta']))
                    mfitArray = np.zeros(imageDim)
                    mfitArray[0:stampPars[3],0:stampPars[1]] = mfxn(*np.indices((stampPars[3], stampPars[1])))

                                   
   
                if keystroke == 'k':
                    print "PYEXAM - Yslice Plotting"
                    pylab.clf()
                    pylab.plot(stampArray.max(0), 'ro')
                    if mfitArray != None:
                        pylab.plot(mfitArray.max(0), 'b', label='moffat fit FWHM='+str(mfitPars['FWHMarcsec']))
                        pylab.axis([0,stampPars[1]-stampPars[0],mfitPars['Bg'],1.3*(mfitPars['Bg']+mfitPars['Peak'])])
                    if gfitArray != None:
                        pylab.plot(gfitArray.max(0), 'g',label='gauss fit FWHM='+str(gfitPars['FWHMarcsec']))
                        pylab.axis([0,stampPars[1]-stampPars[0],gfitPars['Bg'],1.3*(gfitPars['Bg']+gfitPars['Peak'])])
                    pylab.legend()

                if keystroke == 'j':
                    print "PYEXAM - Xslice Plotting"
                    pylab.clf()
                    pylab.plot(stampArray.max(1), 'ro')
                    if mfitArray != None:
                        pylab.plot(mfitArray.max(1), 'b', label='moffat fit FWHM='+str(mfitPars['FWHMarcsec']))
                        pylab.axis([0,stampPars[1]-stampPars[0],mfitPars['Bg'],1.3*(mfitPars['Bg']+mfitPars['Peak'])])
                    if gfitArray != None:
                        pylab.plot(gfitArray.max(1), 'g',label='gauss fit FWHM='+str(gfitPars['FWHMarcsec']))
                        pylab.axis([0,stampPars[1]-stampPars[0],gfitPars['Bg'],1.3*(gfitPars['Bg']+gfitPars['Peak'])])
                    pylab.legend()

                if keystroke == 'r':
                    print "PYEXAM - radial Plotting doesnt work yet!"
                    print "new version"
                    #rad1 = np.sqrt(stampArray.max(1)**2 + stampArray.max(0)**2)
                    pylab.clf()
                    pylab.plot(stampArray.max(1), 'ro')
                    if mfitArray != None:
                        pylab.plot(mfitArray.max(1), 'b', label='moffat fit FWHM='+str(mfitPars['FWHMarcsec']))
                        pylab.axis([0,stampPars[1]-stampPars[0],mfitPars['Bg'],1.3*(mfitPars['Bg']+mfitPars['Peak'])])
                    if gfitArray != None:
                        pylab.plot(gfitArray.max(1), 'g',label='gauss fit FWHM='+str(gfitPars['FWHMarcsec']))
                        pylab.axis([0,stampPars[1]-stampPars[0],gfitPars['Bg'],1.3*(gfitPars['Bg']+gfitPars['Peak'])])
                    pylab.legend()
                    
                    

                if keystroke == 'e':
                    print 'PYEXAM - Plotting contours'
                    pylab.clf()
                    pylab.contour(stampArray)
                                
            elif keystroke == 'q':
                done = 'yes'
                if verbose: print "Done looping through stars"

        if clip:
            if gAllstars:
                if verbose: print "# PYEXAM - performing Gaussian clipping"

                gAllstars, gEllMean, gEllSigma, gFWHMMean, gFWHMSigma = \
                   iqUtil.sigmaClip(gAllstars, outFile, sigma, verbose, niters, garbageStat=garbageStat)
            
                if verbose:
                    print "PYEXAM - Mean Gaussian Ellipticity:", gEllMean
                    print "PYEXAM - Mean Gaussian FWHM (""):", gFWHMMean
                    
            if mAllstars:
                if verbose: print "# PYEXAM - performing Moffat clipping"

                mAllstars, mEllMean, mEllSigma, mFWHMMean, mFWHMSigma = \
                   iqUtil.sigmaClip(mAllstars, outFile, sigma, verbose, niters, garbageStat=garbageStat)

                if verbose:
                    print "PYEXAM - Mean Moffat Ellipticity:", mEllMean
                    print "PYEXAM - Mean Moffat FWHM (""):", mFWHMMean
         

        for star in gAllstars:
            iqUtil.writePars(star, outFile, 'gaussian')

        for star in mAllstars:
            iqUtil.writePars(star, outFile, 'moffat')

 #---------------------------------------------------------------------------

def pyiq (filename, scidata, function, outFile, pixelscale, \
              frame=1, pverbose=True, pymark=True, clip=True, \
              sigma=2.3, niters=4, display=True, imageSigma='default', \
              boxSize=9., residuals=False, saturation=65000,\
              debug=False, xmin=None, xmax=None, ymin=None, ymax=None, garbageStat=False, qa=False):
    '''finds stars and fits gaussian/moffat/both to them

    @param filename: input filename with ".fits[sci,i]" attached
    @type filename: string

    @param scidata: science data array, equal to pyfits.getdata(filename)
    @type scidata: numpy array

    @param function: currently supported are "gauss", "moffat", or "both"
    @type function: string

    @param outFile: opened file object
    @type outFile: string        
    
    @param pixelscale: instrument pixelscale in arcsec/pix
    @type pixelscale: float or int

    @param frame: frame number for marking, only used for multiple science extensions
    @type frame: int

    @param pverbose: warnings will be printed
    @type pverbose: Boolean

    @param pymark: mark displayed images with daofind centers and fwhm circles for fits
    @type pymark: Boolean

    @param clip: Sigma clip outliers in measured FWHM and ellipticities
    @type clip: Boolean

    @param sigma: Sigma threshold for clipping outlier FWHM and ellipticities
    @type clip: float or int

    @param niters: number of times sigma clipping is iterated through parameter list
    @type niters: int

    @param display: Display images?
    @type display: Boolean

    @param imageSigma: sigma of background level, if "default" a star mask is made
                       and std is measured
                       directly from masked image
    @type imageSigma: string

    @param boxSize: aperture size around object x,y to make postage stamp for the
                       fitting function
    @type boxSize: float

    @param residuals:this is problem
    @type residuals: Boolean

    @param saturation: saturation values for daofind inputs
    @type saturation: float or int
    
    @param qa: A flag to use a grid of sub-windows for detecting the sources in 
               the image frames, rather than the entire frame all at once.
    @type qa: Boolean
        '''
    
    if display!=True: 
        pymark=False

    # instantiate the logger object and put into the global variable 
    global log
    if log==None: 
        log = gemLog.getGeminiLog()
    if debug:   
        # update logger object's file handler level to debug
        log.changeLevels(logLevel=log.loglevel(), debug=True)

    # Initialize arrays
    gAllstars = []
    mAllstars = []
    xyArray = []
    if residuals:
        gsubtracted = scidata[:]
        msubtracted = scidata[:]
 
    # Find Objects
    xyArray = iqUtil.pyDaoFind(filename, scidata, pixelscale, 
                               #frame=frame, debug=debug, pymark=pymark, 
                               #display=display, imageSigma=imageSigma, 
                               frame=frame, debug=debug, pymark=False, 
                               display=False, imageSigma=imageSigma, 
                               saturation=saturation, qa=qa)
    
    # Get box size to put around each object
    apertureSize = fit.getApSize(pixelscale, boxSize)
    
    # Remove objects that are too close to each other
    xyArray = iqUtil.removeNeighbors(xyArray, npixapart=apertureSize,
                                   crowded=True, debug=debug)
   
    # Remove objects that are too close to the edge
    xyArray = iqUtil.edgeCheck(scidata,xyArray,npixaway=apertureSize-1,
                             debug=debug, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

    # Label printed value columns
    if debug:
        log.debug('%10s%12s%12s%12s%12s%12s%12s'%('# function', 'Coox',
                'CooY', 'FWHMpix', 'FWHMarcsec', 'Ellip', 'PAdeg'))

     # Label output file columns 
    if frame==1:
        outstring = '%10s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%5s\n'\
        %('# function', 'Cx', 'Cy', 'Bg', 'Peak', 'Wx', 'Wy', 'CooX', 'CooY', \
        'Ellipticity', 'FWHM_pix', 'FWHM_arcsec', 'Theta', 'PA_deg', 'FWHMx', \
         'FWHMy', 'Beta', 'CCD')
        outFile.write(outstring)
        
    # pymark is super slow, mark only the first 10 objects found
    count = 0

    # Loop through objects and fit functions 
    for [xCoo,yCoo] in xyArray:
        log.debug('Starting to loop through the objects found in the images'+
                  ' data frame '+str(frame))
        
        positionCoords= [xCoo,yCoo]
        stampArray, stampPars = fit.makeStamp(scidata, positionCoords, 
                                          apertureSize, outFile, debug=debug)
        # Start time for model fitting
        st = time.time()
        # Call fitprofile to perform model fitting to objects
        gfitPars, gReturnModel, mfitPars, mReturnModel = \
            fit.fitprofile(stampArray, function, positionCoords, 
            outFile, pixelscale, stampPars, debug=debug, frame=frame)
        # End time for model fitting
        et = time.time()
        # Logging the time elapsed during model fitting
        # Below log message was made a debug level as it is repeated for 
        # every object, thus dozens of calls per image, convert back to 
        # fullinfo level in future if deamed more appropriate.
        log.debug('Time to perform gaussian and/or moffat fits: '+str(et-st))

        ######## Gaussian Fits #########
        if gfitPars != None:
            # Below log message was made a debug level as it is repeated for 
            # every object, thus dozens of calls per image, convert back to 
            # status level in future if deamed more appropriate.
            log.debug('Gaussian fits to the stars is being performed')
            
            clipobsPars = ('gauss', gfitPars['CooX'], gfitPars['CooY'],
                           gfitPars['FWHMpix'], gfitPars['FWHMarcsec'],
                           gfitPars['Ellip'], gfitPars['PAdeg'])
            iqUtil.printPars(clipobsPars, debug)
            gAllstars.append(gfitPars)
        
            if pymark:
                if count <= 10:
                    try:
                        iqUtil.iqmark (frame, gfitPars, color=208)
                    except:
                        log.warning('Could not mark stars; is ds9 running?')
                        pymark = False
                else:
                    log.warning('>10 objects found; only first 10 are marked')
                    pymark = False

            if residuals:
                if gAllstars[0] == gfitPars:
                    # Replacing contents of gsubtracted with the gaussian 
                    # subtracted data frame
                    gsubtracted = iqUtil.makeResiduals(scidata, gfitPars, 
                                                       gReturnModel, stampPars)
                else:
                    #$$ Since gsubtracted is set to scidata when initialized
                    #$$ and there are no changes to either gsubtracted or scidata 
                    #$$ above this, this if/else block kinda seems useless??
                    gsubtracted = iqUtil.makeResiduals(gsubtracted, gfitPars, 
                                                       gReturnModel, stampPars)


        #########  Moffat Fits ######### 
        if mfitPars != None:
            # Below log message was made a debug level as it is repeated for 
            # every object, thus dozens of calls per image, convert back to 
            # status level in future if deamed more appropriate.
            log.debug('Moffat fits to the stars is being performed')
            
            clipobsPars = ('moffat', mfitPars['CooX'], mfitPars['CooY'],
                               mfitPars['FWHMpix'], mfitPars['FWHMarcsec'],
                               mfitPars['Ellip'], mfitPars['PAdeg'])
            iqUtil.printPars(clipobsPars, debug)
            mAllstars.append(mfitPars)
                        
            if pymark:
                if count <= 10:
                    try:
                        iqUtil.iqmark (frame, mfitPars, color=206)
                        count+=1
                    except:
                        log.warning('Could not mark stars; is ds9 running?')
                        pymark = False
                else:
                    log.warning('>10 objects found; only first 10 are marked')
                    pymark = False

            if residuals:
                if mAllstars[0] == mfitPars:
                    # Subtract first object
                    msubtracted = iqUtil.makeResiduals(scidata, mfitPars, \
                                               mReturnModel, stampPars)
                else:
                    # Residual subtract already subtracted image
                    msubtracted = iqUtil.makeResiduals(msubtracted, mfitPars, \
                                                mReturnModel, stampPars)
            

        if debug: 
            log.status('Done looping through objects')

    # Write the final residual subtracted image to a fits file
    if residuals:
        resName = gemutil.removeExtension(outFile.name)
        
        if gAllstars != []:
            # Removing previous copy of the file from disk, then 
            # writing a fits file with the gaussian subtracted data frame
            filesystem.deleteFile(resName+'gaussResidual'+str(frame)+'.fits')
            pyfits.writeto(resName+'gaussResidual'+str(frame)+'.fits',\
                gsubtracted)
        if mAllstars != []:
            # Removing previous copy of the file from disk, then 
            # writing a fits file with the moffat subtracted data frame
            filesystem.deleteFile(resName+'moffatResidual'+str(frame)+'.fits')
            pyfits.writeto(resName+'moffatResidual'+str(frame)+'.fits',\
                msubtracted)       

    if clip: # This doesnt get called for GMOS as gemiq does the clipping
             # should be a separate routine in iqUtil as it is used in gemiq
        if gAllstars:
            if debug: 
                log.debug("Performing Gaussian clipping on frame "+str(frame)+
                          ' of file '+filename)
            gAllstars, gEllMean, gEllSigma, gFWHMMean, gFWHMSigma = \
               iqUtil.sigmaClip(gAllstars, outFile, sigma, pverbose, niters, garbageStat)
            for star in gAllstars:
                writePars(star, outFile, 'gaussian')
                clipobsPars = ('gauss', star['CooX'], star['CooY'],
                               star['FWHMpix'], star['FWHMarcsec'],
                               star['Ellip'], star['PAdeg'])
                iqUtil.printPars(clipobsPars, pverbose)

            g1 = "Gauss Ellipticity: "+str(gEllMean)+" +/- "+str(gEllSigma)
            g2 = "Gauss FWHM (arcsec): "+str(gFWHMMean)+" +/- "+str(gFWHMSigma)+'\n'
            if pverbose:
                log.stdinfo(g1+'\n'+g2)
            outFile.write(g1+'\n'+g2)
      
        if mAllstars:
            if debug: 
                log.debug("Performing Moffat clipping on frame "+str(frame)+
                          ' of file '+filename)
            mAllstars, mEllMean, mEllSigma, mFWHMMean, mFWHMSigma = \
                iqUtil.sigmaClip(mAllstars, outFile, sigma, pverbose, niters,
                                  garbageStat=garbageStat)
            for star in mAllstars:
                iqUtil.writePars(star, outFile, 'moffat')
                clipobsPars = ('moffat', star['CooX'], star['CooY'],
                               star['FWHMpix'], star['FWHMarcsec'],
                               star['Ellip'], star['PAdeg'])
                iqUtil.printPars(clipobsPars, pverbose)

            m1 = "Moffat Ellipticity: "+str(mEllMean)+" +/- "+str(mEllSigma)
            m2 = "Moffat FWHM (arcsec): "+str(mFWHMMean)+" +/- "+str(mFWHMSigma)+'\n'
            if pverbose:
                log.stdinfo(m1+'\n'+m2)
            outFile.write(m1+'\n'+m2)
    if debug:
        log.debug("getiq: #gaussian stars/#moffat stars = "+
                  str(len(gAllstars))+ "/"+ str(len(mAllstars)))
    return gAllstars, mAllstars

 #---------------------------------------------------------------------------






