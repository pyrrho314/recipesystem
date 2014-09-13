

Installation
------------

- Hardware requirements: Intel Linux and Mac. NO MS windows OS

- Fluxcal is written in Python and requires the Gemini/Astrodata package

  - Gemini/Astrodata
  - Stsci_python package with Pyraf and pytools

- You can get the module from the SVN repository: Fluxcal_

.. _fluxcal: http://chara.hi.gemini.edu/svn/DRSoftware/gemini_python/trunk/devel/fluxcal/
  - Make sure you have a pathname in your PYTHONPATH variable that points 
      to the directory where you installed FLUXCAL

**BAD PIXEL MASKS**
-------------------

- At this time masking an image using the default BPM files is done
  either as a primitive with *addBPM* or using the fluxcal.py function *addbpm*

   - addBPM
      This primitive uses: Gemini/GMOS/GMOS_BPM_11.fits or Gemini/GMOS/GMOS_BPM_22.fits
      according to the value of CCDSUM. WARNING: These files add artificial objects

   - fluxcal/addbpm
      With the method *addbpm* from fluxcal you can use *better* masks even 
      though some artificial objects could appear at the masked edges. The local
      BPM files are:

      - bpmgmosn11_mosaic.fits.gz   # Mosaic 6218x4608   gmos-n   CDSUM 11
      - bpmgmosn22_123.fits.gz      # MEF 1056x2304      gmos-n   CDSUM 22
      - bpmgmosn22_mosaic.fits.gz   # Mosaic 3108x2304   gmos-n   CDSUM 22
      - bpmgmoss11_mosaic.fits.gz   # Mosaic 6218x4608   gmos-s   CDSUM 11
      - bpmgmoss22_123.fits.gz      # MEF 1056x2304      gmos-s   CDSUM 22
      - bpmgmoss22_mosaic.fits.gz   # Mosaic 3108x2304   gmos-s   CDSUM 22   
      - bpmgmosn11_123.fits.gz   TBD
      - bpmgmoss11_123.fits.gz   TBD
      
      ::

       from astrodata import AstroData
       import fluxcal as fc

       ad=AstroData('mgS20101214S0040.fits')
     
       # Create Fluxcal object. Use BPM in detectsources()
       ff=fc.Fluxcal(ad, addBPM=True,sigma=0.0, threshold=2.5, fwhm=5.5))

       # This step is optional.
       ff.addbpm()            # Mask the images with the appropiate 
                              # BPM according to CCDSUM
                        
       adout=ff.runFC()       # Run fluxcal on ad

       
Running fluxcal
---------------

- FITS file requirements:

 - Image field must be in an IMAGE extension with EXTNAME value 'SCI' and EXTVER of 1.
 - For the Zero point calculation the image header must have
   EXPTIME, AIRMASS keywords. If the input file is not a Gemini file
   then the extinction value is required as one of the input parameters.
 - WCS information to convert x,y to ra,dec.

- Fluxcal will run on all the IMAGE extensions in a FITS file.

  **RESULTS**

  - ZPcorr.log: The output text file with ZP correction values plus others.

  - Input FITS file. The input file will have 2 or more BINTABLE extensions
    with parameters like (x,y) (ra,dec) and magnitudes for each source
    found by the source detection algorithm. The other table extension 
    also contains parameters for each reference start found in the field.
    The EXTNAME for theses are OBJCAT and REFCAT.

::
 
  import fluxcal as fc
  from astrodata import AstroData

  ad = AstroData('rgN20100920S0661.fits')
  # Run the individual scripts
  # Pass threshold=3.5 to detectsources()
  ff = fc.Fluxcal(ad,threshold=3.5)    # Run the individual scripts

  # OR if your want to run detectsources only then

  adout = ff.ds()    # 'ds' is an alias for detectsources.

  # If reference objects were found in the field then the ZP correction
  # was calculated. The results are in the local file "fluxcal.log"
  # OUTPUT FILE IS: 
            zp_rgN20100920S0661.fits

Running individual scripts
--------------------------

You can also run the individuals scripts in another way to allow
you to change some of the available parameters.


::

  from astrodata import AstroData
  import detectSources as ds

  ad = AstroData('gS20101105S0128.fits')
  dd = ds.DetectSources(ad,addBPM=True, sigma=0.0, threshold=2.5, fwhm=5.5)
  adout = runDS() 

  # If you want to save the results please do.
  adout.write('ds_gS20101105S0128.fits')
 
  # this file will contain the OBJCAT table with the source detected
  # parameters.

  # CONTINUE

  import addReferenceCatalogs as ref

  rr = ref.AddReferenceCatalogs(adout)
  adout = rr.getRefs()

  adout.info()     # Now the object contains OBJCAT and REFCAT
  


