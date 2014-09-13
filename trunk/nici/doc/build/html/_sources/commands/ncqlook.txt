ncqlook. Quick look and data quality assesment. 
===============================================

**ncqlook(inputs idir='' odir='' log=True lists=True saturate=5000 nodisplay=False full=False port=5137)**

The ncqlook script produces a quick look analysis of the nici raw files specified in the 'inputs' parameter. It will produce as output a cube of FITS files (512x256xnumber_of_files) in the working directory plus several output files if the 'lists' parameter is kept True. See Parameters for more detail. While is running, each pair of frames are displayed on a ds9 panel.

**Parameters**

* **inputs**
          If left blank then last night NICI raw files resident in the Gemini South repository /net/petrohue/dataflow will be processed. If you want to display data from the repository from a different date, then the format is of the form YYYYMMDD. You can also give a list of files or a unix wild card. See examples. 

* **idir**
          The input directory where the input fits files are located. If left blank, inputs should included directory pathname. 

* **odir**
          The output directory where all the listing and fits files cube will be written. If left blank, they will written in the working directory. 

* **log**
          The default value is True. Will create a log with filename, min-max rms, median value for extension 1 and 2 and values for keywords OBJECT, OBSCLASS, OBSTYPE, MODE, ITIME, NCOADD and optionally 2 numbers representing the core2halo ratio; if one or both of these numbers are missing it means that the algorithm failed to get a meaningful value

* **lists**
          The default value is True. Will create several output files useful 

* **saturate**
          Saturation limit. Default value is 5000. To change give a value.

* **nodisplay**
          Default value True for not displaying current frames.

* **full**
          Default value False. If True it will rebin the frame to 256x256.

**Output files**

* *root_cube.fits*
     FITS file cube.
* *root.log*
    For each frame it contains min-max, and median listing. The values ADI,SDI and 
    ASDI are computed from keywords CRMODE and DICHROIC. The last 4 fields in the 
    log are Exposure time, Ncoads, Core2Halo ratio for red and blue frames.
  * **ADI**
    The frame has this mode if CRMODE is FIXED and DICHROIC has 'Mirror' in the value
    field.
  * **SDI**
    The frame has this mode if CRMODE is FOLLOW and DICHROIC has '50/50' in the value
    field.
  * **ASDI**
    The frame has this mode if CRMODE is FIXED and DICHROIC has '50/50' in the value
    field.
* *root.1_flats*
    Contains calibration files for the ADI mode 
* *root.2_flats*
   Contains calibration files for the ASDI and SDI mode. 
* *root.(adi,sdi,asdi)*
    Contains science object listings. NOTE that these files can have listings of more than one object. You will need to edit these files and create one list per object if you want to use them in ncprepare and ncscience the log file has the necessary information for this. 

**Examples** 
 

1. ncqlook 

    Will do a quick analysis of all the NICI FITS files residing in /net/petrohue/dataflow for the date of last night, displaying each pair of frames on a ds9 frame while a listing of the log file runs on your screen. 

2. ncqlook 20090313 --odir='/tmp' --saturate=3500 --nodisplay

    (Unix command mode)

    List all the NICI fits files from /net/petrohue/dataflow/S20090313S*.fits The output listing will be written in the '/tmp' directory. No display is produced, so ds9 need not be running.

    The output files are:

        * 200903013_cube.fits
        * 200903013.log
        * 200903013.1_flats
        * 200903013.2_flats
        * 200903013.adi
        * 200903013.sdi
        * 200903013.asdi 

3. ncqlook(20090313,odir='/tmp',nodisplay=True)

    This is the syntax for the command in the PYTHON shell. 

4. ncqlook "/data/nici/200903/S2009*.fits" --odir='/tmp' full=True 

    Check all the fits files in the given directory writing the listing
    and cube in the '/tmp' directory. '--full' is the flag to tell
    ncqlook to rebin the frames to 256x256.

