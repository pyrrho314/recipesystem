ncprepare. Find masks center
============================

**ncprepare(inputs oprefix='n' idir='' odir='' fdir='' fsuffix='' dobadpix=True clobber=False logfile='' verbose=False)**

 Ncprepare is a Python script that takes raw NICI data with 2 FITS extensions 
 and calculates the center of each mask -interactively if necessary, adding 
 these values to the header. It will do this after the frames are flat fielded 
 and the blue frame is registered to the red frame coordinate system. The frames 
 are also shifted so that the mask centers are at (512,512). This 
 is a require step before running **ncscience**.  

**Parameters**

* *inputs*
   A input list of FITS files to process. This list can a Unix wildcard pathname, e.g. * .fits, root23[2-9].fits, root??.fits or a @ list, e.g. @file.lis, where 'file.lis' is a text file with a list of FITS files, one file per line or a plain list of FITS filenames separated by commas. 

   **NOTE** If you ran the task **ncmark** with these input files then you should input the output files as input in this script. The script will read the **x and y centers** from the headers.

* *oprefix*
    Default value is ' n'. Is the prefix used for the output filenames. 

* *idir*
    Default is current directory. Directory pathname where the input files reside. 

* *odir*
    Default is current directory. Directory pathname to put the output FITS files. 

* *fdir*
   Directory name where the flats are. The files are: flats_red_<fsuffix>.fits, flats_blue_<fsuffix>.fits, dark_red_<fsuffix>.fits and dark_blue_<fsuffix>.fits. 

* *fsuffix*
   Suffix used by the Calibration files (ncmkflats). If default it will used the
   *suffix* value.

* *dobadpix*
    Default value is True. Correct badpixels the best we can.

* *clobber*
    Default value is False. Set to True to overwrite.

* *logfile* 
     Log filename to hold the script messages. The default name is *gemini.log*

* *verbose* 
    Default value is False. If True the information goes to the terminal as well.


**Mask Centroid notes**

    Mask centroid is done automatically and the 2 FITS header of the output FITS file will have XCEN and YCEN keyword with its coordinates. If the finding algorithm fails then ncprepare will go into "interactive" mode using DS9 to display the frame.

1. Mark the center with left button, then hit 'q' to continue or 's' to skip this frame.
2. The frame is displayed again but at higher resolution. Mark again and press 'q' to continue.

**Examples**

1. ncprepare '*.fits' --odir='/data' --fdir=/data/flats --fsuffix=S20100111

   Prepare all the FITS files in the current directory, find the mask center 
   and update the headers. Write the output files in '/data'. The 'Flats' files 
   are in '/data/flats' and their suffix is 'S20100111'.

2. ncprepare @niciFiles.lis idir='/data/20090312/' odir='/data/reduced' fdir=/data/flats fsuffix=S20100111 clobber=yes (Pyraf mode) 

3. ncprepare @niciFiles.lis --idir='/data/20090312/' --odir='/data/reduced' --fdir=/data/flats --fsuffix=S20100111  --clobber (Unix mode) 

   The input FITS files are in the list file 'niciFiles.lis' as one 
   filename per line. You can put the full pathname of each file in 
   which case do not specified 'idir'. If only filenames are given, 
   then the script will open the FITS files in 'idir'. The *flats* calibration
   directory is in '/data/flats' and the suffix that thos flats have is 
   'S20100111'. The output 
   files are written to 'odir' pathname. Remember that in Unix mode 
   you can get the list of this script by typing 'ncprepare -h'. 

