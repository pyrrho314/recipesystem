ncmark. Manually marks mask centers
====================================

**ncmark(inputs idir='' odir='' oprefix='m' port=5137 logfile='' clobber=False verbose=False)**

    The ncmark script display each frame for the user to manually click on the
    center of the mask. These coodinates are input to a centroid routine for
    finer centering. Uses **ds9**

**Description**

    Sometimes the script **ncprepare** will not correctly compute the
    center of the mask and the manual step of the scripit will also fail. **ncmark**
    is simpler and faster for this purpose. Any small signal in the middle of the
    mask is sufficient to mark the center. The *xcen and ycen* values are written
    to the extension header os the frames. The user marks the center of the mask
    with the left mouse button and hits *q* to quit or *s* to skip the frame.

    **NOTE**
        Make sure you use these output files as input to the *ncprepare* script.

**Parameters**

* **inputs**
    Input filenames. It can be a comma separated list of FITS files, a file
    list containing a list of files (one per line) or a unic wildcard (the string 
    needs to be enclosed in quotes).

* **idir**
    The input directory where the input fits files are located. If left blank, inputs should included directory pathname. 

* **odir**
    The output directory where all the listing and fits files cube will be written. If left blank, they will written in the working directory. 

* **oprefix** 
    Default value is 'm'. Prefix string for the output filenames.

* **logfile**
          A filename that will contain log information from the ncmark.
* **port**
          Default value is 5137. Choose another value if you already have a ds9 application running. Make sure you start ds9 port the '-port #' parameters.

* **clobber**
          Default value is False. If the output file already exists the program will stop. If the value is True, it will delete the file first.

* **verbose**
          Default value is False. The information messages from the program are into the logfile only. If the value is True they are also displayed on the terminal.

**Examples** 
 
1. ncmark '/home/data/S20100109S002*.fits' --odir='/tmp'  -v --logfile='/home/ncmark.log'

    Shows each of the red-blue frames on a ds9 display; the users clicks on the center of the mask with the left mouse button and then hits *q* to accepts the mouse inputs or *s* to skip this frame.

