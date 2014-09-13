#! /usr/bin/env python

from astrodata.AstroData import AstroData, prepOutput
import os
from copy import deepcopy

def ADUToElectron(filelist, odir, oprefix):
    """
     This is a function to convert the ADU counts to electrons
     by multiply the pixel values by the gain.
     Arguments:
       filelist: A python list of FITS filenames
       odir:     Directory pathname for output FITS files
       oprefix:  Prefix for output filenames. Example: If input filename
                 is 'S20100323S0012.fits' and 'oprefix' is 'n', the output 
                 name will be 'nS20100323S0012.fits'
    """

    # Loop through the files in filelist
    for filename in filelist:
        # Open the file as an AstroData object
        adinput = AstroData(filename, mode='readonly')

        # Verify whether the data has already been converted to electrons
        if adinput.phuValue('ELECTRON') != None:
            print "WARNING: File %s has already been converted to electrons"\
                   % filename
            return

        outputname = oprefix + os.path.basename(filename)
        ofile = os.path.join(odir,outputname)
        if os.access(ofile, os.F_OK): os.remove(ofile)

        # Prepare a new output
        #    Propagate PHU and MDF (if applicable) to output.
        #    No pixel extensions yet.
        #    Set output file name.
        #    No overwrite allowed. (default mode for prepOutput)
        #
        # prepOutput copies the adinput PHU and set the name of the new
        # file represented by adout to ofile.

        adout = prepOutput(adinput, ofile)

        # Get the gain values to apply
        # adinput.gain() returns a list, one value for each science extension.
        gain = adinput.gain()

        # Use the deepcopy function to create a true copy and ensure that
        # the original is not modified.s

        adc = deepcopy(adinput) 

        # Multiply each science extension by the gain.
        # Append new extension to already prepared output.
        for extension,g,xn in zip(adc, gain, range(len(gain))):
            extension.data = extension.data * g

            adout.append(data=extension.data, header=extension.header)

        # Update PHU with timestamps
        adout.phuSetKeyValue('ELECTRON', fits_utc(), 
            comment='UT Modified with convertToElectrons')
        adout.phuSetKeyValue('GEM-TLM', fits_utc(), 
            comment='UT Last modification with GEMINI')

        # Write to disk.  The filename was specified when 
        # prepOutput was called.
        adout.write()

        # Close files
        adout.close()
        adc.close()
        adinput.close()

import time
def fits_utc():
   """
     Return a UTC string in FITS format:
     YYYY-MM-DDThh:mm:ss
   """

   gmt = time.gmtime()
   time.asctime(gmt)
   fitsT = '%d-%02d-%02dT%02d:%02d:%02d' % gmt[:6]

   return fitsT

if __name__ == "__main__":

    import optparse

    VERSION = '1.0'

    # Parse input arguments
    usage = 'usage: %prog [options] file1 .. fileN'
    p = optparse.OptionParser(usage=usage, version='v'+VERSION)
    p.add_option('--oprefix', '-p', action='store', dest='prefix', default='elec_',
        help='Prefix for the output files')
    p.add_option('--odir', action='store', default='', help='Output directory pathname')

    (options, args) = p.parse_args()

    ADUToElectron(args, options.odir, options.prefix)

