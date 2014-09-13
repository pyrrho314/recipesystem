#!/usr/bin/env python
#
import optparse
import niciTools as nt
from niciTools import check_dir,robust_sigma,getFileList,nici_noise
from scipy.signal import medfilt2d
import logging as log

import os
from os.path import join
import numpy as np
import pyfits as pf
import time
import sys

command = ''   # Global variable

def ncmkflats(inputs, idir='',odir='',sigma=6, clobber=False,\
                suffix='default', logfile='', verbose=False):
    """
        inputs = ""         NICI flats and darks file list
        (idir = "")         Path for input directory
        (odir = "")         Path for output directory
        (sigma=6)           Set to Nan all pixel above this value from the median
        (clobber = False)   Clobber output file?
        (suffix='default')  If 'default' it will take the rootname of 
                            the first element in the input list.
        (logfile = "")      Logfile
        (verbose = False)   Display progress on monitor


     PARAMETERS:
        Sigma=6  -- All pixels in the flat off (dark)
                 that are way above the median value
                 are set to NaN *** The selection is
                 done in the following way :
                 subtract a low-pass to the flat_off
                 to keep only single 'hot' pixels,
                 and select all pixels more than 6
                 sigma (estimated through
                 robust_sigma) of the median-filtered image.


     DESCRIPTION:
        Produces calibration  given a list of fits files with flats,
        darks and sky exposures

     OUTPUT FILES:
        flat_red_<suffix>.fits
        flat_blue_<suffix>.fits
        sky_red_<suffix>.fits
        sky_blue_<suffix>.fits

        Note: This code has skycube<color,suffix>.fits writeto commented.

     EXAMPLES: 
        See examples in browser:
        file:///<nici_dir>/doc/build/html/ReductionExamples.html
    """

    global command
    log.basicConfig(filename='ncmkflats.log',level=log.DEBUG)

    if odir=='': odir='./'
    if idir=='': idir='./'
    log.info('\n\n START NICI makeFlats +++++++++++++++++++++++++')
    log.info(command)
    log.info( 'inputs_list_len: '+str(len(inputs)))
    log.info('idir: '+idir)
    log.info('odir: '+odir+'\n')

    try:
        odir = check_dir(odir,os.W_OK)
        idir = check_dir(idir,os.R_OK)
    except IOError, strerror:
        log.warning('(%s)' % strerror)
        return

    # We assumed the input directory is specified or is explicit
    # in the input file list.

    # If odir is empty check if we can write in the CWD.
    if odir =='':
        if (not os.access('./', os.W_OK)):
           log.warning("*** Warning: Cannot write output file in the CWD: './'")
        return

    file_ls = getFileList(inputs)

    if len(file_ls) == 0:
        log.warning("\n WARNING: Input list is empty.\n")
        return

    # Find rootname from the first element of 'inputs' if
    # 'suffix' is not 'default'

    if suffix == 'default':
       root = os.path.basename(file_ls[0]).split('.')[0]
       # If root is of the form 20yymmddSnnnn ==> get 20yymmdd
       if len(root) == 14 and root[:3] == 'S20': root = root[:9]
       suffix = root

    # Creare dark and flat lists (same as flats_off or darks)
    dark_ls=[]
    flat_ls=[]
    ti=time.time()
    for file in file_ls:
        hdulist = pf.open(idir+file)
        ah0 = hdulist[0].header

        if ah0['OBSTYPE'] == 'DARK':
            dark_ls.append(file)
        elif ah0['OBSTYPE'] == 'FLAT':
            flat_ls.append(file)
        else:
            continue
	hdulist.close()

    ndarks = len(dark_ls)

    if ndarks == 0:
       log.error("NO FLATS-CLOSED (darks) were encountered in the input list.")
       return
    if len(flat_ls) == 0:
       log.error("NO FLATS were encountered in the input list.")
       return
    cube_dark_r = np.empty([ndarks,1024,1024],dtype=np.float64)
    cube_dark_b = np.empty([ndarks,1024,1024],dtype=np.float64)
    k = 0
    for file in dark_ls:

        hdulist = pf.open(idir+file)
        ah0 = hdulist[0].header
        ah1 = hdulist[1].header
        ah2 = hdulist[2].header

        line = '%s %s %s' % (file,ah0['obstype'],ah0['obsclass'])
        log.info(line)

        # Shutter is closed.
        linearzR = ah1['ITIME_R']*ah1['NCOADD_R']
        linearzB = ah2['ITIME_B']*ah2['NCOADD_B']

        cube_dark_r[k,:,:] = hdulist[1].data/linearzR
        cube_dark_b[k,:,:] = hdulist[2].data/linearzB
        k = k + 1

	hdulist.close()

    t1=time.time()

    dark_r_median = np.median(cube_dark_r,axis=0)
    dark_b_median = np.median(cube_dark_b,axis=0)


    fh = pf.getheader(idir+dark_ls[0])

    # Let's use the cube_dark to calculate cube_sky. The difference is the
    # 'nici_noise' application to each slice.

    cube_sky_r = cube_dark_r
    for i in range(len(dark_ls)):
        cube_sky_r[i,:,:] = nici_noise(cube_dark_r[i])
        
    ofile = join(odir+'sky_red_'+suffix+'.fits')
    log.info("Writing "+ ofile)
    pf.writeto(ofile, np.median(cube_sky_r.astype(np.float32),axis=0),
               header=fh,clobber=clobber)
    #ofile = join(odir,'skycube_red_'+suffix+'.fits')
    #log.info("Writing "+ ofile) 
    #pf.writeto(ofile, cube_sky_r.astype(np.float32),header=fh,clobber=clobber)
    del cube_dark_r
    del cube_sky_r

    cube_sky_b = cube_dark_b
    for i in range(len(dark_ls)):
        cube_sky_b[i,:,:] = nici_noise(cube_dark_b[i])

    ofile = join(odir,'sky_blue_'+suffix+'.fits')
    log.info("Writing "+ ofile)
    pf.writeto(ofile, np.median(cube_sky_b.astype(np.float32),axis=0),
              header=fh,clobber=clobber)

    #ofile = join(odir,'skycube_blue_'+suffix+'.fits')
    #log.info("Writing "+ ofile)
    #pf.writeto(ofile, cube_sky_b.astype(np.float32),header=fh,clobber=clobber)
    del cube_dark_b
    del cube_sky_b

                              # *** here, all pixels in the flat off (dark)
                              # that are way above the median value
                              # are set to NaN *** The selection is
                              # done in the following way :
                              # subtract a low-pass to the flat_off
                              # to keep only single 'hot' pixels,
                              # and select all pixels more than 6
                              # sigma (estimated through
                              # robust_sigma) of the median-filtered image.

    dev_red = dark_r_median - medfilt2d(dark_r_median,7)
    dev_red /= robust_sigma(dev_red)

    dev_blue = dark_b_median - medfilt2d(dark_b_median,7)
    dev_blue /= robust_sigma(dev_blue)

    dark_r_median[np.where(np.abs(dev_red) > sigma)] = np.NaN
    dark_b_median[np.where(np.abs(dev_blue) > sigma)] = np.NaN

    # Now get cube_flat
    nflats = len(flat_ls)
    cube_flat_r = np.empty((nflats,1024,1024),dtype=np.float64)
    cube_flat_b = np.empty((nflats,1024,1024),dtype=np.float64)
    k = 0
    for file in flat_ls:

        hdulist = pf.open(idir+file)
        ah0 = hdulist[0].header
        ah1 = hdulist[1].header
        ah2 = hdulist[2].header

        line = '%s %s %s %s' % (file,ah0['obstype'],ah0['GCALSHUT'],ah0['obsclass'])
        log.info(line)

        linearzR = ah1['ITIME_R']*ah1['NCOADD_R']
        linearzB = ah2['ITIME_B']*ah2['NCOADD_B']
        cube_flat_r[k,:,:] = hdulist[1].data/linearzR
        cube_flat_b[k,:,:] = hdulist[2].data/linearzB
        k = k + 1

	hdulist.close()

    log.info('Number of flats_red:  '+str(len(flat_ls)))
    log.info('Number of flats_blue: '+str(len(flat_ls)))
    log.info('Number of darks_red: '+str(ndarks))
    log.info('Number of darks_blue:'+str(ndarks))
  
    for i in range(nflats):
        cube_flat_r[i,:,:] -= dark_r_median
        cube_flat_r[i,:,:] /= np.median(np.nan_to_num(cube_flat_r[i,:,:]))

    if (nflats > 1):
       flat_r_median = np.median(cube_flat_r,axis=0)
    else:
       flat_r_median = cube_flat_r
    del cube_flat_r

    # set very low or very high values of the flat to 
    # NaN (<5% or >150% of the median value)

    flat_r_median[np.where((flat_r_median > 1.5)|(flat_r_median < 0.05))] = np.nan

    fhon = pf.getheader(idir+flat_ls[0])
    file_red = join(odir,'flat_red_'+suffix+'.fits')
    log.info("Writing "+ file_red)
    pf.writeto(file_red,flat_r_median.astype(np.float32),
         header=fhon,clobber=clobber)

    #FLAT BLUE

    for i in range(nflats):
        cube_flat_b[i,:,:] -= dark_b_median
        cube_flat_b[i,:,:] /= np.median(np.nan_to_num(cube_flat_b[i,:,:]))

    if (nflats > 1):
       flat_b_median = np.median(cube_flat_b,axis=0)
    else:
       flat_b_median = cube_flat_b
    del cube_flat_b
    
    flat_b_median[np.where((flat_b_median > 1.5)|(flat_b_median < 0.05))] = np.nan
    file_blue = join(odir,'flat_blue_'+suffix+'.fits')
    log.info("Writing"+ file_blue)
    pf.writeto(file_blue,flat_b_median.astype(np.float32),
         header=fhon,clobber=clobber)
    
    tt = time.time()-ti
    print 'TIME:','%dm:%.2ds' % (tt/60,tt%60)
    return 'Make flats done'

if __name__ == '__main__':

    # Parse input arguments
    usage = 'usage: %prog file_list [options]'

    p = optparse.OptionParser(usage=usage, version='ncmkflats Sep 2009')
    p.add_option('--idir', action='store', type='string',
                 default='', help='Input directory pathname')
    p.add_option('--odir', action='store', type='string',
                 default='', help='Output directory pathname')
    p.add_option('--clobber', default=False, action='store_true',
                 help='Clobber output files?')
    p.add_option('--sigma', default=6, type='float', action='store',
                 help='Sigmas above median')
    p.add_option('--suffix', default='default', type='string',
                 help='postfix suffix to outfiles')
    p.add_option('--logfile', action='store', type='string', default='',
                 help='logfile name ')
    p.add_option('--verbose','-v', dest='verbose',default=False, 
                 action='store_true', help='toggle on verbose mode')


    (par, args) = p.parse_args()
    ss = sys.argv
    command = ''
    for ar in ss: command += ar+' '    # We want one string not a list.

    iList= args
    if len(iList) == 0:
       print 'options: ', par
       p.print_help()
       exit(0)
    # Generate an input list from idir plus an @ file or a unix template
    #
    inputList = getFileList(iList)
    ncmkflats(inputList,idir=par.idir, odir=par.odir, clobber=par.clobber, 
              suffix=par.suffix, sigma=par.sigma,logfile=par.logfile,
              verbose=par.verbose)

