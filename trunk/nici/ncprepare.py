#!/usr/bin/env python
#
# COMMAND EXAMPLE:
#
#npp.py @in.lis --odir=/tmp/science --idir=/net/petrohue/dataflow --fsuffix=_S20090410 --fdir=/tmp/results/flats -v --clobber

import optparse

import os, sys
from os.path import join
import pyfits as pf
import scipy.ndimage as nd
import niciTools as nt
import numpy as np
import nici_cntrd as nc
import geminiLogger as glog
try:
    import stsci.numdisplay as ndis
except ImportError:
    import numdisplay as ndis

command = ''
def ncprepare (inputs, oprefix='n', idir='', odir='', fdir='', fsuffix='',
              dobadpix=True, clobber=False, logfile='', verbose=False):

    """
        inputs =               Input Fits files on a list
        (oprefix = 'n')        Prefix for output image(s)
        (idir = '')            Path for input raw files
        (odir = '')            Output directory name
        (fdir = '')            Path for calibration files
        (fsuffix = '')         Suffix used in calibration files 
        (dobadpix = True)      Correct bad pixels the best we can 
        (clobber = False)      Replace output file if exits
        (logfile = '')         Logfile name
        (verbose = False)      Display details on terminal?


         DESCRIPTION:
         Prepare science files in 'inputs' list for further analysis 
         with the ncscience script.
         - Shift the mask to frame center  
         - Map blue frame coordinates to the red frame's.

     EXAMPLES:
        PYRAF:
           ncprepare *.fits odir=/tmp

        Python_shell: 
           nicip.ncprepare(flis,idir='/data1/data/tmp',odir='/tmp')
      
        Unix_shell:
           ncprepare.py /data1/data/tmp/*.fits --odir='/tmp'

    """
    global command

    def xy_tran(out):
        """
          Function to transform blue coordinates to red frame coordinates
        """
        xref = out[1] - 990.5897
        yref = out[0] - 37.82109
        x = -0.9985259*xref - 0.0178984*yref
        y = -0.0181331*xref + 0.9991414*yref
        # if we want to recenter as well. Need to think about this.
        #x = x + (512 -409.81)
        #y = y - (512 -650.23)
        return y,x
    
    log = glog.getLogger(name='ncprepare',verbose=verbose,logfile=logfile,
          debug=True)
    log.info("\n START NICI Ncprepare +++++++++++++++++++++++++")
    log.info(command)
    log.info( 'inputs_list_len: '+str(len(inputs)))
    log.info( 'inputdir: '+idir)
    log.info( 'outdir: '+odir)
    log.info( 'flatsdir: '+fdir)

    try:
        odir = nt.check_dir(odir,os.W_OK)
        idir = nt.check_dir(idir,os.R_OK)
        fdir = nt.check_dir(fdir,os.R_OK)
    except IOError, strerror:
        log.warning('(%s)' % strerror)
        #sys.exit(0)
        return

    if (len(oprefix) > 10):
        log.warning( "input prefix is too long (>10 chars), returning")
        return
    
    if len(fsuffix): 
        if not fsuffix.startswith('_'): fsuffix = '_'+fsuffix

    try:
        flat_red  = pf.getdata(fdir + 'flat_red' +fsuffix+'.fits')
        flat_blue = pf.getdata(fdir + 'flat_blue'+fsuffix+'.fits')
        sky_red   = pf.getdata(fdir + 'sky_red'  +fsuffix+'.fits')
        sky_blue  = pf.getdata(fdir + 'sky_blue' +fsuffix+'.fits')
    except IOError, strerror:
        log.error( '(%s, please run ncmkflats.)' % strerror)
        return

    # Expand 'inputs' to a list of files
    inputs = nt.getFileList(inputs)

    ninlen=np.size(inputs) 
    for fname in inputs:
     
        fpath = join(idir,fname.strip())
        hdus = pf.open(fpath)
        if (hdus[0].header).get("INSTRUME") != 'NICI':
            log.info( 'File:'+fname+'is not for NICI. Ignored.' )
            hdus.close()
            continue

        if (hdus[0].header).get("PREPARE"):
            log.info( fpath+" already fixed using ncprepare. -skipped")
            continue

        if (len(hdus) != 3):
            log.info( 'File:'+fname+'does not have 2 image extensions. Ignored' )
            hdus.close()
            continue

        if (hdus[0].header['obsclass']).lower() != 'science':
            log.info( fname+' is not SCIENCE frame. -skipped')
            hdus.close()
            continue

        log.info('%s %s %s %s' % (fname,
                                  hdus[0].header['object'],
                                  hdus[0].header['obstype'], 
                                  hdus[0].header['obsclass']))

        # Get the root name and extension
        root,ext = os.path.basename(fpath).split('.')[0:2]

        outfile = odir+oprefix+root+'.'+ext
        if (os.access(outfile, os.R_OK)):
            if clobber:
               os.remove(outfile)
            else:
	       log.warning('File: '+outfile+ ' already exists.')
	       return

        # Create empty output file
        out = pf.open(outfile,mode='append')
        
        out.append(hdus[0]) 
        # Copy input header
        ophu = out[0].header
        
        date_time = nt.fits_utc()
    
        # Time Stamps
        ophu.update("GEM-TLM",date_time, "UT Last modification with GEMINI")
        ophu.update("PREPARE",date_time, "UT Time stamp for NPREPARE")
        out.flush()


        for ext in [1,2]:       # Do red and Blue

            im = hdus[ext].data       
            hdr = hdus[ext].header       

            im = np.asarray(im,dtype=np.float32)

            # sky subtract flat reduce
            if hdr.has_key('FILTER_R'):
                linearzR = hdr['ITIME_R']*hdr['NCOADD_R']
                #im = (im - sky_red) / flat_red / linearzR
                im = (im/linearzR - sky_red) / flat_red 
            else:
                linearzB = hdr['ITIME_B']*hdr['NCOADD_B']
                #im = (im - sky_blue) / flat_blue / linearzB
                im = (im/linearzB - sky_blue) / flat_blue 

            if dobadpix:               # Aggresive bad pixel removal
                im = remove_badpix(im,log)

            # Create a mask with the Nans.
            bp,im = nt.reset_nans(im)
    
            if ext == 2:
                # register blue frame to red's frame coordinate system
                im = nd.geometric_transform (im,xy_tran)
                #ndis.display(bp,frame=1,zscale=False)
                bp = nd.geometric_transform (bp,xy_tran)
                #ndis.display(bp,frame=2,zscale=False)

            # Get the center of the mask automatically. If not possible
            # nici_cntrd will display the image for the user to click

            xcen,ycen,im = nc.nici_cntrd(im,hdr,center_im=True)
            xco,yco,bp = nc.nici_cntrd(bp,hdr,center_im=True)

            # Now restore the Nans
            im = nt.restore_nans(im,bp)

            #print 'NUMBER OF bp>1',np.sum(bp>1.)
            #print 'NUMBER OF NANS',np.sum(np.isnan(im))

            log.info('%s %s: xcen: %.2f ycen: %.2f' % 
                                      (outfile,
                                      ['red','blue'][ext-1],
                                      xcen,ycen))  

            # Order wcs -in place.
            nt.order_wcs(hdr)
            hdr.update("EXTNAME",'SCI', "SCIENCE extension",after='GCOUNT')
            hdr.__delitem__("EXTVER")
            hdr.update("EXTVER",ext, "SCIENCE ext version",after='EXTNAME')
            pf.append(outfile,im,hdr)
            out.flush()

        out.close()

    return 

def remove_badpix(im,log):
    """
    Aggresive bad pixel removal. (Nan removal)
    If 3 or more of the neighbouring pixels are NOT NaNs, 
    we assume that the value for the NaN pixel is the
    mean of the non-NaNs. If not, leave it as a NaN

    """

    y,x = np.where(~np.isfinite(im))  # find bad pixels (NaNs)
    ny,nx = np.shape(im)

    # we will not attempt to correct badpixels at the edge of the image
    g = np.where((x==0)|(y==0)|(x==(nx-1))|(y==(ny-1)))
    y,x = remove_indx(g[0], y, x)

    im_corr=im.copy()   # im_corr is the 'corrected' image

    for nbad in range(np.size(x)):
        neighbours = im[y[nbad]-1:y[nbad]+2,x[nbad]-1:x[nbad]+2]
                        # if 3 or more of the neighbouring
                        # pixels are NOT NaNs, we assume that
                        # the value for the NaN pixel is the
                        # mean of the non-NaNs. If not, leave
                        # it as a NaN
        if np.sum(np.isfinite(neighbours)) >= 3:
            mm = np.mean(neighbours[np.where(np.isfinite(neighbours))])
            im_corr[y[nbad],x[nbad]] = mm
                        # Just for the fun (yeah right) of it,
                        # print the fraction of NaNs before
                        # and after this correction
    log.info('Badpix: Before and after correction:'+ '%5.2f%1c %5.2f%1c' % \
       (np.sum(np.isfinite(im)==0)/(nx*ny*1.)*100.,'%',\
       np.sum(np.isfinite(im_corr)==0)/(nx*ny*1.)*100.,'%') )
    return im_corr   # replace the image by the corrected one


def remove_indx(idx,x,x1=None,x2=None,x3=None,x4=None):
   """
     Remove from vectors x(i)
     those elements indexed by vector idx.

     The vectors must be of the same size.
   """

   nx = np.size(x)
   ni = np.size(idx)

   m=[]
   m.append(x)
   if x1 != None: m.append(x1)
   if x2 != None: m.append(x2)
   if x3 != None: m.append(x3)
   if x4 != None: m.append(x4)

   m=np.array(m)

   nelems = np.size(x)-np.size(idx)
   j=0
   k=0
   for i in range(nx):
       if idx[k] == i:
          k=min(k+1,ni-1)
          continue
       m[:,j] = m[:,i]
       j=j+1


   # copy to vectors
   x = m[0,:nelems]

   r,c=np.shape(m)
   return m[:r,:nelems]


if __name__ == '__main__':

    # Parse input arguments
    usage = 'Usage: \n       %prog file_list --odir=<output_directory> [options]'

    p = optparse.OptionParser(usage=usage, version='ncprepare_1.0')
    p.add_option('--oprefix', action='store', type='string', default='n', 
                 help='prefixes to use for prepared files')
    p.add_option('--idir', action='store', type='string',
                 default='', help='Input directory pathname')
    p.add_option('--odir', action='store', type='string',
                 default='', help='Output directory pathname')
    p.add_option('--fdir', action='store', type='string',
                 default='', help='directory pathname for calibration files')
    p.add_option('--fsuffix', action='store', type='string', default='n', 
                 help='suffix used in calibration filenames')
    p.add_option('--dobadpix', action='store_true',default=True,
                 help='Correct bad pixels the best we can')
    p.add_option('--clobber', default=False, action='store_true', 
                 help='Clobber output files?')
    p.add_option('--logfile', action='store', type='string', default='',
                 help='logfile name ')
    p.add_option('--verbose','-v', dest='verbose', action='store_true',default=False,
                 help='toggle on verbose mode')



    (par, args) = p.parse_args() 
    ss = sys.argv
    command=''
    for ar in ss: command += ar+' '    # We want one string not a list.

    iList= args 
    if len(iList) == 0:
       print 'options: ', par
       p.print_help()
       sys.exit(0)
    # Generate an input list from idir plus an @ file or a unix template 
    # 
    inputs = nt.getFileList(iList)
    ncprepare (inputs, oprefix=par.oprefix, idir=par.idir, odir=par.odir,
                 fdir=par.fdir, clobber=par.clobber, fsuffix=par.fsuffix,
                 dobadpix=par.dobadpix, logfile=par.logfile, verbose=par.verbose)

