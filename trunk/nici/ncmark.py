#! /usr/bin/env python
import numpy as np
from niciTools import congrid,robust_sigma,gcentroid, getFileList, check_dir
import scipy.ndimage as nd
from peak2halo import peak2halo
from scipy.signal import medfilt2d
import geminiLogger as glog
import pyfits as pf

try:
    import stsci.numdisplay as ndis
except ImportError:
    import numdisplay as ndis
import os
from os.path import join
import sys
import optparse

def ncmark(inputs, idir='', odir='./', oprefix='m', port=5137, logfile='', clobber=False, verbose=False):
    """
    Mark by hand the mask center and update the input header
    with (xcen,ycen) on the output file.
    The SAOIMAGE/ds9 display must be up. If the default port is busy, then it will
    use port 5199, so make sure you start "ds9 -port 5199".
    """

    log = glog.getLogger(name='ncmark',verbose=verbose,logfile=logfile,
          debug=True)
    log.info("\n START NICI Ncprepare +++++++++++++++++++++++++")
    #log.info(command)
    log.info( 'inputs_list_len: '+str(len(inputs)))
    log.info( 'inputdir: '+idir)
    log.info( 'outdir: '+odir)

    try:
        odir = check_dir(odir,os.W_OK)
        idir = check_dir(idir,os.R_OK)
    except IOError, strerror:
        log.warning('(%s)' % strerror)
        #sys.exit(0)
        return


    # Expand 'inputs' to a list of files
    inputs = getFileList(inputs)

    ninlen=np.size(inputs)
    for file in inputs:

        fpath = join(idir,file.strip())
        hdus = pf.open(fpath)
        if (hdus[0].header).get("INSTRUME") != 'NICI':
            log.info( 'File:'+file+'is not for NICI. Ignored.' )
            hdus.close()
            continue

        if (hdus[0].header['obsclass']).lower() != 'science':
            log.info( file+' is not SCIENCE frame. -skipped')
            hdus.close()
            continue

        log.info('%s %s %s %s' % (file,
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
        out.flush()

        for ext in [1,2]:       # Do red and Blue
            hdr =  hdus[ext].header
            im  =  hdus[ext].data
            _findxycenter(hdr, im, ext,log)
            pf.append(outfile,im, hdr)
            out.flush()

def _findxycenter(hdr, im,ext,log):

    def xy_tran(out):
        """
          Function to transform blue coordinates to red frame coordinates
        """
        xref = out[1] - 990.5897
        yref = out[0] - 37.82109
        x = -0.9985259*xref - 0.0178984*yref
        y = -0.0181331*xref + 0.9991414*yref
        return y,x
    
    #hdr = pf.getheader(file)
    #im = pf.getdata(file)
    xcen = hdr.get('xcen')
    ycen = hdr.get('ycen')
    if (xcen == None):
        #ratio,xc,yc = peak2halo('',image=im)
        #xcen = xc
        #ycen = yc
        #if (xcen < 0 or ycen < 0):
        if True:
            if ext == 2:
                # register blue frame to red's frame coordinate system
                im = nd.geometric_transform (im,xy_tran)


            try:
                ndis.display(im,zscale=False)
            except IOError,err:
                sys.stderr.write('\n ***** ERROR: %s Start DS9.\n' % str(err))
                sys.exit(1)
            print " Mark center with left button, then use 'q' to continue, 's' to skip."
            cursor = ndis.readcursor(sample=0)
            cursor = cursor.split()
            if cursor[3] == 's':
                hdr.update("XCEN",-1, "Start mask x-center")
                hdr.update("YCEN",-1, "Start mask y-center")
                updated = True 
                print '\nFrame skipped... ****Make sure not to use it in your science script. ***\n'
                #return updated,xcen,ycen,im
                return xcen,ycen,im
            x1 = float(cursor[0])
            y1 = float(cursor[1])

            xcen,ycen = gcentroid(im, x1,y1,9)
            xcen = xcen[0]
            ycen = ycen[0]
 
            hdr.update("XCEN",xcen, "Start mask x-center")
            hdr.update("YCEN",ycen, "Start mask y-center")
            log.info('                    (%.2f,%.2f) :: ' % (xcen,ycen)+('red','blue')[ext-1])

    return 

if __name__ == '__main__':

    # Parse input arguments
    usage = 'Usage: \n       %prog file_list --odir=<output_directory> [options]'

    p = optparse.OptionParser(usage=usage, version='ncprepare_1.0')
    p.add_option('--idir', action='store', type='string',
                 default='', help='Input directory pathname')
    p.add_option('--odir', action='store', type='string',
                 default='', help='Output directory pathname')
    p.add_option('--oprefix', action='store', type='string',
                 default='m', help='Prefix for output files.')
    p.add_option('--clobber', default=False, action='store_true', 
                 help='Clobber output files?')
    p.add_option('--port', action='store', type='int', default=5137,
                 help="Port number")
    p.add_option('--logfile', action='store', type='string', default='',
                 help='logfile name ')
    p.add_option('--verbose','-v', dest='verbose', action='store_true',default=False,
                 help='toggle on verbose mode')



    (par, args) = p.parse_args() 

    iList= args 
    if len(iList) == 0:
       print 'options: ', par
       p.print_help()
       sys.exit(0)
    # Generate an input list from idir plus an @ file or a unix template 
    # 
    inputs = getFileList(iList)
    ncmark (inputs, idir=par.idir, odir=par.odir, oprefix=par.oprefix, 
                 clobber=par.clobber, port=par.port,  
                 logfile=par.logfile, verbose=par.verbose)

