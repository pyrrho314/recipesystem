#!/usr/bin/env python

import numpy as np
import datetime
import glob
import pyfits as pf
from niciTools import robust_sigma,getFileList,dmstod,parangle,rebin
from scipy import ndimage as nd
import math
#import pysao as sao
try:
    import stsci.numdisplay as ndis
except ImportError:
    import numdisplay as ndis
import matplotlib.mathtext as mathtext
import os, sys, imp
from  peak2halo import peak2halo

import time

def ncqlook(inputs, idir='', odir='./', log=True,lists=True, saturate=5000, nodisplay=False, full=False, port=5137):

    """
        nqlook
             Get NICI data from SBF: /net/petrohue/dataflow repository
             for last night date. Produces a quick display of each
             frame plus log file and lists that can be used with
             niciprepare, nicimkflats and niciscience scripts.
             The output lists are written in the working directory.

        nqlook inputs='20081116' saturate=2000 display=False 
             Get data from 20081116 with a saturation limit at 2000 
             and does not display (for speed)'

        nqlook @inlis idir='/data/myfiles' saturate=2000 display=False 
             Get data  from 'inlis' which is a list of input FITS files
             pathnames -one per line. The actual files are in the directory
             'idir'
             Use a saturation limit at 2000 and does not display (for speed)'
             Even though you have 'inputs' parameter in here, it is not
             taken into account since all the FITS files from the input
             directory will be read.


        NOTE: The files bad_r.fits ans bad_b.fits are read from 
              the nici directory wherever is it installed.

    """ 
    if (inputs == None):
        today = datetime.date.today()
        dd = today.timetuple()
        inputs = str(dd[0]) + '%.2d'%dd[1] + '%.2d'%dd[2]

    if nodisplay == False:
         # NOTE: pyds9 should be installed in ritchie. 
         # Needs to be part of the general installation and
         # an explanatory should be given in the installation Guide
   
         #ds9 = sao.ds9() 
         try:
            pp = 'inet:%d' % port
            ndis.open(imtdev=pp)
         except:
            print "\n ******************  WARNING:  Please start ds9."
            return
         #else:
         #   ndis.close()

    # We need to know were we are. 'nici' is the package name
    fp, pathname, description = imp.find_module('nici')
    pathname += '/'

            #loads the badpixels list for the blue and red channel
    bad_b = pf.getdata(pathname+'bad_b.fits')
    bad_r = pf.getdata(pathname+'bad_r.fits')

    # Chech that we can write in the current directory
    # the output FITS file is written here.
    odir = os.path.join(odir,'')         # Make sure odir has '/'
    if not os.access(odir,os.W_OK):
        print "\nERROR: You do not have write access to \
                this directory (output FITS file)\n"
        return

    if idir != '':
        idir = os.path.join(idir,'')         # Make sure odir has '/'
        if not os.access(idir,os.R_OK):
            print "\nERROR: No access to input directory:",idir
            return

    pathS = '/net/petrohue/dataflow/S'
    froot = 'ncqlook'

    if 'list' not in str(type(inputs)):
        date = str(inputs)
        fits_files = getFileList(date)
        # see if it date format  20yymmdd. Will not work for 20?.fits
        if len(fits_files) == 1:
           dd = fits_files[0] 
           if len(dd) == 8 and dd[:2] == '20':
               print 'Looking for files: ' +pathS+date+'S*.fits'
               fits_files = glob.glob(pathS+date+'S*.fits')
               froot = date
    else:
        # it is already a list:
        fits_files = inputs

    command = 'ncqlook ('+str(inputs)+',idir='+idir+',odir='+odir+',log=' \
             +str(log)+',lists='+str(lists)+',saturate='+str(saturate)\
             +',nodisplay='+str(nodisplay)+',full='+str(full)+')'

    fits_files = np.sort(fits_files)
    if len(fits_files)==0:
        print 'couldn''t find files for :'+pathS +inputs+'S*.fits'
        return

    print '\nfound '+str(len(fits_files))+ ' files.'

    previous_angle = -999 #define a few dummy variables
    previous_time = -999
    previous_ra   = -99
    previous_dec  = -999
    previous_mode = 'NONE'

    nfiles = 0
    nff = len(fits_files)
    cube = np.zeros((nff,256,512),dtype=np.float32)

    hh = 'NAME DATALAB EXPTIME CEN_WAVE FILTER_R FILTER_B (MIN-MAX rms)'+\
         ' [ext1,ext2 median] OBJECT CLASS MODE NCOADD CORE2HALO'
    print command
    print hh
    lg = 0
    if nff > 0 and (log or lists):              # Open outputfiles
       if log:
           lg = open(odir+froot+'.log','w+')
           lg.write('COMMAND: '+command+'\n')
           lg.write(hh+'\n')
       if lists:
           flis_1 = open(odir+froot+'.1_flats','w+')
           flis_2 = open(odir+froot+'.2_flats','w+')
           slis_adi = open(odir+froot+'.adi','w+')
           slis_sdi = open(odir+froot+'.sdi','w+')
           slis_asdi = open(odir+froot+'.asdi','w+')
           slis_other = open(odir+froot+'.other','w+')


    parser = mathtext.MathTextParser("Bitmap")
    for fitsname in fits_files:

        t1=time.time()
        fd = pf.open(idir+fitsname)              # open the Fits file
        fname = os.path.basename(fitsname)

        if len(fd) == 0:
            print fname," is EMPTY."
            continue
        phu = fd[0].header
        inst = phu['instrume']       # instrument name  for this image

        if inst != 'NICI':             #if not NICI, get the next filename
        #    print fname+' .. '+inst
            fd.close()
            continue

        if len(fd) != 3:
           line = ' Nici file '+fitsname+' does not have 2 extensions.'
           print line
           if log:
                lg.write(line+'\n')
           fd.close()
           continue

        nfiles += 1

        im_ext1 = fd[1].data
        hd1 = fd[1].header    # reads extension header 1
        im_ext2 = fd[2].data
        hd2 = fd[2].header    # reads extension header 1
            
        shift_status1=-1 
        shift_status2=-1
                              # extension 1 and blue channel in
                              # extension 2.
                              # -1 if we do not have red channel
                              # 0 is all good
                              # 1 is shifted rows

        if hd1.get('filter_b') == None and hd1.get('filter_r') != None:
            shift_status1 = check_missing_rows (im_ext1,bad_r,log,lg)
        if hd2.get('filter_b') != None and hd1.get('filter_r') == None:
            shift_status2 = check_missing_rows (im_ext2,bad_b,log,lg)
 
        devs = np.zeros(8, dtype=np.float32)     
             # Extracting the RMS of the
             # 8x8-folded image sample. A high value may 
             # indicate readout noise problems
        for sample in range(4):
            x1 = 200+600*(sample % 2)
            x2 = 263+600*(sample % 2)
            y1 = 200+600*(sample/2)
            y2 = 263+600*(sample/2)
            devs[sample] = np.std(fold_88(im_ext1[x1:x2,y1:y2]))

        for sample in range(4):
            x1 = 200+600*(sample % 2)
            x2 = 263+600*(sample % 2)
            y1 = 200+600*(sample/2)
            y2 = 263+600*(sample/2)
            devs[sample+4] = np.std(fold_88(im_ext2[x1:x2,y1:y2]))

        ncoadd_r = hd1['ncoadd_r']
        ncoadd_b = hd2['ncoadd_b']
        itime_r  = hd1['itime_r']
        crmode   = phu.get('crmode')

        if crmode == None:
            line = 'Header keyword CRMODE not found in '+fname+' -skipped.'
            print line
            if log:
                lg.write(line+'\n')
            fd.close()
            continue

        # Start forming the output line.
        obstype = str(phu.get('obstype'))
        oclass = str(phu.get('obsclass'))
        dichr = str(phu.get('dichroic')).lower()
        #adi = ''
        chan = ''
        if 'mirror' in dichr:
            chan = 'B'
        if 'open' in dichr:
            chan = 'R'

        mask = phu.get('fpmw')
        if 'clear' in mask.lower():
            mask = '-NoMASK'
        else:
            mask = ''

        tmp = obstype
        if obstype == 'FLAT':
            tmp = 'FLAT['+phu.get('gcalshut')+']'

        # 1-Channel (ADI)  (Dichroic is Mirror*(100% blue) or Open (100% red))
        # 2-Channel (SDI or ASDI)  (Dichroic is H-50/50*, H-CH4, H/K)
        mode = "OTHER"
        if (obstype == 'FLAT' or obstype == 'DARK'):
            if chan:
                mode = 'ADI-'+chan
            else:
                mode = 'SDI'
        elif (oclass == 'science'):
            if (crmode == 'FIXED'):
                if chan:
                    mode = 'ADI-'+chan+mask
                else:
                    mode = 'ASDI'+mask
            elif (crmode == 'FOLLOW'):  
                if chan:
                    mode = 'Normal'+chan+mask
                else:
                    mode = 'SDI'+mask
            elif '-ENG' in  phu.get('OBSID'):
                mode = '* ENGINEERING DATA *'
            else:
                mode = '***WARNING***: Not ADI, SDI nor ASDI-- '+crmode+' '+dichr

        line = ' '+str(phu.get('object'))+','+ oclass+','+tmp+','+mode

        ti_coad = ' '+str(itime_r)+','+str(ncoadd_r)
        dlab =' '+phu.get('DATALAB')
        exp_time =' %.2f'%(itime_r*ncoadd_r/60.)          # exposure time ... in minutes!
        cenwave =' %.3f'%(float(phu.get('WAVELENG'))/10000.)      # Central wavelenght in microns
        filter_r = ' %10.10s'%phu.get('FILTER_R')
        filter_b = ' %10.10s'%phu.get('FILTER_B')
        mindevs = ' (%8.2f,'%min(devs)
        maxdevs = '%8.2f)'%max(devs)
        medext1 = ' [%8.2f,'%np.median(im_ext1)
        medext2 = '%8.2f]'%np.median(im_ext2)
      
        longl = fname+dlab+exp_time+cenwave+filter_r+filter_b+mindevs+maxdevs+medext1+medext2+line+ti_coad

        if 'ENGINE' in mode:
            print longl
            if log:
                lg.write(longl+'\n')
            continue

        # hh = 'MIN-MAX rms 8x8 boxes: (min - max) MEDIAN: [ext1,ext2],(OBJECT),'\
        #     '(OBSCLASS),(OBSTYPE),(CRMODE),(DICHROIC),(FPMW),(ITIME),(NCOADD)' 
        # Don't print until core2halo value is defined below.
        # print longl
        
        if lists:
            if (obstype == 'FLAT' or obstype == 'DARK'):
                if chan:
                    flis_1.write(fname+'\n')
                else:
                    flis_2.write(fname+'\n')
            elif oclass == 'science':
                md = mode.lower()
                if 'mask' in md:
                    slis_other.write(fname+'\n')
                elif 'adi' in md:
                    slis_adi.write(fname+'\n')
                elif 'asdi' in md:
                    slis_asdi.write(fname+'\n')
                elif 'sdi' in md:
                    slis_sdi.write(fname+'\n')
                else:
                    slis_other.write(fname+'\n')
                
        # extast,hd1,astr ;extract the astrometric structure (astr) from the  
        #                   first extension header
        # cd=astr.cd ;extract the CD matrix from the astrometric structure,  
        #              this matrix defines scale/rotation/shear



        ra_pointing = hd1['crval1']                # pointing RA from astrometry
        dec_pointing = hd1['crval2']               # pointing DEC
        st = np.asfarray(phu['st'].split(':'))
        st = dmstod(phu['st'])                     # sidereal time
        HA = st - ra_pointing/15.                  # defining the hour angle
        harr = HA + np.asfarray([0,-1/60.])
        para_angle = parangle(harr, dec_pointing, -30.240750) 
              # computing the paralactic angle. -30.240750 is the 
              # latitude of Cerro Pachon

        obs_time = (dmstod(phu['LT'])*60.+720) % 1440
              # here time is expressed in minutes after noon

        angular_rate = para_angle[0] - para_angle[1] 
              # computing the rate of change of the paralactic angle
        current_angle = para_angle[0]
        astrom_angular_rate = (current_angle - previous_angle)/ \
                              (obs_time - previous_time)
              # angle change as found by taking the difference in 
              # angle between this frame and the previous one

        telescope_offset = np.sqrt( ((ra_pointing - previous_ra)*\
              math.cos(np.radians(dec_pointing)))**2+\
              (dec_pointing - previous_dec)**2 ) 

        if (crmode == 'FIXED') and (previous_mode == 'FIXED'):
            # we are in ADI mode 
            if (abs(angular_rate - astrom_angular_rate) > 1/10.) \
                and (telescope_offset < 1E-3) and \
                ((obs_time - previous_time) < 3*exp_time): 
                    print '---------- !!!!!!!! ---------------'
                    print 'Expected angular rate : '+str(angular_rate)+\
                          ' deg/minute'
                    print 'Astrometric angular rate : '+\
                          str(astrom_angular_rate)+' deg/minute'


        # CROP THE FRAMES and get core2halo ratios
        crop1=0; crop2=0; c2h_1=-1; c2h_2=-1

        if full:
            out_ext1 = rebin(im_ext1,256,256)
            out_ext2 = rebin(im_ext2,256,256)
        else:
            Ndev = lambda im : (np.amax(im) - np.median(im))/robust_sigma(im) 
            tmp = rebin(im_ext1,32,32)
            if Ndev(tmp) > 50:
                # We might have a peak:
                crop1 = 1
                c2h_1,x,y = peak2halo('',image=im_ext1)
                if x-128 < 0 or x+128 > 1024 or y-128 < 0 or y+128 > 1024:
                    out_ext1 = rebin(im_ext1,256,256)
                    crop1 = 0
                else:
                    #shift frame to center
                    out_ext1 = im_ext1[y-128:y+128,x-128:x+128]
            else:
                out_ext1 = rebin(im_ext1,256,256)

            tmp = rebin(im_ext2,32,32)
            if Ndev(tmp) > 50:
                # We might have a peak:
                crop2 = 1
                c2h_2,x,y = peak2halo('',image=im_ext2)
                if x-128 < 0 or x+128 > 1024 or y-128 < 0 or y+128 > 1024:
                    out_ext2 = rebin(im_ext2,256,256)
                    crop2 = 0
                else:
                    #shift frame to center
                    out_ext2 = im_ext2[y-128:y+128,x-128:x+128]
            else:
                out_ext2 = rebin(im_ext2,256,256)

        if 1:    # Do not process if there is no red flux
            a = out_ext1/ncoadd_r
            saturated_1 = a > saturate
            saturated_1[0:20,0:140] = 0
            sat_red = np.sum(saturated_1)

            #IDL out_ext1>=out_ext1[ (sort(out_ext1))[.01*n_elements(out_ext1)]]
            szx1 = np.size(out_ext1)
            sout = np.sort(out_ext1.ravel())
            limit = sout[.01*szx1]
            xmin = sout[0]; xmax= sout[np.size(sout)-1]

            out_ext1 = np.clip(out_ext1, limit, xmax)
                    # we put lower and upper limits and the 1% and 
                    # 99.5% (OR NCOADD*3500) of the image
            #IDL out_ext1<=( out_ext1[ (sort(out_ext1))[.995*n_elements(out_ext1)]]
            #           < (sxpar(hd1,'ncoadd_r')*saturate) )
            szx1 = np.size(out_ext1)
            if sout[.995*szx1] < ncoadd_r*saturate:
                limit = ncoadd_r*saturate
            else: 
                limit = sout[.995*szx1] 
         
            out_ext1 = np.clip(out_ext1, xmin, limit)

        if 1:
            a = out_ext2/ncoadd_b
            saturated_2 = nd.rotate(a > saturate,5,reshape=False)
            saturated_2[0:20,0:140] = 0
            sat_blue = np.sum(saturated_2)

            szx2 = np.size(out_ext2)
            sout = np.sort(out_ext2.ravel())
            limit = sout[.01*szx2]
            xmin = sout[0]; xmax= sout[np.size(sout)-1]

            out_ext2 = np.clip(out_ext2, limit, xmax)
                    # we put lower and upper limits and the 1% and
                    # 99.5% (OR NCOADD*3500) of the image
            #IDL: out_ext2<=( out_ext2[ (sort(out_ext2))[.995*n_elements(out_ext2)]]
            #           < (sxpar(hd2,'ncoadd_b')*saturate) )
            szx2 = np.size(out_ext2)
            if sout[.995*szx2] < ncoadd_b*saturate:
                limit = ncoadd_b*saturate
            else:
                limit = sout[.995*szx2]

            out_ext2 = np.clip(out_ext2, xmin, limit)

        out = np.zeros((256,512),dtype=np.float32)
              #create the image that will receive the frames to display
        out[:,:256] = ((out_ext1 - np.amin(out_ext1)) / \
                       (np.amax(out_ext1) - np.amin(out_ext1)))[::-1,:]
        a = out_ext2 - np.amin(out_ext2)
        b = np.amax(out_ext2) - np.amin(out_ext2)
        out[:,256:] = nd.rotate( (a/b), 5, reshape=False)[::-1,::-1]

        #extract the root name of the file. 
        name_out = os.path.splitext(os.path.basename(fname))[0]
        name_out += '['+('full','crop')[crop1]+'/'+('full','crop')[crop2]+']'
       
        # Form the text to overlay onto the lower left of the output array.
        rgba, depth1 = parser.to_rgba(name_out, color='gray', fontsize=10, dpi=200)
        c = rgba[::2,::2,3]
        far = np.asarray(c,dtype=np.float32)
        stx = np.shape(far)
        out[:stx[0],:stx[1]] = far[::-1,:]/256


        if not full:
            if sat_red != 0: 
                tmp1=out[:,:256] 
                tmp1[np.where(saturated_1)] = np.nan
                out[:,:256]=tmp1
        if not full:
            if sat_blue != 0:
                tmp2=out[:,256:]   
                tmp2[np.where(saturated_2)] = np.nan
                out[:,256:]=tmp2
        
 
        if nodisplay == False:
            ndis.display(out,z1=-1,z2=1,zscale=False)
            #ds9.view(out)           


        cube[nfiles-1,:,:] = out[::-1,:]

        # Now print 
        if oclass == 'science':
            if chan == 'red': c2h_2=-1
            if chan == 'blue': c2h_1=-1

            if c2h_1 > 0 and c2h_2 > 0:
               longl += ',['+'%.2f,%.2f' % (c2h_1,c2h_2)+']'
            elif c2h_1*c2h_2 == 1:   # both -1
               longl += ',[,]'
            elif c2h_1 > 0:
               longl += ',['+'%.2f,' % (c2h_1)+']'
            else:
               longl += ',['+',%.2f' % (c2h_2)+']'

        if sat_red or sat_blue:
           longl += ' Saturated('+str(sat_red)+','+str(sat_blue)+')'
        print longl
        if log:
            lg.write(longl+'\n')

        previous_angle=current_angle 
             #keep track of some variables for next iteration
        previous_time =obs_time
        previous_ra   =ra_pointing
        previous_dec  =dec_pointing
        previous_mode =crmode
        fd.close()

    if log:
       lg.close()
    if lists:
       flis_1.close()
       flis_2.close()
       slis_adi.close()
       slis_sdi.close()
       slis_asdi.close()
       slis_other.close()

    if (nfiles > 0): 
        pf.writeto(odir+froot+'_cube.fits',cube[:nfiles,::-1,:], clobber=True)
          # write the output cube with cropped/binned images
    else:
        print 'No NICI files found for the specified input: ',inputs
    return

def  check_missing_rows(im,bad,log,lg):
    imrv = im.ravel()
    if abs( np.median(imrv[bad])-np.median(imrv[bad+1]) ) < 50:
        line = 'image is affected by the infamous shifting rows'
        print line
        if log:
            lg.write(line+'\n')
        
        bad_shift=bad.copy()
        top   = np.where(bad > 512*1024.)
        bottom = np.where(bad <= 512*1024.)
        bad_shift[top   ]-=2048
        bad_shift[bottom]+=2048

        if abs( np.median(imrv[bad_shift])-np.median(imrv[bad_shift+1]) ) > 50:
            return 1 #shifted
        line =  'something strange with the file'
        print line
        if log:
            lg.write(line+'\n')
        #return -1
    else:
        return 0 #normal

def  fold_88(im):
     box = np.zeros((8,8), dtype=np.float32)
     sz=np.shape(im)
     y,x = xygen(sz[0],sz[1])
     for i in range(8):
         for j in range(8): 
             box[i,j] = np.median(im[np.where(((x%8) == i) & (y%8 == j))])
     return box

def xygen(ny,nx):
     y = np.arange(nx*ny).reshape(ny,nx)/ny - (ny/2 - 1)
     x = np.rot90 (y)
     return y,x

if __name__ == "__main__":

    import optparse

    # Parse input arguments
    usage = 'usage: %prog wild_card_or_date(yyyymmdd) [options]'

    p = optparse.OptionParser(usage=usage, version='niciprepare_1.0')
    p.add_option('--idir', action='store', type='string',
                 default='', help='Input directory pathname')
    p.add_option('--odir', action='store', type='string',
                 default='./', help='Output directory pathname')
    p.add_option('--saturate', action='store', type='int', default=5000,
                 help='saturate level')
    p.add_option('--nodisplay', action='store_true', default=False, 
                 help="Do not display running images")
    p.add_option('--log', action='store_true', default=True, 
                 help="Write a log file?")
    p.add_option('--lists', action='store_true', default=True, 
                 help="Create list files?")
    p.add_option('--full', action='store_true', dest="full", default=False, 
                 help="Display binned full frames")
    p.add_option('--port', action='store', type='int', default=5137, 
                 help="Port number")

    (par, args) = p.parse_args()

    #inputs= args
    if len(args) == 0:
        inputs=None
    else:
        inputs = args[0]

    ncqlook (inputs, idir=par.idir, odir=par.odir, saturate=par.saturate,nodisplay=par.nodisplay,
            log=par.log, lists=par.lists, full=par.full, port=par.port)

