#!/usr/bin/env python
#
import optparse
import sys
import niciTools as nt

import pyfits as pf
from scipy import ndimage as nd
import scipy.signal
import os
from dist_circle import dist_circle
import numpy as np
import ginterpolate as gi
try:
    from stsci.convolve import boxcar
except ImportError:
    from convolve import boxcar

import time
import logging as log

#TODO   TODO
# the mask size needs to be a variable than can be set via a
# header keyword like ?masksize?
# BADPIX (Correct bad pixels) flag 1: only for extended object. 0: for planet finding

# Global variable to bring the command line to the __init__ function
command=''
def ncscience(inputs, idir='', odir='', central=False, suffix='default',
                 bsize=5, mdfw=11, clobber=False, logfile='', verbose=True):
 
    oc = NcReduce (inputs, idir, odir, central, suffix, bsize, mdfw, clobber,
                  logfile, verbose)

    oc.gen_cubes() 
    oc.write_cubes()
    oc.rot_cubes()
    oc.medfilter()

    oc.loci_sdi()
    oc.loci_medfilter()
    oc.loci_asdi()
    oc.epilog()

class NcReduce:
    """

     Reduce Science frames given as list of fits files.
     inputs: filenames in a python list
    """
    global command
    def __init__(self, inputs, idir='', odir='', central=False, suffix='default', 
                 bsize=5, mdfw=11, clobber=False, logfile='', verbose=True):
        """
        Check input parameters
        @param inputs: Input list file of the form '@in.lis', unix wilcard
                            or a comma separated filenames. Can be .fits or .fits.gz
                            file format.
        @type inputs:  string
        @param inputdir:    Input directory name to prepend to each member of the input list.
        @type inputdir:     string
        @param outputdir:   Directory name where the output files will be written.
        @type outputdir:    string
        @param central=False:
        @type central=False:
        @param suffix='Object':
        @type suffix='Object':
        @param clobber=True:
        @type clobber=True:
        @param logfile='':
        @type logfile='':
        @param verbose=True):

        """
        # Converting the verbose=True to =6 which is how it works in new log
        if verbose:
            verboseNew = 6
        else:
            verboseNew = 0

        log.basicConfig(filename='ncscience.log',level=log.DEBUG)

        self.log = log

        print 'ncc000',inputs
        # We want to write at least the strings below in the logfile,
        log.info("\n START NICI SCIENCE REDUCTION +++++++++++++++++++++++++")
        log.info(command)
        log.info( 'inputs_list_len: '+str(len(inputs)))
        log.info( 'inputdir: '+idir)
        log.info( 'outdir: '+odir)

        if suffix[0] != '_':
           suffix = '_'+suffix
        self.suffix = suffix
        
        try:
            odir = nt.check_dir(odir,os.W_OK)
            idir = nt.check_dir(idir,os.R_OK)
        except IOError, strerror:
            log.warning('(%s)' % strerror)
            sys.exit(0)
            

        self.odir = odir
        self.idir = idir
        self.central = central
        self.clobber = clobber
        self.boxcar_size  = bsize   # Boxcar smoothing size
        self.medfilt_width = mdfw   # Medfiltering width 

        if (not central):
            imsize = (1024,1024)
        else:
            imsize = (512,512)
        self.imsize = imsize

        if idir != '' and not isinstance(inputs, list):
            inputs = idir+inputs

        inputs = nt.getFileList(inputs)
        self.nfiles = len (inputs)

        if (self.nfiles == 0):
            log.error('Number of input files is zero.')
            return
 
        # Let's sort the input list. Usually files are of the form
        # S20091020Sdddd.fits, so it makes sence.
        self.inputs = np.sort(inputs)
        
        # Get the PHU of the 1st element in the list. It will be a 
        # generic header for the output products.
        # TODO: Maybe we need to add WCS info

        self.gphu = pf.getheader(idir+self.inputs[0],0)

        # Calculate the parallactic angles (they will be in order)

        self.pa = nt.parallacticAngle(self.inputs, idir)
        # define l1,l2
        self._get_filters_scaling(self.idir+self.inputs[0])


    def gen_cubes(self):
        self.cubes={}
        self._gen_cube('red')
        self._gen_cube('blue')
  
    def rot_cubes(self):
        if self._file_exists('red','cube_rotate'): return
        rcube = self._read_fits('red','cube')
        rcube = self._rotate_cube(rcube)
        self._write_fits(rcube,'red','cube_rotate')
        fs = lambda cube: np.sum(cube,axis=0) / np.sum(np.isfinite(cube),axis=0)
        self._write_fits(fs(rcube),'red','sumcrunch_rotate')
        self._write_fits(np.median(rcube,axis=0),'red','medcrunch_rotate')
        del rcube

        bcube = self._read_fits('blue','cube')
        bcube = self._rotate_cube(bcube)
        self._write_fits(bcube,'blue','cube_rotate')
        self._write_fits(fs(bcube),'blue','sumcrunch_rotate')
        self._write_fits(np.median(bcube,axis=0),'blue','medcrunch_rotate')
        del bcube
  
    def _gen_cube(self,chan):
        """
          Stack all input science frames into a cube . (cube_red and cube_blue)
          The images have been shifted to a common center with ncprepare. 
        """
         
        log = self.log
        
        if self._file_exists(chan,'cube'): return

        ext = {'red':1,'blue':2}[chan]
        self.cube = np.zeros((self.nfiles,)+self.imsize,dtype=np.float32)

        i = 0
        nfiles = self.nfiles
        for f in self.inputs:
           
            file = self.idir+f
            im = np.asfarray(pf.getdata(file,ext),dtype=np.float32)
        
            # Substract median
            medval = np.median(im[np.where(np.isfinite(im))],axis=None)
            im -= medval

            # Stack on a cube
            if (not self.central):
                self.cube[i,:,:] = im
            else:
                self.cube[i,:,:] = im[512-256:512+256,512-256:512+256]

            log.info('Staking %s[%s] %d/%d ' % (file,chan,i+1,nfiles))
            i += 1
        self.cubes[chan] = self.cube

    def _file_exists (self, chan, name):
        if chan: chan = '_'+chan
        file = self.odir + name + self.suffix +chan+'.fits'
        if os.access(file,os.R_OK) and not self.clobber:
            self.log.info( 'Access: '+file+' already exists.')
            return True
        else: return False

    def write_cubes(self):
   
        if self._file_exists('red','cube'): return

        rcube = self.cubes['red']
        bcube = self.cubes['blue']
        self._write_fits(rcube,'red','cube')
        self._write_fits(bcube,'blue','cube')
        fs = lambda cube: np.sum(cube,axis=0) / np.sum(np.isfinite(cube),axis=0)
        self._write_fits(fs(rcube),'red','sumcrunch')
        self._write_fits(fs(bcube),'blue','sumcrunch')
        self._write_fits(np.median(rcube,axis=0),'red','medcrunch')
        self._write_fits(np.median(bcube,axis=0),'blue','medcrunch')
        del self.cubes

    def _write_fits(self,cube,chan,name):
        """
          write cube 
        """
        log =     self.log
        clobber = self.clobber
        gphu =    self.gphu
        odir =    self.odir
        suffix =  self.suffix
        #cube =    self.cubes[chan]

        if chan: chan = '_'+chan
        file = odir + name + suffix +chan+'.fits'
        if os.access(file,os.R_OK) and self.clobber: os.unlink(file)
        log.info( 'Writing: '+file)
        pf.writeto(file,cube,header=gphu,clobber=clobber)

    def _read_fits(self,chan,name):
        """
          Read FITS
        """
        log =     self.log
        clobber = self.clobber
        gphu =    self.gphu
        odir =    self.odir
        suffix =  self.suffix
        #cube =    self.cubes[chan]

        if chan: chan = '_'+chan
        file = odir + name + suffix +chan+'.fits'
        log.info( 'Reading: '+file)
        data = pf.getdata(file)
        return data
  
    def medfilter(self):

        clobber = self.clobber
        gphu =    self.gphu
        odir =    self.odir
        suffix =  self.suffix

        log = self.log
        #self._get_filters_scaling(self.idir+self.inputs[0])

        if self._file_exists('','cube_sdi'): return

        l1 =      self.l1  
        l2 =      self.l2  

        rcube =  self._read_fits('red','cube')
        bcube =  self._read_fits('blue','cube')
        #rcube = self.cubes['red']
        #bcube = self.cubes['blue']

        sz = np.shape(rcube)
        if not self.central:               # r is a 2d image with pixel values
            r = dist_circle([1024,1024])   # corresponding to the distance to the
        else:                              # center of the image
            r = dist_circle([512,512])

        # Find pixels between 25 and 50 pixels from the center of the image.
        # These values should be modified when using a mask other than 0.32''
        g = np.where((r > 25) & (r < 50))
       
        # The masks data is:  Scale is 0.018 ''/pixels
        #Focal Plane Mask keyword name: FPMW
        # Values are:
        #    Clear
        #    0.90 arcsec
        #    0.65 arcsec
        #    0.46 arcsec
        #    0.32 arcsec radious    18 pixels
        #    0.22 arcsec
        #    Grid
        #    User


        rad = np.zeros(round(np.amax(r))-10,dtype=np.float32)

        cube_diff = np.zeros(np.shape(rcube),dtype=np.float32)

        # Scale rcube or bcube to the shorther wlen (depends on the set
        # of filter used)
        scale = min(l1,l2)/l1
        fn = 0
        if scale == 1.0:
            scale = min(l1,l2)/l2
            fn = 1
        
        ft = lambda i,o: '%.2f %d:%.2d' % (time.time()-i,(time.time()-o)/60,(time.time()-o)%60)
        to = time.time()
        for i in range(sz[0]):
            ti = time.time()
            rim = self._medfilter(rcube[i,:,:],r,rad)
            bim = self._medfilter(bcube[i,:,:],r,rad)
            log.info('medfilter:      ['+str(i+1)+'/'+str(sz[0])+'] ' + ft(ti,to))
            rcube[i,:,:] = rim
            bcube[i,:,:] = bim
            # scale r or b
            ti = time.time()
            sim = self._scale_im((rim,bim)[fn],scale)
            log.info('scale:                 ' +ft(ti,to))
            ti = time.time()
            cube_diff[i,:,:] = self._sdi(((sim,bim),(rim,sim))[fn],self.pa[i],g)
            log.info('sdi (red-blue):        ' +ft(ti,to))

        self._write_fits(cube_diff,'','cube_sdi')
        self._write_fits(np.median(cube_diff,axis=0),'','sdi_medcrunch')

        self._write_fits(rcube,'red','cube_medfilter')
        self._write_fits(bcube,'blue','cube_medfilter')

        self._shift_medfilter(rcube,'red')
        self._shift_medfilter(bcube,'blue')

        rcube = self._rotate_cube(rcube)
        bcube = self._rotate_cube(bcube)
        self._write_fits(np.median(rcube,axis=0),'red','medfilter_medcrunch_rotate')
        self._write_fits(np.median(bcube,axis=0),'blue','medfilter_medcrunch_rotate')

        # Release memory
        del rcube,bcube,cube_diff


    def _sdi(self,(rim,bim),angle,g):
         
        
        l1 = self.l1
        l2 = self.l2

        diff = np.zeros(np.shape(rim),dtype=np.float32)
        bim /= np.sum(rim[g]*bim[g])/np.sum(rim[g]**2)
                        # We project (scalar product)
                        # the pixels of the red frame onto
                        # the blue frame inside the annulus. The amplitude
                        # gives the optimal (in a least-square sense)
                        # amplitude of the blue frame to fit
                        # the red frame. The blue frame is scaled to this
                        # optimal amplitude

        if l2 < l1:
            diff[:,:] = bim - rim
        else:
            diff[:,:] = rim - bim
        
        bp,diff = nt.reset_nans(diff)
        diff = nd.rotate (diff,angle, reshape=False)
        bp = nd.rotate (bp,angle, reshape=False)
        #diff = nd.rotate (diff,-angle, reshape=False)
        #bp = nd.rotate (bp,-angle, reshape=False)
        diff = nt.restore_nans(diff,bp)

        return diff 

    def _medfilter(self,im,r,rad):
        """
           Median filtering of each slice in the cube
        """
        log =     self.log
        pa =      self.pa
        bsize =   self.boxcar_size     # Boxcar smoothing size (def:5)
        mdfw =    self.medfilt_width    # Medfiltering width  (def:11)

        # to accelerate the code, we bin the image and the 'r'
        # image down to 256x256 and do the radial profile on the binned versions
        
        sz = np.shape(im)
        bsy = sz[0]/4
        bsx = sz[1]/4
        rbin = nt.rebin(r,bsy,bsx)

        tmp = np.zeros(self.imsize,dtype=np.float32)

        # COPY: We are not modifying input image
        tmp[:,:] = im

        #NOTE: We are not using reset_num so far

        imbin = nt.rebin(tmp, bsy,bsx)   # bin down to 256x256

        for rr in range(8,np.size(rad)):
            # determine the median profile beyond r=8
            # (irrelevant inside this because of focal plane mask)
            if rr < 50:
                ss = tmp[abs(r-rr) < .5]
                rad[rr] = np.median(ss[np.where(np.isfinite(ss))])
            else:
                ss = imbin[abs(rbin-rr) < 3]
                rad[rr] = np.median(ss[np.where(np.isfinite(ss))])

        # smooth the radial profile beyond r=22
        rad[22:] = boxcar(rad,(bsize,))[22:]
        rad[40:] = boxcar(rad,(bsize,))[40:] # smooth even more

        tmp[:,:] -= gi.ginterpolate(rad,r)

        tmp[:,:] -= scipy.signal.medfilt2d(np.nan_to_num(tmp),mdfw)
        

        #t = time.time()
        #line = '%.2f' % (t-ti)
        #log.info('medfilter['+chan+']: ['+str(i+1)+'/'+str(sz[0])+']' +line)
        return tmp


    def _rotate_cube(self,cube):

        log =  self.log
        pa =   self.pa
        sz = np.shape(cube)
        # Now rotate this cube
        to=time.time()
        for k in range(sz[0]):
            ti=time.time()
            bp,cube[k,:,:] = nt.reset_nans(cube[k,:,:])
            cube[k,:,:] = nd.rotate (cube[k,:,:], pa[k], reshape=False)
            #cube[k,:,:] = nd.rotate (cube[k,:,:], -pa[k], reshape=False)

            #Now rotate the bp matrix 
            #bp = nd.rotate (bp, -pa[k], reshape=False)
            bp = nd.rotate (bp, pa[k], reshape=False)

            # And restore the Nans to cube
            cube[k,:,:] = nt.restore_nans(cube[k,:,:],bp)

            #_form_line('rotating ,k,nk,to,t1,t2,pa[k])
            t = time.time()
            line = 'rotating ['+str(k+1)+'/'+str(sz[0])+']'
            log.info(line+'%.2f %.2f %.2f' % (pa[k],t-ti,t-to))

        return cube
        #fname = self.odir+'medfilter_medcrunch_rotate'+self.suffix+'_'+chan+'.fits'
        #log.info( 'Writing: '+fname)
        #pf.writeto(fname,np.median(cube,axis=0),header=self.gphu,clobber=self.clobber)

    def _get_filters_scaling(self,nici_file):
        """
         Calculates the scaling from red and blue filters used
        """
        # List of filters
        hdr1 = pf.getheader(nici_file,1)
        filter_red = hdr1['filter_r'].split('_')[1]
        hdr2 = pf.getheader(nici_file,2) 
        filter_blue = hdr2['filter_b'].split('_')[1]
        # l1 list
        fr = {'G0724': 1.587, 'G0728': 1.603, 'G0720': 1.628, 'G0735': 1.628,
              'G0742': 1.578, 'G0740': 1.652, 'G0714': 1.701, 'G0710': 2.2718, 
            'G0708': 4.68, 'G0707': 3.78, 'G0706': 2.12, 'G0705': 2.15, 'G0704': 2.20,
                'Block': (-1),
                }
        self.l1 = fr[filter_red]
        # l2 list 
        fb = {'G0722': 1.587, 'G0726': 1.603, 'G0732': 1.628, 'G0743': 1.578,
                'G0737': 1.652, 'G0702': 1.25, 'G0703': 1.65, 'G0712': 1.644,
                'G0709': 2.1239, 'G0711': 2.1686, 'G0713': 1.596, 'Block': (-1),
                }
        self.l2 = fb[filter_blue]

            

    def _scale_im(self,im,scale):

        def _scale_1024(out,scale):
            x = (out[1]-512.)/scale + 512.
            y = (out[0]-512.)/scale + 512.
            return y,x
        def _scale_1600(out,scale):
            x = (out[1]-800.)/scale + 800.
            y = (out[0]-800.)/scale + 800.
            return y,x
        def _scale_512(out, scale):
            x = (out[1]-256.)/scale + 256.
            y = (out[0]-256.)/scale + 256.
            return y,x

        im = np.asarray(im,dtype=np.float32) 
        bp,im = nt.reset_nans(im)
        if (self.central):
            im = nd.geometric_transform(im,_scale_512,extra_arguments=(scale,)) 
            bp = nd.geometric_transform(bp,_scale_512,extra_arguments=(scale,)) 
        else:
            im = nd.geometric_transform(im,_scale_1024,extra_arguments=(scale,))
            bp = nd.geometric_transform(bp,_scale_1024,extra_arguments=(scale,))
        im = nt.restore_nans(im,bp)
        return im

    def _shift_medfilter(self,cube,chan):
        """
          Scaled each sliced of the medfilter_blue file to the scale
          of the corresponding red.
        """

        log =     self.log
        clobber = self.clobber
        gphu =    self.gphu
        odir =    self.odir
        suffix =  self.suffix
        pa =      self.pa
        #cube =    self.cube
        l1 =      self.l1
        l2 =      self.l2
        
        if chan == 'blue':
            scale = min(l1,l2)/l2
        else:
            scale = min(l1,l2)/l1

        # scale:  this is a parameter that depends on
        #        the exact set of filters used

        to = time.time()
        sz = np.shape(cube)
        if scale != 1.0:
            #fname = odir+'cube_medfilter'+ suffix + '_'+ chan +'.fits'
            ## We read the file because the array self.cube has been
            ## rotated
            data = cube.copy()       # Make a local copy
            for i in range(sz[0]):
                data[i,:,:] = self._scale_im(data[i,:,:],scale)
                tt = time.time() - to
                log.info('shift_medfilter['+chan+']: [%d/%d] %d:%.2d' % (i+1,sz[0],tt/60,tt%60))
            self._write_fits(data,chan,'cube_shift_medfilter')
        else:
            self._write_fits(cube,chan,'cube_shift_medfilter')

        #return odir+'cube_shift_medfilter'+suffix+'_'+chan+'.fits'
        return



    def loci_sdi(self):
        log=self.log

        ft = lambda i: '%d:%.2d' % ((time.time()-i)/60,(time.time()-i)%60)
        if self._file_exists('','cube_loci_sdi'): return
        pa = self.pa
        cube =  self._read_fits('','cube_sdi')
        #cube = self.cube_diff
        sz = np.shape(cube)
        bp,cube = nt.reset_nans(cube)
        log.info('Applying loci_subtract...(Can take a while......)')
        ti = time.time()
        cube = loci_subtract(cube, self.pa)
        #bp = loci_subtract(bp, self.pa)
        log.info('                              '+ft(ti))

        for i in range(sz[0]):
            log.info('rotating cube_diff['+str(i+1)+'/'+str(sz[0])+'] %.2f'%-pa[i])
            cube[i,:,:] = nd.rotate (cube[i,:,:],pa[i], reshape=False)
            bp[i,:,:] = nd.rotate (bp[i,:,:],pa[i], reshape=False)
            #cube[i,:,:] = nd.rotate (cube[i,:,:],-pa[i], reshape=False)
            #bp[i,:,:] = nd.rotate (bp[i,:,:],-pa[i], reshape=False)

        cube = nt.restore_nans(cube,bp)
        self._write_fits(cube,'','cube_loci_sdi')
        self._write_fits(np.median(cube,axis=0),'','loci_sdi_medcrunch')
        fs = lambda cube: np.sum(cube,axis=0) / np.sum(np.isfinite(cube),axis=0)
        self._write_fits(fs(cube),'','loci_sdi_sumcrunch')

    def loci_medfilter(self):
        log = self.log

        if self._file_exists('red','cube_loci_medfilter'): return
        ft = lambda i: '%d:%.2d' % ((time.time()-i)/60,(time.time()-i)%60)
        pa = self.pa
        for chan in ['red','blue']:
            cube =  self._read_fits(chan,'cube_medfilter')

            bp,cube = nt.reset_nans(cube)
            log.info('Applying loci method ...(Can take a while......)')
            ti = time.time()
            cube = loci_subtract(cube, pa)
            log.info('                             '+ft(ti))

            sz = np.shape(cube)
            for i in range(sz[0]):
                log.info('rotating cube_medfilter['+str(i+1)+'/'+str(sz[0])+'] %.2f'%-pa[i])
                cube[i,:,:] = nd.rotate (cube[i,:,:],pa[i], reshape=False)
                bp[i,:,:] = nd.rotate (bp[i,:,:],pa[i], reshape=False)
                #cube[i,:,:] = nd.rotate (cube[i,:,:],-pa[i], reshape=False)
                #bp[i,:,:] = nd.rotate (bp[i,:,:],-pa[i], reshape=False)
            cube = nt.restore_nans(cube,bp)

            self._write_fits(cube,chan,'cube_loci_medfilter')
            self._write_fits(np.median(cube,axis=0),chan,'loci_medfilter_medcrunch')
            fs = lambda cube: np.sum(cube,axis=0) / np.sum(np.isfinite(cube),axis=0)
            self._write_fits(fs(cube),chan,'loci_medfilter_sumcrunch')
            del cube,bp

    def loci_asdi(self):
        """
          Combine sdi and adi method on the shift_medfiltered cubes
        """
        log = self.log

        if self._file_exists('','asdi_medcrunch'): return

        ft = lambda i: '%d:%.2d' % ((time.time()-i)/60,(time.time()-i)%60)
        pa = self.pa
        # in case l1,l2 are not defined.
        #self._get_filters_scaling(self.idir+self.inputs[0])
        l1 =      self.l1
        l2 =      self.l2
        
        cube_red = self._read_fits(('blue','red')[l1 < l2],'cube_shift_medfilter')
        cube_blue = self._read_fits(('red','blue')[l1 < l2],'cube_shift_medfilter')

        # Lets not use bpr because we run out of memory
        #bpr,cube_red = nt.reset_nans(cube_red)
        cube_red = np.nan_to_num(cube_red)
        log.info('Applying loci method ...(Can take a while......)')
        ti = time.time()
        # 
        cube_red = loci_subtract(cube_red, pa, np.nan_to_num(cube_blue))
        log.info('                             '+ft(ti))
        del cube_blue

        #cube_red = nt.restore_nans(cube_red,bpr)
        self._write_fits(cube_red,'','cube_asdi')

        #bpr,cube_red = nt.reset_nans(cube_red)
        line = 'rotating combined sdi,adi cube'
        cube_red2=cube_red.copy()
        sz = np.shape(cube_red2)
        for i in range(sz[0]):
            log.info(line+'['+str(i+1)+'/'+str(sz[0])+']')
            cube_red[i,:,:] = nd.rotate (cube_red[i,:,:], pa[i], reshape=False)
            #cube_red[i,:,:] = nd.rotate (cube_red[i,:,:], -pa[i], reshape=False)

        #cube_red = nt.restore_nans(cube_red,bpr)
        self._write_fits(np.median(cube_red,axis=0),'','asdi_medcrunch')
        del cube_red

        for i in range(sz[0]):
            #Do counter rotate (only noise)
            cube_red2[i,:,:] = nd.rotate (cube_red2[i,:,:], -pa[i], reshape=False)
            #cube_red2[i,:,:] = nd.rotate (cube_red2[i,:,:], pa[i], reshape=False)

        self._write_fits(np.median(cube_red2,axis=0),'','asdi_counter_medcrunch')
        del cube_red2
 
    def epilog(self):
        log = self.log
        pa = self.pa
        log.info('total angular range explored : %.2f (%.2f to %.2f)' % (max(pa)-min(pa),min(pa),max(pa)))
    
import time
import os
import copy

def loci_subtract(cube,PA,more_images=None):
   """
	PA: in degrees. Paralactic angle  
 
        Locally Optimized Combination of Images algorithm.
        See Lafreniere et all (2007b) APJ
   """
   PA -= np.amin(PA)
   XtraImages = (more_images is not None)
   cube2 = cube.copy()             # cube2 will contain reference PSF as 
                            # built from the ADI technique
   radeg = 180./np.pi
   sz = np.shape(cube)
   Nslice = sz[0]         # Dimension order is (z,y,x)  
   Nslice_total=Nslice
   if XtraImages:
       Nslice_more_images=(np.shape(more_images))[0]
       PA = np.concatenate((PA,np.zeros(Nslice_more_images)))
       Nslice_total += Nslice_more_images

   Npix_included = 5        # The PSF must move by more than Npix_included 
	                    # to be included in the ADI matrix
                            
                            # Generates an x/y matrix with zero at the center
   # this can be done with meshgrid
   nx = sz[1]
   ny = sz[2]
   y = np.arange(nx*ny).reshape(ny,nx)/nx - nx/2. + 0.5
   x = np.rot90(y) 

   theta = np.arange(nx*ny,dtype=float).reshape(nx,ny)
   theta = np.arctan2(x,y)*radeg

   r = np.sqrt(x**2 + y**2)          # Distance from the center of PSF
   theta += 180                   # just to start at 0

   NRINGS = 11
   NSECTORS = 6
   RINGWD = 30
   for RR in range(NRINGS):          # loop on radial distance of arcs
       r1 = 8 + RR*RINGWD
       r2 = 8 + RINGWD + RR*RINGWD
       Ntheta = np.pi*(r2**2 - r1**2)/Nslice_total
       Ntheta = np.ceil(Ntheta/100.) # we want at least 100 times more pixels 
                                   # than degrees of freedom,
       Ntheta = min(Ntheta,6)      # but no more than 6 sections

       for TT in range(Ntheta):     # loop on angle position of the arcs
          
          theta1 = (360./Ntheta)*TT   # defining radius and angle boundaries 
        			      # of arcs
          theta2 = (360./Ntheta)*(TT+1)

                     # g is the index of pixel within the arc
          g = np.where((r>=r1) & (r<r2) & (theta>=theta1) & (theta<theta2)) 
          ng = np.size(g)/2    # is an (y,x) array of indexes (2)
                     # subsection is a 2D matrix containing in one axis 
                     # the pixels of the arc and second axis the i
                     # Nth slice of the cube
          if ng == 0: break
          subsection = np.empty((Nslice_total,ng),dtype=float)
                     # filling the subsection matrix
          for k in range(Nslice):
              subsection[k,:] = cube[k,:,:][g]
          if XtraImages:
              for i in range(Nslice_more_images):
                  subsection[i+Nslice,:] = more_images[i,:,:][g]

		      # defining a MM, the matrix containing the linear 
                      # projection of Nth slice versus Mth slice of the cube
          MM = np.mat(np.empty((Nslice_total,Nslice_total),dtype=float))
          for ii in range(Nslice_total):
                      # no need to loop from 0 as the matrix is symetrical
              for jj in range(ii,Nslice_total):
                  MM[ii,jj] = np.nansum(subsection[ii,:]*subsection[jj,:])
                  MM[jj,ii] = MM[ii,jj]    
                      
                      # loop on the slices of the cube --NOT the slice 
                      # of the cube with the slices added in more_images
          for nth in range(Nslice):
	      pdiff = abs((PA[nth]-PA)/radeg * (r1+r2)/2 )
                      # defining which slices are included in the 
                      # building of the reference PSF of the Nth slice
              included = np.where( (pdiff >= Npix_included) | (PA == -1))[0]  
                                               # [0] because where
                                               # will produce an array[1,:] 
                                               # rather than [:,]
	      Nincluded = np.size(included)
              if Nincluded == 0: 
                 #print 'size included is zero for nth slice:',nth
                 continue
                      # (MM[included,Nth]) ;amplitude vector of the 
                      # projection of all included slices on the Nth slice
	      amps = (MM[:,included][included,:]).I *(MM[included,nth])
              amps=np.array(amps).flatten()

		      # a 2D image where the psf estimate we are
                      # building will be inserted with the pixel ordered
                      # back to their original position
              delta  = np.empty((nx,ny),dtype=float)
              cdel = delta[g].flatten()
	      for i in range(Nincluded):
	          cdel += amps[i]*subsection[included[i],:]
	          #delta[g] += amps[i]*subsection[included[i],:]
              delta[g] = cdel
              
              cube2[nth,:,:] -= delta
       #print 'adi subtraction R=',r1,'to',r2,'angle=',theta1,'to',theta2

   #-- this's not really ADI ---
   med = np.median(cube,axis=0)   # Median over all the slices
   far = np.where(r >= r2) # in the region outside of the largest arc, 
                   # we put in the cube the median of all slices
   for i in range(Nslice):
       tmp = copy.deepcopy(cube2[i,:,:])
       tmp[far] -= med[far]
       cube2[i,:,:] = tmp

   #----------------------------
   cube = cube2
   return cube

    
## To write TOTAL number of degrees explored


if __name__ == '__main__':
    # Parse input arguments
    usage = 'usage: %prog inputs [options]'

    p = optparse.OptionParser(usage=usage, version='ncscience.0')
    p.add_option('--idir', action='store', type='string',
                 default='', help='Input directory pathname')
    p.add_option('--odir', action='store', type='string',
                 default='', help='Output directory pathname')
    p.add_option('--central','-c', action='store_true',default=False,
                 help='Choose 512x512 output images')
    p.add_option('--suffix', action='store', type='string',
                 default='Object', help='Data set name')
    p.add_option('--bsize', action='store', type='int',
                 default=5, help='Boxcar smoothing size')
    p.add_option('--mdfw', action='store', type='int',
                 default=11, help='Medfiltering width')
    p.add_option('--clobber', default=False, action='store_true',
                 help='Clobber output files?')
    p.add_option('--logfile', action='store', type='string', default='',
                 help='logfile name ')
    p.add_option('--debug','-d', action='store_true',default=False,
                 help='toggle debug messages')
    p.add_option('--verbose','-v', dest='verbose', action='store_true',default=False,
                 help='toggle on verbose mode')

    (par, args) = p.parse_args()
    ss = sys.argv
    command=''
    for ar in ss: command += ar+' '    # We want one string not a list.
    if par.debug:
        par.verbose = True
        print 'options: ', par
        print 'args: ', args
    iList= args
    if len(iList) == 0:
       print 'options: ', par
       p.print_help()
       sys.exit(0)
    # Generate an input list from inputdir plus an @ file or a unix template
    #
    inputs = nt.getFileList(iList)

    print 'ncc0:',iList
    print 'ncc1:',inputs
    ncscience (inputs, idir=par.idir, odir=par.odir, central=par.central,
               suffix=par.suffix, bsize=par.bsize, mdfw=par.mdfw, clobber=par.clobber,
               logfile=par.logfile, verbose=par.verbose)

    #oc = NcReduce (inputs, idir=par.idir,odir=par.odir, central=par.central, 
    #              suffix=par.suffix, bsize=par.bsize, mdfw=par.mdfw, clobber=par.clobber,
    #              logfile=par.logfile, verbose=par.verbose)

    #oc.gen_cubes() 
    #oc.write_cubes()
    #oc.rot_cubes()
    #oc.medfilter()

    #oc.loci_sdi()
    #oc.loci_medfilter()
    #oc.loci_asdi()
    #$$$$$$$$$$$ commenting out next line as it is in above fn ncsciene where oc is actually declared/defined
    #oc.epilog()
    #$$$$$$$$
