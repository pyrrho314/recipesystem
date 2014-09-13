
ncscience. Analysis of science data 
===================================

**ncscience(inputs idir='' odir='' central=False suffix='default' 
bsize=5 mdfw=11 clobber=False logfile='' verbose=False)**

**Parameters**

* *inputs*
   The list of files used in **ncprepare**.  This list can a Unix wildcard pathname, e.g. * .fits, root23[2-9].fits, root??.fits or a @ list, e.g. @file.lis, where 'file.lis' is a text file with a list of FITS files, one file per line or a plain list of FITS filenames separated by commas. 

* *idir*
   Default is current directory. Directory pathname where the input files reside. 

* *odir*
   Default is current directory. Directory pathname to put the output FITS files. 

* *central*
   Default False. Use the whole frame size 1024x1024. If set to True it uses the central area (512x512). 

* *suffix*
   Dataset name. If 'default' it will take the rootname of the first element in the input list.

* *bsize* 
    Default value is 5. This is the boxcar smoothing box size in the *medfiltering* step. 

* *mdfw*
    Default value is 11. This is the median filtering width in the *medfiltering* step. 

* *clobber*
   Default value is False. Set to True to overwrite output files when they exist.

* *logfile*
    Log filename to hold the script messages. The default name is *gemini.log*

* *verbose*
    Default value is False. If True the information goes to the terminal as well.


**Description**

      Ncscience is a collection of python scripts to analyze the science files given in the parameter inputs and produces the following output files in the following order.

::

   GENERATE_CUBES:
      # The output files from 'ncprepare' are stacked into 2 cubes:
      cube_[red,blue].fits
          Stack up red and blue frames.
      medcrunch_[red,blue].fits
          Median reduce through the cube slices.
      sumcrunch_[red,blue].fits
          Sum reduce through the cube slices. The algorithm is: sum(cube)/sum(finite(cube)),
          where 'finite' selects only finite values.
          
   CUBE_ROTATE:
      # Since most of nici exposures are taken with the rotator off, the frames
      # are rotated from one another. To derotate the slices to a common zero
      # angle, we counter rotate with the value of the parallactic angle.

      # Output files are:
      cube_rotate_[red,blue].fits
          Rotated cube to common origin using the parallactic angles. 
      medcrunch_rotate_[red,blue].fits
          Median reduced through the slices of the cube_rotate. 
      sumcrunch_rotate_[red,blue].fits
          Sum reduced of the cube_rotate. (sum(cube)/sum(finite(cube))

   MEDIAN FILTERING:
      # Median filtering of cube slices. This image is the initial cube 
      # minus the median-smoothed image. This is sort-of-an-unsharp-mask 
      # but we use a median filtering and boxcar smoothing. The medfilt width 
      # and the boxcar size are parameters that the user can set. 

      # Input files: The cube_[red,blue].fits

      # Output files are:
      cube_medfilter_[red,blue].fits
          Medfiltered and boxcared subtracted initial cube.
      medfilter_medcrunch_rotate_[red,blue].fits
          Median reduced of rotated cube_medfilter. 
      cube_shift_medfilter_[red,blue].fits
          Scales the two channels to a common 'speckle' size. This is 
          done using the ratio of the central wavelengths of 
          the filter bandpasses. 
      cube_sdi.fits
          Differential imaging of [red-blue] slices from the cube_medfiltered.
      sdi_medcrunch.fits
          Median reduced of cube_sdi.fits 

   LOCI METHOD:
      # The Locally Optimized Combination of Images algorithm described by
      # Lafreniere et al. (ApJ 2007, 660) is used here to construct an estimate
      # of the PSF for every slice and minize speckle noise.
   
      loci_sdi:
          cube_loci_sdi.fits
              LOCI subtraction of cube_sdi plus rotate for each slice. 
          loci_sdi_medcrunch.fits
              Median reduced of the cube_loci_sdi. 
          loci_sdi_sumcrunch.fits
              Sum reduced of the cube_loci_sdi. 

      loci_medfilter:
          #LOCI substraction of the cube_medfilter_[red,blue].fits
          cube_loci_medfilter_[red,blue].fits
              Loci subtraction derotated cube_medfilter_[red,blue].fits
          loci_medfilter_medcrunch_[red,blue].fits
              Median reduced of cube_loci_medfilter. 
          loci_medfilter_sumcrunch.fits
              Sum reduced of cube_loci_medfilter. 

      loci_asdi:
          # LOCI substraction of the combination of sdi and adi methods on the
          # shift_medfiltered cubes.

          cube_asdi.fits
              'Super' Loci subtraction using the red channel as a cube
               and the blue channel as an additional sample of images. 
          asdi_medcrunch.fits
               Rotation and median reduced of the cube_asdi. 
          asdi_counter_medcrunch.fits
               Counter rotation and median reduced of the cube_asdi. 


**Examples**

1. ncscience nS20090312S00[1-3][0-9].fits --odir='/data' --suffix='NiciTest '    (Unix mode)

   Reduce all the matching FITS files . The the flats file located in the given directory and the output files will contain the string 'NiciTest '. 

2. ncscience @ncScience.lis idir='/data/20090312/' odir='/data/reduced' (Pyraf mode)

