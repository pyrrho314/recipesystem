*******
 xytran
*******

The xytran class calculates a set of coefficients to transform the blue frame pixels into the red frame coordinates system. This transformation is of the form:

::

 xr = a + b*xb + c*yb
 yr = d + e*xb + f*yb

To find the coefficients you need to display the Red and Blue frames from a pinhole NICI FITS file, click on control points which you can get by clicking on any point in the red and then the corresponding one in the blue frame. Here is the order of operations.

1. Python

 Start Python, ipython or Pyraf

2. from nici import * 

   Will load the set of nici scripts. Now create an xytran object: (NOTE: you need to have DS9 up first)

 ::
 
   # Creates xytran object
   cp = XYtran(pfiles='S20080812S0052.fits') 

   #The ds9 frame1 is red and frame2 is blue. 
 
   # If you already loaded the images and want to repeat a fit or
   # load only the control points and do the fit then:

   cp = XYtran()

3.  If you already have control points. The files name need to be: 'xyr' and 'xyb', 2 columns separated by blanks

 ::

   cp.loadcp()       # To load the control points

4. cp.getCpoints()     # mark control points

   Please position the cursor in one pinhole in the red frame and hit 'r', move to the corresponding pinhole in the blue frame and hit 'b', continue for at least 6 points. Hit 'q' the finish. 

5. cp.doCfit() # Do the fit

   A listing of the residuals is giving together with a line number. To remove a point from the list that shows a large residual use:

 ::

  cp.delxy(K)

  # Where K is the number in the list you want to remove. 
  # Now repeat 'cp.doCfit()'. A text file 'xycoeff' is created 
  # in your working directory with the coefficients (a..f). 

6. im_blue = cp.transform(im_blue)

   Transform a Blue frame into the coordinate system of the Red. im_blue is the variable containing the blue frame pixels. 


