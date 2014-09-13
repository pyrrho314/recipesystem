.. _linearize:

Resample an image to linear wavelength co-ordinates
---------------------------------------------------

Given that in general an arc image presents distortion, we want to capture this in a function and been able to correct with respect to a given point. 

After :ref:`fittting the arcs <imagefit>`, the method 'fit_image' also calculates an inverse function f(z,y,x) such that we can obtain a pixel value from a given pair (z,y). 

We then generate a set of lambdas with a dispertion value (cdelt = (self.z(nx) - self.z(1))/nx) as (lambdas = (ixx-crpix)*cdelt + crval), where 'ixx' is the array of indices (1..nx) along the dispersion axis.  With the inverse function we obtain the pixel coordinates corresponding to each lambda value. Interpolating the input image values at each of these new pixel coordinates using spline interpolation we linearize the input image.

.. _lin_example:


- Example: Resample an ARC image to linear wavelength co-ordinates

 ::

  from astrodata import AstroData

  from wavecal import Wavecal

  ad = AstroData('gsS20130526S0013.fits')

  # Create a Wavecal object with a GMOS LongSlit file.
  wc = Wavecal(ad)
 
  # Resample the image and output as AstroData object
  adout = wc.resample_image_asAstrodata()
 
.. _lin_sci_example:

- Example: Resample a SCIENCE data array using the resample_imageTo_LinearCoords()

 ::

  # Get the science image ndarray
  science_data = ad['SCI'].data
 
  # Linearize this
  out = wc.resample_imageTo_LinearCoords(science_data)
 
  # Create an AstroData object with the linearized image
  adout = wc.im_astrodata(out)
 
  # See that the SCI header contains the correct WCS
  print adout['SCI'].header.items

