.. _examples:

Examples
--------

.. _example1:

1) Obtain the pixel mapping (wavecal) to wavelength for a GMOS Long Slit image

 ::

  from astrodata import AstroData

  from wavecal import Wavecal

  ad = AstroData('gsS20110310S0137.fits')

  wc = Wavecal(ad)

  # Display parameters values used and fit function
  # coefficients.
  wc.info()

  # Compute the wavelengths from a list of pixel coordinates
  # from the middle row of the image.
  pixels = [200,300,400]

  print wc.z(pixels)

  # Compute the wavelengths from the image ARC peaks
  print wc.z(wc.xpeaks)

  # Plot the features found. 
  wc.plot_features()

  # Or list the features
  # Number, pix, z(pix), user (reference)
  wc.features()


