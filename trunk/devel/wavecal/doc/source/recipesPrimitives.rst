.. _prims:

Recipes
-------

Wavecal functionality is available through primitives function names in a text file
with names 'recipe.<some_name>'

Primitives
----------

There are three primitives: 

- determineWaveCal
   Given a supported Arc Lamp exposure, find the Wavelength Calibration and write the fitting function into a BINTABLE which is appended to the input file with the EXTNAME 'WAVECAL'. The input parameters available via the **reduce** command are:

  - linelist
  - fitfunction
  - fitorder
  - match
  - nsum
  - ntmax
  - minsep
  - debug
 
- approximateWaveCal
   Given a supported ARC exposure finds an approximate wavelength value for the center coordinates updating the WCS.

- wcalResampleToLinearCoords
   Resample the input image to a linear set of coordinates by applying the surface fitting function available in the BINTABLE 'WAVECAL'
