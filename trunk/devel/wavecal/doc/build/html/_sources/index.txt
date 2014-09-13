.. WAVECAL Data Reduction documentation master file, created by
   sphinx-quickstart on Sun Oct  4 20:10:32 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Wavelength Calibration 
============================

Wavelength Calibration finds automatically a mapping from pixel coordinates to wavelengths given an ARC image from a reduced Gemini FITS file. The mapping is possible 
only if you have a reference list containing known wavelengths for the 
spectral lines in the arc image and a rough WCS values in the header. 

The main Wavecal class methods are:

- **wavecal** Finds the mapping between pixel coordinates and wavelengths using a Matching Ratio Algorithm.
- **pix2ref** Given a list of pixel coordinates calculates the associated wavelengths.
- **fit_image** Fit a 3D function f(pix, row, wavelength) to the input image containing an ARC spectrum.
- **im_linearize** Linearize the image applying the surface function evaluator from the method fit_image.
- **im_astrodata** Generate an AstroData object given and ndarray and the input FITS PHU and HDU[‘SCI’].
- **winter** A Winter class method to interactively find a solution.

.. toctree::
   :maxdepth: 2
   :hidden:

   todo
   availability
   inputData
   dataRequirements
   quickExample
   imagefit
   linearize
   examples
   functionality
   wavecal
   winter
   recipesPrimitives

- :ref:`MISSING FUNCTIONALITY <todo>`
- :ref:`Accesing Wavecal <availability>`
- :ref:`Input Data   <inputData>`
- :ref:`Input Data Requirements  <dataRequirements>`
- :ref:`Quick Example   <quickExample>`
- :ref:`Fitting the ARC image  <imagefit>`
- :ref:`Linearize an image  <linearize>`
- :ref:`Primitives and Recipes  <prims>`
- :ref:`Examples <examples>`
 - :ref:`Example 1. Obtain the pixel mapping (wavecal) to wavelength for a GMOS Long Slit image <example1>`
 - :ref:`Example 2. Compute a fitting function to the ARC image. <imagefit>`
 - :ref:`Example 3. Linearize an ARC image. <lin_example>`
 - :ref:`Example 4. Linearize a SCIENCE image. <lin_sci_example>`
 - :ref:`Example 5. If Wavecal fails to converge: Interactive  <winter>`
- :ref:`Work flow.  <functionality>`
- :ref:`Wavecal class <wcal>`
- :ref:`Instrument class <instrument>`
- :ref:`GMOS class <gmos>`
- :ref:`F2 class <f2>`
- :ref:`Gnirs class <gnirs>`
- :ref:`NIRI class <niri>`
- :ref:`NIFS class <nifs>`
- :ref:`MatchLines class <matchlines>`
- :ref:`Winter <winter>`: Interactive line identification.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

