.. _winter:

Winter: Manual line identification program.
-------------------------------------------

If the automatic line identification Wavecal fails to reach a reasonable line matching set, the Winter program can easily find a set. The user clicks on associates peaks in a subplot with reference lines and in a subplot with the input arc peaks.


- The program displays a reference ARC spectrum in wavelength units
  with the same scale as the input ARC spectrum's. This latter spectrum
  should have at least a rough aproximation of a WCS (crpix,crval and cdelt)
  and is displayed as a sesond subplot.

- Looking at the 2 plots the user can visually matched the arc lines, and
  easily make a selection of the correspoding lines in the reference and
  input spectrum. To associate, click a 'u' at the reference line peak
  and 'm' at the line peak in the input spectrum.


- The Winter available keystrokes are: (mouse pointer should be on the plot area)

 ::

  'h'  A list of the available keystrokes.
  'u'  Mark a reference line. Mouse pointer on top of the
       peak in the top subplot. If an "ERROR: Reference line
       not found. Try again." appears you should move the
       pointer closer to the peak. If the error accurs again
       then there is no reference line and you could choose
       another line.
  'm'  Mark the 'u' associated arc peak in the input spectrum
       on the bottom subplot.
  'l'  List the (pix, wavelength) accumulated so far.
  'f'  Fit the above list using a second order polynomial if
       the number of points is 3 or a 4th order Chebyshev polynomial
       if larger.
  'r'  Plot residuals in a new Figure. (Don't know how to delete yet)
  'd'  Delete poins in the residual plot (NOT yet working)
  'q'  To quit the loop.

The fitting function with the mapping is available in the 'z'
fit class object member. Submembers are:

 ::

  z(pix):    Returns the wavelength for pixel coordinates 'pix'.
  z.coeff:   List of fit function coefficients.
  z.fitname: Name of the fitting function used.
  z.order:   Fitting function order.

Example
--------

>>> import winter
>>> reload(winter);wi=winter.Winter(fname)
>>> wi.interactive()
>>> # A plot appears with 2 subplots and the user is in a loop
>>> # until a 'q' is hit.

