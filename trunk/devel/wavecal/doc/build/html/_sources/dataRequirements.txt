.. _dataRequirements:

Requirements
------------

All ARC spectra should be properly prepared with the SCI header containing a rough estimation of the WCS. The **SCI** header should contain (for j 1 and 2):

- CTYPEj
- CRPIXj
- CRVALj
- CDj_j


Example
::

 CTYPE1  = 'LINEAR  '           / R.A. in tangent plane projection               
 CTYPE2  = 'LINEAR  '           / DEC. in tangent plane projection               
 CRPIX1  =                1554. / Ref pix of axis 1                              
 CRPIX2  =                   1. / Ref pix of axis 2                              
 CRVAL1  =                6700. / RA at Ref pix in decimal degrees               
 CRVAL2  =                   1. / DEC at Ref pix in decimal degrees              
 CD1_1   =     -1.3743457180537 / WCS matrix element 1 1                         
 CD2_2   =                   1. / WCS matrix element 2 2                   

