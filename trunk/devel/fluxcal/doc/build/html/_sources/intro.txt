Introduction
____________

  The FLUXCAL program calculates the Zero Point correction from a set
  of reference stars found in the field of the input image. The results are
  appended to 'ZPcorr.log' file in the current working directory.

  Fluxcal appends to the input FITS file a BINTABLE extension with EXTNAME
  OBJCAT with:

  - source ID number (simple serial number unique within this table)
  - Centroid X pixel co-ordinate
  - Centroid Y pixel co-ordinate
  - (Why do we need this? nz) Source flux (in ADUs / whatever units the data are in at this point)
  -  error on flux measurement
  -  If the cataloger measures flux using object and sky appertures, then these values should probably be listed too
  -  RA (assume WCS in image is correct)
  -  Dec (ditto)
  -  RefID (Reference ID of source in the field)
  -  smag (magnitude for the reference source)

  And another BINTABLE REFCAT with:

  - Reference ID from the reference catalog
  - RA (from the Reference Catalog)
  - Dec 
  - X pixel co-ordinate  from the image WCS
  - Y pixel co-ordinate (ditto)
  - smag (magnitude for the filter used)

  Algorithm

  - Read FITS table extension(S) with extname 'OBJCAT'
    containing columns 'refid' and 'smag' with at least
    one data point. These have the reference standard object name
    and magnitudes for the filter used in the exposure.
  - Each x,y position  from OBJCAT containing a non-empty standard
    in 'refid' are used to create a text file 'ZPcoo' in the working
    directory.
  - Call iraf.imexam to get FWHM
  - Call iraf.apphot.phot to get Flux, and iMag (instrumental magnitude)
  - Zero point correction: ( 'k' is  the extinction coefficient)

        ZPCorr = Reference_Magnitude - [iMag - k*(airmass-1)]
  - Image Quality (IQ) is calculated as:

        Seeing = FWHM*pixel_scale

        IQ = [mean(Seeing) / (airmass**0.6)]
  - Append keywords ZP_CORR, ZP_ERROR, IMAGE_IQ

  ERROR Propagation:

  - ZP_ERROR: At this time only the GS 'smith.cat' has error values
    for each magnitude. In the meantime we only calculate the standard
    deviation of the ZPcorr array and assign it to ZP_ERROR.

  OUTPUT FILE:

  - ZPcorr.log is the filename to append the new ZP values:
    "ZPCorr  ZPerror IQ     Sample Filter  Date_time".

