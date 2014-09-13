.. _prog-guide:

*****************************
AstroData Programmer's Guide
*****************************


1. Opening a FITS file.
  
   ad = AstroData(file, mode)

   - file: FITS file pathname
   - mode: 'readonly', 'update', 'append' or 'new'. In 'new' mode, if the 
     file already exists, it will be silently clobbered.

   The 'ad' instance has all the functions and metadata information pertaining
   to the opened FITS file. You can see this list by typing 'dir(ad)'; any particular
   elements' help file can be obtained by typing for example 'help ad.getHeader'

2. Header information in Astrodata

The following methods are available to retrieve information from a MEF file.
Before access them you need to open a MEF file with Astrodata, please see
the 'info' example here.

=================    ================================================
`info`_              List information about the FITS units:
`getHDU`_            Retrieves Header and Data Unit from an extension
`getHeaders`_        Get headers from all units
`getHeader`_         Get one header
`getHeaderValue`_    Returns the value from the given extension's header
`getPHUHeader`_      Returns primary header
`rePHUKeys`_         Returns all keys that matches a Python regular expression
=================    ================================================

**Examples**


.. _info:

info()
  List information about the FITS units.

 ::

        ad = AstroData('nS20100102S0314.fits')
        print ad.info()
        #
        # The output is:
        Filename: nS20100102S0314.fits
        No.    Name         Type      Cards   Dimensions   Format
        0    PRIMARY     PrimaryHDU     169  ()            int16
        1    SCI         ImageHDU        36  (1024, 1024)  float32
        2    SCI         ImageHDU        35  (1024, 1024)  float32

.. _getHDU:

getHDU (extid)
  This function returns the HDU identified by the ``extid`` argument. This
  argument can be an integer or ``(extname, extver)`` tuple.

 ::

        hdu = ad.getHDU(1)     # Get 1st extension hdu pointer.
        header = hdu.header
        data = hdu.data
        phdu = ad.getHDU(0)    # Get PHDU
          

.. _getHeaders:

getHeaders()
  Function returns header member(s) for all extension (except PHU).

 ::

        hdl = ad.getHeaders()
        print hdl[2].items()          # Print all the extension 2 header keywords

.. _getHeader: 

getHeader(extn)
  Function returns header member for SINGLE EXTENSION MEFs (which are those that
  have only one extension plus PHU).

 ::

        hd = ad.getheader(1)
        print hd.get('cd1_1')

.. _getPHUHeader:

getPHUHeader ()
  This function returns PHU header.

 ::

        ad = AstroData('S20100102S0060.fits')
        phu = ad.getPHUHeader()
        print phu['ut']         # Get the UT value from  the PHU.
                                # It is case insensitive
       
.. _getHeaderValue:

getHeaderValue (extn, keyword)
  This function returns the value from the given extension's header.
  extn: identifies which extension int or (EXTNAME, EXTVER) tuple
  keyword: name of header entry to retrieve

 ::

        ad = AstroData('S20100102S0060.fits')
        print ad.getHeaderValue (2, 'cd1_1')   # using the tuple notation
        print ad.getHeaderValue (('SCI',2), 'cd1_1')
       
.. _phuHeader:

phuHeader(keywd)   (alias is phuValue)
   This function returns a header from the primary header unit
   (extension 0 in a MEF).

 ::

        print ad.phuHeader('UT')          # print UT from the PHU
        print ad.phuValue('UT')           # Same

.. _rePHUKeys:

rePHUKeys(re)
   rePHUKeys returns all keys in this dataset's PHU which match the given
   Python regular expression.

 ::

        print ad.rePHUKeys('..CO')      # 2 characters from the beginning
            ['DECOFFSE', 'FOCOFFEN', 'LPCOUNT', 'LPCOADDS']

        print ad.rePHUKeys('\w*ER')      # Any number of character from the beginning
            ['OBSERVER', 'OBSERVAT', 'FILTER_R', 'PRSERVO', 'FILTER_B', 'TTSERVO']

.. _meta-data-info:

Meta data information in Astrodata
===================================

    
 These are the descriptors available for all the Gemini instruments. They present
 a uniform way of calling, but the value returned varies according to the instrument
 you are querying. Some descriptors accept arguments.

 ::

    Example:  ad = AstroData('S20100102S0060.fits')   # Open a FITS file returning AD object
              print ad.airmass()                      # Get the airmass value

==============  =======================================
    airmass     Airmass
    az          Azimuth
    camera      Instrument
    crpa        Cassegrain rotator value (off, on)
    cwave       ND 
    datalab     Data label
    datasec     Data section
    dec         Declination
    detsec      Detecttor section 
    disperser   ND
    el          Elevation
    exptime     Exposure time
    filterid    Filter identification 
    filtername  Filter name
    fpmask      Focal plane mask
    gain        Gain
    progid      Program ID
    instrument  Instrument used
    mdfrow      ND
    nonlinear   ND
    nsciext     Number of Science extensions
    object      Object name
    obsid       Observatin ID
    obsclass    Observation class
    observer    Observer's name
    obsmode     Observation mode
    obstype     Observation type
    obsepoch    Observation epoch
    pixscale    Pixel scale
    progid      Program id
    pupilmask   Pupil mask
    ra          Right ascension 
    rawiq       Raw Image Quality
    rawcc       Raw Cloud Cover
    rawwv       Raw Water Vapour/Transparency
    rawbg       Raw Background
    rawpireq    PI Requirements Met
    rawgemqa    Gemini Quality Assessment
    rdnoise     read out noise
    satlevel    saturatin level
    ssa         SSA name
    telescope   Telescope name
    utdate      UT date
    uttime      UT time
    wdelta      ND
    wrefpix      ND
    xccdbin      ND
    yccdbin      ND
==============  =======================================

.. _AD-utilities:

AstroData Utilities
===================

 These are the utility calls available from the Astrodata object. In order to access these we
 need to have:

 ::


   # Open a GEMINI observation file.
   ad = AstroData('/tmp/existing_gemini_file.fits')

   Example:  ad.discoverTypes()      # displays all the Types associated with this exposure

======================  ==============================================
`discoverTypes`_          List of processing status and typology 
`filename`_               Filename used as Astrodata argument
`getData`_                Get the pixel array associated with the dataset
`setData`_                Set the data portion of a single extension 
`getStatus`_              processing status of the associated data set
`append`_                 Appends more data and header to the associated set
`getTypes`_               Get the associated dataset typology 
`hdulist`_                object returned by Pyfits when opening the FITS file
`checkType`_              Check for a given type
`hdurefcount`_            number of times *getHDUList* has been accessed 
`header`_                 header portion of the associated hdu
`data`_                   data portion of the associated hdu
`close`_                  Closes the AstroData object
`types`_                  List the astrodata types
`typesStatus`_            processing status of the opened file
`countExts`_              number of extensions match a given *EXTNAME* 
`isType`_                 See if the given type belongs to the typology list
`mode`_                   Astrodata opening mode
`write`_                  Write the hdulis to a new FITS file
======================  ==============================================

EXAMPLES
++++++++

.. _discoverTypes:

discoverTypes
   This function provides a list of classifications of both processing status
   and typology which apply to the data encapsulated by this instance, 
   identified by their string names.

 ::

   example: ad.discoverTypes()
                 ['GEMINI', 'NICI_IMAGE', 'NICI', 'GEMINI_SOUTH', 'PREPARED']

.. _filename:

filename
    Filename used as Astrodata argument.  
    Ex:  print ad.filename    # No ()

.. _getData:

getData
    Get the pixel array associated with the dataset. This is for a single extension MEF.
    See *getHDUlist* for a more general example.

.. _getStatus:

getStatus
    Gives the processing status of the associated data set.

.. _append:

append
    Appends more data and header to the associated set
 
 ::

   Example:  hdu = ad.getHDU(1)      # Get data and header from *ad*
                  adf = AstroData('S20100102S0148.fits')
                  adf.append(data=hdu.data, header=hdu.header)

.. _getTypes:

getTypes
    Get the associated dataset typology. Similar to *discoverTypes*

.. _setData:

setData
    Set the data portion of a single extension
        
.. _hdulist:

hdulist
    Is the object returned by Pyfits when opening the FITS file. 

 ::

       # Get the HDU pair for the 1st extension
       hdu = ad.hdulist[1]

.. _checkType:

checkType (type)
    Check for a given type. Return True or False. The possible types are given
    by *getTypes*

 ::

       example:  ad.checkType('NICI')

.. _hdurefcount:

hdurefcount
    Gives the number of times *getHDUList* has been accessed. *relhdul()* should
    be called after each *getHDUList* call.

.. _header:

header
    Get the header portion of the associated hdu.

.. _data:

data
    Get the data portion of the associated hdu.

 ::

       Example: header = ad[2].header       # Get header from extension 2
                data = ad[2].data           # Get data from extension 2

.. _close:

close()
    Closes the AstroData object

.. _types:

types
    List the astrodata types. 

 ::
     
       Example:  ad.types   # No ()
                 ['NICI_FLAT', 'NICI_DARK', 'GEMINI', 'NICI_DARK_OLD', 
                  'NICI_IMAGE', 'NICI', 'GEMINI_SOUTH', 'UNPREPARED']

.. _typesStatus:

typesStatus
    List the processing status of the opened file.

.. _countExts:

countExts (extname)
    Gives the number of extensions in the MEF files that match a given *EXTNAME* keyword value.

 :: 

       Example:  ad.countExts('SCI')    
          Note:  The *EXTNAME* value is case sensitive

.. _isType :

isType (typename)
    Boolean function to determine if the given type matches the list of types.

 ::
  
       Example:   if ad.isType('NICI_FLAT'):
                      print 'Is Nici Flat frame'
                  else:
                      print 'NICI_FLAT is not a type of',ad.filename

.. _mode:

mode
    Tells the mode used by Astrodata when opening the file.

.. _write:

write (filename)
    Write to a new FITS file. The HDUlist has been already created.

 ::

    EXAMPLE:

    ad = AstroData ('/tmp/existing_file.fits')

    # get the hdulist
    hdl = ad.gethdul()

    # Generate some new data
    data = numpy.arange(100*100)

    # Get the 1st extension header
    header = ad.getHeader(1)
    
    # append a keyword
    header.update('newkk',-1.23,comment='new keyw')


    # Replace Header and data for the 1st extension
    hdl[1].header = header      
    hdl[1].data = data

    # Write out to a new file (Error if file exists)
    # NOTE: It will write the same number of extensions as the original 'ad'
    ad.write('/tmp/new.fits')

    # Or equivalently since *hdl* is of type Pyfits we can use the following
    # to overwrite if needed.
    hdl.writeto("/tmp/new.fits",clobber=True)

    # IF you want to write only one Unit, then you will need to 
    # delete the other members of the 'hdulist'

    del hdl[1:3]       # Keep only the PHU
    hdl[0].data = data

    ad.write('/tmp/PHU.fits')      # Notice that we still use the 'ad' descriptor.


