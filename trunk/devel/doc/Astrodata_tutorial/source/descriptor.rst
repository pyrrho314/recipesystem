.. _Make-your-own-descriptor:

Creating a new Instrument descriptor for Astrodata use
======================================================

 To add a new descriptor to AstroData here are the steps:

 ::


     cd $astrodata             # Go to were the astrodata directory is installed in your
                               # machine. If it read-only, then make your own copy.

     cd ../ADCONFIG_Gemini/descriptors    # Astrodata directory is at the same level as
                                          # ADCONFIG_Gemini

     mkdir f2                  # Create new directory for a new intrument 

     cp StandardDescriptorKeyDict.py f2/StandardF2KeyDict.py

     cp nici/calculatorIndex.NICI.py f2/calculatorIndex.F2.py
     
     cp NIRI_RAWDescriptor.py f2/F2_RAWDescriptor.py


     cd f2                                 # cd to f2 and edit the files accordingly

  Edit the new files in f2
  ------------------------
  
  1) Edit StandardF2KeyDict.py   
       This file contains the standard mapping between Astrodata names and F2 
       keyword names. This is a Python dictionary format 'key:value', where 'key'
       is the Astrodata name (format: 'key_f2_<name>') and 'value' is the F2 header keyword name.

  2) Edit calculatorIndex.F2.py
       This is dictionary entry with one pair:

 ::

       calculatorIndex = {"F2_IMAGE":"F2_RAWDescriptor.F2_RAWDescriptorCalc()"}


  3) Edit F2_RAWDescriptor.py

     We are using the ICD document that describe the keyword mapping for
     all Gemini instrument. From the instruction  we made the necessary changes to the
     functions. 

  4) Now add **types**

     cd ../../types

     We should be in a directory: ADCONFIG_Gemini/classifications/types

     mkdir F2

     # Copy the NICI files as an example on how to build new types.

     cp ../nici/* F2


.. _acronyms:

Acronyms
========

- **Unit**  Refer to both header and data portion of any
  extension -including primary unit, of a FITS file 

- **PHDU**  The Primary (Extension 0) Header ``PHU`` and Data ``PDU`` Unit. 
   
- **HDU** Header Data Unit. FITS file reference to header and data 
   portions of a ``Unit`` 

- **HDUList** Pyfits list of FITS descriptors. Each descriptor in the list
  refers to a FITS Unit; e.g. ``phdu = hdulist[0]``

