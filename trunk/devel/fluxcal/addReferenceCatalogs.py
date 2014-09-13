import os
import imp
from copy import deepcopy

import pyfits as pf
import pywcs
import numpy as np
from astrodata import AstroData
from astrodata.adutils import gemLog
from detectSources import DetectSources

class AddReferenceCatalogs(object):
    """ 
      **Syntax**
           ar = addReferenceCatalogs.AddReferenceCatalogs(
                (ad, catalogName=None, logLevel=6, logfile='')
        
      :param ad:
            AstroData object containing a FITS hdulist with at least one
            image extension with EXTNAME='SCI' in the extension header.
      :type ad: AstroData object
      :param catalogName:
            Filename of the catalog. 'smith.cat' and 'gmosn.cat'
            are supported at this time.
      :param logLevel: 
            Verbose level
      :type logLevel: integer,default is 6.
      :param logfile:
             log file name to replace the default value
      :type logfile: String. Default value is 'detectSources.log'
          
      **Description**
        Find standard stars from the current standard catalogs
        gmosn.cat (GN) and smith.cat (GS) that falls into the
        current field of the image extension from input file.
        If any is found creates a FITS table extension to be 
        append to the file.

         ::

                The table description is
                EXTNAME = 'REFCAT'
                EXTVER  = extn      # Extension number of the image.
                REFID   = str       # Unique id from the standard catalog
                RA      = ra        # ra and dec from the catalog
                DEC     = dec
                X       = x         # x,y pixel coordinate calculated with
                Y       = y         # the image WCS
                MAG     = mag       # standard star magnitude
                TUNIT   = filter_name   # filter letter corresponding to
                                      # FILTER2 value in the GMOS PHU

    """
    def __init__(self, ad, catalogName=None, logLevel=6, logfile=''):

        self.ad = ad

        basename = os.path.basename(ad.filename)

        self.outad = deepcopy(ad)
        self.outad.filename = 'ref_'+basename

        if not logfile:
            logfile = 'selectreferences.log'

        self.log = gemLog.getGeminiLog(logName=logfile, logLevel=logLevel)
        log = self.log
        log.defaultCategory(level='ALL', category='selectReferences')

    
        log.info( "\n ****** SELECT REFERENCES for: %s *******"%ad.filename)

    def getRefs(self):
        """ Run the add reference catalog. Actually adding the 
            Bintable to the input ad object.
 
        """

        from pyfits import Column

        log = self.log

        extname = 'REFCAT'

        outad = self.outad

        # Select catalog and format the output data
        usecols,formats,band,delimiter = self.selStdsCatalog()
        refid,ra,dec,fmag = self.readStds(usecols, formats, delimiter)

        # Loop through the SCI extensions
        for scix in outad['SCI']:

            xtver = scix.extver()

            # x,y are the coordinates of the reference stars within the 
            # input image field.
            g,x,y = self.search4standards(ra, dec, xtver)

            log.info("Found %d standards for field in  %s['SCI',%d]"%\
                (len(g[0]),outad.filename,xtver))

            # g: index array with the index of the standards within the field.
            if len(g[0])>0:
                nlines = len(ra)

                # If extension already exists, just update
                if outad[extname,xtver]:
                    log.info('Table already exists,updating values.')
                    tdata = outad[extname,xtver].data
                    theader = outad[extname, xtver].header
                    tdata.field('refid')[:] = refid[g]
                    tdata.field('ra')[:]    = ra[g]
                    tdata.field('dec')[:]   = dec[g]
                    tdata.field('x')[:]     = x
                    tdata.field('y')[:]     = y
                    tdata.field('refmag')[:]  = fmag[g]
                else:
                    c1 = Column (name='refid', format='22A', array=refid[g])
                    c2 = Column (name='ra',    format='E', array=ra[g])
                    c3 = Column (name='dec',   format='E', array=dec[g])
                    c4 = Column (name='x',     format='E', array=x)
                    c5 = Column (name='y',     format='E', array=y)
                    # band:       1-char:  'u','g','r','i' or 'z'
                    c6 = Column (name='refmag',unit=band,format='E',array=fmag[g])
                    colsdef = pf.ColDefs([c1,c2,c3,c4,c5,c6])

                    tbhdu = pf.new_table(colsdef)         # Creates a BINTABLE

                    # pyfits to AstroData
                    tabad = AstroData(tbhdu)

                    # Add or append keywords EXTNAME, EXTVER
                    tabad.rename_ext(extname, xtver)
                
                    outad.append(tabad)
            else:
                log.warning( 'No standard stars were found for this field.')

        return outad

    def selStdsCatalog(self):
        """
          Select standard catalog based on the value of
          the instrument keyword. 
          The values is GMOS_N, the catalog is 'gmosn.cat'
          else 'smith.cat'

          INPUT:
             self.ad: AD object 

          OUTPUT:
              - usecols
              - formats 
              - band:    1-char:  'u','g','r','i' or 'z'
              - delimiter
        """
        filter = self.ad.filter_name() 

        if self.ad.is_type('GMOS_S'): 
            # Read 'smith.cat' catalog of standards for GS.
            # ID,ra,dec,u,ue,un, g,re,gn, r,re,rn, i, ie, in, z,ze,zn
            # 0  1  2   3 4  5   6 7  8   9 10 11  12 13  14 15 16 17

                                        
            self.catName  = 'smith.cat'
            fn = {'u':3,'g':6,'r':9,'i':12,'z':15}     # column position for filter mag
            fne = {'u':4,'g':7,'r':10,'i':13,'z':16}   # column position for filter mag err
            band = filter[0]                   # Get the first filtername letter
                                               # u,g,r,i,z
            usecols = (0,1,2,fn[band],fne[band.lower()])
            formats = ('S22','f4','f4','f4','f4')
            delimiter = ' '

        elif self.ad.is_type('GMOS_N'):
            # gmosn.cat
            # name,Field,RA,DEC,u_mag,v_mag,g_mag,r_mag,i_mag,z_mag,
            # 0    1     2  3   4     5     6     7      8    9
            #             y_mag,j_mag,h_mag,k_mag,lprime_mag,m_mag
            #             10     11    12   13    

            self.catName = 'gmosn.cat'
            fn = {'u':4,'g':6,'r':7,'i':8,'z':9,'y':10,'j':11,'h':12,'k':13}
            band = filter[0]                   # Get the first filtername letter
                                                # u,g,r,i,z
            usecols = (0,2,3,fn[band.lower()])
            #NOTE fmag can have None
            formats = ('S20','f4','f4','S20')
            delimiter = ','
        else:
            raise RuntimeError, 'GMOS file keyword INSTRUME is not GMOS-N nor GMOS-S.'
            return None
            
        return usecols,formats,band,delimiter

    def readStds(self,usecols,formats,delimiter):
        """ Read selected References Catalog (ascii catalogs)

           INPUT:
             - self.catname:  Pathname to input References catalog.
             - usecols:    Columns to read from catalog.
             - formats:     Formats to read columns. EG: ('S20','f4','f4','S20')
                            For more details please see:
                            python: help numpy.loadtxt
            
           OUTPUT:
             # Arrays containing standars info for objects in the field.

             - refid:     Catalog id
             - ra,dec:    ra,dec
             - fmag:      magnitude corrresponding to 
                         the selected band 
            

        """
 
        log = self.log
        fp, pathname, description = imp.find_module('addReferenceCatalogs')

        catName = os.path.join(os.path.dirname(pathname),self.catName)
        log.info( 'Catalog name: '+catName)
        cc = np.loadtxt(catName, usecols=usecols, delimiter=delimiter, 
                         dtype={'names':('ic','ra','dec','fmag'),
                                'formats':formats})
        refid = cc['ic']
        ra =    cc['ra']
        dec =   cc['dec']
        fmag =  cc['fmag']

        # change fmag to float if necessary
        if hasattr(fmag[0],'isupper'):    # See if it is string
            for i in range(len(fmag)):
                if fmag[i] != 'None':
                   fmag[i] = float(fmag[i])
                else:
                   fmag[i] = -999         # Floating flag for None

        return refid,ra,dec,fmag

    def search4standards(self,ra,dec,ext):
        """
          Check if standards given by ra,dec fall in the
          image field determine by its WCS information.

          INPUT:         
            - ra,dec: ra,dec of the whole input catalog.
                      ra is in hour units.
          OUTPUT:
            - g:      Index array of the positions (ra,dec) that 
                      falls in the field.
        """

        scihdr = self.ad['SCI',ext].header

        wcs = pywcs.WCS(scihdr)

        # Convert world coordinates to pixel coordinates
        # The second argument is "origin" -- in this case we're declaring we
        # have 1-based (Fortran-like) coordinates.

        xy = wcs.wcs_sky2pix(zip(ra*15,dec),1)
        xx = xy[:,0]
        yy = xy[:,1]

        naxis1,naxis2 = scihdr['naxis1'],scihdr['naxis2']

        g = np.where((yy>0) & (yy<naxis2) & (xx>0) & (xx<naxis1))

        return g,xx[g],yy[g]
