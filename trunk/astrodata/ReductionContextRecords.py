from datetime import datetime

from astrodata.AstroData import AstroData
from astrodata.generaldata import GeneralData
import os

#------------------------------------------------------------------------------ 
class ReductionContextRecord( object ):
    '''
    The parent record. Contains all members global to all records. (i.e. timestamp)
    '''
    timestamp = None
    
    def __init__( self, timestamp ):
        if timestamp == None:
            timestamp = datetime.now()
        self.timestamp = timestamp
    
  
    
class CalibrationRecord( ReductionContextRecord ):
    '''
    Record for storing all relevant members related to calibration data.
    This is used specifically by the ReductionContext in its calibrations
    member.
    '''
    sciFilename = None
    caltype = None
    filename = None
    source = "all"
    
    def __init__(self, sci_filename, filename, caltype, timestamp = None, source="all"):
        super( CalibrationRecord, self ).__init__( timestamp )
        self.sciFilename = sci_filename
        self.filename = filename
        self.caltype = caltype
        self.source = source
    
    def as_dict(self):
        retd = {"caltype": self.caltype,
                "filename": self.filename,
                "source": self.source
               }
        return retd
        
    def __str__(self):
        rets = """
    sci_filename = %s
    caltype     = %s
    filename    = %s
    timestamp   = %s""" % (  os.path.basename(self.sciFilename), 
                                self.caltype, 
                                os.path.basename(self.filename), 
                                self.timestamp)
        return rets
    
    
class StackableRecord( ReductionContextRecord ):
    '''
    Contains the local cache information for a particular set of stackable data.
    Used in the ReductionContext records stackeep.
    '''
    stkid = None
    filelist = []
    
    def __init__( self, stkid, filelist, timestamp=None ):
        super( StackableRecord, self ).__init__( timestamp )
        self.stkid = stkid
        self.filelist = filelist
    
    def __str__(self):
        rets = """
    stkid     = %s
    filelist  = %s
    timestamp = %s \n""" % ( str(self.stkid), str(self.filelist), self.timestamp )
        return rets  

class FringeRecord( ReductionContextRecord ):
    '''
    Contains the cache information for a set of fringeable data.
    '''
    fringeid = None
    listid = None
    filelist = []
    
    def __init__(self, fringeid, listid, filelist, timestamp=None):
        super( FringeRecord, self ).__init__( timestamp )
        self.fringeid = fringeid
        self.listid = listid
        self.filelist = filelist
        
        
    def __str__(self):
        rets = '''
fringeID   = %s
list_id     = %s
filelist   = %s
timestamp = %s
''' % ( str(self.fringeid), str(self.listid), str(self.filelist), self.timestamp )

    def add(self, image):
        if type(image) == list:
            self.filelist.extend( image )
        else:
            self.filelist.append( image )

##@@FIXME: Because of the nature of how Output -> Input, the name of this record may need to change
## at some point.
class DataObjectRecord( ReductionContextRecord):
    _data_obj = None
    
#    def __init__(self, *args, **args):
#        super(DataObjectRecord, self).__init__(*args, **args) 
        
    def get_data_obj(self):
        return self._data_obj
    def set_data_obj(self, val):
        self._data_obj = val
    def del_data_obj(self):
        self._data_obj = None
        
    data_obj = property(get_data_obj, set_data_obj, del_data_obj)
    
    
    
class AstroDataRecord( DataObjectRecord ):
    '''
    Contains any metadata related to output/input within the ReductionContext.
    This is used specifically in the ReductionContext records inputs and outputs.
    '''
    display_id = None
    filename = None
    ad = None
    parent = None
    
    # making self.ad compatible with DataObjectRecort.data_obj
    def get_data_obj(self):
        return self.ad
    def set_data_obj(self, val):
        self.ad = val
    def del_data_obj(self):
        self.ad = None
    data_obj = property(get_data_obj, set_data_obj, del_data_obj)
    
    def __init__(self, filename, display_id=None, timestamp=None, parent=None, load = True):
        super( AstroDataRecord, self ).__init__( timestamp )
        #print "RCR110:", type(filename), isinstance(filename, AstroData)
        if isinstance(filename, AstroData):
            self.filename = filename.filename
            self.ad = filename
            self.parent = parent #filename.filename
        elif type( filename ) == str:
            self.filename = filename
            if load == True:
                self.ad = AstroData( filename )
            else:
                self.ad = None
            self.parent = parent
        elif type( filename ) == AstroDataRecord:
            adr = filename
            self.display_id = adr.display_id
            self.filename = adr.filename
            self.ad = adr.ad
            self.parent = adr.parent        
            return                      
        else:
            raise "BAD ARGUMENT"
        ##@@TODO: display_id may be obsolete
        self.display_id = display_id
    def load(self):
        self.ad = AstroData(self.filename)
    def is_loaded(self):
        return not (self.ad == None)    
    def __str__(self):
        rets = """
    ad.filename  = %s
    display_id   = %s
    filename     = %s
    timestamp    = %s
    parent       = %s
    astrodata    = %s \n""" % ( self.ad.filename, 
                                str(self.display_id),
                                str(self.filename), 
                                self.timestamp, 
                                self.parent, 
                                str(self.ad) )
        return rets
    
    def data_obj(self):
        return ad
        
class GeneralDataRecord( DataObjectRecord ):
    '''
    Contains any metadata related to output/input within the ReductionContext.
    This is used specifically in the ReductionContext records inputs and outputs.
    '''
    display_id = None
    filename = None
    parent = None
    gd = None
    
    def get_data_obj(self):
        return self.gd
    def set_data_obj(self, val):
        self.gd = val
    def del_data_obj(self):
        self.gd = None
    data_obj = property(get_data_obj, set_data_obj, del_data_obj)
    
    def __init__(self, filename, display_id=None, timestamp=None, parent=None, load = True):
        super( GeneralDataRecord, self ).__init__( timestamp )
        #print "RCR110:", type(filename), isinstance(filename, AstroData)
        if isinstance(filename, GeneralData):
            self.filename = filename.filename
            self.gd = filename
            self.parent = parent #filename.filename
        elif type( filename ) == str:
            self.filename = filename
            if load == True:
                self.gd = GeneralData.create_data_object( filename )
            else:
                self.gd = None
            self.parent = parent
        elif type( filename ) == GeneralDataRecord:
            gdr = filename
            self.display_id = gdr.display_id
            self.filename   = gdr.filename
            self.gd         = gdr.ad
            self.parent     = gdr.parent        
            return                      
        else:
            raise "BAD ARGUMENT"
        ##@@TODO: display_id may be obsolete
        self.display_id = display_id
            
    def load(self):
        self.gd = GeneralData.create_data_object(self.filename)
        print "RCR221: loading %s %s" %(self.filename, self.gd)
        
    def is_loaded(self):
        return not (self.gd == None)    
    def __str__(self):
        rets = """
    gd.filename  = %s
    display_id   = %s
    filename     = %s
    timestamp    = %s
    parent       = %s
    type         = %s \n""" % ( self.gd.filename, 
                                str(self.display_id),
                                str(self.filename), 
                                self.timestamp, 
                                self.parent, 
                                str(self.gd) )
        return rets  
