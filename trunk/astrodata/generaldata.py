import generalclassification
import os

# needs to move to configuration space
_data_object_classes = {
    ".fits": ("astrodata.AstroData", "AstroData"),
    ".json": ("jsondata", "JSONData"),
    ".csv": ("jsondata", "TxtData"),
    ".txt": ("jsondata", "TxtData")
    }

_gen_classification_library = None

class HDUnit(object):
    # abstracts the idea of a unit of data and it's metadata
    pass
    
class GeneralData(object):
    filename = False
    
    @classmethod
    def create_data_object(cls, initarg):
        retv = None
        if (isinstance(initarg, str)):
            # @@ TODO: improve the generalization, externalize to configspace
            for suffix in _data_object_classes:
                instor = _data_object_classes[suffix]
                if (initarg.endswith(suffix)):
                    module, tclass = instor
                    print "gd29:",module, tclass
                    exec("from %s import %s" % instor)
                    retv = eval("%s(initarg)" % tclass) 
        else:
            from astrodata.AstroData import AstroData
            retv = AstroData(initarg)
            
        return retv
    def _get_ad(self):
        return self
    
    ad = property(_get_ad)
                  
    def __getitem__(self, arg):
        return self
    
    def __init__(self, initarg, defer_load = False):
        # initarg will likely be a filename
        pass
    
    def do_write(self):
        # child class does this
        pass
    
    def get_classification_library(self):
        global _gen_classification_library
        if _gen_classification_library == None:
            _gen_classification_library = generalclassification.GeneralClassificationLibrary()
            
        return _gen_classification_library
    
    def close(self):
        print "GeneralData::close called, default parent method, does nothing"
        pass
    
    def get_types(self, prune = False):
        # don't assume any types, subtypes of data 
        from astrodata.AstroDataType import ClassificationLibrary
        cl = ClassificationLibrary.get_classification_library()
        types = cl.discover_types(self)
        return types
        
    types = property(get_types)
    
    def add_suffix(self, suffix = None):
        if not suffix:
            return False
        base = os.path.basename(self.filename)
        dirname = os.path.dirname(self.filename)
        basecore, baseext = os.path.splitext(base)
        
        basecore += suffix
        
        fname = os.path.join(dirname, basecore)+baseext
        return fname
        
    
    def write(self, clobber = False, suffix=None, rename = True):
        oldname = self.filename
        newname = self.add_suffix( suffix)
        if newname:
            if os.path.exists(newname) and not clobber:
                print "gd83: not writing, %s already exists, and clobber==False" % newname
                return False
            
        self.do_write(newname, rename = rename)
        
        if rename and self.filename == oldname:
            print "gd87: rename == True but name was not changed, still %s, expected %s" % (self.filename, newname)
