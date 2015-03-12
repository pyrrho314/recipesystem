# Copyright (C) 2014 Novem LLC, created Craig Allen 2014
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.import generalclassification

import os
from astrodata.Errors import Error
from astrodata.AstroDataType import get_classification_library
from astrodata.adutils import termcolor as tc

cl = get_classification_library() 

class GD_OperationNotPermitted(Error):
    """ GeneralData related exceptions, raised if an operation is not permitted,
        like setting an output_filename property (it's constructed from the
        output_directory and the basename).
    """
    message  = "Operation Not Permitted"

# needs to move to configuration space
# make this a member of GeneralData (?)
_data_object_classes = {
    "tif": ("hbdata","HBRasterData"),
    "fits": ("astrodata.AstroData", "AstroData"),
    "json": ("jsondata", "JSONData"),
    "txt": ("jsondata", "TxtData"),
    "xls": ("pandasdata", "PandasData"),
    "csv": ("pandasdata", "PandasData"),
    "h5": ("pandasdata", "PandasData"),
    "shp": ("hbshapedata", "HBShapeData"),
    "setref": ("jsondata", "ReferenceOnlyData"),
    "dmo": ("cubedata", "CubeDemoData"),
    "dmobson":("cubedata", "CubeDemoData")
    }
_data_object_precedence = [ "tif",
                            "shp",
                            "fits", 
                            "json",
                            "txt",
                            "xls",
                            "csv",
                            "h5",
                            "setref",
                            "dmobson",
                            "dmo"
                          ]
_gen_classification_library = cl

class HDUnit(object):
    # abstracts the idea of a unit of data and it's metadata
    pass
    
class GeneralData(object):
    # filename is a property
    _filename = None
    # output_directory is a property
    _output_directory = None
    _suggested_data_classes = {
        "SHAPEDATA":("hbshapedata", "HBShapeData"),
        "TABLE": ("jsondata", "PandasData"),
        "TXT" : ("jsondata", "TxtData")
        }
    _saved = True
    _changed = False
    
    @classmethod
    def create_data_object(cls, initarg, hint = None):
        #print "gd26: initarg", initarg, hint
        retv = None
        if True: 
            found_instor = None
            # @@ TODO: improve the generalization, externalize to configspace
            if hint:
                instor = hint
                module,tclass = hint
            else:
                if isinstance(initarg, basestring):
                    for suffix in _data_object_precedence:
                        completesuffix = ".%s"%suffix
                        # print "gd71:", completesuffix
                        instor = _data_object_classes[suffix]
                        lowerarg = initarg.lower()
                        # print "gd74:", lowerarg
                        if (lowerarg.endswith(completesuffix)):
                            found_instor = instor
                            break;
                    # if no file matches by extension, setref will be used
            if isinstance(initarg, basestring):
                initarg = os.path.abspath(initarg)
                # print "gd47:", initarg
            if found_instor:
                module, tclass = found_instor
                # print "gd29:",module, tclass
                exec("from %s import %s" % instor)
                retv = eval("%s(initarg)" % tclass)
            else:
                bname,ext = os.path.splitext(initarg)
                raise IOError("generaldata.GeneralData.create_data_object\n" +
                              "DON'T KNOW HOW TO LOAD DATASET\n" +
                              "     %s\n" % initarg +
                              "     unknown suffix: %s" % ext )
        return retv
    
    ############
    #properties
    
    # filename
    debug_filename = False
    def _get_filename(self):
        if self.debug_filename:
            print "gd63: get filename = %s" %self._filename
        return self._filename
        
    def _set_filename(self, name):
        if self.debug_filename:
            if self._filename and name != self._filename:
                print "&*()*&()"*10
                print "gd65: %s -> %s" % (name, self._filename)
                import traceback
                traceback.print_stack()
                print "&*()*&()"*10
        self._filename = name
        if self.debug_filename:
            print "gd73: set filename = %s" %self._filename
        self.prop_put("_data.filename", self._filename)
        return
        
    def _del_filename(self, name):
        if self.debug_filename:
            print "gd76: del filename = %s" %self._filename
        self._set_filename(None)
    filename = property(_get_filename, _set_filename, _del_filename)
    # end filename
    # output directory
    def _get_output_directory(self):
        if self._output_directory == None:
            self._output_directory = os.getcwd()
            
        return self._output_directory
        
    def _set_output_directory(self, val):
        self._output_directory = val
        return
        
    def _del_output_directory(self):
        self._output_directory = None
        return
    output_directory = property(_get_output_directory, 
                                _set_output_directory, 
                                _del_output_directory)
                                
    def _get_output_filename(self):
        return os.join.path(self.output_directory, self.basename)
    def _set_output_filename(self):
        raise GD_OperationNotPermitted("cannot set 'output_filename', it is constructed from 'output_dir' and 'basename'")
    output_filename = property(_get_output_filename, _set_output_filename)
    
    # end output directory
    
    # properties
    #############
        
    # is this still needed?
    def _get_ad(self):
        return self
    
    ad = property(_get_ad)
                  
    def __getitem__(self, arg):
        return self
    
    def __init__(self, initarg, force_load = None):
        # initarg will likely be a filename
        pass
    
    def _basename(self):
        return os.path.basename(self.filename)
    basename = property(_basename)
    
    def _dirname(self):
        return os.path.dirname(self.filename)
    dirname = property(_dirname)
           
    # virtuals... for children
    def do_write(self):
        # child class does this
        pass
        
    # properties
    def prop_exists(self, name):
        return False
    def prop_put(self, name, *args, **argv):
        return None
    def prop_get(self, name):
        return None
    def prop_add(self, name, item):
        return None
    
    def get_classification_library(self):
        global _gen_classification_library
        if _gen_classification_library == None:
            _gen_classification_library = globalClassificationLibrary
            
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
    
    def is_type(self, settype):
        types = self.get_types()
        return settype in types
           
    def add_suffix(self, suffix = None):
        fname = self.with_suffix(suffix = suffix)
        self.filename = fname
        self.is_changed()
        return fname
    
    def with_suffix(self, suffix = None):
        if not suffix:
            return self.filename
        base = os.path.basename(self.filename)
        dirname = os.path.dirname(self.filename)
        basecore, baseext = os.path.splitext(base)
        
        basecore += "_%s" % suffix
        
        fname = os.path.join(dirname, basecore)+baseext
        return fname    
    
    def add_prefix(self, prefix = None):
        fname = self.with_prefix(prefix = prefix)
        self.filename = fname
        self.is_changed()
        return fname
        
    def with_prefix(self, prefix = None):
        
        if not prefix:
            return self.filename
        base = self.basename
        dirname = self.dirname
        fname = os.path.join(dirname, "%s-%s" % (prefix, base))
        # print "with_prefix:", prefix, fname, self.filename
        return fname
        
    def allow_extant_write(self):
        return False
        
    def nativeStorage(self):
        return False
    
    def recommend_data_object(self):
        import gdpgutil as gd
        hints = self.get("type_hints")
        print "gd190: hints=", hints
        if hints:
            modandclass = gd.pick_config(hints, self._suggested_data_classes, "leaves")
            return modandclass
        else: 
            modandclass = gd.pick_config(self, self._suggested_data_classes, "leaves")
            return modandclass
            
    def write(self, suffix=None, ** args):
        # make our filename relative to the output dir
        
        oldname = self.filename
        outname = os.path.join(self.output_directory, self.basename)
        self.filename = outname
        newname = self.filename
        if suffix:
            newname = self.add_suffix(suffix)
        if os.path.exists(newname):
            print "(gd225) %s already exists" % newname
            # check exists_policy
            if self.allow_extant_write():
                # e.g.: SetrefData simply moves the extant out of the way using a ";N" postfix
                print "        but %s allows extant write" % tc.colored(repr(type(self)), attrs=["bold"])
            else:
                raise GD_OperationNotAllowed("General Data does not allow overwriting data by default.")
                
        self.do_write(newname)
        
        self._saved = True
        return True
        
    # file state BEGIN
    def needs_write(self):
        self._saved = False
    def is_changed(self):
        self._changed = True
        self.needs_write()
    # file state END        
    def write_nonexistant(self):
        if not os.path.exists(self.filename):
            self.write()
            return True
        else:
            return False
