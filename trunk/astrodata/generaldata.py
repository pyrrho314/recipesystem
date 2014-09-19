# Copyright (C) 2014 Novem LLC, created Craig Allen 2014
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.import generalclassification

import os
from astrodata.AstroDataType import get_classification_library
cl = get_classification_library() 

# needs to move to configuration space
_data_object_classes = {
    ".fits": ("astrodata.AstroData", "AstroData"),
    ".json": ("jsondata", "JSONData"),
    ".txt": ("jsondata", "TxtData"),
    ".xls": ("jsondata", "PandasData"),
    ".csv": ("jsondata", "PandasData"),
    ".setref": ("jsondata", "SetrefData")
    }

_gen_classification_library = cl

class HDUnit(object):
    # abstracts the idea of a unit of data and it's metadata
    pass
    
class GeneralData(object):
    # filename is a property
    _filename = None
    # output_directory is a property
    _output_directory = None
    
    @classmethod
    def create_data_object(cls, initarg, hint = None):
        #print "gd26: initarg", initarg, hint
        retv = None
        if True: 
            # @@ TODO: improve the generalization, externalize to configspace
            if hint:
                instor = hint
                module,tclass = hint
            else:
                if isinstance(initarg, basestring):
                    for suffix in _data_object_classes:
                        instor = _data_object_classes[suffix]
                        lowerarg = initarg.lower()
                        if (lowerarg.endswith(suffix)):
                            break;
            if isinstance(initarg, basestring):
                initarg = os.path.abspath(initarg)   
                # print "gd47:", initarg     
            module, tclass = instor
            # print "gd29:",module, tclass
            exec("from %s import %s" % instor)
            retv = eval("%s(initarg)" % tclass)
            
        return retv
    
    ############
    #properties
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
        return
        
    def _del_filename(self, name):
        if self.debug_filename:
            print "gd76: del filename = %s" %self._filename
        self._set_filename(None)
    filename = property(_get_filename, _set_filename, _del_filename)
    
    def _get_output_directory(self):
        if self._output_directory == None:
            self.output_directory == os.getcwd()
            raise "catch this"
            
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
        return os.path.dirname(self.dirname)
    dirname = property(_basename)
           
    def do_write(self):
        # child class does this
        pass
    
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
    
    
    
    def add_suffix(self, suffix = None):
        if not suffix:
            return self.filename
        base = os.path.basename(self.filename)
        dirname = os.path.dirname(self.filename)
        basecore, baseext = os.path.splitext(base)
        
        basecore += suffix
        
        fname = os.path.join(dirname, basecore)+baseext
        return fname
        
    def nativeStorage(self):
        return False
        
    def write(self, clobber = False, suffix=None, rename = True):
        # make our filename relative to the output dir
        
        oldname = self.filename
        outname = os.path.join(self.output_directory, self.basename)
        self.filename = outname
        newname = self.add_suffix(suffix)
        if newname:
            if os.path.exists(newname) and not clobber:
                print "gd83: not writing, %s already exists, and clobber==False" % newname
                return False
            
        self.do_write(self.filename, rename = rename)
        
        if rename and self.filename == oldname:
            print "gd87: rename == True but name was not changed, still %s, expected %s" % (self.filename, newname)
