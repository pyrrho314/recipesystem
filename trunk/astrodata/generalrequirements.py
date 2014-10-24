# Copyright (C) 2014 Novem LLC, created Craig Allen 2014
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.import generalclassification

DEBUG=True

from Requirements import Requirement

class FilenameReq(Requirement):
    file_re = None
    
    def __init__(self, re):
        self.file_re = re
    
    def satisfied_by(self, dataset):
        import re
        fmatch = re.match(self.file_re, dataset.filename)
        #print "GR12:", dataset.filename, fmatch
        if fmatch:
            return True
        else:
            return False
        
FILENAME = FilenameReq

class MemberReq(Requirement):
    member_name = None
    def __init__(self, memstr):
        self.member_name = memstr
        
    def satisfied_by(self, dataset):
        mmatch = False
        try: # @@EVAL
            eval("dataset.%s" % self.member_name)
            mmatch = True
        except:
            mmatch = False
        
        if DEBUG:
            print "HASMEMBER %s is %s" %(self.member_name, mmatch)   
        return mmatch
        
HASMEMBER = MemberReq

# class MemberIs
# MEMBERIS

class MemberContains(Requirement):
    member_name = None
    val_check = True
    
    def __init__(self, memname, memval):
        self.member_name = memname
        self.val_check = memval
    
    def satisfied_by(self, dataset):
        mmatch = False
        
        try:
            targ = eval("dataset.%s" % (self.member_name))
            if hasattr(targ, "__contains__"): #isinstance(targ, list):
                #print "supports IN", self.val_check, targ
                mmatch = self.val_check in targ
            else:
                #print "no support IN"
                mmatch = (self.val_check == targ)
            # print "matched=",mmatch
        except:
            mmatch = False
            
        return mmatch
MEMBERCONTAINS = MemberContains

# Properties

class PropertyReq(Requirement):
    prop_name = None
    def __init__(self, memstr):
        self.prop_name = memstr
        
    def satisfied_by(self, dataset):
        mmatch = False
        mmatch = self.prop_exists(self.prop_name)
        return mmatch
        
HASPROP = PropertyReq

class PropertyIsReq(Requirement):
    prop_name = None
    prop_val = None
    def __init__(self, propkey, propval):
        self.prop_key = propkey
        self.prop_val = propval
    def satisfied_by(self, dataset):
        val = dataset.prop_get(self.prop_key)
        # print "gr90:", self.prop_key, self.prop_val, val
        if val == self.prop_val:
            return True
        else:
            return False
        
PROPERTY = PropertyIsReq

class PropContains(Requirement):
    prop_key = None
    prop_val = True
    
    def __init__(self, propkey, propval):
        self.prop_key = propkey
        self.prop_val = propval
    
    def satisfied_by(self, dataset):
        container = self.prop_get(self.prop_key)
        mmatch = False
        try:
            mmatch = self.prop_val in container
        except:
            mmatch = False
        
        return mmatch
        
PROPCONTAINS = PropContains        
       
