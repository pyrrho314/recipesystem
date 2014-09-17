# (c) Novem LLC
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.import generalclassification

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
            
        return mmatch
        
HASMEMBER = MemberReq

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
            if isinstance(targ, list):
                mmatch = self.val_check in targ
            else:
                mmatch = (self.val_check == targ)
        except:
            mmatch = False
        
        return mmatch
        
MEMBERCONTAINS = MemberContains        
       
