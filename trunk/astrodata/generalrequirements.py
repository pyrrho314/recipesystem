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
