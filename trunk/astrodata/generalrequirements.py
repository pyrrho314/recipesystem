from Requirements import Requirement

class FilenameReq(Requirement):
    file_re = None
    
    def __init__(self, re):
        self.file_re = re
    
    def satisfied_by(self, hdulist):
        import re
        fmatch = re.match(self.file_re, hdulist.filename)
        print "GR12:", hdulist.filename, fmatch
        if fmatch:
            return True
        else:
            return False
        
FILENAME = FilenameReq