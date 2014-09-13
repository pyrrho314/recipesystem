import json
import pprint
from astrodata.adutils import ksutil
import generaldata

class JSONData(generaldata.GeneralData):
    json = None
    filename = None
    def __init__(self, initarg, defer_load = False):
        super(JSONData, self) .__init__( initarg)
        if not defer_load:
            self.load(initarg)
            
    def do_write(self, fname, rename = False):
        if rename:
            self.filename = fname
        
        tfile = open(self.filename, "w")
        tfile.write(json.dumps(self.json))
        tfile.close()
        
    
    def load(self, initarg):
        self.filename = initarg
        jsonfile = open(initarg)
        self.json = json.load(jsonfile)
        jsonfile.close()
        #print "jd14:",self.pprint()
        
    def pretty_string(self):
        retstr = ksutil.dict2pretty(self.filename, self.json)
        return retstr
    
    def close(self):
        pass

    def get_types(self, prune = False):
        types = ["JSON"]
        types.append(str(type(self.json)))
        return types
