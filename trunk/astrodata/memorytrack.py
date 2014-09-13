MEMTRACK = False
import os
import json

def memtrack_write(item):
    tpid = os.getpid()
    fname = "memtrack.%d" % os.getpid()
    
    if not os.path.exists(fname):
        print "storing memory info in %s" % fname
        if os.path.exists("memtrack.latest"):
            os.remove("memtrack.latest")
        os.symlink(fname, "memtrack.latest")
        
    mfile = open(fname, "a+");
    
    line = json.dumps(item)+"\n"
    mfile.write(line)
    mfile.close()    

def memtrack(primname = "unknown_point", msg = "", context = None):
    try:
      if MEMTRACK:
        import psutil
        import os
        import json
        import time
        import gc
        
        #gc.collect()
        tpid = os.getpid()
        
        proc = psutil.Process(tpid)
        memi = proc.get_ext_memory_info()
        item = {    "msg": "%s [%s]" %(primname, msg),
                    "primname": primname,
                    "rss": memi.rss,
                    "vms": memi.vms,
                    #"pfaults": memi.pfaults,
                    #"pageins": memi.pageins,
                    "time": time.time()
                }
      memtrack_write(item)
    except:
      pass
    
def memheader(settings):
    if MEMTRACK:
        pass
