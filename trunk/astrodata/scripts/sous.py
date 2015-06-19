#!/usr/bin/env python
import json
from geventwebsocket import WebSocketServer, WebSocketApplication, Resource
from subprocess import Popen, PIPE
from time import sleep
import os
from astrodata import generaldata
from astrodata.generaldata import GeneralData
from astrodata import Lookups
from base64 import b64encode, b64decode
from pylab import *
import numpy as np
dw_info = Lookups.compose_multi_table(  "*/warehouse_settings", 
                                        "warehouse_elements", 
                                        "shelf_addresses", 
                                        "type_shelf_names",
                                        "type_store_precedence"
                                      )
wpack = Lookups.compose_multi_table("*/warehouse_settings", "warehouse_package")
outpacks = []
for pack in wpack:
    outpacks.append(repr(pack));    
    

# type of quick vieew image to use, use name from mimetype (image/<imext>)    
imext = "png"
IMEXT = imext            
                                      
dw_info["warehouse_package"] = outpacks
class EchoApplication(WebSocketApplication):
    _ra_stdin   = None
    _ra_stdout  = None
    _ra_stderr  = None
    _ra_proc    = None
    
    def on_open(self):
        print "Connection opened"

    def on_message(self, message):
        print "message:",type(message),message
        msg = json.loads(message)
        cmd = msg["cmd"]
        if cmd == "depot_msg":
            print "handle depot_msg"
            subcmd = msg["subcmd"] if "subcmd" in msg else None
            ##### DISPLAY
            if subcmd == "display":
                print "client wants a %s" % IMEXT
                fn = msg["options"]["args"][0]
                imgname = "%s.%s" % (fn,imext)
                numnz=-1
                if not os.path.exists(imgname):
                    stats = os.stat(fn)
                    imsize = stats.st_size 
                    a = GeneralData.create_data_object(fn)
                    
                    if imsize > 100000:
                        progct = { "cmd":"nrm_depot",
                                    "subcmd":"display_status",
                                    "status_msg": "creating quick view on server"
                                  }
                        pmsg = json.dumps(progct)
                        self.ws.send(pmsg)
                    
                    nd = a.get_nd(1)
                    if a.data.RasterCount >=3:
                        cd = np.zeros( (nd.shape[0], nd.shape[1], 3), dtype=np.uint8)
                        cd[:,:,0] = nd[:]
                        for i in range(1,a.data.RasterCount-1):
                            xd = a.get_nd(i+i)
                            cd[:,:,i] = xd[:,:]
                        nd = cd

                    if imsize > 100000:
                        progct = { "cmd":"nrm_depot",
                                    "subcmd":"display_status",
                                    "status_msg": "produce %s" % imext
                                  }
                        pmsg = json.dumps(progct)
                        self.ws.send(pmsg)
                        
                    if  nd.shape[0] > 1000:
                        bd = nd[::3,::3]
                    else:
                        bd = nd
                    a = imshow(bd, interpolation = "none", extent=[0,nd.shape[1], nd.shape[0], 0])
                    
                    if imsize > 100000:
                        progct = { "cmd":"nrm_depot",
                                    "subcmd":"display_status",
                                    "status_msg": "transfering"
                                  }
                        pmsg = json.dumps(progct)
                        self.ws.send(pmsg)
                    savefig(imgname, bbox_inches='tight', dpi=32)     
                image = open(imgname)
                done = False
                
                imdata = image.read()
                datastr = b64encode(imdata)
                cmdct = {   
                            "num_nonzero":numnz,
                            "cmd":"nrm_depot",
                            "subcmd":"display",
                            "answering":msg,
                            "data64":"data:image/%s;base64,%s" % (imext,datastr)
                        }
                msg = json.dumps(cmdct)
                self.ws.send(msg)
                return
            elif subcmd == "local_data":
                print "client wants local_data description"
                ldata = {}
                for root,dirs, files in os.walk("."):
                    ldata["root"] = root
                    ldata["dirs"] = dirs
                    ldata["files"] = files
                    datasets = []
                    datasets_ct = {}
                    setrefs = []
                    ldata["datasets"] = datasets
                    print "sous45:", dw_info
                    ldata["datawarehouse"] = dw_info
                    for fil in files:
                        ext = os.path.splitext(fil)[1]
                        
                        #print "sous37: ext", ext
                        ext_type = None
                        if len(ext):
                            ext = ext[1:]
                            # setref pairing
                            # put setrefs in secondary list to check later
                            if ext == "setref":
                                setrefs.append(fil)
                            elif ext in generaldata._data_object_classes:
                                ext_type = ".".join(generaldata._data_object_classes[ext])
                                ds = {
                                        "filename":fil,
                                        "ext_type":ext_type
                                     }
                                imgname = "%s.%s" % (fil,imext)
                                if os.path.exists(imgname):
                                    ds["img_exists"] = True
                                else:
                                    ds["img_exists"] = False
                                datasets.append(ds)
                                # @@ISSUE?: fil cannot already be in dict right?
                                datasets_ct[fil] = ds
                    ## don't recurse into subdirectores
                    break
                
                #setrefs
                for fil in setrefs:
                    rawname = fil[:-7]
                    if rawname in datasets_ct:
                        datasets_ct[rawname]["has_setref"] = True
                        datasets_ct[rawname]["setref_name"] = fil
                        srfile = open(fil)
                        srstr = srfile.read()
                        srfile.close()
                        setrefct = json.loads(srstr)
                        datasets_ct[rawname]["setref"] = setrefct
                
                ldata["cmd"] = "nrm_depot"
                ldata["subcmd"] = "local_data"
                
                mtxt = json.dumps(ldata)
                #print "sous37:",mtxt
                self.ws.send(mtxt)
                return
        elif cmd == "run_recipe":
            # build commands
            cmdargs = ["kit"]
            opts = msg["options"]
            positional = opts["args"]
            del opts["args"]
            
            for key in opts:
                if opts[key] != True:
                    if len(key) == 1:
                        opt = "-%s '%s'" % (key, opts[key])
                    else:
                        opt = "--%s '%s'" % (key, opts[key])
                else:
                    if len(key) == 1:
                        opt = "-%s" %  key
                    else:
                        opt = "--%s" % key  
                cmdargs.append(opt)
           
            cmdargs.extend(positional)
            
            args = cmdargs
            cmdline = " ".join(cmdargs)
            print "sous35: args:", args
            
            proc = Popen(cmdline, shell=True, stdout = PIPE, stderr = PIPE)
            
            print "process = ", proc
            self._ra_proc = proc
            self._ra_stdin = proc.stdin
            self._ra_stdout = proc.stdout
            self._ra_stderr = proc.stderr
            
            done = False
            i = 0
            while not done:
                #sleep(.1)
                proc.stdout.flush()
                #print "reading stdout"
                #buf = proc.stdout.readline().strip()
                #buf = proc.stdout.read(10)
                buf = ""
                #print "reading stderr"
                errbuf = proc.stderr.readline().strip()
                #errbuf = proc.stderr.read(10)
                #
                #print "sous: |%s|%s|%s "%(  buf,errbuf, proc.poll())
                cmdct = {   
                            "cmd":"nrm_log",
                            "stdio": buf,
                            "stderr": errbuf,
                            "ansi":"%s%s" % ( buf, errbuf)
                        }
                msg = json.dumps(cmdct)
                #self.ws.send(buf)
                #self.ws.send(errbuf)
                self.ws.send(msg)
                #print "sous:poll"
                if proc.poll() != None:
                    done = True
                i+=1
                
            print "END CONNECTION"
            return 

    def on_close(self, reason):
        print "closing", reason

WebSocketServer(
    ('', 8000),
    Resource({'/': EchoApplication})
).serve_forever()
# kit -r showInputs -t TABLE 15_01_29.09_08.insta_ctf.nitro.dmo
