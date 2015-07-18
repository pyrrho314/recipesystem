from multiprocessing import Process, Queue
from Queue import Empty
from datetime import datetime, timedelta
import os
import re
import subprocess
import getpass

from time import sleep
from astrodata import Lookups
from astrodata.adutils import ksutil as ks
from astrodata.adutils.dwutil.dwsettings import ingest_sources
from astrodata.adutils.dwutil.dwsettings import warehouse_packages
from astrodata.adutils.dwutil.dwsettings import package_dict

POLL_RESOLUTION = .5 # in seconds
 
class ShelfWatch(object):
    source = None
    frequency = 10
    ingest_package = None
    publish_package = None
    prefix = None
    _last_check = None
    _current_manifests = None
    
    def __init__(self, source = None):
        self._current_manifests = []
        if source:
            self.source = source
            if "elements" in self.source:
                elements = self.source["elements"]
            else:
                elements = {}
            if "shelf_name" not in elements:
                elements["shelf_name"] = source["ingest_shelf"]
            if "frequency" in source:
                self.frequency = source["frequency"]
            if "user" not in source:
                source["user"] = getpass.getuser()
            if "home" not in source:
                if "HOME" in os.environ:
                    source["home"] = os.environ["HOME"]
            print "d_p8:", ks.dict2pretty("source", source)
            in_pt = source["ingest_package_type"]
            in_shelf = source["ingest_shelf"]
            in_package = package_dict[in_pt]()
            self.prefix = in_package.get_store_prefix( elements = elements )
            print "d_p30: source prefix/key:", self.prefix
            self.ingest_package = in_package
            
            if "publish_package_type" in source:
                pub_pt = source["publish_package_type"]
                pub_package = package_dict[pub_pt]()
                self.publish_package = pub_package
                
    def do_check(self):
        ### get files
        now = datetime.now()
        self._last_check = now
        ret = {}
        slist = self.ingest_package.get_store_list(elements = self.source["elements"])
        for filnam in slist:
            basename = os.path.basename(filnam)
            if len(basename)>0:
                if "transferred_files" not in ret:
                    ret["transferred_files"] = []
                indivpkg = self.ingest_package.__class__(storename = filnam)
                indivpkg.deliver_from_warehouse(move = "_ingested")
                ret["transferred_files"].append(filnam)
        
        ## RUN COMMANDS
        # if document is a known type (see warehouse_daemon.py in lookups)
        xfers = None
        if "transferred_files" in ret:
            xfers = ret["transferred_files"]
        if xfers:
            if "commands" in self.source:
                commands = self.source["commands"]
                for command in commands:
                    patt = command["file_pattern"]
                    for filename in xfers:
                        if re.match(patt, filename):
                            for command_line in command["command_lines"]:
                                vargs = {"dataset": os.path.basename(filename)}
                                command_line = command_line.format(**vargs)
                                print "running: %s" % command_line
                                cmdparts = command_line.split()
                                exit_code = subprocess.call(cmdparts)
                                print "   exit_code = %s" % exit_code
        
        if not ret:
            return None
        else:
            return ret
        
    def is_due_for_check(self):
        now = datetime.now()
        if not self._last_check:
            return True
        if now - self._last_check > timedelta(seconds=self.frequency):
            return True
        else:
            return False
            
###### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ######
# NOTE: The process control is done at the module level. Perhaps this should  #
#       be wrapped in a class, but I want a singleton so all processes are    #
#       managed in one place.                                                 #
###### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ######

procs_by_shelf = {}      
queue_by_shelf = {}
source_by_shelf = {}
if ingest_sources:
    for source in ingest_sources:
        source_by_shelf[source["ingest_shelf"]] = source

def start_ingestion_processes():
    global procs_by_shelf, queue_by_shelf
    
    for shelfname in source_by_shelf:
        source = source_by_shelf[shelfname]
        cmd_q = Queue()
        result_q = Queue()
        queues = {"cmd_queue":cmd_q,
                  "result_queue":result_q                
                 }
        proc = Process( target=ingestion_loop,
                        args=(source,), 
                        kwargs= queues
                      )
        procs_by_shelf[shelfname] = proc
        queue_by_shelf[shelfname] = queues
        proc.start()
        
def stop_ingestion_processes():
    global procs_by_shelf, queue_by_shelf
    for shelf_name in queue_by_shelf:
        print "d_p111 Stopping Ingestion Process %s" %shelf_name
        cmd_q = queue_by_shelf[shelf_name]["cmd_queue"]
        cmd_q.put({"command":"stop"})
        
def tend_process_queue():
    global procs_by_shelf, queue_by_shelf
    done = False
    while not done:
        sleep(POLL_RESOLUTION)
        for shelf_name in queue_by_shelf:
            result_q = queue_by_shelf[shelf_name]["result_queue"]
            try:
                mess = result_q.get(False)
            except Empty:
                mess = None
            if mess and ("result" in mess) and mess["result"]:
                print ks.dict2pretty("d_p92 main proc Queue.get: %s" % shelf_name, mess)
    
def ingestion_loop(source, **args):
    """Meant to ingest_report in ingestion_loop_iterator(source):
        sleep(.5)
    """
    watch = ShelfWatch(source)
    cmd_q = args["cmd_queue"]
    result_q = args["result_queue"]
    done = False
    for tx in ingestion_loop_iterator(watch):
        # print ks.dict2pretty("d_p93: tx", tx)
        result_q.put(tx)
        try:
            cmd = cmd_q.get(False)
        except Empty:
            cmd = None
        if cmd :
            print "d_p146:",cmd
            if "command" in cmd:
                #possibly the iterator should handle this
                command = cmd["command"]
                if command == "stop":
                    print "d_p148: stopping ingestion_loop"
                    break;
        sleep(POLL_RESOLUTION)
        
def ingestion_loop_iterator(watch):
    done = False
    while not done:
        result = None
        watchdue = watch.is_due_for_check()
        if watchdue:
            result = watch.do_check()
        #print "d_p58: due?", watchdue
        if result:
            yield {"state":"watching",
               "was_due":watchdue,
               "prefix":watch.prefix,
               "result":result
                  }
        else:
            yield None
    yield {"state":"exiting"}
