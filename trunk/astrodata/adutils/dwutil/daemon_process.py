r"""
daemon_process.py

This module contains functions that help implement the datawarehouse 
daemon mode.

"""
import shutil
from multiprocessing import Process, Queue
from Queue import Empty
from datetime import datetime, timedelta
import os,sys
import re
import subprocess
import getpass
from glob import glob

from time import sleep
from astrodata import Lookups
from astrodata.adutils import ksutil as ks
from astrodata.adutils import termcolor as tc
from astrodata import ConfigSpace
from astrodata.adutils.dwutil.dwsettings import ingest_sources
from astrodata.adutils.dwutil.dwsettings import warehouse_packages
from astrodata.adutils.dwutil.dwsettings import package_dict

POLL_RESOLUTION = .5 # in seconds
def clear_directory(path):
    print "Removing contents of %s" % path
    for fil in os.listdir(path):
        if os.path.islink(fil):
            os.remove(fil)
        elif os.path.isdir(fil):
            shutil.rmtree(fil)
        elif os.path.isfile(fil):
            os.remove(fil)
            
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
            print "Daemon Process __init__: source prefix/key: %s (d_p68)" % self.prefix
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
                            if False: # if command["clean_working_directory"]:
                                # @@WARN: possibly dangerous, Ideally isolate
                                # the daemon with it's own account/permisions
                                rmcontents = os.getcwd()
                                shutil.rmtree(rmcontents)
                                os.mkdir(rmcontent)
                                
                            for command_line in command["command_lines"]:
                                vargs = {"dataset": os.path.basename(filename),
                                         "context": ConfigSpace.get_current_default_context()
                                        }
                                command_line = command_line.format(**vargs)
                                print tc.colored("-"*(len(command_line)+len("running:")+2), None, "on_green")
                                print tc.colored("running:", None, "on_green"), "%s" % command_line
                                print tc.colored("-"*(len(command_line)+len("running:")+2), None, "on_green")
                                cmdparts = command_line.split()
                                # check parts you want to glob
                                convparts = []
                                for part in cmdparts:
                                    if "*" in part:
                                        convparts.extend(glob(part))
                                    else:
                                        convparts.append(part)
                                        
                                exit_code = subprocess.call(convparts)
                                print "   exit_code = %s" % exit_code
                            
                            if command["clean_working_directory"]:
                                # @@WARN: possibly dangerous, Ideally isolate
                                # the daemon with it's own account/permisions
                                rmcontents = os.getcwd()
                                clear_directory(rmcontents)
                                
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
#       managed in one place, which modules perform well at.                  #
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

        local_wd = "workdir_for_shelf-%s" % shelfname
        local_wd = os.path.abspath(local_wd)
        
        queues = {"cmd_queue":cmd_q,
                  "result_queue":result_q,
                  "cwd":local_wd
                 }
        
        if not os.path.exists(local_wd):
            os.mkdir(local_wd)
            
        
        proc = Process( target=ingestion_loop,
                        args=(source,), 
                        kwargs= queues,
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
        try:
            sleep(POLL_RESOLUTION)
            for shelf_name in queue_by_shelf:
                result_q = queue_by_shelf[shelf_name]["result_queue"]
                try:
                    mess = result_q.get(False)
                except Empty:
                    mess = None
                if mess and ("result" in mess) and mess["result"]:
                    print ks.dict2pretty("d_p92 main proc Queue.get: %s" % shelf_name, mess)
        except KeyboardInterrupt:
            print "keyboard interupt, finishing..."
            done = True
            
            
def ingestion_loop(source, **args):
    """Meant to ingest_report in ingestion_loop_iterator(source):
        sleep(.5)
    """
    watch = ShelfWatch(source)
    cmd_q = args["cmd_queue"]
    result_q = args["result_queue"]
    working_dir = args["cwd"]
    os.chdir(working_dir)
    done = False
    try:
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
    except KeyboardInterrupt:
        print "\nIngestion Process Recieved KeyboardInterrupt..."
        return False
        
def ingestion_loop_iterator(watch):
    done = False
    while not done:
        result = None
        watchdue = watch.is_due_for_check()
        if watchdue:
            sys.stdout.write(tc.colored("\rchecking at %s (d_p243)%s" % (datetime.now(), ""), "blue") )  # "\033[J"),
            sys.stdout.flush()
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
