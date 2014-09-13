#!/usr/bin/env python
#
#
#                                                                     QAP Gemini
#
#                                                                        adcc.py
#                                                                        07-2013
# ------------------------------------------------------------------------------
# $Id: adcc.py 4507 2014-01-06 16:56:22Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 4507 $'[11:-2]
__version_date__ = '$Date: 2014-01-06 06:56:22 -1000 (Mon, 06 Jan 2014) $'[7:-2]
# ------------------------------------------------------------------------------
# Updated version of adcc:
#  -- Updated to argparse*
#  -- Updated to externalized ReduceInstanceManager
#  -- Updated PRSProxy version to 0.2 in get_version()
#
# * NOTE: argparse is used if available, but falls back to optparse, which is
#   depracated in 2.7. The get_args() function returns "args" in either case.
#   Interface to args attributes is the same.
#   eg.,
#   >>> args.httpport
#   8777
#   >>> args.reduceport
#   54530
# ------------------------------------------------------------------------------
"""Automated Dataflow Coordination Center"""

import os
import sys
import pickle
import select
import socket

from threading import Thread
from SimpleXMLRPCServer import SimpleXMLRPCServer

from astrodata import AstroData
from astrodata.adutils.reduceutils import prsproxyweb
from astrodata.adutils.reduceutils.prsproxyutil import calibration_search
from astrodata.adutils.reduceutils.reduceInstanceManager import ReduceInstanceManager

# ------------------------------------------------------------------------------
def buildArgParser():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="This is the proxy to PRS functionality, "+\
                            "also invoked locally, e.g. for calibration requests.")

    parser.add_argument("-v", "--verbose", 
                        dest="verbosity", action="store_true",
                        help="increase HTTP client messaging on adcc GET requests.")

    parser.add_argument("-i", "--invoked", 
                        dest = "invoked", action="store_true",
                        help = "Used by processes that invoke prsproxy, so "+\
                        "that PRS proxy knows when to exit. If not present, "+\
                        "the prsproxy registers itself and will only exit by "+\
                        "user control (or by os-level signal).")

    parser.add_argument("--startup-report", 
                        dest = "adccsrn", default=None, 
                        help = "Specify a file name for the adcc startup report")

    parser.add_argument("--preload", 
                        dest = "preload", action="store_true",
                        help = "Useful in proxy mode, where some information "+\
                        "otherwise produced during the first relevant request "+\
                        "is prepared prior to starting the HTTPServer.")

    parser.add_argument("--reload", 
                        dest = "reload", action="store_true",
                        help = "Just like --preload, but uses last, cached "+\
                        "(pickled) directory scan.")

    parser.add_argument("-r", "--reduce-port", 
                        dest = "reduceport", default=54530, type=int,
                        help="Option informs prsproxy of the port on which "+\
                        "reduce listens for xmlrpc commands.")

    parser.add_argument("-p", "--reduce-pid", 
                        dest ="reducepid", default=None, type=int,
                        help = "Option informs prsproxy of the reduce "+\
                        "application's PID.")

    parser.add_argument("-l", "--listen-port", 
                        dest = "listenport", default=53530, type=int,
                        help="prsproxy listener port for the xmlrpc server.")

    parser.add_argument("-w", "--http-port", 
                        dest = "httpport", default=8777, type=int,
                        help="Response port for the web interface. "+\
                        "i.e. http://localhost:<http-port>/")

    args = parser.parse_args()
    return args

# ------------------------------------------------------------------------------
def buildOptParser():
    from optparse import OptionParser

    parser = OptionParser()
    parser.set_description("This is the proxy to PRS functionality, also invoked "+\
                           "locally, e.g. for calibration requests.")

    parser.add_option("-v", "--verbose", 
                      dest="verbosity", action="store_true",
                      help="increase HTTP client messaging on adcc GET requests.")

    parser.add_option("-i", "--invoked", 
                      dest = "invoked", action = "store_true",
                      help = "Used by processes that invoke prsproxy, so "+\
                      "that PRS proxy knows when to exit. If not present, the "+\
                      "prsproxy registers itself and will only exit by user "+\
                      "control (or by os-level signal).")

    parser.add_option("--startup-report", 
                      dest = "adccsrn", default=None, 
                      help = "Specify a file name for the adcc startup report")

    parser.add_option("--preload", 
                      dest="preload", action="store_true",
                      help = "Useful in proxy mode, where some information "+\
                      "otherwise produced during the first relevant request "+\
                      "is prepared prior to starting the HTTPServer.")

    parser.add_option("--reload", 
                      dest="reload", action="store_true",
                      help = "Just like --preload, but uses last, cached "+\
                      "(pickled) directory scan.")

    parser.add_option("-r", "--reduce-port", 
                      dest = "reduceport", default=54530, type="int",
                      help="When invoked by reduce, this is used to inform "+\
                      "the prsproxy of the port on which reduce listens for "+\
                      "xmlrpc commands.")

    parser.add_option("-p", "--reduce-pid", 
                      dest ="reducepid", default=None, type="int",
                      help = "When invoked by reduce, this option is used to "+\
                      "inform the prsproxy of the reduce application's PID.")

    parser.add_option("-l", "--listen-port", 
                      dest = "listenport", default=53530, type="int", 
                      help="prsproxy listener port for the xmlrpc "+\
                      "server.")

    parser.add_option("-w", "--http-port", 
                      dest = "httpport", default=8777, type="int",
                      help="Response port for the web interface. "+\
                      "i.e. http://localhost:<http-port>/")

    args, pos_args = parser.parse_args()

    # No positional arguments to this interface.
    return args

# ------------------------------------------------------------------------------
# class ReduceInstanceManager imported from reduceInstanceManager.py in 
# astrodata.adutils.reduceutils
#
# definition was here.
# ------------------------------- UTILITY FUNCS --------------------------------

def getPersistDir(dirtitle = "adcc"):
    dirs = {"adcc":".adcc"}
    for dirt in dirs.keys():
        if not os.path.exists(dirs[dirt]):
            os.mkdir(dirs[dirt])
    return dirs[dirtitle]

def writeADCCSR(filename, vals=None):
    if filename == None:
        print "adcc93: no filename for sr"
        filename = ".adcc/adccReport"
    print "adcc95: startup report going to", filename
    sr = open(filename, "w+")
    if vals == None:
        sr.write("ADCC ALREADY RUNNING\n")
    else:
        sr.write(repr(vals))
    return

def get_args():
    try:
        args = buildArgParser()
    except ImportError:
        args = buildOptParser()
    return args

def get_version():
    version = [("PRSProxy","0.2")]
    print "prsproxy version:", repr(version)
    return version

# ----------------------------- END UTILITY FUNCS ------------------------------

# ------------------------------- FUTURE FEATURE -------------------------------
if False:
    try:
        from astrodata.FitsStorageFeatures import FitsStorageSetup
    except:
        import traceback
        traceback.print_exc()
    try:
        fss = FitsStorageSetup() # note: uses current working directory!!!
        if not fss.is_setup():
            print """Automated Dataflow Coordination Center:
            The local fits storage database has not been initialized for this
            directory.  This database allows reductions run in the same directory
            to share a common data repository, which can for example be used to
            retrieve best-fit calibrations.
    
            This initialization will only have to be executed one time
            for each working directory. 
            
            please wait...
            """
            fss.setup()
    except:
        msg = "Can't setup Local Fits Storage, some features not available."
        print msg
        print "Error reported:"
        print "-"*len(msg)
        import traceback
        traceback.print_exc()
        print "-"*len(msg)
        print "CONTINUING without Local Fits Storage Database."

# ---------------------------- END FUTURE FEATURE ------------------------------

# begin negotiated startup... we won't run if another adcc owns this directory
# could be done later or in lazy manner, but for now ensure this is present

racefile = ".adcc/adccinfo.py"
args     = get_args()

# caller lock file name
clfn = args.adccsrn
adccdir = getPersistDir()

if os.path.exists(racefile):
    print "ADCC307: adcc already has lockfile"
    from astrodata.Proxies import PRSProxy
    adcc = PRSProxy.get_adcc(check_once = True)
    if adcc == None:
        print "ADCC311: no adcc running, clearing lockfile"
        os.remove(racefile)
    else:
        print "ADCC314: adcc instance found running, halting"
        adcc.unregister()
        writeADCCSR(clfn)
        sys.exit()

# note: try to get a unique port starting at the standard port
findingPort = True
while findingPort:
    try:
        server = SimpleXMLRPCServer(("localhost", args.listenport), allow_none=True)
        print "PRS Proxy listening on port %d..." % args.listenport
        findingPort = False
    except socket.error:
        args.listenport += 1

# write out XMLRPC and HTTP port   
vals = {"xmlrpc_port": args.listenport,
        "http_port"  : args.httpport,
        "pid"        : os.getpid() }

#write racefile and ADCC Startup Report
ports = file(racefile, "w")
ports.write(repr(vals))
ports.close()
writeADCCSR(clfn, vals=vals)

server.register_function(get_version, "get_version")
server.register_function(calibration_search, "calibration_search")

# store the port

rim = ReduceInstanceManager(args.reduceport)
server.register_instance(rim)

if args.preload:
    print "adcc: scanning current directory tree"
    from astrodata.DataSpider import DataSpider
    ds = DataSpider(".")
    dirdict = ds.datasetwalk()
    persistpath = os.path.join(getPersistDir(), "dataspider_cache.pkl")
    ddf = open(persistpath, "w")
    pickle.dump(dirdict, ddf)
    ddf.close()
    print "adcc: done scanning current directory tree"
else:
    ds = None
    dirdict = None
    
if args.reload:
    from astrodata.DataSpider import DataSpider
    print "prsproxy: reloading result of previous scan of directory tree"
    ds = DataSpider(".")
    persistpath = os.path.join(getPersistDir(), "dataspider_cache.pkl")
    ddf = open(persistpath)
    dirdict = pickle.load(ddf)
    ddf.close()

# server.serve_forever(
# start webinterface

webinterface = True

if (webinterface):
    #import multiprocessing
    if ds and dirdict:
        web = Thread(None, prsproxyweb.main, "webface", 
                    kwargs = {"port": args.httpport,
                              "rim" : rim,
                              "dirdict": dirdict,
                              "dataSpider": ds,
                              "verbose": args.verbosity})
    else:
        web = Thread(None, prsproxyweb.main, "webface", 
                    kwargs = {"port":args.httpport,
                              "rim":rim,
                              "verbose": args.verbosity})

    web.start()
    
outerloopdone = False
while True:
    if outerloopdone:
        break
    try:
        while True:
            r, w, x = select.select([server.socket], [], [], 0.5)

            if r:
                server.handle_request() 

            if prsproxyweb.webserverdone:
                print "prsproxy exiting due to command vie http interface"
                print "number of reduce instances abandoned:", rim.numinsts
                outerloopdone = True
                break

            if (args.invoked and rim.finished):
                print "prsproxy exiting, no reduce instances to serve."
                outerloopdone = True
                prsproxyweb.webserverdone = True
                break

    except KeyboardInterrupt:
        if rim.numinsts > 0:
            # note: save reduce pide (pass in register) and 
            #       and check if pids are running!
            print "\nprsproxy: %d instances of reduce running" % rim.numinsts
            #below allows exit anyway
            outerloopdone = True
            prsproxyweb.webserverdone = True
            break

        else:
            print "\nprsproxy: exiting due to Ctrl-C"
            # this directly breaks from the outer loop but outerloopdone for clarity
            outerloopdone = True
            prsproxyweb.webserverdone = True
            # not needed os.kill(os.getpid(), signal.SIGTERM)
            break


if os.path.exists(racefile):
    os.remove(racefile)

sys.exit()
