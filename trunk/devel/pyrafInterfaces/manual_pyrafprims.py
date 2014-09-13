import sys, os


from pyraf import iraf
from astrodata import AstroData

# do showStack instead
from astrodata.RecipeManager import RecipeLibrary, ReductionContext
from astrodata import ReductionObjects # to stuff module level, should use func
from astrodata.ReductionObjects import command_clause
from astrodata.adutils import gemLog
from astrodata import Proxies


# .par file
# name,type,mode,default,min,max,prompt

reductionObject = None
adccpid = None
prs = None

## will be shared
####BEGIN### COPIED FROM REDUCE.PY ####
# GLOBAL/CONSTANTS (could be exported to config file) 
#cachedirs = [".reducecache",
#             ".reducecache/storedcals",
#             ".reducecache/storedcals/storedbiases",
#             ".reducecache/storedcals/storedflats",
#             ".reducecache/storedcals/retrievedbiases",
#             ".reducecache/storedcals/retrievedflats",                        
#             ]
#CALDIR = ".reducecache/storedcals"
#cachedict = {} # constructed below             
#for cachedir in cachedirs:
#    if not os.path.exists(cachedir):                        
#        os.mkdir(cachedir)
#    cachename = os.path.basename(cachedir)
#    if cachename[0] == ".":
#        cachename = cachename[1:]
#    cachedict.update({cachename:cachedir})
    
#calindfile = "./.reducecache/calindex.pkl"
#stkindfile = "./.reducecache/stkindex.pkl" # copied from reduce.py
####END#### COPIED FROM REDUCE.PY ####
log = gemLog.getGeminiLog(logLevel=6)

def mkRO(dataset="", astrotype="", args = None, argv = None):
    global reductionObject, adccpid
    
    if reductionObject == None:
        adccpid = Proxies.start_adcc()
        reduceServer = Proxies.ReduceServer()
        # Playing with module, should use access function
        ReductionObjects.prs = Proxies.PRSProxy.get_adcc(reduce_server=reduceServer)
        
        rl = RecipeLibrary()
        if dataset != "":
            ad = AstroData(dataset)
            ro = rl.retrieve_reduction_object(ad)
        elif astrotype != "":
            ad = None
            ro = rl.retrieve_reduction_object(astrotype = astrotype)

        # using standard command clause supplied in RecipeLibrary module
        ro.register_command_clause(command_clause)
        rc = ReductionContext()
        rc.ro = ro
        # rc.stackFile =  File("stackIndexFile", stkindfile)
        if args:
            rc.add_input(ad)
        
        reductionObject = ro
    else:
        ro = reductionObject
        rc = ro.context

    if args:
        rc.addInputs(args)
    rc.update(argv)
    ro.init(rc)
    
    return reductionObject



def getCalibration( *args, **argv):

    ro = mkRO(astrotype="GEMINI", args=args, argv=argv)
    
    ro.runstep("getCalibration", ro.context)


def setStackable(*args, **argv):
    ro = mkRO(astrotype="GEMINI", args=args, argv=argv)
    
    ro.runstep("setStackable", ro.context)


def showStackable(*args, **argv):
    ro = mkRO(dataset, args, argv)
    
    ro.runstep("showStackable", ro.context)

def addInputs(files, **argv):
    print repr(files), repr(argv)
    rcparmsdict = {"files":files}
    print repr(rcparmsdict)
    ro = mkRO(astrotype="GEMINI", argv=rcparmsdict)
    
    ro.runstep("addInputs", ro.context)
    
def showInputs(*args, **argv):
    #task = showInputs.task
    #parlist =  task.getParList()
        
    ro = mkRO(astrotype="GEMINI", args=args, argv=argv)
    
    ro.runstep("showInputs", ro.context)
    
class TaskManager(object):
    def createTask(self,function):
        fname = function.__name__
        task = iraf.IrafTaskFactory(   taskname=fname,
                            value=parfile, 
                            function=function)
        function.task = task
        par = task.getParList()[0]
        print repr(dir(par))
        
    
import os
currentdir = os.path.dirname(__file__)

parfile = iraf.osfn(os.path.join(currentdir, "showInputs.par"))
t = iraf.IrafTaskFactory(   taskname="showInputs", 
                            value=parfile, 
                            function=showInputs)

parfile = iraf.osfn(os.path.join(currentdir, "showStackable.par"))
t = iraf.IrafTaskFactory(   taskname="showStackable", 
                            value=parfile, 
                            function=showStackable)
    
    
parfile = iraf.osfn(os.path.join(currentdir, "setStackable.par"))
t = iraf.IrafTaskFactory(   taskname="setStackable", 
                            value=parfile, 
                            function=showInputs)

parfile = iraf.osfn(os.path.join(currentdir, "addInputs.par"))
t = iraf.IrafTaskFactory(   taskname="addInputs", 
                            value=parfile, 
                            function=addInputs)
       
