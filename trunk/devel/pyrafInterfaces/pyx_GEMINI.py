import pyraf_prims_GEMINI as ppG
import os

from pyraf import iraf

primslist = ppG.prims
    
class TaskManager(object):
    
    piroModule = None
    currentDir = None

    def __init__(self, pirofile):
        if pirofile[-3:].lower() == ".py":
            pirofile = pirofile[:-3]
            
        pif = None
        exec("import %s as pif" % pirofile)
        reload(pif)
        self.piroModule = pif
        self.currentDir = os.path.dirname(pif.__file__)
    def createTasks(self):
        for prim in self.piroModule.prims:
            pirofunc = eval("self.piroModule.%s" % prim)
            t = self.createTask(pirofunc)
            pirofunc.task = t

    def createTask(self,function):
        fname = function.__name__
        parfile = os.path.join(self.currentDir, fname+".par")
        task = iraf.IrafTaskFactory(   taskname=fname,
                            value=parfile, 
                            function=function)
        # function.task = task
        # par = task.getParList()[0]
        # print repr(dir(par))
        return task

tm = TaskManager("pyraf_prims_GEMINI.py")
tm.createTasks()

