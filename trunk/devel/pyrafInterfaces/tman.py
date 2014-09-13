from pyraf import iraf
import os

class TaskManager(object):
    currentDir = None

    def __init__(self):
        self.currentDir = os.getcwd()

    def createTask(self,function):
        fname = function.__name__
        parfile = os.path.join(self.currentDir, fname+".par")
        task = iraf.IrafTaskFactory(   taskname=fname,
                            value=parfile, 
                            function=function)
        return task


