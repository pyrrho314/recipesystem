
from piro import mkRO
prims = ['addInputs', 'showInputs']

def addInputs(files, **argv):
    print "addInputs running...."
    task = addInputs.task
    
    #need to turn this into a function
    #parlist = task.getParList()
    #print "parlist:", repr(parlist), repr(dir(parlist))
    parmdict = {}
    
    #npar = len(parlist)
    #print "npar=",npar
    #for i in range(0,npar-2):
    #    par = parlist[i]
    #    parmdict.update({par.name:args[i]})
    #print repr(parmdict)
    parmdict = {"files":files}
    ro = mkRO(astrotype="GEMINI",argv=parmdict)
    ro.runstep("addInputs", ro.context)

def showInputs(*args,**argv):
    ro = mkRO(astrotype="GEMINI",args=args, argv=argv)
    ro.runstep("showInputs", ro.context)
