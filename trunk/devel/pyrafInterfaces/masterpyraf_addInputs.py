from piro import mkRO
from tman import TaskManager


def addInputs(args):
    ro = mkRO(astrotype="GEMINI", reuseRO=True)
    # some stuff goes in here to pull out something
    ro.runstep("addInputs", ro.context)


tm = TaskManager()
tm.createTask(addInputs)
