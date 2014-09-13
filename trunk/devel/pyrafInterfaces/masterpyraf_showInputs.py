from piro import mkRO
from tman import TaskManager


def showInputs():
    ro = mkRO(astrotype="GEMINI", reuseRO=True)
    ro.runstep("showInputs", ro.context)


tm = TaskManager()
tm.createTask(showInputs)
