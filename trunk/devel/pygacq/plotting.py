import matplotlib
matplotlib.use('TkAgg')

try:
    import matplotlib.pyplot as plt
except ImportError:
    import ctypes
    ctypes.CDLL("/opt/local/lib/libpng.dylib", mode=ctypes.RTLD_GLOBAL)
    import matplotlib.pyplot as plt

import os
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

import threading
import time

from Queue import Queue

_queue = Queue()

def is_fast_test_mode():
    if "GACQUI" in os.environ:
        uichoice = os.environ["GACQUI"]
        if "fast_test" in uichoice.lower():
            return True

    return False

FAST_TEST_MODE = is_fast_test_mode()

def put_plot(func, args):
    _queue.put((func, args))

def put_poison_pill():
    _queue.put(None)

class QueuePoller(object):
    def __init__(self, queue, parent, interval):
        self._queue = queue
        self._parent = parent
        self._interval = interval
        self._parent.after(self._interval, self.on_timer)
        self._all_figures = []

    def close(self):
        for fig in self._all_figures:
            plt.close(fig)
        self._parent.destroy()

    def on_timer(self):
        if not self._queue.empty():
            obj = self._queue.get()
            if obj is None:
                self.close()
                return

            func, args = obj
            fig = func(*args)
            self._all_figures.append(fig)

        self._parent.after(self._interval, self.on_timer)

def fast_test_mainloop():
    while True:
        obj = _queue.get()
        if obj is None:
            return
    
def tk_mainloop():
    root = Tk.Tk()
    poller = QueuePoller(_queue, root, 100)
    root.withdraw()
    
    Tk.mainloop()

if FAST_TEST_MODE:
    mainloop = fast_test_mainloop
else:
    mainloop = tk_mainloop
