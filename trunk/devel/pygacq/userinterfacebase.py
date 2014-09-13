import collections
import threading
import numpy as np

from gacqlogging import *
import plotting

from Queue import Queue

_WAITING_FOR_UI_EVENT = threading.Event()
_QUEUE = Queue()

class CursorPosition(object):
    def __init__(self, xPos, yPos, frame, key):
        self.xPos = float(xPos)
        self.yPos = float(yPos)
        self.frame = frame
        self.key = key

    def get_position(self):
        return np.array([self.xPos, self.yPos])

    def get_keystroke(self):
        return self.key.lower()

    def get_frame(self):
        return self.frame

    def __repr__(self):
        return ("CursorPosition<" + repr(self.get_position()) + ", " +
                "Frame = " + repr(self.get_frame()) + ", " +
                "Keystroke = " + repr(self.get_keystroke()) + ">")

NEXT_CURSOR_POSITIONS = collections.deque()
def add_cursor_position(xPos, yPos, key):
    NEXT_CURSOR_POSITIONS.append(CursorPosition(xPos, yPos, "100", key))

def clear_cursor_positions():
    NEXT_CURSOR_POSITIONS.clear()

class CursorListener(object):
    def send_cursor_position(self, pos):
        raise RuntimeError("send_cursor_position is not implemented")

    def send_close(self):
        pass

    def has_cached_cursor_position(self):
        return bool(NEXT_CURSOR_POSITIONS)

    def get_cached_cursor_position(self):
        return NEXT_CURSOR_POSITIONS.popleft()

    def get_cursor_position(self):
        if self.has_cached_cursor_position():
            return self.get_cached_cursor_position()

        raise RuntimeError("get_cursor_position not implemented, or a test asked for too many cursor positions")
    
    def read_cursor_positions(self):
        exitcmds = set(['q',
                        'Q'])

        processuicmds = set(['g',
                             'G'])
        while True:
            cPos = self.get_cursor_position()
            if cPos is None:
                self.process_ui_queue()
                continue

            if cPos.get_keystroke()[0] in processuicmds:
                _WAITING_FOR_UI_EVENT.set()
                if _QUEUE.empty():
                    self.process_ui_one_command()
                else:
                    self.process_ui_queue()
                _WAITING_FOR_UI_EVENT.clear()
                continue
            
            if cPos.get_keystroke()[0] in exitcmds:
                self.process_ui_queue()
                self.send_close()
                break

            self.send_cursor_position(cPos)

    def put_ui_command(self, func, args):
        if not _WAITING_FOR_UI_EVENT.is_set():
            warning("Press 'g' in the DS9 window to update DS9 with the point selected in the contour plot")
        _QUEUE.put((func, args))

    def process_ui_one_command(self):
        func, args = _QUEUE.get()
        func(*args)

    def process_ui_queue(self):
        while not _QUEUE.empty():
            self.process_ui_one_command()
            
    def start(self):
        # spawn off a thread to interact with the gui and do the various interactive calculations
        thrd = ReadCursorPositionsThread(listener=self)
        thrd.start()

        # the main thread is then left in charge of the GUI the way Tk likes it
        plotting.mainloop()

        try:
            thrd.join_with_exception()
        except Exception, err:
            error(str(err))
            sys.exit(1)

class ReadCursorPositionsThread(threading.Thread):
    def __init__(self, listener):
        threading.Thread.__init__(self)
        self.listener = listener

    def run(self):
        try:
            try:
                self.listener.read_cursor_positions()
            finally:
                plotting.put_poison_pill()
        except:
            self.exc_info = sys.exc_info()
            raise self.exc_info[1], None, self.exc_info[2]

    def join_with_exception(self):
        self.join()
        if hasattr(self, "exc_info"):
            raise self.exc_info[1], None, self.exc_info[2]

RED = 0
BLUE = 0
GREEN = 0
MAGENTA = 0

def mark_points(points):
    raise RuntimeError("Should be implemented in user interface implementation")

def marker(point, color, mark="+"):
    raise RuntimeError("Should be implemented in user interface implementation")

def circle(center, radius, color):
    raise RuntimeError("Should be implemented in user interface implementation")

def polyline(points, color):
    raise RuntimeError("Should be implemented in user interface implementation")

def undo(*args):
    raise RuntimeError("Should be implemented in user interface implementation")

def display(*args, **kwargs):
    raise RuntimeError("Should be implemented in user interface implementation")

