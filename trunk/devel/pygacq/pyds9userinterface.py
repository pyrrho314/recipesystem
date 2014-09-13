import ds9
import numpy as np
import userinterfacebase as ui
from gacqlogging import *

CursorPosition = ui.CursorPosition
add_cursor_position = ui.add_cursor_position
clear_cursor_positions = ui.clear_cursor_positions

INSTANCE = ds9.ds9()

class CursorListener(ui.CursorListener):
    def get_cursor_position(self):
        _redraw_everything()
        if self.has_cached_cursor_position():
            return self.get_cached_cursor_position()

        cursor = INSTANCE.get("imexam any coordinate image")

        parts  = cursor.split()
        if len(parts) != 3:
            msg = "Did not receive cursor coordinate from DS9, received '%r' instead" % cursor
            error(msg)
            raise RuntimeError(msg)

        key, xpos, ypos = parts 
        frame = 1
        
        return CursorPosition(xpos, ypos, frame, key)

RED = "red"
BLUE = "blue"
GREEN = "green"
MAGENTA = "magenta"

CURRENT_OVERLAYS = set()

def _redraw_everything():
    INSTANCE.set('regions delete all')
    for cmd in CURRENT_OVERLAYS:
        INSTANCE.set(cmd)

def mark_points(points):
    for pnt in points:
        marker(point, RED)
        
def marker(point, color, mark="+"):
    xpos, ypos = point
    cmd = 'regions command {text %f %f #text="%s" font="times 18 bold" color="%s"}' % (xpos, ypos, mark, color)

    CURRENT_OVERLAYS.add(cmd)
    return cmd

def circle(center, radius, color=GREEN):
    xpos, ypos = center
    cmd = 'regions command {circle %f %f %f # color="%s"}' % (xpos, ypos, radius, color)

    CURRENT_OVERLAYS.add(cmd)
    return cmd
             
def polyline(points, color):
    commands = []
    
    prev_pnt = None
    for pnt in points:
        if prev_pnt is not None:
            cmd = 'regions command {line %f %f %f %f # color="%s"}' % (tuple(prev_pnt) + tuple(pnt) + (color, ))
            CURRENT_OVERLAYS.add(cmd)
            commands.append(cmd)
        prev_pnt = pnt

    return commands

def _undo_one(arg):
    if isinstance(arg, str):
        CURRENT_OVERLAYS.remove(arg)
    else:
        for item in arg:
            _undo_one(item)
            

def undo(*args):
    for arg in args:
        _undo_one(arg)

def display(data, *args, **kwargs):
    _clear()
    f32data = np.array(data, dtype=np.float32)
    INSTANCE.set_np2arr(f32data)

def _clear():
    CURRENT_OVERLAYS.clear()
    _redraw_everything()
