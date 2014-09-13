import numdisplay
import userinterfacebase as ui
from gacqlogging import *

class CursorPosition(ui.CursorPosition):
    def get_position(self):
        """ Should theoretically subtract by 1.0, but the rest of gacq seems to depend on 1 based indicing """
        return ui.CursorPosition.get_position(self)

def add_cursor_position(xPos, yPos, key):
    ui.NEXT_CURSOR_POSITIONS.append(CursorPosition(xPos, yPos, "100", key))

clear_cursor_positions = ui.clear_cursor_positions

class CursorListener(ui.CursorListener):
    def get_cursor_position(self):
        if self.has_cached_cursor_position():
            return self.get_cached_cursor_position()

        try:
            cursor = numdisplay.readcursor(sample=0, timeout=0.01)
        except numdisplay.NoDataReady:
            return None

        if cursor == "EOF":
            msg = "EOF returned from DS9, please restart DS9"
            error(msg)
            raise RuntimeError(msg)

        parts  = cursor.split()
        if len(parts) != 4:
            msg = "Did not receive cursor coordinate from DS9, received '%r' instead" % cursor
            error(msg)
            raise RuntimeError(msg)

        return CursorPosition(*parts)

def ds9_round_int(x):
    return int(round(x))

def ds9_round_point(point):
    return map(ds9_round_int, point)

def ds9_round_points(points):
    return [ds9_round_point(crds) for crds in points]

RED = numdisplay.overlay.C_RED
BLUE = numdisplay.overlay.C_BLUE
GREEN = numdisplay.overlay.C_GREEN
MAGENTA = numdisplay.overlay.C_MAGENTA

def mark_points(points):
    points = ds9_round_points(points)
    color = numdisplay.overlay.C_RED
    numdisplay.overlay.point(points=points, color=color)

def marker(point, color, mark="+"):
    xpos, ypos = map(ds9_round_int, point)
    return numdisplay.overlay.marker(x=xpos,
                                     y=ypos,
                                     mark=mark,
                                     size=1,
                                     color=color)

def circle(center, radius, color=GREEN):
    center = ds9_round_point(center)
    return numdisplay.overlay.circle(center=tuple(center),
                                     radius=radius,
                                     color=color)


def polyline(points, color):
    points = ds9_round_points(points)
    return numdisplay.overlay.polyline(points=points, color=color)

def undo(*args):
    numdisplay.overlay.undo(*args)

def display(*args, **kwargs):
    numdisplay.display(*args, **kwargs)
