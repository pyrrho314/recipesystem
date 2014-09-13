import os

import plotting

FAST_TEST_MODE = plotting.is_fast_test_mode()

NUMDISPLAY_MODE = False
if "GACQUI" in os.environ:
    uichoice = os.environ["GACQUI"]
    if "numdisplay" in uichoice.lower():
        NUMDISPLAY_MODE = True

if FAST_TEST_MODE:
    import userinterfacebase as ui
elif NUMDISPLAY_MODE:
    import ds9userinterface as ui
else:
    import pyds9userinterface as ui

CursorPosition = ui.CursorPosition
CursorListener = ui.CursorListener
add_cursor_position = ui.add_cursor_position
clear_cursor_positions = ui.clear_cursor_positions

RED = "RED"
BLUE = "BLUE"
GREEN = "GREEN"
MAGENTA = "MAGENTA"

_colors = ["RED",
           "BLUE",
           "GREEN",
           "MAGENTA"]

__all__ = _colors + ["CursorPosition",
                     "CursorListener",
                     "add_cursor_position",
                     "clear_cursor_positions"]

def _get_function(funcname):
    def newfunc(*args, **kwargs):
        if FAST_TEST_MODE:
            return
        
        if "color" in kwargs:
            kwargs["color"] = getattr(ui, kwargs["color"])

        newargs = []
        for arg in args:
            if arg in _colors:
                arg = getattr(ui, arg)
            newargs.append(arg)

        modfunc = getattr(ui, funcname)
        return modfunc(*newargs, **kwargs)
    return newfunc

current_module = __import__(__name__)
for funcname in ["mark_points",
                 "marker",
                 "circle",
                 "polyline",
                 "clear",
                 "undo",
                 "display"]:
   __all__.append(funcname)
   setattr(current_module, funcname, _get_function(funcname))
