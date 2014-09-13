# Import core Python modules
from __future__ import print_function
import math

# Gemini modules
from gempy.library.astrotools import get_fitted_function

# gacq modules
from gacqlogging import *
from util import get_window
import userinterface as ui
import plotting

# Import scientific modules
import numpy as np
import matplotlib.pyplot as plt

class Selection(object):
    def __init__(self, measured_center, window, fitted_function):
        self.measured_center = measured_center
        self.window = window
        self.fitted_function = fitted_function

    def get_window(self):
        return self.window

    def get_center(self):
        return self.measured_center

    def get_fwhm(self):
        """ return FWHM in pixels """
        return self.fitted_function.get_fwhm()

    def get_rsquared(self):
        return self.fitted_function.get_rsquared()

    def get_fitted_function(self):
        return self.fitted_function

    def __cmp__(self, rhs):
        rhs_rsquared = None
        if rhs is not None:
            rhs_rsquared = rhs.get_rsquared()
        return cmp(self.get_rsquared(), rhs_rsquared)

def plot_selection(selection, pixscale, cursor_listener=None):
    plt.ion()

    # when in debug mode, render into separate windows
    if is_debug_mode(): 
        fig = plt.figure()
        fignum = fig.number
    else:
        fignum = 0
        fig = plt.figure(num=fignum)
        plt.clf()

    fitted_function = selection.get_fitted_function()
    aximage = plt.matshow(fitted_function.get_stamp_data(), fignum=fignum, cmap=plt.cm.gist_earth_r)

    # add a color bar
    fig = aximage.get_figure()
    fig.canvas.set_window_title('GACQuisition')
    fig.colorbar(aximage)

    ax = fig.gca()

    # flip the y axis
    ax.invert_yaxis()

    # put x-axis labels on the bottom
    xaxis = ax.get_xaxis()
    xaxis.tick_bottom()

    # only tick the left side to match
    yaxis = ax.get_yaxis()
    yaxis.tick_left()

    ymin, ymax, xmin, xmax = selection.get_window()
    # set the labels to the original image, convert to 1-based indicing
    xticks = [int(xmin + t + 1) for t in ax.get_xticks()]
    ax.set_xticklabels(xticks)

    yticks = [int(ymin + t + 1) for t in ax.get_yticks()]
    ax.set_yticklabels(yticks)

    plt.contour(fitted_function.get_stamp_data(), cmap=plt.cm.copper)
    #plt.contour(fitted_function.get_model_function(), cmap=plt.cm.copper)

    msg = """
    center = %.1f, %.1f
    fwhm=%.2f arc
    R^2=%.2f"""  % (tuple(selection.get_center()) + 
                    (fitted_function.get_fwhm() * pixscale,
                     fitted_function.get_rsquared()))

    x_center, y_center = fitted_function.get_center() 
    plotted_point = plt.plot([x_center], [y_center], "ro")
    debug("x_width, y_width = %f, %f pix" % fitted_function.get_width())

    plt.text(0.95, 0.05,
             msg,
             fontsize=16,
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=ax.transAxes)

    if cursor_listener is not None:
        def onclick(event):            
            debug('...mouse click=%r, x=%r, y=%r, xdata=%r, ydata=%r' %
                  (event.button, event.x, event.y, event.xdata, event.ydata))

            if event.xdata is None or event.ydata is None:
                return

            while plotted_point:
                point = plotted_point.pop()
                try:
                    point.remove()
                except ValueError as err:
                    pass

            new_point = plt.plot([event.xdata], [event.ydata], "ro")
            plotted_point.extend(new_point)
            plt.draw()

            position = np.array((xmin + event.xdata + 1,
                                 ymin + event.ydata + 1))

            position += cursor_listener.get_offset()

            xpos, ypos = position
            frame = 1
            keystroke = "x"
            cursor = ui.CursorPosition(xpos, ypos, frame, keystroke)
            cursor_listener.put_ui_command(cursor_listener.send_cursor_position, (cursor, ))
    
        fig.canvas.mpl_connect('button_press_event', onclick)
    
    plt.draw()
    return fig

def plot_radial(selection, pixscale):
    plt.ion()

    # when in debug mode, render into separate windows
    if is_debug_mode(): 
        fig = plt.figure()
    else:
        fignum = 1
        fig = plt.figure(num=fignum)
        plt.clf()

    # get the function to plot
    fitted_function = selection.get_fitted_function()
    data = fitted_function.get_stamp_data()
    x_center, y_center = fitted_function.get_center()

    # calculate radial distances
    yindices, xindices = np.indices(data.shape)
    distances = np.sqrt((xindices - x_center)**2 + (yindices - y_center)**2)
    maxdist = distances.max()
    distances *= pixscale

    # flatten 2d space into positive 1d space and plot
    xpoints = distances.flatten()
    ypoints = data.flatten()
    plt.plot(xpoints, ypoints, 'o')

    # construct a 1 dimensional approximation of the 2D moffat
    background = fitted_function.get_background()
    peak = fitted_function.get_peak()
    width = np.mean(np.array(fitted_function.get_width()))
    beta = fitted_function.get_beta()
    def moffat(x):
        return background + peak * (1 + ((x / width)**2.0))**(-beta)

    # plot the 1d moffat
    xcoords = np.linspace(0.0, maxdist)
    plt.plot(xcoords * pixscale, moffat(xcoords), 'red')

    # add some useful information
    fig.canvas.set_window_title('GACQuisition')
    plt.xlabel('arcseconds')
    plt.ylabel('pixel value')
    
    msg = """
    center = %.1f, %.1f
    fwhm=%.2f arc
    R^2=%.2f"""  % (tuple(selection.get_center()) + 
                    (fitted_function.get_fwhm() * pixscale,
                     fitted_function.get_rsquared()))

    ax = fig.gca()
    plt.text(0.95, 0.95,
             msg,
             fontsize=16,
             horizontalalignment='right',
             verticalalignment='top',
             transform=ax.transAxes)

    plt.draw()
    return fig

def assert_int(obj):
    assert isinstance(obj, int)

def get_selection_peak(input_data, curcrds, pixscale,
                       verbose=False, pixel_buffer=None, cursor_listener=None, radial=False):
    seeing_estimate = 0.8
    default_fwhm = seeing_estimate / pixscale

    if pixel_buffer is None:
        pixel_buffer = max(int(round(default_fwhm)), 1)

    best_selection = None
    prev_shape = (None, None)
    
    for i in range(10):
        debug("...using a pixel_buffer =", pixel_buffer)
        xmin, xmax, ymin, ymax = get_window(input_data, curcrds, pixel_buffer)

        stamp_data = input_data[ymin:ymax, xmin:xmax]
        window = (ymin, ymax, xmin, xmax)

        map(assert_int, window)

        debug("...moffat optimization for the data in the window [%i:%i,%i:%i]" % window)

        ydim, xdim = stamp_data.shape
        assert ydim != 0
        assert xdim != 0

        if prev_shape == stamp_data.shape:
            debug("...already attempted to fit a function to a window this size, terminating search.")
            break
        prev_shape = stamp_data.shape

        # fit the subset of data to a moffat function
        fitted_function = get_fitted_function(stamp_data, default_fwhm)
        
        center_x, center_y = fitted_function.get_center()
        # convert to 1-based indicing and original image
        measured_center = np.array([xmin + center_x + 1, ymin + center_y + 1])
        selection = Selection(measured_center, window, fitted_function)

        if is_debug_mode():
            plotting.put_plot(plot_selection, (selection, pixscale, cursor_listener))

        # note, the first time through the loop best_selection is None, this is ok!
        best_selection = max(best_selection, selection)

        if fitted_function.get_success() <= 3:
            debug("...moffat optimization converged.")
            break

        rsquared = fitted_function.get_rsquared()
        if rsquared >= 0.90:
            debug("...moffat optimization did not converge, but R^2 is greater than 0.9 (%.2f)" % rsquared)
            break

        debug("...moffat optimization did not converge.")
        mindim = min(input_data.shape)
        if mindim > 300:
            if pixel_buffer > int(float(mindim) * 0.1):
                debug("...window is already greater than 20% of the field of view, not going to try again.")
                break
        
        debug("...doubling the window size and trying again...")
        pixel_buffer *= 2

    if verbose and (not is_debug_mode() or radial):
        plotting.put_plot(plot_selection, (best_selection, pixscale, cursor_listener))
        if radial:
            plotting.put_plot(plot_radial, (best_selection, pixscale))

    return best_selection

def mark_selection(coords, fwhm):
    mark = ui.marker(coords, color=ui.GREEN)
    if math.isnan(fwhm) or math.isinf(fwhm):
        fwhm = 10.0
    radius = abs(fwhm)
    circle = ui.circle(coords, radius)
    return [mark, circle]

class SelectionCursorListener(ui.CursorListener):
    def __init__(self, acqimage, verbose, mark_center=None,
                 pixel_buffer=None, circle_radius=None, offset=np.array((0.0, 0.0))):
        self.acqimage = acqimage
        self.verbose = verbose
        self.pixel_buffer = pixel_buffer
        self.circle_radius = circle_radius
        self.offset = offset
        self.undo = []
        self.custom_center = None
        self.object_position = None

        if mark_center is not None:
            nx, ny = mark_center
            debug("...marked center nx=", nx, "  ny=", ny)

            self.marked_center = ui.marker(mark_center, ui.MAGENTA)

    def get_offset(self):
        return self.offset

    def has_custom_center(self):
        return self.custom_center is not None

    def get_custom_center(self):
        assert self.has_custom_center()
        return self.custom_center.get_position()

    def get_object_coords(self):
        if self.object_position is None:
            error("No objects selected, I guess the universe has gone cold...")
            raise SystemExit
        return self.object_position.get_position()

    def send_cursor_position(self, pos):
        key = pos.get_keystroke()
        if   key == "c":
            self.handle_custom_center(pos)
            self.custom_center = pos
        elif key == "x":
            pos = self.handle_exact_location(pos)
            self.object_position = pos
        elif key == "e":
            self.handle_plot_fit(pos)
        elif key in ["a", "r"]:
            pos = self.handle_fit(pos)
            self.object_position = pos
        else:
            warning("Unrecognized keystroke '%s', doing nothing..." % key)
    
    def handle_exact_location(self, pos):
        self.undo_previous_selection()
        overlay = ui.marker(pos.get_position(), color=ui.GREEN)
        self.undo = [overlay]
        info("Object selected at (%.1f, %.1f)" % tuple(pos.get_position()))

        xpos, ypos = pos.get_position() - self.get_offset()
        return ui.CursorPosition(xpos, ypos, pos.get_frame(), pos.get_keystroke())

    def handle_custom_center(self, pos):
        self.undo_marked_center()
        self.marked_center = ui.marker(pos.get_position(), color=ui.MAGENTA)

    def handle_plot_fit(self, pos):
        display_plot = self.verbose
        if pos.get_keystroke() in ['e', 'r', 'a']:
            display_plot = True

        radial = pos.get_keystroke() == 'r' 
        selection = get_selection_peak(self.acqimage.get_science_data(),
                                       pos.get_position() - self.get_offset(),
                                       self.acqimage.binned_pixel_scale(),
                                       verbose=display_plot,
                                       pixel_buffer=self.pixel_buffer,
                                       cursor_listener=self,
                                       radial=radial)

        fwhm = selection.get_fwhm() * self.acqimage.binned_pixel_scale()
        args = tuple(selection.get_center()) + (fwhm,)
        info("Object selected at (%.1f, %.1f) with FWHM of %.3f arcseconds" % args)
        
        return selection

    def handle_fit(self, pos):
        self.undo_previous_selection()
        selection = self.handle_plot_fit(pos)

        coords = selection.get_center()
        radius = self.circle_radius
        if radius is None:
            radius = selection.get_fwhm()
        
        self.undo = mark_selection(coords + self.get_offset(), radius)
        
        xpos, ypos = coords 
        return ui.CursorPosition(xpos, ypos, pos.get_frame(), pos.get_keystroke())

    def undo_previous_selection(self):
        for overlay in self.undo:
            ui.undo(overlay)
        self.undo = []

    def undo_marked_center(self):
        if getattr(self, "marked_center", None) is not None:
            ui.undo(self.marked_center)
        self.marked_center = None
        
