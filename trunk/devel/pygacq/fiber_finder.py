import sys
from astrodata import AstroData

import pyfits as pf
import numpy as np

from acquisitionimage import AcquisitionImage
from integratedfieldunit import *
import userinterface as ui

class FiberField(FiberStamp):
    def __init__(self, detsec):
        origin=(0, 0)
        FiberStamp.__init__(self, detsec, detsec.get_data(), origin)


def fiber_finder(fiberfield):
    fiber_center = None
    while fiber_center is None:
        print "== Press 'a' on this fiber, and then 'q' when happy ==" 
        listener = FiberDetectionListener(fiberfield, verbose=True)
        listener.start()
        fiber_center = listener.get_fiber_center()
    return fiber_center


acqimage = AcquisitionImage(sys.argv[1])

if acqimage.focal_plane_mask().upper() != "IFU-2" or acqimage.grating().upper() != "MIRROR":
    print "Should only be run on a IFU-2 flat non-dispersed image"
    sys.exit(1)

fiberfield = FiberField(acqimage.get_full_field_of_view())
ui.display(fiberfield.get_data(), zscale=True)

fname = sys.argv[2]
newfname = fname
if not os.path.exists(fname):
    fname = None
fiberbundles = FiberBundleCollection(fname, acqimage)

try:
    for bundle in fiberbundles.get_bundles():
        if not bundle.has_cached_fiber_positions():
            print "== MDF record for the first fiber for bundle '%s' ==" % bundle.get_block()
            print fiberbundles.get_column_header()
            print fiberbundles.format_record(bundle.get_first_fiber())
            
            fiber_center = fiber_finder(fiberfield)
            bundle.set_first_fiber_center(fiber_center)


            print "== MDF record for the last fiber for bundle '%s' ==" % bundle.get_block()
            print fiberbundles.get_column_header()
            print fiberbundles.format_record(bundle.get_last_fiber())
        
            fiber_center = fiber_finder(fiberfield)
            bundle.set_last_fiber_center(fiber_center)

        positions = bundle.get_fiber_positions()
        ui.mark_points(positions)
        print fiberbundles.get_new_column_header()
        for fiber, pos in zip(bundle.get_fibers(), bundle.get_fiber_positions()):
            print fiberbundles.format_fiber(fiber, pos)

finally:
    fiberbundles.write_new_table(newfname)
