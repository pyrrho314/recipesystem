#!/usr/bin/env python

import sys
from astrodata import mkcalciface
from astrodata.mkcalciface import DD
from datetime import datetime

import glob

dlfilenames = glob.glob("DescriptorsList*.py")

if len(dlfilenames)>1:
    print "MORE THAN ONE FILE, should be just one:"
    print repr(dlfilenames)
    sys.exit()
elif len(dlfilenames)<1:
    print "NO FILES of form DescriptorsList*.py found"
    sys.exit()

dlfilename = dlfilenames[0]

cib = open(dlfilename)
cibsrc = cib.read()
cib.close()
d = {"DescriptorDescriptor":DD, 
     "DD":DD,
     "datetime":datetime
    }
ddlist = eval(cibsrc, d)

dcbody = mkcalciface.mk_calc_iface_body(ddlist)



print dcbody
