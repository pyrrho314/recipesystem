#!/usr/bin/env python

from astrodata import AstroData
from astrodata import memorytrack as mem
import time
import numpy as np
import os, sys
import inspect

import argparse

parser = argparse.ArgumentParser()

options, args =  parser.parse_known_args()

control = {}
for arg in args:
    if arg.startswith("-"):
        control[arg[1:]] = False
    if arg.startswith("+"):
        control[arg[1:]] = True

def lineno():
    """Returns the current line number in our program."""
    return "line #%d" % inspect.currentframe().f_back.f_lineno
    
    
if os.path.exists("memtrack.latest"):
    os.remove("memtrack.latest")

line = 0;

mem.memtrack("testmemprofile", lineno()); line += 1

ad = AstroData("data/N20131223S0243.fits")

mem.memtrack("testmemprofile", lineno()); line += 1

if True:
    datas = []
    for ext in ad:
        mem.memtrack("data pull in %s,%d" % (ext.extname(), ext.extver()), line); line += 1
        #print np.median(ext.data)
        #malloc = np.ones((2000,2000), dtype = np.float32)
        datas.append(ext.data)
    
    
mem.memtrack("before ad.write(..)", lineno());

ad.write("tmp.fits", clobber = True)

mem.memtrack("before ad.close()", lineno());

ad.close()

mem.memtrack("before del(ad)", lineno()); 

del(ad)

mem.memtrack("sleep .25 seconds", lineno()); 
time.sleep(.25)    
mem.memtrack("end sleep .25 seconds", lineno()); 
