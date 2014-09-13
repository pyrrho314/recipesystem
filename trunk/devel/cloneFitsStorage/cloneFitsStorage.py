from subprocess import Popen, call
import os, sys
from CFSconfig import *

#copy *.py from FSDIR to ADDIR

source = os.path.join(FSDIR, "*.py")
dest = ADDIR

import glob
from difflib import Differ
differ = Differ()
for sourcepy in glob.glob(source):
    print "checking..." + os.path.basename(sourcepy) +":",
    destpy = os.path.join(ADDIR, os.path.basename(sourcepy))
    if os.path.exists(destpy):
        # do diff
        s = open(sourcepy)
        d = open(destpy)
        sl = s.readlines()
        dl = d.readlines()
        s.close()
        d.close()
        result = differ.compare(sl,dl)
        difffound = False
        for line in list(result):
            if line[0] == "?" or line[0]== "+" or line[0]=="-":
                difffound = True
                break
    else:
        difffound = True
        
    if difffound:
        print "\n\t\tFiles different, copying..."
        call(["cp", sourcepy, ADDIR])
    else:
        print "Files identical, skipping."
        
    
sourcehtml = os.path.join(FSDIR, "htmldocroot/htmldocs")
destdr   = os.path.join(ADDIR, "htmldocroot")
desthtml = os.path.join(destdr, "htmldocs")

if not os.path.exists(destdr):
    os.mkdir(destdr)
if not os.path.exists(desthtml):
    os.mkdir(desthtml)
doccopy = "cp -r %s/* %s" % (sourcehtml,desthtml)
print "copying htmldocroot/htmldocs..."
print "\t%s" % doccopy
os.system(doccopy)

# make current version file
verfilename = os.path.join(ADDIR, "FSsvnver.ver")
ver = open(verfilename, "w+")

call (["svnversion", FSDIR], stdout = ver)

ver.close()
