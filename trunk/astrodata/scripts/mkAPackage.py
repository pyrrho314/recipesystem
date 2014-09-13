#!/usr/bin/env python

import os, sys
from optparse import OptionParser
from shutil import copytree
import shutil

tempPack = "astrodata_Sample"
tempPackTag = "_Sample"

parser = OptionParser(usage="mkAPackage [options] <template name>")

parser.add_option("-w", dest="whereout", default=".")
(options, args) = parser.parse_args()

if len(args) <1:
    print "You need to give a name for the new template"
    parser.print_help()
    sys.exit(1)
else:
    outpackTag = "_"+args[0]
outdir = os.path.abspath(options.whereout)

dirname = os.path.dirname(__file__)
dirname = os.path.abspath(os.path.join(dirname, "../samples/"))
templatePack = os.path.join(dirname, tempPack)
outPack = os.path.abspath(os.path.join(outdir, "astrodata"+outpackTag))

print ("cloning\n\t%s \nto\n\t%s\n" %(templatePack, outPack))
copytree(templatePack, outPack)
print ("...copied, renaming subdirectories")

movelist = []
print "outpack",outPack

for root, dirs, files in os.walk(outPack):
    if (".svn" in root):
        shutil.rmtree(root)
        # print "skipping svn directory", root

for root, dirs, files in os.walk(outPack):
    print "in %s" % root,
    gotone = False
    for fil in dirs:
        # print "dir:",tempPackTag in fil, fil, outpackTag, fil.replace(tempPackTag, outpackTag)   
        if tempPackTag in fil:
            if gotone == False:
                print ""
            gotone = True
            newfil = fil.replace(tempPackTag, outpackTag)
            print "\tplanning to move %s -> %s" % ( fil, newfil )
            fullfil = os.path.join(root, fil)
            fullnewfil = os.path.join(root, newfil)
            #shutil.move(fullfil, fullnewfil)
            movelist.insert(0,(fullfil,fullnewfil))
    for fil in files:
        # print "dir:",tempPackTag in fil, fil, outpackTag, fil.replace(tempPackTag, outpackTag)   
        if tempPackTag in fil:
            if gotone == False:
                print ""
            gotone = True
            newfil = fil.replace(tempPackTag, outpackTag)
            print "\tplanning to move %s -> %s" % ( fil, newfil )
            fullfil = os.path.join(root, fil)
            fullnewfil = os.path.join(root, newfil)
            #shutil.move(fullfil, fullnewfil)
            movelist.insert(0,(fullfil,fullnewfil))              
    if not gotone:
        print "... nothing to change"
infostr = "Renaming %d Directories" % len(movelist)
print "="*len(infostr)
print infostr
print "-"*len(infostr)

for mov in movelist:
        print "Moving %s" % mov[0]
        print "   --> %s" % mov[1]
        shutil.move(mov[0], mov[1])