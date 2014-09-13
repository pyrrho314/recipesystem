from subprocess import Popen, call
from CFSconfig import *
import os,sys
import difflib

sourcepys = os.path.join(FSDIR, "*.py")
destpys   = os.path.join(ADDIR, "*.py")

sourcestat = open("sourcesvnstat.txt", "w+")
call(" ".join(["svn","status", "--verbose", sourcepys]),shell=True, stdout=sourcestat)
sourcestat.close()
deststat = open("destsvnstat.txt", "w+")
call(" ".join(["svn","status", "--verbose", destpys ]), shell= True, stdout=deststat)
deststat.close()

sourcestat = open("sourcesvnstat.txt")
deststat = open("destsvnstat.txt")

sourcesl = sourcestat.readlines()
destsl   = deststat.readlines()

sourcestat.close()
deststat.close()

def removeUnk(alist):
    newlist = []
    for i in range(0, len(alist)):
        line = alist[i]
        if line[0]!="?":
            pyfile = os.path.basename(line[40:])
            newlist.append(pyfile)
        
    return newlist
sourcesl = removeUnk(sourcesl)
destsl   = removeUnk(destsl)


differ = difflib.Differ()

statusdiff = list(differ.compare(sourcesl, destsl))

missingfile = False
extrafile = False
for diff in statusdiff:
    if diff[0]=="-":
        print "File %s missing from %s" % (diff[2:].strip(),ADDIR)
        print "   This will not be automatically added, use svn add"
        missingfile = True
    if diff[0] == "+":
        print "File %s PRESENT in destination but not in source:"    
        extrafile == True
        
if missingfile:
    print "Svn Add Block to Cut and Paste"
    print "------------------------------"
    print "cd %s" % ADDIR
    for diff in statusdiff:
        if diff[0]=="-":
            print "svn add %s" % diff[2:]
else:
    print "Versions compare well.  NOTE: you must check in the destination version"
    print "svn status", ADDIR
    call (["svn", "status", ADDIR])
    print "----to cut and paste----"
    print "cd %s" % ADDIR
    print "svn ci"
    
    
