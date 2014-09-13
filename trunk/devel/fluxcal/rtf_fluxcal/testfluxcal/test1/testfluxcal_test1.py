#
##############################################################################
#
# PyRAF Test Package - GMOS Set #1
#
# Test script for 'gemini.gmos.fluxcal' (All tests non-interactive.)
#
##############################################################################
#
# Gemini PYTHON package - Test script
#
# Normal-use-test for 'gemini.gmos.fluxcal
#
# Test 1: Run fluxcal
#
#
# Version:  01.11.2011 NZ
#
# Packages that should be loaded before running the script
#    None
#
# Other requirements:
#  The following python modules
#    1. os (loaded with driver)
#    2. shutil (loaded with regress2.py)
#
##############################################################################

# Set directories
#rawdir  = os.environ['fluxcal']+'orig/'
rawdir  = '/home/nzarate/zp/rtf/orig/'
#outdir  = os.environ['fluxcal']+'outdir/'

print ("######################### TEST 1 ############################")
print ("    FLUXCAL  of GMOS IMAGES")
print ("#############################################################")

# Delete files left over from previous run

filelist = ['mgS20101214S0041.fits','mgS20101214S0040.fits','mgS20101105S0128.fits','gS20101214S0041.fits','mgN20100930S0366.fits','gN20100930S0366.fits']
#filelist = open( 'todel', 'r'  )
for file in filelist:
    if os.path.isfile( file.strip() ) is True:
        print "deleting",file
        os.remove( file.strip() )
#filelist.close()

#Create list of input images and copy input images over
#inputlist = open( 'inlist','w+' )
for fi in filelist:
    #inputlist.write(fi) 
    shutil.copy( rawdir + fi, '.' )
#inputlist.close()

#for file in ['mgS20101214S0041.fits','mgS20101214S0040.fits','mgS20101105S0128.fits','gS20101214S0041.fits']:
for file in filelist:

    print "\n ***************************** FLUXCAL for ",file," **************************" 
    ff = fc.Fluxcal(file)
    ff.runFC()

print ("TEST 1: Completed.")

##############################################################################

