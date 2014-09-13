#! /bin/tcsh
# 
# deripple test. This is a set of commands to run the deripple scripts.
#
# Start Python
# 
# import tests
#
setenv deripplepath '/home/nzarate/deripple'

$deripplepath/deripple.py --infile='gn20120419_104_sci.fits' \
      --outfile='out_104_csh_test.fits' --start=550 --end=610 --period=2


