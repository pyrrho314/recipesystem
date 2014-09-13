# 
# Interactive test
#    
# Start python
#    
# Read the input image using pyfits

import pyfits as pf

spectrum = pf.getdata('gn20120419_104_sci.fits')

# Run deripple, returning the ndarray 'outspectrum'
outspectrum = deripple(spectrum, start=550, end=610, period=2)

