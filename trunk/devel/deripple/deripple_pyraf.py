import pyfits as pf
from pyraf import iraf

# Load script deripple from module deripple.py
# The module must be made available to the PYTHON environment
# via PYTHONPATH environment variable.
from deripple import deripple

def _deripple (infile, outfile, start, end, period):
    """
      Script to run deripple under pyraf.

      # Start pyraf
      pyraf
    
      # Set deripplepath to where the deripple.py module is
      # located.
      set deripplepath=/home/nzarate/deripple/

      # Define the pyraf task
      task deripple = deripplepath$deripple.cl
      
      # Load the task into the pyraf environment
      deripple

      # Run the task
      deripple(infile='data/xobj_comb.fits',outfile='out.fits',start=550 ,end=610 ,period=2)
    """
      
    
    infits = pf.open(infile)
    spectrum = infits['SCI'].data
    print 'SPEin',spectrum.shape

    outspectrum = deripple (spectrum,start,end,period)

    # Create a new file with the input PHU and the SCI header
    phu = infits[0]
    pf.writeto(outfile, None, header=infits[0].header, clobber=True)
    pf.append(outfile, outspectrum, header=infits['SCI'].header)

    infits.close()
    print 'SPEout:::',outspectrum.shape
 
parfile = iraf.osfn('deripplepath$deripple.par')
t = iraf.IrafTaskFactory(taskname='deripple', value=parfile, function=_deripple)

