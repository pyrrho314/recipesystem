#! /usr/bin/env python

import numpy as np

def deripple(spectrum, start, end, period):
    """
      Deripple removes a periodic instrument signal from the
      input 'spectrum'.

      Before running this script the user can select about 50 pixels
      of the spectrum which is free from features and shows the ripple
      signal, noting the starting and ending pixel numbers.

      The algorithm is as follows:
      - Within the window  spectrum[start:end] form a list of maximum
        pixel values and another of minimum pixel values. (for a period
        of 2 pixels). For periods >2 make as many lists as the value 
        of period and append the pixel in each corresponding slot.

      - Find the average of each list.
      - With the list of average, make a template periodic feature.
      - Replicate this template as many times as necessary to expand
        to the length of the input spectrum.
      - Divide the input spectrum by the expanded template and
        normalize to the level of the average input spectrum.
      
      spectrum:  1-Dimensional spectrum
      start:     (pixel number) Start of the window within
                 the spectrum. (0-based)
      end:       (pixel number) End of window within the
                 spectrum.  (0-based)
      period:    The period of the ripple in pixels.
    """
    
    if type(spectrum) is not np.ndarray:
        spectrum = np.asarray(spectrum)

    if len(spectrum.shape) != 1:
        raise ValueError('Input spectrum is not 1-Dimension')

    if not ((0<= start < len(spectrum)) and 
            (0< end <= len(spectrum)) and (end > start)):
        raise ValueError("Bad 'start' or 'end' limits")

    # get window section
    window = spectrum[start:end+1]     # If start,end are zero based

    nvals = (end-start)+1 # Number of values in the window.
    p = period            # Alias, otherwise the statements below 
                          # will be very noisy.

    # Make as many lists of pixels as the value of 'period'. Fill in
    # each list with the corresponding pixel value.
    vals = np.zeros(((nvals+p-1)/p,p),np.float32)
    for j,k in enumerate(range(0,len(window)-(p-1),p)):
        for i in range(p):
            vals[j,i] = window[k+i]

    # Now find the average of each layer
    averages = np.zeros(period,np.float32)
    for k in range(period):
        averages[k] = np.mean(vals[:,k])

    # Make a period
    ripple = [k for k in averages]
    
    # Make a copy
    ripples = ripple[:]
    # Replicate this to cover the input spectrum width

    for k in range((len(spectrum)+p-1)/p):
        ripples += ripple
        
    # Make ripples the same length as the spectrum
    rdiff = len(spectrum) - len(ripples)
    if rdiff > 0:
        ripples = ripples + ripple[:diff+1]
    elif rdiff < 0:
        ripples = ripples[:len(spectrum)]

    # Divide both arrays
    dspectrum = spectrum/ripples

    # normalize to spectrum level

    dspectrum = dspectrum*np.mean(ripple)

    return dspectrum

if __name__ == "__main__":
    """
      Running the script under the UNIX shell

      deripple --infile='data/xobj_comb.fits' --outfile='out.fits'
                --start=550 --end=610 --period=2
    """

    import sys
    import optparse
    import pyfits as pf

    # Parse input arguments
    usage = "usage: %prog --infile='in.fits' --outfile='out.fits' --start=550 --end=610 --period=2"

    # EXAMPLE:
    # deripple.py --infile='data/xobj_comb.fits' --outfile='out.fits'
    #               --start=550 --end=610 --period=2

    p = optparse.OptionParser(usage=usage, version='niciprepare_1.0')
    p.add_option('--infile', action='store', type='string',
                 default='', help='Input FITS file')
    p.add_option('--outfile', action='store', type='string',
                 default='./', help='Output FITS file')
    p.add_option('--start', action='store', type='int', 
                 help='start of window')
    p.add_option('--end', action='store', type='int',help="End of window")
    p.add_option('--period', action='store', type='int',help="ripple period")

    (par, args) = p.parse_args()

    if len(sys.argv) == 1:
       print '\n'
       p.print_help()
       sys.exit(0)

    infits = pf.open(par.infile)
    spectrum = infits['SCI'].data

    outspectrum = deripple(spectrum, par.start, par.end, par.period)
     # Create a new file with the input PHU and the SCI header
    phu = infits[0]
    pf.writeto(par.outfile, None, header=infits[0].header, clobber=True)
    pf.append(par.outfile, outspectrum, header=infits['SCI'].header)

    infits.close()

    

    # pf.writeto(par.ourfile,out, header)


