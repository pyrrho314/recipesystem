.. _functionality:

Basic Functionality

- When creating a Wavecal object the following happens:

 - Copies the input parameters to a dictionary 'params'.
 - Gets the input file Astrodata type.
 - Creates an F2, GMOS, GNIRS, NIRI or NIFS object according to the instrument type.
 - The __init__ for the instrument reads the instrument specific parameters from a dictionary 'INSTRUMENT' and checks each value againts the value of 'params', replacing any value into the former that the user has provided.

- Perform the wavelength calibration using the instrument object method 'wavecal()':

 - set_wavecal_data() sets the middle line (self.lpix) from the input image and look for the peaks pixel coordinates (self.xpeaks) up to a maximum of self.ntmax peaks.
 - find_wavesolution() sets and run the Matching ratios algorithm to find  the correct wavelength in the reference line list (self.cuar) to a peak in the arc peaks list (self.xpeaks) 
 - Votes. The matching algorithm does several passes to each triple of arc positions. For each match a vote is added to a matrix votes[peaks,references]. 
 - Fit a function to the array in votes. When all the possible passes are done in the previous step, the pairs (pixel,ref) from the votes matrix containing more than 10 votes are fit using a Legendre function of order 4.
 - Improve the fit. For each entry in the linelist that is within the wavelength range of the input arc spectrum find the pixel position. If an arc peak is found near this peak, get its wavelength via the computed fit function and if this value of closed to the reference, then add it to the list (pix_array, user_array). Fit again.


