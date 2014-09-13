import sys
from astrodata import Errors
from astrodata.adutils import gemLog
from gempy import gemini_tools as gt
from primitives_GMOS import GMOSPrimitives

import registration_functions as rf

class GMOS_IMAGE_alignment_Primitives(GMOSPrimitives):
    """ 
    This is the class of all primitives for the GMOS level of the type 
    hierarchy tree.  It inherits all the primitives to the level above
    , 'GEMINIPrimitives'.

    The list of primitives in this class is: 
    correctWCSToReferenceImage
    alignToReferenceImage
    """
    astrotype = 'GMOS_IMAGE'
    
    def init(self, rc):
        GMOSPrimitives.init(self, rc)
        return rc    

    def alignToReferenceImage(self, rc):
        """
        This primitive applies the transformation encoded in the input images
        WCSs to align them with a reference image, in reference image pixel
        coordinates.  The reference image is taken to be the first image in
        the input list.

        By default, the transformation into the reference frame is done via
        interpolation.  The interpolator parameter specifies the interpolation 
        method.  The options are nearest-neighbor, bilinear, or nth-order 
        spline, with n = 2, 3, 4, or 5.  If interpolator is None, 
        no interpolation is done: the input image is shifted by an integer
        number of pixels, such that the center of the frame matches up as
        well as possible.  The variance plane, if present, is transformed in
        the same way as the science data.  

        The data quality plane, if present, must be handled a little
        differently.  DQ flags are set bit-wise, such that each pixel is the 
        sum of any of the following values: 0=good pixel,
        1=bad pixel (from bad pixel mask), 2=nonlinear, 4=saturated, etc.
        To transform the DQ plane without losing flag information, it is
        unpacked into separate masks, each of which is transformed in the same
        way as the science data.  A pixel is flagged if it had greater than
        1% influence from a bad pixel.  The transformed masks are then added
        back together to generate the transformed DQ plane.
        
        In order not to lose any data, the output image arrays (including the
        reference image's) are expanded with respect to the input image arrays.
        The science and variance data arrays are padded with zeros; the DQ
        plane is padded with ones.  

        The WCS keywords in the headers of the output images are updated
        to reflect the transformation.

        :param interpolator: type of interpolation desired
        :type interpolator: string, possible values are None, 'nearest', 
                            'linear', 'spline2', 'spline3', 'spline4', 
                            or 'spline5'

        :param suffix: string to add on the end of the input filenames to 
                       generate output filenames
        :type suffix: string
        
        """
        log = gemLog.getGeminiLog(logType=rc['logType'],logLevel=rc['logLevel'])
        log.debug(gt.log_message("primitive", "alignToReferenceImage", 
                                 "starting"))
        adoutput_list = rf.align_to_reference_image(
                                         adinput=rc.get_inputs(style='AD'),
                                         interpolator=rc['interpolator'])
        rc.report_output(adoutput_list)
        yield rc

