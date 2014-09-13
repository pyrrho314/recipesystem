import sys
import numpy as np
import pywcs

from astrodata.adutils import gemLog
from gempy import gemini_tools as gt
from astrodata import Errors
from gempy import managers as man


# correct_wcs_to_reference_image and align_to_reference_image functions
# will go into gempy/science/registration

def correct_wcs_to_reference_image(adinput=None, 
                                   method='sources', fallback=None, 
                                   cull_sources=False,
                                   rotate=False, scale=False):

    """
    This function registers images to a reference image by correcting
    the relative error in their world coordinate systems.  The function
    uses points of reference common to the reference image and the
    input images to fit the input WCS to the reference one.  The fit
    is done by a least-squares minimization of the difference between
    the reference points in the input image pixel coordinate system.
    This function is intended to be followed by the align_to_reference_image
    function, which applies the relative transformation encoded in the
    WCS to transform input images into the reference image pixel
    coordinate system.

    The primary registration method is intended to be by direct mapping
    of sources in the image frame to correlated sources in the reference
    frame. This method fails when there are no correlated sources in the
    field, or when the WCSs are very far off to begin with.  As a back-up
    method, the user can try correcting the WCS by the shifts indicated 
    in the POFFSET and QOFFSET header keywords (option fallback='header'), 
    or by hand-selecting common points of reference in an IRAF display
    (option fallback='user').  By default, only the direct method is
    attempted, as it is expected that the relative WCS will generally be
    more correct than either indirect method.  If the user prefers not to
    attempt direct mapping at all, they may set method to either 'user'
    or 'header'.

    In order to use the direct mapping method, sources must have been
    detected in the frame and attached to the AstroData instance in an 
    OBJCAT extension.  This can be accomplished via the detectSources
    primitive.  Running time is optimal, and sometimes the solution is 
    more robust, when there are not too many sources in the OBJCAT.  Try
    running detectSources with threshold=20.  The solution may also be
    more robust if sub-optimal sources are rejected from the set of 
    correlated sources (use option cull_sources=True).  This option may
    substantially increase the running time if there are many sources in
    the OBJCAT.

    It is expected that the relative difference between the WCSs of 
    images to be combined should be quite small, so it may not be necessary
    to allow rotation and scaling degrees of freedom when fitting the image
    WCS to the reference WCS.  However, if it is desired, the options 
    rotate and scale can be used to allow these degrees of freedom.  Note
    that these options refer to rotation/scaling of the WCS itself, not the
    images.  Significant rotation and scaling of the images themselves 
    will generally already be encoded in the WCS, and will be corrected for
    when the images are aligned.

    The WCS keywords in the headers of the output images are updated
    to contain the optimal registration solution.

    Log messages will go to a 'main' type logger object, if it exists.
    or a null logger (ie. no log file, no messages to screen) if it does 
    not.

    :param adinput: images to register. Reference image is assumed to be
                  the first one in the list.  All images must have
                  only one SCI extension.
    :type adinput: AstroData objects, either a single instance or a list

    :param method: method to use to generate reference points. Options
                   are 'sources' to directly map sources from the input image
                   to the reference image, 'user' to select reference
                   points by cursor from an IRAF display, or 'header' to
                   generate reference points from the POFFSET and QOFFSET
                   keywords in the image headers.
    :type method: string, either 'sources', 'user', or 'header'

    :param fallback: back-up method for generating reference points.
                     if the primary method fails.  The 'sources' option
                     cannot be used as the fallback.
    :type fallback: string, either 'user' or 'header'.  

    :param cull_sources: flag to indicate whether sub-optimal sources should
                   be rejected before attempting a direct mapping. If True,
                   sources that are saturated, not well-fit by a Gaussian,
                   too broad, or too elliptical will be eliminated from
                   the list of reference points.
    :type cull_sources: bool
    
    :param rotate: flag to indicate whether the input image WCSs should
                   be allowed to rotate with respect to the reference image
                   WCS
    :type rotate: bool

    :param scale: flag to indicate whether the input image WCSs should
                  be allowed to scale with respect to the reference image
                  WCS.  The same scale factor is applied to all dimensions.
    :type scale: bool

    """

    # Instantiate log
    log = gemLog.getGeminiLog()

    # Ensure that adinput is not None and return
    # a list containing one or more AstroData objects
    adinput = gt.validate_input(adinput=adinput)

    # Keyword to be used for time stamp
    keyword = "REGISTER"

    adoutput_list = []
    try:

        if len(adinput)<2:
            raise Errors.InputError('At least two images must be provided.')
    
        if method is None:
            if fallback is None:
                raise Errors.InputError('Both method and fallback are None; ' +
                                        'not attempting WCS correction.')
            else:
                method = fallback

        n_test = []
        for ad in adinput:

            # Make sure all images have one science extension
            if len(ad['SCI'])!=1:
                raise Errors.InputError('Input images must have only one SCI ' +
                                        'extension.')

            # Get number of objects from OBJCAT
            num_cat = len(ad['OBJCAT'])
            if num_cat==0:
                n_obj = 0
            elif num_cat>1:
                raise Errors.InputError('Input images must have only one ' +
                                        'OBJCAT extension.')
            else:
                n_obj = len(ad['OBJCAT'].data)

            n_test.append(n_obj)


        if n_test[0]==0 and method=='sources':
            log.warning('No objects found in reference image.')
            if fallback is not None:
                log.warning('Only attempting indirect WCS alignment, ' +
                            'via ' + fallback + ' mapping')  
                method=fallback

            else:
                log.warning('WCS can only be corrected indirectly ' +
                            'and fallback method is set to None.  Not ' +
                            'attempting WCS correction.')
                return adinput


        # Reference image is first one supplied
        # (won't be modified)
        reference = adinput[0]
        adoutput_list.append(reference)

        # If no OBJCAT/no sources in reference image, or user choice,
        # use indirect alignment for all images at once
        if method=='header':
            reg_ad = _header_align(reference, adinput[1:])
            adoutput_list.extend(reg_ad)
        elif method=='user':
            reg_ad = _user_align(reference, adinput[1:], rotate, scale)
            adoutput_list.extend(reg_ad)
        elif method!='sources':
            raise Errors.InputError('Did not recognize method' + method)

        # otherwise try to do direct alignment for each image by correlating
        # sources in the reference and input images
        else:

            for i in range(1,len(adinput)):
            
                ad = adinput[i]

                if n_test[i] == 0:
                    log.warning('No objects found in '+ ad.filename)
                    if fallback is not None:
                        log.warning('Only attempting indirect WCS alignment, ' +
                                    'via ' + fallback + ' mapping')  
                        if fallback=='header':
                            adoutput = _header_align(reference, ad)
                        elif fallback=='user':
                            adoutput = _user_align(reference, ad, rotate, scale)
                        else:
                            raise Errors.InputError('Did not recognize ' +
                                                    'fallback method' + fallback)

                    else:
                        log.warning('WCS can only be corrected indirectly '+
                                    'and fallback=None. Not attempting WCS ' +
                                    'correction for ' + ad.filename)
                        continue
                else:
                    log.fullinfo('Number of objects in image %s: %d' %
                                 (ad.filename, n_test[i]))
        
                    log.status('Cross-correlating sources in %s, %s' %
                               (reference.filename, ad.filename))
                    obj_list = _correlate_sources(reference, ad, 
                                                  cull_sources=cull_sources)

                    n_corr = len(obj_list[0])

                    if n_corr==0:
                        log.warning('No correlated sources found.')
                        if fallback is not None:
                            log.warning('Only attempting indirect WCS ' +
                                        'alignment, via ' + fallback + 
                                        ' mapping')

                            if fallback=='header':
                                adoutput = _header_align(reference, ad)
                            elif fallback=='user':
                                adoutput = _user_align(reference, ad, 
                                                       rotate, scale)
                            else:
                                raise Errors.InputError('Did not recognize ' +
                                                        'fallback ' +
                                                        'method' + fallback)
                                              
                        else:
                            log.warning('WCS can only be corrected indirectly '+
                                        'and fallback=None. Not attempting ' +
                                        'WCS correction for ' + ad.filename)
                            continue
                    else:
                        log.fullinfo('Number of correlated sources: %d' % 
                                     n_corr)

                        # Check the fit geometry depending on the 
                        # number of objects
                        if n_corr == 1:
                            log.warning('Too few objects.  Setting ' +
                                        'rotate=False, ' +
                                        'scale=False')
                            rotate=False
                            scale=False

                        log.fullinfo('\nSources used to align frames:')
                        log.fullinfo('  %7s %7s %7s %7s\n%s' % 
                                     (' Ref. x','Ref. y',
                                      'Img. x','Img. y',
                                      '  '+'-'*31))
                        output_obj = zip(obj_list[0],obj_list[1])
                        for obj in output_obj:
                            obj_string = ('  %7.2f %7.2f %7.2f %7.2f' % 
                                          (obj[0][0],obj[0][1],
                                           obj[1][0],obj[1][1]))
                            log.fullinfo(obj_string)
                        log.fullinfo('')

                        adoutput = _align_wcs(reference, ad, [obj_list], 
                                              rotate=rotate, scale=scale)

                gt.mark_history(adinput=adoutput, keyword=keyword)
                adoutput_list.extend(adoutput)

        return adoutput_list

    except:
        # Log the message from the exception
        log.error(repr(sys.exc_info()[1]))
        raise


def align_to_reference_image(adinput, interpolator='linear'):
    """
    This function applies the transformation encoded in the input images
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

    :param adinput: list of images to align.  First image is taken to be
                  the reference image.
    :type adinput: list of AstroData objects

    :param interpolator: type of interpolation desired
    :type interpolator: string, possible values are None, 'nearest', 'linear',
                        'spline2', 'spline3', 'spline4', or 'spline5'

    """

    # instantiate log
    log = gemLog.getGeminiLog()

    # Ensure that adinput is not None and return
    # a list containing one or more AstroData objects
    adinput = gt.validate_input(adinput=adinput)

    # keyword to be used for time stamp
    keyword = "ALIGN"

    adoutput_list = []
    try:
        # check for at least two input images (first one is reference)
        if len(adinput)<2:
            raise Errors.InputError('At least two input images ' +
                                    'must be supplied.')

        # make sure all images have one science extension
        for ad in adinput:
            if len(ad['SCI'])!=1:
                raise Errors.InputError('Input images must have only one ' +
                                        'SCI extension.')

        # load ndimage package if there will be interpolation
        if interpolator=='None':
            interpolator = None
        if interpolator is not None:
            from scipy.ndimage import affine_transform

        # get reference WCS
        reference = adinput[0]
        ref_wcs = pywcs.WCS(reference['SCI'].header)
        ref_shape = reference['SCI'].data.shape
        ref_corners = get_corners(ref_shape)
        naxis = len(ref_shape)

        # first pass: get output image shape required to fit all data in output
        # by transforming corner coordinates of images
        all_corners = [ref_corners]
        xy_img_corners = []
        shifts = []
        for i in range(1,len(adinput)):

            ad = adinput[i]

            img_wcs = pywcs.WCS(ad['SCI'].header)

            img_shape = ad['SCI'].data.shape
            img_corners = get_corners(img_shape)

            xy_corners = [(corner[1],corner[0]) for corner in img_corners]
            xy_img_corners.append(xy_corners)

            if interpolator is None:
                # find shift by transforming center position of field
                # (so that center matches best)
                x1y1 = np.array([img_shape[1]/2.0,img_shape[0]/2.0])
                x2y2 = img_wcs.wcs_sky2pix(ref_wcs.wcs_pix2sky([x1y1],1),1)[0]

                # round shift to nearest integer and flip x and y
                offset = np.roll(np.rint(x2y2-x1y1),1)

                # shift corners of image
                img_corners = [tuple(offset+corner) 
                               for corner in img_corners]
                shifts.append(offset)
            else:
                # transform corners of image via WCS           
                xy_corners = img_wcs.wcs_sky2pix(ref_wcs.wcs_pix2sky(
                                                           xy_corners,0),0)

                img_corners = [(corner[1],corner[0]) for corner in xy_corners]

            all_corners.append(img_corners)


        cenoff = []
        out_shape = []
        for axis in range(naxis):
            # get output shape from corner values
            cvals = [corner[axis] for ic in all_corners for corner in ic]
            out_shape.append(int(max(cvals)-min(cvals)+1))

            # if just shifting, need to set centering shift for reference
            # image from offsets already calculated
            if interpolator is None:
                svals = [shift[axis] for shift in shifts]
                # include a 0 shift for the reference image
                # (in case it's already centered)
                svals.append(0.0)
                cenoff.append(-int(max(svals)))

        out_shape = tuple(out_shape)

        # if not shifting, get offset required to center reference image
        # from the size of the image
        if interpolator is not None:
            incen  = [0.5*(axlen-1) for axlen in ref_shape]
            outcen = [0.5*(axlen-1) for axlen in out_shape]
            cenoff = np.rint(incen) - np.rint(outcen)

        # shift the reference image to keep it in the center of the new array
        # (do the same for VAR and DQ)
        log.fullinfo('Growing reference image to keep all data; ' +
                     'centering data, and updating WCS to account ' +
                     'for shift')
        log.fullinfo('New output shape: '+repr(out_shape))

        ref_corners = [(corner[1]-cenoff[1]+1,corner[0]-cenoff[0]+1) # x,y
                       for corner in ref_corners]
        log.fullinfo('Setting AREA keywords in header to denote original ' +
                     'data area.')
        area_keys = []
        log.fullinfo('AREATYPE = "P4"     / Polygon with 4 vertices')
        area_keys.append(('AREATYPE','P4','Polygon with 4 vertices'))
        for i in range(len(ref_corners)):
            for axis in range(len(ref_corners[i])):
                key_name = 'AREA%i_%i' % (i+1,axis+1)
                key_value = ref_corners[i][axis]
                key_comment = 'Vertex %i, dimension %i' % (i+1,axis+1)
                area_keys.append((key_name,key_value,key_comment))
                log.fullinfo('%-8s = %7.2f  / %s' % 
                             (key_name, key_value,key_comment))

        for ext in reference:
            if ext.extname() not in ['SCI','VAR','DQ']:
                continue

            ref_data = ext.data

            trans_data = np.zeros(out_shape)

            # pad the DQ plane with 1 instead of 0
            if ext.extname()=='DQ':
                trans_data += 1.0

            trans_data[int(-cenoff[0]):int(ref_shape[0]-cenoff[0]),
                       int(-cenoff[1]):int(ref_shape[1]-cenoff[1])] = ref_data

            ext.data = trans_data

            # update the WCS in the reference image to account for the shift
            ext.set_key_value('CRPIX1', ref_wcs.wcs.crpix[0]-cenoff[1])
            ext.set_key_value('CRPIX2', ref_wcs.wcs.crpix[1]-cenoff[0])

            # set area keywords
            for key in area_keys:
                ext.set_key_value(key[0],key[1],key[2])


        # update the WCS in the PHU as well
        reference.phu_set_key_value('CRPIX1', ref_wcs.wcs.crpix[0]-cenoff[1])
        reference.phu_set_key_value('CRPIX2', ref_wcs.wcs.crpix[1]-cenoff[0])

        out_wcs = pywcs.WCS(reference['SCI'].header)

        adoutput_list.append(reference)


        # now transform the data
        for i in range(1,len(adinput)):

            log.status('\nStarting alignment for '+ adinput[i].filename)

            ad = adinput[i]

            sciext = ad['SCI']
            img_wcs = pywcs.WCS(sciext.header)
            img_shape = sciext.data.shape

            if interpolator is None:

                # recalculate shift from new reference wcs
                x1y1 = np.array([img_shape[1]/2.0,img_shape[0]/2.0])
                x2y2 = img_wcs.wcs_sky2pix(out_wcs.wcs_pix2sky([x1y1],1),1)[0]
            
                shift = np.roll(np.rint(x2y2-x1y1),1)
                if np.any(shift>0):
                    log.warning('Shift was calculated to be >0; ' +
                                'interpolator=None '+
                                'may not be appropriate for this data.')
                    shift = np.where(shift>0,0,shift)

                # update PHU WCS keywords
                log.fullinfo('Offsets: '+repr(np.roll(shift,1)))
                log.fullinfo('Updating WCS to track shift in data')
                ad.phu_set_key_value('CRPIX1', img_wcs.wcs.crpix[0]-shift[1])
                ad.phu_set_key_value('CRPIX2', img_wcs.wcs.crpix[1]-shift[0])

            else:
                # get transformation matrix from composite of wcs's
                # matrix = in_sky2pix*out_pix2sky (converts output to input)
                xy_matrix = np.dot(np.linalg.inv(img_wcs.wcs.cd),out_wcs.wcs.cd)

                # switch x and y for compatibility with numpy ordering
                flip_xy = np.roll(np.eye(2),2)
                matrix = np.dot(flip_xy,np.dot(xy_matrix,flip_xy))
                matrix_det = np.linalg.det(matrix)

                # offsets: shift origin of transformation to the reference pixel
                # by subtracting the transformation of the output reference
                # pixel and adding the input reference pixel back in
                refcrpix = np.roll(out_wcs.wcs.crpix,1)
                imgcrpix = np.roll(img_wcs.wcs.crpix,1)
                offset = imgcrpix - np.dot(matrix,refcrpix)

                # then add in the shift of origin due to dithering offset.
                # This is the transform of the reference CRPIX position,
                # minus the original position
                trans_crpix = img_wcs.wcs_sky2pix(out_wcs.wcs_pix2sky(
                                                 [out_wcs.wcs.crpix],1),1)[0]
                trans_crpix = np.roll(trans_crpix,1)
                offset = offset + trans_crpix-imgcrpix

                # Since the transformation really is into the reference
                # WCS coordinate system as near as possible, just set image
                # WCS equal to reference WCS
                log.fullinfo('Offsets: '+repr(np.roll(offset,1)))
                log.fullinfo('Transformation matrix:\n'+repr(matrix))
                log.fullinfo('Updating WCS to match reference WCS')
                ad.phu_set_key_value('CRPIX1', out_wcs.wcs.crpix[0])
                ad.phu_set_key_value('CRPIX2', out_wcs.wcs.crpix[1])
                ad.phu_set_key_value('CRVAL1', out_wcs.wcs.crval[0])
                ad.phu_set_key_value('CRVAL2', out_wcs.wcs.crval[1])
                ad.phu_set_key_value('CD1_1', out_wcs.wcs.cd[0,0])
                ad.phu_set_key_value('CD1_2', out_wcs.wcs.cd[0,1])
                ad.phu_set_key_value('CD2_1', out_wcs.wcs.cd[1,0])
                ad.phu_set_key_value('CD2_2', out_wcs.wcs.cd[1,1])
            

            # transform corners to find new location of original data
            data_corners = out_wcs.wcs_sky2pix(img_wcs.wcs_pix2sky(
                                                  xy_img_corners[i-1],0),1)
            log.fullinfo('Setting AREA keywords in header to denote original ' +
                         'data area.')
            area_keys = []
            log.fullinfo('AREATYPE = "P4"     / Polygon with 4 vertices')
            area_keys.append(('AREATYPE','P4','Polygon with 4 vertices'))
            for i in range(len(data_corners)):
                for axis in range(len(data_corners[i])):
                    key_name = 'AREA%i_%i' % (i+1,axis+1)
                    key_value = data_corners[i][axis]
                    key_comment = 'Vertex %i, dimension %i' % (i+1,axis+1)
                    area_keys.append((key_name,key_value,key_comment))
                    log.fullinfo('%-8s = %7.2f  / %s' % 
                                 (key_name, key_value,key_comment))

            for ext in ad:
                extname = ext.extname()

                if extname not in ['SCI','VAR','DQ']:
                    continue

                log.status('Transforming '+ad.filename+'['+extname+']')

                # Access pixel data
                img_data = ext.data

                if interpolator is None:
                    # just shift the data by an integer number of pixels
                    # (useful for noisy data, also lightning fast)

                    trans_data = np.zeros(out_shape)

                    # pad the DQ plane with 1 instead of 0
                    if extname=='DQ':
                        trans_data += 1

                    trans_data[int(-shift[0]):int(img_shape[0]
                                                  -shift[0]),
                               int(-shift[1]):int(img_shape[1]
                                                  -shift[1])] = img_data
                    
                    matrix_det = 1.0

                    # update the wcs to track the transformation
                    ext.set_key_value('CRPIX1', img_wcs.wcs.crpix[0]-shift[1])
                    ext.set_key_value('CRPIX2', img_wcs.wcs.crpix[1]-shift[0])

                else:
                    # use ndimage to interpolate values

                    # Interpolation method is determined by 
                    # interpolator parameter
                    if interpolator=='nearest':
                        order = 0
                    elif interpolator=='linear':
                        order = 1
                    elif interpolator=='spline2':
                        order = 2
                    elif interpolator=='spline3':
                        order = 3
                    elif interpolator=='spline4':
                        order = 4
                    elif interpolator=='spline5':
                        order = 5
                    else:
                        raise Errors.InputError('Interpolation method ' +
                                                interpolator +
                                                ' not recognized.')


                    if extname=='DQ':

                        # DQ flags are set bit-wise
                        # bit 1: bad pixel (1)
                        # bit 2: nonlinear (2)
                        # bit 3: saturated (4)
                        # A pixel can be 0 (good, no flags), or the sum of
                        # any of the above flags 
                        # (or any others I don't know about)

                        # unpack the DQ data into separate masks
                        # NOTE: this method only works for 8-bit masks!
                        unp = (img_shape[0],img_shape[1],8)
                        unpack_data = np.unpackbits(np.uint8(
                                                      img_data)).reshape(unp)
                    
                        # transform each mask
                        trans_data = np.zeros(out_shape)
                        for j in range(0,8):

                            # skip the transformation if there are no flags set
                            # (but always do the bad pixel mask because it is 
                            # needed to mask the part of the array that was
                            # padded out to match the reference image)
                            if not unpack_data[:,:,j].any() and j!=7:
                                # first bit is j=7 because unpack is backwards 
                                continue

                            mask = np.float32(unpack_data[:,:,j])

                            # if bad pix bit, pad with 1.  Otherwise, pad with 0
                            if j==7:
                                cval = 1
                            else:
                                cval = 0
                            trans_mask = affine_transform(mask, matrix,
                                                         offset=offset,
                                                         output_shape=out_shape,
                                                         order=order, cval=cval)
                            del mask; mask = None

                            # flag any pixels with >1% influence from bad pixel
                            trans_mask = np.where(np.abs(trans_mask)>0.01,
                                                  2**(7-j),0)

                            # add the flags into the overall mask
                            trans_data += trans_mask
                            del trans_data; trans_data = None

                    else: 

                        # transform science and variance data in the same way

                        cval = 0.0
                        trans_data = affine_transform(img_data, matrix,
                                                      offset=offset,
                                                      output_shape=out_shape,
                                                      order=order, cval=cval)

                    
                    # update the wcs
                    ext.set_key_value('CRPIX1', out_wcs.wcs.crpix[0])
                    ext.set_key_value('CRPIX2', out_wcs.wcs.crpix[1])
                    ext.set_key_value('CRVAL1', out_wcs.wcs.crval[0])
                    ext.set_key_value('CRVAL2', out_wcs.wcs.crval[1])
                    ext.set_key_value('CD1_1', out_wcs.wcs.cd[0,0])
                    ext.set_key_value('CD1_2', out_wcs.wcs.cd[0,1])
                    ext.set_key_value('CD2_1', out_wcs.wcs.cd[1,0])
                    ext.set_key_value('CD2_2', out_wcs.wcs.cd[1,1])
           
                    # set area keywords
                    for key in area_keys:
                        ext.set_key_value(key[0],key[1],key[2])


                ext.data = trans_data

            
            # if there was any scaling in the transformation, the
            # pixel size will have changed, and the output should
            # be scaled by the ratio of input pixel size to output
            # pixel size to conserve the total flux in a feature.
            # This factor is the determinant of the transformation
            # matrix.
            if (1.0-matrix_det)>1e-6:
                log.fullinfo('Multiplying by %f to conserve flux' %
                             matrix_det)

                # Allow the arith toolbox to do the multiplication
                # so that variance is handled correctly
                ad.mult(matrix_det)

            # Add time stamp to PHU
            gt.mark_history(adinput=ad, keyword=keyword)
            adoutput_list.append(ad)

        return adoutput_list

    except:
        # Log the message from the exception
        log.error(repr(sys.exc_info()[1]))
        raise





def _cull_sources(ad, img_obj):
    """
    This function takes a list of identified sources in an image, fits
    a Gaussian to each one, and rejects it from the list if it is not
    sufficiently star-like.  The criteria for good sources are that they
    must be fittable by a Gaussian, not be too near the edge of the frame,
    have a peak value below saturation (as defined in the header of
    the image), have ellipticity less than 0.25, and have FWHM less than
    2.4 arcsec.  The return value is a list of the objects that meet these
    criteria, with their positions updated to the fit center.

    :param ad: input image
    :type ad: AstroData instance

    :param img_obj: list of [x,y] positions for sources detected in the
                    input image
    :type img_obj: list
    """

    import scipy.optimize

    if len(ad['SCI'])!=1:
        raise Errors.InputError('Reference image must have only ' +
                                  'one SCI extension.')

    img_data = ad['SCI'].data

    # first guess at background is mean of whole image
    default_bg = img_data.mean()

    # first guess at fwhm is .8 arcsec
    default_fwhm = .8/float(ad.pixel_scale())

    # stamp is 2 times this size on a side
    aperture = default_fwhm 

    # for rejecting saturated sources
    saturation = ad.saturation_level() 

    good_source = []
    for objx,objy in img_obj:

        # array coords start with 0
        objx-=1
        objy-=1

        xlow, xhigh = int(round(objx-aperture)), int(round(objx+aperture)), 
        ylow, yhigh = int(round(objy-aperture)), int(round(objy+aperture)), 

        if (xlow>0 and xhigh<img_data.shape[1] and 
            ylow>0 and yhigh<img_data.shape[0]):
            stamp_data = img_data[ylow:yhigh,xlow:xhigh]
        else:
            # source is too near the edge, skip it
            continue

        # starting values for Gaussian fit
        bg = default_bg
        peak = stamp_data.max()
        x_ctr = (stamp_data.shape[1]-1)/2.0
        y_ctr = (stamp_data.shape[0]-1)/2.0
        x_width = default_fwhm
        y_width = default_fwhm
        theta = 0.

        if peak >= saturation:
            # source is too bright, skip it
            continue

        pars = (bg, peak, x_ctr, y_ctr, x_width, y_width, theta)
    
        # instantiate fit object
        gf = GaussFit(stamp_data)

        # least squares fit of model to data
        try:
            # for scipy versions < 0.9
            new_pars, success = scipy.optimize.leastsq(gf.calc_diff, pars,
                                                       maxfev=1000, 
                                                       warning=False)
        except:
            # for scipy versions >= 0.9
            import warnings
            warnings.simplefilter('ignore')
            new_pars, success = scipy.optimize.leastsq(gf.calc_diff, pars,
                                                       maxfev=1000)


        if success>=4:
            # fit failed, move on
            continue

        (bg, peak, x_ctr, y_ctr, x_width, y_width, theta) = new_pars

        # convert fit parameters to FWHM, ellipticity
        fwhmx = abs(2*np.sqrt(2*np.log(2))*x_width)
        fwhmy = abs(2*np.sqrt(2*np.log(2))*y_width)
        pa = (theta*(180/np.pi))
        pa = pa%360
                
        if fwhmy < fwhmx:
            ellip = 1 - fwhmy/fwhmx
            fwhm = fwhmx
        elif fwhmx < fwhmy:
            ellip = 1 - fwhmx/fwhmy                    
            pa = pa-90 
            fwhm = fwhmy
        else: #fwhmx == fwhmy
            ellip = 0
            fwhm = fwhmx

        if ellip>.25:
            # source not round enough, skip it
            continue

        if fwhm>3*default_fwhm: # ie. 2.4 arcsec -- probably not due to seeing
            # source not pointy enough, skip it
            continue

        # update the position from the fit center
        newx = xlow + x_ctr + 1
        newy = ylow + y_ctr + 1

        good_source.append([newx,newy])

    return good_source


def _correlate_sources(ad1, ad2, delta=None, firstPass=10, cull_sources=False):
    """
    This function takes sources from the OBJCAT extensions in two
    images and attempts to correlate them.  It returns a list of 
    reference source positions and their correlated image source 
    positions.

    :param ad1: reference image
    :type ad1: AstroData instance

    :param ad2: input image
    :type ad2: AstroData instance

    :param delta: maximum distance in pixels to allow a match. If
                  left as None, it will attempt to find an appropriate
                  number (recommended).
    :type delta: float
    
    :param firstPass: estimated maximum distance between correlated
                      sources.  This distance represents the expected
                      mismatch between the WCSs of the input images.
    :type firstPass: float

    :param cull_sources: flag to indicate whether to reject sources that
                   are insufficiently star-like.  If true, will fit
                   a Gaussian to each correlated source, and return
                   the fit center of good sources (rather than the raw
                   OBJCAT position).
    :type cull_sources: bool

    """

    log = gemLog.getGeminiLog()

    # get data and WCS from image 1
    x1 = ad1['OBJCAT'].data.field('x')
    y1 = ad1['OBJCAT'].data.field('y')
    wcs1 = pywcs.WCS(ad1['SCI'].header)

    # get data and WCS from image 2    
    x2 = ad2['OBJCAT'].data.field('x')
    y2 = ad2['OBJCAT'].data.field('y')
    wcs2 = pywcs.WCS(ad2['SCI'].header)

    # convert image 2 data to sky coordinates
    ra2, dec2 = wcs2.wcs_pix2sky(x2,y2,1)

    # convert image 2 sky data to image 1 pixel coordinates
    conv_x2, conv_y2 = wcs1.wcs_sky2pix(ra2,dec2,1)
    

    # find matches
    ind1,ind2 = match_cxy(x1,conv_x2,y1,conv_y2,
                          delta=delta, firstPass=firstPass, log=log)
    
    if len(ind1)!=len(ind2):
        raise Errors.ScienceError('Mismatched arrays returned from match_cxy')

    if len(ind1)<1 or len(ind2)<1:
        return [[],[]]
    else:
        obj_list = [zip(x1[ind1], y1[ind1]),
                    zip(x2[ind2], y2[ind2])]

        if cull_sources:
            log.status('Rejecting non-Gaussian sources')
            obj_list_1 = np.array(_cull_sources(ad1, obj_list[0]))
            obj_list_2 = np.array(_cull_sources(ad2, obj_list[1]))

            x1,y1 = obj_list_1[:,0],obj_list_1[:,1]
            x2,y2 = obj_list_2[:,0],obj_list_2[:,1]

            # re-match sources
            ra2, dec2 = wcs2.wcs_pix2sky(x2,y2,1)
            conv_x2, conv_y2 = wcs1.wcs_sky2pix(ra2,dec2,1)
            ind1,ind2 = match_cxy(x1,conv_x2,y1,conv_y2,
                                  delta=delta, firstPass=firstPass, log=log)

            if len(ind1)!=len(ind2):
                raise Errors.ScienceError('Mismatched arrays returned ' +
                                          'from match_cxy')

            if len(ind1)<1 or len(ind2)<1:
                return [[],[]]
            else:
                obj_list = [zip(x1[ind1], y1[ind1]),
                            zip(x2[ind2], y2[ind2])]

        return obj_list




def _align_wcs(reference, adinput, objIns, rotate=False, scale=False):
    """
    This function fits an input image's WCS to a reference image's WCS
    by minimizing the difference in the input image frame between
    reference points present in both images.

    :param reference: reference image to register other images to. Must
                      have only one SCI extension.
    :type reference: AstroData object
    
    :param adinput: images to register to reference image.  Must have
                  only one SCI extension.
    :type adinput: AstroData objects, either a single instance or a list

    :param objIns: list of object lists, one for each input image
    :type objIns: list of output lists from _correlate_sources

    :param rotate: flag to indicate whether the input image WCSs should
                   be allowed to rotate with respect to the reference image
                   WCS
    :type rotate: bool

    :param scale: flag to indicate whether the input image WCSs should
                  be allowed to scale with respect to the reference image
                  WCS.  The same scale factor is applied to all dimensions.
    :type scale: bool

    """

    log = gemLog.getGeminiLog()

    if not isinstance(adinput,list):
        adinput = [adinput]
    if not isinstance(objIns,list):
        objIns = [objIns]
    if len(objIns) != len(adinput):
        raise Errors.InputError('Argument objIns should have the same ' +
                                'number of elements as adinput')

    import scipy.optimize

    adoutput_list = []
    for i in range(len(adinput)):

        # copy input ad and rename
        ad = adinput[i]

        log.status('Starting WCS adjustment for ' + ad.filename )

        ref_xy, inp_xy = objIns[i]
        ref_xy = np.array(ref_xy)
        inp_xy = np.array(inp_xy)

        ref_wcs = pywcs.WCS(reference['SCI'].header)
        inp_wcs = pywcs.WCS(ad['SCI'].header)

        # convert the reference coordinates to RA/Dec
        ref_radec = ref_wcs.wcs_pix2sky(ref_xy,1)

        # instantiate the alignment object used to fit input
        # WCS to reference WCS
        wcstweak = WCSTweak(inp_wcs, inp_xy, ref_radec, 
                            rotate=rotate, scale=scale)        
        
        # find optimum WCS shift and rotation with
        # starting parameters: dRA, dDec = 0
        # (and dTheta=0 if rotate=True, dMag=1 if scale=True)

        update = False
        if rotate and scale:
            pars = [0,0,0,1]
        elif rotate:
            pars = [0,0,0]
        elif scale:
            pars = [0,0,1]
        else:
            pars = [0,0]

        try:
            # for scipy versions < 0.9
            new_pars,success = scipy.optimize.leastsq(wcstweak.calc_diff, pars,
                                                      warning=False, 
                                                      maxfev=1000)
        except:
            # for scipy versions >= 0.9
            import warnings
            warnings.simplefilter('ignore')
            new_pars,success = scipy.optimize.leastsq(wcstweak.calc_diff, pars,
                                                      maxfev=1000)

        if success<4:
            update = True
            if rotate and scale:
                (dRA, dDec, dTheta, dMag) = new_pars
                log.fullinfo('Best fit dRA, dDec, dTheta, dMag: ' +
                             '%.5f %.5f %.5f %.5f' %
                             (dRA, dDec, dTheta, dMag))
            elif rotate:
                (dRA, dDec, dTheta) = new_pars
                log.fullinfo('Best fit dRA, dDec, dTheta: %.5f %.5f %.5f' %
                             (dRA, dDec, dTheta))
            elif scale:
                (dRA, dDec, dMag) = new_pars
                log.fullinfo('Best fit dRA, dDec, dMag: %.5f %.5f %.5f' %
                             (dRA, dDec, dMag))
            else:
                (dRA, dDec) = new_pars
                log.fullinfo('Best fit dRA, dDec: %.5f %.5f' % (dRA, dDec))
        else:
            log.warning('WCS alignment did not converge. Not updating WCS.')

        # update WCS in ad
        if update:
            log.status('Updating WCS in header')
            ad.phu_set_key_value('CRVAL1', wcstweak.wcs.wcs.crval[0])
            ad.phu_set_key_value('CRVAL2', wcstweak.wcs.wcs.crval[1])
            ad.phu_set_key_value('CD1_1', wcstweak.wcs.wcs.cd[0,0])
            ad.phu_set_key_value('CD1_2', wcstweak.wcs.wcs.cd[0,1])
            ad.phu_set_key_value('CD2_1', wcstweak.wcs.wcs.cd[1,0])
            ad.phu_set_key_value('CD2_2', wcstweak.wcs.wcs.cd[1,1])
            
            for ext in ad:
                if ext.extname() in ['SCI','VAR','DQ']:
                    ext.set_key_value('CRVAL1', wcstweak.wcs.wcs.crval[0])
                    ext.set_key_value('CRVAL2', wcstweak.wcs.wcs.crval[1])
                    ext.set_key_value('CD1_1', wcstweak.wcs.wcs.cd[0,0])
                    ext.set_key_value('CD1_2', wcstweak.wcs.wcs.cd[0,1])
                    ext.set_key_value('CD2_1', wcstweak.wcs.wcs.cd[1,0])
                    ext.set_key_value('CD2_2', wcstweak.wcs.wcs.cd[1,1])


        adoutput_list.append(ad)

    return adoutput_list


def _header_align(reference, adinput):
    """
    This function uses the POFFSET and QOFFSET header keywords
    to get reference points to use in correcting an input WCS to
    a reference WCS.  Positive POFFSET is assumed to mean higher x
    value, and positive QOFFSET is assumed to mean higher y value.
    This function only allows for relative shifts between the images;
    rotations and scales will not be handled properly

    :param reference: reference image to register other images to. Must
                      have only one SCI extension.
    :type reference: AstroData object
    
    :param adinput: images to register to reference image.  Must have
                  only one SCI extension.
    :type adinput: AstroData objects, either a single instance or a list

    """

    log = gemLog.getGeminiLog()

    if not isinstance(adinput,list):
        adinput = [adinput]
 
    # get starting offsets from reference image (first one given)
    pixscale = float(reference.pixel_scale())
    ref_xoff = reference.phu_get_key_value('POFFSET')/pixscale
    ref_yoff = reference.phu_get_key_value('QOFFSET')/pixscale

    # reference position is the center of the reference frame
    data_shape = reference['SCI'].data.shape
    ref_coord = [data_shape[1]/2,data_shape[0]/2]

    log.fullinfo('Pixel scale: %.4f' % pixscale)
    log.fullinfo('Reference offsets: %.4f %.4f' % (ref_xoff, ref_yoff))
    log.fullinfo('Reference coordinates: %.1f %.1f' % 
                 (ref_coord[0], ref_coord[1]))

    objIns = []
    for i in range(len(adinput)):
        ad = adinput[i]
        pixscale = float(ad.pixel_scale())
        xoff = ad.phu_get_key_value('POFFSET')/pixscale
        yoff = ad.phu_get_key_value('QOFFSET')/pixscale

        img_x = xoff-ref_xoff + ref_coord[0]
        img_y = yoff-ref_yoff + ref_coord[1]

        log.fullinfo('For image ' + ad.filename + ':')
        log.fullinfo('   Image offsets: %.4f %.4f' % (xoff, yoff))
        log.fullinfo('   Coordinates to transform: %.4f %.4f' % (img_x, img_y))

        objIns.append(np.array([[ref_coord],[[img_x,img_y]]]))

    adoutput_list = _align_wcs(reference, adinput, objIns, 
                               rotate=False, scale=False)

    return adoutput_list


def _user_align(reference, adinput, rotate, scale):
    """
    This function takes user input to get reference points to use in 
    correcting an input WCS to a reference WCS.  The images are 
    displayed by IRAF and the common points selected by an image cursor.
    If rotation or scaling degrees of freedom are desired, two 
    common points must be selected.  If only shifts are desired, one
    common point must be selected.

    :param reference: reference image to register other images to. Must
                      have only one SCI extension.
    :type reference: AstroData object
    
    :param adinput: images to register to reference image.  Must have
                  only one SCI extension.
    :type adinput: AstroData objects, either a single instance or a list

    :param rotate: flag to indicate whether the input image WCSs should
                   be allowed to rotate with respect to the reference image
                   WCS
    :type rotate: bool

    :param scale: flag to indicate whether the input image WCSs should
                  be allowed to scale with respect to the reference image
                  WCS.  The same scale factor is applied to all dimensions.
    :type scale: bool

    """

    log = gemLog.getGeminiLog()

    # load pyraf modules
    from astrodata.adutils.gemutil import pyrafLoader
    pyraf, gemini, yes, no = pyrafLoader()

    if not isinstance(adinput,list):
        adinput = [adinput]
 
    # start cl manager for iraf display
    all_input = [reference] + adinput
    clm = man.CLManager(imageIns=all_input, funcName='display', log=log)
    tmpfiles = clm.imageInsFiles(type='list')

    # display the reference image
    print " ==> Reference image: " + reference.filename
    pyraf.iraf.display(tmpfiles[0]+'[SCI]', 1)

    if not rotate and not scale:
        # only one object needed for pure shifts
        print "Point to one common object in reference image"
        print "    strike any key"
        words = pyraf.iraf.cl.imcur.split()
        x11 = float(words[0])
        y11 = float(words[1])
        ref_coord = [[x11,y11]]
    else:
        # select two objects for rotation/scaling
        print "Point to first common object in reference image"
        print "    strike any key"
        words = pyraf.iraf.cl.imcur.split()
        x11 = float(words[0])
        y11 = float(words[1])
        print "Point to second common object in reference image"
        print "    strike any key"
        words = pyraf.iraf.cl.imcur.split()
        x12 = float(words[0])
        y12 = float(words[1])
        ref_coord = [[x11,y11],[x12,y12]]


    objIns = []
    for i in range(len(adinput)):
        ad = adinput[i]

        print " ==> Image to be transformed:", ad.filename
        pyraf.iraf.display(tmpfiles[i+1]+'[SCI]', 1)

        if not rotate and not scale:
            print "Point to one common object in image to be transformed"
            print "    coordinates for last image: %.1f, %.1f" % (x11, y11)
            print "    strike any key"
            words = pyraf.iraf.cl.imcur.split()
            x21 = float(words[0])
            y21 = float(words[1])
            img_coord = [[x21, y21]]
        else:
            print "Point to first common object in image to be transformed"
            print "    coordinates for last image: %.1f, %.1f" % (x11, y11)
            print "    strike any key"
            words = pyraf.iraf.cl.imcur.split()
            x21 = float(words[0])
            y21 = float(words[1])

            print "Point to second common object in image to be transformed"
            print "    coordinates for last image: %.1f, %.1f" % (x12, y12)
            print "    strike any key"
            words = pyraf.iraf.cl.imcur.split()
            x22 = float(words[0])
            y22 = float(words[1])
            img_coord = [[x21, y21],[x22,y22]]

        log.fullinfo('Reference coordinates: '+repr(ref_coord))
        log.fullinfo('Coordinates to transform: '+repr(img_coord))
        objIns.append([ref_coord,img_coord])

    # delete temporary files
    clm.finishCL()


    adoutput_list = _align_wcs(reference, adinput, objIns, 
                               rotate=rotate, scale=scale)

    return adoutput_list


# below functions will go into a toolbox

class GaussFit:
    """
    This class provides access to a Gaussian model, intended
    to be fit to a small stamp of data, via a  minimization of the 
    differences between the model and the data.

    Example usage:
    pars = (bg, peak, x_ctr, y_ctr, x_width, y_width, theta)
    gf = GaussFit(stamp_data)
    new_pars, success = scipy.optimize.leastsq(gf.calcDiff, pars, maxfev=1000)
    """


    def __init__(self, stamp_data):
        """
        This instantiates the fitting object.
        
        :param stamp_data: array containing image data, preferably the
                           source to fit plus a little padding
        :type stamp_data: NumPy array
        """
        self.stamp = stamp_data

    def model_gauss2d(self, pars):
        """
        This function returns a Gaussian source in an image array the
        same shape as the stamp_data.  The Gaussian is determined by
        the parameters in pars.
        
        :param pars: Gaussian parameters in this order: background, peak,
                     x-center, y-center, x-width, y-width, position angle
                     (in degrees)
        :type pars: 7-element tuple
        """
        bg, peak, cx, cy, wx, wy, theta = pars
        
        model_fn = lambda y,x: bg + peak*np.exp(-( (((x-cx)*np.cos(theta)
                                                 + (y-cy)*np.sin(theta))
                                                /wx)**2
                                               +(((x-cx)*np.sin(theta)
                                                  -(y-cy)*np.cos(theta))
                                                 /wy)**2)
                                              /2)
        gauss_array = np.fromfunction(model_fn, self.stamp.shape)
        return gauss_array

    def calc_diff(self, pars):
        """
        This function returns an array of the differences between
        the model stamp and the data stamp.  It is intended to be fed
        to an optimization algorithm, such as scipy.optimize.leastsq.

        :param pars: Gaussian parameters in this order: background, peak,
                     x-center, y-center, x-width, y-width, position angle
                     (in degrees)
        :type pars: 7-element tuple
        """
        model = self.model_gauss2d(pars).flatten()
        diff = self.stamp.flatten() - model
        return diff


def get_corners(shape):
    """
    This is a recursive function to calculate the corner indices 
    of an array of the specified shape.

    :param shape: length of the dimensions of the array
    :type shape: tuple of ints, one for each dimension

    """
    if not type(shape)==tuple:
        raise Errors.TypeError('get_corners argument is non-tuple')

    if len(shape)==1:
        corners = [(0,), (shape[0]-1,)]
    else:
        shape_less1 = shape[1:len(shape)]
        corners_less1 = get_corners(shape_less1)
        corners = []
        for corner in corners_less1:
            newcorner = (0,) + corner
            corners.append(newcorner)
            newcorner = (shape[0]-1,) + corner
            corners.append(newcorner)
        
    return corners

def rotate(degs):
    """
    Little helper function to return a basic 2-D rotation matrix.

    :param degs: rotation amount, in degrees
    :type degs: float
    """
    rads = np.radians(degs)
    s = np.sin(rads)
    c = np.cos(rads)
    return np.array([[c,-s],
                     [s,c]])
class WCSTweak:
    """
    This class allows slight tweaking of an image's WCS, to fit to a
    reference WCS, via a minimization of the differences between
    reference points in both images.

    Example usage:
    wcstweak = WCSTweak(inp_wcs, inp_xy, ref_radec)
    pars = [0,0]
    new_pars,success = scipy.optimize.leastsq(wcstweak.calc_diff, pars,
                                              maxfev=1000)
    """

    def __init__(self, wcs, inp, ref, rotate=False, scale=False):
        """
        This instantiates the WCSTweak object.
        
        :param wcs: the input image WCS
        :type wcs: pywcs WCS object

        :param inp: input object position in input pixel frame
        :type inp: NumPy array of [x,y] positions

        :param ref: reference object positions in sky frame (RA/Dec)
        :type ref: NumPy array of [ra,dec] positions

        :param rotate: flag to indicate whether to allow rotation of
                       input WCS with respect to reference WCS
        :type rotate: bool

        :param scale: flag to indicate whether to allow scaling of
                      input WCS with respect to reference WCS
        :type scale: bool
        """
        self.wcs = wcs
        self.inp = inp.flatten() # in input pixel frame
        self.ref = ref           # in ra/dec
        self.rotate = rotate
        self.scale = scale
        self.crval = wcs.wcs.crval.copy()
        self.cd = wcs.wcs.cd.copy()

    def transform_ref(self, pars):
        """
        This function transforms reference RA/Dec into input pixel
        frame, via a WCS tweaked by parameters pars.
        
        :param pars: list of parameters to tweak WCS by. Number of 
                     elements is determined by whether rotation/scaling
                     is allowed.  Order is [dRA, dDec, dTheta, dMag].
                     dTheta and dMag are optional.
        :type pars: list of 2, 3, or 4 elements.
        """
        if self.rotate and self.scale:
            d_ra, d_dec, d_theta, d_mag = pars
            self.wcs.wcs.cd = np.dot(d_mag*rotate(d_theta),self.cd)
        elif self.rotate:
            d_ra, d_dec, d_theta = pars
            self.wcs.wcs.cd = np.dot(rotate(d_theta),self.cd)
        elif self.scale:
            d_ra, d_dec, d_mag = pars
            self.wcs.wcs.cd = d_mag*self.cd
        else:
            d_ra, d_dec = pars

        self.wcs.wcs.crval = self.crval + np.array([d_ra, d_dec])/3600.0

        new_ref = self.wcs.wcs_sky2pix(self.ref,1)
        return new_ref.flatten()

    # calculate residual (called by scipy.optimize.leastsq)
    def calc_diff(self, pars):
        """
        This function returns an array of the differences between the
        input sources and the reference sources in the input pixel frame.
        It is intended to be fed to an optimization algorithm, such as 
        scipy.optimize.leastsq.

        :param pars: list of parameters to tweak WCS by. Number of 
                     elements is determined by whether rotation/scaling
                     is allowed.  Order is [dRA, dDec, dTheta, dMag].
                     dTheta and dMag are optional.
        :type pars: list of 2, 3, or 4 elements.
        """
        new_ref = self.transform_ref(pars)
        diff = self.inp - new_ref
        return diff


def match_cxy (xx, sx, yy, sy, firstPass=50, delta=None, log=None):
    """
    Match reference positions (sx,sy) with those of the 
    object catalog (xx,yy). 
    Select those that are within delta pixels from
    the object positions.
    
    firstPass:  (50) First pass delta radius.

    This matching is a 2 pass algorithm. First pass takes the
    larger delta and adjust the reference positions by the median 
    of the x,y offset. The second pass takes 'delta' value and 
    look for those x,y now closer to the reference positions.
    The units are pixels for the deltas. If you are passing
    degrees, the deltas need to be consistent.


    OUTPUT:
    - obj_index: Index array of the objects matched.
    - ref_index: Index array of the referecences matched.
    """

    # turn to numpy arrays
    xx, sx, yy, sy = map(np.asarray,(xx,sx,yy,sy))

    def getg(xx, sx, yy, sy, deltax=2.5, deltay=2.5):
        """ Return object(xx) and reference(sx) indices of
        common positions.
        OUTPUT
        g:    Indices of the object position common to
        r:    indices of the reference position
        """
        dax=[]; day=[]; g=[]; r=[]
        for k in range(len(sx)):
            gindx,= np.where((abs(xx-sx[k])<deltax) & 
                             (abs(yy-sy[k])<deltay))
            for i in gindx:
                dx = xx[i] - sx[k] 
                dy = yy[i] - sy[k] 

                # if there are multiple matches, keep only the
                # closest one
                if i in g or k in r:
                    if i in g:
                        first_ind = g.index(i)
                    else:
                        first_ind = r.index(k)
                    first_dist = dax[first_ind]**2 + day[first_ind]**2
                    this_dist = dx**2 + dy**2
                    if (first_dist > this_dist):
                        del dax[first_ind]
                        del day[first_ind]
                        del g[first_ind]
                        del r[first_ind]
                        dax.append(dx)
                        day.append(dy)
                        g.append(i)
                        r.append(k)
                else:
                    dax.append(dx)
                    day.append(dy)
                    g.append(i)
                    r.append(k)
                            
        dax,day = map(np.asarray, (dax,day))
        mx = np.median(dax); stdx = np.std(dax)
        my = np.median(day); stdy = np.std(day)
            
        return np.asarray(g),np.asarray(r),mx,my,stdx,stdy 


    # Select only those standards with less than 10 pixels from objects.
    # Get the median values (mx,my) of the differences and add these
    # to the standard positions.

    #NOTE: We are setting a large delta here, we would need to see
    #      median 1st...

    ig,r,mx,my,stdx,stdy = getg(xx,sx,yy,sy, deltax=firstPass,deltay=firstPass)
    log.info('Median differences (x,y):%.2f %.2f, %.2f %.2f' % 
             (mx,my,stdx,stdy)+"[First iteration]")

    if len(r) == 0:
        return ig,r
        
    # Now shift reference position by adding the median of the
    # differences. The standards are now closer to the object positions.
    sx = sx + mx   
    sy = sy + my

    # Select only those that are closer than delta or default(6.5) pixels.
    xxx = xx[ig]; yyy = yy[ig] 
    deltax=delta; deltay=delta
    if delta == None:
        deltax=2*stdx; deltay=2*stdy
    g,r,mx,my,stdx,stdy = getg (xxx, sx, yyy, sy, 
                                deltax=deltax, deltay=deltay)
    log.info('Median differences (x,y):%.2f %.2f %.2f %.2f' %
             (mx,my,stdx,stdy)+"[Second iteration]")

    if g.size == 0:
        indxy,indr=[],[]
    else:
        indxy = ig[g]
        indr = r
        
    return indxy, indr
