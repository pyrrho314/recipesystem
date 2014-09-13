#! /usr/bin/env python

import sys
import math

import numpy as N
import pyfits

import gmos_ifu_data

# examples:
#  imrecon.py N20060131S0019.fits obj.fits object
#  imrecon.py N20060131S0019.fits bkg.fits background

def main (args):

    if len (args) != 3:
        print "syntax:  imrecon.py input output region"
        print "region may be 'object' or 'background'"
        sys.exit()

    if args[2] != 'object' and args[2] != 'background':
        raise RuntimeError, "region '%s' is not supported" % args[2]

    # xxx 5 and 11 are for N20060131S0019.fits xxx
    xref = gmos_ifu_data.XPOS_N + 5
    yref = gmos_ifu_data.YPOS_N + 11

    gmosIfuImage (args[0], args[1], xref, yref, args[2])

def gmosIfuImage (input, output, xref, yref, region_name):
    """Reconstruct an IFU image from a target-acquisition image.

    @param input: name of input FITS file containing target-acq image
    @type input: string
    @param output: name of FITS file to be created
    @type output: string
    @param xref: X pixel location of column of fiber images
    @type xref: float
    @param yref: Y pixel location of the reference fiber
    @type yref: float
    @param region_name: 'object' or 'background'
    @type region_name: string
    """

    fd = pyfits.open (input)

    if region_name == 'object':
        region_list = [gmos_ifu_data.object1, gmos_ifu_data.object2]
        shape = (25, 34)
    elif region_name == 'background':
        region_list = [gmos_ifu_data.background1, gmos_ifu_data.background2]
        shape = (25, 16)

    maskname = getKeyword (fd[0].header, "maskname")
    (nifuslits, sslit) = interpretMaskname (maskname)

    n_all_slices = 0
    flux_list = []
    ext_list = [1, 3]
    for i in range (nifuslits):

        # k indicates the side of the field (which slit)
        # xxx I don't think the logic here is correct; may need to
        # do this differently depending on whether the data are from
        # GMOS-N or GMOS-S. xxx
        if sslit == "red":
            k = 0
            xref_k = xref
        elif sslit == "blue":
            k = 1
            xref_k = xref + gmos_ifu_data.DELTA_XPOS
        elif i == 0:
            k = 0
            xref_k = xref
        else:
            k = 1
            xref_k = xref + gmos_ifu_data.DELTA_XPOS

        image = getImage (fd, ext_list[k])
        region = region_list[k]
        n_slices = len (region)
        flux_list.append (fiberFlux (image, xref_k, yref, region))
        n_all_slices += n_slices

    flux_list.reverse()
    flux = concatFlux (flux_list)

    (hex_x, hex_y) = hexGridXY (n_all_slices)
    (crpix1, crpix2) = adjustZeroPt (hex_x, hex_y, shape)
    ### debugIFU (flux, hex_x, hex_y)

    (nearest, weights) = computeWeights (hex_x, hex_y, shape)
    ofd = fd[0:2]
    ofd[1].data = interpFlux (flux, nearest, weights, shape)

    fd[1].header["crpix1"] = crpix1
    fd[1].header["crpix2"] = crpix2

    ofd.writeto (output)
    ofd.close()
    fd.close()

def debugIFU (flux, hex_x, hex_y):

    fd = open ("debug.txt", "a")
    n = len (flux)
    for k in range (n):
        fd.write ("%10.2f  %8.3f  %8.3f\n" % (flux[k], hex_x[k], hex_y[k]))
    fd.close()

def getImage (fd, ext):
    """Copy out the specified image extension, and subtract bias level.

    @param fd: header/data list for input file
    @type fd: pyfits HDUList object
    @param ext: extension number (we expect ext to be either 1 or 3)
    @type ext: int

    @return: image data with bias subtracted
    @rtype: array
    """

    hdr = fd[ext].header
    data = fd[ext].data

    biassec = getKeyword (hdr, "biassec")
    biassec = biassec.strip()
    if biassec.endswith ("]"):
        biassec = biassec[:-1]
    if biassec.startswith ("["):
        biassec = biassec[1:]
    words = biassec.split (",")
    wx = words[0].split (":")
    wy = words[1].split (":")
    xlow  = int (wx[0]) - 1
    xhigh = int (wx[1]) - 1
    ylow  = int (wy[0]) - 1
    yhigh = int (wy[1]) - 1

    s = "N.median (data[%d:%d,%d:%d].ravel())" % (ylow, yhigh, xlow, xhigh)
    bias_level = eval (s)

    data -= bias_level

    return data

def getKeyword (header, keyword):

    value = header.get (keyword, "N/A")

    return value

def interpretMaskname (maskname):
    """Determine which IFU was used.

    @param maskname: 
    @type maskname: string

    @return: number of slits used (nifuslits) and an indication of which side
        was used (sslit, "red" or "blue", or "" if both slits)
    @rtype: tuple
    """

    sslit = ""
    if maskname == "IFU-1":
        nifuslits = 1
        sslit = "red"
    elif maskname == "IFU-R":
        nifuslits = 1
        sslit = "red"
    elif maskname == "IFU-B":
        nifuslits = 1
        sslit = "blue"
    elif maskname == "IFU-2":
        nifuslits = 2
    elif maskname == "IFU-NS-B":
        nifuslits = 1
        sslit = "blue"
    elif maskname == "IFU-NS-R":
        nifuslits = 1
        sslit = "red"
    elif maskname == "IFU-NS-2":
        nifuslits = 2
    else:
        raise RuntimeError, "maskname '%s' is not valid" % maskname

    return (nifuslits, sslit)

def fiberFlux (image, xref, yref, region, datatype=N.float32):
    """Find the flux (sum of counts) for each fiber.

    @param image: image data (bias subtracted)
    @type image: array
    @param xref: X pixel location of column of fiber images
    @type xref: float
    @param yref: Y pixel location of the reference fiber
    @type yref: float
    @param region: list of tuples, each tuple giving a range of Y pixel
        numbers for a group of fibers
    @type region: list
    @param datatype: data type to use for the output array
    @type datatype: numpy dtype object

    @return: 1-D array with the flux for each fiber
    @rtype: array
    """

    ixref = int (round (xref))
    xlow = ixref - 8
    xhigh = ixref + 10
    dy = gmos_ifu_data.FIBER_Y

    flux = N.zeros (len (region) * gmos_ifu_data.N_FIBERS, dtype=datatype)

    k = 0
    for (y0, y1) in region:
        if y0 < y1:
            y = yref + y0 + dy/2.
            for i in range (gmos_ifu_data.N_FIBERS):
                ylow = int (round (y - 2.5))
                yhigh = int (round (y + dy - 2.5))
                flux[k] = N.sum (image[ylow:yhigh,xlow:xhigh], dtype=N.float64)
                y += dy
                k += 1
        else:
            y = yref + y0 - dy/2.
            for i in range (gmos_ifu_data.N_FIBERS):
                ylow = int (round (y - dy + 2.5))
                yhigh = int (round (y + 2.5))
                flux[k] = N.sum (image[ylow:yhigh,xlow:xhigh], dtype=N.float64)
                y -= dy
                k += 1

    return flux

def concatFlux (flux_list):
    """Concatenate the flux arrays in flux_list.

    @param flux_list: list of flux arrays
    @type flux_list: list

    @return: one array with concatenated input arrays
    @rtype: array

    It is assumed that the elements of flux_list are already in the correct
    order.  If flux_list includes more than one (i.e. two) elements, their
    original order will be the right-hand side followed by the left-hand side,
    the order of the extensions (1 and 3 respectively) in the input image.
    So the calling function should swap the elements of flux_list prior to
    calling this function.
    """

    if len (flux_list) > 1:
        full_length = 0
        for flux in flux_list:
            full_length += len (flux)
        concatenated_flux = N.zeros (full_length, dtype=flux.dtype)
        i = 0
        for flux in flux_list:
            n = len (flux)
            concatenated_flux[i:i+n] = flux
            i += n
    else:
        concatenated_flux = flux_list[0]

    return concatenated_flux


def hexGridXY (n_slices, w=1., datatype=N.float32):
    """Get (x,y) coordinates in the input field of view for each fiber.

    @param n_slices: number of slices over the input image (20 or 40)
    @type n_slices: int
    @param w: width of a hexagonal lenslet (fiber spacing), in pixels
        (width of a hexagon in the smallest dimension)
    @type w: float
    @param datatype: data type for arrays of fiber X and Y coordinates
    @type datatype: numpy dtype

    @return: arrays of X and Y coordinates of fibers
    @rtype: tuple of arrays

    Lensets are packed adjacent to each other vertically and along a
    30-degree diagonal.
    The spacing in the X direction (columns) between slices is w * sqrt (3).
    The spacing in the Y direction (rows) between fibers is w/2.
    """

    sqrt3 = math.sqrt (3.)

    hex_x = N.zeros (n_slices * gmos_ifu_data.N_FIBERS, dtype=datatype)
    hex_y = N.zeros (n_slices * gmos_ifu_data.N_FIBERS, dtype=datatype)

    n_groups = n_slices // 2
    k = 0
    for i in range (n_groups):
        x = i * w * sqrt3
        for j in range (gmos_ifu_data.N_FIBERS):
            y = j * w
            hex_x[k] = x
            hex_y[k] = y
            k += 1
        x += (sqrt3 / 2.)       # next column (second half of a slice)
        for j in range (gmos_ifu_data.N_FIBERS):
            y = j * w + w/2.
            hex_x[k] = x
            hex_y[k] = y
            k += 1

    return (hex_x, hex_y)

def adjustZeroPt (hex_x, hex_y, shape):
    """Shift the coordinates to put the center on a pixel.

    @param hex_x: X positions of fibers; modified in-place
    @type hex_x: 1-D array
    @param hex_y: Y positions of fibers; modified in-place
    @type hex_y: 1-D array
    @param shape: axis lengths of image array to be created
    @type shape: tuple

    @return: reference pixel location (one indexed)
    @rtype: tuple of floats

    When the X and Y coordinates of the fibers are created by hexGridXY,
    the zero point is an arbitrary point.  This function adds offsets to
    the X and Y coordinates so the mean of the X coordinates and the mean
    of the Y coordinates are centered on a pixel.  The pixel numbers in X
    and Y will be returned; these are the FITS keywords CRPIX1 and CRPIX2
    respectively.
    """

    mean_x = hex_x.mean()
    mean_y = hex_y.mean()

    center0 = float (shape[0] // 2)
    center1 = float (shape[1] // 2)

    hex_x += (center1 - mean_x)
    hex_y += (center0 - mean_y)

    return (center1+1., center0+1.)

def computeWeights (hex_x, hex_y, shape):
    """Interpolate the weights at each point in an image.

    @param hex_x: 1-D array of X positions of fibers in field of view
    @type hex_x: array
    @param hex_y: 1-D array of Y positions of fibers in field of view
    @type hex_y: array
    @param shape: axis lengths of image array to be created
    @type shape: tuple

    @return: a tuple containing 'nearest' and 'weights', each of which is
        a list of three-element tuples, with the length of each list equal
        to the number of pixels in the output array which will be created.
        A tuple in 'nearest' is the indices of the three nearest fibers,
        and a tuple in 'weights' is the weights to use for bilinear
        interpolation.
    @rtype: tuple

    The unit of distance for the X and Y positions is the spacing between
    fibers in the Y direction or along a 30-degree diagonal (the closest
    spacing).  This is also the pixel spacing in the output image.
    """

    ny = shape[0]
    nx = shape[1]

    # index numbers and weights for the three nearest fibers, for each pixel
    nearest = []
    weights = []

    n_colinear = 0
    for j in range (ny):
        for i in range (nx):
            klist = []
            for k in range (len (hex_x)):
                dist2 = (hex_x.item(k) - i)**2 + (hex_y.item(k) - j)**2
                # xxx dist2 = (hex_x.item(k) - i)**2 + (hex_y.item(k) - j)**2
                klist.append ((dist2, k))
            klist.sort()
            # index numbers for the three fibers nearest the current pixel
            # (but start with four fibers and select three non-colinear)
            colinear = False
            try:
                (k0, k1, k2) = pickThree (klist[0][1], klist[1][1],
                                      klist[2][1], klist[3][1], hex_x, hex_y)
            except RuntimeError:
                colinear = True
                n_colinear += 1
                (k0, k1, k2) = (0, 0, 0)
                continue
            dx = []
            dy = []
            for k in (k0, k1, k2):
                dx.append (hex_x.item(k) - i)
                dy.append (hex_y.item(k) - j)

            # compute the weights
            if colinear:
                w0 = 0.
                w1 = 0.
                w2 = 0.
            else:
                D = (dx[0] - dx[1]) * (dy[0] - dy[2]) - \
                    (dx[0] - dx[2]) * (dy[0] - dy[1])

                w0 = 1. + \
                    (dx[0] * (dy[2] - dy[1]) - dy[0] * (dx[2] - dx[1])) / D
                w1 = (dx[2] * dy[0] - dx[0] * dy[2]) / D
                w2 = (dx[0] * dy[1] - dx[1] * dy[0]) / D

            nearest.append ((k0, k1, k2))
            weights.append ((w0, w1, w2))

    return (nearest, weights)

def interpFlux (flux, nearest, weights, shape, datatype=N.float32):
    """Interpolate the flux at each point in an image.

    @param flux: 1-D array of fluxes, one for each fiber
    @type flux: array
    @param nearest: list of indices of nearest pixels
    @type nearest: list
    @param weights: list of weights
    @type weights: list
    @param shape: axis lengths of image array to be created
    @type shape: tuple

    @return: 2-D image of the field of view
    @rtype: array

    The unit of distance for the X and Y positions is the spacing between
    fibers in the Y direction or along a 30-degree diagonal (the closest
    spacing).  This is also the pixel spacing in the output image.

    The 'flux' array has one element for each fiber.  In contrast, the
    'nearest' and 'weights' lists have one element for each pixel in the
    output image array.
    """

    outimage = N.zeros (shape, dtype=datatype)
    ny = shape[0]
    nx = shape[1]

    near = N.array (nearest, dtype=N.int32)
    wgt = N.array (weights, dtype=N.float64)

    k0 = near[:,0]
    k1 = near[:,1]
    k2 = near[:,2]
    w0 = wgt[:,0]
    w1 = wgt[:,1]
    w2 = wgt[:,2]

    outimage.ravel()[:] = w0 * flux[k0] + w1 * flux[k1] + w2 * flux[k2]

    return outimage

def pickThree (k0, k1, k2, k3, hex_x, hex_y):
    """Select three of four points close to the current pixel.

    k0 through k3 are indices into lists of rectangular coordinates (and
    fluxes) for the four points closest to the current pixel.  We need
    three points.  k0, k1, k2 would normally be the three to use, but
    they might be colinear, in which case we'll use k0, k1, k3.
    """

    CLOSE_TO_ZERO = 0.1

    dx0 = hex_x[k0]
    dx1 = hex_x[k1]
    dx2 = hex_x[k2]
    dx3 = hex_x[k3]

    dy0 = hex_y[k0]
    dy1 = hex_y[k1]
    dy2 = hex_y[k2]
    dy3 = hex_y[k3]

    d2 = (dx0 - dx1) * (dy0 - dy2) - (dx0 - dx2) * (dy0 - dy1)

    if abs (d2) < CLOSE_TO_ZERO:
        d3 = (dx0 - dx1) * (dy0 - dy3) - (dx0 - dx3) * (dy0 - dy1)
        if abs (d3) < CLOSE_TO_ZERO:
            raise RuntimeError
        return (k0, k1, k3)
    else:
        return (k0, k1, k2)

if __name__ == "__main__":

    main (sys.argv[1:])



