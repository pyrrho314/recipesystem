#! /usr/bin/env python

import os
import time
import copy

import pyfits

import irafutil

"""This file contains the following utilities:

    expandlist (input)
    removeExtension (images)
    appendSuffix (filename, suffix)
    replaceSuffix (filename, suffix)
    gemdate (zone="UT")
    gemhedit (filename=None, extension=0, keyword="", value="", comment="")
    hdrhedit (header=None, keyword="", value="", comment="")
    getkeys (keywords, filename, extension=0, default="not found",
             must_exist=False)
    getkey (keyword, filename, extension=0, default="not found",
            must_exist=False)
    gimverify (image, sci_ext="SCI", dq_ext="DQ")
    fieldsOfTable (input, fields=None)
    joinlists (list1, list2, delim=" ", missing="Missing", shortest=True)
    joinlines (input, delim=" ", missing="Missing", maxchars=161,
               shortest=True)
    printlog (text, logfile=None, verbose=True)
    getdata (filename)
    putdata (filename, data, mef=None, phdr=None, hdr=None, clobber=False)
    convertMEF (filename, template=None)
    copyFile (input, output)
    deleteFile (filename)
"""

# This is used by removeExtension(), appendSuffix(), and replaceSuffix(),
# which are defined in this file.
extensions = [".fits", ".fit", ".pl", ".imh", ".hhh", ".tab"]

gimverify_does_not_exist = 0    # does not exist
gimverify_MEF = 1               # exists and is a MEF file
gimverify_simple_FITS = 2       # exists and is a simple FITS file
gimverify_pl = 3                # exists and is pixel list image
gimverify_imh = 4               # exists and is an imh image
gimverify_hhh = 5               # exists and is an hhh image
gimverify_other = 6             # exists but is not one of the above types

def expandlist (input):
    """Convert a string of comma-separated names to a list of names.

    @param input: one or more names separated by commas;
        a name of the form '@filename' implies that 'filename' is
        the name of a file containing names
    @type input: string

    @return: list of the names in 'input'
    @rtype: list of strings
    """

    filenames = []
    atList (input, filenames)
    return filenames

def atList (input, filenames):
    """Either append the current name, or read contents if it's a file.

    @param input: one or more names (or @name) separated by commas
    @type input: string

    @param filenames: (modified in-place) a list of the names extracted
        from input, or from the contents of input if it's an '@file'
    @type filenames: list
    """

    input = input.strip()
    if not input:
        return

    if input[0] == '@' and input.find (',') < 0:
        # just a single word, and it begins with '@'
        line = irafutil.expandFileName (input[1:])  # expand environment var.
        fd = open (line)
        lines = fd.readlines()
        fd.close()
    else:
        # one or more words, and the first does not begin with '@'
        lines = input.split (',')

    for line in lines:
        line = line.strip()
        if line[0] == '@':
            atList (line, filenames)
        else:
            line = irafutil.expandFileName (line)
            filenames.append (line)

def appendFits (images):
    """Append ".fits" to each name in 'images' that lacks an extension.

    >>> print appendFits ('abc')
    abc.fits
    >>> print appendFits ('abc.fits')
    abc.fits
    >>> print appendFits (['abc', 'xyz.fits'])
    ['abc.fits', 'xyz.fits']

    @param images: a file name or a list of file names
    @type images: a string or a list of strings

    @return: the input file names with ".fits" appended to each, unless
        the name already ended in a recognized extension.
    @rtype: list of strings
    """

    if isinstance (images, str):
        is_a_list = False
        images = [images]
    else:
        is_a_list = True
    modified = []
    for image in images:
        found = False
        # extensions is a list of recognized filename extensions.
        for extn in extensions:
            if image.endswith (extn):
                found = True
                break
        if found:
            modified.append (image)
        else:
            modified.append (image + ".fits")

    if is_a_list:
        return modified
    else:
        return modified[0]

def removeExtension (images):
    """Remove the extension from each file name in the list.

    If a file name does not end with one of the recognized extensions
    (identified by a variable 'extensions' that is global to this file),
    the original file name will be included unchanged in the output list.

    >>> print removeExtension ('abc.fits')
    abc
    >>> print removeExtension (['abc.fits', 'def.imh', 'ghi.fit'])
    ['abc', 'def', 'ghi']

    @param images: a file name or a list of file names
    @type images: a string or a list of strings

    @return: the input file names with extensions removed
    @rtype: list of strings
    """

    if isinstance (images, str):
        is_a_list = False
        images = [images]
    else:
        is_a_list = True
    modified = []
    for image in images:
        found = False
        # extensions is a list of recognized filename extensions.
        for extn in extensions:
            if image.endswith (extn):
                k = image.rfind (extn)
                modified.append (image[:k])
                found = True
                break
        if not found:
            modified.append (image)

    if is_a_list:
        return modified
    else:
        return modified[0]

def appendSuffix (filename, suffix):
    """Append a suffix to the root of the file name.

    If filename does not end with one of the recognized extensions,
    the suffix will just be appended to filename.

    >>> print appendSuffix ('abc.fits', '_flt')
    abc_flt.fits
    >>> print appendSuffix ('abc', '_flt')
    abc_flt
    >>> print appendSuffix ('abc.xyz', '_flt')
    abc.xyz_flt

    @param filename: a file name
    @type filename: string

    @param suffix: the suffix (e.g. "_flt") to append
    @type suffix: string

    @return: the input file name with the suffix included
    @rtype: string
    """

    found = False
    # extensions is a list of recognized filename extensions.
    for extn in extensions:
        if filename.endswith (extn):
            k = filename.rfind (extn)
            newname = filename[:k] + suffix + extn
            found = True
            break

    if not found:
        newname = filename + suffix

    return newname

def replaceSuffix (filename, suffix):
    """Replace the suffix in the file name.

    If filename includes an underscore ("_") character, the slice
    between that point (the rightmost underscore) and the extension will
    be replaced with the specified suffix.  If there is no underscore,
    the suffix will be inserted before the extension.
    If filename does not end with one of the recognized extensions,
    all of the string starting with the rightmost underscore will be
    replaced by the specified suffix, or the suffix will be appended
    if there is no underscore in filename.

    >>> print replaceSuffix ('abc_raw.fits', '_flt')
    abc_flt.fits
    >>> print replaceSuffix ('abc.fits', '_flt')
    abc_flt.fits
    >>> print replaceSuffix ('abc_raw.flub', '_flt')
    abc_flt
    >>> print replaceSuffix ('abc.flub', '_flt')
    abc.flub_flt

    @param filename: a file name
    @type filename: string

    @param suffix: the suffix to replace the existing suffix
    @type suffix: string

    @return: the input file name with the suffix included
    @rtype: string
    """

    found = False
    # extensions is a list of recognized filename extensions.
    for extn in extensions:
        if filename.endswith (extn):
            j = filename.rfind ("_")
            if j >= 0:
                newname = filename[:j] + suffix + extn
            else:
                k = filename.rfind (extn)
                newname = filename[:k] + suffix + extn
            found = True
            break

    if not found:
        j = filename.rfind ("_")
        if j >= 0:
            newname = filename[:j] + suffix
        else:
            newname = filename + suffix

    return newname

def gemdate (zone="UT"):
    """Get the current date and time.

    @param zone: "UT" or "local", to indicate whether the time should be
        UTC (always standard time) or local time (which can be either
        standard or daylight saving time)
    @type zone: string

    @return: date and time formatted as "yyyy-mm-ddThh:mm:ss";
        every value is an integer
    @rtype: string
    """

    if zone == "UT":
        zone_time = time.gmtime (time.time())
    elif zone == "local":
        zone_time = time.localtime (time.time())
    else:
        raise ValueError, "Invalid time zone = %s" % zone

    time_string = time.strftime ("%Y-%m-%dT%H:%M:%S", zone_time)

    return time_string

def gemhedit (filename=None, extension=0, keyword="", value="", comment="",
              delete=False):
    """Modify or delete an existing keyword or add a new keyword to a header.

    @param filename: name of FITS file containing the header to be updated
    @type filename: string

    @param extension: extension number, EXTNAME string,
        or (EXTNAME, EXTVER) tuple
    @type extension: int, string, or tuple

    @param keyword: keyword name
    @type keyword: string

    @param value: value to assign to the keyword
    @type value: string, int, float, or boolean

    @param comment: comment for the keyword; this is only used if the
        keyword is not already present in the header
    @type comment: string

    @param delete: if True and the keyword is present in the header,
        delete the keyword
    @type delete: boolean
    """

    fd = pyfits.open (filename, mode="update")
    header = fd[extension].header

    if delete:
        if header.has_key (keyword.upper()):
            del (header[keyword])
        # else do nothing
    else:
        if header.has_key (keyword.upper()):
            header[keyword] = value
        else:
            header.update (keyword, value=value, comment=comment)

    fd.close()

def hdrhedit (header=None, keyword="", value="", comment=""):
    """Modify an existing keyword or add a new keyword to a header.

    @param header: update the keyword in this header
    @type header: pyfits Header object

    @param keyword: keyword name
    @type keyword: string

    @param value: value to assign to the keyword
    @type value: string, int, float, or boolean

    @param comment: comment for the keyword; this is only used if the
        keyword is not already present in the header
    @type comment: string
    """

    if header.has_key (keyword.upper()):
        header[keyword] = value
    else:
        header.update (keyword, value=value, comment=comment)

def getkeys (keywords, filename, extension=0, default="not found",
             must_exist=False):
    """Get keyword values from a FITS header.

    @param keywords: the names of the keywords to read from the header
    @type keywords: list of strings

    @param filename: name of the FITS file
    @type filename: string

    @param extension: extension number, EXTNAME string,
        or (EXTNAME, EXTVER) tuple
    @type extension: int, string, or tuple

    @param default: value to assign for the value of missing keywords
    @type default: string

    @param must_exist: True if it is an error for any keyword in the
        list to be missing
    @type must_exist: boolean

    @return: tuple of values, one for each keyword in the input keywords list
    @rtype: tuple
    """

    fd = pyfits.open (filename)
    hdr = fd[extension].header
    fd.close()

    results = []
    missing = []
    cardlist = hdr.ascardlist()
    keywords_in_header = cardlist.keys()
    for keyword in keywords:
        keyword = keyword.upper()
        if keyword in keywords_in_header:
            value = hdr[keyword]
        else:
            missing.append (keyword)
            value = default
        results.append (value)

    if must_exist and missing:
        raise RuntimeError, "Missing keywords = %s" % repr (missing)

    return results

def getkey (keyword, filename, extension=0, default="not found",
            must_exist=False):

    keyword = keyword.upper()
    results = getkeys ([keyword], filename=filename, extension=extension,
                       default=default, must_exist=must_exist)

    return results[0]

def gimverify (image, sci_ext="SCI", dq_ext="DQ"):
    """Check whether the specified image exists, and if so, get the type.

    @param image: name of an image; the name must not use wildcard characters
    @type image: string

    @param sci_ext: name or number of the extension for science image data
    @type sci_ext: string or int

    @param dq_ext: name or number of the extension for data quality flags
    @type dq_ext: string or int

    @return: (type, has_dq)
        type is an integer code indicating existence and image type:
            gimverify_does_not_exist - does not exist or no read access
            gimverify_MEF - exists and is a multi-extension FITS (MEF) file;
                for this case we check that the sci extension is an IMAGE
            gimverify_simple_FITS - exists and is a simple FITS file; for this
                case we check that there actually is a data block in the
                primary header/data unit
            gimverify_pl - exists and is an iraf pixel list image
            gimverify_imh - exists and is an imh image
            gimverify_hhh - exists and is an hhh (GEIS) image
            gimverify_other - exists but is not one of the above image types
        has_dq is a boolean flag (always False unless gimverify_MEF):
            true indicates that the file contains a DQ extension (IMAGE)
    @rtype: tuple of int and boolean
    """

    # If the image name includes an extension specification (for FITS or
    # hhh format), strip it off.  Note that if the name includes a wildcard
    # using brackets, this will fail because part of the file name will be
    # chopped off.
    words = image.split ('[')
    if len (words) > 1:
        image = words[0]

    has_dq = False      # this is an initial value than can be reset below
    if not os.access (image, os.R_OK):
        return (gimverify_does_not_exist, has_dq)

    words = image.split ('.')
    if len (words) > 1:
        extension = words[-1]
    else:
        extension = ""

    if extension == "pl":
        type = gimverify_pl
    elif extension == "imh":
        type = gimverify_imh
    elif extension == "hhh":
        type = gimverify_hhh
    elif extension == "fits" or extension == "fit":
        # Find out what type of FITS file this is.
        fd = pyfits.open (image)
        if len (fd) > 1:
            type = gimverify_other      # may be reset below
            try:
                if fd[sci_ext].header["xtension"] == "IMAGE":
                    type = gimverify_MEF
                    try:
                        if fd[dq_ext].header["xtension"] == "IMAGE":
                            has_dq = True
                    except:
                        has_dq = False
            except KeyError:
                type = gimverify_other
        elif fd[0].data is not None:
            type = gimverify_simple_FITS
        else:
            type = gimverify_other
        fd.close()
    else:
        type = gimverify_other

    return (type, has_dq)

def fieldsOfTable (input, fields=None):
    """Extract the specified fields from a list or from a file.

    Blank lines (list elements) and lines that begin with '#' will be
    ignored.  Lines are assumed to contain at least as many words as
    specified in the 'fields' argument.
    Note:  currently does not ignore in-line comments
    Note:  This should return a record array.  xxx

    @param input: either the name of a text file or
        a list of whitespace-separated words
    @type input: either a string (file name) or a list of strings

    @param fields: zero-indexed column numbers, separated by commas or blanks
    @type fields: string, or None

    @return: list of the specified fields from the input
    @rtype: list of strings
    """

    # If the input is a single string, assume it's the name of a file
    # containing the table from which we should extract the columns.
    if isinstance (input, str):
        fd = open (input)
        input = fd.readlines()
        fd.close()

    if fields is None:
        return input

    # this is a string
    elements = fields.replace (",", " ")
    # this is a list of strings
    element_numbers_str = elements.split()
    # this is a list of integers
    element_numbers = [int (x) for x in element_numbers_str]

    result = []
    for line in input:
        words = line.split()
        nwords = len (words)
        # ignore blank lines and comments
        if nwords == 0 or words[0] == "#":
            continue
        output_line = []
        for i in element_numbers:
            if i >= nwords:
                break
            output_line.append (words[i])
        result.append (" ".join (output_line))

    return result

def joinlists (list1, list2, delim=" ", missing="Missing", shortest=True):
    """Join corresponding elements from two input lists.

    This is similar to the iraf.proto.joinlines task, except that the
    input is a pair of lists rather than files (just two input lists),
    and the result is as a list of strings, returned as the function
    value, rather than writing to standard output.  There is no verbose
    mode.  No warnings will be printed.

    @param list1: a list of values
    @type list1: list

    @param list2: another list of values
    @type list2: list

    @param delim: delimiter to separate joined elements
    @type delim: string

    @param missing: string to use for lists with fewer lines,
        if shortest is False
    @type missing: string

    @param shortest: if True, the number of elements in the function
        value will be the smaller of the number of lines in either of
        the input lists;
        if False, the number of elements will be the larger of the
        number lines in either input list
    @type shortest: Boolean

    @return: the contents of the input lists
    @rtype: list of strings
    """

    len1 = len (list1)
    len2 = len (list2)
    min_numlines = min (len1, len2)
    max_numlines = max (len1, len2)

    if min_numlines < max_numlines:
        if shortest:
            numlines = min_numlines
        else:
            numlines = max_numlines
    else:
        numlines = len1

    result = []
    for i in range (numlines):
        if i > len1-1:
            result.append (missing + delim + str (list2[i]))
        elif i > len2-1:
            result.append (str (list1[i]) + delim + missing)
        else:
            result.append (str (list1[i]) + delim + str (list2[i]))

    return result

def joinlines (input, delim=" ", missing="Missing",
               maxchars=161, shortest=True):
    """Join lines from the input list of files.

    This is an implementation of the iraf.proto.joinlines task, with
    the following differences:  The result is as a list of strings,
    returned as the function value, rather than writing to standard
    output.  There is no verbose mode.  No warnings will be printed.

    @param input: names of files, separated by commas (and optional
        whitespace)
    @type input: string

    @param delim: delimiter to separate joined lines
    @type delim: string

    @param missing: string to use for files with fewer lines,
        if shortest is False
    @type missing: string

    @param maxchars: the output strings will be truncated after this length
    @type maxchars: int

    @param shortest: if True, the number of elements in the function
        value will be the smallest number of lines in any input file;
        if False, the number of elements will be the largest number of
        lines in any input file
    @type shortest: Boolean

    @return: the contents of the input files
    @rtype: list of strings
    """

    filenames = input.split (",")
    if not filenames[0]:        # an empty string?
        return filenames

    for i in range (len (filenames)):
        filenames[i] = filenames[i].strip()

    # There will be one element of all_lines for each file in input;
    # all_lines[i] will be a list of the lines (with leading and
    # trailing whitespace removed) of file i from input.
    all_lines = []
    first = True
    for name in filenames:
        fd = open (name)
        lines = fd.readlines()
        fd.close()
        for i in range (len (lines)):
            lines[i] = lines[i].strip()
        all_lines.append (copy.deepcopy (lines))
        numlines = len (lines)
        if first:
            min_numlines = numlines
            max_numlines = numlines
            first = False
        else:
            min_numlines = min (numlines, min_numlines)
            max_numlines = max (numlines, max_numlines)

    if min_numlines < max_numlines:
        if shortest:
            numlines = min_numlines
        else:
            numlines = max_numlines

    if len (all_lines[0]) > numlines:
        result = all_lines[0][0:numlines]
    else:
        result = all_lines[0]
    for k in range (len (result), numlines):
        result.append (missing)

    for i in range (1, len (all_lines)):
        lines = all_lines[i]
        for j in range (len (lines)):
            if j >= numlines:
                break
            result[j] = result[j] + delim + lines[j]
        for j in range (len (lines), numlines):
            result[j] = result[j] + delim + missing

    for j in range (len (result)):
        result[j] = result[j][0:maxchars]

    return result

def printlog (text, logfile=None, verbose=True):
    """Append text to the log file.

    @param text: text string to log
    @type text: string

    @param logfile: name of log file, or None
    @type logfile: string

    @param verbose: if True, then also print to standard output
    @type verbose: boolean
    """

    if logfile == "STDOUT":
        logfile = None
        verbose = True

    if text[0:5] == "ERROR" or text[0:7] == "WARNING":
        verbose = True

    if logfile is not None:
        fd = open (logfile, mode="a")
        fd.write (text + "\n")
        fd.close()

    if verbose:
        print text

def getdata (filename):
    """Read an image from a FITS file.

    The image will be read from the primary header/data unit if the
    primary HDU contains a data block; in this case the third element
    (the extension header) of the function value will be None.
    If the primary HDU consists of only a header, the image will be
    read from the first extension; in this case the extension header
    will also be returned (i.e. will not be None).

    @param filename: name of the FITS file
    @type filename: string

    @return: data, primary header, and either extension header or None
    @rtype: tuple with three elements
    """

    fd = pyfits.open (filename)
    phdr = fd[0].header
    if fd[0].data is not None:
        # take the image from the primary HDU
        hdr = None
        data = fd[0].data
    else:
        # take the image from the first extension HDU
        hdr = fd[1].header
        data = fd[1].data

    return (data, phdr, hdr)

def putdata (filename, data, mef=None, phdr=None, hdr=None, clobber=False):
    """Write an image to a FITS file.

    The boolean mef argument lets you specify explicitly whether the
    file to be written should be simple FITS (primary HDU only) or
    multi-extension FITS (image in the first extension HDU).  If mef is
    None, this will be determined by whether the hdr argument is None.
    If hdr is None, phdr and data will be written to the primary HDU.
    If hdr is not None, only phdr will be written to the primary HDU, and
    hdr and data will be written to the first extension.

    @param filename: name of the FITS file
    @type filename: string

    @param data: image data
    @type data: array

    @param mef: if not None, explicitly specifies whether the output should
        be multi-extension FITS (mef=True) or simple FITS (mef=False)
    @type mef: boolean

    @param phdr: primary header
    @type phdr: pyfits Header object

    @param hdr: extension header or None
    @type hdr: pyfits Header object

    @param clobber: an existing file can be overwritten if clobber=True
    @type clobber: boolean
    """

    if mef is None:
        mef = (hdr is not None)

    if mef:
        # put the image in the first extension
        phdu = pyfits.PrimaryHDU (header=phdr)
        hdu = pyfits.ImageHDU (data, hdr)
        hdulist = pyfits.HDUList ([phdu, hdu])
    else:
        # put the image in the primary HDU
        phdu = pyfits.PrimaryHDU (data, phdr)
        hdulist = pyfits.HDUList (phdu)

    if clobber:
        os.remove (filename)

    hdulist.writeto (filename, output_verify="fix", clobber=clobber)

def convertMEF (filenames, output, extname=["SCI"], template=None):
    """Convert filenames to multi-extension FITS

    @param filenames: names of input FITS file
    @type filenames: string or list of strings

    @param output: name of output FITS file; if this file already exists,
        it will be overwritten
    @type output: string

    @param extname: extension names
    @type extname: string or list of strings

    @param template: name of input file whose primary header should be
        copied to the output primary header
    @type template: string
    """

    # take primary header for output from template
    fd = pyfits.open (template)
    phdu = pyfits.PrimaryHDU (header=fd[0].header)
    hdulist = pyfits.HDUList ([phdu])
    fd.close()

    if isinstance (filenames, str):
        filenames = [filenames]
    if isinstance (extname, str):
        extname = [extname]

    n_files = len (filenames)
    if len (extname) != n_files:
        print "input files and extension names =", filenames, extname
        raise RuntimeError, "lengths of lists must be the same"

    for n in range (n_files):
        fname = filenames[n]
        ename = extname[n]

        # hdr will be None, since 'fname' is a simple FITS file
        (data, phdr, hdr) = getdata (fname)

        # use primary header of 'fname' as extension header
        hdr = phdr.copy()

        # put the image in the extension
        hdr.update ("EXTNAME", ename)
        hdr.update ("EXTVER", 1)
        hdu = pyfits.ImageHDU (data, hdr)
        hdulist.append (hdu)

    if os.access (output, os.F_OK):
        os.remove (output)
    hdulist.writeto (output, output_verify="fix")

def copyFile (input, output):
    """Copy the text file 'input' to the file 'output'."""

    ifd = open (input)
    lines = ifd.readlines()
    ifd.close()
    ofd = open (output, "w")
    ofd.writelines (lines)
    ofd.close()

def deleteFile (filename):
    """Delete a file.

    The file 'filename' will be deleted if it exists, even if it does not
    have write access.

    @param filename: the name of a file
    @type filename: string
    """

    if os.access (filename, os.F_OK):
        os.remove (filename)

def crash (gl=None, errtype=RuntimeError, errmess=""):

    if gl is not None:
        gl.glogprint (loglevel="status", type="error", str=errmess)
        gl.glogclose (fl_success=False)

    if errmess:
        raise errtype, errmess
    else:
        raise errtype


def _test():
    import doctest, gemutil
    return doctest.testmod (gemutil)

if __name__ == "__main__":
    _test()
