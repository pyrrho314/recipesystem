import os, os.path

from pyraf import iraf

def expandFileName(filename):
    """Expand environment variable in a file name.
    If the input file name begins with either a Unix-style or IRAF-style
    environment variable (e.g. $lref/name_dqi.fits or lref$name_dqi.fits
    respectively), this routine expands the variable and returns a complete
    path name for the file.

    @param filename:  a file name, possibly including an environment
        variable
    @type filename:  string

    @return:  the file name with environment variable(s) expanded
    @rtype:  string
    """

    n = filename.find ("$")
    if n == 0:
        # Unix-style file name.
        filename = os.path.expandvars (filename)
        # If filename contains "//", delete one of them.
        double_sep = os.sep + os.sep
        i = filename.find (double_sep)
        if i != -1:
            filename = filename[:i+1] + filename[i+2:]
    elif n > 0:
        # IRAF-style file name.
        temp = iraf.envget(filename[0:n])
        len_temp = len(temp)
        if temp[len_temp-1] == os.sep:
            filename = temp + filename[n+1:]
        else:
            filename = temp + os.sep + filename[n+1:]

    if filename.find ("$") >= 0:
        filename = expandFileName (filename)

    return filename

def imcurXY():
    """Read the IRAF image cursor.

    @return:  the X and Y coordinates of the cursor
    @rtype:  tuple containing two floating-point values
    """

    tmp_cursor = iraf.cl.imcur
    words = tmp_cursor.split()
    return (float(words[0]), float(words[1]))


def splitSection(section):
    """Split an image section into its individual values.

    An image section would have the format '[x11:x12,y11:y12]',
    while alignsec (a pair of image sections) would have the format
    '[x11:x12,y11:y12] [x21:x22,y21:y22]'.  This function can take
    either as input and will return a list of floating-point values,
    four for each image section.

    @param section: a string containing one or more image section
        specifications in IRAF notation
    @type section: string

    @return: the numerical values that are the section limits;
        note that these are one indexed, since they are read directly
        from the image section string
    @rtype: list of float values
    """

    temp = section.replace (",", " ")
    temp = temp.replace (":", " ")
    temp = temp.replace ("[", "")
    temp = temp.replace ("]", "")
    words = temp.split()
    values = [float(x) for x in words]

    return values
