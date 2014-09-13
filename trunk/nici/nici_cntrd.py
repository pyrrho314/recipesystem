import numpy as np
from niciTools import congrid,robust_sigma,gcentroid
import scipy.ndimage as nd
from peak2halo import peak2halo
import scipy.signal.signaltools
try:
    import stsci.numdisplay as ndis
except ImportError:
    import numdisplay as ndis
import sys


def nici_cntrd(im,hdr,center_im=True):
    """
    Read xcen,ycen and update header if necessary
    If the automatic finding of the center mask fails, then
    the interactive session will start. The SAOIMAGE/ds9
    display must be up. If the default port is busy, then it will
    use port 5199, so make sure you start "ds9 -port 5199".
    """

    xcen = hdr.get('xcen')
    ycen = hdr.get('ycen')
    updated = False
    if (xcen == None):
        ratio,xc,yc = peak2halo('',image=im)
        #xcen = xc[0]
        #ycen = yc[0]
        xcen = xc
        ycen = yc
        if (xcen < 0 or ycen < 0):
            try:
                ndis.display(im)
            except IOError,err:
                sys.stderr.write('\n ***** ERROR: %s Start DS9.\n' % str(err))
                sys.exit(1)
            print " Mark center with left button, then use 'q' to continue, 's' to skip."
            cursor = ndis.readcursor(sample=0)
            cursor = cursor.split()
            if cursor[3] == 's':
                hdr.update("XCEN",-1, "Start mask x-center")
                hdr.update("YCEN",-1, "Start mask y-center")
                updated = True 
                print '\nFrame skipped... ****Make sure not to use it in your science script. ***\n'
                #return updated,xcen,ycen,im
                return xcen,ycen,im
            x1 = float(cursor[0])
            y1 = float(cursor[1])


            box = im[y1-64:y1+64,x1-64:x1+64].copy()
            box -= scipy.signal.signaltools.medfilt2d(np.float32(box),11)
            box = box[32:32+64,32:32+64]

            bbx = box * ((box>(-robust_sigma(box)*5)) & \
                         (box <(15*robust_sigma(box))))

            imbb = congrid(bbx,(1024,1024))
            ndis.display(imbb, name='bbx')
            del imbb

            cursor = ndis.readcursor(sample=0)
            cursor = cursor.split()
            x2 = float(cursor[0])
            y2 = float(cursor[1])

            xcen,ycen = gcentroid(box, x2/16., y2/16., 4)

            xcen = (xcen+x1)[0] - 32
            ycen = (ycen+y1)[0] - 32

        hdr.update("XCEN",xcen, "Start mask x-center")
        hdr.update("YCEN",ycen, "Start mask y-center")

        updated = True 
    if center_im:
        # Shift the image. Use ndimage shift function. Make sure
        # the array is float.
        im = np.asarray(np.nan_to_num(im),dtype=np.float32)
        im = nd.shift (im,(512-ycen,512-xcen))     

    #return updated,xcen,ycen,im
    return xcen,ycen,im
