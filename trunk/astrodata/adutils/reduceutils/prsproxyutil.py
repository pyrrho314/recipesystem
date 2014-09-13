#
#                                                                     QAP Gemini
#
#                                                                prsproxyutil.py
#                                                                        07-2013
# ------------------------------------------------------------------------------
# $Id: prsproxyutil.py 4348 2013-08-09 18:47:54Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 4348 $'[11:-2]
__version_date__ = '$Date: 2013-08-09 08:47:54 -1000 (Fri, 09 Aug 2013) $'[7:-2]
# ------------------------------------------------------------------------------
from os.path import join, basename
from xml.dom import minidom
from pprint  import pformat

from astrodata.Lookups import get_lookup_table
# ------------------------------------------------------------------------------
# bad place for this
try:
    CALURL_DICT   = get_lookup_table("Gemini/calurl_dict","calurl_dict")
    UPLOADPROCCAL = CALURL_DICT["UPLOADPROCCAL"]
    _CALMGR       = CALMGR      = CALURL_DICT["CALMGR"]
    _LOCALCALMGR  = LOCALCALMGR = CALURL_DICT["LOCALCALMGR"]
except:
    pass # do without
# ------------------------------------------------------------------------------
CALTYPEDICT = { "arc": "arc",
                "bias": "bias",
                "dark": "dark",
                "flat": "flat",
                "processed_arc":   "processed_arc",
                "processed_bias":   "processed_bias",
                "processed_dark":   "processed_dark",
                "processed_flat":   "processed_flat",
                "processed_fringe": "processed_fringe"}
# -----------------------------------------------------------------------------
RESPONSESTR = """%%%%%%%%%%%%%%%% Request Data BEGIN:
%(sequence)s
%%%%%%%%%% Request Data END

########## Calibration Server Response BEGIN:
%(response)s
########## Calibration Server Response END

########## Nones Report (descriptors that returned None):
%(nones)s
########## Note: all descriptors shown above, scroll up.
        """
# -----------------------------------------------------------------------------

def upload_calibration(filename):
    """Uploads a calibration file to the FITS Store.

    parameters: <string>, the file to be uploaded.
    return:     <void>
    """
    import urllib2

    fn  = basename(filename)
    url = join(UPLOADPROCCAL, fn)
    postdata = open(filename).read()

    try:
        rq = urllib2.Request(url)
        u = urllib2.urlopen(rq, postdata)
        response = u.read()
    except urllib2.HTTPError, error:
        contents = error.read()
        raise
    return


def calibration_search(rq, fullResult=False):
    import urllib, urllib2
    import datetime

    calserv_msg = None
    print "\n@ppu80: calibration_search() ...\n"

    if "descriptors" in rq and "ut_datetime" in rq["descriptors"]:
        utc = rq["descriptors"]["ut_datetime"]
        pyutc = datetime.datetime.strptime(utc.value, "%Y%m%dT%H:%M:%S")
        print "@ppu85: OBS UT Date Time:", pyutc
        rq["descriptors"].update({"ut_datetime":pyutc} )
    
    # if rq has "calurl_dict" use it!!!
    if "calurl_dict" in rq :
        CALMGR = rq["calurl_dict"]["CALMGR"]
        LOCALCALMGR = rq["calurl_dict"]["LOCALCALMGR"]
    else:
        CALMGR = _CALMGR
        LOCALCALMGR = _LOCALCALMGR   
    
    if "source" not in rq:
        source = "central"
    else:
        source = rq["source"]
    
    token = ""             # used for GETs, we're using the post method
    rqurl = None

    print "@ppu104: CALMGR:     ", CALMGR
    print "@ppu105: LOCALCALMGR:", LOCALCALMGR

    if source == "central" or source == "all":
        rqurl = join(CALMGR, CALTYPEDICT[rq['caltype']])
        print "ppu109: CENTRAL SEARCH: rqurl is "+ rqurl

    print "@ppu111: source: ", source
    if source == 'local' or (rqurl == None and source == "all"):
        rqurl = LOCALCALMGR % { "httpport": 8777,
                                "caltype":  CALTYPEDICT[rq['caltype']],
                                } # "tokenstr":tokenstr}
        print "@ppu116: LOCAL SEARCH: rqurl is " + rqurl

    rqurl = rqurl + "/%s" % rq["filename"]
    #print "prs113:", pprint.pformat(rq)
    ### send request
    sequence = [("descriptors", rq["descriptors"]), ("types", rq["types"])]
    postdata = urllib.urlencode(sequence)
    response = "CALIBRATION_NOT_FOUND"

    try:
        print "@ppu126: Request URL: ", rqurl
        calRQ = urllib2.Request(rqurl)
        if source == "local":
            u = urllib2.urlopen(calRQ, postdata)
        else:
            u = urllib2.urlopen(calRQ, postdata)
        response = u.read()
        print "@ppu133:", response
    except urllib2.HTTPError, error:
        calserv_msg = error.read()
        print "ppu136:HTTPError- server returns:", error.read()
        import traceback
        traceback.print_exc()
        return (None, calserv_msg)

    #response = urllib.urlopen(rqurl).read()
    print "prs141:", response

    if fullResult:
        return response

    nones = []
    descripts = rq["descriptors"]
    for desc in descripts:
        if descripts[desc] == None:
            nones.append(desc)

    preerr = RESPONSESTR % { "sequence": pformat(sequence),
                             "response": response.strip(),
                             "nones"   : ", ".join(nones) \
                             if len(nones) > 0 else "No Nones Sent" }

    try:
        dom = minidom.parseString(response)
        calel = dom.getElementsByTagName("calibration")
        calurlel = dom.getElementsByTagName('url')[0].childNodes[0]
        calurlmd5 = dom.getElementsByTagName('md5')[0].childNodes[0]
    except IndexError:
        print "No url for calibration in response, calibration not found"
        return (None, preerr)
    except:
        return (None, preerr)
        
    #print "prs70:", calurlel.data
    
    #@@TODO: test only 
    print "prspu171:", repr(calurlel.data)
    return (calurlel.data, calurlmd5.data)
