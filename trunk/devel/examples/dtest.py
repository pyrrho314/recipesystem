import sys, os
from astrodata import AstroData
from astrodata import Errors
import re

ad = AstroData("../../../test_data/recipedata/N20090703S0163.fits")

ybin = ad.detector_y_bin()
xbin = ad.detector_x_bin()
disp = ad.disperser()

descripts = ["airmass", "amp_read_area", "azimuth", "camera", "cass_rotator_pa",
 "central_wavelength", "coadds", "data_label", "data_section", "dec", "decker",
 "detector_section", "detector_x_bin", "detector_y_bin", "disperser",
 "dispersion", "dispersion_axis", "elevation", "exposure_time", "filter_name",
 "focal_plane_mask", "gain", "gain_setting", "grating", "instrument",
 "local_time", "mdf_row_id", "nod_count", "nod_pixels", "non_linear_level",
 "object", "observation_class", "observation_epoch", "observation_id",
 "observation_type", "pixel_scale", "prism", "program_id", "pupil_mask",
 "qa_state", "ra", "raw_bg", "raw_cc", "raw_iq", "raw_wv", "read_mode",
 "read_noise", "read_speed_setting", "saturation_level", "slit", "telescope",
 "ut_date", "ut_datetime", "ut_time", "wavefront_sensor",
 "wavelength_reference_pixel", "well_depth_setting", "x_offset", "y_offset"]

descripts = ["detector_x_bin","exposure_time", "raw_cc"]      
if True:
    
    ops = ["+","-","*","/","//", "%", "**", "<<",">>", "^", "<", "<=", ">",">=","==",]
    #ops = ["<", "<=", ">",">=","==",]
    exprs = []
    operands = ["'hello'", 10., 10]
    for op in ops:
        for operand in operands:
            expr = "%(lop)s %(op)s %(rop)s" % { "lop": "dval",
                                                "rop": repr(operand),
                                                "op":  op}
            exprs.append(expr)
            expr = "%(lop)s %(op)s %(rop)s" % { "lop": repr(operand),
                                                "rop": "dval",
                                                "op":  op}
            exprs.append(expr)
        expr = "%(lop)s %(op)s %(rop)s" % { "lop": "dval",
                                            "rop": "dval",
                                            "op":  op}
        exprs.append(expr)

    # print repr(exprs)
    for desc in descripts:
        print "CHECKING %s" % desc
        try:
            dval = eval("ad.%s()"%desc)
            print dval
            pydval = dval.pytype(dval)
            print dval,"|", dval.name, "|", repr(dval.pytype)
            for expr in exprs:
                try:
                    le = len(expr)
                    exval = eval(expr)
                    oute = "%s%s%s%s" % ( expr,
                                            " "*(15-le),
                                            " ==> ", 
                                            str(exval))
                    
                    oute += " "*(40-len(oute))
                except Errors.IncompatibleOperand:
                    oute = "IncompatibleOperand: type(dval)= %s " % str(dval.pytype) + expr
                    exval = "FAILED!"
                    print oute
                except TypeError, e:
                    oute = "(%s) TypeError: %s" %(expr,str(e))
                    exval = "FAILED!"
                    
                except Errors.DescriptorValueTypeError, e:
                    oute = "(%s) DVTE: %s" %(expr, str(e))
                    msg = str(e)
                    exval = "FAILED!"
                    print "FAILED %s: %s" %(expr, msg)

                print oute, "||",
                try:
                    controlexpr = re.sub("dval", repr(pydval), expr)
                    cexval = eval(controlexpr)
                    oute = "%s%s%s%s" % ( controlexpr,
                                            " "*(15-le),
                                            " ==> ", 
                                            str(cexval))
                    oute += " "*(60-len(oute))
                except :
                    print "CONTROL EXPR FAILED: ",controlexpr
                    continue
                print oute,
                
                if type(cexval) != type(exval):
                    print "TYPES DIFFER!" + str(type(exval)) + str(type(cexval))
                elif cexval != exval:
                    print "VALUES DIFFER!"  + str(exval) + " != "+str(cexval)
                else:
                    print "Passed Test"

            exec("%s = dval" % desc)
        except KeyError:
            print "FAILED due to KeyError!" #, ad.exception_info
        except Errors.DescriptorTypeError:
            print "FAILED with DescriptorTypeError"
        except Errors.ExistError:
            print "FAILED with ExistError"
        except:
            print "FAILED"
            raise
        print ""

if False:

    v = 1.0/xbin
    y = 1<<xbin # all operators are overloaded
    #print v,y
    data = ad[0].data
    #print data
    d2 = data/xbin.for_numpy() # see craig for what numpy does if you leave of .for_numpy())
    #print d2

    print "\n==Old setup results for units==\n"
    print "centwave m (default) =  7.5e-07 <type 'float'>"
    print "centwave um =  0.75"
    print "centwave nm =  750.0:"
    print "centwave Angstroms =  7500.0"

    ad = AstroData("../../../test_data/gmosspect/N20020810S0117.fits")

    centwave = ad.central_wavelength()
    print "centwave m (default) = ", centwave, str(type(centwave))

    centwave = ad.central_wavelength(asMicrometers=True)
    print "centwave um = ", centwave

    centwave = ad.central_wavelength(asNanometers=True)
    print "centwave nm = ", centwave

    centwave = ad.central_wavelength(asAngstroms=True)
    print "centwave Angstroms = ", centwave


