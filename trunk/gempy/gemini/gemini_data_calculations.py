#
#                                                                     QAP Gemini
#
#                                       gempy.gemini.gemini_data_calculations.py
#                                                                        08-2013
# ------------------------------------------------------------------------------
# $Id: gemini_data_calculations.py 4359 2013-08-22 19:41:34Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 4359 $'[11:-2]
__version_date__ = '$Date: 2013-08-22 09:41:34 -1000 (Thu, 22 Aug 2013) $'[7:-2]
# ------------------------------------------------------------------------------
# The gemini_data_calculations module contains functions that calculate values
# from Gemini data

from time import strptime
from datetime import datetime

import numpy as np

from astrodata import Lookups
from astrodata.gemconstants import SCI
from gempy.gemini import gemini_metadata_utils as gmu
# ------------------------------------------------------------------------------


def get_bias_level(adinput=None, estimate=True):
    if estimate:
        ret_bias_level = _get_bias_level_estimate(adinput=adinput)
    else:
        ret_bias_level = _get_bias_level(adinput=adinput)
    return ret_bias_level

    
def _get_bias_level(adinput=None):
    """
    Determine the bias level value from the science extensions of the input
    AstroData object. The bias level is equal to the median of the overscan
    region.
    """

    # Since this function accesses keywords in the headers of the pixel data
    # extensions, always construct a dictionary where the key of the dictionary
    # is an EXTVER integer
    ret_bias_level = {}
    
    # Get the overscan section value using the appropriate descriptor
    overscan_section_dv = adinput.overscan_section()
    
    # Create a dictionary where the key of the dictionary is an EXTVER integer
    overscan_section_dict = overscan_section_dv.collapse_by_extver()
    
    if not overscan_section_dv.validate_collapse_by_extver(
            overscan_section_dict):
        # The validate_collapse_by_extver function returns False if the values
        # in the dictionary with the same EXTVER are not equal 
        raise Errors.CollapseError()

    if overscan_section_dict is not None:

        # The type of CCD determines the number of contaminated columns in the
        # overscan region. Get the pretty detector name value using the
        # appropriate descriptor.

        detector_name_dv = adinput.detector_name(pretty=True)
        
        if detector_name_dv   == "EEV":       nbiascontam = 4
        elif detector_name_dv == "e2vDD":     nbiascontam = 5
        elif detector_name_dv == "Hamamastu": nbiascontam = 4
        else: nbiascontam = 4
        
        for extver, overscan_section in overscan_section_dict.iteritems():
            
            # Don't include columns at edges
            if overscan_section[0] == 0:
                # Overscan region is on the left
                overscan_section[1] -= nbiascontam
                overscan_section[0] += 1
            else:
                # Overscan region is on the right
                overscan_section[0] += nbiascontam
                overscan_section[1] -= 1
            
            # Extract overscan data. In numpy arrays, y indices come first.
            overscan_data = adinput[SCI,extver].data[
                overscan_section[2]:overscan_section[3],
                overscan_section[0]:overscan_section[1]]
            bias_level = np.median(overscan_data)
            
            # Update the dictionary with the bias level value
            ret_bias_level.update({extver: bias_level})
        
        unique_values = set(ret_bias_level.values())
        if len(unique_values) == 1 and None in unique_values:
            # The bias level was not found for any of the pixel data extensions
            # (all the values in the dictionary are equal to None)
            ret_bias_level = None
    else:
        ret_bias_level = _get_bias_level_estimate(adinput=adinput)
    
    return ret_bias_level


def _get_bias_level_estimate(adinput=None):
    """
    Determine an estiamte of the bias level value from GMOS data.
    """
    # Since this function accesses keywords in the headers of the pixel data
    # extensions, always construct a dictionary where the key of the dictionary
    # is an EXTVER integer
    ret_bias_level = {}
    
    # Get the overscan value and the raw bias level from the header of each
    # pixel data extension as a dictionary where the key of the dictionary is
    # an EXTVER integer
    keyword_value_dict = gmu.get_key_value_dict(adinput=adinput, 
                                                keyword=["OVERSCAN", "RAWBIAS"], 
                                                dict_key_extver=True)
    overscan_value_dict = keyword_value_dict["OVERSCAN"]
    raw_bias_level_dict = keyword_value_dict["RAWBIAS"]
    
    if overscan_value_dict is None:
        
        # If there is no overscan value for any extensions, use the raw bias
        # level value as the value for the bias level
        if raw_bias_level_dict is None:
            
            # If there is no raw bias level value for any extensions, use the
            # static bias levels from the lookup table as the value for the
            # bias level
            ret_bias_level = _get_static_bias_level(adinput=adinput)
        else:
            # Use the raw bias level value as the value for the bias level
            for extver, raw_bias_level in raw_bias_level_dict.iteritems():
                
                if raw_bias_level is None:
                    # If the raw bias level does not exist for a given
                    # extension, use the static bias levels from the lookup
                    # table as the value for the bias level
                    bias_level = _get_static_bias_level_for_ext(
                        adinput=adinput[SCI,extver])
                else:
                    bias_level = raw_bias_level
                
                # Update the dictionary with the bias level value
                ret_bias_level.update({extver: bias_level})
    else:
        for extver, overscan_value in overscan_value_dict.iteritems():
            if overscan_value is None:
                
                # If the overscan value does not exist for a given extension,
                # use the raw bias level value from the header as the value for
                # the bias level
                if raw_bias_level_dict is None:
                    
                    # If there is no raw bias level value for any extensions,
                    # use the static bias levels from the lookup table as the
                    # value for the bias level 
                    bias_level = _get_static_bias_level_for_ext(
                        adinput=adinput[SCI,extver])
                else:
                    raw_bias_level = raw_bias_level_dict[extver]
                    if raw_bias_level is None:
                        
                        # If the raw bias level does not exist for a given 
                        # extension, use the static bias levels from the lookup
                        # table as the value for the bias level 
                        bias_level = _get_static_bias_level_for_ext(
                            adinput=adinput[SCI,extver])
                    else:
                        bias_level = raw_bias_level
            else:
                bias_level = overscan_value
            
            # Update the dictionary with the bias level value
            ret_bias_level.update({ext_name_ver: bias_level})
    
    unique_values = set(ret_bias_level.values())

    if len(unique_values) == 1 and None in unique_values:
        # The bias level was not found for any of the pixel data extensions
        # (all the values in the dictionary are equal to None)
        ret_bias_level = None
        
    return ret_bias_level


def _get_static_bias_level(adinput=None):
    """
    Determine the static bias level value from GMOS data.
    """

    # Since this function accesses keywords in the headers of the pixel data
    # extensions, always construct a dictionary where the key of the dictionary
    # is an EXTVER integer
    static_bias_level = {}
    
    # Get the static bias level lookup table
    gmosampsBias, gmosampsBiasBefore20060831 = Lookups.get_lookup_table(
        "Gemini/GMOS/GMOSAmpTables", "gmosampsBias",
        "gmosampsBiasBefore20060831")
    
    # Get the UT date, read speed setting and gain setting values using the
    # appropriate descriptors
    ut_date_dv = adinput.ut_date()
    read_speed_setting_dv = adinput.read_speed_setting()
    gain_setting_dv = adinput.gain_setting()
    
    # Get the name of the detector amplifier from the header of each pixel data
    # extension as a dictionary
    ampname_dict = gmu.get_key_value_dict(
        adinput=adinput, keyword="AMPNAME", dict_key_extver=True)
    
    if not (ut_date_dv.is_none() and read_speed_setting_dv.is_none() and
            gain_setting_dv.is_none()) and ampname_dict is not None:
        
        # Use as_pytype() to return the values as the default python type
        # rather than an object
        ut_date = str(ut_date_dv)
        read_speed_setting = read_speed_setting_dv.as_pytype()
        
        # Create a gain setting dictionary where the key of the dictionary is
        # an EXTVER integer
        gain_setting_dict = gain_setting_dv.collapse_by_extver()
        
        if not gain_setting_dv.validate_collapse_by_extver(gain_setting_dict):
            # The validate_collapse_by_extver function returns False if the
            # values in the dictionary with the same EXTVER are not equal 
            raise Errors.CollapseError()
        
        obs_ut_date = datetime(*strptime(ut_date, "%Y-%m-%d")[0:6])
        old_ut_date = datetime(2006, 8, 31, 0, 0)
        
        for extver, gain_setting in gain_setting_dict.iteritems():
            ampname  = ampname_dict[extver]
            bias_key = (read_speed_setting, gain_setting, ampname)

            bias_level = None
            if obs_ut_date > old_ut_date:
                if bias_key in gmosampsBias:
                    bias_level = gmosampsBias[bias_key]
            else:
                if bias_key in gmosampsBiasBefore20060831:
                    bias_level = gmosampsBiasBefore20060831[bias_key]
            
            # Update the dictionary with the bias level value
            static_bias_level.update({extver: bias_level})
    
    # if len(static_bias_level) == 1:
    #     # Only one value will be returned
    #     ret_static_bias_level = static_bias_level.values()[0]  #!! Not a dict !
    # else:

    unique_values = set(static_bias_level.values())

    if len(unique_values) == 1 and None in unique_values:
        # The bias level was not found for any of the pixel data extensions
        # (all the values in the dictionary are equal to None)
        ret_static_bias_level = None
    else:
        ret_static_bias_level = static_bias_level
    
    return ret_static_bias_level


def _get_static_bias_level_for_ext(adinput=None):
    """
    Determine the static bias level value from GMOS data.
    """

    # Get the static bias level lookup table
    gmosampsBias, gmosampsBiasBefore20060831 = Lookups.get_lookup_table(
        "Gemini/GMOS/GMOSAmpTables", "gmosampsBias",
        "gmosampsBiasBefore20060831")
    
    # Get the UT date, read speed setting and gain setting values using the
    # appropriate descriptors
    ut_date_dv = adinput.ut_date()
    gain_setting_dv = adinput.gain_setting()
    read_speed_setting_dv = adinput.read_speed_setting()
    
    # Get the name of the detector amplifier from the header of the specific
    # pixel data extension
    ampname = adinput.get_key_value("AMPNAME")
    ret_static_bias_level = None

    if not (ut_date_dv.is_none() and read_speed_setting_dv.is_none() and
            gain_setting.dv.is_none()) and ampname is not None:
        # Use as_pytype() to return the values as the default python type
        # rather than an object
        ut_date = str(adinput.ut_date())
        gain_setting = adinput.gain_setting().as_pytype()
        read_speed_setting = adinput.read_speed_setting().as_pytype()

        obs_ut_date = datetime(*strptime(ut_date, "%Y-%m-%d")[0:6])
        old_ut_date = datetime(2006, 8, 31, 0, 0)
        
        bias_key = (read_speed_setting, gain_setting, ampname)
        ret_static_bias_level = None

        if obs_ut_date > old_ut_date:
            if bias_key in gmosampsBias:
                ret_static_bias_level = gmosampsBias[bias_key]
        else:
            if bias_key in gmosampsBiasBefore20060831:
                ret_static_bias_level = gmosampsBiasBefore20060831[bias_key]
    
    return ret_static_bias_level
