from datetime import datetime
from time import strptime
import math
import numpy as np

from astrodata import Errors
from astrodata import Lookups
from astrodata.Descriptors import DescriptorValue
from astrodata.structuredslice import pixel_exts, bintable_exts
from gempy.gemini import gemini_data_calculations as gdc
from gempy.gemini import gemini_metadata_utils as gmu
import GemCalcUtil

from GMOS_Keywords import GMOS_KeyDict
from GEMINI_Descriptors import GEMINI_DescriptorCalc

class GMOS_DescriptorCalc(GEMINI_DescriptorCalc):
    # Updating the global key dictionary with the local key dictionary
    # associated with this descriptor class
    _update_stdkey_dict = GMOS_KeyDict
    
    def __init__(self):
        GEMINI_DescriptorCalc.__init__(self)
    
    def amp_read_area(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_amp_read_area_dict = {}
        
        # Determine the name of the detector amplifier keyword (ampname) from
        # the global keyword dictionary 
        keyword = self.get_descriptor_key("key_ampname")
        
        # Get the value of the name of the detector amplifier keyword from the
        # header of each pixel data extension as a dictionary where the key of
        # the dictionary is an ("*", EXTVER) tuple
        ampname_dict = gmu.get_key_value_dict(adinput=dataset, keyword=keyword)
        
        if ampname_dict is None:
            # The get_key_value_dict() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        # Get the pretty (1-based indexing) readout area of the CCD using the
        # appropriate descriptor
        detector_section_dv = dataset.detector_section(pretty=True)
        
        if detector_section_dv.is_none():
            # The descriptor functions return None if a value cannot be found
            # and stores the exception info. Re-raise the exception. It will be
            # dealt with by the CalculatorInterface. 
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        # Use as_dict() to return the detector section value as a dictionary
        # where the key of the dictionary is an ("*", EXTVER) tuple 
        detector_section_dict = detector_section_dv.as_dict()
        
        for ext_name_ver, ampname in ampname_dict.iteritems():
            detector_section = detector_section_dict[ext_name_ver]
            
            if ampname is None or detector_section is None:
                amp_read_area = None
            else:
                # Use the composite amp_read_area string as the value
                amp_read_area = "'%s':%s" % (ampname, detector_section)
            
            # Update the dictionary with the amp_read_area value
            ret_amp_read_area_dict.update({ext_name_ver:amp_read_area})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_amp_read_area_dict, name="amp_read_area",
                                 ad=dataset)
        return ret_dv
    
    def array_name(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        array_name_dict = {}
        
        # Determine the name of the array keyword from the global keyword
        # dictionary
        keyword = self.get_descriptor_key("key_array_name")
        
        # Get the value of the name of the array keyword from the header of
        # each pixel data extension as a dictionary where the key of the
        # dictionary is an ("*", EXTVER) tuple
        array_name_dict = gmu.get_key_value_dict(adinput=dataset,
                                                 keyword=keyword)
        if array_name_dict is None:
            # It is possible that the data have been mosaiced, which means that
            # the name of the array keyword no longer exists in the pixel data
            # extensions. Instead, determine the value of the detector name
            # using the appropriate descriptor
            detector_name_dv = dataset.detector_name()
            
            if detector_name_dv.is_none():
                # The descriptor functions return None if a value cannot be
                # found and stores the exception info. Re-raise the exception.
                # It will be dealt with by the CalculatorInterface. 
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            ret_array_name = detector_name_dv
        else:
            ret_array_name = array_name_dict
        
        # Instantiate the return DescriptorValue (DV) object using the newly
        # created dictionary
        ret_dv = DescriptorValue(ret_array_name, name="array_name", ad=dataset)
        
        return ret_dv
    
    def central_wavelength(self, dataset, asMicrometers=False,
                           asNanometers=False, asAngstroms=False, **args):
        # Currently for GMOS data, the central wavelength is recorded in
        # nanometers
        input_units = "nanometers"
        
        # Determine the output units to use
        unit_arg_list = [asMicrometers, asNanometers, asAngstroms]
        
        if unit_arg_list.count(True) == 1:
            # Just one of the unit arguments was set to True. Return the
            # central wavelength in these units.
            if asMicrometers:
                output_units = "micrometers"
            if asNanometers:
                output_units = "nanometers"
            if asAngstroms:
                output_units = "angstroms"
        else:
            # Either none of the unit arguments were set to True or more than
            # one of the unit arguments was set to True. In either case,
            # return the central wavelength in the default units of meters.
            output_units = "meters"
        
        # Determine the central wavelength keyword from the global keyword
        # dictionary
        keyword = self.get_descriptor_key("key_central_wavelength")
        
        # Get the value of the central wavelength keyword from the header of
        # the PHU
        raw_central_wavelength = dataset.phu_get_key_value(keyword)
        
        if raw_central_wavelength is None:
            # The phu_get_key_value() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        else:
            central_wavelength = float(raw_central_wavelength)
        
        # Validate the central wavelength value
        if central_wavelength < 0.0:
            raise Errors.InvalidValueError()
        else:
            # Use the utilities function convert_units to convert the central
            # wavelength value from the input units to the output units
            ret_central_wavelength = GemCalcUtil.convert_units(
                input_units=input_units, input_value=central_wavelength,
                output_units=output_units)
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_central_wavelength,
                                 name="central_wavelength", ad=dataset)
        return ret_dv
    
    def detector_name(self, dataset, pretty=False, **args):
        # Determine the name of the detector keyword from the global keyword
        # dictionary
        keyword1 = self.get_descriptor_key("key_detector_name")
        
        # Get the value of the name of the detector keyword from the header of
        # the PHU
        detector_name = dataset.phu_get_key_value(keyword1)
        
        if detector_name is None:
            # The phu_get_key_value() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        if pretty:
            # Define relationship between the type of the detector and the
            # pretty name of the detector
            pretty_detector_name_dict = {
                "SDSU II CCD": "EEV",
                "SDSU II e2v DD CCD42-90": "e2vDD",
                "S10892-01": "Hamamatsu",
                }
            
            # Determine the type of the detector keyword from the global
            # keyword dictionary
            keyword2 = self.get_descriptor_key("key_detector_type")
            
            # Get the value of the type of the detector keyword from the header
            # of the PHU
            detector_type = dataset.phu_get_key_value(keyword2)
            
            # Return the pretty name of the detector
            if detector_type in pretty_detector_name_dict:
                ret_detector_name = pretty_detector_name_dict[detector_type]
            else:
                raise Errors.TableKeyError()
        else:
            # Return the name of the detectory
            ret_detector_name = detector_name
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_detector_name, name="detector_name",
                                 ad=dataset)
        return ret_dv
    
    def detector_rois_requested(self, dataset, **args):
        # This parses the DETROx GMOS headers and returns a list of ROIs in the
        # form [x1, x2, y1, y2]. These are in physical pixels - should be
        # irrespective of binning. These are 1-based, and inclusive, 2 pixels,
        # starting at pixel 2 would be [2, 3].
        ret_detector_rois_requested_list = []
        
        # Must be single digit ROI number
        for i in range(1, 10):
            x1 = dataset.phu_get_key_value("DETRO%sX" % i)
            xs = dataset.phu_get_key_value("DETRO%sXS" % i)
            y1 = dataset.phu_get_key_value("DETRO%sY" % i)
            ys = dataset.phu_get_key_value("DETRO%sYS" % i)
            
            if x1 is not None:
                # The headers are in the form of a start position and size
                # so make them into start and end pixels here
                xs *= int(dataset.detector_x_bin())
                ys *= int(dataset.detector_y_bin())
                ret_detector_rois_requested_list.append(
                    [x1, x1+xs-1, y1, y1+ys-1])
            else:
                break
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_detector_rois_requested_list,
                                 name="detector_rois_requested", ad=dataset)
        return ret_dv
    
    def detector_roi_setting(self, dataset, **args):
        """
        Attempts to deduce the Name of the ROI, more or less as per the options
        you can select in the OT. Only considers the first ROI.
        
        """
        # Get the lookup table containing the ROI sections
        gmosRoiSettings = Lookups.get_lookup_table("Gemini/GMOS/ROItable",
                                                   "gmosRoiSettings")
        
        ret_detector_roi_setting = "Undefined"
        rois = dataset.detector_rois_requested()
        if rois.is_none():
            # The descriptor functions return None if a value cannot be
            # found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info

        rois = rois.as_list()
        if rois:
            roi = rois[0]
            
            # If we don't recognise it, it's "Custom"
            ret_detector_roi_setting = "Custom"
            
            for s in gmosRoiSettings.keys():
                if roi in gmosRoiSettings[s]:
                    ret_detector_roi_setting = s
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_detector_roi_setting,
                                 name="detector_roi_setting", ad=dataset)
        return ret_dv
    
    def detector_x_bin(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_detector_x_bin_dict = {}
        
        # Determine the ccdsum keyword from the global keyword dictionary 
        keyword = self.get_descriptor_key("key_ccdsum")
        
        # Get the value of the ccdsum keyword from the header of each pixel
        # data extension as a dictionary where the key of the dictionary is an
        # ("*", EXTVER) tuple
        ccdsum_dict = gmu.get_key_value_dict(adinput=dataset, keyword=keyword)
        
        if ccdsum_dict is None:
            # The get_key_value_dict() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        for ext_name_ver, ccdsum in ccdsum_dict.iteritems():
            if ccdsum is None:
                detector_x_bin = None
            else:
                # Use the binning of the x-axis integer as the value
                x_bin, y_bin = ccdsum.split()
                detector_x_bin = int(x_bin)
            
            # Update the dictionary with the binning of the x-axis value
            ret_detector_x_bin_dict.update({ext_name_ver:detector_x_bin})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_detector_x_bin_dict,
                                 name="detector_x_bin", ad=dataset)
        return ret_dv
    
    def detector_y_bin(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_detector_y_bin_dict = {}
        
        # Determine the ccdsum keyword from the global keyword dictionary 
        keyword = self.get_descriptor_key("key_ccdsum")
        
        # Get the value of the ccdsum keyword from the header of each pixel
        # data extension as a dictionary where the key of the dictionary is an
        # ("*", EXTVER) tuple
        ccdsum_dict = gmu.get_key_value_dict(adinput=dataset, keyword=keyword)
        
        if ccdsum_dict is None:
            # The get_key_value_dict() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        for ext_name_ver, ccdsum in ccdsum_dict.iteritems():
            if ccdsum is None:
                detector_y_bin = None
            else:
                # Use the binning of the y-axis integer as the value
                x_bin, y_bin = ccdsum.split()
                detector_y_bin = int(y_bin)
            
            # Update the dictionary with the binning of the y-axis value
            ret_detector_y_bin_dict.update({ext_name_ver:detector_y_bin})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_detector_y_bin_dict,
                                 name="detector_y_bin", ad=dataset)
        return ret_dv
    
    def disperser(self, dataset, stripID=False, pretty=False, **args):
        # Determine the disperser keyword from the global keyword dictionary
        keyword = self.get_descriptor_key("key_disperser")
        
        # Get the value of the disperser keyword from the header of the PHU
        disperser = dataset.phu_get_key_value(keyword)
        
        if disperser is None:
            # The phu_get_key_value() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        if pretty:
            # If pretty=True, use stripID then additionally remove the
            # trailing "+" from the string
            stripID = True
        
        if stripID:
            if pretty:
                # Return the stripped and pretty disperser string
                ret_disperser = gmu.removeComponentID(disperser).strip("+")
            else:
                # Return the stripped disperser string
                ret_disperser = gmu.removeComponentID(disperser)
        else:
            # Return the disperser string
            ret_disperser = str(disperser)
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_disperser, name="disperser", ad=dataset)
        
        return ret_dv
    
    def dispersion(self, dataset, asMicrometers=False, asNanometers=False,
                   asAngstroms=False, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_dispersion_dict = {}
        
        # Currently for GMOS data, the dispersion is recorded in meters (?)
        input_units = "meters"
        
        # Determine the output units to use
        unit_arg_list = [asMicrometers, asNanometers, asAngstroms]
        
        if unit_arg_list.count(True) == 1:
            # Just one of the unit arguments was set to True. Return the
            # dispersion in these units.
            if asMicrometers:
                output_units = "micrometers"
            if asNanometers:
                output_units = "nanometers"
            if asAngstroms:
                output_units = "angstroms"
        else:
            # Either none of the unit arguments were set to True or more than
            # one of the unit arguments was set to True. In either case,
            # return the dispersion in the default units of meters.
            output_units = "meters"
        
        # Determine the dispersion keyword from the global keyword dictionary 
        keyword = self.get_descriptor_key("key_dispersion")
        
        # Get the value of the dispersion keyword from the header of each pixel
        # data extension as a dictionary where the key of the dictionary is an
        # ("*", EXTVER) tuple
        dispersion_dict = gmu.get_key_value_dict(adinput=dataset,
                                                 keyword=keyword)
        if dispersion_dict is None:
            # The get_key_value_dict() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        for ext_name_ver, raw_dispersion in dispersion_dict.iteritems():
            if raw_dispersion is None:
                dispersion = None
            else:
                # Use the utilities function convert_units to convert the
                # dispersion value from the input units to the output units
                dispersion = float(GemCalcUtil.convert_units(
                    input_units=input_units, input_value=float(raw_dispersion),
                    output_units=output_units))
            
            # Update the dictionary with the dispersion value
            ret_dispersion_dict.update({ext_name_ver:dispersion})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_dispersion_dict, name="dispersion",
                                 ad=dataset)
        return ret_dv
    
    def dispersion_axis(self, dataset, **args):
        # The GMOS dispersion axis should always be 1
        ret_dispersion_axis = 1
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_dispersion_axis, name="dispersion_axis",
                                 ad=dataset)
        return ret_dv
    
    def exposure_time(self, dataset, **args):
        # Determine the exposure time keyword from the global keyword
        # dictionary
        keyword = self.get_descriptor_key("key_exposure_time")
        
        # Get the value of the exposure time keyword from the header of the PHU
        exposure_time = dataset.phu_get_key_value(keyword)
        
        if exposure_time is None:
            # The phu_get_key_value() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        # Sanity check for times when the GMOS DC is stoned
        if exposure_time > 10000. or exposure_time < 0.:
            raise Errors.InvalidValueError()
        else:
            # Return the exposure time float
            ret_exposure_time = float(exposure_time)
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_exposure_time, name="exposure_time",
                                 ad=dataset)
        return ret_dv
    
    def focal_plane_mask(self, dataset, stripID=False, pretty=False, **args):
        # Determine the focal plane mask keyword from the global keyword
        # dictionary
        keyword = self.get_descriptor_key("key_focal_plane_mask")
        
        # Get the focal plane mask value from the header of the PHU
        focal_plane_mask = dataset.phu_get_key_value(keyword)
        
        if focal_plane_mask is None:
            # The phu_get_key_value() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        if focal_plane_mask == "None":
            ret_focal_plane_mask = "Imaging"
        else:
            # Return the focal plane mask string
            ret_focal_plane_mask = str(focal_plane_mask)
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_focal_plane_mask, name="focal_plane_mask",
                                 ad=dataset)
        return ret_dv
    
    def gain(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_gain_dict = {}
        
        # If the data have been prepared, take the gain value directly from the
        # appropriate keyword. At some point, a check for the STDHDRSI header
        # keyword should be added, since the function that overwrites the gain
        # keyword also writes the STDHDRSI keyword.
        if "PREPARED" in dataset.types:
            
            # Determine the gain keyword from the global keyword dictionary 
            keyword = self.get_descriptor_key("key_gain")
            
            # Get the value of the gain keyword from the header of each pixel
            # data extension as a dictionary where the key of the dictionary is
            # an ("*", EXTVER) tuple
            gain_dict = gmu.get_key_value_dict(adinput=dataset,
                                               keyword=keyword)
            if gain_dict is None:
                # The get_key_value_dict() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            ret_gain_dict = gain_dict
        
        else:
            # Get the lookup table containing the gain values by amplifier
            gmosampsGain, gmosampsGainBefore20060831 = (
                Lookups.get_lookup_table("Gemini/GMOS/GMOSAmpTables",
                                         "gmosampsGain",
                                         "gmosampsGainBefore20060831"))
            
            # Determine the amplifier integration time keyword (ampinteg) from
            # the global keyword dictionary
            keyword = self.get_descriptor_key("key_ampinteg")
            
            # Get the value of the amplifier integration time keyword from the
            # header of the PHU
            ampinteg = dataset.phu_get_key_value(keyword)
            
            if ampinteg is None:
                # The phu_get_key_value() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            # Get the UT date, gain setting and read speed setting values using
            # the appropriate descriptors
            ut_date_dv = dataset.ut_date()
            gain_setting_dv = dataset.gain_setting()
            read_speed_setting_dv = dataset.read_speed_setting()
            
            if (ut_date_dv.is_none() or gain_setting_dv.is_none() or
                read_speed_setting_dv.is_none()):
                # The descriptor functions return None if a value cannot be
                # found and stores the exception info. Re-raise the exception.
                # It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            # Use as_dict() and as_pytype() to return the values as a
            # dictionary where the key of the dictionary is an ("*", EXTVER)
            # tuple and the default python type, respectively, rather than an
            # object
            ut_date = str(ut_date_dv)
            gain_setting_dict = gain_setting_dv.as_dict()
            read_speed_setting = read_speed_setting_dv.as_pytype()
            
            obs_ut_date = datetime(*strptime(ut_date, "%Y-%m-%d")[0:6])
            old_ut_date = datetime(2006, 8, 31, 0, 0)
            
            # Determine the name of the detector amplifier keyword (ampname)
            # from the global keyword dictionary 
            keyword = self.get_descriptor_key("key_ampname")
            
            # Get the value of the name of the detector amplifier keyword from
            # the header of each pixel data extension as a dictionary where the
            # key of the dictionary is an ("*", EXTVER) tuple
            ampname_dict = gmu.get_key_value_dict(adinput=dataset,
                                                  keyword=keyword)
            if ampname_dict is None:
                # The get_key_value_dict() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            for ext_name_ver, ampname in ampname_dict.iteritems():
                gain_setting = gain_setting_dict[ext_name_ver]
                
                if ampname is None or gain_setting is None:
                    gain = None
                else:
                    gain_key = (read_speed_setting, gain_setting, ampname)
                    
                    if obs_ut_date > old_ut_date:
                        if gain_key in gmosampsGain:
                            gain = gmosampsGain[gain_key]
                        else:
                            raise Errors.TableKeyError()
                    else:
                        if gain_key in gmosampsGainBefore20060831:
                            gain = gmosampsGainBefore20060831[gain_key]
                        else:
                            raise Errors.TableKeyError()
                
                # Update the dictionary with the gain value
                ret_gain_dict.update({ext_name_ver:gain})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_gain_dict, name="gain", ad=dataset)
        
        return ret_dv
    
    def gain_setting(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_gain_setting_dict = {}
        
        # If the data have not been prepared, take the raw gain value directly
        # from the appropriate keyword
        if "PREPARED" not in dataset.types:
            
            # Determine the gain keyword from the global keyword dictionary 
            keyword = self.get_descriptor_key("key_gain")
            
            # Get the value of the gain keyword from the header of each pixel
            # data extension as a dictionary where the key of the dictionary is
            # an ("*", EXTVER) tuple
            gain_dict = gmu.get_key_value_dict(adinput=dataset,
                                               keyword=keyword)
            if gain_dict is None:
                # The get_key_value_dict() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            for ext_name_ver, gain in gain_dict.iteritems():
                if gain is None:
                    gain_setting = None
                elif gain > 3.0:
                    gain_setting = "high"
                else:
                    gain_setting = "low"
                
                # Update the dictionary with the gain setting value
                ret_gain_setting_dict.update({ext_name_ver:gain_setting})
        
        else:
            # If the data have been prepared, take the gain setting value
            # directly from the appropriate keyword in the header of each pixel
            # data extension
            #
            # Determine the gain setting keyword from the global keyword
            # dictionary
            keyword = self.get_descriptor_key("key_gain_setting")
            
            # Get the value of the gain setting keyword from the header of each
            # pixel data extension as a dictionary where the key of the
            # dictionary is an ("*", EXTVER) tuple
            gain_setting_dict = gmu.get_key_value_dict(adinput=dataset,
                                                       keyword=keyword)

            if gain_setting_dict is not None:
                ret_gain_setting_dict = gain_setting_dict
            else:
                # The dataset was not processed using gemini_python. Try to get
                # the gain from the "GAINORIG" keyword in the header of each
                # pixel data extension as a dictionary where the key of the
                # dictionary is an ("*", EXTVER) tuple.
                gain_dict = gmu.get_key_value_dict(adinput=dataset,
                                                   keyword="GAINORIG")
                if gain_dict is None:
                    # Resort to getting the gain using the appropriate
                    # descriptor (this will use the updated gain value). Use
                    # as_dict() to return the gain value as a dictionary where
                    # the key of the dictionary is an ("*", EXTVER) tuple
                    # rather than an object. 
                    gain_dv = dataset.gain()
                    if gain_dv.is_none():
                        # The descriptor functions return None if a value
                        # cannot be found and stores the exception
                        # info. Re-raise the exception. It will be dealt with
                        # by the CalculatorInterface. 
                        if hasattr(dataset, "exception_info"):
                            raise dataset.exception_info
                    gain_dict = gain_dv.as_dict()
                    
                else:
                    for ext_name_ver, gain_orig in gain_dict.iteritems():
                        count_exts = 0
                        count_gain_orig = 0
                        if gain_orig == 1:
                            # If the gain from the "GAINORIG" keyword is equal
                            # to 1, the very original gain was written to
                            # "GAINMULT" 
                            count_gain_orig += 1
                        count_exts += 1
                    
                    if count_exts == count_gain_orig:
                        # The value of "GAINORIG" keyword is equal to 1 in all
                        # pixel data extensions. Try to get the gain from the
                        # "GAINMULT" keyword in the header of each pixel data
                        # extension as a dictionary where the key of the
                        # dictionary is an ("*", EXTVER) tuple.
                        gain_dict = gmu.get_key_value_dict(adinput=dataset,
                                                           keyword="GAINMULT")
                if gain_dict is None:
                    # The get_key_value_dict() function returns None if a value
                    # cannot be found and stores the exception info. Re-raise
                    # the exception. It will be dealt with by the
                    # CalculatorInterface.
                    if hasattr(dataset, "exception_info"):
                        raise dataset.exception_info
                
                for ext_name_ver, gain in gain_dict.iteritems():
                    if gain is None:
                        gain_setting = None
                    elif gain > 3.0:
                        gain_setting = "high"
                    else:
                        gain_setting = "low"
                    
                    # Update the dictionary with the gain setting value
                    ret_gain_setting_dict.update({ext_name_ver:gain_setting})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_gain_setting_dict, name="gain_setting",
                                 ad=dataset)
        return ret_dv
    
    def group_id(self, dataset, **args):
        # For GMOS data, the group id contains the detector_x_bin,
        # detector_y_bin and amp_read_area in addition to the observation id.
        # Get the observation id, the binning of the x-axis and y-axis and the
        # amp_read_area values using the appropriate descriptors.
        observation_id_dv = dataset.observation_id()
        detector_x_bin_dv = dataset.detector_x_bin()
        detector_y_bin_dv = dataset.detector_y_bin()
        amp_read_area_dv = dataset.amp_read_area()
        
        if (observation_id_dv.is_none() or detector_x_bin_dv.is_none() or
            detector_x_bin_dv.is_none() or amp_read_area_dv.is_none()):
            # The descriptor functions return None if a value cannot be found
            # and stores the exception info. Re-raise the exception. It will be
            # dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        # Return the amp_read_area as an ordered list
        amp_read_area_list = amp_read_area_dv.as_list()
        
        # For all data other than data with an AstroData type of GMOS_BIAS, the
        # group id contains the filter_name. Also, for data with an AstroData
        # type of GMOS_BIAS and GMOS_IMAGE_FLAT, the group id does not contain
        # the observation id. 
        if "GMOS_BIAS" in dataset.types:
            ret_group_id = "%s_%s_%s" % (
                detector_x_bin_dv, detector_y_bin_dv, amp_read_area_list)
        else:
            # Get the filter name using the appropriate descriptor
            filter_name_dv = dataset.filter_name(pretty=True)
            
            if filter_name_dv.is_none():
                # The descriptor functions return None if a value cannot be
                # found and stores the exception info. Re-raise the exception.
                # It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            if "GMOS_IMAGE_FLAT" in dataset.types:
                ret_group_id = "%s_%s_%s_%s" % (
                    detector_x_bin_dv, detector_y_bin_dv, filter_name_dv,
                    amp_read_area_list)
            else:
                ret_group_id = "%s_%s_%s_%s_%s" % (
                    observation_id_dv, detector_x_bin_dv, detector_y_bin_dv,
                    filter_name_dv, amp_read_area_list)
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_group_id, name="group_id", ad=dataset)
        
        return ret_dv
    
    def nod_count(self, dataset, **args):
        # The number of nod and shuffle cycles can only be obtained from nod
        # and shuffle data
        if "GMOS_NODANDSHUFFLE" in dataset.types:
            
            # Determine the number of nod and shuffle cycles keyword from the
            # global keyword dictionary
            keyword = self.get_descriptor_key("key_nod_count")
            
            # Get the value of the number of nod and shuffle cycles keyword
            # from the header of the PHU
            nod_count = dataset.phu_get_key_value(keyword)
            
            if nod_count is None:
                # The phu_get_key_value() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            # Return the nod count integer
            ret_nod_count = int(nod_count)
        else:
            raise Errors.DescriptorTypeError()
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_nod_count, name="nod_count", ad=dataset)
        
        return ret_dv
    
    def nod_pixels(self, dataset, **args):
        # The number of pixel rows the charge is shuffled by can only be
        # obtained from nod and shuffle data
        if "GMOS_NODANDSHUFFLE" in dataset.types:
            
            # Determine the number of pixel rows the charge is shuffled by
            # keyword from the global keyword dictionary
            keyword = self.get_descriptor_key("key_nod_pixels")
            
            # Get the value of the number of pixel rows the charge is shuffled
            # by keyword from the header of the PHU
            nod_pixels = dataset.phu_get_key_value(keyword)
            
            if nod_pixels is None:
                # The phu_get_key_value() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            # Return the nod pixels integer
            ret_nod_pixels = int(nod_pixels)
        else:
            raise Errors.DescriptorTypeError()
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_nod_pixels, name="nod_pixels", ad=dataset)
        
        return ret_dv
    
    def nominal_photometric_zeropoint(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_nominal_photometric_zeropoint_dict = {}
        
        # Get the lookup table containing the nominal zeropoints
        table = Lookups.get_lookup_table("Gemini/GMOS/Nominal_Zeropoints",
                                         "nominal_zeropoints")
        
        # Get the values of the gain, detector name and filter name using the
        # appropriate descriptors.
        gain_dv = dataset.gain()
        array_name_dv = dataset.array_name()
        filter_name_dv = dataset.filter_name(pretty=True)
        
        if (gain_dv.is_none() or array_name_dv.is_none() or
            filter_name_dv.is_none()):
            # The descriptor functions return None if a value cannot be found
            # and stores the exception info. Re-raise the exception. It will be
            # dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        # Get the value of the BUNIT keyword from the header of each pixel data
        # extension as a dictionary where the key of the dictionary is an
        # ("*", EXTVER) tuple 
        bunit_dict = gmu.get_key_value_dict(adinput=dataset, keyword="BUNIT")

        # Loop over extvers in one of the DVs as they will all have the same
        # extvers. Then use get_value method to get the value from a given DV
        # for a particular extver
        for extver in gain_dv.ext_vers():
            # key used for bunit dict and for instantiating outuput DV
            ext_name_ver = ("*", extver)
            
            # Obtain the values from the three DVs for the current extver
            array_name = array_name_dv.get_value(extver)
            gain = gain_dv.get_value(extver)
            filter_name = filter_name_dv.get_value(extver)

            if array_name is None:
                nominal_photometric_zeropoint = None
            else:
                # Determine whether data are in ADU or electrons
                if bunit_dict is not None:
                    bunit = bunit_dict[ext_name_ver]
                else:
                    bunit = None
                
                # If bunit is "electron" or None, set the gain factor to 0.0 
                gain_factor = 0.0
                if bunit == "adu":
                    gain_factor = 2.5 * math.log10(gain)

                nominal_zeropoint_key = (array_name, filter_name)
                
                if nominal_zeropoint_key in table:
                    nominal_photometric_zeropoint = (
                        table[nominal_zeropoint_key] - gain_factor)
                else:
                    raise Errors.TableKeyError()
            
            # Update the dictionary with the nominal photometric zeropoint
            # value 
            ret_nominal_photometric_zeropoint_dict.update({
                ext_name_ver: nominal_photometric_zeropoint})

        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_nominal_photometric_zeropoint_dict,
                                 name="nominal_photometric_zeropoint",
                                 ad=dataset)
        return ret_dv
    
    def non_linear_level(self, dataset, **args):
        # Set the non linear level equal to the saturation level for GMOS
        non_linear_level_dv = dataset.saturation_level()
        
        if non_linear_level_dv.is_none():
            # The descriptor functions return None if a value cannot be found
            # and stores the exception info. Re-raise the exception. It will be
            # dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        ret_non_linear_level = non_linear_level_dv
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_non_linear_level, name="non_linear_level",
                                 ad=dataset)
        return ret_dv
    
    def overscan_section(self, dataset, pretty=False, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_overscan_section_dict = {}
        
        # Determine the overscan section keyword from the global keyword
        # dictionary
        keyword = self.get_descriptor_key("key_overscan_section")
        
        # Get the value of the overscan section keyword from the header of each
        # pixel data extension as a dictionary where the key of the dictionary
        # is an ("*", EXTVER) tuple
        overscan_section_dict = gmu.get_key_value_dict(adinput=dataset,
                                                       keyword=keyword)
        if overscan_section_dict is None:
            # The get_key_value_dict() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        dict = overscan_section_dict.iteritems()
        for ext_name_ver, raw_overscan_section in dict:
            if raw_overscan_section is None:
                overscan_section = None
            elif pretty:
                # Use the overscan section string that uses 1-based indexing as
                # the value in the form [x1:x2,y1:y2]
                overscan_section = str(raw_overscan_section)
            else:
                # Use the overscan section list that uses 0-based,
                # non-inclusive indexing as the value in the form
                # [x1, x2, y1, y2] 
                overscan_section = gmu.sectionStrToIntList(
                    raw_overscan_section)
            
            # Update the dictionary with the array section value
            ret_overscan_section_dict.update({ext_name_ver:overscan_section})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_overscan_section_dict,
                                 name="overscan_section", ad=dataset)
        return ret_dv
    
    def pixel_scale(self, dataset, **args):
        # Get the lookup table containing the pixel scale values
        gmosPixelScales = Lookups.get_lookup_table(
            "Gemini/GMOS/GMOSPixelScale", "gmosPixelScales")
        
        # Get the values of the instrument and the binning of the y-axis using
        # the appropriate descriptors
        instrument_dv = dataset.instrument()
        detector_y_bin_dv = dataset.detector_y_bin()
        
        # Determine the detector type keyword from the global keyword
        # dictionary
        keyword = self.get_descriptor_key("key_detector_type")
        
        # Get the value for the detector type keyword from the header of the
        # PHU
        detector_type = dataset.phu_get_key_value(keyword)
        
        if (instrument_dv.is_none() or detector_y_bin_dv.is_none() or
            detector_type is None):
            # The descriptor functions return None if a value cannot be found
            # and stores the exception info. Re-raise the exception. It will be
            # dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        # Use as_pytype() to return the instrument value as the default python
        # type, rather than an object
        instrument = instrument_dv.as_pytype()
        
        # Get the unbinned pixel scale (in arcseconds per unbinned pixel) from
        # the lookup table 
        pixel_scale_key = (instrument, detector_type)
        
        if pixel_scale_key in gmosPixelScales:
            raw_pixel_scale = gmosPixelScales[pixel_scale_key]
        else:
            raise Errors.TableKeyError()
        
        # Return the binned pixel scale value
        ret_pixel_scale = float(detector_y_bin_dv * raw_pixel_scale)
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_pixel_scale, name="pixel_scale",
                                 ad=dataset)
        return ret_dv
    
    def read_mode(self, dataset, **args):
        # For GMOS data, raise an exception if the read_mode descriptor called,
        # since it is not relevant for GMOS data.
        raise Errors.ExistError()
    
    def read_noise(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_read_noise_dict = {}
        
        # If the data have been prepared, take the read noise value directly
        # from the appropriate keyword. At some point, a check for the STDHDRSI
        # header keyword should be added, since the function that overwrites
        # the read noise keyword also writes the STDHDRSI keyword.
        if "PREPARED" in dataset.types:
            
            # Determine the read noise keyword from the global keyword
            # dictionary
            keyword = self.get_descriptor_key("key_read_noise")
            
            # Get the value of the read noise keyword from the header of each
            # pixel data extension as a dictionary where the key of the
            # dictionary is an ("*", EXTVER) tuple
            read_noise_dict = gmu.get_key_value_dict(adinput=dataset,
                                                     keyword=keyword)
            if read_noise_dict is None:
                # The get_key_value_dict() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            for ext_name_ver, raw_read_noise in read_noise_dict.iteritems():
                if raw_read_noise is None:
                    read_noise = None
                else:
                    read_noise = float(raw_read_noise)
                
                # Update the dictionary with the read noise value
                ret_read_noise_dict.update({ext_name_ver:read_noise})
        else:
            
            # Get the lookup table containing the read noise values by
            # amplifier
            gmosampsRdnoise, gmosampsRdnoiseBefore20060831 = (
                Lookups.get_lookup_table("Gemini/GMOS/GMOSAmpTables",
                                         "gmosampsRdnoise",
                                         "gmosampsRdnoiseBefore20060831"))
            
            # Get the UT date, gain setting and read speed setting values using
            # the appropriate descriptors
            ut_date_dv = dataset.ut_date()
            gain_setting_dv = dataset.gain_setting()
            read_speed_setting_dv = dataset.read_speed_setting()
            
            if (ut_date_dv.is_none() or gain_setting_dv.is_none() or
                read_speed_setting_dv.is_none()):
                # The descriptor functions return None if a value cannot be
                # found and stores the exception info. Re-raise the exception.
                # It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            # Use as_dict() and as_pytype() to return the values as a
            # dictionary and the default python type, respectively, rather than
            # an object
            ut_date = str(ut_date_dv)
            gain_setting_dict = gain_setting_dv.as_dict()
            read_speed_setting = read_speed_setting_dv.as_pytype()
            
            obs_ut_date = datetime(*strptime(ut_date, "%Y-%m-%d")[0:6])
            old_ut_date = datetime(2006, 8, 31, 0, 0)
            
            # Determine the name of the detector amplifier keyword (ampname)
            # from the global keyword dictionary 
            keyword = self.get_descriptor_key("key_ampname")
            
            # Get the value of the name of the detector amplifier keyword from
            # the header of each pixel data extension as a dictionary where the
            # key of the dictionary is an ("*", EXTVER) tuple
            ampname_dict = gmu.get_key_value_dict(adinput=dataset,
                                                  keyword=keyword)
            if ampname_dict is None:
                # The get_key_value_dict() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info
            
            for ext_name_ver, ampname in ampname_dict.iteritems():
                gain_setting = gain_setting_dict[ext_name_ver]
                
                if ampname is None or gain_setting is None:
                    read_noise = None
                else:
                    read_noise_key = (
                        read_speed_setting, gain_setting, ampname)
                    
                    if obs_ut_date > old_ut_date:
                        if read_noise_key in gmosampsRdnoise:
                            read_noise = gmosampsRdnoise[read_noise_key]
                        else:
                            raise Errors.TableKeyError()
                    else:
                        if read_noise_key in gmosampsRdnoiseBefore20060831:
                            read_noise = gmosampsRdnoiseBefore20060831[
                                read_noise_key]
                        else:
                            raise Errors.TableKeyError()
                
                # Update the dictionary with the read noise value
                ret_read_noise_dict.update({ext_name_ver:read_noise})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_read_noise_dict, name="read_noise",
                                 ad=dataset)
        return ret_dv
    
    def read_speed_setting(self, dataset, **args):
        # Determine the amplifier integration time keyword (ampinteg) from the
        # global keyword dictionary
        keyword = self.get_descriptor_key("key_ampinteg")
        
        # Get the amplifier integration time from the header of the PHU
        ampinteg = dataset.phu_get_key_value(keyword)
        
        if ampinteg is None:
            # The phu_get_key_value() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        if ampinteg == 1000:
            ret_read_speed_setting = "fast"
        else:
            ret_read_speed_setting = "slow"
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_read_speed_setting,
                                 name="read_speed_setting", ad=dataset)
        return ret_dv
    
    def saturation_level(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_saturation_level_dict = {}
        
        # Get the lookup table containing the saturation values by amplifier
        gmosThresholds = Lookups.get_lookup_table(
            "Gemini/GMOS/GMOSThresholdValues", "gmosThresholds")
        
        # The hard limit for saturation is the controller digitization limit
        controller_limit = 65535
        
        # Determine the name of the bias image and the name of the dark image
        # keywords from the global keyword dictionary 
        keyword1 = self.get_descriptor_key("key_bias_image")
        keyword2 = self.get_descriptor_key("key_dark_image")
        
        # Get the value of the the name of the bias image and the name of the
        # dark image keywords from the header of the PHU 
        bias_image = dataset.phu_get_key_value(keyword1)
        dark_image = dataset.phu_get_key_value(keyword2)
        
        # Determine the name of the detector amplifier (ampname) and the
        # overscan value keywords from the global keyword dictionary
        keyword3 = self.get_descriptor_key("key_ampname")
        keyword4 = self.get_descriptor_key("key_overscan_value")
        
        # Get the value of the name of the detector amplifier, the overscan
        # value and the BUNIT keywords from the header of each pixel data
        # extension as a dictionary, where the key of the dictionary is an
        # ("*", EXTVER) tuple
        keyword_dict = gmu.get_key_value_dict(
            adinput=dataset, keyword=[keyword3, keyword4, "BUNIT"],
            dict_key_extver=True)
        
        ampname_dict = keyword_dict[keyword3]
        
        if ampname_dict is None:
            # The get_key_value_dict() function returns None if a value cannot
            # be found and stores the exception info. Re-raise the exception.
            # It will be dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        overscan_dict = keyword_dict[keyword4]
        bunit_dict = keyword_dict["BUNIT"]
        
        # Get the name of the detector, the gain and the binning of the x-axis
        # and y-axis values using the appropriate descriptors
        detector_name_dv = dataset.detector_name(pretty=True)
        gain_dv = dataset.gain()
        detector_x_bin_dv = dataset.detector_x_bin()
        detector_y_bin_dv = dataset.detector_y_bin()
        
        if (detector_name_dv.is_none() or gain_dv.is_none() or
            detector_x_bin_dv.is_none() or detector_y_bin_dv.is_none()):
            # The descriptor functions return None if a value cannot be found
            # and stores the exception info. Re-raise the exception. It will be
            # dealt with by the CalculatorInterface.
            if hasattr(dataset, "exception_info"):
                raise dataset.exception_info
        
        # Determine the bin factor. If the bin factor is great than 2, the
        # saturation level will be equal to the controller digitization limit.
        bin_factor = detector_x_bin_dv * detector_y_bin_dv
        
        for extver, ampname in ampname_dict.iteritems():
            gain = gain_dv.get_value(extver=extver)
            
            # Determine whether it is required to calculate the bias level
            # (bias level calculations can take some time)
            overscan = None
            if overscan_dict is not None:
                overscan = overscan_dict[extver]
            if overscan is not None:
                # The overscan was subtracted from the data
                data_contains_bias = False
            elif bias_image is not None or dark_image is not None:
                # The overscan was not subtracted from the data, but the bias
                # or dark was subtracted from the data
                data_contains_bias = False
            else:
                # The data still contains a bias level
                data_contains_bias = True
            
            if ((not data_contains_bias) or
                (not detector_name_dv == "EEV" and data_contains_bias and
                 bin_factor <= 2)):
                # Calculate the bias level
                bias_level = gdc.get_bias_level(adinput=dataset, estimate=True)
            
            # Correct the controller limit for bias level and units
            processed_limit = controller_limit
            if not data_contains_bias:
                processed_limit -= bias_level[extver]
            
            # Check units of data (i.e., ADU vs. electrons)
            bunit = None
            if bunit_dict is not None:
                bunit = bunit_dict[extver]
            if bunit == "electron" or bunit == "electrons":
                processed_limit *= gain
            
            if detector_name_dv == "EEV" or bin_factor > 2:
                # For old EEV CCDs, use the detector limit
                saturation = processed_limit
            else:
                # For new GMOS-N CCDs, look up the saturation limit from
                # the table, then correct for binning, units, and bias level
                
                # Get the base saturation value from the lookup table
                if ampname in gmosThresholds:
                    # saturation value is an integer with units ADU
                    saturation = gmosThresholds[ampname]
                else:
                    # This error condition will be hit for all mosaicked
                    # or tiled data
                    raise Errors.TableKeyError()
                
                # Correct the saturation level for binning
                saturation = saturation * bin_factor
                
                # The saturation level is reported in electrons; convert
                # it to ADU if needed
                if bunit != "electron" or bunit != "electrons":
                    saturation = saturation / gain
                
                # The saturation level does not contain the bias; add it
                # in if necessary
                if data_contains_bias:
                    saturation += bias_level[extver]
                
                # Check whether the value is now over the controller limit;
                # if so, set it to the hard limit
                if saturation > processed_limit:
                    saturation = processed_limit
            
            ret_saturation_level_dict.update({extver: saturation})
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_saturation_level_dict,
                                 name="saturation_level", ad=dataset)
        return ret_dv
    
    def wavelength_band(self, dataset, **args):
        if "IMAGE" in dataset.types:
            # If imaging, associate the filter name with a central wavelength
            filter_table = Lookups.get_lookup_table(
                "Gemini/GMOS/GMOSFilterWavelength", "filter_wavelength")
            filter = str(dataset.filter_name(pretty=True))
            if filter in filter_table:
                ctrl_wave = filter_table[filter]
            else:
                raise Errors.TableKeyError()
        else:
            ctrl_wave = dataset.central_wavelength(asMicrometers=True)
        
        min_diff = None
        band = None
        
        for std_band, std_wave in self.std_wavelength_band.items():
            diff = abs(std_wave - ctrl_wave)
            if min_diff is None or diff < min_diff:
                min_diff = diff
                band = std_band
        
        if band is None:
            raise Errors.CalcError()
        
        else:
            ret_wavelength_band = band
        
        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_wavelength_band, name="wavelength_band",
                                 ad=dataset)
        return ret_dv
