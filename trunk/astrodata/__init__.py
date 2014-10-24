# __init__.py for astrodata
# $Rev: 4361 $
__version__ = '0.9.0'

# The following pyfits import previously clobbered astrodata.__version__
# >>> import astrodata
# >>> astrodata.__version__
# '3.1.2'
#
# That is the pyfits version. Below fixes this.

#from pyfits import __version__ as pf_version
#new_pyfits_version = [int(m) for m in pf_version.split('.')[:2]] >= [3,1]
# only works with new_pyfits
new_pyfits_version = True

import data
import datatypes
#import AstroDataType
try:
    import CalculatorInterface
except ImportError:
    pass
try:
    import FITS_Keywords
except ImportError:
    pass
    
import Calculator
import GeminiData
# import RecipeManager
import DataSpider
import Structures
import ReductionContextRecords
import ReductionObjectRequests
# import StackKeeper
import LocalCalibrationService
import gdpgutil
#import tkMonitor
import IDFactory
import LocalCalibrationService
import Errors

import AstroData as astrodata_module
from data import AstroData

# define the constants
from gemconstants import *
