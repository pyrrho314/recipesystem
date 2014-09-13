import sys, os
from astrodata import AstroData
from astrodata import DescriptorUnits as Units
from astrodata.DescriptorUnits import iunits
import GemCalcUtil                

## Intergrating the Units module test



ad = AstroData("../../../test_data/gmosspect/N20020810S0117.fits")


print "\n============= CURRENT setup using gemCalcUtil.py in DV =============\n"
centwave = ad.central_wavelength()
print("CASE 1: Get central_wavelength as Meter (DEFAULT):")
print("    centwave = ad.central_wavelength()")
print("    centwave = %s" % centwave.as_str())

print("\nCASE 2: Get central_wavelength as Micrometer:")
print("    um_centwave = ad.central_wavelength(asMicrometers=True)")
umcentwave = ad.central_wavelength(asMicrometers=True)
print("    um_centwave = %s" % umcentwave.as_str()) 

print("\nCASE 3: Get central_wavelength as Nanometer:")
print("    nm_centwave = ad.central_wavelength(asNanometers=True)")
nmcentwave = ad.central_wavelength(asNanometers=True)
print("    nm_centwave = %s" % nmcentwave.as_str()) 

print("\nCASE 4: Get central_wavelength as Angstroms:")
print("    ang_centwave = ad.central_wavelength(asAngstroms=True)")
angcentwave = ad.central_wavelength(asAngstroms=True)
print("    ang_centwave = %s\n" % angcentwave.as_str()) 
print("centwave.info()")
centwave.info()

print "="*79
print "\n\n","****** Introducing in house DV Unit Conversion: ******\n"
print "First change the 'unit' member var for centwave, this should be "
print "     done in the calculator interface, like how pytype is assigned."
print "     Here we show a simple example by using centwave from the ex. above."
print("\ncentwave.unit = Units.m (look under .unit below now, see the 'm')")
centwave.unit = Units.m
print("centwave.info()")
print("id(centwave) = %d" % id(centwave))
centwave.info()
print("CASE 1: Get central_wavelength as Meter (DEFAULT):")
print("    centwave = ad.central_wavelength()")
print("    centwave = %s" % centwave.as_str())

print("\nCASE 2: Get central_wavelength as Micrometer:")
print("    newDV = centwave.convert_value_to(Units.um)")
newDV = centwave.convert_value_to(Units.um)
print("    newDV = %s" % newDV.as_str())
print("...The centwave DV remains the same...")
print("centwave.info()")
print("id(centwave) = %d" % id(centwave))
centwave.info()
print("...But newDV takes on the converted value...")
print("newDV.info()")
print("id(newDV) = %d" % id(newDV))
newDV.info()
# direct conversion works as well
#print("    um_centwave = %s" % centwave.convert_value_to(Units.um))

print("\nCASE 3: Get central_wavelength as Nanometer:")
print("    nm_centwave = centwave.convert_value_to(Units.nm)")
nm_centwave = centwave.convert_value_to(Units.nm)
print("    nm_centwave = %s" % nm_centwave) 

print("\nCASE 4: Get central_wavelength as Angstroms:")
print("    ang_centwave = centwave.convert_value_to(Units.angstrom)")
ang_centwave = centwave.convert_value_to(Units.angstrom)
print("    ang_centwave = %s\n" % ang_centwave) 

#print("Aliases and nested conversion dict in trunk/astrodata/DescriptorUnits.py)")
#print("iunits = %s (%s)" % (iunits, type(iunits)))
#keys = iunits.keys()
#keys.sort()
#for key in keys:
#    print("\t%-15s: %s" % (key, iunits[key]))


