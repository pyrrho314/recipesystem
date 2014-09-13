from astrodata import AstroData

gmos = AstroData("../../../test_data/gmosspect/N20020810S0117.fits")
niri = AstroData("../../../test_data/niri/spectra/N20100104S0385.fits")
nifs = AstroData("../../../test_data/recipedata/stackable/set001/N20090902S0079.fits")

print gmos.central_wavelength()
print nifs.central_wavelength()
print gmos.central_wavelength()
