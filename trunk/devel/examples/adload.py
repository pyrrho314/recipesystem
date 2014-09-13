from astrodata import AstroData

ad = AstroData("/home/callen/SVN-AD/gemini_python/test_data/recipedata/N20091027S0141.fits")

try:
    filt = ad.filter_name()
except:
    pass
    
