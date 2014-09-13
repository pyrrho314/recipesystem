from astrodata import AstroData

ad = AstroData("../../../test_data/recipedata/N20091027S0141.fits")

# this is a developer script, the way it works is you can set a small bank of variables
# to various descriptor return values (ones that return DescriptorValues!)
# and a series of automated printouts can help you inspect the value and probe its
# behavior

ad.descriptorFormat = "as_dict"
g = ad.detector_x_bin()
h = ad.detector_x_bin()
i = ad.detector_x_bin()
j = ad.detector_x_bin(format="as_dict")
dvsnames = ["g","h","i","j"]

dvs = []
for dvname in dvsnames:
    dv = eval(dvname)
    dv.myname = dvname
    dvs.append(dv)   

for dv in dvs:
    print 'dv "%s" (%s)' % ( dv.myname,dv.name)
    print dv.info()
    
import math

print "using DescriptorValue as a float"
a = math.sin(h)
print "math.sin(h) =",a
a = math.sin(2.0)
print "math.sin(2.0) =",a

print "^"*80

for dv in dvs:
    print "="*20 + dv.myname + "="*20
    print dv

print "()"*80
ad.descriptorFormat = "value"
g = ad.detector_x_bin()
h = ad.gain()
i = ad.gain(format="as_dict")
j = ad.amp_read_area(format="db")
dvsnames = ["g","h","i","j"]
dvs = []
for dvname in dvsnames:
    dv = eval(dvname)
    dv.myname = dvname
    dvs.append(dv)   
    


for dv in dvs:
    print "="*20 + dv.myname + "(%s)" % dv.name + "="*20
    print dv


