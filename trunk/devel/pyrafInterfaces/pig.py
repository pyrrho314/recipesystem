#module header template
mod_header_template="""
from piro import mkRO
"""
func_template="""
def %(primname)s(*args, **argv):
    ro = mkRO(astrotype="GEMINI",args=args,argv=argv)
    ro.runstep("%(primname)s", ro.context)
"""
primsfile = "prims_GEMINI.py"

f = open(primsfile)
primstr = f.read()
f.close()
prims = eval(primstr)
primsblock = ""
for prim in prims:
    primsblock += func_template % {"primname":prim}

primslist = 'prims = %s\n' % repr(prims)
# OUTPUT SECTION OF SCRIPT
outfilename = "pyraf_"+primsfile
outfile = open(outfilename, "w")
outfile.write(mod_header_template)

outfile.write(primslist)
outfile.write(primsblock)
outfile.close()

