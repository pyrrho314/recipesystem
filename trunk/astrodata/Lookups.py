import os
import ConfigSpace
try:
    import pyfits
except:
    pass
    
def get_lookup_table(modname, *lookup, **args):
    """
        get_lookup_table() is used to get lookup table style sets of variables
        from a common facility, allowing the storage in common (global) space
        so that multiple scripts can refer to one lookup table
        without having to manage where this table is stored.  E.g. the Calculator
        (see L{Descriptors}) for NIRI data requires a NIRI lookup table that
        other parts of the package, unrelated to Descriptors, also need to 
        access.  This facility saves these separate components from knowing
        where the configuration is actually stored, or even that other
        parts of the system are relying on it, and ensure that changes will
        affect every part of the system.
        
    @param modname: namespace specifier for the table... in default case this
        is the directory and file name of the module in which the lookup
        table is stored, and the file is pure python.  However, the Lookups
        module can redirect this, using the modname, for example, as a
        key to find the lookup table in a database or elsewhere. Nothing like
        the latter is done at this time, and what is loaded are pure python
        files (e.g. a dict definition) from disk.
    @type modname: string
    @param lookup: name of the lookup table to load
    @type lookup: string
    """
    retval = None
    context = None
    if "context" in args:
        context = args["context"]
        if not context:
            context = "default"
    else:
        context = "default"
    #if not context:
    #    modname = ConfigSpace.lookup_path(modname)
    #else:
    modname = ConfigSpace.lookup_context_path(modname, context=context)
    if not modname or ( not os.path.exists(modname) ):
        return None
    if ".py" in modname:
        f = file(modname)
        g = {}
        l = {}
        
        exec(f,g,l)
        
        #print "L38:",l.keys(),l
        
        f.close()

        if len(lookup) == 1:
            retval = l[lookup[0]]
        elif len(lookup) == 0:
            retval = []
            for key in l:
                retval.append(l[key])
        else:
            retval = []
            for item in lookup:
                retval.append(eval(item))
    elif ".fits" in modname:
        # in this case lookup will have extension ids
        table = pyfits.open(modname)
        if len(lookup) == 1:
            retval = table[lookup[0]]
        else:
            retval = []
            for item in lookup:
                retval.append(table[item])
            
    else:
        raise "this should never happen, tell someone"
    return retval
    
