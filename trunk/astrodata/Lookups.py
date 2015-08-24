import os
import ConfigSpace
from copy import deepcopy

try:
    import pyfits
except:
    pass

def get_lookup_value(modname, *lookup, **args):
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
            context = None
    else:
        context = None
    
    if context == None:
        ConfigSpace.get_current_default_context()
        
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
            if lookup[0] in l:
                retval = l[lookup[0]]
            else:
                retval =  None
        elif len(lookup) == 0:
            retval = []
            for key in l:
                retval.append(l[key])
        else:
            retval = []
            for item in lookup:
                if item in l:
                    retval.append(l[item])
                else:
                    retval.append(None)
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

    if context == None:
        context = ConfigSpace.get_current_default_context()
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
            if lookup[0] in l:
                retval = l[lookup[0]]
            else:
                retval =  [None]
        elif len(lookup) == 0:
            retval = []
            for key in l:
                retval.append(l[key])
        else:
            retval = []
            for item in lookup:
                if item in l:
                    retval.append(l[item])
                else:
                    retval.append(None)
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

def compose_multi_table(lookaddr, *lookups, **args):
    """ Returns a dictionary keyed by the lookups name of a composed combination
    of tables distributed throughout loaded kits.
    """
    retdict = {}
    context = None
    if "context" in args:
        context = args["context"]
    if context == None:
        context = ConfigSpace.get_current_default_context()

    paths = ConfigSpace.lookup_multi_paths(lookaddr, context = context)
    
    if len(paths) == 0:
        return None
    paths.reverse() # so early files override latter, as it's in ADCONFIG path order
    for modname in paths:
        #print "L185: modname = %s" % modname
        f = file(modname)
        g = {}
        l = {}
        
        exec(f,g,l)
        
        #print "L191:",l.keys(),l
        
        f.close()
        
        i = 0
        contributed = False
        for lookup in lookups:
            #print "L197: #%d - lookup = %s" % (i,lookup); i+= 1
            if lookup in l:
                lval = l[lookup]
                valtype = type(lval)
                if valtype == dict:
                    if not lookup in retdict:
                        retdict[lookup] = {}
                    curval = retdict[lookup]
                    curval.update(l[lookup])
                elif valtype == list:
                    if not lookup in retdict:
                        retdict[lookup] = []
                    curval = retdict[lookup]
                    lval.extend(curval)
                    #print "L211 lval = %s" % lval
                    curval = lval
                else:
                    if not lookup in retdict:
                        retdict[lookup] = []
                    curval = retdict[lookup]
                    curval.insert(0,lval)
                retdict[lookup] = curval
                contributed = True
        if not ("_contributors" in retdict):
            retdict["_contributors"] = []
        contribs = retdict["_contributors"]
        if contributed:
            contribs.append(modname)
    return retdict 
   
   
            
        
    
