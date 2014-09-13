class GeneralClassificationLibrary(object):
    _known_types = None
    def __init__(self):
        self._known_types = {}
        
    def get_type_obj(self, typ):
        if (typ in self._known_types):
            gco = self._known_types[typ]
        else:
            gco = GeneralClassificationObj(typ)
            self._known_types[typ] = gco
        
        return gco

    def type_is_child_of(self, childtyp, parenttyp):
        if childtyp in self._known_types:
            ctyp = self._known_types[childtyp]
            return ctyp.child_of(parenttyp)
        else:
            return False
        

class GeneralClassificationObj(object):
    typename = None
    _parent = None
    def __init__(self, typ = None):
        self.typename = typ
        
    def child_of(self, parenttyp):
        return False # parantage not yet supported in GeneralType
    
    def get_parent(self):
        return self._parent
    
    def set_parent(self, val):
        self._parent = val
    
    parent = property(get_parent, set_parent)

_gco = GeneralClassificationObj();
