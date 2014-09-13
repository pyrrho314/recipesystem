import DescriptorUnits as Units
from DescriptorUnits import Unit


import inspect
# functions
def whoami():
    print repr(inspect.stack())
    return inspect.stack()[1][3]
    
def whocalledme():
    return inspect.stack()[2][3]
   
          
class ValueWrapper(object):
    val = None
    descUnits = None
    origval = None
    origUnits = None
    pytype = float
    def __init__(self, val, unit = None, pytype = None):
    
        self.val = val
        self.origval = val
        if unit != None:
            self.descUnits = unit
        else:
            import inspect 
            st = inspect.stack()
            callername = st[1][3]
            callerframe = inspect.stack()[1][0]
            fargs = inspect.getargvalues(callerframe)
            callercalc = fargs[3]["self"]
            unit = eval("callercalc.%s.descUnits" % callername)
        if pytype:
            self.pytype = pytype
        else:
            try:
                pytype = eval("callercalc.%s.pytype" % callername)
            except:
                pytype = float
                
            self.pytype = pytype
        self.descUnits = unit
    
    def convert_value_to(self, newUnits, newType = None):
        retval = self.descUnits.convert(self.val, newUnits)
        if newType:
            retval = newType(retval)
        return retval
    

    def overloaded(self, other):
        other = self.pytype(other)
        funcname = whocalledme()
        if hasattr(other, funcname):
            othertype = type(other)
            return eval("other.%s(othertype(self))" % funcname)
        
    # operands 
    def __add__(self, other):
        return self.overloaded(other)
    def __sub__(self, other):
        return self.overloaded(other)
    def __mul__(self, other):
        return self.overloaded(other)
    def __div__(self, other):
        return self.overloaded(other)
    def __truediv__(self, other):
        return self.overloaded(other)
    def __floordiv__(self, other):
        return self.overloaded(other)
    def __mod__(self, other):
        return self.overloaded(other)
    def __divmod__(self, other):
        return self.overloaded(other)
    def __pow__(self, other):
        return self.overloaded(other)
    def __lshift__(self, other):
        return self.overloaded(other)
    def __rshift__(self, other):
        return self.overloaded(other)
    def __and__(self, other):
        return self.overloaded(other)
    def __xor__(self, other):
        return self.overloaded(other)
    def __or__(self, other):
        return self.overloaded(other)
    # reflected operands
    def __radd__(self, other):
        return self.overloaded(other)
    def __rsub__(self, other):
        return self.overloaded(other)
    def __rmul__(self, other):
        return self.overloaded(other)
    def __rdiv__(self, other):
        return self.overloaded(other)
    def __rtruediv__(self, other):
        return self.overloaded(other)
    def __rfloordiv__(self, other):
        return self.overloaded(other)
    def __rmod__(self, other):
        return self.overloaded(other)
    def __rdivmod__(self, other):
        return self.overloaded(other)
    def __rpow__(self, other):
        return self.overloaded(other)
    def __rlshift__(self, other):
        return self.overloaded(other)
    def __rrshift__(self, other):
        return self.overloaded(other)
    def __rand__(self, other):
        return self.overloaded(other)
    def __rxor__(self, other):
        return self.overloaded(other)
    def __ror__(self, other):
        return self.overloaded(other)
       
    def __str__(self):
        return str(self.val)

    def __int__(self):
        return int(self.val)
        
    def __float__(self):
        return float(self.val)


# descriptor functions
class Calculator(object):

    def eggs(self, color = None):
        if color == "blue":
            return ValueWrapper(640, unit = Units.mile)
        elif color == "green":
            return ValueWrapper(750, unit = Units.nanometers)
        elif color == "red":
            return ValueWrapper(6850)
            return 6850
        else:
            return None
    eggs.descUnits = Units.mile
    eggs.pytype = float

# create the calculator
c = Calculator()
value = c.eggs("blue")

# tests
print "#"*30
print "DESCIPTOR VALUE UNIT CLASS TEST"
print "#"*30
print "\n\nFILTER RETURNS"
print "-"*25
val2 = c.eggs("blue")
print "blue", val2, val2.descUnits, str(type(val2.descUnits))
val2 = c.eggs("green")
print "green", val2, val2.descUnits
val3 = c.eggs("red")
print "red", val3, val3.descUnits

print "\n\nCAST / LEFT RIGHT OPERANDS/ TYPE INHERITANCE"
print "-"*25
print "int(val2) = ", int(val2)
print "val2 + 1 = ", val2+1
print "1 + val2 = ", 1+val2
print "1.0 + val2 = ", 1.0+val2, " (should be float)"

print "\n\nCONVERT UNITS"
print "-"*25
print "val2 = ",val2, val2.descUnits.name
print "val2.convert_value_to(Units.m) =", val2.convert_value_to(Units.m) 
print "val2.convert_value_to(Units.microns) = ",val2.convert_value_to(Units.microns)
print "val2.convert_value_to(Units.nm) = ",val2.convert_value_to(Units.nm)
print "val2.convert_value_to(Units.angstroms) = ",val2.convert_value_to(Units.angstroms)



#print "'hello'+val2"
#print "1",val2




