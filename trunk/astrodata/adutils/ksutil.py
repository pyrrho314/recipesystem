from astrodata import termcolor as tc
from pprint import pprint,pformat

def calc_fulltab(indent):
    tabspc = "-"*4
    fulltab = tabspc*indent
    return fulltab

def dict2pretty(name, var, indent=0, namewidth = None, complete = False, say_type = None):
    #retstr = pformat(var)
    retstr = ""
    fulltab = calc_fulltab(indent)
    tabspc  = calc_fulltab(1)
    _o_namewidth = namewidth
    if not namewidth:
        namewidth = len(str(name))
    if isinstance(var, dict):
        retstr += "\n%(indent)s%(key)s %(type)s:" % {
                                                    "indent":fulltab,
                                                    "key": tc.colored(name, attrs=["bold"]),
                                                    "type": tc.colored(repr(type(var)), attrs=["dark"]),
                                                    "extra": tabspc
                                                 }
        sub_namewidth = maxkeylen(var)
        #print "ks19: sub_nw=", sub_namewidth
        if len(var) == 0:
            retstr += "\n%(indent)s%(tab)s:::empty:::" % {"indent":fulltab,
                                                          "tab":tabspc
                                                          }
        keys = var.keys()
        keys.sort()
        for key in keys:
            value = var[key]
            #print key,value
            newstr = dict2pretty(key, value, indent+1, namewidth = sub_namewidth )
            #print "ks28: indent =", indent, namewidth, _o_namewidth
            #print "ks31: newstr", newstr
            retstr += newstr
    elif isinstance(var, list):
         retstr += "\n%(indent)s%(key)s %(type)s:" % { "indent":fulltab,
                                                    "key": tc.colored(name, attrs=["bold"]),
                                                    "type": tc.colored(repr(type(var)), attrs=["dark"]),
                                               }
         
         listlen = len(var)
         
         if len(var) < 50:
             allstr = True
             reprline = []
             for v in var:
                reprline.append( repr(v))
                if not isinstance(v, basestring):
                    allstr = False
             if allstr:
                oneline = ", ".join(var)
             else:
                oneline = ", ".join(reprline)
             if len(oneline)<120:
                return dict2pretty(name, oneline, indent, say_type = type(var), namewidth = namewidth)
                   
         if listlen > 10 and not complete == True:
            last = listlen - 1
            mid = int(last/2);
            retstr += dict2pretty("[0]", var[0], indent+1, namewidth = namewidth) 
            retstr += dict2pretty("[%d]"%mid, var[mid],indent+1, namewidth = namewidth)
            retstr += dict2pretty("[%d]"%last, var[last], indent+1, namewidth = namewidth)  
         else:
            for i in range(0, listlen):
                key = "[%d]"%i
                value = var[i];
                if hasattr(value, "pretty_string"):
                    retstr += tc.colored("\n[%d]" % i, attrs=["bold"])
                    retstr += value.pretty_string(start_indent = indent+1)
                else:     
                    retstr += dict2pretty(key, value, indent+1, namewidth = namewidth)
    else:
        if say_type:
            stype = repr(say_type)
        else:
            stype = repr(type(var))
            
        if isinstance(var, basestring):
            pvar = var.strip()
        else:
            pvar = repr(var)
        
        
        retstr += "\n%(indent)s%(key)s %(type)s =  %(val)s" % {
                                                        "indent": fulltab,
                                                        "key": tc.colored(
                                                                _pad_str(name, namewidth)
                                                                , attrs=["bold"]),
                                                        "type": tc.colored(stype, attrs=["dark"]),
                                                        "val": pvar
                                                       }
    if indent == 0:
        retstr = retstr.strip()
    return retstr

def maxkeylen (d):
    keys = d.keys()
    kl = 0
    for key in keys:
        keylen = len(key)
        if keylen>kl:
            kl = keylen
    return kl
            
def _pad_str(tstr, length):
    slen = len(tstr)
    pad = " " * (length-slen)
    return tstr+pad

def context_args(context_args):
    """Converts a user supplied list in the format returned by argparse
    e.g. [["a","=","100"],["a=datetime.now()"]]]
    
    The list is turned into a single line, the result split by "=". No equals
    means set the arg to True.  The value is eval-ed.
    """
    d = {}
    # print "r25:", argl, argv
    if context_args:
        from datetime import datetime, date, timedelta
        for item in context_args:
            # print "r27:",item
            expr = "".join(item)
            if  "=" in expr:
                lrval  = expr.split("=")
                lval = lrval[0]
                if len(lrval)>1:
                    theval = eval(lrval[1])
                else:
                    theval = True
                d[lval] = theval
            else:
                if   len(item) == 2:
                    thekey = item[0]
                    theval = eval(item[1])
                    d[thekey] = theval
                elif len(item) == 1:
                    d[item[0]] = True
                    
    return d

# debugging help
# NOTE: uses inspect... might slow down inner loops
def lineno():
    import inspect
    return inspect.currentframe().f_back.f_lineno
    
def called_me(extra_frame = 0):
    import inspect
    (frame, 
    filename, line_number, 
    function_name, lines, 
    index) = inspect.getouterframes(inspect.currentframe())[2+extra_frame]
    fi = {
        "frame":         frame,
        "filename":      filename,
        "line_number":   line_number,
        "function_name": function_name,
        "lines":         lines,
        "index":         index
        }
    return fi
        

