
def dict2pretty(name, var, indent=0, namewidth = None):
    retstr = ""
    tabspc = " "*4
    fulltab = tabspc*indent
    if not namewidth:
        namewidth = len(name)
    if isinstance(var, dict):
        retstr += "\n%(indent)s%(key)s:" % {
                                                    "indent":fulltab,
                                                    "key": name,
                                                    "extra": tabspc
                                                 }
        namewidth = maxkeylen(var)
        if len(var) == 0:
            retstr += "\n%(indent)s%(tab)s:::empty:::" % {"indent":fulltab,
                                                          "tab":tabspc
                                                          }
        keys = var.keys()
        keys.sort()
        for key in keys:
            value = var[key]
            #print key,value
            retstr += dict2pretty(key, value, indent+1, namewidth = namewidth)
    elif isinstance(var, list):
         retstr += "\n%(indent)s%(key)s" % { "indent":fulltab,
                                                    "key": name
                                               }
         
         listlen = len(var)
         
         if len(var) < 50:
             oneline = ", ".join(var)
             if len(oneline)<120:
                return dict2pretty(name, oneline, indent)
            
         if listlen > 10:
            last = listlen - 1
            mid = int(last/2);
            retstr += dict2pretty("[0]", var[0], indent+1, namewidth = 9) 
            retstr += dict2pretty("[%d]"%mid, var[mid],indent+1, namewidth = 9)
            retstr += dict2pretty("[%d]"%last, var[last], indent+1, namewidth = 9)  
         else:
            for i in range(0, listlen):
                key = "[%d]"%i
                value = var[i];
                retstr += dict2pretty(key, value, indent+1, namewidth = 3)
    else:
        retstr += "\n%(indent)s%(key)s = %(val)s" % {"indent": fulltab,
                                                       "key": _pad_str(name, namewidth),
                                                       "val": var
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
        

