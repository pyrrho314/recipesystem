#!/usr/bin/env python
import os
import sys
import re

from inspect import getsourcefile
from astrodata import AstroData
from pprint import pprint
#get color printing started
from astrodata.adutils import terminal
term = terminal
from astrodata.adutils.terminal import TerminalController

REASLSTDOUT = sys.stdout
REALSTDERR = sys.stderr
fstdout = terminal.FilteredStdout()
colorFilter = terminal.ColorFilter(True)
fstdout.addFilter(colorFilter)
sys.stdout = fstdout   
SW = 79

#from optparse import OptionParser
from copy import copy
from astrodata.AstroData import AstroData
from astrodata import RecipeManager 
from astrodata.RecipeManager import RecipeLibrary
from inspect import getsourcefile
from astrodata import Errors

centralRecipeIndex = {}
rl = RecipeLibrary()
colorFilter.on = False
primtypes = RecipeManager.centralPrimitivesIndex.keys()

class PrimInspect():
    """Tool for listing primitives, parameters, and recipes 
    """
    primsdict = None
    name2class = None
    primsdict_kbn = None
    class2instance = None
    primsets = None
    master_dict = None
    allmodules = None
    overrides = None
    
    def __init__(self, use_color=False):
        self.primsdict = {}
        self.name2class = {}
        self.primsdict_kbn = {}
        self.class2instance = {}
        self.master_dict = {}
        self.overrides = {}
        self.allmodules = []
        self.use_color = use_color
        if self.use_color:
            colorFilter.on = True
        self.build_dictionaries()
        primdictkeys = []
        primdictkeys = self.primsdict.keys()
        primdictkeys.sort(self.primsetcmp)
        self.primsets = primdictkeys
        self.master_dict = self.make_master(self.primsets)
    
    def primsetcmp(self,a,b):
        if isinstance(a,type(b)):
            return 1
        elif isinstance(b,type(a)):
            return -1
        else:
            return 0
    
    def build_dictionaries(self):
        self.create_allmodules()
        for primset in self.allmodules:
            pname = primset.__class__.__name__
            self.construct_prims_dict(primset)
            self.construct_primsclass_dict(primset.__class__)
            self.class2instance.update({pname:primset})

    def create_allmodules(self):
        for key in primtypes:
            self.allmodules.extend(rl.retrieve_primitive_set(astrotype=key))

    def construct_prims_dict(self, primset):
        self.primsdict.update({primset:self.get_prim_list(\
            primset.__class__)})
    
    def construct_primsclass_dict(self, startclass):
        if startclass.__name__== "PrimitiveSet":
            return
        self.name2class.update({startclass.__name__:startclass})
        self.primsdict_kbn.update({startclass.__name__:self.get_prim_list(\
            startclass)})
        for base in startclass.__bases__:
            self.construct_primsclass_dict(base)
    
    def get_prim_list(self, cl):
        """get a sorted list of primitive sets, sorted with parents first
        """
        plist = []
        for key in cl.__dict__:
            doappend = True
            fob = eval("cl." + key)
            if not hasattr(fob, "__call__"):
                doappend = False
            if hasattr(fob, "pt_hide"):
                doappend = eval("not fob.pt_hide")
            elif hasattr(cl, "pthide_" + key):
                doappend = eval("not cl.pthide_" + key)
            if key.startswith( "_" ):
                doappend = False   
            if doappend:
                plist.append(key)
        plist.sort()
        return plist
    
    def make_master(self, primset_objects=None):
        """Create master dictionary of adtypes, primitives and parameters.
        Will use this for the future development of rsexplorer
        """
        mdict = {}
        
        # Start with only the AD types
        from pprint import pprint
        for primset in primset_objects:
            pname = primset.__class__.__name__
            mdict.update({pname[:-10]:{'instance':primset}})
            mdict[pname[:-10]].update({'class':self.name2class[pname]})
            mdict[pname[:-10]].update({'primitives':{}})
            
            for prim in self.primsdict_kbn[pname]:
                mdict[pname[:-10]]['primitives'].update({prim:{}})
                if prim in mdict[pname[:-10]]['instance'].param_dict.keys():
                    mdict[pname[:-10]]['primitives'].update(\
                        {prim:mdict[pname[:-10]]['instance'].param_dict[prim]})
            mdict[pname[:-10]].update({'inheritance':{}})
        
        # Add inheritance (includes classes other than AD type names)
        inclass = []
        for adtype in mdict:
            order= []
            inclass.append(mdict[adtype]['class'])

            for clas in mdict[adtype]['class'].__mro__:
                if "Primitives" in clas.__name__:
                    if adtype != clas.__name__[:-10]:
                        order.append(clas.__name__[:-10])
                        mdict[adtype]['inheritance'].update(\
                            {clas.__name__[:-10]:{}})
                        mdict[adtype]['inheritance'][clas.__name__[:-10]].update(\
                            {'class':clas})
                        mdict[adtype]['inheritance'][clas.__name__[:-10]].update(\
                            {'base':clas.__base__.__name__[:-10]})
                        mdict[adtype]['inheritance'][clas.__name__[:-10]].update(\
                            {'primitives':{}})
                        for prim in self.primsdict_kbn[clas.__name__]:
                            mdict[adtype]['inheritance'][clas.__name__[:-10]]\
                                ['primitives'].update({prim:{}})
                            if prim in mdict[adtype]['instance'].param_dict.keys():
                                mdict[adtype]['inheritance'][clas.__name__[:-10]]\
                                    ['primitives'].update({prim:mdict[adtype]\
                                    ['instance'].param_dict[prim]})

            mdict[adtype]['inheritance'].update({'order':order})

        overrides = {}
        overridden = {}
        for adtype in mdict:
            # This will find the overrides for the local prims
            for pkey in mdict[adtype]['primitives'].keys():
                for primset in mdict[adtype]['inheritance']['order']:
                    primdict = mdict[adtype]['inheritance'][primset]['primitives']
                    adpdict = mdict[adtype]['primitives']
                    if pkey in  primdict.keys():
                        #print primset + ":" +  pkey + " overridden by " \
                        #+ adtype + ":" + pkey
                        overrides.update({(adtype, pkey):primset})
        
        #pprint(mdict['GENERAL'])
        self.overrides = overrides
        return mdict
    
    def hides(self, primsetname, prim, instance=None):
        """checks to see if prim hides or is hidden-by another"""
        
        if instance:
            cl = self.name2class[primsetname]
            ps = instance
            psl = rl.retrieve_primitive_set(astrotype= instance.astrotype)
        elif primsetname in self.class2instance:
            ps = self.class2instance[primsetname]
            cl = ps.__class__
            psl = rl.retrieve_primitive_set(astrotype= ps.astrotype)
        else:
            return None
        verb = False # prim == "exit"
        if verb: print "lP200:", len(psl)
        if len(psl)>1:
            before = True
            rets = None
            for ops in psl:
                
                # reason for this comparison: make this work even
                # if retrieve_primitive_set returns new instances...
                if ps.__class__.__name__ == ops.__class__.__name__:
                    before = False
                    if verb : print "lP209: found this by by class name",\
                        repr(ops)
                    continue
                if isinstance(ops, cl):
                    before = False
                    if verb : print "lp213: skipping due to isinstance"
                    continue
                if verb: print "lP215:", repr(ops),repr(dir(ops))
                if hasattr(ops, prim):
                    if verb : print "lP216: hide happens"
                    if before:
                        rets = "${RED}(hidden by " + ops.__class__.__name__ +\
                            ")${NORMAL}"
                    else:
                        rets = '${GREEN}(hides "%s" from %s)${NORMAL}' %\
                            (prim, ops.__class__.__name__)
                    break
            return rets
        return None
    
    def list_recipes(self, pkg="", eng=False, view=None):
        """
        List recipes and sub-recipes or open the recipe to view it.
        """
        retstr = "\n"
        if isinstance(view, str):
            retstr += "="*SW + "\n${RED}RECIPE: %s${NORMAL}\n" % view + "="*SW 
        else:
            retstr += "="*SW  
        cri = RecipeManager.centralRecipeIndex
        if isinstance(view, str):
            for key in cri.keys():
                if key == view:
                    retstr += "\n" + open(cri[key], "rb").read()
        else:
            topkeys = []
            engkeys = []
            subkeys = []
            for key in cri.keys():
                if "Engineering" in cri[key]:
                    engkeys.append(key)
                elif "subrecipes" in cri[key]:
                    subkeys.append(key)
                else:
                    topkeys.append(key)
            topkeys.sort()
            engkeys.sort()
            subkeys.sort()
            pkg = "RECIPES_" + pkg
            retstr += self.list_recipes_str(pkg, topkeys)
            retstr += self.list_recipes_str("Subrecipes", subkeys)
            if eng:
                retstr += self.list_recipes_str("Engineering", engkeys) 
        retstr += "\n\n" + "="*SW
        print(retstr)
    
    def list_recipes_str(self, topdir="", rlist=[]):
        rstr = ""
        rstr += "\n\n${YELLOW}%s${NORMAL}\n" % topdir + "-"*SW
        count = 1
        for r in rlist:
            rstr += "\n    %s. %s" % (count, r)
            count +=1
        return rstr
    
    def list_primitives(self, data=None, adtype=None, info=None, params=None):
        """List primitives and associated parameters
        """
        rstr = ['\n']
        
        adlist = []
        if adtype is None and data is None:
            adlist =  self.master_dict.keys()
            adlist.sort()
        else:
            if data:
                ad = AstroData(data)
                ps = rl.retrieve_primitive_set(dataset=ad)
                adtype = ps[0].astrotype
            adlist.append(adtype)
        for adtype in adlist:
            rstr.append("\n\n" + "="*SW + "\n${RED}" + adtype + " ${NORMAL}\n"\
                        + "="*SW)
            if info:
                clas = self.master_dict[adtype]['class']
                sfull = getsourcefile(clas)
                sdir = os.path.dirname(sfull)
                sfil = os.path.basename(sfull)
                rstr.append("\n${YELLOW}%-13s:%s${NORMAL}" % \
                            ("Class",clas.__name__)) 
                rstr.append("\n${YELLOW}%-13s:%s${NORMAL}" % ("Source",sfull))
                if len(self.master_dict[adtype]["inheritance"]["order"]) == 0:
                    rstr.append("\n${YELLOW}%-13s:%s${NORMAL}" % \
                                ("Inheritance","None"))
                else:
                    rstr.append("\n${YELLOW}%-13s:%s${NORMAL}" % \
                                ("Inheritance","Yes"))
                rstr.append("\n" + "-"*SW)

            primkeys = self.master_dict[adtype]['primitives'].keys()
            count = 0
            primkeys.sort()

            for prim in primkeys:
                count += 1
                rstr.append("\n" + str(count) + ". " + prim)
                if (adtype, prim) in self.overrides.keys():
                    rstr.append("${BLUE} (OVERRIDES " + self.overrides[\
                                (adtype, prim)] + ")${NORMAL}")
                primdict = self.master_dict[adtype]['primitives'][prim]
                if params:
                    rstr.extend(self.add_params(primdict, rstr, info, 8))

            # inheritance
            TW = 4
            sav = ""
            if len(self.master_dict[adtype]['inheritance']['order']) > 0:
                rstr.append("\n\n${RED}        -------- Inheritance by Method Resolution Order --------")
                for primset in self.master_dict[adtype]['inheritance']['order']:
                    rstr.append("\n\n" + " "*TW + "${BLUE}(" + primset + \
                                ")${NORMAL}")
                    primsort = self.master_dict[adtype]['inheritance']\
                                                  [primset]['primitives'].keys()
                    primsort.sort()
                    for prim in primsort:
                        count += 1
                        rstr.append("\n" + " "*TW + str(count) + ". " + prim)
                        if (primset, prim) in self.overrides.keys():
                            rstr.append("${BLUE} (OVERRIDES " + self.overrides\
                                       [(primset, prim)] + ")${NORMAL}")
                        for key in self.overrides.keys():
                            if self.overrides[key] == primset and \
                                prim == key[1] and (key[0] in self.master_dict\
                                [adtype]['inheritance']['order'] or key[0] == \
                                adtype):
                                rstr.append("${RED} (OVERRIDDEN BY " + key[0]\
                                          + ")${NORMAL}")
                                
                        primdict = self.master_dict[adtype]['inheritance']\
                                       [primset]['primitives'][prim]
                        if params:
                            rstr.extend(self.add_params(primdict, rstr, info, 12))

        rstr.append("\n" + "="*SW)
        print("".join(rstr))

    def add_params(self, primdict=None, rstr=None, info=None, indent=None):
        rstr = []
        if len(primdict) > 0:
            for param in primdict:
                if "default" in primdict[param].keys():
                    rstr.append("\n" + " "*indent + "${YELLOW}" + param + ": "\
                                + repr(primdict[param]["default"]) + "${NORMAL}")
                else:
                    rstr.append("\n" + " "*indent + "${YELLOW}" + param + \
                                ": (No default)${NORMAL}")
                         
                if info:
                    for meta in primdict[param]:
                        rstr.append("\n%s%-16s%s" % (" "*(indent + 4), meta,\
                                    ":" + repr(primdict[param][meta])))
        else:
            rstr.append("\n" + " "*indent + "${YELLOW}(No Parameters)${NORMAL}")
        return rstr

    def list_primsets(self, info=None):
        """List the primitive sets (astrodata types), includes inheritance
        and source file.
        """
        rstr = []
        rstr.append("="*SW + "\n")
        count = 1
        for adtype in self.master_dict.keys():
            rstr.append("\n${RED}" + str(count) + ". " + adtype + \
                        " ${NORMAL}\n")
            count += 1
            if info:
                clas = self.master_dict[adtype]['class']
                sfull = getsourcefile(clas)
                sdir = os.path.dirname(sfull)
                sfil = os.path.basename(sfull)
                rstr.append("\n${YELLOW}%-13s:%s${NORMAL}" % \
                           ("Class",clas.__name__)) 
                rstr.append("\n${YELLOW}%-13s:%s${NORMAL}" % ("Source",sfull))
                olist = self.master_dict[adtype]["inheritance"]["order"]
                if len(olist) == 0:
                    rstr.append("\n${YELLOW}%-13s:%s${NORMAL}" % \
                           ("Inheritance","None"))
                elif olist[0] not in self.master_dict.keys():
                    rstr.append("\n${YELLOW}%-13s: %s${NORMAL}" % \
                                ("Inheritance", "Multiple"))
                else:
                    rstr.append("\n${YELLOW}%-13s: %s${NORMAL}" % \
                                ("Inheritance", olist[0]))
                rstr.append("\n${YELLOW}%-13s: %s${NORMAL}" % \
                            ("Local Primitives", str(len(\
                            self.master_dict[adtype]['primitives'].keys()))))
                
                rstr.append("\n" + "-"*SW)
        rstr.append("\n" + "="*SW)
        print("".join(rstr))

        
