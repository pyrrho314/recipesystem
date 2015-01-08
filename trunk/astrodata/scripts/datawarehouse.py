#!/usr/bin/env python
import sys
import getpass
import copy
import os, shutil
from astrodata import Lookups
from astrodata.adutils import ksutil as ks
from astrodata.adutils import termcolor as tc
tc.COLOR_ON = True
from astrodata.generaldata import GeneralData
from datetime import datetime, timedelta

class DataWarehouse:
    stores = None
    def __init__(self):
        self.stores = {}
        
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("datasets", nargs = "*")
    parser.add_argument("--store", default=False, action="store_true")
    parser.add_argument("--info", default=False, action="store_true")
    parser.add_argument("--all", default=False, action="store_true")
    parser.add_argument("--manifest", default=False, action="store_true")
    parser.add_argument("--date_range", default=None)
    parser.add_argument("--shelf", default = "processed_data")
    parser.add_argument("--user", default = None)
    parser.add_argument("--settype", default = None)
    parser.add_argument("--year", default = None, type = int)
    parser.add_argument("--month", default = None, type = int)
    parser.add_argument("--day", default = None, type = int)
    parser.add_argument("--verbose", default = False, action="store_true")
    parser.add_argument("--phrase")
    
    
    package_classes = Lookups.compose_multi_table(
                            "*/warehouse_settings", "warehouse_package")
    
    args = parser.parse_args()
    package_class_list = package_classes["warehouse_package"]
    package_class_struct = package_class_list[0]
    package_type = None
    package_class = None
    for key in package_class_struct:
        package_key = key
        package_class = package_class_struct[key]
        break; # only one supported atm, always first, controlled by path order
    
    if args.info:
        print ks.dict2pretty("contributing files", package_classes["_contributors"])
        for i in range(len(package_class_list)):
            package_def = package_class_list[i]
            print ks.dict2pretty("packager #%d" % i, package_def)
        print "choosing  %s  package class %s" % (tc.colored(key, attrs=["bold"]), tc.colored(package_class, attrs=["dark"]))
        pkg = package_class()
        print ks.dict2pretty("shelf_addresses", pkg.shelf_addresses)
        print ks.dict2pretty("type_shelf_names", pkg.type_shelf_names)
        print ks.dict2pretty("type_store_precedence", pkg.type_store_precedence)
    
    if args.manifest:
        elements = {
                "shelf_name"  : args.shelf,
                }
        pkg = package_class()
        if args.date_range:
            day = datetime.strptime(args.date_range, "%Y%m%d")
            elements.extend({
                "year"        : day.year,
                "month"       : day.month,
                "day"         : day.day,
                "complete_day": day.strftime("%Y%m%d"),
                })

        if args.year:
            elements["year"] = args.year
            daypart = "%4d" % args.year
        if args.month:
            elements["month"] = args.month
            daypart += "%02d" % args.month
        if args.day:
            elements["day"] = args.day
            daypart += "%02d" % args.day
        if not "complete_day" in elements:
            elements["complete_day"] = daypart
        
        if not args.user:
            import getpass
            user = getpass.getuser()
        else:
            user = args.user
        elements["user"] = user    
        
        if args.settype:
            elements["type"] = args.settype
        if args.phrase:
            elements["phrase"] = args.phrase
        pfx = pkg.get_store_prefix( elements = elements )    
        files = pkg.get_store_list( elements = elements )
        print "prefix = ", pfx
        print "num stored items", len(files)
        
        complete = False
        if args.verbose:
            complete = True
        print ks.dict2pretty("files", files, complete = complete)
        sys.exit()
        
        
    datasetnames = args.datasets
    if len(args.datasets)>5:
        if not args.all:
            print tc.colored( "%d datasets, showing first 5, use --all to show all"%len(args.datasets) , "red", "on_white")
            
            datasetnames = args.datasets[:5]
        
    for fname in datasetnames:
        print "  DATASET: %s" % tc.colored(fname, attrs=["bold"])
        setref = GeneralData.create_data_object(fname)
        print "    TYPES: %s " % tc.colored(", ".join( setref.get_types() ) , "blue", "on_white")
        pkg = package_class(setref=setref)
        print "STORE_KEY: %s" % pkg.get_store_path(setref)
        print "STORE_DIR: %s" % os.path.dirname(pkg.get_store_path(setref))
        print ""
        if args.store:
            pkg.transport_to_warehouse()
