#!/usr/bin/env python
import sys
import getpass
import copy
import os, shutil, re
from astrodata import Lookups
from astrodata.adutils import ksutil as ks
from astrodata.adutils import termcolor as tc
tc.COLOR_ON = True
from astrodata.generaldata import GeneralData
from datetime import datetime, timedelta
import subprocess
import glob
class DataWarehouse:
    stores = None
    def __init__(self):
        self.stores = {}
        

def store_datasets(dataset_names):        
    datasetnames = dataset_names
    if len(args.datasets)>5:
        if not args.all:
            print tc.colored( "%d datasets, showing first 5, use --all to show all"%len(args.datasets) , "red", "on_white")
            
            datasetnames = args.datasets[:5]
        
    for fname in datasetnames:
        print "  DATASET: %s" % tc.colored(fname, attrs=["bold"])
        setref = GeneralData.create_data_object(fname)
        print "    TYPES: %s" % tc.colored(", ".join( setref.get_types() ) , "blue", "on_white")
        pkg = package_class(setref=setref)
        print "STORE_KEY: %s" % pkg.get_store_path(setref)
        print "STORE_DIR: %s" % os.path.dirname(pkg.get_store_path(setref))
        print ""
        if args.store:
            pkg.transport_to_warehouse()
                    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("datasets", nargs = "*")
    parser.add_argument("--store", default=False, action="store_true")
    parser.add_argument("--fetch", default=False, action="store_true")
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
    parser.add_argument("--recipe", default = None)
    parser.add_argument("--verbose", default = False, action="store_true")
    parser.add_argument("--region", default = "NPH")
    parser.add_argument("--phrase")
    
    
    package_classes = Lookups.compose_multi_table(
                            "*/warehouse_settings", "warehouse_package")
    
    args = parser.parse_args()
    package_class_list = package_classes["warehouse_package"]
    package_class_struct = package_class_list[0]
    package_type = None
    package_class = None
    
    # some flags imply others
    if args.store:
        args.all = True
    
    if args.region:
        if args.region == "NPH":
            args.region_legacy = "Northern_AOI"
        elif args.region == "SPH":
            args.region_legacy = "Southern_AOI"
        
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
    
    if args.fetch:
        args.manifest = True
    if args.manifest:
        elements = {
                "shelf_name"  : args.shelf,
                }
        pkg = package_class()
        if args.date_range:
            parts = args.date_range.split("-")
            day_start = day_finish = None
            day_start = parts[0]
            if len(parts) > 1:
                day_finish = parts[1]
            
            date_start = datetime.strptime(day_start, "%Y%m%d")
            elements.update({
                "year"        : date_start.year,
                "month"       : date_start.month,
                "day"         : date_start.day,
                "complete_day": date_start.strftime("%Y%m%d"),
                })
            elements.update({
                "year_start"         : date_start.year,
                "month_start"        : date_start.month,
                "day_start"          : date_start.day,
                "complete_day_start" : date_start.strftime("%Y%m%d"),
                })
            
            if day_finish:
                date_finish = datetime.strptime(day_finish, "%Y%m%d")
                elements.update({
                    "year_end"         : date_finish.year,
                    "month_end"        : date_finish.month,
                    "day_end"          : date_finish.day,
                    "complete_day_end" : date_start.strftime("%Y%m%d")
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
            
        if args.region:
            elements["region"] = args.region
            elements["region_legacy"] = args.region_legacy
        pfx = pkg.get_store_prefix( elements = elements )
            
        # declare files one of two ways
        if "complete_day_end" not in elements:
            files = [pkg.get_store_list( elements = elements )]
        else:
            files = []
            date0 = datetime(elements["year_start"], elements["month_start"], elements["day_start"])
            date1 = datetime(elements["year_end"],   elements["month_end"],   elements["day_end"])
            for seekdate in ks.iter_date(date0, date1):
                elements.update({
                "year"        : seekdate.year,
                "month"       : seekdate.month,
                "day"         : seekdate.day,
                "complete_day": seekdate.strftime("%Y%m%d"),
                })
                sfiles = pkg.get_store_list( elements = elements)
                if len(sfiles)>0:
                    files.append(sfiles)
            
        
        print "prefix = %s (dw160)" % pfx
        print "num stored items %d (dw114)" % len(files)
        
        complete = False
        if args.verbose:
            complete = True
        print ks.dict2pretty("files", files, complete = complete)
        if args.fetch:
            print "fetching ... (dw168)"
            for filgroup in files:
                # fetch each from boto
                # get local store name
                recinps = []
                for fil in filgroup:
                    pkg = package_class(storename = fil)
                    pkg.deliver_from_warehouse()
                    if re.match(".*?.tif",fil):
                        recinps.append(pkg.local_path)
                if args.recipe:
                    cmdlist = ["kit", "-r", "%s" % args.recipe, "--invoked"]
                    cmdlist.extend(recinps)
                    exit_code = subprocess.call(cmdlist)
                    print "RECIPE PROCESS EXIT CODE: %s" % exit_code
                if args.store:
                    outdir = os.getcwd()
                    datasets= glob.glob(os.path.join(outdir, "*.tif"))
                    store_datasets(datasets)
                    ## fetch and store means wipe directory! 
                    junk = glob.glob(os.path.join(outdir,"*"))
                    for jfile in junk:
                        if os.path.isdir(jfile):
                            print "IS A DIR:",jfile
                            if jfile.endswith(".tmp"):
                                shutil.rmtree(jfile)
                            else:
                                pass
                        else:
                            os.remove(jfile)
                    
            sys.exit()
if args.store:           
    store_datasets(args.datasets)
