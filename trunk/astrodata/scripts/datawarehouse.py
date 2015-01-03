#!/usr/bin/env python

import getpass
import copy
import os, shutil
from astrodata import Lookups
from astrodata.adutils import ksutil as ks
from astrodata.generaldata import GeneralData

class DataWarehouse:
    stores = None
    def __init__(self):
        self.stores = {}
        
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("datasets", nargs = "*")
    parser.add_argument("--store", default=False, action="store_true")
    package_classes = Lookups.compose_multi_table(
                            "*/warehouse_settings", "warehouse_package")
    
    args = parser.parse_args()
    print ks.dict2pretty("dw25: package_classes", package_classes)
    package_class_struct = package_classes["warehouse_package"][0]
    package_type = None
    package_class = None
    for key in package_class_struct:
        package_key = key
        package_class = package_class_struct[key]
        break; # only one supported atm
    for fname in args.datasets:
        print "DATASET:", fname
        setref = GeneralData.create_data_object(fname)
        print "\ndw35: setref", setref.basename,setref.get_types()
        pkg = package_class(setref=setref)
        print "dw37: store=%s" % pkg.get_store_path(setref)
        pkg.transport_to_warehouse()
