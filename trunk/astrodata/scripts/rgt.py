#!/bin/env python
import argparse

parser = argparse.ArgumentParser("Used to print usefull information from log files.")

parser.add_argument("textfiles", nargs = "+", help = """A log or source file.
                                                     """)
parser.add_argument("-t", "--trace_log", 
                    help="Print Recipe Trace, nested Recipe and Primitive Names", 
                    default=False, action = "store_true")

args = parser.parse_args()

for tfilnam in args.textfiles:
    tfil = open (tfilnam)
    semantype = None
    for line in tfil:
        if "running recipe" in line:
            semantype = "kitchen_log"
            print line,
        if args.trace_log:
            if "PRIMITIVE" in line or "RECIPE" in line:
                print line,
                