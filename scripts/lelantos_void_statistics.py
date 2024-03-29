#!/usr/bin/env python3
from lelantos import tomographic_objects
import argparse


parser = argparse.ArgumentParser(description='Print the statistics of a given void catalog')
parser.add_argument('-i',
                    '--file', help='Void catalog file name', default=None,required=True)
args = vars(parser.parse_args())

if __name__ == '__main__' :
    void_catalog = tomographic_objects.VoidCatalog.init_from_fits(args["file"])
    void_catalog.print_void_statistics()
