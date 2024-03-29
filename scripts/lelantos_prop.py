#!/usr/bin/env python3
from lelantos import tomographic_objects
import argparse


parser = argparse.ArgumentParser(description='Print the property file')
parser.add_argument('-i',
                    '--file', help='Property file name', default=None,required=True)
args = vars(parser.parse_args())

if __name__ == '__main__' :
    prop = tomographic_objects.MapPixelProperty(name = args["file"])
    prop.read()
    prop.print_prop()
