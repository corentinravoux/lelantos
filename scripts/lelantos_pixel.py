#!/usr/bin/env python3
"""Script: print the properties/statistics of a pixel (line-of-sight) file
(lelantos.tomographic_objects.Pixel)."""

from lelantos import tomographic_objects
import argparse


parser = argparse.ArgumentParser(description='Print the pixel file')
parser.add_argument('-i',
                    '--file', help='Pixel file name', default=None,required=True)
args = vars(parser.parse_args())

if __name__ == '__main__' :
    pixel = tomographic_objects.Pixel(name = args["file"])
    pixel.read()
    pixel.print_prop()
