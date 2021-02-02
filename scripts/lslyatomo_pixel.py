#!/usr/bin/env python3
from lslyatomo import tomographic_objects
import argparse


parser = argparse.ArgumentParser(description='Print the pixel file')
parser.add_argument('-f',
                    '--file', help='Pixel file name', default=None,required=True)
args = vars(parser.parse_args())

if __name__ == '__main__' :
    pixel = tomographic_objects.Pixel(name = args["file"])
    pixel.read()
    pixel.print_prop()
