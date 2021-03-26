#!/usr/bin/env python3
from lslyatomo import tomographic_objects
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Cut the DLA catalog')
parser.add_argument('-i',
                    '--input', help='Input DLA catalog',required=True)
parser.add_argument('-o',
                    '--output', help='Output DLA catalog',required=True)
parser.add_argument('-zmin', help='min z coord to cut', default=0.0,required=False)
parser.add_argument('-zmax', help='max z coord to cut', default=np.inf,required=False)
parser.add_argument('-confidencemin', help='min confidence to cut', default=0.0,required=False)
parser.add_argument('-nhimin', help='min nhi to cut', default=0.0,required=False)

args = vars(parser.parse_args())




cat = tomographic_objects.DLACatalog.init_from_fits(args["input"])

print("Number of dlas before cut", cat.coord.shape[0])

mask = cat.cut_catalog_dla(coord_min=float(args["zmin"]),
                           coord_max=float(args["zmax"]),
                           confidence_min=float(args["confidencemin"]),nhi_min=float(args["nhimin"]))

print("Number of dlas after cut", cat.coord.shape[0])

cat.name = args["output"]
cat.write()
