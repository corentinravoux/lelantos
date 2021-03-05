#!/usr/bin/env python3
from lslyatomo import tomographic_objects
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Cut the DLA catalog')
parser.add_argument('-in',
                    '--input', help='Input DLA catalog',required=True)
parser.add_argument('-out',
                    '--output', help='Output DLA catalog',required=True)
parser.add_argument('-ramin', help='min RA coord to cut', default=-np.inf,required=False)
parser.add_argument('-ramax', help='max RA coord to cut', default=np.inf,required=False)
parser.add_argument('-decmin', help='min RA coord to cut', default=-np.inf,required=False)
parser.add_argument('-decmax', help='max DEC coord to cut', default=np.inf,required=False)
parser.add_argument('-zmin', help='min z coord to cut', default=-np.inf,required=False)
parser.add_argument('-zmax', help='max z coord to cut', default=np.inf,required=False)
parser.add_argument('-confidencemin', help='min confidence to cut', default=None,required=False)
parser.add_argument('-nhimin', help='min nhi to cut', default=None,required=False)

args = vars(parser.parse_args())




cat = tomographic_objects.DLACatalog.init_from_fits(args["input"])

print("Number of dlas before cut", cat.coord.shape[0])

mask = cat.cut_catalog(coord_min=(float(args["ramin"]),float(args["decmin"]),float(args["zmin"])),
                       coord_max=(float(args["ramax"]),float(args["decmax"]),float(args["zmax"])),
                       confidence_min=(float(args["confidencemin"]),nhi_min=(float(args["nhimin"]))
cat.apply_mask(mask)

print("Number of dlas after cut", cat.coord.shape[0])

cat.name = args["output"]
cat.write()
