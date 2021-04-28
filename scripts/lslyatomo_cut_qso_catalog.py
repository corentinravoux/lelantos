#!/usr/bin/env python3
from lslyatomo import tomographic_objects
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Cut the QSO catalog')
parser.add_argument('-i',
                    '--input', help='Input QSO catalog',required=True)
parser.add_argument('-o',
                    '--output', help='Output QSO catalog',required=True)
parser.add_argument('-zname',
                    '--redshift-name', help='Pixel file name', default="Z",required=False)
parser.add_argument('-ramin', help='min RA coord to cut', default=-np.inf,required=False)
parser.add_argument('-ramax', help='max RA coord to cut', default=np.inf,required=False)
parser.add_argument('-decmin', help='min RA coord to cut', default=-np.inf,required=False)
parser.add_argument('-decmax', help='max DEC coord to cut', default=np.inf,required=False)
parser.add_argument('-zmin', help='min z coord to cut', default=-np.inf,required=False)
parser.add_argument('-zmax', help='max z coord to cut', default=np.inf,required=False)

args = vars(parser.parse_args())




cat = tomographic_objects.QSOCatalog.init_from_fits(args["input"],redshift_name=args["redshift_name"])

print("Number of quasars before cut", cat.coord.shape[0])

mask = cat.cut_catalog(coord_min=(float(args["ramin"]),float(args["decmin"]),float(args["zmin"])),
                       coord_max=(float(args["ramax"]),float(args["decmax"]),float(args["zmax"])),
                       center_x_coord=True)
cat.apply_mask(mask)
print("Number of quasars after cut", cat.coord.shape[0])

cat.name = args["output"]
cat.write()
