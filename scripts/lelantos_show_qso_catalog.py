#!/usr/bin/env python3
import matplotlib.pyplot as plt
import argparse
from lelantos import tomographic_objects

parser = argparse.ArgumentParser(description='Plot the QSO catalog')
parser.add_argument('-i',
                    '--input', help='Input QSO catalog',required=True)
parser.add_argument('-bins', help='Output QSO catalog', default=100,required=False)
parser.add_argument('-zname',
                    '--redshift-name', help='Pixel file name', default="Z",required=False)
args = vars(parser.parse_args())

catalog = tomographic_objects.QSOCatalog.init_from_fits(args["input"],redshift_name=args["redshift_name"])

RA = catalog.coord[:,0]
RA[RA>180] = RA[RA>180] - 360


DEC = catalog.coord[:,1]
plt.hist2d(RA,DEC,int(args["bins"]))
plt.colorbar()
plt.show()
