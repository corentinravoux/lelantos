from lslyatomo import tomographic_objects
import numpy as np


shape = (450,450,450)
pix_size = 2.19
size = tuple((np.array(shape)-1)*pix_size)

Omega_m = 0.3147
minx = 0.0  
miny = 0.0
minredshift = 2.1
name = "property_file.pickle"
coordinate_transform = "middle"

prop = tomographic_objects.MapPixelProperty.init_false_prop(shape,size,Omega_m,minx,miny,minredshift,coordinate_transform,name=name)
prop.write()
