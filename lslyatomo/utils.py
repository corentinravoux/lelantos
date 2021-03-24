#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

Author: Corentin Ravoux

Description : Void and Over-density finder for Lya Tomographic maps.
Watershed and Simple Spherical techniques are available.
Tested on irene and cobalt (CCRT)
"""



#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################



import numpy as np
import logging, time, os
import fitsio
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.ndimage import map_coordinates
from scipy.optimize import leastsq
import multiprocessing as mp
from multiprocessing.reduction import ForkingPickler, AbstractReducer
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic

#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################

lambdaLy = 1215.673123130217


def mpc_per_pixel(size,shape):
    # return(np.array(size)/(np.array(shape)))
    return(np.array(size)/(np.array(shape)-1))

def pixel_per_mpc(size,shape):
    # return((np.array(shape))/np.array(size))
    return((np.array(shape)-1)/np.array(size))


def get_cosmo_function(Omega_m,Omega_k=0.):
    try:
        from picca import constants
    except:
        import lslyatomo.picca.constants as constants
        print("Picca might be updated, we suggest to install picca independently")
    try:
        Cosmo = constants.cosmo(Omega_m,Ok=Omega_k)
        rcomov = Cosmo.r_comoving
        distang = Cosmo.dm
    except:
        Cosmo = constants.Cosmo(Omega_m,Ok=Omega_k)
        rcomov = Cosmo.get_r_comov
        distang = Cosmo.get_dist_m
    redshift_array = np.linspace(0,5,10000)
    R_array = rcomov(redshift_array)
    Dm_array = rcomov(redshift_array)
    inv_rcomov = interpolate.interp1d(R_array,redshift_array)
    inv_distang = interpolate.interp1d(Dm_array,redshift_array)
    return(rcomov,distang,inv_rcomov,inv_distang)


def return_suplementary_parameters(mode,property=None,zmin=None,zmax=None):
    if(mode == "middle"):
        if(property is not None):
            zmin = property.boundary_sky_coord[0][2]
            zmax = property.boundary_sky_coord[0][1]
            suplementary_parameters = [(zmin + zmax)/2]
        elif((zmin is not None)&(zmin is not None)):
            suplementary_parameters = [(zmin + zmax)/2]
    else:
        suplementary_parameters = None
    return(suplementary_parameters)




def convert_cartesian_to_sky(X,Y,Z,method,inv_rcomov=None,inv_distang=None,distang=None,suplementary_parameters=None):
    """Mpc.h-1 to radians"""
    if(method == "full_angle"):
        (RA,DEC,z)= convert_cartesian_to_sky_full_angle(X,Y,Z,inv_rcomov)
    if(method == "full"):
        (RA,DEC,z)= convert_cartesian_to_sky_full(X,Y,Z,inv_rcomov)
    elif(method == "middle"):
        (RA,DEC,z)= convert_cartesian_to_sky_middle(X,Y,Z,inv_rcomov,distang,suplementary_parameters[0])
    return(RA,DEC,z)

def convert_sky_to_cartesian(RA,DEC,z,method,rcomov=None,distang=None,suplementary_parameters=None):
    """Radians to Mpc.h-1"""
    if(method == "full_angle"):
        (X,Y,Z)= convert_sky_to_cartesian_full_angle(RA,DEC,z,rcomov)
    elif(method == "full"):
        (X,Y,Z)= convert_sky_to_cartesian_full(RA,DEC,z,rcomov)
    elif(method == "middle"):
        (X,Y,Z)= convert_sky_to_cartesian_middle(RA,DEC,z,rcomov,distang,suplementary_parameters[0])
    return(X,Y,Z)

def convert_cartesian_to_sky_full_angle(X,Y,Z,inv_rcomov):
    RA = np.arctan2(X,Z)
    DEC = np.arcsin(Y/np.sqrt(X**2 + Y**2 + Z**2))
    z = inv_rcomov(np.sqrt(X**2 + Y**2 + Z**2))
    return(RA,DEC,z)

def convert_sky_to_cartesian_full_angle(RA,DEC,z,rcomov):
    X = rcomov(z)*np.sin(RA)*np.cos(DEC)
    Y = rcomov(z)*np.sin(DEC)
    Z = rcomov(z)*np.cos(RA)*np.cos(DEC)
    return(X,Y,Z)

def convert_cartesian_to_sky_full(X,Y,Z,inv_rcomov):
    z = inv_rcomov(np.sqrt(X**2 + Y**2 + Z**2))
    RA = X/Z
    DEC = Y/Z
    return(RA,DEC,z)

def convert_sky_to_cartesian_full(RA,DEC,z,rcomov):
    X = rcomov(z)*RA
    Y = rcomov(z)*DEC
    Z = rcomov(z)
    return(X,Y,Z)


def convert_cartesian_to_sky_middle(X,Y,Z,inv_rcomov,distang,middle_z):
    RA = X/distang(middle_z)
    DEC = Y/distang(middle_z)
    z = inv_rcomov(Z)
    return(RA,DEC,z)

def convert_sky_to_cartesian_middle(RA,DEC,z,rcomov,distang,middle_z):
    X = distang(middle_z)*RA
    Y = distang(middle_z)*DEC
    Z = rcomov(z)
    return(X,Y,Z)

def convert_z_cartesian_to_sky_middle(Z,inv_rcomov):
    z = inv_rcomov(Z)
    return(z)

def convert_z_sky_to_cartesian_middle(z,rcomov):
    Z = rcomov(z)
    return(Z)

def get_direction_indexes(direction,rotate):
    if (direction.lower() == "x")|(direction.lower() == "ra"):
        x_index, y_index, index_direction = 2, 1, 0
    elif (direction.lower() == "y")|(direction.lower() == "dec"):
        x_index, y_index, index_direction = 2, 0, 1
    elif (direction.lower() == "z")|(direction.lower() == "redshift"):
        x_index, y_index, index_direction = 1, 0, 2
    if(rotate): x_index,y_index = y_index,x_index
    index_dict = {0:"x",1:"y",2:"z"}
    return(x_index,y_index,index_direction,index_dict)


saclay_mock_lines_per_box = {"box" : 5,"vx" : 5,"vy" : 5,"vz" : 5,"eta_xx":1,"eta_xy":1,"eta_xz":1,"eta_yx":1,"eta_yy":1,"eta_yz":1,"eta_zx":1,"eta_zy":1,"eta_zz":1}

def saclay_mock_box_cosmo_parameters(box_shape,size_cell):
    import cosmolopy.distance as dist
    try :
        from SaclayMocks import constant as saclay_mock_constant
    except :
        from lslyatomo.saclaymocks import constant as saclay_mock_constant
        print("SaclayMocks might be updated, we suggest to install SaclayMocks independently")
    NZ = box_shape[2]
    DZ = size_cell
    LZ = NZ*DZ

    h = saclay_mock_constant.h
    Om = saclay_mock_constant.omega_M_0
    OL = saclay_mock_constant.omega_lambda_0
    Ok = saclay_mock_constant.omega_k_0
    z0 = saclay_mock_constant.z0

    cosmo_fid = {'omega_M_0':Om, 'omega_lambda_0':OL, 'omega_k_0':Ok, 'h':h}
    R_of_z, z_of_R = dist.quick_distance_function(dist.comoving_distance, return_inverse=True, **cosmo_fid)
    R0 = h * R_of_z(z0)
    Rmin = R0 - LZ/2
    Rmax = R0 + LZ/2
    return(R0,z0,R_of_z,z_of_R,Rmin,Rmax,h)

def saclay_mock_center_of_the_box(box_bound):
    ra0_box = (box_bound[0] + box_bound[1])/2
    dec0_box = (box_bound[2] + box_bound[3])/2
    return(ra0_box, dec0_box)


def saclay_mock_coord_dm_map(X,Y,Z,Rmin,size_cell,box_shape,interpolation_method):
    size_cell = size_cell
    center_x = (box_shape[0]-1) / 2
    center_y = (box_shape[1]-1) / 2
    if(interpolation_method.upper() == "NEAREST"):
        n_i = (np.round(center_x + X/size_cell,0)).astype(int)
        n_j = (np.round(center_y + Y/size_cell,0)).astype(int)
        n_k = np.round((Z - Rmin)/size_cell,0).astype(int)
    else:
        n_i = center_x + X/size_cell
        n_j = center_y + Y/size_cell
        n_k = (Z - Rmin)/size_cell
    return(n_i,n_j,n_k)

def saclay_mock_read_box(box_dir,n_x,name_box):
    name = "{}-{}.fits".format(name_box,str(n_x))
    box = fitsio.FITS(os.path.join(box_dir, name))[0][:,:,:][0]
    return(box)

def saclay_mock_get_box(box_dir,box_shape,name_box="box"):
    line_per_box = saclay_mock_lines_per_box[name_box]
    DM_mocks = np.zeros((box_shape[0],box_shape[1],box_shape[2]))
    for i in range(box_shape[0]//line_per_box):
        DM_mocks[i*line_per_box:(i+1)*line_per_box,:,:] = saclay_mock_read_box(box_dir,i,name_box)[:,:]
    return(DM_mocks)


def saclay_mock_sky_to_cartesian(ra,dec,R,ra0,dec0):
    '''
    XYZ of a point P (ra,dec,R) in a frame with
    observer at O, Z along OP, X along ra0, Y along dec0
    angles in radians
    tested that ra,dec, R = box.ComputeRaDecR(R0,ra0,dec0,X,Y,Z)
    x,y,z = box.ComputeXYZ(ra[0],dec[0],R,ra0,dec0)
    print x-X,y-Y,z-R0-Z        prints ~1E-13  for random inputs
    '''
    try :
        from SaclayMocks import box as saclay_mock_box
    except :
        from lslyatomo.saclaymocks import box as saclay_mock_box
        print("SaclayMocks might be updated, we suggest to install SaclayMocks independently")
    X,Y,Z = saclay_mock_box.ComputeXYZ2(ra*(np.pi/180),dec*(np.pi/180),R,ra0*(np.pi/180),dec0*(np.pi/180))
    return(X,Y,Z)




def cut_sky_catalog(ra,dec,z,ramin=None,ramax=None,decmin=None,decmax=None,zmin=None,zmax=None):
    mask = np.full(ra.shape,True)
    if(ramin is not None):
        mask &= ra > np.radians(ramin)
    if(ramax is not None):
        mask &= ra < np.radians(ramax)
    if(decmin is not None):
        mask &= dec > np.radians(decmin)
    if(decmax is not None):
        mask &= dec < np.radians(decmax)
    if(zmin is not None):
        mask &= z > zmin
    if(zmax is not None):
        mask &= z < zmax
    return(mask)

def init_shared_array(shape,full_value=np.inf):
    distance_array = np.full(shape,full_value)
    shared_arr = mp.Array('d', distance_array.flatten())
    del distance_array
    return(shared_arr)

def mp_array_to_numpyarray(mp_arr):
    return np.frombuffer(mp_arr.get_obj())


def bin_ndarray(ndarray, new_shape, operation='mean'):
    """
    From : https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array/29042041
    Bins an ndarray in all axes based on the target shape, by summing or
    averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
    [102 110 118 126 134]
    [182 190 198 206 214]
    [262 270 278 286 294]
    [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg','gauss']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
        elif operation.lower() in ["gauss"]:
            if i!=0 : raise KeyError("gaussian mean is not available for dim higher than 1")
            from scipy import signal
            newndarray = np.zeros(new_shape)
            gaussian_weights = signal.gaussian(int(ndarray.shape[1]),int(ndarray.shape[1])/4)
            for j in range(len(ndarray)):
                newndarray[j] = np.average(ndarray[j],axis=0,weights=gaussian_weights)
            ndarray =newndarray
    return ndarray


def interpolate_map(interpolation_method,map_array,coord):
    if(interpolation_method.upper() == "NEAREST"):
        coord = np.around(coord,decimals=0).astype(int)
        DM_map = map_array[coord[:,:,:,0],coord[:,:,:,1],coord[:,:,:,2]]
    elif(interpolation_method.upper() == "LINEAR"):
        points =  coord.reshape(coord.shape[0]*coord.shape[1]*coord.shape[2],3)
        DM_map = map_coordinates(map_array, np.transpose(points), order=1).reshape((coord.shape[0],coord.shape[1],coord.shape[2]))
    elif(interpolation_method.upper() == "SPLINE"):
        points =  coord.reshape(coord.shape[0]*coord.shape[1]*coord.shape[2],3)
        DM_map = map_coordinates(map_array, np.transpose(points), order=2).reshape((coord.shape[0],coord.shape[1],coord.shape[2]))
    else :
        raise ValueError("Please select NEAREST, LINEAR or SPLINE as interpolation_method")
    return(DM_map)


def interpolate_and_fill_map(interpolation_method,map_array,coord):
    if(interpolation_method.upper() == "NEAREST"):
        coord_nearest = np.around(coord,decimals=0).astype(int)
        map_to_fill = map_array[coord_nearest[:,0],coord_nearest[:,1],coord_nearest[:,2]]
        del coord_nearest
    elif(interpolation_method.upper() == "LINEAR"):
        map_to_fill = map_coordinates(map_array, np.transpose(coord), order=1)
    elif(interpolation_method.upper() == "SPLINE"):
        map_to_fill = map_coordinates(map_array, np.transpose(coord), order=2)
    else :
        raise ValueError("Please select NEAREST, LINEAR or SPLINE as interpolation_method")
    return(map_to_fill)


def gaussian_smoothing(mapdata,sigma):
    gaussian_map = gaussian_filter(mapdata,sigma)
    return(gaussian_map)



class gaussian_fitter_2d(object):

    def __init__(self,inpdata=None):

        self.inpdata = inpdata

    def moments2D(self):
        """ Returns the (amplitude, xcenter, ycenter, xsigma, ysigma, rot, bkg, e) estimated from moments in the 2d input array Data """

        bkg=np.median(np.hstack((self.inpdata[0,:],self.inpdata[-1,:],self.inpdata[:,0],self.inpdata[:,-1])))  #Taking median of the 4 edges points as background
        Data=np.ma.masked_less(self.inpdata-bkg,0)   #Removing the background for calculating moments of pure 2D gaussian
        #We also masked any negative values before measuring moments

        amplitude=Data.max()

        total= float(Data.sum())
        Xcoords,Ycoords= np.indices(Data.shape)

        xcenter= (Xcoords*Data).sum()/total
        ycenter= (Ycoords*Data).sum()/total

        RowCut= Data[int(xcenter),:]  # Cut along the row of data near center of gaussian
        ColumnCut= Data[:,int(ycenter)]  # Cut along the column of data near center of gaussian
        xsigma= np.sqrt(np.ma.sum(ColumnCut* (np.arange(len(ColumnCut))-xcenter)**2)/ColumnCut.sum())
        ysigma= np.sqrt(np.ma.sum(RowCut* (np.arange(len(RowCut))-ycenter)**2)/RowCut.sum())

        #Ellipcity and position angle calculation
        Mxx= np.ma.sum((Xcoords-xcenter)*(Xcoords-xcenter) * Data) /total
        Myy= np.ma.sum((Ycoords-ycenter)*(Ycoords-ycenter) * Data) /total
        Mxy= np.ma.sum((Xcoords-xcenter)*(Ycoords-ycenter) * Data) /total
        e= np.sqrt((Mxx - Myy)**2 + (2*Mxy)**2) / (Mxx + Myy)
        pa= 0.5 * np.arctan(2*Mxy / (Mxx - Myy))
        rot= np.rad2deg(pa)

        return amplitude,xcenter,ycenter,xsigma,ysigma, rot,bkg, e

    def Gaussian2D(self,amplitude, xcenter, ycenter, xsigma, ysigma, rot,bkg):
        """ Returns a 2D Gaussian function with input parameters. rotation input rot should be in degress """
        rot=np.deg2rad(rot)  #Converting to radians
        Xc=xcenter*np.cos(rot) - ycenter*np.sin(rot)  #Centers in rotated coordinates
        Yc=xcenter*np.sin(rot) + ycenter*np.cos(rot)
        #Now lets define the 2D gaussian function
        def Gauss2D(x,y) :
            """ Returns the values of the defined 2d gaussian at x,y """
            xr=x * np.cos(rot) - y * np.sin(rot)  #X position in rotated coordinates
            yr=x * np.sin(rot) + y * np.cos(rot)
            return amplitude*np.exp(-(((xr-Xc)/xsigma)**2 +((yr-Yc)/ysigma)**2)/2) +bkg

        return Gauss2D

    def FitGauss2D(self,ip=None):
        """ Fits 2D gaussian to Data with optional Initial conditions ip=(amplitude, xcenter, ycenter, xsigma, ysigma, rot, bkg)
        Example:
        >>> X,Y=np.indices((40,40),dtype=np.float)
        >>> Data=np.exp(-(((X-25)/5)**2 +((Y-15)/10)**2)/2) + 1
        >>> FitGauss2D(Data)
        (array([  1.00000000e+00,   2.50000000e+01,   1.50000000e+01, 5.00000000e+00,   1.00000000e+01,   2.09859373e-07, 1]), 2)
         """
        if ip is None:   #Estimate the initial parameters form moments and also set rot angle to be 0
            ip=self.moments2D()[:-1]   #Remove ellipticity from the end in parameter list

        Xcoords,Ycoords= np.indices(self.inpdata.shape)
        def errfun(ip):
            dXcoords= Xcoords-ip[1]
            dYcoords= Ycoords-ip[2]
            Weights=np.sqrt(np.square(dXcoords)+np.square(dYcoords)) # Taking radius as the weights for least square fitting
            return np.ravel((self.Gaussian2D(*ip)(*np.indices(self.inpdata.shape)) - self.inpdata)/np.sqrt(Weights))  #Taking a sqrt(weight) here so that while scipy takes square of this array it will become 1/r weight.

        p, success = leastsq(errfun, ip)

        return p,success



def create_log(log_level="info"):
    log = Logger(log_level=log_level)
    log.setup_logging()
    return(log)

def create_report_log(name="Python_Report",log_level="info"):
    log = Logger(name=name,log_level=log_level)
    log.setup_report_logging()
    return(log)

_logging_handler = None

class Logger(object):

    def __init__(self,name="Python_Report",log_level="info"):
        self.name = name
        self.log_level = log_level


    def setup_logging(self):
    	"""
        Taken from https://nbodykit.readthedocs.io/
    	Turn on logging, with the specified level.
    	Parameters
    	----------
    	log_level : 'info', 'debug', 'warning'
    		the logging level to set; logging below this level is ignored
    	"""

    	# This gives:
    	#
    	# [ 000000.43 ]   0: 06-28 14:49  measurestats	INFO	 Nproc = [2, 1, 1]
    	# [ 000000.43 ]   0: 06-28 14:49  measurestats	INFO	 Rmax = 120

    	levels = {"info" : logging.INFO,"debug" : logging.DEBUG,"warning" : logging.WARNING}

    	logger = logging.getLogger();
    	t0 = time.time()


    	class Formatter(logging.Formatter):
    		def format(self, record):
    			s1 = ('[ %09.2f ]: ' % (time.time() - t0))
    			return s1 + logging.Formatter.format(self, record)

    	fmt = Formatter(fmt='%(asctime)s %(name)-15s %(levelname)-8s %(message)s',
    					datefmt='%m-%d %H:%M ')

    	global _logging_handler
    	if _logging_handler is None:
    		_logging_handler = logging.StreamHandler()
    		logger.addHandler(_logging_handler)

    	_logging_handler.setFormatter(fmt)
    	logger.setLevel(levels[self.log_level])



    def setup_report_logging(self):
        levels = {"info" : logging.INFO,"debug" : logging.DEBUG,"warning" : logging.WARNING}
        logging.basicConfig(filename=self.name, filemode='w',level=levels[self.log_level],format='%(asctime)s :: %(levelname)s :: %(message)s')



    @staticmethod
    def add(line,level="info"):
        if(level=="info"):
            logging.info(line)
        if(level=="warning"):
            logging.warning(line)
        if(level=="debug"):
            logging.debug(line)

    @staticmethod
    def add_array_statistics(arr,char):
        if(arr is not None):
            Logger.add(f"Min of {char}: {arr.min()}")
            Logger.add(f"Max of {char}: {arr.max()}")
            Logger.add(f"Mean of {char}: {arr.mean()}")
            Logger.add(f"Standard deviation of {char}: {arr.std()}")


    @staticmethod
    def close():
        logging.shutdown()






def latex_float(float_input,decimals_input="{0:.2g}"):
    """
    example use:
    import matplotlib.pyplot as plt
    plt.figure(),plt.clf()
    plt.plot(np.array([1,2.]),'ko-',label="$P_0="+latex_float(7.63e-5)+'$'),
    plt.legend()
    """
    float_str = decimals_input.format(float_input)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


def return_key(dictionary,string,default_value):
    return(dictionary[string] if string in dictionary.keys() else default_value)



### https://stackoverflow.com/questions/51562221/python-multiprocessing-overflowerrorcannot-serialize-a-bytes-object-larger-t/51568084
### https://stackoverflow.com/questions/47776486/python-struct-error-i-format-requires-2147483648-number-2147483647/47776649#47776649
### Solve the pickle issue OverflowError('cannot serialize a bytes object larger than 4GiB')
### When too large array are considered in multiprocessing
### Add the following lines to you code:
### ctx = mp.get_context()
### ctx.reducer = utils.Pickle4Reducer()
### utils.patch_mp_connection_bpo_17560(log = log)
### You might want to change multiprocessing context with: mp.get_context('fork')



class ForkingPickler4(ForkingPickler):
    def __init__(self, *args):
        if len(args) > 1:
            args[1] = 2
        else:
            args.append(2)
        super().__init__(*args)

    @classmethod
    def dumps(cls, obj, protocol=4):
        return ForkingPickler.dumps(obj, protocol)


def dump(obj, file, protocol=4):
    ForkingPickler4(file, protocol).dump(obj)


class Pickle4Reducer(AbstractReducer):
    ForkingPickler = ForkingPickler4
    register = ForkingPickler4.register
    dump = dump


import functools
import struct
import sys



def patch_mp_connection_bpo_17560(log = None):
    """Apply PR-10305 / bpo-17560 connection send/receive max size update

    See the original issue at https://bugs.python.org/issue17560 and
    https://github.com/python/cpython/pull/10305 for the pull request.

    This only supports Python versions 3.3 - 3.7, this function
    does nothing for Python versions outside of that range.

    """
    patchname = "Multiprocessing connection patch for bpo-17560"
    if not (3, 3) < sys.version_info < (3, 8):
        if log is not None:
            log.add(
            patchname + " not applied, not an applicable Python version: %s",
            sys.version
        )
        return

    from multiprocessing.connection import Connection

    orig_send_bytes = Connection._send_bytes
    orig_recv_bytes = Connection._recv_bytes
    if (
        orig_send_bytes.__code__.co_filename == __file__
        and orig_recv_bytes.__code__.co_filename == __file__
    ):
        if log is not None:
            log.add(patchname + " already applied, skipping")
        return

    @functools.wraps(orig_send_bytes)
    def send_bytes(self, buf):
        n = len(buf)
        if n > 0x7fffffff:
            pre_header = struct.pack("!i", -1)
            header = struct.pack("!Q", n)
            self._send(pre_header)
            self._send(header)
            self._send(buf)
        else:
            orig_send_bytes(self, buf)

    @functools.wraps(orig_recv_bytes)
    def recv_bytes(self, maxsize=None):
        buf = self._recv(4)
        size, = struct.unpack("!i", buf.getvalue())
        if size == -1:
            buf = self._recv(8)
            size, = struct.unpack("!Q", buf.getvalue())
        if maxsize is not None and size > maxsize:
            return None
        return self._recv(size)

    Connection._send_bytes = send_bytes
    Connection._recv_bytes = recv_bytes

    if log is not None:
        log.add(patchname + " applied")





def hist_profile(x, y, bins, range_x,range_y,outlier_insensitive=False):
    w = (y>range_y[0]) & (y<range_y[1])
    if(outlier_insensitive):
        means_result = binned_statistic(x[w], y[w], bins=bins, range=range_x, statistic='median')
        outlier_insensitive_std = lambda x : (np.nanpercentile(x,84.135,axis=0)-np.nanpercentile(x,15.865,axis=0))/2
        std_result = binned_statistic(x[w], y[w], bins=bins, range=range_x, statistic=outlier_insensitive_std)
        nb_entries_result = binned_statistic(x[w], y[w], bins=bins, range=range_x, statistic='count')

    else:
        means_result = binned_statistic(x[w], y[w], bins=bins, range=range_x, statistic='mean')
        std_result = binned_statistic(x[w], y[w], bins=bins, range=range_x, statistic='std')
        nb_entries_result = binned_statistic(x[w], y[w], bins=bins, range=range_x, statistic='count')


    means = means_result.statistic
    std = std_result.statistic
    nb_entries = nb_entries_result.statistic

    errors = std/np.sqrt(nb_entries)

    bin_edges = means_result.bin_edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
    return(bin_centers, means, errors)



#### PLOT ####

def plot_histo(value,value_name,name,**kwargs):
    nb_bins = return_key(kwargs,f"{value_name}_bins",50)
    value_min = return_key(kwargs,f"{value_name}_value_min",None)
    value_max = return_key(kwargs,f"{value_name}_value_max",None)
    alpha = return_key(kwargs,f"{value_name}_alpha",1.0)
    histtype=return_key(kwargs,f"{value_name}_histtype",'bar')
    linestyle=return_key(kwargs,f"{value_name}_linestyle",None)
    ec=return_key(kwargs,f"{value_name}_color",None)

    norm = return_key(kwargs,f"{value_name}_norm",False)
    cumulative = return_key(kwargs,f"{value_name}_cumulative",False)
    log = return_key(kwargs,f"{value_name}_log",False)

    if(norm):
        name = name + "_normalized"
    if(cumulative):
        name = name + "_cumulative"

    if((value_min is None)|(value_min is None)):
        bins = nb_bins
    else:
        bins = np.linspace(value_min,value_max, nb_bins)
    (n, bins, patches) = plt.hist(value, bins, alpha=alpha,histtype=histtype,
                                  linestyle=linestyle,ec=ec,
                                  density=norm,cumulative=cumulative,
                                  log=log)

    return(name, n, bins, patches)

def save_histo(pwd,value,value_name,name,comparison=None,comparison_legend=None,**kwargs):
    xlabel = return_key(kwargs,f"{value_name}_xlabel",value_name)
    ylabel = return_key(kwargs,f"{value_name}_ylabel","#")
    min_lim = return_key(kwargs,f"{value_name}_min_lim",None)
    max_lim = return_key(kwargs,f"{value_name}_max_lim",None)
    plt.figure()
    name_out, n, bins, patches = plot_histo(value,value_name,name,**kwargs)
    if(comparison is not None):
        for i in range(len(comparison)):
            plot_histo(comparison[i],value_name,name,**kwargs)
        if(comparison_legend is not None):
            plt.legend(comparison_legend)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if((min_lim is not None)&(max_lim is not None)):
        plt.xlim([min_lim,max_lim])
    plt.savefig(os.path.join(pwd,f"{name_out}_histo_{value_name}.pdf"), format ="pdf")





def plot_mean_redshift_dependence(value,redshift,value_name,name,**kwargs):
    nb_bins = return_key(kwargs,f"{value_name}_z_bins",50)
    ls=return_key(kwargs,f"{value_name}_linestyle",None)
    color=return_key(kwargs,f"{value_name}_color",None)
    marker=return_key(kwargs,f"{value_name}_marker",'.')
    lambda_obs=return_key(kwargs,f"{value_name}_lambda_obs",False)
    lambda_rest=return_key(kwargs,f"{value_name}_lambda_rest",False)

    outlier_insensitive = return_key(kwargs,f"{value_name}_outlier_insensitive",False)

    if(lambda_obs):
        redshift = (1 + redshift)*lambdaLy

    elif(lambda_rest):
        redshift_qso=return_key(kwargs,f"{value_name}_redshift_qso",None)
        if(redshift_qso is None):
            raise ValueError("""You have chosen to convert your mean redshift
                                dependence plot rest frame wavelength but no
                                QSO redshift was provided""")
        redshift = ((1 + redshift)/(1+ redshift_qso))*lambdaLy

    z_min = return_key(kwargs,f"{value_name}_zmin",np.min(redshift))
    z_max = return_key(kwargs,f"{value_name}_zmax",np.max(redshift))
    range_x = np.array([z_min,z_max])

    bin_centers, means, errors = hist_profile(redshift, value, nb_bins, range_x,
                                              [np.min(value),np.max(value)],
                                              outlier_insensitive=outlier_insensitive)

    plt.errorbar(bin_centers,means, errors,color=color,ls=ls, marker=marker)

    if(outlier_insensitive):
        name = name + "_outlier_insensitive"
        
    return(name)



def save_mean_redshift_dependence(pwd,value,redshift,value_name,name,
                                  comparison=None,
                                  comparison_redshift=None,
                                  comparison_legend=None,**kwargs):
    ylabel = return_key(kwargs,f"{value_name}_xlabel",value_name)
    lambda_obs=return_key(kwargs,f"{value_name}_lambda_obs",False)
    lambda_rest=return_key(kwargs,f"{value_name}_lambda_rest",False)
    if((lambda_obs)&(lambda_rest)):
        raise ValueError(f"""You have chosen to convert your mean redshift
                             dependence plot to both rest frame and observed
                             wavelength. Please choose {value_name}_lambda_obs
                             or {value_name}_lambda_rest""")
    if(lambda_obs):
        default_xlabel = "observed wavelength"
    elif(lambda_rest):
        default_xlabel = "rest frame wavelength"
    else:
        default_xlabel = "redshift"
    xlabel = return_key(kwargs,f"{value_name}_ylabel",default_xlabel)

    plt.figure()
    name_out = plot_mean_redshift_dependence(value,redshift,value_name,name,**kwargs)
    if(comparison is not None):
        for i in range(len(comparison)):
            plot_mean_redshift_dependence(comparison[i],comparison_redshift[i],value_name,name,**kwargs)
        if(comparison_legend is not None):
            plt.legend(comparison_legend)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.savefig(os.path.join(pwd,f"{name_out}_mean_redshift_dependence_{value_name}.pdf"), format ="pdf")



def plot_redshift_dependence(value,redshift,value_name,name,**kwargs):
    z_min = return_key(kwargs,f"{value_name}_zmin",np.min(redshift))
    z_max = return_key(kwargs,f"{value_name}_zmax",np.max(redshift))

    mask = (redshift >= z_min)&(redshift < z_max)
    plt.scatter(redshift[mask],value[mask])
    return(name)



def save_redshift_dependence(pwd,value,redshift,value_name,name,
                             comparison=None,
                             comparison_redshift=None,
                             comparison_legend=None,**kwargs):
    ylabel = return_key(kwargs,f"{value_name}_xlabel",value_name)
    xlabel = return_key(kwargs,f"{value_name}_ylabel","redshift")

    plt.figure()
    name_out = plot_redshift_dependence(value,redshift,value_name,name,**kwargs)
    if(comparison is not None):
        for i in range(len(comparison)):
            plot_redshift_dependence(comparison[i],comparison_redshift[i],value_name,name,**kwargs)
        if(comparison_legend is not None):
            plt.legend(comparison_legend)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.savefig(os.path.join(pwd,f"{name_out}_redshift_dependence_{value_name}.pdf"), format ="pdf")



def plot_ra_dec(ra,dec,name,**kwargs):

    nb_cut = return_key(kwargs,"nb_cut",None)
    plt.figure(figsize=(7,3.5))

    if(nb_cut is not None):
        ramax,ramin,decmax,decmin = np.max(ra),np.min(ra),np.max(dec),np.min(dec)
        interval_ra_array = []
        for cut in range(nb_cut):
            interval_ra_array.append([((cut)/(nb_cut))*(ramax-ramin) + ramin  , ((cut + 1)/(nb_cut))*(ramax-ramin) + ramin])
        for i in range(len(interval_ra_array)):
            plt.plot([interval_ra_array[i][0],interval_ra_array[i][1]],[decmin,decmin],color="orange", linewidth=2)
            plt.plot([interval_ra_array[i][0],interval_ra_array[i][1]],[decmax,decmax],color="orange", linewidth=2)
            plt.plot([interval_ra_array[i][0],interval_ra_array[i][0]],[decmin,decmax],color="orange", linewidth=2)
            plt.plot([interval_ra_array[i][1],interval_ra_array[i][1]],[decmin,decmax],color="orange", linewidth=2)

    plt.plot(ra,dec,'b.', markersize=1.5)
    return(name)


def save_ra_dec(pwd,ra,dec,name,
                comparison_ra=None,comparison_dec=None,
                comparison_legend=None,**kwargs):
    ra_min_lim = return_key(kwargs,"ra_dec_ra_min_lim",np.min(ra))
    ra_max_lim = return_key(kwargs,"ra_dec_ra_max_lim",np.max(ra))
    dec_min_lim = return_key(kwargs,"ra_dec_dec_min_lim",np.min(dec))
    dec_max_lim = return_key(kwargs,"ra_dec_dec_max_lim",np.max(dec))
    deg = return_key(kwargs,"deg",True)
    if(deg):
        label_angle = "deg"
    else:
        label_angle = "rad"
    ylabel = return_key(kwargs,"ra_dec_xlabel",f"DEC [{label_angle}] (J2000)")
    xlabel = return_key(kwargs,"ra_dec_ylabel",f"RA [{label_angle}] (J2000)")

    name_out = plot_ra_dec(ra,dec,name,**kwargs)
    if(comparison_ra is not None):
        for i in range(len(comparison_ra)):
            plot_redshift_dependence(comparison_ra[i],comparison_dec[i],name,**kwargs)
        if(comparison_legend is not None):
            plt.legend(comparison_legend)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([ra_min_lim,ra_max_lim])
    plt.ylim([dec_min_lim,dec_max_lim])
    plt.grid()
    plt.savefig(os.path.join(pwd,f"{name_out}_RA-DEC_diagram.pdf"), format ="pdf")
