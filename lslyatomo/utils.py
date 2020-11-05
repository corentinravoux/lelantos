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
import logging, time
import scipy
import fitsio
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
try :
    from SaclayMocks import box as saclay_mock_box
    from SaclayMocks import constant as saclay_mock_constant
except :
    from lslyatomo.saclaymocks import box as saclay_mock_box
    from lslyatomo.saclaymocks import constant as saclay_mock_constant
    raise Warning("SaclayMocks might be updated, we suggest to install SaclayMocks independently")


try:
    from picca import constants
except:
    import lsstomo.picca.constants as constants
    raise Warning("Picca might be updated, we suggest to install picca independently")


#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################

lambdaLy = 1215.673123130217


def get_cosmo_function(Omega_m,Omega_k=0.):
    Cosmo = constants.cosmo(Omega_m,Ok=Omega_k)
    rcomov = Cosmo.r_comoving
    distang = Cosmo.dm
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
    if(method == "full_angle"):
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

def convert_z_cartesian_to_sky_middle(Z,inv_rcomov,middle_z):
    z = inv_rcomov(Z)
    return(z)

def convert_z_sky_to_cartesian_middle(z,rcomov,middle_z):
    Z = rcomov(z)
    return(Z)




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


def gaussian_smoothing(mapdata,sigma):
    gaussian_map = gaussian_filter(mapdata,sigma)
    return(gaussian_map)



class gaussian_fitter_2d(object):

    def __init__(self,inpdata):

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

        p, success = scipy.optimize.leastsq(errfun, ip)

        return p,success



def create_log(name="Python_Report",log_level="info"):
    log = Logger(name=name,log_level=log_level)
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


saclay_mock_lines_per_box = {"box" : 5,"vx" : 5,"vy" : 5,"vz" : 5,"eta_xx":1,"eta_xy":1,"eta_xz":1,"eta_yx":1,"eta_yy":1,"eta_yz":1,"eta_zx":1,"eta_zy":1,"eta_zz":1}

def saclay_mock_box_cosmo_parameters(box_shape,size_cell):
    import cosmolopy.distance as dist
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
    box = fitsio.FITS(box_dir + "/" + name)[0][:,:,:][0]
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
    X,Y,Z = saclay_mock_box.ComputeXYZ2(ra*(np.pi/180),dec*(np.pi/180),R,ra0*(np.pi/180),dec0*(np.pi/180))
    return(X,Y,Z)
