#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

Author: Corentin Ravoux

Description : routines to extract a DM density box from a Saclay mocks output.
The format of the output box corresponds to the one of a given Tomographic map.
Tested on cori (NERSC)
"""



#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################



import fitsio,os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize
from scipy import interpolate
from scipy.optimize import fsolve
from scipy.ndimage.filters import gaussian_filter
import lsstomo.tomography as tomography
import lsstomo.gaussianfitter as gaussianfitter
from scipy.ndimage import map_coordinates
try :
    from SaclayMocks import box
    from SaclayMocks import constant
except :
    import lsstomo.saclaymocks.box as box
    import lsstomo.saclaymocks.constant as constant
    raise Warning("SaclayMocks might be updated, we suggest to install SaclayMocks independently")
try:
#    import picca.constants as constants_picca
    import lsstomo.picca.constants as constants_picca
except:
    import lsstomo.picca.constants as constants_picca
    raise Warning("Picca might be updated, we suggest to install picca independently")



#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################



class Treat_box():


    def __init__(self,pwd,box_dir,box_shape,size_cell,box_bound,master_file,nb_processors=1,interpolation_method="LINEAR"):

        self.pwd = pwd
        self.box_dir = box_dir
        self.box_shape = box_shape
        self.size_cell = size_cell
        self.box_bound = box_bound
        self.master_file = master_file
        self.nb_processors=nb_processors
        self.lines_per_box = {"vx" : 5,"vy" : 5,"vz" : 5,"eta_xx":1,"eta_xy":1,"eta_xz":1,"eta_yx":1,"eta_yy":1,"eta_yz":1,"eta_zx":1,"eta_zy":1,"eta_zz":1}
        self.interpolation_method = interpolation_method



    def create_Report(self):
        currentdir = os.getcwd()
        os.chdir(self.pwd)
        f = open("Python_Report","w")
        f.write("Report of the Python program\n")
        f.close()
        os.chdir(currentdir)


    def add_Report(self,line):
        currentdir = os.getcwd()
        os.chdir(self.pwd)
        f = open("Python_Report","a")
        f.write(line + "\n")
        f.close()
        os.chdir(currentdir)



    def compute_redshift(self,minredshift,z,rcomov):
        redshift = np.zeros(z.shape)
        for i in range(len(z)):
            redshift[i] = fsolve(self.f,minredshift,args=(z[i],rcomov))
        return(redshift)

    def f(self,redshift,z,rcomov):
        return(rcomov(redshift) - (z))


    def create_my_cosmo_function(self,Om):
        Cosmo = constants_picca.cosmo(Om)
        rcomoving = Cosmo.r_comoving
        rcomov = rcomoving
        distang = Cosmo.dm
        return(distang,rcomov,Cosmo)



    def convert_array_my_cosmo(self,X_tomo_array,Y_tomo_array,Z_tomo_array,distang,rcomov,minredshift,angular_distance_redshift):
        z_array = self.compute_redshift(minredshift,Z_tomo_array,rcomov)
        dist_angular = distang(angular_distance_redshift)
        ra_array = (X_tomo_array / dist_angular)*(180/np.pi)
        dec_array = (Y_tomo_array / dist_angular)*(180/np.pi)
        return(ra_array,dec_array,z_array)




    def get_box_cosmo_parameters(self):
        import cosmolopy.distance as dist
        NZ = self.box_shape[2]
        DZ = self.size_cell
        LZ = NZ*DZ

        h = constant.h
        Om = constant.omega_M_0
        OL = constant.omega_lambda_0
        Ok = constant.omega_k_0
        z0 = constant.z0

        cosmo_fid = {'omega_M_0':Om, 'omega_lambda_0':OL, 'omega_k_0':Ok, 'h':h}
        R_of_z, z_of_R = dist.quick_distance_function(dist.comoving_distance, return_inverse=True, **cosmo_fid)
        R0 = h * R_of_z(z0)
        Rmin = R0 - LZ/2
        Rmax = R0 + LZ/2
        return(R0,z0,R_of_z,z_of_R,Rmin,Rmax,h)


    def get_box_cosmo_parameters_of(self,i_box):
        import cosmolopy.distance as dist
        NZ = self.box_shape[i_box][2]
        DZ = self.size_cell
        LZ = NZ*DZ

        h = constant.h
        Om = constant.omega_M_0
        OL = constant.omega_lambda_0
        Ok = constant.omega_k_0
        z0 = constant.z0

        cosmo_fid = {'omega_M_0':Om, 'omega_lambda_0':OL, 'omega_k_0':Ok, 'h':h}
        R_of_z, z_of_R = dist.quick_distance_function(dist.comoving_distance, return_inverse=True, **cosmo_fid)
        R0 = h * R_of_z(z0)
        Rmin = R0 - LZ/2
        Rmax = R0 + LZ/2
        return(R0,z0,R_of_z,z_of_R,Rmin,Rmax,h)




    def get_box_cosmo_parameters_master(self):
        NZ = self.box_shape[2]
        DZ = self.size_cell
        LZ = NZ*DZ
        Z = fitsio.FITS(self.master_file)[2]["Z"][:]
        R = fitsio.FITS(self.master_file)[2]["DC"][:]
        R_of_z = interpolate.interp1d(Z,R)
        z_of_R = interpolate.interp1d(R,Z)
        h = constant.h
        z0 = constant.z0
        R0 = h * R_of_z(z0)
        Rmin = (R0-LZ/2)
        Rmax = (R0+LZ/2)
        return(R0,z0,R_of_z,z_of_R,Rmin,Rmax,h)




    def get_XYZ_deg(self,ra,dec,R,ra0,dec0):
        '''
        XYZ of a point P (ra,dec,R) in a frame with
        observer at O, Z along OP, X along ra0, Y along dec0
        angles in radians
        tested that ra,dec, R = box.ComputeRaDecR(R0,ra0,dec0,X,Y,Z)
        x,y,z = box.ComputeXYZ(ra[0],dec[0],R,ra0,dec0)
        print x-X,y-Y,z-R0-Z        prints ~1E-13  for random inputs
        '''
        X,Y,Z = box.ComputeXYZ2(ra*(np.pi/180),dec*(np.pi/180),R,ra0*(np.pi/180),dec0*(np.pi/180))
        return(X,Y,Z)






    def create_box_array(self,X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map):
        X_tomo_array = np.linspace(X_tomo_min,X_tomo_max,shape_map[0])
        Y_tomo_array = np.linspace(Y_tomo_min,Y_tomo_max,shape_map[1])
        Z_tomo_array = np.linspace(Z_tomo_min,Z_tomo_max,shape_map[2])
        return(X_tomo_array,Y_tomo_array,Z_tomo_array)


    def get_center_of_the_box(self):
        ra0_box = (self.box_bound[0] + self.box_bound[1])/2
        dec0_box = (self.box_bound[2] + self.box_bound[3])/2
        return(ra0_box, dec0_box)

    def get_center_of_the_box_i(self,i_box):
        ra0_box = (self.box_bound[i_box][0] + self.box_bound[i_box][1])/2
        dec0_box = (self.box_bound[i_box][2] + self.box_bound[i_box][3])/2
        return(ra0_box, dec0_box)

    def read_box(self,n_x):
        name = "box-{}.fits".format(str(n_x))
        box = fitsio.FITS(self.box_dir + "/" + name)[0][:,:,:][0]
        return(box)

    def read_box_prop(self,str_prop,n_x):
        name = "{}-{}.fits".format(str_prop,str(n_x))
        box = fitsio.FITS(self.box_dir + "/" + name)[0][:,:,:][0]
        return(box)

    def get_DM_map_mocks(self):
        DM_mocks = np.zeros((self.box_shape[0],self.box_shape[1],self.box_shape[2]))
        for i in range(self.box_shape[0]):
            DM_mocks[i,:,:] = self.read_box(i)[:,:]
        return(DM_mocks)

    def get_prop_map_mocks(self,str_prop):
        line_per_box = self.lines_per_box[str_prop]
        DM_mocks = np.zeros((self.box_shape[0],self.box_shape[1],self.box_shape[2]))
        for i in range(self.box_shape[0]//line_per_box):
            DM_mocks[i*line_per_box:(i+1)*line_per_box,:,:] = self.read_box_prop(str_prop,i)[:,:]
        return(DM_mocks)

    def read_box_i(self,n_x,i_box):
        name = "box-{}.fits".format(str(n_x))
        box = fitsio.FITS(self.box_dir[i_box] + "/" + name)[0][:,:,:][0]
        return(box)

    def read_box_prop_i(self,str_prop,n_x,i_box):
        name = "{}-{}.fits".format(str_prop,str(n_x))
        box = fitsio.FITS(self.box_dir[i_box] + "/" + name)[0][:,:,:][0]
        return(box)

    def get_DM_map_mocks_i(self,i_box):
        DM_mocks = np.zeros((self.box_shape[i_box][0],self.box_shape[i_box][1],self.box_shape[i_box][2]))
        for i in range(self.box_shape[i_box][0]):
            DM_mocks[i,:,:] = self.read_box_i(i,i_box)[:,:]
        return(DM_mocks)

    def get_prop_map_mocks_i(self,str_prop,i_box):
        line_per_box = self.lines_per_box[str_prop]
        DM_mocks = np.zeros((self.box_shape[i_box][0],self.box_shape[i_box][1],self.box_shape[i_box][2]))
        for i in range(self.box_shape[i_box][0]//line_per_box):
            DM_mocks[i*line_per_box:(i+1)*line_per_box,:,:] = self.read_box_prop_i(str_prop,i,i_box)[:,:]
        return(DM_mocks)


    def gaussianSmoothing(self,mapdata,sigma):
        gaussianMap = gaussian_filter(mapdata,sigma)
        return(gaussianMap)


    def interpolate_dm_map(self,DM_mocks_map,coords_pixels_box_saclay):
        if(self.interpolation_method.upper() == "NEAREST"):
            coords_pixels_box_saclay = coords_pixels_box_saclay.astype(int)
            DM_map = DM_mocks_map[coords_pixels_box_saclay[:,:,:,0],coords_pixels_box_saclay[:,:,:,1],coords_pixels_box_saclay[:,:,:,2]]
        elif(self.interpolation_method.upper() == "LINEAR"):
#            from scipy.interpolate import RegularGridInterpolator
#            x,y,z = np.arange(DM_mocks_map.shape[0]),np.arange(DM_mocks_map.shape[1]),np.arange(DM_mocks_map.shape[2])
#            my_interpolating_function = RegularGridInterpolator((x, y, z), DM_mocks_map,method="linear")
#            points =  coords_pixels_box_saclay.reshape(coords_pixels_box_saclay.shape[0]*coords_pixels_box_saclay.shape[1]*coords_pixels_box_saclay.shape[2],3)
#            DM_map = my_interpolating_function(points).reshape((coords_pixels_box_saclay.shape[0],coords_pixels_box_saclay.shape[1],coords_pixels_box_saclay.shape[2]))
            points =  coords_pixels_box_saclay.reshape(coords_pixels_box_saclay.shape[0]*coords_pixels_box_saclay.shape[1]*coords_pixels_box_saclay.shape[2],3)
            DM_map = map_coordinates(DM_mocks_map, np.transpose(points), order=1).reshape((coords_pixels_box_saclay.shape[0],coords_pixels_box_saclay.shape[1],coords_pixels_box_saclay.shape[2]))
        elif(self.interpolation_method.upper() == "SPLINE"):
            points =  coords_pixels_box_saclay.reshape(coords_pixels_box_saclay.shape[0]*coords_pixels_box_saclay.shape[1]*coords_pixels_box_saclay.shape[2],3)
            DM_map = map_coordinates(DM_mocks_map, np.transpose(points), order=2).reshape((coords_pixels_box_saclay.shape[0],coords_pixels_box_saclay.shape[1],coords_pixels_box_saclay.shape[2]))
        else :
            raise ValueError("Please select NEAREST, LINEAR or SPLINE as interpolation_method")
        return(DM_map)



    def construct_DM_map(self,ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,get_prop=None):
        self.add_Report("Creation of (RA,DEC,R) data matrix")
        coords_ra_dec =np.moveaxis(np.array(np.meshgrid(ra_array,dec_array,h * R_of_z(z_array),indexing='ij')),0,-1)
        self.add_Report("Conversion to (X,Y,Z) coordinates in the Saclay box")
        coords_box_saclay = np.zeros(coords_ra_dec.shape)
        coords_box_saclay[:,:,:,0],coords_box_saclay[:,:,:,1],coords_box_saclay[:,:,:,2] = self.get_XYZ_deg(coords_ra_dec[:,:,:,0],coords_ra_dec[:,:,:,1],coords_ra_dec[:,:,:,2],ra0_box,dec0_box)
        del coords_ra_dec
        self.add_Report("Searching for the nearest pixels of the Saclay box")
        coords_pixels_box_saclay = np.zeros(coords_box_saclay.shape)
        coords_pixels_box_saclay[:,:,:,0],coords_pixels_box_saclay[:,:,:,1],coords_pixels_box_saclay[:,:,:,2] = self.get_coord_DM_map(coords_box_saclay[:,:,:,0],coords_box_saclay[:,:,:,1],coords_box_saclay[:,:,:,2],Rmin)
        del coords_box_saclay
        self.add_Report("Loading of the Saclay map")
        DM_mocks_map = self.get_DM_map_mocks()
        self.add_Report("Creation of the DM map")
#        DM_map = np.zeros(tuple(shape_map_output))
#        DM_map[:,:,:] = DM_mocks_map[coords_pixels_box_saclay[:,:,:,0],coords_pixels_box_saclay[:,:,:,1],coords_pixels_box_saclay[:,:,:,2]]
        DM_map = self.interpolate_dm_map(DM_mocks_map,coords_pixels_box_saclay)
        del DM_mocks_map
        if(get_prop is not None):
            Props_map = []
            for i in range(len(get_prop)):
                self.add_Report("Loading of the Saclay {} map".format(get_prop[i]))
                prop_mocks_map = self.get_prop_map_mocks(get_prop[i])
                self.add_Report("Creation of the {} map".format(get_prop[i]))
#                prop_map = np.zeros(tuple(shape_map_output))
#                prop_map[:,:,:] = prop_mocks_map[coords_pixels_box_saclay[:,:,:,0],coords_pixels_box_saclay[:,:,:,1],coords_pixels_box_saclay[:,:,:,2]]
                prop_map = self.interpolate_dm_map(prop_mocks_map,coords_pixels_box_saclay)
                Props_map.append(prop_map)
                del prop_mocks_map
            Props_map = np.array(Props_map)
        else:
            Props_map = None
        del coords_pixels_box_saclay
        self.add_Report("Multiplying by growth factor at redshift of the LOS")
        return(DM_map,Props_map)



    def interpolate_fill_dm_map(self,DM_map,DM_mocks_map,coords_pixels_box_saclay,mask):
        if(self.interpolation_method.upper() == "NEAREST"):
            coords_pixels_box_saclay = coords_pixels_box_saclay.astype(int)
            DM_map[:,:,:][mask] = DM_mocks_map[coords_pixels_box_saclay[:,:,:,0][mask],coords_pixels_box_saclay[:,:,:,1][mask],coords_pixels_box_saclay[:,:,:,2][mask]]
        elif(self.interpolation_method.upper() == "LINEAR"):
            points =  coords_pixels_box_saclay[mask]
            DM_map[:,:,:][mask] = map_coordinates(DM_mocks_map, np.transpose(points), order=1)
        elif(self.interpolation_method.upper() == "SPLINE"):
            points =  coords_pixels_box_saclay[mask]
            DM_map[:,:,:][mask] = map_coordinates(DM_mocks_map, np.transpose(points), order=2)
        else :
            raise ValueError("Please select NEAREST, LINEAR or SPLINE as interpolation_method")



    def fill_DM_map(self,ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,DM_map,i_box,get_prop=None,Props_map=None):
        self.add_Report("Creation of (RA,DEC,R) data matrix")
        coords_ra_dec =np.moveaxis(np.array(np.meshgrid(ra_array,dec_array,h * R_of_z(z_array),indexing='ij')),0,-1)
        mask1 = (DM_map[:,:,:] == None)
        mask2 = (coords_ra_dec[:,:,:,0]>=self.box_bound[i_box][0])
        mask2 &=(coords_ra_dec[:,:,:,0]<self.box_bound[i_box][1])
        mask2 &=(coords_ra_dec[:,:,:,1]>=self.box_bound[i_box][2])
        mask2 &=(coords_ra_dec[:,:,:,1]<self.box_bound[i_box][3])
        if((len(mask1[mask1==True]) == 0)|(len(mask2[mask2==True])==0)):
            self.add_Report("No need of the Saclay box {}".format(i_box))
            return(DM_map,Props_map)
        mask = mask1&mask2
        del mask1,mask2
        self.add_Report("Conversion to (X,Y,Z) coordinates in the Saclay box {}".format(i_box))
        coords_box_saclay = np.zeros(coords_ra_dec.shape)
        coords_box_saclay[:,:,:,0],coords_box_saclay[:,:,:,1],coords_box_saclay[:,:,:,2] = self.get_XYZ_deg(coords_ra_dec[:,:,:,0],coords_ra_dec[:,:,:,1],coords_ra_dec[:,:,:,2],ra0_box,dec0_box)
        del coords_ra_dec
        self.add_Report("Searching for the nearest pixels of the Saclay box {}".format(i_box))
        coords_pixels_box_saclay = np.zeros(coords_box_saclay.shape)
        coords_pixels_box_saclay[:,:,:,0],coords_pixels_box_saclay[:,:,:,1],coords_pixels_box_saclay[:,:,:,2] = self.get_coord_DM_map_i(coords_box_saclay[:,:,:,0],coords_box_saclay[:,:,:,1],coords_box_saclay[:,:,:,2],Rmin,i_box)
        del coords_box_saclay
        self.add_Report("Loading of the Saclay map {}".format(i_box))
        DM_mocks_map = self.get_DM_map_mocks_i(i_box)
        self.add_Report("Creation of the DM map with box {}".format(i_box))
#        DM_map[:,:,:][mask] = DM_mocks_map[coords_pixels_box_saclay[:,:,:,0][mask],coords_pixels_box_saclay[:,:,:,1][mask],coords_pixels_box_saclay[:,:,:,2][mask]]
        self.interpolate_fill_dm_map(DM_map,DM_mocks_map,coords_pixels_box_saclay,mask)
        del DM_mocks_map
        if(get_prop is not None):
            for i in range(len(get_prop)):
                self.add_Report("Loading of the Saclay {} map {}".format(get_prop[i],i_box))
                prop_mocks_map = self.get_prop_map_mocks_i(get_prop[i],i_box)
                self.add_Report("Creation of the {} map with box {}".format(get_prop[i],i_box))
#                Props_map[i][:,:,:][mask] = prop_mocks_map[coords_pixels_box_saclay[:,:,:,0][mask],coords_pixels_box_saclay[:,:,:,1][mask],coords_pixels_box_saclay[:,:,:,2][mask]]
                self.interpolate_fill_dm_map(Props_map[i],prop_mocks_map,coords_pixels_box_saclay,mask)
                del prop_mocks_map
        del coords_pixels_box_saclay,mask
        return(DM_map,Props_map)



    def construct_DM_LOS(self,ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,Rmin,h):
        self.add_Report("Creation of (RA,DEC,R) data matrix")
        R_array = h * R_of_z(z_array)
        self.add_Report("Conversion to (X,Y,Z) coordinates in the Saclay box")
        X,Y,Z = np.zeros(len(R_array)),np.zeros(len(R_array)),np.zeros(len(R_array))
        X,Y,Z = self.get_XYZ_deg(ra_array,dec_array,R_array,ra0_box,dec0_box)
        self.add_Report("Searching for the nearest pixels of the Saclay box")
        i,j,k = np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int)
        i,j,k = self.get_coord_DM_map(X,Y,Z,Rmin)
        del X,Y,Z
        self.add_Report("Loading of the Saclay map")
        DM_mocks_map = self.get_DM_map_mocks()
        self.add_Report("Creation of the DM los")
        DM_LOS = np.zeros(len(R_array))
        DM_LOS[:] = DM_mocks_map[i[:],j[:],k[:]]
        del DM_mocks_map,i,j,k
        self.add_Report("Multiplying by growth factor at redshift of the LOS")
        DM_LOS = self.multiply_los_by_growth(DM_LOS,z_array)
        return(DM_LOS)

    def fill_DM_LOS(self,ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,Rmin,h,DM_LOS,i_box):
        self.add_Report("Creation of (RA,DEC,R) data matrix")
        R_array = h * R_of_z(z_array)
        mask1 = (DM_LOS[:] == None)
        mask2 = (ra_array[:]>=self.box_bound[i_box][0])
        mask2 &=(ra_array[:]<self.box_bound[i_box][1])
        mask2 &=(dec_array[:]>=self.box_bound[i_box][2])
        mask2 &=(dec_array[:]<self.box_bound[i_box][3])
        if((len(mask1[mask1==True]) == 0)|(len(mask2[mask2==True])==0)):
            self.add_Report("No need of the Saclay box {}".format(i_box))
            return(DM_LOS)
        mask = mask1&mask2
        del mask1,mask2
        self.add_Report("Conversion to (X,Y,Z) coordinates in the Saclay box {}".format(i_box))
        X,Y,Z = np.zeros(len(R_array)),np.zeros(len(R_array)),np.zeros(len(R_array))
        X,Y,Z = self.get_XYZ_deg(ra_array,dec_array,R_array,ra0_box,dec0_box)
        self.add_Report("Searching for the nearest pixels of the Saclay box {}".format(i_box))
        i,j,k = np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int)
        i,j,k = self.get_coord_DM_map_i(X,Y,Z,Rmin,i_box)
        del X,Y,Z
        self.add_Report("Loading of the Saclay map {}".format(i_box))
        DM_mocks_map = self.get_DM_map_mocks_i(i_box)
        self.add_Report("Creation of the DM los with box {}".format(i_box))
        DM_LOS[:][mask] = DM_mocks_map[i[:][mask],j[:][mask],k[:][mask]]
        del DM_mocks_map,i,j,k
        return(DM_LOS)





    def compute_rsd(self,DM_map,Props_map,z_array,Cosmo):
        for i in range(len(z_array)):
            DM_map[:,:,i] = DM_map[:,:,i]  - (((1+z_array[i])* Props_map[0][:,:,i])/Cosmo.hubble(z_array[i]))
        return(DM_map)



    def multiply_by_growth(self,DM_map,z_array):
        Z = fitsio.FITS(self.master_file)[2]["Z"][:]
        G =fitsio.FITS(self.master_file)[2]["G"][:]
        G_of_Z =interpolate.interp1d(Z,G)
        for i in range(len(z_array)):
            DM_map[:,:,i] = DM_map[:,:,i] * G_of_Z(z_array[i])
        return(DM_map)

    def multiply_los_by_growth(self,DM_los,z_array):
        Z = fitsio.FITS(self.master_file)[2]["Z"][:]
        G =fitsio.FITS(self.master_file)[2]["G"][:]
        G_of_Z =interpolate.interp1d(Z,G)
        DM_los = DM_los * G_of_Z(z_array)
        return(DM_los)


    def convert_to_matter_field(self,gaussian_array):
#        sigma_l = constant.sigma_l
        sigma_l = np.std(gaussian_array)
        density_matter = np.exp(gaussian_array - (sigma_l**2/2)) - 1
        return(density_matter)


    def get_coord_DM_map(self,X,Y,Z,Rmin):
        size_cell = self.size_cell
        center_x = (self.box_shape[0]-1) / 2
        center_y = (self.box_shape[1]-1) / 2
        if(self.interpolation_method.upper() == "NEAREST"):
            n_i = (np.round(center_x + X/size_cell,0)).astype(int)
            n_j = (np.round(center_y + Y/size_cell,0)).astype(int)
            n_k = np.round((Z - Rmin)/size_cell,0).astype(int)
        else:
            n_i = center_x + X/size_cell
            n_j = center_y + Y/size_cell
            n_k = (Z - Rmin)/size_cell
        return(n_i,n_j,n_k)

    def get_coord_DM_map_i(self,X,Y,Z,Rmin,i_box):
        size_cell = self.size_cell
        center_x = (self.box_shape[i_box][0]-1) / 2
        center_y = (self.box_shape[i_box][1]-1) / 2
        if(self.interpolation_method.upper() == "NEAREST"):
            n_i = (np.round(center_x + X/size_cell,0)).astype(int)
            n_j = (np.round(center_y + Y/size_cell,0)).astype(int)
            n_k = np.round((Z - Rmin)/size_cell,0).astype(int)
        else:
            n_i = center_x + X/size_cell
            n_j = center_y + Y/size_cell
            n_k = (Z - Rmin)/size_cell
        return(n_i,n_j,n_k)

    def get_rho_DM(self,n_x,n_y,n_z):
        box = self.read_box(n_x)
        delta_rho = box[n_y,n_z]
        return(delta_rho)


    def write_map(self,map_3D,name):
        map_3D.astype(float).ravel().tofile(name)

    def write_los(self,los,name):
        los.astype(float).tofile(name)




    def construct_DM_map_old(self,ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin):
        "obsolete"
        DM_mocks_map = self.get_DM_map_mocks()
        DM_map = np.zeros(tuple(shape_map_output))
        for i in range(len(z_array)):
            R = R_of_z(z_array[i])
            for j in range(len(ra_array)):
                for k in range(len(dec_array)):
                    X,Y,Z = self.get_XYZ_deg(ra_array[j],dec_array[k],R,ra0_box,dec0_box)
                    n_x,n_y,n_z = self.get_coord_DM_map(X,Y,Z,Rmin)
                    DM_map[j,k,i] = DM_mocks_map[n_x,n_y,n_z]
        return(DM_map)


    def construct_DM_map_parallel(self,ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin):
        "obsolete"
        DM_map = np.zeros(tuple(shape_map_output))
        from multiprocessing import Pool
        pool = Pool(self.nb_processors)
        coords = np.zeros((shape_map_output[0],shape_map_output[1],shape_map_output[2],3))
        coords =np.moveaxis(np.array(np.meshgrid(ra_array,dec_array,R_of_z(z_array),indexing='ij')),0,-1)
        self.ra0_box = ra0_box
        self.dec0_box = dec0_box
        self.Rmin = Rmin
        DM_map = pool.map(self.DM_map_parallel,coords)
        del self.ra0_box, self.dec0_box, self.Rmin
        return(DM_map)

    def DM_map_parallel(self,coord):
        "obsolete"
        ra,dec,R = coord[0],coord[1],coord[2]
        X,Y,Z = self.get_XYZ_deg(ra,dec,R,self.ra0_box,self.dec0_box)
        n_x,n_y,n_z = self.get_coord_DM_map(X,Y,Z,self.Rmin)
        delta_rho = self.get_rho_DM(n_x,n_y,n_z)
        return(delta_rho)




########## 1 DM box ##########

    def create_DM_box(self,my_Om,X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output,minredshift,angular_distance_redshift,name,rsd_box=False,growth_multipication=True,matter_field=True):
        """ Main function of the Class"""
        self.create_Report()
        (distang,rcomov,Cosmo) = self.create_my_cosmo_function(my_Om)
        X_tomo_array,Y_tomo_array, Z_tomo_array = self.create_box_array(X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output)
        ra_array,dec_array,z_array = self.convert_array_my_cosmo(X_tomo_array,Y_tomo_array, Z_tomo_array,distang,rcomov,minredshift,angular_distance_redshift)
        (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = self.get_box_cosmo_parameters()
        ra0_box, dec0_box = self.get_center_of_the_box()
        if(rsd_box): get_prop = ["eta_zz"]
        else: get_prop = None
        DM_map,Props_map = self.construct_DM_map(ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,get_prop=get_prop)
        if(growth_multipication): DM_map = self.multiply_by_growth(DM_map,z_array)
        if(matter_field): DM_map = self.convert_to_matter_field(DM_map)
        self.write_map(DM_map,name)
        del DM_map
        if(rsd_box):
            for i in range(len(get_prop)):
                self.write_map(Props_map[i],"{}_{}".format(name,get_prop[i]))
        del Props_map
        self.add_Report("Size and Shape of the reconstructed map :")
        self.add_Report(str((X_tomo_max-X_tomo_min,Y_tomo_max-Y_tomo_min,Z_tomo_max-Z_tomo_min)))
        self.add_Report(str(shape_map_output))



    def create_DM_LOS(self,pixel_name,name,growth_multipication=True,matter_field=True):
        """ Main function of the Class"""
        self.create_Report()
        pixel = np.fromfile(open(pixel_name,'r'),dtype=np.float64)
        pixel = pixel.reshape((len(pixel)//4,4))
        ra_array,dec_array,z_array = np.array(pixel[:,0]),np.array(pixel[:,1]),np.array(pixel[:,2])
        (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = self.get_box_cosmo_parameters()
        ra0_box, dec0_box = self.get_center_of_the_box()
        DM_LOS = self.construct_DM_LOS(ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,Rmin,h)
        self.write_los(DM_LOS,name)


    def create_delta_catalog(self,cat_name,name,name_redshift="Z_QSO_RSD",growth_multipication=True,matter_field=True):
        self.create_Report()
        ra = fitsio.FITS(cat_name)[1]["RA"][:]
        ra[ra>180] = ra[ra>180] -360
        dec = fitsio.FITS(cat_name)[1]["DEC"][:]
        z = fitsio.FITS(cat_name)[1][name_redshift][:]
        (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = self.get_box_cosmo_parameters()
        ra0_box, dec0_box = self.get_center_of_the_box()
        delta_quasars = self.construct_DM_LOS(ra,dec,z,R0,R_of_z,ra0_box,dec0_box,Rmin,h)
        self.write_los(delta_quasars,name)
        self.add_Report("Mean delta extracted : {}".format(np.mean(delta_quasars)))






########## Several box #########



    def create_DM_several_boxes(self,my_Om,X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output,minredshift,angular_distance_redshift,name,rsd_box=False,growth_multipication=True,matter_field=True):
        """ Main function of the Class"""
        self.create_Report()
        (distang,rcomov,Cosmo) = self.create_my_cosmo_function(my_Om)
        X_tomo_array,Y_tomo_array, Z_tomo_array = self.create_box_array(X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output)
        ra_array,dec_array,z_array = self.convert_array_my_cosmo(X_tomo_array,Y_tomo_array, Z_tomo_array,distang,rcomov,minredshift,angular_distance_redshift)
        DM_map = np.full(tuple(shape_map_output),None)
        if(rsd_box):
            get_prop = ["eta_zz","vz"]
            Props_map = np.array([np.full(tuple(shape_map_output),None) for i in range(len(get_prop))])
        else: Props_map,get_prop = None,None
        for i_box in range(len(self.box_dir)):
            (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = self.get_box_cosmo_parameters_of(i_box)
            ra0_box, dec0_box = self.get_center_of_the_box_i(i_box)
            DM_map,Props_map = self.fill_DM_map(ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,DM_map,i_box,get_prop=get_prop,Props_map=Props_map)
        if(len(DM_map[DM_map==None]!=0)):
            self.add_Report("Warning : Not enough boxes to fill the Dark Matter map")
        self.add_Report("Multiplying by growth factor at redshift of the LOS")

        # Good way ? first growth or first density matter conversion ?

        if(growth_multipication): DM_map = self.multiply_by_growth(DM_map,z_array)
        if(matter_field): DM_map = self.convert_to_matter_field(DM_map)
        self.write_map(DM_map,name)
        del DM_map
        if(rsd_box):
            for i in range(len(get_prop)):
                self.write_map(Props_map[i],"{}_{}".format(name,get_prop[i]))
        del Props_map
        self.add_Report("Size and Shape of the reconstructed map :")
        self.add_Report(str((X_tomo_max-X_tomo_min,Y_tomo_max-Y_tomo_min,Z_tomo_max-Z_tomo_min)))
        self.add_Report(str(shape_map_output))



    def create_DM_LOS_several_boxes(self,pixel_name,name,growth_multipication=True,matter_field=True):
        """ Main function of the Class"""
        self.create_Report()
        pixel = np.fromfile(open(pixel_name,'r'),dtype=np.float64)
        pixel = pixel.reshape((len(pixel)//4,4))
        ra_array,dec_array,z_array = np.array(pixel[:,0]),np.array(pixel[:,1]),np.array(pixel[:,2])
        DM_LOS = np.full(z_array.shape,None)
        for i_box in range(len(self.box_dir)):
            (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) =  self.get_box_cosmo_parameters_of(i_box)
            ra0_box, dec0_box = self.get_center_of_the_box_i(i_box)
            DM_LOS = self.fill_DM_LOS(ra_array,dec_array,z_array,R0,R_of_z,ra0_box,dec0_box,Rmin,h,DM_LOS,i_box)
        if(len(DM_LOS[DM_LOS==None]!=0)):
            self.add_Report("Warning : Not enough boxes to fill the Dark Matter map")
        DM_LOS = self.multiply_los_by_growth(DM_LOS,z_array)
        self.write_los(DM_LOS,name)
        self.add_Report("Mean delta extracted : {}".format(np.mean(DM_LOS)))



    def create_delta_catalog_several_boxes(self,cat_name,name,name_redshift="Z_QSO_RSD",growth_multipication=True,matter_field=True):
        self.create_Report()
        ra = fitsio.FITS(cat_name)[1]["RA"][:]
        ra[ra>180] = ra[ra>180] -360
        dec = fitsio.FITS(cat_name)[1]["DEC"][:]
        z = fitsio.FITS(cat_name)[1][name_redshift][:]
        delta_quasars =np.full(z.shape,None)
        for i_box in range(len(self.box_dir)):
            (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) =  self.get_box_cosmo_parameters_of(i_box)
            ra0_box, dec0_box = self.get_center_of_the_box_i(i_box)
            delta_quasars = self.fill_DM_LOS(ra,dec,z,R0,R_of_z,ra0_box,dec0_box,Rmin,h,delta_quasars,i_box)
        if(len(delta_quasars[delta_quasars==None]!=0)):
            self.add_Report("Warning : Not enough boxes to fill the Dark Matter map")
        self.add_Report("Multiplying by growth factor at redshift of the LOS")
        delta_quasars = self.multiply_los_by_growth(delta_quasars,z)
        self.write_los(delta_quasars,name)
        self.add_Report("Mean delta extracted : {}".format(np.mean(delta_quasars)))




######## Main ######




    def create_box(self,my_Om,X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output,minredshift,angular_distance_redshift,name,rsd_box=False,growth_multipication=True,matter_field=True):
        if(type(self.box_dir) == list):
            self.create_DM_several_boxes(my_Om,X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output,minredshift,angular_distance_redshift,name,rsd_box=rsd_box,growth_multipication=growth_multipication,matter_field=matter_field)
        else :
            self.create_DM_box(my_Om,X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output,minredshift,angular_distance_redshift,name,rsd_box=rsd_box,growth_multipication=growth_multipication,matter_field=matter_field)

    def create_LOS(self,pixel_name,name):
        if(type(self.box_dir) == list):
            self.create_DM_LOS_several_boxes(pixel_name,name)
        else :
            self.create_DM_LOS(pixel_name,name)

    def create_catalog(self,cat_name,name,name_redshift="Z_QSO_RSD"):
        if(type(self.box_dir) == list):
            self.create_delta_catalog_several_boxes(cat_name,name,name_redshift=name_redshift)
        else :
            self.create_delta_catalog(cat_name,name,name_redshift=name_redshift)



class Treat_Box(object):

    def __init__(self,pwd,shape_map,box_DM_name,name_map=None,limit=None):

        self.pwd = pwd
        self.name_map = name_map
        self.shape_map = shape_map
        self.box_DM_name = box_DM_name
        self.limit = limit


    def contourplot(self,x,y,ncont=10,colors=None,pltclf=True,binsx=100,binsy=100,log=False) :
        H, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
        # H needs to be rotated and flipped
        H = np.rot90(H)
        H = np.flipud(H)
        xcenters=(xedges[:-1]+xedges[1:])/2.
        ycenters=(yedges[:-1]+yedges[1:])/2.
        if log is True : H[(H>0)]=np.log(H[(H>0)])
        if pltclf : plt.clf()
        plt.contour(xcenters,ycenters,H,ncont,colors=colors,linewidths=2)

    def densityplot(self,x,y,binsx=200,binsy=200,scaleperdeg2=False,maxdensity=0,mindensity=0,pltclf=True,logscale=False) :
        #http://oceanpython.org/2013/02/25/2d-histogram/
        # scaleperdeg2 : for sky map, x and y assumed to be radec
        H, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
        if scaleperdeg2==True :
            delta_ra=abs(xedges[1]-xedges[0])
            delta_dec=abs(yedges[1]-yedges[0])
            cosdelta=np.cos(yedges[:-1]*np.pi/180.) # yedges a (binsy+1) elements
            for k in range(binsy) : H[:,k] = H[:,k]/(delta_ra*delta_dec*cosdelta[k])
        # H needs to be rotated and flipped
        H = np.rot90(H)
        H = np.flipud(H)
        if maxdensity>0 :
            wcut=np.where( (H>maxdensity) )
            H[wcut]=maxdensity
        if mindensity>0 :
            wcut=np.where( (H<mindensity) )
            H[wcut]=0
        Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
        if pltclf : plt.clf()
        if logscale is False :
#            plt.pcolormesh(xedges,yedges,Hmasked, cmap=plt.cm.jet)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            Hmaskednorm = Hmasked/np.max(Hmasked)
            plt.imshow(Hmaskednorm, interpolation='bilinear',cmap='jet',origin="lower",extent = extent,norm=Normalize(vmin=0.0,vmax=1.0))
        else : plt.pcolormesh(xedges,yedges,Hmasked,norm=LogNorm(), cmap=plt.cm.jet)
        cbar = plt.colorbar()
        cbar.set_label("Normalized density")
        if scaleperdeg2==True : plt.title(r'density / deg$^2$')

#    def densityplot(self,x,y,binsx=200,binsy=200,scaleperdeg2=False,maxdensity=0,mindensity=0,pltclf=True,logscale=False) :
#        H, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
#        plt.imshow(H, interpolation='bilinear',cmap='jet_r')#,extent = extentmap )
#        plt.colorbar()

#        plt.imshow(H,norm=LogNorm(), cmap=plt.cm.jet)


    def gaussian_plot(self,X,Y,binsx=200,binsy=200,ncont=4):
        data, xedges,yedges = np.histogram2d(X,Y,bins=(binsx,binsy))
        xcenters=(xedges[:-1]+xedges[1:])/2.
        ycenters=(yedges[:-1]+yedges[1:])/2.
        Fitter = gaussianfitter.gaussian_fitter_2d(data)
        p,success = Fitter.FitGauss2D()
        x,y=np.indices((binsx,binsy),dtype=np.float)
        gauss = Fitter.Gaussian2D(*p)
        data_fitted = gauss(y,x)
        if(type(ncont) != int):
            ncont =np.array(ncont)*p[0]
        plt.contour(xcenters,ycenters,data_fitted.reshape(binsx, binsy),ncont,linewidths=2, colors='w')
        return(p)


    def get_list_boxes(self,gaussian_smoothing=None,cut_redshift_coef=None,cut_pixel=None):
        Tomo = tomography.TreatClamato(self.pwd,self.name_map,self.shape_map,"")
        box_DM = Tomo.readClamatoMapFile_Other(self.box_DM_name,self.shape_map)
        max_x,max_y,max_z =self.shape_map[0],self.shape_map[1],self.shape_map[2]
        if(cut_redshift_coef is not None):
            max_z = self.shape_map[2]//cut_redshift_coef
        if(cut_pixel is not None):
            max_z = cut_pixel
        min_x,min_y,min_z = 0,0,0
        length_list = (max_x-min_x)*(max_y-min_y)*(max_z-min_z)
        box_DM = box_DM[min_x:max_x,min_y:max_y,min_z:max_z]
        if(gaussian_smoothing is not None):
            box_DM = self.gaussianSmoothing(box_DM,gaussian_smoothing)
        list_box = box_DM.reshape(length_list)
        del box_DM
        map_3D = Tomo.readClamatoMapFile()
        map_3D = map_3D[min_x:max_x,min_y:max_y,min_z:max_z]
        list_map = map_3D.reshape(length_list)
        del map_3D, Tomo
        return(list_map,list_box)

    def gaussianSmoothing(self,mapdata,sigma):
        gaussianMap = gaussian_filter(mapdata,sigma)
        return(gaussianMap)

    def plot_gaussian_line(self,gaussian_fit,binsx,binsy):
        xcenter = (gaussian_fit[1] - binsx)*((self.limit[1] - self.limit[0])/binsx)
        ycenter = (gaussian_fit[2] - binsx)*((self.limit[3] - self.limit[2])/binsx)
        xcenter = 0
        ycenter = 0
        angle = gaussian_fit[5]
        pente = 1/np.tan(np.radians(angle))
        x_arrray = np.linspace(self.limit[0],self.limit[1],100)
        y_array = ycenter + pente * (x_arrray - xcenter)
        plt.plot(x_arrray,y_array)


    def plot_scatter_box_DM(self,name,gaussian_smoothing=None,cut_redshift_coef=None,cut_pixel = None,binsx=200,binsy=200,ncont=4,rotate=False):

        (list_map,list_box)= self.get_list_boxes(gaussian_smoothing=gaussian_smoothing,cut_redshift_coef=cut_redshift_coef,cut_pixel=cut_pixel)
        if (gaussian_smoothing is not None):
            if(cut_redshift_coef is not None): name = name + "_smoothing{}_cutpart{}".format(gaussian_smoothing,cut_redshift_coef)
            elif(cut_pixel is not None): name = name + "_smoothing{}_cutpixel{}".format(gaussian_smoothing,cut_pixel)
            else : name = name + "_smoothing{}".format(gaussian_smoothing)
        else:
            if(cut_redshift_coef is not None): name = name + "_cutpart{}".format(cut_redshift_coef)
            elif(cut_pixel is not None): name = name + "_cutpixel{}".format(cut_pixel)
            else : name = name

        #### Analysis ####

        poly_fit = np.polyfit(list_map,list_box,1)
        x_poly = np.linspace(-0.2,0.2,100)
        y_poly = poly_fit[0]*x_poly + poly_fit[1]
        correlation_matrix = np.corrcoef(list_map,list_box)



        #### Plots ####
        xlabel = r"$\delta_{Fmap}$"
        ylabel = r"$\delta_{m}$"
        if rotate : list_box,list_map = list_map, list_box


        plt.figure()
        self.densityplot(list_map,list_box,binsx=binsx,binsy=binsy,scaleperdeg2=False,maxdensity=0,mindensity=0,pltclf=True,logscale=False)
        self.contourplot(list_map,list_box,ncont=ncont,colors="w",pltclf=False,binsx=binsx,binsy=binsy,log=False)
        plt.xlim(self.limit[0:2])
        plt.ylim(self.limit[2:4])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if rotate :
            plt.xlabel(ylabel)
            plt.ylabel(xlabel)
        plt.savefig(name + "_contour.pdf",format="pdf",dpi=100)

        plt.figure()
        self.densityplot(list_map,list_box,binsx=binsx,binsy=binsy,scaleperdeg2=False,maxdensity=0,mindensity=0,pltclf=True,logscale=False)
        gaussian_fit =self.gaussian_plot(list_map,list_box,binsx=binsx,binsy=binsy,ncont = ncont)
#        self.plot_gaussian_line(gaussian_fit,binsx,binsy)
        plt.xlim(self.limit[0:2])
        plt.ylim(self.limit[2:4])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if rotate :
            plt.xlabel(ylabel)
            plt.ylabel(xlabel)
        plt.savefig(name + "_gaussian_fit.pdf",format="pdf")

        plt.figure()
        plt.hist2d(list_map,list_box,300)
        plt.plot(x_poly,y_poly,"r-")
        plt.xlim(self.limit[0:2])
        plt.ylim(self.limit[2:4])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if rotate :
            plt.xlabel(ylabel)
            plt.ylabel(xlabel)
        plt.savefig(name + "_histo2d.pdf",format="pdf")

        del list_box, list_map

        return(poly_fit,gaussian_fit,correlation_matrix)





    def plot_scatter_los_DM(self,pixel_dm,pixel_ra_dec,gaussian_smoothing=None):

        pixel_DM = np.fromfile(open(pixel_dm,"r"),dtype=np.float64)
        pixel_RA_DEC = np.fromfile(open(pixel_ra_dec,"r"),dtype=np.float64)
        pixel_RA_DEC = pixel_RA_DEC.reshape(len(pixel_RA_DEC)//5,5)
        delta_pixel_tomo = pixel_RA_DEC[:,4]
        sigma_pixel_tomo = pixel_RA_DEC[:,3]
        mask = sigma_pixel_tomo < 0.3
        delta_pixel_tomo = delta_pixel_tomo[mask]
        pixel_DM = pixel_DM[mask]
        plt.hist2d(delta_pixel_tomo,pixel_DM,100)
        print(np.corrcoef(delta_pixel_tomo,pixel_DM))




    def plot_DM_field(self,name,map_size,coordSlice,gaussian_smoothing=None,drhomin=-1,drhomax=1,direction ="x"):
        Tomo = tomography.TreatClamato(self.pwd,self.box_DM_name,self.shape_map,"")
        box_map = Tomo.readClamatoMapFile()
        if gaussian_smoothing is not None :
            box_map = self.gaussianSmoothing(box_map,gaussian_smoothing)
        Tomo.printSlice(direction,coordSlice,box_map,map_size,name + "_direction_{}_{}Pixel".format(direction,coordSlice),drhomin,drhomax)
        del Tomo
