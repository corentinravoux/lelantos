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
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize
from scipy import interpolate
from lslyatomo import utils
from lslyatomo import tomographic_objects



#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################



class BoxExtractor():


    def __init__(self,pwd,box_dir,box_shape,size_cell,box_bound,master_file,interpolation_method="LINEAR"):

        self.pwd = pwd
        self.size_cell = size_cell
        self.box_bound = box_bound
        self.box_shape = box_shape
        self.box_dir = box_dir
        self.master_file = master_file
        self.interpolation_method = interpolation_method
        self.log = utils.create_report_log(name=os.path.join(self.pwd,"Python_Report"))





### BOX ###

    def create_box_array(self,X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map):
        X_tomo_array = np.linspace(X_tomo_min,X_tomo_max,shape_map[0])
        Y_tomo_array = np.linspace(Y_tomo_min,Y_tomo_max,shape_map[1])
        Z_tomo_array = np.linspace(Z_tomo_min,Z_tomo_max,shape_map[2])
        return(X_tomo_array,Y_tomo_array,Z_tomo_array)




    def construct_DM_map(self,ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,get_prop=None):
        self.log.add("Creation of (RA,DEC,R) data matrix")
        coords_ra_dec =np.moveaxis(np.array(np.meshgrid(ra_array,dec_array,h * R_of_z(z_array),indexing='ij')),0,-1)
        self.log.add("Conversion to (X,Y,Z) coordinates in the Saclay box")
        coords_box_saclay = np.zeros(coords_ra_dec.shape)
        (coords_box_saclay[:,:,:,0],
         coords_box_saclay[:,:,:,1],
         coords_box_saclay[:,:,:,2]) = utils.saclay_mock_sky_to_cartesian(
                                       coords_ra_dec[:,:,:,0],
                                       coords_ra_dec[:,:,:,1],
                                       coords_ra_dec[:,:,:,2],
                                       ra0_box,dec0_box)
        del coords_ra_dec
        self.log.add("Searching for the nearest pixels of the Saclay box")
        coords_pixels_box_saclay = np.zeros(coords_box_saclay.shape)
        (coords_pixels_box_saclay[:,:,:,0],
         coords_pixels_box_saclay[:,:,:,1],
         coords_pixels_box_saclay[:,:,:,2]) = utils.saclay_mock_coord_dm_map(
                                              coords_box_saclay[:,:,:,0],
                                              coords_box_saclay[:,:,:,1],
                                              coords_box_saclay[:,:,:,2],
                                              Rmin,self.size_cell,
                                              self.box_shape,
                                              self.interpolation_method)
        del coords_box_saclay
        self.log.add("Loading of the Saclay map")
        DM_mocks_map = utils.saclay_mock_get_box(self.box_dir,self.box_shape)
        self.log.add("Creation of the DM map")
        DM_map = utils.interpolate_map(self.interpolation_method,DM_mocks_map,coords_pixels_box_saclay)
        del DM_mocks_map
        if(get_prop is not None):
            Props_map = []
            for i in range(len(get_prop)):
                self.log.add("Loading of the Saclay {} map".format(get_prop[i]))
                prop_mocks_map = utils.saclay_mock_get_box(self.box_dir,self.box_shape,name_box=get_prop[i])
                self.log.add("Creation of the {} map".format(get_prop[i]))
                prop_map = utils.interpolate_map(self.interpolation_method,prop_mocks_map,coords_pixels_box_saclay)
                Props_map.append(prop_map)
                del prop_mocks_map
            Props_map = np.array(Props_map)
        else:
            Props_map = None
        del coords_pixels_box_saclay
        self.log.add("Multiplying by growth factor at redshift of the LOS")
        return(DM_map,Props_map)




    def fill_DM_map(self,ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,DM_map,i_box,get_prop=None,Props_map=None):
        self.log.add("Creation of (RA,DEC,R) data matrix")
        coords_ra_dec =np.moveaxis(np.array(np.meshgrid(ra_array,dec_array,h * R_of_z(z_array),indexing='ij')),0,-1)
        mask1 = (DM_map[:,:,:] == None)
        mask2 = (coords_ra_dec[:,:,:,0]>=np.radians(self.box_bound[i_box][0]))
        mask2 &=(coords_ra_dec[:,:,:,0]<np.radians(self.box_bound[i_box][1]))
        mask2 &=(coords_ra_dec[:,:,:,1]>=np.radians(self.box_bound[i_box][2]))
        mask2 &=(coords_ra_dec[:,:,:,1]<np.radians(self.box_bound[i_box][3]))
        if((len(mask1[mask1==True]) == 0)|(len(mask2[mask2==True])==0)):
            self.log.add("No need of the Saclay box {}".format(i_box))
            return(DM_map,Props_map)
        mask = mask1&mask2
        del mask1,mask2
        self.log.add("Conversion to (X,Y,Z) coordinates in the Saclay box {}".format(i_box))
        coords_box_saclay = np.zeros(coords_ra_dec.shape)
        (coords_box_saclay[:,:,:,0],
         coords_box_saclay[:,:,:,1],
         coords_box_saclay[:,:,:,2]) = utils.saclay_mock_sky_to_cartesian(
                                       coords_ra_dec[:,:,:,0],
                                       coords_ra_dec[:,:,:,1],
                                       coords_ra_dec[:,:,:,2],
                                       ra0_box,dec0_box)
        del coords_ra_dec
        self.log.add("Searching for the nearest pixels of the Saclay box {}".format(i_box))
        coords_pixels_box_saclay = np.zeros(coords_box_saclay.shape)
        coords_pixels_box_saclay[:,:,:,0],coords_pixels_box_saclay[:,:,:,1],coords_pixels_box_saclay[:,:,:,2] = utils.saclay_mock_coord_dm_map(coords_box_saclay[:,:,:,0],coords_box_saclay[:,:,:,1],coords_box_saclay[:,:,:,2],Rmin,self.size_cell,self.box_shape[i_box],self.interpolation_method)
        del coords_box_saclay
        self.log.add("Loading of the Saclay map {}".format(i_box))
        DM_mocks_map = utils.saclay_mock_get_box(self.box_dir[i_box],self.box_shape[i_box])
        self.log.add("Creation of the DM map with box {}".format(i_box))
        DM_map[mask] = utils.interpolate_and_fill_map(self.interpolation_method,DM_mocks_map,coords_pixels_box_saclay[mask])
        del DM_mocks_map
        if(get_prop is not None):
            for i in range(len(get_prop)):
                self.log.add("Loading of the Saclay {} map {}".format(get_prop[i],i_box))
                prop_mocks_map = utils.saclay_mock_get_box(self.box_dir[i_box],self.box_shape[i_box],name_box=get_prop[i])
                self.log.add("Creation of the {} map with box {}".format(get_prop[i],i_box))
                Props_map[i][mask] = utils.interpolate_and_fill_map(self.interpolation_method,prop_mocks_map,coords_pixels_box_saclay[mask])
                del prop_mocks_map
        del coords_pixels_box_saclay,mask,get_prop
        return(DM_map,Props_map)


    def extract_box(self,ra_array,dec_array,z_array,shape_map_output,rsd_box):
        (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = utils.saclay_mock_box_cosmo_parameters(self.box_shape,self.size_cell)
        ra0_box, dec0_box = utils.saclay_mock_center_of_the_box(self.box_bound)
        if(rsd_box): get_prop = ["eta_zz"]
        else: get_prop = None
        DM_map,Props_map = self.construct_DM_map(ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,get_prop=get_prop)
        return(DM_map,Props_map,get_prop)


    def extract_box_multiple(self,ra_array,dec_array,z_array,shape_map_output,rsd_box):
        DM_map = np.full(tuple(shape_map_output),None)
        if(rsd_box):
            get_prop = ["eta_zz","vz"]
            Props_map = np.array([np.full(tuple(shape_map_output),None) for i in range(len(get_prop))])
        else: Props_map,get_prop = None,None
        for i_box in range(len(self.box_dir)):
            (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = utils.saclay_mock_box_cosmo_parameters(self.box_shape[i_box],self.size_cell)
            ra0_box, dec0_box = utils.saclay_mock_center_of_the_box(self.box_bound[i_box])
            DM_map,Props_map = self.fill_DM_map(ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,shape_map_output,Rmin,h,DM_map,i_box,get_prop=get_prop,Props_map=Props_map)
        if(len(DM_map[DM_map==None]!=0)):
            self.log.add("WARNING : Not enough boxes to fill the Dark Matter map")
        self.log.add("Multiplying by growth factor at redshift of the LOS")
        return(DM_map,Props_map,get_prop)



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
        sigma_l = np.std(gaussian_array)
        density_matter = np.exp((gaussian_array - (sigma_l**2/2)).astype(float)) - 1
        return(density_matter)





    def compute_rsd(self,DM_map,Props_map,z_array,Cosmo):
        for i in range(len(z_array)):
            DM_map[:,:,i] = DM_map[:,:,i]  - (((1+z_array[i])* Props_map[0][:,:,i])/Cosmo.hubble(z_array[i]))
        return(DM_map)



### LOS ###



    def construct_DM_LOS(self,ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,Rmin,h):
        self.log.add("Creation of (RA,DEC,R) data matrix")
        R_array = h * R_of_z(z_array)
        self.log.add("Conversion to (X,Y,Z) coordinates in the Saclay box")
        X,Y,Z = np.zeros(len(ra_array)),np.zeros(len(dec_array)),np.zeros(len(R_array))
        X,Y,Z = utils.saclay_mock_sky_to_cartesian(ra_array,dec_array,R_array,ra0_box,dec0_box)
        self.log.add("Searching for the nearest pixels of the Saclay box")
        i,j,k = np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int)
        i,j,k = utils.saclay_mock_coord_dm_map(X,Y,Z,Rmin,self.size_cell,self.box_shape,self.interpolation_method)
        del X,Y,Z
        self.log.add("Loading of the Saclay map")
        DM_mocks_map = utils.saclay_mock_get_box(self.box_dir,self.box_shape)
        self.log.add("Creation of the DM los")
        DM_LOS = np.zeros(len(R_array))
        DM_LOS[:] = DM_mocks_map[i[:],j[:],k[:]]
        del DM_mocks_map,i,j,k
        return(DM_LOS)

    def fill_DM_LOS(self,ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,Rmin,h,DM_LOS,i_box):
        self.log.add("Creation of (RA,DEC,R) data matrix")
        R_array = h * R_of_z(z_array)
        mask1 = (DM_LOS[:] == None)
        mask2 = (ra_array[:]>=np.radians(self.box_bound[i_box][0]))
        mask2 &=(ra_array[:]<np.radians(self.box_bound[i_box][1]))
        mask2 &=(dec_array[:]>=np.radians(self.box_bound[i_box][2]))
        mask2 &=(dec_array[:]<np.radians(self.box_bound[i_box][3]))
        if((len(mask1[mask1==True]) == 0)|(len(mask2[mask2==True])==0)):
            self.log.add("No need of the Saclay box {}".format(i_box))
            return(DM_LOS)
        mask = mask1&mask2
        del mask1,mask2
        self.log.add("Conversion to (X,Y,Z) coordinates in the Saclay box {}".format(i_box))
        X,Y,Z = np.zeros(len(ra_array)),np.zeros(len(dec_array)),np.zeros(len(R_array))
        X,Y,Z = utils.saclay_mock_sky_to_cartesian(ra_array,dec_array,R_array,ra0_box,dec0_box)
        self.log.add("Searching for the nearest pixels of the Saclay box {}".format(i_box))
        i,j,k = np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int),np.zeros(len(R_array)).astype(int)
        i,j,k = utils.saclay_mock_coord_dm_map(X,Y,Z,Rmin,self.size_cell,self.box_shape[i_box],self.interpolation_method)
        del X,Y,Z
        self.log.add("Loading of the Saclay map {}".format(i_box))
        DM_mocks_map = utils.saclay_mock_get_box(self.box_dir[i_box],self.box_shape[i_box])
        self.log.add("Creation of the DM los with box {}".format(i_box))
        DM_LOS[:][mask] = DM_mocks_map[i[:][mask],j[:][mask],k[:][mask]]
        del DM_mocks_map,i,j,k
        return(DM_LOS)




    def extract_los(self,ra_array,dec_array,z_array):
        (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = utils.saclay_mock_box_cosmo_parameters(self.box_shape,self.size_cell)
        ra0_box, dec0_box = utils.saclay_mock_center_of_the_box(self.box_bound)
        los = self.construct_DM_LOS(ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,Rmin,h)
        return(los)


    def extract_los_multiple(self,ra_array,dec_array,z_array):
        los = np.full(z_array.shape,None)
        for i_box in range(len(self.box_dir)):
            (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = utils.saclay_mock_box_cosmo_parameters(self.box_shape[i_box],self.size_cell)
            ra0_box, dec0_box = utils.saclay_mock_center_of_the_box(self.box_bound[i_box])
            los = self.fill_DM_LOS(ra_array,dec_array,z_array,R_of_z,ra0_box,dec0_box,Rmin,h,los,i_box)
        if(len(los[los==None]!=0)):
            self.log.add("Warning : Not enough boxes to fill the Dark Matter map")
        return(los)




### DELTA ###


    def extract_delta_multiple(self,ra,dec,z):
        delta_quasars =np.full(z.shape,None)
        for i_box in range(len(self.box_dir)):
            (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = utils.saclay_mock_box_cosmo_parameters(self.box_shape[i_box],self.size_cell)
            ra0_box, dec0_box = utils.saclay_mock_center_of_the_box(self.box_bound[i_box])
            delta_quasars = self.fill_DM_LOS(ra,dec,z,R_of_z,ra0_box,dec0_box,Rmin,h,delta_quasars,i_box)
        if(len(delta_quasars[delta_quasars==None]!=0)):
            self.log.add("Warning : Not enough boxes to fill the Dark Matter map")
        return(delta_quasars)


    def extract_delta(self,ra,dec,z):
        (R0,z0,R_of_z,z_of_R,Rmin,Rmax,h) = utils.saclay_mock_box_cosmo_parameters(self.box_shape,self.size_cell)
        ra0_box, dec0_box = utils.saclay_mock_center_of_the_box(self.box_bound)
        delta_quasars = self.construct_DM_LOS(ra,dec,z,R_of_z,ra0_box,dec0_box,Rmin,h)
        return(delta_quasars)




######## Main ######



    def create_box(self,map_property_file,name,rsd_box=False,growth_multiplication=True,matter_field=True,gaussian_smoothing=None,shape_map_output=None):
        property_file = tomographic_objects.MapPixelProperty(name=map_property_file)
        property_file.read()
        if(shape_map_output is None):
            shape_map_output = property_file.shape
        X_tomo_min,Y_tomo_min,Z_tomo_min = property_file.boundary_cartesian_coord[0]
        X_tomo_max,Y_tomo_max,Z_tomo_max = property_file.boundary_cartesian_coord[1]
        X_tomo_array,Y_tomo_array, Z_tomo_array = self.create_box_array(X_tomo_min,X_tomo_max,Y_tomo_min,Y_tomo_max,Z_tomo_min,Z_tomo_max,shape_map_output)
        suplementary_parameters = utils.return_suplementary_parameters(property_file.coordinate_transform,property=property_file)
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(property_file.Omega_m)
        ra_array,dec_array,z_array  = utils.convert_cartesian_to_sky(X_tomo_array,
                                                                     Y_tomo_array,
                                                                     Z_tomo_array,
                                                                     property_file.coordinate_transform,
                                                                     inv_rcomov=inv_rcomov,
                                                                     inv_distang=inv_distang,
                                                                     distang=distang,
                                                                     suplementary_parameters=suplementary_parameters)
        del X_tomo_array,Y_tomo_array, Z_tomo_array,suplementary_parameters,rcomov,distang,inv_rcomov,inv_distang
        if(type(self.box_dir) == list):
            (dm_map,prop_maps,prop) = self.extract_box_multiple(ra_array,dec_array,z_array,shape_map_output,rsd_box)
        else :
            (dm_map,prop_maps,prop) = self.extract_box(ra_array,dec_array,z_array,shape_map_output,rsd_box)
        del ra_array,dec_array
        if(growth_multiplication): dm_map = self.multiply_by_growth(dm_map,z_array)
        del z_array
        if(matter_field): dm_map = self.convert_to_matter_field(dm_map)
        if(gaussian_smoothing is not None):
            gaussian_smoothing_pix = gaussian_smoothing*utils.pixel_per_mpc(property_file.size,shape_map_output)
            dm_map = utils.gaussian_smoothing(dm_map,gaussian_smoothing_pix)
        dm_map_object = tomographic_objects.TomographicMap(map_array=np.array(dm_map),name=name)
        del dm_map
        dm_map_object.write()
        del dm_map_object
        if(rsd_box):
            for i in range(len(prop)):
                prop_map_object = tomographic_objects.TomographicMap(map_array=prop_maps[i],name=f"{name}_{prop[i]}")
                prop_map_object.write()
            del prop_map_object
        del prop_maps



    def create_LOS(self,property_file,pixel_name,name,growth_multiplication=True,matter_field=True):
        pixel = tomographic_objects.Pixel(name=pixel_name)
        pixel.read()
        property_file = tomographic_objects.MapPixelProperty(name=property_file)
        property_file.read()
        suplementary_parameters = utils.return_suplementary_parameters(property_file.coordinate_transform,property=property_file)
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(property_file.Omega_m)
        X_array,Y_array,Z_array = np.array(pixel.pixel_array[:,0]),np.array(pixel.pixel_array[:,1]),np.array(pixel.pixel_array[:,2])
        ra_array,dec_array,z_array  = utils.convert_cartesian_to_sky(X_array,
                                                                     Y_array,
                                                                     Z_array,
                                                                     property_file.coordinate_transform,
                                                                     inv_rcomov=inv_rcomov
                                                                     ,inv_distang=inv_distang,
                                                                     distang=distang,
                                                                     suplementary_parameters=suplementary_parameters)
        if(type(self.box_dir) == list):
            los = self.extract_los(ra_array,dec_array,z_array)
        else :
            los = self.extract_los_multiple(ra_array,dec_array,z_array)
        self.log.add("Multiplying by growth factor at redshift of the LOS")
        if(growth_multiplication):
            los = self.multiply_los_by_growth(los,z_array)
        if(matter_field):
            los = self.convert_to_matter_field(los)
        self.log.add("Mean delta extracted : {}".format(np.mean(los)))
        pixel_out = tomographic_objects.Pixel(name=name,pixel_array=los)
        pixel_out.read()


    def create_catalog(self,cat_name,type_catalog,name,growth_multiplication=True,matter_field=True):
        catalog = tomographic_objects.Catalog.init_catalog_from_fits(cat_name,type_catalog)
        ra,dec,z = catalog.coord[:,0],catalog.coord[:,1],catalog.coord[:,2]
        ra[ra>180] = ra[ra>180] -360
        if(type(self.box_dir) == list):
            delta_quasars = self.create_delta_catalog_several_boxes(ra,dec,z)
        else :
            delta_quasars = self.create_delta_catalog(ra,dec,z)
        self.log.add("Multiplying by growth factor at redshift of the LOS")
        if(growth_multiplication):
            delta_quasars = self.multiply_los_by_growth(delta_quasars,z)
        if(matter_field):
            delta_quasars = self.convert_to_matter_field(delta_quasars)
        self.log.add("Mean delta extracted : {}".format(np.mean(delta_quasars)))
        return(delta_quasars)





class BoxPlot(object):

    def __init__(self,pwd,map_name,box_name,property_file,limit=None):

        self.pwd = pwd
        self.map_name = map_name
        self.box_name = box_name
        self.property_file = property_file
        self.limit = limit

    @staticmethod
    def contourplot(x,y,ncont=10,colors=None,pltclf=True,binsx=100,binsy=100,log=False) :
        H, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
        H = np.rot90(H)
        H = np.flipud(H)
        xcenters=(xedges[:-1]+xedges[1:])/2.
        ycenters=(yedges[:-1]+yedges[1:])/2.
        if log is True : H[(H>0)]=np.log(H[(H>0)])
        if pltclf : plt.clf()
        plt.contour(xcenters,ycenters,H,ncont,colors=colors,linewidths=2)

    @staticmethod
    def densityplot(x,y,binsx=200,binsy=200,scaleperdeg2=False,maxdensity=0,mindensity=0,pltclf=True,logscale=False) :
        H, xedges, yedges = np.histogram2d(x,y,bins=(binsx,binsy))
        if scaleperdeg2==True :
            delta_ra=abs(xedges[1]-xedges[0])
            delta_dec=abs(yedges[1]-yedges[0])
            cosdelta=np.cos(yedges[:-1]*np.pi/180.)
            for k in range(binsy) : H[:,k] = H[:,k]/(delta_ra*delta_dec*cosdelta[k])
        H = np.rot90(H)
        H = np.flipud(H)
        if maxdensity>0 :
            wcut=np.where( (H>maxdensity) )
            H[wcut]=maxdensity
        if mindensity>0 :
            wcut=np.where( (H<mindensity) )
            H[wcut]=0
        Hmasked = np.ma.masked_where(H==0,H)
        if pltclf : plt.clf()
        if logscale is False :
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            Hmaskednorm = Hmasked/np.max(Hmasked)
            plt.imshow(Hmaskednorm, interpolation='bilinear',cmap='jet',origin="lower",extent = extent,norm=Normalize(vmin=0.0,vmax=1.0))
        else : plt.pcolormesh(xedges,yedges,Hmasked,norm=LogNorm(), cmap=plt.cm.jet)
        cbar = plt.colorbar()
        cbar.set_label("Normalized density")
        if scaleperdeg2==True : plt.title(r'density / deg$^2$')


    @staticmethod
    def gaussian_plot(X,Y,binsx=200,binsy=200,ncont=4):
        data, xedges,yedges = np.histogram2d(X,Y,bins=(binsx,binsy))
        xcenters=(xedges[:-1]+xedges[1:])/2.
        ycenters=(yedges[:-1]+yedges[1:])/2.
        Fitter = utils.gaussian_fitter_2d(inpdata=data)
        p,success = Fitter.FitGauss2D()
        x,y=np.indices((binsx,binsy),dtype=np.float)
        gauss = Fitter.Gaussian2D(*p)
        data_fitted = gauss(y,x)
        if(type(ncont) != int):
            ncont =np.array(ncont)*p[0]
        plt.contour(xcenters,ycenters,data_fitted.reshape(binsx, binsy),ncont,linewidths=2, colors='w')
        return(p)

    @staticmethod
    def plot_gaussian_line(gaussian_fit,binsx,limit):
        xcenter = (gaussian_fit[1] - binsx)*((limit[1] - limit[0])/binsx)
        ycenter = (gaussian_fit[2] - binsx)*((limit[3] - limit[2])/binsx)
        xcenter = 0
        ycenter = 0
        angle = gaussian_fit[5]
        pente = 1/np.tan(np.radians(angle))
        x_arrray = np.linspace(limit[0],limit[1],100)
        y_array = ycenter + pente * (x_arrray - xcenter)
        plt.plot(x_arrray,y_array)


    @staticmethod
    def plot_scatter_los_DM(self,pixel_dm_name,pixel_name,gaussian_smoothing=None):
        pixel_dm = tomographic_objects.Pixel(name=pixel_dm_name)
        pixel_dm.read()
        pixel = tomographic_objects.Pixel(name=pixel_name)
        pixel.read()
        delta_pixel_tomo = pixel.pixel_array[:,4]
        sigma_pixel_tomo = pixel.pixel_array[:,3]
        mask = sigma_pixel_tomo < 0.3
        delta_pixel_tomo = delta_pixel_tomo[mask]
        pixel_DM = pixel_dm.pixel_array[mask]
        plt.hist2d(delta_pixel_tomo,pixel_DM,100)
        print(np.corrcoef(delta_pixel_tomo,pixel_DM))



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


    def get_list_boxes(self,gaussian_smoothing=None,cut_redshift_coef=None,cut_pixel=None):
        tomography_map = tomographic_objects.TomographicMap.init_from_property_files(self.property_file,name=self.map_name)
        tomography_map.read()
        box = tomographic_objects.TomographicMap.init_from_property_files(self.property_file,name=self.box_name)
        box.read()
        max_x,max_y,max_z =tomography_map.shape[0],tomography_map.shape[1],tomography_map.shape[2]
        if(cut_redshift_coef is not None):
            max_z = tomography_map.shape[2]//cut_redshift_coef
        if(cut_pixel is not None):
            max_z = cut_pixel
        min_x,min_y,min_z = 0,0,0
        length_list = (max_x-min_x)*(max_y-min_y)*(max_z-min_z)
        box_DM = box.map_array[min_x:max_x,min_y:max_y,min_z:max_z]
        if(gaussian_smoothing is not None):
            gaussian_smoothing_pix = gaussian_smoothing*utils.pixel_per_mpc(box.size,box.shape)
            box_DM = utils.gaussian_smoothing(box_DM,gaussian_smoothing_pix)
        list_box = box_DM.reshape(length_list)
        del box_DM, box
        map_3D = tomography_map.map_array[min_x:max_x,min_y:max_y,min_z:max_z]
        list_map = map_3D.reshape(length_list)
        del map_3D, tomography_map
        return(list_map,list_box)
