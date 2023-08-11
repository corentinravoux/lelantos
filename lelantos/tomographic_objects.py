#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Corentin Ravoux

Description :
"""



#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################



import numpy as np
import fitsio
import pickle
import os
from lelantos import utils
from scipy.ndimage.filters import gaussian_filter
import multiprocessing as mp
from scipy.interpolate import interp1d
from picca import data


### Copy to utils ####


def verify_file(file):
    if(not(os.path.isfile(file))):
        raise ValueError("No property file {} was found".format(file))



#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################









#############################################################################
#############################################################################
############################### MAP OBJECTS #####################################
#############################################################################
#############################################################################









class TomographicMap(object):

    def __init__(self,map_array=None,name=None,shape=None,size=None,boundary_cartesian_coord=None,boundary_sky_coord=None,coordinate_transform=None,property_file=None,Omega_m=None):
        self.name = name
        self.shape = shape
        self.size = size
        self.map_array = map_array
        self.boundary_cartesian_coord = boundary_cartesian_coord
        self.boundary_sky_coord = boundary_sky_coord
        self.coordinate_transform = coordinate_transform
        self.property_file = property_file
        self.Omega_m = Omega_m

    @classmethod
    def init_classic(cls,map_array=None,name=None,shape=None,size=None,boundary_cartesian_coord=None,boundary_sky_coord=None,coordinate_transform=None,property_file=None,Omega_m=None):
        if(property_file is not None):
            return(cls.init_from_property_files(property_file,map_array=map_array,name=name))
        else:
            return(cls(map_array=map_array,name=name,shape=shape,size=size,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,property_file=property_file,Omega_m=Omega_m))

    @classmethod
    def init_from_property_files(cls,property_file,map_array=None,name=None):
        Property = MapPixelProperty(name=property_file)
        verify_file(property_file)
        size,shape,boundary_cartesian_coord,boundary_sky_coord,coordinate_transform,Omega_m = None,None,None,None,None,None
        Property.read()
        if(Property.size is not None): size = Property.size
        if(Property.shape is not None): shape = Property.shape
        if(Property.boundary_cartesian_coord is not None): boundary_cartesian_coord = Property.boundary_cartesian_coord
        if(Property.boundary_sky_coord is not None): boundary_sky_coord = Property.boundary_sky_coord
        if(Property.coordinate_transform is not None): coordinate_transform = Property.coordinate_transform
        if(Property.Omega_m is not None): Omega_m = Property.Omega_m
        return(cls(map_array=map_array,name=name,shape=shape,size=size,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,property_file=property_file,Omega_m=Omega_m))


    @classmethod
    def init_by_merging(cls,submap_directory,launching_file_name,name_map,property_file):
        """ Need to generalized (maybe in the task manager class for reading and writting launching files) + naming filename must not be define there"""
        a = pickle.load(open(os.path.join(launching_file_name),"rb"))
        listname,dachshund_params,number_chunks,overlaping = a[0],a[1],a[2],a[3]
        map_chunks,sizemap,shapemap = {},{},{}
        for i in range(len(listname)):
            name = os.path.join(submap_directory,dachshund_params[i]["namemap"])
            sizemap[listname[i]] = (dachshund_params[i]["lx"],dachshund_params[i]["ly"],dachshund_params[i]["lz"])
            shapemap[listname[i]] = (dachshund_params[i]["nx"],dachshund_params[i]["ny"],dachshund_params[i]["nz"])
            submap = TomographicMap(name=name,shape=shapemap[listname[i]],size=sizemap[listname[i]])
            submap.read()
            map_chunks[listname[i]] = submap.map_array
        if overlaping !=0 :
            for i in range(number_chunks[0]):
                for j in range(number_chunks[1]):
                    filename = f'{i:03d}' + f'{j:03d}'
                    sub_map = map_chunks[filename]
                    pixel_to_remove = np.around(utils.pixel_per_mpc(sizemap[filename],shapemap[listname[i]]) * overlaping,decimals=0).astype(int)
                    if ((i==0)&(i == number_chunks[0] - 1)):
                        sub_map = sub_map
                    elif(i==0):
                        sub_map = sub_map[0:len(sub_map)-pixel_to_remove[0],:,:]
                    elif(i == number_chunks[0] - 1):
                        sub_map = sub_map[pixel_to_remove[0]:len(sub_map),:,:]
                    else :
                        sub_map = sub_map[pixel_to_remove[0]:len(sub_map)-pixel_to_remove[0],:,:]
                    if ((j==0)&(j == number_chunks[1] - 1)):
                        sub_map = sub_map
                    elif(j==0) :
                        sub_map = sub_map[:,0:len(sub_map[0])-pixel_to_remove[1],:]
                    elif(j == number_chunks[1] - 1):
                        sub_map = sub_map[:,pixel_to_remove[1]:len(sub_map[0]),:]
                    else :
                        sub_map = sub_map[:,pixel_to_remove[1]:len(sub_map[0])-pixel_to_remove[1],:]
                    map_chunks[filename]=sub_map
        concatenateList = []
        for i in range(number_chunks[0]):
            concatenate = np.concatenate([np.asarray(map_chunks[f'{i:03d}' + f'{j:03d}']) for j in range(number_chunks[1])],axis = 1)
            concatenateList.append(concatenate)
        merged_map_array = np.concatenate([concatenateList[i] for i in range(number_chunks[0])],axis = 0)
        merged_map = cls.init_from_property_files(property_file,map_array=merged_map_array,name=name_map)
        return(merged_map)



    @property
    def mpc_per_pixel(self):
        return(utils.mpc_per_pixel(self.size,self.shape))

    @property
    def pixel_per_mpc(self):
        return(utils.pixel_per_mpc(self.size,self.shape))


    def read(self):
        if(self.name is None):
            raise ValueError("No name was provided for reading map")
        if(self.shape is None):
            raise ValueError("No shape was provided for reading map")
        else:
            with open(self.name,'r') as f:
                self.map_array = np.fromfile(f,dtype=np.float64).reshape(self.shape)



    def write(self):
        if(self.map_array.all() == None):
            raise ValueError("No map array is stored in this map class for writting")
        listmap=np.ravel(self.map_array)
        listmap.tofile(self.name)

    def write_in_vtk(self):
        from pyevtk.hl import gridToVTK
        nx,ny,nz = self.shape
        lx,ly,lz = self.size
        X = np.linspace(0,lx,nx, dtype='float64')
        Y = np.linspace(0,ly,ny, dtype='float64')
        Z = np.linspace(0,lz,nz, dtype='float64')
        x = np.zeros((nx, ny, nz))
        y = np.zeros((nx, ny, nz))
        z = np.zeros((nx, ny, nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x[i,j,k] = X[i]
                    y[i,j,k] = Y[j]
                    z[i,j,k] = Z[k]
        gridToVTK(self.name, x, y, z, pointData = {"DeltaF" : self.map_array})


    def write_property_file(self,property_file_name):
        property_file = MapPixelProperty(name=property_file_name,size=self.size,
                                         shape=self.shape,
                                         boundary_cartesian_coord=self.boundary_cartesian_coord,
                                         boundary_sky_coord=self.boundary_sky_coord,
                                         coordinate_transform=self.coordinate_transform,
                                         Omega_m=self.Omega_m)
        property_file.write()

    def rebin_map(self,new_shape, operation='mean'):
        self.map_array = utils.bin_ndarray(self.map_array, new_shape, operation=operation)
        self.shape = new_shape


    def smooth(self,sigma_mpc):
        sigma = sigma_mpc / self.pixel_per_mpc
        self.map_array = gaussian_filter(self.map_array,sigma)



    def mask_map(self,distance_name,distance):
        distance_map = DistanceMap.init_from_tomographic_map(self,name=distance_name)
        distance_map.read()
        mask = distance_map.get_mask_distance(distance)
        self.map_array[mask] = 0


    def mask_map_from_name(self,distance_map_name,distance):
        distance_map = DistanceMap.init_from_tomographic_map(self,name=distance_map_name)
        mask = distance_map.get_mask_distance(distance)
        self.map_array = np.ma.masked_where(mask,self.map_array)

    def mask_map_from_mask(self,mask):
        self.map_array = np.ma.masked_where(mask,self.map_array)



    def compute_pk3d(self,kmin,kmax,n_k,distance_map=None,criteria_distance_mask=None,log=False):
        """ Never tested properly, might be false in term of mpc per pixel conversion"""
        if((distance_map is not None)&(criteria_distance_mask is not None)):
            self.mask_map_from_name(distance_map,criteria_distance_mask)

        map_3D = self.map_array
        map_fft_3D = np.fft.fftn(map_3D)
        del map_3D
        kx = np.fft.fftfreq(map_fft_3D.shape[0],self.mpc_per_pixel[0])
        ky = np.fft.fftfreq(map_fft_3D.shape[1],self.mpc_per_pixel[1])
        kz = np.fft.fftfreq(map_fft_3D.shape[2],self.mpc_per_pixel[2])
        kx_space , ky_space, kz_space = np.meshgrid(kx,ky,kz)
        del kx,ky,kz
        normalization_factor = (self.mpc_per_pixel[0]*self.mpc_per_pixel[1]*self.mpc_per_pixel[2])/(self.shapeMap[0]*self.shapeMap[1]*self.shapeMap[2])
        power = normalization_factor * np.absolute(map_fft_3D)**2
        norm_k = np.array(map_fft_3D.shape)
        del map_fft_3D
        norm_k = np.sqrt(kx_space[:,:,:]**2 + ky_space[:,:,:]**2 + kz_space[:,:,:]**2)
        del kx_space,ky_space,kz_space
        if(log):
            k_space = np.logspace(kmin,kmax,n_k)
        else :
            k_space = np.linspace(np.min(norm_k),np.max(norm_k),n_k)
        delta_k = (np.max(norm_k)-np.min(norm_k)) / n_k
        Pk_3D = np.zeros(k_space.shape)
        for i in range(len(Pk_3D)) :
            mask = (norm_k < k_space[i] + delta_k)&(norm_k >= k_space[i])
            Pk_3D[i] = np.nanmean(power[mask])
        del power,norm_k, delta_k
        mask2 = Pk_3D != -1
        pk_3D_final = Pk_3D[mask2]
        k_space_final = k_space[mask2]
        del Pk_3D,k_space
        return(k_space_final,pk_3D_final)


    def displace_map(self,property_file_name_to_center,name_out_map,dist_map=None,name_dist_map_out=None,pixel=None,pixel_out=None,qso=None,qso_out=None):
        prop = MapPixelProperty(name= property_file_name_to_center)
        prop.read()
        if(self.shape != prop.shape):
            raise KeyError("Not implemented already for different map shape")
        map_3d = self.readClamatoMapFile()
        new_map = np.zeros(prop.shape)
        if(dist_map is not None):
            dist_map_3d = DistanceMap.init_from_tomographic_map(self,name=dist_map)
            dist_map_3d.read()
            new_dist_map = np.full(prop.shape,np.inf)
        indice = np.transpose(np.indices(prop.shape),axes=(1,2,3,0))
        coord_mpc = indice * prop.mpc_per_pixel
        del indice
        coord_centered = np.zeros(coord_mpc.shape)
        minx,miny,minz = self.boundary_cartesian_coord[0]
        minx2,miny2,minz2 = prop.boundary_cartesian_coord[0]
        coord_centered[:,:,:,0],coord_centered[:,:,:,1],coord_centered[:,:,:,2] = coord_mpc[:,:,:,0] - (minx - minx2),coord_mpc[:,:,:,1] -( miny - miny2),coord_mpc[:,:,:,2] - (minz - minz2)
        del coord_mpc
        indice_centered = np.round(coord_centered / self.mpc_per_pixel,0).astype(int)
        mask =(indice_centered[:,:,:,0] >= 0)&(indice_centered[:,:,:,0] < prop.shape[0])
        mask&=(indice_centered[:,:,:,1] >= 0)&(indice_centered[:,:,:,1] < prop.shape[1])
        mask&=(indice_centered[:,:,:,2] >= 0)&(indice_centered[:,:,:,2] < prop.shape[2])
        new_map[mask] = map_3d[indice_centered[:,:,:,0][mask],indice_centered[:,:,:,1][mask],indice_centered[:,:,:,2][mask]]
        self.map_array = new_map
        self.write()
        if(dist_map is not None):
            new_dist_map[mask] = dist_map_3d[indice_centered[:,:,:,0][mask],indice_centered[:,:,:,1][mask],indice_centered[:,:,:,2][mask]]
            dist_map_3d.map_array = new_dist_map
            dist_map_3d.write()
        del indice_centered,mask
        if(pixel is not None):
            pixel_class = Pixel(name=pixel)
            pixel_class.read()
            pixel_class.pixel_array
            pixel_class.pixel_array[:,0],pixel_class.pixel_array[:,1],pixel_class.pixel_array[:,2] = pixel_class.pixel_array[:,0] + (minx - minx2) ,pixel_class.pixel_array[:,1]  + (miny - miny2),pixel_class.pixel_array[:,2] + (minz - minz2)
            mask = (pixel_class.pixel_array[:,0]>=0)&(pixel_class.pixel_array[:,0]<prop.size[0])
            mask &= (pixel_class.pixel_array[:,1]>=0)&(pixel_class.pixel_array[:,1]<prop.size[1])
            mask &= (pixel_class.pixel_array[:,2]>=0)&(pixel_class.pixel_array[:,2]<prop.size[2])
            pixel_class.name = pixel_out
            pixel_class.write()
        if(qso is not None):
            qso_class = QSOCatalog.init_from_fits(qso)
            qso_class.coord[:,0],qso_class.coord[:,1],qso_class.coord[:,2] = qso_class.coord[:,0] + (minx - minx2) ,qso_class.coord[:,1]  + (miny - miny2),qso_class.coord[:,2] + (minz - minz2)
            mask = (qso_class.coord[:,0]>=0)&(qso_class.coord[:,0]<prop.size[0])
            mask &= (qso_class.coord[:,1]>=0)&(qso_class.coord[:,1]<prop.size[1])
            mask &= (qso_class.coord[:,2]>=0)&(qso_class.coord[:,2]<prop.size[2])
            qso_class.name = qso_out
            qso_class.write()






class DistanceMap(TomographicMap):

    def __init__(self,
                 map_array=None,
                 name=None,
                 shape=None,
                 size=None,
                 boundary_cartesian_coord=None,
                 boundary_sky_coord=None,
                 coordinate_transform=None,
                 property_file=None,
                 Omega_m=None):
        super(DistanceMap,self).__init__(map_array=map_array,name=name,shape=shape,size=size,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,property_file=property_file,Omega_m=Omega_m)


    @classmethod
    def init_from_tomographic_map(cls,tomographic_map,name=None,map_array=None):
        return(cls(map_array=map_array,name=name,shape=tomographic_map.shape,size=tomographic_map.size,boundary_cartesian_coord=tomographic_map.boundary_cartesian_coord,boundary_sky_coord=tomographic_map.boundary_sky_coord,coordinate_transform=tomographic_map.coordinate_transform,property_file=tomographic_map.property_file))


    @classmethod
    def init_by_computing(cls,pixel,tomographic_map,name_map,nb_process=1,radius_local=50):
        cls.log = utils.create_report_log(name=os.path.join("Python_Report"))
        if(tomographic_map.coordinate_transform.lower() == "middle"):
            distance_map = cls.create_mask_aligned_los(pixel,tomographic_map,nb_process=nb_process,radius_local=radius_local)
        return(cls(map_array=distance_map,name=name_map,shape=tomographic_map.shape,size=tomographic_map.size,boundary_cartesian_coord=tomographic_map.boundary_cartesian_coord,boundary_sky_coord=tomographic_map.boundary_sky_coord,coordinate_transform=tomographic_map.coordinate_transform,property_file=tomographic_map.property_file))


    @classmethod
    def create_mask_aligned_los(cls,pixel,tomographic_map,nb_process=1,radius_local=50):
        if((tomographic_map.size is None)|(tomographic_map.shape is None)): raise KeyError("Please provide map size and shape")
        if((pixel.pixel_array.all() is None)): raise KeyError("Please provide pixel array")
        if(nb_process == 1):
            distance_map = cls.create_distance_map_serial(pixel,tomographic_map,radius_local=radius_local)
        else:
            distance_map = cls.create_distance_map_parallel(pixel,tomographic_map,nb_process)
        return(distance_map)



    @classmethod
    def create_distance_map_serial(cls,pixel,tomographic_map,radius_local=50):
        x,y,z = pixel.repack_by_los()
        pixels = [[[x[i],y[i]],z[i]] for i in range(len(x))]
        distance_array = np.full(tomographic_map.shape,np.inf)
        indice = np.transpose(np.indices(tomographic_map.shape),axes=(1,2,3,0))
        indice_mpc = indice * tomographic_map.mpc_per_pixel
        del indice
        for i in range(len(pixels)):
            minimal_coordinates_local = [int(round(((pixels[i][0][0]-radius_local)/tomographic_map.mpc_per_pixel)[0],0)),
                                         int(round(((pixels[i][0][1]-radius_local)/tomographic_map.mpc_per_pixel)[1],0))]
            maximal_coordinates_local = [int(round(((pixels[i][0][0]+radius_local)/tomographic_map.mpc_per_pixel)[0],0)),
                                         int(round(((pixels[i][0][1]+radius_local)/tomographic_map.mpc_per_pixel)[1],0))]
            LOS_range_min_local = distance_array[max(minimal_coordinates_local[0],0):min(tomographic_map.shape[0],maximal_coordinates_local[0]),max(minimal_coordinates_local[1],0):min(tomographic_map.shape[1],maximal_coordinates_local[1]),0:tomographic_map.shape[2]]
            indice_local = indice_mpc[max(minimal_coordinates_local[0],0):min(tomographic_map.shape[0],maximal_coordinates_local[0]),max(minimal_coordinates_local[1],0):min(tomographic_map.shape[1],maximal_coordinates_local[1]),0:tomographic_map.shape[2]] # - np.array([pixels[i][0][0],pixels[i][0][1],0])

            LOS_range_local = np.full(LOS_range_min_local.shape,np.inf)
            mask1 = indice_local[:,:,:,2] < np.min(pixels[i][1])
            LOS_range_local[mask1] = np.sqrt((indice_local[:,:,:,0][mask1] - pixels[i][0][0])**2 + (indice_local[:,:,:,1][mask1] - pixels[i][0][1])**2 + (indice_local[:,:,:,2][mask1] - np.min(pixels[i][1]))**2)
            mask2 = (~mask1)&(indice_local[:,:,:,2] < np.max(pixels[i][1]))
            LOS_range_local[mask2] = np.sqrt((indice_local[:,:,:,0][mask2] - pixels[i][0][0])**2 + (indice_local[:,:,:,1][mask2] - pixels[i][0][1])**2)
            mask3 = ((~mask1)&(~mask2))
            LOS_range_local[mask3] = np.sqrt((indice_local[:,:,:,0][mask3] - pixels[i][0][0])**2 + (indice_local[:,:,:,1][mask3] - pixels[i][0][1])**2 + (indice_local[:,:,:,2][mask3] - np.max(pixels[i][1]))**2)

            LOS_range_min_local[:,:,:] = np.minimum(LOS_range_local,LOS_range_min_local)
            cls.log.add("LOS number {} over {} finished".format(i,len(pixels)))
        mask = distance_array == np.inf
        number_inf = len(mask[mask==True])
        if number_inf !=0 :
            cls.log.add("WARNING : {}% of the new map is still infinity.".format(np.round(number_inf/(tomographic_map.shape[0]*tomographic_map.shape[1]*tomographic_map.shape[2]),2)))
            cls.log.add("Pixels with inf put to the radius_local value. To avoid this warning, increase radius_local")
            distance_array[mask] = radius_local
        return(distance_array)


    @classmethod
    def create_distance_map_parallel(cls,pixel,tomographic_map,nb_process):
        cls.log.add("Creation of global variables")
        x,y,z = pixel.repack_by_los()
        pixels = [[[x[i],y[i]],z[i]] for i in range(len(x))]
        indice = np.transpose(np.indices(tomographic_map.shape),axes=(1,2,3,0))
        global dist_map,shared_arr
        dist_map = indice * tomographic_map.mpc_per_pixel
        del indice
        cls.log.add("Number LOS to treat : {}".format(len(pixels)))
        cls.log.add("Launching of Pool with shared array initialization")
        shared_arr = utils.init_shared_array(tomographic_map.shape)
        pool = mp.Pool(nb_process)
        pool.map(cls.worker_create_distance_map, pixels,tomographic_map.shape)
        cls.log.add("End of Pool")
        cls.log.add("Getting the map from shared array")
        distance_array = utils.mp_array_to_numpyarray(shared_arr).reshape(tomographic_map.shape)
        cls.log.add("Writing of the map")
        return(distance_array)

    @classmethod
    def worker_create_distance_map(cls,pixels,shape):
        cls.log.add("Creation of the dist map for one LOS")
        LOS_range = np.full(shape,np.inf)
        mask1 = (dist_map[:,:,:,2] < np.min(pixels[1]))|(dist_map[:,:,:,2] > np.min(pixels[1]))
        LOS_range[mask1] = np.sqrt((dist_map[:,:,:,0][mask1] - pixels[0][0])**2 + (dist_map[:,:,:,1][mask1] - pixels[0][1])**2 + (dist_map[:,:,:,2][mask1] - np.min(pixels[1]))**2)
        LOS_range[~mask1] = np.sqrt((dist_map[:,:,:,0][~mask1] - pixels[0][0])**2 + (dist_map[:,:,:,1][~mask1] - pixels[0][1])**2)
        cls.log.add("Dist map for one LOS created")
        with shared_arr.get_lock():
            cls.log.add("Lock obtained for a worker")
            LOS_range_min = utils.mp_array_to_numpyarray(shared_arr).reshape(shape)
            LOS_range_min[:,:,:] = np.minimum(LOS_range,LOS_range_min)
            del LOS_range_min
        cls.log.add("Lock released for a worker")
        cls.log.add("LOS treated")
        del LOS_range




    def get_mask_distance(self,distance):
        return(self.map_array > distance)



    def give_redshift_cut(self,distance_criteria,percent_criteria):
        if(self.coordinate_transform !="middle"):
            raise KeyError("Only implemented for middle coordinate transformation")
        mask_distance = self.get_mask_distance(distance_criteria)
        i,percent = 0,0.
        while((percent<percent_criteria)&(i< mask_distance.shape[-1])):
            nb_True = len(mask_distance[:,:,i][mask_distance[:,:,i] == True])
            percent = nb_True/(mask_distance.shape[0]*mask_distance.shape[1])
            i=i+1
        dist_mpc = i * self.mpc_per_pixel[2]
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(self.Omega_m)
        redshift = utils.convert_z_cartesian_to_sky_middle(np.array([dist_mpc]),inv_rcomov)[0]
        return(redshift)








class StackMap(TomographicMap):

    def __init__(self,
                 map_array=None,
                 name=None,
                 shape=None,
                 size=None,
                 boundary_cartesian_coord=None,
                 boundary_sky_coord=None,
                 coordinate_transform=None,
                 property_file=None,
                 Omega_m=None,
                 tomographic_map=None,
                 catalog=None,
                 ellipticity=None,
                 mean_los_distance=None):
        super(StackMap,self).__init__(map_array=map_array,
                                      name=name,
                                      shape=shape,
                                      size=size,
                                      boundary_cartesian_coord=boundary_cartesian_coord,
                                      boundary_sky_coord=boundary_sky_coord,
                                      coordinate_transform=coordinate_transform,
                                      property_file=property_file,
                                      Omega_m=Omega_m)

        self.tomographic_map = tomographic_map
        self.catalog = catalog
        self.ellipticity = ellipticity
        self.mean_los_distance = mean_los_distance

    @classmethod
    def init_by_merging(cls,stack_name,property_stack_name,name,property_name):
        mean_stack = []
        for i in range(len(stack_name)):
            stack = cls.init_stack_by_property_file(property_stack_name[i],name=stack_name)
            stack.read()
            mean_stack.append(stack.map_array)
        stack_class = cls(map_array=np.concatenate(mean_stack,axis=0),name=name,
                          shape=stack.shape,size=stack.size,
                          boundary_cartesian_coord=stack.boundary_cartesian_coord,
                          boundary_sky_coord=stack.boundary_sky_coord,
                          coordinate_transform=stack.coordinate_transform,
                          property_file=property_name)
        if(property_name is not None):
            stack_class.write_property_file(property_name)
        return(stack_class)



    @classmethod
    def init_stack_by_property_file(cls,
                                    property_file_name,
                                    name=None,
                                    tomographic_map_name=None,
                                    property_file_map=None,
                                    catalog_name=None,
                                    type_catalog=None):
        tomographic_map, catalog = None, None
        if(tomographic_map_name is not None)&(property_file_map is not None):
            tomographic_map = TomographicMap.init_from_property_files(property_file_map,name=tomographic_map_name)
        if(catalog_name is not None)&(type_catalog is not None):
            catalog = Catalog.init_catalog_from_fits(catalog_name, type_catalog)
        super(StackMap,cls).__init__()
        stack = cls.init_from_property_files(property_file_name,name=name)
        stack.tomo_map = tomographic_map
        stack.catalog = catalog
        return(stack)



    @classmethod
    def init_by_tomographic_map(cls,
                                tomographic_map,
                                catalog,
                                size_stack,
                                shape_stack,
                                property_file,
                                name=None,
                                interpolation_method="NEAREST",
                                normalized=False,
                                coordinate_convert=None):
        (coord,radius,min_radius) = cls.initiate_stack(coordinate_convert,normalized,tomographic_map,catalog)
        coord_stack = cls.create_init_stack_coordinates(shape_stack,size_stack)
        voids_list,local_maps = [],[]
        for i in range(len(coord)):
            if(normalized):
                coord_stack_local = coord_stack * (radius[i] / min_radius)
            else:
                coord_stack_local = coord_stack
            if(coordinate_convert is not None):
                property = MapPixelProperty(name=tomographic_map.property_file)
                property.read()
                coord_stack_local = coord_stack_local + cls.inv_convert_center(coord[i],tomographic_map.coordinate_transform,coordinate_convert,property)
                coord_stack_local = cls.convert_coord(coord_stack_local,tomographic_map.coordinate_transform,coordinate_convert,property)
            else:
                coord_stack_local = coord_stack_local + coord[i]
            index_stack_local = coord_stack_local/tomographic_map.mpc_per_pixel
            boolean = (np.max(index_stack_local[:,:,:,0])>=tomographic_map.shape[0])
            boolean |=(np.min(index_stack_local[:,:,:,0])<0)
            boolean |=(np.max(index_stack_local[:,:,:,1])>=tomographic_map.shape[1])
            boolean |=(np.min(index_stack_local[:,:,:,1])<0)
            boolean |=(np.max(index_stack_local[:,:,:,2])>=tomographic_map.shape[2])
            boolean |=(np.min(index_stack_local[:,:,:,2])<0)
            if(not(boolean)):
                voids_list.append(coord[i])
                local_maps.append(utils.interpolate_map(interpolation_method,
                                                        tomographic_map.map_array,
                                                        index_stack_local))
        del coord_stack,coord_stack_local,index_stack_local
        local_maps = np.array(local_maps)
        stack = np.mean(local_maps,axis=0)
        del local_maps
        boundary_cartesian_coord = None
        boundary_sky_coord = None
        if(normalized):
            size_stack = size_stack / float(min_radius)
        shape = (2*shape_stack[0]+1,2*shape_stack[1]+1,2*shape_stack[2]+1)
        stack_class =cls(tomographic_map=tomographic_map,
                   catalog=catalog,
                   map_array=stack,
                   name=name,
                   shape=shape,size=(2*size_stack,2*size_stack,2*size_stack),
                   boundary_cartesian_coord=boundary_cartesian_coord,
                   boundary_sky_coord=boundary_sky_coord,
                   coordinate_transform=tomographic_map.coordinate_transform,
                   property_file=property_file)
        stack_class.write_property_file(property_file)
        return(stack_class)



    @classmethod
    def initiate_stack(cls,coordinate_convert,normalized,tomographic_map,catalog):
        if(coordinate_convert is not None):
            if(coordinate_convert.lower() == tomographic_map.coordinate_transform):
                    raise ValueError("Please choose a coordinate transformation different than the map one")
        if(normalized):
            if(catalog.object_type.lower() != "void"):
                raise KeyError("It is not possible to normalize this catalog")
            else:
                coord = catalog.coord[catalog.radius >= normalized]
                radius = catalog.radius[catalog.radius >= normalized]
                min_radius = np.min(radius)
        else:
            coord = catalog.coord
            radius,min_radius = None,None
        return(coord,radius,min_radius)


    @classmethod
    def create_init_stack_coordinates(cls,shape_stack,size_stack):
        x_array = np.linspace(-size_stack,size_stack,2*shape_stack[0]+1)
        y_array = np.linspace(-size_stack,size_stack,2*shape_stack[1]+1)
        z_array = np.linspace(-size_stack,size_stack,2*shape_stack[2]+1)
        return(np.moveaxis(np.array(np.meshgrid(x_array,y_array,z_array,indexing='ij')),0,-1))


    @classmethod
    def convert_coord(cls,coord_stack_local,coordinate_transform,coordinate_convert,property):
        coord_converted = np.zeros(coord_stack_local.shape)
        suplementary_parameters = utils.return_suplementary_parameters(coordinate_convert,property=property)
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(property.Omega_m)
        coord_converted[:,:,:,0],coord_converted[:,:,:,1],coord_converted[:,:,:,2] = utils.convert_cartesian_to_sky(coord_stack_local[:,:,:,0],
                                                                                                                    coord_stack_local[:,:,:,1],
                                                                                                                    coord_stack_local[:,:,:,2],
                                                                                                                    coordinate_convert,
                                                                                                                    inv_rcomov=inv_rcomov,
                                                                                                                    inv_distang=inv_distang,
                                                                                                                    distang=distang,
                                                                                                                    suplementary_parameters=suplementary_parameters)
        suplementary_parameters = utils.return_suplementary_parameters(coordinate_transform,property=property)
        coord_converted[:,:,:,0],coord_converted[:,:,:,1],coord_converted[:,:,:,2] = utils.convert_sky_to_cartesian(coord_converted[:,:,:,0],
                                                                                                                    coord_converted[:,:,:,1],
                                                                                                                    coord_converted[:,:,:,2],
                                                                                                                    coordinate_transform,
                                                                                                                    rcomov=rcomov,
                                                                                                                    distang=distang,
                                                                                                                    suplementary_parameters=suplementary_parameters)
        return(coord_converted)

    @classmethod
    def inv_convert_center(cls,coord,coordinate_transform,coordinate_convert,property):
        suplementary_parameters = utils.return_suplementary_parameters(coordinate_transform,property=property)
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(property.Omega_m)
        X,Y,Z = utils.convert_cartesian_to_sky(coord[0],
                                               coord[1],
                                               coord[2],
                                               coordinate_transform,
                                               inv_rcomov=inv_rcomov,
                                               inv_distang=inv_distang,
                                               distang=distang,
                                               suplementary_parameters=suplementary_parameters)
        suplementary_parameters = utils.return_suplementary_parameters(coordinate_convert,property=property)
        X,Y,Z = utils.convert_sky_to_cartesian(X,Y,Z,
                                               coordinate_convert,
                                               rcomov=rcomov,
                                               distang=distang,
                                               suplementary_parameters=suplementary_parameters)
        return(np.array([X,Y,Z]))




    def compute_distance_to_los(self,pixel_name):
        pixel = Pixel(name=pixel_name)
        pixel.read()
        x,y,z = pixel.repack_by_los()
        diffz=np.zeros(len(self.catalog.coord.shape[0]))
        for i in range(len(diffz)):
            arg_array  = np.argwhere((self.catalog.coord[i,0] == x[:])&(self.catalog.coord[i,0] == y[:]))
            if(len(arg_array)!=0):
                arg = arg_array[0][0]
                diffz[i] = self.catalog.coord[i,2] - z[arg][-1]
        self.mean_los_distance = np.mean(diffz)


    def compute_stack_ellipticity(self,sign=1):
        ellipticity = {}
        sign = np.sign(self.map_array[self.map_array.shape[0]//2,self.map_array.shape[1]//2,self.map_array.shape[2]//2])
        for direction in ["x","y","z"]:
            x_index,y_index = utils.get_direction_indexes(direction,True)[0:2]
            if(direction == "x"):
                Slice = self.map_array[self.map_array.shape[0]//2,:,:]
            elif(direction == "y"):
                Slice = self.map_array[:,self.map_array.shape[1]//2,:]
            elif(direction == "z"):
                Slice = self.map_array[:,:,self.map_array.shape[2]//2]
            gauss = utils.gaussian_fitter_2d(inpdata = sign * Slice)
            p,success = gauss.FitGauss2D()
            angle = p[5]
            # CR - to check, conversion between shape and size need a -1
            p[1] = (p[1] - Slice.shape[0]//2) * self.mpc_per_pixel[x_index]
            p[2] = -(p[2] - Slice.shape[1]//2) * self.mpc_per_pixel[y_index]
            sigma1 = p[3] * np.sqrt( (np.cos(np.radians(angle)) * self.mpc_per_pixel[x_index])**2 + (np.sin(np.radians(angle)) * self.mpc_per_pixel[y_index])**2)
            sigma2 = p[4] * np.sqrt( (np.cos(np.radians(angle)) * self.mpc_per_pixel[y_index])**2 + (np.sin(np.radians(angle)) * self.mpc_per_pixel[x_index])**2)
            if((angle <45)&(angle>-45)):
                sigma1,sigma2 = sigma1,sigma2
                x_index,y_index = x_index,y_index
            elif((angle <135)&(angle>-135)):
                sigma1,sigma2 = sigma1,sigma2
                x_index,y_index = y_index,x_index
            p[3],p[4] = sigma1,sigma2
            p[5] = np.degrees(np.arctan2(np.tan(np.radians(angle))*self.mpc_per_pixel[y_index],self.mpc_per_pixel[x_index]))
            ellipticity[direction] = sigma2/sigma1
            ellipticity[direction + "_gauss"] = p
            if(direction == "x"):
                ellipticity[direction + "_order"] = "sigmay/sigmaz"
            elif(direction == "y"):
                ellipticity[direction + "_order"] = "sigmax/sigmaz"
            elif(direction == "z"):
                ellipticity[direction + "_order"] = "sigmax/sigmay"
        self.ellipticity = ellipticity


    def add_ellipticity_errors(self,ellipticities):
        for direction in ["x","y","z"]:
            sum_ellipticity = 0
            for i in range(len(ellipticities)):
                sum_ellipticity = sum_ellipticity + (ellipticities[direction] - self.ellipticity[direction])**2
            self.ellipticity[direction +"_error"] = np.sqrt(sum_ellipticity/(len(ellipticities)*(len(ellipticities)-1)))


class Pixel(object):

    def __init__(self,pixel_array=None,name=None,boundary_cartesian_coord=None,boundary_sky_coord=None,coordinate_transform=None,property_file=None,Omega_m=None):
        self.name = name
        self.pixel_array = pixel_array
        self.boundary_cartesian_coord = boundary_cartesian_coord
        self.boundary_sky_coord = boundary_sky_coord
        self.coordinate_transform = coordinate_transform
        self.property_file = property_file
        self.Omega_m = Omega_m

        self.z_array = None
        self.dperp_array = None
        self.density_array = None



    @classmethod
    def init_from_property_files(cls,property_file,pixel_array=None,name=None):
        Property = MapPixelProperty(name=property_file)
        verify_file(property_file)
        Property.read()
        boundary_cartesian_coord,boundary_sky_coord,coordinate_transform = None,None,None
        if(Property.boundary_cartesian_coord is not None): boundary_cartesian_coord = Property.boundary_cartesian_coord
        if(Property.boundary_sky_coord is not None): boundary_sky_coord = Property.boundary_sky_coord
        if(Property.coordinate_transform is not None): coordinate_transform = Property.coordinate_transform
        if(Property.Omega_m is not None): Omega_m = Property.Omega_m
        return(cls(pixel_array=pixel_array,name=name,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,property_file=property_file,Omega_m=Omega_m))



    def read(self):
        if(self.name is None):
            raise ValueError("No name was provided for reading pixel")
        else:
            with open(self.name,'r') as f:
                pixel_data = np.fromfile(f,dtype=np.float64)
                self.pixel_array = pixel_data.reshape((len(pixel_data)//5,5))


    def write(self):
        if(self.pixel_array is None):
            raise ValueError("No pixel array is stored in this pixel class for writting")
        listpixel=np.ravel(self.pixel_array)
        listpixel.tofile(self.name)

    def writetxt(self,name_out):
        x,y,z = self.repack_by_los()
        coord = np.array([[x[i],y[i],z[i][0],z[i][-1]] for i in range(len(x))])
        np.savetxt(name_out,coord)



    def print_prop(self):
        log = utils.create_log()
        if(self.name is not None): log.add(f"Arguments of the property file {self.name}")
        if(self.pixel_array is not None):
            log.add_array_statistics(self.pixel_array[:,0],"X")
            log.add_array_statistics(self.pixel_array[:,1],"Y")
            log.add_array_statistics(self.pixel_array[:,2],"Z")
            log.add_array_statistics(self.pixel_array[:,3],"sigma")
            log.add_array_statistics(self.pixel_array[:,4],"delta")
        if(self.boundary_sky_coord is not None): log.add(f"Sky boundaries of the associated map [radians]: {self.boundary_sky_coord}")
        if(self.boundary_cartesian_coord is not None): log.add(f"Cartesian boundaries of the associated map [Mpc.h-1]: {self.boundary_cartesian_coord}")
        if(self.coordinate_transform is not None): log.add(f"Coordinate transformation of the associated map: {self.coordinate_transform}")
        if(self.Omega_m is not None): log.add(f"Value of Omega_m used: {self.Omega_m}")
        log.close()

    def repack_by_los(self):
        coord = self.pixel_array[:,0:3]
        unique_coord = np.unique(coord[:,0:2],axis=0)
        z = []
        for i in range(unique_coord.shape[0]):
            mask = coord[:,0]==unique_coord[i,0]
            mask2 = coord[:,1][mask] == unique_coord[i,1]
            z.append(coord[:,2][mask][mask2])
        return(unique_coord[:,0],unique_coord[:,1],z)



    def compute_mean_separation(self):
        coord_pixels = self.pixel_array[:,0:3]
        coord_distance = coord_pixels[1:,:] - coord_pixels[0:-1,:]
        distance = np.sqrt(coord_distance[:,0]**2 +coord_distance[:,1]**2 + coord_distance[:,2]**2)
        mask_distance = distance < distance.mean()*10
        return(np.mean(distance[mask_distance]))


    @staticmethod
    def compute_mean_distance(x,y,return_array=False):
        dmin = np.zeros(x.shape)
        dist_map = np.zeros(x.shape)
        for i in range(len(x)):
            dist_map = np.sqrt((x-x[i])**2 + (y-y[i])**2)
            dist_map = dist_map[dist_map>0]
            dmin[i] = np.min(dist_map)
        if(len(x)!=0):
            if(return_array):return(dmin)
            else:return(np.mean(dmin))
        else:return(None)



    @staticmethod
    def compute_mean_distance_density_at_z(zpar,x,y,mini_los,maxi_los):
        mask = (zpar > mini_los)&(zpar < maxi_los)
        xintheplane = x[mask]
        yintheplane = y[mask]
        return(Pixel.compute_mean_distance(xintheplane,yintheplane),len(xintheplane))


    @staticmethod
    def return_mean_distance_density(x,y,z,minra,maxra,mindec,maxdec,minz,maxz):
        zpar = np.linspace(1.05*minz,0.95*maxz,100)
        mini_los = np.array([np.min(z[i]) for i in range(len(z))])
        maxi_los = np.array([np.max(z[i]) for i in range(len(z))])
        dperpz,densityz = [],[]
        for i in range(1,len(zpar)):
            dmin,N_los = Pixel.compute_mean_distance_density_at_z(zpar[i],np.array(x),np.array(y),mini_los,maxi_los)
            dperpz.append(dmin)
            densityz.append(N_los/((maxra-minra)*(maxdec-mindec)*((180/np.pi)**2)))
        return(zpar[1:],dperpz,densityz)


    def compute_mean_distance_density(self):
        if(self.coordinate_transform != "middle"):
            raise KeyError("Mean density and separation only implemented for middle coodinate transformation")
        if((self.z_array is None)|(self.dperp_array is None)|(self.density_array is None)):
            x,y,z = self.repack_by_los()
            minra,mindec = self.boundary_sky_coord[0][0:2]
            maxra,maxdec = self.boundary_sky_coord[1][0:2]
            minz = 0.0
            maxz = self.boundary_cartesian_coord[1][2] - self.boundary_cartesian_coord[0][2]
            z_array,self.dperp_array,self.density_array = Pixel.return_mean_distance_density(x,y,z,minra,maxra,mindec,maxdec,minz,maxz)
            (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(self.Omega_m)
            z_array = z_array + self.boundary_cartesian_coord[0][2]
            redshifts = utils.convert_z_cartesian_to_sky_middle(z_array,inv_rcomov)
            self.z_array = redshifts
        return(self.z_array,self.dperp_array,self.density_array)


    def compute_mean_distance_histogram(self,z_value):
        x,y,z = self.repack_by_los()
        mini_los = np.array([np.min(z[i]) for i in range(len(z))])
        maxi_los = np.array([np.max(z[i]) for i in range(len(z))])
        mask = (z_value > mini_los)&(z_value < maxi_los)
        xintheplane = x[mask]
        yintheplane = y[mask]
        return(Pixel.compute_mean_distance(xintheplane,yintheplane,return_array=True))





   # CR - To modify in merge pixels (by an initialization)

    def merge_maps_and_pixels(self,extremum_list,name_map_out,name_pixel_out,size_maps,mask_maps=None,criteria_mask=None,qso_files=None,name_qso_out=None):
        map_out = []
        pixel_out = []
        if(qso_files is not None):qso_out = []
        if((type(self.MapName) is not list)|(type(self.PixelName) is not list)|(type(extremum_list) is not list)): raise KeyError("Give list to merge")
        for i in range(len(self.MapName)):
            map_data = self.readClamatoMapFile_Other(self.MapName[i],self.shapeMap)
            if(mask_maps is not None):
                map_data = self.mask_map_to_3D(map_data,mask_maps[i],criteria_mask)
            map_out.append(map_data)
            pixel_data = self.readClamatoPixelFile_other(self.PixelName[i])
            if(qso_files is not None): qso =  np.transpose(np.array(pickle.load(open(qso_files[i],'rb'))))
            list_max = pickle.load(open(extremum_list[i],'rb'))
            minx,maxx,miny,maxy,minz,maxz,minredshift,maxredshift,minra,maxra,mindec,maxdec = list_max
            if(i==0): minx_out, minra_out = minx,minra
            if(i==len(self.MapName)-1): maxx_out , maxra_out = maxx, maxra
            pixel_data[:,0],pixel_data[:,1],pixel_data[:,2] = pixel_data[:,0] + minx,pixel_data[:,1] + miny,pixel_data[:,2] + minz
            pixel_out.append(pixel_data)
            print(np.min(qso[:,0]),np.max(qso[:,0]))
            if(qso_files is not None):
                qso[:,0],qso[:,1],qso[:,2] = qso[:,0] + minx,qso[:,1] + miny,qso[:,2] + minz
                qso_out.append(qso)
                print(np.min(qso[:,0]),np.max(qso[:,0]))
        list_max_out = (minx_out,maxx_out,miny,maxy,minz,maxz,minredshift,maxredshift,minra_out,maxra_out,mindec,maxdec)
        pickle.dump(list_max_out,open("list_of_maximums_of_data_cube.pickle","wb"))
        size_map_out = (np.sum([size_maps[i][0] for i in range(len(size_maps))]),np.max([size_maps[i][1] for i in range(len(size_maps))]),np.max([size_maps[i][2] for i in range(len(size_maps))]))
        map_out = np.array(np.concatenate(map_out))
        pixel_out = np.array(np.concatenate(pixel_out))
        pixel_out = pixel_out[pixel_out[:,0].argsort()]
        minpixx = np.min(pixel_out[:,0])
        minpixy = np.min(pixel_out[:,1])
        minpixz = np.min(pixel_out[:,2])
        pixel_out[:,0],pixel_out[:,1],pixel_out[:,2] = pixel_out[:,0] - minpixx,pixel_out[:,1] - minpixy,pixel_out[:,2] - minpixz
        qso_out = np.array(np.concatenate(qso_out))
        qso_out = qso_out[qso_out[:,0].argsort()]
        qso_out[:,0],qso_out[:,1],qso_out[:,2] =  qso_out[:,0] - minpixx,qso_out[:,1] - minpixy,qso_out[:,2] - minpixz
        print(np.min(qso_out[:,0]),np.max(qso_out[:,0]),len(qso_out))
        np.savetxt(name_qso_out,qso_out)
        self.create_Report()
        self.add_Report("Writting the map...")
        self.writeClamatoMapFile(name_map_out,map_out)
        self.add_Report("Writting the pixels...")
        self.writeClamatoPixelFile(name_pixel_out,pixel_out)
        self.add_Report("Output map size = {}".format(size_map_out))
        self.add_Report("Output map shape = {}".format(map_out.shape))
        if(self.shapeMap is None): self.shapeMap = map_out.shape
        return(size_map_out)




class MapPixelProperty(object):


    def __init__(self,name=None,size=None,shape=None,boundary_cartesian_coord=None,boundary_sky_coord=None,coordinate_transform=None,Omega_m=None):
        self.name = name
        self.size = size
        self.shape = shape
        self.boundary_cartesian_coord = boundary_cartesian_coord
        self.boundary_sky_coord = boundary_sky_coord
        self.coordinate_transform = coordinate_transform
        self.Omega_m = Omega_m


    @classmethod
    def init_false_prop(cls,shape,size,Omega_m,minx,miny,minredshift,coordinate_transform,name="property_file.pickle"):
        if(coordinate_transform.lower() == "middle"):
            (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(Omega_m)
            minz = utils.convert_z_sky_to_cartesian_middle(np.array([minredshift]),rcomov)[0]
            maxz = minz + size[2]
            maxredshift = utils.convert_z_cartesian_to_sky_middle(np.array([maxz]),inv_rcomov)[0]
            maxx = minx + size[0]
            maxy = miny + size[1]
            suplementary_parameters = utils.return_suplementary_parameters(coordinate_transform,zmin=minredshift,zmax=maxredshift)
            RA,DEC,redshift = utils.convert_cartesian_to_sky(np.array([minx,maxx]),np.array([miny,maxy]),np.array([minz,maxz]),coordinate_transform,inv_rcomov=inv_rcomov,inv_distang=inv_distang,distang=distang,suplementary_parameters=suplementary_parameters)
            minra,maxra = RA[0],RA[1]
            mindec,maxdec = DEC[0],DEC[1]
        else:
            raise NotImplementedError("Only middle coordinate transformation is implemented")
        boundary_cartesian_coord = ((minx,miny,minz),(maxx,maxy,maxz))
        boundary_sky_coord = ((minra,mindec,minredshift),(maxra,maxdec,maxredshift))
        return(cls(name=name,size=size,shape=shape,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,Omega_m=Omega_m))


    @property
    def mpc_per_pixel(self):
        return(utils.mpc_per_pixel(self.size,self.shape))

    @property
    def pixel_per_mpc(self):
        return(utils.pixel_per_mpc(self.size,self.shape))

    def read(self):
        if(self.name is None):
            raise ValueError("No name was given")
        else:
            file = pickle.load(open(self.name,'rb'))
            self.size = file["size"]
            self.shape = file["shape"]
            self.boundary_cartesian_coord = file["boundary_cartesian_coord"]
            self.boundary_sky_coord = file["boundary_sky_coord"]
            self.coordinate_transform = file["coordinate_transform"]
            self.Omega_m = file["Omega_m"]


    def write(self):
        if(self.name is None):
            raise ValueError("No name was given")
        else:
            dict_prop = {"size":self.size ,"shape": self.shape,"boundary_cartesian_coord": self.boundary_cartesian_coord,"boundary_sky_coord": self.boundary_sky_coord,"coordinate_transform": self.coordinate_transform,"Omega_m":self.Omega_m}
            pickle.dump(dict_prop,open(self.name,'wb'))


    def print_prop(self):
        log = utils.create_log()
        if(self.name is not None): log.add(f"Arguments of the property file {self.name}")
        if(self.shape is not None): log.add(f"Shape of the associated map: {self.shape}")
        if(self.size is not None): log.add(f"Size of the associated map: {self.size}")
        if(self.boundary_sky_coord is not None): log.add(f"Sky boundaries of the associated map [radians]: {self.boundary_sky_coord}")
        if(self.boundary_cartesian_coord is not None): log.add(f"Cartesian boundaries of the associated map [Mpc.h-1]: {self.boundary_cartesian_coord}")
        if(self.coordinate_transform is not None): log.add(f"Coordinate transformation of the associated map: {self.coordinate_transform}")
        if(self.Omega_m is not None): log.add(f"Value of Omega_m used: {self.Omega_m}")
        log.close()



#############################################################################
#############################################################################
############################### DELTA #####################################
#############################################################################
#############################################################################




class Delta(object):

    def __init__(self,delta_array=None,name=None,pk1d_type=True):
        self.name = name
        self.delta_array = delta_array
        self.pk1d_type = pk1d_type


    def read(self):
        delta_array = []
        if(self.name is None):
            raise ValueError("No name was provided for reading delta")
        else:
            with fitsio.FITS(self.name) as delta_file:
                for i in range(1,len(delta_file)):
                    delta_array.append(self.read_line(delta_file,i))
                self.delta_array = delta_array
                delta_file.close()

    def read_line(self,delta_file,number_line):
        if(delta_file is None):
            self.read()
        delta = data.Delta.from_fitsio(delta_file[number_line],pk1d_type=self.pk1d_type)
        return(delta)


    def return_params(self,center_ra=True):
        ra, dec, z, delta, sigma,zqso,id = [],[],[],[],[],[],[]
        for i in range(len(self.delta_array)):
            zqso.append(Delta.z_qso(self.delta_array[i]))
            if((Delta.ra(self.delta_array[i]) > np.pi)&(center_ra)):
                ra.append(Delta.ra(self.delta_array[i])-2*np.pi)
            else :
                ra.append(Delta.ra(self.delta_array[i]))
            dec.append(Delta.dec(self.delta_array[i]))
            if(self.delta_array[i].ivar is not None):
                sigma.append(1/np.sqrt(np.asarray(self.delta_array[i].ivar)))
            else :
                sigma.append(np.array([0.0 for i in range(len(self.delta_array[i].delta))]))
            delta.append(self.delta_array[i].delta)
            id.append(Delta.primary_key(self.delta_array[i]))
            z.append(((10**self.delta_array[i].log_lambda / utils.lambdaLy)-1))

        return(np.array(ra),
               np.array(dec),
               z,
               np.array(zqso),
               np.array(id),
               sigma,
               delta)

    # CR - writting need to be adapted to picca
    def write(self):
        if((self.name is None)|(self.delta_array is None)):
            raise ValueError("No delta array or name is stored in this delta class for writting")
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        for i in range(len(self.delta_array)):
            delta =self.delta_array[i]
            fi,head = self.create_fi(delta)
            fits.write(fi,header=head,extname=str(head['THING_ID']))

    def create_fi(self,delta):
        nrows = len(delta.delta)
        head = {}
        if  self.pk1d_type :
            h = np.zeros(nrows, dtype=[('LOGLAM','f8'),('DELTA','f8'),('IVAR','f8'),('DIFF','f8')])
            h['DELTA'] =delta.delta
            h['LOGLAM'] = delta.log_lambda
            h['IVAR'] = delta.ivar
            h['DIFF'] = delta.exposures_diff
            head['MEANSNR'] = delta.mean_snr
            head['MEANRESO'] = delta.mean_reso
            head['MEANZ'] = delta.mean_z
            head['DLL'] = delta.delta_log_lambda
        else :
            h = np.zeros(nrows, dtype=[('LOGLAM','f8'),('DELTA','f8'),('WEIGHT','f8'),('CONT','f8')])
            h['DELTA'] =delta.delta
            h['LOGLAM'] = delta.log_lambda
            h['WEIGHT'] = delta.weights
            h['CONT'] = delta.cont
        head['THING_ID'] = Delta.primary_key(delta)
        head['RA'] = Delta.ra(delta)
        head['DEC'] = Delta.dec(delta)
        head['Z']  = Delta.z_qso(delta)
        head['PLATE'] = delta.plate
        head['MJD'] = delta.mjd
        head['FIBERID'] = delta.fiberid
        return(h,head)

    @staticmethod
    def ra(delta):
        return(delta.ra)

    @staticmethod
    def dec(delta):
        return(delta.dec)

    @staticmethod
    def z_qso(delta):
        return(delta.z_qso)


    @staticmethod
    def primary_key(delta):
        try:
            return(delta.los_id)
        except:
            return(delta.thingid)



#############################################################################
#############################################################################
############################### CATALOGS #####################################
#############################################################################
#############################################################################


# CR - for the cutting routines, add a routine which automaticaly cut additive arrays




class Catalog(object):


    def __init__(self,name=None,coord=None,primary_key=None,catalog_type="sky",coordinate_transform=None,Omega_m=None,boundary_cartesian_coord=None,boundary_sky_coord=None):
        self.name = name
        self.coord = coord
        self.primary_key = primary_key
        self.catalog_type = catalog_type
        self.coordinate_transform = coordinate_transform
        self.Omega_m = Omega_m
        self.boundary_cartesian_coord=boundary_cartesian_coord
        self.boundary_sky_coord=boundary_sky_coord

    @classmethod
    def init_catalog_from_fits(cls,name,type_catalog):
        if(type_catalog.lower() == "qso"):
            return(QSOCatalog.init_from_fits(name))
        if(type_catalog.lower() == "void"):
            return(VoidCatalog.init_from_fits(name))
        if(type_catalog.lower() == "galaxy"):
            return(QSOCatalog.init_from_fits(name))
        if(type_catalog.lower() == "dla"):
            return(DLACatalog.init_from_fits(name))



    @staticmethod
    def load_from_fits(name):
        if(name == None):
            raise ValueError("No name was provided for reading catalog")
        else:
            return(fitsio.FITS(name))

    @staticmethod
    def load_header(catalog):
        coordinate_transform,boundary_cartesian_coord,boundary_sky_coord, Omega_m = None, None, None, None
        if("COORD_TRANSFORM" in catalog[1].read_header()): coordinate_transform = catalog[1].read_header()["COORD_TRANSFORM"]
        if("OMEGA_M" in catalog[1].read_header()): Omega_m = catalog[1].read_header()["OMEGA_M"]
        if("RAMAX" in catalog[1].read_header()):
            ramax = catalog[1].read_header()["RAMAX"]
            ramin = catalog[1].read_header()["RAMIN"]
            decmax = catalog[1].read_header()["DECMAX"]
            decmin = catalog[1].read_header()["DECMIN"]
            redshiftmax = catalog[1].read_header()["REDMAX"]
            redshiftmin = catalog[1].read_header()["REDMIN"]
            boundary_sky_coord = ((ramin,decmin,redshiftmin),(ramax,decmax,redshiftmax))
        if("XMAX" in catalog[1].read_header()):
            xmax = catalog[1].read_header()["XMAX"]
            xmin = catalog[1].read_header()["XMIN"]
            ymax = catalog[1].read_header()["YMAX"]
            ymin = catalog[1].read_header()["YMIN"]
            zmax = catalog[1].read_header()["ZMAX"]
            zmin = catalog[1].read_header()["ZMIN"]
            boundary_cartesian_coord = ((xmin,ymin,zmin),(xmax,ymax,zmax))
        return(coordinate_transform,boundary_cartesian_coord,boundary_sky_coord, Omega_m)


    @staticmethod
    def close(catalog):
        catalog.close()


    def return_header(self):
        head = {}
        if(self.coordinate_transform is not None):
            head["COORD_TRANSFORM"] = self.coordinate_transform
        if(self.Omega_m is not None):
            head["OMEGA_M"] = self.Omega_m
        if(self.boundary_cartesian_coord is not None):
            head["XMAX"] = self.boundary_cartesian_coord[1][0]
            head["XMIN"] = self.boundary_cartesian_coord[0][0]
            head["YMAX"] = self.boundary_cartesian_coord[1][1]
            head["YMIN"] = self.boundary_cartesian_coord[0][1]
            head["ZMAX"] = self.boundary_cartesian_coord[1][2]
            head["ZMIN"] = self.boundary_cartesian_coord[0][2]
        if(self.boundary_sky_coord is not None):
            head["RAMAX"] = self.boundary_sky_coord[1][0]
            head["RAMIN"] = self.boundary_sky_coord[0][0]
            head["DECMAX"] = self.boundary_sky_coord[1][1]
            head["DECMIN"] = self.boundary_sky_coord[0][1]
            head["REDMAX"] = self.boundary_sky_coord[1][2]
            head["REDMIN"] = self.boundary_sky_coord[0][2]
        return(head)



    def convert_to_absolute_coordinates(self,property_file=None):
        if(property_file is not None):
            print("A property file is used instead of the parameters of the catalog")
            prop = MapPixelProperty(name=property_file)
            prop.read()
            boundary_cartesian_coord = prop.boundary_cartesian_coord
            boundary_sky_coord = prop.boundary_sky_coord
        else:
            boundary_cartesian_coord = self.boundary_cartesian_coord
            boundary_sky_coord = self.boundary_sky_coord
        if(self.catalog_type=="cartesian"):
            boundary = boundary_cartesian_coord
        if(self.catalog_type=="sky"):
            boundary = boundary_sky_coord
        self.coord[:,0] = self.coord[:,0]  + boundary[0][0]
        self.coord[:,1] = self.coord[:,1]  + boundary[0][1]
        self.coord[:,2] = self.coord[:,2]  + boundary[0][2]

    def convert_to_normalized_coordinates(self,property_file=None):
        if(property_file is not None):
            print("A property file is used instead of the parameters of the catalog")
            prop = MapPixelProperty(name=property_file)
            prop.read()
            boundary_cartesian_coord = prop.boundary_cartesian_coord
            boundary_sky_coord = prop.boundary_sky_coord
        else:
            boundary_cartesian_coord = self.boundary_cartesian_coord
            boundary_sky_coord = self.boundary_sky_coord
        if(self.catalog_type=="cartesian"):
            boundary = boundary_cartesian_coord
        if(self.catalog_type=="sky"):
            boundary = boundary_sky_coord
        self.coord[:,0] = self.coord[:,0]  - boundary[0][0]
        self.coord[:,1] = self.coord[:,1]  - boundary[0][1]
        self.coord[:,2] = self.coord[:,2]  - boundary[0][2]

    def convert_to_sky(self,property_file=None):
        if(property_file is not None):
            print("A property file is used instead of the parameters of the catalog")
            prop = MapPixelProperty(name=property_file)
            prop.read()
            coordinate_transform = prop.coordinate_transform
            Omega_m = prop.Omega_m
            zmin = prop.boundary_sky_coord[0][2]
            zmax = prop.boundary_sky_coord[1][2]
        else:
            coordinate_transform = self.coordinate_transform
            Omega_m = self.Omega_m
            zmin = self.boundary_sky_coord[0][2]
            zmax = self.boundary_sky_coord[1][2]
        self.convert_to_absolute_coordinates(property_file=property_file)
        if((coordinate_transform.lower()=="middle")|(coordinate_transform.lower()=="full_angle")|(coordinate_transform.lower()=="full")):
            suplementary_parameters = utils.return_suplementary_parameters(coordinate_transform,zmin=zmin,zmax=zmax)
            (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(Omega_m)
            self.coord[:,0],self.coord[:,1],self.coord[:,2] = utils.convert_cartesian_to_sky(self.coord[:,0],self.coord[:,1],self.coord[:,2],coordinate_transform,inv_rcomov=inv_rcomov,inv_distang=inv_distang,distang=distang,suplementary_parameters=suplementary_parameters)
            self.catalog_type = "sky"
        else:
            raise KeyError("Conversion mode not available, please choose between : middle, full or full_angle")

    def convert_to_cartesian(self,property_file=None):
        if(property_file is not None):
            print("A property file is used instead of the parameters of the catalog")
            prop = MapPixelProperty(name=property_file)
            prop.read()
            coordinate_transform = prop.coordinate_transform
            Omega_m = prop.Omega_m
            zmin = prop.boundary_sky_coord[0][2]
            zmax = prop.boundary_sky_coord[1][2]
        else:
            coordinate_transform = self.coordinate_transform
            Omega_m = self.Omega_m
            zmin = self.boundary_sky_coord[0][2]
            zmax = self.boundary_sky_coord[1][2]
        if((coordinate_transform.lower()=="middle")|(coordinate_transform.lower()=="full_angle")|(coordinate_transform.lower()=="full")):
            suplementary_parameters = utils.return_suplementary_parameters(coordinate_transform,zmin=zmin,zmax=zmax)
            (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(Omega_m)
            self.coord[:,0],self.coord[:,1],self.coord[:,2] = utils.convert_sky_to_cartesian(self.coord[:,0],self.coord[:,1],self.coord[:,2],coordinate_transform,rcomov=rcomov,distang=distang,suplementary_parameters=suplementary_parameters)
            self.catalog_type = "cartesian"
        else:
            raise KeyError("Conversion mode not available, please choose between : middle, full or full_angle")
        self.convert_to_normalized_coordinates(property_file=property_file)


    # CR - add center and degree as a catalog parameter (always centered and in radians) + add conversion in the other way
    def convert_coordinates(self,center=False,decenter=False,degree=False,radians=True):
        if(self.catalog_type != "sky"):
            return()
        if(decenter) :
            mask = self.coord[:,0] < 0
            self.coord[:,0][mask] = self.coord[:,0][mask] + 2*np.pi
        if(degree):
            self.coord[:,0] = np.degrees(self.coord[:,0])
            self.coord[:,1] = np.degrees(self.coord[:,1])

    def move_axis(self,moveaxis):
        for i in range(len(moveaxis)):
            self.coord[:,i] = self.coord[:,i] - moveaxis[i]

    @property
    def redshift(self):
        if(self.catalog_type == "cartesian"):
            self.convert_to_sky()
            redshift = self.coord[:,2].copy()
            self.convert_to_cartesian()
        else:
            redshift = self.coord[:,2].copy()
        return(redshift)


    def cut_catalog(self,coord_min=None,coord_max=None,center_x_coord=False):
        if(self.coord.shape[1] != 3):
            raise ValueError("The catalog is not at the good shape for using this function")
        x = self.coord.copy()[:,0]
        y = self.coord.copy()[:,1]
        z = self.coord.copy()[:,2]
        if(center_x_coord):
            mask = x > 180
            x[mask] = x[mask] - 360
        mask_select = np.full(self.coord.shape[0],True)
        if(coord_min is not None):
            mask_select &= (x > coord_min[0]) & (y > coord_min[1]) & (z > coord_min[2])
        if(coord_max is not None):
            mask_select &= (x<coord_max[0])  & (y<coord_max[1]) & (z<coord_max[2])
        return(mask_select)



    def cut_z_catalog(self,z_min=None,z_max=None):
        z = self.coord.copy()
        mask_select = np.full(self.coord.shape,True)
        if(z_min is not None):
            mask_select &= (z > z_min)
        if(z_max is not None):
            mask_select &= (z < z_max)
        return(mask_select)

    def print_statistics(self,close=True):
        log = utils.create_log()
        if(self.name is not None): log.add(f"Arguments of the catalog {self.name}")
        if(self.catalog_type is not None): log.add(f"Type of catalog: {self.catalog_type}")
        if(self.primary_key is not None):
            log.add(f"An example of primary key: {self.primary_key[0]}")
        if(self.coord is not None):
            log.add(f"Number of objects in the catalog: {self.coord.shape[0]}")
            log.add_array_statistics(self.coord[:,0],"coord over the first direction")
            log.add_array_statistics(self.coord[:,1],"coord over the second direction")
            log.add_array_statistics(self.coord[:,2],"coord over the third direction")
        if(close):
            log.close()
        else:
            return(log)

    def convert_to_cross_corr_radec(self):
        if self.catalog_type == "cartesian":
            self.convert_to_sky()
        self.convert_coordinates(decenter=True,degree=True)


    def create_randomized_catalog(self,
                                  name_out,
                                  property_file,
                                  seed=None,
                                  randomized_z_distribution=False):
        prop = MapPixelProperty(name= property_file)
        prop.read()
        if(self.catalog_type == "sky"):
            (minx,miny,minz) = prop.boundary_sky_coord[0]
            (maxx,maxy,maxz) = prop.boundary_sky_coord[1]
        elif(self.catalog_type == "cartesian"):
            (minx,miny,minz) = prop.boundary_cartesian_coord[0]
            (maxx,maxy,maxz) = prop.boundary_cartesian_coord[1]
            maxx, minx = maxx - minx, 0.0
            maxy, miny = maxy - miny, 0.0
            maxz, minz = maxz - minz, 0.0
        if(seed is not None):
            np.random.seed(seed)
        for i in range(len(self.coord)):
            self.coord[i,0] = minx + np.random.rand()*(maxx-minx)
            self.coord[i,1] = miny + np.random.rand()*(maxy-miny)
            if(randomized_z_distribution):
                self.coord[i,2] = minz + np.random.rand()*(maxz-minz)
        self.name = name_out
        self.write()




class QSOCatalog(Catalog):

    def __init__(self,name=None,
                      coord=None,
                      primary_key=None,
                      plate=None,
                      modern_julian_date=None,
                      fiber_id=None,
                      redshift_name="Z",
                      catalog_type="sky",
                      coordinate_transform=None,
                      Omega_m=None,
                      boundary_cartesian_coord=None,
                      boundary_sky_coord=None,
                      weights=None):
        super(QSOCatalog,self).__init__(name=name,
                                        coord=coord,
                                        primary_key=primary_key,
                                        catalog_type=catalog_type,
                                        coordinate_transform=coordinate_transform,
                                        Omega_m=Omega_m,
                                        boundary_cartesian_coord=boundary_cartesian_coord,
                                        boundary_sky_coord=boundary_sky_coord)
        self.redshift_name = redshift_name
        self.plate = plate
        self.weights = weights
        self.modern_julian_date = modern_julian_date
        self.fiber_id = fiber_id
        self.object_type = "qso"


    @classmethod
    def init_from_fits(cls,name,redshift_name="Z"):
        catalog = Catalog.load_from_fits(name)
        (coordinate_transform,
         boundary_cartesian_coord,
         boundary_sky_coord,
         Omega_m) = Catalog.load_header(catalog)
        if("RA" in catalog[1].get_colnames()):
            coord_ra = catalog[1]["RA"][:]
            coord_dec = catalog[1]["DEC"][:]
            coord_z = catalog[1][redshift_name][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            primary_key = catalog[1]["THING_ID"][:]
            plate = catalog[1]["PLATE"][:]
            modern_julian_date = catalog[1]["MJD"][:]
            fiber_id = catalog[1]["FIBERID"][:]
            catalog_type = "sky"
            Catalog.close(catalog)
            return(cls(name=name,
                       coord=coord,
                       primary_key=primary_key,
                       plate=plate,
                       modern_julian_date=modern_julian_date,
                       fiber_id=fiber_id,
                       redshift_name=redshift_name,
                       catalog_type=catalog_type,
                       coordinate_transform=coordinate_transform,
                       Omega_m=Omega_m,
                       boundary_sky_coord=boundary_sky_coord,
                       boundary_cartesian_coord=boundary_cartesian_coord))
        if("X" in catalog[1].get_colnames()):
            coord_ra = catalog[1]["X"][:]
            coord_dec = catalog[1]["Y"][:]
            coord_z = catalog[1]['Z'][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            primary_key = catalog[1]["THING_ID"][:]
            catalog_type = "cartesian"
            Catalog.close(catalog)
            return(cls(name=name,
                       coord=coord,
                       primary_key=primary_key,
                       catalog_type=catalog_type,
                       coordinate_transform=coordinate_transform,
                       Omega_m=Omega_m,
                       boundary_sky_coord=boundary_sky_coord,
                       boundary_cartesian_coord=boundary_cartesian_coord))

    @classmethod
    def init_from_pixel_catalog(cls,
                                quasar_pixels,
                                name=None,
                                coordinate_transform=None,
                                Omega_m=None,
                                boundary_cartesian_coord=None,
                                boundary_sky_coord=None,
                                catalog_type="cartesian"):
        if(name is None): name ="qso_catalog.fits"
        coord_ra = quasar_pixels[:,0]
        coord_dec = quasar_pixels[:,1]
        coord_z = quasar_pixels[:,2]
        coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
        primary_key = quasar_pixels[:,3]
        if(catalog_type =="sky"):
            plate = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(primary_key))]).astype("i8")
            mjd = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(primary_key))]).astype("i8")
            fiber = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(primary_key))]).astype("i8")
            weights = np.asarray([1.0 for i in range(len(primary_key))])
        else:
            plate, mjd, fiber, weights = None, None, None, None
        return(cls(name=name,
                   coord=coord,
                   primary_key=primary_key,
                   catalog_type=catalog_type,
                   coordinate_transform=coordinate_transform,
                   Omega_m=Omega_m,
                   boundary_sky_coord=boundary_sky_coord,
                   boundary_cartesian_coord=boundary_cartesian_coord,
                   plate=plate,
                   modern_julian_date=mjd,
                   fiber_id=fiber,
                   weights=weights))

    def cut_catalog_qso(self,coord_min=None,coord_max=None):
        mask_select = self.cut_catalog(coord_min=coord_min,
                                       coord_max=coord_max,
                                       center_x_coord=True)
        self.apply_mask(mask_select)


    def apply_mask(self,mask):
        if(self.coord is not None):
            self.coord = self.coord[mask]
        if(self.primary_key is not None):
            self.primary_key = self.primary_key[mask]
        if(self.plate is not None):
            self.plate = self.plate[mask]
        if(self.modern_julian_date is not None):
            self.modern_julian_date = self.modern_julian_date[mask]
        if(self.fiber_id is not None):
            self.fiber_id = self.fiber_id[mask]
        if(self.weights is not None):
            self.weights = self.weights[mask]

    def write(self,convert=True):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        nrows = self.coord.shape[0]
        head = self.return_header()
        if(self.catalog_type == "sky"):
            if(convert):
                self.convert_to_cross_corr_radec()
            dtype=[('RA','f8'),
                   ('DEC','f8'),
                   (self.redshift_name,'f8'),
                   ('THING_ID','i8'),
                   ('PLATE','i4'),
                   ('MJD','i4'),
                   ('FIBERID','i2')]
            if(self.weights is not None):
                dtype.append(('WEIGHT','f8'))
            h = np.zeros(nrows, dtype = dtype)
            h['RA'] = self.coord[:,0]
            h['DEC'] = self.coord[:,1]
            h[self.redshift_name] = self.coord[:,2]
            h['THING_ID'] =self.primary_key
            h['PLATE'] =self.plate
            h['MJD'] = self.modern_julian_date
            h['FIBERID'] =self.fiber_id
            if(self.weights is not None):
                h['WEIGHT'] =self.weights
        elif(self.catalog_type == "cartesian"):
            h = np.zeros(nrows, dtype=[('X','f8'),('Y','f8'),('Z','f8'),('THING_ID','i8')])
            h['X'] = self.coord[:,0]
            h['Y'] = self.coord[:,1]
            h['Z'] = self.coord[:,2]
            h['THING_ID'] =self.primary_key
        fits.write(h,header=head)
        fits.close()

    def writetxt(self,name_out,moveaxis=None):
        if(moveaxis is not None):
            self.move_axis(moveaxis)
        np.savetxt(name_out,self.coord)

    def cut_write_catalog(self,name_out,coord_min,coord_max):
        self.read_from_fits()
        self.cut_quasars_catalogs(coord_min,coord_max)
        self.name = name_out
        self.write()



class DLACatalog(Catalog):

    def __init__(self,name=None,coord=None,primary_key=None,z_qso=None,
                      confidence=None,nhi=None,catalog_type="sky",
                      coordinate_transform=None,Omega_m=None,
                      boundary_cartesian_coord=None,boundary_sky_coord=None):
        super(DLACatalog,self).__init__(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type,coordinate_transform=coordinate_transform,Omega_m=Omega_m,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord)
        self.z_qso=z_qso
        self.confidence=confidence
        self.nhi=nhi
        self.object_type = "dla"



    @classmethod
    def init_from_fits(cls,name):
        catalog = Catalog.load_from_fits(name)
        (coordinate_transform,boundary_cartesian_coord,boundary_sky_coord, Omega_m) = Catalog.load_header(catalog)
        if("THING_ID" in catalog["DLACAT"].get_colnames()):
            coord = catalog["DLACAT"]["Z"][:]
            primary_key = catalog["DLACAT"]["THING_ID"][:]
            conf_dla = catalog["DLACAT"]["CONF_DLA"][:]
            nhi_dla = catalog["DLACAT"]["NHI"][:]
            catalog_type = "sky"
            if("Z_QSO" in catalog["DLACAT"].get_colnames()):
                z_qso = catalog["DLACAT"]["Z_QSO"][:]
            else:
                z_qso = None
            Catalog.close(catalog)
            return(cls(name=name,coord=coord,z_qso=z_qso,primary_key=primary_key,confidence=conf_dla,nhi=nhi_dla,catalog_type=catalog_type,coordinate_transform=coordinate_transform,Omega_m=Omega_m,boundary_sky_coord=boundary_sky_coord,boundary_cartesian_coord=boundary_cartesian_coord))
        if("X" in catalog[1].get_colnames()):
            coord_ra = catalog["DLACAT"]["X"][:]
            coord_dec = catalog["DLACAT"]["Y"][:]
            coord_z = catalog["DLACAT"]['Z'][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            z_qso = catalog["DLACAT"]['Z_QSO'][:]
            catalog_type = "cartesian"
            Catalog.close(catalog)
            return(cls(name=name,coord=coord,z_qso=z_qso,catalog_type=catalog_type,coordinate_transform=coordinate_transform,Omega_m=Omega_m,boundary_sky_coord=boundary_sky_coord,boundary_cartesian_coord=boundary_cartesian_coord))

    @classmethod
    def init_from_pixel_catalog(cls,dla_pixels,name=None,coordinate_transform=None,Omega_m=None,boundary_cartesian_coord=None,boundary_sky_coord=None):
        if(name is None): name ="dla_catalog.fits"
        coord_ra = dla_pixels[:,0]
        coord_dec = dla_pixels[:,1]
        coord_z = dla_pixels[:,2]
        coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
        z_qso = dla_pixels[:,3]
        catalog_type = "cartesian"
        return(cls(name=name,coord=coord,z_qso=z_qso,catalog_type=catalog_type,coordinate_transform=coordinate_transform,Omega_m=Omega_m,boundary_sky_coord=boundary_sky_coord,boundary_cartesian_coord=boundary_cartesian_coord))



    def cut_catalog_dla(self,coord_min=None,coord_max=None,confidence_min=None,nhi_min=None):
        if(self.catalog_type == "sky"):
            mask_select = self.cut_z_catalog(z_min=coord_min,z_max=coord_max)
        if(self.catalog_type == "cartesian"):
            mask_select = self.cut_catalog(coord_min=coord_min,coord_max=coord_max,center_x_coord=True) & (self.coord[:,2]!=-1.0)
        self.apply_mask(mask_select,confidence_min=confidence_min,nhi_min=nhi_min)


    def apply_mask(self,mask,confidence_min=None,nhi_min=None):
        if(confidence_min is not None):
            mask &= (self.confidence > confidence_min)
        if(nhi_min is not None):
            mask &= (self.nhi > nhi_min)
        if(self.coord is not None):
            self.coord = self.coord[mask]
        if(self.primary_key is not None):
            self.primary_key = self.primary_key[mask]
        if(self.confidence is not None):
            self.confidence = self.confidence[mask]
        if(self.nhi is not None):
            self.nhi = self.nhi[mask]
        if(self.z_qso is not None):
            self.z_qso = self.z_qso[mask]



    def write(self):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        nrows = self.coord.shape[0]
        head = self.return_header()
        if(self.catalog_type == "sky"):
            if(self.z_qso is not None):
                h = np.zeros(nrows, dtype=[('THING_ID','i8'),('Z','f8'),('CONF_DLA','f8'),('NHI','f8'),('Z_QSO','f8')])
                h['Z_QSO'] =self.z_qso
            else:
                h = np.zeros(nrows, dtype=[('THING_ID','i8'),('Z','f8'),('CONF_DLA','f8'),('NHI','f8')])
            h['THING_ID'] =self.primary_key
            h['Z'] = self.coord
            h['CONF_DLA'] = self.confidence
            h['NHI'] = self.nhi
        elif(self.catalog_type == "cartesian"):
            h = np.zeros(nrows, dtype=[('X','f8'),('Y','f8'),('Z','f8'),('Z_QSO','f8')])
            h['X'] = self.coord[:,1]
            h['Y'] = self.coord[:,2]
            h['Z'] = self.coord[:,3]
            h['Z_QSO'] = self.z_qso
        fits.write(h,header=head,extname="DLACAT")
        fits.close()



    @staticmethod
    def convert_old_format_to_new(input,output):
        file = fitsio.FITS(input)["DLA_CAT"]
        thid = file["THING_ID"][:]
        z = file["Z_DLA"][:]
        nhi = file["NHI_DLA"][:]
        conf = file["CONF_DLA"][:]
        thid_out,z_out,nhi_out,conf_out = [],[],[],[]
        for i in range(len(z)):
            mask = z[i] != -1.
            z_out.append(z[i][mask])
            nhi_out.append(nhi[i][mask])
            conf_out.append(conf[i][mask])
            for j in range(len(z[i][mask])):
                thid_out.append(thid[i])
        dla_cat = DLACatalog(name=output,
                             coord=np.concatenate(z_out),
                             primary_key=np.array(thid_out),
                             confidence=np.concatenate(conf_out),
                             nhi=np.concatenate(nhi_out),
                             catalog_type="sky")
        dla_cat.write()


class GalaxyCatalog(Catalog):

    def __init__(self,name=None,
                      coord=None,
                      primary_key=None,
                      coordinate_transform=None,
                      Omega_m=None,
                      boundary_cartesian_coord=None,
                      boundary_sky_coord=None,
                      confidence=None,
                      standard_deviation=None,
                      magnitude=None,
                      catalog_type="sky"):
        super(GalaxyCatalog,self).__init__(name=name,
                                           coord=coord,
                                           primary_key=primary_key,
                                           catalog_type=catalog_type,
                                           coordinate_transform=coordinate_transform,
                                           Omega_m=Omega_m,
                                           boundary_cartesian_coord=boundary_cartesian_coord,
                                           boundary_sky_coord=boundary_sky_coord)

        self.confidence = confidence
        self.standard_deviation = standard_deviation
        self.magnitude = magnitude
        self.object_type = "galaxy"


    @classmethod
    def init_from_fits(cls,name):
        catalog = Catalog.load_from_fits(name)
        (coordinate_transform,boundary_cartesian_coord,boundary_sky_coord, Omega_m) = Catalog.load_header(catalog)
        if("RA" in catalog[1].get_colnames()):
            coord_ra = catalog[1]["RA"][:]
            coord_dec = catalog[1]["DEC"][:]
            coord_z = catalog[1]['Z'][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            primary_key = catalog[1]["THING_ID"][:]
            confidence = catalog[1]['CONF'][:]
            standard_deviation = catalog[1]["STD"][:]
            magnitude = catalog[1]["MAG"][:]
            catalog_type = "sky"
            Catalog.close(catalog)
            return(cls(name=name,
                       coord=coord,
                       primary_key=primary_key,
                       confidence=confidence,
                       standard_deviation=standard_deviation,
                       magnitude=magnitude,
                       catalog_type=catalog_type))

        elif("X" in catalog[1].get_colnames()):
            coord_ra = catalog[1]["X"][:]
            coord_dec = catalog[1]["Y"][:]
            coord_z = catalog[1]['Z'][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            primary_key = catalog[1]["THING_ID"][:]
            standard_deviation = catalog[1]['STD'][:]
            catalog_type = "cartesian"
            Catalog.close(catalog)
            return(cls(name=name,
                       coord=coord,
                       primary_key=primary_key,
                       standard_deviation=standard_deviation,
                       catalog_type=catalog_type))

    @classmethod
    def init_from_hst_cat(cls,name):
        catalog = Catalog.load_from_fits(name)
        dec = np.array(catalog[1]["dec"][:])
        ra = np.array(catalog[1]["ra"][:])
        photoz = np.array(catalog[1]["PHOTOZ_BEST"][:])
        primary_key = np.array(catalog[1]["THING_ID"][:])
        coord = np.vstack([ra,dec,photoz]).transpose()
        confidence = np.array(catalog[1]["PHOTOZ_CONF_BEST"][:])
        standard_deviation = np.array(catalog[1]["PHOTOZ_STD_BEST"][:])
        magnitude = np.array(catalog[1]["z_cmodel_mag"][:])
        Catalog.close(catalog)
        catalog_type = "sky"
        return(cls(name=name,
                   coord=coord,
                   primary_key=primary_key,
                   confidence=confidence,
                   standard_deviation=standard_deviation,
                   magnitude=magnitude,
                   catalog_type=catalog_type))

    def cut_catalog_galaxy(self,
                           coord_min=None,
                           coord_max=None,
                           standard_deviation_max=None,
                           confidence_min=None,
                           magnitude_max=None,
                           distance_map_name=None,
                           distance_map_prop=None,
                           distance_map_param=None):

        mask_select = self.cut_catalog(coord_min=coord_min,
                                       coord_max=coord_max,
                                       center_x_coord=True)
        string_to_add = ""
        if(distance_map_name is not None):
            if ((distance_map_name is None)|(distance_map_prop is None)|(distance_map_param is None)) is True : raise KeyError("Give a dist map parameter, map and property file please")
            mask_select &= self.cut_distance_map(distance_map_name,
                                                 distance_map_prop,
                                                 distance_map_param)
            string_to_add = string_to_add + f"_cutdistance_{distance_map_param}Mpc"


        self.apply_mask(mask_select,
                        standard_deviation_max=standard_deviation_max,
                        confidence_min=confidence_min,
                        magnitude_max=magnitude_max)
        name_out = f"""{self.name.split(".fits")[0]}{string_to_add}.fits"""
        return(name_out)



    def cut_distance_map(self,
                         distance_map_name,
                         distance_map_prop,
                         distance_map_param):
        distance_map = DistanceMap.init_classic(name=distance_map_name,
                                                property_file=distance_map_prop)
        distance_map.read()
        mask = ~distance_map.get_mask_distance(distance_map_param)
        mask_cut = np.full(self.coord.shape[0],False)
        for i in range(len(self.coord)):
            x_coord = int(round((self.coord[i][0]) / distance_map.mpc_per_pixel[0],0))
            y_coord = int(round((self.coord[i][1]) / distance_map.mpc_per_pixel[1],0))
            z_coord = int(round((self.coord[i][2]) / distance_map.mpc_per_pixel[2],0))
            mask_cut[i] = mask[x_coord,y_coord,z_coord]
        return(mask_cut)


    # CR - change all apply_mask with a list of parameters & getattr(class,str)

    def apply_mask(self,mask,standard_deviation_max=None,confidence_min=None,magnitude_max=None):
        if(standard_deviation_max is not None):
            mask &= self.standard_deviation < standard_deviation_max
        if(magnitude_max is not None):
            mask &= self.magnitude < magnitude_max
        if(confidence_min is not None):
            mask &= self.confidence > confidence_min
        if(self.coord is not None):
            self.coord = self.coord[mask]
        if(self.primary_key is not None):
            self.primary_key = self.primary_key[mask]
        if(self.confidence is not None):
            self.confidence = self.confidence[mask]
        if(self.standard_deviation is not None):
            self.standard_deviation = self.standard_deviation[mask]
        if(self.magnitude is not None):
            self.magnitude = self.magnitude[mask]

    def write(self):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        nrows = self.coord.shape[0]
        head = self.return_header()
        if(self.catalog_type == "sky"):
            h = np.zeros(nrows, dtype=[('RA','f8'),('DEC','f8'),('Z','f8'),('THING_ID','i8'),('CONF','i4'),('STD','i4'),('MAG','i2')])
            h['RA'] = self.coord[:,0]
            h['DEC'] = self.coord[:,1]
            h['Z'] = self.coord[:,2]
            h['THING_ID'] =self.primary_key
            h['CONF'] =self.confidence
            h['STD'] = self.standard_deviation
            h['MAG'] =self.magnitude
        elif(self.catalog_type == "cartesian"):
            h = np.zeros(nrows, dtype=[('X','f8'),('Y','f8'),('Z','f8'),('STD','f8'),('THING_ID','i8')])
            h['X'] = self.coord[:,0]
            h['Y'] = self.coord[:,1]
            h['Z'] = self.coord[:,2]
            h['STD'] = self.standard_deviation
            h['THING_ID'] =self.primary_key
        fits.write(h,header=head)
        fits.close()


class VoidCatalog(Catalog):

    def __init__(self,name=None,coord=None,primary_key=None,radius=None,
                      weights=None,crossing_param=None,central_value=None,
                      los_distance = None,filling_factor=None,
                      mean_value=None,catalog_type="sky",
                      coordinate_transform=None,Omega_m=None,
                      boundary_cartesian_coord=None,boundary_sky_coord=None):
        super(VoidCatalog,self).__init__(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type,coordinate_transform=coordinate_transform,Omega_m=Omega_m,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord)
        self.radius = radius
        self.crossing_param = crossing_param
        self.mean_value = mean_value
        self.central_value = central_value
        self.los_distance = los_distance
        self.filling_factor = filling_factor
        self.weights = weights
        self.object_type = "void"



    @classmethod
    def init_from_fits(cls,name):
        catalog = Catalog.load_from_fits(name)
        (coordinate_transform,boundary_cartesian_coord,boundary_sky_coord, Omega_m) = Catalog.load_header(catalog)
        if("RA" in catalog[1].get_colnames()):
            coord_ra = catalog[1]["RA"][:]
            coord_dec = catalog[1]["DEC"][:]
            coord_z = catalog[1]["Z"][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            catalog_type = "sky"
        elif("X" in catalog[1].get_colnames()):
            coord_ra = catalog[1]["X"][:]
            coord_dec = catalog[1]["Y"][:]
            coord_z = catalog[1]['Z'][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            catalog_type = "cartesian"
        radius = catalog[1]["R"][:]
        primary_key = catalog[1]["THING_ID"][:]
        weights = catalog[1]["WEIGHT"][:]
        crossing_param, central_value, mean_value, filling_factor, los_distance = None, None, None, None, None
        if("CROSSING" in catalog[1].get_colnames()): crossing_param = catalog[1]["CROSSING"][:]
        if("VALUE" in catalog[1].get_colnames()): central_value = catalog[1]["VALUE"][:]
        if("MEAN" in catalog[1].get_colnames()): mean_value = catalog[1]["MEAN"][:]
        if("LOS_DIST" in catalog[1].get_colnames()): los_distance = catalog[1]["LOS_DIST"][:]
        if("FILLING_FACTOR" in catalog[1].read_header()): filling_factor = catalog[1].read_header()["FILLING_FACTOR"]

        Catalog.close(catalog)
        return(cls(name=name,coord=coord,primary_key=primary_key,radius=radius,
                   weights=weights,crossing_param=crossing_param,
                   central_value=central_value,filling_factor=filling_factor,
                   los_distance = los_distance,mean_value=mean_value,
                   catalog_type=catalog_type,
                   coordinate_transform=coordinate_transform,
                   Omega_m=Omega_m,boundary_sky_coord=boundary_sky_coord,
                   boundary_cartesian_coord=boundary_cartesian_coord))


    @classmethod
    def init_from_dictionary(cls,
                             name,
                             radius,
                             coord,
                             catalog_type,
                             coordinate_transform,
                             Omega_m,
                             boundary_cartesian_coord,
                             boundary_sky_coord,
                             other_array=None,
                             other_array_name = None):
        central_value, weights, filling_factor, primary_key, crossing_param, mean_value,los_distance = None,None, None, None, None, None, None
        if(other_array_name is not None):
            if("VALUE" in other_array_name):central_value = other_array[np.argwhere("VALUE" == np.asarray(other_array_name))[0][0]]
            if("MEAN" in other_array_name):mean_value = other_array[np.argwhere("MEAN" == np.asarray(other_array_name))[0][0]]
            if("WEIGHT" in other_array_name):weights = other_array[np.argwhere("WEIGHT" == np.asarray(other_array_name))[0][0]]
            if("FILLING_FACTOR" in other_array_name):filling_factor = other_array[np.argwhere("FILLING_FACTOR" == np.asarray(other_array_name))[0][0]]
            if("THING_ID" in other_array_name):primary_key = other_array[np.argwhere("THING_ID" == np.asarray(other_array_name))[0][0]]
            if("CROSSING" in other_array_name):crossing_param = other_array[np.argwhere("CROSSING" == np.asarray(other_array_name))[0][0]]
            if("LOS_DIST" in other_array_name):los_distance = other_array[np.argwhere("LOS_DIST" == np.asarray(other_array_name))[0][0]]
        return(cls(name=name,
                   coord=coord,
                   primary_key=primary_key,
                   radius=radius,
                   weights=weights,
                   crossing_param=crossing_param,
                   central_value=central_value,
                   filling_factor=filling_factor,
                   los_distance = los_distance,
                   mean_value=mean_value,
                   catalog_type=catalog_type,
                   coordinate_transform=coordinate_transform,
                   Omega_m=Omega_m,
                   boundary_sky_coord=boundary_sky_coord,
                   boundary_cartesian_coord=boundary_cartesian_coord))

    @classmethod
    def init_by_merging(cls,catalog_name,name=None):
        catalog = [VoidCatalog.init_from_fits(name) for name in catalog_name]
        catalog_type = catalog[0].catalog_type
        radius,coord,primary_key,crossing_param,central_value,mean_value,weights,filling_factor,los_distance = None,None,None,None,None,None,None,None,None
        if(catalog[0].radius is not None):radius = np.concatenate([cat.radius for cat in catalog])
        if(catalog[0].coord is not None):coord = np.concatenate([cat.coord for cat in catalog])
        if(catalog[0].primary_key is not None):primary_key = np.concatenate([cat.primary_key for cat in catalog])
        if(catalog[0].crossing_param is not None):crossing_param = np.concatenate([cat.crossing_param for cat in catalog])
        if(catalog[0].central_value is not None):central_value = np.concatenate([cat.central_value for cat in catalog])
        if(catalog[0].mean_value is not None):mean_value = np.concatenate([cat.mean_value for cat in catalog])
        if(catalog[0].los_distance is not None):los_distance = np.concatenate([cat.los_distance for cat in catalog])
        if(catalog[0].weights is not None):weights = np.concatenate([cat.weights for cat in catalog])
        filling_factor_boolean = np.array([catalog[i].filling_factor is not None for i in range(len(catalog))])
        if(len(filling_factor_boolean[filling_factor_boolean==False])==0):
            filling_factor = np.mean([cat.filling_factor for cat in catalog])
        if(catalog[0].coordinate_transform is not None):coordinate_transform = catalog[0].coordinate_transform
        if(catalog[0].Omega_m is not None):Omega_m = catalog[0].Omega_m
        if(catalog[0].boundary_sky_coord is not None):
            minra = np.min([cat.boundary_sky_coord[0][0] for cat in catalog])
            maxra = np.min([cat.boundary_sky_coord[1][0] for cat in catalog])
            mindec = np.min([cat.boundary_sky_coord[0][1] for cat in catalog])
            maxdec = np.min([cat.boundary_sky_coord[1][1] for cat in catalog])
            minredshift = np.min([cat.boundary_sky_coord[0][2] for cat in catalog])
            maxredshift = np.min([cat.boundary_sky_coord[1][2] for cat in catalog])
            boundary_sky_coord = ((minra,mindec,minredshift),(maxra,maxdec,maxredshift))
        if(catalog[0].boundary_cartesian_coord is not None):
            minx = np.min([cat.boundary_cartesian_coord[0][0] for cat in catalog])
            maxx = np.min([cat.boundary_cartesian_coord[1][0] for cat in catalog])
            miny = np.min([cat.boundary_cartesian_coord[0][1] for cat in catalog])
            maxy = np.min([cat.boundary_cartesian_coord[1][1] for cat in catalog])
            minz = np.min([cat.boundary_cartesian_coord[0][2] for cat in catalog])
            maxz = np.min([cat.boundary_cartesian_coord[1][2] for cat in catalog])
            boundary_cartesian_coord = ((minx,miny,minz),(maxx,maxy,maxz))
        return(cls(name=name,coord=coord,primary_key=primary_key,radius=radius,
                   weights=weights,crossing_param=crossing_param,
                   central_value=central_value,filling_factor=filling_factor,
                   los_distance = los_distance,mean_value=mean_value,
                   catalog_type=catalog_type,
                   coordinate_transform=coordinate_transform,
                   Omega_m=Omega_m,boundary_sky_coord=boundary_sky_coord,
                   boundary_cartesian_coord=boundary_cartesian_coord))


    def writetxt(self,name_out,moveaxis=None):
        if(moveaxis is not None):
            self.move_axis(moveaxis)
        coord = np.transpose(np.stack([self.coord[:,0],self.coord[:,1],self.coord[:,2],self.radius]))
        np.savetxt(name_out,coord)


    def write(self,qso_like=False):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        head = self.return_header()
        if(self.catalog_type.lower() == "sky"):
            h = self.create_sky_void_dictionary(qso_like=qso_like)
        elif(self.catalog_type.lower() == "cartesian"):
            h = self.create_cartesian_void_dictionary()
        self.update_dictionary(h,head)
        fits.write(h,header=head)
        fits.close()

    def create_sky_void_dictionary(self,qso_like=False):
        h = {}
        h['RA'] = self.coord[:,0].astype("f8") if self.coord.shape[0] != 0 else np.array([])
        h['DEC'] = self.coord[:,1].astype("f8") if self.coord.shape[0] != 0 else np.array([])
        h['Z'] = self.coord[:,2].astype("f8") if self.coord.shape[0] != 0 else np.array([])
        if(qso_like):
            h['PLATE'] = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(self.radius))]).astype("i8")
            h['MJD'] = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(self.radius))]).astype("i8")
            h['FIBERID'] = np.asarray(['1' + '{0:05d}'.format(i) for i in range(len(self.radius))]).astype("i8")
        return(h)

    def create_cartesian_void_dictionary(self):
        h = {}
        h['X'] = self.coord[:,0].astype("f8") if self.coord.shape[0] != 0 else np.array([])
        h['Y'] = self.coord[:,1].astype("f8") if self.coord.shape[0] != 0 else np.array([])
        h['Z'] = self.coord[:,2].astype("f8") if self.coord.shape[0] != 0 else np.array([])
        return(h)


    def update_dictionary(self,h,head):
        h['R'] = self.radius.astype("f8")
        if(self.primary_key is not None): h['THING_ID'] =np.array(self.primary_key).astype("i8")
        else:h['THING_ID'] = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(self.radius))]).astype("i8")
        if(self.weights is not None): h['WEIGHT'] = np.array(self.weights).astype("f8")
        else:h['WEIGHT'] = np.full(self.radius.shape,1.0)
        if(self.crossing_param  is not None ):
            h["CROSSING"] = np.array(self.crossing_param).astype("f8")
        if(self.central_value is not None):
            h["VALUE"] = np.array(self.central_value).astype("f8")
        if(self.mean_value is not None):
            h["MEAN"] = np.array(self.mean_value).astype("f8")
        if(self.los_distance is not None):
            h["LOS_DIST"] = np.array(self.los_distance).astype("f8")
        if(self.filling_factor is not None):
            head["FILLING_FACTOR"] = self.filling_factor



    def print_void_statistics(self):
        log = self.print_statistics(close=False)
        if(self.radius is not None):
            log.add_array_statistics(self.radius,"void radius")
        if(self.crossing_param is not None):
            log.add_array_statistics(self.crossing_param,"void pixel crossing parameter")
        if(self.central_value is not None):
            log.add_array_statistics(self.central_value,"void central value")
        if(self.mean_value is not None):
            log.add_array_statistics(self.mean_value,"void average value")
        if(self.los_distance is not None):
            log.add_array_statistics(self.los_distance,"void average min distance to los")
        if(self.weights is not None):
            log.add_array_statistics(self.weights,"void weights")
        if(self.filling_factor is not None):
            log.add(f"Filling factor of the void catalog: {self.filling_factor}")
        if(self.coordinate_transform is not None):
            log.add(f"Coordinate transformation of the void catalog: {self.coordinate_transform}")
        if(self.Omega_m is not None):
            log.add(f"Omega_m used for the void catalog: {self.Omega_m}")





    def cut_catalog_void(self,
                         method_cut,
                         coord_min=None,
                         coord_max=None,
                         cut_crossing_param=None,
                         cut_radius=None,
                         distance_map_name=None,
                         distance_map_prop=None,
                         distance_map_param=None,
                         distance_map_percent=None):
        mask_select = self.cut_catalog(coord_min=coord_min,
                                       coord_max=coord_max,
                                       center_x_coord=False)
        string_to_add = ""
        if(type(method_cut) ==list):
            method_cut = tuple(method_cut)
        if(type(method_cut) != tuple):
            method_cut = tuple(str(k.strip()) for k in method_cut.strip().split(','))
        if method_cut == ("ALL",):
            method_cut = ("CROSSING","RADIUS","DIST","BORDER")

        if("CROSSING" in method_cut):
            if (cut_crossing_param is None) : raise KeyError("Crossing parameter is not computed. Please do so with compute_additional_stats()")
            mask_select &= self.cut_crossing_parameter(cut_crossing_param)
            string_to_add = string_to_add + f"_crossing{cut_crossing_param}"
        if("RADIUS" in method_cut):
            if (cut_radius is None) : raise KeyError("Give a radius cutting parameter")
            mask_select &= self.cut_radius(cut_radius)
            string_to_add = string_to_add + f"_cutradius_{cut_radius[0]}rmin_{cut_radius[1]}rmax"
        if("BORDER" in method_cut):
            mask_select &= self.cut_border()
            string_to_add = string_to_add + "_cutborder"
        if("DIST" in method_cut):
            if ((distance_map_name is None)|(distance_map_prop is None)|(distance_map_param is None)|(distance_map_percent is None)) : raise KeyError("Give a dist map parameter, map and property file please")
            mask_select &= self.cut_distance_map(distance_map_name,distance_map_prop,distance_map_param,distance_map_percent)
            string_to_add = string_to_add + f"_cutdistance_{distance_map_param}Mpc_{distance_map_percent}percent"
        self.apply_mask(mask_select)
        name_out = f"""{self.name.split(".fits")[0]}{string_to_add}.fits"""
        return(name_out)


    def apply_mask(self,mask):
        if(self.radius is not None):self.radius = self.radius[mask]
        if(self.coord is not None):self.coord = self.coord[mask]
        if(self.primary_key is not None):self.primary_key = self.primary_key[mask]
        if(self.crossing_param is not None):self.crossing_param = self.crossing_param[mask]
        if(self.central_value is not None):self.central_value = self.central_value[mask]
        if(self.los_distance is not None):self.los_distance = self.los_distance[mask]
        if(self.mean_value is not None):self.mean_value = self.mean_value[mask]
        if(self.weights is not None):self.weights = self.weights[mask]


    def cut_crossing_parameter(self,cut_crossing_param):
        mask = self.crossing_param > cut_crossing_param
        return(mask)

    def cut_radius(self,cut_radius):
        mask = (self.radius >= cut_radius[0])&(self.radius < cut_radius[1])
        return(mask)

    def cut_border(self):
        cut_border_size = np.array(self.boundary_cartesian_coord[1]) - np.array(self.boundary_cartesian_coord[0])
        mask = (self.coord[:,0] - self.radius < 0)|(self.coord[:,0] + self.radius > cut_border_size[0])
        mask |= (self.coord[:,1] - self.radius < 0)|(self.coord[:,1] + self.radius > cut_border_size[1])
        mask |= (self.coord[:,2] - self.radius < 0)|(self.coord[:,2] + self.radius > cut_border_size[2])
        return(~mask)


    def cut_distance_map(self,distance_map_name,distance_map_prop,distance_map_param,distance_map_percent):
        distance_map = DistanceMap.init_classic(name=distance_map_name,property_file=distance_map_prop)
        distance_map.read()
        mask = ~distance_map.get_mask_distance(distance_map_param)
        mask_cut = np.full(self.radius.shape,False)
        for i in range(len(self.coord)):
            x_coord_min  = max(int(round((self.coord[i][0] - self.radius[i]) / distance_map.mpc_per_pixel[0],0)),0)
            x_coord_max  = min(int(round((self.coord[i][0] + self.radius[i]) / distance_map.mpc_per_pixel[0],0)),distance_map.shape[0]-1)
            y_coord_min  = max(int(round((self.coord[i][1] - self.radius[i]) / distance_map.mpc_per_pixel[1],0)),0)
            y_coord_max  = min(int(round((self.coord[i][1] + self.radius[i]) / distance_map.mpc_per_pixel[1],0)),distance_map.shape[1]-1)
            z_coord_min  = max(int(round((self.coord[i][2] - self.radius[i]) / distance_map.mpc_per_pixel[2],0)),0)
            z_coord_max  = min(int(round((self.coord[i][2] + self.radius[i]) / distance_map.mpc_per_pixel[2],0)),distance_map.shape[2]-1)
            sub_mask = mask[x_coord_min:x_coord_max+1,y_coord_min:y_coord_max+1,z_coord_min:z_coord_max+1]
            if(len(sub_mask[sub_mask == True])/len(sub_mask[sub_mask != None]) >= distance_map_percent):
                mask_cut[i] = True
        return(mask_cut)



    def return_array_list(self,other_array_name):
        other_array = []
        if(other_array_name is not None):
            for i in range(len(other_array_name)):
                if(other_array_name[i] == "VALUE"):other_array.append(self.central_value)
                elif(other_array_name[i] == "MEAN"):other_array.append(self.mean_value)
                elif(other_array_name[i] == "WEIGHT"):other_array.append(self.weights)
                elif(other_array_name[i] == "FILLING_FACTOR"):other_array.append(self.filling_factor)
                elif(other_array_name[i] == "THING_ID"):other_array.append(self.primary_key)
                elif(other_array_name[i] == "CROSSING"):other_array.append(self.crossing_param)
                elif(other_array_name[i] == "LOS_DIST"):other_array.append(self.los_distance)
                else: raise KeyError(f"{other_array_name[i]} not available for void catalog")
        return(other_array)


    ### Computing & specific functions ###




    def compute_filling_factor(self):
        if(self.filling_factor is not None):
            return()
        map_size = np.array(self.boundary_cartesian_coord[1]) - np.array(self.boundary_cartesian_coord[0])
        volume_map = map_size[0]*map_size[1]*map_size[2]
        volume_void = np.sum((4/3)*np.pi*(np.array(self.radius))**3)
        self.filling_factor = volume_void/volume_map


    def compute_crossing_criteria(self,pixel_name):
        if(self.crossing_param is not None):
            return()
        pixel = Pixel(name=pixel_name)
        pixel.read()
        coord_pixels = pixel.pixel_array[:,0:3]
        separation_pixels = pixel.compute_mean_separation()
        crossing_param = np.zeros(self.radius.shape)
        for i in range(len(self.radius)):
            diff_pixel = coord_pixels - self.coord[i]
            distance_pixel = np.sqrt(diff_pixel[:,0]**2 + diff_pixel[:,1]**2 +diff_pixel[:,2]**2)
            mask = distance_pixel < self.radius[i]
            length = len(distance_pixel[mask])
            crossing_param[i] = (length * separation_pixels)/self.radius[i]
        self.crossing_param = crossing_param


    def compute_los_distance(self,pixel_name):
        if(self.los_distance is not None):
            return()
        pixel = Pixel(name=pixel_name)
        pixel.read()
        coord_pixels = pixel.pixel_array[:,0:3]
        los_distance = np.full(self.coord.shape[0],np.inf)
        for i in range(len(self.coord)):
            dist_to_pixel = coord_pixels - self.coord[i]
            los_distance[i] = np.min(np.sqrt(dist_to_pixel[:,0]**2 + dist_to_pixel[:,1]**2 + dist_to_pixel[:,2]**2))
        self.los_distance = los_distance



    def get_crossing_qso(self,qso_name):
        qso = QSOCatalog.init_from_fits(qso_name)
        coord_crossing_qso = []
        for i in range(len(self.radius)):
            diff_pixel = qso.coord - self.coord[i]
            distance_pixel = np.sqrt(diff_pixel[:,0]**2 + diff_pixel[:,1]**2)
            mask = distance_pixel < self.radius[i]
            coord_crossing_qso.append(qso.coord[mask])
        return(coord_crossing_qso)




    def compute_cross_corr_parameters(self):
        self.convert_to_sky()
        redshift = self.coord[:,2].copy()
        self.convert_to_cartesian()
        self.convert_to_absolute_coordinates()
        X,Y,Z = self.coord[:,0].copy(),self.coord[:,1].copy(),self.coord[:,2].copy()
        self.convert_to_normalized_coordinates()
        void_coords = np.transpose(np.stack([X,Y,Z,self.weights,redshift]))
        return(void_coords)


    def correct_coordinates(self,method,name_out,inv_g_function,pixel_name,**kwargs):
        if(method == "barycenter"):
            self.correct_coordinates_barycenter(inv_g_function,pixel_name,**kwargs)
        elif(method == "mindist"):
            self.correct_coordinates_mindist(inv_g_function,pixel_name,**kwargs)



    def correct_coordinates_barycenter(self,
                                       inv_g,
                                       pixel_name,
                                       r_max=None,
                                       weight=None,
                                       sigma=None,
                                       A=1,
                                       L_perp=None,
                                       iterate=False,
                                       random_catalog_iterate=None,
                                       pixel_repacked=None):

        pixel = Pixel(name=pixel_name)
        pixel.read()
        if(pixel_repacked is not None):
            (x,y,z) = pickle.load(open(pixel_repacked,"rb"))
        else:
            x,y,z = pixel.repack_by_los()
        if(weight is None):
            function_weight = lambda x : 1
        elif(weight == "gauss"):
            function_weight = lambda x : np.exp(-0.5*(x/sigma)**2)
        elif(weight == "exp"):
            function_weight = lambda x : np.exp(-x/L_perp)
        mini_los = np.array([np.min(z[i]) for i in range(len(z))])
        maxi_los = np.array([np.max(z[i]) for i in range(len(z))])
        for i in range(len(self.coord)):
            mask = (self.coord[i,2] > mini_los)&(self.coord[i,2] < maxi_los)
            if(r_max is not None):
                mask &=(np.sqrt((self.coord[i,0]-x)**2 + (self.coord[i,1]-y)**2)<r_max)
            s_perp_i_x = - ( x[mask] - self.coord[i,0] )
            s_perp_i_y = - ( y[mask] - self.coord[i,1] )
            s_perp_i = np.sqrt(s_perp_i_x**2  + s_perp_i_y **2 )
            r_perp_i = inv_g(s_perp_i)
            r_perp_i_x = s_perp_i_x * (r_perp_i/s_perp_i)
            r_perp_i_y = s_perp_i_y * (r_perp_i/s_perp_i)
            d_r_perp_x = A * np.sum(function_weight(s_perp_i)*(r_perp_i_x - s_perp_i_x)) / np.sum(function_weight(s_perp_i))
            d_r_perp_y = A * np.sum(function_weight(s_perp_i)*(r_perp_i_y - s_perp_i_y)) / np.sum(function_weight(s_perp_i))
            print(f"Void {i+1} over {len(self.coord)}:")
            print(f"    x displacement {d_r_perp_x}")
            print(f"    y displacement {d_r_perp_y}")
            self.coord[i,0] = self.coord[i,0] + d_r_perp_x
            self.coord[i,1] = self.coord[i,1] + d_r_perp_y

    def correct_coordinates_mindist(self,inv_g_min,pixel_name):
        pixel = Pixel(name=pixel_name)
        pixel.read()
        x,y,z = pixel.repack_by_los()
        mini_los = np.array([np.min(z[i]) for i in range(len(z))])
        maxi_los = np.array([np.max(z[i]) for i in range(len(z))])
        for i in range(len(self.coord)):
            mask = (self.coord[i,2] > mini_los)&(self.coord[i,2] < maxi_los)
            s_perp_i_x = x[mask] - self.coord[i,0]
            s_perp_i_y = y[mask] - self.coord[i,1]
            s_perp_i = np.sqrt(s_perp_i_x**2  + s_perp_i_y **2 )
            arg_min = np.argmin(s_perp_i)
            r_perp_i = inv_g_min(s_perp_i[arg_min])
            r_perp_i_x = s_perp_i_x[arg_min] * (r_perp_i/s_perp_i[arg_min])
            r_perp_i_y = s_perp_i_y[arg_min] * (r_perp_i/s_perp_i[arg_min])
            self.coord[i,0] = self.coord[i,0] + r_perp_i_x - s_perp_i_x
            self.coord[i,1] = self.coord[i,1] + r_perp_i_y - s_perp_i_y






    def get_delta_void(self,rmin,rmax,nr,nameout,name_map,map_property_file):
        if(os.path.isfile("{}_rmin{}_rmax{}_nr{}.fits".format(nameout,rmin,rmax,nr))):
            return()
        r_array = np.linspace(rmin,rmax,nr)
        rmask = rmax + 2.
        tomographic_map = TomographicMap.init_from_property_files(map_property_file,name=name_map)
        tomographic_map.read()
        pixels_per_mpc = tomographic_map.pixel_per_mpc
        coord_pixels = np.round(self.coord * pixels_per_mpc,0).astype(int)
        indice_mpc = np.transpose(np.indices(tomographic_map.map_array.shape),axes=(1,2,3,0))/pixels_per_mpc
        delta_array = []
        for i in range(len(self.coord)):
            index = coord_pixels[i]
            number_pixel_maximal_radius = [int(round((rmask*pixels_per_mpc)[0],0)),int(round((rmask*pixels_per_mpc)[1],0)),int(round((rmask*pixels_per_mpc)[2],0))]
            map_local = tomographic_map.map_array[max(index[0]-number_pixel_maximal_radius[0],0):min(tomographic_map.shape[0],index[0]+number_pixel_maximal_radius[0]),
                                                  max(index[1]-number_pixel_maximal_radius[1],0):min(tomographic_map.shape[1],index[1]+number_pixel_maximal_radius[1]),
                                                  max(index[2]-number_pixel_maximal_radius[2],0):min(tomographic_map.shape[2],index[2]+number_pixel_maximal_radius[2])]
            indice_local = indice_mpc[max(index[0]-number_pixel_maximal_radius[0],0):min(tomographic_map.shape[0],index[0]+number_pixel_maximal_radius[0]),
                                      max(index[1]-number_pixel_maximal_radius[1],0):min(tomographic_map.shape[1],index[1]+number_pixel_maximal_radius[1]),
                                      max(index[2]-number_pixel_maximal_radius[2],0):min(tomographic_map.shape[2],index[2]+number_pixel_maximal_radius[2])] - self.coord[i]
            distance_map_local = np.sqrt(indice_local[:,:,:,0]**2 + indice_local[:,:,:,1]**2 + indice_local[:,:,:,2]**2)
            mask = distance_map_local < rmask
            delta_list = map_local[mask]
            r_list = distance_map_local[mask]
            delta = interp1d(r_list,delta_list)
            delta_array.append(delta(r_array))
        delta_array = np.mean(delta_array,axis=0)
        h = fitsio.FITS("{}_rmin{}_rmax{}_nr{}.fits".format(nameout,rmin,rmax,nr),"rw",clobber=True)
        h.write(r_array,extname="R")
        h.write(delta_array,extname="DELTA")
        h.close()





    @staticmethod
    def assessement_matrix_huge_voids(catalog_name,rmin,d_position,d_radius):
        catalog = [VoidCatalog.init_from_fits(name) for name in catalog_name]
        radius_list = [cat.radius for cat in catalog]
        coord_list = [cat.coord for cat in catalog]
        assessement_matrix = np.zeros((len(radius_list),len(radius_list)))
        for i in range(len(radius_list)):
            mask = radius_list[i] > rmin
            ni = len(radius_list[i][mask])
            r_big_i = radius_list[i][mask]
            coord_big_i  = coord_list[i][mask]
            for j in range(len(radius_list)):
                mask = radius_list[j] > rmin
                nj = len(radius_list[j][mask])
                r_big_j =  radius_list[j][mask]
                coord_big_j = coord_list[j][mask]
                nij=0
                for k in range(len(r_big_i)):
                    dist = np.sqrt((coord_big_j[:,0] - coord_big_i[k,0])**2 + (coord_big_j[:,1] - coord_big_i[k,1])**2  + (coord_big_j[:,2] - coord_big_i[k,2])**2 )
                    mask2= (abs(r_big_j - r_big_i[k]) < d_radius)&(dist < d_position)
                    nij = nij + len(r_big_j[mask2])
                assessement_matrix[i,j] = nij /((ni+nj)/2)
        return(assessement_matrix)
