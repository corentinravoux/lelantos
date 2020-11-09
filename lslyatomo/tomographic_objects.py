#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

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
try:
    from picca import data
except:
    import lslyatomo.picca.data as data
    raise Warning("Picca might be updated, we suggest to install picca independently")

from lslyatomo import utils
from scipy.ndimage.filters import gaussian_filter
import multiprocessing as mp
from scipy.interpolate import interp1d

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
        size,shape,boundary_cartesian_coord,boundary_sky_coord,coordinate_transform = None,None,None,None,None
        Property.read()
        if(Property.size is not None): size = Property.size
        if(Property.shape is not None): shape = Property.shape
        if(Property.boundary_cartesian_coord is not None): boundary_cartesian_coord = Property.boundary_cartesian_coord
        if(Property.boundary_sky_coord is not None): boundary_sky_coord = Property.boundary_sky_coord
        if(Property.coordinate_transform is not None): coordinate_transform = Property.coordinate_transform
        if(Property.Omega_m is not None): Omega_m = Property.Omega_m
        return(cls(map_array=map_array,name=name,shape=shape,size=size,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,property_file=property_file,Omega_m=Omega_m))


    @classmethod
    def init_by_merging(cls,launching_file_name,name_map,property_file):
        """ Need to generalized (maybe in the task manager class for reading and writting launching files) + naming filename must not be define there"""
        a = pickle.load(open(os.path.join(launching_file_name),"rb"))
        listname,Dachshundparams,numberChunks,overlaping = a[0],a[1],a[2],a[3]
        mapChunks = {}
        sizemap = {}
        for i in range(len(listname)):
            name = os.path.join("Python_Results",Dachshundparams[i]["namemap"])
            size = (Dachshundparams[i]["lx"],Dachshundparams[i]["ly"],Dachshundparams[i]["lz"])
            shape = (Dachshundparams[i]["nx"],Dachshundparams[i]["ny"],Dachshundparams[i]["nz"])
            submap = TomographicMap(name=name,shape=shape,size=size)
            submap.read()
            mapChunks[listname[i]] = submap.map_array
            sizemap[listname[i]] = (Dachshundparams[i]["lx"],Dachshundparams[i]["ly"])
        lx = 0
        ly = 0
        if overlaping !=0 :
            for i in range(numberChunks[0]):
                for j in range(numberChunks[1]):
                    filename = f'{i:03d}' + f'{j:03d}'
                    Map = mapChunks[filename]
                    nbPixelperMpcx = len(Map)/sizemap[filename][0]
                    nbPixelperMpcy = len(Map[0])/sizemap[filename][1]
                    nbpixelToremovex = int(round(nbPixelperMpcx * overlaping,0))
                    nbpixelToremovey = int(round(nbPixelperMpcy * overlaping,0))
                    if ((i==0)&(i == numberChunks[0] - 1)):
                        Map = Map
                    elif(i==0):
                        Map = Map[0:len(Map)-nbpixelToremovex,:,:]
                    elif(i == numberChunks[0] - 1):
                        Map = Map[nbpixelToremovex:len(Map),:,:]
                    else :
                        Map = Map[nbpixelToremovex:len(Map)-nbpixelToremovex,:,:]
                    if ((j==0)&(j == numberChunks[1] - 1)):
                        Map = Map
                    elif(j==0) :
                        Map = Map[:,0:len(Map[0])-nbpixelToremovey,:]
                    elif(j == numberChunks[1] - 1):
                        Map = Map[:,nbpixelToremovey:len(Map[0]),:]
                    else :
                        Map = Map[:,nbpixelToremovey:len(Map[0])-nbpixelToremovey,:]
                    mapChunks[filename]=Map
            for i in range(numberChunks[0]):
                filename = f'{i:03d}' + f'{0:03d}'
                if ((i==0)&(i == numberChunks[0] - 1)):
                    lx = lx + sizemap[filename][0]
                elif(i==0):
                    lx = lx + sizemap[filename][0] - overlaping
                elif(i == numberChunks[0] - 1):
                    lx = lx + sizemap[filename][0] - overlaping
                else :
                    lx = lx + sizemap[filename][0] - 2*overlaping
            for j in range(numberChunks[1]):
                filename = f'{0:03d}' + f'{j:03d}'
                if ((j==0)&(j == numberChunks[1] - 1)):
                    ly = ly + sizemap[filename][1]
                elif(j==0) :
                    ly = ly + sizemap[filename][1] - overlaping
                elif(j == numberChunks[1] - 1):
                    ly = ly + sizemap[filename][1] - overlaping
                else :
                    ly = ly + sizemap[filename][1] - 2 * overlaping

        else :
            for i in range(numberChunks[0]):
                filename = f'{i:03d}' + f'{0:03d}'
                lx = lx + sizemap[filename][0]
            for j in range(numberChunks[1]):
                filename = f'{0:03d}' + f'{j:03d}'
                ly = ly + sizemap[filename][1]
        concatenateList = []
        for i in range(numberChunks[0]):
            concatenate = np.concatenate([np.asarray(mapChunks[f'{i:03d}' + f'{j:03d}']) for j in range(numberChunks[1])],axis = 1)
            concatenateList.append(concatenate)
        merged_map_array = np.concatenate([concatenateList[i] for i in range(numberChunks[0])],axis = 0)
        map_shape = merged_map_array.shape
        map_size = (lx,ly,Dachshundparams[0]["lz"])
        prop = MapPixelProperty(name=property_file)
        prop.read()
        prop.shape = map_shape
        prop.size = map_size
        prop.write()
        merged_map = cls.init_from_property_files(property_file,map_array=merged_map_array,name=name_map)
        return(merged_map)




    @property
    def mpc_per_pixel(self):
        return(np.array(self.size)/(np.array(self.shape)))
#        return(np.array(self.size)/(np.array(self.shape)-1))

    @property
    def pixel_per_mpc(self):
        return((np.array(self.shape))/np.array(self.size))
        # return((np.array(self.shape)-1)/np.array(self.size))


    def read(self):
        if(self.name is None):
            raise ValueError("No")
        if(self.shape is None):
            raise ValueError("No")
        else:
            with open(self.name,'r') as f:
                self.map_array = np.fromfile(f,dtype=np.float64).reshape(self.shape)



    def write(self):
        if(self.map_array.all() == None):
            raise ValueError("No")
        listmap=np.ravel(self.map_array)
        listmap.tofile(self.name)

    def write_property_file(self,property_file_name):
        property_file = MapPixelProperty(name=property_file_name,size=self.size,shape=self.shape,boundary_cartesian_coord=self.boundary_cartesian_coord,boundary_sky_coord=self.boundary_sky_coord,coordinate_transform=self.coordinate_transform,Omega_m=self.Omega_m)
        property_file.write()

    def rebin_map(self,new_shape, operation='mean'):
        self.map_array = utils.bin_ndarray(self.map_array, new_shape, operation=operation)
        self.shape = new_shape


    def smooth(self,sigma_mpc):
        sigma = sigma_mpc / self.pixel_per_mpc
        self.map_array = gaussian_filter(self.map_array,sigma)



    def mask_map_to_3D(self,mapin,distance_name,distance):
        mapout = self.map_array.copy()
        distance_map = DistanceMap.init_from_tomographic_map(self,name=distance_name)
        distance_map.read()
        mask = distance_map.get_mask_distance(distance)
        mapout[mask] = 0
        return(mapout)


    def mask_map_from_name(self,distance_map_name,distance):
        distance_map = DistanceMap.init_from_tomographic_map(self,name=distance_map_name)
        mask = distance_map.get_mask_distance(distance)
        self.map_array = np.ma.masked_where(mask,self.map_array)

    def mask_map_from_mask(self,mask):
        self.map_array = np.ma.masked_where(mask,self.map_array)




    #### To DO

    def displace_map(self,size_map,max_list,max_list_to_center,size_map_to_center,shape_map_to_center,name_out_map,dist_map=None,name_dist_map_out=None,pixels=None,pixels_out=None,qso=None,qso_out=None):
        if(self.shapeMap != shape_map_to_center):
            raise KeyError("Not implemented already for different map shape")
        map_3d = self.readClamatoMapFile()
        new_map = np.zeros(shape_map_to_center)
        if(dist_map is not None):
            dist_map_3d = self.readClamatoMapFile_Other(dist_map,self.shapeMap)
            new_dist_map = np.full(shape_map_to_center,np.inf)
        indice = np.transpose(np.indices(shape_map_to_center),axes=(1,2,3,0))
        nb_mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap)-1)
        nb_mpc_per_pixels_to_center = np.array(size_map_to_center)/(np.array(shape_map_to_center)-1)
        coord_mpc = indice * nb_mpc_per_pixels_to_center
        del indice
        (minx,maxx,miny,maxy,minz,maxz,minredshift,maxredshift,minra,maxra,mindec,maxdec) = pickle.load(open(max_list,"rb"))
        (minx2,maxx2,miny2,maxy2,minz2,maxz2,minredshift2,maxredshift2,minra2,maxra2,mindec2,maxdec2) = pickle.load(open(max_list_to_center,"rb"))
        coord_centered = np.zeros(coord_mpc.shape)
        coord_centered[:,:,:,0],coord_centered[:,:,:,1],coord_centered[:,:,:,2] = coord_mpc[:,:,:,0] - (minx - minx2),coord_mpc[:,:,:,1] -( miny - miny2),coord_mpc[:,:,:,2] - (minz - minz2)
        if(pixels is not None):
            qso = np.transpose(np.array(pickle.load(open("DataQSOposition.pickle","rb"))))
            pixel_file = self.readClamatoPixelFile_other(pixels)
            pixel_file[:,0],pixel_file[:,1],pixel_file[:,2] = pixel_file[:,0] + (minx - minx2) ,pixel_file[:,1]  + (miny - miny2),pixel_file[:,2] + (minz - minz2)
            qso[:,0],qso[:,1],qso[:,2] = qso[:,0] + (minx - minx2) ,qso[:,1]  + (miny - miny2),qso[:,2] + (minz - minz2)
            mask = (pixel_file[:,0]>=0)&(pixel_file[:,0]<size_map_to_center[0])
            mask &= (pixel_file[:,1]>=0)&(pixel_file[:,1]<size_map_to_center[1])
            mask &= (pixel_file[:,2]>=0)&(pixel_file[:,2]<size_map_to_center[2])
            self.writeClamatoPixelFile(pixels_out,pixel_file[mask])
            mask = (qso[:,0]>=0)&(qso[:,0]<size_map_to_center[0])
            mask &= (qso[:,1]>=0)&(qso[:,1]<size_map_to_center[1])
            mask &= (qso[:,2]>=0)&(qso[:,2]<size_map_to_center[2])
            pickle.dump(np.transpose(qso),open("DataQSOposition.pickle","wb"))
        del coord_mpc
        indice_centered = np.round(coord_centered / nb_mpc_per_pixels,0).astype(int)
        mask =(indice_centered[:,:,:,0] >= 0)&(indice_centered[:,:,:,0] < shape_map_to_center[0])
        mask&=(indice_centered[:,:,:,1] >= 0)&(indice_centered[:,:,:,1] < shape_map_to_center[1])
        mask&=(indice_centered[:,:,:,2] >= 0)&(indice_centered[:,:,:,2] < shape_map_to_center[2])
        new_map[mask] = map_3d[indice_centered[:,:,:,0][mask],indice_centered[:,:,:,1][mask],indice_centered[:,:,:,2][mask]]
        if(dist_map is not None):
            new_dist_map[mask] = dist_map_3d[indice_centered[:,:,:,0][mask],indice_centered[:,:,:,1][mask],indice_centered[:,:,:,2][mask]]
        del indice_centered,mask
        self.writeClamatoMapFile(name_out_map,new_map)
        self.writeClamatoMapFile(name_dist_map_out,new_dist_map)





class DistanceMap(TomographicMap):

    def __init__(self,map_array=None,name=None,shape=None,size=None,boundary_cartesian_coord=None,boundary_sky_coord=None,coordinate_transform=None,property_file=None,Omega_m=None):
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
        mpc_per_pixels = tomographic_map.mpc_per_pixel
        distance_array = np.full(tomographic_map.shape,np.inf)
        indice = np.transpose(np.indices(tomographic_map.shape),axes=(1,2,3,0))
        indice_mpc = indice * mpc_per_pixels
        del indice
        for i in range(len(pixels)):
            minimal_coordinates_local = [int(round(((pixels[i][0][0]-radius_local)/mpc_per_pixels)[0],0)),int(round(((pixels[i][0][1]-radius_local)/mpc_per_pixels)[1],0))]
            maximal_coordinates_local = [int(round(((pixels[i][0][0]+radius_local)/mpc_per_pixels)[0],0)),int(round(((pixels[i][0][1]+radius_local)/mpc_per_pixels)[1],0))]
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
    def create_distance_map_parallel(cls,pixel,tomographic_map,name_dist_los,nb_process):
        cls.log.add("Creation of global variables")
        x,y,z = pixel.repack_by_los()
        pixels = [[[x[i],y[i]],z[i]] for i in range(len(x))]
        mpc_per_pixels = tomographic_map.mpc_per_pixel
        indice = np.transpose(np.indices(tomographic_map.shape),axes=(1,2,3,0))
        global dist_map
        dist_map = indice * mpc_per_pixels
        del indice
        cls.log.add("Number LOS to treat : {}".format(len(pixels)))
        cls.log.add("Launching of Pool with shared array initialization")
        DistanceMap.init_shared_array(tomographic_map.shape)
        pool = mp.Pool(nb_process)
        pool.map(cls.worker_create_distance_map, pixels,tomographic_map.shape)
        cls.log.add("End of Pool")
        cls.log.add("Getting the map from shared array")
        distance_array = DistanceMap.mp_array_to_numpyarray(shared_arr).reshape(tomographic_map.shape)
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
            LOS_range_min = DistanceMap.mp_array_to_numpyarray(shared_arr).reshape(shape)
            LOS_range_min[:,:,:] = np.minimum(LOS_range,LOS_range_min)
            del LOS_range_min
        cls.log.add("Lock released for a worker")
        cls.log.add("LOS treated")
        del LOS_range


    @staticmethod
    def init_shared_array(shape):
        global shared_arr
        distance_array = np.full(shape,np.inf)
        shared_arr = mp.Array('d', distance_array.flatten())
        del distance_array

    @staticmethod
    def mp_array_to_numpyarray(mp_arr):
        return np.frombuffer(mp_arr.get_obj())

        ### TO DO
    def give_redshift_cut(self,mask_los_distance,compute_redshift_out_LOS,sizeMap,redshift_axis=None):
        i,percent = 0,0.
        while((percent<compute_redshift_out_LOS)&(i< mask_los_distance.shape[-1])):
            nb_True = len(mask_los_distance[:,:,i][mask_los_distance[:,:,i] == True])
            percent = nb_True/(mask_los_distance.shape[0]*mask_los_distance.shape[1])
            i=i+1
        dist_mpc = i * (sizeMap[-1]/self.shapeMap[-1])
        if(redshift_axis is not None):
            Om,maxlist_name = redshift_axis
            maxlist = pickle.load(open(maxlist_name,'rb'))
            minz,minredshift = maxlist[4],maxlist[6]
            Cosmo = constants.cosmo(Om)
            self.rcomoving = Cosmo.r_comoving
            redshift = fsolve(self.f,minredshift,args=(dist_mpc + minz))[0]
            print(redshift)



    def get_mask_distance(self,distance):
        return(self.map_array > distance)




class StackMap(TomographicMap):

    def __init__(self,tomo_map=None,catalog=None,map_array=None,name=None,shape=None,size=None,boundary_cartesian_coord=None,boundary_sky_coord=None,coordinate_transform=None,property_file=None,Omega_m=None):
        super(StackMap,self).__init__(map_array=map_array,name=name,shape=shape,size=size,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,property_file=property_file,Omega_m=Omega_m)

        self.tomo_map = tomo_map
        self.catalog = catalog

    @classmethod
    def init_by_tomographic_map(cls,catalog_name,type_catalog,tomographic_map_name,property_file_map,size_stack,name=None):
        tomographic_map = TomographicMap.init_from_property_files(property_file_map,name=tomographic_map_name)
        tomographic_map.read()
        catalog = Catalog.init_catalog_from_fits(catalog_name, type_catalog)
        shape_stack = (np.round(size_stack/tomographic_map.mpc_per_pixel,0)).astype(int)
        index_catalog = (np.round(catalog.coord/tomographic_map.mpc_per_pixel,0)).astype(int)

        mask = (index_catalog[:,0] < shape_stack[0])|(index_catalog[:,0] >= tomographic_map.shape[0] - shape_stack[0])
        mask |=(index_catalog[:,1] < shape_stack[1])|(index_catalog[:,1] >= tomographic_map.shape[1] - shape_stack[1])
        mask |=(index_catalog[:,2] < shape_stack[2])|(index_catalog[:,2] >= tomographic_map.shape[2] - shape_stack[2])

        index_catalog_clean = index_catalog[~mask]
        local_maps = np.zeros((len(index_catalog_clean),2*shape_stack[0]+1,2*shape_stack[1]+1,2*shape_stack[2]+1))
        for i in range(len(index_catalog_clean)):
            local_maps[i] = tomographic_map.map_array[index_catalog_clean[i][0]-shape_stack[0]:index_catalog_clean[i][0]+shape_stack[0]+1,index_catalog_clean[i][1]-shape_stack[1]:index_catalog_clean[i][1]+shape_stack[1]+1,index_catalog_clean[i][2]-shape_stack[2]:index_catalog_clean[i][2]+shape_stack[2]+1]
        stack = np.mean(local_maps,axis=0)


        boundary_cartesian_coord = None
        boundary_sky_coord = None

        return(cls(tomo_map=tomographic_map,catalog=catalog,map_array=stack,name=name,shape=shape_stack,size=size_stack,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=tomographic_map.coordinate_transform,property_file=None))


    @classmethod
    def init_stack_by_property_file(cls,property_file_name,name=None,tomographic_map_name=None,property_file_map=None,catalog_name=None,type_catalog=None):
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



    def stack_voids(self,map_3D,size_map,coord,radius,size_stack,normalized=None,shape="CUBIC",number=None):
        if(normalized is None):
            voids = coord
# shape/size
            number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap)-1)
            number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap))
            shape_stack = (np.round(size_stack/number_Mpc_per_pixels,0)).astype(int)
            index_voids = np.zeros(len(voids))
            index_voids = (np.round(voids/number_Mpc_per_pixels,0)).astype(int)
            mask = (index_voids[:,0] < shape_stack[0])|(index_voids[:,0] >= self.shapeMap[0] - shape_stack[0])
            mask |=(index_voids[:,1] < shape_stack[1])|(index_voids[:,1] >= self.shapeMap[1] - shape_stack[1])
            mask |=(index_voids[:,2] < shape_stack[2])|(index_voids[:,2] >= self.shapeMap[2] - shape_stack[2])
            clean_voids = index_voids[~mask]
            if(number is not None):
                if(number<len(clean_voids)):
                    nb_void = number
                else :
                    nb_void = len(clean_voids)
            else :
                nb_void = len(clean_voids)
            local_maps = np.zeros((nb_void,2*shape_stack[0]+1,2*shape_stack[1]+1,2*shape_stack[2]+1))
            for i in range(nb_void):
                local_maps[i] = map_3D[clean_voids[i][0]-shape_stack[0]:clean_voids[i][0]+shape_stack[0]+1,clean_voids[i][1]-shape_stack[1]:clean_voids[i][1]+shape_stack[1]+1,clean_voids[i][2]-shape_stack[2]:clean_voids[i][2]+shape_stack[2]+1]
            stack = np.mean(local_maps,axis=0)
            return(stack)
        else :
            voids = coord
            radius = radius
            mask_radius = radius >= normalized
            voids = voids[mask_radius]
            radius = radius[mask_radius]
            minimal_radius = np.min(radius)
            size = (size_stack/minimal_radius) * radius
# shape/size
            number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap)-1)
            number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap))
            mask = ((np.round((voids[:,0]-size[:])/number_Mpc_per_pixels[0],0)).astype(int) < 0)|((np.round((voids[:,0] + size[:])/number_Mpc_per_pixels[0],0)).astype(int) >= self.shapeMap[0])
            mask |=((np.round((voids[:,1]-size[:])/number_Mpc_per_pixels[1],0)).astype(int) < 0)|((np.round((voids[:,1] + size[:])/number_Mpc_per_pixels[1],0)).astype(int) >= self.shapeMap[1])
            mask |=((np.round((voids[:,2]-size[:])/number_Mpc_per_pixels[2],0)).astype(int) < 0)|((np.round((voids[:,2] + size[:])/number_Mpc_per_pixels[2],0)).astype(int) >= self.shapeMap[2])
            clean_voids = voids[~mask]
            size = size[~mask]
            shape_stack_min = (np.round(size_stack/number_Mpc_per_pixels,0)).astype(int)
            if(number is not None):
                if(number<len(clean_voids)):
                    nb_void = number
                else :
                    nb_void = len(clean_voids)
            else :
                nb_void = len(clean_voids)
            local_maps = np.zeros((nb_void,2*shape_stack_min[0]+1,2*shape_stack_min[1]+1,2*shape_stack_min[2]+1))
            for i in range(nb_void):
                r_overrmin = (size[i] /size_stack)
                for j in range(2*shape_stack_min[0]+1):
                    for k in range(2*shape_stack_min[1]+1):
                        for l in range(2*shape_stack_min[2]+1):
                            local_maps[i,j,k,l] = map_3D[int(round(clean_voids[i][0]/number_Mpc_per_pixels[0] + r_overrmin * (j-shape_stack_min[0]),0)),int(round(clean_voids[i][1]/number_Mpc_per_pixels[1] + r_overrmin * (k-shape_stack_min[1]),0)),int(round(clean_voids[i][2]/number_Mpc_per_pixels[2] + r_overrmin * (l-shape_stack_min[2]),0))]
            stack = np.mean(local_maps,axis=0)
            return(stack)




    def stack_voids_correction_xyztildestack(self,map_3D,size_map,coord,size_stack,shape="CUBIC"):
        voids = np.array(coord) + np.array([self.minx,self.miny,self.minz])
# shape/size
        number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap)-1)
        number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap))
        shape_stack = (np.round(size_stack/number_Mpc_per_pixels,0)).astype(int)
        x_tilde_array_normalized = np.linspace(-size_stack,size_stack,2*shape_stack[0]+1)
        y_tilde_array_normalized = np.linspace(-size_stack,size_stack,2*shape_stack[1]+1)
        z_tilde_array_normalized = np.linspace(-size_stack,size_stack,2*shape_stack[2]+1)
        coords_tilde_normalized =np.moveaxis(np.array(np.meshgrid(x_tilde_array_normalized,y_tilde_array_normalized,z_tilde_array_normalized,indexing='ij')),0,-1)
        del x_tilde_array_normalized, y_tilde_array_normalized
        clean_voids, local_maps = [],[]
        for i in range(len(voids)):
            center_tilde = list(self.xyz_to_xyztilde(voids[i][0],voids[i][1],voids[i][2]))
            coords_tilde = coords_tilde_normalized + center_tilde
            z_tilde_array = z_tilde_array_normalized + center_tilde[2]
            redshift_array = self.z_to_redshift(z_tilde_array)
            del z_tilde_array
            redshift_3Darray = np.zeros(coords_tilde.shape[0:-1])
            for j in range(len(redshift_3Darray)):
                for k in range(len(redshift_3Darray[j])):
                    redshift_3Darray[j,k,:] = redshift_array
            del redshift_array
            coords = np.zeros(coords_tilde.shape)
            coords[:,:,:,0],coords[:,:,:,1],coords[:,:,:,2] = self.xyztilde_to_xyz(coords_tilde[:,:,:,0],coords_tilde[:,:,:,1],coords_tilde[:,:,:,2],redshift=redshift_3Darray)
            del coords_tilde
            coords_box = coords - np.array([self.minx,self.miny,self.minz])
            del coords
            coords_pixels  = np.round(coords_box/number_Mpc_per_pixels,0).astype(int)
            del coords_box
            boolean = ((np.max(coords_pixels[:,:,:,0])>=self.shapeMap[0])|(np.min(coords_pixels[:,:,:,0])<0))
            boolean |=((np.max(coords_pixels[:,:,:,1])>=self.shapeMap[1])|(np.min(coords_pixels[:,:,:,1])<0))
            boolean |=((np.max(coords_pixels[:,:,:,2])>=self.shapeMap[2])|(np.min(coords_pixels[:,:,:,2])<0))
            if(not(boolean)):
                clean_voids.append(voids[i])
                local_maps.append(map_3D[coords_pixels[:,:,:,0],coords_pixels[:,:,:,1],coords_pixels[:,:,:,2]])
            del coords_pixels
        del coords_tilde_normalized
        local_maps = np.array(local_maps)
        stack = np.mean(local_maps,axis=0)
        del local_maps
        return(stack)


    def stack_voids_correction_xyzstack(self,map_3D,size_map,coord,size_stack,shape="CUBIC"):
        voids = np.array(coord) + np.array([self.minx,self.miny,self.minz])
# shape/size
        number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap)-1)
        number_Mpc_per_pixels = np.array(size_map)/(np.array(self.shapeMap))
        shape_stack = (np.round(size_stack/number_Mpc_per_pixels,0)).astype(int)
        index_voids = np.zeros(len(voids))
        index_voids = (np.round(voids/number_Mpc_per_pixels,0)).astype(int)
        mask = (index_voids[:,0] < shape_stack[0])|(index_voids[:,0] >= self.shapeMap[0] - shape_stack[0])
        mask |=(index_voids[:,1] < shape_stack[1])|(index_voids[:,1] >= self.shapeMap[1] - shape_stack[1])
        mask |=(index_voids[:,2] < shape_stack[2])|(index_voids[:,2] >= self.shapeMap[2] - shape_stack[2])
        clean_voids = index_voids[~mask]
        nb_void = len(clean_voids)
        x_array_normalized = np.linspace(-size_stack,size_stack,2*shape_stack[0]+1)
        y_array_normalized = np.linspace(-size_stack,size_stack,2*shape_stack[1]+1)
        z_array_normalized = np.linspace(-size_stack,size_stack,2*shape_stack[2]+1)
        coords_normalized =np.moveaxis(np.array(np.meshgrid(x_array_normalized,y_array_normalized,z_array_normalized,indexing='ij')),0,-1)
        clean_voids, local_maps = [],[]
        for i in range(len(voids)):
            center_tilde = np.array(self.xyz_to_xyztilde(voids[i][0],voids[i][1],voids[i][2]))
            coords = coords_normalized + voids[i]
            z_array = z_array_normalized + voids[i][2]
            redshift_array = self.z_to_redshift(z_array)
            del z_array
            redshift_3Darray = np.zeros(coords.shape[0:-1])
            for j in range(len(redshift_3Darray)):
                for k in range(len(redshift_3Darray[j])):
                    redshift_3Darray[j,k,:] = redshift_array
            del redshift_array
            coords_tilde = np.zeros(coords.shape)
            coords_tilde[:,:,:,0],coords_tilde[:,:,:,1],coords_tilde[:,:,:,2] = self.xyz_to_xyztilde(coords[:,:,:,0],coords[:,:,:,1],coords[:,:,:,2],redshift=redshift_3Darray)
            coords_box_tilde = coords_tilde - np.array([self.minx,self.miny,self.minz]) + ( voids[i] - center_tilde)
            coords_pixels_tilde  = np.round(coords_box_tilde/number_Mpc_per_pixels,0).astype(int)
            boolean = ((np.max(coords_pixels_tilde[:,:,:,0])>=self.shapeMap[0])|(np.min(coords_pixels_tilde[:,:,:,0])<0))
            boolean |=((np.max(coords_pixels_tilde[:,:,:,1])>=self.shapeMap[1])|(np.min(coords_pixels_tilde[:,:,:,1])<0))
            boolean |=((np.max(coords_pixels_tilde[:,:,:,2])>=self.shapeMap[2])|(np.min(coords_pixels_tilde[:,:,:,2])<0))
            if(not(boolean)):
                clean_voids.append(voids[i])
                local_maps.append(map_3D[coords_pixels_tilde[:,:,:,0],coords_pixels_tilde[:,:,:,1],coords_pixels_tilde[:,:,:,2]])
        for i in range(nb_void):
            local_maps[i] = map_3D[clean_voids[i][0]-shape_stack[0]:clean_voids[i][0]+shape_stack[0]+1,clean_voids[i][1]-shape_stack[1]:clean_voids[i][1]+shape_stack[1]+1,clean_voids[i][2]-shape_stack[2]:clean_voids[i][2]+shape_stack[2]+1]
        stack = np.mean(local_maps,axis=0)
        return(stack)


    def initialize_coordinates_conversion(self,Om,maxlist_name):
        Cosmo = constants.cosmo(Om)
        maxlist = pickle.load(open(maxlist_name,'rb'))
        (minx,maxx,miny,maxy,minz,maxz,minredshift,maxredshift,minra,maxra,mindec,maxdec) = maxlist
        self.minx = minx
        self.maxx = maxx
        self.miny = miny
        self.maxy = maxy
        self.minz = minz
        self.maxz = maxz
        self.minra = minra
        self.maxra = maxra
        self.mindec = mindec
        self.maxdec = maxdec
        self.minredshift = minredshift
        self.meanredshift = (minredshift + maxredshift)/2
        self.rcomoving = Cosmo.r_comoving
        self.redshift_to_dm = Cosmo.dm



    def xyz_to_xyztilde(self,x,y,z,redshift=None):
        z_tilde = z
        if(redshift is None):Z = self.z_to_redshift(z)
        else: Z = redshift
        x_tilde = (self.redshift_to_dm(Z)/self.redshift_to_dm(self.meanredshift)) * (x)
        y_tilde = (self.redshift_to_dm(Z)/self.redshift_to_dm(self.meanredshift)) * (y)
        return(x_tilde,y_tilde,z_tilde)

    def xyztilde_to_xyz(self,x_tilde,y_tilde,z_tilde,redshift=None):
        z = z_tilde
        if(redshift is None):Z = self.z_to_redshift(z)
        else: Z = redshift
        x = (self.redshift_to_dm(self.meanredshift)/self.redshift_to_dm(Z)) * x_tilde
        y = (self.redshift_to_dm(self.meanredshift)/self.redshift_to_dm(Z)) * y_tilde
        return(x,y,z)


    def rcomoving_inverse(self,redshift,z):
        return(self.rcomoving(redshift) - (z))


    def z_to_redshift_scalar(self,z):
        redshift = fsolve(self.rcomoving_inverse,self.minredshift,args=(z))[0]
        return(redshift)

    def z_to_redshift_array3D(self,z):
        redshift = np.zeros(z.shape)
        for i in range(len(redshift)):
            for j in range(len(redshift[i])):
                for k in range(len(redshift[i][j])):
                    redshift[i][j][k] = self.z_to_redshift_scalar(z[i][j][k])
        return(redshift)

    def z_to_redshift_array1D(self,z):
        redshift = np.zeros(z.shape)
        for i in range(len(redshift)):
            redshift[i] = self.z_to_redshift_scalar(z[i])
        return(redshift)

    def z_to_redshift(self,z):
        if((type(z)==np.ndarray)&(z.ndim==3)):
            redshift = self.z_to_redshift_array3D(z)
        elif((type(z)==np.ndarray)&(z.ndim==1)):
            redshift = self.z_to_redshift_array1D(z)
        else:
            redshift = self.z_to_redshift_scalar(z)
        return(redshift)



    def read_stacks(self,name,shape_stack):
        stack = np.fromfile(name)
        return(stack.reshape(shape_stack))

    def merge_stacks(self,list_stacks,shape_stack,name):
        stacks = []
        for i in range(len(list_stacks)):
            stack = self.read_stacks(list_stacks[i],shape_stack)
            stacks.append(stack)
        mean_stack = np.mean(stacks,axis=0)
        self.save_a_stack(mean_stack,name)
        return(mean_stack)

    def save_a_stack(self,stack,name,los_z = None):
        if(los_z is not None):
            np.savetxt("lenght_between_los_and_quasar.txt",[los_z],header="difference in Mpc between end of LOS and quasar center\n")
        stack.tofile(name)





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
            raise ValueError("No")
        else:
            with open(self.name,'r') as f:
                pixel_data = np.fromfile(f,dtype=np.float64)
                self.pixel_array = pixel_data.reshape((len(pixel_data)//5,5))


    def write(self):
        if(self.pixel_array is None):
            raise ValueError("No")
        listpixel=np.ravel(self.pixel_array)
        listpixel.tofile(self.name)


    def repack_by_los(self):
        coord = self.pixel_array
        unique_coord = np.unique(coord[:,0:2],axis=0)
        z = []
        for i in range(unique_coord.shape[0]):
            z.append(coord[:,2][(coord[:,0]==unique_coord[i,0])&(coord[:,1]==unique_coord[i,1])])
        return(unique_coord[:,0],unique_coord[:,1],np.asarray(z))

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
        for i in range(len(zpar)):
            dmin,N_los = Pixel.compute_mean_distance_density_at_z(zpar[i],np.array(x),np.array(y),mini_los,maxi_los)
            dperpz.append(dmin)
            densityz.append(N_los/((maxra-minra)*(maxdec-mindec)*((180/np.pi)**2)))
        return(zpar,dperpz,densityz)


    # CR - to debug
    def compute_mean_distance_density(self):
        if((self.z_array is None)|(self.dperp_array is None)|(self.density_array is None)):
            x,y,z = self.repack_by_los()
            minra,mindec = self.boundary_sky_coord[0][0:2]
            maxra,maxdec = self.boundary_sky_coord[1][0:2]
            minz = 0.0 #self.boundary_cartesian_coord[0][2]
            maxz = 0.0 #self.boundary_cartesian_coord[1][2]
            self.z_array,self.dperp_array,self.density_array = Pixel.return_mean_distance_density(x,y,z,minra,maxra,mindec,maxdec,minz,maxz)
        return(self.z_array,self.dperp_array,self.density_array)


    def compute_mean_distance_histogram(self,z_value):
        x,y,z = self.repack_by_los()
        mini_los = np.array([np.min(z[i]) for i in range(len(z))])
        maxi_los = np.array([np.max(z[i]) for i in range(len(z))])
        mask = (z_value > mini_los)&(z_value < maxi_los)
        xintheplane = x[mask]
        yintheplane = y[mask]
        return(Pixel.compute_mean_distance(xintheplane,yintheplane,return_array=True))





   ### To modify in merge pixels (by an initialization) ###

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




#############################################################################
#############################################################################
############################### DELTA #####################################
#############################################################################
#############################################################################




class Delta(object):

    def __init__(self,delta_file=None,delta_array=None,name=None,pk1d_type=True):
        self.name = name
        self.delta_file = delta_file
        self.delta_array = delta_array
        self.pk1d_type = pk1d_type


    def read(self):
        if(self.name == None):
            raise ValueError("No")
        else:
            self.delta_file = fitsio.FITS(self.name)

    def read_line(self,number_line):
        if(self.delta_file is None):
            self.read_from_fits()
        delta = data.delta.from_fitsio(self.delta_file[number_line],Pk1D_type=self.pk1d_type)
        return(delta)


    def create_fi(self,delta):
        nrows = len(delta.de)
        head = {}
        if  self.pk1d_type :
            h = np.zeros(nrows, dtype=[('LOGLAM','f8'),('DELTA','f8'),('IVAR','f8'),('DIFF','f8')])
            h['DELTA'] =delta.de
            h['LOGLAM'] = delta.ll
            h['IVAR'] = delta.iv
            h['DIFF'] = delta.diff
            head['MEANSNR'] = delta.mean_SNR
            head['MEANRESO'] = delta.mean_reso
            head['MEANZ'] = delta.mean_z
            head['DLL'] = delta.dll
        else :
            h = np.zeros(nrows, dtype=[('LOGLAM','f8'),('DELTA','f8'),('WEIGHT','f8'),('CONT','f8')])
            h['DELTA'] =delta.de
            h['LOGLAM'] = delta.ll
            h['WEIGHT'] = delta.we
            h['CONT'] = delta.co
        head['THING_ID'] = delta.thid
        head['RA'] = delta.ra
        head['DEC'] = delta.dec
        head['Z']  = delta.zqso
        head['PLATE'] = delta.plate
        head['MJD'] = delta.mjd
        head['FIBERID'] = delta.fid
        return(h,head)


    def write_from_delta_list(self):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        for i in range(len(self.delta_array)):
            delta =self.delta_array[i]
            fi,head = self.create_fi(delta)
            fits.write(fi,header=head)

    def close(self):
        self.delta_file.close()




#############################################################################
#############################################################################
############################### CATALOGS #####################################
#############################################################################
#############################################################################


# CR - for the cutting routines, add a routine which automaticaly cut additive arrays






class Catalog(object):


    def __init__(self,name=None,coord=None,primary_key=None,catalog_type="sky"):
        self.name = name
        self.coord = coord
        self.primary_key = primary_key
        self.catalog_type = catalog_type

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
            raise ValueError("No")
        else:
            return(fitsio.FITS(name))

    @staticmethod
    def close(catalog):
        catalog.close()


    def convert_to_absolute_coordinates(self,map_property_file):
        prop = MapPixelProperty(name=map_property_file)
        prop.read()
        if(self.catalog_type=="cartesian"):
            boundary = prop.boundary_cartesian_coord
        if(self.catalog_type=="sky"):
            boundary = prop.boundary_sky_coord
        self.coord[:,0] = self.coord[:,0]  + boundary[0][0]
        self.coord[:,1] = self.coord[:,1]  + boundary[0][1]
        self.coord[:,2] = self.coord[:,2]  + boundary[0][2]

    def convert_to_normalized_coordinates(self,map_property_file):
        prop = MapPixelProperty(name=map_property_file)
        prop.read()
        if(self.catalog_type=="cartesian"):
            boundary = prop.boundary_cartesian_coord
        if(self.catalog_type=="sky"):
            boundary = prop.boundary_sky_coord
        self.coord[:,0] = self.coord[:,0]  - boundary[0][0]
        self.coord[:,1] = self.coord[:,1]  - boundary[0][1]
        self.coord[:,2] = self.coord[:,2]  - boundary[0][2]

    def convert_to_sky(self,map_property_file,mode):
        prop = MapPixelProperty(name=map_property_file)
        prop.read()
        self.convert_to_absolute_coordinates(map_property_file)
        if((mode.lower()=="middle")|(mode.lower()=="full_angle")|(mode.lower()=="full")):
            Omega_m = prop.Omega_m
            suplementary_parameters = utils.return_suplementary_parameters(self.coordinate_transform,property=prop)
            (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(Omega_m)
            self.coord[:,0],self.coord[:,1],self.coord[:,2] = utils.convert_cartesian_to_sky(self.coord[:,0],self.coord[:,1],self.coord[:,2],mode,inv_rcomov=inv_rcomov,inv_distang=inv_distang,distang=distang,suplementary_parameters=suplementary_parameters)
            self.catalog_type = "sky"
        else:
            raise KeyError("Conversion mode not available, please choose between : middle, full or full_angle")

    def convert_to_cartesian(self,map_property_file,mode):
        prop = MapPixelProperty(name=map_property_file)
        prop.read()
        if((mode.lower()=="middle")|(mode.lower()=="full_angle")|(mode.lower()=="full")):
            Omega_m = prop.Omega_m
            suplementary_parameters = utils.return_suplementary_parameters(self.coordinate_transform,property=prop)
            (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(Omega_m)
            self.coord[:,0],self.coord[:,1],self.coord[:,2] = utils.convert_sky_to_cartesian(self.coord[:,0],self.coord[:,1],self.coord[:,2],mode,rcomov=rcomov,distang=distang,suplementary_parameters=suplementary_parameters)
            self.catalog_type = "cartesian"
        else:
            raise KeyError("Conversion mode not available, please choose between : middle, full or full_angle")
        self.convert_to_normalized_coordinates(map_property_file)




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




class QSOCatalog(Catalog):

    def __init__(self,name=None,coord=None,primary_key=None,plate=None,modern_julian_date=None,fiber_id=None,redshift_name="Z",catalog_type="sky"):
        super(QSOCatalog,self).__init__(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type)
        self.redshift_name = redshift_name
        self.plate = plate
        self.modern_julian_date = modern_julian_date
        self.fiber_id = fiber_id


    @classmethod
    def init_from_fits(cls,name,redshift_name="Z"):
        catalog = Catalog.load_from_fits(name)
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
            return(cls(name=name,coord=coord,primary_key=primary_key,plate=plate,modern_julian_date=modern_julian_date,fiber_id=fiber_id,redshift_name=redshift_name,catalog_type=catalog_type))
        if("X" in catalog[1].get_colnames()):
            coord_ra = catalog[1]["X"][:]
            coord_dec = catalog[1]["Y"][:]
            coord_z = catalog[1]['Z'][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            primary_key = catalog[1]["THING_ID"][:]
            catalog_type = "cartesian"
            Catalog.close(catalog)
            return(cls(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type))

    @classmethod
    def init_from_pixel_catalog(cls,quasar_pixels,name=None):
        if(name is None): name ="qso_catalog.fits"
        coord_ra = quasar_pixels[:,0]
        coord_dec = quasar_pixels[:,1]
        coord_z = quasar_pixels[:,2]
        coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
        primary_key = quasar_pixels[:,3]
        catalog_type = "cartesian"
        return(cls(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type))


    def cut_catalog_qso(self,coord_min=None,coord_max=None):
        mask_select = self.cut_catalog(coord_min=coord_min,coord_max=coord_max,center_x_coord=True)
        if(self.coord is not None):
            self.coord = self.coord[mask_select]
        if(self.primary_key is not None):
            self.primary_key = self.primary_key[mask_select]
        if(self.plate is not None):
            self.plate = self.plate[mask_select]
        if(self.modern_julian_date is not None):
            self.modern_julian_date = self.modern_julian_date[mask_select]
        if(self.fiber_id is not None):
            self.fiber_id = self.fiber_id[mask_select]


    def write(self):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        nrows = self.coord.shape[0]
        if(self.catalog_type == "sky"):
            h = np.zeros(nrows, dtype=[('RA','f8'),('DEC','f8'),(self.redshift_name,'f8'),('THING_ID','i8'),('PLATE','i4'),('MJD','i4'),('FIBERID','i2')])
            h['RA'] = self.coord[:,0]
            h['DEC'] = self.coord[:,1]
            h[self.redshift_name] = self.coord[:,2]
            h['THING_ID'] =self.primary_key
            h['PLATE'] =self.plate
            h['MJD'] = self.modern_julian_date
            h['FIBERID'] =self.fiber_id
        elif(self.catalog_type == "cartesian"):
            h = np.zeros(nrows, dtype=[('X','f8'),('Y','f8'),('Z','f8'),('THING_ID','i8')])
            h['X'] = self.coord[:,0]
            h['Y'] = self.coord[:,1]
            h['Z'] = self.coord[:,2]
            h['THING_ID'] =self.primary_key
        fits.write(h)
        fits.close()



    def cut_write_catalog(self,name_out,coord_min,coord_max):
        self.read_from_fits()
        self.cut_quasars_catalogs(coord_min,coord_max)
        self.name = name_out
        self.write()



class DLACatalog(Catalog):

    def __init__(self,name=None,coord=None,primary_key=None,z_qso=None,confidence=None,nhi=None,catalog_type="sky"):
        super(DLACatalog,self).__init__(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type)
        self.z_qso=z_qso
        self.confidence=confidence
        self.nhi=nhi



    @classmethod
    def init_from_fits(cls,name):
        catalog = Catalog.load_from_fits(name)
        if("Z_DLA" in catalog["DLA_CAT"].get_colnames()):
            coord_z = catalog["DLA_CAT"]["Z_DLA"][:]
            z_qso = catalog["DLA_CAT"]["Z_QSO"][:]
            coord = np.vstack([coord_z]).transpose()
            primary_key = catalog["DLA_CAT"]["THING_ID"][:]
            conf_dla = catalog["DLA_CAT"]["CONF_DLA"][:]
            nhi_dla = catalog["DLA_CAT"]["NHI_DLA"][:]
            catalog_type = "sky"
            Catalog.close(catalog)
            return(cls(name=name,coord=coord,primary_key=primary_key,z_qso=z_qso,confidence=conf_dla,nhi=nhi_dla,catalog_type=catalog_type))
        if("X" in catalog[1].get_colnames()):
            coord_ra = catalog["DLA_CAT"]["X"][:]
            coord_dec = catalog["DLA_CAT"]["Y"][:]
            coord_z = catalog["DLA_CAT"]['Z'][:]
            coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
            z_qso = catalog["DLA_CAT"]['Z_QSO'][:]
            catalog_type = "cartesian"
            Catalog.close(catalog)
            return(cls(name=name,coord=coord,z_qso=z_qso,catalog_type=catalog_type))

    @classmethod
    def init_from_pixel_catalog(cls,dla_pixels,name=None):
        if(name is None): name ="dla_catalog.fits"
        coord_ra = dla_pixels[:,0]
        coord_dec = dla_pixels[:,1]
        coord_z = dla_pixels[:,2]
        coord = np.vstack([coord_ra,coord_dec,coord_z]).transpose()
        z_qso = dla_pixels[:,3]
        catalog_type = "cartesian"
        return(cls(name=name,coord=coord,z_qso=z_qso,catalog_type=catalog_type))



    def cut_catalog_dla(self,coord_min=None,coord_max=None,confidence_min=None,nhi_min=None):
        mask_select = self.cut_catalog(coord_min=coord_min,coord_max=coord_max,center_x_coord=True) & (self.coord[:,2]!=-1.0)
        if(confidence_min is not None):
            mask_select &= (self.confidence > confidence_min)
        if(nhi_min is not None):
            mask_select &= (self.nhi > nhi_min)
        if(self.coord is not None):
            self.coord = self.coord[mask_select]
        if(self.primary_key is not None):
            self.primary_key = self.primary_key[mask_select]
        if(self.confidence is not None):
            self.confidence = self.confidence[mask_select]
        if(self.nhi is not None):
            self.nhi = self.nhi[mask_select]
        if(self.z_qso is not None):
            self.z_qso = self.z_qso[mask_select]



    def write(self):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        nrows = len(self.coord_ra)
        if(self.catalog_type == "sky"):
            h = np.zeros(nrows, dtype=[('THING_ID','i8'),('Z_QSO','f8'),('Z_DLA','f8'),('CONF_DLA','f8'),('NHI_DLA','f8')])
            h['THING_ID'] =self.primary_key
            h['Z_DLA'] = self.coord[:,0]
            h['Z_QSO'] = self.z_qso
            h['CONF_DLA'] = self.conf_dla
            h['NHI_DLA'] = self.nhi_dla
        elif(self.catalog_type == "cartesian"):
            h = np.zeros(nrows, dtype=[('X','f8'),('Y','f8'),('Z','f8'),('Z_QSO','f8')])
            h['X'] = self.coord[:,1]
            h['Y'] = self.coord[:,2]
            h['Z'] = self.coord[:,3]
            h['Z_QSO'] = self.z_qso
        fits.write(h)
        fits.close()



class GalaxyCatalog(Catalog):

    def __init__(self,name=None,coord=None,primary_key=None,confidence=None,standard_deviation=None,magnitude=None,catalog_type="sky"):
        super(GalaxyCatalog,self).__init__(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type)

        self.confidence = confidence
        self.standard_deviation = standard_deviation
        self.magnitude = magnitude


    @classmethod
    def init_from_fits(cls,name):
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
        return(cls(name=None,coord=coord,primary_key=primary_key,confidence=confidence,standard_deviation=standard_deviation,magnitude=magnitude,catalog_type=catalog_type))




    def cut_catalog_galaxy(self,coord_min=None,coord_max=None,standard_deviation_max=None,confidence_min=None,magnitude_max=None):
        mask_select = self.cut_catalog(coord_min=coord_min,coord_max=coord_max,center_x_coord=True)
        if(standard_deviation_max is not None):
            mask_select &= self.standard_deviation < standard_deviation_max
        if(magnitude_max is not None):
            mask_select &= self.magnitude < magnitude_max
        if(confidence_min is not None):
            mask_select &= self.confidence > confidence_min
        if(self.coord is not None):
            self.coord = self.coord[mask_select]
        if(self.primary_key is not None):
            self.primary_key = self.primary_key[mask_select]
        if(self.confidence is not None):
            self.confidence = self.confidence[mask_select]
        if(self.standard_deviation is not None):
            self.standard_deviation = self.standard_deviation[mask_select]
        if(self.magnitude is not None):
            self.magnitude = self.magnitude[mask_select]




class VoidCatalog(Catalog):

    def __init__(self,name=None,coord=None,primary_key=None,radius=None,weights=None,crossing_param=None,central_value=None,filling_factor=None,mean_value=None,catalog_type="sky"):
        super(VoidCatalog,self).__init__(name=name,coord=coord,primary_key=primary_key,catalog_type=catalog_type)

        self.radius = radius
        self.crossing_param = crossing_param
        self.central_value = central_value
        self.filling_factor = filling_factor
        self.weights = weights
        self.mean_value = mean_value



    @classmethod
    def init_from_fits(cls,name):
        catalog = Catalog.load_from_fits(name)
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
        crossing_param, central_value, filling_factor = None, None, None
        if("CROSSING" in catalog[1].get_colnames()): crossing_param = catalog[1]["CROSSING"][:]
        if("VALUE" in catalog[1].get_colnames()): central_value = catalog[1]["VALUE"][:]
        if("MEAN" in catalog[1].get_colnames()): mean_value = catalog[1]["MEAN"][:]
        if("FILLING_FACTOR" in catalog[1].read_header()): filling_factor = catalog[1].read_header()["FILLING_FACTOR"]

        Catalog.close(catalog)
        return(cls(name=name,coord=coord,primary_key=primary_key,radius=radius,weights=weights,crossing_param=crossing_param,central_value=central_value,filling_factor=filling_factor,mean_value=mean_value,catalog_type=catalog_type))


    @classmethod
    def init_from_dictionary(cls,name,radius,coord,catalog_type,other_arrays=None,other_array_names = None):
        central_value, weights, filling_factor, primary_key, crossing_param, mean_value = None, None, None, None, None, None
        if(other_array_names is not None):
            if("VALUE" in other_array_names):central_value = other_arrays[np.argwhere("VALUE" == np.asarray(other_array_names))[0][0]]
            if("MEAN" in other_array_names):mean_value = other_arrays[np.argwhere("MEAN" == np.asarray(other_array_names))[0][0]]
            if("WEIGHT" in other_array_names):weights = other_arrays[np.argwhere("WEIGHT" == np.asarray(other_array_names))[0][0]]
            if("FILLING_FACTOR" in other_array_names):filling_factor = other_arrays[np.argwhere("FILLING_FACTOR" == np.asarray(other_array_names))[0][0]]
            if("THING_ID" in other_array_names):primary_key = other_arrays[np.argwhere("THING_ID" == np.asarray(other_array_names))[0][0]]
            if("CROSSING" in other_array_names):crossing_param = other_arrays[np.argwhere("CROSSING" == np.asarray(other_array_names))[0][0]]
        return(cls(name=name,coord=coord,primary_key=primary_key,radius=radius,weights=weights,crossing_param=crossing_param,central_value=central_value,filling_factor=filling_factor,mean_value=mean_value,catalog_type=catalog_type))


    @classmethod
    def init_by_merging(cls,catalog_name,name=None):
        catalog = [VoidCatalog.init_from_fits(name) for name in catalog_name]
        catalog_type = catalog[0].catalog_type
        radius,coord,primary_key,crossing_param,central_value,mean_value,weights,filling_factor = None,None,None,None,None,None,None,None
        if(catalog[0].radius is not None):radius = np.concatenate([cat.radius for cat in catalog])
        if(catalog[0].coord is not None):coord = np.concatenate([cat.coord for cat in catalog])
        if(catalog[0].primary_key is not None):primary_key = np.concatenate([cat.primary_key for cat in catalog])
        if(catalog[0].crossing_param is not None):crossing_param = np.concatenate([cat.crossing_param for cat in catalog])
        if(catalog[0].central_value is not None):central_value = np.concatenate([cat.central_value for cat in catalog])
        if(catalog[0].mean_value is not None):mean_value = np.concatenate([cat.mean_value for cat in catalog])
        if(catalog[0].weights is not None):weights = np.concatenate([cat.weights for cat in catalog])
        if(catalog[0].filling_factor is not None):filling_factor = np.mean([cat.filling_factor for cat in catalog])
        return(cls(name=name,coord=coord,primary_key=primary_key,radius=radius,weights=weights,crossing_param=crossing_param,central_value=central_value,filling_factor=filling_factor,mean_value=mean_value,catalog_type=catalog_type))




    def write(self,qso_like=False):
        fits = fitsio.FITS(self.name,'rw',clobber=True)
        head = {}
        if(self.catalog_type.lower() == "sky"):
            h = self.create_qso_like_dictionary(qso_like=qso_like)
        elif(self.catalog_type.lower() == "cartesian"):
            h = self.create_cartesian_void_dictionary()
        self.update_dictionary(h,head)
        fits.write(h,header=head)
        fits.close()

    def create_sky_void_dictionary(self,qso_like=False):
        h = {}
        h['RA'] = self.coord[:,0].astype("f8")
        h['DEC'] = self.coord[:,1].astype("f8")
        h['Z'] = self.coord[:,2].astype("f8")
        if(qso_like):
            h['PLATE'] = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(self.radius))]).astype("i8")
            h['MJD'] = np.asarray(['1' + '{0:09d}'.format(i) for i in range(len(self.radius))]).astype("i8")
        return(h)

    def create_cartesian_void_dictionary(self):
        h = {}
        h['X'] = self.coord[:,0].astype("f8")
        h['Y'] = self.coord[:,1].astype("f8")
        h['Z'] = self.coord[:,2].astype("f8")
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
        if(self.filling_factor is not None):
            head["FILLING_FACTOR"] = self.filling_factor



    def compute_filling_factor(self,size=None,property_name=None):
        if(size is not None):
            map_size = size
        elif(property_name is not None):
            prop = MapPixelProperty(name=property_name)
            prop.read()
            map_size = prop.size
        volume_map = map_size[0]*map_size[1]*map_size[2]
        volume_void = np.sum((4/3)*np.pi*(np.array(self.radius))**3)
        self.filling_factor = volume_void/volume_map


    def create_crossing_criteria(self,pixel_name):
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

    def get_crossing_qso(self,qso_name):
        qso = QSOCatalog.init_from_fits(qso_name)
        coord_crossing_qso = []
        for i in range(len(self.radius)):
            diff_pixel = qso.coord - self.coord[i]
            distance_pixel = np.sqrt(diff_pixel[:,0]**2 + diff_pixel[:,1]**2)
            mask = distance_pixel < self.radius[i]
            coord_crossing_qso.append(qso.coord[mask])
        return(coord_crossing_qso)



    def cut_catalog_void(self,method_cut,coord_min=None,coord_max=None,cut_crossing_param=None,pixel_name=None,cut_radius=None,distance_map_name=None,distance_map_prop=None,distance_map_param=None,distance_map_percent=None,cut_border_prop=None):
        mask_select = self.cut_catalog(coord_min=coord_min,coord_max=coord_max,center_x_coord=False)
        string_to_add = ""
        if type(method_cut) == str :
            if method_cut == "ALL":
                method_cut = ["CROSSING","RADIUS","DIST","BORDER"]
            else:
                method_cut = [method_cut]

        if("CROSSING" in method_cut):
            if ((cut_crossing_param is None)|(pixel_name is None)) : raise KeyError("Give a crossing parameter and a Pixel file name")
            self.create_crossing_criteria(pixel_name)
            mask_select &= self.cut_crossing_parameter(cut_crossing_param)
            string_to_add = string_to_add + f"_crossing{cut_crossing_param}"
        if("RADIUS" in method_cut):
            if (cut_radius is None) : raise KeyError("Give a radius cutting parameter")
            mask_select &= self.cut_radius(cut_radius)
            string_to_add = string_to_add + f"_cutradius{cut_radius}"
        if("BORDER" in method_cut):
            if cut_border_prop is None : raise KeyError("Give a property file for border cutting")
            mask_select &= self.cut_border(cut_border_prop)
            string_to_add = string_to_add + "_cutborder"
        if("DIST" in method_cut):
            if ((distance_map_name is None)|(distance_map_prop is None)|(distance_map_param is None)|(distance_map_percent is None)) is None : raise KeyError("Give a dist map parameter, map and property file please")
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
        if(self.mean_value is not None):self.mean_value = self.mean_value[mask]
        if(self.weights is not None):self.weights = self.weights[mask]


    def cut_crossing_parameter(self,cut_crossing_param):
        mask = self.crossing_param > cut_crossing_param
        return(mask)

    def cut_radius(self,cut_radius):
        mask = (self.radius >= cut_radius[0])&(self.radius < cut_radius[1])
        return(mask)

    def cut_border(self,cut_border_prop):
        prop = MapPixelProperty(name=cut_border_prop)
        prop.read()
        cut_border_size = prop.size
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