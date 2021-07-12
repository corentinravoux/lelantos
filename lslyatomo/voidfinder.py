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



import os
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from scipy.stats import ks_2samp
from functools import partial
from lslyatomo import tomographic_objects,utils
from scipy.optimize import curve_fit






def create_merged_catalog(pwd,list_catalog_name,merged_catalog_name):
    void_merged = tomographic_objects.VoidCatalog.init_by_merging(list_catalog_name,name=os.path.join(pwd,merged_catalog_name))
    void_merged.write()


def compute_additional_stats(catalog_name,pixel_name):
    void = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    void.compute_filling_factor()
    void.compute_crossing_criteria(pixel_name)
    void.compute_los_distance(pixel_name)
    void.write()



def cut_catalog(pwd,
                catalog_name,
                method_cut,
                coord_min=None,
                coord_max=None,
                cut_crossing_param=None,
                cut_radius=None,
                distance_map_name=None,
                distance_map_prop=None,
                distance_map_param=None,
                distance_map_percent=None):
    void_cut = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    cut_catalog_name = void_cut.cut_catalog_void(method_cut,
                                                 coord_min=coord_min,
                                                 coord_max=coord_max,
                                                 cut_crossing_param=cut_crossing_param,
                                                 cut_radius=cut_radius,
                                                 distance_map_name=distance_map_name,
                                                 distance_map_prop=distance_map_prop,
                                                 distance_map_param=distance_map_param,
                                                 distance_map_percent=distance_map_percent)
    void_cut.name = os.path.join(pwd,cut_catalog_name)
    void_cut.write()
    return(void_cut.name)

def correct_void_coordinates(pwd,catalog_name,method,name_out,inv_g_function,pixel_name,**kwargs):
    void = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    void.correct_coordinates(method,name_out,inv_g_function,pixel_name,**kwargs)
    void.name = os.path.join(pwd,name_out)
    void.write()



def get_crossing_qso(catalog_name,qso_name):
    void = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    qso = void.get_crossing_qso(qso_name)
    return(qso)


def create_qso_like_catalog(catalog_name):
    void = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    void.name = f"""{void.name.split(".fits")[0]}_qso_like.fits"""
    void.convert_to_cross_corr_radec()
    void.write(qso_like=True)

def qso_to_3d(catalog_name,new_name,moveaxis=None):
    qso = tomographic_objects.QSOCatalog.init_from_fits(catalog_name)
    qso.writetxt(new_name,moveaxis=moveaxis)

def void_to_3d(catalog_name,new_name,moveaxis=None):
    void = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    void.writetxt(new_name,moveaxis=moveaxis)





#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################



class VoidFinder(object):

    def __init__(self,pwd,map_name,params_void_finder,map_shape=None,map_size=None,map_property_file=None,number_core=1,find_cluster=False,split_map=None,split_overlap=None,restart=False,delete_option="CLUSTERS"):
        self.pwd = pwd
        self.params_void_finder = params_void_finder
        self.number_core = number_core
        self.find_cluster=find_cluster
        self.split_map = split_map
        self.split_overlap = split_overlap
        self.delete_option = delete_option
        self.restart = restart

        self.map_name = map_name
        self.map_shape = map_shape
        self.map_size = map_size
        self.map_property_file = map_property_file
        self.map_mpc_per_pixel = None
        self.map_coordinate_transform = None
        self.map_Omega_m = None
        self.map_boundary_cartesian_coord = None
        self.map_boundary_sky_coord = None

        log_name = f"void_finder_report_{self.get_name_catalog()}.txt"
        self.log = utils.create_report_log(name=os.path.join(self.pwd,log_name))

        self.save_temporary_file = False
        self.temporary_file_name = "{}" + f"_temporary_void_catalog_{self.get_name_catalog()}.fits"

    def initialize_finder(self):
        tomographic_map = tomographic_objects.TomographicMap.init_classic(name=self.map_name,
                                                                          shape=self.map_shape,
                                                                          size=self.map_size,
                                                                          property_file=self.map_property_file)
        tomographic_map.read()
        if(self.map_shape is None): self.map_shape = tomographic_map.shape
        if(self.map_size is None): self.map_size = tomographic_map.size
        self.map_mpc_per_pixel = tomographic_map.mpc_per_pixel
        self.map_coordinate_transform = tomographic_map.coordinate_transform
        self.map_Omega_m = tomographic_map.Omega_m
        self.map_boundary_cartesian_coord = tomographic_map.boundary_cartesian_coord
        self.map_boundary_sky_coord = tomographic_map.boundary_sky_coord
        map_array = tomographic_map.map_array
        del tomographic_map
        return(map_array)



    def find_voids(self):
        map_array = self.initialize_finder()
        if(self.split_map is None):
            return(self.find_voids_single_map(map_array))
        else :
            return(self.find_voids_map_split(map_array))


    def find_voids_single_map(self,map_array):
        if(self.params_void_finder["method"].upper() == "WATERSHED"):
            (radius, coord, other_array, other_array_name) = self.find_voids_watershed((self.map_name,
                                                                                       self.map_mpc_per_pixel,
                                                                                       map_array))
        elif(self.params_void_finder["method"].upper() == "SPHERICAL"):
            (radius, coord, other_array, other_array_name) = self.find_voids_sphere((self.map_name,
                                                                                    self.map_mpc_per_pixel,
                                                                                    map_array))
        else :
            raise ValueError("The method_void chosen is not implemented, try : WATERSHED or SPHERICAL")
        del map_array
        name = self.save_voids(radius, coord,other_array,other_array_name,
                               self.map_coordinate_transform,self.map_Omega_m,
                               self.map_boundary_cartesian_coord,
                               self.map_boundary_sky_coord)
        return(name)


    def find_voids_map_split(self,map_array):
        self.save_temporary_file = True
        self.log.add(f"Splitting of the map in chunks: {self.split_map}")
        self.log.add("Saving of temporary files activated")
        map_chunks = self.split_map_in_chunks(map_array)
        del map_array
        list_index_map_chunks = [f'{i:03d}' + f'{j:03d}' for j in range(self.split_map[1]) for i in range(self.split_map[0])]
        if(self.restart):
            self.log.add("Restarting of the void finder, searching for temporary files...")
            other_array_name_restart =["VALUE","MEAN"] # CR - might change in the future the treatment of other_array_name variable
            (map_chunks,list_index_map_chunks_restart) = self.restart_calculation(map_chunks,list_index_map_chunks,other_array_name_restart)
            list_index_map_chunks = list_index_map_chunks_restart
        if(self.params_void_finder["method"] == "WATERSHED"):
            if(self.number_core > 1):
                pool = mp.Pool(self.number_core)
                list_map_name = [map_chunks[list_index_map_chunks[i]]["map_name"] for i in range(len(list_index_map_chunks))]
                list_map_mpc_per_pixel = [map_chunks[list_index_map_chunks[i]]["map_mpc_per_pixel"] for i in range(len(list_index_map_chunks))]
                list_map_array = [map_chunks[list_index_map_chunks[i]]["map_array"] for i in range(len(list_index_map_chunks))]
                out_finder = pool.map(self.find_voids_watershed,zip(list_map_name,list_map_mpc_per_pixel,list_map_array))
                del list_map_name,list_map_mpc_per_pixel,list_map_array
                for i in range(len(list_index_map_chunks)):
                    map_chunks[list_index_map_chunks[i]]["radius"] = out_finder[i][0]
                    map_chunks[list_index_map_chunks[i]]["coord"] = out_finder[i][1]
                    map_chunks[list_index_map_chunks[i]]["other_array"] = out_finder[i][2]
                    map_chunks[list_index_map_chunks[i]]["other_array_name"] = out_finder[i][3]
            else:
                for i in range(len(list_index_map_chunks)):
                    out_finder = self.find_voids_watershed((map_chunks[list_index_map_chunks[i]]["map_name"],
                                                           map_chunks[list_index_map_chunks[i]]["map_mpc_per_pixel"],
                                                           map_chunks[list_index_map_chunks[i]]["map_array"]))
                    map_chunks[list_index_map_chunks[i]]["radius"] = out_finder[0]
                    map_chunks[list_index_map_chunks[i]]["coord"] = out_finder[1]
                    map_chunks[list_index_map_chunks[i]]["other_array"] = out_finder[2]
                    map_chunks[list_index_map_chunks[i]]["other_array_name"] = out_finder[3]
        elif(self.params_void_finder["method"] == "SPHERICAL"):
            for i in range(len(list_index_map_chunks)):
                out_finder = self.find_voids_sphere((map_chunks[list_index_map_chunks[i]]["map_name"],
                                                    map_chunks[list_index_map_chunks[i]]["map_mpc_per_pixel"],
                                                    map_chunks[list_index_map_chunks[i]]["map_array"]))
                map_chunks[list_index_map_chunks[i]]["radius"] = out_finder[0]
                map_chunks[list_index_map_chunks[i]]["coord"] = out_finder[1]
                map_chunks[list_index_map_chunks[i]]["other_array"] = out_finder[2]
                map_chunks[list_index_map_chunks[i]]["other_array_name"] = out_finder[3]
        else :
            raise ValueError("The method_void chosen is not implemented, try : WATERSHED or SPHERICAL")
        self.log.add("End of the split finding procedure. Merging catalogs")
        (radius, coord,other_array,other_array_name) = self.merge_chunks(map_chunks)
        coord_clean, radius_clean, other_array_clean = self.delete_voids(self.map_mpc_per_pixel,radius,coord,other_array=other_array,mpc=True)
        name = self.save_voids(radius_clean, coord_clean,other_array_clean,other_array_name,
                               self.map_coordinate_transform,self.map_Omega_m,
                               self.map_boundary_cartesian_coord,
                               self.map_boundary_sky_coord)
        self.delete_temporary_files(map_chunks,list_index_map_chunks)
        del map_chunks,list_index_map_chunks
        return(name)


    def split_map_in_chunks(self,map_array):
        subIntervalx = self.map_shape[0]//self.split_map[0]
        subIntervaly = self.map_shape[1]//self.split_map[1]
        if(self.split_overlap is None):
            overlaping_x, overlaping_y = 0,0
        elif(type(self.split_overlap) == int):
            overlaping_x, overlaping_y = self.split_overlap,self.split_overlap
        else :
            overlaping_x, overlaping_y = self.split_overlap[0],self.split_overlap[1]
        map_chunks = {}
        for i in range(self.split_map[0]):
            for j in range(self.split_map[1]):
                map_chunks[f'{i:03d}' + f'{j:03d}']={}
                if((i==self.split_map[0]-1)&(i==0)):
                    pixel_x_interval = [0, self.map_shape[0]]
                elif i == 0 :
                    pixel_x_interval = [0, subIntervalx + overlaping_x]
                elif i == self.split_map[0]-1 :
                    pixel_x_interval = [i*subIntervalx - overlaping_x, self.map_shape[0]]
                else:
                    pixel_x_interval = [i*subIntervalx - overlaping_x, (i+1)*subIntervalx + overlaping_x]
                if((j==self.split_map[1]-1)&(j==0)):
                    pixel_y_interval = [0, self.map_shape[1]]
                elif j == 0 :
                    pixel_y_interval = [0 , subIntervaly + overlaping_y]
                elif j == self.split_map[1]-1 :
                    pixel_y_interval = [ j*subIntervaly - overlaping_y , self.map_shape[1]]
                else:
                    pixel_y_interval = [ j*subIntervaly - overlaping_y, (j+1)*subIntervaly + overlaping_y]
                min_x_interval = pixel_x_interval[0]*self.map_mpc_per_pixel[0]
                min_y_interval = pixel_y_interval[0]*self.map_mpc_per_pixel[1]

                map_shape = (pixel_x_interval[1]-pixel_x_interval[0],
                             pixel_y_interval[1]-pixel_y_interval[0],
                             self.map_shape[2])
                map_size = utils.get_map_size(map_shape, self.map_mpc_per_pixel)
                map_chunks[f'{i:03d}' + f'{j:03d}']["map_name"]=f"{self.map_name}_{i:03d}{j:03d}"
                map_chunks[f'{i:03d}' + f'{j:03d}']["map_mpc_per_pixel"]= self.map_mpc_per_pixel
                map_chunks[f'{i:03d}' + f'{j:03d}']["map_array"]=map_array[pixel_x_interval[0]:pixel_x_interval[1],
                                                                           pixel_y_interval[0]:pixel_y_interval[1],:]
                map_chunks[f'{i:03d}' + f'{j:03d}']["map_size"]=map_size
                map_chunks[f'{i:03d}' + f'{j:03d}']["map_min"]=(min_x_interval,min_y_interval,0)
        return(map_chunks)

    def merge_chunks(self,map_chunks):
        radius_to_contatenate = []
        coord_to_contatenate = []
        if(self.split_overlap is not None):
            if(type(self.split_overlap) == int):
                overlaping_x, overlaping_y = self.split_overlap,self.split_overlap
            else :
                overlaping_x, overlaping_y = self.split_overlap[0],self.split_overlap[1]
        other_array_name = map_chunks[list(map_chunks.keys())[0]]["other_array_name"]
        other_array = [[] for i in range(len(other_array_name))]
        for i in range(self.split_map[0]):
            for j in range(self.split_map[1]):
                coord_chunks = map_chunks[f'{i:03d}' + f'{j:03d}']["coord"]
                mpc_per_pixel = map_chunks[f'{i:03d}' + f'{j:03d}']["map_mpc_per_pixel"]
                if(coord_chunks.shape[0] !=0):
                    if(self.split_overlap is not None):
                        if((i==self.split_map[0]-1)&(i==0)):
                            pixel_x_interval = [0.,0.]
                        elif i == 0 :
                            pixel_x_interval = [0, overlaping_x]
                        elif i == self.split_map[0]-1 :
                            pixel_x_interval = [overlaping_x, 0]
                        else:
                            pixel_x_interval = [overlaping_x, overlaping_x]
                        if((j==self.split_map[1]-1)&(j==0)):
                            pixel_y_interval = [0,0]
                        elif j == 0 :
                            pixel_y_interval = [0, overlaping_y]
                        elif j == self.split_map[1]-1 :
                            pixel_y_interval = [overlaping_y, 0]
                        else:
                            pixel_y_interval = [overlaping_y, overlaping_y]
                        mask = (coord_chunks[:,0] > pixel_x_interval[0]*mpc_per_pixel[0])
                        mask &= (coord_chunks[:,0] < map_chunks[f'{i:03d}' + f'{j:03d}']["map_size"][0] - (pixel_x_interval[1]*mpc_per_pixel[0]))
                        mask &= (coord_chunks[:,1] > pixel_y_interval[0]*mpc_per_pixel[1])
                        mask &= (coord_chunks[:,1] < map_chunks[f'{i:03d}' + f'{j:03d}']["map_size"][1] - (pixel_y_interval[1]*mpc_per_pixel[1]))
                    else:
                        mask = np.full(coord_chunks.shape[0],True)
                    radius_to_contatenate.append(map_chunks[f'{i:03d}' + f'{j:03d}']["radius"][mask])
                    for k in range(len(other_array_name)):
                        other_array[k].append(np.array(map_chunks[f'{i:03d}' + f'{j:03d}']["other_array"][k])[mask])
                    coord_to_contatenate.append((map_chunks[f'{i:03d}' + f'{j:03d}']["coord"]
                                                 + np.array(map_chunks[f'{i:03d}' + f'{j:03d}']["map_min"]))[mask])
        if(len(radius_to_contatenate) == 0):
            radius = np.empty(0)
            coord = np.empty(0)
        else:
            radius = np.concatenate(radius_to_contatenate,axis=0)
            coord = np.concatenate(coord_to_contatenate,axis=0)
        for k in range(len(other_array_name)):
            other_array[k] = np.concatenate(other_array[k],axis=0)
        return(radius, coord,other_array,other_array_name)




    def find_voids_watershed(self,parameters):
        map_name,map_mpc_per_pixel,map_array = parameters
        self.log.add(f"Beginning of the Watershed finding for the map {map_name}")
        if(self.find_cluster):
            mask = map_array < self.params_void_finder["threshold"]
        else :
            mask = map_array > self.params_void_finder["threshold"]
        index_under_density = np.argwhere(mask)
        self.log.add(f"Number of pixels for the map {map_name} = {len(index_under_density)}")
        map_under_density = map_array[mask]
        del map_array
        cluster_map,clusters = self.create_watershed_clusters(index_under_density)
        self.log.add(f"Pixel clusters created for the map {map_name}")
        centers = np.zeros((len(clusters),3))
        radius_shed = np.zeros(len(clusters))
        delta_max = np.zeros((len(clusters)))
        delta_mean = np.zeros((len(clusters)))
        volume_cell = map_mpc_per_pixel[0]*map_mpc_per_pixel[1]*map_mpc_per_pixel[2]
        mask_clust = None
        for i in range(len(clusters)):
            mask_clust = cluster_map == clusters[i]
            if(self.find_cluster):
                arg_center = np.argmin(map_under_density[mask_clust])
            else :
                arg_center = np.argmax(map_under_density[mask_clust])
            delta_max[i] = map_under_density[mask_clust][arg_center]
            delta_mean[i] = np.mean(map_under_density[mask_clust])
            centers[i] = index_under_density[mask_clust][arg_center]
            volume_shed = len(map_under_density[mask_clust]) * volume_cell
            radius_shed[i] = ((3 * volume_shed)/(4*np.pi))**(1/3)
        self.log.add(f"Computation of radius and center finished for the map {map_name}")
        mask_radius = radius_shed > self.params_void_finder["minimal_radius"]
        radius = radius_shed[mask_radius]
        coord = centers[mask_radius]
        delta_max = delta_max[mask_radius]
        delta_mean = delta_mean[mask_radius]
        self.log.add(f"Masking of low radius done for the map {map_name}")
        coord_clean, radius_clean, other_array_clean = self.delete_voids(map_mpc_per_pixel,radius,coord,other_array=[delta_max,delta_mean])
        other_array_name =["VALUE","MEAN"]
        coord_clean = self.convert_to_Mpc(map_mpc_per_pixel,coord_clean)
        del mask,mask_clust,mask_radius,cluster_map,clusters,map_under_density,centers,index_under_density,radius_shed
        self.log.add(f"End of the Watershed finding for the map {map_name}")
        if(self.save_temporary_file):
            self.save_temporary_catalog(map_name,radius_clean, coord_clean, other_array_clean, other_array_name)
        return(radius_clean, coord_clean, other_array_clean, other_array_name)






    def create_watershed_clusters(self,indices):
        cluster_map = np.zeros(indices[:,0].shape,dtype = np.int64)
        mask_clusters = cluster_map == 0
        cluster_number = 0
        while(len(cluster_map[mask_clusters])!=0):
            indice_normalized = (indices - indices[np.argwhere(mask_clusters)[0][0]])
            dist_index = np.sqrt(indice_normalized[:,0]**2 + indice_normalized[:,1]**2 + indice_normalized[:,2]**2)
            mask_dist = dist_index <= self.params_void_finder["dist_clusters"]
            clust = np.unique(cluster_map[mask_dist])
            clust = clust[clust!=0]
            if(len(clust)==0):
                cluster_number+=1
                cluster_map[mask_dist]=cluster_number
            else :
                if(len(clust) == 1):
                    mask_clust = mask_dist & (cluster_map == 0)
                    cluster_map[mask_clust]=clust[0]
                else :
                    cluster_number+=1
                    for c in clust:
                        mask_clust = cluster_map == c
                        cluster_map[mask_clust] = cluster_number
            mask_clusters = cluster_map == 0
        clusters = np.unique(cluster_map)
        return(cluster_map,clusters)



    def create_watershed_clusters2(self,indices):
        from sklearn.cluster import AgglomerativeClustering
        linkage="simple"
        cluster_map = AgglomerativeClustering(distance_threshold=1.5,n_clusters=None,linkage=linkage).fit(indices).labels_
        clusters = np.unique(cluster_map)
        return(cluster_map,clusters)

    def create_watershed_clusters3(self,indices):
        from scipy.cluster.hierarchy import fclusterdata
        cluster_map = fclusterdata(indices, t=1,criterion="distance")
        clusters = np.unique(cluster_map)
        return(cluster_map,clusters)



    def find_voids_sphere(self,parameters):
        global map_array_spherical_voidfinder
        map_name,map_mpc_per_pixel,map_array_spherical_voidfinder = parameters
        self.log.add(f"Beginning of the Simple spherical finding for the map {map_name}")
        maximal_radius = np.around(self.params_void_finder["maximal_radius"]/map_mpc_per_pixel,decimals=0).astype(int)
        global indice_spherical_voidfinder
        indice_spherical_voidfinder = np.transpose(np.indices(map_array_spherical_voidfinder.shape),axes=(1,2,3,0))

        if(self.find_cluster):
            mask = map_array_spherical_voidfinder < self.params_void_finder["threshold"]
        else :
            mask = map_array_spherical_voidfinder > self.params_void_finder["threshold"]
        coord = np.argwhere(mask)
        del mask
        self.log.add(f"Number of pixels for the map {map_name} = {coord.shape[0]}")

        if(self.number_core > 1):
            self.log.add(f"{self.number_core} processes used, start of multiprocessing routines")
            self.log.add(f"Start of pool for the map {map_name}")
            func = partial(self.find_the_sphere,map_mpc_per_pixel,maximal_radius)
            with mp.Pool(self.number_core) as pool:
                pool_results = np.array(pool.map(func,coord))
            if(pool_results.shape[0] !=0):
                radius , mean_value = pool_results[:,0], pool_results[:,1]
            else:
                radius , mean_value = np.array([]),np.array([])
            self.log.add(f"End of pool for the map {map_name}")
        else :
            self.log.add("Start of serial calculation")
            radius = np.zeros(coord.shape[0])
            mean_value = np.zeros(coord.shape[0])
            for i in range(len(coord)):
                radius[i],mean_value[i] = self.find_the_sphere(map_mpc_per_pixel,
                                                               maximal_radius,
                                                               coord[i])
        del indice_spherical_voidfinder,maximal_radius
        mask = radius == 0
        coord = coord[~mask]
        radius = radius[~mask]
        mean_value = mean_value[~mask]
        coord_clean, radius_clean, other_array_clean = self.delete_voids(map_mpc_per_pixel,
                                                                   radius,coord,
                                                                   other_array=[mean_value])
        nearest_coord = np.round(coord_clean,0).astype(int)
        other_array_clean.append(map_array_spherical_voidfinder[nearest_coord[:,0],nearest_coord[:,1],nearest_coord[:,2]])
        coord_clean = self.convert_to_Mpc(map_mpc_per_pixel,coord_clean)
        del map_array_spherical_voidfinder,coord,mask,radius,nearest_coord,map_mpc_per_pixel
        other_array_name=["MEAN","VALUE"]
        self.log.add(f"End of the Simple spherical finding for the map {map_name}")
        if(self.save_temporary_file):
            self.save_temporary_catalog(map_name,radius_clean, coord_clean, other_array_clean, other_array_name)
        return(radius_clean, coord_clean,other_array_clean,other_array_name)




    def find_the_sphere(self,mpc_per_pixel,maximal_radius,coord):
        map_local = map_array_spherical_voidfinder[max(coord[0]-maximal_radius[0],0):min(map_array_spherical_voidfinder.shape[0],coord[0]+maximal_radius[0]),
                                                   max(coord[1]-maximal_radius[1],0):min(map_array_spherical_voidfinder.shape[1],coord[1]+maximal_radius[1]),
                                                   max(coord[2]-maximal_radius[2],0):min(map_array_spherical_voidfinder.shape[2],coord[2]+maximal_radius[2])]
        distance_map = (indice_spherical_voidfinder[max(coord[0]-maximal_radius[0],0):min(map_array_spherical_voidfinder.shape[0],coord[0]+maximal_radius[0]),
                                                    max(coord[1]-maximal_radius[1],0):min(map_array_spherical_voidfinder.shape[1],coord[1]+maximal_radius[1]),
                                                    max(coord[2]-maximal_radius[2],0):min(map_array_spherical_voidfinder.shape[2],coord[2]+maximal_radius[2])]- coord)*mpc_per_pixel
        distance_map = np.sqrt(distance_map[:,:,:,0]**2 + distance_map[:,:,:,1]**2 + distance_map[:,:,:,2]**2)
        radius = self.params_void_finder["minimal_radius"]
        if(self.find_cluster):
            boolean = np.mean(map_local[distance_map < radius]) < self.params_void_finder["average"]
        else:
            boolean = np.mean(map_local[distance_map < radius]) > self.params_void_finder["average"]
        while(boolean&(radius<self.params_void_finder["maximal_radius"])):
            radius = radius + self.params_void_finder["radius_step"]
            mean_value = np.mean(map_local[distance_map < radius])
            if(self.find_cluster):
                boolean = mean_value < self.params_void_finder["average"]
            else:
                boolean = mean_value > self.params_void_finder["average"]
        del distance_map,map_local,boolean
        if((radius>=self.params_void_finder["maximal_radius"])|(radius<=self.params_void_finder["minimal_radius"])):
            radius, mean_value= 0, 0
        return(radius,mean_value)



    def delete_voids(self,mpc_per_pixel,radius,coord,other_array=None,mpc=False):
        if(self.delete_option == "CLUSTERS"):
            new_coord, new_radius, new_other_array = self.delete_overlapers_clusters(mpc_per_pixel,radius,coord,other_array=other_array,mpc=mpc)
        elif(self.delete_option == "ITERATION"):
            new_coord, new_radius, new_other_array = self.iterate_overlap_deletion(mpc_per_pixel,radius,coord,other_array=other_array,mpc=mpc)
        elif(self.delete_option == "NONE"):
            new_coord, new_radius, new_other_array = coord, radius, other_array
        else:
            raise ValueError("The delete_option chosen is not implemented, try : CLUSTERS, ITERATION or NONE")
        if(other_array is not None):
            return(new_coord, new_radius, new_other_array)
        else:
            return(new_coord, new_radius)


    def iterate_overlap_deletion(self,mpc_per_pixel,radius,coord,other_array=None,mpc=False):
        # CR - Weird, to optimize
        if other_array is not None:
            others_arrays_copies =[]
            for i in range(len(other_array)):
                others_arrays_copies.append(other_array[i].copy())
        else:
            others_arrays_copies =None
        radius_copy = radius.copy()
        coord_copy = coord.copy()
        nb_voids_delete = len(coord)
        new_coord,new_radius,new_others_arrays = coord,radius,other_array
        while(nb_voids_delete>0):
            new_coord,new_radius,nb_voids_delete,new_others_arrays = self.delete_overlapers(mpc_per_pixel,radius_copy,coord_copy,other_array=others_arrays_copies,mpc=mpc)
            coord_copy = new_coord
            radius_copy = new_radius
            others_arrays_copies = new_others_arrays
        if(other_array is not None):
            return(new_coord,new_radius,new_others_arrays)
        else :
            return(new_coord,new_radius,None)




    def delete_overlapers(self,mpc_per_pixel,radius,coord,other_array=None,mpc=False):
        if(other_array is not None):
            other_array_clean = [[] for i in range(len(other_array))]
        else :
            other_array_clean = None
        coord_cleaned = []
        radius_cleaned = []
        radius_delete=radius.copy()
        mask = radius_delete > 0
        while(len(radius_delete[mask])!=0):
            coord_left=coord[mask]
            radius_left=radius_delete[mask]
            if(other_array_clean is not None):
                other_array_left= []
                for i in range(len(other_array)):
                    other_array_left.append(other_array[i][mask])
            else:
                other_array_left=None
            rad=radius_left[0]
            index = coord_left[0]
            sum_rad = (radius_left + rad)
            if(mpc):indice_normalized = (coord_left - index)
            else:indice_normalized = (coord_left - index)*mpc_per_pixel
            dist_rad = np.sqrt(indice_normalized[:,0]**2 + indice_normalized[:,1]**2 + indice_normalized[:,2]**2)
            mask2 = sum_rad >= dist_rad
            maxi = np.argwhere(radius_left[mask2] == np.amax(radius_left[mask2])).flatten().tolist()
            coord_cleaned.append(np.mean(coord_left[mask2][maxi],axis=0))
            radius_cleaned.append(np.mean(radius_left[mask2][maxi],axis=0))
            if(other_array_clean is not None):
                for i in range(len(other_array)):
                    other_array_clean[i].append(np.mean(other_array_left[i][mask2][maxi]))
            radius_left[mask2]=0
            radius_delete[mask]=radius_left
            mask = radius_delete > 0
        del mask,radius_delete,radius_left,coord_left,indice_normalized,sum_rad,dist_rad
        if(other_array_clean is not None):
            for i in range(len(other_array_clean)):
                other_array_clean[i] = np.array(other_array_clean[i])
        return(np.array(coord_cleaned),np.array(radius_cleaned),len(coord)-len(coord_cleaned),other_array_clean)



    def find_overlapers(self,mpc_per_pixel,index,rad,radius,coord,mpc=False):
        sum_rad = (radius + rad)
        if(mpc):indice_normalized = (coord - index)
        else:indice_normalized = (coord - index)*mpc_per_pixel
        dist_rad = np.sqrt(indice_normalized[:,0]**2 + indice_normalized[:,1]**2 + indice_normalized[:,2]**2)
        overlapers_mask = sum_rad >= dist_rad
        del sum_rad,dist_rad,indice_normalized
        return(overlapers_mask)



    def create_overlaper_map(self,mpc_per_pixel,radius,coord,mpc=False):
        cluster_map = np.zeros(radius.shape,dtype = np.int64)
        mask_clusters = cluster_map == 0
        cluster_number = 0
        while(len(cluster_map[mask_clusters])!=0):
            arg = np.argwhere(mask_clusters)[0][0]
            rad = radius[arg]
            index = coord[arg]
            overlapers_mask = self.find_overlapers(mpc_per_pixel,index,rad,radius,coord,mpc=mpc)
            clust = np.array(list(set(cluster_map[overlapers_mask])))
            clust = clust[clust!=0]
            if(len(clust)==0):
                cluster_number+=1
                cluster_map[overlapers_mask]=cluster_number
            else :
                if(len(clust) == 1):
                    mask_clust = overlapers_mask & (cluster_map == 0)
                    cluster_map[mask_clust]=clust[0]
                else :
                    cluster_number+=1
                    for c in clust:
                        mask_clust = cluster_map == c
                        cluster_map[mask_clust] = cluster_number
            mask_clusters = cluster_map == 0
        clusters = np.unique(cluster_map)
        return(cluster_map,clusters)



    def delete_overlapers_clusters(self,mpc_per_pixel,radius,coord,other_array=None,mpc=False):
        if(other_array is not None):
            other_array_clean = [[] for i in range(len(other_array))]
        else :
            other_array_clean = None
        coord_cleaned = []
        radius_cleaned = []
        (cluster_map,clusters)=self.create_overlaper_map(mpc_per_pixel,radius,coord,mpc=mpc)
        for c in clusters :
            mask = cluster_map == c
            radius_cluster = radius[mask]
            coord_cluster = coord[mask]
            if(other_array is not None):
                other_array_cluster= []
                for i in range(len(other_array)):
                    other_array_cluster.append(other_array[i][mask])
            maxi = np.argwhere(radius_cluster == np.amax(radius_cluster)).flatten().tolist()
            coord_cleaned.append(np.mean(coord_cluster[maxi],axis=0))
            radius_cleaned.append(np.mean(radius_cluster[maxi],axis=0))
            if(other_array is not None):
                for i in range(len(other_array)):
                    other_array_clean[i].append(np.mean(other_array_cluster[i][maxi]))
        if(other_array is not None):
            return(np.array(coord_cleaned),np.array(radius_cleaned),other_array_clean)
        else :
            return(np.array(coord_cleaned),np.array(radius_cleaned),None)





    def convert_to_Mpc(self,mpc_per_pixel,coord):
        if(coord.shape[0] != 0):
            coord = coord * np.array(mpc_per_pixel)
        return(coord)


    def save_temporary_catalog(self,map_name,radius, coord, other_array, other_array_name):
        dict_void = {"R" : radius, "COORD" : coord}
        for i in range(len(other_array)):
            dict_void[other_array_name[i]] = other_array[i]
        name = self.temporary_file_name.format(map_name)
        void = tomographic_objects.VoidCatalog.init_from_dictionary(name,radius,coord,"cartesian",
                                                                    self.map_coordinate_transform,
                                                                    self.map_Omega_m,
                                                                    self.map_boundary_cartesian_coord,
                                                                    self.map_boundary_sky_coord,
                                                                    other_array=other_array,
                                                                    other_array_name = other_array_name)
        void.write()

    def delete_temporary_files(self,map_chunks,list_index_map_chunks):
        for i in range(len(list_index_map_chunks)):
            map_name = map_chunks[list_index_map_chunks[i]]["map_name"]
            file_name = self.temporary_file_name.format(map_name)
            if(os.path.isfile(file_name)):
                os.remove(file_name)

    def restart_calculation(self,map_chunks,list_index_map_chunks,other_array_name):
        list_index_map_chunks_restart = []
        for i in range(len(list_index_map_chunks)):
            map_name = map_chunks[list_index_map_chunks[i]]["map_name"]
            tmp_file_name = self.temporary_file_name.format(map_name)
            if(os.path.isfile(tmp_file_name)):
                self.log.add(f"Temporary file {tmp_file_name} found and added to the calculation")
                tmp_catalog = tomographic_objects.VoidCatalog.init_from_fits(tmp_file_name)
                map_chunks[list_index_map_chunks[i]]["radius"] = tmp_catalog.radius
                map_chunks[list_index_map_chunks[i]]["coord"] = tmp_catalog.coord
                map_chunks[list_index_map_chunks[i]]["other_array"] = tmp_catalog.return_array_list(other_array_name)
                map_chunks[list_index_map_chunks[i]]["other_array_name"] = other_array_name
            else:
                list_index_map_chunks_restart.append(list_index_map_chunks[i])
        return(map_chunks,list_index_map_chunks_restart)

    def save_voids(self,radius,coord,other_array,other_array_name,coordinate_transform,Omega_m,boundary_cartesian_coord,boundary_sky_coord):
        dict_void = {"R" : radius, "COORD" : coord}
        for i in range(len(other_array)):
            dict_void[other_array_name[i]] = other_array[i]
        name = os.path.join(self.pwd,f"Catalog_{self.get_name_catalog()}.fits")
        void = tomographic_objects.VoidCatalog.init_from_dictionary(name,radius,coord,"cartesian",
                                                                    coordinate_transform,
                                                                    Omega_m,boundary_cartesian_coord,
                                                                    boundary_sky_coord,
                                                                    other_array=other_array,
                                                                    other_array_name = other_array_name)
        void.write()
        return(name)

    def get_name_catalog(self):
        if (self.find_cluster):
            name_out= "Clusters"
        else:
            name_out= "Voids"
        if(self.params_void_finder["method"]=="SPHERICAL"):
            name = f"""{name_out}_{self.params_void_finder["method"]}_{self.params_void_finder["threshold"]}threshold_{self.params_void_finder["average"]}average_{self.params_void_finder["minimal_radius"]}rmin_{self.delete_option}_deletion"""
        elif(self.params_void_finder["method"]=="WATERSHED"):
            name = f"""{name_out}_{self.params_void_finder["method"]}_{self.params_void_finder["threshold"]}threshold_{self.params_void_finder["dist_clusters"]}dist_clusters_{self.params_void_finder["minimal_radius"]}rmin_{self.delete_option}_deletion"""
        else :
            raise ValueError("The method_void chosen is not implemented, try : WATERSHED or SPHERICAL")
        return(name)





class PlotVoid(object):

    def __init__(self,pwd,void_catalog):
        self.pwd = pwd
        self.void = tomographic_objects.Catalog.init_catalog_from_fits(void_catalog, "void")

    def load_catalog(self,comparison,value_name):
        comparison_redshift,comparison_value = None,None
        if(comparison is not None):
            comparison_value = []
            comparison_redshift = []
            for i in range(len(comparison)):
                catalog = tomographic_objects.Catalog.init_catalog_from_fits(comparison[i], "void")
                comparison_value.append(getattr(catalog,value_name))
                comparison_redshift.append(catalog.redshift)
        value = getattr(self.void,value_name)
        return(value,comparison_value,comparison_redshift)


    def plot_histo(self,value_name,name,comparison=None,
                   comparison_legend=None,loaded_value=None,
                   **kwargs):
        if(loaded_value is None):
            (value,comparison_value,comparison_redshift) = self.load_catalog(comparison,value_name)
        else:
            value , comparison_value = loaded_value, comparison
        utils.save_histo(self.pwd,value,value_name,
                         name,comparison=comparison_value,
                         comparison_legend=comparison_legend,
                         **kwargs)

    def plot_mean_redshift_dependence(self,value_name,name,
                                      comparison=None,
                                      comparison_redshift=None,
                                      comparison_legend=None,
                                      loaded_value=None,
                                      **kwargs):
        if(loaded_value is None):
            (value,comparison_value,comparison_redshift) = self.load_catalog(comparison,value_name)
        else:
            value , comparison_value = loaded_value, comparison
        redshift = self.void.redshift
        utils.save_mean_redshift_dependence(self.pwd,value,redshift,
                                            value_name,name,
                                            comparison=comparison_value,
                                            comparison_redshift=comparison_redshift,
                                            comparison_legend=None,
                                            **kwargs)

    def plot_redshift_dependence(self,value_name,name,
                                 comparison=None,
                                 comparison_redshift=None,
                                 comparison_legend=None,
                                 loaded_value=None,
                                 **kwargs):
        if(loaded_value is None):
            (value,comparison_value,comparison_redshift) = self.load_catalog(comparison,value_name)
        else:
            value , comparison_value = loaded_value, comparison
        redshift = self.void.redshift
        utils.save_redshift_dependence(self.pwd,value,redshift,
                                            value_name,name,
                                            comparison=comparison_value,
                                            comparison_redshift=comparison_redshift,
                                            comparison_legend=None,
                                            **kwargs)



    def plot(self,value_names,name,
             comparison=None,comparison_legend=None,
             histo=True,mean_z_dependence=True,
             z_dependence=True,
             **kwargs):
        for value_name in value_names:
            (value,comparison_value,comparison_redshift) = self.load_catalog(comparison,value_name)
            if(histo):
                self.plot_histo(value_name,name,
                                comparison=comparison_value,
                                comparison_legend=comparison_legend,
                                loaded_value=value,
                                **kwargs)
            if(mean_z_dependence)&(value_name!="redshift"):
                self.plot_mean_redshift_dependence(value_name,name,
                                                   comparison=comparison_value,
                                                   comparison_redshift=comparison_redshift,
                                                   comparison_legend=comparison_legend,
                                                   loaded_value=value,
                                                   **kwargs)
            if(z_dependence)&(value_name!="redshift"):
                self.plot_redshift_dependence(value_name,name,
                                              comparison=comparison_value,
                                              comparison_redshift=comparison_redshift,
                                              comparison_legend=comparison_legend,
                                              loaded_value=value,
                                              **kwargs)



    def compute_ks_stat(self,comparison):
        (radius,comparison) = self.load_catalog(comparison,"radius")
        KS_stat, p_value = ks_2samp(radius,comparison[0])
        return(KS_stat,p_value)



    def plot_radius_histo_fit_expo(self,expo_fit_rmin=0,comparison=None,**kwargs):
        fit_function = lambda x,a,b : np.exp(a*x+b)

        (radius,comparison) = self.load_catalog(comparison,"radius")
        (name, n, bins, patches) = utils.plot_histo(radius,"radius","",**kwargs)
        bin_center = (bins[1:] + bins[0:-1]) / 2
        mask = bin_center>expo_fit_rmin
        n,bins_fit = n[mask],bin_center[mask]
        fit = curve_fit(fit_function,bins_fit,n)
        plt.plot(bins_fit,fit_function(bins_fit,*fit[0]),'r')
        perr = np.sqrt(np.diag(fit[1]))

        perr_comparison = []
        if(comparison is not None):
            for i in range(len(comparison)):
                (name, n, bins, patches) = utils.plot_histo(comparison[i],"radius","",**kwargs)
                bin_center = (bins[1:] + bins[0:-1]) / 2
                mask = bin_center>expo_fit_rmin
                n,bins_fit = n[mask],bin_center[mask]
                fit = curve_fit(fit_function,bins_fit,n)
                plt.plot(bins_fit,fit_function(bins_fit,*fit[0]),'b')
                perr_comparison.append(np.sqrt(np.diag(fit[1])))
        return(perr,perr_comparison)
