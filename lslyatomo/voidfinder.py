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
from multiprocessing import Pool
from scipy.stats import ks_2samp
from functools import partial
from lslyatomo import tomographic_objects,utils
from scipy.optimize import curve_fit



def create_merged_catalog(pwd,list_catalog_name,merged_catalog_name):
    void_merged = tomographic_objects.VoidCatalog.init_by_merging(list_catalog_name,name=os.path.join(pwd,merged_catalog_name))
    void_merged.write()

def cut_catalog(pwd,catalog_name,method_cut,cut_crossing_param=None,
                pixel_name=None,cut_radius=None,distance_map_name=None,
                distance_map_prop=None,distance_map_param=None,
                distance_map_percent=None,cut_border_prop=None):
    void_cut = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    cut_catalog_name = void_cut.cut_catalog_void(method_cut,cut_crossing_param=cut_crossing_param,
                                                 pixel_name=pixel_name,cut_radius=cut_radius,
                                                 distance_map_name=distance_map_name,
                                                 distance_map_prop=distance_map_prop,
                                                 distance_map_param=distance_map_param,
                                                 distance_map_percent=distance_map_percent,
                                                 cut_border_prop=cut_border_prop)
    void_cut.name = os.path.join(pwd,cut_catalog_name)
    void_cut.write()

def compute_filling_factor(catalog_name,property_name):
    void = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    void.compute_filling_factor(property_name=property_name)
    void.write()

def get_crossing_qso(catalog_name,qso_name):
    void = tomographic_objects.VoidCatalog.init_from_fits(catalog_name)
    qso = void.get_crossing_qso(qso_name)
    return(qso)



#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################



class VoidFinder(object):

    def __init__(self,pwd,map_name,params_void_finder,map_shape=None,map_size=None,property_file=None,number_core=1,restart=False,find_cluster=False,split_map=None,split_overlap=None,delete_option="CLUSTERS"):
        self.pwd = pwd
        self.map_name = map_name
        self.params_void_finder = params_void_finder
        self.number_core = number_core
        self.restart = restart
        self.find_cluster=find_cluster
        self.split_map = split_map
        self.split_overlap = split_overlap
        self.delete_option = delete_option

        self.tomographic_map = tomographic_objects.TomographicMap.init_classic(name=map_name,shape=map_shape,size=map_size,property_file=property_file)
        self.tomographic_map.read()

        log_name = f"void_finder_report_{self.get_name_catalog()}.txt"
        self.log = utils.create_report_log(name=os.path.join(self.pwd,log_name))


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


    def find_voids(self):
        if(self.split_map is None):
            return(self.find_voids_single_map())
        else :
            return(self.find_voids_map_split())

    def find_voids_single_map(self):
        if(self.params_void_finder["method"] == "WATERSHED"):
            (radius, coord_Mpc,other_arrays,other_array_names) = self.find_voids_watershed(self.tomographic_map)
        elif(self.params_void_finder["method"] == "SPHERICAL"):
            (radius, coord_Mpc,other_arrays,other_array_names) = self.find_voids_sphere(self.tomographic_map)
        else :
            raise ValueError("The method_void chosen is not implemented, try : WATERSHED or SPHERICAL")
        self.save_voids(radius, coord_Mpc,other_arrays,other_array_names,self.tomographic_map.coordinate_transform,self.tomographic_map.Omega_m,self.tomographic_map.boundary_cartesian_coord,self.tomographic_map.boundary_sky_coord)
        return(radius, coord_Mpc)


    def find_voids_map_split(self):
        Chunks = self.split_map_in_chunks()
        if(self.params_void_finder["method"] == "WATERSHED"):
            if(self.number_core > 1):
                list_chunks = list(Chunks.items())
                list_map_chunks = [tomographic_objects.TomographicMap.init_classic(name=Chunks[str(i) + str(j)]["map_name"],shape=Chunks[str(i) + str(j)]["map_shape"],size=Chunks[str(i) + str(j)]["map_size"],map_array=Chunks[str(i) + str(j)]["map"]) for i in range(self.split_map[0]) for j in range(self.split_map[1])]
                pool = Pool(self.number_core)
                map_output = pool.map(self.find_voids_watershed,list_map_chunks)
                for i in range(len(list_chunks)):
                    Chunks[list_chunks[i][0]]["radius"],Chunks[list_chunks[i][0]]["coord_Mpc"],Chunks[list_chunks[i][0]]["other_arrays"],Chunks[list_chunks[i][0]]["other_array_names"]= map_output[i]
                del list_chunks,map_output
            else:
                for i in range(self.split_map[0]):
                    for j in range(self.split_map[1]):
                        map_chunks = tomographic_objects.TomographicMap.init_classic(name=Chunks[str(i) + str(j)]["map_name"],shape=Chunks[str(i) + str(j)]["map_shape"],size=Chunks[str(i) + str(j)]["map_size"],map_array=Chunks[str(i) + str(j)]["map"])
                        Chunks[str(i) + str(j)]["radius"],Chunks[str(i) + str(j)]["coord_Mpc"],Chunks[str(i) + str(j)]["other_arrays"],Chunks[str(i) + str(j)]["other_array_names"]  = self.find_voids_watershed(map_chunks)
        elif(self.params_void_finder["method"] == "SPHERICAL"):
            for i in range(self.split_map[0]):
                for j in range(self.split_map[1]):
                    map_chunks = tomographic_objects.TomographicMap.init_classic(name=Chunks[str(i) + str(j)]["map_name"],shape=Chunks[str(i) + str(j)]["map_shape"],size=Chunks[str(i) + str(j)]["map_size"],map_array=Chunks[str(i) + str(j)]["map"])
                    Chunks[str(i) + str(j)]["radius"],Chunks[str(i) + str(j)]["coord_Mpc"],Chunks[str(i) + str(j)]["other_arrays"],Chunks[str(i) + str(j)]["other_array_names"] = self.find_voids_sphere(map_chunks)
        else :
            raise ValueError("The method_void chosen is not implemented, try : WATERSHED or SPHERICAL")
        (radius, coord_Mpc,other_arrays,other_array_names) = self.merge_chunks(Chunks)
        new_coord_Mpc, new_radius, new_other_arrays = self.delete_voids(self.tomographic_map,radius,coord_Mpc,other_arrays=other_arrays,mpc=True)
        self.save_voids(new_radius, new_coord_Mpc,new_other_arrays,other_array_names,self.tomographic_map.coordinate_transform,self.tomographic_map.Omega_m,self.tomographic_map.boundary_cartesian_coord,self.tomographic_map.boundary_sky_coord)
        return(new_radius, new_coord_Mpc)


    def split_map_in_chunks(self):
        map_3D = self.tomographic_map.map_array
        number_Mpc_per_pixels = self.tomographic_map.mpc_per_pixel
        pixels_x = self.tomographic_map.shape[0]
        pixels_y = self.tomographic_map.shape[1]
        subIntervalx = pixels_x//self.split_map[0]
        subIntervaly = pixels_y//self.split_map[1]
        if(self.split_overlap is None):
            overlaping_x,overlaping_y = 0,0
        else :
            overlaping_x = int(np.round(self.split_overlap/number_Mpc_per_pixels[0],0))
            overlaping_y = int(np.round(self.split_overlap/number_Mpc_per_pixels[1],0))
        Chunks = {}
        for i in range(self.split_map[0]):
            for j in range(self.split_map[1]):
                Chunks[str(i) + str(j)]={}
                if((i==self.split_map[0]-1)&(i==0)):
                    pixel_x_interval = [i*subIntervalx, (i+1)*subIntervalx]
                elif i == 0 :
                    pixel_x_interval = [i*subIntervalx, (i+1)*subIntervalx + overlaping_x]
                elif i == self.split_map[0]-1 :
                    pixel_x_interval = [i*subIntervalx - overlaping_x, self.tomographic_map.shape[0]]
                else:
                    pixel_x_interval = [i*subIntervalx - overlaping_x, (i+1)*subIntervalx + overlaping_x]
                if((j==self.split_map[1]-1)&(j==0)):
                    pixel_y_interval = [  j*subIntervaly, (j+1)*subIntervaly]
                elif j == 0 :
                    pixel_y_interval = [  j*subIntervaly, (j+1)*subIntervaly + overlaping_y]
                elif j == self.split_map[1]-1 :
                    pixel_y_interval = [ j*subIntervaly - overlaping_y , self.tomographic_map.shape[1]]
                else:
                    pixel_y_interval = [  j*subIntervaly -overlaping_y, (j+1)*subIntervaly + overlaping_y]
                size_x_interval = np.array(pixel_x_interval)*number_Mpc_per_pixels[0]
                size_y_interval = np.array(pixel_y_interval)*number_Mpc_per_pixels[1]
                Chunks[str(i) + str(j)]["map"]=map_3D[pixel_x_interval[0]:pixel_x_interval[1],pixel_y_interval[0]:pixel_y_interval[1],:]
                Chunks[str(i) + str(j)]["map_size"]=(size_x_interval[1]-size_x_interval[0],size_y_interval[1]-size_y_interval[0],self.tomographic_map.size[2])
                Chunks[str(i) + str(j)]["map_min"]=(size_x_interval[0],size_y_interval[0],0)
                Chunks[str(i) + str(j)]["map_shape"]=(pixel_x_interval[1]-pixel_x_interval[0],pixel_y_interval[1]-pixel_y_interval[0],self.tomographic_map.shape[2])
                Chunks[str(i) + str(j)]["map_name"]="{}_{}_{}".format(self.map_name,i,j)
        return(Chunks)


    def merge_chunks(self,Chunks):
        radius_to_contatenate = []
        coord_Mpc_to_contatenate = []
        other_array_names = Chunks[list(Chunks.keys())[0]]["other_array_names"]
        other_arrays_to_contatenate = [[] for i in range(len(other_array_names))]
        for i in range(self.split_map[0]):
            for j in range(self.split_map[1]):
                coord_mpc = Chunks[str(i) + str(j)]["coord_Mpc"]
                if(coord_mpc.shape[0] !=0):
                    if(self.split_overlap is not None):
                        if((i==self.split_map[0]-1)&(i==0)):
                            pixel_x_interval = [0,0]
                        elif i == 0 :
                            pixel_x_interval = [0, self.split_overlap]
                        elif i == self.split_map[0]-1 :
                            pixel_x_interval = [ self.split_overlap,0]
                        else:
                            pixel_x_interval = [ self.split_overlap,self.split_overlap]
                        if((j==self.split_map[1]-1)&(j==0)):
                            pixel_y_interval = [0,0]
                        elif j == 0 :
                            pixel_y_interval = [0, self.split_overlap]
                        elif j == self.split_map[1]-1 :
                            pixel_y_interval =  [ self.split_overlap,0]
                        else:
                            pixel_y_interval = [ self.split_overlap,self.split_overlap]
                        mask = (coord_mpc[:,0] > pixel_x_interval[0])
                        mask &= (coord_mpc[:,0] < Chunks[str(i) + str(j)]["map_size"][0] - pixel_x_interval[1])
                        mask &= (coord_mpc[:,1] > pixel_y_interval[0])
                        mask &= (coord_mpc[:,1] < Chunks[str(i) + str(j)]["map_size"][1] - pixel_y_interval[1])
                    else:
                        mask = np.full(coord_mpc.shape[0],True)
                    radius_to_contatenate.append(Chunks[str(i) + str(j)]["radius"][mask])
                    for k in range(len(other_array_names)):
                        other_arrays_to_contatenate[k].append(np.array(Chunks[str(i) + str(j)]["other_arrays"][k])[mask])
                    coord_Mpc_to_contatenate.append((Chunks[str(i) + str(j)]["coord_Mpc"] + np.array(Chunks[str(i) + str(j)]["map_min"]))[mask])
        if(len(radius_to_contatenate) == 0):
            radius = np.empty(0)
            coord_Mpc = np.empty(0)
        else:
            radius = np.concatenate(radius_to_contatenate,axis=0)
            coord_Mpc = np.concatenate(coord_Mpc_to_contatenate,axis=0)
        for k in range(len(other_array_names)):
            other_arrays_to_contatenate[k] = np.concatenate(other_arrays_to_contatenate[k],axis=0)
        return(radius, coord_Mpc,other_arrays_to_contatenate,other_array_names)



    def find_voids_watershed(self,tomographic_map):
        self.log.add("Beginning of the Watershed finding for the map {}".format(tomographic_map.name))
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        global map_3D
        map_3D = tomographic_map.map_array
        if(self.find_cluster):
            mask = map_3D < self.params_void_finder["threshold"]
        else :
            mask = map_3D > self.params_void_finder["threshold"]
        index_under_density = np.argwhere(mask)
        self.log.add("Number of pixels for the map {} = {}".format(tomographic_map.name,len(index_under_density)))
        map_under_density = map_3D[mask]
        cluster_map,clusters = self.create_watershed_clusters(index_under_density)
        self.log.add("Pixel clusters created for the map {}".format(tomographic_map.name))
        centers = np.zeros((len(clusters),3))
        radius_shed = np.zeros(len(clusters))
        delta_max = np.zeros((len(clusters)))
        delta_mean = np.zeros((len(clusters)))
        volume_cell = number_Mpc_per_pixels[0]*number_Mpc_per_pixels[1]*number_Mpc_per_pixels[2]
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
        self.log.add("Computation of radius and center finished for the map {}".format(tomographic_map.name))
        mask_radius = radius_shed > self.params_void_finder["minimal_radius"]
        radius = radius_shed[mask_radius]
        coord = centers[mask_radius]
        delta_max = delta_max[mask_radius]
        delta_mean = delta_mean[mask_radius]
        self.log.add("Masking of low radius done for the map {}".format(tomographic_map.name))
        new_coord, new_radius, new_other_arrays = self.delete_voids(tomographic_map,radius,coord,other_arrays=[delta_max,delta_mean])
        other_array_names =["VALUE","MEAN"]
        new_coord_Mpc = self.convert_to_Mpc(tomographic_map,new_coord,new_radius)
        del map_3D,mask,mask_clust,mask_radius,cluster_map,clusters,map_under_density,centers,index_under_density,radius_shed
        self.log.add("End of the Watershed finding for the map {}".format(tomographic_map.name))
        return(new_radius, new_coord_Mpc,new_other_arrays,other_array_names)






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




    def find_voids_sphere(self,tomographic_map):
        self.log.add("Beginning of the Simple spherical finding for the map {}".format(tomographic_map.name))
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        number_pixel_maximal_radius = [int(round((self.params_void_finder["maximal_radius"]/number_Mpc_per_pixels)[0],0)),int(round((self.params_void_finder["maximal_radius"]/number_Mpc_per_pixels)[1],0)),int(round((self.params_void_finder["maximal_radius"]/number_Mpc_per_pixels)[2],0))]
        global map_3D
        map_3D = tomographic_map.map_array
        if(self.find_cluster):
            mask = map_3D < self.params_void_finder["threshold"]
        else :
            mask = map_3D > self.params_void_finder["threshold"]
        coord = np.argwhere(mask)
        del mask
        global indice
        indice = np.transpose(np.indices(map_3D.shape),axes=(1,2,3,0))
        radius = np.full(len(coord),-1)
        self.log.add("Number of pixels for the map {} = {}".format(tomographic_map.name,len(coord)))
        if(self.restart):
            radius = self.restart_calculation(radius,coord)
        mask = radius < 0
        radius_to_compute = radius[mask]
        coord_to_compute = coord[mask]
        mean_value = np.zeros(len(radius_to_compute))
        if(self.number_core > 1):
            self.log.add("Start of pool for the map {}".format(tomographic_map.name))
            if(self.find_cluster):
                func = partial(self.find_the_sphere_cluster,number_Mpc_per_pixels,number_pixel_maximal_radius)
            else:
                func = partial(self.find_the_sphere,number_Mpc_per_pixels,number_pixel_maximal_radius)
            pool = Pool(self.number_core)
            out_pool = np.array(pool.map(func,coord_to_compute))
            if(len(out_pool) !=0): radius_to_compute , mean_value = out_pool[:,0], out_pool[:,1]
            else: radius_to_compute , mean_value = np.array([]),np.array([])
            self.log.add("End of pool for the map {}".format(tomographic_map.name))
        else :
            if(self.find_cluster):
                for i in range(len(radius_to_compute)):
                    radius_to_compute[i],mean_value[i] = self.find_the_sphere_cluster(number_Mpc_per_pixels,number_pixel_maximal_radius,coord_to_compute[i])
            else:
                for i in range(len(radius_to_compute)):
                    radius_to_compute[i],mean_value[i] = self.find_the_sphere(number_Mpc_per_pixels,number_pixel_maximal_radius,coord_to_compute[i])
        radius[mask] = radius_to_compute
        del mask,radius_to_compute,coord_to_compute
        mask = radius == 0
        coord = coord[~mask]
        radius = radius[~mask]
        mean_value = mean_value[~mask]
        new_coord, new_radius, new_other_arrays = self.delete_voids(tomographic_map,radius,coord,other_arrays=[mean_value])
        nearest_coord = np.round(new_coord,0).astype(int)
        new_other_arrays.append(map_3D[nearest_coord[:,0],nearest_coord[:,1],nearest_coord[:,2]])
        new_coord_Mpc = self.convert_to_Mpc(tomographic_map,new_coord,new_radius)
        other_array_names=["MEAN","VALUE"]
        del indice,coord,map_3D,mask,radius
        self.log.add("End of the Simple spherical finding for the map {}".format(tomographic_map.name))
        return(new_radius, new_coord_Mpc,new_other_arrays,other_array_names)


    def find_the_sphere(self,number_Mpc_per_pixels,number_pixel_maximal_radius,index):
        map_local = map_3D[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]
        indice_local = (indice[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]- index)*number_Mpc_per_pixels
        distance_map = np.sqrt(indice_local[:,:,:,0]**2 + indice_local[:,:,:,1]**2 + indice_local[:,:,:,2]**2)
        del indice_local
        rayon = self.params_void_finder["minimal_radius"]
        boolean = np.mean(map_local[distance_map < rayon]) > self.params_void_finder["average"]
        while(boolean&(rayon<self.params_void_finder["maximal_radius"])):
            rayon = rayon + self.params_void_finder["radius_step"]
            mean_value = np.mean(map_local[distance_map < rayon])
            boolean = mean_value > self.params_void_finder["average"]
        del distance_map,map_local,boolean
        if((rayon>=self.params_void_finder["maximal_radius"])|(rayon<=self.params_void_finder["minimal_radius"])):
            return(0,0)
        else:
            return(rayon,mean_value)


    def find_the_sphere_cluster(self,number_Mpc_per_pixels,number_pixel_maximal_radius,index):
        map_local = map_3D[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]
        indice_local = (indice[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]- index)*number_Mpc_per_pixels
        distance_map = np.sqrt(indice_local[:,:,:,0]**2 + indice_local[:,:,:,1]**2 + indice_local[:,:,:,2]**2)
        del indice_local
        rayon = self.params_void_finder["minimal_radius"]
        boolean = np.mean(map_local[distance_map < rayon]) < self.params_void_finder["average"]
        while(boolean&(rayon<self.params_void_finder["maximal_radius"])):
            rayon = rayon + self.params_void_finder["radius_step"]
            mean_value = np.mean(map_local[distance_map < rayon])
            boolean = mean_value < self.params_void_finder["average"]
        del distance_map,map_local,boolean
        if((rayon>=self.params_void_finder["maximal_radius"])|(rayon<=self.params_void_finder["minimal_radius"])):
            return(0,0)
        else:
            return(rayon,mean_value)


### Overlapper deletion algorithms


    def delete_voids(self,tomographic_map,radius,coord,other_arrays=None,mpc=False):
        if(self.delete_option == "CLUSTERS"):
            new_coord, new_radius, new_other_arrays = self.delete_overlapers_clusters(tomographic_map,radius,coord,other_arrays=other_arrays,mpc=mpc)
        elif(self.delete_option == "ITERATION"):
            new_coord, new_radius, new_other_arrays = self.iterate_overlap_deletion(tomographic_map,radius,coord,other_arrays=other_arrays,mpc=mpc)
        elif(self.delete_option == "NONE"):
            True
        else:
            raise ValueError("The delete_option chosen is not implemented, try : CLUSTERS, ITERATION or NONE")
        if(other_arrays is not None):
            return(new_coord, new_radius, new_other_arrays)
        else:
            return(new_coord, new_radius)


    def iterate_overlap_deletion(self,tomographic_map,radius,coord,other_arrays=None,mpc=False):
        # CR - Weird, to optimize
        if other_arrays is not None:
            others_arrays_copies =[]
            for i in range(len(other_arrays)):
                others_arrays_copies.append(other_arrays[i].copy())
        else:
            others_arrays_copies =None
        radius_copy = radius.copy()
        coord_copy = coord.copy()
        nb_voids_delete = len(coord)
        new_coord,new_radius,new_others_arrays = coord,radius,other_arrays
        while(nb_voids_delete>0):
            new_coord,new_radius,nb_voids_delete,new_others_arrays = self.delete_overlapers(tomographic_map,radius_copy,coord_copy,other_arrays=others_arrays_copies,mpc=mpc)
            coord_copy = new_coord
            radius_copy = new_radius
            others_arrays_copies = new_others_arrays
        if(other_arrays is not None):
            return(new_coord,new_radius,new_others_arrays)
        else :
            return(new_coord,new_radius,None)




    def delete_overlapers(self,tomographic_map,radius,coord,other_arrays=None,mpc=False):
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        if(other_arrays is not None):
            other_arrays_clean = [[] for i in range(len(other_arrays))]
        else :
            other_arrays_clean = None
        coord_cleaned = []
        radius_cleaned = []
        radius_delete=radius.copy()
        mask = radius_delete > 0
        while(len(radius_delete[mask])!=0):
            coord_left=coord[mask]
            radius_left=radius_delete[mask]
            if(other_arrays_clean is not None):
                other_arrays_left= []
                for i in range(len(other_arrays)):
                    other_arrays_left.append(other_arrays[i][mask])
            else:
                other_arrays_left=None
            rad=radius_left[0]
            index = coord_left[0]
            sum_rad = (radius_left + rad)
            if(mpc):indice_normalized = (coord_left - index)
            else:indice_normalized = (coord_left - index)*number_Mpc_per_pixels
            dist_rad = np.sqrt(indice_normalized[:,0]**2 + indice_normalized[:,1]**2 + indice_normalized[:,2]**2)
            mask2 = sum_rad >= dist_rad
            maxi = np.argwhere(radius_left[mask2] == np.amax(radius_left[mask2])).flatten().tolist()
            coord_cleaned.append(np.mean(coord_left[mask2][maxi],axis=0))
            radius_cleaned.append(np.mean(radius_left[mask2][maxi],axis=0))
            if(other_arrays_clean is not None):
                for i in range(len(other_arrays)):
                    other_arrays_clean[i].append(np.mean(other_arrays_left[i][mask2][maxi]))
            radius_left[mask2]=0
            radius_delete[mask]=radius_left
            mask = radius_delete > 0
        del mask,radius_delete,radius_left,coord_left,indice_normalized,sum_rad,dist_rad
        if(other_arrays_clean is not None):
            for i in range(len(other_arrays_clean)):
                other_arrays_clean[i] = np.array(other_arrays_clean[i])
        return(np.array(coord_cleaned),np.array(radius_cleaned),len(coord)-len(coord_cleaned),other_arrays_clean)



    def find_overlapers(self,tomographic_map,index,rad,radius,coord,mpc=False):
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        sum_rad = (radius + rad)
        if(mpc):indice_normalized = (coord - index)
        else:indice_normalized = (coord - index)*number_Mpc_per_pixels
        dist_rad = np.sqrt(indice_normalized[:,0]**2 + indice_normalized[:,1]**2 + indice_normalized[:,2]**2)
        overlapers_mask = sum_rad >= dist_rad
        del sum_rad,dist_rad,indice_normalized
        return(overlapers_mask)



    def create_overlaper_map(self,tomographic_map,radius,coord,mpc=False):
        cluster_map = np.zeros(radius.shape,dtype = np.int64)
        mask_clusters = cluster_map == 0
        cluster_number = 0
        while(len(cluster_map[mask_clusters])!=0):
            arg = np.argwhere(mask_clusters)[0][0]
            rad = radius[arg]
            index = coord[arg]
            overlapers_mask = self.find_overlapers(tomographic_map,index,rad,radius,coord,mpc=mpc)
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



    def delete_overlapers_clusters(self,tomographic_map,radius,coord,other_arrays=None,mpc=False):
        if(other_arrays is not None):
            other_arrays_clean = [[] for i in range(len(other_arrays))]
        else :
            other_arrays_clean = None
        coord_cleaned = []
        radius_cleaned = []
        (cluster_map,clusters)=self.create_overlaper_map(tomographic_map,radius,coord,mpc=mpc)
        for c in clusters :
            mask = cluster_map == c
            radius_cluster = radius[mask]
            coord_cluster = coord[mask]
            if(other_arrays is not None):
                other_arrays_cluster= []
                for i in range(len(other_arrays)):
                    other_arrays_cluster.append(other_arrays[i][mask])
            maxi = np.argwhere(radius_cluster == np.amax(radius_cluster)).flatten().tolist()
            coord_cleaned.append(np.mean(coord_cluster[maxi],axis=0))
            radius_cleaned.append(np.mean(radius_cluster[maxi],axis=0))
            if(other_arrays is not None):
                for i in range(len(other_arrays)):
                    other_arrays_clean[i].append(np.mean(other_arrays_cluster[i][maxi]))
        if(other_arrays is not None):
            return(np.array(coord_cleaned),np.array(radius_cleaned),other_arrays_clean)
        else :
            return(np.array(coord_cleaned),np.array(radius_cleaned),None)



### Other algorithms



    def restart_calculation(self,radius,coord):
        file = open(self.restart,"r")
        line_file = file.readlines()
        file.close()
        for i in range(len(line_file)):
            line = line_file[i].strip().split()
            coord_line = line_file[i].split('[')[-1].split(']')[0]
            if((line[5]=="Sphere")&(line[6]=="found")):
                index = np.array([int(coord_line.split()[0]),
                                  int(coord_line.split()[1]),
                                  int(coord_line.split()[2])])
                mask = coord == index
                arg = np.argwhere(mask)
                del mask
                radius[arg]=float(line[10].split(",")[0])
        del line, line_file
        self.log.add("restart coordinates")
        mask = radius >= 0
        self.log.add("Sphere found : R = " + str(radius[mask]) + ", Coord = " + str(coord[mask]))
        return(radius)



    def convert_to_Mpc(self,tomographic_map,new_coord,new_radius):
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        if(new_coord.shape[0] != 0):
            coord_Mpc = np.zeros(new_coord.shape)
            coord_Mpc = new_coord * np.array(number_Mpc_per_pixels)
        else:
            coord_Mpc = new_coord
        return(coord_Mpc)



    def save_voids(self,radius,coord,other_arrays,other_array_names,coordinate_transform,Omega_m,boundary_cartesian_coord,boundary_sky_coord):
        dict_void = {"R" : radius, "COORD" : coord}
        for i in range(len(other_arrays)):
            dict_void[other_array_names[i]] = other_arrays[i]
        name = os.path.join(self.pwd,f"Catalog_{self.get_name_catalog()}.fits")
        void = tomographic_objects.VoidCatalog.init_from_dictionary(name,radius,coord,"cartesian",coordinate_transform,Omega_m,boundary_cartesian_coord,boundary_sky_coord,other_arrays=other_arrays,other_array_names = other_array_names)
        void.write()







class PlotVoid(object):

    def __init__(self,pwd,void_catalog,nb_bins=30):
        self.pwd = pwd
        self.void = tomographic_objects.Catalog.init_catalog_from_fits(void_catalog, "void")
        self.nb_bins = nb_bins



    def plot_histo_radius(self,rmin,rmax,name,norm=False):
        plt.figure()
        bins = np.linspace(rmin,rmax, self.nb_bins)
        plt.hist(self.void.radius, bins, alpha=1, label='x',density=norm)
        plt.grid()
        plt.ylabel("Number of voids")
        plt.xlabel("Radius of the void in Mpc.h-1")
        plt.savefig(os.path.join(self.pwd,f"{name}_histo_radius.pdf"), format ="pdf")
        plt.show()
        plt.close()


    def plot_radius_redshift(self,name):
        redshift = self.void.redshift
        plt.figure()
        plt.plot(redshift,self.void.radius,"b.")
        plt.xlabel("Redshift")
        plt.ylabel("Radius size")
        plt.grid()
        plt.savefig(os.path.join(self.pwd,f"{name}_radius_redshift.pdf"), format ="pdf")
        plt.show()
        plt.close()

    def plot_meanradius_redshift(self,name):
        redshift_catalog = self.void.redshift
        maxredshift = np.max(redshift_catalog)
        minredshift = np.min(redshift_catalog)
        redshift = np.linspace(minredshift,maxredshift,self.nb_bins)
        mean_radius = []
        for i in range(self.nb_bins):
            mask = (redshift_catalog < minredshift + (i+1/self.nb_bins)*((maxredshift - minredshift)))
            mask &= (redshift_catalog > minredshift + (i/self.nb_bins)*((maxredshift - minredshift)))
            mean_radius.append(np.mean(self.void.radius[mask]))
        plt.figure()
        plt.plot(redshift,mean_radius,"b.")
        plt.xlabel("Redshift")
        plt.ylabel("Radius size")
        plt.grid()
        plt.savefig(os.path.join(self.pwd,f"{name}_meanradius_redshift.pdf"), format ="pdf")
        plt.show()
        plt.close()


    def plot_histo_redshift(self,name):
        redshift = self.void.redshift
        plt.figure()
        plt.hist(redshift,self.nb_bins)
        plt.xlabel("Redshift")
        plt.ylabel("Number of voids")
        plt.grid()
        plt.savefig(os.path.join(self.pwd,f"{name}_histo_redshift.pdf"), format ="pdf")
        plt.show()
        plt.close()


    def load_and_plot_comparison(self,catalog_name,name,legend,rmin,rmax,norm=False,log_scale=True,factor_add=None,expo_fit_rmin=None,other_catalogs=None):
        if(other_catalogs is not None):
            other_radius = []
            for i in range(len(other_catalogs)):
                other_void = tomographic_objects.Catalog.init_catalog_from_fits(other_catalogs[i], "void")
                if(factor_add is not None):
                    other_radius.apped(other_void.radius*factor_add[i])
                else:
                    other_radius.apped(other_void.radius)
        else:
            other_radius = None
        catalog = tomographic_objects.Catalog.init_catalog_from_fits(catalog_name, "void")
        KS_stat, p_value = ks_2samp(self.void.radius,catalog.radius)
        self.plot_void_histogram_comparison(catalog.radius,name,legend,rmin,rmax,norm=norm,log_scale=log_scale,expo_fit_rmin=expo_fit_rmin,other_radius=other_radius)
        return(KS_stat, p_value)


    def plot_void_histogram_comparison(self,radius2,name,legend,rmin,rmax,norm=False,log_scale=True,expo_fit_rmin=None,other_radius=None):
        plt.figure()
        bins = np.linspace(rmin,rmax, self.nb_bins)
        (n1, bins1, patches1)  = plt.hist(self.void.radius, bins, alpha=0.5, label=legend[0],density=norm)
        (n2, bins2, patches2)  = plt.hist(radius2, bins, alpha=0.5, label=legend[0],density=norm)
        if(other_radius is not None):
            for i in range(len(other_radius)):
                (n_other, bins_other, patches_other)  = plt.hist(other_radius[i], bins, label=legend[i+1],histtype='step',linestyle='dashed',ec="k")
        if(expo_fit_rmin is not None):
            bin_center = (bins1[1:] + bins1[0:-1]) / 2
            mask = bin_center>expo_fit_rmin
            n1,n2,bins_fit = n1[mask],n2[mask],bin_center[mask]
            fit_function = lambda x,a,b : np.exp(a*x+b)
            fit = curve_fit(fit_function,bins_fit,n1)
            fit2 = curve_fit(fit_function,bins_fit,n2)
            plt.plot(bins_fit,fit_function(bins_fit,*fit[0]),'r')
            plt.plot(bins_fit,fit_function(bins_fit,*fit2[0]),'b')
            perr = np.sqrt(np.diag(fit[1]))
            perr2 = np.sqrt(np.diag(fit2[1]))
        else:
            perr,perr2 = None,None
        if(log_scale):
            plt.yscale("log")
        plt.grid()
        plt.xlabel("Void radius in Mpc.h" + r"$^{-1}$")
        if(log_scale):
            plt.savefig(os.path.join(self.pwd,"histogram_voids_radius_"+ name + "log.pdf"),format="pdf")
        else :
            plt.savefig(os.path.join(self.pwd,"histogram_voids_radius_"+ name + "lin.pdf"),format="pdf")
        plt.show()
        plt.close()
        return(perr,perr2)
