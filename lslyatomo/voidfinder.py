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



import os,pickle,logging
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.stats import ks_2samp
from functools import partial
from lslyatomo import tomographic_objects



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
        self.create_Report()
        
    def create_Report(self):
        currentdir = os.getcwd()
        os.chdir(self.pwd)
        logging.basicConfig(filename='Python_Report',level=logging.INFO,format='%(asctime)s :: %(levelname)s :: %(message)s')
        os.chdir(currentdir)  
        
    def add_Report(self,line):
        logging.info(line)


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
        self.save_voids(radius, coord_Mpc,other_arrays,other_array_names)
        return(radius, coord_Mpc)


    def find_voids_map_split(self):
        Chunks = self.split_map_in_chunks()
        if(self.params_void_finder["method"] == "WATERSHED"):
            if(self.number_core !=1):
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
        self.save_voids(new_radius, new_coord_Mpc,new_other_arrays,other_array_names)
        return(new_radius, new_coord_Mpc)

    def cutInChunks(self,Om,numberChunks,overlaping,shuffle=None):
        if overlaping == None :
            overlaping = 0.0
        (x,y,z,deltas,sigmas,maxlist,zqso,redshift,redshift_qso,zdlas) = self.getdata(Om)
        if(shuffle is not None):
            if(shuffle["method"] == "delta") :            
                deltas = self.shuffle_deltas(deltas)
            elif(shuffle["method"] == "radec"):
                x,y = self.shuffle_radec(x,y)
            elif(shuffle["method"] == "delta_sigma") :            
                deltas,sigmas = self.shuffle_deltas_sigmas(deltas,sigmas)
            self.createDachshundInput(Om,shuffle["name"],x,y,z,deltas,sigmas)
        pickle.dump([x,y,zqso], open("DataQSOposition.pickle","wb"))
        if(self.dla_catalog is not None):
            dla_array = []    
            for i in range(len(zdlas)):
                for j in range(len(zdlas[i])):
                    dla_array.append([x[i],y[i],zdlas[i][j]])
            dla_array = np.transpose(np.array(dla_array))
            pickle.dump(dla_array, open("DataDLAposition.pickle","wb"))
        minx,maxx,miny,maxy= tuple(maxlist[0:4])
        intervalx = (maxx - minx)
        intervaly = (maxy - miny)
        subIntervalx = intervalx/numberChunks[0]
        subIntervaly = intervaly/numberChunks[1]
        Chunks = {}
        for i in range(numberChunks[0]):
            for j in range(numberChunks[1]):
                Chunks[str(i) + str(j)]=[[],[],[],[],[]]
        for i in range(numberChunks[0]):
            for j in range(numberChunks[1]):
                if((i==numberChunks[0]-1)&(i==0)):
                    intervalxChunk = [i*subIntervalx, (i+1)*subIntervalx]
                elif i == 0 :
                    intervalxChunk = [i*subIntervalx, (i+1)*subIntervalx + overlaping]
                elif i == numberChunks[0]-1 :
                    intervalxChunk = [i*subIntervalx - overlaping, maxx - minx]
                else:
                    intervalxChunk = [i*subIntervalx - overlaping, (i+1)*subIntervalx + overlaping]
                if((j==numberChunks[1]-1)&(j==0)):
                    intervalyChunk = [  j*subIntervaly, (j+1)*subIntervaly]
                elif j == 0 :
                    intervalyChunk = [  j*subIntervaly, (j+1)*subIntervaly + overlaping]
                elif j == numberChunks[1]-1 :
                    intervalyChunk = [ j*subIntervaly - overlaping , maxy - miny]
                else:
                    intervalyChunk = [  j*subIntervaly -overlaping, (j+1)*subIntervaly + overlaping]
                for k in range(len(x)):
                    if((x[k]<intervalxChunk[1])&(x[k]>=intervalxChunk[0])&(y[k]<intervalyChunk[1])&(y[k]>=intervalyChunk[0])):
                        Chunks[str(i) + str(j)][0].append(x[k])
                        Chunks[str(i) + str(j)][1].append(y[k])
                        Chunks[str(i) + str(j)][2].append(z[k])
                        Chunks[str(i) + str(j)][3].append(deltas[k])
                        Chunks[str(i) + str(j)][4].append(sigmas[k])
                Chunks[str(i) + str(j)].append([intervalxChunk[0],intervalxChunk[1],intervalyChunk[0],intervalyChunk[1],maxlist[4],maxlist[5]])
                Chunks[str(i) + str(j)][0] = np.asarray(Chunks[str(i) + str(j)][0]) - intervalxChunk[0]
                Chunks[str(i) + str(j)][1] = np.asarray(Chunks[str(i) + str(j)][1]) - intervalyChunk[0]
        Chunks["overlaping"]=overlaping
        return(Chunks)


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
        radius = np.concatenate(radius_to_contatenate,axis=0)
        coord_Mpc = np.concatenate(coord_Mpc_to_contatenate,axis=0)
        for k in range(len(other_array_names)):
            other_arrays_to_contatenate[k] = np.concatenate(other_arrays_to_contatenate[k],axis=0)
        return(radius, coord_Mpc,other_arrays_to_contatenate,other_array_names)

        
                
    def find_voids_watershed(self,tomographic_map):
        self.add_Report("Beginning of the Watershed finding for the map {}".format(tomographic_map.name))
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        global map_3D
        map_3D = self.tomographic_map.map_array
        if(self.find_cluster):
            mask = map_3D < self.params_void_finder["threshold"]
        else :
            mask = map_3D > self.params_void_finder["threshold"]
        index_under_density = np.argwhere(mask)
        self.add_Report("Number of pixels for the map {} = {}".format(tomographic_map.name,len(index_under_density)))
        map_under_density = map_3D[mask]
        cluster_map,clusters = self.create_watershed_clusters(index_under_density)
        self.add_Report("Pixel clusters created for the map {}".format(tomographic_map.name))
        centers = np.zeros((len(clusters),3))
        radius_shed = np.zeros(len(clusters))
        delta_max = np.zeros((len(clusters)))
        volume_cell = number_Mpc_per_pixels[0]*number_Mpc_per_pixels[1]*number_Mpc_per_pixels[2]
        for i in range(len(clusters)):
            mask_clust = cluster_map == clusters[i]
            if(self.find_cluster):
                arg_center = np.argmin(map_under_density[mask_clust])
            else :
                arg_center = np.argmax(map_under_density[mask_clust])
            delta_max[i] = map_under_density[mask_clust][arg_center]
            centers[i] = index_under_density[mask_clust][arg_center]
            volume_shed = len(map_under_density[mask_clust]) * volume_cell
            radius_shed[i] = ((3 * volume_shed)/(4*np.pi))**(1/3)
        self.add_Report("Computation of radius and center finished for the map {}".format(tomographic_map.name))
        mask_radius = radius_shed > self.params_void_finder["minimal_radius"]
        radius = radius_shed[mask_radius]
        coord = centers[mask_radius]
        delta_max = delta_max[mask_radius]
        self.add_Report("Masking of low radius done for the map {}".format(tomographic_map.name))
        new_coord, new_radius, new_other_arrays = self.delete_voids(tomographic_map,radius,coord,other_arrays=[delta_max])
        other_array_names =["VALUE"]
        new_coord_Mpc = self.convert_to_Mpc(tomographic_map,new_coord,new_radius)
        del map_3D,mask,mask_clust,mask_radius,cluster_map,clusters,map_under_density,centers,index_under_density,radius_shed
        self.add_Report("End of the Watershed finding for the map {}".format(tomographic_map.name))
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
        self.add_Report("Beginning of the Simple spherical finding for the map {}".format(tomographic_map.name))
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixels
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
        self.add_Report("Number of pixels for the map {} = {}".format(tomographic_map.name,len(coord)))
        if(self.restart):
            radius = self.restart_calculation(radius,coord)            
        mask = radius < 0
        radius_to_compute = radius[mask]
        coord_to_compute = coord[mask]
        if(self.number_core != 0):
            self.add_Report("Start of pool for the map {}".format(tomographic_map.name))
            if(self.find_cluster):
                func = partial(self.find_the_sphere_cluster,number_Mpc_per_pixels,number_pixel_maximal_radius)
            else:
                func = partial(self.find_the_sphere,number_Mpc_per_pixels,number_pixel_maximal_radius)
            pool = Pool(self.number_core)
            radius_to_compute = pool.map(func,coord_to_compute)
            self.add_Report("End of pool for the map {}".format(tomographic_map.name))
        else :
            if(self.find_cluster):
                for i in range(len(radius_to_compute)):
                    radius_to_compute[i] = self.find_the_sphere_cluster(number_Mpc_per_pixels,number_pixel_maximal_radius,coord_to_compute[i])
            else:
                for i in range(len(radius_to_compute)):
                    radius_to_compute[i] = self.find_the_sphere(number_Mpc_per_pixels,number_pixel_maximal_radius,coord_to_compute[i])                
        radius[mask] = radius_to_compute
        del mask,radius_to_compute,coord_to_compute
        mask = radius == 0
        coord = coord[~mask]
        radius = radius[~mask]
        new_coord, new_radius = self.delete_voids(tomographic_map,radius,coord)
        new_coord_Mpc = self.convert_to_Mpc(tomographic_map,new_coord,new_radius)
        other_arrays,other_array_names=[],[]
        del indice,coord,map_3D,mask,radius
        self.add_Report("End of the Simple spherical finding for the map {}".format(tomographic_map.name))
        return(new_radius, new_coord_Mpc,other_arrays,other_array_names)


    def find_the_sphere(self,number_Mpc_per_pixels,number_pixel_maximal_radius,index):
        map_local = map_3D[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]
        indice_local = (indice[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]- index)*number_Mpc_per_pixels
        distance_map = np.sqrt(indice_local[:,:,:,0]**2 + indice_local[:,:,:,1]**2 + indice_local[:,:,:,2]**2)
        del indice_local
        rayon = self.params_void_finder["minimal_radius"]
        boolean = np.mean(map_local[distance_map < rayon]) > self.params_void_finder["average"]
        while(boolean&(rayon<self.params_void_finder["maximal_radius"])):
            rayon = rayon + self.params_void_finder["radius_step"]
            boolean = np.mean(map_local[distance_map < rayon]) > self.params_void_finder["average"]
        del distance_map,map_local,boolean
        if((rayon>=self.params_void_finder["maximal_radius"])|(rayon<=self.params_void_finder["minimal_radius"])):return(0)
        return(rayon)


    def find_the_sphere_cluster(self,number_Mpc_per_pixels,number_pixel_maximal_radius,index):
        map_local = map_3D[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]
        indice_local = (indice[max(index[0]-number_pixel_maximal_radius[0],0):min(map_3D.shape[0],index[0]+number_pixel_maximal_radius[0]),max(index[1]-number_pixel_maximal_radius[1],0):min(map_3D.shape[1],index[1]+number_pixel_maximal_radius[1]),max(index[2]-number_pixel_maximal_radius[2],0):min(map_3D.shape[2],index[2]+number_pixel_maximal_radius[2])]- index)*number_Mpc_per_pixels
        distance_map = np.sqrt(indice_local[:,:,:,0]**2 + indice_local[:,:,:,1]**2 + indice_local[:,:,:,2]**2)
        del indice_local
        rayon = self.params_void_finder["minimal_radius"]
        boolean = np.mean(map_local[distance_map < rayon]) < self.params_void_finder["average"]
        while(boolean&(rayon<self.params_void_finder["maximal_radius"])):
            rayon = rayon + self.params_void_finder["radius_step"]
            boolean = np.mean(map_local[distance_map < rayon]) < self.params_void_finder["average"]
        del distance_map,map_local,boolean
        if((rayon>=self.params_void_finder["maximal_radius"])|(rayon<=self.params_void_finder["minimal_radius"])):return(0)
        return(rayon)



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
        if other_arrays is not None:
            others_arrays_copies =[]
            for i in range(len(other_arrays)):
                others_arrays_copies.append(other_arrays[i].copy())
        else:
            others_arrays_copies =None
        radius_copy = radius.copy()
        coord_copy = coord.copy()
        nb_voids_delete = len(coord)
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
                index = np.array([int(coord_line.split()[0]),int(coord_line.split()[1]),int(coord_line.split()[2])])
                mask = coord == index
                arg = np.argwhere(mask)
                del mask
                radius[arg]=float(line[10].split(",")[0])
        del line, line_file
        self.add_Report("restart coordinates")
        mask = radius >= 0
        self.add_Report("Sphere found : R = " + str(radius[mask]) + ", Coord = " + str(coord[mask]))        
        return(radius)
        


    def convert_to_Mpc(self,tomographic_map,new_coord,new_radius):
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        if(new_coord.shape != 0):
            coord_Mpc = np.zeros(new_coord.shape)
            coord_Mpc = new_coord * np.array(number_Mpc_per_pixels)
        else:
            coord_Mpc = new_coord
        return(coord_Mpc)



    def save_voids(self,radius,coord,other_arrays,other_array_names):
        dict_void = {"R" : radius, "COORD" : coord}
        for i in range(len(other_arrays)):
            dict_void[other_array_names[i]] = other_arrays[i]
        if (self.find_cluster):
            name_out= "Clusters"
        else:
            name_out= "Voids"
        if(self.params_void_finder["method"]=="SPHERICAL"):
            name = "Dictionary_{}_{}_{}threshold_{}average_{}rmin_{}_deletion".format(name_out,self.params_void_finder["method"],self.params_void_finder["threshold"],self.params_void_finder["average"],self.params_void_finder["minimal_radius"],self.delete_option)        
        elif(self.params_void_finder["method"]=="WATERSHED"):
            name = "Dictionary_{}_{}_{}threshold_{}dist_clusters_{}rmin_{}_deletion".format(name_out,self.params_void_finder["method"],self.params_void_finder["threshold"],self.params_void_finder["dist_clusters"],self.params_void_finder["minimal_radius"],self.delete_option)        
        else :
            raise ValueError("The method_void chosen is not implemented, try : WATERSHED or SPHERICAL")  
        void = tomographic_objects.VoidCatalog.init_from_dictionary(f"{name}.fits",radius,coord,"cartesian",other_arrays=other_arrays,other_array_names = other_array_names)
        void.write()



        







class VoidPlot(object):
    
    def __init__(self,pwd):
        self.pwd = pwd

    

    
    def create_void_histogram_multiple(self,radius,nb_bins,rmin,rmax,name,norm=False):
        plt.figure()
        bins = np.linspace(rmin,rmax, nb_bins)
        (n, bins, patches) = plt.hist(radius, bins, alpha=1, label='x',density=norm)
        n = np.mean(n,axis=0)
        pickle.dump([n,bins],open(name,'wb'))

    
    def plot_radius_vs_redshift(self,name_dict,name):
        dict_void = self.load_dictionary(name_dict)
        radius = dict_void["radius"]
        z = dict_void["coord"][:,2]
        plt.figure()
        plt.plot(z,radius,"b.")
        plt.xlabel("Redshift")
        plt.ylabel("Radius size")
        plt.grid()
        plt.savefig(name + ".pdf", format ="pdf")
        plt.close()

    def plot_meanradius_vs_redshift(self,name_dict,name,nb_bins):
        dict_void = self.load_dictionary(name_dict)
        radius = dict_void["radius"]
        z = dict_void["coord"][:,2]
        maxredshift = np.max(z)
        minredshift = np.min(z)
        redshift = np.linspace(minredshift,maxredshift,nb_bins)
        mean_radius = []
        for i in range(nb_bins):
            mask = (z < minredshift + (i+1/nb_bins)*((maxredshift - minredshift))) & (z > minredshift + (i/nb_bins)*((maxredshift - minredshift)))
            mean_radius.append(np.mean(radius[mask]))
        plt.figure()
        plt.plot(redshift,mean_radius,"b.")
        plt.xlabel("Redshift")
        plt.ylabel("Radius size")
        plt.grid()
        plt.savefig(name + ".pdf", format ="pdf")
        plt.close()


    def plot_histo_redshift(self,name_dict,name,nb_bins):
        dict_void = self.load_dictionary(name_dict)
        z = dict_void["coord"][:,2]
        plt.figure()
        plt.hist(z,nb_bins)
        plt.xlabel("Redshift")
        plt.ylabel("Number of voids")
        plt.grid()
        plt.savefig(name + ".pdf", format ="pdf")
        plt.close()

    def plot_void_histogram(self,radius,nb_bins,name,norm=False):
        plt.figure()
        bins = np.linspace(np.min(radius),np.max(radius), nb_bins)
        plt.hist(radius, bins, alpha=1, label='x',density=norm)
        plt.grid()
        plt.ylabel("Number of voids")
        plt.xlabel("Radius of the void in Mpc.h-1")
        plt.savefig("histogram_voids_radius_"+ name + ".pdf",format="pdf")

        

    def load_and_plot_comparison(self,name1,name2,nb_bins,name,legend,rmin,rmax,norm=False,hist1=False,hist2=False,log_scale=True,supplementary_void=None,factor_add=None,legend_supp=None,expo_fit=False):
        bins1,bins2 = False,False
        if(not hist1):
            (radius1, coord1,filling_factor1) = self.load_void_dictionary(name1)
        else :
            radius1 = pickle.load(open(name1,'rb'))[0]
            bins1 = pickle.load(open(name1,'rb'))[1]
            assert(len(bins1) == nb_bins)
        if(not hist2):
            (radius2, coord2,filling_factor2) = self.load_void_dictionary(name2)
        else :
            radius2 = pickle.load(open(name2,'rb'))[0]
            bins2 = pickle.load(open(name2,'rb'))[1]
            assert(len(bins2) == nb_bins)
        if (supplementary_void is not None):
            radius_supp = []
            for i in range(len(supplementary_void)):
                (radius, coord,filling_factor) = self.load_void_dictionary(supplementary_void[i])
                if(factor_add is not None):
                    radius = list(radius)*factor_add[i]
                radius_supp.append(radius)
        else :
            radius_supp = None
        KS_stat, p_value = ks_2samp(radius1,radius2)
        self.plot_void_histogram_comparison(radius1,radius2,nb_bins,name,legend,rmin,rmax,norm=norm,hist1=hist1,hist2=hist2,log_scale=log_scale,radius_supp=radius_supp,legend_supp=legend_supp,expo_fit=expo_fit)
        return(KS_stat, p_value)

    def load_and_plot(self,name_dict,nb_bins,name,norm=False):
        (radius, coord,filling_factor) = self.load_void_dictionary(name_dict)
        self.plot_void_histogram(radius,nb_bins,name,norm=norm)


    def expo(self,x,a,b):
        return(np.exp(a*x+b))
        
    def plot_void_histogram_comparison(self,radius1,radius2,nb_bins,name,legend,rmin,rmax,norm=False,hist1=False,hist2=False,log_scale=True,radius_supp=None,legend_supp=None,expo_fit=False):
        plt.figure()
        bins = np.linspace(rmin,rmax, nb_bins)
        if(not hist1):
            (n1, bins1, patches1)  = plt.hist(radius1, bins, alpha=0.5, label=legend[0],density=norm)
        else :            
            (n1, bins1, patches1)  = plt.hist(bins[:-1],bins,weights=radius1, alpha=0.5, label=legend[0],density=norm)
        if(not hist2):
            (n2, bins2, patches2)  = plt.hist(radius2, bins, alpha=0.5, label=legend[1],density=norm)
        else :
            (n2, bins2, patches2)  = plt.hist(bins[:-1],weights=radius2, alpha=0.5, label=legend[0],density=norm)
        if(radius_supp is not None):
            for i in range(len(radius_supp)):
                (n, bins, patches)  = plt.hist(radius_supp[i], bins, label=legend_supp[i],histtype='step',linestyle='dashed',ec="k")
        if(expo_fit):
            from scipy.optimize import curve_fit
            bin_center = (bins1[1:] + bins1[0:-1]) / 2
            mask = bin_center>14
            n1,n2,bins_fit = n1[mask],n2[mask],bin_center[mask]
            fit = curve_fit(self.expo,bins_fit,n1)
            fit2 = curve_fit(self.expo,bins_fit,n2)
            plt.plot(bins_fit,self.expo(bins_fit,*fit[0]),'r')
            plt.plot(bins_fit,self.expo(bins_fit,*fit2[0]),'b')    
            perr = np.sqrt(np.diag(fit[1]))
            perr2 = np.sqrt(np.diag(fit2[1]))
            print(perr,perr2)
        if(log_scale):
            plt.yscale("log")
        plt.grid()
        if(legend_supp is not None):
            plt.legend(legend + legend_supp)
        else :
            plt.legend(legend )
        plt.xlabel("Void radius in Mpc.h" + r"$^{-1}$")
        if(log_scale):
            plt.savefig("histogram_voids_radius_"+ name + "log.pdf",format="pdf")
        else :
            plt.savefig("histogram_voids_radius_"+ name + "lin.pdf",format="pdf")


    def plot_void_histogram_subset(self,radius,coord_Mpc,n,name,mapsize):
        plt.figure()
        sub_interval = mapsize[0]/n
        intervales = [[i*sub_interval,(i+1)*sub_interval] for i in range(n)]
        legend = ["Ra22-33","Ra33-44"]
        color = ['red','blue']
        for i in range(len(intervales)) :
            mask = (coord_Mpc[:,0] > intervales[i][0]) &  (coord_Mpc[:,0] < intervales[i][1])
            plt.hist(radius[mask],20,fill=False, edgecolor = color[i])
            plt.grid()
            plt.title("Histogram of the void radius")
        plt.legend(legend)
        plt.savefig("histogram_voids_radius_"+ name +".pdf",format="pdf")


    def write_qso_to_mask(self,cross,nameout):
        file = open(nameout,"w")
        file.write("#Proto-cluster number & RA QSO (deg) & DEC QSO (deg) & redshift QSO\n")
        for i in range(len(cross)):
            for j in range(len(cross[i])):
                file.write("{} , {} , {} , {}\n".format(i+1,cross[i][j][0],cross[i][j][1],cross[i][j][2]))
        file.close()




