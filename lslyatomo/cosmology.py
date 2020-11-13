#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

Author: Corentin Ravoux

Description : Classes to return a Dachshund input based on data that can be
treated via picca software.
"""



#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################



import os,pickle,glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic
from random import sample
from matplotlib.lines import Line2D
from lslyatomo import utils
from lslyatomo import tomographic_objects




#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################





def get_deltas(namefile,center_ra=True):
    """ Extract delta properties """
    ra, dec, z, deltas, sigmas,zqsos,ids = [],[],[],[],[],[],[]
    for i in range(len(namefile)):
        delta_tomo = tomographic_objects.Delta(name=namefile[i],pk1d_type=True)
        delta_tomo.read()
        for j in range(1,len(delta_tomo.delta_file)):
            delta = delta_tomo.read_line(j)
            zqsos.append(delta.zqso)
            if((delta.ra > np.pi)&(center_ra)):ra.append(delta.ra-2*np.pi)
            else : ra.append(delta.ra)
            dec.append(delta.dec)
            if(delta.iv is not None): sigmas.append(1/np.sqrt(np.asarray(delta.iv)))
            else : sigmas.append(np.array([0.0 for i in range(len(delta.de))]))
            deltas.append(delta.de)
            ids.append(delta.thid)
            z.append(((10**delta.ll / utils.lambdaLy)-1))
        delta_tomo.close()
    return(np.array(ra),np.array(dec),np.asarray(z),np.array(zqsos),np.array(ids),np.asarray(sigmas),np.asarray(deltas))


def get_merged_multiple_exposure_deltas(namefile):
    """ Merge deltas with repeated observation"""
    # Pack LOS by Id in the dict Deltas
    ra, dec, z, deltas, sigmas,zqso,Idlist = [],[],[],[],[],[],[],[]
    for i in range(len(namefile)):
        delta_tomo = tomographic_objects.Delta(name=namefile[i],pk1d_type=True)
        delta_tomo.read()
        for j in range(1,len(delta_tomo.delta_file)):
            Delta = delta_tomo.read_line(j)
            Idlist.append(Delta.thid)
        delta_tomo.close()
    Idlist = list(set(Idlist))
    Deltas = {}
    for i in range(len(Idlist)):
        Deltas[Idlist[i]] = []
        for j in range(len(namefile)):
            delta_tomo = tomographic_objects.Delta(name=namefile[j],pk1d_type=True)
            delta_tomo.read()
            for k in range(1,len(delta_tomo.delta_file)):
                Delta = delta_tomo.read_line(k)
                if(Delta.thid == Idlist[i]):
                    Deltas[Idlist[i]].append(Delta)
            delta_tomo.close()
    # For each pack of LOS
    for i in range(len(Idlist)):

        # Get the data
        zqso.append(Deltas[Idlist[i]][0].zqso)
        ra.append(Deltas[Idlist[i]][0].ra)
        dec.append(Deltas[Idlist[i]][0].dec)
        listsigmas,listz,listdelta = [],[],[]
        for j in range(len(Deltas[Idlist[i]])):
            listsigmas.append(1/np.sqrt(np.asarray(Deltas[Idlist[i]][j].iv)))
            listdelta.append(Deltas[Idlist[i]][j].de)
            listz.append(((10**np.asarray(Deltas[Idlist[i]][j].ll) / utils.lambdaLy)-1))

        # Get the list of common elements in the pack to have minimum and maximum redshifts
        zmin,zmax = 0, 10**10
        zcommon = listz[0]
        for j in range(1,len(listz)):
            zcommon = list(set(zcommon).intersection(listz[j]))
        zmin = np.min(zcommon)
        zmax = np.max(zcommon)

        # Deleting pixels at the beginning and at the end of LOSs
        for j in range(len(listz)):
            lineToDeleteFirst = 0
            k = 0
            while(listz[j][k] != zmin):
                lineToDeleteFirst = lineToDeleteFirst + 1
                k = k + 1
            lineToDeleteLast = 0
            k = -1
            while(listz[j][k] != zmax):
                lineToDeleteLast = lineToDeleteLast + 1
                k = k -1
            listz[j] = listz[j][lineToDeleteFirst:len(listz[j])-lineToDeleteLast]
            listsigmas[j] = listsigmas[j][lineToDeleteFirst:len(listsigmas[j])-lineToDeleteLast]
            listdelta[j] = listdelta[j][lineToDeleteFirst:len(listdelta[j])-lineToDeleteLast]

        # Ensuring that all LOS have the same lenght
        lenlists = []
        for j in range(len(listz)):
            lenlists.append(len(listz[j]))

        # Selection of the pixels to delete in case of a missing pixel along one LOS + Deletion
        while(np.max(lenlists)!=np.min(lenlists)):
            eltTodelete = [[] for j in range(len(listz))]
            mi = 10**10
            mins = []
            for j in range(len(lenlists)):
                if(lenlists[j] < mi):
                    mi = lenlists[j]
                    mins =[j]
                elif (lenlists[j] == mi):
                    mins.append(j)
            for j in range(len(mins)):
                for k in range(len(listz)):
                    for m in range(len(listz[k])):
                        if((k!=mins[j])&(np.isin([listz[k][m]],listz[mins[j]]) == False)):
                            eltTodelete[k].append(m)
            newlistz, newlistsigma, newlistdelta = [[] for n in range(len(listz))],[[] for n in range(len(listz))],[[] for n in range(len(listz))]
            for j in range(len(eltTodelete)):
                eltTodeletej = list(set(eltTodelete[j]))
                for k in range(len(listz[j])):
                    if (np.isin([k], eltTodeletej)==False):
                        newlistz[j].append(listz[j][k])
                        newlistsigma[j].append(listsigmas[j][k])
                        newlistdelta[j].append(listdelta[j][k])
            listz = newlistz
            listsigmas = newlistsigma
            listdelta = newlistdelta
            lenlists = []
            for j in range(len(listz)):
                lenlists.append(len(listz[j]))


        # Weighted merging of the LOSs
        zmerged, sigmamerged, deltamerged = listz[0],[],[]
        for k in range(len(listz[0])):
            sigma = 0
            delta = 0
            sumsigma = 0
            for m in range(len(listz)):
                sumsigma = sumsigma + 1/(listsigmas[m][k]**2)
            for m in range(len(listz)):
                delta = delta + listdelta[m][k] / ((listsigmas[m][k]**2)*sumsigma)

                sigma = sigma + 1 / ((listsigmas[m][k]**2)*(sumsigma**2))
            sigma = np.sqrt(sigma/len(listz))
            sigmamerged.append(sigma)
            deltamerged.append(delta)
        deltas.append(deltamerged)
        z.append(zmerged)
        sigmas.append(sigmamerged)
    return(np.array(ra),np.array(dec),np.asarray(z),np.array(zqso),np.array(Idlist),np.asarray(sigmas),np.asarray(deltas))











#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################








class DeltaConverter():

    def __init__(self,pwd,Omega_m,delta_path,coordinate_transform,plot_pixel_properties,software,return_qso_catalog=None,return_dla_catalog=None,dla_catalog=None,return_sky_catalogs=False,repeat=False):

        self.pwd = pwd
        self.delta_path = delta_path
        self.Omega_m = Omega_m
        self.coordinate_transform = coordinate_transform
        self.plot_pixel_properties = plot_pixel_properties
        self.return_qso_catalog = return_qso_catalog
        self.return_dla_catalog = return_dla_catalog
        self.dla_catalog = dla_catalog
        self.return_sky_catalogs = return_sky_catalogs
        self.repeat = repeat
        self.software = software



    def transform_delta_to_pixel_file(self,rebin=None,shuffle=None,sigma_min=None,sigma_max=None,z_cut_min=None,z_cut_max=None,dec_cut_min=None,dec_cut_max=None,ra_cut_min=None,ra_cut_max=None):
        namefile = glob.glob(os.path.join(self.delta_path,"delta-*"))
        properties_map_pixels = {}
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(self.Omega_m)
        if(self.repeat):
            (ra,dec,z,zqso,ids,sigmas,deltas) = get_merged_multiple_exposure_deltas(namefile)
        else:
            (ra,dec,z,zqso,ids,sigmas,deltas)  = get_deltas(namefile)

        if(shuffle is not None):
            (ra,dec,deltas,sigmas) = self.shuffle_data(shuffle,ra,dec,deltas,sigmas)

        if(rebin is not None):
            (z,deltas,sigmas) = self.rebin_data(z,deltas,sigmas,rebin)


        if(self.return_dla_catalog is not None):
            zdlas,z_qso_dlas = [],[]
            if(self.dla_catalog is None):
                raise KeyError("Please give a DLA catalog name or turn off the return_dla_catalog option")
            dla_catalog = tomographic_objects.DLACatalog.init_from_fits(self.dla_catalog)
            for i in range(len(ids)):
                mask = (dla_catalog.primary_key == ids[i])
                if(z_cut_min is not None):
                    mask &=(dla_catalog.coord_z > z_cut_min)
                if(z_cut_max is not None):
                    mask &=(dla_catalog.coord_z < z_cut_max)
                zdlas.append(dla_catalog.coord_z[mask])
                z_qso_dlas.append(dla_catalog.z_qso[mask])




        sky_deltas = np.array([[ra[i],dec[i],z[i][j],sigmas[i][j],deltas[i][j]] for i in range(len(ra)) for j in range(len(z[i]))])
        sky_deltas = sky_deltas[utils.cut_sky_catalog(sky_deltas[:,0],sky_deltas[:,1],sky_deltas[:,2],ramin=ra_cut_min,ramax=ra_cut_max,decmin=dec_cut_min,decmax=dec_cut_max,zmin=z_cut_min,zmax=z_cut_max)]

        suplementary_parameters = utils.return_suplementary_parameters(self.coordinate_transform,zmin=np.min(sky_deltas[:,2]),zmax=np.max(sky_deltas[:,2]))
        cartesian_deltas = np.zeros(sky_deltas.shape)
        cartesian_deltas[:,0],cartesian_deltas[:,1],cartesian_deltas[:,2] = utils.convert_sky_to_cartesian(sky_deltas[:,0],sky_deltas[:,1],sky_deltas[:,2],self.coordinate_transform,rcomov=rcomov,distang=distang,suplementary_parameters=suplementary_parameters)
        cartesian_deltas[:,3],cartesian_deltas[:,4] = sky_deltas[:,3],sky_deltas[:,4]
        properties_map_pixels["minx"],properties_map_pixels["miny"],properties_map_pixels["minz"] = np.min(cartesian_deltas[:,0]),np.min(cartesian_deltas[:,1]),np.min(cartesian_deltas[:,2])
        cartesian_deltas = cartesian_deltas - np.array([properties_map_pixels["minx"],properties_map_pixels["miny"],properties_map_pixels["minz"],0,0])


        if(sigma_min is not None):
            cartesian_deltas[:,3][cartesian_deltas[:,3] < sigma_min] = sigma_min
        if(sigma_max is not None):
            cartesian_deltas[:,3][cartesian_deltas[:,3] > sigma_max] = sigma_max


        if(self.return_qso_catalog is not None):
            sky_qso_catalog = np.array([[ra[i],dec[i],zqso[i],ids[i]] for i in range(len(ra))])
            sky_qso_catalog= sky_qso_catalog[utils.cut_sky_catalog(sky_qso_catalog[:,0],sky_qso_catalog[:,1],sky_qso_catalog[:,2],ramin=ra_cut_min,ramax=ra_cut_max,decmin=dec_cut_min,decmax=dec_cut_max,zmin=z_cut_min,zmax=z_cut_max),:]
            cartesian_qso_catalog = np.zeros(sky_qso_catalog.shape)
            cartesian_qso_catalog[:,0],cartesian_qso_catalog[:,1],cartesian_qso_catalog[:,2] = utils.convert_sky_to_cartesian(sky_qso_catalog[:,0],sky_qso_catalog[:,1],sky_qso_catalog[:,2],self.coordinate_transform,rcomov=rcomov,distang=distang,suplementary_parameters=suplementary_parameters)
            cartesian_qso_catalog[:,3] = sky_qso_catalog[:,3]
            cartesian_qso_catalog = cartesian_qso_catalog - np.array([properties_map_pixels["minx"],properties_map_pixels["miny"],properties_map_pixels["minz"],0])

        else:
            sky_qso_catalog,cartesian_qso_catalog=None,None

        if(self.return_dla_catalog is not None):
            sky_dla_catalog = np.array([[ra[i],dec[i],zdlas[i][j],z_qso_dlas[i][j]] for i in range(len(ra)) for j in range(len(zdlas[i]))])
            sky_dla_catalog = sky_dla_catalog[utils.cut_sky_catalog(sky_dla_catalog[:,0],sky_dla_catalog[:,1],sky_dla_catalog[:,2],ramin=ra_cut_min,ramax=ra_cut_max,decmin=dec_cut_min,decmax=dec_cut_max,zmin=z_cut_min,zmax=z_cut_max)]
            cartesian_dla_catalog = np.zeros(sky_dla_catalog.shape)
            cartesian_dla_catalog[:,0],cartesian_dla_catalog[:,1],cartesian_dla_catalog[:,2] = utils.convert_sky_to_cartesian(sky_dla_catalog[:,0],sky_dla_catalog[:,1],sky_dla_catalog[:,2],self.coordinate_transform,rcomov=rcomov,distang=distang,suplementary_parameters=suplementary_parameters)
            cartesian_qso_catalog[:,3] = sky_qso_catalog[:,3]
            cartesian_dla_catalog = cartesian_dla_catalog - np.array([properties_map_pixels["minx"],properties_map_pixels["miny"],properties_map_pixels["minz"],0])
        else:
            sky_dla_catalog,cartesian_dla_catalog=None,None


        return(cartesian_deltas,cartesian_qso_catalog,cartesian_dla_catalog,sky_deltas,sky_qso_catalog,sky_dla_catalog,properties_map_pixels)



    def shuffle_data(self,shuffle,ra,dec,deltas,sigmas):
        if(shuffle == "radec"):
            ra,dec = self.shuffle_arrays(ra,dec)
        elif(shuffle == "deltasigma") :
            deltas,sigmas = self.shuffle_deltas_sigmas(deltas,sigmas)
        return(ra,dec,deltas,sigmas)


    def shuffle_deltas_sigmas(self,deltas,sigmas):
        len_deltas = [len(deltas[i]) for i in range(len(deltas))]
        delta_list, sigma_list = [],[]
        for i in range(len(deltas)):
            delta_list = delta_list + [deltas[i][j] for j in range(len(deltas[i]))]
            sigma_list = sigma_list + [sigmas[i][j] for j in range(len(sigmas[i]))]
        delta_list ,sigma_list = self.shuffle_arrays(delta_list,sigma_list)
        deltas_shuffle,sigmas_shuffle = [],[]
        for i in range(len(len_deltas)):
            i_begin = int(np.sum(len_deltas[:(i)]))
            i_end = np.sum(len_deltas[:(i+1)])
            deltas_shuffle.append(delta_list[i_begin:i_end])
            sigmas_shuffle.append(sigma_list[i_begin:i_end])
        return(deltas_shuffle,sigmas_shuffle)


    def shuffle_arrays(self,x,y):
        return(np.random.permutation(x), np.random.permutation(y))


    def rebin_data(self,z,deltas,sigmas,bin_pixel,method="gauss"):
        for i in range(len(z)):
            if(len(z[i])>1):
                if(len(z[i])<=bin_pixel):
                    z[i]=[np.mean(z[i])]
                    deltas[i] = [np.mean(deltas[i])]
                    sigmas[i] = [np.mean(sigmas[i])]
                else :
                    new_shape = len(z[i])//bin_pixel
                    first_coord = len(z[i]) - (len(z[i])//bin_pixel)*bin_pixel
                    if(first_coord==0):
                        z[i] = utils.bin_ndarray(np.array(z[i])[:],[new_shape],operation=method)
                        deltas[i] = utils.bin_ndarray(np.array(deltas[i])[:],[new_shape],operation=method)
                        sigmas[i] = utils.bin_ndarray(np.array(sigmas[i])[:],[new_shape],operation=method)
                    else :
                        z[i] = np.concatenate([[np.mean(np.array(z[i])[:first_coord])],utils.bin_ndarray(np.array(z[i])[first_coord:],[new_shape],operation=method)])
                        deltas[i] = np.concatenate([[np.mean(np.array(deltas[i])[:first_coord])],utils.bin_ndarray(np.array(deltas[i])[first_coord:],[new_shape],operation=method)])
                        sigmas[i] = np.concatenate([[np.mean(np.array(sigmas[i])[:first_coord])],utils.bin_ndarray(np.array(sigmas[i])[first_coord:],[new_shape],operation=method)])
        return(z,deltas,sigmas)



    def create_input_files(self,coordinates_to_write,properties,name_pixel,sky_coordinates=None,create_launcher=None):
        if(self.software.lower() == "dachshund"):
            self.create_dachshund_input_files(coordinates_to_write,properties,name_pixel,sky_coordinates=sky_coordinates,create_launcher=create_launcher)


    def create_dachshund_input_files(self,coordinates_to_write,properties,name_pixel,sky_coordinates=None,create_launcher=None):
        pixel = tomographic_objects.Pixel(name=name_pixel,pixel_array=coordinates_to_write)
        pixel.write()
        if(create_launcher is not None):
            self.create_dachshund_launcher(np.max(coordinates_to_write[:,0]),np.max(coordinates_to_write[:,1]),np.max(coordinates_to_write[:,2]),len(coordinates_to_write),properties["shape"][0],properties["shape"][1],properties["shape"][2],properties["sigma_f"],properties["lperp"],properties["lpar"],properties["name_pixel"],properties["name_map"],create_launcher)



    def create_dachshund_launcher(self,lx,ly,lz,npix,nx,ny,nz,sigmaf,lperp,lpar,namepixel,namemap,nameinput):
        f = open(f'{nameinput}.cfg',"w")
        f.write("#lx, ly, lz: the domain size in each direction.\n")
        f.write("#num_pixels: the *total* number of pixels.\n")
        f.write("#map_nx, map_ny, map_nz: the number of map points. The map points are arbitrary but for now these n's are used to setup a uniform grid across the domain given above.\n")
        f.write("#corr_var_s: the signal cov prefactor sigma_f^2\n")
        f.write("#corr_l_perp: the signal cov perp scale.\n")
        f.write("#corr_l_para: the signal cov para scale.\n")
        f.write("#pcg_max_iter: the PCG max number of iterations. 100 should be good.\n")
        f.write("#pcg_tol: the PCG stopping tolerance. I found 1.0e-3 is good enough. Set it very small if you want the most accurate map.\n")
        f.write("lx = {}\n".format(lx))
        f.write("ly = {}\n".format(ly))
        f.write("lz = {}\n".format(lz))
        f.write("\n")
        f.write("# From output of GEN_DACH_INPUT.PRO\n")
        f.write("num_pixels = {}\n".format(npix))
        f.write("\n")
        f.write("map_nx = {}\n".format(nx))
        f.write("map_ny = {}\n".format(ny))
        f.write("map_nz = {}\n".format(nz))
        f.write("\n")
        f.write("corr_var_s = {}\n".format(sigmaf))
        f.write("corr_l_perp = {}\n".format(lperp))
        f.write("corr_l_para = {}\n".format(lpar))
        f.write("\n")
        f.write("pcg_max_iter = 500\n")
        f.write("pcg_tol = 1.0e-3\n")
        f.write("#pcg_step_r = 1\n")
        f.write("\n")
        f.write("option_map_covar = 0\n")
        f.write("option_noise_covar = 0\n")
        f.write("pixel_data_path = {}\n".format(namepixel))
        f.write("map_path = {}\n".format(namemap))
        f.close()


    def create_dachshund_map_pixel_property_file(self,name_out,cartesian_coordinates,sky_coordinates,shape,properties_map_pixels):
        size = (np.max(cartesian_coordinates[:,0]),np.max(cartesian_coordinates[:,1]),np.max(cartesian_coordinates[:,2]))
        coordinate_transform = self.coordinate_transform
        boundary_cartesian_coord = ((properties_map_pixels["minx"],properties_map_pixels["miny"],properties_map_pixels["minz"]),(properties_map_pixels["minx"]+size[0],properties_map_pixels["miny"]+size[1],properties_map_pixels["minz"]+size[2]))
        boundary_sky_coord = ((np.min(sky_coordinates[:,0]),np.min(sky_coordinates[:,1]),np.min(sky_coordinates[:,2])),(np.max(sky_coordinates[:,0]),np.max(sky_coordinates[:,1]),np.max(sky_coordinates[:,2])))
        property_file = tomographic_objects.MapPixelProperty(name=name_out,size=size,shape=shape,boundary_cartesian_coord=boundary_cartesian_coord,boundary_sky_coord=boundary_sky_coord,coordinate_transform=coordinate_transform,Omega_m=self.Omega_m)
        property_file.write()



    def create_serial_input(self,nameout,properties,cartesian_deltas,sky_deltas,property_file_name,properties_map_pixels):
        self.create_input_files(cartesian_deltas,properties,properties["name_pixel"],sky_deltas,create_launcher=nameout)
        if(self.return_sky_catalogs):
            self.create_input_files(sky_deltas,properties,"{}_sky_coordinates".format(properties["name_pixel"]),sky_deltas,create_launcher=None)
        self.create_dachshund_map_pixel_property_file(property_file_name,cartesian_deltas,sky_deltas,properties["shape"],properties_map_pixels)




    def cut_in_chunks(self,cartesian_deltas,number_chunks,overlaping,shape_sub_map):
        if (overlaping is None) :
            overlaping = 0.0
        minx,maxx,miny,maxy= np.min(cartesian_deltas[:,0]),np.max(cartesian_deltas[:,0]),np.min(cartesian_deltas[:,1]),np.max(cartesian_deltas[:,1])
        intervalx = (maxx - minx)
        intervaly = (maxy - miny)
        subIntervalx = intervalx/number_chunks[0]
        subIntervaly = intervaly/number_chunks[1]
        Chunks = {}
        shape_x =  number_chunks[0]*shape_sub_map[0]
        shape_y =  number_chunks[1]*shape_sub_map[1]
        remove_shape_x,remove_shape_y = 0,0
        for i in range(number_chunks[0]):
            for j in range(number_chunks[1]):
                filename = f'{i:03d}' + f'{j:03d}'
                Chunks[filename]={}
                if((i==number_chunks[0]-1)&(i==0)):
                    intervalxChunk = [i*subIntervalx, (i+1)*subIntervalx]
                elif i == 0 :
                    intervalxChunk = [i*subIntervalx, (i+1)*subIntervalx + overlaping]
                elif i == number_chunks[0]-1 :
                    intervalxChunk = [i*subIntervalx - overlaping, intervalx]
                else:
                    intervalxChunk = [i*subIntervalx - overlaping, (i+1)*subIntervalx + overlaping]
                if((j==number_chunks[1]-1)&(j==0)):
                    intervalyChunk = [  j*subIntervaly, (j+1)*subIntervaly]
                elif j == 0 :
                    intervalyChunk = [  j*subIntervaly, (j+1)*subIntervaly + overlaping]
                elif j == number_chunks[1]-1 :
                    intervalyChunk = [ j*subIntervaly - overlaping , intervaly]
                else:
                    intervalyChunk = [  j*subIntervaly -overlaping, (j+1)*subIntervaly + overlaping]
                mask = (cartesian_deltas[:,0] < intervalxChunk[1])&(cartesian_deltas[:,0] >= intervalxChunk[0])
                mask &= (cartesian_deltas[:,1] < intervalyChunk[1])&(cartesian_deltas[:,1] >= intervalyChunk[0])
                chunks_deltas = []
                chunks_deltas.append(cartesian_deltas[:,0][mask] - intervalxChunk[0])
                chunks_deltas.append(cartesian_deltas[:,1][mask] - intervalyChunk[0])
                chunks_deltas.append(cartesian_deltas[:,2][mask] )
                chunks_deltas.append(cartesian_deltas[:,3][mask] )
                chunks_deltas.append(cartesian_deltas[:,4][mask] )
                chunks_deltas = np.transpose(np.stack(chunks_deltas))
                Chunks[filename]["coord"] = chunks_deltas
                Chunks[filename]["limits"] = [intervalxChunk[0],intervalxChunk[1],intervalyChunk[0],intervalyChunk[1],np.min(cartesian_deltas[:,2]),np.max(cartesian_deltas[:,2])]
                nbPixelperMpcx = shape_sub_map[0]/(intervalxChunk[1] - intervalxChunk[0])
                nbPixelperMpcy = shape_sub_map[1]/(intervalyChunk[1] - intervalyChunk[0])
                nbpixelToremovex = int(round(nbPixelperMpcx * overlaping,0))
                nbpixelToremovey = int(round(nbPixelperMpcy * overlaping,0))
                if(number_chunks[0] !=1):
                    if ((i==0)|(i == number_chunks[0] - 1)):
                        remove_shape_x = remove_shape_x + nbpixelToremovex
                    else:
                        remove_shape_x = remove_shape_x + 2*nbpixelToremovex
                if(number_chunks[1] !=1):
                    if ((j==0)|(j == number_chunks[1] - 1)):
                        remove_shape_y = remove_shape_y + nbpixelToremovey
                    else:
                        remove_shape_y = remove_shape_y + 2*nbpixelToremovey
        shape_x = shape_x - remove_shape_x//number_chunks[1]
        shape_y = shape_y - remove_shape_y//number_chunks[0]
        Chunks["overlaping"]=overlaping
        return(Chunks,(shape_x,shape_y))



    def create_parallel_input(self,nameout,properties,cartesian_deltas,sky_deltas,number_chunks,overlaping,shape_sub_map,property_file_name,properties_map_pixels):
        self.create_input_files(cartesian_deltas,properties,properties["name_pixel"],sky_deltas,create_launcher=nameout)
        chunks ,shape= self.cut_in_chunks(cartesian_deltas,number_chunks,overlaping,shape_sub_map)
        self.create_dachshund_map_pixel_property_file(property_file_name,cartesian_deltas,sky_deltas,(shape[0],shape[1],shape_sub_map[2]),properties_map_pixels)
        filename = []
        Dachshundparams = []
        for i in range(len(list(chunks.keys()))):
            key = list(chunks.keys())[i]
            if key != 'overlaping' :
                filename.append(key)
                self.create_input_files(chunks[key]["coord"],properties,"{}_{}".format(properties["name_pixel"],key))
                Dachshundparams.append({})
                Dachshundparams[i]["maxx"]=chunks[key]["limits"][1]
                Dachshundparams[i]["minx"]=chunks[key]["limits"][0]
                Dachshundparams[i]["maxy"]=chunks[key]["limits"][3]
                Dachshundparams[i]["miny"]=chunks[key]["limits"][2]
                Dachshundparams[i]["maxz"]=chunks[key]["limits"][5]
                Dachshundparams[i]["minz"]=chunks[key]["limits"][4]
                Dachshundparams[i]["lx"]=chunks[key]["limits"][1] - chunks[key]["limits"][0]
                Dachshundparams[i]["ly"]=chunks[key]["limits"][3] - chunks[key]["limits"][2]
                Dachshundparams[i]["lz"]=chunks[key]["limits"][5] - chunks[key]["limits"][4]
                Dachshundparams[i]["npix"]=len(chunks[key]["coord"])
                Dachshundparams[i]["nx"]=shape_sub_map[0]
                Dachshundparams[i]["ny"]=shape_sub_map[1]
                Dachshundparams[i]["nz"]=shape_sub_map[2]
                Dachshundparams[i]["sigmaf"]=properties["sigma_f"]
                Dachshundparams[i]["lperp"]=properties["lperp"]
                Dachshundparams[i]["lpar"]=properties["lpar"]
                Dachshundparams[i]["namepixel"]="{}_{}".format(properties["name_pixel"],key)
                Dachshundparams[i]["namemap"]="map_{}_{}".format(properties["name_pixel"],key)
                Dachshundparams[i]["nameinput"]= "input_{}.cfg".format(key)
        pickle.dump([filename,Dachshundparams,number_chunks,overlaping],open("data_launch_dachshund.pickle","wb"))


    ### Modification : dissociate map properties and launcher properties

    #### Main routine





    def transform_delta(self,mode,nameout,properties,rebin=False,shuffle=None,sigma_min=None,sigma_max=None,z_cut_min=None,z_cut_max=None,dec_cut_min=None,dec_cut_max=None,ra_cut_min=None,ra_cut_max=None,number_chunks=None,overlaping=None,shape_sub_map=None,property_file_name="property_file.pickle"):
        (cartesian_deltas,cartesian_qso_catalog,cartesian_dla_catalog,sky_deltas,sky_qso_catalog,sky_dla_catalog,properties_map_pixels) = self.transform_delta_to_pixel_file(rebin=rebin,shuffle=shuffle,sigma_min=sigma_min,sigma_max=sigma_max,z_cut_min=z_cut_min,z_cut_max=z_cut_max,dec_cut_min=dec_cut_min,dec_cut_max=dec_cut_max,ra_cut_min=ra_cut_min,ra_cut_max=ra_cut_max)
        if(mode.lower() == "serial"):
            self.create_serial_input(nameout,properties,cartesian_deltas,sky_deltas,property_file_name,properties_map_pixels)
        elif(mode.lower() == "parallel"):
            self.create_parallel_input(nameout,properties,cartesian_deltas,sky_deltas,number_chunks,overlaping,shape_sub_map,property_file_name,properties_map_pixels)
        else:
            raise KeyError("Please choose a mode between serial and parallel")
        self.create_additional_catalogs(cartesian_qso_catalog,cartesian_dla_catalog,sky_qso_catalog,sky_dla_catalog)
        if(self.plot_pixel_properties):
            pixel = tomographic_objects.Pixel.init_from_property_files(property_file_name,pixel_array=cartesian_deltas,name=properties["name_pixel"])
            pixel_analyzer = PixelAnalizer(self.pwd,pixel=pixel)
            pixel_analyzer.analyze_pixels(False,True,name_dperp=nameout,coupled_plot=True)

    def create_additional_catalogs(self,cartesian_qso_catalog,cartesian_dla_catalog,sky_qso_catalog,sky_dla_catalog):
        if(self.return_dla_catalog is not None):
            dla_catalog = tomographic_objects.DLACatalog.init_from_pixel_catalog(cartesian_dla_catalog,name=self.return_dla_catalog)
            dla_catalog.write()
            if(self.return_sky_catalogs):
                dla_catalog = tomographic_objects.DLACatalog.init_from_pixel_catalog(sky_dla_catalog,name=f"{self.return_dla_catalog}_sky_coordinates")
                dla_catalog.write()
        if(self.return_qso_catalog is not None):
            quasar_catalog = tomographic_objects.QSOCatalog.init_from_pixel_catalog(cartesian_qso_catalog,name=self.return_qso_catalog)
            quasar_catalog.write()
            if(self.return_sky_catalogs):
                quasar_catalog = tomographic_objects.QSOCatalog.init_from_pixel_catalog(sky_qso_catalog,name=f"{self.return_qso_catalog}_sky_coordinates")
                quasar_catalog.write()





class PixelAnalizer(object):

    ### Include in pixel class instead of the coord


    def __init__(self,pwd,pixel=None):

        self.pwd = pwd
        self.pixel = pixel




    @staticmethod
    def write_density_file(z,dperp,name):
        pickle.dump([z,dperp],open(name,"wb"))

    @staticmethod
    def write_dperp_file(z,dperp,name):
        pickle.dump([z,dperp],open(name,"wb"))


    @staticmethod
    def read_density_file(name):
        a = pickle.dump(open(name,"wb"))
        return(a[0],a[1])

    @staticmethod
    def read_dperp_file(name):
        a = pickle.load(open(name,"wb"))
        return(a[0],a[1])




    @staticmethod
    def plot_histogram_mean_distance(zpar,dmin,name_histo,nb_bins=50):
        plt.figure()
        plt.hist(dmin,nb_bins)
        plt.xlabel("minimal distance histogram at Z={}".format(zpar))
        plt.savefig("{}_at_Z{}.pdf".format(name_histo,zpar),format="pdf")




    @staticmethod
    def plot_mean_distance_density(zpar,dperpz,densityz,nameout,coupled_plot=False,comparison=False,dperp_comparison=None,density_comparison=None,zpar_comparison=None,legend=None):
        if(coupled_plot is False):

            plt.figure()
            plt.xlabel("Redshift")
            plt.ylabel("Mean los separation [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]")
            plt.grid()
            plt.plot(zpar,dperpz)
            if(comparison):
                plt.plot(zpar_comparison,dperp_comparison)
                plt.legend(legend)
            plt.savefig("separation_{}.pdf".format(nameout), format = "pdf")

            plt.figure()
            plt.xlabel("Redshift")
            plt.ylabel("Density [" + r"$\mathrm{deg^{-2}}$" + "]")
            plt.grid()
            plt.plot(zpar,densityz)
            if(comparison):
                plt.plot(zpar_comparison,density_comparison)
                plt.legend(legend)
            plt.savefig("density_{}.pdf".format(nameout), format = "pdf")

        else:
            fig, ax1 = plt.subplots()

            color = 'tab:blue'
            ax1.set_xlabel("Redshift z")
            ax1.set_ylabel("Mean separation [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]", color=color)
            ax1.plot(zpar,dperpz, color=color)
            if(comparison):
                ax1.plot(zpar_comparison,dperp_comparison, color=color,linestyle ="--")
            ax1.tick_params(axis='y', labelcolor=color)

            ax2 = ax1.twinx()

            color = 'tab:orange'
            ax2.set_ylabel("Density [" + r"$\mathrm{deg^{-2}}$" + "]", color=color)
            ax2.plot(zpar,densityz, color=color)
            if(comparison):
                ax2.plot(zpar_comparison,density_comparison, color=color,linestyle ="--")
            ax2.tick_params(axis='y', labelcolor=color)

            fig.tight_layout()
            if(comparison):
                legend_elements = [Line2D([0], [0], color='k', lw=1, label=legend[0]),Line2D([0], [0], color='k',linestyle ="--", lw=1, label=legend[1])]
                ax1.legend(handles=legend_elements,loc = "upper center")
            plt.savefig("{}.pdf".format(nameout),format ="pdf",dpi=200)


### Improved on lsstomo :
        # if(coupled_plot is False):

        #     plt.figure()
        #     plt.plot(stripe_dperp[0],stripe_dperp[1])
        #     plt.plot(saclay_dperp[0],saclay_dperp[1])
        #     if(otherdperp is not None):
        #         for i in range(len(dperp_other)):
        #             plt.plot(dperp_other[i][0],dperp_other[i][1])
        #     plt.ylabel("Mean los separation [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]")
        #     plt.xlabel("Redshift")
        #     plt.legend(legend)
        #     plt.grid()
        #     plt.savefig("separation_{}.pdf".format(nameout), format = "pdf")


        #     plt.figure()
        #     plt.plot(stripe_density[0],stripe_density[1])
        #     plt.plot(saclay_density[0],saclay_density[1])
        #     if(otherdensity is not None):
        #         for i in range(len(otherdensity)):
        #             plt.plot(density_other[i][0],density_other[i][1])
        #     plt.ylabel("Density [" + r"$\mathrm{deg^{-2}}$" + "]")
        #     plt.xlabel("Redshift")
        #     plt.legend(legend)
        #     plt.grid()
        #     plt.savefig("density_{}.pdf".format(nameout), format = "pdf")

        # else :
        #     fig, ax1 = plt.subplots()
        #     line = ["dotted","dashdot","densely dashdotdotted"]
        #     color = 'tab:blue'
        #     ax1.set_xlabel("Redshift")
        #     ax1.set_ylabel("Mean separation [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]", color=color)
        #     ax1.plot(stripe_dperp[0],stripe_dperp[1], color=color)
        #     ax1.plot(saclay_dperp[0],saclay_dperp[1], color=color,linestyle ="--")
        #     if(otherdperp is not None):
        #         for i in range(len(dperp_other)):
        #             ax1.plot(dperp_other[i][0],dperp_other[i][1], color=color,linestyle =line[i])
        #     ax1.tick_params(axis='y', labelcolor=color)

        #     ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        #     color = 'tab:orange'
        #     ax2.set_ylabel("Density [" + r"$\mathrm{deg^{-2}}$" + "]", color=color)  # we already handled the x-label with ax1
        #     ax2.plot(stripe_density[0],stripe_density[1], color=color)
        #     ax2.plot(saclay_density[0],saclay_density[1], color=color,linestyle ="--")
        #     if(otherdensity is not None):
        #         for i in range(len(otherdensity)):
        #             ax2.plot(density_other[i][0],density_other[i][1], color=color,linestyle =line[i])
        #     ax2.tick_params(axis='y', labelcolor=color)

        #     fig.tight_layout()  # otherwise the right y-label is slightly clipped
        #     legend_elements = [Line2D([0], [0], color='k', lw=1, label=legend[0]),Line2D([0], [0], color='k',linestyle ="--", lw=1, label=legend[1])]
        #     if(otherdperp is not None):
        #         for i in range(len(dperp_other)):
        #             legend_elements.append(Line2D([0], [0], color='k',linestyle =line[i], lw=1, label=legend[i+2]))
        #     ax1.legend(handles=legend_elements,loc = "upper center")
        #     plt.savefig("{}.pdf".format(nameout),format ="pdf",dpi=200)








    @staticmethod
    def plot_histo_mean_distance_comparison(dpername1,dpername2,densityname1,densityname2,nameout,legend,coupled_plot=False):
        if(type(dpername1) is str):
            zpar, dperp = PixelAnalizer.read_dperp_file(dpername1)
            zpar_comparison, dperp_comparison = PixelAnalizer.read_dperp_file(dpername2)
            zpar, density = PixelAnalizer.read_density_file(densityname1)
            zpar_comparison, density_comparison = PixelAnalizer.read_density_file(densityname2)
        else :
            zpar, dperp  = np.mean([PixelAnalizer.read_dperp_file(dpername1[i]) for i in range(len(dpername1))],axis=0)
            zpar_comparison, dperp_comparison = np.mean([PixelAnalizer.read_dperp_file(dpername2[i]) for i in range(len(dpername2))],axis=0)
            zpar, density = np.mean([PixelAnalizer.read_density_file(densityname1[i]) for i in range(len(densityname1))],axis=0)
            zpar_comparison, density_comparison = np.mean([PixelAnalizer.read_density_file(densityname2[i]) for i in range(len(densityname2))],axis=0)
        PixelAnalizer.plot_mean_distance_density(zpar,dperp,density,nameout,coupled_plot=False,comparison=True,dperp_comparison=dperp_comparison,density_comparison=density_comparison,zpar_comparison=zpar_comparison,legend=legend)



    def compute_plot_histo_mean_distance(self,zpar,name_histo):
        dmin = self.pixel.compute_mean_distance_histogram(zpar)
        PixelAnalizer.plot_histogram_mean_distance(zpar,dmin,name_histo,nb_bins=50)


    def compute_plot_mean_distance_density(self,nameout,coupled=False,plot=True):
        (zpar,dperpz,densityz) = self.pixel.compute_mean_distance_density()
        if(plot):
            PixelAnalizer.write_dperp_file(zpar,dperpz,f'{nameout}_dperp_file')
            PixelAnalizer.write_density_file(zpar,densityz,f'{nameout}_density_file')
            PixelAnalizer.plot_mean_distance_density(zpar,dperpz,densityz,nameout,coupled_plot=coupled)
        return(zpar,dperpz,densityz)






    def analyze_pixels(self,compute_histo,compute_mean_distance_density,histo_zpar=None,name_histo="histogram_dperp",name_dperp="density_mean_distance",coupled_plot=False):
        if(compute_histo):
            self.compute_plot_histo_mean_distance(histo_zpar,name_histo)
        if(compute_mean_distance_density):
            self.compute_plot_mean_distance_density(name_dperp,coupled=coupled_plot)








class DeltaAnalyzer(object):

    def __init__(self,pwd,delta_path,center_ra=True,z_cut_min=None,z_cut_max=None,dec_cut_min=None,dec_cut_max=None,ra_cut_min=None,ra_cut_max=None,degree=True):
        self.pwd = pwd
        self.delta_path = delta_path
        self.center_ra = center_ra
        self.z_cut_min = z_cut_min
        self.z_cut_max = z_cut_max
        self.dec_cut_min = dec_cut_min
        self.dec_cut_max = dec_cut_max
        self.ra_cut_min = ra_cut_min
        self.ra_cut_max = ra_cut_max
        self.degree=degree


    def get_ra_dec(self):
        """ Obtain arrays of RA and DEC coordinates from a list or a name of a delta file in pickle, fits or ascii format"""
        namefile = glob.glob(os.path.join(self.delta_path,"delta-*.fits"))
        (ra,dec,z,zqso,ids,sigmas,deltas)  = get_deltas(namefile,center_ra=self.center_ra)
        pixel_coord = np.array([[ra[i],dec[i],z[i][j],sigmas[i][j],deltas[i][j],zqso[i],ids[i]] for i in range(len(ra)) for j in range(len(z[i]))])
        pixel_coord = pixel_coord[utils.cut_sky_catalog(pixel_coord[:,0],pixel_coord[:,1],pixel_coord[:,2],ramin=self.ra_cut_min,ramax=self.ra_cut_max,decmin=self.dec_cut_min,decmax=self.dec_cut_max,zmin=self.z_cut_min,zmax=self.z_cut_max)]
        z,zqso,ids,sigmas,deltas = pixel_coord[:,2],pixel_coord[:,5],pixel_coord[:,6],pixel_coord[:,3],pixel_coord[:,4]
        ra = np.array([pixel_coord[:,0][index] for index in sorted(np.unique(pixel_coord[:,0], return_index=True)[1])])
        dec = np.array([pixel_coord[:,1][index] for index in sorted(np.unique(pixel_coord[:,1], return_index=True)[1])])
        if(self.degree): ra,dec = np.degrees(ra),np.degrees(dec)
        return(ra,dec,z,zqso,ids,sigmas,deltas)




    @staticmethod
    def plot_snr_diagram(sigmas,deltas,plot_name,nb_bin=100):
        signal_noise = abs((deltas + 1)/sigmas)
        plt.hist(signal_noise,nb_bin,density=True)
        plt.xlabel("signal-to-noise ratio")
        plt.savefig(f"{plot_name}.pdf",format="pdf")


    @staticmethod
    def plot_ra_dec_diagram(ra,dec,plot_name,nb_cut=None,deg=True):
        plt.figure(figsize=(7,3.5))
        plt.plot(ra,dec,'b.', markersize=1.5)
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
        if(deg):
            plt.xlabel("RA [deg] (J2000)")
            plt.ylabel("DEC [deg] (J2000)")
        else:
            plt.xlabel("RA [rad] (J2000)")
            plt.ylabel("DEC [rad] (J2000)")
        plt.grid()
        plt.savefig("RA-DEC_diagram_{}.pdf".format(plot_name),format = "pdf")



    @staticmethod
    def plot_los_density(ra,dec,plot_name,nb_interval=20,different_sign_region=False):
        ra_interval = np.linspace(np.min(ra),np.max(ra),nb_interval)
        ra_size = abs((np.max(ra)-np.min(ra))/nb_interval)
        ra_array = []
        for i in range(len(ra_interval)-1):
            ra_array.append((ra_interval[i] + ra_interval[i+1])/2)
        ra_array=np.asarray(ra_array)
        density_array = np.zeros(ra_array.shape)
        density_array_plus = np.zeros(ra_array.shape)
        density_array_minus = np.zeros(ra_array.shape)
        for i in range(len(ra_interval)-1):
            mask =  (ra > ra_interval[i]) & (ra < ra_interval[i+1])
            density_array[i] = len(ra[mask])
            density_array_plus[i] = len(ra[mask & (dec >= 0)])
            density_array_minus[i] = len(ra[mask & (dec < 0)])
        maxdec = np.max(dec)
        mindec = np.min(dec)
        plt.figure()
        plt.plot(ra_array,density_array/abs(ra_size*(maxdec-mindec)))
        plt.title("LOS density in function of RA")
        plt.grid()
        plt.savefig(f"los_density_{plot_name}.pdf",format = "pdf")
        if(different_sign_region):
            plt.figure()
            plt.plot(ra_array,density_array_plus/abs(ra_size*maxdec))
            plt.title("LOS density in function of RA for DEC >= 0")
            plt.grid()
            plt.savefig(f"los_density_dec_positive_{plot_name}.pdf",format = "pdf")
            plt.figure()
            plt.plot(ra_array,density_array_minus/abs(ra_size*mindec))
            plt.title("LOS density in function of RA for DEC < 0")
            plt.grid()
            plt.savefig(f"los_density_dec_negative_{plot_name}.pdf",format = "pdf")



    @staticmethod
    def plot_los_histogram(ra,dec,plot_name,nb_bins=20,different_sign_region=False):
        plt.figure()
        plt.hist(ra,nb_bins)
        plt.title("LOS histogram in function of RA")
        plt.grid()
        plt.savefig(f"los_histogram_{plot_name}.pdf",format = "pdf")
        if(different_sign_region):
            plt.figure()
            plt.hist(ra[dec >= 0],nb_bins)
            plt.title("LOS histogram in function of RA for DEC >= 0")
            plt.grid()
            plt.savefig(f"los_histogram_dec_positive_{plot_name}.pdf",format = "pdf")
            plt.figure()
            plt.hist(ra[dec < 0],nb_bins)
            plt.title("LOS histogram in function of RA for DEC < 0")
            plt.grid()
            plt.savefig(f"los_histogram_dec_negative_{plot_name}.pdf",format = "pdf")



    @staticmethod
    def plot_histogram_deltas(deltas,plot_name,nb_bins=200,delta_min=-2,delta_max=2):
        plt.figure()
        data, bins , patches = plt.hist(deltas,nb_bins)
        plt.xlim([delta_min,delta_max])
        bin_centers= np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
        fit_function = lambda x, A, mu, sigma : A * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2))
        popt, pcov = curve_fit(fit_function, xdata=bin_centers, ydata=data, p0=[1, 0.0, 0.1])
        x = np.linspace(min(bins),max(bins),1000)
        y = fit_function(x, *popt)
        plt.plot(x,y,'r--',linewidth=2)
        mu,sigma = popt[1],popt[2]
        plt.text(0.8 ,2* np.max(data)/6,"mu = " +str(round(mu,8)))
        plt.text(0.8,1.5* np.max(data)/6,"sigma = " + str(round(sigma,8)))
        plt.savefig(f"delta_map_histo_{plot_name}.pdf",format="pdf")






    @staticmethod
    def plot_delta_binned_stat(redshift,deltas,redshift_qso,plot_name, nb_bins=50,error_bar=True) :

        list_RF = ((1 + redshift)/(1+ redshift_qso))*utils.lambdaLy
        list_lobs = (1 + redshift)*utils.lambdaLy

        plt.figure()
        bin_centers, means, errors = DeltaAnalyzer.hist_profile(list_lobs,deltas,nb_bins, (np.min(list_lobs),np.max(list_lobs)), (np.min(deltas),np.max(deltas)))
        if(error_bar) :plt.errorbar(x=bin_centers, y=means, yerr=errors, linestyle='none', marker='.',color='blue', label="dr8")
        else :plt.plot(bin_centers, means, linestyle='none', marker='.',color='blue', label="dr8")
        plt.ylabel("flux contrast pixels")
        plt.xlabel("observed wavelength")
        plt.savefig(f"delta_binned_statistics_obs_frame_{plot_name}.pdf",format="pdf")

        plt.figure()
        bin_centers, means, errors = DeltaAnalyzer.hist_profile(list_RF,deltas,nb_bins, (np.min(list_RF),np.max(list_RF)), (np.min(deltas),np.max(deltas)))
        if(error_bar) :plt.errorbar(x=bin_centers, y=means, yerr=errors, linestyle='none', marker='.',color='blue', label="dr8")
        else :plt.plot(bin_centers, means, linestyle='none', marker='.',color='blue', label="dr8")
        plt.ylabel("flux contrast pixels")
        plt.xlabel("rest frame wavelength")
        plt.savefig(f"delta_binned_statistics_rest_frame_{plot_name}.pdf",format="pdf")



    @staticmethod
    def hist_profile(x, y, bins, range_x, range_y) :
        w = (y>range_y[0]) & (y<range_y[1])
        means_result = binned_statistic(x[w], [y[w], y[w]**2], bins=bins, range=range_x, statistic='mean')
        nb_entries_result = binned_statistic(x[w], y[w], bins=bins, range=range_x, statistic='count')
        means, means2 = means_result.statistic
        nb_entries = nb_entries_result.statistic
        errors = np.sqrt(means2 - means**2)/np.sqrt(nb_entries)
        bin_edges = means_result.bin_edges
        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
        return bin_centers, means, errors



    ### COMPARE TWO DELTAs FOLDER ###


    def compare_deltas(self,delta_path2,log_scale=False):
        " Print different properties of two pickled deltas files to compare them"


        (ra,dec,redshift,zqso,ids,sigmas,delta)  = self.get_ra_dec()
        comp = DeltaAnalyzer(self.pwd,delta_path2,center_ra=self.center_ra,z_cut_min=self.z_cut_min,z_cut_max=self.z_cut_max,dec_cut_min=self.dec_cut_min,dec_cut_max=self.dec_cut_max,ra_cut_min=self.ra_cut_min,ra_cut_max=self.ra_cut_max,degree=self.degree)
        (ra_comp,dec_comp,redshift_comp,zqso_comp,ids_comp,sigmas_comp,delta_comp)  = comp.get_ra_dec()

        # Transform to list
        sigmas = np.concatenate(sigmas)
        sigmas_comp = np.concatenate(sigmas_comp)
        delta = np.concatenate(delta)
        delta_comp = np.concatenate(delta_comp)
        redshift = np.concatenate(redshift)
        redshift_comp = np.concatenate(redshift_comp)



        # Histogram of sigmas
        plt.figure()
        bins = np.linspace(np.min(sigmas),np.max(sigmas), 100)
        plt.hist(sigmas, bins, alpha=0.5, label="mocks",density=True,log=log_scale)
        plt.hist(sigmas_comp, bins, alpha=0.5, label="dataBOSS",density=True,log=log_scale)
        plt.xlabel("Noise on flux contrast " +r"$\sigma_{\delta_{F}}$")
        plt.legend(["Mocks","Stripe 82 data"])
        plt.savefig("histogram_sigmas_normalized_cutz.pdf",format="pdf")


        # Mean sigma in function of redshift
        plt.figure()
        z = np.linspace(np.min(redshift),np.max(redshift),50)
        sigma_mean = np.zeros(len(z)-1)
        sigma_mean_comp = np.zeros(len(z)-1)
        for i in range(len(z)-1):
            mask2 = (redshift < z[i+1])&(redshift >= z[i])
            mask3 = (redshift_comp < z[i+1])&(redshift_comp >= z[i])
            sigma_mean[i]=np.mean(sigmas[mask2])
            sigma_mean_comp[i]=np.mean(sigmas_comp[mask3])
        plt.plot(z[:-1],sigma_mean)
        plt.plot(z[:-1],sigma_mean_comp)
        plt.savefig("sigmas_redshift_bin_DR16_mocks.pdf",format="pdf")

        # Histogram of deltas
        plt.figure()
        bins = np.linspace(np.min(delta),np.max(delta), 100)
        plt.hist(delta, bins, alpha=0.5, label="mocks",density=True,log=log_scale)
        plt.hist(delta_comp, bins, alpha=0.5, label="dataBOSS",density=True,log=log_scale)
        plt.legend(["mocks","dataBOSS"])
        plt.savefig("histo_deltas_1D_mocks_DR16.pdf",format="pdf")


        # Print scalar statistical data
        print("Redshift interval for 1 =",np.max(redshift),np.min(redshift))
        print("Redshift interval for 2 =",np.max(redshift_comp),np.min(redshift_comp))
        print("Maximal sigma for 1 =",np.max(sigmas))
        print("Maximal sigma for 2 =",np.max(sigmas_comp))
        print("Mean sigma for 1 =",np.mean(sigmas))
        print("Mean sigma for 2 =",np.mean(sigmas_comp))
        print("Median sigma for 1 =",np.median(sigmas))
        print("Median sigma for 2 =",np.median(sigmas_comp))
        print("Mean delta for 1 =",np.mean(delta))
        print("Mean delta for 1 =",np.mean(delta_comp))





    #### MAIN ROUTINES ####



    def analyze_deltas(self,plot_name,plot_ra_dec=False,plot_histo_ra_dec=False,plot_density_ra_dec=False,plot_histo_delta=False,plot_snr=False,plot_binned_delta=False,nb_cut_ra_dec=None,nb_interval_histo_ra_dec = 20,nb_interval_density_ra_dec = 20,nb_bins_histo_ra_dec = 20,different_sign_region = False,nb_bins_histo_delta=200,delta_min_histo_delta=-2,delta_max_histo_delta=2,nb_bins_snr=100, nb_bins_binned_stat=50,error_bar_binned_stat=True):
        (ra,dec,z,zqso,ids,sigmas,deltas)=self.get_ra_dec()

        if(plot_ra_dec):DeltaAnalyzer.plot_ra_dec_diagram(ra,dec,plot_name,nb_cut=nb_cut_ra_dec,deg=self.degree)
        if(plot_density_ra_dec):DeltaAnalyzer.plot_los_density(ra,dec,plot_name,nb_interval=nb_interval_density_ra_dec,different_sign_region=different_sign_region)
        if(plot_histo_ra_dec):DeltaAnalyzer.plot_los_histogram(ra,dec,plot_name,nb_bins=nb_bins_histo_ra_dec,different_sign_region=different_sign_region)

        if(plot_histo_delta):DeltaAnalyzer.plot_histogram_deltas(deltas,plot_name,nb_bins=nb_bins_histo_delta,delta_min=delta_min_histo_delta,delta_max=delta_max_histo_delta)
        if(plot_snr):DeltaAnalyzer.plot_snr_diagram(sigmas,deltas,plot_name,nb_bin=nb_bins_snr)
        if(plot_binned_delta):DeltaAnalyzer.plot_delta_binned_stat(z,deltas,zqso,plot_name, nb_bins=nb_bins_binned_stat,error_bar=error_bar_binned_stat)









class DeltaModifier(object):

    def __init__(self,pwd,delta_path,center_ra=True,z_cut_min=None,z_cut_max=None,dec_cut_min=None,dec_cut_max=None,ra_cut_min=None,ra_cut_max=None,degree=True):
        self.pwd = pwd
        self.delta_path = delta_path
        self.center_ra = center_ra
        self.dec_cut_min = dec_cut_min
        self.dec_cut_max = dec_cut_max
        self.ra_cut_min = ra_cut_min
        self.ra_cut_max = ra_cut_max



    def modify_deltas(self,name_out,number_cut,random_density_parameter=False,number_repeat=1,iterative_selection_parameters=None):
        deltas = self.get_new_healpix(number_cut,random_density_parameter=random_density_parameter,number_repeat=number_repeat,iterative_selection_parameters=None)
        self.save_deltas(deltas,name_out,number_cut)


    def save_deltas(self,deltas,name_out,number_cut):
        for cut in range(number_cut):
            delta = tomographic_objects.Delta(name=f"{name_out}_{cut}.fits",delta_file=deltas[cut])
            delta.write_from_delta_list()


    def get_new_healpix(self,number_cut,random_density_parameter=False,number_repeat=1,iterative_selection_parameters=None):
        ramax,ramin,decmax,decmin = self.ra_cut_min,self.ra_cut_max,self.dec_cut_min,self.dec_cut_max
        if(random_density_parameter is None):
            deltas = self.create_healpix(ramax,ramin,decmax,decmin,number_cut,random=None,return_len_ra=False)
        else :
            if (number_repeat == 1):
                deltas = self.create_healpix(ramax,ramin,decmax,decmin,number_cut,random=random_density_parameter,return_len_ra=False)
            else :
                if(iterative_selection_parameters is None):return KeyError("Please dictionary parameter for the iterative selection")
                deltas = self.iterate_healpix_creation(ramax,ramin,decmax,decmin,number_cut,random_density_parameter,number_repeat,iterative_selection_parameters)
        return(deltas)


    def create_healpix(self,ramax,ramin,decmax,decmin,number_cut,random=None,return_len_ra=False):
        namefile = glob.glob(os.path.join(self.delta_path,"delta-*.fits"))
        deltas ={}
        ra_array = []
        dec_array = []
        for cut in range(number_cut):
            deltas[cut]=[]
        for i in range(len(namefile)) :
            delta_tomo = tomographic_objects.Delta(name=namefile[i],pk1d_type=True)
            delta_tomo.read()
            for j in range(1,len(delta_tomo.delta_file)):
                delta = delta_tomo.read_line(j)
                if(self.center_ra_to_zero):
                    if(delta.ra *180 /np.pi  > 180):
                        ra =((delta.ra * 180 / np.pi)-360)
                    else:
                        ra = ((delta.ra * 180 / np.pi))
                else :
                    ra = ((delta.ra * 180 / np.pi))
                dec = (delta.dec* 180 / np.pi)
                if((ra>ramin)&(ra<ramax)&(dec>decmin)&(dec<decmax)):
                    ra_array.append(ra)
                    dec_array.append(dec)
                    for cut in range(number_cut):
                        interval_ra =  ((cut)/(number_cut))*(ramax-ramin) + ramin  , ((cut + 1)/(number_cut))*(ramax-ramin) + ramin
                        if((ra>interval_ra[0])&(ra<=interval_ra[1])):
                            if(self.center_ra_to_zero):
                                delta.ra = delta.ra - 2*np.pi
                            deltas[cut].append(delta)
            delta_tomo.close()
        if (random is not None):
            deltas = self.randomize_choice_of_los(deltas,random,ramax,ramin,decmax,decmin,number_cut,len(ra_array))
        if(return_len_ra):
            return(deltas,len(ra_array))
        else :
            return(deltas)


    def randomize_choice_of_los(self,deltas_dict,random,ramax,ramin,decmax,decmin,number_cut,number_ra):
        deltas = deltas_dict.copy()
        density = number_ra/((ramax-ramin)*(decmax-decmin))
        utils.Logger.add("density before random choice ="  + str(density))
        random_cut = random/density
        ra_random = []
        dec_random = []
        for cut in range(number_cut) :
            number_of_delta_to_select = int(round(len(deltas[cut])*random_cut,0))
            deltas[cut] = sample(deltas[cut],number_of_delta_to_select)
            for i in range(len(deltas[cut])):
                ra_random.append(deltas[cut][i].ra*180 /np.pi)
                dec_random.append(deltas[cut][i].dec*180 /np.pi)
        density = len(ra_random)/((ramax-ramin)*(decmax-decmin))
        utils.Logger.add("density after random choice ="  + str(density))
        return(deltas)



    def iterate_healpix_creation(self,ramax,ramin,decmax,decmin,number_cut,random_density_parameter,number_repeat,iterative_selection_parameters,property_file_name):
        """ To optimize & test"""
        density_names = iterative_selection_parameters["density_names"]
        dperp_names = iterative_selection_parameters["separation_names"]
        zmin = iterative_selection_parameters["zmin"]
        zmax = iterative_selection_parameters["zmax"]
        Om = iterative_selection_parameters["Om"]
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(Om)
        suplementary_parameters = iterative_selection_parameters["suplementary_parameters"]
        coordinate_transform = iterative_selection_parameters["coordinate_transform"]
        density_ref ={}
        dperp_ref = {}
        for cut in range(number_cut):
            density_ref[cut] = PixelAnalizer.read_density_file(density_names[cut])[1]
            dperp_ref[cut] = PixelAnalizer.read_dperp_file(dperp_names[cut])[1]
        min_density, min_dperp = np.inf,np.inf
        delta_to_keep = {}
        utils.Logger.add("Beginning of the random iteration selection")
        deltas_dict,number_ra = self.create_healpix(ramax,ramin,decmax,decmin,number_cut,random=False,return_len_ra=True)
        deltas_random={}
        for i in range(number_repeat):
            utils.Logger.add("repeat "+ str(i))
            utils.Logger.add("deltas "+ str(i) + " computed")
            diff_dperp , diff_density = [], []
            deltas_random = self.randomize_choice_of_los(deltas_dict,random_density_parameter,ramax,ramin,decmax,decmin,number_cut,number_ra)
            for cut in range(number_cut):
                delta_file = tomographic_objects.Delta(delta_file=None,delta_array=deltas_random[cut],name="delta_" + str(cut) + "_to_test.pickle")
                delta_file.write_from_delta_list()
                namefile = "delta_" + str(cut) + "_to_test.pickle"
                (ra,dec,z,zqso,ids,sigmas,deltas) = get_deltas(namefile)
                sky_deltas = np.array([[ra[i],dec[i],z[i][j],sigmas[i][j],deltas[i][j]] for i in range(len(ra)) for j in range(len(z[i]))])
                sky_deltas = sky_deltas[utils.cut_sky_catalog(sky_deltas[:,0],sky_deltas[:,1],sky_deltas[:,2],ramin=ramin,ramax=ramax,decmin=decmin,decmax=decmax,zmin=zmin,zmax=zmax)]
                cartesian_deltas = np.zeros(sky_deltas.shape)
                cartesian_deltas[:,0],cartesian_deltas[:,1],cartesian_deltas[:,2] = utils.convert_sky_to_cartesian(sky_deltas[:,0],sky_deltas[:,1],sky_deltas[:,2],coordinate_transform,rcomov=rcomov,distang=distang,suplementary_parameters=suplementary_parameters)
                pixel = tomographic_objects.Pixel.init_from_property_files(property_file_name,pixel_array=cartesian_deltas,name=None)
                pixel_analyzer = PixelAnalizer(self.pwd,pixel=pixel)
                (zpar,dperpz,densityz) = pixel_analyzer.compute_plot_mean_distance_density("",plot=False)
                diff_dperp.append(np.mean(abs(np.array(dperpz) -np.array(dperp_ref[cut]))))
                diff_density.append(np.mean(abs(np.array(densityz) -np.array(density_ref[cut]))))
            utils.Logger.add("Mean difference in term LOS density : {}".format(np.mean(diff_density)))
            utils.Logger.add("Mean difference in term of Mean LOS separation : {}".format(np.mean(diff_dperp)))
            if((np.mean(diff_density)<min_density)&(np.mean(diff_dperp)<min_dperp)):
                utils.Logger.add("Better at the repeat " + str(i))
                min_density = np.mean(diff_density)
                min_dperp = np.mean(diff_dperp)
                delta_to_keep = deltas_random
        for cut in range(number_cut):
            os.remove("delta_" + str(cut) + "_to_test.pickle")
        utils.Logger.add("End of the random iteration selection")
        return(delta_to_keep)
