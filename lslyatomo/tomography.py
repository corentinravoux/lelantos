#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

Author: Corentin Ravoux

Description : Classes used to treat input and output of the Dachshund software.
Can be used to plot the slices of a Tomographic map or to return a Paraview
file for 3D visualization. Can be also used to compare a simulation to a data
analysis.
"""



#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################


import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lslyatomo import tomographic_objects
from lslyatomo import utils



def create_merged_map(launching_file_name,map_name,property_file):
    map_merged = tomographic_objects.TomographicMap.init_by_merging(launching_file_name,map_name,property_file)
    map_merged.write()

def rebin_map(map_name,property_file,new_shape,new_name):
    map_class = tomographic_objects.TomographicMap.init_from_property_files(property_file,name=map_name)
    map_class.read()
    map_class.rebin_map(new_shape)
    map_class.name = new_name
    map_class.write()

def create_distance_map(map_name,pixel_file,property_file,nb_process=1,radius_local=50):
    pixel = tomographic_objects.Pixel(name=pixel_file)
    pixel.read()
    tomographic_map = tomographic_objects.TomographicMap.init_from_property_files(property_file)
    distance_map = tomographic_objects.DistanceMap.init_by_computing(pixel,tomographic_map,map_name,nb_process=nb_process,radius_local=radius_local)
    distance_map.write()





#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################







class TomographyPlot(object):

    def __init__(self,pwd,map_name=None,map_shape=None,pixel_name=None,property_file=None,**kwargs):
        self.pwd = pwd
        self.map_name = map_name
        self.map_shape = map_shape
        self.pixel_name = pixel_name
        self.property_file = property_file
        self.kwargs = kwargs



    def load_tomographic_objects(self,qso=None,void=None,galaxy=None,distance_mask=None,cut_plot=None):
        tomographic_map = tomographic_objects.TomographicMap.init_classic(name=self.map_name,shape=self.map_shape,property_file=self.property_file)
        tomographic_map.read()

        pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map = None,None,None,None,None
        if(self.pixel_name is not None):
            if(self.property_file is not None):
                pixel = tomographic_objects.Pixel.init_from_property_files(self.property_file,name=self.pixel_name)
            else:
                pixel = tomographic_objects.Pixel(name=self.pixel_name)
            pixel.read()
        if(qso is not None):
            quasar_catalog = tomographic_objects.QSOCatalog.init_from_fits(qso)
        if(void is not None):
            void_catalog = tomographic_objects.VoidCatalog.init_from_fits(void)
        if(galaxy is not None):
            void_catalog = tomographic_objects.GalaxyCatalog.init_from_fits(galaxy)
        if(distance_mask is not None):
            dist_map = tomographic_objects.DistanceMap.init_from_tomographic_map(tomographic_map,name=distance_mask)
            dist_map.read()

        if(cut_plot is not None):
            self.cut_objects(cut_plot,tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map)

        return(tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map)




    def cut_objects(self,cut_plot,tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map):

        size_map = np.array(cut_plot)*tomographic_map.size
        shape_map = np.round(np.array(cut_plot)*tomographic_map.shape,0).astype(int)
        tomographic_map.shape = shape_map
        tomographic_map.size = size_map
        tomographic_map.map_array = tomographic_map.map_array[0:shape_map[0],0:shape_map[1],0:shape_map[2]]

        if(dist_map is not None):
            dist_map.shape = shape_map
            dist_map.size = size_map
            dist_map.map_array = dist_map.map_array[0:shape_map[0],0:shape_map[1],0:shape_map[2]]


        if(pixel is not None):
            mask = (pixel.pixel_array[:,0] < size_map[0])&(pixel.pixel_array[:,1] < size_map[1])&(pixel.pixel_array[:,2] < size_map[2])
            pixel.pixel_array = pixel.pixel_array[mask]

        if(quasar_catalog is not None):
            mask = (quasar_catalog.coord[:,0] < size_map[0])&(quasar_catalog.coord[:,1] < size_map[1])&(quasar_catalog.coord[:,2] < size_map[2])
            quasar_catalog.apply_mask(mask)

        if(void_catalog is not None):
            mask = (void_catalog.coord[:,0] < size_map[0])&(void_catalog.coord[:,1] < size_map[1])&(void_catalog.coord[:,2] < size_map[2])
            void_catalog.apply_mask(mask)

        if(galaxy_catalog is not None):
            mask = (galaxy_catalog.coord[:,0] < size_map[0])&(galaxy_catalog.coord[:,1] < size_map[1])&(galaxy_catalog.coord[:,2] < size_map[2])
            galaxy_catalog.apply_mask(mask)

    def mask_tomographic_objects(self,void_catalog,dist_map,tomographic_map,criteria_distance_mask = None,minimal_void_crossing = None):
        mask_void = None
        if(void_catalog is not None):
            if(minimal_void_crossing is not None):
                crossing_param = void_catalog.crossing_param
                mask_void = crossing_param > minimal_void_crossing
        if(dist_map is not None):
            mask_distance = dist_map.get_mask_distance(criteria_distance_mask)
            tomographic_map.mask_map_from_mask(mask_distance)
        return(mask_void)

    def delta_histogram(self,listdeltas,nb_bins,norm=True,gauss_fit=True,alpha=1):
        data, bins , patches = plt.hist(listdeltas,nb_bins,density=norm,alpha=alpha)
        if(gauss_fit):
            bin_centers= np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
            fit_function = lambda x,A,mu,sigma : A * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2))
            popt, pcov = curve_fit(fit_function, xdata=bin_centers, ydata=data, p0=[1, 0.0, 0.1])
            x = np.linspace(min(bins),max(bins),1000)
            y = fit_function(x, *popt)
            mu,sigma = popt[1],popt[2]
            print(mu,sigma)
            plt.plot(x,y,'r--',linewidth=2)



    def select_pixel_in(self,center_mpc,space_mpc,pixel,index_direction):
        mask = pixel.pixel_array[:,index_direction] < center_mpc + space_mpc/2
        mask &= pixel.pixel_array[:,index_direction] >= center_mpc - space_mpc/2
        return(mask)

    def select_qso_in(self,center_mpc,space_mpc,quasar_catalog,index_direction):
        mask = quasar_catalog.coord[:,index_direction] < center_mpc + space_mpc/2
        mask &= quasar_catalog.coord[:,index_direction] >= center_mpc - space_mpc/2
        return(mask)

    def select_void_in(self,center_mpc,space_mpc,void,index_direction):
        mask = void.coord[:,index_direction] - void.radius < center_mpc
        mask &= void.coord[:,index_direction] + void.radius >= center_mpc
        return(mask)

    def select_galaxy_in(self,center_mpc,space_mpc,galaxy,index_direction):
        mask = galaxy.coord[:,index_direction] < center_mpc + space_mpc/2
        mask &= galaxy.coord[:,index_direction] >= center_mpc - space_mpc/2
        return(mask)

    @staticmethod
    def get_direction_informations(direction,rotate,size_map):
        if (direction.lower() == "x")|(direction.lower() == "ra"):
            x_index, y_index, index_direction = 2, 1, 0
        elif (direction.lower() == "y")|(direction.lower() == "dec"):
            x_index, y_index, index_direction = 2, 0, 1
        elif (direction.lower() == "z")|(direction.lower() == "redshift"):
            x_index, y_index, index_direction = 1, 0, 2
        if(rotate): x_index,y_index = y_index,x_index
        index_dict = {0:"x",1:"y",2:"z"}


        if(rotate):extentmap = [0,size_map[x_index],0,size_map[y_index]]
        else:extentmap = [0,size_map[x_index],size_map[y_index],0]

        xlab = f"Relative comoving distance in the {index_dict[x_index]} direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
        ylab = f"Relative comoving distance in the {index_dict[y_index]} direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"

        return(x_index, y_index, index_direction, extentmap, xlab, ylab)



    def print_one_slice(self,name,tomographic_map,direction,space_mpc,center_mpc,deltamin,deltamax,pixel=None,quasar_catalog=None,void_catalog=None,mask_void=None,galaxy_catalog=None,rotate = False,hd=False,redshift_axis=False,xlim=None,ylim=None):

        size_map = tomographic_map.size
        pixel_per_mpc = tomographic_map.pixel_per_mpc

        (x_index, y_index, index_direction, extentmap, xlab, ylab) = TomographyPlot.get_direction_informations(direction,rotate,size_map)

        center_pix = int(round(center_mpc*pixel_per_mpc[index_direction],0))
        if direction == "x":
            map_slice = tomographic_map.map_array[center_pix,:,:]
            if(rotate): map_slice = np.transpose(np.flip(tomographic_map.map_array[center_pix,:,:],axis=1))
        elif direction == "y" :
            map_slice = tomographic_map.map_array[:,center_pix,:]
            if(rotate): map_slice = np.transpose(np.flip(tomographic_map.map_array[:,center_pix,:],axis=1))
        elif direction == "z" :
            map_slice = tomographic_map.map_array[:,:,center_pix]
            if(rotate): map_slice = np.transpose(np.flip(tomographic_map.map_array[:,:,center_pix],axis=1))


        pixel_in, pixel_bis_in, qso_in, qso_bis_in, void_in, void_bis_in, galaxy_in = None, None, None, None, None, None, None
        if(pixel is not None):
            mask_pixel_in = self.select_pixel_in(center_mpc,space_mpc,pixel,index_direction)
            pixel_in = pixel.pixel_array[mask_pixel_in]
            mask_pixel_bis_in = (self.select_pixel_in(center_mpc,2*space_mpc,pixel,index_direction))&(~mask_pixel_in)
            pixel_bis_in = pixel.pixel_array[mask_pixel_bis_in]
        if(quasar_catalog is not None):
            mask_qso_in = self.select_qso_in(center_mpc,space_mpc,quasar_catalog,index_direction)
            qso_in = quasar_catalog.coord[mask_qso_in]
            mask_qso_bis_in = (self.select_qso_in(center_mpc,2*space_mpc,quasar_catalog,index_direction))&(~mask_qso_in)
            qso_bis_in = quasar_catalog.coord[mask_qso_bis_in]
        if(void_catalog is not None):
            mask_void_in = self.select_void_in(center_mpc,space_mpc,void_catalog,index_direction)
            void_coord_in = void_catalog.coord[mask_void_in]
            void_radius_in = void_catalog.radius[mask_void_in]
            effective_radius = np.sqrt(void_radius_in**2 - (void_coord_in[:,index_direction] - center_mpc)**2)
            void_in = np.transpose(np.vstack([void_coord_in[:,0],void_coord_in[:,1],void_coord_in[:,2],effective_radius]))
            if(mask_void is not None):
                void_bis_in = void_in[~mask_void]
                void_in = void_in[mask_void]
        if(galaxy_catalog is not None):
            mask_galaxy_in = self.select_void_in(center_mpc,space_mpc,void_catalog,index_direction)
            galaxy_in = np.transpose(np.vstack([galaxy_catalog.coord[:,0],galaxy_catalog.coord[:,1],galaxy_catalog.coord[:,2],galaxy_catalog.error_z]))[mask_galaxy_in]

        name_plot = f"{name}_direction_{direction}_mpc_{center_mpc}"
        TomographyPlot.plot_slice(map_slice,extentmap,xlab,ylab,name_plot,deltamin,deltamax,x_index,y_index,pixel_in=pixel_in,pixel_bis_in=pixel_bis_in,qso_in=qso_in,qso_bis_in=qso_bis_in,void_in=void_in,void_bis_in=void_bis_in,galaxy_in=galaxy_in,redshift_axis=redshift_axis,tomographic_map=tomographic_map,rotate=rotate,hd=hd,xlim=xlim,ylim=ylim,**self.kwargs)



    @staticmethod
    def plot_slice(map_slice,extentmap,xlab,ylab,name,deltamin,deltamax,x_index,y_index,pixel_in=None,pixel_bis_in=None,qso_in=None,qso_bis_in=None,void_in=None,void_bis_in=None,galaxy_in=None,redshift_axis=False,tomographic_map=None,rotate=False,hd=False,xlim=None,ylim=None,**kwargs):
        plt.figure()
        fig = plt.gcf()
        ax = plt.gca()
        size = fig.get_size_inches()
        fig.set_size_inches(1.75*size)
        plt.xlabel(xlab)
        plt.ylabel(ylab)

        im = TomographyPlot.add_elements(map_slice,extentmap,deltamin,deltamax,x_index,y_index,pixel_in=pixel_in,pixel_bis_in=pixel_bis_in,qso_in=qso_in,qso_bis_in=qso_bis_in,void_in=void_in,void_bis_in=void_bis_in,galaxy_in=galaxy_in,**kwargs)

        if(rotate):
            cbar = plt.colorbar(im,ax=ax, orientation='horizontal', fraction=.1)
        else:
            cbar = plt.colorbar(im,ax=ax, fraction=.1)
        cbar.set_label("Flux contrast " + r"$\delta_{Fmap}$")

        if(redshift_axis):
            TomographyPlot.add_reshift_axe(tomographic_map,rotate=rotate,**kwargs)

        if(xlim is not None):plt.xlim(xlim)
        if(ylim is not None):plt.ylim(ylim)

        if(hd):plt.savefig(f"{name}.pdf",format="pdf", dpi=200)
        else:plt.savefig(f"{name}.pdf",format="pdf")

        plt.show()
        plt.close()


    @staticmethod
    def add_elements(map_slice,extentmap,deltamin,deltamax,x_index,y_index,pixel_in=None,pixel_bis_in=None,qso_in=None,qso_bis_in=None,void_in=None,void_bis_in=None,galaxy_in=None,**kwargs):
        im =plt.imshow(map_slice, interpolation='bilinear',cmap='jet_r',vmin = deltamin, vmax = deltamax,extent = extentmap)
        if(pixel_in is not None):
            plt.scatter(pixel_in[:,x_index],pixel_in[:,y_index],linewidth=utils.return_key(kwargs,"linewidth_pixel",1),marker=utils.return_key(kwargs,"marker_pixel","s"),color="k")
        if(pixel_bis_in is not None):
            plt.scatter(pixel_bis_in[:,x_index],pixel_bis_in[:,y_index]
                        ,linewidth=0.5*utils.return_key(kwargs,"linewidth_pixel",1)
                        ,marker=utils.return_key(kwargs,"marker_pixel",".")
                        ,color=utils.return_key(kwargs,"grey_pixel","0.5")
                        ,alpha=utils.return_key(kwargs,"transparency_pixel",0.5))
        if(qso_in is not None):
            plt.plot(qso_in[:,x_index],qso_in[:,y_index], "k*",linewidth=1)
        if(qso_bis_in is not None):
            plt.plot(qso_bis_in[:,x_index],qso_bis_in[:,y_index], "k*",linewidth=1,fillstyle='none')
        if(void_in is not None):
            for i in range(len(void_in)):
                circle = plt.Circle((void_in[i,x_index],void_in[i,y_index]),void_in[i,3],fill=False,color=utils.return_key(kwargs,"color_voids","r"))
                plt.gcf().gca().add_artist(circle)
        if(void_bis_in is not None):
            for i in range(len(void_bis_in)):
                circle = plt.Circle((void_bis_in[i,x_index],void_bis_in[i,y_index]),void_bis_in[i,3],fill=False,color=utils.return_key(kwargs,"color_voids_outside","k"))
                plt.gcf().gca().add_artist(circle)
        if(galaxy_in):
            plt.plot(galaxy_in[:,x_index],galaxy_in[:,y_index],"rx")
            plt.errorbar(galaxy_in[:,x_index],galaxy_in[:,y_index],xerr=2*galaxy_in[:,3], capsize = 0.01, ecolor = 'red',fmt = 'none')
        return(im)


    @staticmethod
    def add_reshift_axe(tomographic_map,rotate=False,**kwargs):
        if(tomographic_map.coordinate_transform != "middle"):
            return()

        fig = plt.gcf()
        ax1 = fig.axes[0]
        ax2 = fig.add_axes(ax1.get_position(), frameon=False)


        if(rotate):
            bounds = ax1.get_ybound()+ tomographic_map.boundary_cartesian_coord[0][2]
        else:
            bounds = ax1.get_xbound()+ tomographic_map.boundary_cartesian_coord[0][2]

        z_array = np.linspace(bounds[0],bounds[1],1000)
        (rcomov,distang,inv_rcomov,inv_distang) = utils.get_cosmo_function(tomographic_map.Omega_m)
        redshifts = utils.convert_z_cartesian_to_sky_middle(z_array,inv_rcomov)
        redshift_to_plot = np.unique(np.around(redshifts,decimals=1))
        tick_position = np.array([np.argmin(np.abs(redshifts - redshift_to_plot[i])) for i in range(len(redshift_to_plot))])/(len(redshifts) - 1)


        ax1.tick_params(labelbottom='on',labelleft='on', bottom='on', left='on')
        if(rotate):
            ax2.set_yticks(tick_position)
            ax2.set_yticklabels(redshift_to_plot)
            ax2.tick_params(labelright='on', right='on',labelbottom=None,labelleft=None, bottom=None, left=None)
            ax2.set_ylabel('Redshift')
        else:
            ax2.set_xticks(tick_position)
            ax2.set_xticklabels(redshift_to_plot)
            ax2.tick_params(labeltop='on', top='on',labelbottom=None,labelleft=None, bottom=None, left=None)
            ax2.set_xlabel('Redshift')


        position_redshift_axe = utils.return_key(kwargs, "position_redshift_axe", "other")
        outward_redshift_axe = utils.return_key(kwargs, "outward_redshift_axe", 50)
        if(position_redshift_axe == "other"):
            if(rotate):
                ax2.yaxis.set_label_position('right')
                ax2.yaxis.set_ticks_position('right')
            else:
                ax2.xaxis.set_label_position('top')
                ax2.xaxis.set_ticks_position('top')
        elif(position_redshift_axe == "same"):
            if(rotate):
                ax2.yaxis.set_label_position('left')
                ax2.yaxis.set_ticks_position('left')
                ax2.spines['left'].set_position(('outward', outward_redshift_axe)) # put redshift axis at the bottom
            else:
                ax2.xaxis.set_label_position('bottom')
                ax2.xaxis.set_ticks_position('bottom')
                ax2.spines['bottom'].set_position(('outward', outward_redshift_axe))





    def plot(self,name,direction,space,center_mpc,deltamin,deltamax,qso=None,void=None,galaxy=None,distance_mask = None,criteria_distance_mask = None,rotate = False,hd=False,minimal_void_crossing = None,redshift_axis=False,cut_plot=None):
        tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map = self.load_tomographic_objects(qso=qso,void=void,galaxy=galaxy,distance_mask = distance_mask,cut_plot=cut_plot)
        mask_void = self.mask_tomographic_objects(void_catalog,dist_map,tomographic_map,criteria_distance_mask = criteria_distance_mask, minimal_void_crossing = minimal_void_crossing)
        self.print_one_slice(name,tomographic_map,direction,space,center_mpc,deltamin,deltamax,pixel=pixel,quasar_catalog=quasar_catalog,void_catalog=void_catalog,mask_void=mask_void,galaxy_catalog=galaxy_catalog,rotate = rotate,hd=hd,redshift_axis=redshift_axis)



    def plot_all_slice(self,name,direction,space,deltamin,deltamax,qso=None,void=None,galaxy=None,distance_mask = None,criteria_distance_mask = None,rotate = False,hd=False,minimal_void_crossing = None,redshift_axis=False,cut_plot=None):
        tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map = self.load_tomographic_objects(qso=qso,void=void,galaxy=galaxy,distance_mask = distance_mask,cut_plot=cut_plot)
        mask_void = self.mask_tomographic_objects(void_catalog,dist_map,tomographic_map,criteria_distance_mask = criteria_distance_mask, minimal_void_crossing = minimal_void_crossing)
        center_mpc = space/2
        index_dict = {"x":0,"y":1,"z":2,"ra":0,"dec":1,"redshift":2}
        while(center_mpc + space/2 <tomographic_map.size[index_dict[direction]]):
            self.print_one_slice(name,tomographic_map,direction,space,center_mpc,deltamin,deltamax,pixel=pixel,quasar_catalog=quasar_catalog,void_catalog=void_catalog,mask_void=mask_void,galaxy_catalog=galaxy_catalog,rotate = rotate,hd=hd,redshift_axis=redshift_axis)
            center_mpc += space
            print(center_mpc,tomographic_map.size[index_dict[direction]])


    def plot_integrate_image(self,zmin,zmax,name,deltamin,deltamax,rotate=False,hd=False,cut_plot=None):
        tomographic_map = self.load_tomographic_objects(cut_plot=cut_plot)[0]
        (x_index, y_index, index_direction, extentmap, xlab, ylab) = TomographyPlot.get_direction_informations("z",rotate,tomographic_map.size)
        map_data = tomographic_map.map_array
        i_pix_end = int(round((zmax*tomographic_map.pixel_per_mpc[index_direction]),0))
        integrated_map = np.mean(map_data[:,:,0:i_pix_end],axis=2)
        name = f"{name}_integrated_map"
        TomographyPlot.plot_slice(integrated_map,extentmap,xlab,ylab,name,deltamin,deltamax,x_index,y_index,rotate=rotate,hd=hd,**self.kwargs)
        return(integrated_map)


    def plot_catalog_centered_maps(self,direction,name,space,deltamin,deltamax,void,nb_plot,radius_centered,qso=None,rotate = False, hd=False):
        load = self.load_tomographic_objects(void=void,qso=qso)
        tomographic_map,pixel,void_catalog = load[0],load[1],load[3]
        if qso is not None: qso = load[2]
        arg = np.array(void_catalog.radius).argsort()[-nb_plot:][::-1]
        coords = np.array(void_catalog.coord)[arg]
        (x_index, y_index, index_direction, extentmap, xlab, ylab) = TomographyPlot.get_direction_informations(direction,rotate,tomographic_map.size)
        for i in range(len(coords)):
            xlim = [coords[i][x_index]-radius_centered,coords[i][x_index]+radius_centered]
            ylim = [coords[i][y_index]-radius_centered,coords[i][y_index]+radius_centered]
            center_mpc = coords[i][index_direction]
            name = "{}_number{}_radius{}".format(name,i,np.array(void_catalog.radius)[arg][i])
            self.print_one_slice(name,tomographic_map,direction,space,center_mpc,deltamin,deltamax,pixel=pixel,quasar_catalog=qso,void_catalog=void_catalog,rotate = rotate,hd=hd,xlim=xlim,ylim=ylim)



    def plot_delta_histogram(self,name,nb_bins,gauss_fit=True,norm=True,distance_mask=None,criteria_distance_mask=None,log_scale=True,cut_plot=None):
        tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map = self.load_tomographic_objects(distance_mask = distance_mask,cut_plot=cut_plot)
        self.mask_tomographic_objects(void_catalog,dist_map,tomographic_map,criteria_distance_mask = criteria_distance_mask)
        listdeltas = tomographic_map.map_array.ravel()
        plt.figure()
        self.delta_histogram(listdeltas,nb_bins,norm=norm,gauss_fit=gauss_fit)
        if(norm):
            plt.title("Histogram of deltas normalized")
        else:
            plt.title("Histogram of deltas")
        if(log_scale):
            plt.yscale("log")
        plt.grid()
        plt.savefig("{}.pdf".format(name),format = "pdf")




    def plot_delta_histogram_comparison(self,name,name_second_map,nb_bins,legend,gauss_fit=True,norm=True,distance_mask=None,distance_mask2=None,criteria_distance_mask=None,log_scale=True,cut_plot=None):

        tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map = self.load_tomographic_objects(distance_mask = distance_mask,cut_plot=cut_plot)
        self.mask_tomographic_objects(void_catalog,dist_map,tomographic_map,criteria_distance_mask = criteria_distance_mask)
        listdeltas = tomographic_map.map_array.ravel()

        tomo_plot = TomographyPlot(self.pwd,map_name=name_second_map,map_shape=self.map_shape,pixel_name=self.pixel_name,property_file=self.property_file)
        tomographic_map2,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map2 = tomo_plot.load_tomographic_objects(distance_mask = distance_mask2,cut_plot=cut_plot)
        tomo_plot.mask_tomographic_objects(void_catalog,dist_map2,tomographic_map2,criteria_distance_mask = criteria_distance_mask)
        listdeltas2 = tomographic_map2.map_array.ravel()

        plt.figure()
        self.delta_histogram(listdeltas,nb_bins,norm=norm,gauss_fit=gauss_fit,alpha=0.5)
        self.delta_histogram(listdeltas2,nb_bins,norm=norm,gauss_fit=gauss_fit,alpha=0.5)
        if(norm):
            plt.title("Histogram of deltas normalized")
        else:
            plt.title("Histogram of deltas")
        if(log_scale):
            plt.yscale("log")
        plt.legend(legend)
        plt.xlabel("$\delta_{Fmap}$")
        plt.grid()
        plt.savefig("{}.pdf".format(name),format = "pdf")



    def compare_two_map(self,name,name_second_map,distance_mask,distance_second_mask,dist_extremum,bin_dist,legend,shuffle_map=None,cut_plot=None):

        tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map = self.load_tomographic_objects(distance_mask = distance_mask,cut_plot=cut_plot)

        tomo_plot = TomographyPlot(self.pwd,map_name=name_second_map,map_shape=self.map_shape,pixel_name=self.pixel_name,property_file=self.property_file)
        tomographic_map2,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map2 = tomo_plot.load_tomographic_objects(distance_mask = distance_second_mask,cut_plot=cut_plot)

        corr = []
        if(shuffle_map is not None):
            tomo_plot_shuffle = TomographyPlot(self.pwd,map_name=shuffle_map,map_shape=self.map_shape,pixel_name=self.pixel_name,property_file=self.property_file)
            tomographic_shuffle= tomo_plot_shuffle.load_tomographic_objects(cut_plot=cut_plot)[0]
            corr_shuffle = []

        dist_range = np.linspace(dist_extremum[0],dist_extremum[1],bin_dist)
        for i in range(len(dist_range)):
            mask = (dist_map < dist_range[i])&(dist_map2<dist_range[i])
            corr.append(np.corrcoef(tomographic_map.map_array[mask],tomographic_map2.map_array[mask])[0][1])
            if(shuffle_map is not None):
                corr_shuffle.append(np.corrcoef(tomographic_map.map_array[mask],tomographic_shuffle.map_array[mask])[0][1])
        plt.plot(dist_range,corr)
        plt.plot(dist_range,corr_shuffle)
        plt.legend(legend)
        plt.xlabel("Distance to the nearest los [" + r"$\mathrm{h^{-1}Mpc}$" + "]")
        plt.ylabel("Correlation coefficient")
        plt.grid()
        plt.savefig("{}.pdf".format(name),format = "pdf")
        return(corr,dist_range)



    # CR - move to tomographic_objects ?
    def compute_Pk3D(self,size_map,n_k,kmin,kmax,log=False,distance_mask=None,criteria_distance_mask=None,cut_plot=None):
        tomographic_map,pixel,quasar_catalog,void_catalog,galaxy_catalog,dist_map = self.load_tomographic_objects(distance_mask = distance_mask,cut_plot=cut_plot)
        self.mask_tomographic_objects(void_catalog,dist_map,tomographic_map,criteria_distance_mask = criteria_distance_mask)

        map_3D = tomographic_map.map_array
        map_fft_3D = np.fft.fftn(map_3D)
        number_Mpc_per_pixels = tomographic_map.mpc_per_pixel
        kx = np.fft.fftfreq(map_fft_3D.shape[0],number_Mpc_per_pixels[0])
        ky = np.fft.fftfreq(map_fft_3D.shape[1],number_Mpc_per_pixels[1])
        kz = np.fft.fftfreq(map_fft_3D.shape[2],number_Mpc_per_pixels[2])
        kx_space , ky_space, kz_space = np.meshgrid(kx,ky,kz)
        normalization_factor = (number_Mpc_per_pixels[0]*number_Mpc_per_pixels[1]*number_Mpc_per_pixels[2])/(self.shapeMap[0]*self.shapeMap[1]*self.shapeMap[2])
        power = normalization_factor * np.absolute(map_fft_3D)**2
        norm_k = np.array(map_fft_3D.shape)
        norm_k = np.sqrt(kx_space[:,:,:]**2 + ky_space[:,:,:]**2 + kz_space[:,:,:]**2)
        if(log):
            k_space = np.logspace(kmin,kmax,n_k)
        else :
            k_space = np.linspace(np.min(norm_k),np.max(norm_k),n_k)
        delta_k = (np.max(norm_k)-np.min(norm_k)) / n_k
        Pk_3D = np.zeros(k_space.shape)
        for i in range(len(Pk_3D)) :
            mask = (norm_k < k_space[i] + delta_k)&(norm_k >= k_space[i])
#            if(len(power[mask])==0):
#                Pk_3D[i] = -1
#            else :
            Pk_3D[i] = np.mean(power[mask])
        mask2 = Pk_3D != -1
        pk_3D_final = Pk_3D[mask2]
        k_space_final = k_space[mask2]
        del kx,ky,kz,kx_space,ky_space,kz_space,power,norm_k,Pk_3D,delta_k,k_space,map_3D,map_fft_3D
        if(log):
            plt.semilogx(k_space_final,pk_3D_final)
        else :
            plt.plot(k_space_final,pk_3D_final)
        plt.grid()
        plt.xlabel("Comoving wavevector in h.Mpc-1")
        plt.ylabel("Pk 3D")
        plt.savefig("Pk_3D.pdf",format = "pdf")
        return(pk_3D_final,k_space_final)






class TomographyStack(object):

    def __init__(self,pwd,map_name,catalog_name,type_catalog,property_file_stack,size_stack,name_stack,map_shape=None,map_size=None,property_file=None):
        self.pwd = pwd
        self.size_stack = size_stack
        self.name_stack = name_stack
        self.property_file_stack = property_file_stack

        self.tomographic_map = tomographic_objects.TomographicMap.init_classic(name=map_name,shape=map_shape,size=map_size,property_file=property_file)
        self.tomographic_map.read()
        self.catalog = tomographic_objects.Catalog.init_catalog_from_fits(catalog_name, type_catalog)


    def stack(self):
        stack = tomographic_objects.StackMap.init_by_tomographic_map(self.tomographic_map,self.catalog,self.size_stack,name=self.name_stack)
        stack.write_property_file(self.property_file_stack)
        stack.write()
        self.stack = stack

    # CR - to test
    def compute_mean_distance_to_los(self,pixel_name):
        pixel = tomographic_objects.Pixel.init_from_property_files(self.tomographic_map.property_file,name=pixel_name)
        pixel.read()
        x,y,z = pixel.repack_by_los()
        diffz=np.zeros(len(self.catalog.coord.shape[0]))
        for i in range(len(diffz)):
            arg_array  = np.argwhere((self.catalog.coord[i,0] == x[:])&(self.catalog.coord[i,0] == y[:]))
            if(len(arg_array)!=0):
                arg = arg_array[0][0]
                diffz[i] = self.catalog.coord[i,2] - z[arg][-1]
        return(np.mean(diffz))




    #
    # @staticmethod
    # def plot_stack(self,stack,name,deltamin,deltamax,los_quasar=None):
    #
    #     for direction_normale in ["x","y","z"]:
    #         plot_stack(self,stack,name,deltamin,deltamax,direction,los_quasar=None)
    #


    @staticmethod
    def plot_stack(stack,name_plot,deltamin,deltamax,direction,los_quasar=None,rotate=False,**kwargs):
        size_map = stack.size
        (x_index, y_index, index_direction, extentmap, xlab, ylab) = TomographyPlot.get_direction_informations(direction,rotate,(size_map,size_map,size_map))
        if direction == "x":
            stack_slice = stack.map_array[stack.shape[0]//2,:,:]
            if(rotate): stack_slice = np.transpose(np.flip(stack_slice,axis=1))
        elif direction == "y" :
            stack_slice = stack.map_array[:,stack.shape[1]//2,:]
            if(rotate): stack_slice = np.transpose(np.flip(stack_slice,axis=1))
        elif direction == "z" :
            stack_slice = stack.map_array[:,:,stack.shape[2]//2]
            if(rotate): stack_slice = np.transpose(np.flip(stack_slice,axis=1))
        # TomographyPlot.plot_slice(stack_slice,extentmap,xlab,ylab,name_plot,deltamin,deltamax,x_index,y_index,**kwargs)
        plt.imshow(stack_slice, interpolation='bilinear',cmap='jet_r',vmin = deltamin, vmax = deltamax,extent = extentmap)

        # plt.figure()
        # if(direction == "x"):
        #     Slice = stack[stack.shape[0]//2,:,:]
        #     xlab = "Comoving distance in the z direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
        #     ylab = "Comoving distance in the y direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
        # elif(direction == "y"):
        #     Slice = stack[:,stack.shape[1]//2,:]
        #     xlab = "Comoving distance in the z direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
        #     ylab = "Comoving distance in the x direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
        # elif(direction == "z"):
        #     Slice = stack[:,:,stack.shape[2]//2]
        #     xlab = "Comoving distance in the y direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
        #     ylab = "Comoving distance in the x direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
        # extentmap = [-stack.size,+stack.size,-stack.size,+stack.size]
        # plt.imshow(Slice, interpolation='bilinear',cmap='jet_r',vmin = deltamin, vmax = deltamax,extent = extentmap)
        # plt.plot([0],[0],'k*')
        # plt.xlabel(xlab)
        # plt.ylabel(ylab)
        # plt.colorbar()




        # if(los_quasar is not None):
        #     if (los_quasar < stack.size):
        #         if((direction_normale=="x")|(direction_normale=="y")):
        #             plt.plot([-stack.size,-los_quasar],[0,0],"r-")
        # plt.savefig("stack_"+ name +"{}.pdf".format(direction_normale),format="pdf")


    ##### ELLIPTICITY ####



    def compute_stack_ellipticity(self,stack,Mpcperpixels,deltamin,deltamax,name,size_stack,shape_stack,signe=1,ncont=4,ticks=True,length_between_los_and_quasar=None):
        elipticity = {}
        for direction in ["x","y","z"]:
            plt.figure()
            if(direction == "x"):
                Slice = np.transpose(stack[stack.shape[0]//2,:,:])
                xlab = "Comoving distance in the y direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
                ylab = "Comoving distance in the z direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
            elif(direction == "y"):
                Slice = np.transpose(stack[:,stack.shape[1]//2,:])
                xlab = "Comoving distance in the x direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
                ylab = "Comoving distance in the z direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
            elif(direction == "z"):
                Slice = np.transpose(stack[:,:,stack.shape[2]//2])
                xlab = "Comoving distance in the x direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
                ylab = "Comoving distance in the y direction [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]"
            extentmap = [-size_stack,+size_stack,-size_stack,+size_stack]
            binsx,binsy = Slice.shape[0],Slice.shape[1]
            gauss = utils.gaussian_fitter_2d(signe * Slice)
            p,success = gauss.FitGauss2D()
            x,y=np.indices((binsx,binsy),dtype=np.float)
            if(direction == "x"):
                angle = p[5]
                if((angle <45)&(angle>-45)):
                    sigmay,sigmaz,iy,iz = p[4],p[3],1,2
                elif((angle <135)&(angle>-135)):
                    sigmay,sigmaz,iy,iz = p[3],p[4],2,1
                else :
                    sigmay,sigmaz,iy,iz = p[4],p[3],1,2
                sigmay = (np.cos(np.radians(angle))**2 * Mpcperpixels[iy] + (np.sin(np.radians(angle))**2)* Mpcperpixels[iz])*sigmay
                sigmaz = (np.cos(np.radians(angle))**2 * Mpcperpixels[iz] + (np.sin(np.radians(angle))**2)* Mpcperpixels[iy])*sigmaz
                if((angle <45)&(angle>-45)):
                    p[4],p[3] = sigmay,sigmaz
                elif((angle <135)&(angle>-135)):
                    p[3],p[4] = sigmay,sigmaz
                else :
                    p[4],p[3] = sigmay,sigmaz
                xcenter,ycenter = (x - shape_stack[2]//2) * Mpcperpixels[2],-(y - shape_stack[1]//2) * Mpcperpixels[1]
                elipticity["x direction, y over z"] = sigmay/sigmaz
                elname ="sigmay/sigmaz = " + str(np.round(sigmay/sigmaz,2))
            elif(direction == "y"):
                angle = p[5]
                if((angle <45)&(angle>-45)):
                    sigmax,sigmaz,ix,iz = p[4],p[3],0,2
                elif((angle <135)&(angle>-135)):
                    sigmax,sigmaz,ix,iz = p[3],p[4],2,0
                else :
                    sigmax,sigmaz,ix,iz = p[4],p[3],0,2
                sigmax = (np.cos(np.radians(angle))**2 * Mpcperpixels[ix] + (np.sin(np.radians(angle))**2)* Mpcperpixels[iz])*sigmax
                sigmaz = (np.cos(np.radians(angle))**2 * Mpcperpixels[iz] + (np.sin(np.radians(angle))**2)* Mpcperpixels[ix])*sigmaz
                if((angle <45)&(angle>-45)):
                    p[4],p[3] = sigmax,sigmaz
                elif((angle <135)&(angle>-135)):
                    p[3],p[4] = sigmax,sigmaz
                else :
                    p[4],p[3] = sigmax,sigmaz
                xcenter,ycenter = (x - shape_stack[0]//2) * Mpcperpixels[0],-(y - shape_stack[2]//2) * Mpcperpixels[2]
                elipticity["y direction, x over z"] = sigmax/sigmaz
                elname ="sigmax/sigmaz = " + str(np.round(sigmax/sigmaz,2))
            elif(direction == "z"):
                angle = p[5]
                if((angle <45)&(angle>-45)):
                    sigmax,sigmay,ix,iy = p[4],p[3],0,1
                elif((angle <135)&(angle>-135)):
                    sigmax,sigmay,ix,iy = p[3],p[4],1,0
                else :
                    sigmax,sigmay,ix,iy = p[4],p[3],0,1
                sigmax = (np.cos(np.radians(angle))**2 * Mpcperpixels[ix] + (np.sin(np.radians(angle))**2)* Mpcperpixels[iy])*sigmax
                sigmay = (np.cos(np.radians(angle))**2 * Mpcperpixels[iy] + (np.sin(np.radians(angle))**2)* Mpcperpixels[ix])*sigmay
                if((angle <45)&(angle>-45)):
                    p[4],p[3] = sigmax,sigmay
                elif((angle <135)&(angle>-135)):
                    p[3],p[4] = sigmax,sigmay
                else :
                    p[4],p[3] = sigmax,sigmay
                xcenter,ycenter = (x - shape_stack[0]//2) * Mpcperpixels[0], -(y - shape_stack[1]//2) * Mpcperpixels[1]
                elipticity["z direction, x over y"] = sigmax/sigmay
                elname ="sigmax/sigmay = " + str(np.round(sigmax/sigmay,2))
            gaussian = gauss.Gaussian2D(*p)
            data_fitted = gaussian(y,x)
            if(ncont is not None): plt.contour(xcenter,ycenter,data_fitted.reshape(binsx, binsy),ncont,linewidths=2, colors='w')
            plt.imshow(Slice, interpolation='bilinear',cmap='jet_r',vmin = deltamin, vmax = deltamax,extent=extentmap)
            cbar = plt.colorbar()
            cbar.set_label("Flux contrast " + r"$\delta_{Fmap}$")
            if(length_between_los_and_quasar is not None):
                if (length_between_los_and_quasar < size_stack):
                    if((direction=="x")|(direction=="y")):
                        plt.plot([0,0],[-size_stack,-length_between_los_and_quasar],"r-")
            if(ticks): plt.text(-0.9*size_stack,0.9*size_stack,elname)
            plt.xlabel(xlab)
            plt.ylabel(ylab)
            plt.savefig(name + "_gaussian_fit_direction_" + direction + ".pdf",format="pdf")
        return(elipticity,p)




    def create_and_plot_stack_voids_or_cluster(self,mapsize,clusters,size_stack,deltamin,deltamax,nameout,normalized=None,number=None,coordinate_correction=None):
        map_3d = self.readClamatoMapFile()
        dict_void_or_cluster =pickle.load(open(clusters,"rb"))
        coord = dict_void_or_cluster["coord"]
        radius = dict_void_or_cluster["radius"]
        if(coordinate_correction is not None):
            Om,maxlist_name = coordinate_correction
            self.initialize_coordinates_conversion(Om,maxlist_name)
            stack = self.stack_voids_correction_xyztildestack(map_3d,mapsize,coord,size_stack)
        else :
            stack = self.stack_voids(map_3d,mapsize,coord,radius,size_stack,normalized=normalized,number=number)
        self.plot_stack(stack,nameout,deltamin,deltamax,size_stack)
        self.save_a_stack(stack,"stack_{}.bin".format(nameout))
        print("shape of stack {}".format(stack.shape))


    def create_and_plot_stack_quasars(self,mapsize,zqso_file,size_stack,deltamin,deltamax,nameout,coordinate_correction=None):
        qso = np.transpose(np.asarray(pickle.load(open(zqso_file,"rb"))))
        map_3d = self.readClamatoMapFile()
        if(coordinate_correction is not None):
            Om,maxlist_name = coordinate_correction
            self.initialize_coordinates_conversion(Om,maxlist_name)
            stack = self.stack_voids_correction_xyztildestack(map_3d,mapsize,qso,size_stack)
            diffz = None
        else :
            stack,diffz = self.stack_quasars(map_3d,mapsize,qso,size_stack)
        self.plot_stack(stack,nameout,deltamin,deltamax,size_stack,los_quasar=diffz)
        self.save_a_stack(stack,"stack_{}.bin".format(nameout),los_z = diffz)
        return(stack)

    def create_and_plot_stack_galaxy(self,mapsize,galaxy_name_file,galaxy_type_file,galaxy_zmin,galaxy_zmax,galaxy_confidence,galaxy_std,galaxy_omega_m,galaxy_mag_max,size_stack,deltamin,deltamax,nameout):
        import Picca
        map_3d = self.readClamatoMapFile()
        galaxy = Picca.Treat_HSC(self.pwd,galaxy_name_file,galaxy_type_file,galaxy_zmin,galaxy_zmax,galaxy_confidence,galaxy_std,galaxy_omega_m)
        dict_galaxy = galaxy.get_galaxy(maxlist="list_of_maximums_of_data_cube.pickle",mag_max=galaxy_mag_max)
        stack = self.stack_galaxy(map_3d,mapsize,dict_galaxy,size_stack)
        self.plot_stack(stack,nameout,deltamin,deltamax,size_stack)
        self.save_a_stack(stack,"stack_{}.bin".format(nameout))


    def merge_and_compute_ellipticity(self,list_stacks,sizemap,name_stack_merged,deltamin,deltamax,size_stack,shape_stack,signe = 1, ncont = 4,ticks=True,length_between_los_and_quasar=None):
        Mpcperpixels = np.array(sizemap)/(np.array(self.shapeMap)-1)
        stack = self.merge_stacks(list_stacks,shape_stack,name_stack_merged)
        nameout = name_stack_merged.split(".bin")[0]
        (elipticity,p) = self.compute_stack_ellipticity(stack,Mpcperpixels,deltamin,deltamax,nameout,size_stack,shape_stack,signe=signe,ncont=ncont,ticks=ticks,length_between_los_and_quasar=length_between_los_and_quasar)
        return(elipticity)


    def compute_jack_knife_ellipticity(self,snap,list_stacks,list_stacks_snap,sizemap,name_stack_merged,name_mean_stack_snap,deltamin,deltamax,size_stack,shape_stack,signe = 1, ncont = 4,ticks=True,length_between_los_and_quasar=None):
        ellipticities = {}
        for i in range(len(snap)):
            list_stacks_knife = list_stacks_snap.copy()
            list_stacks_knife.remove(list_stacks_snap[i])
            (elipticity) = self.merge_and_compute_ellipticity(list_stacks_knife,sizemap,name_mean_stack_snap[i],deltamin,deltamax,size_stack,shape_stack,signe = signe, ncont = ncont)
            ellipticities[snap[i]] = elipticity
        ellipticities["full"] = self.merge_and_compute_ellipticity(list_stacks,sizemap,name_stack_merged,deltamin,deltamax,size_stack,shape_stack,signe = signe, ncont = ncont)
        sigma ={}
        for el in ["x direction, y over z","y direction, x over z","z direction, x over y"]:
            mean = ellipticities["full"][el]
            Sum = 0
            for s in snap :
                Sum = Sum + (ellipticities[s][el] - mean)**2
            sigma[el] = np.sqrt(Sum/(len(snap)*(len(snap)-1)))
        return(ellipticities,sigma)
