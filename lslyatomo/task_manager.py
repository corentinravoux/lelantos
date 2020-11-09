#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

Author: Corentin Ravoux

Description : Task manager developed to launch several Tomography jobs.
Tested on Cobalt and Irene clusters.
"""



#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################


import os,time,pickle
from subprocess import call
from lslyatomo import utils


#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################




class Machine(object):
    
    def __init__(self,ending_str,error_str,wait_check):
        self.ending_str = ending_str
        self.error_str = error_str
        self.wait_check = wait_check


    def is_finished(self,f):
        for i in range(len(f)):
            if f[i].strip() == self.ending_str :
                return(True)
        return(False)

    def gives_error(self,f):
        for i in range(len(f)):
            for j in range(len(self.error_str)):
                if f[i].strip().split()[0] == self.error_str[j] :
                    return(True)
        return(False)




class Nersc(Machine):
    
    def __init__(self):
        ending_str = "Execution Sum Up"
        error_str = ["srun:"]
        wait_check = True
        super(Nersc,self).__init__(ending_str,error_str,wait_check)

        self.project_name = "desi"
        self.launcher_name = "Tomography_start.sl"
        self.out_file_name = "{}.out"
        self.error_file_name = "{}.err"


    def create_launcher(self,dir_path,software_command_line,software_name):
        f = open(os.path.join(dir_path,self.launcher_name),"w")
        f.write("#!/bin/bash -l\n")
        f.write("#SBATCH -N 1\n")
        f.write("#SBATCH -C haswell\n")
        f.write("#SBATCH -q regular \n")
        f.write("#SBATCH -J {software_name}" + "\n")
        f.write("#SBATCH -t 06:00:00 \n")
        f.write("#SBATCH -L project \n")
        f.write("#SBATCH -A {self.project_name}" + "\n")
        f.write("#SBATCH -n 1\n")
        f.write("#SBATCH -c 1\n")
        f.write(f"#SBATCH -o {self.out_file_name.format(software_name)}" + "\n")
        f.write(f"#SBATCH -e {self.error_file_name.format(software_name)}" + "\n")
        f.write("\n")
        f.write(f"srun {software_command_line}" +"\n")



    def launch(self,dir_path=None,software_command_line=None,software_name=None):
        call(["sbatch",os.path.join(dir_path,self.launcher_name)])




class Irene(Machine):
    
    def __init__(self):
        ending_str = "Execution Sum Up"
        error_str = ["srun:"]
        wait_check = True
        super(Irene,self).__init__(ending_str,error_str,wait_check)

        self.project_name = "gen10586"
        self.launcher_name = "Tomography_start.sl"
        self.out_file_name = "{}.out"
        self.error_file_name = "{}.err"


    def create_launcher(self,dir_path,software_command_line,software_name):
        f = open(os.path.join(dir_path,self.launcher_name),"w")
        f.write("#!/bin/bash -l\n")
        f.write("\n")
        f.write(f"#MSUB -r {software_name}" + "\n")
        f.write("#MSUB -T 60000  \n")
        f.write("#MSUB -q skylake \n")
        f.write(f"#MSUB -o {self.out_file_name.format(software_name)}" + "\n")
        f.write(f"#MSUB -e {self.error_file_name.format(software_name)}" + "\n")
        f.write("#MSUB -m scratch,work  \n")
        f.write("#MSUB -n 1\n")
        f.write(f"#MSUB -A {self.project_name}" + "\n")
        f.write("\n")
        f.write("\n")
        f.write(f"ccc_mprun {software_command_line}" +"\n")
        f.close()



    def launch(self,dir_path=None,software_command_line=None,software_name=None):
        call(["ccc_msub",os.path.join(dir_path,self.launcher_name)])



class PersonalComputer(Machine):
    
    def __init__(self,ending_str="",error_str=[]):
        wait_check = False
        super(PersonalComputer,self).__init__(ending_str,error_str,wait_check)

        self.out_file_name = "{}.out"
        self.error_file_name = "{}.err"


    def create_launcher(self,dir_path,command_line,software_name):
        return()
        
        
    def launch(self,dir_path=None,software_command_line=None,software_name=None):
        out_name = os.path.join(dir_path,self.out_file_name.format(software_name))
        err_name = os.path.join(dir_path,self.error_file_name.format(software_name))
        call(software_command_line.split(), stdout=open(out_name,'w'), stderr=open(err_name,'w'))




class TomographySoftware(object):
    
    def __init__(self,exec_file):
        self.exec_file = exec_file
        
   

class Borg(TomographySoftware):
    
    name = "borg"
    
    def __init__(self):
        
        source_path = os.path.dirname(os.path.realpath(__file__))
        name_exec = "borg.exe"
        exec_file = os.path.join(source_path,'exec',name_exec)

        super(Borg,self).__init__(exec_file)
        self.command_line = self.exec_file + " {}"    


    def create_input(self,dir_path,launcher_params,name):
        return()



class Dachshund(TomographySoftware):
    
    name = "dachshund"
    
    def __init__(self):
        
        source_path = os.path.dirname(os.path.realpath(__file__))
        name_exec = "dachshund.exe"
        exec_file = os.path.join(source_path,'exec',name_exec)
        
        super(Dachshund,self).__init__(exec_file)
        self.command_line = self.exec_file + " {}"    



    def create_input(self,dir_path,launcher_params,name):
        lx,ly,lz,npix,nx,ny,nz,sigmaf,lperp,lpar,namepixel,namemap = launcher_params["lx"],launcher_params["ly"],launcher_params["lz"],launcher_params["npix"],launcher_params["nx"],launcher_params["ny"],launcher_params["nz"],launcher_params["sigmaf"],launcher_params["lperp"],launcher_params["lpar"],launcher_params["namepixel"],launcher_params["namemap"]
        f = open(name,"w")
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
        f.write("pixel_data_path = {}\n".format(os.path.join(dir_path,namepixel)))
        f.write("map_path = {}\n".format(os.path.join(dir_path,namemap)))
        f.close()






class TomographyManager(object):
    
    available_software = ("dachshund","borg")
    available_machine = ("irene","pc")
    
    def __init__(self,pwd,software,machine,name_pixel,launch_file,**kwargs):
        self.pwd = pwd
        self.name_pixel = name_pixel
        self.launch_file = launch_file
        
        
        self.log = utils.create_report_log(name=os.path.join(self.pwd,"Python_Report"))
        self.software = self.init_sofware(software)
        self.machine = self.init_machine(machine,kwargs=kwargs)
   


    def init_sofware(self,software):
        if(software.lower() == "dachshund"):
            return(Dachshund())
        elif(software.lower() == "borg"):
            return(Borg())
        else: return KeyError(f"The software {software} is not available, please choose in {TomographyManager.available_software}")
            
            
            
    def init_machine(self,machine,**kwargs):
        ### Deal with the kwargs
        if(machine.lower() == "irene"):
            return(Irene())
        if(machine.lower() == "nersc"):
            return(Nersc())
        elif(machine.lower() == "pc"):
            return(PersonalComputer())
        else: return KeyError(f"The machine {machine} is not available, please choose in {TomographyManager.available_machine}")



    @staticmethod
    def create_python_dir(pwd):
        if (os.path.isdir(os.path.join(pwd,"Tmp")) == False):
            os.mkdir(os.path.join(pwd,"Tmp"))
        if (os.path.isdir(os.path.join(pwd,"Python_Results")) == False):
            os.mkdir(os.path.join(pwd,"Python_Results"))
            
   
    @staticmethod
    def create_dir(pwd,dirnames):
        for i in range(len(dirnames)):
            dir_path = os.path.join(pwd,"Tmp",dirnames[i])
            if (os.path.isdir(dir_path) == False):
                os.mkdir(dir_path)

         
    @staticmethod
    def create_file(name):
        if(os.path.isfile(name) == False):
            f = open(name,"w")
            f.write("wait")
            f.close()
            




    def is_finished(self,file_lines):
        return(self.machine.is_finished(file_lines))
            
    def gives_error(self,file_lines):
        return(self.machine.gives_error(file_lines))


    def wait_until_finished(self,pwd,listname):
        list_finished = [False for i in range(len(listname))]
        out_file_name = self.machine.out_file_name.format(self.software.name)
        while(self.all_finished(list_finished) == False):
            for i in range(len(listname)):
                out_file = os.path.join(pwd,'Tmp',listname[i],out_file_name)
                TomographyManager.create_file(out_file)
                file = open(out_file,"r")
                file_lines = file.readlines()
                file.close()
                list_finished[i] = self.is_finished(file_lines)
            time.sleep(10)
        self.log.add("All calculation are finished")

    def check_errors(self,pwd,listname):
        list_error = [False for i in range(len(listname))]
        error_file_name = self.machine.error_file_name.format(self.software.name)
        self.log.add("List of errors :")
        self.log.add("")
        for i in range(len(listname)):
            error_file = os.path.join(pwd,'Tmp',listname[i],error_file_name)
            TomographyManager.create_file(error_file)
            file = open(error_file,"r")
            file_lines = file.readlines()
            file.close()
            list_error[i] = self.gives_error(file_lines)
            self.log.add(f'The file {listname[i]} gave an error : {list_error[i]}')
    

    def all_finished(self,listFinished):
        allFinished = True
        for i in range(len(listFinished)):
            if(listFinished[i]==False):
                allFinished = False
        return(allFinished)





    def create_launcher(self,dir_path,launcher_param):
        launcher_name = launcher_param["nameinput"]
        self.software.create_input(dir_path,launcher_param,os.path.join(dir_path,launcher_name))
        self.machine.create_launcher(dir_path,self.software.command_line.format(os.path.join(dir_path,launcher_name)),self.software.name)
        
    def launch(self,dir_path,launcher_param):
        launcher_name = launcher_param["nameinput"]
        self.machine.launch(dir_path = dir_path ,software_command_line= self.software.command_line.format(os.path.join(dir_path,launcher_name)),software_name=self.software.name)


    def treat_launch(self,pwd,listname):
        wait_check = self.machine.wait_check
        if(wait_check):
            self.wait_until_finished(pwd,listname)
        self.check_errors(pwd,listname)




    ### Main routine


    def launch_all(self):
        launching_file = pickle.load(open(self.launch_file,"rb")) 
        listname,launcher_params = launching_file[0],launching_file[1]
        TomographyManager.create_python_dir(self.pwd)
        TomographyManager.create_dir(self.pwd,listname)
        for i in range(len(listname)):
            dir_path = os.path.join(self.pwd,"Tmp",listname[i])
            call(["cp",os.path.join(self.pwd,f"{self.name_pixel}_{listname[i]}"),dir_path])
            self.create_launcher(dir_path,launcher_params[i])
            self.launch(dir_path,launcher_params[i])
            self.log.add("Launch of the input " + str(launcher_params[i]["nameinput"]))
        time.sleep(5)
        self.treat_launch(self.pwd,listname)

    
    
    def copy(self):
        launching_file = pickle.load(open(self.launch_file,"rb")) 
        listname,launcher_params = launching_file[0],launching_file[1]
        for i in range(len(listname)):
            call(["cp",self.pwd + "/Tmp/" + listname[i] +"/" + launcher_params[i]["namepixel"],self.pwd + "/Python_Results/"])
            call(["cp",self.pwd + "/Tmp/" + listname[i] +"/" + launcher_params[i]["namemap"],self.pwd + "/Python_Results/"])
    
    
    
    
    def launch_all_and_copy(self):
        self.launch_all()
        self.copy()

        
