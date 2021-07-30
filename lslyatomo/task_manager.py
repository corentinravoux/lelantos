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


import os,time,pickle,shutil
from subprocess import call
from lslyatomo import utils
import multiprocessing as mp
from functools import partial



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


    def load_cluster_optional_arguments(self):
        N = utils.return_key(self.kwargs,"N",1)
        n = utils.return_key(self.kwargs,"n",1)
        c = utils.return_key(self.kwargs,"c",1)
        return(N,n,c)



class Nersc(Machine):

    def __init__(self,**kwargs):
        ending_str = "Execution Sum Up"
        error_str = ["srun:"]
        wait_check = True
        super(Nersc,self).__init__(ending_str,error_str,wait_check)

        self.project_name = "desi"
        self.launcher_name = "Tomography_start.sl"
        self.out_file_name = "{}.out"
        self.error_file_name = "{}.err"
        self.mem_per_cpu = 1952
        self.launch_mode = utils.return_key(kwargs,"mode","separated")
        self.kwargs = kwargs


    def create_launcher(self,pwd,dir_paths,software_command_lines,software_name):
        if(self.launch_mode == "separated"):
            self.create_launcher_separated(dir_paths,software_command_lines,software_name)
        elif(self.launch_mode == "unified"):
            self.create_launcher_unified(pwd,dir_paths,software_command_lines,software_name)



    def create_launcher_separated(self,dir_paths,software_command_lines,software_name):
        queue = utils.return_key(self.kwargs,"queue","regular")
        partition = utils.return_key(self.kwargs,"partition","knl")
        (N,n,c) = self.load_cluster_optional_arguments()
        time = utils.return_key(self.kwargs,"time","06:00:00")
        for i in range(len(dir_paths)):
            f = open(os.path.join(dir_paths[i],self.launcher_name),"w")
            f.write("#!/bin/bash -l\n")
            f.write(f"#SBATCH -N {N}" + "\n")
            f.write(f"#SBATCH -n {n}" + "\n")
            f.write(f"#SBATCH -c {c}" + "\n")
            f.write(f"#SBATCH -C {partition}" + "\n")
            f.write(f"#SBATCH -q {queue}" + "\n")
            f.write(f"#SBATCH -J {software_name}" + "\n")
            f.write(f"#SBATCH -t {time}" + "\n")
            f.write("#SBATCH -L project \n")
            f.write(f"#SBATCH -A {self.project_name}" + "\n")
            f.write(f"#SBATCH -o {os.path.join(dir_paths[i],self.out_file_name.format(software_name))}" + "\n")
            f.write(f"#SBATCH -e {os.path.join(dir_paths[i],self.error_file_name.format(software_name))}" + "\n")
            f.write("\n")
            f.write(f"export OMP_NUM_THREADS={c}" + "\n")
            f.write("\n")
            f.write(f"""srun {software_command_lines[i]}"""+" \n")

    def create_launcher_unified(self,pwd,dir_paths,software_command_lines,software_name):
        queue = utils.return_key(self.kwargs,"queue","regular")
        partition = utils.return_key(self.kwargs,"partition","knl")
        (N,n,c) = self.load_cluster_optional_arguments()
        time = utils.return_key(self.kwargs,"time","06:00:00")
        f = open(os.path.join(pwd,self.launcher_name),"w")
        f.write("#!/bin/bash -l\n")
        f.write(f"#SBATCH -N {N}" + "\n")
        f.write(f"#SBATCH -C {partition}" + "\n")
        f.write(f"#SBATCH -q {queue}" + "\n")
        f.write(f"#SBATCH -J {software_name}" + "\n")
        f.write(f"#SBATCH -t {time}" + "\n")
        f.write("#SBATCH -L project \n")
        f.write(f"#SBATCH -A {self.project_name}" + "\n")
        f.write(f"#SBATCH -o {os.path.join(pwd,self.launcher_name)}.out" + "\n")
        f.write(f"#SBATCH -e {os.path.join(pwd,self.launcher_name)}.err" + "\n")
        f.write("\n")
        for i in range(len(dir_paths)):
            f.write(f"""srun -n {n} -c {c} --mem {int(c*self.mem_per_cpu)}"""
                    +f""" -o {os.path.join(dir_paths[i],self.out_file_name.format(software_name))}"""
                    +f""" -e {os.path.join(dir_paths[i],self.error_file_name.format(software_name))}"""
                    +f""" {software_command_lines[i]}"""
                    +" & \n \n")
        f.write("wait")


    def launch_unified(self,pwd):
        call(["sbatch",os.path.join(pwd,self.launcher_name)])

    def launch_separated(self,dir_paths):
        for i in range(len(dir_paths)):
            call(["sbatch",os.path.join(dir_paths[i],self.launcher_name)])


    def launch(self,dir_paths=None,pwd=None,software_command_lines=None,software_name=None):
        if(self.launch_mode == "separated"):
            self.launch_separated(dir_paths)
        elif(self.launch_mode == "unified"):
            self.launch_unified(pwd)



class Irene(Machine):

    def __init__(self,**kwargs):
        ending_str = ""
        error_str = ["srun:"]
        wait_check = True
        super(Irene,self).__init__(ending_str,error_str,wait_check)

        self.project_name = "gen12028"
        self.launcher_name = "Tomography_start.sl"
        self.out_file_name = "{}.out"
        self.error_file_name = "{}.err"
        self.kwargs = kwargs


    def create_launcher(self,pwd,dir_paths,software_command_lines,software_name):
        (N,n,c) = self.load_cluster_optional_arguments()
        partition = utils.return_key(self.kwargs,"partition","rome")
        time = utils.return_key(self.kwargs,"time","60000")
        for i in range(len(dir_paths)):
            f = open(os.path.join(dir_paths[i],self.launcher_name),"w")
            f.write("#!/bin/bash -l\n")
            f.write("\n")
            f.write(f"#MSUB -r {software_name}" + "\n")
            f.write(f"#MSUB -T {time}" + "\n")
            f.write(f"#MSUB -q {partition}" + "\n")
            f.write(f"#MSUB -o {os.path.join(dir_paths[i],self.out_file_name.format(software_name))}" + "\n")
            f.write(f"#MSUB -e {os.path.join(dir_paths[i],self.error_file_name.format(software_name))}" + "\n")
            f.write("#MSUB -m scratch,work  \n")
            f.write(f"#MSUB -N {N}" + "\n")
            f.write(f"#MSUB -n {n}" + "\n")
            f.write(f"#MSUB -c {c}" + "\n")
            f.write(f"#MSUB -A {self.project_name}" + "\n")
            f.write("\n")
            f.write("export OMP_NUM_THREADS=2\n")
            f.write("\n")
            f.write(f"ccc_mprun {software_command_lines[i]}" +"\n")
            f.close()



    def launch(self,dir_paths=None,pwd=None,software_command_lines=None,software_name=None):
        for i in range(len(dir_paths)):
            call(["ccc_msub",os.path.join(dir_paths[i],self.launcher_name)])



class Bash(Machine):

    def __init__(self,ending_str="",error_str=[],**kwargs):
        wait_check = False
        super(Bash,self).__init__(ending_str,error_str,wait_check)

        self.out_file_name = "{}.out"
        self.error_file_name = "{}.err"
        self.kwargs = kwargs


    def create_launcher(self,pwd,dir_path,command_line,software_name):
        return()


    def launch(self,dir_paths=None,pwd=None,software_command_lines=None,software_name=None):

        number_process = utils.return_key(self.kwargs,"n",1)
        if(number_process == 1):
            self.launch_serial(dir_paths=dir_paths,
                               pwd=pwd,
                               software_command_lines=software_command_lines,
                               software_name=software_name)
        else:
            self.launch_parallel(number_process,
                                 dir_paths=dir_paths,
                                 pwd=pwd,
                                 software_command_lines=software_command_lines,
                                 software_name=software_name)


    def launch_parallel(self,number_process,dir_paths=None,pwd=None,software_command_lines=None,software_name=None):
        list_launch = [[dir_paths[i],software_command_lines[i]] for i in range(len(dir_paths))]
        func = partial(self.launch_single,software_name)
        with mp.Pool(number_process) as pool:
            pool.map(func,list_launch)


    def launch_single(self,software_name,param_launch):
        dir_paths,software_command_lines = param_launch[0], param_launch[1]
        out_name = os.path.join(dir_paths,self.out_file_name.format(software_name))
        err_name = os.path.join(dir_paths,self.error_file_name.format(software_name))
        print(f"launch of {software_command_lines}")
        call(software_command_lines.split(), stdout=open(out_name,'w'), stderr=open(err_name,'w'))
        print(f"end of {software_command_lines}")


    def launch_serial(self,dir_paths=None,pwd=None,software_command_lines=None,software_name=None):
        for i in range(len(dir_paths)):
            out_name = os.path.join(dir_paths[i],self.out_file_name.format(software_name))
            err_name = os.path.join(dir_paths[i],self.error_file_name.format(software_name))
            print(f"launch of {software_command_lines[i]}")
            call(software_command_lines[i].split(), stdout=open(out_name,'w'), stderr=open(err_name,'w'))






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
        self.ending_str = "Total time"



    def create_input(self,dir_paths,launcher_params,launcher_names):
        for i in range(len(dir_paths)):

            name = os.path.join(dir_paths[i],launcher_names[i])
            (lx,ly,lz,npix,
             nx,ny,nz,sigmaf,
             lperp,lpar,
             namepixel,namemap) = (launcher_params[i]["lx"],launcher_params[i]["ly"],
                                   launcher_params[i]["lz"],launcher_params[i]["npix"],
                                   launcher_params[i]["nx"],launcher_params[i]["ny"],
                                   launcher_params[i]["nz"],launcher_params[i]["sigmaf"],
                                   launcher_params[i]["lperp"],launcher_params[i]["lpar"],
                                   launcher_params[i]["namepixel"],launcher_params[i]["namemap"])
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
            f.write("pixel_data_path = {}\n".format(os.path.join(dir_paths[i],namepixel)))
            f.write("map_path = {}\n".format(os.path.join(dir_paths[i],namemap)))
            f.close()






class TomographyManager(object):

    available_software = ("dachshund","borg")
    available_machine = ("irene","nersc","bash")

    def __init__(self,
                 pwd,
                 software,
                 machine,
                 name_pixel,
                 launch_file,
                 symlink_folder = None,
                 **kwargs):
        self.pwd = pwd
        if(symlink_folder is not None):
            self.link_tomography_folder(symlink_folder)

        self.name_pixel = name_pixel
        self.launch_file = launch_file

        self.log = utils.create_report_log(name=os.path.join(self.pwd,"Python_Report"))
        self.software = self.init_sofware(software)
        self.machine = self.init_machine(machine,**kwargs)
        if(self.machine.ending_str == "")&(self.software.ending_str is not None):
            self.machine.ending_str = self.software.ending_str



    def init_sofware(self,software):
        if(software.lower() == "dachshund"):
            return(Dachshund())
        elif(software.lower() == "borg"):
            return(Borg())
        else: return KeyError(f"The software {software} is not available, please choose in {TomographyManager.available_software}")



    def init_machine(self,machine,**kwargs):
        if(machine.lower() == "irene"):
            return(Irene(**kwargs))
        if(machine.lower() == "nersc"):
            return(Nersc(**kwargs))
        elif(machine.lower() == "bash"):
            return(Bash(**kwargs))
        else: return KeyError(f"The machine {machine} is not available, please choose in {TomographyManager.available_machine}")



    @staticmethod
    def create_python_dir(pwd):
        if (os.path.isdir(os.path.join(pwd,"Tmp")) == False):
            os.mkdir(os.path.join(pwd,"Tmp"))



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


    # CR - make it software dependent
    def copy_files(self,dir_paths,listname):
        for i in range(len(listname)):
            call(["cp",os.path.join(self.pwd,f"{self.name_pixel}_{listname[i]}"),dir_paths[i]])

    def create_machine_launcher(self,dir_paths,launcher_params):
        launcher_names = [launcher_params[i]["nameinput"] for i in range(len(launcher_params))]
        command_lines = [self.software.command_line.format(os.path.join(dir_paths[i],launcher_names[i])) for i in range(len(dir_paths))]
        self.software.create_input(dir_paths,launcher_params,launcher_names)
        self.machine.create_launcher(self.pwd,dir_paths,command_lines,self.software.name)

    def launch(self,dir_paths,launcher_params):
        launcher_names = [launcher_params[i]["nameinput"] for i in range(len(launcher_params))]
        command_lines = [self.software.command_line.format(os.path.join(dir_paths[i],launcher_names[i])) for i in range(len(dir_paths))]
        self.machine.launch(dir_paths = dir_paths, pwd=self.pwd ,software_command_lines= command_lines,software_name=self.software.name)


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
        dir_paths= [os.path.join(self.pwd,"Tmp",listname[i]) for i in range(len(listname))]
        self.copy_files(dir_paths,listname)
        self.create_machine_launcher(dir_paths,launcher_params)
        self.launch(dir_paths,launcher_params)
        for i in range(len(listname)):
            self.log.add("Launch of the input " + str(launcher_params[i]["nameinput"]))
        time.sleep(5)
        self.treat_launch(self.pwd,listname)



    def copy(self):
        launching_file = pickle.load(open(self.launch_file,"rb"))
        listname,launcher_params = launching_file[0],launching_file[1]
        for i in range(len(listname)):
            call(["cp",self.pwd + "/Tmp/" + listname[i] +"/" + launcher_params[i]["namemap"],self.pwd])


    def remove_tmp(self):
        shutil.rmtree(os.path.join(self.pwd,"Tmp"), ignore_errors=True)

    def link_tomography_folder(self,symlink_folder):
        os.makedirs(symlink_folder,exist_ok=True)
        if(os.path.isdir(self.pwd)):
            try:
                os.remove(self.pwd)
            except:
                shutil.rmtree(self.pwd)
        os.symlink(symlink_folder,self.pwd,target_is_directory =True)
