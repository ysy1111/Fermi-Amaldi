import numpy as np
import pandas as pd
import subprocess as sbp
import os
import utils
from pandas_ods_reader import read_ods
from math import sqrt
from error_calculator import *
import time
from skopt.space import Real, Integer
def d3_input_generator(type = 'GGA',func = 'BLYP'):
    tpath = '/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_{type}/'.format(type = type)
    os.chdir(tpath)
    sbp.run('rm -r *',shell=True) #prepare an empty folder

    cmd = 'cp -r /home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_GGAd3/* ' + tpath
    sbp.run(cmd, shell=True)
    os.chdir(tpath)
    list1 = os.listdir('.')
    #print(list1)
    for f in list1:
        print(f)
        os.chdir(f)
        files = f+'.nw'
        cmd = "sed -i '/XC becke88 lyp/c\XC {func}' {file}".format(func = func, file=files)
        sbp.run(cmd, shell=True)
        os.chdir('..')

def input_generator(type = 'hybridGGA',func = 'B3LYP'):
    tpath = '/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_{type}/'.format(type = type)
    os.chdir(tpath)
    sbp.run('rm -r *',shell=True) #prepare an empty folder

    cmd = 'cp -r /home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_GGA/* ' + tpath
    sbp.run(cmd, shell=True)
    os.chdir(tpath)
    list1 = os.listdir('.')
    #print(list1)
    for f in list1:
        print(f)
        os.chdir(f)
        files = f+'.nw'
        cmd = "sed -i '/XC becke88 lyp/c\XC {func}' {file}".format(func = func, file=files)
        sbp.run(cmd, shell=True)
        os.chdir('..')
#d3_input_generator(type = 'hybridmGGAd3', func = 'M06-2X')
#input_generator(type = 'hybridGGA', func = 'B3LYP')

def create_sbatch_to_origin_input(hybrid = False,type = 'hybridGGA'): #file_edit_helper
    '''
    #!/bin/bash
    #SBATCH --mem-per-cpu=2G             # memory, roughly 2 times %mem defined in the input name.com file
    #SBATCH --job-name= 
    #SBATCH --time=0-00:10           # time (DD-HH:MM)
    #SBATCH --ntasks=8               # number of MPI processes
    #SBATCH --account=def-yawang
    module load StdEnv/2020
    mpirun /home/ysy1111/nwchem-7.0.2-release/bin/LINUX64/nwchem .nw > .out
#for hybrid functional:  mpirun /home/ysy1111/nwchem-origin/nwchem-7.0.2-release/bin/LINUX64/nwchem .nw > .out
    '''
    target_dir = '/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_{type}/'.format(type = type)
    src_sh = '/home/ysy1111/MyProj/FAhybrid/submit.sh'
    os.chdir(target_dir)
    list1 = os.listdir('.')
    for f in list1:
        os.chdir(f)           
        sys = f
        input_file = '{0}.nw'.format(sys)
        output_file = '{0}.out'.format(sys)
        log_file = '{0}.log'.format(sys)

        sh_name = sys +'.sh'
        if os.path.exists(sh_name):
            sbp.run('rm {0}'.format(sh_name),shell = True)
        fin = open(src_sh,'r')
        fout = open(sh_name,'w')
        for i,line in enumerate(fin):
            if i == 2:
                new_line = line.split('\n')[0] + f +'\n'
                #print(i,line,new_line)
                fout.write(new_line)
            elif i == 7:
                if hybrid:
                    new_line = 'mpirun /home/ysy1111/nwchem-origin/nwchem-7.0.2-release/bin/LINUX64/nwchem ' + input_file + ' > ' + output_file + ' 2>' + log_file
                else:
                    new_line = line + input_file + ' > ' + output_file + ' 2>' + log_file
                #print(i,line,new_line)
                fout.write(new_line)
            else:
                #print(i,line)
                fout.write(line)

        os.chdir('..')

#create_sbatch_to_origin_input(True,'hybridmGGAd3')

def clean_up(subfolder): #preprocess: CLEAN UP the job area~
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    sbp.run('rm -r *',shell=True)

def search_space(fitting_mod):
    if fitting_mod == 'XC_only': # shrink exchange the range to 1 
        space  = [Real(0.9, 1.1, name='a0'),
             Real(0, 2, name='a1')
           ]
    elif fitting_mod == 'Ex_FA':
        space  = [Real(0.8, 1.2, name='a0'),
             Real(-0.2, 0.2, name='a2')
           ]
    elif fitting_mod =='Exchange_only':
        space  = [Real(0.7, 1.3, name='a0')]
    elif fitting_mod =='Correlation_only':
        space  = [Real(0, 2, name='a1')]
    elif fitting_mod =='FA_only':
        space  = [Real(-0.2, 0.2, name='a2')]
    elif fitting_mod == 'All':
        space  = [Real(0.9, 1.1, name='a0'),
            Real(0, 2, name='a1'),
            Real(-0.1, 0.1, name='a2')
          ]
    return space

def copy_files_to_scratch(subfolder,basic_method): #STEP1: COPY all the input files to the work folder
    #print(subfolder,basic_method)
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    cmd = 'cp -r /home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_{type}/* '.format(type = basic_method) + tpath
    sbp.run(cmd,shell=True)

def edit_input_files_para(subfolder,x,basic_method,fitting_mod): #STEP1.1: CHANGE the parameter
    if fitting_mod == 'XC_only':
        a0,a1 = x
    elif fitting_mod =='Exchange_only':
        a0 = x[0]
    elif fitting_mod =='Correlation_only':
        a1 = x[0]
    elif fitting_mod =='FA_only':
        a2 = x[0]
    elif fitting_mod == 'Ex_FA':
        a0,a2 = x
    elif fitting_mod == 'All':
        a0,a1,a2 = x
        
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    list1 = os.listdir('.')

    for f in list1:
        os.chdir(f)
        files = f + '.nw'
        if basic_method == 'LDA': # hybrid with LDA
            if fitting_mod == 'XC_only' or fitting_mod == 'All':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater {a0} vwn_5 {a1}' {file}".format(file=files,a0=a0,a1=a1)
            elif fitting_mod =='Exchange_only' or fitting_mod =='Ex_FA':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater {a0} vwn_5 ' {file}".format(file=files,a0=a0)
            elif fitting_mod =='Correlation_only':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater vwn_5 {a1}' {file}".format(file=files,a1=a1)
            elif fitting_mod == 'FA_only':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater {a0} vwn_5 ' {file}".format(file=files,a0 = 1.056541) # use the correct a0 value
        elif basic_method == 'LDA_rpa': # hybrid with LDA
            if fitting_mod == 'XC_only' or fitting_mod == 'All':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater {a0} vwn_1_rpa {a1}' {file}".format(file=files,a0=a0,a1=a1)
            elif fitting_mod =='Exchange_only' or fitting_mod =='Ex_FA':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater {a0} vwn_1_rpa ' {file}".format(file=files,a0=a0)
            elif fitting_mod =='Correlation_only':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater vwn_1_rpa {a1}' {file}".format(file=files,a1=a1)
            elif fitting_mod == 'FA_only':
                cmd = "sed -i '/XC slater vwn_5/c\XC slater {a0} vwn_1_rpa ' {file}".format(file=files,a0 = 1.056541) # use the correct a0 value
        elif basic_method == 'GGA': # hybrid with GGA
            if fitting_mod == 'XC_only' or fitting_mod == 'All':
                cmd = "sed -i '/XC becke88 lyp/c\XC becke88 {a0} lyp {a1}' {file}".format(file=files,a0=a0,a1=a1)
            elif fitting_mod =='Exchange_only' or fitting_mod =='Ex_FA':
                cmd = "sed -i '/XC becke88 lyp/c\XC becke88 {a0} lyp' {file}".format(file=files,a0=a0)
            elif fitting_mod =='Correlation_only':
                cmd = "sed -i '/XC becke88 lyp/c\XC becke88 lyp {a1}' {file}".format(file=files,a1=a1)
        else:
            cmd = "" # empty command
        sbp.run(cmd, shell=True)

        os.chdir('..')

# copy_files_to_scratch('job1','LDA')    
# edit_input_files_para('job1',[1,1,1],'LDA')


def create_a2txt(x,subfolder,fitting_mod):  #STEP2: CREATE the a2.txt file for each system
    if fitting_mod == 'FA_only' or fitting_mod == 'All' or fitting_mod =='Ex_FA':
        a2 = x[-1]
    else:
        a2 = 0
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    list1 = os.listdir('.')
    for f in list1:
        os.chdir(f)
        f = open('a2.txt','w')
        f.write(str(a2))
        f.close()
        os.chdir('..')


def run_job(subfolder): #STEP2: RUN all jobs
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    list1 = os.listdir('.')
    for f in list1:
        os.chdir(f)
        sh_name = f + '.sh'
        log_file = f + '.log'
        #print(sh_name)
        cmd = 'sbatch {0} > /dev/null 2>{1}'.format(sh_name,log_file)
        sbp.run(cmd,shell=True) # submit the job
        os.chdir('..')
        
#run_job('FA_scaling_job')
def check_running_jobs(origin_l): #STEP3: check whether the all jobs have been finished
    '''
    Now check whether all the jobs have been finished or not 
    '''
    start = time.perf_counter()
    flag = True
    while(flag):
        time.sleep(60)
        #print(flag)
        p1 = sbp.run('sq',capture_output=True,shell=True)
        output = str(p1.stdout).split()
        if len(output) == origin_l:
            flag=False
    end = time.perf_counter()
    print(f'All the jobs have been finished in {round(end-start)} second(s)')

    '''
    Check whether the submitted jobs were run or not (check the existance of log and out files)
    '''


def error_compute(subfolder,basic_method,flag_t=False,flag_a=False,flag_p=False,flag_i=False,flag_e=False,verbose = 1): #STEP4: COMPUTE the error
    #ratio_LDA = [1, 0.0544, 0.0064, 0.0331, 0.0486]
    #ratio_GGA = [1, 0.1492, 0.0320, 0.3144, 0.4369]
    ratio = [1,1,1,1,1]
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    # x,nx --> sum(square error) in kcal/mol**2, number of successful job
    if flag_t:
        t,nt = total_energy_error(tpath,verbose)
    else:
        t=0
        nt=1

    if flag_a:
        a,na  = atomization_energy_error(tpath,verbose)
    else:
        a=0
        na=1

    if flag_p:
        p,np  = proton_affinity_error(tpath,verbose)
    else:
        p=0
        np=1

    if flag_e:
        e,ne  = electron_affinity_error(tpath,verbose)
    else:
        e=0
        ne=1

    if flag_i:
        i,ni  = ionization_energy_error(tpath,verbose)
    else:
        i=0
        ni=1
        
    ave_t = sqrt(t/nt) # Square-Root-Error
    ave_a = sqrt(a/na)
    ave_p = sqrt(p/np)
    ave_e = sqrt(e/ne)
    ave_i = sqrt(i/ni)
    turn_on_properties = flag_a + flag_p + flag_t + flag_i + flag_e
    # if basic_method == 'LDA' or basic_method == 'LDA_rpa':
    #     all_ave = (ave_t + ave_a/ratio[1] + ave_p/ratio_LDA[2] + ave_e/ratio_LDA[3] + ave_i/ratio_LDA[4])/turn_on_properties
    # elif basic_method == 'GGA':
    #     all_ave = (ave_t + ave_a/ratio[1] + ave_p/ratio_GGA[2] + ave_e/ratio_GGA[3] + ave_i/ratio_GGA[4])/turn_on_properties
    all_ave = (ave_t + ave_a/ratio[1] + ave_p/ratio[2] + ave_e/ratio[3] + ave_i/ratio[4])/turn_on_properties   
    print("""
    - total_energy_error=%.8f kcal/mol
    - atomization_energy_error=%.8f kcal/mol
    - proton_affinity_error=%.8f kcal/mol
    - electron_affinity_error=%.8f kcal/mol
    - ionization_energy_error=%.8f kcal/mol
    - weighted_average_error=%.8f
    """ % (ave_t, ave_a, ave_p, ave_e, ave_i, all_ave))

    return all_ave, ave_t, ave_a, ave_p, ave_e, ave_i

def details_record(x,all_ave, ave_t, ave_a, ave_p, ave_e, ave_i,trial_tag,fitting_mod):
    
    os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/callbacks')
    f = open(trial_tag,'a')
    if fitting_mod == 'XC_only':
        f.write("""parameters:
            - a0=%.6f
            - a1=%.6f
            """ % (x[0], x[1],))
    elif fitting_mod =='Exchange_only':
        f.write("""parameters:
            - a0=%.6f
            """ % (x[0]))
    elif fitting_mod =='Correlation_only':
        f.write("""parameters:
            - a1=%.6f
            """ % (x[0]))
    elif fitting_mod =='FA_only':
        f.write("""parameters:
            - a2=%.6f
            """ % (x[0]))
    elif fitting_mod =='Ex_FA':
        f.write("""parameters:
            - a0=%.6f
            - a2=%.6f
            """ % (x[0], x[1],))
    elif fitting_mod == 'All':
        f.write("""parameters:
            - a0=%.6f
            - a1=%.6f
            - a2=%.6f
            """ % (x[0], x[1], x[2]))
    f.write("""
    - total_energy_error=%.6f kcal/mol
    - atomization_energy_error=%.6f kcal/mol
    - proton_affinity_error=%.6f kcal/mol
    - electron_affinity_error=%.6f kcal/mol
    - ionization_energy_error=%.6f kcal/mol
    - ave_err=%.6f kcal/mol
    """ % (ave_t, ave_a, ave_p, ave_e, ave_i, all_ave))
    f.close()
#error_compute('FA_scaling_job','GGA')

def printFinalResult(fitting_mod,res_gp):
    if fitting_mod == 'XC_only':
        print("""Best parameters:
            - a0=%.6f
            - a1=%.6f
            """ % (res_gp.x[0], res_gp.x[1],))
    elif fitting_mod =='Exchange_only':
        print("""Best parameters:
            - a0=%.6f
            """ % (res_gp.x[0]))
    elif fitting_mod =='Correlation_only':
        print("""Best parameters:
            - a1=%.6f
            """ % (res_gp.x[0]))
    elif fitting_mod =='FA_only':
        print("""Best parameters:
            - a2=%.6f
            """ % (res_gp.x[0]))
    if fitting_mod == 'Ex_FA':
        print("""Best parameters:
            - a0=%.6f
            - a2=%.6f
            """ % (res_gp.x[0], res_gp.x[1],))
    elif fitting_mod == 'All':
        print("""Best parameters:
            - a0=%.6f
            - a1=%.6f
            - a2=%.6f
            """ % (res_gp.x[0], res_gp.x[1], res_gp.x[2]))
    print("Best score=%.12f" % res_gp.fun)


def buildTotalMatrix(subfolder,flag_t=False,flag_a=False,flag_p=False,flag_i=False,flag_e=False,fitted_energy = []):
    X = []
    y = []
    if flag_a:
        X_AE,y_AE = AE_matrix_builder(subfolder,fitted_energy)
        X += X_AE
        y += y_AE
    if flag_p: 
        X_PA,y_PA = PA_matrix_builder(subfolder,fitted_energy)
        X += X_PA
        y += y_PA
    if flag_i: 
        X_IE,y_IE = IE_matrix_builder(subfolder,fitted_energy)
        X += X_IE
        y += y_IE 
       
    return X,y