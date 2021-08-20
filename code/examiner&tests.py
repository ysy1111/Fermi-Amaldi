import numpy as np
import pandas as pd
import subprocess as sbp
import os
import utils
from pandas_ods_reader import read_ods
import matplotlib.pyplot as plt
import re
def examine_copy_from_local():
    
    os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/all_system_input_files/')
    list1 = os.listdir('.')
    print(len(list1))
    os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/doc')
    files = ['Atomization energy.ods', 'Electron affinities.ods', 'Ionization energy.ods', 'Proton affinities.ods', 'Total energy.ods']
    df1 = read_ods(files[0], 1)
    df2 = read_ods(files[1], 1)
    df3 = read_ods(files[2], 1)
    df4 = read_ods(files[3], 1)
    df5 = read_ods(files[4], 1)
    df  = pd.read_csv(r'all_system.csv')
    all_list = list(df['SYSTEM'])
    te_list = list(df5['SYSTEM'])
    at_list = list(df1['SYSTEM'])
    ea_list = list(df2['SYSTEM1']) + list(df2['SYSTEM2'])
    ie_list = list(df3['SYSTEM1']) + list(df3['SYSTEM2'])
    pa_list = list(df4['SYSTEM1']) + list(df4['SYSTEM2'])
    ans = []
    print('1')
    for s in at_list:
        if s not in list1:
            print(s)
            
    print('2')        
    for s in ea_list:
        if s not in list1:
            print(s)
            
    print('3')     
    for s in ie_list:
        if s not in list1:
            print(s)
            
    print('4')       
    for s in pa_list:
        if s not in list1:
            print(s)
            
    print('5')        
    for s in te_list:
        if s not in list1:
            print(s)
    
    for s in list1:
        if s not in all_list:
            print(s)
    os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/all_system_input_files/')
    for f in list1:
        os.chdir(f)
        files = os.listdir('.')
        if utils._atom_number_(f) == 1:
            if len(files) != 1:
                print(f)
        else:
            if len(files) != 2:
                print(f)
        os.chdir('..')

def map_plot():
    x = np.arange(1, 10)
    y = x.reshape(-1, 1)
    h = x * y

    cs = plt.contourf(h, levels=[10, 30, 50],
        colors=['#808080', '#A0A0A0', '#C0C0C0'], extend='both')
    cs.cmap.set_over('red')
    cs.cmap.set_under('blue')
    cs.changed()
    plt.show()

def energy_per_electron():
    # os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/all_system_input_files_LDA/')
    # list1 = os.listdir('.')
    # print(len(list1))
    os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/doc')
    files = ['Atomization energy.ods', 'Electron affinities.ods', 'Ionization energy.ods', 'Proton affinities.ods', 'Total energy.ods']
    df1 = read_ods(files[0], 1)
    df2 = read_ods(files[1], 1)
    df3 = read_ods(files[2], 1)
    df4 = read_ods(files[3], 1)
    df5 = read_ods(files[4], 1)
    systems = df5['SYSTEM']
    energy = df5['G4-energy']
    ave_energy = []
    for i,sys in enumerate(systems):
        en = utils._electron_number_(sys)
        ave_energy.append(energy[i]/en)
    print(max(ave_energy),min(ave_energy),len(ave_energy))
#energy_per_electron()
#map_plot()

def copy_file_for_mp2_geom_opt():
    subfolder = 'job2'
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    cmd = 'cp -r /home/ysy1111/MyProj/FAhybrid/Gn/all_system_input_files_LDA/* ' + tpath
    sbp.run(cmd,shell=True)

def search_mult(Myfile):

    with open(Myfile) as f:
        data = f.read()
        match = re.search(r'mult [0-9]', data, re.DOTALL)
        if match:
            ans = match.group().split()
            return ans[-1]
    return None

def writeFile(Myfile,s,mult):
    multTrans = {'1':'singlet','2':'doublet','3':'triplet','4':'quartet','5':'quintet'} 
    if s[-1] == '+':
        ion_mode = 1
    elif s[-1] == '-':
        ion_mode = 2
    else:
        ion_mode = 0   
    f = open(Myfile,'w')
    f.write('start ' + s + '\n')
    f.write("title '" + s + "'\n")
    if ion_mode == 1:
        f.write('charge 1.0\n')
    elif ion_mode == 2:
        f.write('charge -1.0\n')
    f.write('geometry units au \n')
    f.write('  load ' + s + '.xyz\n')
    f.write('end\n')
    f.write('basis\n')
    f.write('  * library 6-31G(2df,p) \n')
    f.write('end \n')
    f.write('scf\n')
    f.write('  uhf\n')
    f.write('  ' + multTrans[mult] + '\n')
    f.write('  maxiter 299\n')
    f.write('end\n')
    f.write('driver\n')
    f.write('  maxiter 299\n')
    f.write('  xyz ' + s + '\n')
    f.write('end\n')
    f.write('task mp2 optimize')
    f.close()
# scf
# uhf
# DOUBLET
# end
# driver
# maxiter 99
# xyz BeH-mp2
# end
# task mp2 optimize
# task mp2 frequencies 

def modify_input_for_geom_opt():
    subfolder = 'job2'
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    list1 = os.listdir('.')
    
    for f in list1:
        atomNum = utils._atom_number_(f)
        if atomNum > 1:
            os.chdir(f)
            fileOrigin = f + '.nw'
            fileNew = f + '_mp2.nw'
            mult = search_mult(fileOrigin)
            writeFile(fileNew,f,mult)
            #spb.run(cmd,shell=True,)
            os.chdir('..')

def file_process():
    subfolder = 'job2'
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    list1 = os.listdir('.')
    for f in list1:
        atomNum = utils._atom_number_(f)
        if atomNum > 1:
            os.chdir(f)
            cmd = 'rm {0}.nw'.format(f)
            sbp.run(cmd,shell=True)
            cmd = 'mv {0}_mp2.nw {0}.nw'.format(f)
            sbp.run(cmd,shell=True)
            os.chdir('..')

def modified_sh_time():
    tpath = "/home/ysy1111/MyProj/FAhybrid/Gn/all_system_input_files_LDA"
    os.chdir(tpath)
    list1 = os.listdir('.')
    for f in list1:
        atomNum = utils._atom_number_(f)
        if atomNum > 1:
            os.chdir(f)
            sh_name = f + '.sh'
            cmd = "sed -i '/#SBATCH --time=0-00:40           # time (DD-HH:MM)/c\#SBATCH --time=0-03:00           # time (DD-HH:MM)' {file}".format(file = sh_name)
            sbp.run(cmd,shell=True)
            os.chdir('..')

#modified_sh_time()


def run_geom_opt():
    subfolder = 'job2'
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    list1 = os.listdir('.')
    for f in list1:
        atomNum = utils._atom_number_(f)
        if atomNum > 1:
            os.chdir(f)
            sh_name = f + '.sh'
            log_name = f + '.log'
            cmd = 'sbatch {0} > /dev/null 2>{1}'.format(sh_name,log_name)
            sbp.run(cmd,shell=True)
            os.chdir('..')

def identify_the_failed_job():
    subfolder = 'job2'
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)
    list1 = os.listdir('.')
    failed_sys = []
    suc_sys = []
    for f in list1:
        atomNum = utils._atom_number_(f)
        if atomNum > 1:
            os.chdir(f)
            out_name = f + '.out'
            if os.path.exists(out_name): # first make sure the output file exist
                pp=sbp.run('tail -n 1 {0}'.format(out_name),capture_output=True, text=True,shell=True)
                words_list = pp.stdout.split() 
                if words_list[:3] != ['Total', 'times', 'cpu:']: #verify the out file ends up with 'Total times cpu: xxx'
                    failed_sys.append(f)
                else:
                    suc_sys.append(f)
            os.chdir('..')
    print(failed_sys)
    print(suc_sys)

def change_basis():


    target_folder = "/home/ysy1111/MyProj/FAhybrid/Gn/all_system_input_files_LDA"

    os.chdir(target_folder)
    list1 = os.listdir('.')
    for s in list1:
        file_path = './' + s + '/' + s + '.nw'
        # f = open(file_path,'w')
        cmd = "sed -i 's/.*library.*/  * library 6-311+G(2d,p) /' {file}".format(file=file_path)
        sbp.run(cmd,shell=True,capture_output=True)

def input_for_AE_IE_PA():
    tpath = '/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_GGAd3'


    os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/doc')
    files = ['Atomization energy.ods', 'Electron affinities.ods', 'Ionization energy.ods', 'Proton affinities.ods', 'Total energy.ods']
    dfAE = read_ods(files[0], 1)
    dfIE = read_ods(files[2], 1)
    dfPA = read_ods(files[3], 1)
    sys_AE = set(np.array(dfAE['SYSTEM']))
    sys_IE1 = set(np.array(dfIE['SYSTEM1']))
    sys_IE2 = set(np.array(dfIE['SYSTEM2']))
    sys_PA1 = set(np.array(dfPA['SYSTEM1']))
    sys_PA2 = set(np.array(dfPA['SYSTEM2']))
    set_Atom = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar'}
    total_set = set_Atom.union(sys_AE)
    total_set = total_set.union(sys_IE1)
    total_set = total_set.union(sys_IE2)
    total_set = total_set.union(sys_PA1)
    total_set = total_set.union(sys_PA2)
    for sys in total_set:
        cmd = 'cp -r /home/ysy1111/MyProj/FAhybrid/Gn/all_system_input_files_GGA/{file} '.format(file = sys) + tpath
        sbp.run(cmd,shell=True)

def edit_input_add_d3():
    os.chdir('/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_GGAd3')
    list1 = os.listdir('.')
    for sys in list1:
        os.chdir(sys)
        cmd = "sed -i '/decomp/a\ \ disp vdw 3' {}.nw".format(sys)
        sbp.run(cmd,shell = True)
        os.chdir('..')

def edit_input_mGGA():
    target_dir = '/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_mGGA'
    source_dir = '/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_GGAd3'
    os.chdir(target_dir)
    list1 = os.listdir('.')
    for l in list1:
        os.chdir(l)
        sys = "".join([l,'.nw'])
        cmd = "sed -i '/XC becke88 lyp/c\XC m06-L' {file}".format(file=sys)
        sbp.run(cmd, shell=True)
        os.chdir('..')

def reEditClocedShellInput(type):
    tpath = '/home/ysy1111/MyProj/FAhybrid/Gn/sub_all_system_input_files_{type}/'.format(type = type)
    os.chdir(tpath)
    list1 = os.listdir('.')
    # list1 = ['C']
    for f in list1:
        os.chdir(f)
        files = f+'.nw'
        cmd = "grep 'mult' {file}".format(file = files)
        p1 = sbp.run(cmd, capture_output = True, shell=True)
        mult = int(p1.stdout.split()[1])
        if mult == 1:
            cmd = "sed -i '/odft/d' ./{}".format(files)
            sbp.run(cmd, shell=True)
            cmd = "sed -i '/mult/d' ./{}".format(files)
            sbp.run(cmd, shell=True)
        os.chdir('..')


reEditClocedShellInput('hybridmGGAd3')
#edit_input_mGGA()
# edit_input_add_d3()

#change_basis()
# copy_file_for_mp2_geom_opt()
# modify_input_for_geom_opt()
# file_process()
# break_sym()
# run_geom_opt()