import numpy as np
import pandas as pd
import subprocess as sbp
import os
import utils
from pandas_ods_reader import read_ods
from sys import exit
from math import sqrt
# from main import trial_tag
hartree_to_kcal = 627.503
kcal_to_hartree = 1/hartree_to_kcal
cm_to_kcal = 1.0/349.757
trial_tag = 'BLYP_Ex_Ec_FA_AE_IE_PA'
class energyDecomp:
    def __init__(self):
        self.E_tot = 0
        self.E_1e = 0
        self.E_coul = 0
        self.E_x = 0
        self.E_c = 0
        self.E_nuc = 0
        self.E_FA = 0
        self.E_zp = 0
    
    def setValue(self,part,value):
        if part == 'E_tot':
            self.E_tot = value
        elif part == 'E_1e':
            self.E_1e = value
        elif part == 'E_coul':
            self.E_coul = value
        elif part == 'E_x':
            self.E_x = value
        elif part == 'E_c':
            self.E_c = value
        elif part == 'E_nuc':
            self.E_nuc = value
        elif part == 'E_FA':
            self.E_FA = value
        elif part == 'E_zp':
            self.E_zp = value

    def print(self):
        print("""
        - Total DFT energy         = %.6f
        - One electron energy      = %.6f
        - Coulomb energy           = %.6f
        - Exchange energy          = %.6f
        - Correlation energy       = %.6f
        - Nuclear repulsion energy = %.6f
        - Fermi-Amaldi energy      = %.6f
        - Zero-point correction    = %.6f
        """ % (self.E_tot, self.E_1e, self.E_coul, self.E_x, self.E_c, self.E_nuc, self.E_FA, self.E_zp))
    
    def plus(self,atom):
        self.E_tot  += atom.E_tot
        self.E_1e   += atom.E_1e
        self.E_coul += atom.E_coul
        self.E_x    += atom.E_x
        self.E_c    += atom.E_c
        self.E_nuc  += atom.E_nuc
        self.E_FA   += atom.E_FA
        self.E_zp   += atom.E_zp
        return self
    
    def minus(self,atom):
        self.E_tot  -= atom.E_tot
        self.E_1e   -= atom.E_1e
        self.E_coul -= atom.E_coul
        self.E_x    -= atom.E_x
        self.E_c    -= atom.E_c
        self.E_nuc  -= atom.E_nuc
        self.E_FA   -= atom.E_FA
        self.E_zp   -= atom.E_zp
        return self

    def scale2kcalmol(self):
        self.E_tot  *= hartree_to_kcal
        self.E_1e   *= hartree_to_kcal
        self.E_coul *= hartree_to_kcal
        self.E_x    *= hartree_to_kcal
        self.E_c    *= hartree_to_kcal
        self.E_nuc  *= hartree_to_kcal
        self.E_FA   *= hartree_to_kcal
        self.E_zp   *= hartree_to_kcal
        return self

    def mult(self,value):
        self.E_tot  *= value
        self.E_1e   *= value
        self.E_coul *= value
        self.E_x    *= value
        self.E_c    *= value
        self.E_nuc  *= value
        self.E_FA   *= value
        self.E_zp   *= value
        return self
    
    def energyCheck(self):
        print("E_tot - \sum E_part = %.6f" % (self.E_tot - self.E_1e - self.E_coul - self.E_x - self.E_c - self.E_nuc))

def read_energies(s):

    Energy = energyDecomp()
    if s == 'H+':
        return None
    os.chdir(s)
    out_file = s + '.out'
    if os.path.exists(out_file): # first make sure the output file exist
        pp=sbp.run('tail -n 1 {0}'.format(out_file),capture_output=True, text=True,shell=True)
        words_list = pp.stdout.split() 
        if words_list[:3] == ['Total', 'times', 'cpu:']: #verify the out file ends up with 'Total times cpu: xxx'
            cmd = 'cat {0}.out | grep "Total DFT energy" | head -n 3 | tail -n 1'.format(s)
            p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
            Energy.setValue('E_tot',float(p1.stdout.split()[-1]))

            cmd = 'cat {0}.out | grep "One electron energy" | head -n 3 | tail -n 1'.format(s)
            p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
            Energy.E_1e = float(p1.stdout.split()[-1])

            cmd = 'cat {0}.out | grep "Coulomb energy" | head -n 3 | tail -n 1'.format(s)
            p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
            Energy.setValue('E_coul', float(p1.stdout.split()[-1]))

            N = utils._electron_number_(s)
            Energy.setValue('E_FA', - Energy.E_coul / N) 

            cmd = 'cat {0}.out | grep "Exchange energy" | head -n 3 | tail -n 1'.format(s)
            p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
            Energy.setValue('E_x', float(p1.stdout.split()[-1]))

            cmd = 'cat {0}.out | grep "Correlation energy" | head -n 3 | tail -n 1'.format(s)
            p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
            Energy.setValue('E_c', float(p1.stdout.split()[-1]))

            cmd = 'cat {0}.out | grep "Nuclear repulsion energy" | head -n 3 | tail -n 1'.format(s)
            p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
            Energy.setValue('E_nuc', float(p1.stdout.split()[-1]))

            if utils._atom_number_(s) > 1:
                cmd = 'cat -n {0}.out | grep "Zero-Point correction to Energy" | tail -1 '.format(s)
                p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
                Energy.setValue('E_zp', float(p1.stdout.split()[-2]))
            #print(calc_energy)
        else:
            print('Failed system: ',s)
            os.chdir('..')
            return None
    else:
        print('Failed system: ',s)
        os.chdir('..')
        return None
    os.chdir('..')
    return Energy


def data_loading(file_name):
    if file_name not in ['Atomization energy','Electron affinities','Ionization energy','Proton affinities','Total energy']:
        print("Wrong file")
        exit()

    file_path = '/home/ysy1111/MyProj/FAhybrid/Gn/doc/' + file_name + '.ods'
    df = read_ods(file_path, 1)
    return df

def read_total_energy(s):
    if s == 'H+':
        return 0,0
    os.chdir(s)
    out_file = s + '.out'
    if os.path.exists(out_file): # first make sure the output file exist
        pp=sbp.run('tail -n 1 {0}'.format(out_file),capture_output=True, text=True,shell=True)
        words_list = pp.stdout.split() 
        if words_list[:3] == ['Total', 'times', 'cpu:']: #verify the out file ends up with 'Total times cpu: xxx'
            cmd = 'cat {0}.out | grep "Total DFT energy" | head -n 3 | tail -n 1'.format(s)
            p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
            
            calc_energy = float(p1.stdout.split()[-1])
            zpe = 0
            if utils._atom_number_(s) > 1:
                cmd = 'cat -n {0}.out | grep "Zero-Point correction to Energy" | tail -1 '.format(s)
                p1 = sbp.run(cmd, capture_output=True, text=True,shell=True)
                zpe = np.float(p1.stdout.split()[-2])
            #print(calc_energy)
        else:
            print('Failed system: ',s)
            os.chdir('..')
            return None,None
    else:
        print('Failed system: ',s)
        os.chdir('..')
        return None,None
    os.chdir('..')
    return calc_energy,zpe

def total_energy_error(data_path,verbose):
    

    #get_the_info
    df = data_loading('Total energy')
    #print(df)
    SYS = np.array(list(df['SYSTEM']))
    expt_value = np.array(list(df["G4-energy"]))
    calc_value = np.array([])
    error_value = np.array([])
    expt_value_s = np.array([])
    sys = np.array([])
    successful_job_counter = 0
    TOTAL_ERROR = 0
    os.chdir(data_path)

    for i,s in enumerate(SYS):
        Molecule_energy,zpe = read_total_energy(s)
        if not Molecule_energy:
            continue
        # Molecule_energy *= hartree_to_kcal
        # if(abs(Molecule_energy - expt_value[i])>3):
        #     print(s,Molecule_energy,expt_value[i])
        Molecule_energy += zpe
        TOTAL_ERROR += (Molecule_energy - expt_value[i])**2
        sys = np.append(sys,s)
        calc_value = np.append(calc_value,Molecule_energy)
        error_value = np.append(error_value,Molecule_energy - expt_value[i])
        expt_value_s = np.append(expt_value_s,expt_value[i])
        successful_job_counter += 1
    if verbose == 2:
        space = {'sys':sys,'expt_value':expt_value_s,'calc_value':calc_value,'error':error_value}
        df = pd.DataFrame(space,columns = ['sys','expt_value','calc_value','error'])
        path = "/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}_totalEnergyErrorDetails.csv".format(trial_tag)
        df.to_csv(path, index = False, header=True)
    return TOTAL_ERROR*(hartree_to_kcal**2),successful_job_counter


# tpath = '/home/ysy1111/scratch/FA_scaling_job/'
# sq_error,nj = total_energy_error(tpath)
# print(nj,sqrt(sq_error/nj))


def atomization_energy_error(data_path,verbose):
    #get_the_info
    df = data_loading('Atomization energy')
    SYS = np.array(list(df['SYSTEM']))
    expt_value = np.array(list(df["expt.(kcal/mol)"]))
    # ZPE_cm = np.array(list(df["ZPE(cm-1)"]))
    # ZPE_kcal = ZPE_cm * cm_to_kcal

    calc_value = np.array([])
    sys = np.array([])
    error_value = np.array([])
    expt_value_s = np.array([])

    successful_job_counter = 0
    TOTAL_ERROR = 0
    os.chdir(data_path)

    for i,s in enumerate(SYS):
        Unreadable_atom = False
        # get the energy of atoms
        atom_list = utils._element_list_(s) # get the atom info of a molecule(ion),i.e., C2H3(+) -> {'H':3,'C':2}
        atom_energy_sum = 0
        for el in atom_list:
            single_atom_energy,zpe = read_total_energy(el)
            if not single_atom_energy:
                Unreadable_atom = True
                break #go to the next molecule
            atom_energy_sum += atom_list[el]*single_atom_energy
        if Unreadable_atom:
            continue
        atom_energy_sum *= hartree_to_kcal # Unify Unit to kcal/mol
        #get the energy of the molucule
        Molecule_energy,zpe = read_total_energy(s)
        if not Molecule_energy:
            continue
        Molecule_energy += zpe    
        Molecule_energy *= hartree_to_kcal

        #calc_atomization_energy = atom_energy_sum - Molecule_energy - ZPE_kcal[i]
        calc_atomization_energy = atom_energy_sum - Molecule_energy
        #print(s,calc_atomization_energy,expt_value[i],(calc_atomization_energy - expt_value[i]))
        TOTAL_ERROR += (calc_atomization_energy - expt_value[i])**2
        sys = np.append(sys,s)
        calc_value = np.append(calc_value,calc_atomization_energy)
        error_value = np.append(error_value,calc_atomization_energy - expt_value[i])
        expt_value_s = np.append(expt_value_s,expt_value[i])
        successful_job_counter += 1
    if verbose == 2:
        space = {'sys':sys,'expt_value':expt_value_s,'calc_value':calc_value,'error':error_value}
        df = pd.DataFrame(space,columns = ['sys','expt_value','calc_value','error'])
        path = "/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}_AtomizationEnergyErrorDetails.csv".format(trial_tag)
        df.to_csv(path, index = False, header=True)
    return TOTAL_ERROR,successful_job_counter

# tpath = '/home/ysy1111/scratch/FA_scaling_job/'
# sq_error,nj = atomization_energy_error(tpath)
# print(nj,sqrt(sq_error/nj))        
    
def electron_affinity_error(data_path,verbose):
    #get_the_info
    df = data_loading('Electron affinities')
    #print(df)
    SYS1 = np.array(list(df['SYSTEM1']))
    SYS2 = np.array(list(df['SYSTEM2']))
    expt_value = -np.array(list(df["expt.(kcal/mol)"]))
    ZPE1_cm = np.array(list(df["ZPE1(cm-1)"]))
    ZPE1_kcal = ZPE1_cm * cm_to_kcal
    ZPE2_cm = np.array(list(df["ZPE2(cm-1)"]))
    ZPE2_kcal = ZPE2_cm * cm_to_kcal

    calc_value = np.array([])
    error_value = np.array([])
    sys = np.array([])
    expt_value_s = np.array([])

    successful_job_counter = 0
    TOTAL_ERROR = 0
    os.chdir(data_path)
    for i,s in enumerate(SYS1):
        Molecule_energy_1,zpe1 = read_total_energy(s)
        if not Molecule_energy_1:
            continue
        Molecule_energy_2,zpe2 = read_total_energy(SYS2[i])
        if not Molecule_energy_2:
            continue
        Molecule_energy_1 += zpe1   
        Molecule_energy_2 += zpe2 
        Molecule_energy_1 *= hartree_to_kcal
        Molecule_energy_2 *= hartree_to_kcal
        # ZPE_1 = ZPE1_kcal[i]
        # ZPE_2 = ZPE2_kcal[i]
        # calc_pa =  -(Molecule_energy_1 + ZPE_1 - Molecule_energy_2 - ZPE_2)
        calc_pa =  -(Molecule_energy_1 - Molecule_energy_2)
        TOTAL_ERROR += (expt_value[i] - calc_pa)**2
        successful_job_counter += 1
        sys = np.append(sys,s)
        calc_value = np.append(calc_value,calc_pa)
        error_value = np.append(error_value,calc_pa - expt_value[i])
        expt_value_s = np.append(expt_value_s,expt_value[i])
        #print(s,calc_pa,expt_value[i])
        if verbose == 2:
            space = {'sys':sys,'expt_value':expt_value_s,'calc_value':calc_value,'error':error_value}
            df = pd.DataFrame(space,columns = ['sys','expt_value','calc_value','error'])
            path = "/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}_ElectronAffinityErrorDetails.csv".format(trial_tag)
            df.to_csv(path, index = False, header=True)
    return TOTAL_ERROR,successful_job_counter

# tpath = '/home/ysy1111/scratch/FA_scaling_job/'
# sq_error,nj = electron_affinity_error(tpath)    
# print(nj,sqrt(sq_error/nj)) 

def proton_affinity_error(data_path,verbose):
    #get_the_info
    df = data_loading('Proton affinities')
    #print(df)
    SYS1 = np.array(list(df['SYSTEM1']))
    SYS2 = np.array(list(df['SYSTEM2']))
    expt_value = np.array(list(df["expt.(kcal/mol)"]))
    ZPE1_cm = np.array(list(df["ZPE1(cm-1)"]))
    ZPE1_kcal = ZPE1_cm * cm_to_kcal
    ZPE2_cm = np.array(list(df["ZPE2(cm-1)"]))
    ZPE2_kcal = ZPE2_cm * cm_to_kcal

    calc_value = np.array([])
    error_value = np.array([])
    sys = np.array([])
    expt_value_s = np.array([])

    successful_job_counter = 0
    TOTAL_ERROR = 0
    os.chdir(data_path)
    for i,s in enumerate(SYS1):
        Molecule_energy_1,zpe1 = read_total_energy(s)
        if not Molecule_energy_1:
            continue
        Molecule_energy_2,zpe2 = read_total_energy(SYS2[i])
        if not Molecule_energy_2:
            continue
        Molecule_energy_1 += zpe1   
        Molecule_energy_2 += zpe2 
        Molecule_energy_1 *= hartree_to_kcal
        Molecule_energy_2 *= hartree_to_kcal
        # ZPE_1 = ZPE1_kcal[i]
        # ZPE_2 = ZPE2_kcal[i]
        # calc_pa =  (Molecule_energy_1 + ZPE_1 - Molecule_energy_2 - ZPE_2)
        calc_pa =  (Molecule_energy_1 - Molecule_energy_2)
        TOTAL_ERROR += (expt_value[i] - calc_pa)**2
        sys = np.append(sys,s)
        calc_value = np.append(calc_value,calc_pa)
        error_value = np.append(error_value,calc_pa - expt_value[i])
        expt_value_s = np.append(expt_value_s,expt_value[i])
        successful_job_counter += 1
    
    if verbose == 2:
        space = {'sys':sys,'expt_value':expt_value_s,'calc_value':calc_value,'error':error_value}
        df = pd.DataFrame(space,columns = ['sys','expt_value','calc_value','error'])
        path = "/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}_ProtonAffinityErrorDetails.csv".format(trial_tag)
        df.to_csv(path, index = False, header=True)
    return TOTAL_ERROR,successful_job_counter
    


# tpath = '/home/ysy1111/scratch/FA_scaling_job/'
# sq_error,nj = proton_affinity_error(tpath) 
# print(nj,sqrt(sq_error/nj)) 
def ionization_energy_error(data_path,verbose):
    #get_the_info
    df = data_loading('Ionization energy')
    #print(df)
    SYS1 = np.array(list(df['SYSTEM1']))
    SYS2 = np.array(list(df['SYSTEM2']))
    expt_value = np.array(list(df["expt.(kcal/mol)"]))

    sys = np.array([])
    calc_value = np.array([])
    error_value = np.array([])
    expt_value_s = np.array([])

    successful_job_counter = 0
    TOTAL_ERROR = 0
    os.chdir(data_path)
    #=====The ignored systems are the ones that Gaussian and NWChem can not agree. 
    #=====It is believed that NWChem does not give the correct total energy for these systems.
    ignored_system = ['CH3F','SC','CN','C4H4O']
    for i,s in enumerate(SYS1):
        if s in ignored_system:
            continue
        Molecule_energy_1,zpe1 = read_total_energy(s)
        if not Molecule_energy_1:
            continue
        Molecule_energy_2,zpe2 = read_total_energy(SYS2[i])
        if not Molecule_energy_2:
            continue
        Molecule_energy_1 += zpe1   
        Molecule_energy_2 += zpe2 
        Molecule_energy_1 *= hartree_to_kcal
        Molecule_energy_2 *= hartree_to_kcal
        calc_pa =  -(Molecule_energy_1 - Molecule_energy_2)
        TOTAL_ERROR += (expt_value[i] - calc_pa)**2
        sys = np.append(sys,s)
        calc_value = np.append(calc_value,calc_pa)
        error_value = np.append(error_value,calc_pa - expt_value[i])
        expt_value_s = np.append(expt_value_s,expt_value[i])
        successful_job_counter += 1
        # if abs(expt_value[i] - calc_pa)>40:
        #     print(s,Molecule_energy_1*kcal_to_hartree,Molecule_energy_2*kcal_to_hartree)
        #     print(calc_pa,expt_value[i])
    if verbose == 2:
        space = {'sys':sys,'expt_value':expt_value_s,'calc_value':calc_value,'error':error_value}
        df = pd.DataFrame(space,columns = ['sys','expt_value','calc_value','error'])
        path = "/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}_ionizationEnergyErrorDetails.csv".format(trial_tag)
        df.to_csv(path, index = False, header=True)
    return TOTAL_ERROR,successful_job_counter

# tpath = '/home/ysy1111/scratch/FA_scaling_job/'
# sq_error,nj = ionization_energy_error(tpath) 
# print(nj,sqrt(sq_error/nj)) 

'''
This part of the code is used to build the matrix for linear-regression.
    E_tot -- Total DFT energy 
    E_1e -- One electron energy 
    E_coul -- Coulomb energy 
    E_x -- Exchange energy : In my implementation, FA energy is combined into Exchange energy. 
                      Thus any existing Exchange energy should be (Exchange energy - 1/N * Coulomb energy).
                      But for linear fitting, FA contribution is turned off. directly calculate it from Coulomb
    E_c -- Correlation energy 
    E_nuc -- Nuclear repulsion energy 
    E_zp -- Zero-Point correction to Energy
    E_FA -- Fermi-Amaldi energy = 1/N * Coulomb energy
 ---1. For atomization energy fitting matrix:
        \sigma_atoms (E_1e+E_coul+E_c+E_nuc) - molecule (E_1e+E_coul+E_c+E_nuc+E_zp) 
         + a0 * [\sigma_atoms ( E_x ) - molecule ( E_x )]
         + a2 * [\sigma_atoms ( E_FA ) - molecule ( E_FA )] = E_exp
        
                                                                                                 [a0]
    --> [\sigma_atoms ( E_x ) - molecule ( E_x )    \sigma_atoms ( E_FA ) - molecule ( E_FA )] * [  ] = 
                                                                                                 [a2]
        [E_exp - (\sigma_atoms (E_1e+E_coul+E_c+E_nuc) - molecule (E_1e+E_coul+E_c+E_nuc+E_zp) )]
 ---2. For ionization energy fitting matrix:
        molecule+ (E_1e+E_coul+E_c+E_nuc+E_zp)  - molecule (E_1e+E_coul+E_c+E_nuc+E_zp) 
         + a0 * [molecule+ ( E_x ) - molecule ( E_x )]
         + a2 * [molecule+ ( E_FA ) - molecule ( E_FA )] = E_exp
                                                                                            [a0]         
    ---> [molecule+ ( E_x ) - molecule ( E_x )    molecule+ ( E_FA ) - molecule ( E_FA ) ] *[  ] =
                                                                                            [a2]

         [E_exp - (molecule+ (E_1e+E_coul+E_c+E_nuc+E_zp)  - molecule (E_1e+E_coul+E_c+E_nuc+E_zp))]
 ---3. Proton affinity is the same with ionization energy                                                                                   
'''
def AE_matrix_builder(subfolder,fitted_properties = []):
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)

    df = data_loading('Atomization energy')
    SYS = np.array(list(df['SYSTEM']))
    expt_value = np.array(list(df["expt.(kcal/mol)"]))

    successful_job_counter = 0

    X = []
    y = []
    #SYS = ["O2"]
    for i,s in enumerate(SYS):
        atomEnergySum = energyDecomp()
        Unreadable_atom = False
        # get the energy of atoms
        atom_list = utils._element_list_(s) # get the atom info of a molecule(ion),i.e., C2H3(+) -> {'H':3,'C':2}
        for atom in atom_list:
            atomEnergy = read_energies(atom)
            #======
            # atomEnergy.print()
            #======
            if not atomEnergy:
                Unreadable_atom = True
                break #go to the next molecule
            atomEnergySum.plus(atomEnergy.mult(atom_list[atom]))
        if Unreadable_atom:
            continue

        moleculeEnergy = read_energies(s)
        # moleculeEnergy.print()

        if not moleculeEnergy:
            continue
        
        energiesDiff = atomEnergySum.minus(moleculeEnergy)
        #======
        # energiesDiff.print()
        # energiesDiff.energyCheck()
        #======

        energiesDiff.scale2kcalmol()
        Energy_list = {'E_1e': energiesDiff.E_1e, 'E_coul': energiesDiff.E_coul, 'E_x': energiesDiff.E_x, 'E_c': energiesDiff.E_c, 'E_nuc': energiesDiff.E_nuc, 'E_zp': energiesDiff.E_zp, 'E_FA': energiesDiff.E_FA}
        minus_FA = 0 if 'E_FA' in fitted_properties else Energy_list['E_FA']
        E_part = sum([Energy_list[e] for e in Energy_list]) - sum([Energy_list[e] for e in fitted_properties]) - minus_FA
        #print("diff = ", E_part - (energiesDiff.E_1e + energiesDiff.E_coul + energiesDiff.E_nuc + energiesDiff.E_zp))
        #currentY = expt_value[i] - (energiesDiff.E_1e + energiesDiff.E_coul + energiesDiff.E_c + energiesDiff.E_nuc + energiesDiff.E_zp)
        currentY = expt_value[i] - E_part
        y.append(currentY)
        currentX = [Energy_list[e] for e in fitted_properties]
        #currentX = [energiesDiff.E_x, energiesDiff.E_FA]
        #print(alter_X, energiesDiff.E_x, energiesDiff.E_c)
        X.append(currentX)
        #print(s, currentX)

        error = expt_value[i] - (energiesDiff.E_tot + energiesDiff.E_zp)
        #print(i, s,error)

        successful_job_counter += 1
    print("Atomization Energy Matrix Building Complete!")
    return X,y
#AE_matrix_builder('job_GGA',fitted_properties = ['E_x','E_c'])

def IE_matrix_builder(subfolder,fitted_properties=[]):
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)

    df = data_loading('Ionization energy')
    SYS1 = np.array(list(df['SYSTEM1']))
    SYS2 = np.array(list(df['SYSTEM2']))
    expt_value = np.array(list(df["expt.(kcal/mol)"]))
    # print(SYS1)
    # print(SYS2)
    successful_job_counter = 0

    X = []
    y = []
    ignored_system = ['CH3F','SC','CN','C4H4O']
    # SYS1 = ["Li"]
    # SYS2 = ["Li+"]
    for i in range(len(SYS1)):
        if SYS1[i] in ignored_system:
            continue
        moleculeEnergy1 = read_energies(SYS1[i])
        # #===
        # moleculeEnergy1.print()
        # #===
        if not moleculeEnergy1:
            continue
        moleculeEnergy2 = read_energies(SYS2[i])
        # #===
        # moleculeEnergy2.print()
        # #===
        if not moleculeEnergy2:
            continue
        
        energiesDiff = moleculeEnergy2.minus(moleculeEnergy1)
        # # ======
        # energiesDiff.print()
        # # ======

        energiesDiff.scale2kcalmol()
        Energy_list = {'E_1e': energiesDiff.E_1e, 'E_coul': energiesDiff.E_coul, 'E_x': energiesDiff.E_x, 'E_c': energiesDiff.E_c, 'E_nuc': energiesDiff.E_nuc, 'E_zp': energiesDiff.E_zp, 'E_FA': energiesDiff.E_FA}
        minus_FA = 0 if 'E_FA' in fitted_properties else Energy_list['E_FA']
        E_part = sum([Energy_list[e] for e in Energy_list]) - sum([Energy_list[e] for e in fitted_properties]) - minus_FA
        #currentY = expt_value[i] - (energiesDiff.E_1e + energiesDiff.E_coul + energiesDiff.E_c + energiesDiff.E_nuc + energiesDiff.E_zp)
        currentY = expt_value[i] - E_part
        y.append(currentY)
        
        #currentX = [energiesDiff.E_x, energiesDiff.E_FA]
        #print(s, currentX)
        currentX = [Energy_list[e] for e in fitted_properties]
        X.append(currentX)

        error = expt_value[i] - (energiesDiff.E_tot + energiesDiff.E_zp)
        #print(i, SYS1[i], SYS2[i], error)

        successful_job_counter += 1
    print("Ionization Energy Matrix Building Complete!")
    return X,y

def PA_matrix_builder(subfolder,fitted_properties=[]):
    tpath = '/home/ysy1111/scratch/' + subfolder + '/'
    os.chdir(tpath)

    df = data_loading('Proton affinities')
    SYS1 = np.array(list(df['SYSTEM1']))
    SYS2 = np.array(list(df['SYSTEM2']))
    expt_value = np.array(list(df["expt.(kcal/mol)"]))
    # print(SYS1)
    # print(SYS2)
    successful_job_counter = 0

    X = []
    y = []
    # SYS1 = ["Li"]
    # SYS2 = ["Li+"]
    ignored_system = []
    for i in range(len(SYS1)):
        if SYS1[i] in ignored_system:
            continue
        moleculeEnergy1 = read_energies(SYS1[i])
        # #===
        # moleculeEnergy1.print()
        # #===
        if not moleculeEnergy1:
            continue
        moleculeEnergy2 = read_energies(SYS2[i])
        # #===
        # moleculeEnergy2.print()
        # #===
        if not moleculeEnergy2:
            continue
        
        energiesDiff = moleculeEnergy1.minus(moleculeEnergy2)
        # # ======
        # energiesDiff.print()
        # # ======

        energiesDiff.scale2kcalmol()
        Energy_list = {'E_1e': energiesDiff.E_1e, 'E_coul': energiesDiff.E_coul, 'E_x': energiesDiff.E_x, 'E_c': energiesDiff.E_c, 'E_nuc': energiesDiff.E_nuc, 'E_zp': energiesDiff.E_zp, 'E_FA': energiesDiff.E_FA}
        minus_FA = 0 if 'E_FA' in fitted_properties else Energy_list['E_FA']
        E_part = sum([Energy_list[e] for e in Energy_list]) - sum([Energy_list[e] for e in fitted_properties]) - minus_FA
        #currentY = expt_value[i] - (energiesDiff.E_1e + energiesDiff.E_coul + energiesDiff.E_c + energiesDiff.E_nuc + energiesDiff.E_zp)
        currentY = expt_value[i] - E_part
        y.append(currentY)
        currentX = [Energy_list[e] for e in fitted_properties]
        #currentX = [energiesDiff.E_x, energiesDiff.E_FA]
        #print(s, currentX)
        X.append(currentX)
        
        #error = expt_value[i] - (energiesDiff.E_tot + energiesDiff.E_zp)
        #print(i, SYS1[i], SYS2[i], error)

        successful_job_counter += 1
    print("Proton Affinity Matrix Building Complete!")
    return X,y

# PA_matrix_builder('job_GGA')