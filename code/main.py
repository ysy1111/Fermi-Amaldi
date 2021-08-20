import numpy as np
import pandas as pd
import subprocess as sbp
from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args
import utils
from skopt.plots import plot_convergence
import os
import skopt
from skopt import load
import matplotlib.pyplot as plt
from tools import *
import time
import pickle
from sklearn.linear_model import LinearRegression
#job id#
basic_method = 'GGA'  #LDA/GGA/LDA_rpa  GGA contains GGA-d3 (check tools copy_files_to_scratch)
subfolder = 'job1' # job1  job2  job_GGA  job_GGAd3  job_hybridmGGAd3  job_LDA  job_mGGAd3
#restored_tag = 'Gn_all_pure_LDA_repara'
#'BLYP_Ex1.01743487e,EFA5.35548303e-04_Ex_FA_AE_IE_PA'
trial_tag = 'BLYP_Ex_Ec_FA_AE_IE_PA' #label the fitting job 
initial_points = [[1,1,0]] # only for start()
fitting_mod = ['Exchange_only','Correlation_only','FA_only','XC_only','Ex_FA','All']
fm = 'All' # control the set of parameters that need to be fitted

p_origin = sbp.run('sq -u ysy1111',capture_output=True,shell=True)
output = str(p_origin.stdout).split()
origin_l = len(output)    

def f(x):
    print(x)
    clean_up(subfolder)
    copy_files_to_scratch(subfolder, basic_method)
    edit_input_files_para(subfolder, x, basic_method, fm)
    create_a2txt(x, subfolder, fm)
    run_job(subfolder)
    check_running_jobs(origin_l)
    all_ave, ave_t, ave_a, ave_p, ave_e, ave_i = error_compute(subfolder,basic_method,flag_t=False,flag_a=True,flag_p=True,flag_i=True,flag_e=False,verbose = 2)
    details_record(x,all_ave, ave_t, ave_a, ave_p, ave_e, ave_i,trial_tag,fm)
    return all_ave

# f([1,0])

def start():
    space = search_space(fm)
    use_named_args(space)
    checkpoint_callback = skopt.callbacks.CheckpointSaver("/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}.pkl".format(trial_tag))


    res_gp = gp_minimize(func=f, dimensions=space, n_calls=55, n_random_starts = 15, random_state=50,x0 = initial_points,callback = [checkpoint_callback], kappa=1)

    printFinalResult(fm,res_gp)


def restart():
    
    space = search_space(fm)
    #==========================================================================!!!!!!
    #res = load("/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/Gn_all_pure_LDA_repara.pkl")
    res = load("/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}.pkl".format(trial_tag))  # should be the same with checkpoint file
    #==========================================================================!!!!!!
    x0 = res.x_iters
    y0 = res.func_vals
    # print(type(x0),type(y0))
    x0.append([1.01743487, 5.35548303e-04])
    y0 = np.append(y0,[4.649331],axis = 0)
    # print(type(x0),type(y0))
    #===============================================
    # x0_filter = []
    # y0_filter = []
    # for i,x in enumerate(x0):
    #     if x[0]>=0.9 and x[0]<=1.1:
    #         x0_filter.append([x[0]])
    #         y0_filter.append(y0[i])
    # x0 = x0_filter
    # y0 = y0_filter
    # print(x0)
    # filter the x0 and y0
    # only for Gn_all_pure_LDA_repara
    #================================
    checkpoint_callback = skopt.callbacks.CheckpointSaver("/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}.pkl".format(trial_tag))
    
    res_gp = gp_minimize(func=f, dimensions=space, n_random_starts = 2, x0 = x0, y0=y0, n_calls=30, callback = [checkpoint_callback], kappa=1)


    printFinalResult(fm,res_gp)

if __name__ == '__main__':
    # x = [1,0]
    # x = [0.9395924, 0.00861975] # LDA linear regression para E_x E_FA
    # x = [1.01743487, 5.35548303e-04] # GGA linear regression para E_x E_FA
    # x = [1.0223, 1.5660e-05] #GGA + FA,Bayesion
    # x = [1.02604396, 0.94771026] # BLYP linear regression: E_x E_c
    x = [ 1.02653263e+00,  9.43178060e-01, -1.25555222e-04]

    f(x)

    # error_compute(subfolder,basic_method,flag_t=False,flag_a=True,flag_p=True,flag_i=True,flag_e=False,verbose = 2)

    # start()

    #restart()

    # res = load("/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}.pkl".format(trial_tag))
    # x0 = np.array(res.x_iters)
    # y0 = np.array(res.func_vals)
    # n_var = len(x0[0])
    # print(y0)
    # if fm == 'XC_only':
    #     a0 = x0[:,0]
    #     a1 = x0[:,1]
    #     para_error = {'a0':a0,'a1':a1,'LOG_AVE_RMSE':np.log(y0)}
    #     df = pd.DataFrame(para_error,columns = ['a0','a1','LOG_AVE_RMSE'])
    # elif fm =='Exchange_only':
    #     a0 = x0[:,0]
    #     para_error = {'a0':a0,'LOG_AVE_RMSE':np.log(y0)}
    #     df = pd.DataFrame(para_error,columns = ['a0','LOG_AVE_RMSE'])
    # elif fm == 'Ex_FA':
    #     a0 = x0[:,0]
    #     a2 = x0[:,1]
    #     para_error = {'a0':a0,'a2':a2,'LOG_AVE_RMSE':np.log(y0)}
    #     df = pd.DataFrame(para_error,columns = ['a0','a2','LOG_AVE_RMSE'])
    # elif fm =='Correlation_only':
    #     a1 = x0[:,0]
    #     para_error = {'a1':a1,'LOG_AVE_RMSE':np.log(y0)}
    #     df = pd.DataFrame(para_error,columns = ['a1','LOG_AVE_RMSE'])
    # elif fm =='FA_only':
    #     a2 = x0[:,0]
    #     para_error = {'a2':a2,'LOG_AVE_RMSE':np.log(y0)}
    #     df = pd.DataFrame(para_error,columns = ['a2','LOG_AVE_RMSE'])
    # elif fm == 'All':
    #     a0 = x0[:,0]
    #     a1 = x0[:,1]
    #     a2 = x0[:,2]
    #     para_error = {'a0':a0,'a1':a1,'a2':a2,'LOG_AVE_RMSE':np.log(y0)}
    #     df = pd.DataFrame(para_error,columns = ['a0','a1','a2','LOG_AVE_RMSE'])

    # #remember to change the file name!!!!!!!!!!!!!!!!!!!!!!!
    # path = "/home/ysy1111/MyProj/FAhybrid/Gn/callbacks/{}.csv".format(trial_tag)
    # df.to_csv(path, index = False, header=True)
    # print(df)
    # print('eval points: %d' % len(x0))
    # printFinalResult(fm,res)

    # X,y = buildTotalMatrix(subfolder,flag_t=False,flag_a=True,flag_p=True,flag_i=True,flag_e=False,fitted_energy = ['E_x','E_c','E_FA'])
    # reg = LinearRegression(fit_intercept = False).fit(X, y)
    # print(reg.coef_)