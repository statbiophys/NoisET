import os
import shutil
import math
import time
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm.notebook import trange, tqdm

# Library functions to generate TCR repertoires


##===================================Initial-Distributions==========================================
def rho_counts_theo_minus_x(A, B, N_0):
    
    # I am disretizing the logspace with nfbins = 100000
    
    Cmin = 1
    freq_dtype = 'float32'

    N_cells = int(1e10)
    S_c = -(A+B/2)*N_cells/(N_0-1)
    
    alpha = -2*A/B
    
    nbins_1 = 100000
    
    logcountvec = np.linspace(np.log10(Cmin),np.log10(N_0), nbins_1)
    log_countvec_minus = np.array(np.log(np.power(10,logcountvec)) ,dtype=freq_dtype).flatten() 
    log_rho_minus = np.log(-(S_c/A))+ np.log(1-np.exp(-alpha*log_countvec_minus))
    
    N_clones_1 = -(S_c/A)*(np.log(N_0) - (1/alpha)*(1 - N_0**(-alpha)))
    
    
    return log_rho_minus, log_countvec_minus, N_clones_1

def rho_counts_theo_plus_x(A, B, N_0):
    
    # I am disretizing the logspace with nfbins = 100000, I can put a better discretization than for the minus
    # distribution
    
    Cmax = int(1e10)
    #Cmax = np.inf
    freq_dtype = 'float32'

    N_cells = int(1e10)
    S_c = -(A+B/2)*N_cells/(N_0 -1)
    
    alpha = -2*A/B
    
    nbins_2 = 100000
    
    logcountvec = np.linspace(np.log10(N_0),np.log10(Cmax), nbins_2 )
    log_countvec_plus = np.array(np.log(np.power(10,logcountvec)) ,dtype=freq_dtype).flatten() 
    log_rho_plus = np.log(N_0**alpha-1) + np.log(-(S_c/A)) -(alpha)*log_countvec_plus
    
    N_clones_2 = -(S_c/(A*alpha))*(1 - N_0**(-alpha))
    

    return log_rho_plus, log_countvec_plus, N_clones_2


def get_distsample(pmf,Nsamp, dtype='uint32'):
    '''
    generates Nsamp index samples of dtype (e.g. uint16 handles up to 65535 indices) from discrete probability mass function pmf.
    Handles multi-dimensional domain. N.B. Output is sorted.
    '''
    #assert np.sum(pmf)==1, "cmf not normalized!"
    
    shape = np.shape(pmf)
    sortindex = np.argsort(pmf, axis=None)#uses flattened array
    pmf = pmf.flatten()
    pmf = pmf[sortindex]
    cmf = np.cumsum(pmf)
   #print('cumulative distribution is equal to: ' + str(cmf[-1]))
    choice = np.random.uniform(high = cmf[-1], size = int(float(Nsamp)))
    index = np.searchsorted(cmf, choice)
    index = sortindex[index]
    index = np.unravel_index(index, shape)
    index = np.transpose(np.vstack(index))
    sampled_inds = np.array(index[np.argsort(index[:,0])], dtype=dtype)
    return sampled_inds

#======================================Propagator===============================================================


def gaussian_matrix(x_vec, x_i_vec_unique, A, B, t):
    
    
    x_vec_reshaped = np.reshape(x_vec, (len(x_vec), 1))
    ones_vec = np.ones((len(x_i_vec_unique), 1))
    M = np.multiply(ones_vec, x_vec_reshaped.T)
    x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))
    
    return (1/np.sqrt(2*np.pi*B*t))*np.exp((-1/(2*B*t))*(M - x_i_unique_reshaped - A*t)**2)

def gaussian_adsorption_matrix(x_vec, x_i_vec_unique, A, B, t):
    
    a = 0
    gauss = gaussian_matrix(x_vec, x_i_vec_unique, A, B, t)
    gauss_a = gaussian_matrix(x_vec, 2*a-x_i_vec_unique, A, B, t)
    x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))
    return gauss - np.exp((A*(a-x_i_unique_reshaped))/(B/2)) * gauss_a

def extinction_vector(x_i, A, B, t): 
    
    nbins = 2000
    eps = 1e-20
    #eps = 0
    x_vec = np.linspace(eps, np.max(x_i) - A*t + 3*np.sqrt(B*t), nbins)
    
    x_i_sorted = np.sort(x_i)
    
    xiind_vals, xi_start_ind, xi_counts=np.unique(x_i_sorted, return_counts=True,return_index=True)
    Prop_Matrix = gaussian_adsorption_matrix(x_vec, xiind_vals, A, B, t)
    
    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    p_ext = 1 - integ
    
    p_ext_new = np.zeros((len(x_i)))
    for it,xiind in enumerate(xiind_vals):
        p_ext_new[xi_start_ind[it]:xi_start_ind[it]+xi_counts[it]] = p_ext[it]
        
    test = np.random.uniform(0,1, size = (len(p_ext_new))) > p_ext_new
    results_extinction = test.astype(int)
    
    return results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext

#--------------------------------------------Source-term-no-frequency-dependency-------------------------------------------------------

## Functions for the source term
def gaussian_matrix_time(x_vec, x_i_scal, A, B, tvec_unique):
    
    x_vec_reshaped = np.reshape(x_vec, (len(x_vec), 1))
    ones_vec = np.ones((len(tvec_unique), 1))
    M = np.multiply(ones_vec, x_vec_reshaped.T)
    tvec_unique_reshaped = np.reshape(tvec_unique, (len(tvec_unique), 1))
    #x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))
    
    return (1/np.sqrt(2*np.pi*B*tvec_unique_reshaped))*np.exp((-1/(2*B*tvec_unique_reshaped))*(M - x_i_scal - A*tvec_unique_reshaped)**2)

def gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tvec_unique):
    
    a = 0
    gauss = gaussian_matrix_time(x_vec, x_i_scal, A, B, tvec_unique)
    gauss_a = gaussian_matrix_time(x_vec, 2*a-x_i_scal, A, B, tvec_unique)
    
    return gauss - np.exp((A*(a-x_i_scal))/(B/2)) * gauss_a

def Prop_Matrix_source( A, B, tvec): 
    
    #st = time.time()
    nbins = 2000
    N_0 = 40
    x_i_scal = np.log(N_0)
    t = np.max(tvec)
    x_vec = np.linspace(0, x_i_scal - A*t + 2*np.sqrt(B*t), nbins)
    
    tvec_sorted = np.sort(tvec)
    
    tiind_vals, ti_start_ind, ti_counts=np.unique(tvec_sorted, return_counts=True,return_index=True)
    Prop_Matrix = gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tiind_vals)
    
    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    
    
    #print(time.time() - st)
    return Prop_Matrix, x_vec, tiind_vals, ti_start_ind, ti_counts, integ

def extinction_vector_source(A, B, tvec): 
    
    #st = time.time()
    nbins = 2000
    N_0 = 40
    x_i_scal = np.log(N_0)
    t = np.max(tvec)
    x_vec = np.linspace(0, x_i_scal - A*t + 2*np.sqrt(B*t), nbins)
    
    tvec_sorted = np.sort(tvec)
    
    tiind_vals, ti_start_ind, ti_counts=np.unique(tvec_sorted, return_counts=True,return_index=True)
    Prop_Matrix = gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tiind_vals)
    
    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    p_ext = 1 - integ
    
    p_ext_new = np.zeros((len(tvec)))
    for it,tiind in enumerate(tiind_vals):
        p_ext_new[ti_start_ind[it]:ti_start_ind[it]+ti_counts[it]] = p_ext[it]
        
    test = np.random.uniform(0,1, size = (len(p_ext_new))) > p_ext_new
    results_extinction = test.astype(int)
    
    #print(time.time() - st)
    return results_extinction, Prop_Matrix, x_vec, tiind_vals, ti_start_ind, ti_counts, p_ext

#---------------------------------Linear-dependency-of-frequencies--------------------------------------

# Linear case A(f) = A_0 + A_1f 
# B(f) = B_0 + B_1f
def gaussian_matrix_linear_freq(x_vec, x_i_vec_unique, A_0, A_1, B_0, B_1, t):
    
    
    x_vec_reshaped = np.reshape(x_vec, (len(x_vec), 1))
    ones_vec = np.ones((len(x_i_vec_unique), 1))
    M = np.multiply(ones_vec, x_vec_reshaped.T)
    x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))
    
    B = B_0 + B_1*np.exp(x_i_unique_reshaped)
    A = A_0 + A_1*np.exp(x_i_unique_reshaped)
    
    return (1/np.sqrt(2*np.pi*B*t))*np.exp((-1/(2*B*t))*(M - x_i_unique_reshaped - A*t)**2)

def gaussian_adsorption_matrix_linear_freq(x_vec, x_i_vec_unique, A_0, A_1, B_0, B_1, t):
    
    a = 0
    gauss = gaussian_matrix_linear_freq(x_vec, x_i_vec_unique, A_0, A_1, B_0, B_1, t)
    gauss_a = gaussian_matrix_linear_freq(x_vec, 2*a-x_i_vec_unique, A_0, A_1, B_0, B_1, t)
    x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))
    A = A_0 + A_1*np.exp(x_i_unique_reshaped)
    B = B_0 + B_1*np.exp(x_i_unique_reshaped)
    return gauss - np.exp((A*(a-x_i_unique_reshaped))/(B/2)) * gauss_a

def extinction_vector_linear_freq(x_i, A_0, A_1, B_0, B_1, t): 
    
    nbins = 2000    
    x_i_sorted = np.sort(x_i)
    
    xiind_vals, xi_start_ind, xi_counts=np.unique(x_i_sorted, return_counts=True,return_index=True)
    xiind_vals_reshaped = np.reshape(xiind_vals, (len(xiind_vals), 1))
    A = A_0 + A_1*np.exp(xiind_vals_reshaped)
    B = B_0 + B_1*np.exp(xiind_vals_reshaped)
    
    x_vec = np.linspace(1e-20, np.max(x_i) - A_0 + 3*np.sqrt(B_0*t), nbins)
    print(x_vec[-1])
    
    Prop_Matrix = gaussian_adsorption_matrix_linear_freq(x_vec, xiind_vals, A_0, A_1, B_0, B_1, t)
    
    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    p_ext = 1 - integ
    
    p_ext_new = np.zeros((len(x_i)))
    for it,xiind in enumerate(xiind_vals):
        p_ext_new[xi_start_ind[it]:xi_start_ind[it]+xi_counts[it]] = p_ext[it]
        
    test = np.random.uniform(0,1, size = (len(p_ext_new))) > p_ext_new
    results_extinction = test.astype(int)
    

    return results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext

## source term for the frequencies dependency


def gaussian_matrix_time_linear_freq(x_vec, x_i_scal, A_0, A_1, B_0, B_1, tvec_unique):
    
    x_vec_reshaped = np.reshape(x_vec, (len(x_vec), 1))
    ones_vec = np.ones((len(tvec_unique), 1))
    M = np.multiply(ones_vec, x_vec_reshaped.T)
    tvec_unique_reshaped = np.reshape(tvec_unique, (len(tvec_unique), 1))
    #x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))

    A = A_0 + A_1*np.exp(x_i_scal)
    B = B_0 + B_1*np.exp(x_i_scal)

    
    return (1/np.sqrt(2*np.pi*B*tvec_unique_reshaped))*np.exp((-1/(2*B*tvec_unique_reshaped))*(M - x_i_scal - A*tvec_unique_reshaped)**2)

def gaussian_adsorption_matrix_time_linear_freq(x_vec, x_i_scal, A_0, A_1, B_0, B_1, tvec_unique):
    
    a = 0

    A = A_0 + A_1*np.exp(x_i_scal)
    B = B_0 + B_1*np.exp(x_i_scal)

    gauss = gaussian_matrix_time(x_vec, x_i_scal, A, B, tvec_unique)
    gauss_a = gaussian_matrix_time(x_vec, 2*a-x_i_scal, A, B, tvec_unique)
    
    return gauss - np.exp((A*(a-x_i_scal))/(B/2)) * gauss_a

def Prop_Matrix_source_linear_freq( A_0, A_1, B_0, B_1, tvec): 
    
    #st = time.time()
    nbins = 2000
    N_0 = 40
    x_i_scal = np.log(N_0)
    t = np.max(tvec)

    A = A_0 + A_1*np.exp(x_i_scal)
    B = B_0 + B_1*np.exp(x_i_scal)


    x_vec = np.linspace(0, x_i_scal - A*t + 2*np.sqrt(B*t), nbins)
    
    tvec_sorted = np.sort(tvec)
    
    tiind_vals, ti_start_ind, ti_counts=np.unique(tvec_sorted, return_counts=True,return_index=True)
    Prop_Matrix = gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tiind_vals)
    
    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    
    
    #print(time.time() - st)
    return Prop_Matrix, x_vec, tiind_vals, ti_start_ind, ti_counts, integ

def extinction_vector_source_linear_freq(A_0, A_1, B_0, B_1, tvec): 
    
    #st = time.time()
    nbins = 2000
    N_0 = 40
    x_i_scal = np.log(N_0)
    t = np.max(tvec)
    A = A_0 + A_1*np.exp(x_i_scal)
    B = B_0 + B_1*np.exp(x_i_scal)

    x_vec = np.linspace(0, x_i_scal - A*t + 2*np.sqrt(B*t), nbins)
    
    tvec_sorted = np.sort(tvec)
    
    tiind_vals, ti_start_ind, ti_counts=np.unique(tvec_sorted, return_counts=True,return_index=True)
    Prop_Matrix = gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tiind_vals)
    
    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    p_ext = 1 - integ
    
    p_ext_new = np.zeros((len(tvec)))
    for it,tiind in enumerate(tiind_vals):
        p_ext_new[ti_start_ind[it]:ti_start_ind[it]+ti_counts[it]] = p_ext[it]
        
    test = np.random.uniform(0,1, size = (len(p_ext_new))) > p_ext_new
    results_extinction = test.astype(int)
    
    #print(time.time() - st)
    return results_extinction, Prop_Matrix, x_vec, tiind_vals, ti_start_ind, ti_counts, p_ext



#-----------------------------------------------------------------------------------------------------------------

def generator_diffusion_linear_freq(i, A_0, A_1, B_0, B_1, N_0, t):
    
    '''
    Generate steady state distributions of log-populations/ populations for one individual with following parameters:
    A_0 :
    B_0 :
    A_1 : 
    B_1 : 
    N_0 :
    
    Second-step of this code is to use Euler-Mayurana scheme to simulate stochastic simulations from initial 
    time to 2 years afterwards, it is important to have small enough dt for the absorbing 
    boundary condition to be satisfied. 

    '''
    
    eps = 1e-20
    
    ## Choose initial size of the immune system to be 1e10 (for a mouse)
    N_cells = int(1e10)
    
    #Parameters for the repertoire generation
    alpha_rho = -1 + 2*A_0/B_0
    N_ext = 1
    freq_dtype = 'float32' 

    print(np.random.random())
    
    #==========================generate the steady state distribution===============================
    
    #for counts < N0: 
    logrhofvec,logfvec, N_clones_1 = rho_counts_theo_minus_x(A_0, B_0, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_1,dtype='uint32').flatten()
    print("generation population smaller than N_0: check")
    
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_minus = np.sum(counts_generated)
    print(str(C_f_minus) + ' N_cells smaller than N_0')
    log_cminus_generated = logcvec_generated
    logrhofvec_1,logfvec_1 = logrhofvec,logfvec
    print(str(N_clones_1) + ' N_clones_1')
    
    #for counts > N0:
    
    logrhofvec,logfvec, N_clones_2 = rho_counts_theo_plus_x(A_0, B_0, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_2,dtype='uint32').flatten()
    print("generation population larger than N_0: check")
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_plus = np.sum(counts_generated)
    print(str(C_f_plus) + ' N_cells larger than N_0')
    log_cplus_generated = logcvec_generated
    logrhofvec_2,logfvec_2 = logrhofvec,logfvec
    print(str(N_clones_2) + ' N_clones_2')
    
    #===================================================
    
    N_clones = int(N_clones_1 + N_clones_2)
    print('N_clones= ' + str(N_clones))

    S_c = - (A_0 + B_0/2)*(N_cells/(N_0-1))
    print('N_clones_theory= ' + str(-(S_c/A_0)*np.log(N_0)))
    
    
    x_i = np.concatenate((log_cminus_generated, log_cplus_generated), axis = None)
    
    N_total_cells_generated = np.sum(np.exp(x_i))
    print("N_total_cells_generated/N_total_cells:" + str(N_total_cells_generated/N_cells))
    
    
    
    results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext = extinction_vector_linear_freq(x_i, A_0, A_1, B_0, B_1, t)
    #x_vec = np.linspace(0, 30*B*t, 2000)
    dx=np.asarray(np.diff(x_vec)/2., dtype='float32')
    
    x_i_noext= x_i[np.where(results_extinction ==1)]
    
    #print('NO_EXTINCT'+ str(len(x_i_noext)))
    
    
    #xiind_vals, xi_start_ind, xi_counts=np.unique(x_i,return_counts=True,return_index=True)
    #Prop_Matrix = gaussian_adsorption_matrix(x_vec, xiind_vals, A, B, t)
    
    
    x_f = np.zeros((len(x_i)))
    
    for i in range(len(xiind_vals)): 
        
        #Prop_adsorp = Prop_Matrix[i,:]
        #print(np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1]))
        
        if (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1])) < 1e-7:
            pass
        
        else:
        
            Prop_adsorp = Prop_Matrix[i,:] / (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1]))


            #print(np.dot(dx, Prop_adsorp[1:] + Prop_adsorp[:-1]))
            integ = Prop_adsorp[np.newaxis,:]
            #print(np.dot(dx[np.newaxis,:], integ[i,1:] + integ[i, :-1]))
            #integ_bis = integ/np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])))
            f_samples_inds =get_distsample(np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), xi_counts[i],dtype='uint32').flatten()

            x_f[xi_start_ind[i]:xi_start_ind[i]+xi_counts[i]] = x_vec[f_samples_inds]
    
    x_f = np.multiply(x_f,results_extinction)
    
    
    #x_f[x_f == 0] = -eps
    x_f[x_f == 0] = -np.inf
    
    N_extinction = np.sum(1- results_extinction)
    N_extinction_bis = len(x_f[x_f == -np.inf])
    
    print('Number of extinction= ' + str(N_extinction))
    print('Number of extinction bis= ' + str(N_extinction_bis))
    sim_ext = (N_extinction/len(results_extinction))*100
    theo_ext = (-A_0/np.log(N_0))*100
    print('simulations % of extinction= ' + str((N_extinction/len(results_extinction))*100/t) + '%')
    print('theoretical % of extinction= ' + str((-A_0/np.log(N_0))*100) + '%')

    N_source = S_c*t

    print('Number of insertions= ' +str(N_source))

    N_source = int(N_source)

    eps = 1e-8

    nbins_time = 10000
    time_vec_span = np.linspace(eps, t, nbins_time)
    time_vec = np.random.choice(time_vec_span, N_source)
    time_vec = np.sort(time_vec)
    
    results_extinction_source, Prop_Matrix_source, x_vec_source, tiind_vals, ti_start_ind, ti_counts, p_ext_source = extinction_vector_source_linear_freq(A_0, A_1, B_0, B_1, time_vec)

    #if np.sum(1-results_extinction_source)/len(x_i_LB) < 1e-2:
    #    pass

    dx_source=np.asarray(np.diff(x_vec_source)/2., dtype='float32')



    x_source_LB = np.zeros((N_source))
    for i in range(len(tiind_vals)): 
        
        if (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1])) < 1e-7:
            pass
        
        else:
            Prop_adsorp_s = Prop_Matrix_source[i,:]
            Prop_adsorp_s = Prop_Matrix_source[i,:] / (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1]))


            integ = Prop_adsorp_s[np.newaxis,:]
            f_samples_inds_s = get_distsample(np.asarray((dx_source[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), ti_counts[i],dtype='uint32').flatten()

            x_source_LB[ti_start_ind[i]:ti_start_ind[i]+ti_counts[i]] = x_vec_source[f_samples_inds_s]
            
        
    x_source_LB = np.multiply(x_source_LB, results_extinction_source)
    
    x_source_LB[x_source_LB == 0] = -np.inf


            

    
    return x_i, x_f, Prop_Matrix, p_ext, results_extinction, time_vec, results_extinction_source, x_source_LB

#---------------------------------------------------------------------------------------------------------------



def generator_diffusion_LB(i, A, B, N_0, t):
    
    '''
    Generate steady state distributions of log-populations/ populations for one individual with following parameters:
    A :
    B : 
    N_0 :
    
    Second-step of this code is to use Euler-Mayurana scheme to simulate stochastic simulations from initial 
    time to 2 years afterwards, it is important to have small enough dt for the absorbing 
    boundary condition to be satisfied. 

    '''
    
    eps = 1e-20
    
    ## Choose initial size of the immune system to be 1e10 (for a mouse)
    N_cells = int(1e10)
    
    #Parameters for the repertoire generation
    alpha_rho = -1 + 2*A/B
    N_ext = 1
    freq_dtype = 'float32' 
    
    #==========================generate the steady state distribution===============================
    
    #for counts < N0: 
    logrhofvec,logfvec, N_clones_1 = rho_counts_theo_minus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_1,dtype='uint32').flatten()
    #print("generation population smaller than N_0: check")
    
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_minus = np.sum(counts_generated)
    print(str(C_f_minus) + ' cells smaller than N_0')
    log_cminus_generated = logcvec_generated
    logrhofvec_1,logfvec_1 = logrhofvec,logfvec
    print(str(N_clones_1) + ' N_clones_1')
    
    #for counts > N0:
    
    logrhofvec,logfvec, N_clones_2 = rho_counts_theo_plus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_2,dtype='uint32').flatten()
    #print("generation population larger than N_0: check")
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_plus = np.sum(counts_generated)
    print(str(C_f_plus) + ' N_cells larger than N_0')
    log_cplus_generated = logcvec_generated
    logrhofvec_2,logfvec_2 = logrhofvec,logfvec
    print(str(N_clones_2) + ' N_clones_2')
    
    #===================================================
    
    N_clones = int(N_clones_1 + N_clones_2)
    print('N_clones= ' + str(N_clones))

    S_c = - (A + B/2)*(N_cells/(N_0-1))
    print('N_clones_theory= ' + str(-(S_c/A)*np.log(N_0)))
    
    
    x_i = np.concatenate((log_cminus_generated, log_cplus_generated), axis = None)
    
    N_total_cells_generated = np.sum(np.exp(x_i))
    print("N_total_cells_generated/N_total_cells:" + str(N_total_cells_generated/N_cells))
    
    
    
    results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext = extinction_vector(x_i, A, B, t)
    #x_vec = np.linspace(0, 30*B*t, 2000)
    dx=np.asarray(np.diff(x_vec)/2., dtype='float32')
    
    x_i_noext= x_i[np.where(results_extinction ==1)]
    
    #print('NO_EXTINCT'+ str(len(x_i_noext)))
    
    
    #xiind_vals, xi_start_ind, xi_counts=np.unique(x_i,return_counts=True,return_index=True)
    #Prop_Matrix = gaussian_adsorption_matrix(x_vec, xiind_vals, A, B, t)
    
    
    x_f = np.zeros((len(x_i)))
    
    for i in range(len(xiind_vals)): 
        
        
        if (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1])) < 1e-7:
            pass
        
        else:
        
            Prop_adsorp = Prop_Matrix[i,:] / (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1]))


            #print(np.dot(dx, Prop_adsorp[1:] + Prop_adsorp[:-1]))
            integ = Prop_adsorp[np.newaxis,:]
            #print(np.dot(dx[np.newaxis,:], integ[i,1:] + integ[i, :-1]))
            #integ_bis = integ/np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])))
            f_samples_inds =get_distsample(np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), xi_counts[i],dtype='uint32').flatten()

            x_f[xi_start_ind[i]:xi_start_ind[i]+xi_counts[i]] = x_vec[f_samples_inds]
    
    x_f = np.multiply(x_f,results_extinction)
    
    
    #x_f[x_f == 0] = -eps
    x_f[x_f == 0] = -np.inf
    
    N_extinction = np.sum(1- results_extinction)
    N_extinction = len(x_f[x_f == -np.inf])
    
    print('Number of extinction= ' + str(N_extinction))
    sim_ext = (N_extinction/len(results_extinction))*100
    theo_ext = (-A/np.log(N_0))*100
    print('simulations % of extinction= ' + str((N_extinction/len(results_extinction))*100/t) + '%')
    print('theoretical % of extinction= ' + str((-A/np.log(N_0))*100) + '%')
        

    #Source term

    N_source = S_c*t

    print('Number of insertions= ' +str(N_source))

    N_source = int(N_source)

    eps = 1e-8
    time_vec_span = np.linspace(eps, t, 5000)
    time_vec = np.random.choice(time_vec_span, N_source)
    time_vec = np.sort(time_vec)
    
    results_extinction_source, Prop_Matrix_source, x_vec_source, tiind_vals, ti_start_ind, ti_counts, p_ext_source = extinction_vector_source(A, B, time_vec)

    #if np.sum(1-results_extinction_source)/len(x_i_LB) < 1e-2:
    #    pass

    dx_source=np.asarray(np.diff(x_vec_source)/2., dtype='float32')



    x_source_LB = np.zeros((N_source))
    for i in range(len(tiind_vals)): 
        
        if (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1])) < 1e-7:
            pass
        
        else:
            Prop_adsorp_s = Prop_Matrix_source[i,:]
            Prop_adsorp_s = Prop_Matrix_source[i,:] / (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1]))


            integ = Prop_adsorp_s[np.newaxis,:]
            f_samples_inds_s = get_distsample(np.asarray((dx_source[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), ti_counts[i],dtype='uint32').flatten()

            x_source_LB[ti_start_ind[i]:ti_start_ind[i]+ti_counts[i]] = x_vec_source[f_samples_inds_s]
            
        
    x_source_LB = np.multiply(x_source_LB, results_extinction_source)
    
    x_source_LB[x_source_LB == 0] = -np.inf


    
    return x_i, x_f, Prop_Matrix, p_ext, results_extinction, time_vec, results_extinction_source, x_source_LB


def generator_diffusion_EM(B, A, N_0, t, nt):
    
    '''
    Generate steady state distributions of log-populations/ populations for one individual with following parameters:
    A :
    B : 
    N_0 :
    
    Second-step of this code is to use Euler-Mayurana scheme to simulate stochastic simulations from initial 
    time to 2 years afterwards, it is important to have small enough dt for the absorbing 
    boundary condition to be satisfied. 

    '''
    
    ## Choose initial size of the immune system to be 1e10 (for a mouse)
    N_cells = int(1e10)
    
    #Parameters for the repertoire generation
    alpha_rho = -1 + 2*A/B
    N_ext = 1
    tvec = np.linspace(0, t, nt)
    dt = np.diff(tvec)[0] # time discretization : every week
    freq_dtype = 'float32' 
    
    #==========================generate the steady state distribution===============================
    
    #for counts < N0: 
    logrhofvec,logfvec, N_clones_1 = rho_counts_theo_minus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_1,dtype='uint32').flatten()
    print("generation population smaller than N_0: check")
    
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_minus = np.sum(counts_generated)
    print(str(C_f_minus) + ' N_cells smaller than N_0')
    log_cminus_generated = logcvec_generated
    logrhofvec_1,logfvec_1 = logrhofvec,logfvec
    print(str(N_clones_1) + ' N_clones_1')
    
    #for counts > N0:
    
    logrhofvec,logfvec, N_clones_2 = rho_counts_theo_plus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_2,dtype='uint32').flatten()
    print("generation population larger than N_0: check")
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_plus = np.sum(counts_generated)
    print(str(C_f_plus) + ' N_cells larger than N_0')
    log_cplus_generated = logcvec_generated
    logrhofvec_2,logfvec_2 = logrhofvec,logfvec
    print(str(N_clones_2) + ' N_clones_2')
    
    #===================================================
    
    N_clones = int(N_clones_1 + N_clones_2)
    print('N_clones=' + str(N_clones))

    S_c = - (A + B/2)*(N_cells/(N_0-1))
    print('N_clones_theory=' + str(-(S_c/A)*np.log(N_0)))
    
    
    N_c = int(S_c*dt)
    N_total_new = N_c*len(np.diff(tvec))
    N_total = int(N_clones + 1.5*N_total_new)
    
    print('number of new-clones added to the repertoire in dt ' + str(N_c))
    print('number of new-clones added to the repertoire in ' + str(t) + ' years ' + str(N_total_new))
          
    New_clones_logfreq = -np.inf*np.ones((int(1.5*N_total_new))) 
    log_nvec_generated = np.concatenate((log_cminus_generated, log_cplus_generated, New_clones_logfreq), axis = None)
    
    print("C_min = " + str(np.exp(np.min(log_cminus_generated))))

    N_total_cells_generated = np.sum(np.exp(log_nvec_generated))
    print("N_total_cells_generated/N_total_cells:" + str(N_total_cells_generated/N_cells))
    
    
    #================EM-scheme===========================
    
    x = log_nvec_generated
    print('shape of x is ' + str(np.shape(x)))
    N_clones_not_extinct_total = [N_clones]
    N_clones_no_source = [N_clones]
    
    N_cells_traj = [N_total_cells_generated]
    N_cells_no_source = [N_total_cells_generated]
    
    
    N_clones_bis_total = N_clones
    
    for i in range(len(tvec)-1):
        
        
        # 1/ Diffusion of the existing clones from t -> t + dt
        x[:N_clones_bis_total ] = x[:N_clones_bis_total] + A*(dt) + np.sqrt(B*dt)*np.random.randn(N_clones_bis_total) 
        
        # 2/ Extinction of all clones at time t + dt that have a log-size < log(Next)
        x[x<= np.log(N_ext)] = -np.inf
        #N_clones_extinct.append(len(x[x == -np.inf]))
        
        # 3/ Introduction of N_c_poisson clones at time t + dt of log-size x_int = log(Nint)
        N_c_poisson = np.random.poisson(N_c, 1)
        x[N_clones_bis_total:N_clones_bis_total + int(N_c_poisson)] = np.log(N_0)
        
        # 4/ Keep track of the total number of clones and the total number of cells
        N_clones_not_extinct_total.append(len(x[x != -np.inf]))
        N_cells_traj.append(np.sum(np.exp(x)))

        
        # 5/ Update of the total number of clones that have existed since intial time t_0
        N_clones_bis_total = N_clones_bis_total + int(N_c_poisson)
          
        x_bis = x[:N_clones]
        N_cells_no_source.append(np.sum(np.exp(x_bis)))
        N_clones_no_source.append(len(x_bis[x_bis != -np.inf]))
            
          
    N_clones_not_extinct_total = np.reshape(np.array(N_clones_not_extinct_total), (len(tvec)))
    N_cells_traj = np.reshape(np.array(N_cells_traj), (len(tvec)))
          
    N_clones_no_source = np.reshape(np.array(N_clones_no_source), (len(tvec)))
    N_cells_no_source= np.reshape(np.array(N_cells_no_source), (len(tvec)))
    
    print('% extinction EM' + str(((N_clones_no_source[0]-N_clones_no_source[-1])/N_clones_no_source[0])*100))
    print('theoretical % of extinction EM= ' + str((-A/np.log(N_0))*100) + '%')
    
    sim_ext = ((N_clones_no_source[0]-N_clones_no_source[-1])/N_clones_no_source[0])*100
    theo_ext = (-A/np.log(N_0))*100
          
    log_nvec_generated = np.concatenate((log_cminus_generated, log_cplus_generated, New_clones_logfreq), axis = None)
          
    return log_nvec_generated, x, N_cells_traj, N_clones_not_extinct_total, N_cells_no_source, N_clones_no_source


def generator_diffusion_EM_LB(B, A, N_0, t, nt):
    
    '''
    Generate steady state distributions of log-populations/ populations for one individual with following parameters:
    A :
    B : 
    N_0 :
    
    Second-step of this code is to use Euler-Mayurana scheme to simulate stochastic simulations from initial 
    time to 2 years afterwards, it is important to have small enough dt for the absorbing 
    boundary condition to be satisfied. 

    '''
    
    ## Choose initial size of the immune system to be 1e10 (for a mouse)
    N_cells = int(1e10)

    #print(i)
    
    #Parameters for the repertoire generation
    alpha_rho = -1 + 2*A/B
    N_ext = 1
    tvec = np.linspace(0, t, nt)
    dt = np.diff(tvec)[0] # time discretization : every week
    freq_dtype = 'float32' 
    
    #==========================generate the steady state distribution===============================
    
    #for counts < N0: 
    logrhofvec,logfvec, N_clones_1 = rho_counts_theo_minus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_1,dtype='uint32').flatten()
    #print("generation population smaller than N_0: check")
    
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_minus = np.sum(counts_generated)
    print(str(C_f_minus) + ' N_cells smaller than N_0')
    log_cminus_generated = logcvec_generated
    logrhofvec_1,logfvec_1 = logrhofvec,logfvec
    print(str(N_clones_1) + ' N_clones_1')
    
    #for counts > N0:
    
    logrhofvec,logfvec, N_clones_2 = rho_counts_theo_plus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_2,dtype='uint32').flatten()
    #print("generation population larger than N_0: check")
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_plus = np.sum(counts_generated)
    print(str(C_f_plus) + ' N_cells larger than N_0')
    log_cplus_generated = logcvec_generated
    logrhofvec_2,logfvec_2 = logrhofvec,logfvec
    print(str(N_clones_2) + ' N_clones_2')
    
    #===================================================
    
    N_clones = int(N_clones_1 + N_clones_2)
    print('N_clones=' + str(N_clones))

    S_c = - (A + B/2)*(N_cells/(N_0-1))
    print('N_clones_theory=' + str(-(S_c/A)*np.log(N_0)))
    
    
    N_c = int(S_c*dt)
    N_total_new = N_c*len(np.diff(tvec))
    N_total = int(N_clones + 1.5*N_total_new)
    
    #print('number of new-clones added to the repertoire in dt ' + str(N_c))
    #print('number of new-clones added to the repertoire in ' + str(t) + ' years ' + str(N_total_new))
          
    New_clones_logfreq = -np.inf*np.ones((int(1.5*N_total_new))) 
    log_nvec_generated = np.concatenate((log_cminus_generated, log_cplus_generated, New_clones_logfreq), axis = None)
    
    #print("C_min = " + str(np.exp(np.min(log_cminus_generated))))

    N_total_cells_generated = np.sum(np.exp(log_nvec_generated))
    print("EM-N_total_cells_generated/N_total_cells:" + str(N_total_cells_generated/N_cells))

    #================LB-scheme===========================

    x_i_LB = np.concatenate((log_cminus_generated, log_cplus_generated), axis = None)
    
    N_total_cells_generated = np.sum(np.exp(x_i_LB))
    print("LB-N_total_cells_generated/N_total_cells:" + str(N_total_cells_generated/N_cells))
    
    
    
    results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext = extinction_vector(x_i_LB, A, B, t)
    #x_vec = np.linspace(0, 30*B*t, 2000)
    dx=np.asarray(np.diff(x_vec)/2., dtype='float32')
    
    x_i_noext= x_i_LB[np.where(results_extinction ==1)]
    
    x_f_LB = np.zeros((len(x_i_LB)))
    
    for i in range(len(xiind_vals)): 
        
        
        if (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1])) < 1e-7:
            pass
        
        else:
        
            Prop_adsorp = Prop_Matrix[i,:] / (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1]))


            #print(np.dot(dx, Prop_adsorp[1:] + Prop_adsorp[:-1]))
            integ = Prop_adsorp[np.newaxis,:]
            #print(np.dot(dx[np.newaxis,:], integ[i,1:] + integ[i, :-1]))
            #integ_bis = integ/np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])))
            f_samples_inds =get_distsample(np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), xi_counts[i],dtype='uint32').flatten()

            x_f_LB[xi_start_ind[i]:xi_start_ind[i]+xi_counts[i]] = x_vec[f_samples_inds]
    
    x_f_LB = np.multiply(x_f_LB,results_extinction)
    
    
    #x_f[x_f == 0] = -eps
    x_f_LB[x_f_LB == 0] = -np.inf
    
    #N_extinction = np.sum(1- results_extinction)
    N_extinction = len(x_f_LB[x_f_LB == -np.inf])
    
    print('Number of extinctions= ' + str(N_extinction))

    N_source = S_c*t

    print('Number of insertions= ' +str(N_source))

    N_source = int(N_source)

    eps = 1e-8

    time_vec_span = np.linspace(eps, t, 5000)
    time_vec = np.random.choice(time_vec_span, N_source)
    time_vec = np.sort(time_vec)
    
    results_extinction_source, Prop_Matrix_source, x_vec_source, tiind_vals, ti_start_ind, ti_counts, p_ext_source = extinction_vector_source(A, B, time_vec)

    #if np.sum(1-results_extinction_source)/len(x_i_LB) < 1e-2:
    #    pass

    dx_source=np.asarray(np.diff(x_vec_source)/2., dtype='float32')



    x_source_LB = np.zeros((N_source))
    for i in range(len(tiind_vals)): 
        
        if (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1])) < 1e-7:
            pass
        
        else:
            Prop_adsorp_s = Prop_Matrix_source[i,:]
            Prop_adsorp_s = Prop_Matrix_source[i,:] / (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1]))


            integ = Prop_adsorp_s[np.newaxis,:]
            f_samples_inds_s = get_distsample(np.asarray((dx_source[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), ti_counts[i],dtype='uint32').flatten()

            x_source_LB[ti_start_ind[i]:ti_start_ind[i]+ti_counts[i]] = x_vec_source[f_samples_inds_s]
            
        
    x_source_LB = np.multiply(x_source_LB, results_extinction_source)
    
    x_source_LB[x_source_LB == 0] = -np.inf




    sim_ext = (N_extinction/len(results_extinction))*100
    theo_ext = (-A/np.log(N_0))*100
    print('simulations % of extinction LB= ' + str((N_extinction/len(results_extinction))*100/t) + '%')
    print('theoretical % of extinction LB= ' + str((-A/np.log(N_0))*100) + '%')
        
    
    
    #================EM-scheme===========================
    
    x = log_nvec_generated
    N_clones_not_extinct_total = [N_clones]
    N_clones_no_source = [N_clones]
    
    N_cells_traj = [N_total_cells_generated]
    N_cells_no_source = [N_total_cells_generated]
    
    
    N_clones_bis_total = N_clones
    
    for i in range(len(tvec)-1):
        
        
        # 1/ Diffusion of the existing clones from t -> t + dt
        x[:N_clones_bis_total ] = x[:N_clones_bis_total] + A*(dt) + np.sqrt(B*dt)*np.random.randn(N_clones_bis_total) 
        
        # 2/ Extinction of all clones at time t + dt that have a log-size < log(Next)
        x[x<= np.log(N_ext)] = -np.inf
        #N_clones_extinct.append(len(x[x == -np.inf]))
        
        # 3/ Introduction of N_c_poisson clones at time t + dt of log-size x_int = log(Nint)
        N_c_poisson = np.random.poisson(N_c, 1)
        x[N_clones_bis_total:N_clones_bis_total + int(N_c_poisson)] = np.log(N_0)
        
        # 4/ Keep track of the total number of clones and the total number of cells
        N_clones_not_extinct_total.append(len(x[x != -np.inf]))
        N_cells_traj.append(np.sum(np.exp(x)))

        
        # 5/ Update of the total number of clones that have existed since intial time t_0
        N_clones_bis_total = N_clones_bis_total + int(N_c_poisson)
          
        x_bis = x[:N_clones]
        N_cells_no_source.append(np.sum(np.exp(x_bis)))
        N_clones_no_source.append(len(x_bis[x_bis != -np.inf]))
            
          
    N_clones_not_extinct_total = np.reshape(np.array(N_clones_not_extinct_total), (len(tvec)))
    N_cells_traj = np.reshape(np.array(N_cells_traj), (len(tvec)))
          
    N_clones_no_source = np.reshape(np.array(N_clones_no_source), (len(tvec)))
    N_cells_no_source= np.reshape(np.array(N_cells_no_source), (len(tvec)))
    
    print('simulations % of extinction EM= ' + str(((N_clones_no_source[0]-N_clones_no_source[-1])/N_clones_no_source[0])*100))
    print('theoretical % of extinction EM= ' + str((-A/np.log(N_0))*100*t) + '%')
    
    sim_ext = ((N_clones_no_source[0]-N_clones_no_source[-1])/N_clones_no_source[0])*100
    theo_ext = (-A/np.log(N_0))*100
          
    log_nvec_generated = np.concatenate((log_cminus_generated, log_cplus_generated, New_clones_logfreq), axis = None)
          
    return log_nvec_generated, x, N_cells_traj, N_clones_not_extinct_total, N_cells_no_source, N_clones_no_source, x_i_LB, x_f_LB, x_source_LB #N_cells_final_LB


def generator_diffusion_trajectories(B, A, N_0, t):
    
    
    
    '''
    Generate steady state distributions of log-populations/ populations for one individual with following parameters:
    A :
    B : 
    N_0 :
    t : t is a vector here !!!
    
    Second-step of this code is to use Euler-Mayurana scheme to simulate stochastic simulations from initial 
    time to 2 years afterwards, it is important to have small enough dt for the absorbing 
    boundary condition to be satisfied. 

    '''
    
    ## Choose initial size of the immune system to be 1e10 (for a mouse)
    N_cells = int(1e10)
    
    #Parameters for the repertoire generation
    alpha_rho = -1 + 2*A/B
    N_ext = 1
    freq_dtype = 'float32' 
    eps = 1e-8
    
    #==========================generate the steady state distribution===============================
    
    #for counts < N0: 
    logrhofvec,logfvec, N_clones_1 = rho_counts_theo_minus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_1,dtype='uint32').flatten()
    print("generation population smaller than N_0: check")
    
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_minus = np.sum(counts_generated)
    print(str(C_f_minus) + ' N_cells smaller than N_0')
    log_cminus_generated = logcvec_generated
    logrhofvec_1,logfvec_1 = logrhofvec,logfvec
    print(str(N_clones_1) + ' N_clones_1')
    
    #for counts > N0:
    
    logrhofvec,logfvec, N_clones_2 = rho_counts_theo_plus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_2,dtype='uint32').flatten()
    print("generation population larger than N_0: check")
    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_plus = np.sum(counts_generated)
    print(str(C_f_plus) + ' N_cells larger than N_0')
    log_cplus_generated = logcvec_generated
    logrhofvec_2,logfvec_2 = logrhofvec,logfvec
    print(str(N_clones_2) + ' N_clones_2')
    
    #===================================================
    
    N_clones = int(N_clones_1 + N_clones_2)
    print('N_clones=' + str(N_clones))

    S_c = - (A + B/2)*(N_cells/(N_0-1))
    print('N_clones_theory=' + str(-(S_c/A)*np.log(N_0)))
    
    
    x_i = np.concatenate((log_cminus_generated, log_cplus_generated), axis = None)
    N_total_cells_generated = np.sum(np.exp(x_i))
    print("LB-N_total_cells_generated/N_total_cells:" + str(N_total_cells_generated/N_cells))
    
    initial_log_f = x_i
    
    st = time.time()
    delta_t_vec = np.diff(t)
    x_f = np.zeros((len(x_i), len(delta_t_vec)))
    N_cells_total = [np.sum(np.exp(x_i))]
    
    
    
    for k in tqdm(range(len(delta_t_vec))):
        
        
       
        results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext = extinction_vector(initial_log_f, A, B, delta_t_vec[k])
        print('TIME BETWEEN 2 GENERATIONS ' + str(delta_t_vec[k]) + 'YEAR(S)')
        
        print('TOTAL EXTINCTION IN ' + str(delta_t_vec[k]) + 'YEAR(S) :' + str(np.sum(1- results_extinction)))
    #x_vec = np.linspace(0, 30*B*t, 2000)
        dx = np.asarray(np.diff(x_vec)/2., dtype='float32')

        print(len(xiind_vals))
    
        #x_i_noext= x_i[np.where(results_extinction ==1)]
    
    #print('NO_EXTINCT'+ str(len(x_i_noext)))
    
    
    #xiind_vals, xi_start_ind, xi_counts=np.unique(x_i,return_counts=True,return_index=True)
    #Prop_Matrix = gaussian_adsorption_matrix(x_vec, xiind_vals, A, B, t)
    
    
    
        for i in tqdm(range(len(xiind_vals))): 
        
            Prop_adsorp = Prop_Matrix[i,:]
            integ= Prop_adsorp[np.newaxis,:]
            f_samples_inds =get_distsample(np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), xi_counts[i],dtype='uint32').flatten()

            x_f[xi_start_ind[i]:xi_start_ind[i]+xi_counts[i], k] = x_vec[f_samples_inds]
            
        
            
    
        x_f_vec1D = np.multiply(x_f[:,k], results_extinction)
        x_f_vec1D[x_f_vec1D == 0] = -np.inf

        N_extinction_bis = len(x_f_vec1D[x_f_vec1D == -np.inf])
    
        
        print('Number of extinction bis= ' + str(N_extinction_bis))
        x_f[:,k] = x_f_vec1D
        
        x_f_vec1D = np.sort(x_f_vec1D.flatten())
        
        
        initial_log_f = x_f_vec1D


        N_source = int(S_c*delta_t_vec[k])
        print('Number of insertions= ' +str(N_source))
        time_vec_span = np.linspace(eps, delta_t_vec[k], 5000)
        time_vec = np.random.choice(time_vec_span, N_source)
        time_vec = np.sort(time_vec)
    
        results_extinction_source, Prop_Matrix_source, x_vec_source, tiind_vals, ti_start_ind, ti_counts, p_ext_source = extinction_vector_source(A, B, time_vec)
        dx_source=np.asarray(np.diff(x_vec_source)/2., dtype='float32')
        x_source_LB = np.zeros((N_source))
        for i in range(len(tiind_vals)): 
        
            if (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1])) < 1e-7:
                pass
        
            else:
                Prop_adsorp_s = Prop_Matrix_source[i,:]
                Prop_adsorp_s = Prop_Matrix_source[i,:] / (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1]))


                integ = Prop_adsorp_s[np.newaxis,:]
                f_samples_inds_s = get_distsample(np.asarray((dx_source[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), ti_counts[i],dtype='uint32').flatten()

                x_source_LB[ti_start_ind[i]:ti_start_ind[i]+ti_counts[i]] = x_vec_source[f_samples_inds_s]
            
        
        x_source_LB = np.multiply(x_source_LB, results_extinction_source)
        x_source_LB[x_source_LB == 0] = -np.inf

        N_cells_k = np.sum(np.exp(initial_log_f)) + np.sum(np.exp(x_source_LB))
        N_cells_total.append(N_cells_k)
        
    print("--- %s seconds ---" % (time.time() - st))
    
    
    
    return x_i, x_f, Prop_Matrix, p_ext, results_extinction, N_cells_total

def experimental_sampling_diffusion_Poisson(NreadsI, NreadsII, x_0, x_2, t, N_cell_0, N_cell_2):
    
    
    #----------------------------Counts generation --------------------------------------------
    
    ##Initial condition
    N_total_0 = len(x_0[x_0 != -np.inf])
    x_0_bis = x_0[x_0 != -np.inf]
    
    print('Number of clones at initial time ' + str(N_total_0))
    
    N_total_2 = len(x_2[x_2 != -np.inf])
    x_2_bis = x_2[x_2 != -np.inf]
    
    
    print('Number of clones after ' + str(t) + ' year(s) ' +  str(N_total_2))
    
    #N_total = min(N_total_0, N_total_2)
    assert len(x_0) == len(x_2)
    N_total = len(x_0)
    
    x_2_final = x_2[:N_total]
    

    f_vec_initial = np.exp(x_0)/N_cell_0
    m=float(NreadsI)*f_vec_initial
    n_counts_day_0 = np.random.poisson(m, size =(1, int(N_total)))
    n_counts_day_0 = n_counts_day_0[0,:]
    
    #print('done')
    
    #Final condition
    f_vec_end = np.exp(x_2_final)/N_cell_2
    m=float(NreadsII)*f_vec_end
    #print(m)
    print('MEAN N : ' + str(np.mean(m)))
    n_counts_day_1 = np.random.poisson(m, size =(1, int(N_total)))
    print(n_counts_day_1)
    n_counts_day_1 = n_counts_day_1[0,:]
    
    #print('done')
    
    
    #n_counts_day_0_bis = np.random.choice(n_counts_day_0, N_total)
    #n_counts_day_1_bis = np.random.choice

    #-------------------------------Creation of the data set-------------------------------------
    
    obs=np.logical_or(n_counts_day_0>0, n_counts_day_1>0)
    n1_samples=n_counts_day_0[obs]
    n2_samples=n_counts_day_1[obs]
    pair_samples_df= pd.DataFrame({'Clone_count_1':n1_samples,'Clone_count_2':n2_samples})
    
    pair_samples_df['Clone_frequency_1'] = pair_samples_df['Clone_count_1'] / np.sum(pair_samples_df['Clone_count_1'])
    pair_samples_df['Clone_frequency_2'] = pair_samples_df['Clone_count_2'] / np.sum(pair_samples_df['Clone_count_2'])
    
    
    return pair_samples_df

def experimental_sampling_diffusion_NegBin(NreadsI, NreadsII, paras, x_0, x_2, N_cell_0, N_cell_2):
    
    
    #----------------------------Counts generation --------------------------------------------
    
    ##Initial condition
    N_total_0 = len(x_0[x_0 != -np.inf])
    x_0_bis = x_0[x_0 != -np.inf]
    
    print('Number of clones at initial time ' + str(N_total_0))
    
    N_total_2 = len(x_2[x_2 != -np.inf])
    x_2_bis = x_2[x_2 != -np.inf]
    
    print('Number of clones after 2 years ' + str(N_total_2))
    
    #N_total = min(N_total_0, N_total_2)
    assert len(x_0) == len(x_2)
    N_total = len(x_0)
    

    f_vec_initial = np.exp(x_0)/N_cell_0
    m=float(NreadsI)*f_vec_initial
    print(m)
 
    beta_mv=paras[0]
    alpha_mv=paras[1]
   
    v=m+beta_mv*np.power(m,alpha_mv)

    pvec=1-m/v
    nvec=m*m/v/pvec

    pvec = np.nan_to_num(pvec, nan=0.0)
    nvec = np.nan_to_num(nvec, nan=1e-30)

    n_counts_day_0 = np.random.negative_binomial(nvec, 1-pvec, size =(1, int(N_total)))
    n_counts_day_0 = n_counts_day_0[0,:]
    print(n_counts_day_0)
    
    
    #Final condition
    f_vec_end = np.exp(x_2)/N_cell_2
    m_end=float(NreadsII)*f_vec_end
    print(m_end)

    v_end=m_end+beta_mv*np.power(m_end,alpha_mv)
    pvec_end=1-m_end/v_end
    nvec_end=m_end*m_end/v_end/pvec_end

    pvec_end = np.nan_to_num(pvec_end, nan=0.0)
    nvec_end = np.nan_to_num(nvec_end, nan=1e-30)


    n_counts_day_1 = np.random.negative_binomial(nvec_end, 1-pvec_end, size =(1, int(N_total)))
    n_counts_day_1 = n_counts_day_1[0,:]
    print(n_counts_day_1)
    
    #print('done')
    
    
    #n_counts_day_0_bis = np.random.choice(n_counts_day_0, N_total)
    #n_counts_day_1_bis = np.random.choice

    #-------------------------------Creation of the data set-------------------------------------
    
    obs=np.logical_or(n_counts_day_0>0, n_counts_day_1>0)
    n1_samples=n_counts_day_0[obs]
    n2_samples=n_counts_day_1[obs]
    pair_samples_df= pd.DataFrame({'Clone_count_1':n1_samples,'Clone_count_2':n2_samples})
    
    pair_samples_df['Clone_frequency_1'] = pair_samples_df['Clone_count_1'] / np.sum(pair_samples_df['Clone_count_1'])
    pair_samples_df['Clone_frequency_2'] = pair_samples_df['Clone_count_2'] / np.sum(pair_samples_df['Clone_count_2'])
    
    
    return pair_samples_df


def synthetic_data_generator(B, A, NreadsI, NreadsII, i, t, N_0):



    # EM generation

    #N_0=40
    nt = 100

    #A = -0.25
    #NreadsI = int(1e6)
    #NreadsII = int(1e6)
    alpha = -1 +2*A/B

    st_1 = time.time()

    #1/ Generate log-population at initial time from steady-state distribution + GBM diffusion trajectories for 2 years
    log_nvec_generated, x, N_cells_traj, N_clones_not_extinct_total, N_cells_no_source, N_clones_no_source = generator_diffusion_EM(B, A, N_0, t, nt)

    #print("--- %s seconds ---" % (time.time() - st_1))
    #print('EM' + str(t) + ' year generation done')

    N_clones_ini_EM = N_clones_no_source[0]

    #print('N_clones_ini_EM= ' + str(N_clones_ini_EM) )

    x_i_EM_bis = log_nvec_generated[:N_clones_ini_EM]
    x_f_EM_bis = x[:N_clones_ini_EM]

    # create useful dataframes

    ## traj

    #print('ok')

    dic = {'#clones with source' : N_clones_not_extinct_total, '#clones without any source' : N_clones_no_source, '#cells with source' : N_cells_traj, '#cells without any source' : N_cells_no_source}

    table = pd.DataFrame(data=dic)

    table.to_csv('trajectories_nt_' + str(nt) + str(A) + str(B) +'_alpha_.csv', sep= '\t')


    #4/ Generate data-frame for diffusion
    #N_cells_day_0 = N_cells_traj[0]
    #N_cells_day_1 = N_cells_traj[1] ## LINE TO CHECK 


    N_cells_day_0 = int(1e10)
    N_cells_day_1 = int(1e10)

    df_diffusion_EM  = experimental_sampling_diffusion(NreadsI, NreadsII, x_i_EM_bis, x_f_EM_bis,t,  N_cells_day_0, N_cells_day_1)


    df_diffusion_EM.to_csv('EM_synthetic_sym_Poisson_nt_' + str(nt)+ '_N0_40_A_' + str(A) + '_B' + str(round(B, 2)) + '_' + str(t) + str(i) + '_years' + str(i)  + '.csv' , sep= '\t')


    # LB generation

    x_i_LB, x_f_LB, Prop_Matrix, p_ext, results_extinction  = generator_diffusion_LB(B, A, N_0, t)
    #N_cells_day_0_LB, N_cells_day_1_LB = np.sum(np.exp(x_i_LB)), np.sum(np.exp(x_i_LB))
    N_cells_day_0_LB, N_cells_day_1_LB = int(1e10), int(1e10)
    print('LB' + str(t) + ' year generation done')

        
    print('Number of generated cells LB:' + str(N_cells_day_0_LB)) 

    df_diffusion_LB  = experimental_sampling_diffusion(NreadsI, NreadsII, x_i_LB, x_f_LB, t, N_cells_day_0_LB, N_cells_day_1_LB)

    df_diffusion_LB.to_csv('LB_synthetic_sym_Poisson_N0_40_A_' + str(A) + '_B' + str(round(B, 2)) +'_' +  str(t)  + '_years' + str(i)  + '.csv' , sep= '\t')

def synthetic_data_generator_EM_LB(B, A, i,  NreadsI, NreadsII, t, N_0, method):


    if NreadsI == NreadsII:
        key_sym = 'sym'

    else:
        key_sym = 'asym'

    #Name of the directory

    dirName = 'Synthetic_data_comparison_' + method + '_key_sym_' + str(t) +  '_year(s)'     
    os.makedirs(dirName, exist_ok=True) 

    # EM generation

    #N_0=40
    nt = 100
    paras = [0.7, 1.1]
    #A = -0.25
    #NreadsI = int(1e6)
    #NreadsII = int(1e6)
    alpha = -1 +2*A/B

    st_1 = time.time()

    #1/ Generate log-population at initial time from steady-state distribution + GBM diffusion trajectories for 2 years
    log_nvec_generated, x, N_cells_traj, N_clones_not_extinct_total, N_cells_no_source, N_clones_no_source, x_i_LB, x_f_LB, x_source_LB = generator_diffusion_EM_LB(B, A, N_0, t, nt)

    print("--- %s seconds ---" % (time.time() - st_1))
    #print('EM' + str(t) + ' year generation done')

    N_clones_ini_EM = N_clones_no_source[0]

    print('N_clones_ini_EM= ' + str(N_clones_ini_EM) )

    x_i_EM_bis = log_nvec_generated[:N_clones_ini_EM]
    x_f_EM_bis = x[:N_clones_ini_EM]

    # create useful dataframes

    ## traj

    #print('ok')

    dic = {'#clones with source' : N_clones_not_extinct_total, '#clones without any source' : N_clones_no_source, '#cells with source' : N_cells_traj, '#cells without any source' : N_cells_no_source}

    table = pd.DataFrame(data=dic)

    table.to_csv('trajectories_nt_' + str(nt) + str(A) + str(B) +'_alpha_.csv', sep= '\t')


    #4/ Generate data-frame for diffusion
    N_cells_day_0 = N_cells_traj[0]
    N_cells_day_1 = N_cells_traj[-1] ## LINE TO CHECK 
    print('EM - NUMBER OF CELLS AT INITIAL TIME')
    print(N_cells_day_0)

    print('EM - NUMBER OF CELLS AT FINAL TIME')
    print(N_cells_day_1)


    #5/ EM data-frame generations
    if method == 'negative_binomial':

        df_diffusion_EM  = experimental_sampling_diffusion_NegBin(NreadsI, NreadsII, paras, x_i_EM_bis, x_f_EM_bis, N_cells_day_0, N_cells_day_1)
        df_diffusion_EM.to_csv(dirName + '/EM_synthetic_' + key_sym + '_NegBin_nt_' + str(nt)+ '_N0_'+ str(N_0) + '_A_' + str(A) + '_B_' + str(B) + '_' + str(t) + '_years_' + str(i)  + '.csv' , sep= '\t')

    elif method == 'poisson':

        df_diffusion_EM  = experimental_sampling_diffusion_Poisson(NreadsI, NreadsII, x_i_EM_bis, x_f_EM_bis,t,  N_cells_day_0, N_cells_day_1)
        df_diffusion_EM.to_csv(dirName + '/EM_synthetic_' + key_sym + '_Poisson_nt_' + str(nt)+ '_N0_'+ str(N_0) + '_A_' + str(A) + '_B_' + str(B) + '_' + str(t) + '_years_' + str(i)  + '.csv' , sep= '\t')


    #6/ LB generation

    #x_i_LB, x_f_LB, Prop_Matrix, p_ext, results_extinction  = generator_diffusion_LB(B, A, N_0, t)
    N_cells_day_0_LB, N_cells_day_1_LB = np.sum(np.exp(x_i_LB)), np.sum(np.exp(x_f_LB)) + np.sum(np.exp(x_source_LB))  #N_cells_final_LB
    print('LB - NUMBER OF CELLS AT INITIAL TIME')
    print(N_cells_day_0_LB)

    print('LB - NUMBER OF CELLS AT FINAL TIME')
    print(N_cells_day_1_LB)


    if method == 'negative_binomial':

        df_diffusion_LB  = experimental_sampling_diffusion_NegBin(NreadsI, NreadsII, paras, x_i_LB, x_f_LB, N_cells_day_0_LB, N_cells_day_1_LB)
        df_diffusion_LB.to_csv(dirName + '/LB_synthetic_' + key_sym + '_NegBin_' + 'N0_'+ str(N_0) + '_A_' + str(A) + '_B_' + str(B) +'_' +  str(t)  + '_years_' + str(i)  + '.csv' , sep= '\t')

    elif method == 'poisson': 

        df_diffusion_LB  = experimental_sampling_diffusion_Poisson(NreadsI, NreadsII, x_i_LB, x_f_LB, t, N_cells_day_0_LB, N_cells_day_1_LB)
        df_diffusion_LB.to_csv(dirName + '/LB_synthetic_' + key_sym + '_Poisson_' + 'N0_'+ str(N_0) + '_A_' + str(A) + '_B_' + str(B) +'_' +  str(t)  + '_years_' + str(i)  + '.csv' , sep= '\t')




def synthetic_data_generator_LB_linear_freq(i, A_0, A_1, B_0, B_1,  NreadsI, NreadsII, t, N_0, method):


    print('N_iter ' + str(i)) 
    np.random.seed(i)
    if NreadsI == NreadsII:
        key_sym = '_sym_'

    else:
        key_sym = '_asym_'

    #Name of the directory

    dirName = 'Synthetic_data_linear_freq_' + method + key_sym + str(t) + '_A_0_' + str(A_0) + '_B_0_' + str(B_0) + '_A_1_' + str(A_1) + str(t) + '_year(s)'   
    print('repository name' + dirName)  
    os.makedirs(dirName, exist_ok=True) 

    # EM generation

    #N_0=40
    paras = [0.7, 1.1]
    #A = -0.25
    #NreadsI = int(1e6)
    #NreadsII = int(1e6)
    alpha = -1 +2*A_0/B_0
    print('alpha : ' + str(alpha))

    st_1 = time.time()

    #1/ Generate log-population at initial time from steady-state distribution + GBM diffusion trajectories for 2 years
    x_i_LB, x_f_LB, Prop_Matrix_LB, p_ext_LB, results_extinction_LB, time_vec_LB, results_extinction_source_LB, x_source_LB = generator_diffusion_linear_freq(i, A_0, A_1, B_0, B_1, N_0, t)

    print("--- %s seconds ---" % (time.time() - st_1))
    
    #6/ LB generation

    #x_i_LB, x_f_LB, Prop_Matrix, p_ext, results_extinction  = generator_diffusion_LB(B, A, N_0, t)
    N_cells_day_0_LB, N_cells_day_1_LB = np.sum(np.exp(x_i_LB)), np.sum(np.exp(x_f_LB)) + np.sum(np.exp(x_source_LB))  #N_cells_final_LB
    print('LB - NUMBER OF CELLS AT INITIAL TIME')
    print(N_cells_day_0_LB)

    print('LB - NUMBER OF CELLS AT FINAL TIME')
    print(N_cells_day_1_LB)

    print('SHAPE_X_I ' +  str(np.shape(x_i_LB)))
    print('SHAPE_X_F ' +  str(np.shape(x_f_LB)))


    if method == 'negative_binomial':

        df_diffusion_LB  = experimental_sampling_diffusion_NegBin(NreadsI, NreadsII, paras, x_i_LB, x_f_LB, N_cells_day_0_LB, N_cells_day_1_LB)
        df_diffusion_LB.to_csv(dirName + '/LB_synthetic' + key_sym + 'NegBin_' + 'N0_'+ str(N_0) + '_A_0_' + str(A_0) + '_A_1_' + str(A_1) + '_B_0_' + str(B_0) + '_B_1' + str(B_1) +'_' +  str(t)  + '_years_' + str(i)  + '.csv' , sep= '\t')

    elif method == 'poisson': 

        df_diffusion_LB  = experimental_sampling_diffusion_Poisson(NreadsI, NreadsII, x_i_LB, x_f_LB, t, N_cells_day_0_LB, N_cells_day_1_LB)
        df_diffusion_LB.to_csv(dirName + '/LB_synthetic' + key_sym + 'Poisson_' + 'N0_'+ str(N_0) + '_A_0_' + str(A_0) + '_A_1_' + str(A_1) + '_B_0_' + str(B_0) + '_B_1' + str(B_1) +  str(t)  + '_years_' + str(i)  + '.csv' , sep= '\t')


def synthetic_data_generator_traj(i, A, B, NreadsI, NreadsII, t, N_0, method):

    """
    t : is a vector
    """


    print('N_iter ' + str(i)) 
    np.random.seed(i)
    if NreadsI == NreadsII:
        key_sym = '_sym_'

    else:
        key_sym = '_asym_'

    paras = [0.7, 1.1]

    #Name of the directory

    dirName = 'Synthetic_data_traj_' + method + key_sym + 'A_' + str(A) + '_B_' + str(B) + '_' + str(t[-1]) + '_year(s)'   
    print('repository name' + dirName)  
    os.makedirs(dirName, exist_ok=True) 

    
    st_1 = time.time()

    #1/ Generate log-population at initial time from steady-state distribution + GBM diffusion trajectories for 2 years
    x_i, x_f, Prop_Matrix, p_ext, results_extinction, N_cells_total = generator_diffusion_trajectories(B, A, N_0, t)

    x_i = np.reshape(x_i, (len(x_i), 1))

    print(str(np.shape(x_i)) + 'SHAPE_X_I') 
    print(str(np.shape(x_f)) + 'SHAPE_X_F')

    X_pop = np.concatenate((x_i, x_f), axis = 1)
    print('SHAPE X_pop' + str(np.shape(X_pop)))

    print("--- %s seconds ---" % (time.time() - st_1))
    
    #2/ Experimental sampling generation

    delta_t_vec = np.diff(t)
    
    for k in tqdm(range(len(delta_t_vec))):
        N_cells_day_0_LB, N_cells_day_1_LB = N_cells_total[k], N_cells_total[k+1]  #N_cells_final_LB
        print('LB - NUMBER OF CELLS AT TIME ' + str(k))
        print(N_cells_day_0_LB)

        print('LB - NUMBER OF CELLS AT TIME ' +str(k+1))
        print(N_cells_day_1_LB)

        x_i_LB = X_pop[:, k]
        x_f_LB = X_pop[:, k+1]

        print(x_i_LB == x_f_LB)

        print('SHAPE_X_I ' +  str(np.shape(x_i_LB)))
        print('SHAPE_X_F ' +  str(np.shape(x_f_LB)))


        if method == 'negative_binomial':

            df_diffusion_LB  = experimental_sampling_diffusion_NegBin(NreadsI, NreadsII, paras, x_i_LB, x_f_LB, N_cells_day_0_LB, N_cells_day_1_LB)
            df_diffusion_LB.to_csv(dirName + '/LB_synthetic' + key_sym + 'NegBin_' + 'N0_'+ str(N_0) + '_A_' + str(A) + '_B_' + str(B) + '_step_' + str(k) + '_' + str(delta_t_vec[k])  + '_years_' + str(i)  + '.csv' , sep= '\t')

        elif method == 'poisson': 

            df_diffusion_LB  = experimental_sampling_diffusion_Poisson(NreadsI, NreadsII, x_i_LB, x_f_LB, delta_t_vec[k], N_cells_day_0_LB, N_cells_day_1_LB)
            df_diffusion_LB.to_csv(dirName + '/LB_synthetic' + key_sym + 'Poisson_' + 'N0_'+ str(N_0) + '_A_' + str(A) + '_B_' + str(B) + '_step_' + str(k) + '_' + str(delta_t_vec[k])  + '_years_' + str(i)  + '.csv' , sep= '\t')


def synthetic_data_generator_LB(i, A, B, NreadsI, NreadsII, t, N_0, method):


    print('N_iter ' + str(i)) 
    np.random.seed(i)
    if NreadsI == NreadsII:
        key_sym = '_sym_'

    else:
        key_sym = '_asym_'

    #Name of the directory

    name_B = str(B)
    if len(name_B) == 3:
        name_B = name_B + '0'
    else:
        pass

    dirName = 'Synthetic_data_' + method + key_sym + 'A_' + str(A) + '_B_' + name_B + '_' + str(t) + '_year(s)'   
    print('repository name' + dirName)  
    os.makedirs(dirName, exist_ok=True) 

    # EM generation

    #N_0=40
    paras = [0.7, 1.1]
    #A = -0.25
    #NreadsI = int(1e6)
    #NreadsII = int(1e6)
    alpha = -1 +2*A/B
    print('alpha : ' + str(alpha))

    st_1 = time.time()

    #1/ Generate log-population at initial time from steady-state distribution + GBM diffusion trajectories for 2 years
    x_i_LB, x_f_LB, Prop_Matrix_LB, p_ext_LB, results_extinction_LB, time_vec_LB, results_extinction_source_LB, x_source_LB = generator_diffusion_LB(i, A, B, N_0, t)

    print("--- %s seconds ---" % (time.time() - st_1))
    
    #6/ LB generation

    #x_i_LB, x_f_LB, Prop_Matrix, p_ext, results_extinction  = generator_diffusion_LB(B, A, N_0, t)
    N_cells_day_0_LB, N_cells_day_1_LB = np.sum(np.exp(x_i_LB)), np.sum(np.exp(x_f_LB)) + np.sum(np.exp(x_source_LB))  #N_cells_final_LB
    print('LB - NUMBER OF CELLS AT INITIAL TIME')
    print(N_cells_day_0_LB)

    print('LB - NUMBER OF CELLS AT FINAL TIME')
    print(N_cells_day_1_LB)

    print('SHAPE_X_I ' +  str(np.shape(x_i_LB)))
    print('SHAPE_X_F ' +  str(np.shape(x_f_LB)))


    if method == 'negative_binomial':

        df_diffusion_LB  = experimental_sampling_diffusion_NegBin(NreadsI, NreadsII, paras, x_i_LB, x_f_LB, N_cells_day_0_LB, N_cells_day_1_LB)
        df_diffusion_LB.to_csv(dirName + '/LB_synthetic' + key_sym + 'NegBin_' + 'N0_'+ str(N_0) + '_A_' + str(A) + '_B_' + str(B) + '_' +  str(t)  + '_years_' + str(i)  + '.csv' , sep= '\t')

    elif method == 'poisson': 

        df_diffusion_LB  = experimental_sampling_diffusion_Poisson(NreadsI, NreadsII, x_i_LB, x_f_LB, t, N_cells_day_0_LB, N_cells_day_1_LB)
        df_diffusion_LB.to_csv(dirName + '/LB_synthetic' + key_sym + 'Poisson_' + 'N0_'+ str(N_0) + '_A_' + str(A) +  '_B_' + str(B) + '_' +  str(t)  + '_years_' + str(i)  + '.csv' , sep= '\t')

