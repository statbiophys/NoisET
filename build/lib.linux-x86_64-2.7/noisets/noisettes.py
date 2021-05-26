
"""
Functions library for NoisET - construction of noisettes package

Copyright (C) 2021 Meriem Bensouda Koraichi
   This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""


# Import python libraries
import os
import time
import math
from copy import deepcopy
from decimal import Decimal
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import nbinom
from scipy.stats import poisson
from scipy.stats import rv_discrete
from datetime import datetime, date
from scipy.optimize import minimize

import utils as ut


#===============================Data-Pre-Processing===================================

class Data_Process():

    def __init__(self, path, filename1, filename2, id_labels, count_labels, 
                 delimiters=['\t','\t']):
        """
        Initialize the tables of counts of the two samples. The tables can be tsv files
        or gz-compressed tsv files. 
        Arguments:
            path: string variable specifying the path to the data-set directory
            filename1: string specifying the name of the first sample associated 
                to the first sample
            filename2: string specifying the name of the second sample associated 
                to the second sample
            id_labels: pair of strings which says the coulmn of the table in the two
                files containing the identifier of the count (e.g. the id of the clone 
                CDR3 seuqnece)
            count_labels: pair of strings which says the coulmn of the table in the two
                files containing the counts
        """
        self.path = path
        self.filename1 = filename1
        self.filename2 = filename2
        self.count_labels = count_labels
        self.id_labels = id_labels
        self.delimiters = delimiters
    

    def import_data(self):
        """
        Import the tables with the parameters specified in the constructor
        
        """

        mincount = 0
        maxcount = np.inf
        
        #headerline=0 #line number of headerline
        #newnames=['Clone_fraction','Clone_count','ntCDR3','AACDR3']
        count_label = 'counts'
        id_label = 'id'

        if self.filename1[-2:] == 'gz':
            chunk1=pd.read_csv(self.path+self.filename1, delimiter=self.delimiters[0], 
                               compression = 'gzip')
        else:
            chunk1=pd.read_csv(self.path+self.filename1, delimiter=self.delimiters[0])
        chunk1 = chunk1.loc[:,[self.id_labels[0], self.count_labels[0]]]

        if self.filename2[-2:] == 'gz':
            chunk2=pd.read_csv(self.path+self.filename2, delimiter=self.delimiters[1], 
                               compression = 'gzip')
        else:
            chunk2=pd.read_csv(self.path+self.filename2, delimiter=self.delimiters[1])
        chunk2 = chunk2.loc[:,[self.id_labels[1], self.count_labels[1]]]

        chunk1 = chunk1.rename(columns={self.id_labels[0]: id_label, 
                                        self.count_labels[0]: count_label})
        chunk2 = chunk2.rename(columns={self.id_labels[1]: id_label, 
                                        self.count_labels[1]: count_label})

        suffixes=('1','2')
        mergedFrame=pd.merge(chunk1, chunk2, on=id_label, suffixes=suffixes, how='outer')

        for labelit in suffixes:
            mergedFrame.loc[:,count_label+labelit].fillna(int(0),inplace=True)
            mergedFrame.loc[:,count_label+labelit].astype(int)

        # Check for duplicates??

        #def dummy(x):
        #    val=x[0]
        #    if pd.isnull(val):
        #        val=x[1]    
        #    return val

        # Ask Meriem what's happening next:
        #mergedFrame.loc[:,newnames[3]+suffixes[0]]=mergedFrame.loc[:,[newnames[3]+suffixes[0],newnames[3]+suffixes[1]]].apply(dummy,axis=1) #assigns AA sequence to clones, creates duplicates
        #mergedFrame.drop(newnames[3]+suffixes[1], 1,inplace=True) #removes duplicates
        #mergedFrame.rename(columns = {newnames[3]+suffixes[0]:newnames[3]}, inplace = True)
        #mergedFrame=mergedFrame[[newname+suffix for newname in newnames[:2] for suffix in suffixes]+[newnames[2],newnames[3]]]
        
        # TODO: ask in input mincount and maxcount
        #filterout=((mergedFrame.Clone_count_1<mincount) & (mergedFrame.Clone_count_2==0)) | ((mergedFrame.Clone_count_2<mincount) & (mergedFrame.Clone_count_1==0)) #has effect only if mincount>0
        #number_clones=len(mergedFrame)
        #return number_clones,mergedFrame.loc[((mergedFrame.Clone_count_1<=maxcount) & (mergedFrame.Clone_count_2<=maxcount)) & ~filterout]
        
        return mergedFrame


#===============================Noise-Model=========================================

class Noise_Model:


    """ Noise_Model:
    creation of an oject + methods to learn null noise model"""


    def __init__(self, data, count_labels=['counts1', 'counts2'], verbose=True):
        """
        Compute the noise model given the data used as input.
        - pandas.DataFrame having two colums with label specified in count_labels.
        - a 2-dim list which contains all the pairs [n1, n2] of counts.
        """

        df = ut.data_parser(data, count_labels)
        self.verbose = verbose
        self.sparserep = ut.get_sparserep(df)


    def learn_null_model(self, noise_model, init_paras=[],  print_to_path=None):
        """
        Performs constrained maximization of null model likelihood for finding the 
        noise-model parameters of the data.
        Arguments:
            noise_model: 0 for a negative binomial poisson, 1 for a negative binomial
                2 for a poisson.
            init_paras: initial guesses of the parameters. If empty list a default init 
                condition is chosen
            print_to_path: if a string, it prints the parameters to the specified path
        """

        # constraint type 1 gives only low error modes, see paper for details.
        constr_type = 1

        # Choice of the model:
        # Parameters initialization depending on the model 
        if noise_model < 1:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'm_total', 'fmin']
            if len(init_paras)==0:
                init_paras = np.array([-2.0, 1.5, 1.2, 6.7, -9.7])
        elif noise_model == 1:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'fmin']
            if len(init_paras)==0:
                init_paras = np.array([-2.0, 0.45, 1.1, -10.0])
        else:
            parameter_labels = ['alph_rho', 'fmin']
            if len(init_paras)==0:
                init_paras = np.array([-2.1,  -9.5])

        assert len(parameter_labels) == len(init_paras), "number of model and initial paras differ!"

        condict = {'type': 'eq', 'fun': self._nullmodel_constr_fn, 'args': (noise_model, constr_type)}


        partialobjfunc = partial(self._get_Pn1n2, noise_model=noise_model)
        nullfunctol = 1e-6
        nullmaxiter = 200

        if self.verbose:
            header = ['Iter'] + parameter_labels
            print(''.join(['{' + str(it) + ':9s} ' for it in range(len(init_paras) + 1)]).format(*header))
            
        global curr_iter
        curr_iter = 1
        callbackp = partial(self._callback, nparas=len(init_paras), noise_model= noise_model)
        outstruct = minimize(partialobjfunc, init_paras, method='SLSQP', callback=callbackp, constraints=condict,
                        options={'ftol': nullfunctol, 'disp': True, 'maxiter': nullmaxiter})
            
        constr_value = self._nullmodel_constr_fn(outstruct.x, noise_model, constr_type)

        if noise_model < 1:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'm_total', 'fmin']
        elif noise_model == 1:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'fmin']
        else:
            parameter_labels = ['alph_rho', 'fmin']
        self.params = {'label' : parameter_labels, 'value': outstruct.x}
        df = pd.DataFrame(data = self.params)

        if type(print_to_path) == str:
            df.to_csv(print_to_path, sep = '\t')

        #np.save('nullpara' + str(noise_model), outstruct.x)

        return outstruct, constr_value



    #-----------------------------Null-Model-optimization--------------------------
        
    def _get_Pn1n2(self, paras, noise_model):

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,NreadsI,NreadsII = self.sparserep
            
        nfbins = 1200
        freq_dtype = float

        # Parameters

        alpha = paras[0]
        fmin = np.power(10,paras[-1])

        # 
        logrhofvec, logfvec, normconst = ut.get_rhof(alpha,nfbins,fmin,freq_dtype)

        # 

        logfvec_tmp=deepcopy(logfvec)

        logPn1_f = ut.get_logPn_f(unicountvals_1, NreadsI,logfvec_tmp, noise_model, paras)
        logPn2_f = ut.get_logPn_f(unicountvals_2, NreadsII,logfvec_tmp, noise_model, paras)

        # for the trapezoid integral methods

        dlogfby2=np.diff(logfvec)/2

        # Compute P(0,0) for the normalization constraint
        integ = np.exp(logrhofvec + logPn2_f[:, 0] + logPn1_f[:, 0] + logfvec)
        Pn0n0 = np.dot(dlogfby2, integ[1:] + integ[:-1])

        #print("computing P(n1,n2)")
        Pn1n2 = np.zeros(len(sparse_rep_counts))  # 1D representation
        for it, (ind1, ind2) in enumerate(zip(indn1, indn2)):
            integ = np.exp(logPn1_f[:, ind1] + logrhofvec + logPn2_f[:, ind2] + logfvec)
            Pn1n2[it] = np.dot(dlogfby2, integ[1:] + integ[:-1])
        Pn1n2 /= 1. - Pn0n0  # renormalize
        return -np.dot(sparse_rep_counts, np.where(Pn1n2 > 0, np.log(Pn1n2), 0)) / float(np.sum(sparse_rep_counts))


    def _callback(self, paras, nparas, noise_model):
        '''prints iteration info. called by scipy.minimize'''

        global curr_iter
        #curr_iter = 0
        global Loss_function 
        if self.verbose:
            print(''.join(['{0:d} ']+['{'+str(it)+':3.6f} ' for it in range(1,len(paras)+1)]).format(*([curr_iter]+list(paras))))
            #print ('{' + str(len(paras)+1) + ':3.6f}'.format( [self.get_Pn1n2(paras, sparse_rep, acq_model_type)]))
            Loss_function = self._get_Pn1n2(paras, noise_model)
            print(Loss_function)
        curr_iter += 1
        

    # Constraints for the Null-Model, no filtered 
    def _nullmodel_constr_fn(self, paras, noise_model, constr_type):
            
        '''
        returns either or both of the two level-set functions: log<f>-log(1/N), with N=Nclones/(1-P(0,0)) and log(Z_f), with Z_f=N<f>_{n+n'=0} + sum_i^Nclones <f>_{f|n,n'}
        '''

    # Choice of the model: 

        indn1, indn2, sparse_rep_counts, unicountvals_1, unicountvals_2, NreadsI, NreadsII = self.sparserep

        #Variables that would be chosen in the future by the user 
        nfbins = 1200
        freq_dtype = float

        alpha = paras[0]  # power law exponent
        fmin = np.power(10, paras[-1]) # true minimal frequency 

        logrhofvec, logfvec, normconst = ut.get_rhof(alpha,nfbins,fmin,freq_dtype)
        dlogfby2 = np.diff(logfvec) / 2.  # 1/2 comes from trapezoid integration below

        integ = np.exp(logrhofvec + 2 * logfvec)
        avgf_ps = np.dot(dlogfby2, integ[:-1] + integ[1:])

        logPn1_f = ut.get_logPn_f(unicountvals_1, NreadsI, logfvec, noise_model, paras)
        logPn2_f = ut.get_logPn_f(unicountvals_2, NreadsII, logfvec, noise_model, paras)

        integ = np.exp(logPn1_f[:, 0] + logPn2_f[:, 0] + logrhofvec + logfvec)
        Pn0n0 = np.dot(dlogfby2, integ[1:] + integ[:-1])
        logPnng0 = np.log(1 - Pn0n0)
        avgf_null_pair = np.exp(logPnng0 - np.log(np.sum(sparse_rep_counts)))

        C1 = np.log(avgf_ps) - np.log(avgf_null_pair)

        integ = np.exp(logPn1_f[:, 0] + logPn2_f[:, 0] + logrhofvec + 2 * logfvec)
        log_avgf_n0n0 = np.log(np.dot(dlogfby2, integ[1:] + integ[:-1]))

        integ = np.exp(logPn1_f[:, indn1] + logPn2_f[:, indn2] + logrhofvec[:, np.newaxis] + logfvec[:, np.newaxis])
        log_Pn1n2 = np.log(np.sum(dlogfby2[:, np.newaxis] * (integ[1:, :] + integ[:-1, :]), axis=0))
        integ = np.exp(np.log(integ) + logfvec[:, np.newaxis])
        tmp = deepcopy(log_Pn1n2)
        tmp[tmp == -np.Inf] = np.Inf  # since subtracted in next line
        avgf_n1n2 = np.exp(np.log(np.sum(dlogfby2[:, np.newaxis] * (integ[1:, :] + integ[:-1, :]), axis=0)) - tmp)
        log_sumavgf = np.log(np.dot(sparse_rep_counts, avgf_n1n2))

        logNclones = np.log(np.sum(sparse_rep_counts)) - logPnng0
        Z = np.exp(logNclones + np.log(Pn0n0) + log_avgf_n0n0) + np.exp(log_sumavgf)

        C2 = np.log(Z)
        
        # print('C1:'+str(C1)+' C2:'+str(C2))
        if constr_type == 0:
            return C1
        elif constr_type == 1:
            return C2
        else:
            return C1, C2



#============================================Differential expression =============================================================

class Expansion_Model:
    
    """
    Explain Methods for this class
    """


    def __init__(self, data, count_labels=['counts1', 'counts2'], 
                 id_label='id', verbose=True):
        """
        Compute the noise model given the data used as input.
        - pandas.DataFrame having two colums called "counts1" and "counts2".
        - a 2-dim list which contains all the pairs [n1, n2] of counts.
        """

        df = ut.data_parser(data, count_labels, id_label)
        self.verbose= verbose
        self.sparserep = ut.get_sparserep(df)
        self.df = df


    def _get_Ps(self, alp, sbar, smax, stp):
        '''
        generates symmetric exponential distribution over log fold change
        with effect size sbar and nonresponding fraction 1-alp at s=0.
        computed over discrete range of s from -smax to smax in steps of size stp
        '''
        lamb=-stp/sbar
        smaxt=round(smax/stp)
        s_zeroind=int(smaxt)
        Z=2*(np.exp((smaxt+1)*lamb)-1)/(np.exp(lamb)-1)-1
        Ps=alp*np.exp(lamb*np.fabs(np.arange(-smaxt,smaxt+1)))/Z
        Ps[s_zeroind]+=(1-alp)
        return Ps

    def _callbackFdiffexpr(self, Xi): #case dependent
        '''prints iteration info. called scipy.minimize'''
        
        if self.verbose:
            print('{0: 3.6f}   {1: 3.6f}   '.format(Xi[0], Xi[1])+'\n')   
    

    def learning_dynamics_expansion_polished(self, paras_1, paras_2,  noise_model):
        """
        Different uses for this function to explain, explain the use of NreadsItrue and NreadsIItrue
        paras : 
        sparse_rep :
        noise_model : 
        not_filtered : True if all the repertoire is used to infer the dynamics paramers, false if data is filtered
        NreadsItrue : the total number of reads in the original sample at time 1
        NreadsII true: the total number of reads in the original sample at time 2
        time_unit : day, months or years
        time_1 : indication of the first sample extraction time
        time_2 : indication of the second sample extraction time
        """

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,NreadsI,NreadsII = self.sparserep

        alpha_rho = paras_1[0]
        fmin = np.power(10,paras_1[-1])
        freq_dtype = 'float64'
        nfbins = 1200 #Accuracy of the integration

        logrhofvec, logfvec, normconst = ut.get_rhof(alpha_rho, nfbins, fmin, freq_dtype)

        #Definition of svec
        smax = 25.0     #maximum absolute logfold change value
        s_step = 0.1
        s_0 = -1
        
        s_step_old= s_step
        logf_step= logfvec[1] - logfvec[0] #use natural log here since f2 increments in increments in exp().  
        f2s_step= int(round(s_step/logf_step)) #rounded number of f-steps in one s-step
        s_step= float(f2s_step)*logf_step
        smax= s_step*(smax/s_step_old)
        svec= s_step*np.arange(0,int(round(smax/s_step)+1))   
        svec= np.append(-svec[1:][::-1],svec)

        smaxind=(len(svec)-1)/2
        f2s_step=int(round(s_step/logf_step)) #rounded number of f-steps in one s-step
        logfmin=logfvec[0 ]-f2s_step*smaxind*logf_step
        logfmax=logfvec[-1]+f2s_step*smaxind*logf_step
        
        logfvecwide = np.linspace(logfmin,logfmax,len(logfvec)+2*smaxind*f2s_step) #a wider domain for the second frequency f2=f1*exp(s)
            
        # Compute P(n1|f) and P(n2|f), each in an iteration of the following loop

        for it in range(2):
            if it == 0:
                unicounts=unicountvals_1
                logfvec_tmp=deepcopy(logfvec)
                Nreads = NreadsI
                paras = paras_1
            else:
                unicounts=unicountvals_2
                logfvec_tmp=deepcopy(logfvecwide) #contains s-shift for sampled data method
                Nreads = NreadsII
                paras = paras_2
            if it == 0:
                logPn1_f = ut.get_logPn_f( unicounts, Nreads, logfvec_tmp, noise_model, paras)

            else:
                logPn2_f = ut.get_logPn_f(unicounts, Nreads, logfvec_tmp, noise_model, paras)

        #for the trapezoid method
        dlogfby2=np.diff(logfvec)/2 

        # Computing P(n1,n2|f,s)
        Pn1n2_s=np.zeros((len(svec), len(unicountvals_1), len(unicountvals_2))) 

        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1,indn2)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec )
                Pn1n2_s[s_it, n1_it, n2_it] = np.dot(dlogfby2,integ[1:] + integ[:-1])
            
    
        Pn0n0_s = np.zeros(svec.shape)
        for s_it,s in enumerate(svec):    
            integ=np.exp(logPn1_f[:,0]+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),0]+logrhofvec+logfvec)
            Pn0n0_s[s_it]=np.dot(dlogfby2,integ[1:]+integ[:-1])
            
    
        N_obs = np.sum(sparse_rep_counts)
        if self.verbose:
            print("N_obs: " + str(N_obs))
    
            
        def cost(PARAS):

            alp = PARAS[0]
            sbar = PARAS[1]

            Ps = _get_Ps(self,alp,sbar,smax,s_step)
            Pn0n0=np.dot(Pn0n0_s,Ps)
            Pn1n2_ps=np.sum(Pn1n2_s*Ps[:,np.newaxis,np.newaxis],0)
            Pn1n2_ps/=1-Pn0n0
            if self.verbose:
                print(Pn0n0)

       

            Energy = - np.dot(sparse_rep_counts/float(N_obs),np.where(Pn1n2_ps[indn1,indn2]>0,np.log(Pn1n2_ps[indn1,indn2]),0)) 
                
            return Energy

    #--------------------------Compute-the-grid-----------------------------------------
        
        if self.verbose:
            print('Calculation Surface : \n')
        st = time.time()

        npoints = 20 #to be chosen by the user 
        alpvec = np.logspace(-3,np.log10(0.99), npoints)
        sbarvec = np.linspace(0.01,5, npoints)

        LSurface =np.zeros((len(sbarvec),len(alpvec)))
        for i in range(len(sbarvec)):
            for j in range(len(alpvec)):
                LSurface[i, j]=  - cost([alpvec[j], sbarvec[i]])
        
        alpmesh, sbarmesh = np.meshgrid(alpvec, sbarvec)
        a,b = np.where(LSurface == np.max(LSurface))
        if self.verbose:
            print("--- %s seconds ---" % (time.time() - st))
    
    
    #------------------------------Optimization----------------------------------------------
        
        optA = alpmesh[a[0],b[0]]
        optB = sbarmesh[a[0],b[0]]

        if self.verbose:
            print('polish parameter estimate from '+ str(optA)+' '+str(optB))
        initparas=(optA,optB)  
    

        outstruct = minimize(cost, initparas, method='SLSQP', callback=_callbackFdiffexpr, tol=1e-6,options={'ftol':1e-8 ,'disp': True,'maxiter':300})

        return outstruct.x, Pn1n2_s, Pn0n0_s, svec


    def learning_dynamics_expansion(self, paras_1, paras_2, noise_model, display_plot=False):
        """
        Different uses for this function to explain, explain the use of NreadsItrue and NreadsIItrue
        paras : 
        sparse_rep :
        noise_model: 
        """

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,NreadsI,NreadsII = self.sparserep

        alpha_rho = paras_1[0]
        fmin = np.power(10,paras_1[-1])
        freq_dtype = 'float64'
        nfbins = 1200 #Accuracy of the integration

        logrhofvec, logfvec, normconst = ut.get_rhof(alpha_rho, nfbins, fmin, freq_dtype)

        #Definition of svec
        smax = 25.0     #maximum absolute logfold change value
        s_step = 0.1
        s_0 = -1
        
        s_step_old= s_step
        logf_step= logfvec[1] - logfvec[0] #use natural log here since f2 increments in increments in exp().  
        f2s_step= int(round(s_step/logf_step)) #rounded number of f-steps in one s-step
        s_step= float(f2s_step)*logf_step
        smax= s_step*(smax/s_step_old)
        svec= s_step*np.arange(0,int(round(smax/s_step)+1))   
        svec= np.append(-svec[1:][::-1],svec)

        smaxind=(len(svec)-1)/2
        f2s_step=int(round(s_step/logf_step)) #rounded number of f-steps in one s-step
        logfmin=logfvec[0 ]-f2s_step*smaxind*logf_step
        logfmax=logfvec[-1]+f2s_step*smaxind*logf_step
        
        logfvecwide = np.linspace(logfmin,logfmax,int(len(logfvec)+2*smaxind*f2s_step)) #a wider domain for the second frequency f2=f1*exp(s)
            
        # Compute P(n1|f) and P(n2|f), each in an iteration of the following loop

        for it in range(2):
            if it == 0:
                unicounts=unicountvals_1
                logfvec_tmp=deepcopy(logfvec)
                Nreads = NreadsI
                paras = paras_1
            else:
                unicounts=unicountvals_2
                logfvec_tmp=deepcopy(logfvecwide) #contains s-shift for sampled data method
                Nreads = NreadsII
                paras = paras_2
            if it == 0:
                logPn1_f = ut.get_logPn_f(unicounts, Nreads, logfvec_tmp, noise_model, paras)

            else:
                logPn2_f = ut.get_logPn_f(unicounts, Nreads, logfvec_tmp, noise_model, paras)

        #for the trapezoid method
        dlogfby2=np.diff(logfvec)/2 

        # Computing P(n1,n2|f,s)
        Pn1n2_s=np.zeros((len(svec), len(unicountvals_1), len(unicountvals_2))) 

        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1,indn2)):
                integ = np.exp(logrhofvec+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+logPn1_f[:,n1_it]+ logfvec )
                Pn1n2_s[s_it, n1_it, n2_it] = np.dot(dlogfby2,integ[1:] + integ[:-1])
            
    
        Pn0n0_s = np.zeros(svec.shape)
        for s_it,s in enumerate(svec):    
            integ=np.exp(logPn1_f[:,0]+logPn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),0]+logrhofvec+logfvec)
            Pn0n0_s[s_it]=np.dot(dlogfby2,integ[1:]+integ[:-1])
            
   
        N_obs = np.sum(sparse_rep_counts)
        if self.verbose:
            print("N_obs: " + str(N_obs))
    
            
        def cost(PARAS):

            alp = PARAS[0]
            sbar = PARAS[1]

            Ps = self._get_Ps(alp,sbar,smax,s_step)
            Pn0n0=np.dot(Pn0n0_s,Ps)
            Pn1n2_ps=np.sum(Pn1n2_s*Ps[:,np.newaxis,np.newaxis],0)
            Pn1n2_ps/=1-Pn0n0
            #print(Pn0n0)

       

            Energy = - np.dot(sparse_rep_counts/float(N_obs),np.where(Pn1n2_ps[indn1,indn2]>0,np.log(Pn1n2_ps[indn1,indn2]),0)) 
                
            return Energy

    #--------------------------Compute-the-grid-----------------------------------------
        if self.verbose:
            print('Calculation Surface : \n')
        st = time.time()

        npoints = 50 #to be chosen by the user 
        alpvec = np.logspace(-3,np.log10(0.99), npoints)
        sbarvec = np.linspace(0.01,5, npoints)

        LSurface =np.zeros((len(sbarvec),len(alpvec)))
        for i in range(len(sbarvec)):
            for j in range(len(alpvec)):
                LSurface[i, j]=  - cost([alpvec[j], sbarvec[i]])
        
        alpmesh, sbarmesh = np.meshgrid(alpvec, sbarvec)
        a,b = np.where(LSurface == np.max(LSurface))
        if self.verbose:
            print("--- %s seconds ---" % (time.time() - st))
    
    #---------------------------Plot-the-grid-------------------------------------------
        if display_plot:

            fig, ax =plt.subplots(1, figsize=(10,8))

         
            a,b = np.where(LSurface == np.max(LSurface))

            ax.contour(alpmesh, sbarmesh, LSurface, linewidths=1, colors='k', linestyles = 'solid')
            plt.contourf(alpmesh, sbarmesh, LSurface, 20, cmap = 'viridis', alpha= 0.8)

            xmax = alpmesh[a[0],b[0]]
            ymax = sbarmesh[a[0],b[0]]
            text= r"$ alpha={:.3f}, s={:.3f} $".format(xmax, ymax)
            bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
            arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=80")
            kw = dict(xycoords='data',textcoords="axes fraction",
                arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
            plt.annotate(text, xy=(xmax, ymax), xytext=(0.94,0.96), **kw)
            plt.xlabel(r'$ \alpha, \ size \ of \ the \ repertoire \ that \ answers \ to \ the \ vaccine $') 
            plt.ylabel(r'$ s_{bar}, \ characteristic \ expansion \ decrease $')
            plt.xscale('log')
            plt.yscale('log')
            plt.grid()
            plt.title(r'$Grid \ Search \ graph \ for \ \alpha \ and \ s_{bar} \ parameters. $')
            plt.colorbar()

        return LSurface, Pn1n2_s, Pn0n0_s, svec
 

    def build_table(self, svec, Ps, Pn1n2_s, Pn0n0_s, pthresh, smedthresh):
        '''
        takes learned diffexpr model, Pn1n2_s*Ps, computes posteriors over (n1,n2) pairs, and writes to file a table of data with clones as rows and columns as measures of thier posteriors 
        print_expanded=True orders table as ascending by , else descending
        pthresh is the threshold in 'p-value'-like (null hypo) probability, 1-P(s>0|n1_i,n2_i), where i is the row (i.e. the clone) n.b. lower null prob implies larger probability of expansion
        smedthresh is the threshold on the posterior median, below which clones are discarded
        '''

        indn1_d,indn2_d,sparse_rep_counts, unicountvals_1_d, unicountvals_2_d, NreadsI, NreadsII = self.sparserep
        Psn1n2_ps=Pn1n2_s*Ps[:,np.newaxis,np.newaxis] 
    
        #compute marginal likelihood (neglect renormalization , since it cancels in conditional below) 
        Pn1n2_ps=np.sum(Psn1n2_ps,0)

        Ps_n1n2ps=Pn1n2_s*Ps[:,np.newaxis,np.newaxis]/Pn1n2_ps[np.newaxis,:,:]
        #compute cdf to get p-value to threshold on to reduce output size
        cdfPs_n1n2ps=np.cumsum(Ps_n1n2ps,0)
    

        def dummy(row, cdfPs_n1n2ps, unicountvals_1_d, unicountvals_2_d):
            '''
            when applied to dataframe, generates 'p-value'-like (null hypo) probability, 1-P(s>0|n1_i,n2_i), where i is the row (i.e. the clone)
            '''
            return cdfPs_n1n2ps[np.argmin(np.fabs(svec)),row['counts1']==unicountvals_1_d,row['counts2']==unicountvals_2_d][0]
        
        dummy_part=partial(dummy,cdfPs_n1n2ps=cdfPs_n1n2ps,unicountvals_1_d=unicountvals_1_d,unicountvals_2_d=unicountvals_2_d)
    
        cdflabel=r'$1-P(s>0)$'
        subset = self.df
        subset[cdflabel]=subset.apply(dummy_part, axis=1)
        subset=subset[subset[cdflabel]<pthresh].reset_index(drop=True)

        #go from clone count pair (n1,n2) to index in unicountvals_1_d and unicountvals_2_d
        data_pairs_ind_1=np.zeros((len(subset),),dtype=int)
        data_pairs_ind_2=np.zeros((len(subset),),dtype=int)
        for it in range(len(subset)):
            data_pairs_ind_1[it]=np.where(int(subset.iloc[it].counts1)==unicountvals_1_d)[0]
            data_pairs_ind_2[it]=np.where(int(subset.iloc[it].counts2)==unicountvals_2_d)[0]   
        #posteriors over data clones
        Ps_n1n2ps_datpairs=Ps_n1n2ps[:,data_pairs_ind_1,data_pairs_ind_2]
    
        #compute posterior metrics
        mean_est=np.zeros((len(subset),))
        max_est= np.zeros((len(subset),))
        slowvec= np.zeros((len(subset),))
        smedvec= np.zeros((len(subset),))
        shighvec=np.zeros((len(subset),))
        pval=0.025 #double-sided comparison statistical test
        pvalvec=[pval,0.5,1-pval] #bound criteria defining slow, smed, and shigh, respectively
        for it,column in enumerate(np.transpose(Ps_n1n2ps_datpairs)):
            mean_est[it]=np.sum(svec*column)
            max_est[it]=svec[np.argmax(column)]
            forwardcmf=np.cumsum(column)
            backwardcmf=np.cumsum(column[::-1])[::-1]
            inds=np.where((forwardcmf[:-1]<pvalvec[0]) & (forwardcmf[1:]>=pvalvec[0]))[0]
            slowvec[it]=np.mean(svec[inds+np.ones((len(inds),),dtype=int)])  #use mean in case there are two values
            inds=np.where((forwardcmf>=pvalvec[1]) & (backwardcmf>=pvalvec[1]))[0]
            smedvec[it]=np.mean(svec[inds])
            inds=np.where((forwardcmf[:-1]<pvalvec[2]) & (forwardcmf[1:]>=pvalvec[2]))[0]
            shighvec[it]=np.mean(svec[inds+np.ones((len(inds),),dtype=int)])
    
        colnames=(r'$\bar{s}$',r'$s_{max}$',r'$s_{3,high}$',r'$s_{2,med}$',r'$s_{1,low}$')
        for it,coldata in enumerate((mean_est,max_est,shighvec,smedvec,slowvec)):
            subset.insert(0,colnames[it],coldata)
        #oldcolnames=( 'AACDR3',  'ntCDR3', 'Clone_count_1', 'Clone_count_2', 'Clone_fraction_1', 'Clone_fraction_2')
        #newcolnames=('CDR3_AA', 'CDR3_nt',        r'$n_1$',        r'$n_2$',           r'$f_1$',           r'$f_2$')
        #subset=subset.rename(columns=dict(zip(oldcolnames, newcolnames)))
    
        #select only clones whose posterior median pass the given threshold
        return subset[subset[r'$s_{2,med}$']>smedthresh].sort_values(by=cdflabel,ascending=True)
        
        # Ask Meriem about contracted 

        # print("writing to: "+outpath)
        # if print_expanded:
        #     subset=subset.sort_values(by=cdflabel,ascending=True)
        #     strout='expanded'
        # else:
        #     subset=subset.sort_values(by=cdflabel,ascending=False)
        #     strout='contracted'
        # subset.to_csv(outpath+'top_'+strout+'.csv',sep='\t',index=False)



    def expansion_table(self, paras_1, paras_2, noise_model, pval_threshold, smed_threshold):

        L_surface, Pn1n2_s_d, Pn0n0_s_d, svec = self.learning_dynamics_expansion(paras_1, paras_2, noise_model)
        npoints= 50 # same as in learning_dynamics_expansion
        smax = 25.0     
        s_step = 0.1
        alpvec = np.logspace(-3,np.log10(0.99), npoints)
        sbarvec = np.linspace(0.01,5, npoints)
        maxinds=np.unravel_index(np.argmax(L_surface),np.shape(L_surface))
        optsbar=sbarvec[maxinds[0]]
        optalp=alpvec[maxinds[1]]
        optPs= self._get_Ps(optalp,optsbar,smax,s_step)

        return self.build_table(svec, optPs, Pn1n2_s_d, Pn0n0_s_d, pval_threshold, smed_threshold)


#============================================Generate Synthetic Data =============================================================

class Generator:


    # def get_rhof(self, alpha_rho, fmin, freq_nbins=800, freq_dtype='float64'):

    #     '''
    #     generates power law (power is alpha_rho) clone frequency distribution over 
    #     freq_nbins discrete logarithmically spaced frequences between fmin and 1 of dtype freq_dtype
    #     Outputs log probabilities obtained at log frequencies'''
    #     fmax=1e0
    #     logfvec=np.linspace(np.log10(fmin),np.log10(fmax),freq_nbins)
    #     logfvec=np.array(np.log(np.power(10,logfvec)) ,dtype=freq_dtype).flatten()  
    #     logrhovec=logfvec*alpha_rho
    #     integ=np.exp(logrhovec+logfvec,dtype=freq_dtype)
    #     normconst=np.log(np.dot(np.diff(logfvec)/2.,integ[1:]+integ[:-1]))
    #     logrhovec-=normconst 
    #     return logrhovec,logfvec

    def get_distsample(self, pmf, Nsamp, dtype='uint32'):
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
        choice = np.random.uniform(high = cmf[-1], size = int(float(Nsamp)))
        index = np.searchsorted(cmf, choice)
        index = sortindex[index]
        index = np.unravel_index(index, shape)
        index = np.transpose(np.vstack(index))
        sampled_inds = np.array(index[np.argsort(index[:,0])],dtype=dtype)
        return sampled_inds

    
    def gen_synthetic_data_Null(self, paras, noise_model, NreadsI, NreadsII, Nsamp):
        '''
        outputs an array of observed clone frequencies and corresponding dataframe of pair counts
        for a null model learned from a dataset pair with NreadsI and NreadsII number of reads, respectively.
        Crucial for RAM efficiency, sampling is conditioned on being observed in each of the three (n,0), (0,n'), and n,n'>0 conditions
        so that only Nsamp clones need to be sampled, rather than the N clones in the repertoire.
        Note that no explicit normalization is applied. It is assumed that the values in paras are consistent with N<f>=1 
        (e.g. were obtained through the learning done in this package).
        '''

    
        alpha = paras[0] #power law exponent
        fmin=np.power(10,paras[-1])
        if noise_model<1:
            m_total=float(np.power(10, paras[3])) 
            r_c1=NreadsI/m_total
            r_c2=NreadsII/m_total
            r_cvec=[r_c1,r_c2]
        if noise_model<2:
            beta_mv= paras[1]
            alpha_mv=paras[2]
    
        logrhofvec,logfvec,n = ut.get_rhof(alpha, 800, fmin)
        fvec=np.exp(logfvec)
        dlogf=np.diff(logfvec)/2.
    
        #generate measurement model distribution, Pn_f
        Pn_f=np.empty((len(logfvec),),dtype=object) #len(logfvec) samplers
    
        #get value at n=0 to use for conditioning on n>0 (and get full Pn_f here if noise_model=1,2)
        m_max=1e3 #conditioned on n=0, so no edge effects
    
        Nreadsvec=(NreadsI,NreadsII)
        for it in range(2):
            Pn_f=np.empty((len(fvec),),dtype=object)
            if noise_model==2:
                m1vec=Nreadsvec[it]*fvec
                for find,m1 in enumerate(m1vec):
                    Pn_f[find]=poisson(m1)
                logPn0_f=-m1vec
            elif noise_model==1:
                m1=Nreadsvec[it]*fvec
                v1=m1+beta_mv*np.power(m1,alpha_mv)
                p=1-m1/v1
                n=m1*m1/v1/p
                for find,(n,p) in enumerate(zip(n,p)):
                    Pn_f[find]=nbinom(n,1-p)
                Pn0_f=np.asarray([Pn_find.pmf(0) for Pn_find in Pn_f])
                logPn0_f=np.log(Pn0_f)
            
            elif noise_model==0:
                m1=m_total*fvec
                v1=m1+beta_mv*np.power(m1,alpha_mv)
                p=1-m1/v1
                n=m1*m1/v1/p
                Pn0_f=np.zeros((len(fvec),))
                for find in range(len(Pn0_f)):
                    nbtmp=nbinom(n[find],1-p[find]).pmf(np.arange(m_max+1))
                    ptmp=poisson(r_cvec[it]*np.arange(m_max+1)).pmf(0)
                    Pn0_f[find]=np.sum(np.exp(np.log(nbtmp)+np.log(ptmp)))
                logPn0_f=np.log(Pn0_f)
            else:
                print('acq_model is 0,1,or 2 only')
            
            if it==0:
                Pn1_f=Pn_f
                logPn10_f=logPn0_f
            else:
                Pn2_f=Pn_f
                logPn20_f=logPn0_f

        #3-quadrant q|f conditional distribution (qx0:n1>0,n2=0;q0x:n1=0,n2>0;qxx:n1,n2>0)
        logPqx0_f=np.log(1-np.exp(logPn10_f))+logPn20_f
        logPq0x_f=logPn10_f+np.log(1-np.exp(logPn20_f))
        logPqxx_f=np.log(1-np.exp(logPn10_f))+np.log(1-np.exp(logPn20_f))
        #3-quadrant q,f joint distribution
        logPfqx0=logPqx0_f+logrhofvec
        logPfq0x=logPq0x_f+logrhofvec
        logPfqxx=logPqxx_f+logrhofvec
        #3-quadrant q marginal distribution 
        Pqx0=np.trapz(np.exp(logPfqx0+logfvec),x=logfvec)
        Pq0x=np.trapz(np.exp(logPfq0x+logfvec),x=logfvec)
        Pqxx=np.trapz(np.exp(logPfqxx+logfvec),x=logfvec)
    
        #3 quadrant conditional f|q distribution
        Pf_qx0=np.where(Pqx0>0,np.exp(logPfqx0-np.log(Pqx0)),0)
        Pf_q0x=np.where(Pq0x>0,np.exp(logPfq0x-np.log(Pq0x)),0)
        Pf_qxx=np.where(Pqxx>0,np.exp(logPfqxx-np.log(Pqxx)),0)
    
        #3-quadrant q marginal distribution
        newPqZ=Pqx0 + Pq0x + Pqxx
        Pqx0/=newPqZ
        Pq0x/=newPqZ
        Pqxx/=newPqZ

        Pfqx0=np.exp(logPfqx0)
        Pfq0x=np.exp(logPfq0x)
        Pfqxx=np.exp(logPfqxx)
    
        print('Model probs: '+str(Pqx0)+' '+str(Pq0x)+' '+str(Pqxx))

        #get samples 
        num_samples=Nsamp
        q_samples=np.random.choice(range(3), num_samples, p=(Pqx0,Pq0x,Pqxx))
        vals,counts=np.unique(q_samples,return_counts=True)
        num_qx0=counts[0]
        num_q0x=counts[1]
        num_qxx=counts[2]
        print('q samples: '+str(sum(counts))+' '+str(num_qx0)+' '+str(num_q0x)+' '+str(num_qxx))
        print('q sampled probs: '+str(num_qx0/float(sum(counts)))+' '+str(num_q0x/float(sum(counts)))+' '+str(num_qxx/float(sum(counts))))
    
        #x0
        integ=np.exp(np.log(Pf_qx0)+logfvec)
        f_samples_inds= self.get_distsample(dlogf*(integ[1:] + integ[:-1]),num_qx0).flatten()
        f_sorted_inds=np.argsort(f_samples_inds)
        f_samples_inds=f_samples_inds[f_sorted_inds] 
        qx0_f_samples=fvec[f_samples_inds]
        find_vals,f_start_ind,f_counts=np.unique(f_samples_inds,return_counts=True,return_index=True)
        qx0_samples=np.zeros((num_qx0,))
        if noise_model<1:
            qx0_m_samples=np.zeros((num_qx0,))
            #conditioning on n>0 applies an m-dependent factor to Pm_f, which can't be incorporated into the ppf method used for noise_model 1 and 2. 
            #We handle that here by using a custom finite range sampler, which has the drawback of having to define an upper limit. 
            #This works so long as n_max/r_c<<m_max, so depends on highest counts in data (n_max). My data had max counts of 1e3-1e4.
            #Alternatively, could define a custom scipy RV class by defining it's PMF, but has to be array-compatible which requires care. 
            m_samp_max=int(1e5) 
            mvec=np.arange(m_samp_max)   
    
        for it,find in enumerate(find_vals):
            if noise_model==0:      
                m1=m_total*fvec[find]
                v1=m1+beta_mv*np.power(m1,alpha_mv)
                p=1-m1/v1
                n=m1*m1/v1/p
                Pm1_f=nbinom(n,1-p)
            
                Pm1_f_adj=np.exp(np.log(1-np.exp(-r_c1*mvec))+np.log(Pm1_f.pmf(mvec))-np.log((1-np.power(np.exp(r_c1+np.log(1-p))/(np.exp(r_c1)-p),n)))) #adds m-dependent factor due to conditioning on n>0...
                Pm1_f_adj_obj=rv_discrete(name='nbinom_adj',values=(mvec,Pm1_f_adj/np.sum(Pm1_f_adj)))
                qx0_m_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pm1_f_adj_obj.rvs(size=f_counts[it])
            
                mvals,minds,m_counts=np.unique(qx0_m_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]],return_inverse=True,return_counts=True)
                for mit,m in enumerate(mvals):
                    Pn1_m1=poisson(r_c1*m)
                    samples=np.random.random(size=m_counts[mit]) * (1-Pn1_m1.cdf(0)) + Pn1_m1.cdf(0)
                    qx0_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]][minds==mit]=Pn1_m1.ppf(samples)
 
        
            elif noise_model>0:
                samples=np.random.random(size=f_counts[it]) * (1-Pn1_f[find].cdf(0)) + Pn1_f[find].cdf(0)
                qx0_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pn1_f[find].ppf(samples)
            else:
                print('acq_model is 0,1, or 2 only')
        qx0_pair_samples=np.hstack((qx0_samples[:,np.newaxis],np.zeros((num_qx0,1)))) 
    
        #0x
        integ=np.exp(np.log(Pf_q0x)+logfvec)
        f_samples_inds=self.get_distsample(dlogf*(integ[1:] + integ[:-1]),num_q0x).flatten()
        f_sorted_inds=np.argsort(f_samples_inds)
        f_samples_inds=f_samples_inds[f_sorted_inds] 
        q0x_f_samples=fvec[f_samples_inds]
        find_vals,f_start_ind,f_counts=np.unique(f_samples_inds,return_counts=True,return_index=True)
        q0x_samples=np.zeros((num_q0x,))
        if noise_model<1:
            q0x_m_samples=np.zeros((num_q0x,))
        for it,find in enumerate(find_vals):
            if noise_model==0:
                m2=m_total*fvec[find]
                v2=m2+beta_mv*np.power(m2,alpha_mv)
                p=1-m2/v2
                n=m2*m2/v2/p
                Pm2_f=nbinom(n,1-p)
            
                Pm2_f_adj=np.exp(np.log(1-np.exp(-r_c2*mvec))+np.log(Pm2_f.pmf(mvec))-np.log((1-np.power(np.exp(r_c2+np.log(1-p))/(np.exp(r_c2)-p),n)))) #adds m-dependent factor due to conditioning on n>0...
                Pm2_f_adj_obj=rv_discrete(name='nbinom_adj',values=(mvec,Pm2_f_adj/np.sum(Pm2_f_adj)))
                q0x_m_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pm2_f_adj_obj.rvs(size=f_counts[it])

                mvals,minds,m_counts=np.unique(q0x_m_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]],return_inverse=True,return_counts=True)
                for mit,m in enumerate(mvals):
                    Pn2_m2=poisson(r_c2*m)
                    samples=np.random.random(size=m_counts[mit]) * (1-Pn2_m2.cdf(0)) + Pn2_m2.cdf(0)
                    q0x_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]][minds==mit]=Pn2_m2.ppf(samples)
        
                       
        
            elif noise_model > 0:
                samples=np.random.random(size=f_counts[it]) * (1-Pn2_f[find].cdf(0)) + Pn2_f[find].cdf(0)
                q0x_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pn2_f[find].ppf(samples)
            else:
                print('acq_model is 0,1,or 2 only')
        q0x_pair_samples=np.hstack((np.zeros((num_q0x,1)),q0x_samples[:,np.newaxis]))
    
        #qxx
        integ=np.exp(np.log(Pf_qxx)+logfvec)
        f_samples_inds=self.get_distsample(dlogf*(integ[1:] + integ[:-1]),num_qxx).flatten()        
        f_sorted_inds=np.argsort(f_samples_inds)
        f_samples_inds=f_samples_inds[f_sorted_inds] 
        qxx_f_samples=fvec[f_samples_inds]
        find_vals,f_start_ind,f_counts=np.unique(f_samples_inds,return_counts=True,return_index=True)
        qxx_n1_samples=np.zeros((num_qxx,))
        qxx_n2_samples=np.zeros((num_qxx,))
        if noise_model<1:
            qxx_m1_samples=np.zeros((num_qxx,))
            qxx_m2_samples=np.zeros((num_qxx,))
        for it,find in enumerate(find_vals):
            if noise_model==0:
                m1=m_total*fvec[find]
                v1=m1+beta_mv*np.power(m1,alpha_mv)
                p=1-m1/v1
                n=m1*m1/v1/p
                Pm1_f=nbinom(n,1-p)
            
                Pm1_f_adj=np.exp(np.log(1-np.exp(-r_c1*mvec))+np.log(Pm1_f.pmf(mvec))-np.log((1-np.power(np.exp(r_c1+np.log(1-p))/(np.exp(r_c1)-p),n)))) #adds m-dependent factor due to conditioning on n>0...
                if np.sum(Pm1_f_adj)==0:
                    qxx_m1_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=1
                else:
                    Pm1_f_adj_obj=rv_discrete(name='nbinom_adj',values=(mvec,Pm1_f_adj/np.sum(Pm1_f_adj)))
                    qxx_m1_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pm1_f_adj_obj.rvs(size=f_counts[it])

                mvals,minds,m_counts=np.unique(qxx_m1_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]],return_inverse=True,return_counts=True)
                for mit,m in enumerate(mvals):
                    Pn1_m1=poisson(r_c1*m)
                    samples=np.random.random(size=m_counts[mit]) * (1-Pn1_m1.cdf(0)) + Pn1_m1.cdf(0)
                    qxx_n1_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]][minds==mit]=Pn1_m1.ppf(samples)
                
                m2=m_total*fvec[find]
                v2=m2+beta_mv*np.power(m2,alpha_mv)
                p=1-m2/v2
                n=m2*m2/v2/p
                Pm2_f=nbinom(n,1-p)
            
                Pm2_f_adj=np.exp(np.log(1-np.exp(-r_c2*mvec))+np.log(Pm2_f.pmf(mvec))-np.log((1-np.power(np.exp(r_c2+np.log(1-p))/(np.exp(r_c2)-p),n)))) #adds m-dependent factor due to conditioning on n>0...
                if np.sum(Pm1_f_adj)==0:
                    qxx_m2_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=1
                else:
                    Pm2_f_adj_obj=rv_discrete(name='nbinom_adj',values=(mvec,Pm2_f_adj/np.sum(Pm2_f_adj)))
                    qxx_m2_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pm2_f_adj_obj.rvs(size=f_counts[it])

                mvals,minds,m_counts=np.unique(qxx_m2_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]],return_inverse=True,return_counts=True)
                for mit,m in enumerate(mvals):
                    Pn2_m2=poisson(r_c2*m)
                    samples=np.random.random(size=m_counts[mit]) * (1-Pn2_m2.cdf(0)) + Pn2_m2.cdf(0)
                    qxx_n2_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]][minds==mit]=Pn2_m2.ppf(samples)    

                          
            elif noise_model>0:
                samples=np.random.random(size=f_counts[it]) * (1-Pn1_f[find].cdf(0)) + Pn1_f[find].cdf(0)
                qxx_n1_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pn1_f[find].ppf(samples)
                samples=np.random.random(size=f_counts[it]) * (1-Pn2_f[find].cdf(0)) + Pn2_f[find].cdf(0)
                qxx_n2_samples[f_start_ind[it]:f_start_ind[it]+f_counts[it]]=Pn2_f[find].ppf(samples)
            else:
                print('acq_model is 0,1, or 2 only')
            
        qxx_pair_samples=np.hstack((qxx_n1_samples[:,np.newaxis],qxx_n2_samples[:,np.newaxis]))
    
        pair_samples=np.vstack((q0x_pair_samples,qx0_pair_samples,qxx_pair_samples))
        f_samples=np.concatenate((q0x_f_samples,qx0_f_samples,qxx_f_samples))
        output_m_samples=False
        if noise_model<1 and output_m_samples:                
            m1_samples=np.concatenate((q0x_m1_samples,qx0_m1_samples,qxx_m1_samples))
            m2_samples=np.concatenate((q0x_m2_samples,qx0_m2_samples,qxx_m2_samples))
    
        pair_samples_df=pd.DataFrame({'Clone_count_1':pair_samples[:,0],'Clone_count_2':pair_samples[:,1]})

        pair_samples_df['Clone_fraction_1'] = pair_samples_df['Clone_count_1']/np.sum(pair_samples_df['Clone_count_1'])
        pair_samples_df['Clone_fraction_2'] = pair_samples_df['Clone_count_2']/np.sum(pair_samples_df['Clone_count_2'])
    
        return f_samples,pair_samples_df
