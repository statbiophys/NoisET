# Import of the functional libraries
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
from scipy.optimize import minimize
from scipy.stats import nbinom

#from tqdm.notebook import trange, tqdm

#-----------------------------Transform-the-data------------------------------------

# ask the user the path to go and find the data
# transform the initial dataframes into dataframes that are usable in all the code

#filename = ask the user check how to do online

def change_data(filename,Patient, n, time, countsname, frequenciesname, nucleotidesname, aminoacidname, sequencestatusstate ):
    """ 
    Funtion that creates a neww data frame that is usable by all the codes that follow.

    Input : 
    filename :
    Patient : 
    n : 
    time :
    countsname : 
    frequenciesname :
    nucleotidesname :
    aminoacidname :
    sequencestatusstate : 

    Output :
    new_df : 



    """

    df = pd.read_table(filename, Patient, n, time )
    if sequencestatusstate : 
        # for Harlan Robins data-set, if not don't forget to extract only the "in-frame" data
        df = df.where(df['sequenceStatus'] == 'In')
    else :
        pass

    df = df[[countsname, frequenciesname, nucleotidesname, aminoacidname]]
    df = df.dropna()
    df = df.rename(columns = {nucleotidesname: 'N. Seq. CDR3', aminoacidname: 'AA. Seq. CDR3', countsname: 'Clone count', frequenciesname:'Clone fraction' })
    df.to_csv(Patient + '_'+ time + '_F' + n +'_.txt', index=None, sep = '\t')

    new_df = pd.read_csv(Patient + '_'+ time + '_F' + n +'_.txt', sep = '\t')

    return new_df

#===============================Data-Pre-Processing===================================

class Data_Process:

    """Explain this class of methods """

    def import_data(self,path,filename1,filename2,mincount,maxcount,colnames1,colnames2):
        '''
        Reads in Yellow fever data from two datasets and merges based on nt sequence.
        Outputs dataframe of pair counts for all clones.
        Considers clones with counts between mincount and maxcount
        Uses specified column names and headerline in stored fasta file.
        '''
    
        headerline=0 #line number of headerline
        newnames=['Clone_fraction','Clone_count','ntCDR3','AACDR3']    
        with open(path+filename1, 'r') as f:
            F1Frame_chunk=pd.read_csv(f,delimiter='\t',usecols=colnames1,header=headerline)[colnames1]
        with open(path+filename2, 'r') as f:
            F2Frame_chunk=pd.read_csv(f,delimiter='\t',usecols=colnames2,header=headerline)[colnames2]
        F1Frame_chunk.columns=newnames
        F2Frame_chunk.columns=newnames
        suffixes=('_1','_2')
        mergedFrame=pd.merge(F1Frame_chunk,F2Frame_chunk,on=newnames[2],suffixes=suffixes,how='outer')
        for nameit in [0,1]:
            for labelit in suffixes:
                mergedFrame.loc[:,newnames[nameit]+labelit].fillna(int(0),inplace=True)
                if nameit==1:
                    mergedFrame.loc[:,newnames[nameit]+labelit].astype(int)
        def dummy(x):
            val=x[0]
            if pd.isnull(val):
                val=x[1]    
            return val
        mergedFrame.loc[:,newnames[3]+suffixes[0]]=mergedFrame.loc[:,[newnames[3]+suffixes[0],newnames[3]+suffixes[1]]].apply(dummy,axis=1) #assigns AA sequence to clones, creates duplicates
        mergedFrame.drop(newnames[3]+suffixes[1], 1,inplace=True) #removes duplicates
        mergedFrame.rename(columns = {newnames[3]+suffixes[0]:newnames[3]}, inplace = True)
        mergedFrame=mergedFrame[[newname+suffix for newname in newnames[:2] for suffix in suffixes]+[newnames[2],newnames[3]]]
        filterout=((mergedFrame.Clone_count_1<mincount) & (mergedFrame.Clone_count_2==0)) | ((mergedFrame.Clone_count_2<mincount) & (mergedFrame.Clone_count_1==0)) #has effect only if mincount>0
        number_clones=len(mergedFrame)
        return number_clones,mergedFrame.loc[((mergedFrame.Clone_count_1<=maxcount) & (mergedFrame.Clone_count_2<=maxcount)) & ~filterout]
    
    def get_sparserep(self,counts):


        '''
        Tranforms {(n1,n2)} data stored in pandas dataframe to a sparse 1D representation.
        unicountvals_1(2) are the unique values of n1(2).
        sparse_rep_counts gives the counts of unique pairs.
        indn1(2) is the index of unicountvals_1(2) giving the value of n1(2) in that unique pair.
        len(indn1)=len(indn2)=len(sparse_rep_counts)
        '''


        counts['paircount'] = 1  # gives a weight of 1 to each observed clone

        clone_counts = counts.groupby(['Clone_count_1', 'Clone_count_2']).sum()
        sparse_rep_counts = np.asarray(clone_counts.values.flatten(), dtype=int)
        clonecountpair_vals = clone_counts.index.values
        indn1 = np.asarray([clonecountpair_vals[it][0] for it in range(len(sparse_rep_counts))], dtype=int)
        indn2 = np.asarray([clonecountpair_vals[it][1] for it in range(len(sparse_rep_counts))], dtype=int)
        NreadsI = np.sum(counts['Clone_count_1'])
        NreadsII = np.sum(counts['Clone_count_2'])

        unicountvals_1, indn1 = np.unique(indn1, return_inverse=True)
        unicountvals_2, indn2 = np.unique(indn2, return_inverse=True)

        return indn1, indn2, sparse_rep_counts, unicountvals_1, unicountvals_2, NreadsI, NreadsII
    
    
    def filter_df(self, n_limit_1, n_limit_2, flag, path,filename1,filename2,mincount,maxcount,colnames1,colnames2):
        
        """Function that filters the initial data ensemble
        Explain the use of the flag, and that df is a global variable / to be improved """
        
        if flag :
            n, df = self.import_data(path,filename1,filename2,mincount,maxcount,colnames1,colnames2)
        else :
            pass 
        df_bis = df[(df['Clone_count_1'] > n_limit_1) & (df['Clone_count_1'] < n_limit_2) & (df['Clone_count_2'] != 0) ]
        
        return df_bis
    
    def scatter_plot(self, n_limit_1, n_limit_2, time_unit, time_1, time_2, Patient, path,filename1,filename2,mincount,maxcount,colnames1,colnames2): 
        
        """Function that plots the scatter-plot of the wanted initial data frame 
        
        n_limit_1 : 
        n_limit_2 : 
        time_unit : 
        time_1 :
        time_2 :
        Patient : """
        
        
        
        X = np.arange(1e5)
        
        n, df = self.import_data(path,filename1,filename2,mincount,maxcount,colnames1,colnames2)
        df_bis = df[(df['Clone_count_1'] > n_limit_1) & (df['Clone_count_1'] < n_limit_2) & (df['Clone_count_2'] != 0) ]

        plt.figure(figsize=(10,8))
        plt.style.use("seaborn-whitegrid")
        plt.scatter(df_bis['Clone_fraction_1'], df_bis['Clone_fraction_2'], c='none', edgecolor='DarkBlue')
        plt.plot(X,X,c='tomato')
        plt.xlabel(r'$clone \ frequency \ ' + time_unit + '_{' +  time_1 +  '}$', fontsize=20 )
        plt.ylabel(r'$clone \ frequency \ ' + time_unit + '_{' +  time_2 +  '}$', fontsize=20  )
        plt.xscale('log')
        plt.yscale('log')
        plt.axis([8e-7, 2e-2,  8e-7, 2e-2])
        plt.title(r'$' + Patient + '\ between \  ' + time_unit + '_{' +  time_1 +  '} ' + ' \ and \ ' + time_unit + '_{' +  time_2 +  '}' + ' \ ' + str(n_limit_1) + ' < n_1 < ' + str(n_limit_2) + '$', fontsize=20)
        plt.show()

# Put the data characteristics here also 
# for all the filenames construct a Jaccard Index (for the Jaccard MAtrix)
# In the future link the responding clones (Alice - read - precise tracking vaccine in YFV // revaccination )
        

#===============================Noise-Model=========================================

class RNASeq_NoiseModel:

    """ Explain this class of methods"""

    def NegBinPar(self,m,v,mvec): 
        '''
        Same as NegBinParMtr, but for m and v being scalars.
        Assumes m>0.
        Output is (len(mvec),) array
        '''
        mmax=mvec[-1]
        p = 1-m/v
        r = m*m/v/p
        NBvec=np.arange(mmax+1,dtype=float)   
        NBvec[1:]=np.log((NBvec[1:]+r-1)/NBvec[1:]*p) #vectorization won't help unfortuneately here since log needs to be over array
        NBvec[0]=r*math.log(m/v)
        NBvec=np.exp(np.cumsum(NBvec)[mvec]) #save a bit here
        return NBvec

    def NegBinParMtr(self,m,v,nvec): #speed up only insofar as the log and exp are called once on array instead of multiple times on rows
        ''' 
        computes NegBin probabilities over the ordered (but possibly discontiguous) vector (nvec) 
        for mean/variance combinations given by the mean (m) and variance (v) vectors. 
        Note that m<v for negative binomial.
        Output is (len(m),len(nvec)) array
        '''
        nmax=nvec[-1]
        p = 1-m/v
        r = m*m/v/p
        NBvec=np.arange(nmax+1,dtype=float)
        NBvec=np.log((NBvec+r[:,np.newaxis]-1)*(p[:,np.newaxis]/NBvec))
        NBvec[:,0]=r*np.log(m/v) #handle NBvec[0]=0, treated specially when m[0]=0, see below
        NBvec=np.exp(np.cumsum(NBvec,axis=1)) #save a bit here
        if m[0]==0:
            NBvec[0,:]=0.
            NBvec[0,0]=1.
        NBvec=NBvec[:,nvec]
        return NBvec

    def PoisPar(self, Mvec,unicountvals):
        #assert Mvec[0]==0, "first element needs to be zero"
        nmax=unicountvals[-1]
        nlen=len(unicountvals)
        mlen=len(Mvec)
        Nvec=unicountvals
        logNvec=-np.insert(np.cumsum(np.log(np.arange(1,nmax+1))),0,0.)[unicountvals] #avoid n=0 nans  
        Nmtr=np.exp(Nvec[np.newaxis,:]*np.log(Mvec)[:,np.newaxis]+logNvec[np.newaxis,:]-Mvec[:,np.newaxis]) # np.log(Mvec) throws warning: since log(0)=-inf
        if Mvec[0]==0:
            Nmtr[0,:]=np.zeros((nlen,)) #when m=0, n=0, and so get rid of nans from log(0)
            Nmtr[0,0]=1. #handled belowacq_model_type
        if unicountvals[0]==0: #if n=0 included get rid of nans from log(0)
            Nmtr[:,0]=np.exp(-Mvec)
        return Nmtr

    def get_rhof(self,alpha_rho, nfbins,fmin,freq_dtype):
        '''
        generates power law (power is alpha_rho) clone frequency distribution over 
        freq_nbins discrete logarithmically spaced frequences between fmin and 1 of dtype freq_dtype
        Outputs log probabilities obtained at log frequencies'''
        fmax=1e0
        logfvec=np.linspace(np.log10(fmin),np.log10(fmax), nfbins)
        logfvec=np.array(np.log(np.power(10,logfvec)) ,dtype=freq_dtype).flatten()  
        logrhovec=logfvec*alpha_rho
        integ=np.exp(logrhovec+logfvec,dtype=freq_dtype)
        normconst=np.log(np.dot(np.diff(logfvec)/2.,integ[1:]+integ[:-1]))
        logrhovec-=normconst 
        return logrhovec,logfvec, normconst


    def get_logPn_f(self,unicounts,Nreads,logfvec,model,paras):

        """"""


         # Choice of the model:

        if model == "NB+Poiss":
            acq_model_type = 0
        elif model == "Poiss+NB":
            acq_model_type = 1
        elif model == "NBonly":
            acq_model_type = 2
        elif model =="Poissonly":
            acq_model_type = 3
    
        if acq_model_type<2:

            m_total=float(np.power(10, paras[3])) 
            r_c=Nreads/m_total
        if acq_model_type<3:

            beta_mv= paras[1]
            alpha_mv=paras[2]
        
        if acq_model_type<2: #for models that include cell counts
            #compute parametrized range (mean-sigma,mean+5*sigma) of m values (number of cells) conditioned on n values (reads) appearing in the data only 
            nsigma=5.
            nmin=300.
            #for each n, get actual range of m to compute around n-dependent mean m
            m_low =np.zeros((len(unicounts),),dtype=int)
            m_high=np.zeros((len(unicounts),),dtype=int)
            for nit,n in enumerate(unicounts):
                mean_m=n/r_c
                dev=nsigma*np.sqrt(mean_m)
                m_low[nit] =int(mean_m-  dev) if (mean_m>dev**2) else 0                         
                m_high[nit]=int(mean_m+5*dev) if (      n>nmin) else int(10*nmin/r_c)
            m_cellmax=np.max(m_high)
            #across n, collect all in-range m
            mvec_bool=np.zeros((m_cellmax+1,),dtype=bool) #cheap bool
            nvec=range(len(unicounts))
            for nit in nvec:
                mvec_bool[m_low[nit]:m_high[nit]+1]=True  #mask vector
            mvec=np.arange(m_cellmax+1)[mvec_bool]                
            #transform to in-range index
            for nit in nvec:
                m_low[nit]=np.where(m_low[nit]==mvec)[0][0]
                m_high[nit]=np.where(m_high[nit]==mvec)[0][0]

        Pn_f=np.zeros((len(logfvec),len(unicounts)))
        if acq_model_type==0:

            mean_m=m_total*np.exp(logfvec)
            var_m=mean_m+beta_mv*np.power(mean_m,alpha_mv)
            Poisvec = self.PoisPar(mvec*r_c,unicounts)
            for f_it in range(len(logfvec)):
                NBvec=self.NegBinPar(mean_m[f_it],var_m[f_it],mvec)
                for n_it,n in enumerate(unicounts):
                    Pn_f[f_it,n_it]=np.dot(NBvec[m_low[n_it]:m_high[n_it]+1],Poisvec[m_low[n_it]:m_high[n_it]+1,n_it]) 
        elif acq_model_type==1:

            Poisvec= self.PoisPar(m_total*np.exp(logfvec),mvec)
            mean_n=r_c*mvec
            NBmtr= self.NegBinParMtr(mean_n,mean_n+beta_mv*np.power(mean_m,alpha_mv),unicounts)
            for f_it in range(len(logfvec)):
                Poisvectmp=Poisvec[f_it,:]
                for n_it,n in enumerate(unicounts):
                    Pn_f[f_it,n_it]=np.dot(Poisvectmp[m_low[n_it]:m_high[n_it]+1],NBmtr[m_low[n_it]:m_high[n_it]+1,n_it]) 
        elif acq_model_type==2:

            mean_n=Nreads*np.exp(logfvec)
            var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
            Pn_f = self.NegBinParMtr(mean_n,var_n,unicounts)
        elif acq_model_type==3:

            mean_n=Nreads*np.exp(logfvec)
            Pn_f= self.PoisPar(mean_n,unicounts)
        else:
            print('acq_model is 0,1,2, or 3 only')

        return np.log(Pn_f)

    #-----------------------------Null-Model-optimization--------------------------
    
    def get_Pn1n2_notfiltered(self,paras, sparse_rep, model):

        """
        For the whole frequencies spectrum !! 
        """
        # Choice of the model:

        if model == "NB+Poiss":
            acq_model_type = 0
        elif model == "Poiss+NB":
            acq_model_type = 1
        elif model == "NBonly":
            acq_model_type = 2
        elif model =="Poissonly":
            acq_model_type = 3

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,NreadsI,NreadsII = sparse_rep
        
        nfbins = 1200
        freq_dtype = float

        # Parameters

        alpha = paras[0]
        fmin = np.power(10,paras[3])

        # 
        logrhofvec, logfvec, normconst = self.get_rhof(alpha,nfbins,fmin,freq_dtype)

        # 

        logfvec_tmp=deepcopy(logfvec)

        logPn1_f = self.get_logPn_f(unicountvals_1, NreadsI,logfvec_tmp,acq_model_type,paras)
        logPn2_f = self.get_logPn_f(unicountvals_2, NreadsII,logfvec_tmp,acq_model_type,paras)

        # for the trapezoid integral methods

        dlogfby2=np.diff(logfvec)/2
        #print("computing P(n1,n2)")
        Pn1n2 = np.zeros((len(sparse_rep_counts)))
        for it, (n1_it, n2_it) in enumerate(zip(indn1, indn2)) :
            integ = np.exp(logrhofvec+logPn2_f[:,n2_it]+logPn1_f[:,n1_it]+ logfvec )
            Pn1n2[it] = np.dot(dlogfby2,integ[1:] + integ[:-1])
         
        #Compute P(0,0) for the normalization constraint
        integ = np.exp(logrhofvec+logPn2_f[:,0]+logPn1_f[:,0]+ logfvec)
        Pn0n0 = np.dot(dlogfby2,integ[1:] + integ[:-1])

        N_samp = np.sum(sparse_rep_counts)
        C = N_samp*np.log(1-Pn0n0)

        return - np.dot(sparse_rep_counts, np.log(Pn1n2)) + C


    def callback(self,paras, nparas):
        '''prints iteration info. called by scipy.minimize'''

        global curr_iter
        print(''.join(['{0:d} ']+['{'+str(it)+':3.6f} ' for it in range(1,len(paras)+1)]).format(*([curr_iter]+list(paras))))
        curr_iter += 1


    # Constraints for the Null-Model, no filtered 
    def nullmodel_constr_fn(self,paras, sparse_rep, model, constr_type):
        
        '''
        returns either or both of the two level-set functions: log<f>-log(1/N), with N=Nclones/(1-P(0,0)) and log(Z_f), with Z_f=N<f>_{n+n'=0} + sum_i^Nclones <f>_{f|n,n'}
        '''

         # Choice of the model:

        if model == "NB+Poiss":
            acq_model_type = 0
        elif model == "Poiss+NB":
            acq_model_type = 1
        elif model == "NBonly":
            acq_model_type = 2
        elif model =="Poissonly":
            acq_model_type = 3

        indn1, indn2, sparse_rep_counts, unicountvals_1, unicountvals_2, NreadsI, NreadsII = sparse_rep

        #Variables that would be chosen in the future by the user 
        nfbins = 1200
        freq_dtype = float

        alpha = paras[0]  # power law exponent
        fmin = np.power(10, paras[-1]) # true minimal frequency 

        logrhofvec, logfvec, normconst = self.get_rhof(alpha,nfbins,fmin,freq_dtype)
        dlogfby2 = np.diff(logfvec) / 2.  # 1/2 comes from trapezoid integration below

        integ = np.exp(logrhofvec + 2 * logfvec)
        avgf_ps = np.dot(dlogfby2, integ[:-1] + integ[1:])

        logPn1_f = self.get_logPn_f(unicountvals_1, NreadsI, logfvec, acq_model_type, paras)
        logPn2_f = self.get_logPn_f(unicountvals_2, NreadsII, logfvec, acq_model_type, paras)

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


    
        # Null-Model optimization learning 

    def learn_null_model(self, sparse_rep, model, init_paras, constr_type=1):  # constraint type 1 gives only low error modes, see paper for details.
        '''
        performs constrained maximization of null model likelihood
        '''
        

         # Choice of the model:

        if model == "NB+Poiss":
            acq_model_type = 0
        elif model == "Poiss+NB":
            acq_model_type = 1
        elif model == "NBonly":
            acq_model_type = 2
        elif model =="Poissonly":
            acq_model_type = 3


        # Parameters initialization depending on the model 
        if acq_model_type < 2:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'm_total', 'fmin']
        elif acq_model_type == 2:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'fmin']
        else:
            parameter_labels = ['alph_rho', 'fmin']

        assert len(parameter_labels) == len(init_paras), "number of model and initial paras differ!"

        condict = {'type': 'eq', 'fun': self.nullmodel_constr_fn, 'args': (sparse_rep, acq_model_type, constr_type)}


        partialobjfunc = partial(self.get_Pn1n2_notfiltered, sparse_rep=sparse_rep, acq_model_type=acq_model_type)
        nullfunctol = 1e-6
        nullmaxiter = 200
        header = ['Iter'] + parameter_labels
        print(''.join(['{' + str(it) + ':9s} ' for it in range(len(init_paras) + 1)]).format(*header))
        
        global curr_iter
        curr_iter = 1
        callbackp = partial(self.callback, nparas=len(init_paras))
        outstruct = minimize(partialobjfunc, init_paras, method='SLSQP', callback=callbackp, constraints=condict,
                         options={'ftol': nullfunctol, 'disp': True, 'maxiter': nullmaxiter})
        
        constr_value = self.nullmodel_constr_fn(outstruct.x, sparse_rep, acq_model_type, constr_type)

        print(outstruct)
        return outstruct, constr_value


#=========================DIFFEXPR=================================

class Propagator:

    """
    Definition of all the propagators 
    """


