  #from tqdm.notebook import trange, tqdm

    #-----------------------------Transform-the-data------------------------------------

    # ask the user the path to go and find the data
    # transform the initial dataframes into dataframes that are usable in all the code

    #filename = ask the user check how to do online

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


#-----------------------------Transform-the-data------------------------------------

# ask the user the path to go and find the data
# transform the initial dataframes into dataframes that are usable in the entire code

#filename = ask the user check how to do online

def change_data(filename, Patient, n, time, countsname, frequenciesname, nucleotidesname, aminoacidname, sequencestatusstate):
    """ 
    Funtion that creates a new data frame that is usable by all the codes that follow.

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

class Data_Process(object):

    """Explain this class of methods 
    path :
    filename1 :
    filename2 : 
    mincount : 
    maxcount : 
    colnames1 : 
    colnames2 : """

    def __init__(self, path, filename1, filename2, mincount, maxcount, colnames1,  colnames2):
        self.path = path
        self.filename1 = filename1
        self.filename2 = filename2
        self.mincount = mincount
        self.maxcount = maxcount
        self.colnames1 = colnames1
        self.colnames2 = colnames2
    

    def import_data(self):
        '''
        Reads in Yellow fever data from two datasets and merges based on nt sequence.
        Outputs dataframe of pair counts for all clones.
        Considers clones with counts between mincount and maxcount
        Uses specified column names and headerline in stored fasta file.
        
        number_clones:
        df:
        '''
        
        headerline=0 #line number of headerline
        newnames=['Clone_fraction','Clone_count','ntCDR3','AACDR3']    
        with open(self.path+ self.filename1, 'r') as f:
            F1Frame_chunk=pd.read_csv(f,delimiter='\t',usecols=self.colnames1,header=headerline)[self.colnames1]
        with open(self.path+self.filename2, 'r') as f:
            F2Frame_chunk=pd.read_csv(f,delimiter='\t',usecols=self.colnames2,header=headerline)[self.colnames2]
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
        filterout=((mergedFrame.Clone_count_1<self.mincount) & (mergedFrame.Clone_count_2==0)) | ((mergedFrame.Clone_count_2<self.mincount) & (mergedFrame.Clone_count_1==0)) #has effect only if mincount>0
        number_clones=len(mergedFrame)
        return number_clones,mergedFrame.loc[((mergedFrame.Clone_count_1<=self.maxcount) & (mergedFrame.Clone_count_2<=self.maxcount)) & ~filterout]
        

           
    def filter_df(self, n_limit_1, n_limit_2, flag):
            
        """Function that filters the initial data ensemble
        Explain the use of the flag, and that df is a global variable / to be improved 
        n_limit_1 :
        n_limit_2:
        flag:
        
        df_bis"""
            
        if flag :
            n, df = self.import_data()
        else :
            pass 
        df_bis = df[(df['Clone_count_1'] > n_limit_1) & (df['Clone_count_1'] < n_limit_2) & (df['Clone_count_2'] != 0) ]
            
        return df_bis

    def get_sparserep(self, df): 
        """
        Tranforms {(n1,n2)} data stored in pandas dataframe to a sparse 1D representation.
        unicountvals_1(2) are the unique values of n1(2).
        sparse_rep_counts gives the counts of unique pairs.
        ndn1(2) is the index of unicountvals_1(2) giving the value of n1(2) in that unique pair.
        len(indn1)=len(indn2)=len(sparse_rep_counts)"""
        
        counts = df.loc[:,['Clone_count_1', 'Clone_count_2']]
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

    def scatter_plot(self, n_limit_1, n_limit_2, time_unit, time_1, time_2, Patient): 
            
        """Function that plots the scatter-plot of the wanted initial data frame 
            
        n_limit_1 : 
        n_limit_2 : 
        time_unit : 
        time_1 :
        time_2 :
        Patient : """
            
            
            
        X = np.arange(1e5)
            
        n, df = self.import_data()
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
    
    def frequencies_dist_inters(self, n_limit_1, n_limit_2, time_unit, time_1, time_2, Patient):

        """Be careful, it is not the true frequencies distributions of the two time points but the frequencies distributions of 
        the input data (persistent clones between two time points)
        note: if you want the true frequencies distributions, ...
        n_limit_1 : 
        n_limit_2 :
        time_unit : 
        time_1 :
        time_2 : 
        Patient : 
        """

        n, df = self.import_data()
        df_bis = df[(df['Clone_count_1'] > n_limit_1) & (df['Clone_count_1'] < n_limit_2) & (df['Clone_count_2'] != 0) ]

        f_detect1 = np.min(df_bis['Clone_fraction_1'])
        f_detect2 = np.min(df_bis['Clone_fraction_2'])

        plt.figure(figsize=(10,10))

        plt.hist(df_bis['Clone_fraction_1'],bins = 1000, log = True, color= 'skyblue', label = r'$Data \ ' + time_unit + '_{' +  time_1 +  '}'  + '$')
        plt.axvline(f_detect1, color = 'DarkBlue', linestyle = 'dashed', linewidth = 6, label = r'$fdetect'+ '_{' +  time_1 +  '} = '+ str(round(f_detect1,8)) + '$' )
        plt.hist(df_bis['Clone_fraction_2'],bins = 1000, log = True, color= 'salmon', label = r'$Data \ ' + time_unit + '_{' +  time_2 +  '}' + '$')
        plt.axvline(f_detect2, color = 'DarkRed', linestyle = 'dashed', linewidth = 6, label = r'$fdetect'+ '_{' +  time_2 +  '} = '+ str(round(f_detect2,8)) + '$'  )
        plt.xscale('log')
        plt.xlabel(r'$Frequencies$', fontsize = 20)
        plt.legend()
        plt.title(r'$Frequencies\ histogram \ for \ ' + Patient + '$', fontsize = 20)
        plt.show()

        x_1 = np.log(f_detect1)
        x_2 = np.log(f_detect2)

        return x_1, x_2

    def cumulative_plots(self, n_limit_1, n_limit_2, time_unit, time_1, time_2, Patient):

        """
        n_limit_1 : 
        n_limit_2 :
        time_unit : 
        time_1 :
        time_2 : 
        Patient : 

        """


        n, df = self.import_data()
        df_bis = df[(df['Clone_count_1'] > n_limit_1) & (df['Clone_count_1'] < n_limit_2) & (df['Clone_count_2'] != 0) ]


        def cumulative_frequencies_date1(df):
            df_sorted = df.sort_values(by = ['Clone_fraction_1'], ascending = True)
            n , c = df_sorted.shape
            X =     np.array(df_sorted['Clone_fraction_1'])
            Y = np.zeros((n,1))
            for i in tqdm(range(n)):
                Bool = X >= X[i]
                Y[i] = np.sum(Bool)
    
            return Y, X

        def cumulative_frequencies_date2(df):
            
            df_sorted = df.sort_values(by = ['Clone_fraction_2'], ascending = True)
            n , c = df_sorted.shape
            X = np.array(df_sorted['Clone_fraction_2'])
            Y = np.zeros((n,1))
            for i in tqdm(range(n)):
                Bool = X >= X[i]
                Y[i] = np.sum(Bool)
    
            return Y, X

        
        Y_1, X_1 = cumulative_frequencies_date1(df_bis)
        Y_2, X_2 = cumulative_frequencies_date2(df_bis)
        Y_1_bis = np.reshape(Y_1, len(Y_1))
        X_1_bis = np.reshape(X_1, len(X_1))

        Y_2_bis = np.reshape(Y_2, len(Y_2))
        X_2_bis = np.reshape(X_2, len(X_2))

        slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(np.log(X_1_bis), np.log(Y_1_bis) )
        y_1 = intercept_1 + slope_1*np.array(np.log(X_1_bis))

        slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(np.log(X_2_bis), np.log(Y_2_bis) )
        y_2 = intercept_2 + slope_2*np.array(np.log(X_2_bis))

        plt.figure(figsize = (12,10))
        
        plt.plot(X_1_bis , Y_1_bis, color='DarkBlue', linestyle='dashed', linewidth=2, markersize=12, label = r'$Cumulative \  frequencies \  distribution \ ' + Patient + ' \ '  + time_unit  +  time_1 + '$')
        plt.plot(X_2_bis , Y_2_bis, color='lightblue', linestyle='dashed', linewidth=2, markersize=12, label = r'$Cumulative \  frequencies \ distribution  \ ' + Patient + ' \ ' + time_unit  +  time_2 + '$')
        plt.plot(X_1_bis , np.exp(y_1), color='red', linestyle='solid', linewidth=2, markersize=12, label = r'$fit \ ( \alpha + 1)  = ' +  str(round(slope_1, 3)) + '$' )
        plt.plot(X_2_bis , np.exp(y_2), color='DarkRed', linestyle='solid', linewidth=2, markersize=12, label = r'$fit \ ( \alpha + 1)  = ' +  str(round(slope_2, 3)) + '$')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$frequencies$', fontsize=20 )
        plt.ylabel(r'$Cumulative \ distritution$', fontsize=20 )
        plt.title(r'$Cumulative \ frequencies \ distributions \ ' + Patient + '\ between \  ' + time_unit + '_{' +  time_1 +  '} ' + ' \ and \ ' + time_unit + '_{' +  time_2 +  '}' + ' \ ' + str(n_limit_1) + ' < n_1 < ' + str(n_limit_2) + '$', fontsize=20)
        plt.axis([4e-7, 1e-2, 0, 1e5])
        plt.legend()
        plt.show()

        return y_1






    
    

    # Put the data characteristics here also 
    # for all the filenames construct a Jaccard Index (for the Jaccard MAtrix)
    # In the future link the responding clones (Alice - read - precise tracking vaccine in YFV // revaccination )
            

#===============================Noise-Model=========================================

class RNAUMIs_NoiseModel:


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


    def get_logPn_f(self,unicounts,Nreads,logfvec, acq_model_type, paras):

        """"""


        # Choice of the model:
        
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
        
    def get_Pn1n2(self,paras, sparse_rep, acq_model_type):

        """
        
        """
        # Choice of the model:


        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,NreadsI,NreadsII = sparse_rep
            
        nfbins = 1200
        freq_dtype = float

        # Parameters

        alpha = paras[0]
        fmin = np.power(10,paras[-1])

        # 
        logrhofvec, logfvec, normconst = self.get_rhof(alpha,nfbins,fmin,freq_dtype)

        # 

        logfvec_tmp=deepcopy(logfvec)

        logPn1_f = self.get_logPn_f(unicountvals_1, NreadsI,logfvec_tmp,acq_model_type,paras)
        logPn2_f = self.get_logPn_f(unicountvals_2, NreadsII,logfvec_tmp,acq_model_type,paras)

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

    


    def callback(self,paras, nparas, sparse_rep, acq_model_type):
        '''prints iteration info. called by scipy.minimize'''

        global curr_iter
        #curr_iter = 0
        global Loss_function 
        print(''.join(['{0:d} ']+['{'+str(it)+':3.6f} ' for it in range(1,len(paras)+1)]).format(*([curr_iter]+list(paras))))
        #print ('{' + str(len(paras)+1) + ':3.6f}'.format( [self.get_Pn1n2(paras, sparse_rep, acq_model_type)]))
        Loss_function = self.get_Pn1n2(paras, sparse_rep, acq_model_type)
        print(Loss_function)
        curr_iter += 1
        


    # Constraints for the Null-Model, no filtered 
    def nullmodel_constr_fn(self,paras, sparse_rep, acq_model_type, constr_type):
            
        '''
        returns either or both of the two level-set functions: log<f>-log(1/N), with N=Nclones/(1-P(0,0)) and log(Z_f), with Z_f=N<f>_{n+n'=0} + sum_i^Nclones <f>_{f|n,n'}
        '''

    # Choice of the model: 

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

    def learn_null_model(self, sparse_rep, acq_model_type, init_paras, constr_type=1, display_loss_function = False):  # constraint type 1 gives only low error modes, see paper for details.
        '''
        performs constrained maximization of null model likelihood
        '''
            

        # Choice of the model:
        # Parameters initialization depending on the model 
        if acq_model_type < 2:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'm_total', 'fmin']
        elif acq_model_type == 2:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'fmin']
        else:
            parameter_labels = ['alph_rho', 'fmin']

        assert len(parameter_labels) == len(init_paras), "number of model and initial paras differ!"

        condict = {'type': 'eq', 'fun': self.nullmodel_constr_fn, 'args': (sparse_rep, acq_model_type, constr_type)}


        partialobjfunc = partial(self.get_Pn1n2, sparse_rep=sparse_rep, acq_model_type=acq_model_type)
        nullfunctol = 1e-6
        nullmaxiter = 200
        header = ['Iter'] + parameter_labels
        print(''.join(['{' + str(it) + ':9s} ' for it in range(len(init_paras) + 1)]).format(*header))
            
        global curr_iter
        curr_iter = 1
        callbackp = partial(self.callback, nparas=len(init_paras), sparse_rep = sparse_rep, acq_model_type= acq_model_type)
        outstruct = minimize(partialobjfunc, init_paras, method='SLSQP', callback=callbackp, constraints=condict,
                        options={'ftol': nullfunctol, 'disp': True, 'maxiter': nullmaxiter})
            
        constr_value = self.nullmodel_constr_fn(outstruct.x, sparse_rep, acq_model_type, constr_type)

        print(outstruct)

        if display_loss_function :

            plt.figure()

            plt.style.use('seaborn-whitegrid')

            n_iter = np.arange(len(Loss_function))
            plt.plot(n_iter, Loss_function, marker = 'o', markersize = 5, linestyle = '--', linewidth =1, color = 'DarkBlue')
            plt.xlabel(r'$Number \ of \ iterations$')
            plt.ylabel(r'$Loss \ function$')


        return outstruct, constr_value, Loss_function

#============================================Differential expression =============================================================

class Expansion_Model:
    
    """
    """

    
    
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
        return logrhovec,logfvec

    
    def get_logPn_f(self,unicounts,Nreads,logfvec, acq_model_type, paras):

        """"""


        # Choice of the model:
        
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

    def get_Ps(self, alp,sbar,smax,stp):
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

    def callbackFdiffexpr(self, Xi): #case dependent
        '''prints iteration info. called scipy.minimize'''
               
        print('{0: 3.6f}   {1: 3.6f}   '.format(Xi[0], Xi[1])+'\n')   
    

    def learning_dynamics_expansion_polished(self, paras_1, paras_2, sparse_rep, acq_model_type):
        """
        Different uses for this function to explain, explain the use of NreadsItrue and NreadsIItrue
        paras : 
        sparse_rep :
        acq_model_type : 
        not_filtered : True if all the repertoire is used to infer the dynamics paramers, false if data is filtered
        NreadsItrue : the total number of reads in the original sample at time 1
        NreadsII true: the total number of reads in the original sample at time 2
        time_unit : day, months or years
        time_1 : indication of the first sample extraction time
        time_2 : indication of the second sample extraction time
        """

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,NreadsI,NreadsII = sparse_rep

        alpha_rho = paras_1[0]
        fmin = np.power(10,paras_1[-1])
        freq_dtype = 'float64'
        nfbins = 1200 #Accuracy of the integration


        logrhofvec, logfvec = get_rhof(self, alpha_rho, nfbins, fmin, freq_dtype)

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
                logPn1_f = get_logPn_f(self, unicounts, Nreads, logfvec_tmp, acq_model_type, paras)

            else:
                logPn2_f = get_logPn_f(self, unicounts, Nreads, logfvec_tmp, acq_model_type, paras)

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
        print("N_obs: " + str(N_obs))
    
            
        def cost(PARAS):

            alp = PARAS[0]
            sbar = PARAS[1]

            Ps = get_Ps(self,alp,sbar,smax,s_step)
            Pn0n0=np.dot(Pn0n0_s,Ps)
            Pn1n2_ps=np.sum(Pn1n2_s*Ps[:,np.newaxis,np.newaxis],0)
            Pn1n2_ps/=1-Pn0n0
            print(Pn0n0)

       

            Energy = - np.dot(sparse_rep_counts/float(N_obs),np.where(Pn1n2_ps[indn1,indn2]>0,np.log(Pn1n2_ps[indn1,indn2]),0)) 
                
            return Energy

    #--------------------------Compute-the-grid-----------------------------------------
        
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
        print("--- %s seconds ---" % (time.time() - st))
    
    
    #------------------------------Optimization----------------------------------------------
        
        optA = alpmesh[a[0],b[0]]
        optB = sbarmesh[a[0],b[0]]
                  
        print('polish parameter estimate from '+ str(optA)+' '+str(optB))
        initparas=(optA,optB)  
    

        outstruct = minimize(cost, initparas, method='SLSQP', callback=callbackFdiffexpr, tol=1e-6,options={'ftol':1e-8 ,'disp': True,'maxiter':300})

        return outstruct.x, Pn1n2_s, Pn0n0_s, svec

    def learning_dynamics_expansion(self, paras_1, paras_2, sparse_rep, acq_model_type, display_plot=False):
        """
        Different uses for this function to explain, explain the use of NreadsItrue and NreadsIItrue
        paras : 
        sparse_rep :
        acq_model_type : 
        """

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,NreadsI,NreadsII = sparse_rep

        alpha_rho = paras_1[0]
        fmin = np.power(10,paras_1[-1])
        freq_dtype = 'float64'
        nfbins = 1200 #Accuracy of the integration


        logrhofvec, logfvec = self.get_rhof(alpha_rho, nfbins, fmin, freq_dtype)

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
                logPn1_f = self.get_logPn_f(unicounts, Nreads, logfvec_tmp, acq_model_type, paras)

            else:
                logPn2_f = self.get_logPn_f(unicounts, Nreads, logfvec_tmp, acq_model_type, paras)

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
        print("N_obs: " + str(N_obs))
    
            
        def cost(PARAS):

            alp = PARAS[0]
            sbar = PARAS[1]

            Ps = self.get_Ps(alp,sbar,smax,s_step)
            Pn0n0=np.dot(Pn0n0_s,Ps)
            Pn1n2_ps=np.sum(Pn1n2_s*Ps[:,np.newaxis,np.newaxis],0)
            Pn1n2_ps/=1-Pn0n0
            #print(Pn0n0)

       

            Energy = - np.dot(sparse_rep_counts/float(N_obs),np.where(Pn1n2_ps[indn1,indn2]>0,np.log(Pn1n2_ps[indn1,indn2]),0)) 
                
            return Energy

    #--------------------------Compute-the-grid-----------------------------------------
        
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
 

    def save_table(self, outpath, svec, Ps,Pn1n2_s, Pn0n0_s,  subset, unicountvals_1_d, unicountvals_2_d, indn1_d, indn2_d, print_expanded, pthresh, smedthresh):
        '''
        takes learned diffexpr model, Pn1n2_s*Ps, computes posteriors over (n1,n2) pairs, and writes to file a table of data with clones as rows and columns as measures of thier posteriors 
        print_expanded=True orders table as ascending by , else descending
        pthresh is the threshold in 'p-value'-like (null hypo) probability, 1-P(s>0|n1_i,n2_i), where i is the row (i.e. the clone) n.b. lower null prob implies larger probability of expansion
        smedthresh is the threshold on the posterior median, below which clones are discarded
        '''

        Psn1n2_ps=Pn1n2_s*Ps[:,np.newaxis,np.newaxis] 
    
        #compute marginal likelihood (neglect renormalization , since it cancels in conditional below) 
        Pn1n2_ps=np.sum(Psn1n2_ps,0)

        Ps_n1n2ps=Pn1n2_s*Ps[:,np.newaxis,np.newaxis]/Pn1n2_ps[np.newaxis,:,:]
        #compute cdf to get p-value to threshold on to reduce output size
        cdfPs_n1n2ps=np.cumsum(Ps_n1n2ps,0)
    

        def dummy(row,cdfPs_n1n2ps,unicountvals_1_d,unicountvals_2_d):
            '''
            when applied to dataframe, generates 'p-value'-like (null hypo) probability, 1-P(s>0|n1_i,n2_i), where i is the row (i.e. the clone)
            '''
            return cdfPs_n1n2ps[np.argmin(np.fabs(svec)),row['Clone_count_1']==unicountvals_1_d,row['Clone_count_2']==unicountvals_2_d][0]
        dummy_part=partial(dummy,cdfPs_n1n2ps=cdfPs_n1n2ps,unicountvals_1_d=unicountvals_1_d,unicountvals_2_d=unicountvals_2_d)
    
        cdflabel=r'$1-P(s>0)$'
        subset[cdflabel]=subset.apply(dummy_part, axis=1)
        subset=subset[subset[cdflabel]<pthresh].reset_index(drop=True)

        #go from clone count pair (n1,n2) to index in unicountvals_1_d and unicountvals_2_d
        data_pairs_ind_1=np.zeros((len(subset),),dtype=int)
        data_pairs_ind_2=np.zeros((len(subset),),dtype=int)
        for it in range(len(subset)):
            data_pairs_ind_1[it]=np.where(int(subset.iloc[it].Clone_count_1)==unicountvals_1_d)[0]
            data_pairs_ind_2[it]=np.where(int(subset.iloc[it].Clone_count_2)==unicountvals_2_d)[0]   
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
        oldcolnames=( 'AACDR3',  'ntCDR3', 'Clone_count_1', 'Clone_count_2', 'Clone_fraction_1', 'Clone_fraction_2')
        newcolnames=('CDR3_AA', 'CDR3_nt',        r'$n_1$',        r'$n_2$',           r'$f_1$',           r'$f_2$')
        subset=subset.rename(columns=dict(zip(oldcolnames, newcolnames)))
    
        #select only clones whose posterior median pass the given threshold
        subset=subset[subset[r'$s_{2,med}$']>smedthresh]
    
        print("writing to: "+outpath)
        if print_expanded:
            subset=subset.sort_values(by=cdflabel,ascending=True)
            strout='expanded'
        else:
            subset=subset.sort_values(by=cdflabel,ascending=False)
            strout='contracted'
        subset.to_csv(outpath+'top_'+strout+'.csv',sep='\t',index=False)



    def expansion_table(self, outpath, paras_1, paras_2, sparse_rep, df, acq_model_type, pval_threshold, smed_threshold):

        """
        
        """

        L_surface, Pn1n2_s_d, Pn0n0_s_d, svec = self.learning_dynamics_expansion(paras_1, paras_2, sparse_rep, acq_model_type)
        npoints= 50 # same as in learning_dynamics_expansion
        smax = 25.0     
        s_step = 0.1
        alpvec = np.logspace(-3,np.log10(0.99), npoints)
        sbarvec = np.linspace(0.01,5, npoints)
        maxinds=np.unravel_index(np.argmax(L_surface),np.shape(L_surface))
        optsbar=sbarvec[maxinds[0]]
        optalp=alpvec[maxinds[1]]
        optPs= self.get_Ps(optalp,optsbar,smax,s_step)
        pval_expanded = True

        indn1,indn2,sparse_rep_counts, unicountvals_1, unicountvals_2, NreadsI, NreadsII = sparse_rep

        self.save_table(outpath, svec, optPs, Pn1n2_s_d, Pn0n0_s_d,  df, unicountvals_1, unicountvals_2, indn1, indn2, pval_expanded, pval_threshold, smed_threshold)





    