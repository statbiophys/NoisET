
"""
Utilities for the classes in noisettes.py

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

import pandas as pd
import numpy as np
import math


def data_parser(data, count_labels, id_label=None):
    """
    Load the data and convert it to the proper DataFrame. Data can be a 2d array or
    a dataframe containing count_labels[0] and count_labels[1] columns
    """

    # Dataframe parsing
    if type(data) == pd.core.frame.DataFrame:

        if count_labels[0] not in data or count_labels[1] not in data:
            raise Exception("The DataFrame does not contains columns with the \
                specifeid count names")

        if type(id_label) == str:
            if id_label not in data:
                raise Exception("The DataFrame does not contains columns with \
                    the specifeid id name")
            data = data.rename(columns={id_label: 'id'})

        data = data.rename(columns={count_labels[0]: 'counts1', count_labels[1]: 'counts2'})

    # 2d list parsing
    elif type(data) == list or type(data) == np.ndarray:
        data = np.array(data)
        if len(data.shape) != 2 or data.shape[1] != 2:
            raise Exception("Invalid list dimension")
        data = pd.DataFrame(data, columns=['counts1', 'counts2'])

        # In the case of lists the ids are the indexes of the list
        data['id'] = np.arange(len(data))

    else:
        raise Exception("Invalid data format. It should be a DataFrame or a list")

    return data


def get_sparserep(df): 
    """
    Tranforms {(n1,n2)} data stored in pandas dataframe to a sparse 1D representation.
    unicountvals_1(2) are the unique values of n1(2).
    sparse_rep_counts gives the counts of unique pairs.
    ndn1(2) is the index of unicountvals_1(2) giving the value of n1(2) in that unique pair.
    len(indn1)=len(indn2)=len(sparse_rep_counts)
    """
    
    counts = df.loc[:,['counts1', 'counts2']]
    counts['paircount'] = 1  # gives a weight of 1 to each observed clone

    clone_counts = counts.groupby(['counts1', 'counts2']).sum()
    sparse_rep_counts = np.asarray(clone_counts.values.flatten(), dtype=int)
    clonecountpair_vals = clone_counts.index.values
    indn1 = np.asarray([clonecountpair_vals[it][0] for it in range(len(sparse_rep_counts))], dtype=int)
    indn2 = np.asarray([clonecountpair_vals[it][1] for it in range(len(sparse_rep_counts))], dtype=int)
    NreadsI = np.sum(counts['counts1'])
    NreadsII = np.sum(counts['counts2'])

    unicountvals_1, indn1 = np.unique(indn1, return_inverse=True)
    unicountvals_2, indn2 = np.unique(indn2, return_inverse=True)

    return indn1, indn2, sparse_rep_counts, unicountvals_1, unicountvals_2, NreadsI, NreadsII


def NegBinPar(m,v,mvec): 
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


def NegBinParMtr(m,v,nvec): #speed up only insofar as the log and exp are called once on array instead of multiple times on rows
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


def PoisPar( Mvec,unicountvals):
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


def get_rhof(alpha_rho, nfbins, fmin, freq_dtype='float64'):
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
    return logrhovec, logfvec, normconst


def get_logPn_f(unicounts, Nreads, logfvec, noise_model, paras):

    # Choice of the model:
    
    if noise_model<1:

        m_total=float(np.power(10, paras[3])) 
        r_c=Nreads/m_total
    if noise_model<2:

        beta_mv= paras[1]
        alpha_mv=paras[2]
        
    if noise_model<1: #for models that include cell counts
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
    if noise_model==0:

        mean_m=m_total*np.exp(logfvec)
        var_m=mean_m+beta_mv*np.power(mean_m,alpha_mv)
        Poisvec = PoisPar(mvec*r_c,unicounts)
        for f_it in range(len(logfvec)):
            NBvec=NegBinPar(mean_m[f_it],var_m[f_it],mvec)
            for n_it,n in enumerate(unicounts):
                Pn_f[f_it,n_it]=np.dot(NBvec[m_low[n_it]:m_high[n_it]+1],Poisvec[m_low[n_it]:m_high[n_it]+1,n_it]) 
    
    elif noise_model==1:

        mean_n=Nreads*np.exp(logfvec)
        var_n=mean_n+beta_mv*np.power(mean_n,alpha_mv)
        Pn_f = NegBinParMtr(mean_n,var_n,unicounts)
    elif noise_model==2:

        mean_n=Nreads*np.exp(logfvec)
        Pn_f= PoisPar(mean_n,unicounts)
    else:
        print('acq_model is 0,1, or 2 only')

    return np.log(Pn_f)
