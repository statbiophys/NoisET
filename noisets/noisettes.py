"""
# NoisET<sup>*</sup>  NOIse sampling learning & Expansion detection of T-cell receptors using Bayesian inference.

High-throughput sequencing of T- and B-cell receptors makes it possible to track immune
repertoires across time, in different tissues, in acute and chronic diseases or in healthy individuals. However
quantitative comparison between repertoires is confounded by variability in the read count of each receptor
clonotype due to sampling, library preparation, and expression noise. We present an easy-to-use python
package NoisET that implements and generalizes a previously developed Bayesian method in [Puelma Touzel et al, 2020](<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007873&rev=2>). It can be used
to learn experimental noise models for repertoire sequencing from replicates, and to detect responding
clones following a stimulus. The package was tested on different repertoire sequencing technologies and
datasets. NoisET package is desribed [here](<https://arxiv.org/abs/2102.03568>).

<sup>* NoisET should be pronounced as "noisettes" (ie hazelnuts in French).</sup>

Functions library for NoisET - construction of noisettes package
Copyright (C) 2021 Meriem Bensouda Koraichi.
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

# Installation

Python 3

NoisET is a python /3.6 software. It is available on PyPI and can be downloaded and installed through pip:

```console
$ pip install noisets
```
Watch out, Data pre-processing, diversity estimates and generation of neutral TCR clonal dynamics is not possible yet with installation with pip. Use only the sudo command below.

To install NoisET and try the tutorial dusplayed in this github: gitclone the file in your working environment.
Using the terminal, go to NoisET directory and write the following command :

```console
$ sudo python setup.py install
```

If you do not have the following python libraries (that are useful to use NoisET) : numpy, pandas, matplotlib, seaborn, scipy, scikit-learn, please do the following commands, to try first to install the dependencies separately: :
 ```
python -m pip install -U pip
python -m pip install -U matplotlib
pip install numpy
pip install pandas
pip install matplotlib
pip install seaborn
pip install -U scikit-learn

 ```
# Documentation

## Command lines with terminal

A tutorial is available at https://github.com/mbensouda/NoisET_tutorial .
Three commands are available to use :
- `noiset-noise` To infer Null noise model: NoisET first function (1)
- `noiset-nullgenerator` To qualitatively check consistency of NoisET first function
- `noiset-detection` To detect responding clones to a stimulus: NoisET second function (2)

All options are described typing one of the previous commands + `--help`or `-h`. Options are also described in the following READme.

## 1/ Infer noise model

To infer null noise model: NoisET first function (1), use the command `noiset-noise`

At the command prompt, type:
```console
$ noiset-noise --path 'DATA_REPO/' --f1 'FILENAME1_X_REP1' --f2 'FILENAME2_X_REP2' --(noisemodel)
```

Several options are needed to learn noise model from two replicate samples associated to one individual at a specific time point:

#### 1/ Data information:

- `--path 'PATHTODATA'`: set path to data file
- `--f1 'FILENAME1_X_REP1'`: filename for individual X replicate 1
- `--f2 'FILENAME2_X_REP2'`: filename for individual X replicate 2

If your TCR CDR3 clonal populations features (ie clonal fractions, clonal counts, clonal nucleotide CDR3 sequences and clonal amino acid sequences) have different column names than: ('Clone fraction', 'Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3), you can specify the name directly by using:

- `--specify`
- `--freq 'frequency'` : Column label associated to clonal fraction
- `--counts 'counts'`:  Column label associated to clonal count
- `--ntCDR3 'ntCDR3'`:  Column label associated to clonal CDR3 nucleotides sequence
- `--AACDR3 'AACDR3'`:  Column label associated to clonal CDR3 amino acid sequence

#### 2/ Choice of noise model: (parameters meaning described in Methods section)
- `--NBPoisson`: Negative Binomial + Poisson Noise Model - 5 parameters
- `--NB`: Negative Binomial - 4 parameters
- `--Poisson`: Poisson - 2 parameters

#### 3/ Example:
At the command prompt, type:
```console
$ noiset-noise --path 'data_examples/' --f1 'Q1_0_F1_.txt.gz' --f2 'Q1_0_F2_.txt.gz' --NB
```
This command line will learn four parameters associated to negative binomial null noise Model `--NB` for individual Q1 at day 0.
A '.txt' file is created in the working directory: it is a 5/4/2 parameters data-set regarding on NBP/NB/Poisson noise model. In this example, it is a four parameters table (already created in data_examples repository).
You can run previous examples using data (Q1 day 0/ day15) provided in the data_examples folder - data from [Precise tracking of vaccine-responding T cell clones reveals convergent and personalized response in identical twins, Pogorelyy et al, PNAS](https://www.pnas.org/content/115/50/12704)

#### 4/ Example with `--specify`:

At the command prompt, type:
```console
$ noiset-noise --path 'data_examples/' --f1 'replicate_1_1.tsv.gz' --f2 'replicate_1_2.tsv.gz' --specify --freq 'frequencyCount' --counts 'count' --ntCDR3 'nucleotide' --AACDR3 'aminoAcid' --NB
```
As previously this command enables us to learn four parameters associated to negative binomial null noise model `--NB` for one individual in cohort produced in [Model to improve specificity for identification of clinically-relevant expanded T cells in peripheral blood, Rytlewski et al, PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0213684).

## 2/ Generate synthetic data from null model learning:

To qualitatively check consistency of NoisET first function (1) with experiments or for other reasons, it can be useful to generates synthetic replicates from the null model (described in Methods section).
One can also generalte healthy RepSeq samples dynamics using the noise model which has been learned in a first step anf giving the time-scale dynamics of turnover of the repertoire as defined in https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. Check [here](<https://github.com/statbiophys/NoisET/blob/master/NoisET%20example%20-%20Null%20model%20learning%20.ipynb>).

To generate synthetic TCR RepSeq data replicates having chosen sampling noise characteristics, use the command `noiset-nullgenerator`

 ```console
 $ noiset-nullgenerator --(noise-model) --nullpara 'NULLPARAS' --nreads_1 float --nreads_2 float --Nclones float --output 'SYNTHETICDATA'
 ```

#### 1/ Choice of noise model:
The user must chose one of the three possible models for the probability that a TCR has <strong> an empirical count n </strong> knowing that its  <strong> true frequency is f </strong>, P(n|f): a Poisson distribution `--Poisson`, a negative binomial distribution `--NB`, or a two-step model combining Negative-Binomial and a Poisson distribution `--NBP`. n is the empirical clone size and  depends on the experimental protocol.
For each P(n|f), a set of parameters is learned.

- `--NBPoisson`: Negative Binomial + Poisson Noise Model - 5 parameters 5 parameters described in [Puelma Touzel et al, 2020](<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007873&rev=2>): power-law exponent of clonotypes frequencies distributions `'alph_rho'`, minimum of clonotype frequencies distribution `'fmin'`, `'beta'` and `'alpha'`, parameters of negative binomial distribution constraining mean and variance of P(m|f) distribution (m being the number of cells associated to a clonotype in the experiemental sample), and `'m_total'` the total number of cells in the sample of interest..
- `--NB`: Negative Binomial - 4 parameters: power-law of the clonotypes frequencies distributions (same ansatz than in [Puelma Touzel et al, 2020](<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007873&rev=2>) `'alph_rho'`, minimum of clonotype frequencies distribution `'fmin'`, `'beta'` and `'alpha'`, parameters of negative binomial distribution constraining mean and variance of P(n|f) distribution. <em> NB(fNreads, fNreads + betafNreads<sup>alpha</sup>) </em>. (Nreads is the total number of reads in the sample of interest.)
- `--Poisson`: Poisson - 2 parameters power-law of the clonotypes frequencies distributions (same ansatz than in [Puelma Touzel et al, 2020](<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007873&rev=2>)`'alph_rho'` and minimum of clonotype frequencies distribution `'fmin'`. P(n|f) is a Poisson distribution of parameter <em> fNreads </em>. (Nreads is the total number of reads in the sample of interest.)

#### 2/ Specify learned noise parameters:
- `--nullpara 'PATHTOFOLDER/NULLPARAS.txt'`: parameters learned thanks to NoisET function (1) \
!!! Make sure to match correctly the noise model and the null parameter file content : 5 parameters for `--NBP`, 4 parameters for `--NB`and 2 parameters
for `--Poisson`.

#### 3/ Sequencing properties of data:
- `--nreads_1 NNNN`: total number  of reads in first replicate - it should match the actual data. In the example below, it is the sum of 'Clone count' in 'Q1_0_F1_.txt.gz'.
- `--Nreads2 NNNN`: total number  of reads in second replicate - it should match the actual data. In the example below, it is the sum of 'Clone count' in 'Q1_0_F2_.txt.gz'.
- `--Nclones NNNN`: total number of clones in union of two replicates - it should match the actual data. In the example below, it is the number of clones present in both replicates : 'Q1_0_F1_.txt.gz' and 'Q1_0_F2_.txt.gz'.

#### 4/ Output file
`--output 'SYNTHETICDATA'`: name of the output file where you can find the synthetic data set.

At the command prompt, type
 ```console
 $ noiset-nullgenerator --NB --nullpara 'data_examples/nullpara1.txt' --nreads_1 829578 --nreads_2 954389 --Nclones 776247 --output 'test'
 ```
 Running this line, you create a 'synthetic_test.csv' file with four columns : 'Clone_count_1', 'Clone_count_2', 'Clone_fraction_1', 'Clone_fraction_2', resctively synthetic read counts and frequencies that you would have found in an experimental sample of same learned parameters 'nullpara1.txt', 'nreads_1', 'nreads_2' and 'Nclones'.

## 3/ Detect responding clones:

Detects responding clones to a stimulus: NoisET second function (2)

To detect responding clones from two RepSeq data at time_1 and time_2, use the command `noiset-detection`

```console
$ noiset-detection --(noisemodel)  --nullpara1 'FILEFORPARAS1' --nullpara2 'FILEFORPARAS1' --path 'REPO/' --f1 'FILENAME_TIME_1' --f2 'FILENAME_TIME_2' --pval float --smedthresh float --output 'DETECTIONDATA'
```
Several options are needed to learn noise model from two replicate samples associated to one individual at a specific time point:

#### 1/ Choice of noise model:
- `--NBPoisson`: Negative Binomial + Poisson Noise Model - 5 parameters
- `--NB`: Negative Binomial - 4 parameters
- `--Poisson`: Poisson - 2 parameters

#### 2/ Specify learned parameters for both time points:
(they can be the same for both time points if replicates are not available but to use carefully as mentioned in [ARTICLE])
- `--nullpara1 'PATH/FOLDER/NULLPARAS1.txt'`: parameters learned thanks to NoisET function (1) for time 1
- `--nullpara2 'PATH/FOLDER/NULLPARAS2.txt'`: parameters learned thanks to NoisET function (1) for time 2

!!! Make sure to match correctly the noise model and the null parameters file content : 5 parameters for `--NBP`, 4 parameters for `--NB`and 2 parameters
for `--Poisson`.

#### 3/ Data information:

- `--path 'PATHTODATA'`: set path to data file
- `--f1 'FILENAME1_X_time1'`: filename for individual X time 1
- `--f2 'FILENAME2_X_time2'`: filename for individual X time 2

If your TCR CDR3 clonal populations features (ie clonal fractions, clonal counts, clonal nucleotides CDR3 sequences and clonal amino acids sequences) have different column names than: ('Clone fraction', 'Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3), you can specify the name by using:

- `--specify`
- `--freq 'frequency'` : Column label associated to clonal fraction
- `--counts 'counts'`:  Column label associated to clonal count
- `--ntCDR3 'ntCDR3'`:  Column label associated to clonal CDR3 nucleotides sequence
- `--AACDR3 'AACDR3'`:  Column label associated to clonal CDR3 amino acid sequence

#### 4/ Detection thresholds: (More details in Methods section).
- `--pval XXX` : p-value threshold for the expansion/contraction - use 0.05 as a default value.
- `--smedthresh XXX` : log fold change median threshold for the expansion/contraction - use 0 as a default value.

#### 5/ Output file
`--output 'DETECTIONDATA'`: name of the output file (.csv) where you can find a list of the putative responding clones with statistics features. (More details in Methods section).


At the command prompt, type
```console
$ noiset-detection --NB  --nullpara1 'data_examples/nullpara1.txt' --nullpara2 'data_examples/nullpara1.txt' --path 'data_examples/' --f1 'Q1_0_F1_.txt.gz' --f2 'Q1_15_F1_.txt.gz' --pval 0.05 --smedthresh 0 --output 'detection'
```

Ouput: table containing all putative detected clones with statistics features about logfold-change variable <em> s </em>: more theoretical description [Puelma Touzel et al, 2020](<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007873&rev=2>).

## Python package

"""
from copy import deepcopy
from datetime import datetime, date
from decimal import Decimal
from functools import partial
import math
from multiprocessing import Pool, cpu_count
import os
import shutil
import time
from typing import Callable, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import cm, colors, colorbar
import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
import seaborn as sns
from scipy.special import gammaln
from scipy.stats import nbinom
from scipy.stats import poisson
from scipy.stats import rv_discrete
from scipy.optimize import minimize, OptimizeResult
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering

###===================================TOOLS-TO-GENERATE-NEUTRAL-TCR-REP-SEQ-TRAJECTORIES=====================================================
#  Library functions to generate TCR repertoires
##------------------------Initial-Distributions------------------------
def _rho_counts_theo_minus_x(A, B, N_0):
    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user 

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

def _rho_counts_theo_plus_x(A, B, N_0):
    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user 
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


def _get_distsample(pmf,Nsamp, dtype='uint32'):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user 
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

##------------------------Propagator------------------------
def _gaussian_matrix(x_vec, x_i_vec_unique, A, B, t):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user
    x_vec_reshaped = np.reshape(x_vec, (len(x_vec), 1))
    ones_vec = np.ones((len(x_i_vec_unique), 1))
    M = np.multiply(ones_vec, x_vec_reshaped.T)
    x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))

    return (1/np.sqrt(2*np.pi*B*t))*np.exp((-1/(2*B*t))*(M - x_i_unique_reshaped - A*t)**2)

def _gaussian_adsorption_matrix(x_vec, x_i_vec_unique, A, B, t):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user
    a = 0
    gauss = _gaussian_matrix(x_vec, x_i_vec_unique, A, B, t)
    gauss_a = _gaussian_matrix(x_vec, 2*a-x_i_vec_unique, A, B, t)
    x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))
    return gauss - np.exp((A*(a-x_i_unique_reshaped))/(B/2)) * gauss_a

def _extinction_vector(x_i, A, B, t):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user
    nbins = 2000
    eps = 1e-20
    #eps = 0
    x_vec = np.linspace(eps, np.max(x_i) - A*t + 3*np.sqrt(B*t), nbins)

    x_i_sorted = np.sort(x_i)

    xiind_vals, xi_start_ind, xi_counts=np.unique(x_i_sorted, return_counts=True,return_index=True)
    Prop_Matrix = _gaussian_adsorption_matrix(x_vec, xiind_vals, A, B, t)

    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    p_ext = 1 - integ

    p_ext_new = np.zeros((len(x_i)))
    for it,xiind in enumerate(xiind_vals):
        p_ext_new[xi_start_ind[it]:xi_start_ind[it]+xi_counts[it]] = p_ext[it]

    test = np.random.uniform(0,1, size = (len(p_ext_new))) > p_ext_new
    results_extinction = test.astype(int)

    return results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext

#------------------------Source-term-no-frequency-dependency------------------------

def _gaussian_matrix_time(x_vec, x_i_scal, A, B, tvec_unique):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user

    x_vec_reshaped = np.reshape(x_vec, (len(x_vec), 1))
    ones_vec = np.ones((len(tvec_unique), 1))
    M = np.multiply(ones_vec, x_vec_reshaped.T)
    tvec_unique_reshaped = np.reshape(tvec_unique, (len(tvec_unique), 1))
    #x_i_unique_reshaped = np.reshape(x_i_vec_unique, (len(x_i_vec_unique), 1))

    return (1/np.sqrt(2*np.pi*B*tvec_unique_reshaped))*np.exp((-1/(2*B*tvec_unique_reshaped))*(M - x_i_scal - A*tvec_unique_reshaped)**2)

def _gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tvec_unique):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user

    a = 0
    gauss = _gaussian_matrix_time(x_vec, x_i_scal, A, B, tvec_unique)
    gauss_a = _gaussian_matrix_time(x_vec, 2*a-x_i_scal, A, B, tvec_unique)

    return gauss - np.exp((A*(a-x_i_scal))/(B/2)) * gauss_a

def _Prop_Matrix_source( A, B, tvec):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user

    nbins = 2000
    N_0 = 40
    x_i_scal = np.log(N_0)
    t = np.max(tvec)
    x_vec = np.linspace(0, x_i_scal - A*t + 2*np.sqrt(B*t), nbins)

    tvec_sorted = np.sort(tvec)

    tiind_vals, ti_start_ind, ti_counts=np.unique(tvec_sorted, return_counts=True,return_index=True)
    Prop_Matrix = _gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tiind_vals)

    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)

    return Prop_Matrix, x_vec, tiind_vals, ti_start_ind, ti_counts, integ

def _extinction_vector_source(A, B, tvec):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user

    nbins = 2000
    N_0 = 40
    x_i_scal = np.log(N_0)
    t = np.max(tvec)
    x_vec = np.linspace(0, x_i_scal - A*t + 2*np.sqrt(B*t), nbins)

    tvec_sorted = np.sort(tvec)

    tiind_vals, ti_start_ind, ti_counts=np.unique(tvec_sorted, return_counts=True,return_index=True)
    Prop_Matrix = _gaussian_adsorption_matrix_time(x_vec, x_i_scal, A, B, tiind_vals)

    dx =np.asarray(np.diff(x_vec)/2., dtype='float32')
    integ = np.sum(dx*(Prop_Matrix[:, 1:]+Prop_Matrix[:, :-1]), axis = 1)
    p_ext = 1 - integ

    p_ext_new = np.zeros((len(tvec)))
    for it,tiind in enumerate(tiind_vals):
        p_ext_new[ti_start_ind[it]:ti_start_ind[it]+ti_counts[it]] = p_ext[it]

    test = np.random.uniform(0,1, size = (len(p_ext_new))) > p_ext_new
    results_extinction = test.astype(int)


    return results_extinction, Prop_Matrix, x_vec, tiind_vals, ti_start_ind, ti_counts, p_ext

##------------------------Function-to-generate-in-silico-Rep-Seq-samples------------------------

def _generator_diffusion_LB(A, B, N_0, t):

    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user

    eps = 1e-20

    ## Choose initial size of the immune system to be 1e10 (for a mouse)
    N_cells = int(1e10)

    #Parameters for the repertoire generation
    alpha_rho = -1 + 2*A/B
    N_ext = 1
    freq_dtype = 'float32'

    #==========================generate the steady state distribution===============================

    #for counts < N0: 
    logrhofvec,logfvec, N_clones_1 = _rho_counts_theo_minus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=_get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), N_clones_1,dtype='uint32').flatten()
    #print("generation population smaller than N_0: check")

    logcvec_generated = logfvec[f_samples_inds]
    counts_generated = np.exp(logcvec_generated)
    C_f_minus = np.sum(counts_generated)
    print(str(C_f_minus) + ' cells smaller than N_0')
    log_cminus_generated = logcvec_generated
    logrhofvec_1,logfvec_1 = logrhofvec,logfvec
    print(str(N_clones_1) + ' N_clones_1')

    #for counts > N0:

    logrhofvec,logfvec, N_clones_2 = _rho_counts_theo_plus_x(A, B, N_0)
    dlogfby2=np.asarray(np.diff(logfvec)/2., dtype='float32')
    integ=np.exp(logrhofvec[np.newaxis,:])
    f_samples_inds=_get_distsample(np.asarray((dlogfby2[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'),N_clones_2,dtype='uint32').flatten()
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



    results_extinction, Prop_Matrix, x_vec, xiind_vals, xi_start_ind, xi_counts, p_ext = _extinction_vector(x_i, A, B, t)
    #x_vec = np.linspace(0, 30*B*t, 2000)
    dx=np.asarray(np.diff(x_vec)/2., dtype='float32')

    x_i_noext= x_i[np.where(results_extinction ==1)]
    x_f = np.zeros((len(x_i)))

    for i in range(len(xiind_vals)):
        if (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1])) < 1e-7:
            pass
        else:

            Prop_adsorp = Prop_Matrix[i,:] / (np.dot(dx, Prop_Matrix[i,1:] + Prop_Matrix[i, :-1]))

            integ = Prop_adsorp[np.newaxis,:]
            f_samples_inds = _get_distsample(np.asarray((dx[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), xi_counts[i],dtype='uint32').flatten()

            x_f[xi_start_ind[i]:xi_start_ind[i]+xi_counts[i]] = x_vec[f_samples_inds]

    x_f = np.multiply(x_f,results_extinction)


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

    results_extinction_source, Prop_Matrix_source, x_vec_source, tiind_vals, ti_start_ind, ti_counts, p_ext_source = _extinction_vector_source(A, B, time_vec)

    dx_source=np.asarray(np.diff(x_vec_source)/2., dtype='float32')

    x_source_LB = np.zeros((N_source))
    for i in range(len(tiind_vals)):

        if (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1])) < 1e-7:
            pass

        else:
            Prop_adsorp_s = Prop_Matrix_source[i,:]
            Prop_adsorp_s = Prop_Matrix_source[i,:] / (np.dot(dx_source, Prop_Matrix_source[i,1:] + Prop_Matrix_source[i, :-1]))


            integ = Prop_adsorp_s[np.newaxis,:]
            f_samples_inds_s = _get_distsample(np.asarray((dx_source[np.newaxis,:]*(integ[:,1:]+integ[:,:-1])).flatten(),dtype='float32'), ti_counts[i],dtype='uint32').flatten()

            x_source_LB[ti_start_ind[i]:ti_start_ind[i]+ti_counts[i]] = x_vec_source[f_samples_inds_s]


    x_source_LB = np.multiply(x_source_LB, results_extinction_source)

    x_source_LB[x_source_LB == 0] = -np.inf



    return x_i, x_f, Prop_Matrix, p_ext, results_extinction, time_vec, results_extinction_source, x_source_LB

def _experimental_sampling_diffusion_Poisson(nreads_1, nreads_2, x_0, x_2, t, N_cell_0, N_cell_2):


    # This function has been made to generate TCR clonal frequencies distribution from the theoretical model described in the paper 
    # https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1. 
    # this function is not made for a NoisET user


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
    m=float(nreads_1)*f_vec_initial
    n_counts_day_0 = np.random.poisson(m, size =(1, int(N_total)))
    n_counts_day_0 = n_counts_day_0[0,:]

    #print('done')

    #Final condition
    f_vec_end = np.exp(x_2_final)/N_cell_2
    m=float(nreads_2)*f_vec_end
    #print(m)
    print('MEAN N : ' + str(np.mean(m)))
    n_counts_day_1 = np.random.poisson(m, size =(1, int(N_total)))
    print(n_counts_day_1)
    n_counts_day_1 = n_counts_day_1[0,:]


    #-------------------------------Creation of the data set-------------------------------------

    obs=np.logical_or(n_counts_day_0>0, n_counts_day_1>0)
    n1_samples=n_counts_day_0[obs]
    n2_samples=n_counts_day_1[obs]
    pair_samples_df= pd.DataFrame({'Clone_count_1':n1_samples,'Clone_count_2':n2_samples})

    pair_samples_df['Clone_frequency_1'] = pair_samples_df['Clone_count_1'] / np.sum(pair_samples_df['Clone_count_1'])
    pair_samples_df['Clone_frequency_2'] = pair_samples_df['Clone_count_2'] / np.sum(pair_samples_df['Clone_count_2'])


    return pair_samples_df

def _experimental_sampling_diffusion_NegBin(nreads_1, nreads_2, paras, x_0, x_2, N_cell_0, N_cell_2):


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
    m=float(nreads_1)*f_vec_initial
    print(m)

    beta_mv=paras[1]
    alpha_mv=paras[2]

    v=m+beta_mv*np.power(m,alpha_mv)

    pvec=1-m/v
    nvec=m*m/v/pvec

    pvec = np.nan_to_num(pvec, nan=0.0)
    nvec = np.nan_to_num(nvec, nan=1e-30)

    print(pvec)
    print(1-pvec)
    print(np.sum(pvec>=1))
    print(nvec)

    n_counts_day_0 = np.random.negative_binomial(nvec, 1-pvec, size =(1, int(N_total)))
    n_counts_day_0 = n_counts_day_0[0,:]
    print(n_counts_day_0)


    #Final condition
    f_vec_end = np.exp(x_2)/N_cell_2
    m_end=float(nreads_2)*f_vec_end
    print(m_end)

    v_end=m_end+beta_mv*np.power(m_end,alpha_mv)
    pvec_end=1-m_end/v_end
    nvec_end=m_end*m_end/v_end/pvec_end

    pvec_end = np.nan_to_num(pvec_end, nan=0.0)
    nvec_end = np.nan_to_num(nvec_end, nan=1e-30)


    n_counts_day_1 = np.random.negative_binomial(nvec_end, 1-pvec_end, size =(1, int(N_total)))
    n_counts_day_1 = n_counts_day_1[0,:]
    print(n_counts_day_1)


    #-------------------------------Creation of the data set-------------------------------------

    obs=np.logical_or(n_counts_day_0>0, n_counts_day_1>0)
    n1_samples=n_counts_day_0[obs]
    n2_samples=n_counts_day_1[obs]
    pair_samples_df= pd.DataFrame({'Clone_count_1':n1_samples,'Clone_count_2':n2_samples})

    pair_samples_df['Clone_frequency_1'] = pair_samples_df['Clone_count_1'] / np.sum(pair_samples_df['Clone_count_1'])
    pair_samples_df['Clone_frequency_2'] = pair_samples_df['Clone_count_2'] / np.sum(pair_samples_df['Clone_count_2'])


    return pair_samples_df

#==========================================================================================================================


#===============================Longitudinal-Data-Pre-Processing===================================

class longitudinal_analysis():

    """
    This class provides some tool to inspect and compute some simple statistics on longitudinal data associated with
    one individual (it is independent of the NoisET software).

    ...

    Attributes
    ----------
    clone_count_label : str
        label in the clonotype tables indicating the clonotype count
    seq_label : str
        label in the clonotype tables indicating the sequence of the receptor
    clones : dict of pandas.DataFrame
        dictionary containing the clonotype tables as pandas frames. The keys are
        strings "patient_time", replicated are merged. Created in the initalization
    times : list of float
        ordered times of the imported tables. Created in the initialization
    unique_clones : list of str
        list of all the unique clonotype sequences in all the time points
    time_occurrence : list of int
        number of time points in which each clonotype appears. The index
        refers to the clonotype in the unique_clones list

    Methods
    -------

    compute_clone_time_occurrence()
        It creates two new attribues: the list of uniqe clonotypes in all the dataset
        "self.unique_clones" and the time occurrence of each of them "self.time_occurrence".
        the time occurrence is the number of time points in which the clone appears.

    plot_hist_persistence(figsize=(12,10))
        It plots the distribution of time occurrence of the unique clonotypes

    top_clones_set(n_top_clones)
        Compute the set of top clones as the union of the "n_top_clones" most abundant
        clonotype in each time point

    build_traj_frame(top_clones_set)
        Compute the set of top clones as the union of the "n_top_clones" most abundant
        clonotype in each time point

    plot_trajectories(n_top_clones, colormap='viridis', figsize=(12,10))
        Function to plot the trajectories of the first "n_top_clones". Colors of the
        trajectories represent the cumulative frequency in all the time points.

    PCA_traj(n_top_clones, nclus=4)
        Perform PCA over the normalized trajectories of n_top_clones TCR clones.
        The normalization consists in dividing the whole trajectory by its maximum value.
        After PCA the trajectories are clustered in the two principal componets space
        with a hierarchical clustering algorithm.

    plot_PCA2(n_top_clones, nclus=4, colormap='viridis', figsize=(12,10))
        Plotting the trajectories in the space of their two principal components and
        clustering them as in "PCA_traj".

    plot_PCA_clusters_traj(n_top_clones, nclus=4, colormap='viridis', figsize=(12,10))
        Plotting the trajectories grouped by PCA clusters
    """




    def __init__(self, patient, data_folder, sequence_label='N. Seq. CDR3', clone_count_label='Clone count',
                 replicate1_label='_F1', replicate2_label='_F2', separator='\t'):
        """
        Import all the clonotypes of a given patient and store them in the dictionary "self.clones".
        It also creates the list of times "self.times". During this process the replicates at the
        same time points are merged together.
        The names of the tables containing TCR should be structured as "patient_time_replicate.csv".
        Those tables should be cvs files compressed in a zip archive (see the example notebook).

        Parameters
        ----------
        patient : str
            The ID of the patient
        data_folder : str
            folder name containing the csv files listing the T-cell receptors
        separator : str
            separator symbol in the csv tables
        """

        self.clone_count_label = clone_count_label
        self.seq_label = sequence_label
        self.unique_clones = None
        self.time_occurrence = None
        self.times = []
        clones_repl = dict()

        # Iteration over all the file in the folder for importing each table
        for file_name in os.listdir(data_folder):
        # If the name before the underscore corresponds to the chosen patient..
            if file_name.split('_')[0] == patient:
                # Import the table
                frame = pd.read_csv(data_folder+file_name, sep='\t', compression=dict(method='zip'))
                # Store it in a dictionary where the key contains the patient, the time
                # and the replicate.
                clones_repl[file_name[:-10]] = frame
                # Reading the time from the name and storing it
                self.times.append(int(file_name.split('_')[1]))
                print('Clonotypes',file_name[:-10],'imported')

        # Sorting the unique times
        self.times = np.sort(list(set(self.times)))
        self.clones = self._merge_replicates(patient, clones_repl, replicate1_label, replicate2_label)


    def _merge_replicates(self, patient, clones_repl, repl1_label, repl2_label):

        clones_merged = dict()

        # Iteration over the times
        for it, t in enumerate(self.times):
            # Building the ids correponding at 1st and 2nd replicate at given time point
            id_F1 = patient + '_' + str(t) + repl1_label
            id_F2 = patient + '_' + str(t) + repl2_label
            # Below all the rows of one table are appended to the rows of the other
            merged_replicates = clones_repl[id_F1].merge(clones_repl[id_F2], how='outer')
            # But there are common clonotypes that now appear in two different rows 
            # (one for the first and one for the second replicate)! 
            # Below we collapse those common sequences and the counts of the two are summed 
            merged_replicates = merged_replicates.groupby(self.seq_label, as_index=False).agg({self.clone_count_label:sum})
            depth = merged_replicates[self.clone_count_label].sum()
            merged_replicates['Clone freq'] = merged_replicates[self.clone_count_label] / depth
            merged_replicates = merged_replicates.sort_values('Clone freq', ascending=False)
            # The merged table is then added to the dictionary
            clones_merged[patient + '_' + str(t)] = merged_replicates

        return clones_merged


    def compute_clone_time_occurrence(self):

        """
        It creates two new attribues: the list of uniqe clonotypes in all the dataset
        "self.unique_clones" and the time occurrence of each of them "self.time_occurrence".
        the time occurrence is the number of time points in which the clone appears.
        """

        all_clones = np.array([])
        for id_, cl in self.clones.items():
            all_clones = np.append(all_clones, cl[self.seq_label].values)

        # The following function returns the list of unique clonotypes and the number of
        # repetitions for each of them. 
        # Note that the number of repetitions is exactly the time occurrence
        self.unique_clones, self.time_occurrence = np.unique(all_clones, return_counts=True)


    def plot_hist_persistence(self, figsize=(12,10)):

        """
        It plots the distribution of time occurrence of the unique clonotypes

        Parameters
        ----------
        figsize : tuple
            width, height in inches

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            axes where to draw the plot
        fig : matplotlib.figure.Figure
            matplotlib figure
        """

        if type(self.unique_clones) != np.ndarray:
            self.compute_clone_time_occurrence()

        fig, ax = plt.subplots(figsize=figsize)

        plt.rc('xtick', labelsize = 30)
        plt.rc('ytick', labelsize = 30)

        ax.set_yscale('log')
        ax.set_xlabel('Time occurrence', fontsize = 30)
        ax.set_ylabel('Counts', fontsize = 30)
        ax.hist(self.time_occurrence, bins=np.arange(1,len(self.times)+2)-0.5, rwidth=0.6)

        return fig, ax


    def top_clones_set(self, n_top_clones):

        """
        Compute the set of top clones as the union of the "n_top_clones" most abundant
        clonotype in each time point

        Parameters
        ----------
        n_top_clones : int
            number of most abundant clontypes in each time point

        Returns
        -------
        top_clones : set of str
            set of top clones
        """

        top_clones = set()
        for id_, cl in self.clones.items():
            top_clones_at_time = cl.sort_values(self.clone_count_label, ascending=False)[:n_top_clones]
            top_clones = top_clones.union(top_clones_at_time[self.seq_label].values)
        return top_clones


    def build_traj_frame(self, clone_set):

        """
        This builds a dataframe containing the frequency at all the time points for each
        of the clonotypes specified in clone_set.
        The dataframe has also a field that contains the cumulative frequency.

        Parameters
        ----------
        clones_set : iterable of str
            list of clonotypes whose temporal trajectory is drawn

        Returns
        -------
        traj_frame : pandas.DataFrame
            dataframe containing the frequency at all the time points
        """

        traj_frame = pd.DataFrame(index=clone_set)
        traj_frame['Clone cumul freq'] = 0

        for id_, cl in self.clones.items():

            # Getting the time from the index of clones_merged
            t = id_.split('_')[1]
            # Selecting the clonotypes that are both in the frame at the given time 
            # point and in the list of top_clones_set
            top_clones_at_time = clone_set.intersection(set(cl[self.seq_label]))
            # Creating a sub-dataframe containing only the clone in top_clones_at_time
            clones_at_time = cl.set_index(self.seq_label).loc[top_clones_at_time]
            # Creating a new column in the trajectory frames for the counts at that time
            traj_frame['t'+str(t)] = traj_frame.index.map(clones_at_time['Clone freq'].to_dict())
            # The clonotypes not present at that time are NaN. Below we convert NaN in 0s
            traj_frame = traj_frame.fillna(0)
            # The cumulative count for each clonotype is updated
            traj_frame['Clone cumul freq'] += traj_frame['t'+str(t)]

        return traj_frame



    # Plot clonal trajectories


    def plot_trajectories(self, n_top_clones, colormap='viridis', figsize=(12,10)):

        """
        Function to plot the trajectories of the first "n_top_clones". Colors of the
        trajectories represent the cumulative frequency in all the time points.

        Parameters
        ----------
        n_top_clones : int
            number of most abundant clontypes in each time point

        colormap  : str
            colors of the trajectories

        figsize : tuple
            width, height in inches

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            axes where to draw the plot
        fig : matplotlib.figure.Figure
            matplotlib figure
        """

        cmap = cm.get_cmap(colormap)
        top_clones = self.top_clones_set(n_top_clones)
        traj_frame = self.build_traj_frame(top_clones)

        fig, ax = plt.subplots(figsize=figsize)
        plt.rc('xtick', labelsize = 30)
        plt.rc('ytick', labelsize = 30)
        ax.set_yscale('log')
        ax.set_xlabel('time', fontsize = 25)
        ax.set_ylabel('frequency', fontsize = 25)

        log_counts = np.log10(traj_frame['Clone cumul freq'].values)
        max_log_count = max(log_counts)
        min_log_count = min(log_counts)

        for id_, row in traj_frame.iterrows():
            traj = row.drop(['Clone cumul freq']).to_numpy()
            log_count = np.log10(row['Clone cumul freq'])
            norm_log_count = (log_count-min_log_count)/(max_log_count-min_log_count)
            plt.plot(self.times, traj, c=cmap(norm_log_count))


        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min(log_counts), vmax=max(log_counts)))
        cb = plt.colorbar(sm)
        cb.set_label('Log10 cumulative frequency', fontsize = 25)

        return fig, ax


    def PCA_traj(self, n_top_clones, nclus=4):

        """
        Perform PCA over the normalized trajectories of n_top_clones TCR clones.
        The normalization consists in dividing the whole trajectory by its maximum value.
        After PCA the trajectories are clustered in the two principal componets space
        with a hierarchical clustering algorithm.

        Parameters
        ----------
        n_top_clones : int
            number of most abundant clontypes in each time point to consider in the PCA

        nclus : float
            number of clusters


        Returns
        -------
        pca : sklearn.decomposition._pca.PCA
            object containing the result of the principal component analysis

        clustering : sklearn.cluster._agglomerative.AgglomerativeClustering
            object containing the result of the hierarchical clustering
        """

        #Getting the top n_top_clones clonotypes at each time point
        top_clones = self.top_clones_set(n_top_clones)
        #Building a trajectory dataframe
        traj_frame = self.build_traj_frame(top_clones)

        #Converting it in a numpy matrix
        traj_matrix = traj_frame.drop(['Clone cumul freq'], axis = 1).to_numpy()

        # Normalize each trajectory by its maximum
        norm_traj_matrix = traj_matrix/np.max(traj_matrix, axis=1)[:, np.newaxis]

        pca = PCA(n_components =2).fit(norm_traj_matrix.T)
        clustering = AgglomerativeClustering(n_clusters = nclus)
        clustering = clustering.fit(pca.components_.T)

        return pca, clustering


    def plot_PCA2(self, n_top_clones, nclus=4, colormap='viridis', figsize=(12,10)):

        """
        Plotting the trajectories in the space of their two principal components and
        clustering them as in "PCA_traj".

        Parameters
        ----------
        n_top_clones : int
            number of most abundant clontypes in each time point to consider in the PCA

        nclus : float
            number of clusters

        colormap : str
            colormap indicating the different clusters

        figsize : tuple
            width, height in inches

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            axes where to draw the plot
        fig : matplotlib.figure.Figure
            matplotlib figure
        """


        cmap = cm.get_cmap(colormap)
        pca, clustering = self.PCA_traj(n_top_clones, nclus)

        fig, ax = plt.subplots(figsize=figsize)
        ax.set_title('PCA components (%i trajs)' %pca.n_features_, fontsize = 25)
        ax.set_xlabel('First component (expl var: %3.2f)'%pca.explained_variance_ratio_[0], fontsize = 25)
        ax.set_ylabel('Second component (expl var: %3.2f)'%pca.explained_variance_ratio_[1], fontsize = 25)
        for c_ind in range(clustering.n_clusters):
            x = pca.components_[0][clustering.labels_ == c_ind]
            y = pca.components_[1][clustering.labels_ == c_ind]
            ax.scatter(x, y, alpha=0.2, color=cmap(c_ind/clustering.n_clusters))

        return fig, ax


    def plot_PCA_clusters_traj(self, n_top_clones, nclus=4, colormap='viridis', figsize=(12,10)):

        """
        Plotting the trajectories grouped by PCA clusters

        Parameters
        ----------
        n_top_clones : int
            number of most abundant clontypes in each time point to consider in the PCA

        nclus : float
            number of clusters

        colormap : str
            colormap indicating the different clusters

        figsize : tuple
            width, height in inches


        Returns
        -------
        axs : tuple of matplotlib.axes._subplots.AxesSubplot
            axis where to draw the plot
        fig : matplotlib.figure.Figure
            matplotlib figure
        """

        cmap = cm.get_cmap(colormap)
        pca, clustering = self.PCA_traj(n_top_clones, nclus)

        n_cl = clustering.n_clusters

        #Getting the top n_top_clones clonotypes at each time point
        top_clones = self.top_clones_set(n_top_clones)
        #Building a trajectory dataframe
        traj_frame = self.build_traj_frame(top_clones)

        #Converting it in a numpy matrix
        traj_matrix = traj_frame.drop(['Clone cumul freq'], axis=1).to_numpy()

        # Normalize each trajectory by its maximum
        norm_traj_matrix = traj_matrix/np.max(traj_matrix, axis=1)[:, np.newaxis]

        fig, axs = plt.subplots(2, n_cl, figsize=(5*n_cl, 12))
        for cl in range(n_cl):
            trajs = norm_traj_matrix[clustering.labels_ == cl]
            axs[0][cl].set_xlabel('Time', fontsize = 15)
            axs[0][cl].set_ylabel('Normalized frequency', fontsize = 15)
            axs[1][cl].set_xlabel('Time', fontsize = 15)
            axs[1][cl].set_ylabel('Normalized frequency', fontsize = 15)
            for traj in trajs:
                axs[0][cl].plot(self.times, traj, alpha=0.2, color=cmap(cl/n_cl))
            axs[1][cl].set_ylim(0,1)
            axs[1][cl].errorbar(self.times, np.mean(trajs, axis=0), 
                                yerr=np.std(trajs, axis=0), lw=3, color=cmap(cl/n_cl))
            #axs[1][cl].fill_between(times, np.quantile(trajs, 0.75, axis=0), np.quantile(trajs, 0.25, axis=0), color=colors[cl])

        plt.tight_layout()
        return fig, axs

#===============================Data-Pre-Processing===================================

class Data_Process():

    """
    ## TODO in the future, merge this class with other classes.

    A class used to represent longitudinal RepSeq data and pre-analysis of the longitudinal data associated with
    one individual.

    ...

    Attributes
    ----------
    path : str
        the name of the path to get access to the data files to use for our analysis
    filename1 : str
        the name of the file of the RepSeq sample which can be the first replicate when deciphering the experimental noise
        or the first time point RepSeq sample when analysing responding clones to a stimulus between two time points.
    filename2 : str
        the name of the file of the RepSeq sample which can be the second replicate when deciphering the experimental noise
        or the second time point RepSeq sample when analysing responding clones to a stimulus between two time points.
    colnames1 : str
        list of columns names of data-set - first sample
    colnames2 : str
        list of columns names of data-set - second sample


    Methods
    -------

    import_data() :
        to import and merged two RepSeq samples and build a unique data-frame with frequencies and abundances of all TCR clones present in the
        union of both samples.


    """

    def __init__(self, path, filename1, filename2, colnames1,  colnames2):

        self.path = path
        self.filename1 = filename1
        self.filename2 = filename2
        self.colnames1 = colnames1
        self.colnames2 = colnames2


    def import_data(self):
        """
        TOFILL

        Parameters
        ----------
        NONE


        Returns
        -------
        number_clones
            numpy array, number of clones in the data frame which is the union of the two RepSeq used as entries of the function

        df
            pandas data-frame which is the data-frame containing the informations labeled in colnames vector string
            for both RepSeq samples taken as input.

        """

        mincount = 0
        maxcount = np.inf

        headerline=0 #line number of headerline
        newnames=['Clone_fraction','Clone_count','ntCDR3','AACDR3']

        if self.filename1[-2:] == 'gz':
            F1Frame_chunk=pd.read_csv(self.path + self.filename1, delimiter='\t',usecols=self.colnames1,header=headerline, compression = 'gzip')[self.colnames1]
        else:
            F1Frame_chunk=pd.read_csv(self.path + self.filename1, delimiter='\t',usecols=self.colnames1,header=headerline)[self.colnames1]

        if self.filename2[-2:] == 'gz':
            F2Frame_chunk=pd.read_csv(self.path + self.filename2, delimiter='\t',usecols=self.colnames2,header=headerline, compression = 'gzip')[self.colnames2]

        else:
            F2Frame_chunk=pd.read_csv(self.path + self.filename2, delimiter='\t',usecols=self.colnames2,header=headerline)[self.colnames2]

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
        mergedFrame.drop(newnames[3]+suffixes[1], axis = 1,inplace=True) #removes duplicates
        mergedFrame.rename(columns = {newnames[3]+suffixes[0]:newnames[3]}, inplace = True)
        mergedFrame=mergedFrame[[newname+suffix for newname in newnames[:2] for suffix in suffixes]+[newnames[2],newnames[3]]]
        filterout=((mergedFrame.Clone_count_1<mincount) & (mergedFrame.Clone_count_2==0)) | ((mergedFrame.Clone_count_2<mincount) & (mergedFrame.Clone_count_1==0)) #has effect only if mincount>0
        number_clones=len(mergedFrame)
        return number_clones,mergedFrame.loc[((mergedFrame.Clone_count_1<=maxcount) & (mergedFrame.Clone_count_2<=maxcount)) & ~filterout]


def get_sparserep(df):
    """
    Tranforms {(n1,n2)} data stored in pandas dataframe to a sparse 1D representation.
    unicountvals_1(2) are the unique values of n1(2).
    sparse_rep_counts gives the counts of unique pairs.
    ndn1(2) is the index of unicountvals_1(2) giving the value of n1(2) in that unique pair.
    len(indn1)=len(indn2)=len(sparse_rep_counts)


    Parameters
    ----------
    df : pandas data frame
        data-frame which is the output of the method .import_data() for one Data_Process instance.
        these data-frame should give the list of TCR clones present in two replicates RepSeq samples
        associated to their clone frequencies and clone abundances in the first and second replicate.


    Returns
    -------
    indn1
        numpy array list of indexes of all values of unicountvals_1

    indn2
        numpy array list of indexes of all values of unicountvals_2

    sparse_rep_counts
        numpy array, # of clones having the read counts pair {(n1,n2)}

    unicountvals_1
        numpy array list of unique counts values present in the first sample in df[clone_count_1]

    unicountvals_2
        numpy array list of unique counts values present in the second sample in df[clone_count_2]

    Nreads1
        float, total number of counts/reads in the first sample referred in df by "_1"

    Nreads2
        float, total number of counts/reads in the second sample referred in df by "_2"

    """
    counts = df.loc[:,['Clone_count_1', 'Clone_count_2']]
    counts['paircount'] = 1  # gives a weight of 1 to each observed clone

    clone_counts = counts.groupby(['Clone_count_1', 'Clone_count_2']).sum()
    sparse_rep_counts = np.asarray(clone_counts.values.flatten(), dtype=int)
    clonecountpair_vals = clone_counts.index.values
    indn1 = np.asarray([clonecountpair_vals[it][0] for it in range(len(sparse_rep_counts))], dtype=int)
    indn2 = np.asarray([clonecountpair_vals[it][1] for it in range(len(sparse_rep_counts))], dtype=int)
    nreads_1 = np.sum(counts['Clone_count_1'])
    nreads_2 = np.sum(counts['Clone_count_2'])

    unicountvals_1, indn1 = np.unique(indn1, return_inverse=True)
    unicountvals_2, indn2 = np.unique(indn2, return_inverse=True)

    return indn1, indn2, sparse_rep_counts, unicountvals_1, unicountvals_2, nreads_1, nreads_2

def _pois_log_likelihoods(x: ArrayLike,
                          mu: ArrayLike
                         ) -> Union[np.float64, np.ndarray]:
    log_likelihoods = -mu + x * np.log(mu) - gammaln(x + 1)
    log_likelihoods[np.isnan(log_likelihoods)] = 0
    return log_likelihoods

def _nbinom_log_likelihoods(x: ArrayLike,
                            r: ArrayLike,
                            log_p: ArrayLike,
                            log_1mp: ArrayLike
                           ) -> Union[np.float64, np.ndarray]:
    ln_binom = gammaln(x + r) - gammaln(r) - gammaln(x + 1)
    return ln_binom + x * log_p + r * log_1mp

def _get_rhof(alpha_rho: float,
              fmin: float,
              nfbins: int = 1200,
              freq_dtype: type = np.float64
             ) -> Tuple[np.ndarray]:
    """
    generates power law (power is alpha_rho) clone frequency distribution over
    freq_nbins discrete logarithmically spaced frequences between fmin and 1 of dtype freq_dtype
    Outputs log probabilities obtained at log frequencies
    """
    logfvec = np.linspace(fmin, 0, nfbins)
    logfvec = logfvec * np.log(10)
    d_logfvec_div_2 = np.diff(logfvec) / 2
    logrhovec = logfvec * alpha_rho
    integ = np.exp(logrhovec + logfvec, dtype=freq_dtype)
    normconst = np.log(np.dot(d_logfvec_div_2, integ[1:] + integ[:-1]))
    logrhovec -= normconst
    return logrhovec, logfvec, normconst, d_logfvec_div_2

def _log_pn_f_0(unicounts: np.ndarray,
                nreads: int,
                logfvec: np.ndarray,
                paras: np.ndarray
               ) -> np.ndarray:
    m_total = 10.**paras[3]
    r_c = nreads / m_total

    nsigma = 5.
    nmin = 300.
    m_low = np.zeros(len(unicounts), dtype=np.int64)
    m_high = np.zeros(len(unicounts), dtype=np.int64)

    mean_m = unicounts / r_c
    dev = nsigma * np.sqrt(mean_m)
    m_low = (mean_m - dev).astype(np.int64)
    m_high = (mean_m + 5 * dev).astype(np.int64)
    m_low[mean_m <= dev**2] = 0
    m_high[unicounts <= nmin] = 10 * nmin / r_c
    m_cellmax=np.max(m_high)
    mvec=np.arange(m_cellmax + 1)

    pois_probs = np.exp(_pois_log_likelihoods(unicounts, (mvec * r_c)[:, None]))
    # Remove rows which have probabilities below machine precision to avoid
    # unnecessary calculations.
    keep = np.any(pois_probs, axis=1)
    pois_probs = pois_probs[keep]

    beta_mv, alpha_mv = paras[1], paras[2]

    mean_m = m_total * np.exp(logfvec)
    log_mean = np.log(mean_m)
    log_mean_multiplier = np.log(beta_mv) + (alpha_mv - 1) * log_mean
    log_1mp = -np.logaddexp(0, log_mean_multiplier)
    log_p = log_mean_multiplier + log_1mp
    r = np.exp(log_mean - log_mean_multiplier)
    nbinom_probs = np.exp(_nbinom_log_likelihoods(mvec[keep], r[:, None],
                                                  log_p[:, None], log_1mp[:, None]))
    return np.log(nbinom_probs @ pois_probs)

def _log_pn_f_1(unicounts: np.ndarray,
                nreads: int,
                logfvec: np.ndarray,
                paras: np.ndarray
               ) -> np.ndarray:
    beta_mv, alpha_mv = paras[1], paras[2]
    mean_n = nreads * np.exp(logfvec)
    log_mean = np.log(mean_n)
    log_mean_multiplier = np.log(beta_mv) + (alpha_mv - 1) * log_mean
    log_1mp = -np.logaddexp(0, log_mean_multiplier)
    log_p = log_mean_multiplier + log_1mp
    r = np.exp(log_mean - log_mean_multiplier)
    return _nbinom_log_likelihoods(unicounts, r[:, None], log_p[:, None], log_1mp[:, None])

def _log_pn_f_2(unicounts: np.ndarray,
                nreads: int,
                logfvec: np.ndarray,
                *args
               ) -> np.ndarray:
    mean_n = nreads * np.exp(logfvec)
    return _pois_log_likelihoods(unicounts, mean_n[:, None])

def _get_pn1n2(paras: np.ndarray,
               sparse_rep: Tuple[np.ndarray, int],
               log_pn_f_func: Callable,
               nfbins: int = 1200,
               freq_dtype: type = np.float64
              ) -> np.ndarray:
    """
    Tool to compute likelihood of the noise model. It is not useful for the user.
    """
    alpha = paras[0]
    fmin = paras[-1]

    logrhofvec, logfvec, normconst, dlogfby2 = _get_rhof(alpha, fmin, nfbins, freq_dtype)
    logfvec_tmp = deepcopy(logfvec)

    (indn1, indn2, _, unicountvals_1, unicountvals_2,
     nreads_1, nreads_2) = sparse_rep

    log_pn1_f = log_pn_f_func(unicountvals_1, nreads_1, logfvec_tmp, paras)
    log_pn2_f = log_pn_f_func(unicountvals_2, nreads_2, logfvec_tmp, paras)

    # Compute P(0,0) for the normalization constraint
    integ = np.exp(logrhofvec + log_pn2_f[:, 0] + log_pn1_f[:, 0] + logfvec)
    pn0n0 = np.dot(dlogfby2, integ[1:] + integ[:-1])

    integ = np.exp(log_pn1_f[:, indn1] + log_pn2_f[:, indn2] + (logfvec + logrhofvec)[:, None])
    pn1n2 = np.dot(dlogfby2, integ[1:] + integ[:-1])

    pn1n2 /= 1. - pn0n0  # renormalize
    return pn1n2

def _avg_log_likelihood_pn1n2(paras: np.ndarray,
                              sparse_rep: Tuple[np.ndarray, int],
                              log_pn_f_func: Callable,
                              nfbins: int = 1200,
                              freq_dtype: Union[str, type] = np.float64
                             ) -> np.float64:
    _, _, sparse_rep_counts, _, _, _, _ = sparse_rep
    pn1n2 = _get_pn1n2(paras, sparse_rep, log_pn_f_func, nfbins, freq_dtype)
    log_pn1n2 = np.where(pn1n2 > 0, np.log(pn1n2), 0)
    return -np.dot(sparse_rep_counts, log_pn1n2) / np.sum(sparse_rep_counts)

#===============================Noise-Model===================================================================
### Noise Model

class NoiseModel():

    """
    A class used to build an object associated to methods in order to learn the experimental noise from same day
    biological RepSeq samples.

    ...

    Methods
    -------

    get_sparserep(df) :
        get sparse representation of the abundances / frequencies of the TCR clones present in both RepSeq samples of interest.
        this changes the data input to fasten the algorithm

    learn_null_model(df, noise_model, init_paras,  output_dir = None, filename = None, display_loss_function = False) :
        function to optimize the likelihood associated to the experimental noise model and get the associated parameters.

    diversity_estimate(df, paras, noise_model) :
        function to get the estimation of diversity from the noise model information.
    """
    def __init__(self,
                 num_frequency_bins: int = 1200,
                 freq_dtype: Union[str, type] = np.float64
                ) -> None:
        self.nfbins = num_frequency_bins
        self.freq_dtype = freq_dtype

    def get_sparserep(self,
                      df: pd.DataFrame
                     ) -> Tuple[np.ndarray, int]:
        return get_sparserep(df)

    def _callback(self,
                  paras: np.ndarray,
                  nparas: int,
                  sparse_rep: Tuple[np.ndarray, int],
                  noise_model: int,
                  log_pn_f_func: Callable
                 ) -> None:
        '''prints iteration info. called by scipy.minimize. Not useful for the user.'''

        global curr_iter
        #curr_iter = 0
        global Loss_function
        print(''.join(['{0:d} ']+['{'+str(it)+':3.6f} ' for it in range(1,len(paras)+1)]).format(*([curr_iter]+list(paras))))
        #print ('{' + str(len(paras)+1) + ':3.6f}'.format( [self.get_Pn1n2(paras, sparse_rep, acq_model_type)]))
        Loss_function = _avg_log_likelihood_pn1n2(paras, sparse_rep, log_pn_f_func, self.nfbins, self.freq_dtype)
        print(Loss_function)
        curr_iter += 1

    # Constraints for the Null-Model, no filtered 
    def _nullmodel_constr_fn(self,
                             paras: np.ndarray,
                             sparse_rep: Tuple[np.ndarray, int],
                             log_pn_f_func: Callable,
                             constr_type: int
                            ) -> Union[np.float64, Tuple[np.float64]]:
        """
        returns either or both of the two level-set functions: log<f>-log(1/N), with N=Nclones/(1-P(0,0)) and log(Z_f), with Z_f=N<f>_{n+n'=0} + sum_i^Nclones <f>_{f|n,n'}
        not useful for the user
        """
        indn1, indn2, sparse_rep_counts, unicountvals_1, unicountvals_2, nreads_1, nreads_2 = sparse_rep

        alpha = paras[0]  # power law exponent
        fmin = paras[-1] # true minimal frequency 

        logrhofvec, logfvec, normconst, dlogfby2 = _get_rhof(alpha, fmin, self.nfbins, self.freq_dtype)

        integ = np.exp(logrhofvec + 2 * logfvec)
        avgf_ps = np.dot(dlogfby2, integ[:-1] + integ[1:])

        log_pn1_f = log_pn_f_func(unicountvals_1, nreads_1, logfvec, paras)
        log_pn2_f = log_pn_f_func(unicountvals_2, nreads_2, logfvec, paras)

        integ = np.exp(log_pn1_f[:, 0] + log_pn2_f[:, 0] + logrhofvec + logfvec)
        Pn0n0 = np.dot(dlogfby2, integ[1:] + integ[:-1])
        logPnng0 = np.log(1 - Pn0n0)
        avgf_null_pair = np.exp(logPnng0 - np.log(np.sum(sparse_rep_counts)))

        C1 = np.log(avgf_ps) - np.log(avgf_null_pair)

        integ = np.exp(log_pn1_f[:, 0] + log_pn2_f[:, 0] + logrhofvec + 2 * logfvec)
        log_avgf_n0n0 = np.log(np.dot(dlogfby2, integ[1:] + integ[:-1]))

        log_integ = log_pn1_f[:, indn1] + log_pn2_f[:, indn2] + (logrhofvec + logfvec)[:, None]
        integ = np.exp(log_integ)
        log_Pn1n2 = np.log(np.einsum('i,ij->j', dlogfby2, integ[1:] + integ[:-1]))
        integ = np.exp(log_integ + logfvec[:, None])
        tmp = deepcopy(log_Pn1n2)
        tmp[tmp == -np.Inf] = np.inf  # since subtracted in next line
        avgf_n1n2 = np.exp(np.log(np.einsum('i,ij->j', dlogfby2, integ[1:] + integ[:-1])) - tmp)
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

    def learn_null_model(self,
                         df: pd.DataFrame,
                         noise_model: int,
                         init_paras: np.ndarray,
                         output_dir: str = None,
                         filename: str = None,
                         constr_type: int = 1,
                         tol: float = 1e-6,
                         maxiter: int = 200,
                        ) -> Tuple[OptimizeResult, float]:
        """
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame which is the output of the method .import_data() for one Data_Process instance.
            These data-frame should give the list of TCR clones present in two replicates RepSeq samples
            associated to their clone frequencies and clone abundances in the first and second replicate.
        noise_model: int
            Choice of noise model.
            Options are 0, 1, or 2.
        init_paras: numpy.ndarray
            Initial vector of parameters to start the optimization of the model.
        output_dir : str, optional
            The output directory name to which the values of the parameters will be saved.
        constr_type : int, default 1
            Specify which constraint to use.
            Constraint type 1 gives only low error modes, see paper for details.

        Returns
        -------
        outstruct : scipy.optimize.OptimizeResult
            The optimization result which has attributes ``x``, the parameters
            which give the solution, and ``fun``, the value of the solution.
        constr_value : float
            The value of the constraint at the solution.

        """
        sparse_rep = self.get_sparserep(df)

        if noise_model == 0:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'm_total', 'fmin']
            bounds = [(-np.inf, 0), (0, np.inf), (-np.inf, np.inf), (0, np.inf), (-np.inf, 0)]
            log_pn_f_func = _log_pn_f_0
        elif noise_model == 1:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'fmin']
            bounds = [(-np.inf, 0), (0, np.inf), (-np.inf, np.inf), (-np.inf, 0)]
            log_pn_f_func = _log_pn_f_1
        elif noise_model == 2:
            parameter_labels = ['alph_rho', 'fmin']
            bounds = [(-10, 0), (-20, 0)]
            log_pn_f_func = _log_pn_f_2
        else:
            raise ValueError('noise_model must be 0, 1, or 2.')

        assert len(parameter_labels) == len(init_paras), "number of model and initial paras differ!"

        condict = {'type': 'eq',
                   'fun': self._nullmodel_constr_fn,
                   'args': (sparse_rep, log_pn_f_func, constr_type)}

        partialobjfunc = partial(_avg_log_likelihood_pn1n2, sparse_rep=sparse_rep,
                                 log_pn_f_func=log_pn_f_func,
                                 nfbins=self.nfbins, freq_dtype=self.freq_dtype)

        header = ['Iter'] + parameter_labels
        print(''.join(['{' + str(it) + ':9s} ' for it in range(len(init_paras) + 1)]).format(*header))

        global curr_iter
        curr_iter = 1
        options = {'ftol': tol, 'disp': True, 'maxiter': maxiter}
        callbackp = partial(self._callback, nparas=len(init_paras),
                            sparse_rep=sparse_rep, noise_model=noise_model,
                            log_pn_f_func=log_pn_f_func)
        outstruct = minimize(partialobjfunc, init_paras, method='SLSQP',
                             callback=callbackp, constraints=condict,
                             options=options, bounds=bounds)

        constr_value = self._nullmodel_constr_fn(outstruct.x, sparse_rep,
                                                 log_pn_f_func, constr_type)

        if noise_model < 1:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'm_total', 'fmin']
            d = {'label' : parameter_labels, 'value': outstruct.x}
            df = pd.DataFrame(data = d)
        elif noise_model == 1:
            parameter_labels = ['alph_rho', 'beta', 'alpha', 'fmin']
            d = {'label' : parameter_labels, 'value': outstruct.x}
            df = pd.DataFrame(data = d)
        else:
            parameter_labels = ['alph_rho', 'fmin']
            d = {'label' : parameter_labels, 'value': outstruct.x}
            df = pd.DataFrame(data = d)

        if (output_dir == None) & (filename == None):
            df.to_csv('nullpara' + str(noise_model)+ '.txt', sep = '\t')

        elif (output_dir != None) & (filename == None):
            df.to_csv(output_dir + '/nullpara' + str(noise_model)+ '.txt', sep = '\t')

        else :
            df.to_csv(output_dir + '/' + filename + '.txt', sep = '\t')

        return outstruct, constr_value

    def diversity_estimate(self,
                           df: pd.DataFrame,
                           paras: np.ndarray,
                           noise_model: int
                          ) -> int:
        """
        Estimate diversity of the individual repertoire from the experimental noise learning step.

        Parameters
        ----------
        df : data-frame
            The data-frame which has been used to learn the noise model
        paras : numpy array
            vector containing the noise parameters
        noise_model : int
            choice of noise model

        Returns
        -------
        diversity_estimate
            float, diversity estimate from the noise model inference.

        """
        if noise_model == 0:
            log_pn_f_func = _log_pn_f_0
        elif noise_model == 1:
            log_pn_f_func = _log_pn_f_1
        elif noise_model == 2:
            log_pn_f_func = _log_pn_f_2
        else:
            raise ValueError('noise_model must be 0, 1, or 2.')

        sparse_rep = self.get_sparserep(df)

        _, _, sparse_rep_counts, _, _, nreads_1, nreads_2 = sparse_rep

        alpha = paras[0]
        fmin = paras[-1]

        logrhofvec, logfvec, _, dlogfby2 = _get_rhof(alpha, fmin, self.nfbins, self.freq_dtype)

        logfvec_tmp = deepcopy(logfvec)

        zero_arr = np.zeros(1, dtype=np.int64)
        log_pn1_f = log_pn_f_func(zero_arr, nreads_1, logfvec_tmp, paras).ravel()
        log_pn2_f = log_pn_f_func(zero_arr, nreads_2, logfvec_tmp, paras).ravel()

        # Compute P(0,0) for the normalization constraint
        integ = np.exp(logrhofvec + log_pn1_f + log_pn2_f + logfvec)
        Pn0n0 = np.dot(dlogfby2, integ[1:] + integ[:-1])

        N_obs = np.sum(sparse_rep_counts)

        return int(N_obs / (1 - Pn0n0))

#============================================Differential expression =============================================================

class Expansion_Model():

    """
    A class used to build an object associated to methods in order to select significant expanding or
    contracting clones from RepSeq samples taken at two different time points.

    ...

    Methods
    -------
    expansion_table(outpath, paras_1, paras_2, df, noise_model, pval_threshold, smed_threshold):
        generate the table of clones that have been significantly detected to be responsive to an acute stimuli.
    """
    def _get_Ps(self, alp,sbar,smax,stp):
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

    def _get_Ps_mesh(self, alpmesh, sbarmesh, smax, s_step):
        '''
        generates symmetric exponential distribution over log fold change
        with effect size sbar and nonresponding fraction 1-alp at s=0.
        computed over discrete range of s from -smax to smax in steps of size stp
        '''
        lamb = -s_step / sbarmesh
        smaxt = round(smax / s_step)
        s_zeroind = int(smaxt)
        Z = 2 * (np.exp((smaxt + 1) * lamb) - 1) / (np.exp(lamb) - 1) - 1
        Ps = alpmesh[:, :, None] * np.exp(lamb[:, :, None] * np.fabs(np.arange(-smaxt, smaxt + 1))) / Z[:, :, None]
        Ps[:, :, s_zeroind] += (1 - alpmesh)
        return Ps

    def _callbackFdiffexpr(self, Xi): #case dependent
        '''prints iteration info. called scipy.minimize'''

        print('{0: 3.6f}   {1: 3.6f}   '.format(Xi[0], Xi[1])+'\n')


    def _learning_dynamics_expansion_polished(self, df, paras_1, paras_2,  noise_model):
        """
        function to infer the expansion mode parameters - not usable by the user.
        """

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,nreads_1,nreads_2 = self.get_sparserep(df)

        alpha_rho = paras_1[0]
        fmin = np.power(10,paras_1[-1])
        freq_dtype = 'float64'
        nfbins = 1200 #Accuracy of the integration


        logrhofvec, logfvec, _, _ = _get_rhof(alpha_rho, fmin, nfbins, freq_dtype)

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
                Nreads = nreads_1
                paras = paras_1
            else:
                unicounts=unicountvals_2
                logfvec_tmp=deepcopy(logfvecwide) #contains s-shift for sampled data method
                Nreads = nreads_2
                paras = paras_2
            if it == 0:
                log_pn1_f = self._get_logPn_f( unicounts, Nreads, logfvec_tmp, noise_model, paras)

            else:
                log_pn2_f = self._get_logPn_f(unicounts, Nreads, logfvec_tmp, noise_model, paras)

        #for the trapezoid method
        dlogfby2=np.diff(logfvec)/2

        # Computing P(n1,n2|f,s)
        Pn1n2_s=np.zeros((len(svec), len(unicountvals_1), len(unicountvals_2)))

        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1,indn2)):
                integ = np.exp(logrhofvec+log_pn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+log_pn1_f[:,n1_it]+ logfvec )
                Pn1n2_s[s_it, n1_it, n2_it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        Pn0n0_s = np.zeros(svec.shape)
        for s_it,s in enumerate(svec):
            integ=np.exp(log_pn1_f[:,0]+log_pn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),0]+logrhofvec+logfvec)
            Pn0n0_s[s_it]=np.dot(dlogfby2,integ[1:]+integ[:-1])


        N_obs = np.sum(sparse_rep_counts)
        print("N_obs: " + str(N_obs))


        def cost(PARAS):

            alp = PARAS[0]
            sbar = PARAS[1]

            Ps = _get_Ps(self,alp,sbar,smax,s_step)
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


        outstruct = minimize(cost, initparas, method='SLSQP', callback=_callbackFdiffexpr, tol=1e-6,options={'ftol':1e-8 ,'disp': True,'maxiter':300})

        return outstruct.x, Pn1n2_s, Pn0n0_s, svec

    def _learning_dynamics_expansion(self, sparse_rep, paras_1, paras_2, noise_model, display_plot=False):
        """
        function to infer the expansion mode parameters - not usable by the user.
        """

        indn1,indn2,sparse_rep_counts,unicountvals_1,unicountvals_2,nreads_1,nreads_2 = sparse_rep

        alpha_rho = paras_1[0]
        fmin = np.power(10,paras_1[-1])
        freq_dtype = 'float64'
        nfbins = 1200 #Accuracy of the integration


        logrhofvec, logfvec = self._get_rhof(alpha_rho, nfbins, fmin, freq_dtype)

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
                Nreads = nreads_1
                paras = paras_1
            else:
                unicounts=unicountvals_2
                logfvec_tmp=deepcopy(logfvecwide) #contains s-shift for sampled data method
                Nreads = nreads_2
                paras = paras_2
            if it == 0:
                log_pn1_f = self._get_logPn_f(unicounts, Nreads, logfvec_tmp, noise_model, paras)

            else:
                log_pn2_f = self._get_logPn_f(unicounts, Nreads, logfvec_tmp, noise_model, paras)

        #for the trapezoid method
        dlogfby2=np.diff(logfvec)/2

        # Computing P(n1,n2|f,s)
        Pn1n2_s=np.zeros((len(svec), len(unicountvals_1), len(unicountvals_2)))

        for s_it,s in enumerate(svec):
            for it,(n1_it, n2_it) in enumerate(zip(indn1,indn2)):
                integ = np.exp(logrhofvec+log_pn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),n2_it]+log_pn1_f[:,n1_it]+ logfvec )
                Pn1n2_s[s_it, n1_it, n2_it] = np.dot(dlogfby2,integ[1:] + integ[:-1])


        Pn0n0_s = np.zeros(svec.shape)
        for s_it,s in enumerate(svec):
            integ=np.exp(log_pn1_f[:,0]+log_pn2_f[f2s_step*s_it:(f2s_step*s_it+len(logfvec)),0]+logrhofvec+logfvec)
            Pn0n0_s[s_it]=np.dot(dlogfby2,integ[1:]+integ[:-1])


        N_obs = np.sum(sparse_rep_counts)
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

        def create_landscape(alpha_mesh, sbar_mesh):
            Ps = self._get_Ps(alpha_mesh, sbar_mesh)
            Pn0n0 = np.einsum('ijk,k->ij', p_s, Pn0n0_s)
            Pn1n2_ps = np.einsum('ijk,klm', Ps, Pn1n2_s)
            Pn1n2_ps /= 1 - Pn0n0[:, :, None, None]
            log_Pn1n2_ps = np.where(Pn1n2_ps[:, :, indn1, indn2] > 0,
                                    np.log(Pn1n2_ps[:, :, indn1, indn2]),
                                    0)
            frequencies = sparse_rep_counts / float(N_obs)
            energy = np.einsum('ijk,k->ij', log_Pn1n2_ps, frequencies)
            return energy

    #--------------------------Compute-the-grid-----------------------------------------

        print('Calculation Surface : \n')
#        st = time.time()

        npoints = 50 #to be chosen by the user 
        alpvec = np.logspace(-3,np.log10(0.99), npoints)
        sbarvec = np.linspace(0.01,5, npoints)

#        LSurface =np.zeros((len(sbarvec),len(alpvec)))
#        for i in range(len(sbarvec)):
#            for j in range(len(alpvec)):
#                LSurface[i, j]=  - cost([alpvec[j], sbarvec[i]])
#
#        alpmesh, sbarmesh = np.meshgrid(alpvec, sbarvec)
#        a,b = np.where(LSurface == np.max(LSurface))
#        print("--- %s seconds ---" % (time.time() - st))

        st = time.time()
        alpmesh, sbarmesh = np.meshgrid(alpvec, sbarvec)
        LSurface = create_landscape(alpvec, sbarvec)
        print(f'--- {time.time() - st} seconds ---')

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


    def _save_table(self, outpath, svec, Ps,Pn1n2_s, Pn0n0_s,  subset, unicountvals_1_d, unicountvals_2_d, indn1_d, indn2_d, print_expanded, pthresh, smedthresh):
        '''
        takes learned diffexpr model, Pn1n2_s*Ps, computes posteriors over (n1,n2) pairs, and writes to file a table of data with clones as rows and columns as measures of thier posteriors
        print_expanded=True orders table as ascending by , else descending
        pthresh is the threshold in 'p-value'-like (null hypo) probability, 1-P(s>0|n1_i,n2_i), where i is the row (i.e. the clone) n.b. lower null prob implies larger probability of expansion
        smedthresh is the threshold on the posterior median, below which clones are discarded

        not usable by the user.
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

    def expansion_table(self, outpath, paras_1, paras_2, df, noise_model, pval_threshold, smed_threshold):

        '''
        generate the table of clones that have been significantly detected to be responsive to an acute stimuli.

        Parameters
        ----------
        outpath  : str
            Name of the directory where to store the output table
        paras_1  : numpy array
            parameters of the noise model that has been learned at time_1
        paras_2  : numpy array
            parameters of the noise model that has been learned at time_2
        df       : pandas dataframe
            pandas dataframe merging the two RepSeq data at time_1 and time_2

        noise_model : int
            choice of noise model 0: Poisson, 1: negative Binomial, 2: negative Binomial + Poisson

        pval_threshold : float
            P-value threshold to detect and discriminate if a TCR clone has expanded

        smed_threshold : float
            median of the log-fold change threshold to detect if a TCR clone has expanded

        Returns
        -------
        data-frame - csv file
            the output is a csv file of columns : $s_{1-low}$, $s_{2-med}$, $s_{3-high}$, $s_{max}$, $\bar{s}$, $f_1$, $f_2$, $n_1$, $n_2$, 'CDR3_nt', 'CDR3_AA' and '$p$-value'
        '''

        sparse_rep = self.get_sparserep(df)
        L_surface, Pn1n2_s_d, Pn0n0_s_d, svec = self._learning_dynamics_expansion(sparse_rep, paras_1, paras_2, noise_model)
        npoints= 50 # same as in learning_dynamics_expansion
        smax = 25.0
        s_step = 0.1
        alpvec = np.logspace(-3,np.log10(0.99), npoints)
        sbarvec = np.linspace(0.01,5, npoints)
        maxinds=np.unravel_index(np.argmax(L_surface),np.shape(L_surface))
        optsbar=sbarvec[maxinds[0]]
        optalp=alpvec[maxinds[1]]
        optPs= self._get_Ps(optalp,optsbar,smax,s_step)
        pval_expanded = True

        indn1,indn2,sparse_rep_counts, unicountvals_1, unicountvals_2, nreads_1, nreads_2 = sparse_rep

        self._save_table(outpath, svec, optPs, Pn1n2_s_d, Pn0n0_s_d,  df, unicountvals_1, unicountvals_2, indn1, indn2, pval_expanded, pval_threshold, smed_threshold)


#============================================Generate Synthetic Data =============================================================

class Generator:

    """
    A class used to build an object to generate in-Silico (synthetic) RepSeq samples, in the case of replicates at
    the same day and in the case of having 2 samples generated at an initial time for the first one and some time after (months, years)
    for the second one using the geometric Brownian motion model decribed in https://www.biorxiv.org/content/10.1101/2022.05.01.490247v1.

    ...

    Methods
    -------

    gen_synthetic_data_Null(paras, noise_model, nreads_1,nreads_2,Nsamp):
        generate in-silico same day RepSeq replicates.

    generate_trajectories(tau, theta, method, paras_1, paras_2, t_ime, filename, nreads_1 = '1e6', nreads_2 = '1e6'):
        generate in-silico t_ime apart RepSeq samples.
    """

    def _get_rhof(self, alpha_rho, fmin, freq_nbins=800, freq_dtype='float64'):

        '''
        generates power law (power is alpha_rho) clone frequency distribution over
        freq_nbins discrete logarithmically spaced frequences between fmin and 1 of dtype freq_dtype
        Outputs log probabilities obtained at log frequencies'''
        fmax=1e0
        logfvec=np.linspace(np.log10(fmin),np.log10(fmax),freq_nbins)
        logfvec=np.array(np.log(np.power(10,logfvec)) ,dtype=freq_dtype).flatten()
        logrhovec=logfvec*alpha_rho
        integ=np.exp(logrhovec+logfvec,dtype=freq_dtype)
        normconst=np.log(np.dot(np.diff(logfvec)/2.,integ[1:]+integ[:-1]))
        logrhovec-=normconst
        return logrhovec,logfvec

    def _get_distsample(self, pmf,Nsamp,dtype='uint32'):
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


    def gen_synthetic_data_Null(self, paras, noise_model, nreads_1,nreads_2,Nsamp):
        '''
        outputs an array of observed clone frequencies and corresponding dataframe of pair counts
        for a null model learned from a dataset pair with nreads_1 and nreads_2 number of reads, respectively.
        Crucial for RAM efficiency, sampling is conditioned on being observed in each of the three (n,0), (0,n'), and n,n'>0 conditions
        so that only Nsamp clones need to be sampled, rather than the N clones in the repertoire.
        Note that no explicit normalization is applied. It is assumed that the values in paras are consistent with N<f>=1 
        (e.g. were obtained through the learning done in this package).
        '''


        alpha = paras[0] #power law exponent
        fmin=np.power(10,paras[-1])
        if noise_model<1:
            m_total=float(np.power(10, paras[3]))
            r_c1=nreads_1/m_total
            r_c2=nreads_2/m_total
            r_cvec=[r_c1,r_c2]
        if noise_model<2:
            beta_mv= paras[1]
            alpha_mv=paras[2]

        logrhofvec,logfvec = self._get_rhof(alpha,fmin)
        fvec=np.exp(logfvec)
        dlogf=np.diff(logfvec)/2.

        #generate measurement model distribution, Pn_f
        Pn_f=np.empty((len(logfvec),),dtype=object) #len(logfvec) samplers

        #get value at n=0 to use for conditioning on n>0 (and get full Pn_f here if noise_model=1,2)
        m_max=1e3 #conditioned on n=0, so no edge effects

        Nreadsvec=(nreads_1,nreads_2)
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
        f_samples_inds= self._get_distsample(dlogf*(integ[1:] + integ[:-1]),num_qx0).flatten()
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
        f_samples_inds=self._get_distsample(dlogf*(integ[1:] + integ[:-1]),num_q0x).flatten()
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
        f_samples_inds=self._get_distsample(dlogf*(integ[1:] + integ[:-1]),num_qxx).flatten()
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


    def generate_trajectories(self, tau, theta, method, paras_1, paras_2, t_ime, filename, nreads_1 = '1e6', nreads_2 = '1e6'):


        """
        generate in-silico t_ime apart RepSeq samples.

        Parameters
        ----------
        paras_1  : numpy array
            parameters of the noise model that has been learnt at time_1
        paras_2  : numpy array
            parameters of the noise model that has been learnt at time_2
        method   : str
            'negative_binomial' or 'poisson'
        tau      : float
            first time-scale parameter of the dynamics
        theta    : float
            second time-scale parameter of the dynamics
        t_ime    : float
            number of years between both synthetic sampling (between time_1 and time_2)
        filename : str
            name of the file in which the dataframe is stored

        Returns
        -------
        data-frame - csv file
            the output is a csv file of columns : 'Clone_count_1' (at time_1) 'Clone_count_2' (at time_2) and the frequency counterparts 'Clone_frequency_1' and 'Clone_frequency_2'
        """

        np.seterr(divide = 'ignore')
        np.warnings.filterwarnings('ignore')

        method = 'negative_binomial'


        # Synthetic data generation

        print('execution starting...')

        st = time.time()

        #Values of the parameters
        A = -1/tau
        B = 1/theta
        N_0 = 40
        nreads_1 = float(NreadsI)
        nreads_2 = float(nreads_1I)

        t = float(t_ime)

        if nreads_1 == nreads_2:
            key_sym = '_sym_'

        else:
            key_sym = '_asym_'

        # Name of the directory


        dirName = 'output'
        os.makedirs(dirName, exist_ok=True)

        paras = paras_1 #Just put a and b of the negative binomiale distribution [0.7, 1.1]
        alpha = -1 +2*A/B
        #print('alpha : ' + str(alpha))

        #1/ Generate log-population at initial time from steady-state distribution + GBM diffusion trajectories for 2 years
        x_i_LB, x_f_LB, Prop_Matrix_LB, p_ext_LB, results_extinction_LB, time_vec_LB, results_extinction_source_LB, x_source_LB = _generator_diffusion_LB(A, B, N_0, t)

        #x_i_LB, x_f_LB, Prop_Matrix, p_ext, results_extinction  = generator_diffusion_LB(B, A, N_0, t)
        N_cells_day_0_LB, N_cells_day_1_LB = np.sum(np.exp(x_i_LB)), np.sum(np.exp(x_f_LB)) + np.sum(np.exp(x_source_LB))  #N_cells_final_LB
        print('NUMBER OF CELLS AT INITIAL TIME')
        print(N_cells_day_0_LB)

        print('NUMBER OF CELLS AT FINAL TIME')
        print(N_cells_day_1_LB)

        #print('SHAPE_X_I ' +  str(np.shape(x_i_LB)))
        #print('SHAPE_X_F ' +  str(np.shape(x_f_LB)))


        if method == 'negative_binomial':

            df_diffusion_LB  = _experimental_sampling_diffusion_NegBin(nreads_1, nreads_2, paras, x_i_LB, x_f_LB, N_cells_day_0_LB, N_cells_day_1_LB)
            df_diffusion_LB.to_csv(filename + '.csv' , sep= '\t')

        elif method == 'poisson':

            df_diffusion_LB  = _experimental_sampling_diffusion_Poisson(nreads_1, nreads_2, x_i_LB, x_f_LB, t, N_cells_day_0_LB, N_cells_day_1_LB)
            df_diffusion_LB.to_csv(filename + '.csv' , sep= '\t')
