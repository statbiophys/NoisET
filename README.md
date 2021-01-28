# NoisET  NOIse sampling learning & Expansion detection of T-cell receptors using Bayesian inference.
(NoisET should be pronounced like "noisettes" (ie nuts in French)).
High-throughput sequencing of T- and B-cell receptors makes it possible to track immune
repertoires across time, in different tissues, in acute and chronic diseases or in healthy individuals. However
quantitative comparison between repertoires is confounded by variability in the read count of each receptor
clonotype due to sampling, library preparation, and expression noise. We present an easy-to-use python
package NoisET that implements and generalizes a previously developed Bayesian method in [Puelma Touzel et al, 2020](<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007873&rev=2>). It can be used
to learn experimental noise models for repertoire sequencing from replicates, and to detect responding
clones following a stimulus. The package was tested on different repertoire sequencing technologies and
datasets.

NoisET package is desribed in  <https://arxiv.org/pdf/1912.08304.pdf> #change the arxiv file to NoisET file

----------------------------------------------------------------------------------------------------------------------------

# Installation

Python 3 

To install NoisET, gitclone the file in your working environment. 
Using terminal, go inside NoisET directory and write the following command : 

```console
$ sudo python setup.py install
```
NoisET can be installed using python pip. (not yet allowed!)

```console
$ pip install noisets
```

# Documentation

## 1/ Infer Noise Model 

### A/ Command Line

To Infer Null Model noise: NoisET first function (1), use the command `noiset-noise`
Several options are needed to learn noise model from two replicates samples associated to one individual at a specific time point:

#### 1/ Data information:

`--path 'PATHTODATA'`: set path to data file \
`--f1 'FILENAME1_X_REP1'`: filename for individual X replicate 1 \
`--f2 'FILENAME2_X_REP2'`: filename for individual X replicate 2 

If your TCRs CDR3 clonal populations features (ie clonal fractions, clonal counts, clonal nucleotides CDR3 sequences and clonal amino acids sequences) have different column names than (respectively): ('Clone fraction', 'Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3), you can specify it by using: 

`--specify` \
`--freq 'frequency'` : Column label associated to clonal fraction \
`--counts 'counts'`:  Column label associated to clonal count  \
`--ntCDR3 'ntCDR3'`:  Column label associated to clonal CDR3 nucleotides sequence  \
`--AACDR3 'AACDR3'`:  Column label associated to clonal CDR3 amino acid sequence


#### 2/ Choice of noise model: (described in Methods section)
`--NBPoisson`: Negative Binomial + Poisson Noise Model - 5 parameters \
`--NB`: Negative Binomial - 4 parameters  \
`--Poisson`: Poisson - 2 parameters 

#### 3/ Example:
At the command prompt, type:
```console
$ noiset-noise --path '../data_examples/' --f1 'Q1_0_F1_.txt' --f2 'Q1_0_F2_.txt' --NB
```
This command line will learn four parameters associated to Negative Binomial Noise Model `--NB` for individual Q1 of the associated study at day 0.
A .npy file is created in the working directory: it is a 5/4/2 parameters vector depending of the choise of NBP/NB/Poisson noise model.
You can run previous example using data provided in data_examples folder - data from [Precise tracking of vaccine-responding T cell clones reveals convergent and personalized response in identical twins, Pogorelyy et al, PNAS](https://www.pnas.org/content/115/50/12704) 

### B/ Python Package 

For Python users, it is possible to use NoisET as a package importing it as mentioned before. A jupyter notebook explaining the use of all the functions of interest is provided: NoisET example - Null model learning.ipynb
```python 
import noisets
from noisets import noisettes as ns
```
You can download Jupyter notebook and modify it with your own PATHTODATA / datafile specificities.

## 2/ Generate synthetic data for null model learning:

To check qualitatively consistency of NoisET first function (1) with experiments or for other reasons, it can be useful to generates synthetic replicates from null model (described in Methods section).

### A/ Command Line

#### 1/ Choice of noise model:
As before `--NBP`, `--NB`, or `--Poisson`.

#### 2/ Specify learnt parameters:
`--nullpara 'PATH/FOLDER/NULLPARAS.npy'`: parameters learnt thanks to NoisET function (1) \
!!! Watch out to match correctly the noise model and the null parameters file content : 5 parameters for `--NBP`, 4 parameters for `--NB`and 2 parameters
for `--Poisson`. 

#### 3/ Sequencing properties of data:
`--NreadsI NNNN`: total number  of reads in first replicate - it should match the actual data. \
`--Nreads2 NNNN`: total number  of reads in second replicate - it should match the actual data. \
`--Nclones NNNN`: total number of clones in union of two replicates - it should match the actual data.

### 4/ Output file
`--filename 'SYNTHETICDATA'`: name of the output file where you can find the synthetic data set. 

At the command prompt, type 
 ```console
 $ noiset-nullgenerator --NB --nullpara '../data_examples/parameters_1.npy' --NreadsI 25000 --NreadsII 25000 --Nclones 20000 --filename 'test'  
 ```
 ### B/ Python Package 

For Python users, it is possible to use NoisET as a package importing it as mentioned before. A jupyter notebook explaining the use of all the functions of interest is provided: NoisET example - Null model learning.ipynb
```python 
import noisets
from noisets import noisettes as ns
```
You can download Jupyter notebook and modify it with your own PATHTODATA / datafile specificities - vizualization tools are also provided.


 ## 3/ Detect Responding clones:
 
To detect responding clones to a stimulus: NoisET second function (2)

### A/ Command Line

#### 1/ Choice of noise model:
As before `--NBP`, `--NB`, or `--Poisson`.

#### 2/ Specify learnt parameters for both time points:
(they can be the same for both time points if replicates are not available but to use carefully as mentioned in [ARTICLE]) \
`--nullpara1 'PATH/FOLDER/NULLPARAS1.npy'`: parameters learnt thanks to NoisET function (1) for time 1 \
`--nullpara2 'PATH/FOLDER/NULLPARAS2.npy'`: parameters learnt thanks to NoisET function (1) for time 2  

!!! Watch out to match correctly the noise model and the null parameters file content : 5 parameters for `--NBP`, 4 parameters for `--NB`and 2 parameters
for `--Poisson`. 

#### 3/ Data information:

`--path 'PATHTODATA'`: set path to data file \
`--f1 'FILENAME1_X_time1'`: filename for individual X time 1 \
`--f2 'FILENAME2_X_time2'`: filename for individual X time 2 

If your TCRs CDR3 clonal populations features (ie clonal fractions, clonal counts, clonal nucleotides CDR3 sequences and clonal amino acids sequences) have different column names than (respectively): ('Clone fraction', 'Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3), you can specify it by using: 

`--specify` \
`--freq 'frequency'` : Column label associated to clonal fraction \
`--counts 'counts'`:  Column label associated to clonal count  \
`--ntCDR3 'ntCDR3'`:  Column label associated to clonal CDR3 nucleotides sequence  \
`--AACDR3 'AACDR3'`:  Column label associated to clonal CDR3 amino acid sequence

#### 4/ Detection thresholds: (More details in Methods section).
`--pval XXX` : p-value threshold for the expansion/contraction - use 0.05 as default value. \
`--smedthresh XXX` : log fold change mediane threshold for the expansion/contraction - use 0 as default value. 


At the command prompt, type 
```console
$ noiset-detection --NB  --nullpara1 '../data_examples/parameters_1.npy' --nullpara2 '../data_examples/parameters_1.npy' --path '../data_examples/' --f1 'Q1_0_F1_.txt' --f2 'Q1_15_F1_.txt' --pval 0.05 --smedthresh 0 
```
 ### B/ Python Package 

For Python users, it is possible to use NoisET as a package importing it as mentioned before. A jupyter notebook explaining the use of all the functions of interest is provided: NoisET example - detection responding clones.ipynb
```python 
import noisets
from noisets import noisettes as ns
```
You can download Jupyter notebook and modify it with your own PATHTODATA / datafile specificities - vizualization tools are also provided.

# Methods
