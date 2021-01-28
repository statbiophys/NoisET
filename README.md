# NoisET  NOIse sampling learning & Expansion detection of T-cell receptors using Bayesian inference.
(NoisET should be pronounced like "noisettes" (aka nuts in French)).
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
`--freq 'frequency'` : Column label associated to the clonal fraction \
`--counts 'counts'`:  Column label associated to the clonal count  \
`--ntCDR3 'ntCDR3'`:  Column label associated to the clonal CDR3 nucleotides sequence  \
`--AACDR3 'AACDR3'`:  Column label associated to the clonal CDR3 amino acid sequence


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

 At the command prompt, type 
 ```$ noiset-nullgenerator --NB --NreadsI 25000 --NreadsII 25000 --Nclones 20000 --filename 'test' --nullpara '../data_examples/parameters_1.npy'  
 ```
 
 

To detect responding clones to a stimulus: NoisET second function (2)
At the command prompt, type 
```$ noiset-detection --NB --freq 'Clone fraction' --counts 'Clone count' --ntCDR3 'N. Seq CDR3' --AACDR3 'AA. Seq. CDR3' --path '../data_examples/' --f1 'Q1_0_F1_.txt' --f2 'Q1_15_F1_.txt' --nullpara1 'parameters_1.npy' --nullpara2 'parameters_1.npy' --pval 0.05 --smedthresh 0  `


## First function : Noise Model

</br>

### Infer the model
First function of NoisET is to infer statistical null model of sequence counts and variability, using RepSeq experiments.
In the main code file noisettes.py, this function is associated to the class Experimental_Noise_Model. In the notebook [Null Model Learning],
replicates of patient $S_1$ from [Pogorelyy PNASref] are used to extract sampling noise dispersion. An other example using replicates from another technology () is also displayed in the notebook.

!! default input : initial parameters, maxcount mincount

To learn the null model, the user should choose one of  the fourth noise model : noise_model 
noise_model:
- 0 : NB + Poisson / 5 parameters
- 1 : NB / 4 parameters
- 2 : Poisson / 2 parameters

### Check the model



## Second function: Differential Expression: 

To detect responding clonotypes, the user provides, in addtion to the two datasets to be compared, two sets of experimental noise parameters learnt at both times
using the first function. 





