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

## Installation

To install NoisET, git clone the file in your working environment. 
Using terminal, go inside NoisET directory and write the following command : 

```console
sudo python setup.py install
```

Once done, you will be able to import NoisET as introduced in the NoisET examples files.
```python 
import noisets
```
NoisET can be installed using python pip. (not yet allowed!)

```console
pip install noisets
```

## Command Lines:

To Infer Null Model noise: NoisET first function (1)
At the command prompt, type ` (env) machine: user$ noiset-noise --NB --freq 'Clone fraction' --counts 'Clone count' --ntCDR3 'N. Seq CDR3' --AACDR3 'AA. Seq. CDR3' --path 'data_examples/' --f1 'Q1_0_F1_.txt' --f2 'Q1_0_F2_.txt' `

To generate biological replicates - to check consistency with experiments 
At the command prompt, type ` (env) machine: user$ noiset-nullgenerator --NB --NreadsI 25000 --NreadsII 25000 --Nclones 20000 --filename 'test' --nullpara 'data_examples/parameters_1.npy'  `

To detect responding clones to a stimulus: NoisET second function (2)
At the command prompt, type ` (env) machine: user$ noiset-detection --NB --freq 'Clone fraction' --counts 'Clone count' --ntCDR3 'N. Seq CDR3' --AACDR3 'AA. Seq. CDR3' --path 'data_examples/' --f1 'Q1_0_F1_.txt' --f2 'Q1_15_F1_.txt' --nullpara1 'parameters_1.npy' --nullpara2 'parameters_1.npy' --pval 0.05 --smedthresh 0  `


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





