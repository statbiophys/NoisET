# NoisET  NOIse sampling learning & Expansion detection of T-cell receptors using Bayesian inference.
(NoisET should be pronounced like "noisettes" (aka nuts in French)).
NoisET is an easy-to-use Python package allowing to assess sampling noise [ref] and detect responding Tcells clones to a given stimulus using Bayesian inference 
and RepSeq longitudinal data.
Statistical Model used in this package are described in <https://arxiv.org/pdf/1912.08304.pdf>.
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





