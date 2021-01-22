# NoisET  NOIse sampling learning & Expansion detection of T-cell receptors using Bayesian inference.
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

Once done, you will be able to import noisets as introduced in the noisets examples files.
```python 
import noisets
```
NoisET can be installed using python pip. (not yet allowed!)

```console
pip install noisets
```

## Noise Model:

</br>

(Put a filter + Visualization )

Not the same class of functions:

Choose the model choosing a keyword:

One class for RNA sequences Noise Models:
- four cases. 
- output: the learning of the parameters of the model.
- 1-step: Poisson (counting) / Negative Binomial (counting)
- 2-steps : NegativeBinomial + Poisson / Poisson + NegativeBinomial 


One class for gDNA sequences (motivation Harlan Robins Adaptive Data - look at the current noise of the data for these replicates)
- limiting branching Process 


## Differential Expression: 

- Max’s propagator  (choosing one of the null model + implementing the parameters that have been learnt before)
- Diffusion propagator (A, B, fluctuating fitness per clone) / check the frequency dependency.) 
- Birth-and-Death propagator;

Application:
Clonal Expansion:
- look at the differences between these 3 methods to be able to discriminate responding clones / also comparison with edgeR

TCR dynamics as a short time scale for the bulk repertoire:
- Find a way to learn the likelihood estimator for a trajectory that is not Markovian (propagator of a trajectory) for the 3 dynamics model.
-  Same for a long-time scale. (In the article put other parameters CD4/CD8, naive/memory).



