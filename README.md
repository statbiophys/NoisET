# NoisEts

<https://arxiv.org/pdf/1912.08304.pdf>

----------------------------------------------------------------------------------------------------------------------------

## Installation

To install Noisets, git clone the file ijn your working environment. 
Using terminal, go inside the NoisEts directory and write the following command : 

```console
sudo python setup.py install
```

Once done, you will be able to import noisets as introduced in the noisets examples files.

NoisEttes can be installed using python pip.

```console
pip install NoisEttes
```

## Noise Model:

</br>

(Put a filter + Visualization )

Not the same class of functions:

Choose the model choosing a keyword:

One class for RNA sequences Noise Models:
- fourth cases. 
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


```python 
import NoisEttes
```

