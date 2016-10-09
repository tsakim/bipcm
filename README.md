# Partial Bipartite Configuration Model with one-sided constraints for Python

## About
The module contains a Python implementation of the partial Bipartite Configuration
Model with one-sided degree constraints (pBiCM), which can be used as a
statistical null model for undirected and binary bipartite networks (see
reference \[1\]).

Given the biadjacency matrix of a bipartite network, the probabilities of links
between nodes of different bipartite network layers are calculated for the
corresponding ensemble average graph. 

To address the question of node similarity within one bipartite layer, it is
possible to perform a statistical validation of the number of common nearest
neighbors and to calculate the p-values of the corresponding Lambda-motifs, as
described by Saracco et al. \[1\].
 
## Author 
Mika Straka

### This Version
The newest version can be found on
[https://github.com/tsakim/bicmi](https://github.com/tsakim/bicmi).

## Dependencies
* [ctypes](https://docs.python.org/2/library/ctypes.html)
* [multiprocessing](https://docs.python.org/2/library/multiprocessing.html)
* [scipy](https://www.scipy.org/)
* [numpy](www.numpy.org)
* [doctest](https://docs.python.org/2/library/doctest.html)
* [poibin](https://github.com/tsakim/poibin) Module for the Poisson Binomial
  probability distribution 

## Usage
Be `input_mat` a two-dimensional binary numpy array describing the biadjacency
matrix of an undirected bipartite network. The nodes of the two distinct
bipartite layers are ordered along the columns and rows, respectively. In the
algorithm, the two layers are identified by the boolean values `True` for the
row-nodes and `False` for the column-nodes.

Import the module
```python
from src.bicmi import BiCMa
```
and initialize the partial Bipartite Configuration Model with one sided
constraint for the matrix `input_mat` with 
```python
cma = BiCMa(bin_mat=td, constraint=<constraint>)
```
where `<constraint> == True` constrains the degrees of the row-nodes and
`<constraint> == False` the degrees of the column nodes.

In order to analyze the similarity of the row-layer nodes and to save the
p-values of the corresponding Lambda-motifs in the folder `bicmi/output/`,
use
```python
cma.lambda_motifs_main(bip_set=True, write=True, filename=<filename>, delim='\t') ```
```
For the column-layer nodes, use
```python
cma.lambda_motifs_main(bip_set=False, write=True, filename=<filename>, delim='\t')
```
`bip_set` defines the bipartite node set for which the p-values should be
saved. The default name of the ouput file is
`bicmi_pval_constr_<constraint>proj_<bip_set>.csv`

### NB: Main folder
Note that saving the files requires the name of the main directory,
which contains the folder `src` and thus the file `src/bicmi.py`.
If the name of the main directory is *not* the default `bicmi`, the BiCMa
instance has to be initialized as 
```python
cma = BiCMa(bin_mat=td, constraint=<constraint>, main_dir=<main directory name>)
```

## Testing
The methods have been implemented using the doctest module. To run the tests,
execute 
```
$ python -m doctest bicmi_tests.txt
```
in the folder `src` in the command line. If you want to run the tests in
verbose mode, use 
```
$ python -m doctest -v bicmi_tests.txt
```
Note that `bicmi.py` and `bicmi_tests.txt` have to be in the same directory.

## References
* \[1\] [Saracco, Di Clemente, Gabrielli, Squartini - Randomizing bipartite networks:
the case of the World Trade Web](http://www.nature.com/articles/srep10595)
---
Copyright (c) 2015-2016 Mika J. Straka 
