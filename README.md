# Bipartite Partial Configuration Model for Python

The Bipartite Partial Configuration Model (BiPCM) is a statistical null model
for binary bipartite networks. It offers an unbiased method of
analyzing node similarities and obtaining statistically validated monopartite
projections \[Saracco2016\].

The BiPCM is related to the [Bipartite Configuration Model
(BiCM)](https://github.com/tsakim/bicm) \[Saracco2015\], but imposes only
constraints on the degrees of one bipartite node layer. It belongs to a series
of entropy-based null model for binary biparite networks, see also

* [BiCM](https://github.com/tsakim/bicm)
* [BiRG](https://github.com/tsakim/birg)

Please consult the original articles for details about the underlying methods
and applications to user-movie and international trade databases
\[Saracco2016, Straka2016\].
 
## Author 

Mika J. Straka

## Version and Documentation

The newest version of the module can be found on
[https://github.com/tsakim/bipcm](https://github.com/tsakim/bipcm).

The complete documentation is available at
[http://bipcm.readthedocs.io/](http://bipcm.readthedocs.io/) and in the file
`docs/BiPCM_manual.pdf`
<<<<<<< HEAD

## How to cite

If you use the `bipcm` module, please cite its location on Github
[https://github.com/tsakim/bipcm](https://github.com/tsakim/bipcm) and the
original article \[Saracco2016\]. 

=======

## How to cite

If you use the `bipcm` module, please cite its location on Github
[https://github.com/tsakim/bipcm](https://github.com/tsakim/bipcm) and the
original article \[Saracco2016\]. 

>>>>>>> d18f8bc617dea777e57d4d9d8b306ce6e96d5150
### References

\[Saracco2015\] [F. Saracco, R. Di Clemente, A. Gabrielli, T. Squartini, Randomizing bipartite networks: the case of the World Trade Web, Scientific Reports 5, 10595 (2015)](http://www.nature.com/articles/srep10595).

\[Saracco2016\] [F. Saracco, M. J. Straka, R. Di Clemente, A. Gabrielli, G. Caldarelli, T. Squartini, Inferring monopartite projections of bipartite networks: an entropy-based approach, arXiv preprint arXiv:1607.02481](https://arxiv.org/abs/1607.02481)

\[Straka2016\] [M. J. Straka, F. Saracco, G. Caldarelli, Product Similarities in International Trade from Entropy-based Null Models, Complex Networks 2016, 130-132 (11 2016), ISBN 978-2-9557050-1-8](http://www.complexnetworks.org/BookOfAbstractCNA16.pdf)

<!---
It is also possible to calculate the probabilites of two nodes (i, j) having 0, 1, ..., M nearest neighbors in common, where M is the number of nodes in the opposite bipartite layer. To do so for the row-nodes, execute 
```python
cma.save_lambda_probdist(bip_set=True, write=True, filename=<filename>, delim='\t')
```
For the column-nodes, use `bip_set=False`.The default name of the ouput file is `pbicm_lambdaprob_contr_<constraint>_layer_<bip_set>.csv`.

Note that, if `<constraint> == <bip_set>`, the output file contains N(N - 1) rows for the node pairs (1, 2), (1, 3), ..., (1, N), (2, 3), ..., (N - 1, N), and M + 1 columns contain the probabilities of having 0, 1, ..., M common neighbors.
If `<constraint> != <bip_set>`, the file contains only one row containing the probabilities for two nodes (i, j) sharing 0, 1, ..., M + 1 common neighbors, since the probability distribution is the same for all node pairs in the bipartite nodes set `bip_set`.  

-->

---
Copyright (c) 2015-2017 Mika J. Straka 
