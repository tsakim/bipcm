.. _overview:

Overview
================================================================================

The ``bipcm`` module is an implementation of the Bipartite Partial
Configuration Model (BiPCM) as described in the article [Saracco2016]_. The
BiPCM can be used as a statistical null model to analyze the similarity of
nodes in undirected bipartite networks. The similarity criterion is based on
the number of common neighbors of nodes, which is expressed in terms of
:math:`\Lambda`-motifs in the original article [Saracco2016]_. Subsequently,
one can obtain unbiased statistically validated monopartite projections of the
original bipartite network.

The construction of the BiPCM, just like the related `BiCM
<https://github.com/tsakim/bipcm>`_ and `BiRG
<https://github.com/tsakim/birg>`_ models, is based on the generation of a
grandcanonical ensemble of bipartite graphs subject to certain constraints. The
constraints can be of different types. For instance, in the Bipartite Random
Graph (BiRG) the total number of edges is fixed. In the case of the BiCM
the average degrees of all the nodes are constrained. In the BiPCM, on the
other hand, only the degrees of one bipartite layer are constrained. 

The average graph of the ensemble can be calculated analytically using the
entropy-maximization principle and provides a statistical null model, which can
be used for establishing statistically significant node similarities. In
general, they are referred to as entropy-based null models. For more
information and a detailed explanation of the underlying methods, please refer
to [Saracco2016]_.  

By using the ``bipcm`` module, the user can obtain the BiPCM null model which
corresponds to the input matrix representing an undirected bipartite network.
To address the question of node similarity, the p-values of the observed
numbers of common neighbors can be calculated and used for statistical
verification. For an illustration and further details, please refer to
[Saracco2016]_ and [Straka2016]_.

Dependencies
--------------------------------------------------------------------------------

``bipcm`` is written in `Python 2.7` and uses the following modules:

* `poibin <https://github.com/tsakim/poibin>`_ Module for the Poisson Binomial
  probability distribution 
* `scipy <https://www.scipy.org/>`_
* `numpy <http://www.numpy.org>`_
* `doctest <https://docs.python.org/2/library/doctest.html>`_ For unit testing

