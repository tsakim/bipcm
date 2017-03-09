Bipartite Partial Configuration Model - Documentation
================================================================================

The Bipartite Partial Configuration Model (BiPCM) is a statistical null model
for binary bipartite networks. It offers an unbiased method for  
analyzing node similarities and obtaining statistically validated monopartite
projections [Saracco2016]_.

The BiPCM is related to the `Bipartite Configuration Model (BiCM)
<https://github.com/tsakim/bicm>`_ [Saracco2015]_, but imposes only constraints
on the degrees of one bipartite node layer. It belongs to a series of
entropy-based null models for binary biparite networks, see also

* `BiCM <https://github.com/tsakim/bicm>`_ - Bipartite Configuration Model
* `BiRG <https://github.com/tsakim/birg>`_ - Bipartite Random Graph

Please consult the original articles for details about the underlying methods
and applications to user-movie and international trade databases
[Saracco2016]_, [Straka2016]_.

An example case is illustrated in the :ref:`tutorial`.

How to cite
--------------------------------------------------------------------------------

If you use the ``bipcm`` module, please cite its `location on Github
<https://github.com/tsakim/bipcm>`_ and the original article [Saracco2016]_. 

References
````````````````````````````````````````````````````````````````````````````````

.. [Saracco2015] `F. Saracco, R. Di Clemente, A. Gabrielli, T. Squartini, Randomizing bipartite networks: the case of the World Trade Web, Scientific Reports 5, 10595 (2015) <http://www.nature.com/articles/srep10595>`_

.. [Saracco2016] `F. Saracco, M. J. Straka, R. Di Clemente, A. Gabrielli, G. Caldarelli, T. Squartini, Inferring monopartite projections of bipartite networks: an entropy-based approach, arXiv preprint arXiv:1607.02481 <https://arxiv.org/abs/1607.02481>`_

.. [Straka2016] `M. J. Straka, F. Saracco, G. Caldarelli, Product Similarities in International Trade from Entropy-based Null Models, Complex Networks 2016, 130-132 (11 2016), ISBN 978-2-9557050-1-8 <http://www.complexnetworks.org/BookOfAbstractCNA16.pdf>`_

Getting Started
================================================================================

.. toctree::
   :maxdepth: 2

   ./source/overview
   ./source/quickstart
   ./source/tutorial
   ./source/testing
   ./source/src
   ./source/license
   ./source/contact

Indices and tables
================================================================================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

