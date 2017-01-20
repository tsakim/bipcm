BiPCM Quickstart
================================================================================

The calculation of the p-values of node similarities with the ``bipcm`` module
is straightforward as shown below. The validated node similarities can be used
to obtain an unbiased monopartite projection of the bipartite network, as
illustrated in [Saracco2016]_.

For more detailed explanations of the methods, please refer to [Saracco2016]_,
the :ref:`tutorial` and the :ref:`api`.

Calculating the p-values of the node similarities
--------------------------------------------------------------------------------

Be ``mat`` a two-dimensional binary NumPy array, which describes the
`biadjacency matrix
<https://en.wikipedia.org/w/index.php?title=Adjacency_matrix&oldid=751840428#Adjacency_matrix_of_a_bipartite_graph>`_
of an undirected bipartite network. The nodes of the two bipartite layers are
ordered along the columns and rows, respectively. In the algorithm, the two
layers are identified by the boolean values ``True`` for the **row-nodes** and
``False`` for the **column-nodes**.

Import the module and initialize the Bipartite Partial Configuration Model:: 

    >>> from src.bipcm import BiPCM
    >>> pcm = BiPCM(bin_mat=mat, constraint=<bool>)

The parameter ``constraint`` specifies whether the degrees of the
row-nodes (``constraint = True``) or the degrees of the column-nodes
(``constraint = False``) should be constrained.  

In order to analyze the similarity of the row-layer nodes and to save the
p-values of the corresponding :math:`\Lambda`-motifs, i.e. of the number of
shared neighbors [Saracco2016]_, use::

    >>> pcm.lambda_motifs_main(bip_set=True, filename=<filename>)

For the column-layer nodes, use::

    >>> pcm.lambda_motifs_main(bip_set=False, filename=<filename>)

``bip_set`` selects the bipartite node set for which the p-values should be
calculated and saved. The filename *<filename>* should contain a relative path
declaration. The default name of the output file is
*pval_constr_<constraint>_proj_<bip_set>.csv*, where *<constraint>* and
*<bip_set>* are either *rows* or *columns* depending on the degree constraint
and the parameter choice in ``lambda_motifs_main``. By default, the values in
the file are separated by tabs, which can be changed using the ``delim``
keyword. 

Subsequently, the p-values can be used to perform a multiple hypotheses testing
and to obtain statistically validated monopartite projections [Saracco2016]_.

If the p-values should not be saved but returned by 
``lambda_motifs_main``, use::

    >>> pcm.lambda_motifs_main(bip_set=True, write=False)

.. By default, the file is saved in a human-readable CSV format. The information can also be saved as a binary NumPy file ``.npy`` by using::

..    >>> cm.save_matrix(cm.adj_matrix, filename=<filename>, binary=True)

