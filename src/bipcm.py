# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 10:00:21 2016

Module:
    bipcm - Bipartite Partial Configuration Model

Author:
    Mika Straka

Description:
    Implementation of the Bipartite Partial Configuration Model with one-sided
    constraint (BiPCM) for binary undirected bipartite networks [Saracco2016]_.

    Given the biadjacency matrix of a bipartite graph in the form of a binary
    array as input, the module allows the user to calculate the biadjacency
    matrix of the ensemble average graph :math:`<G>^*` of the BiPCM null model,
    where the matrix entries correspond to the link probabilities
    :math:`<G>^*_{rc} = p_{rc}` between nodes of the two distinct bipartite
    node sets. Subsequently, one can calculate the p-values of the node
    similarities for nodes in the same bipartite layer [Saracco2016]_.
    The BiPCM derives from the Bipartite Configuration Model as presented in
    [Saracco2015]_, in which the degrees of both bipartite layers are fixed
    on average. For the BiPCM, the constraints are relaxed and only the degrees
    of one bipartite layer are imposed [Saracco2016]_.

Usage:
    Be ``mat`` a two-dimensional binary NumPy array. The nodes of the two
    bipartite layers are ordered along the columns and rows, respectively. In
    the algorithm, the two layers are identified by the boolean values ``True``
    for the **row-nodes** and ``False`` for the **column-nodes**.

    Import the module and initialize the Bipartite Partial Configuration Model::

        >>> from src.bipcm import BiPCM
        >>> cma = BiPCM(bin_mat=mat, constraint=<constraint>)

    where ``<constraint> = True`` constrains the degrees of the row-nodes and
    ``<constraint> = False`` the degrees of the column nodes.

    In order to analyze the similarity of the row-layer nodes and to save the
    p-values of the corresponding :math:`\\Lambda`-motifs, use::

        >>> cma.lambda_motifs_main(bip_set=True, filename=<filename>)

    For the column-layer nodes, use::

        >>> cma.lambda_motifs_main(bip_set=False, filename=<filename>)

    ``bip_set`` selects the bipartite node set for which the p-values should be
    calculated and saved. The filename *<filename>* should contain a relative
    path declaration. The default name of the output file is
    *pval_constr_<constraint>_proj_<bip_set>.csv*, where *<contraint>* and
    *<bip_set>* are either *rows* or *columns* depending on the degree
    constraints and the parameter choice in ``lambda_motifs_main``. By default,
    the values in the file are separated by tabs, which can be changed using the
    ``delim`` keyword.

Reference:
    [Saracco2015] F. Saracco, R. Di Clemente, A. Gabrielli, T. Squartini,
    Randomizing bipartite networks: the case of the World Trade Web, Scientific
    Reports 5, 10595 (2015)

    [Saracco2016] F. Saracco, M. J. Straka, R. Di Clemente, A. Gabrielli, G.
    Caldarelli, T. Squartini, Inferring monopartite projections of bipartite
    networks: an entropy-based approach, arXiv preprint arXiv:1607.02481
"""

import numpy as np
from scipy.stats import binom
from poibin.poibin import PoiBin


class BiPCM:
    """Bipartite Partial Configuration Model for binary bipartite newtworks

    This class implements the Bipartite Configuration Model (BiCM), which can
    be used as a null model for the analysis of undirected and binary bipartite
    networks. The class provides methods to calculate the biadjacency matrix of
    the null model and to quantify node similarities in terms of p-values.
    """

    def __init__(self, bin_mat, constraint):
        """Initialize the parameters of the BiPCM.

        :param bin_mat: binary input matrix describing the biadjacency matrix
                of a bipartite graph with the nodes of one layer along the rows
                and the nodes of the other layer along the columns.
        :type bin_mat: numpy.array
        :param constraint: constraints either the degrees of the row-nodes
            (``True``), or the degrees of the column-nodes columns (``False``)
        :type constraint: bool
        """
        self.bin_mat = np.array(bin_mat)
        self.check_constraint(constraint)
        self.const_set = constraint
        self.check_input_matrix_is_binary()
        [self.num_rows, self.num_columns] = self.bin_mat.shape
        self.eprob_seq = self.get_edge_prob_seq()

# ------------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------------

    def check_input_matrix_is_binary(self):
        """Check that the entries of the input matrix are 0 or 1."""
        assert np.all(np.logical_or(self.bin_mat == 0, self.bin_mat == 1)), \
            "Input matrix is not binary."

    @staticmethod
    def check_constraint(constraint):
        """Check that the constraint parameter is either ``True`` of ``False``.

        :param constraint: constraint on the degrees of either the row-nodes
            (``True``) or column-nodes (``False``)
        :type constraint: bool
        """
        assert constraint in [False, True], \
            "Constraint has to be True or False."

    def set_degree_seq(self):
        """Set the degree sequence [degrees row-nodes, degrees column-nodes]."""
        dseq = np.empty(self.num_rows + self.num_columns)
        dseq[self.num_rows:] = np.squeeze(np.sum(self.bin_mat, axis=0))
        dseq[:self.num_rows] = np.squeeze(np.sum(self.bin_mat, axis=1))
        assert dseq.size == (self.num_rows + self.num_columns)
        return dseq

    def get_edge_prob_seq(self):
        """Return an array with the link probabilities of the BiPCM null model.

        In the first part of the array, the row degrees are fixed. In the
        second part, the column degrees are fixed.

        :returns: array of link probabilities
        :rtype: numpy.array
        """
        dseq = self.set_degree_seq()
        eprob_seq = np.zeros(dseq.size)
        eprob_seq[:self.num_rows] = dseq[:self.num_rows] / \
                                         self.num_columns
        eprob_seq[self.num_rows:] = dseq[self.num_rows:] / \
                                         self.num_rows
        return eprob_seq

# ------------------------------------------------------------------------------
# Total log-likelihood of the observed Lambda motifs in the input matrix
# ------------------------------------------------------------------------------

    def lambda_loglike(self, bip_set=False):
        """Return the log-likelihood of the number of :math:`\\Lambda`-motifs.

        The total log-likelihood of the number of observed
        :math:`\\Lambda`-motifs in the input matrix is calculated according to
        the BiPCM null model.

        :param bip_set: analyze :math:`\\Lambda`-motifs of row-nodes (``True``)
            or column-nodes (``False``)
        :type bip_set: bool
        """
        plam_mat = self.get_plambda_matrix()
        nlam_mat = self.get_lambda_motif_matrix(self.bin_mat, bip_set)
        p_mat = self.get_proj_pmat(plam_mat, nlam_mat, bip_set)
        logp = np.log(p_mat[np.triu_indices_from(p_mat, k=1)])
        loglike = logp.sum()
        print loglike
        return loglike

    def get_proj_pmat(self, plam_mat, nlam_mat, bip_set=False):
        """Return the probabilities of the observed :math:`\\Lambda`-motifs.

        The probabilities of the :math:`\\Lambda`-motifs between the nodes
        specified by ``bip_set`` in the input matrix are calculated and
        returned.

        If the node set ``bip_set`` is the same as the  constrained one, the
        :math:`\\Lambda`-motifs follow a Binomial probability distribution.
        Otherwise, all the node pairs follow the same Poisson Binomial
        distribution.

        The probability mass function is given by :math:`pmf(k) = Pr(X = k)`.

        .. note::
            The lower triangular part including the diagonal is set to 0 since
            the matrix is symmetric.

        :param plam_mat: matrix of Lambda motif probabilities
        :type plam_mat: numpy.array
        :param nlam_mat: matrix of observed number of Lambda motifs
        :type nlam_mat: numpy.array
        :param bip_set: select row-nodes (``True``) or column-nodes (``False``)
        :type bip_set: bool
        :returns: matrix containing the probabilities of the
            :math:`\\Lambda`-motifs
        :rtype: numpy.array
        """
        if bip_set:
            m = self.num_columns
        elif not bip_set:
            m = self.num_rows
        else:
            errmsg = "'" + str(bip_set) + "' " + 'not supported.'
            raise NameError(errmsg)
        n = nlam_mat.shape[0]
        pmat = np.zeros(nlam_mat.shape)
        if bip_set != self.const_set:
            pb = PoiBin(plam_mat[np.diag_indices_from(plam_mat)])
            for i in xrange(n):
                pmat[i, i + 1:] = pb.pmf(nlam_mat[i, i + 1:])
        elif bip_set == self.const_set:
            # if the sets correspond, the matrix dimensions should be the same
            assert plam_mat.shape[0] == nlam_mat.shape[0]
            for i in xrange(n):
                for j in xrange(i + 1, n):
                    bn = binom(m, plam_mat[i, j])
                    pmat[i, j] = bn.pmf(nlam_mat[i, j])
        return pmat

# ------------------------------------------------------------------------------
# Lambda motifs
# ------------------------------------------------------------------------------

    def lambda_motifs_main(self, bip_set=False, write=True, filename=None,
                           delim='\t'):
        """Calculate and save the p-values of the :math:`\\Lambda`-motifs.

        For each node couple in the bipartite layer specified by ``bip_set``,
        :math:`\\Lambda`-motifs and calculate the corresponding p-value.

        :param bip_set: select row-nodes (``True``) or column-nodes (``False``)
        :type bip_set: bool
        :param write: if ``True``, the pvalues are saved in the specified file
        :type write: bool
        :param filename: name of the file which will contain the p-values,
            default is *pval_constr_<contraint>_proj_<rows OR columns>.csv*
        :type filename: str
        :param delim: delimiter between entries in file, default is tab
        :type delim: str
        """
        plam_mat = self.get_plambda_matrix()
        nlam_mat = self.get_lambda_motif_matrix(self.bin_mat, bip_set)
        pval_mat = self.get_lambda_pvalues(plam_mat, nlam_mat, bip_set)
        if write:
            if filename is None:
                if bip_set:
                    b = 'rows'
                elif not bip_set:
                    b = 'columns'
                else:
                    errmsg = "'" + str(bip_set) + "' " + 'not supported.'
                    raise NameError(errmsg)
                if self.const_set:
                    constr = "rows"
                else:
                    constr = "columns"
                fname = 'pval_constr_' + constr + '_proj_' + b + '.csv'
            else:
                fname = filename
            self.save_matrix(pval_mat, filename=fname, delim=delim)
        else:
            return pval_mat

    def get_plambda_matrix(self):
        """Return the :math:`\\Lambda`-motif probability matrix.

        Return a square matrix of Lambda probabilities for the nodes given
        the degree constraints on the node set self.const_set.

        .. note::
            The lower triangular part excluding the diagonal is set to 0 since
            the matrix is symmetric.
        """
        if self.const_set:
            ps = self.eprob_seq[:self.num_rows]
        elif not self.const_set:
            ps = self.eprob_seq[self.num_rows:]
        mat = np.zeros((ps.size, ps.size))
        for i in xrange(ps.size):
            mat[i, i:] = ps[i] * ps[i:]
        return mat

    def get_lambda_motif_matrix(self, mm, bip_set=False):
        """Return the number of :math:`\\Lambda`-motifs as found in ``mm``.

        Given the binary input matrix ``mm``, count the number of
        :math:`\\Lambda`-motifs between node couples of the bipartite layer
        specified by ``bip_set``.

        :param mm: binary matrix
        :type mm: numpy.array
        :param bip_set: selects row-nodes (``True``) or column-nodes (``False``)
        :type bip_set: bool
        :returns: square matrix of observed :math:`\\Lambda`-motifs
        :rtype: numpy.array
        """
        if bip_set:
            l2_mat = np.dot(mm, np.transpose(mm))
            assert l2_mat.shape == (self.num_rows, self.num_rows)
        elif not bip_set:
            l2_mat = np.dot(np.transpose(mm), mm)
            assert l2_mat.shape == (self.num_columns, self.num_columns)
        else:
            errmsg = str(bip_set) + 'not supported.'
            raise NameError(errmsg)
        # set diagonal to zero
        di = np.diag_indices_from(l2_mat)
        l2_mat[di] = 0
        return l2_mat.astype(int)

    def get_lambda_pvalues(self, plam_mat, nlam_mat, bip_set=False):
        """Return the p-values for the :math:`\\Lambda`-motifs in ``nlam_mat``.

        Calculate the p-values for the numbers of observed
        :math:`\\Lambda`-motifs as given in the parameter ``nlam_mat`` for the
        bipartite node layer ``bip_set``. The probabilities for the single
        :math:`\\Lambda`-motifs are given in ``plam_mat``.

        If ``bip_set`` corresponds to the constrained bipartite node set, the
        :math:`\\Lambda`-motifs follow a Binomial probability distribution.
        Otherwise, all the node pairs follow the same Poisson Binomial
        probability distribution.

        .. note::

            :math:`pval(k) = Pr(X >= k) = 1 - Pr(X < k) = 1 - cdf(k) + pmf(k)`

            The lower triangular part (including the diagonal) of the returned
            matrixis set to zero.

        :param plam_mat: matrix of :math:`\\Lambda`-motif probabilities
        :type plam_mat: numpy.array
        :param nlam_mat: matrix of observed number of Lambda motifs
        :type nlam_mat: numpy.array
        :param bip_set: selects row-nodes (``True``) or column-nodes (``False``)
        :type bip_set: bool
        :returns: matrix of the p-values for the :math:`\\Lambda`-motifs
        :rtype: numpy.array
        """
        if bip_set:
            m = self.num_columns
        elif not bip_set:
            m = self.num_rows
        else:
            errmsg = "'" + str(bip_set) + "' " + 'not supported.'
            raise NameError(errmsg)

        n = nlam_mat.shape[0]
        pval_mat = np.zeros(nlam_mat.shape)

        if bip_set != self.const_set:
            pb = PoiBin(plam_mat[np.diag_indices_from(plam_mat)])
            for i in xrange(n):
                pval_mat[i, i + 1:] = pb.pval(nlam_mat[i, i + 1:])
        elif bip_set == self.const_set:
            # if the sets correspond, the matrix dimensions should be the same
            assert plam_mat.shape[0] == nlam_mat.shape[0]
            for i in xrange(n):
                for j in xrange(i + 1, n):
                    bn = binom(m, plam_mat[i, j])
                    pval_mat[i, j] = 1. - bn.cdf(nlam_mat[i, j]) \
                                     + bn.pmf(nlam_mat[i, j])
        return pval_mat

# ------------------------------------------------------------------------------
# Probability distributions for Lambda values
# ------------------------------------------------------------------------------

#    def save_lambda_probdist(self, bip_set=False, write=True, filename=None,
#                             delim='\t', binary=True):
#        """Obtain and save the probabilities of all possible values of the
#        Lambda motifs for the node set defined by bip_set. The matrix can
#        either be saved as human-readable ASCII or as a binary Numpy file.
#
#        :param bip_set: analyze row-nodes (True) or column-nodes (False)
#        :type bip_set: bool
#        :param write: if True, the pvalues are saved in an external file
#        :type write: bool
#        :param filename: name of the output file, default name is
#            bipcm_pval_constr_<contraint>_proj_<rows OR columns>.csv
#        :param delim: delimiter to use if file is saved as .csv
#        :param binary: if true, save as binary .npy file. Otherwise as .csv
#                        file
#        """
#        plam_mat = self.get_plambda_matrix()
#        lambdaprobs_mat = self.get_lambda_probdist(plam_mat, bip_set)
#        if write:
#            if filename is None:
#                if bip_set:
#                    b = 'rows'
#                elif not bip_set:
#                    b = 'columns'
#                else:
#                    errmsg = "'" + str(bip_set) + "' " + 'not supported.'
#                    raise NameError(errmsg)
#                if self.const_set:
#                    constr = "rows"
#                else:
#                    constr = "columns"
#                fname = 'bipcm_lambdaprob_constr_' + constr + '_layer_' + b
#            else:
#                fname = filename
#            if binary:
#                fname += '.npy'
#            else:
#                fname += '.csv'
#            self.save_matrix(lambdaprobs_mat, filename=fname, delim=delim,
#                             binary=binary)
#        else:
#            return lambdaprobs_mat
#
#    def get_lambda_probdist(self, plam_mat, bip_set=False):
#        """Return a matrix which contains the probabilities for each node pair
#        (i, j) in the bipartite node set bip_set of observing every possible
#        number of common neighbors in the opposite bipartite layer.
#
#        :param plam_mat: matrix of Lambda motif probabilities
#        :type plam_mat: numpy.array
#        :param bip_set: selects row-nodes (True) or column-nodes (False)
#        :type bip_set: bool
#        """
#        if bip_set:
#            n = self.num_rows
#            m = self.num_columns
#        elif not bip_set:
#            n = self.num_columns
#            m = self.num_rows
#        else:
#            errmsg = "'" + str(bip_set) + "' " + 'not supported.'
#            raise NameError(errmsg)
#
#        lambda_values = np.arange(m + 1)
#
#        if bip_set != self.const_set:
#            lambdaprobs_mat = np.zeros([1, m + 1])
#            pb = PoiBin(plam_mat[np.diag_indices_from(plam_mat)])
#            lambdaprobs_mat[0, :] = pb.pmf(lambda_values)
#        elif bip_set == self.const_set:
#            lambdaprobs_mat = np.zeros([n * (n - 1) / 2, m + 1])
#            # if the sets correspond, the matrix dimensions should be the same
#            for i in xrange(n):
#                for j in xrange(i + 1, n):
#                    bn = binom(m, plam_mat[i, j])
#                    k = self.triumat2flat_idx(i, j, n)
#                    lambdaprobs_mat[k, :] = bn.pmf(lambda_values)
#        else:
#            errmsg = "'" + str(bip_set) + "' " + 'not supported.'
#            raise NameError(errmsg)
#        return lambdaprobs_mat

# ------------------------------------------------------------------------------
# Auxiliary methods
# ------------------------------------------------------------------------------

    @staticmethod
    def triumat2flat_idx(i, j, n):
        """Convert an matrix index couple to a flattened array index.

        Given a square matrix of dimension :math:`n` and an index couple
        :math:`(i, j)` *of the upper triangular part* of the matrix, the
        function returns the index which the matrix element would have in a
        flattened array.

        .. note::
            * :math:`i \\in [0, ..., n - 1]`
            * :math:`j \\in [i + 1, ..., n - 1]`
            * returned index :math:`\\in [0,\\, n (n - 1) / 2 - 1]`

        :param i: row index
        :type i: int
        :param j: column index
        :type j: int
        :param n: dimension of the square matrix
        :type n: int
        :returns: flattened array index
        :rtype: int
        """
        return int((i + 1) * n - (i + 2) * (i + 1) / 2. - (n - (j + 1)) - 1)

#    @staticmethod
#    def triumat2flat_idx(idx_i, idx_j, n):
#        """Convert index couple (i, j) into index in one-dimensional array index
#        for the upper triangular part of a square matrix with dimension n. I.e.
#        idx_i runs from 0, ..., n and idx_j runs from idx_i + 1, ..., n
#
#        NB: the returned indices start from 0.
#        """
#        return int((idx_i + 1) * n - (idx_i + 2) * (idx_i + 1) / 2.
#                   - (n - (idx_j + 1)) - 1)

#    @staticmethod
#    def get_main_dir(main_dir_name='bipcm'):
#        """Return the absolute path to the main directory which contains the
#        folders "src" and "output".
#        Note that the default directory name is "bipcm".
#
#        :param main_dir_name: name of the main directory of the program.
#        :type main_dir_name: string
#        """
#        s = os.getcwd()
#        dirpath = s[:s.index(main_dir_name) + len(main_dir_name) + 1]
#        return dirpath

#    def save_matrix(self, mat, filename, delim='\t', binary=False):
#        """Save the input matrix in a csv-file.
#
#        :param mat: two-dimensional matrix
#        :param filename: name of the output file
#        :param delim: delimiter between values in file.
#        :param binary: if true, save as binary .npy file. Otherwise as .csv
#                        file
#        """
#        fname = ''.join([self.main_dir, '/output/', filename])
#        if binary:
#            np.save(fname, mat)
#        else:
#            np.savetxt(fname, mat, delimiter=delim)

    @staticmethod
    def save_matrix(mat, filename, delim='\t', binary=False):
        """Save the matrix ``mat`` in the file ``filename``.

        The matrix can either be saved as a binary NumPy ``.npy`` file or as a
        human-readable CSV file.

        .. note:: The relative path has to be provided in the filname, e.g.
                *../data/pvalue_matrix.csv*

        :param mat: two-dimensional matrix
        :type mat: numpy.array
        :param filename: name of the output file
        :type filename: str
        :param delim: delimiter between values in file
        :type delim: str
        :param binary: if ``True``, save as binary ``.npy``, otherwise as a
            CSV file
        :type binary: bool
        """
        if binary:
            np.save(filename, mat)
        else:
            np.savetxt(filename, mat, delimiter=delim)

################################################################################
# Main
################################################################################

if __name__ == "__main__":
    pass
