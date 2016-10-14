# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 10:00:21 2016

Module:
    bipcm - Bipartite Configuration Model with one-sided constraints

Author:
    Mika Straka

Description:
    Implementation of the Bipartite Partial Configuration Model with one-sided
    constraint (BiPCM).

    Given the biadjacency matrix of a bipartite graph in the form of a binary
    array as input, the module calculates the link probabilities for the
    corresponding ensemble average graph <G>^*, where the degrees of one
    bipartite node set are fixed and used as constraints. The user can
    calculate and save the p-values of the Lambda motifs of the two distinct
    bipartite node sets, which correspond to the row and column indices of the
    biadjacency matrix.

Usage:
    Be <td> a binary matrix in the form of an numpy array. The nodes of the two
    distinct bipartite layers are ordered along the columns and rows. In the
    following, rows and columns are referred to by the booleans "True" and
    "False" and the nodes are called "row-nodes" and "column-nodes",
    respectively.

    Initialize the Bipartite Partial Configuration Model with one sided
    constraint for the matrix <td> with

        $ cma = BiPCM(bin_mat=td, constraint=<constraint>)

    where <constraint> == True constrains the degrees of the row-nodes and
    <constraint> == False the degrees of the column nodes.

    In order to analyze the similarity of the row-layer nodes and to save the
    p-values of the corresponding Lambda-motifs in the folder "bipcm/output/",
    use

        $ cma.lambda_motifs_main(bip_set=True, write=True, filename=<filename>, 
                                delim='\t')

    For the column-layer nodes, use

        $ cma.lambda_motifs_main(bip_set=False, write=True, filename=<filename>,
                                delim='\t')

   "bip_set" defines the bipartite node set for which the p-values should be
   saved. The default name of the ouput file is
        'bipcm_pval_constr_<constraint>proj_<bip_set>.csv'

NB Main folder
    Note that saving the files requires the name of the main directory
    which contains the folder "src" and itself contains the file bicm.py.
    If the folder name is NOT the default "bipcm", the BiCM instance has to be
    initialized as

        $ cma = BiPCM(bin_mat=td, constraint=<constraint>,
                        main_dir=<main directory name>)
"""

import os
import numpy as np
from scipy.stats import binom
from poibin.poibin import PoiBin


class BiPCM:
    """Create the Bipartite Configuration model with one-sided constraint for
    the input matrix and analyze the Lambda motifs.
    """

    def __init__(self, bin_mat, constraint, main_dir='bipcm'):
        """Initialize the parameters of the BiPCM.

        :param bin_mat: binary input matrix describing the biadjacency matrix
                of a bipartite graph with the nodes of one layer along the rows
                and the nodes of the other layer along the columns.
        :type bin_mat: np.array
        :param constraint: constraints either the degrees of the row-nodes
            (True), or the degrees of the column-nodes columns (False)
        :type constraint: bool
        :param main_dir: directory containing the src/ and output/ folders
        """
        self.bin_mat = np.array(bin_mat)
        self.check_constraint(constraint)
        self.const_set = constraint
        self.check_input_matrix_is_binary()
        [self.num_rows, self.num_columns] = self.bin_mat.shape
        self.eprob_seq = self.get_edge_prob_seq()
        self.main_dir = self.get_main_dir()

# ------------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------------

    def check_input_matrix_is_binary(self):
        """Check that the input matrix is binary, i.e. entries are all either
        0 or 1.
        """
        assert np.all(np.logical_or(self.bin_mat == 0, self.bin_mat == 1)), \
            "Input matrix is not binary."

    def check_constraint(self, constr):
        """Check that the constraint is either True of False."""
        assert constr in [False, True], \
            "Constraint has to be True or False."

    def set_degree_seq(self):
        """Set the degree sequence [degrees row-nodes, degrees column-nodes]."""
        dseq = np.empty(self.num_rows + self.num_columns)
        dseq[self.num_rows:] = np.squeeze(np.sum(self.bin_mat, axis=0))
        dseq[:self.num_rows] = np.squeeze(np.sum(self.bin_mat, axis=1))
        assert dseq.size == (self.num_rows + self.num_columns)
        return dseq

    def get_edge_prob_seq(self):
        """Return an array with the link probabilities. In the first part of the
        array, the row degrees are fixed. In the second part, the column
        degrees are fixed.

        :return eprob_seq: array of link probabilities
        :rtype eprob_seq: np.array
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
        """Calculate the total log-likelihood of the number of observed Lambda
        motifs in the input matrix according to the null model.

        :param bip_set: analyze row-nodes (True) or column-nodes (False)
        :type bip_set: bool
        :param write: if True, the pvalues are saved in an external file
        :type write: bool
        """
        plam_mat = self.get_plambda_matrix()
        nlam_mat = self.get_lambda_motif_matrix(bip_set)
        p_mat = self.get_proj_pmat(plam_mat, nlam_mat, bip_set)
        logp = np.log(p_mat[np.triu_indices_from(p_mat, k=1)])
        loglike = logp.sum()
        print loglike
        return loglike

    def get_proj_pmat(self, plam_mat, nlam_mat, bip_set=False):
        """Return a matrix of probabilities of the observed Lambda motifs in
        the input matrix of the given bipartite node set.

        If the node set is the constrained one, the Lambda motifs follow a
        Binomial probability distribution.
        Otherwise, all the node pairs follow the same Poisson Binomial
        probability distribution.

        NB:
            pmf(k) = Pr(X = k)

        The lower triangular part (including the diagonal) is set to zero.

        :param plam_mat: matrix of Lambda motif probabilities
        :type plam_mat: np.array
        :param nlam_mat: matrix of observed number of Lambda motifs
        :type nlam_mat: np.array
        :param bip_set: selects row-nodes (True) or column-nodes (False)
        :type bip_set: bool
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
        """Obtain and save the p-values of the Lambda motifs observed in the
        binary input matrix for the node set defined by bip_set.

        :param bip_set: analyze row-nodes (True) or column-nodes (False)
        :type bip_set: bool
        :param write: if True, the pvalues are saved in an external file
        :type write: bool
        """
        plam_mat = self.get_plambda_matrix()
        nlam_mat = self.get_lambda_motif_matrix(bip_set)
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
                elif not self.const_set:
                    constr = "columns"
                fname = 'bipcm_pval_' + 'constr_' + constr + 'proj_' +  b  + '.csv'
            else:
                fname = filename
            self.save_matrix(pval_mat, filename=fname, delim=delim)
        else:
            return pval_mat

    def get_plambda_matrix(self):
        """Return a square matrix of Lambda probabilities for the nodes given
        the degree constraints on the node set self.const_set.

        NB:
            The lower triangular part excludint the diagonal is set to 0 since
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

    def get_lambda_motif_matrix(self, bip_set=False):
        """Return a matrix of observed Lambda motifs, i.e.

            A_{ij} = N(Lambda_{ij}).

        :param bip_set: selects row-nodes (True, fix rows) or
                        column-nodes (False, fix columns)
        :type bip_set: bool
        """
        if bip_set:
            l2_mat = np.dot(self.bin_mat, np.transpose(self.bin_mat))
            assert l2_mat.shape == (self.num_rows, self.num_rows)
        elif not bip_set:
            l2_mat = np.dot(np.transpose(self.bin_mat), self.bin_mat)
            assert l2_mat.shape == (self.num_columns, self.num_columns)
        else:
            errmsg = str(bip_set) + 'not supported.'
            raise NameError(errmsg)
        # set diagonal to zero
        di = np.diag_indices_from(l2_mat)
        l2_mat[di] = 0
        return l2_mat.astype(int)

    def get_lambda_pvalues(self, plam_mat, nlam_mat, bip_set=False):
        """Return a matrix of p-values for the observed Lambda motifs in the
        input matrix of the given bipartite node set.

        If the node set is the constrained one, the Lambda motifs follow a
        Binomial probability distribution.
        Otherwise, all the node pairs follow the same Poisson Binomial
        probability distribution.

        NB:
            pval(k) = Pr(X >= k) = 1 - Pr(X < k) = 1 - cdf(k) + pmf(k)

        The lower triangular part (including the diagonal) is set to zero.

        :param plam_mat: matrix of Lambda motif probabilities
        :type plam_mat: np.array
        :param nlam_mat: matrix of observed number of Lambda motifs
        :type nlam_mat: np.array
        :param bip_set: selects row-nodes (True) or column-nodes (False)
        :type bip_set: bool
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
# Auxiliary methods
# ------------------------------------------------------------------------------

    @staticmethod
    def get_main_dir(main_dir_name='bipcm'):
        """Return the absolute path to the main directory which contains the
        folders "src" and "output".
        Note that the default directory name is "bipcm".

        :param main_dir_name: name of the main directory of the program.
        :type main_dir_name: string
        """
        s = os.getcwd()
        dirpath = s[:s.index(main_dir_name) + len(main_dir_name) + 1]
        return dirpath

    def save_matrix(self, mat, filename, delim='\t'):
        """Save the input matrix in a csv-file.

        :param mat: two-dimensional matrix
        :param filename: name of the output file
        :param delim: delimiter between values in file.
        """
        fname = ''.join([self.main_dir, '/output/', filename])
        np.savetxt(fname, mat, delimiter=delim)

################################################################################
# Main
################################################################################

if __name__ == "__main__":
    pass
