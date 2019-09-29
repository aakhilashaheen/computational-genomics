#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
#from subprocess import Popen, PIPE

#def make_arg_parser():
#    parser = argparse.ArgumentParser(prog='template.py',
#                          version="%prog 1.0",
#                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument("-q","--query",
#                      default=argparse.SUPPRESS,
#                      required=True,
#                      help="Path to query fasta [required]") 
#    parser.add_argument("-r","--ref",
#                      default=argparse.SUPPRESS,
#                      required=True,
#                      help="Path to reference fasta [required]")
#    parser.add_argument("-t","--taxonomy",
#                      default=None,
#                      required=True,
#                      help="Path to taxonomy file [required]")
#    parser.add_argument("-o","--output",
#                      default=None,
#                      required=True,
#                      help="Path to output file [required]")
#    parser.add_argument("-c","--command",
#                      default='./burst',
#                      help="Path to BURST command") 
#    parser.add_argument("-V","--verbose",
#                      action="store_true",
#                      help="Verbose output")
#    return parser

# Runs plain old needleman wunsch algorithm
def needleman_wunsch(query, ref):
    n = len(ref)
    m = len(query)
    matrix = np.zeros((n+1, m+1))
    MISMATCH_PENALITY = -3
    MATCH_BOOST = 1
    GAP_PENALITY = -2
    matrix[0][0] = 0
    for i in range(1, m+1):
        matrix[0][i] = matrix[0][i-1] + GAP_PENALITY
    for i in range(1, n+1):
        matrix[i][0] = matrix[i-1][0] + GAP_PENALITY
    for i in range (1, n+1):
        for j in range(1, m+1):
            diag = matrix[i-1][j-1]+MATCH_BOOST if ref[i-1] == query[j-1] else  matrix[i-1][j-1]+MISMATCH_PENALITY
            top = matrix[i-1][j]+ GAP_PENALITY
            left = matrix[i][j-1]+ GAP_PENALITY
            matrix[i][j] = max([diag, top, left])



    print(repr(matrix))
            



#if __name__ == '__main__':
#    parser = make_arg_parser()
#    args = parser.parse_args()
#    return_value, stdout, stderr = run_burst(args.query, args.ref, args.taxonomy, args.output)
#    print(return_value)
#    print(stdout)
#    print(stderr)
    ### YOUR CODE HERE ###


def __main__():
    needleman_wunsch("ACAA", "ACTGA")

__main__()
