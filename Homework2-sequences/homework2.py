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
    rows = len(ref)
    columns = len(query)
    matrix = np.zeros((rows+1, columns+1))
    MISMATCH_PENALITY = -3
    MATCH_BOOST = 1
    GAP_PENALITY = -2
    matrix[0,0] = 0
    for i in range(1, columns+1):
        matrix[0,i] = matrix[0,i-1] + GAP_PENALITY
    for i in range(1, rows+1):
        matrix[i,0] = matrix[i-1, 0] + GAP_PENALITY
    print(repr(matrix))
    for i in range (1, rows+1):
        for j in range(1, columns+1):
            diag = matrix[i-1, j-1]+MATCH_BOOST if ref[i-1] == query[j-1] else  matrix[i-1, j-1]+MISMATCH_PENALITY
            top = matrix[i-1, j]+ GAP_PENALITY
            left = matrix[i, j-1]+ GAP_PENALITY
            matrix[i, j] = max([diag, top, left])

    
    #To get the alignments
    global_score = matrix[rows, columns]
    i = rows
    j = columns
    alignment1 = ""
    alignment2 = ""
    
    while i > 0 and j > 0:
        diag = matrix[i-1, j-1]
        up = matrix[i-1, j]
        left = matrix[i, j-1]
        maximum = max([diag, up, left])

        if maximum == diag:
            alignment1 =  ref[i-1] + alignment1
            alignment2 =  query[j-1] + alignment2
            i -=1
            j -=1

        elif maximum == left:
            alignment1 = '-' + alignment1
            alignment2 = query[j-1] + alignment2
            j -=1

        else:
            alignment1 = ref[i-1] + alignment1
            alignment2 = '-' + alignment2
            i -=1
    
    while j > 0:
        alignment1 = '-' + alignment1
        alignment2 = query[j-1] + alignment2
        j -=1

    while i > 0:
        alignment1 = ref[i-1] + alignment1
        alignment2 = '-' + alignment2
        i -=1

    print(alignment1)
    print(alignment2)
    print(global_score)
    



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
    needleman_wunsch("AAAA", "AAAA")

__main__()
