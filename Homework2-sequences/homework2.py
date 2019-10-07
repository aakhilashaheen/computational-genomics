#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
import matplotlib.pyplot as plt
#from subprocess import Popen, PIPE
import random

def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework2.py',
                          version="%prog 1.0",
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]")
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-m","--matches",
                      default=None,
                      required=False,
                      help="matches.txt file [optional]")
    return parser

#This implements the basic needleman wunsch algorithm approach with dynamic programming
def needleman_wunsch(query, ref):
    rows = len(ref)
    columns = len(query)
    #Append an extra row and column for the algorithm 
    matrix = np.zeros((rows+1, columns+1))
    MISMATCH_PENALITY = -3
    MATCH_BOOST = 1
    GAP_PENALITY = -2
    #The first index is set as 0
    matrix[0,0] = 0

    #Filling up the initial starting values
    for i in range(1, columns+1):
        matrix[0,i] = matrix[0,i-1] + GAP_PENALITY
    for i in range(1, rows+1):
        matrix[i,0] = matrix[i-1, 0] + GAP_PENALITY
    
    #Starting the algorithm for other cells in the matrix
    for i in range (1, rows+1):
        for j in range(1, columns+1):
            diag = matrix[i-1, j-1]+MATCH_BOOST if ref[i-1] == query[j-1] else  matrix[i-1, j-1]+MISMATCH_PENALITY
            top = matrix[i-1, j]+ GAP_PENALITY
            left = matrix[i, j-1]+ GAP_PENALITY
            matrix[i, j] = max([diag, top, left])
        
    
    #To get the alignments
    score = matrix[rows, columns]
    i = rows
    j = columns
    alignment1 = ""
    alignment2 = ""
    
    while i > 0 and j > 0:
        diag = matrix[i-1, j-1]
        up = matrix[i-1, j]
        left = matrix[i, j-1]

        if diag + MATCH_BOOST == matrix[i][j] or diag + MISMATCH_PENALITY == matrix[i][j]:
            alignment1 =  ref[i-1] + alignment1
            alignment2 =  query[j-1] + alignment2
            i -=1
            j -=1

        elif left + GAP_PENALITY == matrix[i][j]:
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
        
    #print("Reference:")
    #print(ref)
    #print("Query:")
    #print(query)
    #print("Alignment1 :")
    #print(alignment1)
    #print("Alignment2:")
    #print(alignment2)
    #print("Score:")
    #print(score)
    
    return score, alignment1, alignment2

#Anchored needleman wunsch algorithm wrapper
def anchored_needleman_wunsch(query, ref, matches):
    total_score = 0
    #First match index
    start1, end1, start2, end2 = matches[0]
    temp_query = query[0:start1-1]
    #print("NW")
    #print(temp_query)
    temp_ref = ref[0:start2-1]
    #print(temp_ref)
    score, alignment1, alignment2 = needleman_wunsch(temp_query, temp_ref)
    total_score = score
    align1 = alignment1
    align2 = alignment2
    prevEnd1 = 0
    prevEnd2 = 0
    findPrevNW = 0
    for key in matches:
        #Accounting for sequences before the first match 
        if prevEnd1 is not 0 and prevEnd2 is not 0:
            findPrevNW = 1
        start1, end1, start2, end2 = matches[key]
        temp_query = query[start1-1:end1]
        align1 = align1 + temp_query
        #print("Matches")
        #print(temp_query)
        temp_ref = ref[start2-1: end2]
        align2 = align2 + temp_ref
        #print(temp_ref)
        score = sum(temp_query[i] == temp_ref[i] for i in range(len(temp_query)))
        score = score + sum(temp_query[i] != temp_ref[i] for i in range(len(temp_query))) * -3
        total_score = total_score + score 
        if findPrevNW == 1:
            temp_query = query[prevEnd1 : end1]
            #print("NW")
            #print(temp_query)
            temp_ref = ref[prevEnd2 : end2]
            #print(temp_ref)
            score, alignment1, alignment2 = needleman_wunsch(temp_query, temp_ref)
            align1 = align1 + alignment1
            align2 = align2 + alignment2
            total_score = total_score + score
        prevEnd1 = end1
        prevEnd2 = end2
    
    #Accounting for sequences after the last match
    temp_query = query[prevEnd1:]
    temp_ref = ref[prevEnd2:]
    #print("NW")
    #print(temp_query)
    #print(temp_ref)
    score, alignment1, alignment2 = needleman_wunsch(temp_query, temp_ref)
    total_score = total_score + score
    align1 = align1 + alignment1
    align2 = align2 + alignment2
    print("Total Score:")
    print(total_score)
    print("Alignment 1 :")
    print(align1)
    print("Alignment 2 :")
    print(align2)

    return total_score

#Reads the ref and query fasta files
def readFasta(fileName):
    with open(fileName) as f:
        sequence = f.read().splitlines()
        #print(sequence)
    return "".join(sequence[1:])

#Reads the matches fasta files
def readMatchesFasta(fileName):
    indices = {}
    with open(fileName) as f:
        content = f.readlines()
    i = 0
    for line in content:
        values = []
        values = line.split("\t")
        values[3] = values[3].split("\r")[0]
        for j in range(len(values)):
            values[j] = int(values[j])
        indices[i] = values
        i +=1
    return indices   

#Runs the standard needleman wunsch with shuffle 10000 times
def run_10000_perm(query, ref, original_score):
    scores =[]
    for i in range(10000):
        query_temp = list(query)
        ref_temp = list(query)
        random.shuffle(query_temp)
        query_temp = ''.join(query_temp)
        random.shuffle(ref_temp)
        ref_temp = ''.join(ref_temp)
        score, alignment1, alignment2 = needleman_wunsch(query_temp, ref_temp)
        scores.append(score)
    scores_np = np.array(scores)
    #Plots the required histograms
    plt.hist(scores_np, bins='auto')
    plt.axvline(x=original_score, color='r', linestyle='dashed', linewidth=1)
    plt.show()
    return

def __main__():
    parser = make_arg_parser()
    args = parser.parse_args()
    query = readFasta(args.query)
    ref = readFasta(args.ref)
    if not args.matches:
        original_score, alignment1, alignment2 = needleman_wunsch(query=query,ref=ref)
        #print(original_score)
        #print(alignment1)
        #print(alignment2)
        run_10000_perm(query, ref, original_score)
    else:
        matches = readMatchesFasta(args.matches)
        anchored_needleman_wunsch(query=query, ref=ref, matches=matches)
__main__()
