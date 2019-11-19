import numpy as np
import utils
import sys
import collections
import matplotlib.pyplot as plt
import random

class Bonus:

  def bonusCalc(self, filename, regions):
    dictionary = utils.read_fasta_file(filename)
    # Generating 100 random sequences
    n_keys = dictionary.keys()
    rand_100 = random.sample(n_keys, 100)
    rand_seqs = {}
    whole_seqs = {}

    for key in rand_100:
        rand_seqs[key] = dictionary[key]
    
    # Fill whole sequences without gaps
    for i in range(1474):
        gaps = 0
        for key, val in rand_seqs.items():
            if val[i] == '-':
                gaps += 1
        if gaps == 100:
            continue
        for key, val in rand_seqs.items():
            if key in whole_seqs.keys():
                whole_seqs[key] += val[i]
            else:
                whole_seqs[key] = val[i]
    utils.write_fasta_file(whole_seqs, 'whole_seqs.fna')

    #get regions 2 and 4
    reg1 = regions[0]
    reg4 = regions[3]
    reg1_dict = {}
    reg4_dict = {}
    for key, value in rand_seqs.items():
        reg1_dict[key] = value[reg1[0]:reg1[1]+1]
    utils.write_fasta_file(reg1_dict, 'reg1.fna')
    for key, value in rand_seqs.items():
        reg4_dict[key] = value[reg4[0]:reg4[1]+1]
    utils.write_fasta_file(reg4_dict, 'reg4.fna')


