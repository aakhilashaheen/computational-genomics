import numpy as np
import utils
import sys
import collections
import matplotlib.pyplot as plt
import random

class Bonus:

  def bonusCalc(self, filename, regions):
    # dictionary = utils.read_fasta_file(filename)
    # # Generating 100 random sequences
    # n_keys = dictionary.keys()
    # rand_100 = random.sample(n_keys, 100)
    # rand_seqs = {}
    # for key in rand_100:
    #     rand_seqs[key] = dictionary[key]

    # #get regions 2 and 4
    # reg2 = regions[1]
    # reg4 = regions[3]
    # reg2_dict = {}
    # reg4_dict = {}
    # for key, value in rand_seqs.items():
    #     reg2_dict[key] = value[reg2[0]:reg2[1]+1]
    # utils.write_fasta_file(reg2_dict, 'FASTA/reg1.fna')
    # for key, value in rand_seqs.items():
    #     reg4_dict[key] = value[reg4[0]:reg4[1]+1]
    # utils.write_fasta_file(reg4_dict, 'FASTA/reg4.fna')
    # utils.write_fasta_file(rand_seqs, 'FASTA/100-random.fna')
    #select 100 sequences
    dictionary = utils.read_fasta_file(filename)
    first100 = {key: dictionary[key] for key in dictionary.keys()[:100]}
    first100_no_gaps = collections.defaultdict(str)
    for i in xrange(1474):
        # skip the gaps
        gap_count = 0
        for key, value in first100.items():
            if value[i] == '-':
                gap_count += 1
        # if its all gaps then ignore it
        if gap_count == 100:
            continue
        for key, value in first100.items():
            first100_no_gaps[key] += value[i]
    #get regions 1 and 4
    reg1 = regions[0]
    reg4 = regions[3]
    reg1_dict = {}
    reg4_dict = {}
    for key, value in first100_no_gaps.items():
        reg1_dict[key] = value[reg1[0]:reg1[1]+1]
    utils.write_fasta_file(reg1_dict, 'FASTA/reg1.fna')
    for key, value in first100_no_gaps.items():
        reg4_dict[key] = value[reg4[0]:reg4[1]+1]
    utils.write_fasta_file(reg4_dict, 'FASTA/reg4.fna')
    utils.write_fasta_file(first100_no_gaps, 'FASTA/100-no-gaps.fna')