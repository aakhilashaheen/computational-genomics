import numpy as np
import utils
import sys
import collections


'''
calculate the variability
count up the # of each base there is (excluding gaps)
find the base that has the most count and divide by the total count for that position
'''
def get_variabilities(dictionary):
    all_variabilites = []
    for i in range(1474):
        variability_per_region = {}
        total_count = 0
        for key, value in dictionary.items():
            #print(len(value))
            if value[i] == '-':
                continue
            if value[i] in variability_per_region:
                variability_per_region[value[i]] += 1
            else:
                variability_per_region[value[i]] = 1
            total_count+=1
        if len(variability_per_region) > 0:
            # get the maximum count and divide by total count
            conserved = max(variability_per_region, key=variability_per_region.get)
            all_variabilities.append(variability_per_region[conserved] / float(total_count))
        else:
            all_variabilities.append(0)
    return all_variabilities

def main(filename):
    alignments = utils.read_fasta_file(filename)
    variabilities = get_variabilities(alignments)
    utils.write_variabilities(variabilities)
    print(len(variabilities))

if __name__ == '__main__':
    main(sys.argv[1])
