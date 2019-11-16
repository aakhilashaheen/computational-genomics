import numpy as np
import utils
import sys
import collections
import matplotlib.pyplot as plt

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
            all_variabilites.append(variability_per_region[conserved] / float(total_count))
        else:
            all_variabilites.append(0)
    return all_variabilites


def plot_variabilities(variabilities):
    x = np.arange(1, len(variabilities)+1)
    y = np.array(variabilities)
    # plt.plot(x,y)
    # plt.show()


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def plot_smooth_variabilites(variabilities):
    y = np.array(variabilities)
    y = moving_average(np.array(y), n=20)
    x = np.arange(1, len(variabilities)+1)
    plt.plot(x,y)
    plt.show()    

def main(filename):
    alignments = utils.read_fasta_file(filename)
    variabilities = get_variabilities(alignments)
    utils.write_variabilities(variabilities)
    print(len(variabilities))
    plot_variabilities(variabilities)
    plot_smooth_variabilites(variabilities)

if __name__ == '__main__':
    main(sys.argv[1])
