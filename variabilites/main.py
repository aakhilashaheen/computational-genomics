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
        for key, value in dictionary.items():
            #print(key, value[i])
            if value[i] == '-':
                continue
            if value[i] in variability_per_region:
                variability_per_region[value[i]] += 1
            else:
                variability_per_region[value[i]] = 1
        
        if len(variability_per_region) > 0:
            total_count = 0
            for count in variability_per_region.values():
                total_count += count
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


def running_mean(x, N=10):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def plot_smooth_variabilites(variabilities):
    y = np.array(running_mean(variabilities, N=40))
    x = np.arange(1, len(y)+1)
    plt.plot(x,y)
    plt.show()    

def get_variable_regions(variabilities):
    # some example data
    threshold = 0.75
    values = np.array(variabilities)
    x = range(len(values))

    # split it up
    above_threshold = np.maximum(values - threshold, 0)
    below_threshold = np.minimum(values, threshold)

    # and plot it
    fig, ax = plt.subplots()
    ax.bar(x, below_threshold, 0.35, color="g")
    ax.bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)

    # horizontal line indicating the threshold
    ax.plot([0., 9.5], [threshold, threshold], "k--")

    fig.savefig("look-ma_a-threshold-plot.png")       


def main(filename):
    alignments = utils.read_fasta_file(filename)
    variabilities = get_variabilities(alignments)
    utils.write_variabilities(variabilities)
    #plot_variabilities(variabilities)
    #plot_smooth_variabilites(variabilities)
    get_variable_regions(variabilities)

if __name__ == '__main__':
    main(sys.argv[1])
