import numpy as np
import utils
import sys
import collections
import matplotlib.pyplot as plt
import bonus
'''
calculate the identities
count up the # of each base there is (excluding gaps)
find the base that has the most count and divide by the total count for that position
'''
def get_conserved(dictionary):
    all_conserved = []

    for i in range(1474):
        total_count = len(dictionary.keys())
        conserved_per_region = {}
        for key, value in dictionary.items():
            if value[i] == '-':
                continue
            if value[i] in conserved_per_region:
                conserved_per_region[value[i]] += 1
            else:
                conserved_per_region[value[i]] = 1

        if len(conserved_per_region) > 0:
            # get the maximum count and divide by total count
            conserved = max(conserved_per_region, key=conserved_per_region.get)
            all_conserved.append(conserved_per_region[conserved] / float(total_count))
        else:
            all_conserved.append(0)
    return all_conserved

# Calculates the running mean for smoothing
def running_mean(variabilies, N):
    sum = 0
    result = list( 0 for x in variabilies)
 
    for i in range( 0, N ):
        sum = sum + variabilies[i]
        result[i] = sum / (i+1)
 
    for i in range( N, len(variabilies) ):
        sum = sum - variabilies[i-N] + variabilies[i]
        result[i] = sum / N
 
    return result

# Gets the smoothing for variabilites and plots it
def plot_smooth_variabilites(variabilities):
    y = np.array(running_mean(variabilities, N=30))
    x = np.arange(1, len(y)+1)
    plt.plot(x,y)
    plt.show()    
    return y

# Identifying the variable regions using a window size and threshold for density
# The density estimation is set using experiments.
def get_variable_regions(variabilities):
    window_size = 40
    threshold = 0.81
    prevInterval = []
    total = len(variabilities)
    i = 0
    regions = []
    while i < total-window_size:
        window = variabilities[i: i+window_size]     
        count = 0
        for j in window:
            if j <= threshold:
                count+=1
        if count == 40:
            if len(regions) == 0:
                regions.append([i, i+window_size])
            else:
                lastRegion = regions[len(regions)-1]
                if lastRegion[1] >= i:
                    #Merge Regions
                    popped = regions.pop()
                    regions.append([popped[0], i+window_size])
                else:
                    regions.append([i,i+window_size])
        i+=1

    return regions

#Plots the variabilites with region markers
def plot_variabilities_with_regions(variabilities):
    y = np.array(running_mean(variabilities, N=30))
    x = np.arange(1, len(y)+1)
    plt.plot(x,y)
    for x1, x2 in get_variable_regions(y):
        plt.plot([x1, x2], [.6, .6], 'r-', lw=2)
    plt.show()    


def main(filename):
    # Reads the fasta file and gets the alignments
    alignments = utils.read_fasta_file(filename)

    # Finds the conserved percentages
    conserved = get_conserved(alignments)

    #Writes the identities to a file
    utils.write_identities(conserved)

    # Plots the smoothened vriabilites
    y = plot_smooth_variabilites(conserved)

    # Gets the regions with high variablilty
    regions = get_variable_regions(y)

    # Writes the regions to a file
    utils.write_regions(regions)

    # Plots the variabilities with region markers
    plot_variabilities_with_regions(conserved)

    #Bonus calculations
    Bonus_instance = bonus.Bonus()
    Bonus_instance.bonusCalc(filename, regions)

if __name__ == '__main__':
    main(sys.argv[1])
