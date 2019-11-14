# This file contains the utils that could be helpful in general

# Function that reads the given fasta file and gets the sequences.
def read_fasta_file(filename):
    alignments = {}
    with open(filename) as f:
        lines = f.read().splitlines()
    curr = None
    for index, line in enumerate(lines):
        if index % 2 == 0:
            curr = line[1:]
        else:
            alignments[curr] = line
    return alignments

def write_variabilities(variabilities):
    with open('variability.txt', 'w') as f:
        for variability in variabilities:
            #print(variability)
            f.write(str(variability) + '\n')

