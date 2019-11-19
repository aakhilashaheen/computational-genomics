# This file contains the utils that could be helpful in general

# Function that reads the given fasta file and gets the sequences.
def read_fasta_file(filename):
    sequences = {}
    with open(filename) as f:
        lines = f.read().splitlines()
    curr = None
    for index, line in enumerate(lines):
        if index % 2 == 0:
            curr = line[1:]
        else:
            sequences[curr] = line
    return sequences

# Function that writes the identities to a file.
def write_identities(conserved):
    with open('identities.txt', 'w') as f:
        for c in conserved:
            #print(variability)
            f.write(str(c) + '\n')

'''
write fasta file (bonus) given a dictionary of id->sequence
'''
def write_fasta_file(dictionary, filename):
    with open(filename, 'w') as f:
        for key, value in dictionary.items():
            f.write('>' + key + '\n')
            f.write(value + '\n')
