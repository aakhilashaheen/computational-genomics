import sys
# function to calculate the difference between two input id_sequences
# simply counts mismatches and return the percentage as a float
def distance(seq1,seq2):
    total = len(seq1)
    mismatch = 0
    for i in xrange(total):
        if seq1[i] != seq2[i]:
            mismatch += 1
    if mismatch == 0:
        return 0
    return mismatch / float(total)


# Fucntion to get the overall distance matrix for all the sequences.
def get_distance_matrix(ids, sequences):
    total = len(ids)
    matrix = [[0] * total for _ in range(total)]
    for i, id1 in enumerate(ids):
        for j, id2 in enumerate(ids):
            matrix[i][j] = distance(sequences[id1], sequences[id2])
    return matrix


# writes the distance matrix to file
def write_distance_matrix(ids, matrix):
    num_ids = len(ids)
    with open('genetic_distances.txt', 'w') as f:
    	# First append the row of ids.
        f.write('\t' + '\t'.join(ids) + '\n')
        for i in range(num_ids):
            f.write(ids[i] + '\t' + '\t'.join(map(str, matrix[i])) + '\n')



# Function to read the fast files with sequences. This returns the sequence,id map
# and the id list for downstream use
def read_fasta_file(filename):
    sequences = {}
    ids = []
    with open(filename) as f:
        content = f.read().splitlines()
    currId = None
    for index, line in enumerate(content):
        if index % 2 == 0:
            currId = line[1:]
            ids.append(currId)
        else:
            sequences[currId] = line
    return ids, sequences





def main(filename):
    # read the fna file and return the ids and id sequence mappings
    ids, sequences = read_fasta_file(filename)

    # generate the difference matrix
    distances = get_distance_matrix(ids, sequences)

    # write the distances.txt file (distance matrix)
    write_distance_matrix(ids, distances)

if __name__ == '__main__':
	main(sys.argv[1])