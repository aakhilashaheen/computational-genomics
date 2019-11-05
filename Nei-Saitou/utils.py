def calculateQMatrix(distanceMatrix):
   '''
   Calculates the entire Q matrix.

   Receive: Distance matrix.
   Return: Q matrix.
   '''

   q = numpy.zeros((len(distanceMatrix),len(distanceMatrix)), int)

   for i in range(1, len(distanceMatrix)):
      for j in range(i):
         q[i][j] = calculateQMatrix(distanceMatrix, i, j)

   return q

def minQValue(QMatrix):
	pass


def calculateNewDistanceMatrix(ids, sequences):
	pass


# function to calculate the difference between two input id_sequences
# simply counts mismatches and return the percentage as a float
def distance(seq1,seq2):
    total = len(seq1)
    mismatch = 0
    for i in range(total):
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

