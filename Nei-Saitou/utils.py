
# Function to calculate the distance between two sequences as a helper 
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
def get_distMatrix(ids, sequences):
    total = len(ids)
    matrix = [[0] * total for _ in range(total)]
    for i, id1 in enumerate(ids):
        for j, id2 in enumerate(ids):
            matrix[i][j] = distance(sequences[id1], sequences[id2])
    return matrix


# writes the distance matrix to a file
def write_distMatrix(ids, matrix):
    num_ids = len(ids)
    with open('genetic_distances.txt', 'w') as f:
    	# First append the row of ids.
        f.write('\t' + '\t'.join(ids) + '\n')
        for i in range(num_ids):
            f.write(ids[i] + '\t' + '\t'.join(map(str, matrix[i])) + '\n')


# Calculates the Q matrix from the distance matrix and returns the two closest tips and their distance
def construct_Q_matrix(distMatrix):
    N = len(distMatrix)
    # Matrix to store the Q values
    Q_matrix = [[0] * N for _ in range(N)]

    # We also want to return the two closest tips
    mini = minj = 0
    minVal = float('inf')

    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            Q_matrix[i][j] = (N - 2) * distMatrix[i][j] - sum(distMatrix[i]) - sum(distMatrix[j])
            # keep track of the minimum value and indices in the Q matrix
            if Q_matrix[i][j] < minVal:
                mini = i
                minj = j
                minVal = Q_matrix[i][j]
    return mini, minj, minVal


# Calculates the edge lengths to the u node and returns the lengths.
def calculate_edge_lengths(distMatrix, mini, minj):
    N = len(distMatrix)
    # Distances to the new internal node
    edge_i = (0.5) * distMatrix[mini][minj] + (0.5 / (N - 2))* (sum(distMatrix[mini]) - sum(distMatrix[minj]))
    # distance from minj node to u
    edge_j = distMatrix[mini][minj] - edge_i
    return edge_i, edge_j


def calculate_new_distMatrix(distMatrix, mini, minj):
    N = len(distMatrix)
    copyMatrix = [[0] * (N + 1) for _ in range(N + 1)]

    # copy everything over
    for i in range(N):
        for j in range(N):
            copyMatrix[i][j] = distMatrix[i][j]
    # update the distances to the new node
    for k in range(N):
        copyMatrix[N][k] = (0.5) * (distMatrix[mini][k] + distMatrix[minj][k] - distMatrix[mini][minj])
        copyMatrix[k][N] = copyMatrix[N][k]

    # shrink the tmp matrix to create the new distance matrix
    new_distMatrix = [[0] * (N - 1) for _ in range(N - 1)]

    # where to move values from tmp values to the new_distMatrix
    keep_i = keep_j = 0
    for i in range(N + 1):
        # Replacing these two with a new node
        if i == mini or i == minj:
            continue
        keep_j = 0
        for j in range(N + 1):
            # Replacing these two with the new node
            if j == mini or j == minj:
                continue
            new_distMatrix[keep_i][keep_j] = copyMatrix[i][j]
            keep_j += 1
        keep_i += 1

    return new_distMatrix