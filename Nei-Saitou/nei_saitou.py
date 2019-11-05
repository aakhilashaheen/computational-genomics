
import node 
def nei_saitou(ids, distMatrix, nodeIdCounter, child_parent, id_index_map):	
    # nodeId starts at 120 and decrements
    # keeps track of child parent relationships for constructing the tree later
    # keeps track of id -> index mappings for O(1) lookup

    #Starting the index for internal nodes
    if len(distMatrix) <= 2:
    # add the relationship into our dictionary and also set the root
        if ids[0] in id_index_map:
            node_1 = id_index_map[ids[0]] 
        else:
        	node_1 = ids[0]
       	if ids[1] in id_index_map:
        	node_2 = id_index_map[ids[1]]  
        else: 
        	node_2 = ids[1]
        child_parent[node_1] = (node_2, distance_matrix[0][1])
        root = node_2
        return root
    # Find the Q matrix
    min_i, min_j, min_val = QMatrixCalc(ids, distMatrix)
    #Calculate the edge distance to tips
    d_i, d_j = edgeLengths(distMatrix, min_i, min_j)
    #Create new nodes and join them to an internal node
    # store child parent relations
    # we want to get the node id. ranging from 1-120
    child_1 = id_index_map[ids[min_i]] if ids[min_i] in id_index_map else ids[min_i]
    child_2 = id_index_map[ids[min_j]] if ids[min_j] in id_index_map else ids[min_j]
    child_parent[child_1] = (str(nodeIdCounter), q_i_u)
    child_parent[child_2] = (str(nodeIdCounter), q_j_u)
    # remove the used ids from the id list and also add the new nodeid
    ids = [e for e in ids if e not in (ids[min_i], ids[min_j])] + [str(nodeIdCounter)]
    # create a tmp distance matrix with new node at the end (index N)
    # used to calculate distances to the new node and then used to create the new smaller
    # distance matrix
    tmp_matrix = [[0] * (N + 1) for _ in range(N + 1)]
    # copy everything over
    for i in range(N):
        for j in xrange(N):
            tmp_matrix[i][j] = distance_matrix[i][j]
    # update the distances to the new node
    for k in range(N):
        # based on equations in the book
        tmp_matrix[N][k] = 1 / 2.0 * (distance_matrix[min_i][k] + distance_matrix[min_j][k] - distance_matrix[min_i][min_j])
        tmp_matrix[k][N] = tmp_matrix[N][k]

    # shrink the tmp matrix to create the new distance matrix
    new_distance_matrix = [[0] * (N - 1) for _ in range(N - 1)]

    # where to move values from tmp values to the new_distance_matrix
    i_pointer = j_pointer = 0
    for i in range(N + 1):
        # if it is the minimum indices, skip to the next iteration
        if i == min_i or i == min_j:
            continue
        j_pointer = 0
        for j in range(N + 1):
            if j == min_i or j == min_j:
                continue
            # set the corresponding values in the new distance matrix
            new_distance_matrix[i_pointer][j_pointer] = tmp_matrix[i][j]
            j_pointer += 1
        i_pointer += 1

    # recursively call with the new distance matrix
    nei_saitou(ids, new_distance_matrix, nodeIdCounter-1)
# Function that takes in the distance amtrix and outputs the Q matrix according 
# to the article in Wikipedia on neighbour joining
def QMatrixCalc(distMatrix):
    N = len(distMatrix)
    QMatrix = [[0] * N for _ in range(N)]
    min_i = 0 
    min_j = 0 
    minVal = 0.0
    for i in range(N):
        for j in range(N):
            # We do not have to calculate these
            if i == j:
                continue
            QMatrix[i][j] = (N - 2) * distMatrix[i][j] - sum(distMatrix[i]) - sum(distMatrix[j])
                # Get the minimum Q value in the matrix
            if QMatrix[i][j] < minVal:
            min_i = i
            min_j = j
            minVal = Q_matrix[i][j]
    return min_i, min_j, minVal


def edgeLengths(distMatrix, min_i, min_j):
    d_i = 0.5 * distMatrix[min_i][min_j] + 0.5/(2*(n-2))(sum(distMatrix[i]) - sum(distMatrix[j]))
    d_j = distMatrix[min_i][min_j] - d_i
    return d_i, d_j

def build_phylogeny(ids, distances):
    child_parent = {}
    id_index_map = {}
    for i, id in enumerate(ids):
        id_index_map[id] = str(i + 1)
    root = nei_saitou(ids, distances, 120, child_parent, id_index_map)
    # construct tree from the child_parent relations
    rootNode = TreeNode(root)
    # queue for breadth first search
    queue = [rootNode]

    # bfs to construct the tree
    # will end when all nodes have been traversed
    while len(queue) > 0:
        tmp = []
        for node in queue:
            # for every parent node taht matches the current node, create
            # a new child node and add it to the queue and add it as a child
            for child, (parent, distance) in child_parent.items():
                if parent == node.id:
                    childNode = TreeNode(child)
                    node.add_child(childNode, distance)
                    tmp.append(childNode)
        # delete the node in the dictionary so we dont reuse it
        for node in tmp:
            del child_parent[node.id]
        # the queue is set to the next level
        queue = tmp

    # return the root of the tree
    return rootNode